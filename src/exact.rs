#![forbid(unsafe_code)]

//! Exact arithmetic operations via arbitrary-precision rational numbers.
//!
//! This module is only compiled when the `"exact"` Cargo feature is enabled.
//!
//! # Architecture
//!
//! ## Determinants
//!
//! All determinant methods (`det_exact`, `det_exact_f64`,
//! `det_exact_rounded_f64`, and `det_sign_exact`) share the same integer-scaled
//! determinant core. Each f64 entry is decomposed via `f64_decompose` into
//! `mantissa × 2^exponent`, then all entries are scaled to a common `BigInt`
//! matrix (shifting by `e - e_min`). D≤4 uses direct integer expansions; larger
//! matrices use fraction-free Bareiss elimination entirely in `BigInt`
//! arithmetic — no `BigRational`, no GCD, no denominator tracking. The result
//! is `(det_int, total_exp)` where `det = det_int × 2^(D × e_min)`. `det_exact`
//! wraps this with `bigint_exp_to_bigrational` to reconstruct a reduced
//! `BigRational`; `det_exact_f64` converts the same pair only when the exact
//! value is representable as finite binary64; `det_exact_rounded_f64` rounds
//! the same exact value to finite binary64; and `det_sign_exact` reads the sign
//! directly from `det_int` (the scale factor is always positive).
//!
//! `det_sign_exact` adds a two-stage adaptive-precision optimisation inspired
//! by Shewchuk's robust geometric predicates:
//!
//! 1. **Fast filter (D ≤ 4)**: compute `det_direct()` and a conservative error
//!    bound. If `|det| > bound`, the f64 sign is provably correct — return
//!    immediately without allocating.
//! 2. **Exact fallback**: run integer-only Bareiss for a guaranteed-correct
//!    sign.
//!
//! ## Linear system solve
//!
//! `solve_exact`, `solve_exact_f64`, and `solve_exact_rounded_f64` solve
//! `A x = b` with a hybrid algorithm that reuses the integer-only Bareiss core
//! used for determinants.  Matrix and RHS entries are decomposed via
//! `f64_decompose` into `mantissa × 2^exponent`, scaled to a shared
//! base `2^e_min`, and assembled into a `BigInt` augmented system
//! `(A | b)`.  Forward elimination runs entirely in `BigInt` with
//! fraction-free Bareiss updates — no `BigRational`, no GCD
//! normalisation in the `O(D³)` phase.  Once the system is upper
//! triangular, back-substitution is performed in `BigRational`, where
//! fractions are inherent; this phase is only `O(D²)` so the rational
//! overhead is modest.  First-non-zero pivoting is used throughout;
//! since all arithmetic is exact, any non-zero pivot gives the correct
//! result (no numerical stability concern).  Every finite `f64` is
//! exactly representable as a rational, so the result is provably correct.
//! `solve_exact_f64` returns `Vector<D>` only when every exact component is
//! exactly representable as finite binary64; `solve_exact_rounded_f64` returns
//! the exact components rounded to finite binary64.
//!
//! ## f64 → integer decomposition
//!
//! Both the determinant and solve paths share a single conversion
//! primitive, `f64_decompose`, which extracts `(mantissa, exponent,
//! sign)` from the IEEE 754 binary64 bit representation (\[9\]).  The
//! determinant path combines those components into a `BigInt` matrix
//! (for Bareiss) and a `2^(D × e_min)` scale factor, while the solve
//! path builds a `BigInt` augmented system and lifts the
//! upper-triangular result into `BigRational` for back-substitution.
//! See Goldberg \[10\] for background on floating-point representation
//! and exact rational reconstruction.  Reference numbers refer to
//! `REFERENCES.md`.
//!
//! ## Validation
//!
//! Public `Matrix` / `Vector` values are finite by construction before exact
//! methods reach the integer-Bareiss core. The decomposition helpers for those
//! domain types can then call `f64_decompose` without repeating stored-entry
//! validation; `f64_decompose` itself is therefore never called with non-finite
//! input from the public API.

use core::hint::cold_path;
use core::mem::take;
use core::num::NonZeroU64;
use std::array::from_fn;

use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::ToPrimitive;

use crate::matrix::Matrix;
use crate::vector::Vector;
use crate::{LaError, UnrepresentableReason};

const F64_SIGNIFICAND_BITS: i64 = 53;
const F64_FRACTION_BITS: i64 = 52;
const F64_MIN_BINARY_EXPONENT: i64 = -1074;
const F64_MIN_NORMAL_EXPONENT: i64 = -1022;
const F64_MAX_BINARY_EXPONENT: i64 = 1023;
const F64_EXPONENT_BIAS: i64 = 1023;
const F64_FRACTION_MASK: u64 = (1u64 << 52) - 1;

/// Decompose a finite `f64` into its IEEE 754 components.
///
/// Returns `None` for ±0.0, or `Some((mantissa, exponent, is_negative))` with a
/// non-zero mantissa where the value is exactly
/// `(-1)^is_negative × mantissa × 2^exponent` and `mantissa` is odd (trailing
/// zeros stripped).  See `REFERENCES.md` \[9-10\].
///
/// # Errors
/// Returns [`LaError::NonFinite`] if `x` is NaN or infinite.
const fn f64_decompose(x: f64) -> Result<Option<(NonZeroU64, i32, bool)>, LaError> {
    let bits = x.to_bits();
    let biased_exp = ((bits >> 52) & 0x7FF) as i32;
    let fraction = bits & 0x000F_FFFF_FFFF_FFFF;

    // ±0.0
    if biased_exp == 0 && fraction == 0 {
        return Ok(None);
    }

    if biased_exp == 0x7FF {
        cold_path();
        return Err(LaError::non_finite_at(0));
    }

    let (mantissa, raw_exp) = if biased_exp == 0 {
        // Subnormal: (-1)^s × 0.fraction × 2^(-1022)
        //          = (-1)^s × fraction × 2^(-1074)
        (fraction, -1074_i32)
    } else {
        // Normal: (-1)^s × 1.fraction × 2^(biased_exp - 1023)
        //       = (-1)^s × (2^52 | fraction) × 2^(biased_exp - 1075)
        ((1u64 << 52) | fraction, biased_exp - 1075)
    };

    // Strip trailing zeros so the mantissa is odd.
    let tz = mantissa.trailing_zeros();
    let mantissa = mantissa >> tz;
    let Some(mantissa) = NonZeroU64::new(mantissa) else {
        cold_path();
        return Ok(None);
    };
    let exponent = raw_exp + tz.cast_signed();
    let is_negative = bits >> 63 != 0;

    Ok(Some((mantissa, exponent, is_negative)))
}

/// Convert a `BigInt × 2^exp` pair to a reduced `BigRational`.
///
/// When `exp < 0` (denominator is `2^(-exp)`), shared factors of 2 are
/// stripped from `value` to keep the fraction in lowest terms without a
/// full GCD computation.
fn bigint_exp_to_bigrational(mut value: BigInt, mut exp: i32) -> BigRational {
    if value == BigInt::from(0) {
        return BigRational::from_integer(BigInt::from(0));
    }

    // Strip shared powers of 2 between value and the 2^(-exp) denominator.
    if exp < 0
        && let Some(tz) = value.trailing_zeros()
    {
        let exp_abs = exp.unsigned_abs();
        let reduce = tz.min(u64::from(exp_abs));
        value >>= reduce;
        #[allow(clippy::cast_possible_truncation)]
        let reduce = reduce as u32;
        let remaining_abs = exp_abs - reduce;
        exp = match remaining_abs {
            0 => 0,
            2_147_483_648 => i32::MIN,
            value => -value.cast_signed(),
        };
    }

    if exp >= 0 {
        BigRational::new_raw(value << exp.cast_unsigned(), BigInt::from(1u32))
    } else {
        BigRational::new_raw(value, BigInt::from(1u32) << exp.unsigned_abs())
    }
}

/// Convert an exact rational result to `f64` only when the conversion is exact.
///
/// This supports the strict `*_exact_f64` public APIs by accepting only dyadic
/// rational values that fit in finite binary64. The optional `index` is attached
/// to [`LaError::Unrepresentable`] for vector-valued solve components.
///
/// # Errors
/// Returns [`LaError::Unrepresentable`] with
/// [`UnrepresentableReason::RequiresRounding`] when the rational denominator is
/// not a power of two and the rounded value would still be finite.
fn exact_rational_to_finite_f64(exact: &BigRational, index: Option<usize>) -> Result<f64, LaError> {
    let numerator = exact.numer();
    if numerator.sign() == Sign::NoSign {
        return Ok(0.0);
    }

    let denominator = exact.denom();
    let Some(denominator_exp) = denominator.trailing_zeros() else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            rounded_rational_unrepresentable_reason(exact),
        ));
    };

    if denominator.bits().checked_sub(1) != Some(denominator_exp) {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            rounded_rational_unrepresentable_reason(exact),
        ));
    }

    let Ok(denominator_exp) = i32::try_from(denominator_exp) else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            rounded_rational_unrepresentable_reason(exact),
        ));
    };

    bigint_exp_to_finite_f64(numerator.clone(), -denominator_exp, index)
}

/// Classify a failed exact-rational-to-`f64` conversion by the rounded result.
///
/// Strict exact conversion has already failed when this helper is called. It
/// preserves the [`UnrepresentableReason`] recovery contract: callers may retry
/// with a rounded API only when that rounded result would still be finite.
fn rounded_rational_unrepresentable_reason(exact: &BigRational) -> UnrepresentableReason {
    match exact.to_f64() {
        Some(value) if value.is_finite() => UnrepresentableReason::RequiresRounding,
        _ => UnrepresentableReason::NotFinite,
    }
}

/// Convert an exact rational result to a rounded finite `f64`.
fn exact_rational_to_rounded_f64(
    exact: &BigRational,
    index: Option<usize>,
) -> Result<f64, LaError> {
    let Some(value) = exact.to_f64() else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    };
    if value.is_finite() {
        Ok(value)
    } else {
        cold_path();
        Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ))
    }
}

/// Convert a `BigInt × 2^exp` pair to an exactly represented finite `f64`.
///
/// This avoids allocating a [`BigRational`] when determinant and solve paths
/// already have an integer significand plus binary exponent. The optional
/// `index` is forwarded to [`LaError::Unrepresentable`] for vector-valued solve
/// components; determinant callers pass `None`.
///
/// # Errors
/// Returns [`LaError::Unrepresentable`] with
/// [`UnrepresentableReason::RequiresRounding`] when the exact nonzero value
/// would need rounding or underflows below the smallest positive subnormal.
///
/// Returns [`LaError::Unrepresentable`] with
/// [`UnrepresentableReason::NotFinite`] when the exact value cannot be
/// represented by any finite `f64`.
fn bigint_exp_to_finite_f64(
    mut value: BigInt,
    exp: i32,
    index: Option<usize>,
) -> Result<f64, LaError> {
    if value == BigInt::from(0) {
        return Ok(0.0);
    }

    let is_negative = value.sign() == Sign::Minus;
    if is_negative {
        value = -value;
    }

    let mut exp = i64::from(exp);
    if let Some(tz) = value.trailing_zeros() {
        value >>= tz;
        let Ok(tz) = i64::try_from(tz) else {
            cold_path();
            return Err(LaError::unrepresentable(
                index,
                UnrepresentableReason::NotFinite,
            ));
        };
        let Some(updated_exp) = exp.checked_add(tz) else {
            cold_path();
            return Err(LaError::unrepresentable(
                index,
                UnrepresentableReason::NotFinite,
            ));
        };
        exp = updated_exp;
    }

    let bit_len = value.bits();
    let Ok(bit_len) = i64::try_from(bit_len) else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    };
    let Some(top_bit_exp) = exp.checked_add(bit_len - 1) else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    };
    if exp < F64_MIN_BINARY_EXPONENT {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::RequiresRounding,
        ));
    }
    if top_bit_exp > F64_MAX_BINARY_EXPONENT {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    }
    if bit_len > F64_SIGNIFICAND_BITS {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::RequiresRounding,
        ));
    }

    let Some(mantissa) = value.to_u64() else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    };
    let sign = if is_negative { 1u64 << 63 } else { 0 };

    if top_bit_exp < F64_MIN_NORMAL_EXPONENT {
        let Ok(shift) = u32::try_from(exp - F64_MIN_BINARY_EXPONENT) else {
            cold_path();
            return Err(LaError::unrepresentable(
                index,
                UnrepresentableReason::RequiresRounding,
            ));
        };
        Ok(f64::from_bits(sign | (mantissa << shift)))
    } else {
        let Ok(biased_exp) = u64::try_from(top_bit_exp + F64_EXPONENT_BIAS) else {
            cold_path();
            return Err(LaError::unrepresentable(
                index,
                UnrepresentableReason::NotFinite,
            ));
        };
        let Ok(shift) = u32::try_from(F64_FRACTION_BITS - (bit_len - 1)) else {
            cold_path();
            return Err(LaError::unrepresentable(
                index,
                UnrepresentableReason::RequiresRounding,
            ));
        };
        let significand = mantissa << shift;
        Ok(f64::from_bits(
            sign | (biased_exp << F64_FRACTION_BITS) | (significand & F64_FRACTION_MASK),
        ))
    }
}

/// Convert a `BigInt × 2^exp` determinant pair to a rounded finite `f64`.
fn bigint_exp_to_rounded_f64(value: BigInt, exp: i32) -> Result<f64, LaError> {
    let exact = bigint_exp_to_bigrational(value, exp);
    exact_rational_to_rounded_f64(&exact, None)
}

// -----------------------------------------------------------------------
// Shared integer-Bareiss primitives
// -----------------------------------------------------------------------
//
// Both `bareiss_det_int` (determinants) and `gauss_solve` (linear system
// solve) follow the same pipeline: decompose every f64 entry into
// `(mantissa, exponent, is_negative)`, track the minimum exponent across
// non-zero entries, scale each entry by `2^(exp − e_min)` to build a
// fully-integer `BigInt` matrix, and run Bareiss fraction-free forward
// elimination.  The helpers below factor out each stage so the two
// callers differ only in post-processing (± sign for det, back-sub for
// solve) and in whether they carry a RHS through the elimination.

/// Decomposed finite f64 in the form `(-1)^is_negative · mantissa · 2^exponent`.
///
/// `Zero` represents ±0.0. Non-zero entries carry a [`NonZeroU64`] mantissa, so
/// the exact-arithmetic paths cannot accidentally combine an absent mantissa
/// with active exponent/sign fields after decomposition.
#[derive(Clone, Copy, Default)]
enum Component {
    #[default]
    Zero,
    NonZero {
        mantissa: NonZeroU64,
        exponent: i32,
        is_negative: bool,
    },
}

/// Decompose every entry of a finite `D×D` matrix via `f64_decompose`.
///
/// Returns the per-entry components and the minimum exponent across non-zero
/// entries. If every entry is zero, the exponent is `i32::MAX`.
fn decompose_finite_matrix<const D: usize>(
    m: &Matrix<D>,
) -> Result<([[Component; D]; D], i32), LaError> {
    let mut components = [[Component::default(); D]; D];
    let mut e_min = i32::MAX;
    for (r, row) in m.rows().iter().enumerate() {
        for (c, &entry) in row.iter().enumerate() {
            if let Some((mantissa, exponent, is_negative)) =
                f64_decompose(entry).map_err(|_| LaError::non_finite_cell(r, c))?
            {
                components[r][c] = Component::NonZero {
                    mantissa,
                    exponent,
                    is_negative,
                };
                e_min = e_min.min(exponent);
            }
        }
    }
    Ok((components, e_min))
}

/// Decompose every entry of a finite length-`D` vector via `f64_decompose`.
///
/// Returns the per-entry components and the minimum exponent across non-zero
/// entries. If every entry is zero, the exponent is `i32::MAX`.
fn decompose_finite_vec<const D: usize>(v: &Vector<D>) -> Result<([Component; D], i32), LaError> {
    let mut components = [Component::default(); D];
    let mut e_min = i32::MAX;
    let data = v.as_array();
    for (i, &entry) in data.iter().enumerate() {
        if let Some((mantissa, exponent, is_negative)) =
            f64_decompose(entry).map_err(|_| LaError::non_finite_at(i))?
        {
            components[i] = Component::NonZero {
                mantissa,
                exponent,
                is_negative,
            };
            e_min = e_min.min(exponent);
        }
    }
    Ok((components, e_min))
}

/// Convert a single decomposed component to its scaled `BigInt`
/// representation: `(±mantissa) << (exp − e_min)`.
#[inline]
fn component_to_bigint(c: Component, e_min: i32) -> BigInt {
    match c {
        Component::Zero => BigInt::from(0),
        Component::NonZero {
            mantissa,
            exponent,
            is_negative,
        } => {
            let v = BigInt::from(mantissa.get()) << (exponent - e_min).cast_unsigned();
            if is_negative { -v } else { v }
        }
    }
}

/// Build a `D×D` integer matrix from a component table, scaled to the
/// shared base `2^e_min`.
fn build_bigint_matrix<const D: usize>(
    components: &[[Component; D]; D],
    e_min: i32,
) -> [[BigInt; D]; D] {
    from_fn(|r| from_fn(|c| component_to_bigint(components[r][c], e_min)))
}

/// Build a length-`D` integer vector from a component array, scaled to
/// the shared base `2^e_min`.
fn build_bigint_vec<const D: usize>(components: &[Component; D], e_min: i32) -> [BigInt; D] {
    from_fn(|i| component_to_bigint(components[i], e_min))
}

/// Compute a 2×2 determinant from a scaled integer matrix.
#[inline]
fn det2_bigint<const D: usize>(a: &[[BigInt; D]; D]) -> BigInt {
    &a[0][0] * &a[1][1] - &a[0][1] * &a[1][0]
}

/// Compute a 3×3 determinant from scaled integer entries.
#[inline]
#[allow(clippy::too_many_arguments)]
fn det3_bigint_entries(
    a00: &BigInt,
    a01: &BigInt,
    a02: &BigInt,
    a10: &BigInt,
    a11: &BigInt,
    a12: &BigInt,
    a20: &BigInt,
    a21: &BigInt,
    a22: &BigInt,
) -> BigInt {
    let m00 = a11 * a22 - a12 * a21;
    let m01 = a10 * a22 - a12 * a20;
    let m02 = a10 * a21 - a11 * a20;
    a00 * m00 - a01 * m01 + a02 * m02
}

/// Compute a 3×3 determinant from a scaled integer matrix.
#[inline]
fn det3_bigint<const D: usize>(a: &[[BigInt; D]; D]) -> BigInt {
    det3_bigint_entries(
        &a[0][0], &a[0][1], &a[0][2], &a[1][0], &a[1][1], &a[1][2], &a[2][0], &a[2][1], &a[2][2],
    )
}

/// Compute a 4×4 determinant from a scaled integer matrix.
#[inline]
fn det4_bigint<const D: usize>(a: &[[BigInt; D]; D]) -> BigInt {
    let mut det = BigInt::from(0);

    if a[0][0].sign() != Sign::NoSign {
        let c00 = det3_bigint_entries(
            &a[1][1], &a[1][2], &a[1][3], &a[2][1], &a[2][2], &a[2][3], &a[3][1], &a[3][2],
            &a[3][3],
        );
        det += &a[0][0] * c00;
    }
    if a[0][1].sign() != Sign::NoSign {
        let c01 = det3_bigint_entries(
            &a[1][0], &a[1][2], &a[1][3], &a[2][0], &a[2][2], &a[2][3], &a[3][0], &a[3][2],
            &a[3][3],
        );
        det -= &a[0][1] * c01;
    }
    if a[0][2].sign() != Sign::NoSign {
        let c02 = det3_bigint_entries(
            &a[1][0], &a[1][1], &a[1][3], &a[2][0], &a[2][1], &a[2][3], &a[3][0], &a[3][1],
            &a[3][3],
        );
        det += &a[0][2] * c02;
    }
    if a[0][3].sign() != Sign::NoSign {
        let c03 = det3_bigint_entries(
            &a[1][0], &a[1][1], &a[1][2], &a[2][0], &a[2][1], &a[2][2], &a[3][0], &a[3][1],
            &a[3][2],
        );
        det -= &a[0][3] * c03;
    }

    det
}

/// Outcome of a Bareiss forward-elimination pass.
#[derive(Debug)]
enum BareissResult {
    /// Elimination completed; `sign` is `±1` based on the parity of row
    /// swaps (relevant for determinants; solves discard it).
    Upper { sign: i8 },
    /// Column `pivot_col` has no non-zero pivot at or below its diagonal.
    Singular { pivot_col: usize },
}

/// Run Bareiss fraction-free forward elimination on the `D×D` integer
/// matrix `a`, optionally augmented with a length-`D` RHS vector.
///
/// When `rhs` is `Some`, row swaps and the inner-loop Bareiss update are
/// mirrored on the RHS (treating it as column `D+1` of an augmented
/// system).  On return, `a` is upper triangular and the last pivot lives
/// in `a[D-1][D-1]`.
///
/// First-non-zero pivoting is used: since all arithmetic is exact, any
/// non-zero pivot is valid — no tolerance is required.
fn bareiss_forward_eliminate<const D: usize>(
    a: &mut [[BigInt; D]; D],
    mut rhs: Option<&mut [BigInt; D]>,
) -> BareissResult {
    let zero = BigInt::from(0);
    let mut prev_pivot = BigInt::from(1);
    let mut sign: i8 = 1;

    for k in 0..D {
        // First-non-zero pivot search.
        if a[k][k] == zero {
            let mut found = false;
            for i in (k + 1)..D {
                if a[i][k] != zero {
                    a.swap(k, i);
                    if let Some(r) = &mut rhs {
                        r.swap(k, i);
                    }
                    sign = -sign;
                    found = true;
                    break;
                }
            }
            if !found {
                cold_path();
                return BareissResult::Singular { pivot_col: k };
            }
        }

        // Elimination.  The Bareiss update reads the current `a[i][k]`
        // in both the inner `j`-loop and the RHS update, so zero it only
        // *after* those reads.
        for i in (k + 1)..D {
            for j in (k + 1)..D {
                a[i][j] = (&a[k][k] * &a[i][j] - &a[i][k] * &a[k][j]) / &prev_pivot;
            }
            if let Some(r) = &mut rhs {
                r[i] = (&a[k][k] * &r[i] - &a[i][k] * &r[k]) / &prev_pivot;
            }
            a[i][k].clone_from(&zero);
        }

        prev_pivot.clone_from(&a[k][k]);
    }

    // Post-conditions (debug builds only): `a` is upper triangular with
    // non-zero pivots.  These catch future regressions in the inner-loop
    // update or pivot-search logic without runtime cost in release.
    // Indexed iteration is clearer than iterator chains here because the
    // checks read disjoint cells across rows and columns at each step.
    #[cfg(debug_assertions)]
    #[allow(clippy::needless_range_loop)]
    for k in 0..D {
        assert_ne!(a[k][k], zero, "pivot at ({k}, {k}) must be non-zero");
        for i in (k + 1)..D {
            assert_eq!(a[i][k], zero, "sub-diagonal at ({i}, {k}) must be zero");
        }
    }

    BareissResult::Upper { sign }
}

/// Compute the determinant scale exponent `D × e_min`.
///
/// This centralizes the scale-overflow classification used by public exact
/// determinant APIs: [`Matrix::det_exact`], [`Matrix::det_exact_f64`], and
/// [`Matrix::det_sign_exact`] all surface failures from this helper as
/// [`LaError::DeterminantScaleOverflow`].
///
/// # Errors
/// Returns [`LaError::DeterminantScaleOverflow`] if `D` cannot fit in the
/// internal `i32` exponent multiplier or if `D × e_min` overflows `i32`.
fn determinant_scale_exp<const D: usize>(e_min: i32) -> Result<i32, LaError> {
    let Ok(d_i32) = i32::try_from(D) else {
        cold_path();
        return Err(LaError::determinant_scale_overflow(D, e_min));
    };
    let Some(total_exp) = e_min.checked_mul(d_i32) else {
        cold_path();
        return Err(LaError::determinant_scale_overflow(D, e_min));
    };
    Ok(total_exp)
}

/// Compute the exact determinant from integer-scaled entries.
///
/// Returns `(det_int, scale_exp)` where the true determinant is
/// `det_int × 2^scale_exp`.  Since the scale factor `2^scale_exp` is always
/// positive, `det_int.sign()` gives the sign of the determinant directly.
///
/// All arithmetic is in `BigInt` — no `BigRational`, no GCD, no denominator
/// tracking.  Each f64 entry is decomposed into `mantissa × 2^exponent` and
/// scaled to a common base `2^e_min` so every entry becomes an integer. D≤4
/// uses direct determinant expansions; larger matrices use Bareiss elimination
/// whose inner-loop division is exact (guaranteed by the algorithm).
///
fn bareiss_det_int_finite<const D: usize>(m: &Matrix<D>) -> Result<(BigInt, i32), LaError> {
    // D == 0 has no `a[D-1][D-1]` to read; shortcut to the empty-product
    // determinant.
    if D == 0 {
        return Ok((BigInt::from(1), 0));
    }

    let (components, e_min) = decompose_finite_matrix(m)?;

    // All entries are zero → singular (det = 0).
    if e_min == i32::MAX {
        return Ok((BigInt::from(0), 0));
    }

    let mut a = build_bigint_matrix(&components, e_min);
    let det_int = match D {
        1 => a[0][0].clone(),
        2 => det2_bigint(&a),
        3 => det3_bigint(&a),
        4 => det4_bigint(&a),
        _ => {
            let sign = match bareiss_forward_eliminate(&mut a, None) {
                BareissResult::Upper { sign } => sign,
                BareissResult::Singular { .. } => {
                    cold_path();
                    return Ok((BigInt::from(0), 0));
                }
            };

            if sign < 0 {
                -&a[D - 1][D - 1]
            } else {
                a[D - 1][D - 1].clone()
            }
        }
    };

    let total_exp = determinant_scale_exp::<D>(e_min)?;
    Ok((det_int, total_exp))
}

/// Compute the exact determinant of a `D×D` matrix using integer-only Bareiss
/// elimination and return the result as a `BigRational`.
fn bareiss_det_finite<const D: usize>(m: &Matrix<D>) -> Result<BigRational, LaError> {
    let (det_int, total_exp) = bareiss_det_int_finite(m)?;
    Ok(bigint_exp_to_bigrational(det_int, total_exp))
}

/// Solve `A x = b` exactly after matrix and RHS finiteness has been proven.
///
/// Public [`Matrix`] / [`Vector`] values are finite by construction before
/// reaching this helper, so decomposition can proceed without rediscovering
/// stored NaN/∞ entries.
///
/// # Errors
/// Returns [`LaError::Singular`] if the matrix is exactly singular.
fn gauss_solve_finite<const D: usize>(
    m: &Matrix<D>,
    b: &Vector<D>,
) -> Result<[BigRational; D], LaError> {
    let (m_components, m_e_min) = decompose_finite_matrix(m)?;
    let (b_components, b_e_min) = decompose_finite_vec(b)?;
    gauss_solve_components(m_components, m_e_min, b_components, b_e_min)
}

/// Solve an exact integer-scaled augmented system from decomposed components.
///
/// Forward elimination runs in [`BigInt`] using fraction-free Bareiss updates
/// \[7\]. This is exact arithmetic, so there is no floating-point conditioning or
/// roundoff error in the elimination itself; ill-conditioned inputs can still
/// produce large exact numerators and denominators in the final solution. The
/// elimination phase performs `O(D³)` integer operations and Bareiss exact
/// division controls intermediate integer growth compared with naive fraction
/// arithmetic. The resulting upper-triangular system is then lifted into
/// [`BigRational`] for back-substitution, limiting rational arithmetic to the
/// `O(D²)` phase.
///
/// # Errors
/// Returns [`LaError::Singular`] if the matrix component table represents an
/// exactly singular matrix.
fn gauss_solve_components<const D: usize>(
    m_components: [[Component; D]; D],
    m_e_min: i32,
    b_components: [Component; D],
    b_e_min: i32,
) -> Result<[BigRational; D], LaError> {
    let mut e_min = m_e_min.min(b_e_min);
    if e_min == i32::MAX {
        e_min = 0;
    }

    let mut a = build_bigint_matrix(&m_components, e_min);
    let mut rhs = build_bigint_vec(&b_components, e_min);

    match bareiss_forward_eliminate(&mut a, Some(&mut rhs)) {
        BareissResult::Upper { .. } => {}
        BareissResult::Singular { pivot_col } => {
            cold_path();
            return Err(LaError::Singular { pivot_col });
        }
    }

    let mut x: [BigRational; D] = from_fn(|_| BigRational::from_integer(BigInt::from(0)));
    for i in (0..D).rev() {
        let mut sum = BigRational::from_integer(take(&mut rhs[i]));
        for j in (i + 1)..D {
            let a_ij = BigRational::from_integer(take(&mut a[i][j]));
            sum -= &a_ij * &x[j];
        }
        let a_ii = BigRational::from_integer(take(&mut a[i][i]));
        x[i] = sum / &a_ii;
    }

    Ok(x)
}

/// Exact determinant for a finite-by-construction matrix.
///
/// This is the private implementation target for [`Matrix::det_exact`]. Keeping
/// the helper separate from the public method keeps the exact core focused on
/// the Bareiss computation while relying on the public [`Matrix`] finite-storage
/// invariant.
///
/// # Errors
/// Returns [`LaError::DeterminantScaleOverflow`] if determinant scaling
/// overflows the internal exponent representation.
#[inline]
fn det_exact_finite<const D: usize>(m: &Matrix<D>) -> Result<BigRational, LaError> {
    bareiss_det_finite(m)
}

/// Exact determinant converted to finite `f64` without rounding.
///
/// This preserves the strict contract of [`Matrix::det_exact_f64`]: if the exact
/// determinant is not representable as a finite binary64 value, callers receive
/// a typed [`LaError::Unrepresentable`] instead of a rounded result.
///
/// # Errors
/// Returns [`LaError::DeterminantScaleOverflow`] if determinant scaling
/// overflows the internal exponent representation.
///
/// Returns [`LaError::Unrepresentable`] with
/// [`UnrepresentableReason::RequiresRounding`] when the exact determinant is
/// finite but not exactly representable as binary64, or
/// [`UnrepresentableReason::NotFinite`] when no finite `f64` can represent it.
#[inline]
fn det_exact_f64_finite<const D: usize>(m: &Matrix<D>) -> Result<f64, LaError> {
    let (det_int, total_exp) = bareiss_det_int_finite(m)?;
    bigint_exp_to_finite_f64(det_int, total_exp, None)
}

/// Exact determinant rounded to finite `f64`.
///
/// This is the intentionally lossy counterpart to [`det_exact_f64_finite`] and
/// the private implementation target for [`Matrix::det_exact_rounded_f64`].
///
/// # Errors
/// Returns [`LaError::DeterminantScaleOverflow`] if determinant scaling
/// overflows the internal exponent representation.
///
/// Returns [`LaError::Unrepresentable`] with
/// [`UnrepresentableReason::NotFinite`] if the rounded result would be NaN or
/// infinity.
#[inline]
fn det_exact_rounded_f64_finite<const D: usize>(m: &Matrix<D>) -> Result<f64, LaError> {
    let (det_int, total_exp) = bareiss_det_int_finite(m)?;
    bigint_exp_to_rounded_f64(det_int, total_exp)
}

/// Exact linear solve for finite inputs.
///
/// This is the private implementation target for [`Matrix::solve_exact`].
/// Public [`Matrix`] and [`Vector`] values are finite by construction, so this
/// helper can focus on the exact Bareiss/rational solve.
///
/// # Errors
/// Returns [`LaError::Singular`] if the matrix is exactly singular.
#[inline]
fn solve_exact_finite<const D: usize>(
    m: &Matrix<D>,
    b: Vector<D>,
) -> Result<[BigRational; D], LaError> {
    gauss_solve_finite(m, &b)
}

/// Exact linear solve converted to finite `f64` components without rounding.
///
/// This preserves the strict contract of [`Matrix::solve_exact_f64`]: each exact
/// component must already be representable as finite binary64.
///
/// # Errors
/// Returns [`LaError::Singular`] if the matrix is exactly singular.
///
/// Returns [`LaError::Unrepresentable`] with the failing component index when an
/// exact solution component requires rounding or cannot be represented as a
/// finite `f64`.
#[inline]
fn solve_exact_f64_finite<const D: usize>(
    m: &Matrix<D>,
    b: Vector<D>,
) -> Result<Vector<D>, LaError> {
    let exact = solve_exact_finite(m, b)?;
    let mut result = [0.0f64; D];
    for (i, val) in exact.iter().enumerate() {
        result[i] = exact_rational_to_finite_f64(val, Some(i))?;
    }
    Ok(Vector::new_unchecked(result))
}

/// Exact linear solve rounded to finite `f64` components.
///
/// This is the intentionally lossy counterpart to [`solve_exact_f64_finite`] and
/// the private implementation target for [`Matrix::solve_exact_rounded_f64`].
///
/// # Errors
/// Returns [`LaError::Singular`] if the matrix is exactly singular.
///
/// Returns [`LaError::Unrepresentable`] with
/// [`UnrepresentableReason::NotFinite`] if any rounded component would be NaN or
/// infinity.
#[inline]
fn solve_exact_rounded_f64_finite<const D: usize>(
    m: &Matrix<D>,
    b: Vector<D>,
) -> Result<Vector<D>, LaError> {
    let exact = solve_exact_finite(m, b)?;
    let mut result = [0.0f64; D];
    for (i, val) in exact.iter().enumerate() {
        result[i] = exact_rational_to_rounded_f64(val, Some(i))?;
    }
    Ok(Vector::new_unchecked(result))
}

/// Exact determinant sign for an already finite matrix.
///
/// The fast `f64` filter may reject overflowed scalar intermediates as
/// inconclusive, then fall back to integer Bareiss sign computation.
///
/// # Errors
/// Returns [`LaError::NonFinite`] if a direct determinant or error-bound
/// computation detects a non-finite condition that is not an inconclusive scalar
/// overflow.
///
/// Returns [`LaError::DeterminantScaleOverflow`] if determinant scaling
/// overflows the internal exponent representation.
#[inline]
fn det_sign_exact_finite<const D: usize>(m: &Matrix<D>) -> Result<i8, LaError> {
    match (m.det_direct(), m.det_errbound()) {
        (Ok(Some(det_f64)), Ok(Some(err))) => {
            if det_f64 > err {
                return Ok(1);
            }
            if det_f64 < -err {
                return Ok(-1);
            }
        }
        (Err(LaError::NonFinite { row: None, .. }), _)
        | (_, Err(LaError::NonFinite { row: None, .. })) => {}
        (Err(err), _) | (_, Err(err)) => return Err(err),
        _ => {}
    }

    cold_path();
    let (det_int, _) = bareiss_det_int_finite(m)?;
    Ok(match det_int.sign() {
        Sign::Plus => 1,
        Sign::Minus => -1,
        Sign::NoSign => 0,
    })
}

impl<const D: usize> Matrix<D> {
    /// Exact determinant using arbitrary-precision rational arithmetic.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Returns the determinant as an exact [`BigRational`] value. Every finite
    /// `f64` is exactly representable as a rational, so the conversion is
    /// lossless and the result is provably correct.
    ///
    /// # When to use
    ///
    /// Use this when you need the exact determinant *value* — for example,
    /// exact volume computation or distinguishing truly-degenerate simplices
    /// from near-degenerate ones.  If you only need the *sign*, prefer
    /// [`det_sign_exact`](Self::det_sign_exact) which has a fast f64 filter.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let det = m.det_exact()?;
    /// // det = 1*4 - 2*3 = -2  (exact)
    /// assert_eq!(det, BigRational::from_integer((-2).into()));
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::DeterminantScaleOverflow`] if determinant scaling
    /// overflows the internal exponent representation.
    #[inline]
    pub fn det_exact(&self) -> Result<BigRational, LaError> {
        det_exact_finite(self)
    }

    /// Exact determinant converted to `f64`.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Computes the exact determinant with the same integer Bareiss core used by
    /// [`det_exact`](Self::det_exact), then converts the exact scaled integer
    /// result to `f64` only if the result is exactly representable as a finite
    /// binary64 value.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let det = m.det_exact_f64()?;
    /// assert!((det - (-2.0)).abs() <= f64::EPSILON);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::DeterminantScaleOverflow`] if determinant scaling
    /// overflows the internal exponent representation.
    ///
    /// Returns [`LaError::Unrepresentable`] if the exact determinant cannot be
    /// represented exactly as a finite `f64`.
    #[inline]
    pub fn det_exact_f64(&self) -> Result<f64, LaError> {
        det_exact_f64_finite(self)
    }

    /// Exact determinant rounded to `f64`.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Computes the exact determinant with the same integer Bareiss core used by
    /// [`det_exact`](Self::det_exact), then rounds the exact value to a finite
    /// binary64 value. Unlike [`det_exact_f64`](Self::det_exact_f64), this method
    /// is intentionally lossy and may round non-dyadic or underflowing nonzero
    /// exact determinants.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([
    ///     [1.0 + f64::EPSILON, 0.0],
    ///     [0.0, 1.0 - f64::EPSILON],
    /// ])?;
    ///
    /// assert_eq!(
    ///     m.det_exact_f64(),
    ///     Err(LaError::Unrepresentable {
    ///         index: None,
    ///         reason: UnrepresentableReason::RequiresRounding,
    ///     })
    /// );
    /// assert_eq!(m.det_exact_rounded_f64()?.to_bits(), 1.0f64.to_bits());
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::DeterminantScaleOverflow`] if determinant scaling
    /// overflows the internal exponent representation.
    ///
    /// Returns [`LaError::Unrepresentable`] if the rounded determinant would be
    /// NaN or infinite.
    #[inline]
    pub fn det_exact_rounded_f64(&self) -> Result<f64, LaError> {
        det_exact_rounded_f64_finite(self)
    }

    /// Exact linear system solve using hybrid integer/rational arithmetic.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Solves `A x = b` where `A` is `self` and `b` is the given vector.
    /// Returns the exact solution as `[BigRational; D]`.  Every finite `f64`
    /// is exactly representable as a rational, so the conversion is lossless
    /// and the result is provably correct.
    ///
    /// # When to use
    ///
    /// Use this when you need a provably correct solution — for example,
    /// exact circumcenter computation for near-degenerate simplices where
    /// f64 arithmetic may produce wildly wrong results.
    ///
    /// # Algorithm
    ///
    /// Matrix and RHS entries are decomposed via IEEE 754 bit extraction and
    /// scaled to a shared power-of-two base so the augmented system `(A | b)`
    /// becomes integer-valued.  Forward elimination runs entirely in `BigInt`
    /// with fraction-free Bareiss updates — no `BigRational`, no GCD, no
    /// denominator tracking in the `O(D³)` phase.  Only the upper-triangular
    /// result is lifted into `BigRational` for back-substitution (the `O(D²)`
    /// phase where fractions are inherent).  First-non-zero pivoting is used
    /// throughout; since all arithmetic is exact, any non-zero pivot yields
    /// the correct answer (no numerical-stability concerns).
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // A x = b  where A = [[1,2],[3,4]], b = [5, 11]  →  x = [1, 2]
    /// let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let b = Vector::<2>::try_new([5.0, 11.0])?;
    /// let x = a.solve_exact(b)?;
    /// assert_eq!(x[0], BigRational::from_integer(1.into()));
    /// assert_eq!(x[1], BigRational::from_integer(2.into()));
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if the matrix is exactly singular.
    #[inline]
    pub fn solve_exact(&self, b: Vector<D>) -> Result<[BigRational; D], LaError> {
        solve_exact_finite(self, b)
    }

    /// Exact linear system solve converted to `f64`.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Computes the exact [`BigRational`] solution via
    /// [`solve_exact`](Self::solve_exact) and converts each component to `f64`
    /// only if that component is exactly representable as a finite binary64
    /// value.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let b = Vector::<2>::try_new([5.0, 11.0])?;
    /// let x = a.solve_exact_f64(b)?.into_array();
    /// assert!((x[0] - 1.0).abs() <= f64::EPSILON);
    /// assert!((x[1] - 2.0).abs() <= f64::EPSILON);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if the matrix is exactly singular.
    /// Returns [`LaError::Unrepresentable`] if any component of the exact solution
    /// cannot be represented exactly as a finite `f64`.
    #[inline]
    pub fn solve_exact_f64(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        solve_exact_f64_finite(self, b)
    }

    /// Exact linear system solve rounded to `f64`.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Computes the exact [`BigRational`] solution via
    /// [`solve_exact`](Self::solve_exact) and rounds each component to a finite
    /// binary64 value. Unlike [`solve_exact_f64`](Self::solve_exact_f64), this
    /// method is intentionally lossy and may round non-dyadic or underflowing
    /// nonzero exact components.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<1>::try_from_rows([[3.0]])?;
    /// let b = Vector::<1>::try_new([1.0])?;
    ///
    /// assert_eq!(
    ///     a.solve_exact_f64(b),
    ///     Err(LaError::Unrepresentable {
    ///         index: Some(0),
    ///         reason: UnrepresentableReason::RequiresRounding,
    ///     })
    /// );
    /// assert_eq!(a.solve_exact_rounded_f64(b)?.into_array(), [1.0 / 3.0]);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if the matrix is exactly singular.
    /// Returns [`LaError::Unrepresentable`] if any rounded component would be
    /// NaN or infinite.
    #[inline]
    pub fn solve_exact_rounded_f64(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        solve_exact_rounded_f64_finite(self, b)
    }

    /// Exact determinant sign using adaptive-precision arithmetic.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Returns `1` if `det > 0`, `-1` if `det < 0`, and `0` if `det == 0` (singular).
    ///
    /// For D ≤ 4, a fast f64 filter is tried first: `det_direct()` is compared
    /// against a conservative error bound derived from the matrix permanent.
    /// If the f64 result clearly exceeds the bound, the sign is returned
    /// immediately without allocating.  Otherwise (and always for D ≥ 5),
    /// integer-only Bareiss elimination (`bareiss_det_int`) computes the exact
    /// sign without constructing any `BigRational` values.
    ///
    /// # When to use
    ///
    /// Use this when the sign of the determinant must be correct regardless of
    /// floating-point conditioning (e.g. geometric predicates on near-degenerate
    /// configurations).  For well-conditioned matrices the fast filter resolves
    /// the sign without touching `BigRational`, so the overhead is minimal.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<3>::try_from_rows([
    ///     [1.0, 2.0, 3.0],
    ///     [4.0, 5.0, 6.0],
    ///     [7.0, 8.0, 9.0],
    /// ])?;
    /// // This matrix is singular (row 3 = row 1 + row 2 in exact arithmetic).
    /// assert_eq!(m.det_sign_exact()?, 0);
    ///
    /// assert_eq!(Matrix::<3>::identity().det_sign_exact()?, 1);
    /// # Ok::<(), LaError>(())
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::DeterminantScaleOverflow`] if determinant scaling
    /// overflows the internal exponent representation.
    #[inline]
    pub fn det_sign_exact(&self) -> Result<i8, LaError> {
        det_sign_exact_finite(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_SINGULAR_TOL;

    use core::assert_matches;
    use num_traits::Signed;
    use pastey::paste;
    use std::array::from_fn;

    // -----------------------------------------------------------------------
    // Test helpers
    // -----------------------------------------------------------------------

    /// Build an exact `BigRational` from an `f64` via IEEE 754 bit decomposition.
    ///
    /// Thin wrapper over [`f64_decompose`] that packs the mantissa/exponent
    /// pair into a fully-formed `BigRational` of the form `±m · 2^e`.  The
    /// production code paths (`bareiss_det_int`, `gauss_solve`) instead
    /// decompose every entry into a shared-scale `BigInt` matrix, which
    /// avoids per-entry GCD work in the elimination loops — so this helper
    /// is not used by them and lives here to keep test assertions concise
    /// (e.g. `assert_eq!(x[0], f64_to_bigrational(3.0))`).
    ///
    /// See `REFERENCES.md` \[9-10\] for the IEEE 754 standard and Goldberg's
    /// survey of floating-point representation.
    ///
    /// # Panics
    /// Panics if `x` is NaN or infinite.
    fn f64_to_bigrational(x: f64) -> BigRational {
        let Some((mantissa, exponent, is_negative)) =
            f64_decompose(x).expect("test helper requires finite f64 input")
        else {
            return BigRational::from_integer(BigInt::from(0));
        };

        let numer = if is_negative {
            -BigInt::from(mantissa.get())
        } else {
            BigInt::from(mantissa.get())
        };

        if exponent >= 0 {
            BigRational::new_raw(numer << exponent.cast_unsigned(), BigInt::from(1u32))
        } else {
            BigRational::new_raw(numer, BigInt::from(1u32) << (-exponent).cast_unsigned())
        }
    }

    // -----------------------------------------------------------------------
    // Macro-generated per-dimension tests (D=2..5)
    // -----------------------------------------------------------------------

    macro_rules! gen_internal_matrix_exact_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<internal_matrix_exact_paths_reuse_validation_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let b = Vector::<$d>::new([1.0; $d]);

                    assert_eq!(
                        det_exact_finite(&a).unwrap(),
                        BigRational::from_integer(BigInt::from(1))
                    );
                    assert!((det_exact_f64_finite(&a).unwrap() - 1.0).abs() <= f64::EPSILON);
                    assert_eq!(det_sign_exact_finite(&a).unwrap(), 1);

                    let exact = solve_exact_finite(&a, b).unwrap();
                    for value in exact {
                        assert_eq!(value, BigRational::from_integer(BigInt::from(1)));
                    }

                    let exact_f64 = a.solve_exact_f64(b).unwrap();
                    for value in exact_f64.into_array() {
                        assert!((value - 1.0).abs() <= f64::EPSILON);
                    }
                }
            }
        };
    }

    gen_internal_matrix_exact_tests!(2);
    gen_internal_matrix_exact_tests!(3);
    gen_internal_matrix_exact_tests!(4);
    gen_internal_matrix_exact_tests!(5);

    macro_rules! gen_det_exact_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<det_exact_identity_ $d d>]() {
                    let det = Matrix::<$d>::identity().det_exact().unwrap();
                    assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
                }
            }
        };
    }

    gen_det_exact_tests!(2);
    gen_det_exact_tests!(3);
    gen_det_exact_tests!(4);
    gen_det_exact_tests!(5);

    macro_rules! gen_det_exact_f64_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<det_exact_f64_identity_ $d d>]() {
                    let det = Matrix::<$d>::identity().det_exact_f64().unwrap();
                    assert!((det - 1.0).abs() <= f64::EPSILON);
                }
            }
        };
    }

    gen_det_exact_f64_tests!(2);
    gen_det_exact_f64_tests!(3);
    gen_det_exact_f64_tests!(4);
    gen_det_exact_f64_tests!(5);

    /// For D ≤ 4, `det_exact_f64` should agree with `det_direct` on matrices
    /// whose exact determinant is representable in f64.
    macro_rules! gen_det_exact_f64_agrees_with_det_direct {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<det_exact_f64_agrees_with_det_direct_ $d d>]() {
                    // Power-of-two diagonal entries make the determinant
                    // exactly representable in binary64.
                    let mut rows = [[0.0f64; $d]; $d];
                    let mut value = 2.0;
                    for (i, row) in rows.iter_mut().enumerate() {
                        row[i] = value;
                        value *= 2.0;
                    }
                    let m = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let exact = m.det_exact_f64().unwrap();
                    let direct = m.det_direct().unwrap().unwrap();
                    assert_eq!(exact.to_bits(), direct.to_bits());
                }
            }
        };
    }

    gen_det_exact_f64_agrees_with_det_direct!(2);
    gen_det_exact_f64_agrees_with_det_direct!(3);
    gen_det_exact_f64_agrees_with_det_direct!(4);

    #[test]
    fn det_sign_exact_d0_is_positive() {
        assert_eq!(Matrix::<0>::zero().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_d1_positive() {
        let m = Matrix::<1>::try_from_rows([[42.0]]).unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_d1_negative() {
        let m = Matrix::<1>::try_from_rows([[-3.5]]).unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_d1_zero() {
        let m = Matrix::<1>::try_from_rows([[0.0]]).unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), 0);
    }

    #[test]
    fn det_sign_exact_identity_2d() {
        assert_eq!(Matrix::<2>::identity().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_identity_3d() {
        assert_eq!(Matrix::<3>::identity().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_identity_4d() {
        assert_eq!(Matrix::<4>::identity().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_identity_5d() {
        assert_eq!(Matrix::<5>::identity().det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_singular_duplicate_rows() {
        let m = Matrix::<3>::try_from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [1.0, 2.0, 3.0], // duplicate of row 0
        ])
        .unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), 0);
    }

    #[test]
    fn det_sign_exact_singular_linear_combination() {
        // Row 2 = row 0 + row 1 in exact arithmetic.
        let m = Matrix::<3>::try_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [5.0, 7.0, 9.0]])
            .unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), 0);
    }

    #[test]
    fn det_sign_exact_negative_det_row_swap() {
        // Swapping two rows of the identity negates the determinant.
        let m = Matrix::<3>::try_from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_negative_det_known() {
        // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
        let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_agrees_with_det_for_spd() {
        // SPD matrix → positive determinant.
        let m = Matrix::<3>::try_from_rows([[4.0, 2.0, 0.0], [2.0, 5.0, 1.0], [0.0, 1.0, 3.0]])
            .unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), 1);
        assert!(m.det().unwrap() > 0.0);
    }

    /// Near-singular matrix with an exact perturbation.
    ///
    /// The base matrix `[[1,2,3],[4,5,6],[7,8,9]]` is exactly singular (rows in
    /// arithmetic progression).  Adding `2^-50` to entry (0,0) makes
    /// `det = 2^-50 × cofactor(0,0) = 2^-50 × (5×9 − 6×8) = −3 × 2^-50 < 0`.
    /// Both f64 `det_direct()` and `det_sign_exact()` should agree here.
    #[test]
    fn det_sign_exact_near_singular_perturbation() {
        let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
        let m = Matrix::<3>::try_from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ])
        .unwrap();
        // Exact: det = perturbation × (5×9 − 6×8) = perturbation × (−3) < 0.
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    /// For D ≤ 4, well-conditioned matrices should hit the fast filter
    /// and never allocate `BigRational`.  We can't directly observe this,
    /// but we verify correctness for a range of known signs.
    #[test]
    fn det_sign_exact_fast_filter_positive_4x4() {
        let m = Matrix::<4>::try_from_rows([
            [2.0, 1.0, 0.0, 0.0],
            [1.0, 3.0, 1.0, 0.0],
            [0.0, 1.0, 4.0, 1.0],
            [0.0, 0.0, 1.0, 5.0],
        ])
        .unwrap();
        // SPD tridiagonal → positive det.
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_fast_filter_negative_4x4() {
        // Swap rows 0 and 1 of the above → negate det.
        let m = Matrix::<4>::try_from_rows([
            [1.0, 3.0, 1.0, 0.0],
            [2.0, 1.0, 0.0, 0.0],
            [0.0, 1.0, 4.0, 1.0],
            [0.0, 0.0, 1.0, 5.0],
        ])
        .unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_subnormal_entries() {
        // Subnormal f64 values should convert losslessly.
        let tiny = 5e-324_f64; // smallest positive subnormal
        assert!(tiny.is_subnormal());

        let m = Matrix::<2>::try_from_rows([[tiny, 0.0], [0.0, tiny]]).unwrap();
        // det = tiny^2 > 0
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_returns_err_on_nan() {
        assert_eq!(
            Matrix::<2>::try_from_rows([[f64::NAN, 0.0], [0.0, 1.0]]),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_sign_exact_returns_err_on_infinity() {
        assert_eq!(
            Matrix::<2>::try_from_rows([[f64::INFINITY, 0.0], [0.0, 1.0]]),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn exact_public_methods_reject_unchecked_nonfinite_matrix_before_computation() {
        let m = Matrix::<2>::from_rows_unchecked([[1.0, 0.0], [f64::NAN, 1.0]]);
        let b = Vector::<2>::new([1.0, 1.0]);
        let expected = Err(LaError::NonFinite {
            row: Some(1),
            col: 0,
        });

        assert_eq!(m.det_exact().map(|_| ()), expected);
        assert_eq!(m.det_exact_f64().map(|_| ()), expected);
        assert_eq!(m.det_exact_rounded_f64().map(|_| ()), expected);
        assert_eq!(m.det_sign_exact().map(|_| ()), expected);
        assert_eq!(m.solve_exact(b).map(|_| ()), expected);
        assert_eq!(m.solve_exact_f64(b).map(|_| ()), expected);
        assert_eq!(m.solve_exact_rounded_f64(b).map(|_| ()), expected);
    }

    #[test]
    fn exact_solve_public_methods_reject_unchecked_nonfinite_rhs_before_computation() {
        let a = Matrix::<2>::identity();
        let b = Vector::<2>::new_unchecked([1.0, f64::INFINITY]);
        let expected = Err(LaError::NonFinite { row: None, col: 1 });

        assert_eq!(a.solve_exact(b).map(|_| ()), expected);
        assert_eq!(a.solve_exact_f64(b).map(|_| ()), expected);
        assert_eq!(a.solve_exact_rounded_f64(b).map(|_| ()), expected);
    }

    #[test]
    fn det_sign_exact_returns_err_on_nan_5x5() {
        // D ≥ 5 bypasses the fast filter, exercising the bareiss_det path.
        let mut m = Matrix::<5>::identity();
        assert_eq!(
            m.set(2, 3, f64::NAN),
            Err(LaError::NonFinite {
                row: Some(2),
                col: 3
            })
        );
    }

    #[test]
    fn det_sign_exact_returns_err_on_infinity_5x5() {
        let mut m = Matrix::<5>::identity();
        assert_eq!(
            m.set(0, 0, f64::INFINITY),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_sign_exact_pivot_needed_5x5() {
        // D ≥ 5 skips the fast filter → exercises Bareiss pivoting.
        // Permutation matrix with a single swap (rows 0↔1) → det = −1.
        let m = Matrix::<5>::try_from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ])
        .unwrap();
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_5x5_known() {
        // det of a permutation matrix with two swaps = +1 (even permutation).
        let m = Matrix::<5>::try_from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ])
        .unwrap();
        // Two transpositions → even permutation → det = +1
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    // -----------------------------------------------------------------------
    // Direct tests for internal helpers (coverage of private functions)
    // -----------------------------------------------------------------------

    #[test]
    fn det_errbound_d0_is_zero() {
        assert_eq!(Matrix::<0>::zero().det_errbound(), Ok(Some(0.0)));
    }

    #[test]
    fn det_errbound_d1_is_zero() {
        assert_eq!(
            Matrix::<1>::try_from_rows([[42.0]]).unwrap().det_errbound(),
            Ok(Some(0.0))
        );
    }

    #[test]
    fn det_errbound_d2_positive() {
        let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        let bound = m.det_errbound().unwrap().unwrap();
        assert!(bound > 0.0);
        // bound = ERR_COEFF_2 * (|1*4| + |2*3|) = ERR_COEFF_2 * 10
        assert!(crate::ERR_COEFF_2.mul_add(-10.0, bound).abs() < 1e-30);
    }

    #[test]
    fn det_errbound_d3_positive() {
        let m = Matrix::<3>::identity();
        let bound = m.det_errbound().unwrap().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d3_non_identity() {
        // Non-identity matrix to exercise all code paths in D=3 case
        let m = Matrix::<3>::try_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]])
            .unwrap();
        let bound = m.det_errbound().unwrap().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d4_positive() {
        let m = Matrix::<4>::identity();
        let bound = m.det_errbound().unwrap().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d4_non_identity() {
        // Non-identity matrix to exercise all code paths in D=4 case
        let m = Matrix::<4>::try_from_rows([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 3.0, 0.0],
            [0.0, 0.0, 0.0, 4.0],
        ])
        .unwrap();
        let bound = m.det_errbound().unwrap().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d5_is_none() {
        assert_eq!(Matrix::<5>::identity().det_errbound(), Ok(None));
    }

    // -----------------------------------------------------------------------
    // f64_decompose tests
    // -----------------------------------------------------------------------

    #[test]
    fn f64_decompose_zero() {
        assert!(f64_decompose(0.0).unwrap().is_none());
        assert!(f64_decompose(-0.0).unwrap().is_none());
    }

    #[test]
    fn f64_decompose_one() {
        let (mant, exp, neg) = f64_decompose(1.0).unwrap().unwrap();
        assert_eq!(mant.get(), 1);
        assert_eq!(exp, 0);
        assert!(!neg);
    }

    #[test]
    fn f64_decompose_negative() {
        let (mant, exp, neg) = f64_decompose(-3.5).unwrap().unwrap();
        // -3.5 = -7 × 2^(-1), mantissa is 7 (odd after stripping)
        assert_eq!(mant.get(), 7);
        assert_eq!(exp, -1);
        assert!(neg);
    }

    #[test]
    fn f64_decompose_subnormal() {
        let tiny = 5e-324_f64;
        assert!(tiny.is_subnormal());
        let (mant, exp, neg) = f64_decompose(tiny).unwrap().unwrap();
        assert_eq!(mant.get(), 1);
        assert_eq!(exp, -1074);
        assert!(!neg);
    }

    #[test]
    fn f64_decompose_power_of_two() {
        let (mant, exp, neg) = f64_decompose(1024.0).unwrap().unwrap();
        assert_eq!(mant.get(), 1);
        assert_eq!(exp, 10); // 1024 = 2^10
        assert!(!neg);
    }

    #[test]
    fn f64_decompose_rejects_nan() {
        assert_eq!(
            f64_decompose(f64::NAN),
            Err(LaError::NonFinite { row: None, col: 0 })
        );
    }

    #[test]
    fn component_to_bigint_distinguishes_zero_from_nonzero_mantissa() {
        assert_eq!(
            component_to_bigint(Component::default(), -10),
            BigInt::from(0)
        );

        let positive = Component::NonZero {
            mantissa: NonZeroU64::new(3).unwrap(),
            exponent: 4,
            is_negative: false,
        };
        assert_eq!(component_to_bigint(positive, 1), BigInt::from(24));

        let negative = Component::NonZero {
            mantissa: NonZeroU64::new(5).unwrap(),
            exponent: 3,
            is_negative: true,
        };
        assert_eq!(component_to_bigint(negative, 1), BigInt::from(-20));
    }

    #[test]
    fn determinant_scale_exp_multiplies_dimension_and_min_exponent() {
        assert_eq!(determinant_scale_exp::<4>(-1074), Ok(-4296));
    }

    #[test]
    fn determinant_scale_exp_rejects_dimension_too_large_for_i32() {
        assert_eq!(
            determinant_scale_exp::<{ i32::MAX as usize + 1 }>(-1074),
            Err(LaError::DeterminantScaleOverflow {
                dim: i32::MAX as usize + 1,
                min_exponent: -1074,
            })
        );
    }

    #[test]
    fn determinant_scale_exp_rejects_exponent_product_overflow() {
        assert_eq!(
            determinant_scale_exp::<3_000_000>(-1074),
            Err(LaError::DeterminantScaleOverflow {
                dim: 3_000_000,
                min_exponent: -1074,
            })
        );
    }

    // -----------------------------------------------------------------------
    // bareiss_det_int tests
    // -----------------------------------------------------------------------

    #[test]
    fn bareiss_det_int_d0() {
        let m = Matrix::<0>::zero();
        let (det, exp) = bareiss_det_int_finite(&m).unwrap();
        assert_eq!(det, BigInt::from(1));
        assert_eq!(exp, 0);
    }

    /// Table-driven coverage of the D=1 fast-path: each 1×1 matrix
    /// decomposes to `(±mant, exp)` directly.  Includes an integer, zero,
    /// a negative fractional, and a positive fractional case — the
    /// combinations that exercise the sign handling, the all-zero early
    /// return, trailing-zero stripping, and negative exponent scaling.
    #[test]
    fn bareiss_det_int_d1_cases() {
        let cases: &[(f64, i64, i32)] = &[
            // (input, expected_det_int, expected_exp)
            (7.0, 7, 0),    // integer → (7, 0)
            (0.0, 0, 0),    // all-zero early return → (0, 0)
            (-3.5, -7, -1), // -3.5 = -7 × 2^(-1)
            (0.5, 1, -1),   // 0.5  =  1 × 2^(-1)
        ];
        for &(input, expected_det_int, expected_exp) in cases {
            let m = Matrix::<1>::try_from_rows([[input]]).unwrap();
            let (det, exp) = bareiss_det_int_finite(&m).unwrap();
            assert_eq!(
                det,
                BigInt::from(expected_det_int),
                "det_int for input={input}"
            );
            assert_eq!(exp, expected_exp, "exp for input={input}");
        }
    }

    #[test]
    fn bareiss_det_int_d2_known() {
        // det([[1,2],[3,4]]) = -2
        let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        let (det_int, total_exp) = bareiss_det_int_finite(&m).unwrap();
        // Reconstruct and verify.
        let det = bigint_exp_to_bigrational(det_int, total_exp);
        assert_eq!(det, BigRational::from_integer(BigInt::from(-2)));
    }

    #[test]
    fn bareiss_det_int_all_zeros() {
        let m = Matrix::<3>::zero();
        let (det, _) = bareiss_det_int_finite(&m).unwrap();
        assert_eq!(det, BigInt::from(0));
    }

    #[test]
    fn bareiss_det_int_sign_matches_det_sign_exact() {
        // The sign of det_int should match det_sign_exact for various matrices.
        let m = Matrix::<3>::try_from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let (det_int, _) = bareiss_det_int_finite(&m).unwrap();
        assert_eq!(det_int.sign(), Sign::Minus); // det = -1
    }

    #[test]
    fn bareiss_det_int_fractional_entries() {
        // Entries with negative exponents: 0.5 = 1×2^(-1), 0.25 = 1×2^(-2).
        // det([[0.5, 0.25], [1.0, 1.0]]) = 0.5×1.0 − 0.25×1.0 = 0.25
        let m = Matrix::<2>::try_from_rows([[0.5, 0.25], [1.0, 1.0]]).unwrap();
        let (det_int, total_exp) = bareiss_det_int_finite(&m).unwrap();
        let det = bigint_exp_to_bigrational(det_int, total_exp);
        assert_eq!(det, BigRational::new(BigInt::from(1), BigInt::from(4)));
    }

    #[test]
    fn bareiss_det_int_d3_with_pivoting() {
        // Zero on diagonal → exercises pivot swap inside bareiss_det_int.
        let m = Matrix::<3>::try_from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let (det_int, total_exp) = bareiss_det_int_finite(&m).unwrap();
        let det = bigint_exp_to_bigrational(det_int, total_exp);
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    /// Per AGENTS.md: dimension-generic tests must cover D=2–5.
    macro_rules! gen_bareiss_det_int_identity_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<bareiss_det_int_identity_ $d d>]() {
                    let m = Matrix::<$d>::identity();
                    let (det_int, total_exp) = bareiss_det_int_finite(&m).unwrap();
                    let det = bigint_exp_to_bigrational(det_int, total_exp);
                    assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
                }
            }
        };
    }

    gen_bareiss_det_int_identity_tests!(2);
    gen_bareiss_det_int_identity_tests!(3);
    gen_bareiss_det_int_identity_tests!(4);
    gen_bareiss_det_int_identity_tests!(5);

    // -----------------------------------------------------------------------
    // bigint_exp_to_bigrational tests
    // -----------------------------------------------------------------------

    #[test]
    fn bigint_exp_to_bigrational_zero() {
        let r = bigint_exp_to_bigrational(BigInt::from(0), -50);
        assert_eq!(r, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn bigint_exp_to_bigrational_positive_exp() {
        // 3 × 2^2 = 12
        let r = bigint_exp_to_bigrational(BigInt::from(3), 2);
        assert_eq!(r, BigRational::from_integer(BigInt::from(12)));
    }

    #[test]
    fn bigint_exp_to_bigrational_negative_exp_reduced() {
        // 6 × 2^(-2) = 6/4 → reduced to 3/2 (strip one shared factor of 2)
        let r = bigint_exp_to_bigrational(BigInt::from(6), -2);
        assert_eq!(*r.numer(), BigInt::from(3));
        assert_eq!(*r.denom(), BigInt::from(2));
    }

    #[test]
    fn bigint_exp_to_bigrational_negative_exp_reduces_to_integer() {
        // 8 × 2^(-3) = 1 after stripping every denominator factor.
        let r = bigint_exp_to_bigrational(BigInt::from(8), -3);
        assert_eq!(r, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn bigint_exp_to_bigrational_negative_exp_already_odd() {
        // 3 × 2^(-2) = 3/4 (already in lowest terms since 3 is odd)
        let r = bigint_exp_to_bigrational(BigInt::from(3), -2);
        assert_eq!(*r.numer(), BigInt::from(3));
        assert_eq!(*r.denom(), BigInt::from(4));
    }

    #[test]
    fn bigint_exp_to_bigrational_negative_value() {
        // -5 × 2^1 = -10
        let r = bigint_exp_to_bigrational(BigInt::from(-5), 1);
        assert_eq!(r, BigRational::from_integer(BigInt::from(-10)));
    }

    #[test]
    fn bigint_exp_to_bigrational_negative_value_with_denominator() {
        // -3 × 2^(-2) = -3/4
        let r = bigint_exp_to_bigrational(BigInt::from(-3), -2);
        assert_eq!(*r.numer(), BigInt::from(-3));
        assert_eq!(*r.denom(), BigInt::from(4));
    }

    // -----------------------------------------------------------------------
    // bareiss_det (wrapper) tests
    // -----------------------------------------------------------------------

    #[test]
    fn bareiss_det_d1_returns_entry() {
        let det = Matrix::<1>::try_from_rows([[7.0]])
            .unwrap()
            .det_exact()
            .unwrap();
        assert_eq!(det, f64_to_bigrational(7.0));
    }

    #[test]
    fn bareiss_det_d3_with_pivoting() {
        // First column has zero on diagonal → exercises pivot swap + break.
        let m = Matrix::<3>::try_from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let det = m.det_exact().unwrap();
        // det of this permutation matrix = -1
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    #[test]
    fn bareiss_det_singular_all_zeros_in_column() {
        // Column 1 is all zeros below diagonal after elimination → singular.
        let m = Matrix::<3>::try_from_rows([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn det_sign_exact_overflow_determinant_finite_entries() {
        // Entries near f64::MAX are finite, but the f64 determinant overflows
        // to infinity.  The fast filter should be skipped and Bareiss should
        // compute the correct positive sign.
        let big = f64::MAX / 2.0;
        assert!(big.is_finite());
        let m = Matrix::<3>::try_from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]])
            .unwrap();
        // det = big^2 > 0
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    // -----------------------------------------------------------------------
    // det_exact: dimension-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn det_exact_d0_is_one() {
        let det = Matrix::<0>::zero().det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn det_exact_known_2x2() {
        // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
        let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(-2)));
    }

    #[test]
    fn det_exact_known_dense_4x4() {
        let m = Matrix::<4>::try_from_rows([
            [4.0, 1.0, 3.0, 2.0],
            [0.0, 5.0, 2.0, 1.0],
            [7.0, 2.0, 6.0, 3.0],
            [1.0, 8.0, 4.0, 9.0],
        ])
        .unwrap();

        assert_eq!(
            m.det_exact(),
            Ok(BigRational::from_integer(BigInt::from(92)))
        );
    }

    #[test]
    fn det_exact_singular_returns_zero() {
        // Rows in arithmetic progression → exactly singular.
        let m = Matrix::<3>::try_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
            .unwrap();
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn det_exact_near_singular_perturbation() {
        // Same 2^-50 perturbation case: exact det = -3 × 2^-50.
        let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
        let m = Matrix::<3>::try_from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ])
        .unwrap();
        let det = m.det_exact().unwrap();
        // det should be exactly -3 × 2^-50.
        let expected = BigRational::new(BigInt::from(-3), BigInt::from(1u64 << 50));
        assert_eq!(det, expected);
    }

    #[test]
    fn det_exact_5x5_permutation() {
        // Single swap (rows 0↔1) → det = -1.
        let m = Matrix::<5>::try_from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ])
        .unwrap();
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    // -----------------------------------------------------------------------
    // det_exact_f64: dimension-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn det_exact_f64_known_2x2() {
        let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        let det = m.det_exact_f64().unwrap();
        assert!((det - (-2.0)).abs() <= f64::EPSILON);
    }

    #[test]
    fn det_exact_f64_overflow_returns_err() {
        // Entries near f64::MAX produce a determinant too large for f64.
        let big = f64::MAX / 2.0;
        let m = Matrix::<3>::try_from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]])
            .unwrap();
        // det = big^2, which overflows f64.
        assert_eq!(
            m.det_exact_f64(),
            Err(LaError::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::NotFinite,
            })
        );
    }

    #[test]
    fn det_exact_rounded_f64_overflow_returns_err() {
        let big = f64::MAX / 2.0;
        let m = Matrix::<3>::try_from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]])
            .unwrap();

        assert_eq!(
            m.det_exact_rounded_f64(),
            Err(LaError::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::NotFinite,
            })
        );
    }

    #[test]
    fn det_exact_f64_underflow_returns_err_for_nonzero_exact_result() {
        let tiny = f64::from_bits(1);
        let m = Matrix::<2>::try_from_rows([[tiny, 0.0], [0.0, tiny]]).unwrap();

        assert!(m.det_exact().unwrap().is_positive());
        assert_eq!(
            m.det_exact_f64(),
            Err(LaError::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::RequiresRounding,
            })
        );
    }

    #[test]
    fn det_exact_f64_rejects_inexact_rounding() {
        let m = Matrix::<2>::try_from_rows([[1.0 + f64::EPSILON, 0.0], [0.0, 1.0 - f64::EPSILON]])
            .unwrap();

        assert_eq!(
            m.det_exact(),
            Ok(BigRational::new(
                (BigInt::from(1_u128) << 104_u32) - BigInt::from(1),
                BigInt::from(1_u128 << 104),
            ))
        );
        assert_eq!(
            m.det_exact_f64(),
            Err(LaError::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::RequiresRounding,
            })
        );
    }

    #[test]
    fn det_exact_f64_accepts_min_positive_subnormal() {
        let tiny = f64::from_bits(1);
        let m = Matrix::<1>::try_from_rows([[tiny]]).unwrap();

        assert_eq!(m.det_exact_f64().unwrap().to_bits(), tiny.to_bits());
    }

    #[test]
    fn det_exact_f64_accepts_max_finite_binary64() {
        let m = Matrix::<1>::try_from_rows([[f64::MAX]]).unwrap();

        assert_eq!(m.det_exact_f64().unwrap().to_bits(), f64::MAX.to_bits());
    }

    #[test]
    fn det_exact_rounded_f64_rounds_inexact_result() {
        let m = Matrix::<2>::try_from_rows([[1.0 + f64::EPSILON, 0.0], [0.0, 1.0 - f64::EPSILON]])
            .unwrap();

        assert_eq!(
            m.det_exact_rounded_f64().unwrap().to_bits(),
            1.0f64.to_bits()
        );
    }

    // -----------------------------------------------------------------------
    // solve_exact: macro-generated per-dimension tests (D=2..5)
    // -----------------------------------------------------------------------

    /// Helper: build an arbitrary RHS vector for dimension `$d`.
    fn arbitrary_rhs<const D: usize>() -> Vector<D> {
        let values = [1.0, -2.5, 3.0, 0.25, -4.0];
        let mut arr = [0.0f64; D];
        for (dst, src) in arr.iter_mut().zip(values.iter()) {
            *dst = *src;
        }
        Vector::<D>::new(arr)
    }

    macro_rules! gen_solve_exact_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_identity_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let b = arbitrary_rhs::<$d>();
                    let x = a.solve_exact(b).unwrap();
                    for (i, xi) in x.iter().enumerate() {
                        assert_eq!(*xi, f64_to_bigrational(b.as_array()[i]));
                    }
                }

                #[test]
                fn [<solve_exact_singular_ $d d>]() {
                    // Zero matrix is singular.
                    let a = Matrix::<$d>::zero();
                    let b = arbitrary_rhs::<$d>();
                    assert_eq!(a.solve_exact(b), Err(LaError::Singular { pivot_col: 0 }));
                }
            }
        };
    }

    gen_solve_exact_tests!(2);
    gen_solve_exact_tests!(3);
    gen_solve_exact_tests!(4);
    gen_solve_exact_tests!(5);

    macro_rules! gen_solve_exact_f64_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_f64_identity_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let b = arbitrary_rhs::<$d>();
                    let x = a.solve_exact_f64(b).unwrap().into_array();
                    for i in 0..$d {
                        assert!((x[i] - b.as_array()[i]).abs() <= f64::EPSILON);
                    }
                }
            }
        };
    }

    gen_solve_exact_f64_tests!(2);
    gen_solve_exact_f64_tests!(3);
    gen_solve_exact_f64_tests!(4);
    gen_solve_exact_f64_tests!(5);

    /// For D ≤ 4, `solve_exact_f64` should agree with `Lu::solve` on
    /// well-conditioned matrices.
    macro_rules! gen_solve_exact_f64_agrees_with_lu {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_f64_agrees_with_lu_ $d d>]() {
                    // Diagonally dominant integer matrix with an exactly
                    // representable target solution.  The exact result can
                    // therefore be parsed into f64 without rounding.
                    let mut rows = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            rows[r][c] = if r == c {
                                f64::from($d) + 1.0
                            } else {
                                1.0
                            };
                        }
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let x_true = {
                        let mut arr = [0.0f64; $d];
                        for (dst, src) in arr.iter_mut().zip([1.0, -2.0, 3.0, -4.0, 5.0]) {
                            *dst = src;
                        }
                        arr
                    };
                    let mut b_arr = [0.0f64; $d];
                    for i in 0..$d {
                        let mut sum = 0.0;
                        for j in 0..$d {
                            sum = rows[i][j].mul_add(x_true[j], sum);
                        }
                        b_arr[i] = sum;
                    }
                    let b = Vector::<$d>::new(b_arr);
                    let exact = a.solve_exact_f64(b).unwrap().into_array();
                    let lu_sol = a.lu(DEFAULT_SINGULAR_TOL).unwrap()
                        .solve(b).unwrap().into_array();
                    for i in 0..$d {
                        assert_eq!(exact[i].to_bits(), x_true[i].to_bits());
                        let eps = lu_sol[i].abs().mul_add(1e-12, 1e-12);
                        assert!((exact[i] - lu_sol[i]).abs() <= eps);
                    }
                }
            }
        };
    }

    gen_solve_exact_f64_agrees_with_lu!(2);
    gen_solve_exact_f64_agrees_with_lu!(3);
    gen_solve_exact_f64_agrees_with_lu!(4);
    gen_solve_exact_f64_agrees_with_lu!(5);

    /// Round-trip: for a well-conditioned integer matrix `A` and integer
    /// target `x0`, solving `A x = A x0` must return `x0` exactly.  All
    /// intermediate values stay small enough that `A * x0` is exactly
    /// representable in `f64`, so the round-trip is a precise equality
    /// check on the hybrid BigInt/BigRational path.
    macro_rules! gen_solve_exact_roundtrip_tests {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)]
                fn [<solve_exact_roundtrip_ $d d>]() {
                    // A = D * I + J (diag = D+1, off-diag = 1).  Invertible
                    // for any D >= 1 and cheap to multiply by hand.
                    let mut rows = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            rows[r][c] = if r == c {
                                f64::from($d) + 1.0
                            } else {
                                1.0
                            };
                        }
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();

                    // x0 = [1, 2, ..., D].
                    let mut x0 = [0.0f64; $d];
                    for i in 0..$d {
                        x0[i] = (i + 1) as f64;
                    }

                    // b = A * x0 computed in f64.  With small integers the
                    // multiply-add sequence is exact.
                    let mut b_arr = [0.0f64; $d];
                    for r in 0..$d {
                        let mut sum = 0.0_f64;
                        for c in 0..$d {
                            sum = rows[r][c].mul_add(x0[c], sum);
                        }
                        b_arr[r] = sum;
                    }
                    let b = Vector::<$d>::new(b_arr);

                    let x = a.solve_exact(b).unwrap();
                    for i in 0..$d {
                        assert_eq!(x[i], f64_to_bigrational(x0[i]));
                    }
                }
            }
        };
    }

    gen_solve_exact_roundtrip_tests!(2);
    gen_solve_exact_roundtrip_tests!(3);
    gen_solve_exact_roundtrip_tests!(4);
    gen_solve_exact_roundtrip_tests!(5);

    // -----------------------------------------------------------------------
    // solve_exact: dimension-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn solve_exact_d0_returns_empty() {
        let a = Matrix::<0>::zero();
        let b = Vector::<0>::zero();
        let x = a.solve_exact(b).unwrap();
        assert!(x.is_empty());
    }

    #[test]
    fn solve_exact_known_2x2() {
        // [[1,2],[3,4]] x = [5, 11] → x = [1, 2]
        let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        let b = Vector::<2>::new([5.0, 11.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], BigRational::from_integer(BigInt::from(1)));
        assert_eq!(x[1], BigRational::from_integer(BigInt::from(2)));
    }

    #[test]
    fn solve_exact_pivoting_needed() {
        // First column has zero on diagonal → pivot swap required.
        let a = Matrix::<3>::try_from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let b = Vector::<3>::new([2.0, 3.0, 4.0]);
        let x = a.solve_exact(b).unwrap();
        // x = [3, 2, 4]
        assert_eq!(x[0], f64_to_bigrational(3.0));
        assert_eq!(x[1], f64_to_bigrational(2.0));
        assert_eq!(x[2], f64_to_bigrational(4.0));
    }

    #[test]
    fn solve_exact_fractional_result() {
        // [[2, 1], [1, 3]] x = [1, 1] → x = [2/5, 1/5]
        let a = Matrix::<2>::try_from_rows([[2.0, 1.0], [1.0, 3.0]]).unwrap();
        let b = Vector::<2>::new([1.0, 1.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], BigRational::new(BigInt::from(2), BigInt::from(5)));
        assert_eq!(x[1], BigRational::new(BigInt::from(1), BigInt::from(5)));
    }

    #[test]
    fn solve_exact_singular_duplicate_rows() {
        let a = Matrix::<3>::try_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [1.0, 2.0, 3.0]])
            .unwrap();
        let b = Vector::<3>::new([1.0, 2.0, 3.0]);
        assert_matches!(a.solve_exact(b), Err(LaError::Singular { .. }));
    }

    #[test]
    fn solve_exact_5x5_permutation() {
        // Permutation matrix (swap rows 0↔1): P x = b → x = P^T b.
        let a = Matrix::<5>::try_from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ])
        .unwrap();
        let b = Vector::<5>::new([10.0, 20.0, 30.0, 40.0, 50.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], f64_to_bigrational(20.0));
        assert_eq!(x[1], f64_to_bigrational(10.0));
        assert_eq!(x[2], f64_to_bigrational(30.0));
        assert_eq!(x[3], f64_to_bigrational(40.0));
        assert_eq!(x[4], f64_to_bigrational(50.0));
    }

    /// Entries near `f64::MAX / 2` are finite but their product would
    /// overflow to ±∞ in pure f64 arithmetic.  The `BigInt` augmented-system
    /// path computes the correct solution without any overflow.  The D×D
    /// case uses a diagonal matrix with `big` on every diagonal and a RHS
    /// of `[big, …, big, 0]`, giving the known solution `[1, …, 1, 0]`.
    macro_rules! gen_solve_exact_large_finite_entries_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_large_finite_entries_ $d d>]() {
                    let big = f64::MAX / 2.0;
                    assert!(big.is_finite());
                    // D×D diagonal matrix with `big` on the diagonal.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = big;
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    // RHS = [big, …, big, 0] → x = [1, …, 1, 0].
                    let mut b_arr = [big; $d];
                    b_arr[$d - 1] = 0.0;
                    let b = Vector::<$d>::new(b_arr);
                    let x = a.solve_exact(b).unwrap();
                    for i in 0..($d - 1) {
                        assert_eq!(x[i], BigRational::from_integer(BigInt::from(1)));
                    }
                    assert_eq!(x[$d - 1], BigRational::from_integer(BigInt::from(0)));
                }
            }
        };
    }

    gen_solve_exact_large_finite_entries_tests!(2);
    gen_solve_exact_large_finite_entries_tests!(3);
    gen_solve_exact_large_finite_entries_tests!(4);
    gen_solve_exact_large_finite_entries_tests!(5);

    /// Matrix and RHS entries span many orders of magnitude (from
    /// `f64::MIN_POSITIVE` up through `1e100`).  This exercises the
    /// shared `e_min` scaling: even the largest shift keeps every entry a
    /// representable `BigInt`.  The D×D case alternates `huge`/`tiny`
    /// along the diagonal with a matching RHS, giving `x = [1, …, 1]`.
    macro_rules! gen_solve_exact_mixed_magnitude_entries_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_mixed_magnitude_entries_ $d d>]() {
                    let tiny = f64::MIN_POSITIVE; // 2^-1022, smallest normal
                    let huge = 1.0e100_f64;
                    // Alternate huge/tiny along the diagonal.
                    let mut rows = [[0.0f64; $d]; $d];
                    let mut b_arr = [0.0f64; $d];
                    for i in 0..$d {
                        let val = if i % 2 == 0 { huge } else { tiny };
                        rows[i][i] = val;
                        b_arr[i] = val;
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let b = Vector::<$d>::new(b_arr);
                    let x = a.solve_exact(b).unwrap();
                    for i in 0..$d {
                        assert_eq!(x[i], BigRational::from_integer(BigInt::from(1)));
                    }
                }
            }
        };
    }

    gen_solve_exact_mixed_magnitude_entries_tests!(2);
    gen_solve_exact_mixed_magnitude_entries_tests!(3);
    gen_solve_exact_mixed_magnitude_entries_tests!(4);
    gen_solve_exact_mixed_magnitude_entries_tests!(5);

    /// Subnormal RHS entries must survive the decomposition and
    /// back-substitution paths unchanged.  The D×D case uses the identity
    /// matrix and RHS `[1·tiny, 2·tiny, …, D·tiny]`; each entry remains a
    /// valid subnormal f64 (integer multiples of `2^-1074` fit in the
    /// 52-bit subnormal mantissa for the small integers used here).
    macro_rules! gen_solve_exact_subnormal_rhs_tests {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)]
                fn [<solve_exact_subnormal_rhs_ $d d>]() {
                    let tiny = 5e-324_f64; // smallest positive subnormal
                    assert!(tiny.is_subnormal());
                    let a = Matrix::<$d>::identity();
                    // b[i] = (i+1) · tiny — each entry remains a valid subnormal.
                    let mut b_arr = [0.0f64; $d];
                    for i in 0..$d {
                        b_arr[i] = (i + 1) as f64 * tiny;
                        assert!(b_arr[i].is_subnormal());
                    }
                    let b = Vector::<$d>::new(b_arr);
                    let x = a.solve_exact(b).unwrap();
                    for i in 0..$d {
                        assert_eq!(x[i], f64_to_bigrational((i + 1) as f64 * tiny));
                    }
                }
            }
        };
    }

    gen_solve_exact_subnormal_rhs_tests!(2);
    gen_solve_exact_subnormal_rhs_tests!(3);
    gen_solve_exact_subnormal_rhs_tests!(4);
    gen_solve_exact_subnormal_rhs_tests!(5);

    /// Pivoting path with a zero top-left entry forces a row swap in the
    /// `BigInt` forward-elimination loop and propagates it to the RHS.
    /// Combined with a fractional solution, this exercises the
    /// `BigRational` back-substitution after integer forward elimination.
    ///
    /// The 2×2 block `[[0, 1], [2, 1]]` with rhs `[3, 4]` (→ `x = [1/2, 3]`)
    /// is embedded into the top-left of a D×D identity matrix.  Remaining
    /// rows contribute pass-through equalities `x[i] = b[i]`, so the same
    /// fractional solution appears at indices 0 and 1 regardless of D.
    macro_rules! gen_solve_exact_pivot_swap_fractional_tests {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)]
                // `2..$d` is empty when D=2 (no padded rows); that is the
                // intended behaviour of the macro, not a bug.
                #[allow(clippy::reversed_empty_ranges)]
                fn [<solve_exact_pivot_swap_with_fractional_result_ $d d>]() {
                    // Top-left 2×2: A = [[0, 1], [2, 1]].  After swap:
                    // [[2, 1], [0, 1]], rhs = [4, 3] → x[1] = 3, x[0] = 1/2.
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][1] = 1.0;
                    rows[1][0] = 2.0;
                    rows[1][1] = 1.0;
                    // Identity padding for the remaining rows.
                    for i in 2..$d {
                        rows[i][i] = 1.0;
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    // b = [3, 4, 12, 13, …]; padded entries are arbitrary
                    // finite integers so the identity block gives x[i] = b[i].
                    let mut b_arr = [0.0f64; $d];
                    b_arr[0] = 3.0;
                    b_arr[1] = 4.0;
                    for i in 2..$d {
                        b_arr[i] = (i + 10) as f64;
                    }
                    let b = Vector::<$d>::new(b_arr);
                    let x = a.solve_exact(b).unwrap();
                    assert_eq!(x[0], BigRational::new(BigInt::from(1), BigInt::from(2)));
                    assert_eq!(x[1], BigRational::from_integer(BigInt::from(3)));
                    for i in 2..$d {
                        assert_eq!(x[i], f64_to_bigrational((i + 10) as f64));
                    }
                }
            }
        };
    }

    gen_solve_exact_pivot_swap_fractional_tests!(2);
    gen_solve_exact_pivot_swap_fractional_tests!(3);
    gen_solve_exact_pivot_swap_fractional_tests!(4);
    gen_solve_exact_pivot_swap_fractional_tests!(5);

    /// Mid-elimination pivot swap: the 3×3 block
    /// `[[1, 2, 3], [0, 0, 4], [0, 5, 6]]` has a non-zero pivot at k=0 but
    /// a zero pivot at k=1, so the swap happens *during* forward
    /// elimination rather than at the start.  With rhs `[6, 7, 8]` the
    /// exact solution is `[7/4, -1/2, 7/4]`.  For D > 3 the block is
    /// embedded into the top-left of a D×D identity matrix so the same
    /// fractional solution appears in `x[0..3]` and `x[i] = b[i]` for
    /// `i >= 3`.
    macro_rules! gen_solve_exact_mid_pivot_swap_tests {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)]
                // `3..$d` is empty when D=3 (no padded rows); that is the
                // intended behaviour of the macro, not a bug.
                #[allow(clippy::reversed_empty_ranges)]
                fn [<solve_exact_mid_pivot_swap_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = 1.0; rows[0][1] = 2.0; rows[0][2] = 3.0;
                    // rows[1][0..2] are zero; rows[1][2] = 4.
                    rows[1][2] = 4.0;
                    rows[2][1] = 5.0; rows[2][2] = 6.0;
                    // Identity padding for the remaining rows.
                    for i in 3..$d {
                        rows[i][i] = 1.0;
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let mut b_arr = [0.0f64; $d];
                    b_arr[0] = 6.0;
                    b_arr[1] = 7.0;
                    b_arr[2] = 8.0;
                    for i in 3..$d {
                        b_arr[i] = (i + 10) as f64;
                    }
                    let b = Vector::<$d>::new(b_arr);
                    let x = a.solve_exact(b).unwrap();
                    // x[0..3] = [7/4, -1/2, 7/4].
                    assert_eq!(x[0], BigRational::new(BigInt::from(7), BigInt::from(4)));
                    assert_eq!(x[1], BigRational::new(BigInt::from(-1), BigInt::from(2)));
                    assert_eq!(x[2], BigRational::new(BigInt::from(7), BigInt::from(4)));
                    for i in 3..$d {
                        assert_eq!(x[i], f64_to_bigrational((i + 10) as f64));
                    }
                }
            }
        };
    }

    gen_solve_exact_mid_pivot_swap_tests!(3);
    gen_solve_exact_mid_pivot_swap_tests!(4);
    gen_solve_exact_mid_pivot_swap_tests!(5);

    /// Rank-deficient singular: the last column is identically zero and the
    /// leading `(D-1)×(D-1)` block is full rank, so every intermediate
    /// pivot is non-zero and the singularity surfaces only at the final
    /// column.  The matrix is identity in the top-left `(D-1)×(D-1)` with
    /// a row of ones as the last row (and an all-zero last column), so the
    /// rank is exactly `D-1`.  `solve_exact` must return
    /// `LaError::Singular { pivot_col: D - 1 }`.
    macro_rules! gen_solve_exact_singular_rank_deficient_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_singular_rank_deficient_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..($d - 1) {
                        rows[i][i] = 1.0;
                        rows[$d - 1][i] = 1.0;
                    }
                    // Last column is left all-zero → rank exactly D-1.
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let b = Vector::<$d>::new([1.0; $d]);
                    assert_eq!(
                        a.solve_exact(b),
                        Err(LaError::Singular { pivot_col: $d - 1 })
                    );
                }
            }
        };
    }

    gen_solve_exact_singular_rank_deficient_tests!(2);
    gen_solve_exact_singular_rank_deficient_tests!(3);
    gen_solve_exact_singular_rank_deficient_tests!(4);
    gen_solve_exact_singular_rank_deficient_tests!(5);

    /// Zero RHS with a non-singular matrix.  Every Bareiss update reads
    /// `rhs[k]` and `rhs[i]`, both initialised to zero; every update
    /// produces zero; back-substitution therefore yields `x = 0`
    /// regardless of the matrix entries.  This exercises the
    /// back-substitution `mem::take` path against an all-zero `rhs`.
    macro_rules! gen_solve_exact_zero_rhs_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_zero_rhs_ $d d>]() {
                    // A = D*I + J (diagonally dominant, invertible).
                    let mut rows = [[1.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = f64::from($d) + 1.0;
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let b = Vector::<$d>::zero();
                    let x = a.solve_exact(b).unwrap();
                    for xi in &x {
                        assert_eq!(*xi, BigRational::from_integer(BigInt::from(0)));
                    }
                }
            }
        };
    }

    gen_solve_exact_zero_rhs_tests!(2);
    gen_solve_exact_zero_rhs_tests!(3);
    gen_solve_exact_zero_rhs_tests!(4);
    gen_solve_exact_zero_rhs_tests!(5);

    // -----------------------------------------------------------------------
    // Adversarial-input coverage mirroring `benches/exact.rs`
    // -----------------------------------------------------------------------
    //
    // These tests pin the behaviour of the extreme-input benchmark groups
    // (`exact_near_singular_3x3`, `exact_large_entries_3x3`,
    // `exact_hilbert_{4x4,5x5}`) so a regression would be caught even
    // when benchmarks are not running.

    /// Multiply `A · x` entirely in `BigRational`, using `f64_to_bigrational`
    /// to lift each matrix entry.  Used by residual assertions for inputs
    /// whose exact solution has no closed form we can easily type out.
    fn bigrational_matvec<const D: usize>(a: &Matrix<D>, x: &[BigRational; D]) -> [BigRational; D] {
        from_fn(|i| {
            let mut sum = BigRational::from_integer(BigInt::from(0));
            for (aij, xj) in a.rows()[i].iter().zip(x.iter()) {
                sum += f64_to_bigrational(*aij) * xj;
            }
            sum
        })
    }

    fn hilbert<const D: usize>() -> Matrix<D> {
        let rows = from_fn(|r| from_fn(|c| 1.0 / f64::from(u32::try_from(r + c + 1).unwrap())));
        Matrix::<D>::try_from_rows(rows).unwrap()
    }

    /// Near-singular 3×3 solve (matches the `exact_near_singular_3x3`
    /// bench).  With `A = [[1+2^-50, 2, 3], [4, 5, 6], [7, 8, 9]]` and
    /// `x0 = [1, 1, 1]`, `A · x0 = [6 + 2^-50, 15, 24]`; every component is
    /// exactly representable in `f64` (`6` has ulp `2^-50` at its exponent).
    /// `solve_exact` must recover `x0` exactly — the fractional denominator
    /// introduced by `det(A) = -3 × 2^-50` cancels cleanly against the
    /// augmented RHS.
    #[test]
    fn solve_exact_near_singular_3x3_integer_x0() {
        let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
        let a = Matrix::<3>::try_from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ])
        .unwrap();
        let b = Vector::<3>::new([6.0 + perturbation, 15.0, 24.0]);
        let x = a.solve_exact(b).unwrap();
        let one = BigRational::from_integer(BigInt::from(1));
        assert_eq!(x[0], one);
        assert_eq!(x[1], one);
        assert_eq!(x[2], one);
    }

    #[test]
    fn solve_exact_f64_near_singular_benchmark_rhs_is_rejection_path() {
        let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
        let a = Matrix::<3>::try_from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ])
        .unwrap();
        let b = Vector::<3>::new([1.0, 2.0, 3.0]);

        assert_eq!(
            a.solve_exact_f64(b),
            Err(LaError::Unrepresentable {
                index: Some(2),
                reason: UnrepresentableReason::RequiresRounding,
            })
        );
        assert_eq!(
            a.solve_exact_rounded_f64(b).unwrap().into_array()[2].to_bits(),
            (1.0f64 / 3.0).to_bits()
        );
    }

    /// Large-entry 3×3 solve (matches the `exact_large_entries_3x3`
    /// bench).  `A = big · I + (1 - I)` with `big = f64::MAX / 2` and
    /// `b = [big, 1, 1] = A · [1, 0, 0]`.  The `BigInt` augmented system
    /// sees entries of ~1023 bits on the diagonal and unit entries
    /// elsewhere; Bareiss elimination still produces the exact integer
    /// solution `[1, 0, 0]`.
    #[test]
    fn solve_exact_large_entries_3x3_unit_vector() {
        let big = f64::MAX / 2.0;
        assert!(big.is_finite());
        let a = Matrix::<3>::try_from_rows([[big, 1.0, 1.0], [1.0, big, 1.0], [1.0, 1.0, big]])
            .unwrap();
        let b = Vector::<3>::new([big, 1.0, 1.0]);
        let x = a.solve_exact(b).unwrap();
        let zero = BigRational::from_integer(BigInt::from(0));
        let one = BigRational::from_integer(BigInt::from(1));
        assert_eq!(x[0], one);
        assert_eq!(x[1], zero);
        assert_eq!(x[2], zero);
    }

    /// Determinant of the large-entry 3×3 is roughly `big^3`, which
    /// overflows `f64`. `det_direct()` therefore reports a computed
    /// [`LaError::NonFinite`], the fast filter inside `det_sign_exact`
    /// treats that as inconclusive, and the Bareiss fallback resolves the
    /// positive sign correctly. `det_exact_f64` must report `Unrepresentable`.
    #[test]
    fn det_sign_exact_large_entries_3x3_positive() {
        let big = f64::MAX / 2.0;
        let a = Matrix::<3>::try_from_rows([[big, 1.0, 1.0], [1.0, big, 1.0], [1.0, 1.0, big]])
            .unwrap();
        // Fast filter is inconclusive (big^3 overflows f64 to +∞), so
        // this exercises the Bareiss cold path.
        assert_matches!(a.det_direct(), Err(LaError::NonFinite { row: None, .. }));
        assert_eq!(a.det_sign_exact().unwrap(), 1);
        // Cross-validate: the exact `BigRational` determinant must agree
        // on sign with `det_sign_exact`, and `det_exact_f64` must reject the
        // conversion (the value is representable in BigRational but far exceeds f64).
        assert!(a.det_exact().unwrap().is_positive());
        assert_eq!(
            a.det_exact_f64(),
            Err(LaError::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::NotFinite,
            })
        );
    }

    /// Hilbert matrices are symmetric positive-definite, so
    /// `det_sign_exact` must return `1` for every D.  For D=2..=4 the
    /// fast f64 filter resolves the positive sign without falling
    /// through (Hilbert's determinant is tiny but still well above the
    /// `det_errbound` cushion); for D=5 the filter is skipped entirely
    /// and the Bareiss path handles inputs whose `(mantissa, exponent)`
    /// pairs all differ.
    macro_rules! gen_det_sign_exact_hilbert_positive_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<det_sign_exact_hilbert_positive_ $d d>]() {
                    let h = hilbert::<$d>();
                    assert_eq!(h.det_sign_exact().unwrap(), 1);
                }
            }
        };
    }

    gen_det_sign_exact_hilbert_positive_tests!(2);
    gen_det_sign_exact_hilbert_positive_tests!(3);
    gen_det_sign_exact_hilbert_positive_tests!(4);
    gen_det_sign_exact_hilbert_positive_tests!(5);

    macro_rules! gen_solve_exact_f64_hilbert_benchmark_rhs_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_f64_hilbert_benchmark_rhs_is_rejection_path_ $d d>]() {
                    let h = hilbert::<$d>();
                    let b = Vector::<$d>::new([1.0; $d]);

                    assert_matches!(
                        h.solve_exact_f64(b),
                        Err(LaError::Unrepresentable {
                            reason: UnrepresentableReason::RequiresRounding,
                            ..
                        })
                    );
                    assert_matches!(h.solve_exact_rounded_f64(b), Ok(_));
                }
            }
        };
    }

    gen_solve_exact_f64_hilbert_benchmark_rhs_tests!(4);
    gen_solve_exact_f64_hilbert_benchmark_rhs_tests!(5);

    /// `solve_exact` on a Hilbert matrix must produce a solution whose
    /// residual `A · x - b` is *exactly* zero in `BigRational` arithmetic.
    /// Hilbert entries (`1/3`, `1/5`, `1/6`, `1/7`, …) are non-terminating
    /// in binary, so this is a stronger test than the
    /// `gen_solve_exact_roundtrip_tests` construction (which requires the
    /// RHS to be representable as an exact `f64` product).
    macro_rules! gen_solve_exact_hilbert_residual_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<solve_exact_hilbert_residual_ $d d>]() {
                    let h = hilbert::<$d>();
                    // Use a non-trivial RHS with both positive and negative
                    // entries to avoid accidental structural cancellation.
                    let mut b_arr = [0.0f64; $d];
                    for i in 0usize..$d {
                        let sign = if i.is_multiple_of(2) { 1.0 } else { -1.0 };
                        b_arr[i] = sign * f64::from(u32::try_from(i + 1).unwrap());
                    }
                    let b = Vector::<$d>::new(b_arr);
                    let x = h.solve_exact(b).unwrap();
                    let ax = bigrational_matvec(&h, &x);
                    for i in 0..$d {
                        assert_eq!(ax[i], f64_to_bigrational(b_arr[i]));
                    }
                }
            }
        };
    }

    gen_solve_exact_hilbert_residual_tests!(2);
    gen_solve_exact_hilbert_residual_tests!(3);
    gen_solve_exact_hilbert_residual_tests!(4);
    gen_solve_exact_hilbert_residual_tests!(5);

    // -----------------------------------------------------------------------
    // solve_exact_f64: dimension-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn solve_exact_f64_known_2x2() {
        let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        let b = Vector::<2>::new([5.0, 11.0]);
        let x = a.solve_exact_f64(b).unwrap().into_array();
        assert!((x[0] - 1.0).abs() <= f64::EPSILON);
        assert!((x[1] - 2.0).abs() <= f64::EPSILON);
    }

    #[test]
    fn solve_exact_f64_overflow_returns_err() {
        // [[1/big, 0], [0, 1/big]] x = [big, big] → x = [big², big²],
        // which overflows f64.
        let big = f64::MAX / 2.0;
        let a = Matrix::<2>::try_from_rows([[1.0 / big, 0.0], [0.0, 1.0 / big]]).unwrap();
        let b = Vector::<2>::new([big, big]);
        assert_eq!(
            a.solve_exact_f64(b),
            Err(LaError::Unrepresentable {
                index: Some(0),
                reason: UnrepresentableReason::NotFinite,
            })
        );
    }

    #[test]
    fn solve_exact_f64_huge_non_dyadic_component_returns_not_finite() {
        let a = Matrix::<1>::try_from_rows([[3.0 * f64::MIN_POSITIVE]]).unwrap();
        let b = Vector::<1>::new([f64::MAX]);

        assert_eq!(
            a.solve_exact_f64(b),
            Err(LaError::Unrepresentable {
                index: Some(0),
                reason: UnrepresentableReason::NotFinite,
            })
        );
    }

    #[test]
    fn solve_exact_rounded_f64_overflow_returns_err() {
        let big = f64::MAX / 2.0;
        let a = Matrix::<2>::try_from_rows([[1.0 / big, 0.0], [0.0, 1.0 / big]]).unwrap();
        let b = Vector::<2>::new([big, big]);

        assert_eq!(
            a.solve_exact_rounded_f64(b),
            Err(LaError::Unrepresentable {
                index: Some(0),
                reason: UnrepresentableReason::NotFinite,
            })
        );
    }

    #[test]
    fn solve_exact_f64_underflow_returns_err_for_nonzero_exact_component() {
        let tiny = f64::from_bits(1);
        let a = Matrix::<1>::try_from_rows([[2.0]]).unwrap();
        let b = Vector::<1>::new([tiny]);

        assert_eq!(
            a.solve_exact_f64(b),
            Err(LaError::Unrepresentable {
                index: Some(0),
                reason: UnrepresentableReason::RequiresRounding,
            })
        );
    }

    #[test]
    fn solve_exact_f64_accepts_smallest_subnormal_result() {
        let tiny = f64::from_bits(1);
        let a = Matrix::<1>::identity();
        let b = Vector::<1>::new([tiny]);

        assert_eq!(
            a.solve_exact_f64(b).unwrap().into_array()[0].to_bits(),
            tiny.to_bits()
        );
    }

    #[test]
    fn solve_exact_f64_rejects_non_dyadic_component() {
        let a = Matrix::<1>::try_from_rows([[3.0]]).unwrap();
        let b = Vector::<1>::new([1.0]);

        assert_eq!(
            a.solve_exact_f64(b),
            Err(LaError::Unrepresentable {
                index: Some(0),
                reason: UnrepresentableReason::RequiresRounding,
            })
        );
    }

    #[test]
    fn solve_exact_rounded_f64_rounds_non_dyadic_component() {
        let a = Matrix::<1>::try_from_rows([[3.0]]).unwrap();
        let b = Vector::<1>::new([1.0]);

        assert_eq!(
            a.solve_exact_rounded_f64(b).unwrap().into_array()[0].to_bits(),
            (1.0f64 / 3.0).to_bits()
        );
    }

    // -----------------------------------------------------------------------
    // exact solve boundary tests
    // -----------------------------------------------------------------------

    #[test]
    fn gauss_solve_d1() {
        let a = Matrix::<1>::try_from_rows([[2.0]]).unwrap();
        let b = Vector::<1>::new([6.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], f64_to_bigrational(3.0));
    }

    #[test]
    fn gauss_solve_singular_column_all_zero() {
        let a = Matrix::<3>::try_from_rows([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let b = Vector::<3>::new([1.0, 2.0, 3.0]);
        assert_eq!(a.solve_exact(b), Err(LaError::Singular { pivot_col: 1 }));
    }

    // -----------------------------------------------------------------------
    // f64_to_bigrational tests
    // -----------------------------------------------------------------------

    #[test]
    fn f64_to_bigrational_positive_zero() {
        let r = f64_to_bigrational(0.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn f64_to_bigrational_negative_zero() {
        let r = f64_to_bigrational(-0.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn f64_to_bigrational_one() {
        let r = f64_to_bigrational(1.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn f64_to_bigrational_negative_one() {
        let r = f64_to_bigrational(-1.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(-1)));
    }

    #[test]
    fn f64_to_bigrational_half() {
        let r = f64_to_bigrational(0.5);
        assert_eq!(r, BigRational::new(BigInt::from(1), BigInt::from(2)));
    }

    #[test]
    fn f64_to_bigrational_quarter() {
        let r = f64_to_bigrational(0.25);
        assert_eq!(r, BigRational::new(BigInt::from(1), BigInt::from(4)));
    }

    #[test]
    fn f64_to_bigrational_negative_three_and_a_half() {
        // -3.5 = -7/2
        let r = f64_to_bigrational(-3.5);
        assert_eq!(r, BigRational::new(BigInt::from(-7), BigInt::from(2)));
    }

    #[test]
    fn f64_to_bigrational_integer() {
        let r = f64_to_bigrational(42.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(42)));
    }

    #[test]
    fn f64_to_bigrational_power_of_two() {
        let r = f64_to_bigrational(1024.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(1024)));
    }

    #[test]
    fn f64_to_bigrational_subnormal() {
        let tiny = 5e-324_f64; // smallest positive subnormal
        assert!(tiny.is_subnormal());
        let r = f64_to_bigrational(tiny);
        // 5e-324 = 1 × 2^(-1074)
        assert_eq!(
            r,
            BigRational::new(BigInt::from(1), BigInt::from(1u32) << 1074u32)
        );
    }

    #[test]
    fn f64_to_bigrational_already_lowest_terms() {
        // 0.5 should produce numer=1, denom=2 (already reduced).
        let r = f64_to_bigrational(0.5);
        assert_eq!(*r.numer(), BigInt::from(1));
        assert_eq!(*r.denom(), BigInt::from(2));
    }

    #[test]
    fn f64_to_bigrational_round_trip() {
        // -0.0 is excluded: it maps to BigRational(0) which round-trips
        // to +0.0 (correct; tested separately in f64_to_bigrational_negative_zero).
        let values = [
            0.0,
            1.0,
            -1.0,
            0.5,
            0.25,
            0.1,
            42.0,
            -3.5,
            1e10,
            1e-10,
            f64::MAX / 2.0,
            f64::MIN_POSITIVE,
            5e-324,
        ];
        for &v in &values {
            let r = f64_to_bigrational(v);
            let back = r.to_f64().expect("round-trip to_f64 failed");
            assert!(
                v.to_bits() == back.to_bits(),
                "round-trip failed for {v}: got {back}"
            );
        }
    }

    #[test]
    fn f64_decompose_rejects_nonfinite_inputs() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                f64_decompose(value),
                Err(LaError::NonFinite { row: None, col: 0 })
            );
        }
    }
}
