//! Exact arithmetic operations via arbitrary-precision rational numbers.
//!
//! This module is only compiled when the `"exact"` Cargo feature is enabled.
//!
//! # Architecture
//!
//! ## Determinants
//!
//! All three determinant methods (`det_exact`, `det_exact_f64`, `det_sign_exact`)
//! share the same integer-only Bareiss core (`bareiss_det_int`).  Each f64
//! entry is decomposed via `f64_decompose` into `mantissa × 2^exponent`,
//! all entries are scaled to a common `BigInt` matrix (shifting by
//! `e - e_min`), and Bareiss elimination runs entirely in `BigInt`
//! arithmetic — no `BigRational`, no GCD, no denominator tracking.
//! The result is `(det_int, total_exp)` where `det = det_int × 2^(D × e_min)`.
//! `bareiss_det` wraps this with `bigint_exp_to_bigrational` to reconstruct
//! a reduced `BigRational`; `det_sign_exact` reads the sign directly from
//! `det_int` (the scale factor is always positive).
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
//! `solve_exact` and `solve_exact_f64` solve `A x = b` with a hybrid
//! algorithm that reuses the integer-only Bareiss core used for
//! determinants.  Matrix and RHS entries are decomposed via
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
//! exactly representable as a rational, so the result is provably
//! correct.
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
//! `decompose_matrix` / `decompose_vec` fold an `is_finite()` check
//! into the same pass that decomposes each entry, returning
//! `Err(LaError::NonFinite { row, col })` on the first NaN or ±∞
//! encountered.  This error is propagated through `bareiss_det_int`,
//! `bareiss_det`, and `gauss_solve` via the `?` operator, so every
//! public entry point that reaches the integer-Bareiss core is
//! automatically validated — `f64_decompose` itself is therefore
//! never called with non-finite input from the public API.

use core::hint::cold_path;
use core::mem::take;
use std::array::from_fn;

use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::ToPrimitive;

use crate::LaError;
use crate::matrix::Matrix;
use crate::vector::Vector;

/// Decompose a finite `f64` into its IEEE 754 components.
///
/// Returns `None` for ±0.0, or `Some((mantissa, exponent, is_negative))` where
/// the value is exactly `(-1)^is_negative × mantissa × 2^exponent` and
/// `mantissa` is odd (trailing zeros stripped).  See `REFERENCES.md` \[9-10\].
///
/// # Panics
/// Panics if `x` is NaN or infinite.
fn f64_decompose(x: f64) -> Option<(u64, i32, bool)> {
    let bits = x.to_bits();
    let biased_exp = ((bits >> 52) & 0x7FF) as i32;
    let fraction = bits & 0x000F_FFFF_FFFF_FFFF;

    // ±0.0
    if biased_exp == 0 && fraction == 0 {
        return None;
    }

    // NaN / Inf — callers must validate finiteness before reaching here.
    assert!(biased_exp != 0x7FF, "non-finite f64 in exact conversion");

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
    let exponent = raw_exp + tz.cast_signed();
    let is_negative = bits >> 63 != 0;

    Some((mantissa, exponent, is_negative))
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
        let reduce = tz.min(u64::from((-exp).cast_unsigned()));
        value >>= reduce;
        exp += i32::try_from(reduce).expect("reduce ≤ -exp which fits in i32");
    }

    if exp >= 0 {
        BigRational::new_raw(value << exp.cast_unsigned(), BigInt::from(1u32))
    } else {
        BigRational::new_raw(value, BigInt::from(1u32) << (-exp).cast_unsigned())
    }
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
/// Zero entries have `mantissa == 0`; the other fields are unused in that
/// case.  `Default` yields such a zero component, which is what the
/// per-entry initialiser in `decompose_matrix` / `decompose_vec` produces
/// for ±0.0 cells.
#[derive(Clone, Copy, Default)]
struct Component {
    mantissa: u64,
    exponent: i32,
    is_negative: bool,
}

/// Decompose every entry of a `D×D` matrix via `f64_decompose`,
/// validating finiteness in the same pass.  Returns the per-entry
/// components and the minimum exponent across non-zero entries.  If
/// every entry is zero, the exponent is `i32::MAX`.
///
/// # Errors
/// Returns [`LaError::NonFinite`] with `row: Some(r), col: c` pointing
/// at the first non-finite entry encountered (row-major order).
fn decompose_matrix<const D: usize>(m: &Matrix<D>) -> Result<([[Component; D]; D], i32), LaError> {
    let mut components = [[Component::default(); D]; D];
    let mut e_min = i32::MAX;
    for (r, row) in m.rows.iter().enumerate() {
        for (c, &entry) in row.iter().enumerate() {
            if !entry.is_finite() {
                cold_path();
                return Err(LaError::non_finite_cell(r, c));
            }
            if let Some((mantissa, exponent, is_negative)) = f64_decompose(entry) {
                components[r][c] = Component {
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

/// Decompose every entry of a length-`D` vector via `f64_decompose`,
/// validating finiteness in the same pass.  Returns the per-entry
/// components and the minimum exponent across non-zero entries.  If
/// every entry is zero, the exponent is `i32::MAX`.
///
/// # Errors
/// Returns [`LaError::NonFinite`] with `row: None, col: i` pointing at
/// the first non-finite entry encountered.
fn decompose_vec<const D: usize>(v: &Vector<D>) -> Result<([Component; D], i32), LaError> {
    let mut components = [Component::default(); D];
    let mut e_min = i32::MAX;
    for (i, &entry) in v.data.iter().enumerate() {
        if !entry.is_finite() {
            cold_path();
            return Err(LaError::non_finite_at(i));
        }
        if let Some((mantissa, exponent, is_negative)) = f64_decompose(entry) {
            components[i] = Component {
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
/// representation: `(±mantissa) << (exp − e_min)`.  Zero components map
/// to `BigInt::from(0)`.
#[inline]
fn component_to_bigint(c: Component, e_min: i32) -> BigInt {
    if c.mantissa == 0 {
        BigInt::from(0)
    } else {
        let v = BigInt::from(c.mantissa) << (c.exponent - e_min).cast_unsigned();
        if c.is_negative { -v } else { v }
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

/// Compute the exact determinant using integer-only Bareiss elimination.
///
/// Returns `(det_int, scale_exp)` where the true determinant is
/// `det_int × 2^scale_exp`.  Since the scale factor `2^scale_exp` is always
/// positive, `det_int.sign()` gives the sign of the determinant directly.
///
/// All arithmetic is in `BigInt` — no `BigRational`, no GCD, no denominator
/// tracking.  Each f64 entry is decomposed into `mantissa × 2^exponent` and
/// scaled to a common base `2^e_min` so every entry becomes an integer.
/// The Bareiss inner-loop division is exact (guaranteed by the algorithm).
///
/// # Errors
/// Returns [`LaError::NonFinite`] (propagated from `decompose_matrix`) if
/// any matrix entry is NaN or infinite.
fn bareiss_det_int<const D: usize>(m: &Matrix<D>) -> Result<(BigInt, i32), LaError> {
    // D == 0 has no `a[D-1][D-1]` to read; shortcut to the empty-product
    // determinant.
    if D == 0 {
        return Ok((BigInt::from(1), 0));
    }

    let (components, e_min) = decompose_matrix(m)?;

    // All entries are zero → singular (det = 0).
    if e_min == i32::MAX {
        return Ok((BigInt::from(0), 0));
    }

    let mut a = build_bigint_matrix(&components, e_min);
    let sign = match bareiss_forward_eliminate(&mut a, None) {
        BareissResult::Upper { sign } => sign,
        BareissResult::Singular { .. } => {
            cold_path();
            return Ok((BigInt::from(0), 0));
        }
    };

    let det_int = if sign < 0 {
        -&a[D - 1][D - 1]
    } else {
        a[D - 1][D - 1].clone()
    };

    // det(original) = det_int × 2^(D × e_min)
    let d_i32 = i32::try_from(D).expect("dimension exceeds i32");
    let total_exp = e_min
        .checked_mul(d_i32)
        .expect("exponent overflow in bareiss_det_int");

    Ok((det_int, total_exp))
}

/// Compute the exact determinant of a `D×D` matrix using integer-only Bareiss
/// elimination and return the result as a `BigRational`.
///
/// # Errors
/// Returns [`LaError::NonFinite`] if any matrix entry is NaN or infinite.
fn bareiss_det<const D: usize>(m: &Matrix<D>) -> Result<BigRational, LaError> {
    let (det_int, total_exp) = bareiss_det_int(m)?;
    Ok(bigint_exp_to_bigrational(det_int, total_exp))
}

/// Solve `A x = b` using a hybrid BigInt/BigRational algorithm.
///
/// Forward elimination runs entirely in `BigInt` using fraction-free
/// Bareiss updates on the augmented system `(A | b)`: every f64 entry
/// is decomposed into `mantissa × 2^exponent` and scaled to a shared
/// base `2^e_min` so both the matrix and the RHS become integer.
/// Because the same power-of-two scaling is applied to both sides of
/// `A x = b`, the solution is unchanged.  Row swaps also swap the RHS
/// row; no sign tracking is needed (pivot permutations do not affect
/// the solution of a linear system).
///
/// After forward elimination, the upper-triangular `BigInt` system and
/// its RHS are lifted into `BigRational` for back-substitution, where
/// fractions are inherent.  This keeps the expensive `O(D³)` phase
/// GCD-free and limits `BigRational` work to the cheaper `O(D²)` phase.
///
/// Returns the exact solution as `[BigRational; D]`.
///
/// # Errors
/// Returns [`LaError::NonFinite`] (propagated from `decompose_matrix` /
/// `decompose_vec`) if any matrix or vector entry is NaN or infinite.
/// The matrix is validated before the vector, matching public-API order.
/// Returns [`LaError::Singular`] if the matrix is exactly singular.
fn gauss_solve<const D: usize>(m: &Matrix<D>, b: &Vector<D>) -> Result<[BigRational; D], LaError> {
    // Decompose both matrix and RHS (validating finiteness in one pass);
    // the shared minimum exponent makes every entry of the augmented
    // system an integer after scaling.
    let (m_components, m_e_min) = decompose_matrix(m)?;
    let (b_components, b_e_min) = decompose_vec(b)?;
    let mut e_min = m_e_min.min(b_e_min);

    // All matrix + RHS entries are zero.  For `D > 0` this surfaces as
    // singular inside forward elimination; for `D == 0` the elimination
    // loop body is empty and we return `Ok([])` without touching e_min.
    // Pick any finite value so the shift computation is well-defined (the
    // resulting BigInts are all zero either way).
    if e_min == i32::MAX {
        e_min = 0;
    }

    let mut a = build_bigint_matrix(&m_components, e_min);
    let mut rhs = build_bigint_vec(&b_components, e_min);

    if let BareissResult::Singular { pivot_col } = bareiss_forward_eliminate(&mut a, Some(&mut rhs))
    {
        cold_path();
        return Err(LaError::Singular { pivot_col });
    }

    // Back-substitution in `BigRational`.  Only the upper triangle of `a`
    // and the transformed `rhs` are read, each exactly once — so we
    // `mem::take` instead of `clone` to avoid a per-entry allocation.
    // `BigInt::default()` is the zero value and does not allocate.
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

impl<const D: usize> Matrix<D> {
    /// Exact determinant using arbitrary-precision rational arithmetic.
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
    /// let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let det = m.det_exact().unwrap();
    /// // det = 1*4 - 2*3 = -2  (exact)
    /// assert_eq!(det, BigRational::from_integer((-2).into()));
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix entry is NaN or infinite.
    #[inline]
    pub fn det_exact(&self) -> Result<BigRational, LaError> {
        bareiss_det(self)
    }

    /// Exact determinant converted to `f64`.
    ///
    /// Computes the exact [`BigRational`] determinant via [`det_exact`](Self::det_exact)
    /// and converts it to the nearest `f64`.  This is useful when you want the
    /// most accurate f64 determinant possible without committing to `BigRational`
    /// in your downstream code.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let det = m.det_exact_f64().unwrap();
    /// assert!((det - (-2.0)).abs() <= f64::EPSILON);
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix entry is NaN or infinite.
    /// Returns [`LaError::Overflow`] if the exact determinant is too large to
    /// represent as a finite `f64`.
    #[inline]
    pub fn det_exact_f64(&self) -> Result<f64, LaError> {
        let exact = self.det_exact()?;
        let val = exact.to_f64().unwrap_or(f64::INFINITY);
        if val.is_finite() {
            Ok(val)
        } else {
            cold_path();
            Err(LaError::Overflow { index: None })
        }
    }

    /// Exact linear system solve using hybrid integer/rational arithmetic.
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
    /// // A x = b  where A = [[1,2],[3,4]], b = [5, 11]  →  x = [1, 2]
    /// let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let b = Vector::<2>::new([5.0, 11.0]);
    /// let x = a.solve_exact(b).unwrap();
    /// assert_eq!(x[0], BigRational::from_integer(1.into()));
    /// assert_eq!(x[1], BigRational::from_integer(2.into()));
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix or vector entry is NaN or
    /// infinite.
    /// Returns [`LaError::Singular`] if the matrix is exactly singular.
    #[inline]
    pub fn solve_exact(&self, b: Vector<D>) -> Result<[BigRational; D], LaError> {
        gauss_solve(self, &b)
    }

    /// Exact linear system solve converted to `f64`.
    ///
    /// Computes the exact [`BigRational`] solution via
    /// [`solve_exact`](Self::solve_exact) and converts each component to the
    /// nearest `f64`.  This is useful when you want the most accurate f64
    /// solution possible without committing to `BigRational` downstream.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let b = Vector::<2>::new([5.0, 11.0]);
    /// let x = a.solve_exact_f64(b).unwrap().into_array();
    /// assert!((x[0] - 1.0).abs() <= f64::EPSILON);
    /// assert!((x[1] - 2.0).abs() <= f64::EPSILON);
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix or vector entry is NaN or
    /// infinite.
    /// Returns [`LaError::Singular`] if the matrix is exactly singular.
    /// Returns [`LaError::Overflow`] if any component of the exact solution is
    /// too large to represent as a finite `f64`.
    #[inline]
    pub fn solve_exact_f64(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        let exact = self.solve_exact(b)?;
        let mut result = [0.0f64; D];
        for (i, val) in exact.iter().enumerate() {
            let f = val.to_f64().unwrap_or(f64::INFINITY);
            if !f.is_finite() {
                cold_path();
                return Err(LaError::Overflow { index: Some(i) });
            }
            result[i] = f;
        }
        Ok(Vector::new(result))
    }

    /// Exact determinant sign using adaptive-precision arithmetic.
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
    /// let m = Matrix::<3>::from_rows([
    ///     [1.0, 2.0, 3.0],
    ///     [4.0, 5.0, 6.0],
    ///     [7.0, 8.0, 9.0],
    /// ]);
    /// // This matrix is singular (row 3 = row 1 + row 2 in exact arithmetic).
    /// assert_eq!(m.det_sign_exact().unwrap(), 0);
    ///
    /// assert_eq!(Matrix::<3>::identity().det_sign_exact().unwrap(), 1);
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if any matrix entry is NaN or infinite.
    #[inline]
    pub fn det_sign_exact(&self) -> Result<i8, LaError> {
        // Stage 1: f64 fast filter for D ≤ 4.
        //
        // When entries are large (e.g. near f64::MAX) the determinant can
        // overflow to infinity even though every individual entry is finite.
        // In that case the fast filter is inconclusive; fall through to the
        // exact Bareiss path.  For NaN/±∞ entries IEEE 754 propagates
        // non-finite through `det_direct()`, the `det_f64.is_finite()`
        // guard fails, and we also fall through — validation then happens
        // inside `bareiss_det_int` via `decompose_matrix`.
        match self.det_direct() {
            Some(det_f64)
                if let Some(err) = self.det_errbound()
                    && det_f64.is_finite() =>
            {
                if det_f64 > err {
                    return Ok(1);
                }
                if det_f64 < -err {
                    return Ok(-1);
                }
            }
            _ => {}
        }

        // Stage 2: integer Bareiss fallback — the 2^(D×e_min) scale factor
        // is always positive, so det_int.sign() == det(A).sign().  This is
        // the cold path: the fast filter resolves the vast majority of
        // well-conditioned calls without allocating.  `bareiss_det_int`
        // validates finiteness via `decompose_matrix`, so NaN/±∞ inputs
        // surface here as `Err(LaError::NonFinite)`.
        cold_path();
        let (det_int, _) = bareiss_det_int(self)?;
        Ok(match det_int.sign() {
            Sign::Plus => 1,
            Sign::Minus => -1,
            Sign::NoSign => 0,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_PIVOT_TOL;

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
        let Some((mantissa, exponent, is_negative)) = f64_decompose(x) else {
            return BigRational::from_integer(BigInt::from(0));
        };

        let numer = if is_negative {
            -BigInt::from(mantissa)
        } else {
            BigInt::from(mantissa)
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

    macro_rules! gen_det_exact_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<det_exact_identity_ $d d>]() {
                    let det = Matrix::<$d>::identity().det_exact().unwrap();
                    assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
                }

                #[test]
                fn [<det_exact_err_on_nan_ $d d>]() {
                    let mut m = Matrix::<$d>::identity();
                    m.set(0, 0, f64::NAN);
                    assert_eq!(m.det_exact(), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }

                #[test]
                fn [<det_exact_err_on_inf_ $d d>]() {
                    let mut m = Matrix::<$d>::identity();
                    m.set(0, 0, f64::INFINITY);
                    assert_eq!(m.det_exact(), Err(LaError::NonFinite { row: Some(0), col: 0 }));
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

                #[test]
                fn [<det_exact_f64_err_on_nan_ $d d>]() {
                    let mut m = Matrix::<$d>::identity();
                    m.set(0, 0, f64::NAN);
                    assert_eq!(m.det_exact_f64(), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }
            }
        };
    }

    gen_det_exact_f64_tests!(2);
    gen_det_exact_f64_tests!(3);
    gen_det_exact_f64_tests!(4);
    gen_det_exact_f64_tests!(5);

    /// For D ≤ 4, `det_exact_f64` should agree with `det_direct` on
    /// well-conditioned matrices.
    macro_rules! gen_det_exact_f64_agrees_with_det_direct {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)]
                fn [<det_exact_f64_agrees_with_det_direct_ $d d>]() {
                    // Diagonally dominant → well-conditioned.
                    let mut rows = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            rows[r][c] = if r == c {
                                (r as f64) + f64::from($d) + 1.0
                            } else {
                                0.1 / ((r + c + 1) as f64)
                            };
                        }
                    }
                    let m = Matrix::<$d>::from_rows(rows);
                    let exact = m.det_exact_f64().unwrap();
                    let direct = m.det_direct().unwrap();
                    let eps = direct.abs().mul_add(1e-12, 1e-12);
                    assert!((exact - direct).abs() <= eps);
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
        let m = Matrix::<1>::from_rows([[42.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_d1_negative() {
        let m = Matrix::<1>::from_rows([[-3.5]]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_d1_zero() {
        let m = Matrix::<1>::from_rows([[0.0]]);
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
        let m = Matrix::<3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [1.0, 2.0, 3.0], // duplicate of row 0
        ]);
        assert_eq!(m.det_sign_exact().unwrap(), 0);
    }

    #[test]
    fn det_sign_exact_singular_linear_combination() {
        // Row 2 = row 0 + row 1 in exact arithmetic.
        let m = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [5.0, 7.0, 9.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), 0);
    }

    #[test]
    fn det_sign_exact_negative_det_row_swap() {
        // Swapping two rows of the identity negates the determinant.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_negative_det_known() {
        // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_agrees_with_det_for_spd() {
        // SPD matrix → positive determinant.
        let m = Matrix::<3>::from_rows([[4.0, 2.0, 0.0], [2.0, 5.0, 1.0], [0.0, 1.0, 3.0]]);
        assert_eq!(m.det_sign_exact().unwrap(), 1);
        assert!(m.det(DEFAULT_PIVOT_TOL).unwrap() > 0.0);
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
        let m = Matrix::<3>::from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);
        // Exact: det = perturbation × (5×9 − 6×8) = perturbation × (−3) < 0.
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    /// For D ≤ 4, well-conditioned matrices should hit the fast filter
    /// and never allocate `BigRational`.  We can't directly observe this,
    /// but we verify correctness for a range of known signs.
    #[test]
    fn det_sign_exact_fast_filter_positive_4x4() {
        let m = Matrix::<4>::from_rows([
            [2.0, 1.0, 0.0, 0.0],
            [1.0, 3.0, 1.0, 0.0],
            [0.0, 1.0, 4.0, 1.0],
            [0.0, 0.0, 1.0, 5.0],
        ]);
        // SPD tridiagonal → positive det.
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_fast_filter_negative_4x4() {
        // Swap rows 0 and 1 of the above → negate det.
        let m = Matrix::<4>::from_rows([
            [1.0, 3.0, 1.0, 0.0],
            [2.0, 1.0, 0.0, 0.0],
            [0.0, 1.0, 4.0, 1.0],
            [0.0, 0.0, 1.0, 5.0],
        ]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_subnormal_entries() {
        // Subnormal f64 values should convert losslessly.
        let tiny = 5e-324_f64; // smallest positive subnormal
        assert!(tiny.is_subnormal());

        let m = Matrix::<2>::from_rows([[tiny, 0.0], [0.0, tiny]]);
        // det = tiny^2 > 0
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    #[test]
    fn det_sign_exact_returns_err_on_nan() {
        let m = Matrix::<2>::from_rows([[f64::NAN, 0.0], [0.0, 1.0]]);
        assert_eq!(
            m.det_sign_exact(),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_sign_exact_returns_err_on_infinity() {
        let m = Matrix::<2>::from_rows([[f64::INFINITY, 0.0], [0.0, 1.0]]);
        assert_eq!(
            m.det_sign_exact(),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_sign_exact_returns_err_on_nan_5x5() {
        // D ≥ 5 bypasses the fast filter, exercising the bareiss_det path.
        let mut m = Matrix::<5>::identity();
        m.set(2, 3, f64::NAN);
        assert_eq!(
            m.det_sign_exact(),
            Err(LaError::NonFinite {
                row: Some(2),
                col: 3
            })
        );
    }

    #[test]
    fn det_sign_exact_returns_err_on_infinity_5x5() {
        let mut m = Matrix::<5>::identity();
        m.set(0, 0, f64::INFINITY);
        assert_eq!(
            m.det_sign_exact(),
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
        let m = Matrix::<5>::from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        assert_eq!(m.det_sign_exact().unwrap(), -1);
    }

    #[test]
    fn det_sign_exact_5x5_known() {
        // det of a permutation matrix with two swaps = +1 (even permutation).
        let m = Matrix::<5>::from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        // Two transpositions → even permutation → det = +1
        assert_eq!(m.det_sign_exact().unwrap(), 1);
    }

    // -----------------------------------------------------------------------
    // Direct tests for internal helpers (coverage of private functions)
    // -----------------------------------------------------------------------

    #[test]
    fn det_errbound_d0_is_zero() {
        assert_eq!(Matrix::<0>::zero().det_errbound(), Some(0.0));
    }

    #[test]
    fn det_errbound_d1_is_zero() {
        assert_eq!(Matrix::<1>::from_rows([[42.0]]).det_errbound(), Some(0.0));
    }

    #[test]
    fn det_errbound_d2_positive() {
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
        // bound = ERR_COEFF_2 * (|1*4| + |2*3|) = ERR_COEFF_2 * 10
        assert!(crate::ERR_COEFF_2.mul_add(-10.0, bound).abs() < 1e-30);
    }

    #[test]
    fn det_errbound_d3_positive() {
        let m = Matrix::<3>::identity();
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d3_non_identity() {
        // Non-identity matrix to exercise all code paths in D=3 case
        let m = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]]);
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d4_positive() {
        let m = Matrix::<4>::identity();
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d4_non_identity() {
        // Non-identity matrix to exercise all code paths in D=4 case
        let m = Matrix::<4>::from_rows([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 3.0, 0.0],
            [0.0, 0.0, 0.0, 4.0],
        ]);
        let bound = m.det_errbound().unwrap();
        assert!(bound > 0.0);
    }

    #[test]
    fn det_errbound_d5_is_none() {
        assert_eq!(Matrix::<5>::identity().det_errbound(), None);
    }

    // -----------------------------------------------------------------------
    // f64_decompose tests
    // -----------------------------------------------------------------------

    #[test]
    fn f64_decompose_zero() {
        assert!(f64_decompose(0.0).is_none());
        assert!(f64_decompose(-0.0).is_none());
    }

    #[test]
    fn f64_decompose_one() {
        let (mant, exp, neg) = f64_decompose(1.0).unwrap();
        assert_eq!(mant, 1);
        assert_eq!(exp, 0);
        assert!(!neg);
    }

    #[test]
    fn f64_decompose_negative() {
        let (mant, exp, neg) = f64_decompose(-3.5).unwrap();
        // -3.5 = -7 × 2^(-1), mantissa is 7 (odd after stripping)
        assert_eq!(mant, 7);
        assert_eq!(exp, -1);
        assert!(neg);
    }

    #[test]
    fn f64_decompose_subnormal() {
        let tiny = 5e-324_f64;
        assert!(tiny.is_subnormal());
        let (mant, exp, neg) = f64_decompose(tiny).unwrap();
        assert_eq!(mant, 1);
        assert_eq!(exp, -1074);
        assert!(!neg);
    }

    #[test]
    fn f64_decompose_power_of_two() {
        let (mant, exp, neg) = f64_decompose(1024.0).unwrap();
        assert_eq!(mant, 1);
        assert_eq!(exp, 10); // 1024 = 2^10
        assert!(!neg);
    }

    #[test]
    #[should_panic(expected = "non-finite f64 in exact conversion")]
    fn f64_decompose_panics_on_nan() {
        f64_decompose(f64::NAN);
    }

    // -----------------------------------------------------------------------
    // bareiss_det_int tests
    // -----------------------------------------------------------------------

    #[test]
    fn bareiss_det_int_d0() {
        let (det, exp) = bareiss_det_int(&Matrix::<0>::zero()).unwrap();
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
            let (det, exp) = bareiss_det_int(&Matrix::<1>::from_rows([[input]])).unwrap();
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
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let (det_int, total_exp) = bareiss_det_int(&m).unwrap();
        // Reconstruct and verify.
        let det = bigint_exp_to_bigrational(det_int, total_exp);
        assert_eq!(det, BigRational::from_integer(BigInt::from(-2)));
    }

    #[test]
    fn bareiss_det_int_all_zeros() {
        let (det, _) = bareiss_det_int(&Matrix::<3>::zero()).unwrap();
        assert_eq!(det, BigInt::from(0));
    }

    #[test]
    fn bareiss_det_int_sign_matches_det_sign_exact() {
        // The sign of det_int should match det_sign_exact for various matrices.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let (det_int, _) = bareiss_det_int(&m).unwrap();
        assert_eq!(det_int.sign(), Sign::Minus); // det = -1
    }

    #[test]
    fn bareiss_det_int_fractional_entries() {
        // Entries with negative exponents: 0.5 = 1×2^(-1), 0.25 = 1×2^(-2).
        // det([[0.5, 0.25], [1.0, 1.0]]) = 0.5×1.0 − 0.25×1.0 = 0.25
        let m = Matrix::<2>::from_rows([[0.5, 0.25], [1.0, 1.0]]);
        let (det_int, total_exp) = bareiss_det_int(&m).unwrap();
        let det = bigint_exp_to_bigrational(det_int, total_exp);
        assert_eq!(det, BigRational::new(BigInt::from(1), BigInt::from(4)));
    }

    #[test]
    fn bareiss_det_int_d3_with_pivoting() {
        // Zero on diagonal → exercises pivot swap inside bareiss_det_int.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let (det_int, total_exp) = bareiss_det_int(&m).unwrap();
        let det = bigint_exp_to_bigrational(det_int, total_exp);
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    /// Non-finite matrix entries surface as `LaError::NonFinite` with the
    /// row/col of the first offending entry.
    #[test]
    fn bareiss_det_int_errs_on_nan() {
        let mut m = Matrix::<3>::identity();
        m.set(1, 2, f64::NAN);
        assert_eq!(
            bareiss_det_int(&m),
            Err(LaError::NonFinite {
                row: Some(1),
                col: 2
            })
        );
    }

    #[test]
    fn bareiss_det_int_errs_on_inf() {
        let mut m = Matrix::<2>::identity();
        m.set(0, 0, f64::INFINITY);
        assert_eq!(
            bareiss_det_int(&m),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    /// Per AGENTS.md: dimension-generic tests must cover D=2–5.
    macro_rules! gen_bareiss_det_int_identity_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<bareiss_det_int_identity_ $d d>]() {
                    let (det_int, total_exp) =
                        bareiss_det_int(&Matrix::<$d>::identity()).unwrap();
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
    fn bareiss_det_d0_is_one() {
        let det = bareiss_det(&Matrix::<0>::zero()).unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn bareiss_det_d1_returns_entry() {
        let det = bareiss_det(&Matrix::<1>::from_rows([[7.0]])).unwrap();
        assert_eq!(det, f64_to_bigrational(7.0));
    }

    #[test]
    fn bareiss_det_d3_with_pivoting() {
        // First column has zero on diagonal → exercises pivot swap + break.
        let m = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let det = bareiss_det(&m).unwrap();
        // det of this permutation matrix = -1
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    #[test]
    fn bareiss_det_singular_all_zeros_in_column() {
        // Column 1 is all zeros below diagonal after elimination → singular.
        let m = Matrix::<3>::from_rows([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let det = bareiss_det(&m).unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn det_sign_exact_overflow_determinant_finite_entries() {
        // Entries near f64::MAX are finite, but the f64 determinant overflows
        // to infinity.  The fast filter should be skipped and Bareiss should
        // compute the correct positive sign.
        let big = f64::MAX / 2.0;
        assert!(big.is_finite());
        let m = Matrix::<3>::from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]]);
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
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(-2)));
    }

    #[test]
    fn det_exact_singular_returns_zero() {
        // Rows in arithmetic progression → exactly singular.
        let m = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn det_exact_near_singular_perturbation() {
        // Same 2^-50 perturbation case: exact det = -3 × 2^-50.
        let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
        let m = Matrix::<3>::from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);
        let det = m.det_exact().unwrap();
        // det should be exactly -3 × 2^-50.
        let expected = BigRational::new(BigInt::from(-3), BigInt::from(1u64 << 50));
        assert_eq!(det, expected);
    }

    #[test]
    fn det_exact_5x5_permutation() {
        // Single swap (rows 0↔1) → det = -1.
        let m = Matrix::<5>::from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    // -----------------------------------------------------------------------
    // det_exact_f64: dimension-specific tests
    // -----------------------------------------------------------------------

    #[test]
    fn det_exact_f64_known_2x2() {
        let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let det = m.det_exact_f64().unwrap();
        assert!((det - (-2.0)).abs() <= f64::EPSILON);
    }

    #[test]
    fn det_exact_f64_overflow_returns_err() {
        // Entries near f64::MAX produce a determinant too large for f64.
        let big = f64::MAX / 2.0;
        let m = Matrix::<3>::from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]]);
        // det = big^2, which overflows f64.
        assert_eq!(m.det_exact_f64(), Err(LaError::Overflow { index: None }));
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
                        assert_eq!(*xi, f64_to_bigrational(b.data[i]));
                    }
                }

                #[test]
                fn [<solve_exact_err_on_nan_matrix_ $d d>]() {
                    let mut a = Matrix::<$d>::identity();
                    a.set(0, 0, f64::NAN);
                    let b = arbitrary_rhs::<$d>();
                    assert_eq!(a.solve_exact(b), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }

                #[test]
                fn [<solve_exact_err_on_inf_matrix_ $d d>]() {
                    let mut a = Matrix::<$d>::identity();
                    a.set(0, 0, f64::INFINITY);
                    let b = arbitrary_rhs::<$d>();
                    assert_eq!(a.solve_exact(b), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }

                #[test]
                fn [<solve_exact_err_on_nan_vector_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let mut b_arr = [1.0f64; $d];
                    b_arr[0] = f64::NAN;
                    let b = Vector::<$d>::new(b_arr);
                    assert_eq!(a.solve_exact(b), Err(LaError::NonFinite { row: None, col: 0 }));
                }

                #[test]
                fn [<solve_exact_err_on_inf_vector_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let mut b_arr = [1.0f64; $d];
                    b_arr[$d - 1] = f64::INFINITY;
                    let b = Vector::<$d>::new(b_arr);
                    assert_eq!(a.solve_exact(b), Err(LaError::NonFinite { row: None, col: $d - 1 }));
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
                        assert!((x[i] - b.data[i]).abs() <= f64::EPSILON);
                    }
                }

                #[test]
                fn [<solve_exact_f64_err_on_nan_ $d d>]() {
                    let mut a = Matrix::<$d>::identity();
                    a.set(0, 0, f64::NAN);
                    let b = arbitrary_rhs::<$d>();
                    assert_eq!(a.solve_exact_f64(b), Err(LaError::NonFinite { row: Some(0), col: 0 }));
                }
            }
        };
    }

    gen_solve_exact_f64_tests!(2);
    gen_solve_exact_f64_tests!(3);
    gen_solve_exact_f64_tests!(4);
    gen_solve_exact_f64_tests!(5);

    /// For D ≤ 4, `solve_exact_f64` should agree with `Lu::solve_vec` on
    /// well-conditioned matrices.
    macro_rules! gen_solve_exact_f64_agrees_with_lu {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)]
                fn [<solve_exact_f64_agrees_with_lu_ $d d>]() {
                    // Diagonally dominant → well-conditioned.
                    let mut rows = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            rows[r][c] = if r == c {
                                (r as f64) + f64::from($d) + 1.0
                            } else {
                                0.1 / ((r + c + 1) as f64)
                            };
                        }
                    }
                    let a = Matrix::<$d>::from_rows(rows);
                    let b = arbitrary_rhs::<$d>();
                    let exact = a.solve_exact_f64(b).unwrap().into_array();
                    let lu_sol = a.lu(DEFAULT_PIVOT_TOL).unwrap()
                        .solve_vec(b).unwrap().into_array();
                    for i in 0..$d {
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
                    let a = Matrix::<$d>::from_rows(rows);

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
                            sum += rows[r][c] * x0[c];
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
        let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
        let b = Vector::<2>::new([5.0, 11.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], BigRational::from_integer(BigInt::from(1)));
        assert_eq!(x[1], BigRational::from_integer(BigInt::from(2)));
    }

    #[test]
    fn solve_exact_pivoting_needed() {
        // First column has zero on diagonal → pivot swap required.
        let a = Matrix::<3>::from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
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
        let a = Matrix::<2>::from_rows([[2.0, 1.0], [1.0, 3.0]]);
        let b = Vector::<2>::new([1.0, 1.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], BigRational::new(BigInt::from(2), BigInt::from(5)));
        assert_eq!(x[1], BigRational::new(BigInt::from(1), BigInt::from(5)));
    }

    #[test]
    fn solve_exact_singular_duplicate_rows() {
        let a = Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [1.0, 2.0, 3.0]]);
        let b = Vector::<3>::new([1.0, 2.0, 3.0]);
        assert!(matches!(a.solve_exact(b), Err(LaError::Singular { .. })));
    }

    #[test]
    fn solve_exact_5x5_permutation() {
        // Permutation matrix (swap rows 0↔1): P x = b → x = P^T b.
        let a = Matrix::<5>::from_rows([
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);
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
                    let a = Matrix::<$d>::from_rows(rows);
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
                    let a = Matrix::<$d>::from_rows(rows);
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
                    let a = Matrix::<$d>::from_rows(rows);
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
                    let a = Matrix::<$d>::from_rows(rows);
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
                    let a = Matrix::<$d>::from_rows(rows);
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
                    let a = Matrix::<$d>::from_rows(rows);
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
            for (aij, xj) in a.rows[i].iter().zip(x.iter()) {
                sum += f64_to_bigrational(*aij) * xj;
            }
            sum
        })
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
        let a = Matrix::<3>::from_rows([
            [1.0 + perturbation, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]);
        let b = Vector::<3>::new([6.0 + perturbation, 15.0, 24.0]);
        let x = a.solve_exact(b).unwrap();
        let one = BigRational::from_integer(BigInt::from(1));
        assert_eq!(x[0], one);
        assert_eq!(x[1], one);
        assert_eq!(x[2], one);
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
        let a = Matrix::<3>::from_rows([[big, 1.0, 1.0], [1.0, big, 1.0], [1.0, 1.0, big]]);
        let b = Vector::<3>::new([big, 1.0, 1.0]);
        let x = a.solve_exact(b).unwrap();
        let zero = BigRational::from_integer(BigInt::from(0));
        let one = BigRational::from_integer(BigInt::from(1));
        assert_eq!(x[0], one);
        assert_eq!(x[1], zero);
        assert_eq!(x[2], zero);
    }

    /// Determinant of the large-entry 3×3 is roughly `big^3`, which
    /// overflows `f64`.  `det_direct()` therefore returns `±∞`, the fast
    /// filter inside `det_sign_exact` falls through on the `is_finite()`
    /// guard, and the Bareiss fallback resolves the positive sign
    /// correctly.  `det_exact_f64` must report `Overflow`.
    #[test]
    fn det_sign_exact_large_entries_3x3_positive() {
        let big = f64::MAX / 2.0;
        let a = Matrix::<3>::from_rows([[big, 1.0, 1.0], [1.0, big, 1.0], [1.0, 1.0, big]]);
        // Fast filter is inconclusive (big^3 overflows f64 to +∞), so
        // this exercises the Bareiss cold path.
        assert!(!a.det_direct().is_some_and(f64::is_finite));
        assert_eq!(a.det_sign_exact().unwrap(), 1);
        // Cross-validate: the exact `BigRational` determinant must agree
        // on sign with `det_sign_exact`, and `det_exact_f64` must overflow
        // (the value is representable in BigRational but far exceeds f64).
        assert!(a.det_exact().unwrap().is_positive());
        assert_eq!(a.det_exact_f64(), Err(LaError::Overflow { index: None }));
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
                #[allow(clippy::cast_precision_loss)]
                fn [<det_sign_exact_hilbert_positive_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            rows[r][c] = 1.0 / ((r + c + 1) as f64);
                        }
                    }
                    let h = Matrix::<$d>::from_rows(rows);
                    assert_eq!(h.det_sign_exact().unwrap(), 1);
                }
            }
        };
    }

    gen_det_sign_exact_hilbert_positive_tests!(2);
    gen_det_sign_exact_hilbert_positive_tests!(3);
    gen_det_sign_exact_hilbert_positive_tests!(4);
    gen_det_sign_exact_hilbert_positive_tests!(5);

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
                #[allow(clippy::cast_precision_loss)]
                fn [<solve_exact_hilbert_residual_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            rows[r][c] = 1.0 / ((r + c + 1) as f64);
                        }
                    }
                    let h = Matrix::<$d>::from_rows(rows);
                    // Use a non-trivial RHS with both positive and negative
                    // entries to avoid accidental structural cancellation.
                    let mut b_arr = [0.0f64; $d];
                    for i in 0usize..$d {
                        let sign = if i.is_multiple_of(2) { 1.0 } else { -1.0 };
                        b_arr[i] = sign * ((i + 1) as f64);
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
        let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
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
        let a = Matrix::<2>::from_rows([[1.0 / big, 0.0], [0.0, 1.0 / big]]);
        let b = Vector::<2>::new([big, big]);
        assert_eq!(
            a.solve_exact_f64(b),
            Err(LaError::Overflow { index: Some(0) })
        );
    }

    // -----------------------------------------------------------------------
    // gauss_solve: internal helper tests
    // -----------------------------------------------------------------------

    #[test]
    fn gauss_solve_d0_returns_empty() {
        let a = Matrix::<0>::zero();
        let b = Vector::<0>::zero();
        assert_eq!(gauss_solve(&a, &b).unwrap().len(), 0);
    }

    #[test]
    fn gauss_solve_d1() {
        let a = Matrix::<1>::from_rows([[2.0]]);
        let b = Vector::<1>::new([6.0]);
        let x = gauss_solve(&a, &b).unwrap();
        assert_eq!(x[0], f64_to_bigrational(3.0));
    }

    #[test]
    fn gauss_solve_singular_column_all_zero() {
        let a = Matrix::<3>::from_rows([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let b = Vector::<3>::new([1.0, 2.0, 3.0]);
        assert_eq!(gauss_solve(&a, &b), Err(LaError::Singular { pivot_col: 1 }));
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
    #[should_panic(expected = "non-finite f64 in exact conversion")]
    fn f64_to_bigrational_panics_on_nan() {
        f64_to_bigrational(f64::NAN);
    }

    #[test]
    #[should_panic(expected = "non-finite f64 in exact conversion")]
    fn f64_to_bigrational_panics_on_inf() {
        f64_to_bigrational(f64::INFINITY);
    }

    #[test]
    #[should_panic(expected = "non-finite f64 in exact conversion")]
    fn f64_to_bigrational_panics_on_neg_inf() {
        f64_to_bigrational(f64::NEG_INFINITY);
    }
}
