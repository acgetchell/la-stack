#![forbid(unsafe_code)]

//! Exact arithmetic operations via arbitrary-precision rational numbers.
//!
//! This module is only compiled when the `"exact"` Cargo feature is enabled.
//! Exactness begins with the finite binary64 values already stored in
//! [`Matrix`] and [`Vector`]: each value is lifted losslessly to a rational.
//! These APIs cannot recover information rounded away before construction.
//!
//! # Architecture
//!
//! ## Determinants
//!
//! All determinant methods (`det_exact`, `det_exact_f64`,
//! `det_exact_rounded_f64`, and `det_sign_exact`) share the same integer-scaled
//! determinant core. Each proven-finite f64 entry is decomposed via
//! `decompose_proven_finite_f64` into `mantissa × 2^exponent`, then all entries
//! are scaled to a common `BigInt`
//! matrix (shifting by `e - e_min`). D≤4 uses direct integer expansions; larger
//! matrices use fraction-free Bareiss elimination \[7\] entirely in `BigInt`
//! arithmetic — no `BigRational`, no GCD, no denominator tracking. The result
//! is `(det_int, total_exp)` where `det = det_int × 2^(D × e_min)`. `det_exact`
//! wraps this with `big_int_exp_to_big_rational` to reconstruct a reduced
//! `BigRational`; `det_exact_f64` converts the same pair only when the exact
//! value is representable as finite binary64; `det_exact_rounded_f64` rounds
//! the same exact value to finite binary64; and `det_sign_exact` reads the sign
//! directly from `det_int` (the scale factor is always positive).
//!
//! `det_sign_exact` adds a two-stage adaptive-precision optimisation inspired
//! by Shewchuk's robust geometric predicates \[8\]:
//!
//! 1. **Fast filter (D ≤ 4)**: compute `det_direct()` and a conservative error
//!    bound. If `|det| > bound`, the f64 sign is provably correct — return
//!    immediately without allocating.
//! 2. **Exact fallback**: evaluate the scaled `BigInt` matrix directly for
//!    D ≤ 4 or with Bareiss elimination for D ≥ 5, yielding a
//!    guaranteed-correct sign.
//!
//! ## Linear system solve
//!
//! `solve_exact`, `solve_exact_f64`, and `solve_exact_rounded_f64` solve
//! `A x = b` with a hybrid algorithm that shares the determinant path's exact
//! integer scaling and then applies Bareiss elimination to the augmented
//! system. Matrix and RHS entries are decomposed via
//! `decompose_proven_finite_f64` into `mantissa × 2^exponent`. Each side is
//! independently scaled to its own minimum exponent before assembly into a
//! `BigInt` augmented system `(A | b)`; the resulting solution is adjusted by
//! the exact power-of-two ratio between those scales. This avoids inflating one
//! side's integers merely because the other side has much smaller entries.
//! Forward elimination runs entirely in `BigInt` with
//! fraction-free Bareiss updates \[7\] — no `BigRational`, no GCD
//! normalisation in the `O(D³)` phase.  Once the system is upper
//! triangular, back-substitution is performed in `BigRational`, where
//! fractions are inherent; this phase is only `O(D²)` so the rational
//! overhead is modest.  First-non-zero pivoting is used throughout;
//! since all arithmetic is exact, any non-zero pivot gives the correct
//! result (no numerical stability concern). Every finite `f64` is exactly
//! representable as a rational, so the result is exact for the stored inputs.
//! `solve_exact_f64` returns `Vector<D>` only when every exact component is
//! exactly representable as finite binary64; `solve_exact_rounded_f64` returns
//! the exact components rounded to finite binary64.
//!
//! ## f64 → integer decomposition
//!
//! Both the determinant and solve paths share a single conversion
//! primitive, `decompose_proven_finite_f64`, which parses the IEEE 754 binary64 bit
//! representation into a proof-bearing component (\[9\]). The
//! determinant path combines those components into a `BigInt` matrix for
//! direct expansion or Bareiss elimination and a `2^(D × e_min)` scale factor,
//! while the solve
//! path builds a `BigInt` augmented system and lifts the
//! upper-triangular result into `BigRational` for back-substitution.
//! See Goldberg \[10\] for background on floating-point representation
//! and conversion. Reference numbers refer to
//! `REFERENCES.md`.
//!
//! ## Validation
//!
//! Public `Matrix` / `Vector` values are finite by construction before exact
//! methods reach the integer-scaled exact core. The decomposition helpers consume
//! that proof without repeating stored-entry validation; a fallible raw-f64
//! decomposition remains only to test rejection at the primitive boundary.

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

/// The exact sign of a determinant.
///
/// Available with the `exact` Cargo feature.
///
/// This type makes the three possible outcomes explicit instead of exposing a
/// raw integer that could contain values other than −1, 0, or +1.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// let sign = Matrix::<2>::identity().det_sign_exact();
/// assert_eq!(sign, DeterminantSign::Positive);
/// assert_eq!(sign.as_i8(), 1);
/// ```
#[must_use]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DeterminantSign {
    /// The determinant is strictly negative.
    Negative,
    /// The determinant is exactly zero.
    Zero,
    /// The determinant is strictly positive.
    Positive,
}

impl DeterminantSign {
    /// Return the conventional numeric sign −1, 0, or +1.
    #[inline]
    #[must_use]
    pub const fn as_i8(self) -> i8 {
        match self {
            Self::Negative => -1,
            Self::Zero => 0,
            Self::Positive => 1,
        }
    }
}

/// Convert an already-computed exact result to finite binary64 output.
///
/// This extension trait is implemented for [`BigRational`] determinants and
/// `[BigRational; D]` exact solutions. It lets callers retain the exact value,
/// try the strict no-rounding contract, and recover with explicit rounding
/// without repeating determinant evaluation or linear-system elimination.
/// [`BigRational::new_raw`] values are interpreted by their mathematical
/// quotient: denominator signs and common factors do not change the result. A
/// zero denominator is rejected as [`UnrepresentableReason::NotFinite`].
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let matrix = Matrix::<2>::try_from_rows([
///     [1.0 + f64::EPSILON, 0.0],
///     [0.0, 1.0 - f64::EPSILON],
/// ])?;
/// let exact = matrix.det_exact()?;
/// let rounded = match exact.try_to_f64() {
///     Ok(value) => value,
///     Err(error) if error.requires_rounding() => exact.to_rounded_f64()?,
///     Err(error) => return Err(error),
/// };
/// assert_eq!(rounded.to_bits(), 1.0_f64.to_bits());
///
/// let system = Matrix::<1>::try_from_rows([[3.0]])?;
/// let rhs = Vector::<1>::try_new([3.0])?;
/// let exact_solution = system.solve_exact(rhs)?;
/// assert_eq!(exact_solution.try_to_f64()?.into_array(), [1.0]);
/// # Ok(())
/// # }
/// ```
pub trait ExactF64Conversion {
    /// Finite binary64 output produced by the conversion.
    type Output;

    /// Convert only when every exact value already has an exact finite
    /// binary64 representation.
    ///
    /// The candidate conversion follows IEEE 754 round-to-nearest,
    /// ties-to-even, but this strict method returns it only when no rounding is
    /// required.
    ///
    /// # Errors
    /// Returns [`LaError::Unrepresentable`] with
    /// [`UnrepresentableReason::RequiresRounding`] when finite binary64 output
    /// would require rounding, or [`UnrepresentableReason::NotFinite`] when
    /// rounding cannot produce finite output. Exact solution errors include the
    /// first failing component index.
    fn try_to_f64(&self) -> Result<Self::Output, LaError>;

    /// Round the exact value to finite binary64 output.
    ///
    /// Rounding follows IEEE 754 round-to-nearest, ties-to-even.
    ///
    /// # Errors
    /// Returns [`LaError::Unrepresentable`] with
    /// [`UnrepresentableReason::NotFinite`] when rounding cannot produce finite
    /// output. Exact solution errors include the first failing component index.
    fn to_rounded_f64(&self) -> Result<Self::Output, LaError>;
}

const F64_SIGNIFICAND_BITS: i64 = 53;
const F64_FRACTION_BITS: i64 = 52;
const F64_MIN_BINARY_EXPONENT: i64 = -1074;
const F64_MIN_NORMAL_EXPONENT: i64 = -1022;
const F64_MAX_BINARY_EXPONENT: i64 = 1023;
const F64_EXPONENT_BIAS: i64 = 1023;
const F64_FRACTION_MASK: u64 = (1u64 << 52) - 1;

/// Decompose an `f64` whose finiteness has already been proven into its IEEE
/// 754 components.
///
/// This helper is total for every bit pattern so proof-bearing callers never
/// need to recover from a second finiteness check. Its result is meaningful as
/// an exact real value only when `x` is finite.
const fn decompose_proven_finite_f64(x: f64) -> Component {
    let bits = x.to_bits();
    let biased_exp = ((bits >> 52) & 0x7FF) as i32;
    let fraction = bits & 0x000F_FFFF_FFFF_FFFF;

    // ±0.0
    if biased_exp == 0 && fraction == 0 {
        return Component::Zero;
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

    // Strip trailing zeros so the mantissa is odd. The zero bit patterns
    // returned above are the only finite values with a zero mantissa.
    let tz = mantissa.trailing_zeros();
    let Some(mantissa) = NonZeroU64::new(mantissa >> tz) else {
        return Component::Zero;
    };

    Component::NonZero {
        mantissa,
        exponent: raw_exp + tz.cast_signed(),
        is_negative: bits >> 63 != 0,
    }
}

/// Parse an arbitrary `f64` into its exact IEEE 754 components.
///
/// Returns [`Component::Zero`] for ±0.0, or [`Component::NonZero`] with a
/// non-zero mantissa where the value is exactly
/// `(-1)^is_negative × mantissa × 2^exponent` and `mantissa` is odd (trailing
/// zeros stripped).  See `REFERENCES.md` \[9-10\].
///
/// # Errors
/// Returns [`LaError::NonFinite`] if `x` is NaN or infinite.
#[cfg(test)]
const fn decompose_f64(x: f64) -> Result<Component, LaError> {
    let bits = x.to_bits();
    let biased_exp = ((bits >> 52) & 0x7FF) as i32;

    if biased_exp == 0x7FF {
        cold_path();
        return Err(LaError::non_finite_input_scalar());
    }

    Ok(decompose_proven_finite_f64(x))
}

/// Convert a [`BigInt`] × `2^exp` pair to a reduced [`BigRational`].
///
/// When `exp < 0` (denominator is `2^(-exp)`), shared factors of 2 are
/// stripped from `value` to keep the fraction in lowest terms without a
/// full GCD computation.
fn big_int_exp_to_big_rational(mut value: BigInt, mut exp: i32) -> BigRational {
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
        let remaining_abs = u64::from(exp_abs) - reduce;
        exp = negative_exponent_from_magnitude(remaining_abs);
    }

    if exp >= 0 {
        BigRational::new_raw(value << exp.cast_unsigned(), BigInt::from(1u32))
    } else {
        BigRational::new_raw(value, BigInt::from(1u32) << exp.unsigned_abs())
    }
}

/// Reconstruct a non-positive `i32` exponent from its unsigned magnitude.
///
/// The exact-conversion path can produce magnitude 2^31 for `i32::MIN`, which
/// is one greater than `i32::MAX` and therefore cannot be converted before
/// negation. Magnitudes derived from an `i32` never exceed that boundary.
#[inline]
fn negative_exponent_from_magnitude(magnitude: u64) -> i32 {
    if magnitude == u64::from(i32::MIN.unsigned_abs()) {
        return i32::MIN;
    }

    let Ok(value) = i32::try_from(magnitude) else {
        cold_path();
        unreachable!("negative exponent magnitude exceeds the i32 domain");
    };
    -value
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
    if exact.denom().sign() == Sign::NoSign {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    }

    if exact.numer().sign() == Sign::NoSign {
        return Ok(0.0);
    }

    let denominator = exact.denom();
    if denominator.sign() == Sign::Plus
        && let Some(denominator_exp) = positive_power_of_two_exponent(denominator)
        && let Ok(denominator_exp) = i32::try_from(denominator_exp)
    {
        return big_int_exp_ref_to_finite_f64(exact.numer(), -denominator_exp, index, || {
            rounded_rational_unrepresentable_reason(exact)
        });
    }

    // `BigRational::new_raw` can expose a negative denominator or uncancelled
    // common factors. Normalization is necessary before deciding whether the
    // mathematical quotient is dyadic. Canonical dyadic values take the
    // borrowed fast path above and do not clone.
    let reduced = exact.reduced();
    reduced_rational_to_finite_f64(&reduced, index)
}

/// Return `k` exactly when `value` is the positive integer `2^k`.
fn positive_power_of_two_exponent(value: &BigInt) -> Option<u64> {
    if value.sign() != Sign::Plus {
        return None;
    }

    let exponent = value.trailing_zeros()?;
    (value.bits().checked_sub(1) == Some(exponent)).then_some(exponent)
}

/// Strictly convert a reduced rational with a positive denominator.
fn reduced_rational_to_finite_f64(
    exact: &BigRational,
    index: Option<usize>,
) -> Result<f64, LaError> {
    let Some(denominator_exp) = positive_power_of_two_exponent(exact.denom()) else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            rounded_rational_unrepresentable_reason(exact),
        ));
    };
    let Ok(denominator_exp) = i32::try_from(denominator_exp) else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            rounded_rational_unrepresentable_reason(exact),
        ));
    };

    big_int_exp_ref_to_finite_f64(exact.numer(), -denominator_exp, index, || {
        rounded_rational_unrepresentable_reason(exact)
    })
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

/// Convert an exact rational result to a rounded finite `f64` using IEEE 754
/// round-to-nearest, ties-to-even.
fn exact_rational_to_rounded_f64(
    exact: &BigRational,
    index: Option<usize>,
) -> Result<f64, LaError> {
    if exact.denom().sign() == Sign::NoSign {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    }
    if exact.numer().sign() == Sign::NoSign {
        return Ok(0.0);
    }

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

impl ExactF64Conversion for BigRational {
    type Output = f64;

    #[inline]
    fn try_to_f64(&self) -> Result<Self::Output, LaError> {
        exact_rational_to_finite_f64(self, None)
    }

    #[inline]
    fn to_rounded_f64(&self) -> Result<Self::Output, LaError> {
        exact_rational_to_rounded_f64(self, None)
    }
}

impl<const D: usize> ExactF64Conversion for [BigRational; D] {
    type Output = Vector<D>;

    #[inline]
    fn try_to_f64(&self) -> Result<Self::Output, LaError> {
        let mut result = [0.0; D];
        for (index, value) in self.iter().enumerate() {
            result[index] = exact_rational_to_finite_f64(value, Some(index))?;
        }
        Vector::try_new(result)
    }

    #[inline]
    fn to_rounded_f64(&self) -> Result<Self::Output, LaError> {
        let mut result = [0.0; D];
        for (index, value) in self.iter().enumerate() {
            result[index] = exact_rational_to_rounded_f64(value, Some(index))?;
        }
        Vector::try_new(result)
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
fn shifted_magnitude_to_u64(value: &BigInt, shift: u64) -> Option<u64> {
    let word_bits = u64::from(u64::BITS);
    let word_index = usize::try_from(shift / word_bits).ok()?;
    let bit_shift = u32::try_from(shift % word_bits).ok()?;
    let mut digits = value.iter_u64_digits().skip(word_index);
    let low = digits.next()? >> bit_shift;
    if bit_shift == 0 {
        Some(low)
    } else {
        let high = digits.next().unwrap_or(0) << (u64::BITS - bit_shift);
        Some(low | high)
    }
}

/// Return whether a bit is set in the magnitude of `value`.
fn magnitude_bit_is_set(value: &BigInt, bit: u64) -> bool {
    let word_bits = u64::from(u64::BITS);
    let Ok(word_index) = usize::try_from(bit / word_bits) else {
        return false;
    };
    let bit_index = u32::try_from(bit % word_bits).unwrap_or(0);
    value
        .iter_u64_digits()
        .nth(word_index)
        .is_some_and(|word| word & (1_u64 << bit_index) != 0)
}

/// Return whether any magnitude bit below `exclusive_end` is set.
fn magnitude_has_lower_bits(value: &BigInt, exclusive_end: u64) -> bool {
    let word_bits = u64::from(u64::BITS);
    let Ok(full_words) = usize::try_from(exclusive_end / word_bits) else {
        return value.sign() != Sign::NoSign;
    };
    let partial_bits = u32::try_from(exclusive_end % word_bits).unwrap_or(0);
    let mut digits = value.iter_u64_digits();

    for _ in 0..full_words {
        if digits.next().unwrap_or(0) != 0 {
            return true;
        }
    }

    if partial_bits == 0 {
        false
    } else {
        let mask = (1_u64 << partial_bits) - 1;
        digits.next().is_some_and(|word| word & mask != 0)
    }
}

/// Right-shift a magnitude and round the retained integer to nearest-even.
fn rounded_shifted_magnitude_to_u64(value: &BigInt, shift: u64) -> Option<u64> {
    if shift > value.bits() {
        return Some(0);
    }
    let retained = shifted_magnitude_to_u64(value, shift).unwrap_or(0);
    if shift == 0 {
        return Some(retained);
    }

    let guard_bit = shift - 1;
    let increment = magnitude_bit_is_set(value, guard_bit)
        && (magnitude_has_lower_bits(value, guard_bit) || retained & 1 != 0);
    retained.checked_add(u64::from(increment))
}

/// Classify an inexact integer conversion, evaluating rounding only in the
/// maximum exponent bin where it can overflow to infinity.
#[inline]
fn inexact_big_int_reason(
    top_bit_exp: i64,
    rounded_reason: impl FnOnce() -> UnrepresentableReason,
) -> UnrepresentableReason {
    if top_bit_exp < F64_MAX_BINARY_EXPONENT {
        UnrepresentableReason::RequiresRounding
    } else {
        rounded_reason()
    }
}

/// Round a `BigInt × 2^exp` pair directly to finite binary64.
///
/// The implementation reads only the magnitude bits needed for the binary64
/// significand and rounding decision. It therefore avoids constructing a
/// potentially enormous [`BigRational`] denominator for very negative
/// exponents.
fn big_int_exp_ref_to_rounded_f64(
    value: &BigInt,
    exp: i32,
    index: Option<usize>,
) -> Result<f64, LaError> {
    if value.sign() == Sign::NoSign {
        return Ok(0.0);
    }

    let sign = if value.sign() == Sign::Minus {
        1_u64 << 63
    } else {
        0
    };
    let Ok(bit_len) = i64::try_from(value.bits()) else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    };
    let Some(mut top_bit_exp) = i64::from(exp).checked_add(bit_len - 1) else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    };
    if top_bit_exp > F64_MAX_BINARY_EXPONENT {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    }

    if top_bit_exp >= F64_MIN_NORMAL_EXPONENT {
        let mut significand = if bit_len <= F64_SIGNIFICAND_BITS {
            let Some(magnitude) = shifted_magnitude_to_u64(value, 0) else {
                cold_path();
                unreachable!("nonzero integer must expose magnitude digits");
            };
            let shift = u32::try_from(F64_SIGNIFICAND_BITS - bit_len)
                .unwrap_or_else(|_| unreachable!("normal significand shift must fit u32"));
            magnitude
                .checked_shl(shift)
                .unwrap_or_else(|| unreachable!("normal significand must fit u64"))
        } else {
            let shift = u64::try_from(bit_len - F64_SIGNIFICAND_BITS)
                .unwrap_or_else(|_| unreachable!("positive significand shift must fit u64"));
            rounded_shifted_magnitude_to_u64(value, shift)
                .unwrap_or_else(|| unreachable!("rounded binary64 significand must fit u64"))
        };

        if significand == 1_u64 << F64_SIGNIFICAND_BITS {
            significand >>= 1;
            top_bit_exp += 1;
        }
        if top_bit_exp > F64_MAX_BINARY_EXPONENT {
            cold_path();
            return Err(LaError::unrepresentable(
                index,
                UnrepresentableReason::NotFinite,
            ));
        }

        let biased_exp = u64::try_from(top_bit_exp + F64_EXPONENT_BIAS)
            .unwrap_or_else(|_| unreachable!("normal exponent must be positive"));
        return Ok(f64::from_bits(
            sign | (biased_exp << F64_FRACTION_BITS) | (significand & F64_FRACTION_MASK),
        ));
    }

    let subnormal_shift = i64::from(exp) - F64_MIN_BINARY_EXPONENT;
    let significand = if subnormal_shift >= 0 {
        let Some(magnitude) = shifted_magnitude_to_u64(value, 0) else {
            cold_path();
            unreachable!("nonzero integer must expose magnitude digits");
        };
        let shift = u32::try_from(subnormal_shift)
            .unwrap_or_else(|_| unreachable!("subnormal left shift must fit u32"));
        magnitude
            .checked_shl(shift)
            .unwrap_or_else(|| unreachable!("subnormal significand must fit u64"))
    } else {
        let shift = u64::try_from(-subnormal_shift)
            .unwrap_or_else(|_| unreachable!("subnormal right shift must fit u64"));
        rounded_shifted_magnitude_to_u64(value, shift)
            .unwrap_or_else(|| unreachable!("rounded subnormal significand must fit u64"))
    };

    if significand == 1_u64 << F64_FRACTION_BITS {
        return Ok(f64::from_bits(sign | (1_u64 << F64_FRACTION_BITS)));
    }
    Ok(f64::from_bits(sign | significand))
}

/// Borrowed core for exact integer-and-exponent conversion.
///
/// The normalized significand is read directly from the [`BigInt`] digits, so
/// successful strict conversion does not clone an already-computed exact
/// result. `rounded_reason` is evaluated only when finite output would require
/// rounding.
fn big_int_exp_ref_to_finite_f64(
    value: &BigInt,
    exp: i32,
    index: Option<usize>,
    rounded_reason: impl FnOnce() -> UnrepresentableReason,
) -> Result<f64, LaError> {
    if value.sign() == Sign::NoSign {
        return Ok(0.0);
    }

    let is_negative = value.sign() == Sign::Minus;
    let mut exp = i64::from(exp);
    let Some(trailing_zeros) = value.trailing_zeros() else {
        cold_path();
        unreachable!("nonzero integer must have a least-significant set bit");
    };
    let Ok(trailing_zeros_i64) = i64::try_from(trailing_zeros) else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    };
    let Some(updated_exp) = exp.checked_add(trailing_zeros_i64) else {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    };
    exp = updated_exp;

    let Some(bit_len) = value.bits().checked_sub(trailing_zeros) else {
        cold_path();
        unreachable!("trailing-zero count cannot exceed integer bit length");
    };
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
    if top_bit_exp > F64_MAX_BINARY_EXPONENT {
        cold_path();
        return Err(LaError::unrepresentable(
            index,
            UnrepresentableReason::NotFinite,
        ));
    }
    if exp < F64_MIN_BINARY_EXPONENT {
        cold_path();
        // A low least-significant exponent normally rounds to a finite value,
        // but a very wide integer in the maximum exponent bin can round up to
        // infinity.
        let reason = inexact_big_int_reason(top_bit_exp, rounded_reason);
        return Err(LaError::unrepresentable(index, reason));
    }
    if bit_len > F64_SIGNIFICAND_BITS {
        cold_path();
        // Rounding can overflow only when the exact value already occupies the
        // maximum binary64 exponent bin. Avoid the full rounding calculation
        // for every ordinary inexact conversion.
        let reason = inexact_big_int_reason(top_bit_exp, rounded_reason);
        return Err(LaError::unrepresentable(index, reason));
    }

    let Some(mantissa) = shifted_magnitude_to_u64(value, trailing_zeros) else {
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

fn big_int_exp_to_finite_f64(
    value: &BigInt,
    exp: i32,
    index: Option<usize>,
) -> Result<f64, LaError> {
    big_int_exp_ref_to_finite_f64(value, exp, index, || {
        match big_int_exp_ref_to_rounded_f64(value, exp, index) {
            Ok(_) => UnrepresentableReason::RequiresRounding,
            Err(_) => UnrepresentableReason::NotFinite,
        }
    })
}

/// Convert a `BigInt × 2^exp` determinant pair to a rounded finite `f64`.
fn big_int_exp_to_rounded_f64(value: &BigInt, exp: i32) -> Result<f64, LaError> {
    big_int_exp_ref_to_rounded_f64(value, exp, None)
}

// -----------------------------------------------------------------------
// Shared integer-scaling and Bareiss primitives
// -----------------------------------------------------------------------
//
// Both `exact_det_int_finite` (determinants) and `bareiss_solve_finite` (linear
// systems) parse every f64 entry into a proof-bearing component, track the
// minimum exponent across non-zero entries, and scale each entry by
// `2^(exp − e_min)`. Determinants then use direct expansions for D≤4 and
// fraction-free Bareiss elimination for D≥5; solves scale the matrix and RHS
// independently, use Bareiss elimination on the augmented system, and restore
// their exact power-of-two scale ratio after rational back-substitution.

/// Decomposed finite f64 in the form `(-1)^is_negative · mantissa · 2^exponent`.
///
/// `Zero` represents ±0.0. Non-zero entries carry a [`NonZeroU64`] mantissa, so
/// the exact-arithmetic paths cannot accidentally combine an absent mantissa
/// with active exponent/sign fields after decomposition.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
enum Component {
    #[default]
    Zero,
    NonZero {
        mantissa: NonZeroU64,
        exponent: i32,
        is_negative: bool,
    },
}

impl Component {
    /// Return the exponent carried by a non-zero component.
    const fn exponent(self) -> Option<i32> {
        match self {
            Self::Zero => None,
            Self::NonZero { exponent, .. } => Some(exponent),
        }
    }
}

mod decomposition {
    use super::Component;

    /// A component collection paired with its proven minimum non-zero exponent.
    ///
    /// `None` represents an all-zero collection; no sentinel exponent is stored.
    /// Private fields prevent callers from supplying components and their proof
    /// independently.
    #[derive(Clone, Debug, Eq, PartialEq)]
    pub(super) struct Decomposed<T> {
        components: T,
        min_exponent: Option<i32>,
    }

    impl<T> Decomposed<T> {
        /// Borrow the parsed components.
        pub(super) const fn components(&self) -> &T {
            &self.components
        }

        /// Return the minimum exponent, or `None` when every component is zero.
        pub(super) const fn min_exponent(&self) -> Option<i32> {
            self.min_exponent
        }
    }

    impl<const D: usize> Decomposed<[Component; D]> {
        /// Derive a vector decomposition and its proof together.
        pub(super) fn from_vector_components(components: [Component; D]) -> Self {
            let min_exponent = components
                .iter()
                .filter_map(|component| component.exponent())
                .min();
            Self {
                components,
                min_exponent,
            }
        }
    }

    impl<const D: usize> Decomposed<[[Component; D]; D]> {
        /// Derive a matrix decomposition and its proof together.
        pub(super) fn from_matrix_components(components: [[Component; D]; D]) -> Self {
            let min_exponent = components
                .iter()
                .flatten()
                .filter_map(|component| component.exponent())
                .min();
            Self {
                components,
                min_exponent,
            }
        }
    }

    /// A scaling exponent derived from a component collection's minimum.
    ///
    /// The private field prevents raw construction outside this proof-owning
    /// module.
    #[derive(Clone, Copy, Debug, Eq, PartialEq)]
    pub(super) struct ScaleExponent {
        value: i32,
    }

    impl ScaleExponent {
        /// Canonical scale for an empty or all-zero component collection.
        pub(super) const ZERO: Self = Self { value: 0 };

        /// Select a scale for one decomposed collection.
        pub(super) const fn for_decomposed<T>(decomposed: &Decomposed<T>) -> Self {
            let value = match decomposed.min_exponent() {
                Some(exponent) => exponent,
                None => 0,
            };
            Self { value }
        }

        /// Select the lower of two already-derived collection scales.
        pub(super) const fn min(self, other: Self) -> Self {
            if self.value < other.value {
                self
            } else {
                other
            }
        }

        /// Return the proven collection exponent.
        pub(super) const fn get(self) -> i32 {
            self.value
        }

        /// Compute a non-negative shift from this proven collection exponent.
        ///
        /// # Panics
        /// Panics only if a private decomposition invariant is broken and an entry
        /// exponent is lower than the common minimum.
        pub(super) fn shift_for(self, exponent: i32) -> u32 {
            let Some(shift) = exponent.checked_sub(self.value) else {
                unreachable!("finite f64 exponent difference cannot overflow");
            };
            let Ok(shift) = u32::try_from(shift) else {
                unreachable!("scale exponent cannot exceed a component exponent");
            };
            shift
        }
    }
}

use decomposition::{Decomposed, ScaleExponent};

/// Decompose a matrix whose finite-storage invariant has already been proven.
fn decompose_proven_finite_matrix<const D: usize>(
    m: &Matrix<D>,
) -> Decomposed<[[Component; D]; D]> {
    let components =
        from_fn(|row| from_fn(|col| decompose_proven_finite_f64(m.as_rows()[row][col])));
    Decomposed::from_matrix_components(components)
}

/// Decompose a vector whose finite-storage invariant has already been proven.
fn decompose_proven_finite_vector<const D: usize>(v: &Vector<D>) -> Decomposed<[Component; D]> {
    let components = from_fn(|index| decompose_proven_finite_f64(v.as_array()[index]));
    Decomposed::from_vector_components(components)
}

/// Convert a single decomposed component to its scaled `BigInt`
/// representation: `(±mantissa) << (exp − e_min)`.
#[inline]
fn component_to_big_int(component: Component, scale: ScaleExponent) -> BigInt {
    match component {
        Component::Zero => BigInt::from(0),
        Component::NonZero {
            mantissa,
            exponent,
            is_negative,
        } => {
            let value = BigInt::from(mantissa.get()) << scale.shift_for(exponent);
            if is_negative { -value } else { value }
        }
    }
}

/// Build a `D×D` integer matrix from components, scaled to a shared base.
fn build_big_int_matrix<const D: usize>(
    components: &[[Component; D]; D],
    scale: ScaleExponent,
) -> [[BigInt; D]; D] {
    from_fn(|row| from_fn(|col| component_to_big_int(components[row][col], scale)))
}

/// Build a length-`D` integer vector from components, scaled to a shared base.
fn build_big_int_vec<const D: usize>(
    components: &[Component; D],
    scale: ScaleExponent,
) -> [BigInt; D] {
    from_fn(|index| component_to_big_int(components[index], scale))
}

/// Compute a 2×2 determinant from a scaled integer matrix.
#[inline]
fn det2_big_int<const D: usize>(a: &[[BigInt; D]; D]) -> BigInt {
    &a[0][0] * &a[1][1] - &a[0][1] * &a[1][0]
}

/// Compute an exact 3×3 determinant from borrowed scaled-integer entries.
///
/// This fixed-shape kernel serves both the direct D=3 determinant path and the
/// 3×3 minors used by the D=4 expansion. Borrowing entries avoids cloning
/// [`BigInt`] values while keeping every operation in exact integer arithmetic.
#[inline]
fn det3_big_int_entries(a: [[&BigInt; 3]; 3]) -> BigInt {
    let m00 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    let m01 = a[1][0] * a[2][2] - a[1][2] * a[2][0];
    let m02 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    a[0][0] * m00 - a[0][1] * m01 + a[0][2] * m02
}

/// Compute a 3×3 determinant from a scaled integer matrix.
#[inline]
fn det3_big_int<const D: usize>(a: &[[BigInt; D]; D]) -> BigInt {
    det3_big_int_entries([
        [&a[0][0], &a[0][1], &a[0][2]],
        [&a[1][0], &a[1][1], &a[1][2]],
        [&a[2][0], &a[2][1], &a[2][2]],
    ])
}

/// Compute a 4×4 determinant from a scaled integer matrix.
#[inline]
fn det4_big_int<const D: usize>(a: &[[BigInt; D]; D]) -> BigInt {
    let mut det = BigInt::from(0);

    if a[0][0].sign() != Sign::NoSign {
        let c00 = det3_big_int_entries([
            [&a[1][1], &a[1][2], &a[1][3]],
            [&a[2][1], &a[2][2], &a[2][3]],
            [&a[3][1], &a[3][2], &a[3][3]],
        ]);
        det += &a[0][0] * c00;
    }
    if a[0][1].sign() != Sign::NoSign {
        let c01 = det3_big_int_entries([
            [&a[1][0], &a[1][2], &a[1][3]],
            [&a[2][0], &a[2][2], &a[2][3]],
            [&a[3][0], &a[3][2], &a[3][3]],
        ]);
        det -= &a[0][1] * c01;
    }
    if a[0][2].sign() != Sign::NoSign {
        let c02 = det3_big_int_entries([
            [&a[1][0], &a[1][1], &a[1][3]],
            [&a[2][0], &a[2][1], &a[2][3]],
            [&a[3][0], &a[3][1], &a[3][3]],
        ]);
        det += &a[0][2] * c02;
    }
    if a[0][3].sign() != Sign::NoSign {
        let c03 = det3_big_int_entries([
            [&a[1][0], &a[1][1], &a[1][2]],
            [&a[2][0], &a[2][1], &a[2][2]],
            [&a[3][0], &a[3][1], &a[3][2]],
        ]);
        det -= &a[0][3] * c03;
    }

    det
}

/// Outcome of a Bareiss forward-elimination pass.
#[derive(Debug)]
enum BareissResult {
    /// Elimination completed; `odd_swaps` records the parity of row
    /// swaps (relevant for determinants; solves discard it).
    Upper { odd_swaps: bool },
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
    let mut odd_swaps = false;

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
                    odd_swaps = !odd_swaps;
                    found = true;
                    break;
                }
            }
            if !found {
                cold_path();
                return BareissResult::Singular { pivot_col: k };
            }
        }

        // The final pivot has now been proven non-zero. There are no rows or
        // columns left to eliminate, and `prev_pivot` would never be read again.
        if k + 1 == D {
            break;
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
    #[cfg(debug_assertions)]
    for (k, row) in a.iter().enumerate() {
        assert_ne!(row[k], zero, "pivot at ({k}, {k}) must be non-zero");
        for (i, lower_row) in a.iter().enumerate().skip(k + 1) {
            assert_eq!(
                lower_row[k], zero,
                "sub-diagonal at ({i}, {k}) must be zero"
            );
        }
    }

    BareissResult::Upper { odd_swaps }
}

/// Compute the determinant scale exponent `D × e_min`.
///
/// This centralizes the scale-overflow classification used by exact
/// determinant value APIs. Sign-only evaluation deliberately bypasses this
/// bookkeeping because a positive binary scale cannot change determinant sign.
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

/// Compute the determinant integer and its shared per-entry scale.
///
/// Returns `(det_int, scale)` where the true determinant is
/// `det_int × 2^(D × scale)`. Since that scale factor is always positive,
/// callers interested only in the sign do not need to form `D × scale`.
///
/// All arithmetic is in `BigInt` — no `BigRational`, no GCD, no denominator
/// tracking.  Each f64 entry is decomposed into `mantissa × 2^exponent` and
/// scaled to a common base `2^e_min` so every entry becomes an integer. D≤4
/// uses direct determinant expansions; larger matrices use Bareiss elimination
/// whose inner-loop division is exact (guaranteed by the algorithm).
///
fn scaled_det_int_finite<const D: usize>(m: &Matrix<D>) -> (BigInt, ScaleExponent) {
    let decomposed = decompose_proven_finite_matrix(m);
    scaled_det_int_decomposed(&decomposed)
}

/// Compute a determinant integer from a proof-bearing component table.
fn scaled_det_int_decomposed<const D: usize>(
    decomposed: &Decomposed<[[Component; D]; D]>,
) -> (BigInt, ScaleExponent) {
    // D == 0 has no `a[D-1][D-1]` to read; shortcut to the empty-product
    // determinant.
    if D == 0 {
        return (BigInt::from(1), ScaleExponent::ZERO);
    }

    if decomposed.min_exponent().is_none() {
        return (BigInt::from(0), ScaleExponent::ZERO);
    }
    let scale = ScaleExponent::for_decomposed(decomposed);
    let mut a = build_big_int_matrix(decomposed.components(), scale);
    let det_int = match D {
        1 => take(&mut a[0][0]),
        2 => det2_big_int(&a),
        3 => det3_big_int(&a),
        4 => det4_big_int(&a),
        _ => {
            let odd_swaps = match bareiss_forward_eliminate(&mut a, None) {
                BareissResult::Upper { odd_swaps } => odd_swaps,
                BareissResult::Singular { .. } => {
                    cold_path();
                    return (BigInt::from(0), ScaleExponent::ZERO);
                }
            };

            let det = take(&mut a[D - 1][D - 1]);
            if odd_swaps { -det } else { det }
        }
    };

    (det_int, scale)
}

/// Compute the exact determinant as an integer plus one total binary scale.
///
/// Zero determinants use exponent zero because their value is independent of
/// scale. Non-zero determinants validate `D × e_min` for the value-producing
/// exact APIs; sign-only callers use [`scaled_det_int_finite`] directly.
fn exact_det_int_finite<const D: usize>(m: &Matrix<D>) -> Result<(BigInt, i32), LaError> {
    let (det_int, scale) = scaled_det_int_finite(m);
    if det_int.sign() == Sign::NoSign {
        return Ok((det_int, 0));
    }
    let total_exp = determinant_scale_exp::<D>(scale.get())?;
    Ok((det_int, total_exp))
}

/// Compute the exact determinant of a `D×D` matrix using direct `BigInt`
/// expansions for D≤4 or integer-only Bareiss elimination for D≥5, then return
/// the result as a `BigRational`.
fn exact_det_finite<const D: usize>(m: &Matrix<D>) -> Result<BigRational, LaError> {
    let (det_int, total_exp) = exact_det_int_finite(m)?;
    Ok(big_int_exp_to_big_rational(det_int, total_exp))
}

/// Solve `A x = b` exactly after matrix and RHS finiteness has been proven.
///
/// Public [`Matrix`] / [`Vector`] values are finite by construction before
/// reaching this helper, so decomposition can proceed without rediscovering
/// stored NaN/∞ entries.
///
/// # Errors
/// Returns [`LaError::Singular`] if the matrix is exactly singular.
fn bareiss_solve_finite<const D: usize>(
    m: &Matrix<D>,
    b: &Vector<D>,
) -> Result<[BigRational; D], LaError> {
    let matrix = decompose_proven_finite_matrix(m);
    let rhs = decompose_proven_finite_vector(b);
    bareiss_solve_components(&matrix, &rhs)
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
fn bareiss_solve_components<const D: usize>(
    matrix: &Decomposed<[[Component; D]; D]>,
    rhs: &Decomposed<[Component; D]>,
) -> Result<[BigRational; D], LaError> {
    const MAX_SHARED_SCALE_GAP_BITS: u32 = 64;

    let independent_matrix_scale = ScaleExponent::for_decomposed(matrix);
    let independent_rhs_scale = ScaleExponent::for_decomposed(rhs);
    let scale_gap = independent_matrix_scale
        .get()
        .abs_diff(independent_rhs_scale.get());
    let (matrix_scale, rhs_scale) = if scale_gap <= MAX_SHARED_SCALE_GAP_BITS {
        let shared = independent_matrix_scale.min(independent_rhs_scale);
        (shared, shared)
    } else {
        (independent_matrix_scale, independent_rhs_scale)
    };
    let mut a = build_big_int_matrix(matrix.components(), matrix_scale);
    let mut rhs = build_big_int_vec(rhs.components(), rhs_scale);

    match bareiss_forward_eliminate(&mut a, Some(&mut rhs)) {
        BareissResult::Upper { .. } => {}
        BareissResult::Singular { pivot_col } => {
            cold_path();
            return Err(LaError::singular_exact(pivot_col));
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

    let solution_scale_exp = rhs_scale
        .get()
        .checked_sub(matrix_scale.get())
        .unwrap_or_else(|| unreachable!("finite f64 scale difference cannot overflow i32"));
    if solution_scale_exp != 0 {
        let solution_scale = big_int_exp_to_big_rational(BigInt::from(1_u8), solution_scale_exp);
        for component in &mut x {
            *component *= &solution_scale;
        }
    }

    Ok(x)
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
    let (det_int, total_exp) = exact_det_int_finite(m)?;
    big_int_exp_to_finite_f64(&det_int, total_exp, None)
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
/// [`UnrepresentableReason::NotFinite`] if rounding cannot produce a finite `f64`.
#[inline]
fn det_exact_rounded_f64_finite<const D: usize>(m: &Matrix<D>) -> Result<f64, LaError> {
    let (det_int, total_exp) = exact_det_int_finite(m)?;
    big_int_exp_to_rounded_f64(&det_int, total_exp)
}

/// Exact determinant sign for an already finite matrix.
///
/// The fast `f64` filter treats overflowed or underflow-sensitive scalar
/// intermediates as inconclusive, then falls back to exact integer sign
/// computation: direct expansion for D≤4 or Bareiss elimination for D≥5.
///
#[inline]
fn det_sign_exact_finite<const D: usize>(m: &Matrix<D>) -> DeterminantSign {
    if let Ok(Some(estimate)) = m.det_direct_with_errbound() {
        let det_f64 = estimate.determinant();
        let error_bound = estimate.absolute_error_bound();
        if det_f64 > error_bound {
            return DeterminantSign::Positive;
        }
        if det_f64 < -error_bound {
            return DeterminantSign::Negative;
        }
    }

    cold_path();
    let decomposed = decompose_proven_finite_matrix(m);
    let (det_int, _) = scaled_det_int_decomposed(&decomposed);
    match det_int.sign() {
        Sign::Plus => DeterminantSign::Positive,
        Sign::Minus => DeterminantSign::Negative,
        Sign::NoSign => DeterminantSign::Zero,
    }
}

impl<const D: usize> Matrix<D> {
    /// Exact determinant using arbitrary-precision rational arithmetic.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Returns the determinant as an exact [`BigRational`] value. Every finite
    /// `f64` is exactly representable as a rational, so the conversion is
    /// lossless and the result is exact for the stored binary64 entries. It
    /// cannot recover precision lost before matrix construction.
    ///
    /// # When to use
    ///
    /// Use this when you need the exact determinant *value* — for example,
    /// volume computation over stored coordinates or distinguishing simplices
    /// that are exactly degenerate at those coordinates from near-degenerate
    /// ones. If you only need the *sign*, prefer
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
        exact_det_finite(self)
    }

    /// Exact determinant converted to `f64`.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Computes the exact determinant with the same integer-scaled core used by
    /// [`det_exact`](Self::det_exact), then converts the exact scaled integer
    /// result to `f64` only if the result is exactly representable as a finite
    /// binary64 value. The candidate conversion follows IEEE 754
    /// round-to-nearest, ties-to-even, but is returned only when no rounding is
    /// required.
    ///
    /// When callers also need the exact value or may recover with explicit
    /// rounding, compute [`det_exact`](Self::det_exact) once and use
    /// [`ExactF64Conversion`] on the returned [`BigRational`].
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
    /// Computes the exact determinant with the same integer-scaled core used by
    /// [`det_exact`](Self::det_exact), then rounds the exact value to a finite
    /// binary64 value using IEEE 754 round-to-nearest, ties-to-even. Unlike
    /// [`det_exact_f64`](Self::det_exact_f64), this method is intentionally lossy
    /// and may round non-dyadic or underflowing nonzero exact determinants.
    ///
    /// # Examples
    /// ```
    /// use core::assert_matches;
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([
    ///     [1.0 + f64::EPSILON, 0.0],
    ///     [0.0, 1.0 - f64::EPSILON],
    /// ])?;
    ///
    /// assert_matches!(
    ///     m.det_exact_f64(),
    ///     Err(LaError::Unrepresentable {
    ///         index: None,
    ///         reason: UnrepresentableReason::RequiresRounding,
    ///         ..
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
    /// Returns [`LaError::Unrepresentable`] if rounding cannot produce a finite `f64`.
    #[inline]
    pub fn det_exact_rounded_f64(&self) -> Result<f64, LaError> {
        det_exact_rounded_f64_finite(self)
    }

    /// Exact linear system solve using hybrid integer/rational arithmetic.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Solves `A x = b` where `A` is `self` and `b` is the given vector.
    /// Returns the exact solution as `[BigRational; D]`. Every finite `f64` is
    /// exactly representable as a rational, so the conversion is lossless and
    /// the result is exact for the stored binary64 entries. It cannot recover
    /// precision lost before matrix or vector construction.
    ///
    /// # When to use
    ///
    /// Use this when you need a solution exact for the stored inputs — for
    /// example, circumcenter computation over stored coordinates for
    /// near-degenerate simplices where f64 arithmetic may produce wildly wrong
    /// results.
    ///
    /// # Algorithm
    ///
    /// Matrix and RHS entries are decomposed via IEEE 754 bit extraction and
    /// independently scaled to their own power-of-two bases so both sides of
    /// the augmented system `(A | b)` become integer-valued without needless
    /// cross-side shifts. After solving that integer system, the exact
    /// power-of-two ratio between the RHS and matrix scales is restored.
    /// Forward elimination runs entirely in `BigInt`
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
        bareiss_solve_finite(self, &b)
    }

    /// Exact linear system solve converted to `f64`.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Computes the exact [`BigRational`] solution via
    /// [`solve_exact`](Self::solve_exact) and converts each component to `f64`
    /// only if that component is exactly representable as a finite binary64
    /// value. The candidate conversion follows IEEE 754 round-to-nearest,
    /// ties-to-even, but is returned only when no rounding is required.
    ///
    /// When callers also need the exact solution or may recover with explicit
    /// rounding, compute [`solve_exact`](Self::solve_exact) once and use
    /// [`ExactF64Conversion`] on the returned array.
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
        self.solve_exact(b)?.try_to_f64()
    }

    /// Exact linear system solve rounded to `f64`.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Computes the exact [`BigRational`] solution via
    /// [`solve_exact`](Self::solve_exact) and rounds each component to a finite
    /// binary64 value using IEEE 754 round-to-nearest, ties-to-even. Unlike
    /// [`solve_exact_f64`](Self::solve_exact_f64), this method is intentionally
    /// lossy and may round non-dyadic or underflowing nonzero exact components.
    ///
    /// # Examples
    /// ```
    /// use core::assert_matches;
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<1>::try_from_rows([[3.0]])?;
    /// let b = Vector::<1>::try_new([1.0])?;
    ///
    /// assert_matches!(
    ///     a.solve_exact_f64(b),
    ///     Err(LaError::Unrepresentable {
    ///         index: Some(0),
    ///         reason: UnrepresentableReason::RequiresRounding,
    ///         ..
    ///     })
    /// );
    /// assert_eq!(a.solve_exact_rounded_f64(b)?.into_array(), [1.0 / 3.0]);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if the matrix is exactly singular.
    /// Returns [`LaError::Unrepresentable`] if rounding any component cannot
    /// produce a finite `f64`.
    #[inline]
    pub fn solve_exact_rounded_f64(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        self.solve_exact(b)?.to_rounded_f64()
    }

    /// Exact determinant sign using adaptive-precision arithmetic.
    ///
    /// Requires the `exact` Cargo feature.
    ///
    /// Returns [`DeterminantSign::Positive`], [`DeterminantSign::Negative`], or
    /// [`DeterminantSign::Zero`] according to the determinant that is exact for
    /// the stored binary64 entries. This cannot recover precision lost before
    /// matrix construction.
    ///
    /// For D ≤ 4, a fast f64 filter is tried first: `det_direct()` is compared
    /// against a conservative error bound derived from the matrix permanent.
    /// If the f64 result clearly exceeds the bound, the sign is returned
    /// immediately without allocating. Otherwise, exact integer arithmetic
    /// computes the sign without constructing any `BigRational` values: direct
    /// `BigInt` expansions for D ≤ 4 and Bareiss elimination for D ≥ 5.
    ///
    /// # When to use
    ///
    /// Use this when the sign of the determinant over the stored entries must be
    /// correct regardless of floating-point conditioning (e.g. geometric
    /// predicates on near-degenerate stored coordinates). For well-conditioned
    /// matrices the fast filter resolves the sign without touching
    /// `BigRational`, so the overhead is minimal.
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
    /// assert_eq!(m.det_sign_exact(), DeterminantSign::Zero);
    ///
    /// assert_eq!(Matrix::<3>::identity().det_sign_exact(), DeterminantSign::Positive);
    /// # Ok::<(), LaError>(())
    /// ```
    #[inline]
    pub fn det_sign_exact(&self) -> DeterminantSign {
        det_sign_exact_finite(self)
    }
}

#[cfg(test)]
mod tests {
    use core::assert_matches;
    use std::array::from_fn;

    use num_traits::Signed;
    use pastey::paste;
    use proptest::prelude::*;

    use super::*;
    use crate::{
        ArithmeticOperation, DEFAULT_SINGULAR_TOL, NonFiniteLocation, NonFiniteOrigin,
        SingularityReason,
    };

    // -----------------------------------------------------------------------
    // Test helpers
    // -----------------------------------------------------------------------

    /// Build an exact `BigRational` from an `f64` via IEEE 754 bit decomposition.
    ///
    /// Thin wrapper over [`decompose_f64`] that packs the mantissa/exponent
    /// pair into a fully-formed `BigRational` of the form `±m · 2^e`.  The
    /// production code paths (`exact_det_int_finite`, `bareiss_solve_finite`) instead
    /// decompose entries into scaled `BigInt` collections, which avoids
    /// per-entry GCD work in the elimination loops — so this helper
    /// is not used by them and lives here to keep test assertions concise
    /// (e.g. `assert_eq!(x[0], f64_to_big_rational(3.0))`).
    ///
    /// See `REFERENCES.md` \[9-10\] for the IEEE 754 standard and Goldberg's
    /// survey of floating-point representation.
    ///
    /// # Panics
    /// Panics if `x` is NaN or infinite.
    fn f64_to_big_rational(x: f64) -> BigRational {
        let component = decompose_f64(x).expect("test helper requires finite f64 input");
        let Component::NonZero {
            mantissa,
            exponent,
            is_negative,
        } = component
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

    fn assert_non_finite_input_scalar<T>(result: &Result<T, LaError>) {
        let Err(error) = result else {
            panic!("expected a non-finite scalar-input error");
        };
        assert!(matches!(
            *error,
            LaError::NonFinite {
                location: NonFiniteLocation::Scalar,
                origin: NonFiniteOrigin::Input,
                ..
            }
        ));
    }

    fn assert_unrepresentable<T>(
        result: &Result<T, LaError>,
        expected_index: Option<usize>,
        expected_reason: UnrepresentableReason,
    ) {
        let Err(error) = result else {
            panic!("expected an exact-to-f64 conversion error");
        };
        assert!(matches!(
            *error,
            LaError::Unrepresentable { index, reason, .. }
                if index == expected_index && reason == expected_reason
        ));
    }

    // -----------------------------------------------------------------------
    // Macro-generated per-dimension tests (D=2..5)
    // -----------------------------------------------------------------------

    macro_rules! gen_exact_identity_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<exact_identity_paths_ $d d>]() {
                    let matrix = Matrix::<$d>::identity();
                    let one = BigRational::from_integer(BigInt::from(1));

                    assert_eq!(matrix.det_exact().unwrap(), one);
                    assert_eq!(matrix.det_exact_f64().unwrap().to_bits(), 1.0_f64.to_bits());
                    assert_eq!(
                        matrix.det_exact_rounded_f64().unwrap().to_bits(),
                        1.0_f64.to_bits()
                    );
                    assert_eq!(matrix.det_sign_exact(), DeterminantSign::Positive);
                }
            }
        };
    }

    gen_exact_identity_tests!(2);
    gen_exact_identity_tests!(3);
    gen_exact_identity_tests!(4);
    gen_exact_identity_tests!(5);

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
        assert_eq!(
            Matrix::<0>::zero().det_sign_exact(),
            DeterminantSign::Positive
        );
    }

    #[test]
    fn det_sign_exact_d1_positive() {
        let m = Matrix::<1>::try_from_rows([[42.0]]).unwrap();
        assert_eq!(m.det_sign_exact(), DeterminantSign::Positive);
    }

    #[test]
    fn det_sign_exact_d1_negative() {
        let m = Matrix::<1>::try_from_rows([[-3.5]]).unwrap();
        assert_eq!(m.det_sign_exact(), DeterminantSign::Negative);
    }

    #[test]
    fn det_sign_exact_d1_zero() {
        let m = Matrix::<1>::try_from_rows([[0.0]]).unwrap();
        assert_eq!(m.det_sign_exact(), DeterminantSign::Zero);
    }

    #[test]
    fn det_sign_exact_singular_duplicate_rows() {
        let m = Matrix::<3>::try_from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [1.0, 2.0, 3.0], // duplicate of row 0
        ])
        .unwrap();
        assert_eq!(m.det_sign_exact(), DeterminantSign::Zero);
    }

    #[test]
    fn det_sign_exact_singular_linear_combination() {
        // Row 2 = row 0 + row 1 in exact arithmetic.
        let m = Matrix::<3>::try_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [5.0, 7.0, 9.0]])
            .unwrap();
        assert_eq!(m.det_sign_exact(), DeterminantSign::Zero);
    }

    #[test]
    fn det_sign_exact_negative_det_row_swap() {
        // Swapping two rows of the identity negates the determinant.
        let m = Matrix::<3>::try_from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        assert_eq!(m.det_sign_exact(), DeterminantSign::Negative);
    }

    #[test]
    fn det_sign_exact_negative_det_known() {
        // det([[1,2],[3,4]]) = 1*4 - 2*3 = -2
        let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        assert_eq!(m.det_sign_exact(), DeterminantSign::Negative);
    }

    #[test]
    fn det_sign_exact_agrees_with_det_for_spd() {
        // SPD matrix → positive determinant.
        let m = Matrix::<3>::try_from_rows([[4.0, 2.0, 0.0], [2.0, 5.0, 1.0], [0.0, 1.0, 3.0]])
            .unwrap();
        assert_eq!(m.det_sign_exact(), DeterminantSign::Positive);
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
        assert_eq!(m.det_sign_exact(), DeterminantSign::Negative);
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
        assert_eq!(m.det_sign_exact(), DeterminantSign::Positive);
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
        assert_eq!(m.det_sign_exact(), DeterminantSign::Negative);
    }

    #[test]
    fn det_sign_exact_subnormal_entries() {
        // Subnormal f64 values should convert losslessly.
        let tiny = 5e-324_f64; // smallest positive subnormal
        assert!(tiny.is_subnormal());

        let m = Matrix::<2>::try_from_rows([[tiny, 0.0], [0.0, tiny]]).unwrap();
        // det = tiny^2 > 0
        assert_eq!(m.det_sign_exact(), DeterminantSign::Positive);
    }

    #[test]
    fn det_sign_exact_falls_back_when_subnormal_rounding_reverses_direct_sign() {
        let scale = 2.0_f64.powi(-360);
        let matrix = Matrix::<3>::try_from_rows([
            [-5.0 * scale, 3.0 * scale, 6.0 * scale],
            [0.0, -7.0 * scale, -7.0 * scale],
            [2.0 * scale, -3.0 * scale, -4.0 * scale],
        ])
        .unwrap();

        assert_eq!(
            matrix.det_direct().unwrap().unwrap().to_bits(),
            (-f64::from_bits(1)).to_bits()
        );
        assert_eq!(matrix.det_errbound(), Ok(None));
        assert!(matrix.det_exact().unwrap().is_positive());
        assert_eq!(matrix.det_sign_exact(), DeterminantSign::Positive);
    }

    #[test]
    fn det_sign_exact_falls_back_for_bit_exact_underflow_counterexample() {
        let matrix = Matrix::<3>::try_from_rows([
            [
                f64::from_bits(9_218_868_437_227_405_311),
                f64::from_bits(13_830_554_455_654_793_216),
                0.0,
            ],
            [
                f64::from_bits(6_790_500_848_393_242_208),
                f64::from_bits(2_184_621_143_747_520_227),
                f64::from_bits(2_187_555_472_467_513_745),
            ],
            [
                0.0,
                f64::from_bits(2_184_859_204_554_904_434),
                f64::from_bits(2_184_762_736_385_916_910),
            ],
        ])
        .unwrap();

        assert!(matrix.det_direct().unwrap().unwrap().is_sign_positive());
        assert_eq!(matrix.det_errbound(), Ok(None));
        assert!(matrix.det_exact().unwrap().is_negative());
        assert_eq!(matrix.det_sign_exact(), DeterminantSign::Negative);
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
        assert_eq!(m.det_sign_exact(), DeterminantSign::Negative);
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
        assert_eq!(m.det_sign_exact(), DeterminantSign::Positive);
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
    fn det_errbound_d3_non_identity() {
        // Non-identity matrix to exercise all code paths in D=3 case
        let m = Matrix::<3>::try_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]])
            .unwrap();
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

    // -----------------------------------------------------------------------
    // decompose_f64 tests
    // -----------------------------------------------------------------------

    #[test]
    fn decompose_f64_zero() {
        assert_eq!(decompose_f64(0.0), Ok(Component::Zero));
        assert_eq!(decompose_f64(-0.0), Ok(Component::Zero));
    }

    #[test]
    fn decompose_f64_one() {
        assert_eq!(
            decompose_f64(1.0),
            Ok(Component::NonZero {
                mantissa: NonZeroU64::new(1).unwrap(),
                exponent: 0,
                is_negative: false,
            })
        );
    }

    #[test]
    fn decompose_f64_negative() {
        assert_eq!(
            decompose_f64(-3.5),
            Ok(Component::NonZero {
                mantissa: NonZeroU64::new(7).unwrap(),
                exponent: -1,
                is_negative: true,
            })
        );
    }

    #[test]
    fn decompose_f64_subnormal() {
        let tiny = f64::from_bits(1);
        assert!(tiny.is_subnormal());
        assert_eq!(
            decompose_f64(tiny),
            Ok(Component::NonZero {
                mantissa: NonZeroU64::new(1).unwrap(),
                exponent: -1074,
                is_negative: false,
            })
        );
    }

    #[test]
    fn decompose_f64_normalizes_mixed_subnormal_mantissa() {
        let value = f64::from_bits(0x000C_0000_0000_0000);
        assert!(value.is_subnormal());
        assert_eq!(
            decompose_f64(value),
            Ok(Component::NonZero {
                mantissa: NonZeroU64::new(3).unwrap(),
                exponent: -1024,
                is_negative: false,
            })
        );
    }

    #[test]
    fn decompose_f64_power_of_two() {
        assert_eq!(
            decompose_f64(1024.0),
            Ok(Component::NonZero {
                mantissa: NonZeroU64::new(1).unwrap(),
                exponent: 10,
                is_negative: false,
            })
        );
    }

    #[test]
    fn decompose_f64_rejects_nan() {
        assert_non_finite_input_scalar(&decompose_f64(f64::NAN));
    }

    proptest! {
        #[test]
        fn finite_f64_round_trips_through_exact_decomposition(bits in any::<u64>()) {
            let value = f64::from_bits(bits);
            prop_assume!(value.is_finite());

            if let Ok(Component::NonZero { mantissa, .. }) = decompose_f64(value) {
                prop_assert_eq!(mantissa.get() & 1, 1);
            }

            let exact = f64_to_big_rational(value);
            let reconstructed = exact_rational_to_finite_f64(&exact, None);

            prop_assert_eq!(reconstructed, Ok(value));
        }
    }

    #[test]
    fn wide_low_exponent_value_reports_non_finite_rounded_result() {
        // (2^2099 - 1) × 2^-1075 lies just below 2^1024 and rounds to +∞.
        let value = (BigInt::from(1_u8) << 2099_u32) - BigInt::from(1_u8);
        let result = big_int_exp_to_finite_f64(&value, -1075, None);

        assert!(!result.as_ref().unwrap_err().requires_rounding());
        assert_unrepresentable(&result, None, UnrepresentableReason::NotFinite);
    }

    #[test]
    fn direct_big_int_rounding_handles_extreme_negative_exponent_without_large_denominator() {
        let positive = big_int_exp_ref_to_rounded_f64(&BigInt::from(1_u8), i32::MIN, None).unwrap();
        let negative =
            big_int_exp_ref_to_rounded_f64(&BigInt::from(-1_i8), i32::MIN, None).unwrap();

        assert_eq!(positive.to_bits(), 0.0_f64.to_bits());
        assert_eq!(negative.to_bits(), (-0.0_f64).to_bits());
    }

    proptest! {
        #[test]
        fn direct_big_int_rounding_matches_rational_oracle(
            value in any::<i128>(),
            exp in -1200_i32..=1200_i32,
        ) {
            let value = BigInt::from(value);
            let direct = big_int_exp_ref_to_rounded_f64(&value, exp, None);
            let exact = big_int_exp_to_big_rational(value, exp);
            let oracle = exact_rational_to_rounded_f64(&exact, None);

            match (direct, oracle) {
                (Ok(actual), Ok(expected)) => {
                    prop_assert_eq!(actual.to_bits(), expected.to_bits());
                }
                (Err(actual), Err(expected)) => {
                    prop_assert_eq!(actual, expected);
                }
                (actual, expected) => {
                    prop_assert_eq!(actual, expected);
                }
            }
        }
    }

    #[test]
    fn component_to_big_int_distinguishes_zero_from_nonzero_mantissa() {
        let baseline = Component::NonZero {
            mantissa: NonZeroU64::new(1).unwrap(),
            exponent: 1,
            is_negative: false,
        };
        let positive = Component::NonZero {
            mantissa: NonZeroU64::new(3).unwrap(),
            exponent: 4,
            is_negative: false,
        };
        let negative = Component::NonZero {
            mantissa: NonZeroU64::new(5).unwrap(),
            exponent: 3,
            is_negative: true,
        };

        let decomposed =
            Decomposed::from_vector_components([Component::Zero, baseline, positive, negative]);
        let scale = ScaleExponent::for_decomposed(&decomposed);

        assert_eq!(
            component_to_big_int(Component::Zero, scale),
            BigInt::from(0)
        );
        assert_eq!(component_to_big_int(positive, scale), BigInt::from(24));
        assert_eq!(component_to_big_int(negative, scale), BigInt::from(-20));
    }

    #[test]
    fn decomposed_all_zero_uses_no_sentinel_exponent() {
        let decomposed = decompose_proven_finite_matrix(&Matrix::<2>::zero());
        assert_eq!(decomposed.min_exponent(), None);

        let scale = ScaleExponent::for_decomposed(&decomposed);
        assert_eq!(scale, ScaleExponent::ZERO);
        assert_eq!(scale.get(), 0);
        assert_eq!(
            build_big_int_matrix(decomposed.components(), scale),
            [
                [BigInt::from(0), BigInt::from(0)],
                [BigInt::from(0), BigInt::from(0)]
            ]
        );
    }

    #[test]
    fn matrix_and_rhs_scales_are_derived_independently() {
        let tiny = f64::from_bits(1);
        let matrix = Matrix::<2>::try_from_rows([[f64::MAX, 0.0], [0.0, 1.0]]).unwrap();
        let rhs = Vector::<2>::try_new([tiny, 0.0]).unwrap();
        let matrix = decompose_proven_finite_matrix(&matrix);
        let rhs = decompose_proven_finite_vector(&rhs);

        assert_eq!(matrix.min_exponent(), Some(0));
        assert_eq!(rhs.min_exponent(), Some(-1074));

        let matrix_scale = ScaleExponent::for_decomposed(&matrix);
        let rhs_scale = ScaleExponent::for_decomposed(&rhs);
        assert_eq!(matrix_scale.get(), 0);
        assert_eq!(matrix_scale.shift_for(0), 0);
        assert_eq!(rhs_scale.get(), -1074);
        assert_eq!(rhs_scale.shift_for(-1074), 0);
    }

    proptest! {
        #[test]
        fn derived_scale_yields_nonnegative_shifts(bits in any::<[u64; 4]>()) {
            let values = bits.map(f64::from_bits);
            prop_assume!(values.iter().all(|value| value.is_finite()));
            let matrix = Matrix::<2>::try_from_rows([
                [values[0], values[1]],
                [values[2], values[3]],
            ]).unwrap();
            let decomposed = decompose_proven_finite_matrix(&matrix);
            let scale = ScaleExponent::for_decomposed(&decomposed);

            for component in decomposed.components().iter().flatten() {
                if let Some(exponent) = component.exponent() {
                    prop_assert!(exponent >= scale.get());
                    prop_assert_eq!(
                        scale.shift_for(exponent),
                        u32::try_from(exponent - scale.get()).unwrap(),
                    );
                }
            }
        }
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

    #[test]
    fn negative_exponent_from_magnitude_covers_i32_domain_boundaries() {
        assert_eq!(negative_exponent_from_magnitude(0), 0);
        assert_eq!(negative_exponent_from_magnitude(1), -1);
        assert_eq!(
            negative_exponent_from_magnitude(i32::MAX.cast_unsigned().into()),
            -i32::MAX
        );
        assert_eq!(
            negative_exponent_from_magnitude(i32::MIN.unsigned_abs().into()),
            i32::MIN
        );
    }

    #[test]
    #[should_panic(expected = "negative exponent magnitude exceeds the i32 domain")]
    fn negative_exponent_from_magnitude_rejects_values_above_i32_domain() {
        let _ = negative_exponent_from_magnitude(u64::from(i32::MIN.unsigned_abs()) + 1);
    }

    // -----------------------------------------------------------------------
    // Exact scaled-integer determinant tests
    // -----------------------------------------------------------------------

    #[test]
    fn exact_det_int_d0() {
        let m = Matrix::<0>::zero();
        let (det, exp) = exact_det_int_finite(&m).unwrap();
        assert_eq!(det, BigInt::from(1));
        assert_eq!(exp, 0);
    }

    /// Table-driven coverage of the D=1 fast-path: each 1×1 matrix
    /// decomposes to `(±mant, exp)` directly.  Includes an integer, zero,
    /// a negative fractional, and a positive fractional case — the
    /// combinations that exercise the sign handling, the all-zero early
    /// return, trailing-zero stripping, and negative exponent scaling.
    #[test]
    fn exact_det_int_d1_cases() {
        let cases: &[(f64, i64, i32)] = &[
            // (input, expected_det_int, expected_exp)
            (7.0, 7, 0),    // integer → (7, 0)
            (0.0, 0, 0),    // all-zero early return → (0, 0)
            (-3.5, -7, -1), // -3.5 = -7 × 2^(-1)
            (0.5, 1, -1),   // 0.5  =  1 × 2^(-1)
        ];
        for &(input, expected_det_int, expected_exp) in cases {
            let m = Matrix::<1>::try_from_rows([[input]]).unwrap();
            let (det, exp) = exact_det_int_finite(&m).unwrap();
            assert_eq!(
                det,
                BigInt::from(expected_det_int),
                "det_int for input={input}"
            );
            assert_eq!(exp, expected_exp, "exp for input={input}");
        }
    }

    #[test]
    fn exact_det_int_d2_known() {
        // det([[1,2],[3,4]]) = -2
        let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        let (det_int, total_exp) = exact_det_int_finite(&m).unwrap();
        // Reconstruct and verify.
        let det = big_int_exp_to_big_rational(det_int, total_exp);
        assert_eq!(det, BigRational::from_integer(BigInt::from(-2)));
    }

    #[test]
    fn exact_det_int_all_zeros() {
        let m = Matrix::<3>::zero();
        let (det, _) = exact_det_int_finite(&m).unwrap();
        assert_eq!(det, BigInt::from(0));
    }

    #[test]
    fn exact_det_int_fractional_entries() {
        // Entries with negative exponents: 0.5 = 1×2^(-1), 0.25 = 1×2^(-2).
        // det([[0.5, 0.25], [1.0, 1.0]]) = 0.5×1.0 − 0.25×1.0 = 0.25
        let m = Matrix::<2>::try_from_rows([[0.5, 0.25], [1.0, 1.0]]).unwrap();
        let (det_int, total_exp) = exact_det_int_finite(&m).unwrap();
        let det = big_int_exp_to_big_rational(det_int, total_exp);
        assert_eq!(det, BigRational::new(BigInt::from(1), BigInt::from(4)));
    }

    #[test]
    fn exact_det_int_d3_direct_expansion_handles_zero_diagonal() {
        // A zero diagonal entry does not require pivoting in the direct D=3 expansion.
        let m = Matrix::<3>::try_from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let (det_int, total_exp) = exact_det_int_finite(&m).unwrap();
        let det = big_int_exp_to_big_rational(det_int, total_exp);
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    // -----------------------------------------------------------------------
    // big_int_exp_to_big_rational tests
    // -----------------------------------------------------------------------

    #[test]
    fn big_int_exp_to_big_rational_zero() {
        let r = big_int_exp_to_big_rational(BigInt::from(0), -50);
        assert_eq!(r, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn big_int_exp_to_big_rational_positive_exp() {
        // 3 × 2^2 = 12
        let r = big_int_exp_to_big_rational(BigInt::from(3), 2);
        assert_eq!(r, BigRational::from_integer(BigInt::from(12)));
    }

    #[test]
    fn big_int_exp_to_big_rational_negative_exp_reduced() {
        // 6 × 2^(-2) = 6/4 → reduced to 3/2 (strip one shared factor of 2)
        let r = big_int_exp_to_big_rational(BigInt::from(6), -2);
        assert_eq!(*r.numer(), BigInt::from(3));
        assert_eq!(*r.denom(), BigInt::from(2));
    }

    #[test]
    fn big_int_exp_to_big_rational_negative_exp_reduces_to_integer() {
        // 8 × 2^(-3) = 1 after stripping every denominator factor.
        let r = big_int_exp_to_big_rational(BigInt::from(8), -3);
        assert_eq!(r, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn big_int_exp_to_big_rational_negative_exp_already_odd() {
        // 3 × 2^(-2) = 3/4 (already in lowest terms since 3 is odd)
        let r = big_int_exp_to_big_rational(BigInt::from(3), -2);
        assert_eq!(*r.numer(), BigInt::from(3));
        assert_eq!(*r.denom(), BigInt::from(4));
    }

    #[test]
    fn big_int_exp_to_big_rational_negative_value() {
        // -5 × 2^1 = -10
        let r = big_int_exp_to_big_rational(BigInt::from(-5), 1);
        assert_eq!(r, BigRational::from_integer(BigInt::from(-10)));
    }

    #[test]
    fn big_int_exp_to_big_rational_negative_value_with_denominator() {
        // -3 × 2^(-2) = -3/4
        let r = big_int_exp_to_big_rational(BigInt::from(-3), -2);
        assert_eq!(*r.numer(), BigInt::from(-3));
        assert_eq!(*r.denom(), BigInt::from(4));
    }

    // -----------------------------------------------------------------------
    // Public exact determinant wrapper tests
    // -----------------------------------------------------------------------

    #[test]
    fn det_exact_d1_returns_entry() {
        let det = Matrix::<1>::try_from_rows([[7.0]])
            .unwrap()
            .det_exact()
            .unwrap();
        assert_eq!(det, f64_to_big_rational(7.0));
    }

    #[test]
    fn det_exact_d3_direct_expansion_handles_zero_diagonal() {
        // Direct D=3 expansion handles a zero diagonal entry without pivoting.
        let m = Matrix::<3>::try_from_rows([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let det = m.det_exact().unwrap();
        // det of this permutation matrix = -1
        assert_eq!(det, BigRational::from_integer(BigInt::from(-1)));
    }

    #[test]
    fn det_exact_d3_singular_zero_column_returns_zero() {
        // A zero column makes the direct D=3 determinant exactly zero.
        let m = Matrix::<3>::try_from_rows([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let det = m.det_exact().unwrap();
        assert_eq!(det, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn det_sign_exact_overflow_determinant_finite_entries() {
        // Entries near f64::MAX are finite, but the f64 determinant overflows
        // to infinity. The fast filter is inconclusive and the direct `BigInt`
        // expansion computes the correct positive sign.
        let big = f64::MAX / 2.0;
        assert!(big.is_finite());
        let m = Matrix::<3>::try_from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]])
            .unwrap();
        // det = big^2 > 0
        assert_eq!(m.det_sign_exact(), DeterminantSign::Positive);
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
        assert_unrepresentable(&m.det_exact_f64(), None, UnrepresentableReason::NotFinite);
    }

    #[test]
    fn det_exact_rounded_f64_overflow_returns_err() {
        let big = f64::MAX / 2.0;
        let m = Matrix::<3>::try_from_rows([[0.0, 0.0, 1.0], [big, 0.0, 1.0], [0.0, big, 1.0]])
            .unwrap();

        assert_unrepresentable(
            &m.det_exact_rounded_f64(),
            None,
            UnrepresentableReason::NotFinite,
        );
    }

    #[test]
    fn det_exact_f64_underflow_returns_err_for_nonzero_exact_result() {
        let tiny = f64::from_bits(1);
        let m = Matrix::<2>::try_from_rows([[tiny, 0.0], [0.0, tiny]]).unwrap();

        assert!(m.det_exact().unwrap().is_positive());
        assert_unrepresentable(
            &m.det_exact_f64(),
            None,
            UnrepresentableReason::RequiresRounding,
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
        assert_unrepresentable(
            &m.det_exact_f64(),
            None,
            UnrepresentableReason::RequiresRounding,
        );
    }

    #[test]
    fn det_exact_f64_accepts_max_finite_binary64() {
        let m = Matrix::<1>::try_from_rows([[f64::MAX]]).unwrap();

        assert_eq!(m.det_exact_f64().unwrap().to_bits(), f64::MAX.to_bits());
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
                fn [<solve_exact_identity_paths_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let b = arbitrary_rhs::<$d>();
                    let exact = a.solve_exact(b).unwrap();
                    let strict_f64 = a.solve_exact_f64(b).unwrap().into_array();

                    for i in 0..$d {
                        assert_eq!(exact[i], f64_to_big_rational(b.as_array()[i]));
                        assert_eq!(strict_f64[i].to_bits(), b.as_array()[i].to_bits());
                    }
                }

                #[test]
                fn [<solve_exact_singular_ $d d>]() {
                    // Zero matrix is singular.
                    let a = Matrix::<$d>::zero();
                    let b = arbitrary_rhs::<$d>();
                    assert_matches!(
                        a.solve_exact(b),
                        Err(LaError::Singular {
                            pivot_col: 0,
                            reason: SingularityReason::Exact,
                            ..
                        })
                    );
                }
            }
        };
    }

    gen_solve_exact_tests!(2);
    gen_solve_exact_tests!(3);
    gen_solve_exact_tests!(4);
    gen_solve_exact_tests!(5);

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
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "dimensions and indices are at most five and exactly representable as f64"
                )]
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
                        assert_eq!(x[i], f64_to_big_rational(x0[i]));
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
        assert_eq!(x[0], f64_to_big_rational(3.0));
        assert_eq!(x[1], f64_to_big_rational(2.0));
        assert_eq!(x[2], f64_to_big_rational(4.0));
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
        assert_matches!(
            a.solve_exact(b),
            Err(LaError::Singular {
                reason: SingularityReason::Exact,
                ..
            })
        );
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
        assert_eq!(x[0], f64_to_big_rational(20.0));
        assert_eq!(x[1], f64_to_big_rational(10.0));
        assert_eq!(x[2], f64_to_big_rational(30.0));
        assert_eq!(x[3], f64_to_big_rational(40.0));
        assert_eq!(x[4], f64_to_big_rational(50.0));
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
    /// `f64::MIN_POSITIVE` up through `1e100`). This exercises each
    /// collection's independently derived minimum exponent: even the largest
    /// within-collection shift remains a representable `BigInt`. The D×D case
    /// alternates `huge`/`tiny` along the diagonal with a matching RHS, giving
    /// `x = [1, …, 1]`.
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

    #[test]
    fn solve_exact_restores_independent_matrix_and_rhs_scales() {
        let large = 2.0_f64.powi(500);
        let tiny = 2.0_f64.powi(-1000);

        let large_matrix = Matrix::<2>::try_from_rows([[large, 0.0], [0.0, large]]).unwrap();
        let tiny_rhs = Vector::<2>::new([tiny, -2.0 * tiny]);
        let small_solution = large_matrix.solve_exact(tiny_rhs).unwrap();
        assert_eq!(
            small_solution[0],
            BigRational::new(BigInt::from(1_u8), BigInt::from(1_u8) << 1500_u32)
        );
        assert_eq!(
            small_solution[1],
            BigRational::new(BigInt::from(-1_i8), BigInt::from(1_u8) << 1499_u32)
        );

        let tiny_matrix = Matrix::<1>::try_from_rows([[tiny]]).unwrap();
        let large_rhs = Vector::<1>::new([large]);
        let large_solution = tiny_matrix.solve_exact(large_rhs).unwrap();
        assert_eq!(
            large_solution[0],
            BigRational::from_integer(BigInt::from(1_u8) << 1500_u32)
        );
    }

    /// Subnormal RHS entries must survive the decomposition and
    /// back-substitution paths unchanged.  The D×D case uses the identity
    /// matrix and RHS `[1·tiny, 2·tiny, …, D·tiny]`; each entry remains a
    /// valid subnormal f64 (integer multiples of `2^-1074` fit in the
    /// 52-bit subnormal mantissa for the small integers used here).
    macro_rules! gen_solve_exact_subnormal_rhs_tests {
        ($d:literal) => {
            paste! {
                #[test]
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "indices are at most five and exactly representable as f64"
                )]
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
                        assert_eq!(x[i], f64_to_big_rational((i + 1) as f64 * tiny));
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
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "indices and test offsets are small integers exactly representable as f64"
                )]
                fn [<solve_exact_pivot_swap_with_fractional_result_ $d d>]() {
                    // Top-left 2×2: A = [[0, 1], [2, 1]].  After swap:
                    // [[2, 1], [0, 1]], rhs = [4, 3] → x[1] = 3, x[0] = 1/2.
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][1] = 1.0;
                    rows[1][0] = 2.0;
                    rows[1][1] = 1.0;
                    // Identity padding for the remaining rows.
                    for (i, row) in rows.iter_mut().enumerate().skip(2) {
                        row[i] = 1.0;
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    // b = [3, 4, 12, 13, …]; padded entries are arbitrary
                    // finite integers so the identity block gives x[i] = b[i].
                    let mut b_arr = [0.0f64; $d];
                    b_arr[0] = 3.0;
                    b_arr[1] = 4.0;
                    for (i, value) in b_arr.iter_mut().enumerate().skip(2) {
                        *value = (i + 10) as f64;
                    }
                    let b = Vector::<$d>::new(b_arr);
                    let x = a.solve_exact(b).unwrap();
                    assert_eq!(x[0], BigRational::new(BigInt::from(1), BigInt::from(2)));
                    assert_eq!(x[1], BigRational::from_integer(BigInt::from(3)));
                    for (i, value) in x.iter().enumerate().skip(2) {
                        assert_eq!(value, &f64_to_big_rational((i + 10) as f64));
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
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "indices and test offsets are small integers exactly representable as f64"
                )]
                fn [<solve_exact_mid_pivot_swap_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = 1.0; rows[0][1] = 2.0; rows[0][2] = 3.0;
                    // rows[1][0..2] are zero; rows[1][2] = 4.
                    rows[1][2] = 4.0;
                    rows[2][1] = 5.0; rows[2][2] = 6.0;
                    // Identity padding for the remaining rows.
                    for (i, row) in rows.iter_mut().enumerate().skip(3) {
                        row[i] = 1.0;
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let mut b_arr = [0.0f64; $d];
                    b_arr[0] = 6.0;
                    b_arr[1] = 7.0;
                    b_arr[2] = 8.0;
                    for (i, value) in b_arr.iter_mut().enumerate().skip(3) {
                        *value = (i + 10) as f64;
                    }
                    let b = Vector::<$d>::new(b_arr);
                    let x = a.solve_exact(b).unwrap();
                    // x[0..3] = [7/4, -1/2, 7/4].
                    assert_eq!(x[0], BigRational::new(BigInt::from(7), BigInt::from(4)));
                    assert_eq!(x[1], BigRational::new(BigInt::from(-1), BigInt::from(2)));
                    assert_eq!(x[2], BigRational::new(BigInt::from(7), BigInt::from(4)));
                    for (i, value) in x.iter().enumerate().skip(3) {
                        assert_eq!(value, &f64_to_big_rational((i + 10) as f64));
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
    /// exact singularity at `pivot_col = D - 1`.
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
                    assert_matches!(
                        a.solve_exact(b),
                        Err(LaError::Singular {
                            pivot_col,
                            reason: SingularityReason::Exact,
                            ..
                        }) if pivot_col == $d - 1
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

    /// Multiply `A · x` entirely in `BigRational`, using `f64_to_big_rational`
    /// to lift each matrix entry.  Used by residual assertions for inputs
    /// whose exact solution has no closed form we can easily type out.
    fn big_rational_matvec<const D: usize>(
        a: &Matrix<D>,
        x: &[BigRational; D],
    ) -> [BigRational; D] {
        from_fn(|i| {
            let mut sum = BigRational::from_integer(BigInt::from(0));
            for (aij, xj) in a.as_rows()[i].iter().zip(x.iter()) {
                sum += f64_to_big_rational(*aij) * xj;
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
    /// treats that as inconclusive, and the direct `BigInt` fallback resolves
    /// the positive sign correctly. `det_exact_f64` must report `Unrepresentable`.
    #[test]
    fn det_sign_exact_large_entries_3x3_positive() {
        let big = f64::MAX / 2.0;
        let a = Matrix::<3>::try_from_rows([[big, 1.0, 1.0], [1.0, big, 1.0], [1.0, 1.0, big]])
            .unwrap();
        // Fast filter is inconclusive (big^3 overflows f64 to +∞), so
        // this exercises the direct `BigInt` cold path.
        assert_matches!(
            a.det_direct(),
            Err(LaError::NonFinite {
                location: NonFiniteLocation::Scalar,
                origin: NonFiniteOrigin::Computation {
                    operation: ArithmeticOperation::Determinant,
                    ..
                },
                ..
            })
        );
        assert_eq!(a.det_sign_exact(), DeterminantSign::Positive);
        // Cross-validate: the exact `BigRational` determinant must agree
        // on sign with `det_sign_exact`, and `det_exact_f64` must reject the
        // conversion (the value is representable in BigRational but far exceeds f64).
        assert!(a.det_exact().unwrap().is_positive());
        assert_unrepresentable(&a.det_exact_f64(), None, UnrepresentableReason::NotFinite);
    }

    /// Hilbert matrices are symmetric positive-definite, so
    /// `det_sign_exact` must return [`DeterminantSign::Positive`] for every D.
    /// For D=2..=4 the
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
                    assert_eq!(h.det_sign_exact(), DeterminantSign::Positive);
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
                    let ax = big_rational_matvec(&h, &x);
                    for i in 0..$d {
                        assert_eq!(ax[i], f64_to_big_rational(b_arr[i]));
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
        assert_unrepresentable(
            &a.solve_exact_f64(b),
            Some(0),
            UnrepresentableReason::NotFinite,
        );
    }

    #[test]
    fn solve_exact_f64_huge_non_dyadic_component_returns_not_finite() {
        let a = Matrix::<1>::try_from_rows([[3.0 * f64::MIN_POSITIVE]]).unwrap();
        let b = Vector::<1>::new([f64::MAX]);

        assert_unrepresentable(
            &a.solve_exact_f64(b),
            Some(0),
            UnrepresentableReason::NotFinite,
        );
    }

    #[test]
    fn solve_exact_rounded_f64_overflow_returns_err() {
        let big = f64::MAX / 2.0;
        let a = Matrix::<2>::try_from_rows([[1.0 / big, 0.0], [0.0, 1.0 / big]]).unwrap();
        let b = Vector::<2>::new([big, big]);

        assert_unrepresentable(
            &a.solve_exact_rounded_f64(b),
            Some(0),
            UnrepresentableReason::NotFinite,
        );
    }

    #[test]
    fn solve_exact_f64_underflow_returns_err_for_nonzero_exact_component() {
        let tiny = f64::from_bits(1);
        let a = Matrix::<1>::try_from_rows([[2.0]]).unwrap();
        let b = Vector::<1>::new([tiny]);

        assert_unrepresentable(
            &a.solve_exact_f64(b),
            Some(0),
            UnrepresentableReason::RequiresRounding,
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

    // -----------------------------------------------------------------------
    // exact solve boundary tests
    // -----------------------------------------------------------------------

    #[test]
    fn bareiss_solve_d1() {
        let a = Matrix::<1>::try_from_rows([[2.0]]).unwrap();
        let b = Vector::<1>::new([6.0]);
        let x = a.solve_exact(b).unwrap();
        assert_eq!(x[0], f64_to_big_rational(3.0));
    }

    #[test]
    fn bareiss_solve_singular_column_all_zero() {
        let a = Matrix::<3>::try_from_rows([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let b = Vector::<3>::new([1.0, 2.0, 3.0]);
        assert_matches!(
            a.solve_exact(b),
            Err(LaError::Singular {
                pivot_col: 1,
                reason: SingularityReason::Exact,
                ..
            })
        );
    }

    // -----------------------------------------------------------------------
    // f64_to_big_rational tests
    // -----------------------------------------------------------------------

    #[test]
    fn f64_to_big_rational_positive_zero() {
        let r = f64_to_big_rational(0.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn f64_to_big_rational_negative_zero() {
        let r = f64_to_big_rational(-0.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(0)));
    }

    #[test]
    fn f64_to_big_rational_one() {
        let r = f64_to_big_rational(1.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(1)));
    }

    #[test]
    fn f64_to_big_rational_negative_one() {
        let r = f64_to_big_rational(-1.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(-1)));
    }

    #[test]
    fn f64_to_big_rational_half() {
        let r = f64_to_big_rational(0.5);
        assert_eq!(r, BigRational::new(BigInt::from(1), BigInt::from(2)));
    }

    #[test]
    fn f64_to_big_rational_quarter() {
        let r = f64_to_big_rational(0.25);
        assert_eq!(r, BigRational::new(BigInt::from(1), BigInt::from(4)));
    }

    #[test]
    fn f64_to_big_rational_negative_three_and_a_half() {
        // -3.5 = -7/2
        let r = f64_to_big_rational(-3.5);
        assert_eq!(r, BigRational::new(BigInt::from(-7), BigInt::from(2)));
    }

    #[test]
    fn f64_to_big_rational_integer() {
        let r = f64_to_big_rational(42.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(42)));
    }

    #[test]
    fn f64_to_big_rational_power_of_two() {
        let r = f64_to_big_rational(1024.0);
        assert_eq!(r, BigRational::from_integer(BigInt::from(1024)));
    }

    #[test]
    fn f64_to_big_rational_subnormal() {
        let tiny = 5e-324_f64; // smallest positive subnormal
        assert!(tiny.is_subnormal());
        let r = f64_to_big_rational(tiny);
        // 5e-324 = 1 × 2^(-1074)
        assert_eq!(
            r,
            BigRational::new(BigInt::from(1), BigInt::from(1u32) << 1074u32)
        );
    }

    #[test]
    fn f64_to_big_rational_already_lowest_terms() {
        // 0.5 should produce numer=1, denom=2 (already reduced).
        let r = f64_to_big_rational(0.5);
        assert_eq!(*r.numer(), BigInt::from(1));
        assert_eq!(*r.denom(), BigInt::from(2));
    }

    #[test]
    fn f64_to_big_rational_round_trip() {
        // -0.0 is excluded: it maps to BigRational(0) which round-trips
        // to +0.0 (correct; tested separately in f64_to_big_rational_negative_zero).
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
            let r = f64_to_big_rational(v);
            let back = r.to_f64().expect("round-trip to_f64 failed");
            assert_eq!(
                v.to_bits(),
                back.to_bits(),
                "round-trip failed for {v}: got {back}"
            );
        }
    }

    #[test]
    fn decompose_f64_rejects_non_finite_inputs() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_non_finite_input_scalar(&decompose_f64(value));
        }
    }
}
