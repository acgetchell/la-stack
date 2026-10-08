#![forbid(unsafe_code)]

//! Scale-safe unsigned angles from compensated exterior products and dot products.

use crate::rounding::two_sum_error;
use crate::{ArithmeticOperation, LaError, NonFiniteLocation, NonFiniteOrigin, VectorOperand};

/// Unsigned angle between two finite, nonzero vectors, in radians in `[0, π]`.
///
/// Borrows equal-length coordinate slices without allocating. Every positive
/// ambient length is supported, including runtime `D + 1` spherical coordinates;
/// no matrix dispatch limit applies. Magnitudes may differ arbitrarily and the
/// input norms need not be representable. Signed zeros are treated as zero.
/// Identical inputs return positive zero, and opposite inputs return `π`.
///
/// This is a rounded numerical operation, not an exact parallelism predicate or
/// a correctly rounded result. No certified absolute error bound is provided.
/// An `atan2` of the exterior-product norm and dot product avoids the endpoint
/// sensitivity of `acos`. Power-of-two scaling and compensated products preserve
/// small differences without first rounding normalized directions. Work is
/// quadratic in the ambient length, with constant auxiliary storage; the intended
/// performance scope is small dimensions. An angle below the binary64 output
/// range can round to zero. See the
/// [derivation and limitations](https://github.com/acgetchell/la-stack/blob/main/docs/mathematical_basis.md#unsigned-vector-angles)
/// and `REFERENCES.md` \[18–19\].
///
/// For finite fixed-size storage, use [`crate::Vector::angle`], which shares the
/// numerical implementation without repeating coordinate validation.
///
/// # Errors
/// Validation proceeds in this order:
/// - [`LaError::DimensionMismatch`] for unequal lengths (including one empty).
/// - [`LaError::NonFinite`] for the first non-finite coordinate, checking the
///   left operand before the right; its input location includes operand/index.
/// - [`LaError::EmptyVector`] for two empty inputs.
/// - [`LaError::ZeroVector`] for an all-zero operand, left before right.
///
/// A non-finite computed result is reported as [`LaError::NonFinite`] with
/// [`ArithmeticOperation::VectorAngle`] computation provenance.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let tiny = f64::from_bits(1);
/// assert_eq!(angle_between(&[1.0, 0.0, 0.0], &[1.0, tiny, 0.0])?, tiny);
/// # Ok(())
/// # }
/// ```
#[inline]
pub fn angle_between(left: &[f64], right: &[f64]) -> Result<f64, LaError> {
    if left.len() != right.len() {
        return Err(LaError::DimensionMismatch {
            left: left.len(),
            right: right.len(),
        });
    }
    let left_scale = input_scale(left, VectorOperand::Left)?;
    let right_scale = input_scale(right, VectorOperand::Right)?;
    angle_with_scales(left, right, left_scale, right_scale)
}

fn input_scale(values: &[f64], operand: VectorOperand) -> Result<f64, LaError> {
    let mut scale = 0.0_f64;
    for (index, &value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(LaError::NonFinite {
                location: NonFiniteLocation::VectorOperandEntry { operand, index },
                origin: NonFiniteOrigin::Input,
            });
        }
        scale = scale.max(value.abs());
    }
    Ok(scale)
}

/// Fixed-size callers already carry finite-coordinate and equal-length proofs.
#[inline]
pub(crate) fn angle_finite<const D: usize>(
    left: &[f64; D],
    right: &[f64; D],
) -> Result<f64, LaError> {
    let maximum = |values: &[f64; D]| values.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    angle_with_scales(left, right, maximum(left), maximum(right))
}

/// Scale a direction so its largest magnitude is in `[2^500, 2^501)`.
///
/// Direct `value / maximum` can discard a half-subnormal component whose
/// difference from the other direction contributes a representable full angle.
/// Two power-of-two multiplications cover maxima down to `2^-1074`, preserving
/// significands unless a component far below the angular output range underflows.
/// Underflow here affects only components below roughly `2^-1574` relative to
/// the maximum, far below the output range even after a slice-length reduction.
struct DirectionScale {
    first: f64,
    second: f64,
}

impl DirectionScale {
    /// The caller has proved `maximum` finite and strictly positive.
    fn new(mut maximum: f64) -> Self {
        let adjustment = if maximum < f64::MIN_POSITIVE {
            maximum *= f64::from_bits((1023 + 52) << 52);
            52
        } else {
            0
        };
        let bits = maximum.to_bits();
        let exponent = (bits >> 52).cast_signed() - 1023 - adjustment;
        let shift = 500 - exponent;
        let first_shift = shift.min(1023);
        // first_shift is in [-523, 1023]; the remaining shift is in [0, 551].
        Self {
            first: f64::from_bits((first_shift + 1023).cast_unsigned() << 52),
            second: f64::from_bits((shift - first_shift + 1023).cast_unsigned() << 52),
        }
    }

    #[inline]
    fn apply(&self, value: f64) -> f64 {
        (value * self.first) * self.second
    }
}

/// Compensated `a*b - c*d`, symmetric under exchanging the two products.
///
/// Both products are bounded by `2^1002` after direction scaling. FMA captures
/// their residuals (reference [18]); `FastTwoSum` captures subtraction roundoff
/// (reference [17]). Relevant tiny minors remain normal at this enlarged scale.
/// This is an approximate minor, not an exact determinant predicate.
#[inline]
fn difference_of_products(a: f64, b: f64, c: f64, d: f64) -> f64 {
    let left = a * b;
    let right = c * d;
    let left_error = a.mul_add(b, -left);
    let right_error = c.mul_add(d, -right);
    let difference = left - right;
    let difference_error = two_sum_error(left, -right, difference);
    difference + ((left_error - right_error) + difference_error)
}

/// Finite equal-length inputs and their maximum magnitudes share one kernel.
fn angle_with_scales(
    left: &[f64],
    right: &[f64],
    left_scale: f64,
    right_scale: f64,
) -> Result<f64, LaError> {
    const REDUCE: f64 = f64::from_bits((1023 - 500) << 52);

    if left.is_empty() {
        return Err(LaError::EmptyVector);
    }
    for (scale, operand) in [
        (left_scale, VectorOperand::Left),
        (right_scale, VectorOperand::Right),
    ] {
        if scale == 0.0 {
            return Err(LaError::ZeroVector { operand });
        }
    }

    let left_scale = DirectionScale::new(left_scale);
    let right_scale = DirectionScale::new(right_scale);
    let mut dot = 0.0_f64;
    let mut dot_error = 0.0_f64;
    let mut exterior_norm = 0.0_f64;
    for (i, (&left_i, &right_i)) in left.iter().zip(right).enumerate() {
        let left_i = left_scale.apply(left_i);
        let right_i = right_scale.apply(right_i);
        // Reduce products before accumulating: even a usize-sized sum and
        // exterior norm retain ample overflow headroom on 64-bit targets.
        let product = (left_i * right_i) * REDUCE;
        let next_dot = dot + product;
        dot_error += two_sum_error(dot, product, next_dot);
        dot = next_dot;
        for (&left_j, &right_j) in left[i + 1..].iter().zip(&right[i + 1..]) {
            let left_j = left_scale.apply(left_j);
            let right_j = right_scale.apply(right_j);
            let minor = difference_of_products(left_i, right_j, left_j, right_i);
            exterior_norm = exterior_norm.hypot(minor * REDUCE);
        }
    }

    // Lagrange's identity: ||left ∧ right|| = ||left|| ||right|| sin(θ).
    // The dot product has the same norm factor times cos(θ). atan2 cancels it
    // without forming the norm product or halving a subnormal angle.
    let angle = exterior_norm.atan2(dot + dot_error);
    if angle.is_finite() {
        Ok(angle)
    } else {
        Err(LaError::non_finite_computation_scalar(
            ArithmeticOperation::VectorAngle,
        ))
    }
}
