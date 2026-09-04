#![forbid(unsafe_code)]

//! Outward-rounded intervals and fixed-size interval determinant signs.

use crate::rounding::{compare_product_with_rounded, two_sum_error};
use crate::{ArithmeticOperation, IntervalBound, IntervalOperand, LaError, Matrix};

/// Largest dimension supported by [`IntervalMatrix::det`] and
/// [`IntervalMatrix::det_sign`].
///
/// A subset-DP determinant needs `2^D` partial intervals. The implementation
/// reserves 128 entries inline, covering the geometry-oriented D ≤ 7 scope
/// without heap allocation.
pub const MAX_INTERVAL_MATRIX_DIM: usize = 7;

/// A closed finite binary64 interval `[lower, upper]`.
///
/// Construction keeps both endpoints finite and ordered. Arithmetic rounds
/// outward, so every successful result contains the exact-real result of the
/// corresponding operation on all represented inputs. Both IEEE-754 signed
/// zeros are accepted and canonicalized to `+0.0`; subnormal bounds are
/// retained.
///
/// This is a deliberately small proof-bearing surface, not a general-purpose
/// interval arithmetic package. Division is intentionally absent.
///
/// # Examples
/// ```
/// use la_stack::{Interval, LaError};
///
/// # fn main() -> Result<(), LaError> {
/// let difference = Interval::try_from_subtraction(1.0, 0.1)?;
/// let square = difference.try_square()?;
/// assert!(difference.lower() < difference.upper());
/// assert!(square.contains((1.0_f64 - 0.1).powi(2)));
/// # Ok(())
/// # }
/// ```
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interval {
    lower: f64,
    upper: f64,
}

/// Sign evidence from an outward-rounded interval determinant.
///
/// `Positive`, `Negative`, and `Zero` are proofs about every determinant
/// represented by the interval matrix. `Inconclusive` means the computed
/// enclosure overlaps zero and must not be interpreted as exact singularity.
#[must_use]
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntervalDeterminantSign {
    /// The determinant interval is strictly greater than zero.
    Positive,
    /// The determinant interval is strictly less than zero.
    Negative,
    /// The determinant interval is exactly `[0, 0]`.
    Zero,
    /// The determinant interval contains zero and at least one nonzero value.
    Inconclusive,
}

/// Fixed-size square matrix of outward-rounded [`Interval`] entries.
///
/// Storage is the inline array `[[Interval; D]; D]`. Determinants use a
/// division-free Leibniz subset DP through D=7, so zero-containing pivot
/// intervals never require a special case and no heap allocation occurs.
///
/// # Examples
/// ```
/// use la_stack::{IntervalDeterminantSign, IntervalMatrix, LaError};
///
/// # fn main() -> Result<(), LaError> {
/// let matrix = IntervalMatrix::<3>::try_from_point_rows([
///     [0.0, 1.0, 0.0],
///     [1.0, 0.0, 0.0],
///     [0.0, 0.0, 1.0],
/// ])?;
/// assert_eq!(matrix.det_sign()?, IntervalDeterminantSign::Negative);
/// # Ok(())
/// # }
/// ```
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IntervalMatrix<const D: usize> {
    rows: [[Interval; D]; D],
}

/// Canonicalize either signed representation of real zero to `+0.0`.
#[inline]
const fn canonical_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

/// Turn a finite rounded sum into the tight adjacent-float enclosure implied by
/// its exact `TwoSum` residual.
#[inline]
const fn rounded_add_bounds(
    left: f64,
    right: f64,
    operation: ArithmeticOperation,
) -> Result<(f64, f64), LaError> {
    let rounded = left + right;
    if !rounded.is_finite() {
        return Err(LaError::interval_range_exhausted(operation));
    }

    let error = two_sum_error(left, right, rounded);
    if !error.is_finite() {
        return Err(LaError::non_finite_computation_scalar(operation));
    }
    let (lower, upper) = if error < 0.0 {
        (rounded.next_down(), rounded)
    } else if error > 0.0 {
        (rounded, rounded.next_up())
    } else {
        (rounded, rounded)
    };
    if !lower.is_finite() || !upper.is_finite() {
        return Err(LaError::interval_range_exhausted(operation));
    }

    Ok((canonical_zero(lower), canonical_zero(upper)))
}

/// Turn a finite rounded product into the tight adjacent-float enclosure of the
/// exact binary64-input product.
#[inline]
const fn rounded_product_bounds(
    left: f64,
    right: f64,
    operation: ArithmeticOperation,
) -> Result<(f64, f64), LaError> {
    if left == 0.0 || right == 0.0 {
        return Ok((0.0, 0.0));
    }

    let rounded = left * right;
    if !rounded.is_finite() {
        return Err(LaError::interval_range_exhausted(operation));
    }

    let relation = compare_product_with_rounded(left, right, rounded);
    let (lower, upper) = if relation < 0 {
        (rounded.next_down(), rounded)
    } else if relation > 0 {
        (rounded, rounded.next_up())
    } else {
        (rounded, rounded)
    };
    if !lower.is_finite() || !upper.is_finite() {
        return Err(LaError::interval_range_exhausted(operation));
    }

    Ok((canonical_zero(lower), canonical_zero(upper)))
}

impl Interval {
    /// Exact real zero.
    pub const ZERO: Self = Self {
        lower: 0.0,
        upper: 0.0,
    };

    /// Exact real one.
    pub const ONE: Self = Self {
        lower: 1.0,
        upper: 1.0,
    };

    /// Construct a closed interval from finite ordered bounds.
    ///
    /// Signed zero endpoints are canonicalized to `+0.0`.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when either endpoint is NaN or infinity.
    /// Returns [`LaError::InvertedInterval`] when `lower > upper`.
    #[inline]
    pub const fn try_new(lower: f64, upper: f64) -> Result<Self, LaError> {
        if !lower.is_finite() {
            return Err(LaError::non_finite_input_interval_bound(
                IntervalBound::Lower,
            ));
        }
        if !upper.is_finite() {
            return Err(LaError::non_finite_input_interval_bound(
                IntervalBound::Upper,
            ));
        }
        if lower > upper {
            return Err(LaError::inverted_interval(lower, upper));
        }
        Ok(Self::new_unchecked(lower, upper))
    }

    /// Construct a point interval from a finite binary64 value.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when `value` is NaN or infinity.
    #[inline]
    pub const fn point(value: f64) -> Result<Self, LaError> {
        match Self::try_new(value, value) {
            Ok(interval) => Ok(interval),
            Err(LaError::NonFinite { .. }) => Err(LaError::non_finite_input_scalar()),
            Err(error) => Err(error),
        }
    }

    /// Enclose the exact-real subtraction of two finite binary64 inputs.
    ///
    /// Unlike subtracting first and then calling [`point`](Self::point), this
    /// method preserves the rounding uncertainty introduced by the subtraction.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] for a non-finite input, preserving whether
    /// it was the left or right operand. Returns
    /// [`LaError::IntervalRangeExhausted`] when the exact difference has no
    /// finite binary64 enclosure.
    #[inline]
    pub const fn try_from_subtraction(left: f64, right: f64) -> Result<Self, LaError> {
        if !left.is_finite() {
            return Err(LaError::non_finite_input_interval_operand(
                IntervalOperand::Left,
            ));
        }
        if !right.is_finite() {
            return Err(LaError::non_finite_input_interval_operand(
                IntervalOperand::Right,
            ));
        }
        match rounded_add_bounds(left, -right, ArithmeticOperation::IntervalSubtraction) {
            Ok((lower, upper)) => Ok(Self::new_unchecked(lower, upper)),
            Err(error) => Err(error),
        }
    }

    /// Return the finite lower bound.
    #[inline]
    #[must_use]
    pub const fn lower(self) -> f64 {
        self.lower
    }

    /// Return the finite upper bound.
    #[inline]
    #[must_use]
    pub const fn upper(self) -> f64 {
        self.upper
    }

    /// Return whether this interval contains the finite `value`.
    #[inline]
    #[must_use]
    pub const fn contains(self, value: f64) -> bool {
        value.is_finite() && self.lower <= value && value <= self.upper
    }

    /// Add two intervals with outward rounding.
    ///
    /// # Errors
    /// Returns [`LaError::IntervalRangeExhausted`] when the exact result range
    /// has no finite binary64 enclosure.
    #[inline]
    pub const fn try_add(&self, other: &Self) -> Result<Self, LaError> {
        self.try_add_for(other, ArithmeticOperation::IntervalAddition)
    }

    /// Multiply two intervals with outward rounding.
    ///
    /// # Errors
    /// Returns [`LaError::IntervalRangeExhausted`] when the exact result range
    /// has no finite binary64 enclosure.
    #[inline]
    pub const fn try_mul(&self, other: &Self) -> Result<Self, LaError> {
        self.try_mul_for(other, ArithmeticOperation::IntervalMultiplication)
    }

    /// Negate an interval exactly by swapping and negating its endpoints.
    #[inline]
    pub const fn negate(&self) -> Self {
        Self::new_unchecked(-self.upper, -self.lower)
    }

    /// Square an interval with outward rounding.
    ///
    /// An interval spanning zero has exact lower bound zero. The upper bound is
    /// the outward-rounded square of the endpoint with greatest magnitude.
    ///
    /// # Errors
    /// Returns [`LaError::IntervalRangeExhausted`] when the exact square range
    /// has no finite binary64 enclosure.
    #[inline]
    pub const fn try_square(&self) -> Result<Self, LaError> {
        let operation = ArithmeticOperation::IntervalSquare;
        let left_square = match rounded_product_bounds(self.lower, self.lower, operation) {
            Ok(bounds) => bounds,
            Err(error) => return Err(error),
        };
        let right_square = match rounded_product_bounds(self.upper, self.upper, operation) {
            Ok(bounds) => bounds,
            Err(error) => return Err(error),
        };
        let lower = if self.lower <= 0.0 && self.upper >= 0.0 {
            0.0
        } else if left_square.0 < right_square.0 {
            left_square.0
        } else {
            right_square.0
        };
        let upper = if left_square.1 > right_square.1 {
            left_square.1
        } else {
            right_square.1
        };
        Ok(Self::new_unchecked(lower, upper))
    }

    /// Construct an interval after its finite ordered-bound invariant is known.
    #[inline]
    const fn new_unchecked(lower: f64, upper: f64) -> Self {
        Self {
            lower: canonical_zero(lower),
            upper: canonical_zero(upper),
        }
    }

    /// Add while attributing range failure to the owning public operation.
    #[inline]
    const fn try_add_for(
        &self,
        other: &Self,
        operation: ArithmeticOperation,
    ) -> Result<Self, LaError> {
        if self.is_zero() {
            return Ok(*other);
        }
        if other.is_zero() {
            return Ok(*self);
        }

        let lower = match rounded_add_bounds(self.lower, other.lower, operation) {
            Ok((lower, _)) => lower,
            Err(error) => return Err(error),
        };
        let upper = match rounded_add_bounds(self.upper, other.upper, operation) {
            Ok((_, upper)) => upper,
            Err(error) => return Err(error),
        };
        Ok(Self::new_unchecked(lower, upper))
    }

    /// Multiply while attributing range failure to the owning public operation.
    #[inline]
    const fn try_mul_for(
        &self,
        other: &Self,
        operation: ArithmeticOperation,
    ) -> Result<Self, LaError> {
        if self.is_zero() || other.is_zero() {
            return Ok(Self::ZERO);
        }
        if self.is_one() {
            return Ok(*other);
        }
        if other.is_one() {
            return Ok(*self);
        }
        if self.is_point() && other.is_point() {
            return match rounded_product_bounds(self.lower, other.lower, operation) {
                Ok((lower, upper)) => Ok(Self::new_unchecked(lower, upper)),
                Err(error) => Err(error),
            };
        }

        self.try_mul_by_sign(other, operation)
    }

    /// Select only the endpoint products that can attain each range extremum.
    #[inline]
    const fn try_mul_by_sign(
        &self,
        other: &Self,
        operation: ArithmeticOperation,
    ) -> Result<Self, LaError> {
        let self_nonnegative = self.lower >= 0.0;
        let self_nonpositive = self.upper <= 0.0;
        let other_nonnegative = other.lower >= 0.0;
        let other_nonpositive = other.upper <= 0.0;

        if self_nonnegative {
            if other_nonnegative {
                return Self::try_product_extrema(
                    (self.lower, other.lower),
                    (self.upper, other.upper),
                    operation,
                );
            }
            if other_nonpositive {
                return Self::try_product_extrema(
                    (self.upper, other.lower),
                    (self.lower, other.upper),
                    operation,
                );
            }
            return Self::try_product_extrema(
                (self.upper, other.lower),
                (self.upper, other.upper),
                operation,
            );
        }
        if self_nonpositive {
            if other_nonnegative {
                return Self::try_product_extrema(
                    (self.lower, other.upper),
                    (self.upper, other.lower),
                    operation,
                );
            }
            if other_nonpositive {
                return Self::try_product_extrema(
                    (self.upper, other.upper),
                    (self.lower, other.lower),
                    operation,
                );
            }
            return Self::try_product_extrema(
                (self.lower, other.upper),
                (self.lower, other.lower),
                operation,
            );
        }
        if other_nonnegative {
            return Self::try_product_extrema(
                (self.lower, other.upper),
                (self.upper, other.upper),
                operation,
            );
        }
        if other_nonpositive {
            return Self::try_product_extrema(
                (self.upper, other.lower),
                (self.lower, other.lower),
                operation,
            );
        }

        let lower_left = match rounded_product_bounds(self.lower, other.upper, operation) {
            Ok(bounds) => bounds,
            Err(error) => return Err(error),
        };
        let lower_right = match rounded_product_bounds(self.upper, other.lower, operation) {
            Ok(bounds) => bounds,
            Err(error) => return Err(error),
        };
        let upper_left = match rounded_product_bounds(self.lower, other.lower, operation) {
            Ok(bounds) => bounds,
            Err(error) => return Err(error),
        };
        let upper_right = match rounded_product_bounds(self.upper, other.upper, operation) {
            Ok(bounds) => bounds,
            Err(error) => return Err(error),
        };
        let lower = if lower_left.0 < lower_right.0 {
            lower_left.0
        } else {
            lower_right.0
        };
        let upper = if upper_left.1 > upper_right.1 {
            upper_left.1
        } else {
            upper_right.1
        };
        Ok(Self::new_unchecked(lower, upper))
    }

    /// Enclose the selected exact lower and upper product extrema.
    #[inline]
    const fn try_product_extrema(
        lower_factors: (f64, f64),
        upper_factors: (f64, f64),
        operation: ArithmeticOperation,
    ) -> Result<Self, LaError> {
        let lower = match rounded_product_bounds(lower_factors.0, lower_factors.1, operation) {
            Ok((lower, _)) => lower,
            Err(error) => return Err(error),
        };
        let upper = match rounded_product_bounds(upper_factors.0, upper_factors.1, operation) {
            Ok((_, upper)) => upper,
            Err(error) => return Err(error),
        };
        Ok(Self::new_unchecked(lower, upper))
    }

    /// Return whether this interval is exactly real zero.
    #[inline]
    const fn is_zero(&self) -> bool {
        self.lower == 0.0 && self.upper == 0.0
    }

    /// Return whether this interval is exactly real one.
    #[inline]
    const fn is_one(&self) -> bool {
        self.lower.to_bits() == 1.0_f64.to_bits() && self.upper.to_bits() == 1.0_f64.to_bits()
    }

    /// Return whether this interval contains one binary64 point.
    #[inline]
    const fn is_point(&self) -> bool {
        self.lower.to_bits() == self.upper.to_bits()
    }
}

impl Default for Interval {
    #[inline]
    fn default() -> Self {
        Self::ZERO
    }
}

impl<const D: usize> IntervalMatrix<D> {
    /// Construct an interval matrix from already-validated interval rows.
    #[inline]
    pub const fn from_rows(rows: [[Interval; D]; D]) -> Self {
        Self { rows }
    }

    /// Lift finite binary64 rows into point intervals.
    ///
    /// This preserves the stored binary64 values exactly; it does not recover
    /// uncertainty from arithmetic performed before this call.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with matrix coordinates for the first NaN
    /// or infinity in row-major order.
    #[inline]
    pub const fn try_from_point_rows(rows: [[f64; D]; D]) -> Result<Self, LaError> {
        let mut intervals = [[Interval::ZERO; D]; D];
        let mut row = 0;
        while row < D {
            let mut column = 0;
            while column < D {
                let value = rows[row][column];
                if !value.is_finite() {
                    return Err(LaError::non_finite_input_matrix(row, column));
                }
                intervals[row][column] = Interval::new_unchecked(value, value);
                column += 1;
            }
            row += 1;
        }
        Ok(Self::from_rows(intervals))
    }

    /// Lift a finite [`Matrix`] into point intervals.
    ///
    /// Earlier rounded expression construction is not enclosed; use interval
    /// operations while constructing derived coefficients when that uncertainty
    /// belongs in the proof.
    #[inline]
    pub const fn from_matrix(matrix: &Matrix<D>) -> Self {
        let matrix_rows = matrix.as_rows();
        let mut intervals = [[Interval::ZERO; D]; D];
        let mut row = 0;
        while row < D {
            let mut column = 0;
            while column < D {
                let value = matrix_rows[row][column];
                intervals[row][column] = Interval::new_unchecked(value, value);
                column += 1;
            }
            row += 1;
        }
        Self::from_rows(intervals)
    }

    /// All-zero interval matrix.
    #[inline]
    pub const fn zero() -> Self {
        Self::from_rows([[Interval::ZERO; D]; D])
    }

    /// Identity interval matrix.
    #[inline]
    pub const fn identity() -> Self {
        let mut matrix = Self::zero();
        let mut index = 0;
        while index < D {
            matrix.rows[index][index] = Interval::ONE;
            index += 1;
        }
        matrix
    }

    /// Borrow the row-major interval storage.
    #[inline]
    pub const fn as_rows(&self) -> &[[Interval; D]; D] {
        &self.rows
    }

    /// Consume this matrix and return its row-major interval storage.
    #[inline]
    pub const fn into_rows(self) -> [[Interval; D]; D] {
        self.rows
    }

    /// Get an interval entry with bounds checking.
    #[inline]
    #[must_use]
    pub const fn get(&self, row: usize, column: usize) -> Option<Interval> {
        if row < D && column < D {
            Some(self.rows[row][column])
        } else {
            None
        }
    }

    /// Get an interval entry while preserving index context on failure.
    ///
    /// # Errors
    /// Returns [`LaError::IndexOutOfBounds`] when either index is not `< D`.
    #[inline]
    pub const fn try_get(&self, row: usize, column: usize) -> Result<Interval, LaError> {
        if row < D && column < D {
            Ok(self.rows[row][column])
        } else {
            Err(LaError::index_out_of_bounds(row, column, D))
        }
    }

    /// Set an interval entry with bounds checking.
    ///
    /// Validation is unnecessary for the value because [`Interval`] already
    /// carries the finite ordered-bound proof.
    ///
    /// # Errors
    /// Returns [`LaError::IndexOutOfBounds`] when either index is not `< D`.
    #[inline]
    pub const fn set(&mut self, row: usize, column: usize, value: Interval) -> Result<(), LaError> {
        if row >= D || column >= D {
            return Err(LaError::index_out_of_bounds(row, column, D));
        }
        self.rows[row][column] = value;
        Ok(())
    }

    /// Enclose the determinant with division-free subset dynamic programming.
    ///
    /// For each column subset, the DP stores the determinant interval of the
    /// leading rows and those columns. This evaluates the Leibniz expansion in
    /// `D × 2^(D-1)` products and additions without choosing or dividing by a
    /// pivot. The returned interval therefore encloses every exact-real
    /// determinant represented by the input intervals, subject only to an
    /// explicit range failure.
    ///
    /// The D=0 determinant follows the empty-product convention and is `[1, 1]`.
    ///
    /// # Errors
    /// Returns [`LaError::UnsupportedDimension`] for D>7. Returns
    /// [`LaError::IntervalRangeExhausted`] with interval-determinant provenance
    /// when an exact intermediate has no finite binary64 enclosure; callers can
    /// then proceed to an exact or higher-range fallback.
    #[inline]
    pub const fn det(&self) -> Result<Interval, LaError> {
        if D > MAX_INTERVAL_MATRIX_DIM {
            return Err(LaError::unsupported_dimension(D, MAX_INTERVAL_MATRIX_DIM));
        }

        let state_count = 1_usize << D;
        let mut partials = [Interval::ZERO; 1 << MAX_INTERVAL_MATRIX_DIM];
        partials[0] = Interval::ONE;
        let operation = ArithmeticOperation::IntervalDeterminant;

        let mut subset = 1;
        while subset < state_count {
            let row = subset.count_ones() as usize - 1;
            let mut sum = Interval::ZERO;
            let mut column = 0;
            while column < D {
                let column_bit = 1_usize << column;
                if subset & column_bit != 0 {
                    let previous = subset ^ column_bit;
                    let mut term =
                        match partials[previous].try_mul_for(&self.rows[row][column], operation) {
                            Ok(term) => term,
                            Err(error) => return Err(error),
                        };
                    let columns_after = (subset >> (column + 1)).count_ones();
                    if !columns_after.is_multiple_of(2) {
                        term = term.negate();
                    }
                    sum = match sum.try_add_for(&term, operation) {
                        Ok(next_sum) => next_sum,
                        Err(error) => return Err(error),
                    };
                }
                column += 1;
            }
            partials[subset] = sum;
            subset += 1;
        }

        Ok(partials[state_count - 1])
    }

    /// Return proof-bearing determinant sign evidence.
    ///
    /// An interval strictly on one side of zero proves that sign. Only the
    /// singleton interval `[0, 0]` proves `Zero`; every other overlap with zero
    /// is [`IntervalDeterminantSign::Inconclusive`].
    ///
    /// # Errors
    /// Propagates the dimension and arithmetic range failures from
    /// [`det`](Self::det).
    #[inline]
    pub const fn det_sign(&self) -> Result<IntervalDeterminantSign, LaError> {
        let determinant = match self.det() {
            Ok(determinant) => determinant,
            Err(error) => return Err(error),
        };
        if determinant.lower > 0.0 {
            Ok(IntervalDeterminantSign::Positive)
        } else if determinant.upper < 0.0 {
            Ok(IntervalDeterminantSign::Negative)
        } else if determinant.lower == 0.0 && determinant.upper == 0.0 {
            Ok(IntervalDeterminantSign::Zero)
        } else {
            Ok(IntervalDeterminantSign::Inconclusive)
        }
    }
}

impl<const D: usize> Default for IntervalMatrix<D> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(test)]
mod tests {
    use core::assert_matches;

    use pastey::paste;

    use super::*;
    use crate::{IntervalBound, IntervalOperand, NonFiniteLocation, NonFiniteOrigin};

    #[test]
    fn point_and_bounds_enforce_interval_invariants() {
        assert_eq!(Interval::point(-0.0).unwrap().lower().to_bits(), 0);
        assert_eq!(Interval::try_new(-0.0, 0.0).unwrap(), Interval::ZERO);
        assert_matches!(
            Interval::point(f64::NAN),
            Err(LaError::NonFinite {
                location: NonFiniteLocation::Scalar,
                origin: NonFiniteOrigin::Input,
                ..
            })
        );
        assert_matches!(
            Interval::try_new(2.0, 1.0),
            Err(LaError::InvertedInterval {
                lower: 2.0,
                upper: 1.0,
                ..
            })
        );
    }

    #[test]
    fn constructors_preserve_non_finite_input_locations() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                Interval::try_new(value, 0.0),
                Err(LaError::NonFinite {
                    location: NonFiniteLocation::IntervalBound {
                        bound: IntervalBound::Lower,
                    },
                    origin: NonFiniteOrigin::Input,
                })
            );
            assert_eq!(
                Interval::try_new(0.0, value),
                Err(LaError::NonFinite {
                    location: NonFiniteLocation::IntervalBound {
                        bound: IntervalBound::Upper,
                    },
                    origin: NonFiniteOrigin::Input,
                })
            );
            assert_eq!(
                Interval::try_from_subtraction(value, 0.0),
                Err(LaError::NonFinite {
                    location: NonFiniteLocation::IntervalOperand {
                        operand: IntervalOperand::Left,
                    },
                    origin: NonFiniteOrigin::Input,
                })
            );
            assert_eq!(
                Interval::try_from_subtraction(0.0, value),
                Err(LaError::NonFinite {
                    location: NonFiniteLocation::IntervalOperand {
                        operand: IntervalOperand::Right,
                    },
                    origin: NonFiniteOrigin::Input,
                })
            );
        }

        assert_eq!(
            Interval::try_new(f64::NAN, f64::INFINITY),
            Err(LaError::NonFinite {
                location: NonFiniteLocation::IntervalBound {
                    bound: IntervalBound::Lower,
                },
                origin: NonFiniteOrigin::Input,
            })
        );
        assert_eq!(
            Interval::try_from_subtraction(f64::NAN, f64::INFINITY),
            Err(LaError::NonFinite {
                location: NonFiniteLocation::IntervalOperand {
                    operand: IntervalOperand::Left,
                },
                origin: NonFiniteOrigin::Input,
            })
        );

        let rows = [[0.0, f64::NAN], [f64::INFINITY, 0.0]];
        assert_eq!(
            IntervalMatrix::<2>::try_from_point_rows(rows),
            Err(LaError::NonFinite {
                location: NonFiniteLocation::MatrixCell { row: 0, col: 1 },
                origin: NonFiniteOrigin::Input,
            })
        );
    }

    #[test]
    fn exact_operations_remain_point_intervals() -> Result<(), LaError> {
        let one = Interval::point(1.0)?;
        let two = Interval::point(2.0)?;
        assert_eq!(one.try_add(&two)?, Interval::point(3.0)?);
        assert_eq!(two.try_mul(&two)?, Interval::point(4.0)?);
        assert_eq!(Interval::try_from_subtraction(3.0, 2.0)?, one);
        assert_eq!(
            Interval::try_new(-2.0, -1.0)?.negate(),
            Interval::try_new(1.0, 2.0)?
        );
        Ok(())
    }

    #[test]
    fn inexact_operations_expand_only_in_the_required_direction() -> Result<(), LaError> {
        let subtraction = Interval::try_from_subtraction(1.0, 0.1)?;
        let rounded_subtraction = 1.0 - 0.1;
        assert!(subtraction.contains(rounded_subtraction));
        assert!(subtraction.lower() < subtraction.upper());

        let product = Interval::point(0.1)?.try_mul(&Interval::point(0.2)?)?;
        assert!(product.contains(0.1 * 0.2));
        assert!(product.lower() < product.upper());

        let below_one = 1.0 - f64::EPSILON;
        let above_one = 1.0 + f64::EPSILON;
        let binade_boundary = Interval::point(below_one)?.try_mul(&Interval::point(above_one)?)?;
        assert_eq!(
            binade_boundary,
            Interval::try_new(1.0_f64.next_down(), 1.0)?
        );
        Ok(())
    }

    #[test]
    fn cancellation_preserves_an_exact_ulp_difference() -> Result<(), LaError> {
        let next = 1.0_f64.next_up();
        let difference = Interval::try_from_subtraction(next, 1.0)?;
        assert_eq!(difference, Interval::point(f64::EPSILON)?);
        Ok(())
    }

    #[test]
    fn underflowed_product_still_encloses_the_positive_exact_result() -> Result<(), LaError> {
        let least_subnormal = f64::from_bits(1);
        let product = Interval::point(least_subnormal)?.try_mul(&Interval::point(0.5)?)?;
        assert_eq!(product, Interval::try_new(0.0, least_subnormal)?);
        Ok(())
    }

    #[test]
    fn range_failure_preserves_interval_operation() -> Result<(), LaError> {
        let error = Interval::point(f64::MAX)?
            .try_mul(&Interval::point(2.0)?)
            .unwrap_err();
        assert_eq!(
            error,
            LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalMultiplication,
            }
        );
        Ok(())
    }

    #[test]
    fn rounded_maximum_detects_exact_sum_beyond_finite_range() -> Result<(), LaError> {
        let maximum = Interval::point(f64::MAX)?;
        let tiny = Interval::point(f64::MIN_POSITIVE)?;
        assert_eq!(
            maximum.try_add(&maximum),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalAddition,
            })
        );
        assert_eq!(
            maximum.try_add(&tiny),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalAddition,
            })
        );
        let nonnegative = Interval::try_new(0.0, f64::MAX)?;
        assert_eq!(
            nonnegative.try_add(&nonnegative),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalAddition,
            })
        );

        let finite_difference = maximum.try_add(&tiny.negate())?;
        assert_eq!(finite_difference.upper().to_bits(), f64::MAX.to_bits());
        assert!(finite_difference.lower() < finite_difference.upper());
        Ok(())
    }

    #[test]
    fn subtraction_and_square_preserve_distinct_range_operations() -> Result<(), LaError> {
        assert_eq!(
            Interval::try_from_subtraction(f64::MAX, -f64::MIN_POSITIVE),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalSubtraction,
            })
        );
        assert_eq!(
            Interval::point(f64::MAX)?.try_square(),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalSquare,
            })
        );
        assert_eq!(
            Interval::try_new(-1.0, f64::MAX)?.try_square(),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalSquare,
            })
        );
        Ok(())
    }

    #[test]
    fn determinant_overflow_reports_interval_determinant_range_failure() -> Result<(), LaError> {
        let matrix = IntervalMatrix::<2>::try_from_point_rows([[f64::MAX, 0.0], [0.0, 2.0]])?;
        assert_eq!(
            matrix.det(),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalDeterminant,
            })
        );

        let accumulating =
            IntervalMatrix::<2>::try_from_point_rows([[f64::MAX, f64::MAX], [-1.0, 1.0]])?;
        assert_eq!(
            accumulating.det(),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalDeterminant,
            })
        );
        assert_eq!(
            accumulating.det_sign(),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalDeterminant,
            })
        );
        Ok(())
    }

    #[test]
    fn determinant_reports_intermediate_exhaustion_before_exact_cancellation() -> Result<(), LaError>
    {
        let matrix = IntervalMatrix::<2>::try_from_point_rows([[f64::MAX, f64::MAX], [2.0, 2.0]])?;
        assert_eq!(
            matrix.det(),
            Err(LaError::IntervalRangeExhausted {
                operation: ArithmeticOperation::IntervalDeterminant,
            })
        );
        Ok(())
    }

    #[test]
    fn square_spanning_zero_has_exact_zero_lower_bound() -> Result<(), LaError> {
        let square = Interval::try_new(-2.0, 3.0)?.try_square()?;
        assert_eq!(square, Interval::try_new(0.0, 9.0)?);
        Ok(())
    }

    #[test]
    fn multiplication_selects_correct_extrema_in_every_sign_quadrant() -> Result<(), LaError> {
        for (left, right, expected) in [
            ((2.0, 3.0), (4.0, 5.0), (8.0, 15.0)),
            ((2.0, 3.0), (-5.0, -4.0), (-15.0, -8.0)),
            ((2.0, 3.0), (-5.0, 4.0), (-15.0, 12.0)),
            ((-3.0, -2.0), (4.0, 5.0), (-15.0, -8.0)),
            ((-3.0, -2.0), (-5.0, -4.0), (8.0, 15.0)),
            ((-3.0, -2.0), (-5.0, 4.0), (-12.0, 15.0)),
            ((-3.0, 2.0), (4.0, 5.0), (-15.0, 10.0)),
            ((-3.0, 2.0), (-5.0, -4.0), (-10.0, 15.0)),
            ((-3.0, 2.0), (-5.0, 4.0), (-12.0, 15.0)),
        ] {
            let product = Interval::try_new(left.0, left.1)?
                .try_mul(&Interval::try_new(right.0, right.1)?)?;
            assert_eq!(product, Interval::try_new(expected.0, expected.1)?);
        }

        assert_eq!(
            Interval::ZERO.try_mul(&Interval::try_new(-f64::MAX, f64::MAX)?)?,
            Interval::ZERO
        );
        Ok(())
    }

    #[test]
    fn multiplication_rejects_unrepresentable_selected_extrema() -> Result<(), LaError> {
        let half_maximum = f64::MAX / 2.0;
        for (left, right) in [
            ((half_maximum, f64::MAX), (-2.0, -1.0)),
            ((half_maximum, f64::MAX), (1.0, 2.0)),
            ((-f64::MAX, 1.0), (-1.0, 2.0)),
            ((-1.0, f64::MAX), (-2.0, 1.0)),
            ((-f64::MAX, 1.0), (-2.0, 1.0)),
            ((-1.0, f64::MAX), (-1.0, 2.0)),
        ] {
            let result =
                Interval::try_new(left.0, left.1)?.try_mul(&Interval::try_new(right.0, right.1)?);
            assert_eq!(
                result,
                Err(LaError::IntervalRangeExhausted {
                    operation: ArithmeticOperation::IntervalMultiplication,
                }),
                "left={left:?}, right={right:?}"
            );
        }
        Ok(())
    }

    macro_rules! gen_interval_identity_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<interval_identity_sign_is_positive_ $d d>]() {
                    let matrix = IntervalMatrix::<$d>::identity();
                    assert_eq!(matrix.det(), Ok(Interval::ONE));
                    assert_eq!(
                        matrix.det_sign(),
                        Ok(IntervalDeterminantSign::Positive)
                    );
                }
            }
        };
    }

    gen_interval_identity_tests!(2);
    gen_interval_identity_tests!(3);
    gen_interval_identity_tests!(4);
    gen_interval_identity_tests!(5);
    gen_interval_identity_tests!(6);
    gen_interval_identity_tests!(7);

    #[test]
    fn determinant_sign_handles_row_swap_and_exact_singularity() -> Result<(), LaError> {
        let swapped = IntervalMatrix::<3>::try_from_point_rows([
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ])?;
        assert_eq!(swapped.det_sign()?, IntervalDeterminantSign::Negative);

        let singular = IntervalMatrix::<3>::try_from_point_rows([
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [0.0, 0.0, 1.0],
        ])?;
        assert_eq!(singular.det_sign()?, IntervalDeterminantSign::Zero);
        Ok(())
    }

    #[test]
    fn wide_determinant_interval_is_inconclusive() -> Result<(), LaError> {
        let matrix = IntervalMatrix::<2>::from_rows([
            [Interval::ONE, Interval::ZERO],
            [Interval::ZERO, Interval::try_new(-1.0, 1.0)?],
        ]);
        assert_eq!(matrix.det_sign()?, IntervalDeterminantSign::Inconclusive);
        Ok(())
    }

    #[test]
    fn determinant_rejects_dimensions_above_supported_stack_dp() {
        assert_matches!(
            IntervalMatrix::<8>::identity().det(),
            Err(LaError::UnsupportedDimension {
                requested: 8,
                max: MAX_INTERVAL_MATRIX_DIM,
                ..
            })
        );
    }

    #[test]
    fn matrix_accessors_preserve_validated_storage() -> Result<(), LaError> {
        let source = Matrix::<2>::identity();
        let mut intervals = IntervalMatrix::from_matrix(&source);
        let value = Interval::try_new(2.0, 3.0)?;
        intervals.set(0, 1, value)?;
        assert_eq!(intervals.get(0, 1), Some(value));
        assert_eq!(intervals.get(2, 0), None);
        assert_eq!(intervals.try_get(0, 1)?, value);
        assert_matches!(
            intervals.try_get(2, 0),
            Err(LaError::IndexOutOfBounds {
                row: 2,
                col: 0,
                dim: 2,
                ..
            })
        );
        assert_eq!(intervals.as_rows()[0][1], value);
        assert_eq!(intervals.into_rows()[0][1], value);
        Ok(())
    }

    #[test]
    fn rejected_matrix_set_is_failure_atomic() -> Result<(), LaError> {
        let mut matrix = IntervalMatrix::<2>::identity();
        let before = matrix;
        let value = Interval::try_new(2.0, 3.0)?;

        assert_eq!(
            matrix.set(2, 0, value),
            Err(LaError::IndexOutOfBounds {
                row: 2,
                col: 0,
                dim: 2,
            })
        );
        assert_eq!(matrix, before);
        assert_eq!(
            matrix.set(0, 2, value),
            Err(LaError::IndexOutOfBounds {
                row: 0,
                col: 2,
                dim: 2,
            })
        );
        assert_eq!(matrix, before);
        Ok(())
    }
}
