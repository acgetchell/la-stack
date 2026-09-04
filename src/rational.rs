#![forbid(unsafe_code)]

//! Exact-input fixed-size matrices and vectors.
//!
//! [`RationalMatrix`] and [`RationalVector`] preserve caller-supplied
//! [`BigRational`] coefficients without a binary64 round trip. Determinants and
//! solves clear denominators with a positive scale per row, then reuse the
//! crate's fraction-free [`BigInt`] Bareiss backend. The positive row scales
//! preserve determinant sign; determinant values divide by their product; and
//! solves apply the same row scale to the matrix and right-hand side.

use std::array::from_fn;

use num_bigint::{BigInt, Sign};
use num_rational::BigRational;

use crate::exact::{det_big_int, solve_big_int};
use crate::{DeterminantSign, ExactF64Conversion, LaError, Vector};

/// Exact rational square matrix with compile-time dimension `D`.
///
/// Construction validates that every denominator is non-zero and canonicalizes
/// every entry to lowest terms with a positive denominator. The private storage
/// then carries those invariants, so determinant and solve methods do not repeat
/// input validation. Unlike [`crate::Matrix`], entries are already exact rational
/// values rather than finite binary64 values interpreted exactly.
///
/// Direct field construction is intentionally unavailable:
///
/// ```compile_fail
/// use la_stack::{BigRational, RationalMatrix};
///
/// let _ = RationalMatrix::<1> {
///     rows: [[BigRational::from_integer(1.into())]],
/// };
/// ```
#[derive(Clone, Debug, Eq, PartialEq)]
#[must_use]
pub struct RationalMatrix<const D: usize> {
    rows: [[BigRational; D]; D],
}

/// Exact rational vector with compile-time dimension `D`.
///
/// Construction validates that every denominator is non-zero and canonicalizes
/// every entry to lowest terms with a positive denominator. Solutions returned
/// by [`RationalMatrix::solve`] and [`crate::Matrix::solve_exact`] also use this
/// type, making any later conversion to [`Vector`] explicit through
/// [`ExactF64Conversion`].
#[derive(Clone, Debug, Eq, PartialEq)]
#[must_use]
pub struct RationalVector<const D: usize> {
    data: [BigRational; D],
}

impl<const D: usize> RationalMatrix<D> {
    /// Try to create an exact matrix from row-major rational storage.
    ///
    /// Raw, non-reduced [`BigRational::new_raw`] values and negative
    /// denominators are accepted, interpreted as their mathematical quotient,
    /// and stored canonically. A raw zero denominator is not a rational value
    /// and is rejected at this construction boundary.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let matrix = RationalMatrix::<2>::try_from_rows([
    ///     [
    ///         BigRational::new(1.into(), 3.into()),
    ///         BigRational::from_integer(2.into()),
    ///     ],
    ///     [
    ///         BigRational::from_integer(1.into()),
    ///         BigRational::new(5.into(), 2.into()),
    ///     ],
    /// ])?;
    /// assert_eq!(
    ///     matrix.det(),
    ///     BigRational::new((-7).into(), 6.into())
    /// );
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] at the first matrix cell whose raw
    /// rational denominator is zero.
    pub fn try_from_rows(rows: [[BigRational; D]; D]) -> Result<Self, LaError> {
        for (row_index, row) in rows.iter().enumerate() {
            for (col_index, value) in row.iter().enumerate() {
                if value.denom().sign() == Sign::NoSign {
                    return Err(LaError::non_finite_input_matrix(row_index, col_index));
                }
            }
        }
        Ok(Self {
            rows: rows.map(|row| row.map(canonicalize_rational)),
        })
    }

    /// Try to create an exact matrix by evaluating a function at every cell.
    ///
    /// The function is evaluated once per cell in row-major order.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let diagonal = RationalMatrix::<3>::try_from_fn(|row, col| {
    ///     BigRational::from_integer(u8::from(row == col).into())
    /// })?;
    /// assert_eq!(diagonal.det_sign(), DeterminantSign::Positive);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] at the first generated cell whose raw
    /// rational denominator is zero.
    pub fn try_from_fn(
        mut make_entry: impl FnMut(usize, usize) -> BigRational,
    ) -> Result<Self, LaError> {
        let rows = from_fn(|row| from_fn(|col| make_entry(row, col)));
        Self::try_from_rows(rows)
    }

    /// Return the all-zero exact matrix.
    pub fn zero() -> Self {
        Self {
            rows: from_fn(|_| from_fn(|_| BigRational::from_integer(BigInt::from(0)))),
        }
    }

    /// Borrow the row-major exact storage.
    #[must_use]
    pub const fn as_rows(&self) -> &[[BigRational; D]; D] {
        &self.rows
    }

    /// Consume the matrix and return its row-major exact storage.
    #[must_use]
    pub fn into_rows(self) -> [[BigRational; D]; D] {
        self.rows
    }

    /// Borrow one entry, returning `None` for an out-of-bounds index.
    #[must_use]
    pub fn get(&self, row: usize, col: usize) -> Option<&BigRational> {
        self.rows.get(row)?.get(col)
    }

    /// Replace one exact entry while preserving the canonical non-zero-
    /// denominator invariant.
    ///
    /// # Errors
    /// Returns [`LaError::IndexOutOfBounds`] when `(row, col)` lies outside the
    /// matrix, or [`LaError::NonFinite`] when `value` has a raw zero
    /// denominator.
    pub fn set(&mut self, row: usize, col: usize, value: BigRational) -> Result<(), LaError> {
        if row >= D || col >= D {
            return Err(LaError::index_out_of_bounds(row, col, D));
        }
        if value.denom().sign() == Sign::NoSign {
            return Err(LaError::non_finite_input_matrix(row, col));
        }
        self.rows[row][col] = canonicalize_rational(value);
        Ok(())
    }

    /// Return the provably exact determinant sign.
    ///
    /// This path clears denominators and reads the sign of the resulting
    /// integer determinant. It does not construct a rational determinant.
    /// For D=0, the empty-product determinant has positive sign.
    pub fn det_sign(&self) -> DeterminantSign {
        let (integer_rows, _) = self.integer_rows();
        match det_big_int(integer_rows).sign() {
            Sign::Minus => DeterminantSign::Negative,
            Sign::NoSign => DeterminantSign::Zero,
            Sign::Plus => DeterminantSign::Positive,
        }
    }

    /// Return the exact determinant.
    ///
    /// Denominators are cleared independently per row. If row `i` uses
    /// positive scale `sᵢ`, the integer determinant is divided by `∏ᵢ sᵢ`.
    /// For D=0, this returns the empty-product determinant `1`.
    #[must_use]
    pub fn det(&self) -> BigRational {
        let (integer_rows, row_scales) = self.integer_rows();
        let determinant_denominator = row_scales.iter().product();
        BigRational::new(det_big_int(integer_rows), determinant_denominator)
    }

    /// Solve `A x = b` exactly.
    ///
    /// Each augmented row is multiplied by one positive common denominator,
    /// then fraction-free Bareiss forward elimination runs in [`BigInt`]. Only
    /// the `O(D²)` back-substitution phase constructs [`BigRational`] values.
    /// For D=0, the empty matrix and vector have the unique empty solution.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let zero = BigRational::from_integer(0.into());
    /// let one = BigRational::from_integer(1.into());
    /// let matrix = RationalMatrix::<2>::try_from_rows([
    ///     [BigRational::new(1.into(), 2.into()), zero.clone()],
    ///     [zero, BigRational::new(1.into(), 3.into())],
    /// ])?;
    /// let rhs = RationalVector::try_new([one.clone(), one])?;
    ///
    /// let solution = matrix.solve(&rhs)?.try_to_f64()?.into_array();
    /// assert_eq!(solution, [2.0, 3.0]);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] with exact-singularity metadata when a
    /// pivot column contains no non-zero entry.
    pub fn solve(&self, rhs: &RationalVector<D>) -> Result<RationalVector<D>, LaError> {
        let row_scales: [BigInt; D] = from_fn(|row| {
            common_denominator(
                self.rows[row]
                    .iter()
                    .chain(core::iter::once(&rhs.data[row])),
            )
        });
        let integer_rows =
            from_fn(|row| from_fn(|col| integer_at_scale(&self.rows[row][col], &row_scales[row])));
        let integer_rhs = from_fn(|row| integer_at_scale(&rhs.data[row], &row_scales[row]));
        solve_big_int(integer_rows, integer_rhs).map(RationalVector::from_canonical_array)
    }

    /// Clear matrix denominators with one positive common denominator per row.
    fn integer_rows(&self) -> ([[BigInt; D]; D], [BigInt; D]) {
        let row_scales: [BigInt; D] = from_fn(|row| common_denominator(self.rows[row].iter()));
        let integer_rows =
            from_fn(|row| from_fn(|col| integer_at_scale(&self.rows[row][col], &row_scales[row])));
        (integer_rows, row_scales)
    }
}

impl<const D: usize> RationalVector<D> {
    /// Wrap values produced by `BigRational` arithmetic, which preserves the
    /// canonical representation established at public input boundaries.
    pub(crate) const fn from_canonical_array(data: [BigRational; D]) -> Self {
        Self { data }
    }

    /// Try to create an exact vector from rational storage.
    ///
    /// Raw, non-reduced values and negative denominators are accepted and
    /// stored canonically. A raw zero denominator is rejected.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let rhs = RationalVector::<2>::try_new([
    ///     BigRational::new(1.into(), 2.into()),
    ///     BigRational::from_integer(3.into()),
    /// ])?;
    /// assert_eq!(rhs.try_to_f64()?.into_array(), [0.5, 3.0]);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] at the first vector entry whose raw
    /// rational denominator is zero.
    pub fn try_new(data: [BigRational; D]) -> Result<Self, LaError> {
        for (index, value) in data.iter().enumerate() {
            if value.denom().sign() == Sign::NoSign {
                return Err(LaError::non_finite_input_vector(index));
            }
        }
        Ok(Self {
            data: data.map(canonicalize_rational),
        })
    }

    /// Try to create an exact vector by evaluating a function at every index.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] at the first generated entry whose raw
    /// rational denominator is zero.
    pub fn try_from_fn(make_entry: impl FnMut(usize) -> BigRational) -> Result<Self, LaError> {
        Self::try_new(from_fn(make_entry))
    }

    /// Return the all-zero exact vector.
    pub fn zero() -> Self {
        Self {
            data: from_fn(|_| BigRational::from_integer(BigInt::from(0))),
        }
    }

    /// Borrow the exact backing array.
    #[must_use]
    pub const fn as_array(&self) -> &[BigRational; D] {
        &self.data
    }

    /// Consume the vector and return its exact backing array.
    #[must_use]
    pub fn into_array(self) -> [BigRational; D] {
        self.data
    }

    /// Borrow one entry, returning `None` for an out-of-bounds index.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<&BigRational> {
        self.data.get(index)
    }
}

impl<const D: usize> ExactF64Conversion for RationalVector<D> {
    type Output = Vector<D>;

    fn try_to_f64(&self) -> Result<Self::Output, LaError> {
        self.data.try_to_f64()
    }

    fn to_rounded_f64(&self) -> Result<Self::Output, LaError> {
        self.data.to_rounded_f64()
    }
}

/// Reduce one validated rational and make its denominator positive before
/// publishing it through the exact-input storage types.
fn canonicalize_rational(value: BigRational) -> BigRational {
    let (numerator, denominator) = value.into_raw();
    BigRational::new(numerator, denominator)
}

/// Return a positive least common multiple of all raw denominator magnitudes.
fn common_denominator<'a>(values: impl Iterator<Item = &'a BigRational>) -> BigInt {
    values.fold(BigInt::from(1), |scale, value| {
        least_common_multiple(scale, denominator_magnitude(value))
    })
}

/// Return the positive magnitude of a validated non-zero denominator.
fn denominator_magnitude(value: &BigRational) -> BigInt {
    match value.denom().sign() {
        Sign::Minus => -value.denom(),
        Sign::Plus => value.denom().clone(),
        Sign::NoSign => unreachable!("RationalMatrix and RationalVector validate denominators"),
    }
}

/// Convert a rational to an integer using a positive divisible scale.
fn integer_at_scale(value: &BigRational, scale: &BigInt) -> BigInt {
    let denominator = denominator_magnitude(value);
    let multiplier = scale / denominator;
    let numerator = match value.denom().sign() {
        Sign::Minus => -value.numer(),
        Sign::Plus => value.numer().clone(),
        Sign::NoSign => unreachable!("RationalMatrix and RationalVector validate denominators"),
    };
    numerator * multiplier
}

/// Return the positive least common multiple of two positive integers.
fn least_common_multiple(lhs: BigInt, rhs: BigInt) -> BigInt {
    let gcd = greatest_common_divisor(lhs.clone(), rhs.clone());
    (lhs / gcd) * rhs
}

/// Euclidean greatest common divisor for positive integers.
fn greatest_common_divisor(mut lhs: BigInt, mut rhs: BigInt) -> BigInt {
    let zero = BigInt::from(0);
    while rhs != zero {
        let remainder = &lhs % &rhs;
        lhs = rhs;
        rhs = remainder;
    }
    lhs
}

#[cfg(test)]
mod tests {
    use pastey::paste;

    use super::*;
    use crate::{NonFiniteLocation, NonFiniteOrigin, SingularityReason, UnrepresentableReason};

    fn ratio(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn exact_input_preserves_non_binary64_coefficients() {
        let matrix = RationalMatrix::<2>::try_from_rows([
            [ratio(1, 3), ratio(1, 10)],
            [ratio(2, 7), ratio(3, 11)],
        ])
        .unwrap();

        assert_eq!(matrix.det(), ratio(24, 385));
        assert_eq!(matrix.det_sign(), DeterminantSign::Positive);
    }

    #[test]
    fn determinant_handles_pivoting_and_singularity() {
        let pivoting = RationalMatrix::<3>::try_from_rows([
            [ratio(0, 1), ratio(1, 2), ratio(0, 1)],
            [ratio(2, 3), ratio(0, 1), ratio(0, 1)],
            [ratio(0, 1), ratio(0, 1), ratio(3, 5)],
        ])
        .unwrap();
        assert_eq!(pivoting.det(), ratio(-1, 5));
        assert_eq!(pivoting.det_sign(), DeterminantSign::Negative);

        let singular = RationalMatrix::<3>::try_from_rows([
            [ratio(1, 2), ratio(1, 3), ratio(1, 5)],
            [ratio(1, 1), ratio(2, 3), ratio(2, 5)],
            [ratio(0, 1), ratio(1, 1), ratio(1, 1)],
        ])
        .unwrap();
        assert_eq!(singular.det(), ratio(0, 1));
        assert_eq!(singular.det_sign(), DeterminantSign::Zero);
    }

    #[test]
    fn row_clearing_accepts_raw_signs_factors_and_dyadic_exponents() {
        let matrix = RationalMatrix::<2>::try_from_rows([
            [
                BigRational::new_raw(BigInt::from(6), BigInt::from(-8)),
                BigRational::new_raw(BigInt::from(3), BigInt::from(1_u8) << 80_u32),
            ],
            [
                BigRational::new_raw(BigInt::from(-10), BigInt::from(-20)),
                BigRational::new_raw(BigInt::from(14), BigInt::from(21)),
            ],
        ])
        .unwrap();

        let expected =
            ratio(-1, 2) - BigRational::new(BigInt::from(3), BigInt::from(1_u8) << 81_u32);
        assert_eq!(matrix.det(), expected);
        assert_eq!(matrix.det_sign(), DeterminantSign::Negative);
        assert_eq!(matrix.as_rows()[0][0], ratio(-3, 4));
        assert_eq!(matrix.as_rows()[1][0], ratio(1, 2));
        assert_eq!(matrix.as_rows()[1][1], ratio(2, 3));
    }

    #[test]
    fn construction_and_set_store_equal_quotients_identically() {
        let huge_factor = BigInt::from(1_u8) << 256_u32;
        let raw = RationalMatrix::<2>::try_from_rows([
            [
                BigRational::new_raw(huge_factor.clone(), &huge_factor * 2_u8),
                BigRational::new_raw(BigInt::from(0), -&huge_factor),
            ],
            [ratio(0, 1), ratio(1, 1)],
        ])
        .unwrap();
        let canonical = RationalMatrix::<2>::try_from_rows([
            [ratio(1, 2), ratio(0, 1)],
            [ratio(0, 1), ratio(1, 1)],
        ])
        .unwrap();

        assert_eq!(raw, canonical);
        assert_eq!(raw.as_rows()[0][1].denom(), &BigInt::from(1));

        let mut updated = RationalMatrix::<2>::zero();
        updated
            .set(
                0,
                0,
                BigRational::new_raw(&huge_factor * 3_u8, -&huge_factor * 6_u8),
            )
            .unwrap();
        assert_eq!(updated.get(0, 0), Some(&ratio(-1, 2)));
    }

    #[test]
    fn exact_solve_scales_matrix_and_rhs_together() {
        let matrix = RationalMatrix::<3>::try_from_rows([
            [ratio(0, 1), ratio(1, 3), ratio(1, 5)],
            [ratio(2, 7), ratio(1, 11), ratio(0, 1)],
            [ratio(1, 13), ratio(0, 1), ratio(3, 2)],
        ])
        .unwrap();
        let expected = [ratio(2, 3), ratio(-4, 5), ratio(7, 9)];
        let rhs = RationalVector::try_from_fn(|row| {
            matrix.as_rows()[row]
                .iter()
                .zip(expected.iter())
                .map(|(coefficient, solution)| coefficient * solution)
                .sum()
        })
        .unwrap();

        assert_eq!(matrix.solve(&rhs).unwrap().into_array(), expected);
    }

    #[test]
    fn exact_solve_accepts_signed_unreduced_matrix_and_rhs_values() {
        let matrix = RationalMatrix::<2>::try_from_rows([
            [
                BigRational::new_raw(BigInt::from(10), BigInt::from(20)),
                BigRational::new_raw(BigInt::from(-5), BigInt::from(-15)),
            ],
            [
                BigRational::new_raw(BigInt::from(-6), BigInt::from(-15)),
                BigRational::new_raw(BigInt::from(9), BigInt::from(21)),
            ],
        ])
        .unwrap();
        let rhs = RationalVector::try_new([
            BigRational::new_raw(BigInt::from(0), BigInt::from(-99)),
            BigRational::new_raw(BigInt::from(34), BigInt::from(-70)),
        ])
        .unwrap();

        let solution = matrix.solve(&rhs).unwrap();
        assert_eq!(solution.as_array(), &[ratio(2, 1), ratio(-3, 1)]);
        assert_eq!(
            &matrix.as_rows()[0][0] * &solution.as_array()[0]
                + &matrix.as_rows()[0][1] * &solution.as_array()[1],
            rhs.as_array()[0]
        );
        assert_eq!(
            &matrix.as_rows()[1][0] * &solution.as_array()[0]
                + &matrix.as_rows()[1][1] * &solution.as_array()[1],
            rhs.as_array()[1]
        );
    }

    #[test]
    fn singular_solve_preserves_exact_pivot_metadata() {
        let matrix = RationalMatrix::<2>::try_from_rows([
            [ratio(1, 2), ratio(1, 3)],
            [ratio(1, 1), ratio(2, 3)],
        ])
        .unwrap();
        let rhs = RationalVector::<2>::zero();

        assert!(matches!(
            matrix.solve(&rhs),
            Err(LaError::Singular {
                pivot_col: 1,
                reason: SingularityReason::Exact,
                ..
            })
        ));
    }

    #[test]
    fn constructors_reject_raw_zero_denominators_with_locations() {
        let matrix_error = RationalMatrix::<1>::try_from_rows([[BigRational::new_raw(
            BigInt::from(1),
            BigInt::from(0),
        )]])
        .unwrap_err();
        assert!(matches!(
            matrix_error,
            LaError::NonFinite {
                location: NonFiniteLocation::MatrixCell { row: 0, col: 0, .. },
                origin: NonFiniteOrigin::Input,
                ..
            }
        ));

        let vector_error =
            RationalVector::<1>::try_new([BigRational::new_raw(BigInt::from(1), BigInt::from(0))])
                .unwrap_err();
        assert!(matches!(
            vector_error,
            LaError::NonFinite {
                location: NonFiniteLocation::VectorEntry { index: 0, .. },
                origin: NonFiniteOrigin::Input,
                ..
            }
        ));
    }

    macro_rules! gen_rejected_set_is_failure_atomic_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<rejected_set_is_failure_atomic_ $d d>]() {
                    let mut matrix = RationalMatrix::<$d>::try_from_fn(|row, col| {
                        BigRational::from_integer(BigInt::from(row * $d + col + 1))
                    })
                    .unwrap();
                    let original = matrix.clone();

                    assert!(matches!(
                        matrix.set($d, 0, ratio(1, 2)),
                        Err(LaError::IndexOutOfBounds {
                            row: $d,
                            col: 0,
                            dim: $d,
                            ..
                        })
                    ));
                    assert_eq!(matrix, original);

                    let zero_denominator =
                        BigRational::new_raw(BigInt::from(1), BigInt::from(0));
                    assert!(matches!(
                        matrix.set(0, 1, zero_denominator),
                        Err(LaError::NonFinite {
                            location: NonFiniteLocation::MatrixCell { row: 0, col: 1, .. },
                            origin: NonFiniteOrigin::Input,
                            ..
                        })
                    ));
                    assert_eq!(matrix, original);
                }
            }
        };
    }

    gen_rejected_set_is_failure_atomic_tests!(2);
    gen_rejected_set_is_failure_atomic_tests!(3);
    gen_rejected_set_is_failure_atomic_tests!(4);
    gen_rejected_set_is_failure_atomic_tests!(5);

    #[test]
    fn zero_dimension_uses_empty_determinant_and_unique_solve() {
        let matrix = RationalMatrix::<0>::zero();
        assert_eq!(matrix.det(), ratio(1, 1));
        assert_eq!(matrix.det_sign(), DeterminantSign::Positive);
        assert_eq!(
            matrix.solve(&RationalVector::zero()).unwrap().into_array(),
            []
        );
    }

    #[test]
    fn rational_vector_conversion_is_explicit() {
        let vector = RationalVector::<2>::try_new([ratio(1, 2), ratio(1, 3)]).unwrap();
        assert!(matches!(
            vector.try_to_f64(),
            Err(LaError::Unrepresentable {
                index: Some(1),
                reason: UnrepresentableReason::RequiresRounding,
                ..
            })
        ));
        let rounded = vector.to_rounded_f64().unwrap().into_array();
        assert_eq!(rounded[0].to_bits(), 0.5_f64.to_bits());
        assert_eq!(rounded[1].to_bits(), (1.0_f64 / 3.0).to_bits());
    }

    macro_rules! gen_pivoting_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<pivoting_determinant_and_solve_ $d d>]() {
                    let matrix = RationalMatrix::<$d>::try_from_fn(|row, col| {
                        let is_permutation_entry = (row == 0 && col == 1)
                            || (row == 1 && col == 0)
                            || (row >= 2 && row == col);
                        BigRational::from_integer(BigInt::from(u8::from(is_permutation_entry)))
                    })
                    .unwrap();
                    let expected = from_fn(|index| {
                        BigRational::new(BigInt::from(1), BigInt::from(index + 2))
                    });
                    let rhs = RationalVector::try_from_fn(|row| {
                        matrix.as_rows()[row]
                            .iter()
                            .zip(expected.iter())
                            .map(|(coefficient, component)| coefficient * component)
                            .sum()
                    })
                    .unwrap();

                    assert_eq!(matrix.det_sign(), DeterminantSign::Negative);
                    assert_eq!(matrix.det(), ratio(-1, 1));
                    assert_eq!(matrix.solve(&rhs).unwrap().into_array(), expected);
                }
            }
        };
    }

    gen_pivoting_tests!(2);
    gen_pivoting_tests!(3);
    gen_pivoting_tests!(4);
    gen_pivoting_tests!(5);
    gen_pivoting_tests!(6);
    gen_pivoting_tests!(7);
    gen_pivoting_tests!(8);
}
