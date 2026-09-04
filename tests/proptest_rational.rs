#![forbid(unsafe_code)]

//! Independent property checks for exact rational-input matrices through D=8.

#![cfg(feature = "exact")]

use std::array::from_fn;

use pastey::paste;
use proptest::prelude::*;

use la_stack::prelude::*;

#[path = "common/proptest_config.rs"]
mod proptest_config;
use proptest_config::with_default_cases;

fn zero() -> BigRational {
    BigRational::from_integer(BigInt::from(0))
}

fn rational(numerator: i16, denominator: u8) -> BigRational {
    BigRational::new(
        BigInt::from(numerator),
        BigInt::from(u32::from(denominator)),
    )
}

fn determinant_sign(value: &BigRational) -> DeterminantSign {
    if value.is_positive() {
        DeterminantSign::Positive
    } else if value.is_negative() {
        DeterminantSign::Negative
    } else {
        DeterminantSign::Zero
    }
}

/// Straightforward rational Gaussian elimination, independent of the
/// production denominator-clearing and integer Bareiss backend.
fn rational_determinant_gaussian<const D: usize>(mut rows: [[BigRational; D]; D]) -> BigRational {
    let zero = zero();
    let mut determinant = BigRational::from_integer(BigInt::from(1));
    let mut odd_swaps = false;

    for pivot_col in 0..D {
        let Some(pivot_row) = (pivot_col..D).find(|&row| rows[row][pivot_col] != zero) else {
            return zero;
        };
        if pivot_row != pivot_col {
            rows.swap(pivot_col, pivot_row);
            odd_swaps = !odd_swaps;
        }

        let pivot = rows[pivot_col][pivot_col].clone();
        let pivot_entries = rows[pivot_col].clone();
        determinant *= &pivot;
        for row_entries in rows.iter_mut().skip(pivot_col + 1) {
            let factor = &row_entries[pivot_col] / &pivot;
            for (entry, pivot_entry) in row_entries
                .iter_mut()
                .zip(pivot_entries.iter())
                .skip(pivot_col + 1)
            {
                *entry -= &factor * pivot_entry;
            }
            row_entries[pivot_col] = zero.clone();
        }
    }

    if odd_swaps { -determinant } else { determinant }
}

fn rational_matvec<const D: usize>(
    rows: &[[BigRational; D]; D],
    vector: &[BigRational; D],
) -> [BigRational; D] {
    from_fn(|row| {
        rows[row]
            .iter()
            .zip(vector.iter())
            .map(|(coefficient, component)| coefficient * component)
            .sum()
    })
}

macro_rules! gen_rational_properties {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(with_default_cases(12))]

                #[test]
                fn [<rational_det_sign_value_and_solve_agree_with_reference_ $d d>](
                    entries in proptest::collection::vec((-5_i16..=5_i16, 1_u8..=9_u8), $d * $d),
                    solution_entries in proptest::collection::vec((-5_i16..=5_i16, 1_u8..=9_u8), $d),
                ) {
                    let rows: [[BigRational; $d]; $d] = from_fn(|row| {
                        from_fn(|col| {
                            let (numerator, denominator) = entries[row * $d + col];
                            rational(numerator, denominator)
                        })
                    });
                    let expected_determinant = rational_determinant_gaussian(rows.clone());
                    let matrix = RationalMatrix::try_from_rows(rows.clone()).unwrap();

                    prop_assert_eq!(matrix.det(), expected_determinant.clone());
                    prop_assert_eq!(matrix.det_sign(), determinant_sign(&expected_determinant));

                    if expected_determinant != zero() {
                        let expected_solution: [BigRational; $d] = from_fn(|index| {
                            let (numerator, denominator) = solution_entries[index];
                            rational(numerator, denominator)
                        });
                        let rhs = RationalVector::try_new(rational_matvec(&rows, &expected_solution)).unwrap();
                        let actual_solution = matrix.solve(&rhs).unwrap();

                        prop_assert_eq!(actual_solution.as_array(), &expected_solution);
                    }
                }

                #[test]
                fn [<rational_duplicate_row_is_singular_ $d d>](
                    entries in proptest::collection::vec((-5_i16..=5_i16, 1_u8..=9_u8), $d * $d),
                ) {
                    let mut rows: [[BigRational; $d]; $d] = from_fn(|row| {
                        from_fn(|col| {
                            let (numerator, denominator) = entries[row * $d + col];
                            rational(numerator, denominator)
                        })
                    });
                    rows[$d - 1] = rows[0].clone();
                    let matrix = RationalMatrix::try_from_rows(rows).unwrap();

                    prop_assert_eq!(matrix.det_sign(), DeterminantSign::Zero);
                    prop_assert_eq!(matrix.det(), zero());
                    prop_assert!(matches!(
                        matrix.solve(&RationalVector::zero()),
                        Err(LaError::Singular {
                            reason: SingularityReason::Exact,
                            ..
                        })
                    ), "duplicate rows must produce exact singularity");
                }
            }
        }
    };
}

gen_rational_properties!(2);
gen_rational_properties!(3);
gen_rational_properties!(4);
gen_rational_properties!(5);
gen_rational_properties!(6);
gen_rational_properties!(7);
gen_rational_properties!(8);
