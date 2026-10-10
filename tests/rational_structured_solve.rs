//! Exact solve regressions for triangular recognition and the general fallback.

#![cfg(feature = "exact")]
#![forbid(unsafe_code)]

use std::array::from_fn;

use la_stack::{BigInt, BigRational, LaError, RationalMatrix, RationalVector};
use num_traits::Zero;
use pastey::paste;
use proptest::{collection, prelude::*};

#[path = "../benches/common/bench_utils.rs"]
mod bench_utils;
#[path = "common/proptest_config.rs"]
mod proptest_config;
#[path = "../benches/common/rational_solve.rs"]
pub mod rational_solve;

use proptest_config::with_default_cases;
use rational_solve::legacy_solve;

fn matvec<const D: usize>(rows: &[[BigRational; D]; D], x: &[BigRational; D]) -> [BigRational; D] {
    from_fn(|row| rows[row].iter().zip(x).map(|(a, x)| a * x).sum())
}

fn check_solution<const D: usize>(rows: &[[BigRational; D]; D], expected: &[BigRational; D]) {
    let matrix = RationalMatrix::try_from_rows(rows.clone()).unwrap();
    for x in [expected, &from_fn(|_| BigRational::zero())] {
        let b = matvec(rows, x);
        let reference = legacy_solve(
            rows.clone().into_iter().map(Vec::from).collect(),
            b.to_vec(),
        )
        .unwrap();
        assert_eq!(reference, x);
        let rhs = RationalVector::try_new(b.clone()).unwrap();
        let solution = matrix.solve(&rhs).unwrap();
        assert_eq!(solution.as_array(), x);
        assert_eq!(solution.as_array().as_slice(), reference);
        assert_eq!(matvec(rows, solution.as_array()), b);
        // Both inputs remain borrowed, including across repeated solves.
        assert_eq!(matrix.as_rows(), rows);
        assert_eq!(rhs.as_array(), &b);
    }
}

fn check_wide_and_fallback<const D: usize>() {
    let expected: [BigRational; D] =
        from_fn(|col| BigRational::new(BigInt::from(col + 1), (BigInt::from(1) << 257) - 1));
    for lower in [false, true] {
        let mut rows = from_fn(|row| {
            from_fn(|col| {
                if (lower && col > row) || (!lower && col < row) {
                    return BigRational::zero();
                }
                let sign = if (row + col).is_multiple_of(2) { 1 } else { -1 };
                BigRational::new(
                    ((BigInt::from(1) << (256 + row)) + BigInt::from(col + 1)) * sign,
                    (BigInt::from(1) << (300 + col)) - BigInt::from(2 * row + 1),
                )
            })
        });
        for _ in 0..D.max(1) {
            check_solution(&rows, &expected);
            if D > 0 {
                rows.rotate_left(1);
            }
        }
        rows.reverse();
        check_solution(&rows, &expected);
    }
    // Every row has two edge coefficients for D>=2, so neither triangular
    // recognition can succeed. Strict diagonal dominance proves invertibility.
    let rows = from_fn(|row| {
        from_fn(|col| {
            let numerator = if row == col { 2 * D + 1 } else { 1 };
            BigRational::new(BigInt::from(numerator), BigInt::from(2 * row + 3))
        })
    });
    check_solution(&rows, &expected);
}

fn check_singular<const D: usize>() {
    if D == 0 {
        assert_eq!(
            RationalMatrix::<D>::zero().solve(&RationalVector::zero()),
            Ok(RationalVector::zero())
        );
    }
    for pivot in 0..D {
        let mut rows: [[BigRational; D]; D] = from_fn(|row| {
            from_fn(|col| {
                if row == col && row != pivot {
                    BigRational::from_integer(1.into())
                } else {
                    BigRational::zero()
                }
            })
        });
        for _ in 0..D {
            let matrix = RationalMatrix::try_from_rows(rows.clone()).unwrap();
            assert_eq!(
                matrix.solve(&RationalVector::zero()),
                Err(LaError::singular_exact(pivot))
            );
            rows.rotate_left(1);
        }
    }
    if D >= 2 {
        let mut rows: [[BigRational; D]; D] = from_fn(|row| {
            from_fn(|col| BigRational::from_integer(BigInt::from(usize::from(row == col))))
        });
        let (first, remaining) = rows.split_at_mut(1);
        remaining[D - 2].clone_from(&first[0]);
        for _ in 0..D {
            let matrix = RationalMatrix::try_from_rows(rows.clone()).unwrap();
            assert_eq!(
                matrix.solve(&RationalVector::zero()),
                Err(LaError::singular_exact(D - 1))
            );
            rows.rotate_left(1);
        }
    }
}

fn check_random_triangular<const D: usize>(
    entries: &[(i16, u8)],
    components: &[(i16, u8)],
    permutation_keys: &[u64],
    lower: bool,
) {
    let expected = from_fn(|col| {
        BigRational::new(
            BigInt::from(components[col].0),
            BigInt::from(components[col].1),
        )
    });
    let mut order: [usize; D] = from_fn(|row| row);
    order.sort_by_key(|&row| permutation_keys[row]);
    let rows = from_fn(|index| {
        let row = order[index];
        from_fn(|col| {
            if (lower && col > row) || (!lower && col < row) {
                return BigRational::zero();
            }
            let (mut numerator, denominator) = entries[row * D + col];
            if row == col && numerator == 0 {
                numerator = 1;
            }
            BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
        })
    });
    check_solution::<D>(&rows, &expected);
}

macro_rules! dimension_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<wide_structured_and_fallback_ $d d>]() { check_wide_and_fallback::<$d>(); }

            #[test]
            fn [<singular_metadata_ $d d>]() { check_singular::<$d>(); }

            proptest! {
                #![proptest_config(with_default_cases(24))]
                #[test]
                fn [<permuted_triangular_ $d d>](
                    entries in collection::vec((-9_i16..=9, 1_u8..=17), $d * $d),
                    components in collection::vec((-9_i16..=9, 1_u8..=17), $d),
                    permutation_keys in collection::vec(any::<u64>(), $d),
                    lower in any::<bool>(),
                ) {
                    check_random_triangular::<$d>(&entries, &components, &permutation_keys, lower);
                }
            }
        }
    };
}

dimension_tests!(0);
dimension_tests!(1);
dimension_tests!(2);
dimension_tests!(3);
dimension_tests!(4);
dimension_tests!(5);
dimension_tests!(6);
dimension_tests!(7);
dimension_tests!(8);

#[test]
fn rank_deficiency_with_nonzero_rows_preserves_pivot_metadata() {
    for (entries, pivot) in [
        ([[1, 2, 3], [0, 0, 4], [0, 0, 5]], 1),
        ([[1, 0, 0], [2, 0, 0], [3, 4, 5]], 2),
        ([[1, 2, 3], [2, 4, 6], [4, 5, 6]], 2),
    ] {
        let mut rows = entries.map(|row| row.map(|value| BigRational::from_integer(value.into())));
        for _ in 0..3 {
            let matrix = RationalMatrix::try_from_rows(rows.clone()).unwrap();
            assert_eq!(
                matrix.solve(&RationalVector::zero()),
                Err(LaError::singular_exact(pivot))
            );
            rows.rotate_left(1);
        }
    }
}
