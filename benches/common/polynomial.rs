#![forbid(unsafe_code)]

//! A consumer-shaped, maximally coincident lifted matrix and factorial oracle.

use std::array::from_fn;
use std::cmp::Ordering;
use std::collections::BTreeMap;

use la_stack::{
    BigRational, PolynomialDeterminantBudget, SparsePolynomial, polynomial_determinant,
};

use crate::bench_utils::OrAbort;

pub(crate) const BUDGET: PolynomialDeterminantBudget = PolynomialDeterminantBudget {
    max_subset_transitions: 1_024,
    max_term_products: 55_631_284,
};

/// A fixture checked against independent permutation enumeration before timing.
pub(crate) struct ValidatedPolynomialInput<const N: usize> {
    pub(crate) matrix: [[SparsePolynomial; N]; N],
}

/// All unperturbed points are at the origin. A coordinate is `t^(3^variable)`,
/// its lift is the sum of coordinate squares, and the last column is one.
/// Each permutation chooses one squared coordinate in its lifted row. Base-three
/// digits (0, 1, 2) identify that row, its axis and every other coordinate choice
/// uniquely, so the `N! * (N-2)` monomials are distinct and have coefficient ±1.
pub(crate) fn coincident_lifted<const N: usize>() -> ValidatedPolynomialInput<N> {
    assert!((3..=8).contains(&N));
    let axes = N - 2;
    let powers: [[u128; N]; N] = from_fn(|row| {
        from_fn(|col| {
            if col < axes {
                3_u128.pow(u32::try_from(row * axes + col).or_abort("bounded variable index"))
            } else {
                0
            }
        })
    });
    let matrix = from_fn(|row| {
        from_fn(|col| {
            let exponents = match col.cmp(&axes) {
                Ordering::Less => vec![powers[row][col]],
                Ordering::Equal => (0..axes).map(|axis| 2 * powers[row][axis]).collect(),
                Ordering::Greater => vec![0],
            };
            SparsePolynomial::try_from_terms(
                exponents
                    .into_iter()
                    .map(|e| (e, BigRational::from_integer(1.into()))),
            )
            .or_abort("integer polynomial coefficients")
        })
    });
    let mut expected = BTreeMap::new();
    enumerate(&powers, &mut from_fn(|i| i), 0, &mut expected);
    assert_eq!(expected.len(), (1..=N).product::<usize>() * axes);
    let determinant =
        polynomial_determinant(&matrix, BUDGET).or_abort("lifted determinant validation");
    assert_eq!(determinant.terms(), &expected);
    ValidatedPolynomialInput { matrix }
}

fn enumerate<const N: usize>(
    powers: &[[u128; N]; N],
    columns: &mut [usize; N],
    row: usize,
    result: &mut BTreeMap<u128, BigRational>,
) {
    if row < N {
        for col in row..N {
            columns.swap(row, col);
            enumerate(powers, columns, row + 1, result);
            columns.swap(row, col);
        }
        return;
    }
    let inversions: usize = (0..N)
        .map(|i| (i + 1..N).filter(|&j| columns[i] > columns[j]).count())
        .sum();
    for axis in 0..N - 2 {
        let exponent = columns
            .iter()
            .enumerate()
            .map(|(i, &j)| {
                if j == N - 2 {
                    2 * powers[i][axis]
                } else {
                    powers[i][j]
                }
            })
            .sum();
        let coefficient =
            BigRational::from_integer(if inversions.is_multiple_of(2) { 1 } else { -1 }.into());
        assert_eq!(
            result.insert(exponent, coefficient),
            None,
            "unique base-three encoding"
        );
    }
}
