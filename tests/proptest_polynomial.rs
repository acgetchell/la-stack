#![forbid(unsafe_code)]

//! Determinant multilinearity for arbitrary canonical polynomial entries.

#![cfg(feature = "exact")]

use std::array::from_fn;

use la_stack::prelude::*;
use pastey::paste;
use proptest::prelude::*;

#[path = "common/proptest_config.rs"]
mod proptest_config;

const BUDGET: PolynomialDeterminantBudget = PolynomialDeterminantBudget {
    max_subset_transitions: 1_024,
    max_term_products: 1_000_000,
};

fn row_scaling<const D: usize>(terms: &[(u8, i16, i16)], scale: i16) {
    let mut matrix: [[SparsePolynomial; D]; D] = from_fn(|row| {
        from_fn(|col| {
            let start = 2 * (row * D + col);
            SparsePolynomial::try_from_terms(
                terms[start..start + 2]
                    .iter()
                    .map(|&(e, n, d)| (u128::from(e), BigRational::new(n.into(), d.into()))),
            )
            .unwrap()
        })
    });
    let determinant = polynomial_determinant(&matrix, BUDGET).unwrap();
    let factor = BigRational::from_integer(scale.into());
    for entry in &mut matrix[0] {
        *entry =
            SparsePolynomial::try_from_terms(entry.terms().iter().map(|(&e, c)| (e, c * &factor)))
                .unwrap();
    }
    let expected = SparsePolynomial::try_from_terms(
        determinant.terms().iter().map(|(&e, c)| (e, c * &factor)),
    )
    .unwrap();
    assert_eq!(polynomial_determinant(&matrix, BUDGET).unwrap(), expected);
}

macro_rules! dimension_properties {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(proptest_config::with_default_cases(32))]
                #[test]
                fn [<row_scaling_ $d>](
                    terms in prop::collection::vec((0_u8..6, -4_i16..5, 1_i16..7), 2 * $d * $d),
                    scale in -3_i16..4,
                ) {
                    row_scaling::<$d>(&terms, scale);
                }
            }
        }
    };
}
dimension_properties!(2);
dimension_properties!(3);
dimension_properties!(4);
dimension_properties!(5);
