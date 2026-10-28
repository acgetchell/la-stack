#![forbid(unsafe_code)]

//! Preserve cancellation and handle bounded exact-polynomial evaluation.

use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    let one = SparsePolynomial::try_from_terms([(0, BigRational::from_integer(1.into()))])?;
    let perturbed = SparsePolynomial::try_from_terms([
        (0, BigRational::from_integer(1.into())),
        (3, BigRational::new(1.into(), 2.into())),
    ])?;
    let matrix = [[one.clone(), one.clone()], [one, perturbed]];
    let budget = PolynomialDeterminantBudget {
        max_subset_transitions: 4,
        max_term_products: 5,
    };
    let determinant = polynomial_determinant(&matrix, budget)?;
    assert_eq!(
        determinant.leading_term(),
        Some((3, &BigRational::new(1.into(), 2.into())))
    );
    println!("Exact determinant terms: {:?}", determinant.terms());

    let exhausted = polynomial_determinant(
        &matrix,
        PolynomialDeterminantBudget {
            max_term_products: 4,
            ..budget
        },
    );
    match exhausted {
        Err(LaError::PolynomialBudgetExceeded {
            work: PolynomialWork::TermProducts,
            limit: 4,
            completed: 4,
            requested: 1,
            ..
        }) => {}
        other => panic!("expected the explicit budget failure, got {other:?}"),
    }
    Ok(())
}
