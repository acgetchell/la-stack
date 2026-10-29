#![forbid(unsafe_code)]

//! Independent permutation expansions and public polynomial contracts.

#![cfg(feature = "exact")]

use core::assert_matches;
use std::array::from_fn;
use std::collections::BTreeMap;

use la_stack::prelude::*;
use num_traits::{One, Zero};
use pastey::paste;

const BUDGET: PolynomialDeterminantBudget = PolynomialDeterminantBudget {
    max_subset_transitions: 1_024,
    max_term_products: 1_000_000,
};

fn rational(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(numerator.into(), denominator.into())
}

fn polynomial(terms: &[(u128, i64)]) -> SparsePolynomial {
    SparsePolynomial::try_from_terms(terms.iter().map(|&(e, c)| (e, rational(c, 1))))
        .expect("integer coefficients")
}

// A dense coefficient vector and factorial permutation enumeration deliberately
// avoid the production sparse convolution and column-subset recurrence.
fn permutation_determinant<const D: usize>(
    matrix: &[[Vec<BigRational>; D]; D],
) -> Vec<BigRational> {
    fn visit<const D: usize>(
        matrix: &[[Vec<BigRational>; D]; D],
        columns: &mut [usize; D],
        row: usize,
        result: &mut Vec<BigRational>,
    ) {
        if row < D {
            for next in row..D {
                columns.swap(row, next);
                visit(matrix, columns, row + 1, result);
                columns.swap(row, next);
            }
            return;
        }
        let inversions: usize = (0..D)
            .map(|i| (i + 1..D).filter(|&j| columns[i] > columns[j]).count())
            .sum();
        let mut product = vec![rational(
            if inversions.is_multiple_of(2) { 1 } else { -1 },
            1,
        )];
        for (i, &j) in columns.iter().enumerate() {
            let entry = &matrix[i][j];
            let mut next = vec![BigRational::zero(); product.len() + entry.len() - 1];
            for (a, x) in product.iter().enumerate() {
                for (b, y) in entry.iter().enumerate() {
                    next[a + b] += x * y;
                }
            }
            product = next;
        }
        result.resize(result.len().max(product.len()), BigRational::zero());
        for (target, coefficient) in result.iter_mut().zip(product) {
            *target += coefficient;
        }
    }
    let mut result = Vec::new();
    visit(matrix, &mut from_fn(|i| i), 0, &mut result);
    result
}

fn dense_terms(coefficients: &[BigRational]) -> BTreeMap<u128, BigRational> {
    coefficients
        .iter()
        .enumerate()
        .filter(|(_, coefficient)| !coefficient.is_zero())
        .map(|(exponent, coefficient)| (exponent as u128, coefficient.clone()))
        .collect()
}

fn check_permutations<const D: usize>() {
    for seed in 0..6 {
        let dense: [[Vec<BigRational>; D]; D] = from_fn(|row| {
            from_fn(|col| {
                (0..3)
                    .map(|e| {
                        let value =
                            (row * 7 + col * 11 + e * 13 + row * col * (seed + 1) + seed) % 9;
                        rational(
                            i64::try_from(value).unwrap() - 4,
                            i64::try_from(e + col + 1).unwrap(),
                        )
                    })
                    .collect()
            })
        });
        let mut matrix = dense.clone().map(|row| {
            row.map(|coefficients| {
                SparsePolynomial::try_from_terms(
                    coefficients
                        .into_iter()
                        .enumerate()
                        .map(|(e, c)| (e as u128, c)),
                )
                .unwrap()
            })
        });
        let expected = dense_terms(&permutation_determinant(&dense));
        let result = polynomial_determinant(&matrix, BUDGET).unwrap();
        assert_eq!(result.terms(), &expected, "order {D}, seed {seed}");
        matrix.swap(0, D - 1);
        let negative: BTreeMap<_, _> = expected.iter().map(|(&e, c)| (e, -c)).collect();
        assert_eq!(
            polynomial_determinant(&matrix, BUDGET).unwrap().terms(),
            &negative
        );
        let (first, rest) = matrix.split_at_mut(1);
        rest[0].clone_from(&first[0]);
        assert_eq!(
            polynomial_determinant(&matrix, BUDGET).unwrap(),
            SparsePolynomial::zero()
        );
    }
}

// det(J + diag(t^e_i)) = product(t^e_i) + sum_i product_{j != i}(t^e_j).
// Multilinearity leaves either all diagonal columns or exactly one J column;
// two J columns are identical and contribute zero. This checks every coefficient
// with independent analytical evidence through the largest supported order.
fn check_degenerate<const D: usize>() {
    let exponents: [u128; D] = from_fn(|i| 3_u128.pow(u32::try_from(i).unwrap()));
    let matrix: [[SparsePolynomial; D]; D] = from_fn(|row| {
        from_fn(|col| {
            if row == col {
                polynomial(&[(0, 1), (exponents[row], 1)])
            } else {
                polynomial(&[(0, 1)])
            }
        })
    });
    let total: u128 = exponents.iter().sum();
    let expected: BTreeMap<_, _> = std::iter::once((total, BigRational::one()))
        .chain(exponents.iter().map(|e| (total - e, BigRational::one())))
        .collect();
    let result = polynomial_determinant(&matrix, BUDGET).unwrap();
    assert_eq!(result.terms(), &expected);
    assert_eq!(result.leading_term().unwrap().0, total - exponents[D - 1]);
}

macro_rules! dimension_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<permutation_oracle_ $d>]() { check_permutations::<$d>(); }
            #[test]
            fn [<degenerate_ $d>]() { check_degenerate::<$d>(); }
        }
    };
}
dimension_tests!(2);
dimension_tests!(3);
dimension_tests!(4);
dimension_tests!(5);

#[test]
fn degenerate_largest_orders() {
    check_degenerate::<6>();
    check_degenerate::<7>();
    check_degenerate::<8>();
}

#[test]
fn canonicalization_merges_cancels_and_orders() {
    const ZERO: SparsePolynomial = SparsePolynomial::zero();
    let result = SparsePolynomial::try_from_terms([
        (7, BigRational::new_raw(2.into(), (-4).into())),
        (2, rational(3, 5)),
        (7, rational(1, 2)),
        (1, BigRational::new_raw(0.into(), (-3).into())),
        (2, rational(-1, 10)),
        (9, BigRational::new_raw((-6).into(), (-8).into())),
    ])
    .unwrap();
    assert_eq!(
        result.terms(),
        &BTreeMap::from([(2, rational(1, 2)), (9, rational(3, 4))])
    );
    assert_eq!(result.leading_term(), Some((2, &rational(1, 2))));
    assert!(!result.is_zero());
    assert_eq!(ZERO, SparsePolynomial::default());
    assert!(ZERO.is_zero());
    assert_eq!(ZERO.leading_term(), None);
}

#[test]
fn invalid_raw_coefficient_is_not_hidden_by_zero_removal() {
    for numerator in [0, 1] {
        let result = SparsePolynomial::try_from_terms([
            (4, rational(1, 1)),
            (7, BigRational::new_raw(numerator.into(), 0.into())),
        ]);
        assert_matches!(
            result,
            Err(LaError::NonFinite {
                location: NonFiniteLocation::PolynomialTerm {
                    index: 1,
                    exponent: 7,
                    ..
                },
                origin: NonFiniteOrigin::Input,
                ..
            })
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "non-finite input value at polynomial term 1 (exponent 7)"
        );
    }
}

#[test]
fn empty_zero_and_unsupported_orders_are_distinct() {
    let none = PolynomialDeterminantBudget {
        max_subset_transitions: 0,
        max_term_products: 0,
    };
    assert_eq!(
        polynomial_determinant::<0>(&[], none).unwrap(),
        polynomial(&[(0, 1)])
    );
    let zero = [[SparsePolynomial::zero()]];
    assert_eq!(
        polynomial_determinant(
            &zero,
            PolynomialDeterminantBudget {
                max_subset_transitions: 1,
                ..none
            }
        )
        .unwrap()
        .leading_term(),
        None
    );
    let too_large: [[SparsePolynomial; 9]; 9] = from_fn(|_| from_fn(|_| SparsePolynomial::zero()));
    assert_eq!(
        polynomial_determinant(&too_large, none),
        Err(LaError::unsupported_dimension(9, MAX_POLYNOMIAL_MATRIX_DIM))
    );
}

#[test]
fn early_coefficients_cancel_before_leading_term_is_selected() {
    // (1+t)(1-t+t^2+t^5) - 1 = t^3+t^5+t^6.
    let matrix = [
        [polynomial(&[(0, 1), (1, 1)]), polynomial(&[(0, 1)])],
        [
            polynomial(&[(0, 1)]),
            polynomial(&[(0, 1), (1, -1), (2, 1), (5, 1)]),
        ],
    ];
    assert_eq!(
        polynomial_determinant(&matrix, BUDGET).unwrap(),
        polynomial(&[(3, 1), (5, 1), (6, 1)])
    );
}

#[test]
fn budgets_are_inclusive_deterministic_and_do_not_mutate_inputs() {
    let one = polynomial(&[(0, 1)]);
    let matrix = [
        [one.clone(), one.clone()],
        [one, polynomial(&[(0, 1), (3, 1)])],
    ];
    let original = matrix.clone();
    let exact = PolynomialDeterminantBudget {
        max_subset_transitions: 4,
        max_term_products: 5,
    };
    for (budget, expected) in [
        (
            PolynomialDeterminantBudget {
                max_subset_transitions: 3,
                ..exact
            },
            (PolynomialWork::SubsetTransitions, 3, 0, 4),
        ),
        (
            PolynomialDeterminantBudget {
                max_term_products: 4,
                ..exact
            },
            (PolynomialWork::TermProducts, 4, 4, 1),
        ),
        (
            PolynomialDeterminantBudget {
                max_term_products: 0,
                ..exact
            },
            (PolynomialWork::TermProducts, 0, 0, 1),
        ),
    ] {
        match polynomial_determinant(&matrix, budget) {
            Err(LaError::PolynomialBudgetExceeded {
                work,
                limit,
                completed,
                requested,
                ..
            }) => {
                assert_eq!((work, limit, completed, requested), expected);
            }
            other => panic!("expected budget exhaustion, got {other:?}"),
        }
        assert_eq!(matrix, original);
    }
    assert_eq!(
        polynomial_determinant(&matrix, exact).unwrap(),
        polynomial(&[(3, 1)])
    );
    let zero: [[SparsePolynomial; 8]; 8] = from_fn(|_| from_fn(|_| SparsePolynomial::zero()));
    assert_matches!(
        polynomial_determinant(
            &zero,
            PolynomialDeterminantBudget {
                max_subset_transitions: 1_023,
                max_term_products: usize::MAX,
            }
        ),
        Err(LaError::PolynomialBudgetExceeded {
            work: PolynomialWork::SubsetTransitions,
            limit: 1_023,
            completed: 0,
            requested: 1_024,
            ..
        })
    );
    assert!(
        polynomial_determinant(
            &zero,
            PolynomialDeterminantBudget {
                max_subset_transitions: 1_024,
                max_term_products: 0,
            }
        )
        .unwrap()
        .is_zero()
    );
}

#[test]
fn exponent_boundaries_and_error_precedence() {
    let largest = polynomial(&[(u128::MAX, 1)]);
    assert_eq!(
        polynomial_determinant(&[[largest.clone()]], BUDGET).unwrap(),
        largest
    );
    let matrix = [
        [largest, SparsePolynomial::zero()],
        [SparsePolynomial::zero(), polynomial(&[(1, 1)])],
    ];
    assert_matches!(
        polynomial_determinant(&matrix, BUDGET),
        Err(LaError::PolynomialExponentOverflow {
            row: 1,
            col: 1,
            left: u128::MAX,
            right: 1,
            ..
        })
    );
    assert_matches!(
        polynomial_determinant(
            &matrix,
            PolynomialDeterminantBudget {
                max_term_products: 1,
                ..BUDGET
            }
        ),
        Err(LaError::PolynomialBudgetExceeded {
            work: PolynomialWork::TermProducts,
            limit: 1,
            completed: 1,
            requested: 1,
            ..
        })
    );
    let boundary = [
        [polynomial(&[(u128::MAX - 1, 1)]), SparsePolynomial::zero()],
        [SparsePolynomial::zero(), polynomial(&[(1, 1)])],
    ];
    assert_eq!(
        polynomial_determinant(&boundary, BUDGET).unwrap(),
        polynomial(&[(u128::MAX, 1)])
    );
}

#[test]
fn coefficients_remain_exact_beyond_binary64_range() {
    let wide = BigInt::one() << 2048;
    let a = BigRational::new(&wide + 1, 3.into());
    let b = BigRational::new(5.into(), &wide - 1);
    let matrix = [
        [
            SparsePolynomial::try_from_terms([(0, a.clone()), (7, b.clone())]).unwrap(),
            SparsePolynomial::zero(),
        ],
        [
            SparsePolynomial::zero(),
            SparsePolynomial::try_from_terms([(0, a.clone()), (7, -b.clone())]).unwrap(),
        ],
    ];
    assert_eq!(
        polynomial_determinant(&matrix, BUDGET).unwrap().terms(),
        &BTreeMap::from([(0, &a * &a), (14, -(&b * &b))])
    );
}
