#![forbid(unsafe_code)]

//! Bounded exact determinants over sparse univariate rational polynomials.
//!
//! The subset recurrence reorganizes the Leibniz identity (REFERENCES.md \[12\]).
//! See the [invariant and resource contract](https://github.com/acgetchell/la-stack/blob/main/docs/mathematical_basis.md#exact-sparse-polynomial-determinants).

use std::collections::{BTreeMap, btree_map::Entry};

use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::rational::canonicalize_rational;
use crate::{LaError, NonFiniteLocation, NonFiniteOrigin, PolynomialWork};

/// Largest supported order of [`polynomial_determinant`].
pub const MAX_POLYNOMIAL_MATRIX_DIM: usize = 8;

/// An exact polynomial `Σ coefficient[e] * t^e` with nonnegative `u128` exponents.
///
/// Terms are stored in increasing exponent order. Coefficients are reduced
/// rational numbers with positive denominators; zero coefficients are absent.
/// Duplicate input exponents are summed exactly, including cancellation.
/// [`Self::leading_term`] means the *least* exponent with a nonzero coefficient,
/// appropriate to an infinitesimal `t → 0+`, rather than the highest degree.
///
/// Storage and arbitrary-precision coefficients allocate on the heap. Construction
/// consumes owned terms; queries borrow them; cloning makes a deep copy. The
/// determinant budget does not cover construction. Multivariate encodings,
/// perturbation-variable assignment and interpretation of coefficient signs
/// belong to callers, who must prove their encoding preserves the desired order
/// and exponent addition. No multivariate ordering is inferred here.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
#[must_use]
pub struct SparsePolynomial {
    terms: BTreeMap<u128, BigRational>,
}

impl SparsePolynomial {
    /// Construct the zero polynomial without allocating.
    pub const fn zero() -> Self {
        Self {
            terms: BTreeMap::new(),
        }
    }

    /// Consume terms, normalize rational coefficients and combine equal exponents.
    ///
    /// Input order is immaterial to a successful result. Raw unreduced or
    /// negative-denominator rationals are accepted. Validation precedes zero
    /// removal, so even `0/0` is rejected. The iterator is consumed only up to
    /// the first invalid coefficient; this constructor has no work budget.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with input origin and a
    /// [`NonFiniteLocation::PolynomialTerm`] for a zero denominator.
    pub fn try_from_terms(
        terms: impl IntoIterator<Item = (u128, BigRational)>,
    ) -> Result<Self, LaError> {
        let mut polynomial = Self::zero();
        for (index, (exponent, coefficient)) in terms.into_iter().enumerate() {
            if coefficient.denom().is_zero() {
                return Err(LaError::NonFinite {
                    location: NonFiniteLocation::PolynomialTerm { index, exponent },
                    origin: NonFiniteOrigin::Input,
                });
            }
            polynomial.add_term(exponent, canonicalize_rational(coefficient));
        }
        Ok(polynomial)
    }

    /// Borrow all nonzero terms in increasing exponent order.
    #[must_use]
    pub const fn terms(&self) -> &BTreeMap<u128, BigRational> {
        &self.terms
    }

    /// Whether every coefficient is exactly zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Borrow the first nonzero coefficient in increasing exponent order.
    ///
    /// `None` proves this is the zero polynomial. A failed determinant evaluation
    /// never returns a polynomial, so it cannot be confused with this result.
    #[must_use]
    pub fn leading_term(&self) -> Option<(u128, &BigRational)> {
        self.terms
            .first_key_value()
            .map(|(&exponent, coefficient)| (exponent, coefficient))
    }

    fn add_term(&mut self, exponent: u128, coefficient: BigRational) {
        if coefficient.is_zero() {
            return;
        }
        match self.terms.entry(exponent) {
            Entry::Vacant(entry) => {
                entry.insert(coefficient);
            }
            Entry::Occupied(mut entry) => {
                *entry.get_mut() += coefficient;
                if entry.get().is_zero() {
                    entry.remove();
                }
            }
        }
    }
}

/// Caller-selected work limits for one [`polynomial_determinant`] invocation.
///
/// No default policy is supplied. All limits, including zero, are valid.
/// These count combinatorial work, not elapsed time, bytes or integer bit
/// operations: a single rational product can be expensive for wide inputs.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[must_use]
pub struct PolynomialDeterminantBudget {
    /// Maximum structural subset transitions. The complete `D * 2^(D-1)`
    /// envelope is checked before allocation, even for an all-zero matrix.
    /// Order zero requires zero transitions; order eight requires 1,024.
    pub max_subset_transitions: usize,
    /// Maximum individual nonzero coefficient-pair multiplications. Each is
    /// charged before exponent addition or rational arithmetic; zero partials
    /// and zero entries cost no products. Cancellation does not refund work.
    pub max_term_products: usize,
}

/// Compute the complete exact determinant polynomial of a borrowed square array.
///
/// Supports orders `0..=8`; the empty determinant is the constant one. Singular
/// polynomial matrices return the zero polynomial successfully. All cancellation
/// terms are retained, so [`SparsePolynomial::leading_term`] is valid only after
/// the complete calculation succeeds. There is no early sign fallback, truncation,
/// division, rounding or conditioning-dependent numerical error.
///
/// The kernel allocates `2^D` sparse states (at most 256) and consumes each partial
/// once all its predecessors have contributed. Inputs are neither cloned nor
/// changed. At most one new stored term is created per charged product, plus the
/// initial constant one. Rational sizes, allocator overhead, and caller-owned
/// inputs are outside these combinatorial bounds; allocation failure follows
/// Rust's allocator behavior. Limits do not promise a byte or wall-clock bound.
///
/// # Errors
/// In order: [`LaError::UnsupportedDimension`] for `D > 8`,
/// [`LaError::PolynomialBudgetExceeded`] for insufficient structural transitions,
/// then product-budget or [`LaError::PolynomialExponentOverflow`] failures during
/// evaluation. An overflowing intermediate exponent is an error even if later
/// contributions could cancel it. Errors discard partial results and leave the
/// input reusable. Traversal is increasing subset mask, column, left exponent,
/// then right exponent, making budget accounting deterministic.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
/// # fn main() -> Result<(), LaError> {
/// let one = SparsePolynomial::try_from_terms([(0, BigRational::from_integer(1.into()))])?;
/// let perturbed = SparsePolynomial::try_from_terms([
///     (0, BigRational::from_integer(1.into())),
///     (3, BigRational::new(1.into(), 2.into())),
/// ])?;
/// let matrix = [[one.clone(), one.clone()], [one, perturbed]];
/// let determinant = polynomial_determinant(&matrix, PolynomialDeterminantBudget {
///     max_subset_transitions: 4,
///     max_term_products: 5,
/// })?;
/// assert_eq!(determinant.leading_term(), Some((3, &BigRational::new(1.into(), 2.into()))));
/// # Ok(())
/// # }
/// ```
pub fn polynomial_determinant<const D: usize>(
    matrix: &[[SparsePolynomial; D]; D],
    budget: PolynomialDeterminantBudget,
) -> Result<SparsePolynomial, LaError> {
    if D > MAX_POLYNOMIAL_MATRIX_DIM {
        return Err(LaError::unsupported_dimension(D, MAX_POLYNOMIAL_MATRIX_DIM));
    }
    let state_count = 1usize << D;
    let transitions = D * (state_count / 2);
    if transitions > budget.max_subset_transitions {
        return Err(LaError::PolynomialBudgetExceeded {
            work: PolynomialWork::SubsetTransitions,
            limit: budget.max_subset_transitions,
            completed: 0,
            requested: transitions,
        });
    }
    let mut partials = vec![SparsePolynomial::zero(); state_count];
    partials[0].add_term(0, BigRational::one());
    let mut products = 0;
    for mask in 0..state_count - 1 {
        let row = mask.count_ones() as usize;
        let partial = std::mem::take(&mut partials[mask]);
        if partial.is_zero() {
            continue;
        }
        for (col, entry) in matrix[row].iter().enumerate() {
            if mask & (1 << col) != 0 {
                continue;
            }
            let negative = (mask >> (col + 1)).count_ones() % 2 != 0;
            let next = &mut partials[mask | (1 << col)];
            for (&left, left_coefficient) in &partial.terms {
                for (&right, right_coefficient) in &entry.terms {
                    if products == budget.max_term_products {
                        return Err(LaError::PolynomialBudgetExceeded {
                            work: PolynomialWork::TermProducts,
                            limit: budget.max_term_products,
                            completed: products,
                            requested: 1,
                        });
                    }
                    products += 1;
                    let exponent =
                        left.checked_add(right)
                            .ok_or(LaError::PolynomialExponentOverflow {
                                row,
                                col,
                                left,
                                right,
                            })?;
                    let coefficient = left_coefficient * right_coefficient;
                    next.add_term(exponent, if negative { -coefficient } else { coefficient });
                }
            }
        }
    }
    Ok(std::mem::take(&mut partials[state_count - 1]))
}
