#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![doc = include_str!("../README.md")]

#[cfg(doc)]
mod readme_doctests {
    //! Executable versions of README examples.
    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// // This system requires pivoting (a[0][0] = 0), so it's a good LU demo.
    /// let a = Matrix::<5>::from_rows([
    ///     [0.0, 1.0, 1.0, 1.0, 1.0],
    ///     [1.0, 0.0, 1.0, 1.0, 1.0],
    ///     [1.0, 1.0, 0.0, 1.0, 1.0],
    ///     [1.0, 1.0, 1.0, 0.0, 1.0],
    ///     [1.0, 1.0, 1.0, 1.0, 0.0],
    /// ]);
    ///
    /// let b = Vector::<5>::new([14.0, 13.0, 12.0, 11.0, 10.0]);
    ///
    /// let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();
    /// let x = lu.solve_vec(b).unwrap().into_array();
    ///
    /// // Floating-point rounding is expected; compare with a tolerance.
    /// let expected = [1.0, 2.0, 3.0, 4.0, 5.0];
    /// for (x_i, e_i) in x.iter().zip(expected.iter()) {
    ///     assert!((*x_i - *e_i).abs() <= 1e-12);
    /// }
    /// ```
    fn solve_5x5_example() {}

    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// // This matrix is symmetric positive-definite (A = L*L^T) so LDLT works without pivoting.
    /// let a = Matrix::<5>::from_rows([
    ///     [1.0, 1.0, 0.0, 0.0, 0.0],
    ///     [1.0, 2.0, 1.0, 0.0, 0.0],
    ///     [0.0, 1.0, 2.0, 1.0, 0.0],
    ///     [0.0, 0.0, 1.0, 2.0, 1.0],
    ///     [0.0, 0.0, 0.0, 1.0, 2.0],
    /// ]);
    ///
    /// let det = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap().det();
    /// assert!((det - 1.0).abs() <= 1e-12);
    /// ```
    fn det_5x5_ldlt_example() {}
}

#[cfg(feature = "exact")]
mod exact;
#[cfg(feature = "exact")]
pub use num_bigint::BigInt;
#[cfg(feature = "exact")]
pub use num_rational::BigRational;
#[cfg(feature = "exact")]
pub use num_traits::{FromPrimitive, Signed, ToPrimitive};

mod ldlt;
mod lu;
mod matrix;
mod vector;

use core::fmt;

// ---------------------------------------------------------------------------
// Error-bound constants for `Matrix::det_errbound()`.
//
// For `D ∈ {2, 3, 4}`, `Matrix::det_direct()` evaluates the Leibniz expansion
// of the determinant as a tree of f64 multiplies and fused multiply-adds
// (FMAs).  Following Shewchuk's error-analysis methodology (REFERENCES.md
// [8]), the absolute error of that computation is bounded by
//
//     |det_direct(A) - det_exact(A)|  ≤  ERR_COEFF_D · p(|A|)
//
// where `p(|A|)` is the **absolute Leibniz sum**
//
//     p(|A|) = Σ_σ ∏ᵢ |A[i, σ(i)]|,
//
// i.e. the same cofactor-expansion tree as `det_direct` but with each
// entry replaced by its magnitude.  Note that `p(|A|)` is *not* the
// combinatorial matrix permanent — the name "permanent" appears in the
// source for brevity and to match the cited literature.
//
// Each constant has the shape `a · EPS + b · EPS²`: the linear term bounds
// the first-order rounding and the quadratic term absorbs the interaction
// of errors in nested FMAs.  The coefficients `a` and `b` are conservative
// over-estimates derived from the longest dependency chain of `det_direct`
// at that dimension.
//
// These constants are NOT feature-gated — they rely only on f64 arithmetic
// and are useful for adaptive-precision logic even without the `exact`
// feature.  Most callers should prefer `Matrix::det_errbound()`, which
// applies these constants to the actual matrix; the raw constants are
// exposed for advanced use cases (composing the bound with a pre-reduced
// permanent, rolling a custom adaptive filter, etc.).  See
// `Matrix::det_sign_exact()` (behind the `exact` feature) for the
// reference adaptive-filter that consumes these internally.
// ---------------------------------------------------------------------------

const EPS: f64 = f64::EPSILON; // 2^-52

/// Absolute error coefficient for [`Matrix::<2>::det_direct`](crate::Matrix::det_direct).
///
/// For any 2×2 matrix `A = [[a, b], [c, d]]` with finite f64 entries,
///
/// ```text
/// |A.det_direct() - det_exact(A)|  ≤  ERR_COEFF_2 · (|a·d| + |b·c|)
/// ```
///
/// `det_direct` evaluates `a·d - b·c` as one multiply followed by one FMA
/// (2 rounding events); the linear `3·EPS` term bounds those roundings
/// and the quadratic `16·EPS²` term is a conservative cushion for their
/// interaction.  Derivation follows Shewchuk's framework; see
/// `REFERENCES.md` \[8\].
///
/// Prefer [`Matrix::det_errbound`](crate::Matrix::det_errbound) unless
/// you already have the absolute-Leibniz sum available; see
/// [`Matrix::det_sign_exact`](crate::Matrix::det_sign_exact) (under the
/// `exact` feature) for the reference adaptive-precision filter.
///
/// # Example
/// ```
/// use la_stack::prelude::*;
///
/// let m = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
/// let det = m.det_direct().unwrap();
/// // Compute the bound from the raw constant for illustration; most
/// // callers would just use `m.det_errbound().unwrap()` instead.
/// let p = (1.0_f64 * 4.0).abs() + (2.0_f64 * 3.0).abs();
/// let bound = ERR_COEFF_2 * p;
/// if det.abs() > bound {
///     // The f64 sign is provably correct without exact arithmetic.
/// }
/// ```
pub const ERR_COEFF_2: f64 = 3.0 * EPS + 16.0 * EPS * EPS;

/// Absolute error coefficient for [`Matrix::<3>::det_direct`](crate::Matrix::det_direct).
///
/// For any 3×3 matrix `A` with finite f64 entries,
///
/// ```text
/// |A.det_direct() - det_exact(A)|  ≤  ERR_COEFF_3 · p(|A|)
/// ```
///
/// where `p(|A|)` is the absolute Leibniz sum (the same cofactor
/// expansion as `det_direct` but with `|·|` at every leaf).
/// `det_direct` for D=3 uses three 2×2 FMA minors combined by a nested
/// FMA, yielding the `8·EPS + 64·EPS²` bound.  See `REFERENCES.md`
/// \[8\] for the Shewchuk framework these bounds follow.
///
/// Prefer [`Matrix::det_errbound`](crate::Matrix::det_errbound) over this
/// constant for typical use; see [`ERR_COEFF_2`] for a worked example.
pub const ERR_COEFF_3: f64 = 8.0 * EPS + 64.0 * EPS * EPS;

/// Absolute error coefficient for [`Matrix::<4>::det_direct`](crate::Matrix::det_direct).
///
/// For any 4×4 matrix `A` with finite f64 entries,
///
/// ```text
/// |A.det_direct() - det_exact(A)|  ≤  ERR_COEFF_4 · p(|A|)
/// ```
///
/// where `p(|A|)` is the absolute Leibniz sum.  `det_direct` for D=4
/// hoists six 2×2 minors, combines them into four 3×3 cofactors, then
/// reduces those with an FMA row combination, yielding the
/// `12·EPS + 128·EPS²` bound.  See `REFERENCES.md` \[8\] for the
/// Shewchuk framework these bounds follow.
///
/// Prefer [`Matrix::det_errbound`](crate::Matrix::det_errbound) over this
/// constant for typical use; see [`ERR_COEFF_2`] for a worked example.
pub const ERR_COEFF_4: f64 = 12.0 * EPS + 128.0 * EPS * EPS;

/// Default absolute threshold used for singularity/degeneracy detection.
///
/// This is intentionally conservative for geometric predicates and small systems.
///
/// Conceptually, this is an absolute bound for deciding when a scalar should be treated
/// as "numerically zero" (e.g. LU pivots, LDLT diagonal entries).
pub const DEFAULT_SINGULAR_TOL: f64 = 1e-12;

/// Default absolute pivot magnitude threshold used for LU pivot selection / singularity detection.
///
/// This name is kept for backwards compatibility; prefer [`DEFAULT_SINGULAR_TOL`] when the
/// tolerance is not specifically about pivot selection.
pub const DEFAULT_PIVOT_TOL: f64 = DEFAULT_SINGULAR_TOL;

/// Linear algebra errors.
///
/// This enum is `#[non_exhaustive]` — downstream `match` arms must include a
/// wildcard (`_`) pattern to compile, allowing new variants to be added in
/// future minor releases without breaking existing code.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum LaError {
    /// The matrix is (numerically) singular.
    Singular {
        /// The factorization column/step where a suitable pivot/diagonal could not be found.
        pivot_col: usize,
    },
    /// A non-finite value (NaN/∞) was encountered.
    ///
    /// The `(row, col)` coordinate follows a consistent convention across the crate:
    ///
    /// - `row: Some(r), col: c` — a *stored* matrix cell at `(r, c)` is non-finite.
    ///   Used by `Matrix::det`, `Lu::factor`, `Ldlt::factor`, and the `solve_vec`
    ///   paths when they detect a corrupt stored factor (only reachable via
    ///   direct struct construction; `factor` itself rejects such inputs).
    /// - `row: None, col: c` — the non-finite value is either a *vector input*
    ///   entry at index `c`, or a *computed intermediate* at step `c`
    ///   (e.g. an accumulator that overflowed during forward/back substitution).
    NonFinite {
        /// Row of the non-finite entry for a stored matrix cell, or `None` for
        /// a vector-input entry or a computed intermediate. See the variant
        /// docs for the full convention.
        row: Option<usize>,
        /// Column index (stored cell), vector index, or factorization/solve
        /// step where the non-finite value was detected.
        col: usize,
    },
    /// The exact result overflows the target representation (e.g. `f64`).
    ///
    /// Returned by `Matrix::det_exact_f64` and `Matrix::solve_exact_f64`
    /// (requires `exact` feature) when an exact value is too large to
    /// represent as a finite `f64`.
    Overflow {
        /// For vector results (e.g. `solve_exact_f64`), the index of the
        /// component that overflowed.  `None` for scalar results.
        index: Option<usize>,
    },
}

impl LaError {
    /// Construct a [`LaError::NonFinite`] pinpointing a stored matrix cell at `(row, col)`.
    ///
    /// Use this for non-finite values read from a stored `Matrix` entry or
    /// factorization cell.  The resulting error has `row: Some(row), col`,
    /// matching the stored-cell convention documented on
    /// [`NonFinite`](Self::NonFinite).  For vector-input entries or computed
    /// intermediates, use [`non_finite_at`](Self::non_finite_at).
    #[inline]
    #[must_use]
    pub const fn non_finite_cell(row: usize, col: usize) -> Self {
        Self::NonFinite {
            row: Some(row),
            col,
        }
    }

    /// Construct a [`LaError::NonFinite`] pinpointing a vector-input entry or
    /// computed-intermediate step at index `col`.
    ///
    /// Use this for non-finite values in a `Vector` input or an accumulator
    /// that overflowed during forward/back substitution.  The resulting error
    /// has `row: None, col`, matching the vector/intermediate convention
    /// documented on [`NonFinite`](Self::NonFinite).  For stored matrix cells,
    /// use [`non_finite_cell`](Self::non_finite_cell).
    #[inline]
    #[must_use]
    pub const fn non_finite_at(col: usize) -> Self {
        Self::NonFinite { row: None, col }
    }
}

impl fmt::Display for LaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::Singular { pivot_col } => {
                write!(f, "singular matrix at pivot column {pivot_col}")
            }
            Self::NonFinite { row: Some(r), col } => {
                write!(f, "non-finite value at ({r}, {col})")
            }
            Self::NonFinite { row: None, col } => {
                write!(f, "non-finite value at index {col}")
            }
            Self::Overflow { index: Some(i) } => {
                write!(
                    f,
                    "exact result overflows the target representation at index {i}"
                )
            }
            Self::Overflow { index: None } => {
                write!(f, "exact result overflows the target representation")
            }
        }
    }
}

impl std::error::Error for LaError {}

pub use ldlt::Ldlt;
pub use lu::Lu;
pub use matrix::Matrix;
pub use vector::Vector;

/// Common imports for ergonomic usage.
///
/// This prelude re-exports the primary types and constants: [`Matrix`], [`Vector`], [`Lu`],
/// [`Ldlt`], [`LaError`], [`DEFAULT_PIVOT_TOL`], [`DEFAULT_SINGULAR_TOL`], and the determinant
/// error bound coefficients [`ERR_COEFF_2`], [`ERR_COEFF_3`], and [`ERR_COEFF_4`].
///
/// When the `exact` feature is enabled, [`BigInt`] and [`BigRational`]
/// are also re-exported so callers can construct exact values (e.g. as
/// the expected result of `Matrix::det_exact`) without adding
/// `num-bigint` / `num-rational` to their own dependencies.  The most
/// commonly needed `num-traits` items are re-exported alongside them:
/// [`FromPrimitive`] for `BigRational::from_f64` / `from_i64`,
/// [`ToPrimitive`] for `BigRational::to_f64` / `to_i64`, and [`Signed`]
/// for `.is_positive()` / `.is_negative()` / `.abs()`.
pub mod prelude {
    pub use crate::{
        DEFAULT_PIVOT_TOL, DEFAULT_SINGULAR_TOL, ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4, LaError,
        Ldlt, Lu, Matrix, Vector,
    };

    #[cfg(feature = "exact")]
    pub use crate::{BigInt, BigRational, FromPrimitive, Signed, ToPrimitive};
}

#[cfg(test)]
mod tests {
    use super::*;

    use approx::assert_abs_diff_eq;

    #[test]
    fn default_singular_tol_is_expected() {
        assert_abs_diff_eq!(DEFAULT_SINGULAR_TOL, 1e-12, epsilon = 0.0);
        assert_abs_diff_eq!(DEFAULT_PIVOT_TOL, DEFAULT_SINGULAR_TOL, epsilon = 0.0);
    }

    #[test]
    fn laerror_display_formats_singular() {
        let err = LaError::Singular { pivot_col: 3 };
        assert_eq!(err.to_string(), "singular matrix at pivot column 3");
    }

    #[test]
    fn laerror_display_formats_nonfinite_with_row() {
        let err = LaError::NonFinite {
            row: Some(1),
            col: 2,
        };
        assert_eq!(err.to_string(), "non-finite value at (1, 2)");
    }

    #[test]
    fn laerror_display_formats_nonfinite_without_row() {
        let err = LaError::NonFinite { row: None, col: 3 };
        assert_eq!(err.to_string(), "non-finite value at index 3");
    }

    #[test]
    fn laerror_display_formats_overflow() {
        let err = LaError::Overflow { index: None };
        assert_eq!(
            err.to_string(),
            "exact result overflows the target representation"
        );
    }

    #[test]
    fn laerror_display_formats_overflow_with_index() {
        let err = LaError::Overflow { index: Some(2) };
        assert_eq!(
            err.to_string(),
            "exact result overflows the target representation at index 2"
        );
    }

    #[test]
    fn laerror_is_std_error_with_no_source() {
        let err = LaError::Singular { pivot_col: 0 };
        let e: &dyn std::error::Error = &err;
        assert!(e.source().is_none());
    }

    #[test]
    fn prelude_reexports_compile_and_work() {
        use crate::prelude::*;

        // Use the items so we know they are in scope and usable.
        let m = Matrix::<2>::identity();
        let v = Vector::<2>::new([1.0, 2.0]);
        let _ = m.lu(DEFAULT_PIVOT_TOL).unwrap().solve_vec(v).unwrap();
        let _ = m.ldlt(DEFAULT_SINGULAR_TOL).unwrap().solve_vec(v).unwrap();
    }

    /// Exercise every exact-feature re-export via the prelude so a future
    /// refactor that drops one (e.g. removing `Signed` from the prelude
    /// list) fails to compile rather than silently breaking downstream.
    #[cfg(feature = "exact")]
    #[test]
    fn prelude_exact_reexports_compile_and_work() {
        use crate::prelude::*;

        // `BigInt` and `BigRational` constructors.
        let n = BigInt::from(7);
        let r = BigRational::from_integer(n.clone());
        assert_eq!(*r.numer(), n);

        // `FromPrimitive::from_f64` / `from_i64` on `BigRational`.
        let half = BigRational::from_f64(0.5).unwrap();
        let two = BigRational::from_i64(2).unwrap();
        assert_eq!(
            half.clone() + half.clone(),
            BigRational::from_integer(BigInt::from(1))
        );

        // `Signed::is_positive` / `is_negative` / `abs`.
        assert!(half.is_positive());
        assert!(!half.is_negative());
        let neg = -half.clone();
        assert!(neg.is_negative());
        assert_eq!(neg.abs(), half);

        // `ToPrimitive::to_f64` / `to_i64`.
        assert!((half.to_f64().unwrap() - 0.5).abs() <= f64::EPSILON);
        assert_eq!(two.to_i64().unwrap(), 2);
    }
}
