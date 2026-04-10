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
pub use num_rational::BigRational;

mod ldlt;
mod lu;
mod matrix;
mod vector;

use core::fmt;

// ---------------------------------------------------------------------------
// Error-bound constants for determinant error analysis.
//
// These constants bound the absolute error of `det_direct()` relative to the
// *permanent* (sum of absolute products in the Leibniz expansion). The
// constants are conservative over-estimates following Shewchuk's methodology.
//
// These are NOT feature-gated because they use pure f64 arithmetic and are
// useful for adaptive-precision logic even without the `exact` feature.
// ---------------------------------------------------------------------------

const EPS: f64 = f64::EPSILON; // 2^-52

/// Error coefficient for D=2 determinant error bound.
///
/// Accounts for one f64 multiply + one FMA → 2 rounding events.
/// Used in computing the absolute error bound for 2×2 determinants.
pub const ERR_COEFF_2: f64 = 3.0 * EPS + 16.0 * EPS * EPS;

/// Error coefficient for D=3 determinant error bound.
///
/// Accounts for three 2×2 FMA minors + nested FMA combination.
/// Used in computing the absolute error bound for 3×3 determinants.
pub const ERR_COEFF_3: f64 = 8.0 * EPS + 64.0 * EPS * EPS;

/// Error coefficient for D=4 determinant error bound.
///
/// Accounts for six hoisted 2×2 minors → four 3×3 cofactors → FMA row combination.
/// Used in computing the absolute error bound for 4×4 determinants.
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
    NonFinite {
        /// Row of the non-finite entry (for matrix inputs), or `None` when
        /// the error originates from a vector input or a computed intermediate.
        row: Option<usize>,
        /// Column index (for matrix inputs), vector index, or factorization
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
/// When the `exact` feature is enabled, `BigRational` is also
/// re-exported for use with `Matrix::det_exact`.
pub mod prelude {
    pub use crate::{
        DEFAULT_PIVOT_TOL, DEFAULT_SINGULAR_TOL, ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4, LaError,
        Ldlt, Lu, Matrix, Vector,
    };

    #[cfg(feature = "exact")]
    pub use crate::BigRational;
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
}
