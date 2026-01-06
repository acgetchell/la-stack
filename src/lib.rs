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

mod ldlt;
mod lu;
mod matrix;
mod vector;

use core::fmt;

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
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LaError {
    /// The matrix is (numerically) singular.
    Singular {
        /// The factorization column/step where a suitable pivot/diagonal could not be found.
        pivot_col: usize,
    },
    /// A non-finite value (NaN/∞) was encountered.
    NonFinite {
        /// The column being processed when non-finite values were detected.
        pivot_col: usize,
    },
}

impl fmt::Display for LaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::Singular { pivot_col } => {
                write!(f, "singular matrix at pivot column {pivot_col}")
            }
            Self::NonFinite { pivot_col } => {
                write!(
                    f,
                    "non-finite value encountered at pivot column {pivot_col}"
                )
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
/// [`Ldlt`], [`LaError`], [`DEFAULT_PIVOT_TOL`], and [`DEFAULT_SINGULAR_TOL`].
pub mod prelude {
    pub use crate::{DEFAULT_PIVOT_TOL, DEFAULT_SINGULAR_TOL, LaError, Ldlt, Lu, Matrix, Vector};
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
    fn laerror_display_formats_nonfinite() {
        let err = LaError::NonFinite { pivot_col: 2 };
        assert_eq!(
            err.to_string(),
            "non-finite value encountered at pivot column 2"
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
