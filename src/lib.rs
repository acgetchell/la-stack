#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Small, stack-allocated linear algebra for fixed dimensions.
//!
//! This crate intentionally focuses on a minimal, explicit API surface:
//! - const-generic sizes (no dynamic dimensions)
//! - stack-backed `Copy` types
//! - explicit algorithms (LU, solve, det)

mod lu;
mod matrix;
mod vector;

use core::fmt;

/// Default absolute pivot tolerance used for singularity detection.
///
/// This is intentionally conservative for geometric predicates and small systems.
pub const DEFAULT_PIVOT_TOL: f64 = 1e-12;

/// Linear algebra errors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LaError {
    /// The matrix is (numerically) singular.
    Singular {
        /// The column where a suitable pivot could not be found.
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

pub use lu::Lu;
pub use matrix::Matrix;
pub use vector::Vector;
