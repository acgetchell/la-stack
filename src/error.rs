//! Error types and helpers for linear algebra operations.

use core::fmt;

use crate::Tolerance;

/// Reason an exact result cannot satisfy an exact-to-`f64` conversion contract.
///
/// `RequiresRounding` is recoverable when the caller is willing to opt into a
/// rounded exact-to-`f64` API. `NotFinite` means even the rounded result would
/// not be a finite `f64`.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// let err = LaError::unrepresentable(None, UnrepresentableReason::RequiresRounding);
/// assert!(err.requires_rounding());
///
/// let err = LaError::unrepresentable(None, UnrepresentableReason::NotFinite);
/// assert_eq!(
///     err.unrepresentable_reason(),
///     Some(UnrepresentableReason::NotFinite)
/// );
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum UnrepresentableReason {
    /// A finite `f64` exists only after rounding, but the requested conversion
    /// requires an exact binary64 representation.
    RequiresRounding,
    /// The exact value would convert to NaN or infinity rather than a finite
    /// `f64`.
    NotFinite,
}

/// Linear algebra errors.
///
/// This enum is `#[non_exhaustive]` — downstream `match` arms must include a
/// wildcard (`_`) pattern to compile, allowing new variants to be added in
/// future minor releases without breaking existing code.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// let err = LaError::unrepresentable(None, UnrepresentableReason::RequiresRounding);
/// assert!(err.requires_rounding());
///
/// match LaError::unsupported_dimension(8, MAX_STACK_MATRIX_DISPATCH_DIM) {
///     LaError::UnsupportedDimension { requested, max } => {
///         assert_eq!((requested, max), (8, MAX_STACK_MATRIX_DISPATCH_DIM));
///     }
///     _ => unreachable!("constructor returns the requested variant"),
/// }
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
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
    /// - `row: Some(r), col: c` — the non-finite value is tied to a matrix/factor
    ///   cell at `(r, c)`, either because a stored input/factor cell is already
    ///   non-finite or because factorization computed a non-finite value for
    ///   that cell before storing it.
    /// - `row: None, col: c` — the non-finite value is tied to a vector entry,
    ///   determinant product, solve accumulator, or other scalar/intermediate
    ///   that has no matrix row coordinate.
    NonFinite {
        /// Row of the non-finite entry for a stored matrix cell, or `None` for
        /// a vector-input entry or a computed intermediate. See the variant
        /// docs for the full convention.
        row: Option<usize>,
        /// Column index (stored cell), vector index, or factorization/solve
        /// step where the non-finite value was detected.
        col: usize,
    },
    /// An exact result cannot satisfy the requested finite `f64` conversion.
    ///
    /// Returned by [`Matrix::det_exact_f64`](crate::Matrix::det_exact_f64) and
    /// [`Matrix::solve_exact_f64`](crate::Matrix::solve_exact_f64) (requires the
    /// `exact` feature) when the exact rational value is too large, too small,
    /// or would require rounding in binary64. Also returned by the rounded
    /// exact-to-`f64` APIs when the rounded result would be NaN or infinite.
    Unrepresentable {
        /// For vector results (e.g. `solve_exact_f64`), the index of the
        /// component that failed conversion. `None` for scalar results.
        index: Option<usize>,
        /// Why the requested conversion cannot return a finite `f64`.
        reason: UnrepresentableReason,
    },
    /// Exact determinant scaling overflowed the internal exponent representation.
    DeterminantScaleOverflow {
        /// Matrix dimension `D`.
        dim: usize,
        /// Minimum decomposed f64 exponent among non-zero matrix entries.
        min_exponent: i32,
    },
    /// A requested runtime matrix dimension has no stack-dispatch arm.
    UnsupportedDimension {
        /// Runtime dimension requested by the caller.
        requested: usize,
        /// Largest runtime dimension supported by the dispatch helper.
        max: usize,
    },
    /// A matrix index is outside the `D×D` storage domain.
    IndexOutOfBounds {
        /// Requested row index.
        row: usize,
        /// Requested column index.
        col: usize,
        /// Matrix dimension `D`; valid row and column indices are `< dim`.
        dim: usize,
    },
    /// A tolerance value is not finite and non-negative.
    InvalidTolerance {
        /// Raw tolerance supplied by the caller.
        value: f64,
    },
    /// A matrix required to be symmetric has an asymmetric off-diagonal pair.
    Asymmetric {
        /// Row index of the first asymmetric pair.
        row: usize,
        /// Column index of the first asymmetric pair.
        col: usize,
        /// Matrix dimension `D`.
        dim: usize,
    },
    /// A symmetric matrix failed the positive-semidefinite LDLT domain check.
    NotPositiveSemidefinite {
        /// Factorization column/step where a negative LDLT diagonal was found.
        pivot_col: usize,
        /// Negative diagonal value observed at that step.
        value: f64,
    },
}

impl LaError {
    /// Construct a [`LaError::NonFinite`] pinpointing a stored matrix cell at `(row, col)`.
    ///
    /// Use this for non-finite values read from a stored [`Matrix`](crate::Matrix)
    /// entry or factorization cell, and for non-finite factorization updates
    /// that would be stored at `(row, col)` if accepted.  The resulting error has
    /// `row: Some(row), col`, matching the matrix/factor-cell convention
    /// documented on [`NonFinite`](Self::NonFinite).  For vector-input entries
    /// or scalar intermediates without a matrix row coordinate, use
    /// [`non_finite_at`](Self::non_finite_at).
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(
    ///     LaError::non_finite_cell(1, 2),
    ///     LaError::NonFinite {
    ///         row: Some(1),
    ///         col: 2,
    ///     }
    /// );
    /// ```
    #[inline]
    #[must_use]
    pub const fn non_finite_cell(row: usize, col: usize) -> Self {
        Self::NonFinite {
            row: Some(row),
            col,
        }
    }

    /// Construct a [`LaError::NonFinite`] pinpointing a vector-input entry or
    /// computed scalar/intermediate at index `col`.
    ///
    /// Use this for non-finite values in a [`Vector`](crate::Vector) input,
    /// determinant scalar, tolerance-scale accumulator, or solve accumulator
    /// that overflowed during forward/back substitution.  The resulting error
    /// has `row: None, col`, matching the vector/scalar-intermediate convention
    /// documented on [`NonFinite`](Self::NonFinite).  For stored matrix cells or
    /// computed factorization updates tied to a matrix cell, use
    /// [`non_finite_cell`](Self::non_finite_cell).
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(
    ///     LaError::non_finite_at(2),
    ///     LaError::NonFinite { row: None, col: 2 }
    /// );
    /// ```
    #[inline]
    #[must_use]
    pub const fn non_finite_at(col: usize) -> Self {
        Self::NonFinite { row: None, col }
    }

    /// Construct a [`LaError::Unrepresentable`] for exact-to-`f64` conversion.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(
    ///     LaError::unrepresentable(Some(2), UnrepresentableReason::RequiresRounding),
    ///     LaError::Unrepresentable {
    ///         index: Some(2),
    ///         reason: UnrepresentableReason::RequiresRounding,
    ///     }
    /// );
    /// ```
    #[inline]
    #[must_use]
    pub const fn unrepresentable(index: Option<usize>, reason: UnrepresentableReason) -> Self {
        Self::Unrepresentable { index, reason }
    }

    /// Return the reason for an exact-to-`f64` conversion failure.
    ///
    /// This is a concise alternative to matching the full
    /// [`LaError::Unrepresentable`] variant when callers only need the
    /// conversion reason.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let err = LaError::unrepresentable(None, UnrepresentableReason::RequiresRounding);
    /// assert_eq!(
    ///     err.unrepresentable_reason(),
    ///     Some(UnrepresentableReason::RequiresRounding)
    /// );
    /// assert_eq!(LaError::Singular { pivot_col: 0 }.unrepresentable_reason(), None);
    /// ```
    #[inline]
    #[must_use]
    pub const fn unrepresentable_reason(&self) -> Option<UnrepresentableReason> {
        match self {
            Self::Unrepresentable { reason, .. } => Some(*reason),
            _ => None,
        }
    }

    /// Return `true` when strict exact-to-`f64` conversion only failed because
    /// rounding would be required.
    ///
    /// This is useful at the call site that wants to retry with an explicit
    /// rounded exact-to-`f64` API while still propagating non-finite conversion
    /// failures.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let err = LaError::unrepresentable(None, UnrepresentableReason::RequiresRounding);
    /// assert!(err.requires_rounding());
    ///
    /// let err = LaError::unrepresentable(None, UnrepresentableReason::NotFinite);
    /// assert!(!err.requires_rounding());
    /// ```
    #[inline]
    #[must_use]
    pub const fn requires_rounding(&self) -> bool {
        matches!(
            self,
            Self::Unrepresentable {
                reason: UnrepresentableReason::RequiresRounding,
                ..
            }
        )
    }

    /// Construct a [`LaError::DeterminantScaleOverflow`] for exact determinant scaling.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(
    ///     LaError::determinant_scale_overflow(3, -1074),
    ///     LaError::DeterminantScaleOverflow {
    ///         dim: 3,
    ///         min_exponent: -1074,
    ///     }
    /// );
    /// ```
    #[inline]
    #[must_use]
    pub const fn determinant_scale_overflow(dim: usize, min_exponent: i32) -> Self {
        Self::DeterminantScaleOverflow { dim, min_exponent }
    }

    /// Construct a [`LaError::UnsupportedDimension`] for runtime stack dispatch.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(
    ///     LaError::unsupported_dimension(8, MAX_STACK_MATRIX_DISPATCH_DIM),
    ///     LaError::UnsupportedDimension {
    ///         requested: 8,
    ///         max: MAX_STACK_MATRIX_DISPATCH_DIM,
    ///     }
    /// );
    /// ```
    #[inline]
    #[must_use]
    pub const fn unsupported_dimension(requested: usize, max: usize) -> Self {
        Self::UnsupportedDimension { requested, max }
    }

    /// Construct a [`LaError::IndexOutOfBounds`] for a `D×D` matrix index.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(
    ///     LaError::index_out_of_bounds(2, 0, 2),
    ///     LaError::IndexOutOfBounds {
    ///         row: 2,
    ///         col: 0,
    ///         dim: 2,
    ///     }
    /// );
    /// ```
    #[inline]
    #[must_use]
    pub const fn index_out_of_bounds(row: usize, col: usize, dim: usize) -> Self {
        Self::IndexOutOfBounds { row, col, dim }
    }

    /// Construct a [`LaError::InvalidTolerance`] for a raw tolerance value.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(
    ///     LaError::invalid_tolerance(-1.0),
    ///     LaError::InvalidTolerance { value: -1.0 }
    /// );
    /// ```
    #[inline]
    #[must_use]
    pub const fn invalid_tolerance(value: f64) -> Self {
        Self::InvalidTolerance { value }
    }

    /// Construct a [`LaError::Asymmetric`] for a `D×D` matrix.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(
    ///     LaError::asymmetric(0, 1, 3),
    ///     LaError::Asymmetric {
    ///         row: 0,
    ///         col: 1,
    ///         dim: 3,
    ///     }
    /// );
    /// ```
    #[inline]
    #[must_use]
    pub const fn asymmetric(row: usize, col: usize, dim: usize) -> Self {
        Self::Asymmetric { row, col, dim }
    }

    /// Construct a [`LaError::NotPositiveSemidefinite`] for LDLT factorization.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(
    ///     LaError::not_positive_semidefinite(1, -3.0),
    ///     LaError::NotPositiveSemidefinite {
    ///         pivot_col: 1,
    ///         value: -3.0,
    ///     }
    /// );
    /// ```
    #[inline]
    #[must_use]
    pub const fn not_positive_semidefinite(pivot_col: usize, value: f64) -> Self {
        Self::NotPositiveSemidefinite { pivot_col, value }
    }

    /// Parse a raw tolerance into a finite, non-negative [`Tolerance`].
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// assert_eq!(LaError::validate_tolerance(1e-12)?.get(), 1e-12);
    ///
    /// let raw = 0.0;
    /// let tol = LaError::validate_tolerance(raw)?;
    /// let _lu = Matrix::<2>::identity().lu(tol)?;
    ///
    /// assert_eq!(
    ///     LaError::validate_tolerance(-1.0),
    ///     Err(LaError::InvalidTolerance { value: -1.0 })
    /// );
    /// # Ok::<(), LaError>(())
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::InvalidTolerance`] when `value` is NaN, infinite, or
    /// negative.
    #[inline]
    pub const fn validate_tolerance(value: f64) -> Result<Tolerance, Self> {
        Tolerance::new(value)
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
            Self::Unrepresentable {
                index: Some(i),
                reason: UnrepresentableReason::RequiresRounding,
            } => write!(
                f,
                "exact result requires rounding to fit finite f64 at index {i}"
            ),
            Self::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::RequiresRounding,
            } => write!(f, "exact result requires rounding to fit finite f64"),
            Self::Unrepresentable {
                index: Some(i),
                reason: UnrepresentableReason::NotFinite,
            } => write!(f, "exact result does not round to finite f64 at index {i}"),
            Self::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::NotFinite,
            } => write!(f, "exact result does not round to finite f64"),
            Self::DeterminantScaleOverflow { dim, min_exponent } => {
                write!(
                    f,
                    "exact determinant scale exponent overflows for dimension {dim} with minimum entry exponent {min_exponent}"
                )
            }
            Self::UnsupportedDimension { requested, max } => {
                write!(
                    f,
                    "unsupported matrix dimension {requested}; maximum stack-dispatch dimension is {max}"
                )
            }
            Self::IndexOutOfBounds { row, col, dim } => {
                write!(
                    f,
                    "matrix index ({row}, {col}) is out of bounds for dimension {dim}"
                )
            }
            Self::InvalidTolerance { value } => {
                write!(f, "invalid tolerance {value}; expected finite value >= 0")
            }
            Self::Asymmetric { row, col, dim } => {
                write!(
                    f,
                    "matrix is not symmetric for dimension {dim}: asymmetric pair ({row}, {col})"
                )
            }
            Self::NotPositiveSemidefinite { pivot_col, value } => {
                write!(
                    f,
                    "matrix is not positive semidefinite at LDLT pivot column {pivot_col}: diagonal value {value} < 0"
                )
            }
        }
    }
}

impl std::error::Error for LaError {}

#[cfg(test)]
mod tests {
    use super::*;

    use core::assert_matches;

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
    fn laerror_display_formats_unrepresentable_requires_rounding() {
        let err = LaError::Unrepresentable {
            index: None,
            reason: UnrepresentableReason::RequiresRounding,
        };
        assert_eq!(
            err.to_string(),
            "exact result requires rounding to fit finite f64"
        );
    }

    #[test]
    fn laerror_display_formats_unrepresentable_requires_rounding_with_index() {
        let err = LaError::Unrepresentable {
            index: Some(2),
            reason: UnrepresentableReason::RequiresRounding,
        };
        assert_eq!(
            err.to_string(),
            "exact result requires rounding to fit finite f64 at index 2"
        );
    }

    #[test]
    fn laerror_display_formats_unrepresentable_not_finite() {
        let err = LaError::Unrepresentable {
            index: None,
            reason: UnrepresentableReason::NotFinite,
        };
        assert_eq!(err.to_string(), "exact result does not round to finite f64");
    }

    #[test]
    fn laerror_display_formats_unrepresentable_not_finite_with_index() {
        let err = LaError::Unrepresentable {
            index: Some(2),
            reason: UnrepresentableReason::NotFinite,
        };
        assert_eq!(
            err.to_string(),
            "exact result does not round to finite f64 at index 2"
        );
    }

    #[test]
    fn laerror_unrepresentable_reason_reports_typed_reason() {
        let rounding = LaError::Unrepresentable {
            index: Some(2),
            reason: UnrepresentableReason::RequiresRounding,
        };
        let not_finite = LaError::Unrepresentable {
            index: None,
            reason: UnrepresentableReason::NotFinite,
        };

        assert_eq!(
            rounding.unrepresentable_reason(),
            Some(UnrepresentableReason::RequiresRounding)
        );
        assert_eq!(
            not_finite.unrepresentable_reason(),
            Some(UnrepresentableReason::NotFinite)
        );
        assert_eq!(
            LaError::Singular { pivot_col: 0 }.unrepresentable_reason(),
            None
        );
    }

    #[test]
    fn laerror_requires_rounding_only_matches_rounding_reason() {
        assert!(
            LaError::Unrepresentable {
                index: Some(2),
                reason: UnrepresentableReason::RequiresRounding,
            }
            .requires_rounding()
        );
        assert!(
            !LaError::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::NotFinite,
            }
            .requires_rounding()
        );
        assert!(!LaError::Singular { pivot_col: 0 }.requires_rounding());
    }

    #[test]
    fn laerror_display_formats_determinant_scale_overflow() {
        let err = LaError::DeterminantScaleOverflow {
            dim: 3,
            min_exponent: -1074,
        };
        assert_eq!(
            err.to_string(),
            "exact determinant scale exponent overflows for dimension 3 with minimum entry exponent -1074"
        );
    }

    #[test]
    fn laerror_display_formats_unsupported_dimension() {
        let err = LaError::UnsupportedDimension {
            requested: 8,
            max: crate::MAX_STACK_MATRIX_DISPATCH_DIM,
        };
        assert_eq!(
            err.to_string(),
            "unsupported matrix dimension 8; maximum stack-dispatch dimension is 7"
        );
    }

    #[test]
    fn laerror_display_formats_index_out_of_bounds() {
        let err = LaError::IndexOutOfBounds {
            row: 3,
            col: 0,
            dim: 3,
        };
        assert_eq!(
            err.to_string(),
            "matrix index (3, 0) is out of bounds for dimension 3"
        );
    }

    #[test]
    fn laerror_display_formats_invalid_tolerance() {
        let err = LaError::InvalidTolerance { value: -1.0 };
        assert_eq!(
            err.to_string(),
            "invalid tolerance -1; expected finite value >= 0"
        );
    }

    #[test]
    fn validate_tolerance_matches_tolerance_new() {
        for value in [0.0, 1e-12, f64::MAX] {
            assert_eq!(LaError::validate_tolerance(value), Tolerance::new(value));
        }

        assert_eq!(
            LaError::validate_tolerance(-1.0),
            Err(LaError::InvalidTolerance { value: -1.0 })
        );
        assert_matches!(
            LaError::validate_tolerance(f64::NAN),
            Err(LaError::InvalidTolerance { value }) if value.is_nan()
        );
        assert_eq!(
            LaError::validate_tolerance(f64::INFINITY),
            Err(LaError::InvalidTolerance {
                value: f64::INFINITY,
            })
        );
        assert_eq!(
            LaError::validate_tolerance(f64::NEG_INFINITY),
            Err(LaError::InvalidTolerance {
                value: f64::NEG_INFINITY,
            })
        );
    }

    #[test]
    fn laerror_display_formats_asymmetric() {
        let err = LaError::Asymmetric {
            row: 0,
            col: 2,
            dim: 3,
        };
        assert_eq!(
            err.to_string(),
            "matrix is not symmetric for dimension 3: asymmetric pair (0, 2)"
        );
    }

    #[test]
    fn laerror_display_formats_not_positive_semidefinite() {
        let err = LaError::NotPositiveSemidefinite {
            pivot_col: 1,
            value: -3.0,
        };
        assert_eq!(
            err.to_string(),
            "matrix is not positive semidefinite at LDLT pivot column 1: diagonal value -3 < 0"
        );
    }

    #[test]
    fn laerror_is_std_error_with_no_source() {
        let err = LaError::Singular { pivot_col: 0 };
        let e: &dyn std::error::Error = &err;
        assert!(e.source().is_none());
    }
}
