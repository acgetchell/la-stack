#![forbid(unsafe_code)]

//! Typed error categories and helpers for linear algebra operations.

use core::fmt;

/// Floating-point operation that produced a non-finite intermediate or result.
///
/// # Examples
/// ```
/// use la_stack::ArithmeticOperation;
///
/// assert_eq!(ArithmeticOperation::LuSolve.to_string(), "LU solve");
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum ArithmeticOperation {
    /// Matrix infinity-norm calculation.
    MatrixInfinityNorm,
    /// Matrix symmetry validation.
    SymmetryCheck,
    /// LU factorization.
    LuFactorization,
    /// LDLT factorization.
    LdltFactorization,
    /// Forward or backward substitution with an LU factorization.
    LuSolve,
    /// Forward, diagonal, or backward substitution with an LDLT factorization.
    LdltSolve,
    /// Determinant calculation.
    Determinant,
    /// Determinant error-bound calculation.
    DeterminantErrorBound,
    /// Vector dot-product calculation.
    VectorDotProduct,
    /// Vector squared-norm calculation.
    VectorSquaredNorm,
}

impl fmt::Display for ArithmeticOperation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::MatrixInfinityNorm => "matrix infinity norm",
            Self::SymmetryCheck => "symmetry check",
            Self::LuFactorization => "LU factorization",
            Self::LdltFactorization => "LDLT factorization",
            Self::LuSolve => "LU solve",
            Self::LdltSolve => "LDLT solve",
            Self::Determinant => "determinant",
            Self::DeterminantErrorBound => "determinant error bound",
            Self::VectorDotProduct => "vector dot product",
            Self::VectorSquaredNorm => "vector squared norm",
        })
    }
}

/// Factorization whose pivot policy rejected a matrix as numerically singular.
///
/// # Examples
/// ```
/// use la_stack::FactorizationKind;
///
/// assert_eq!(FactorizationKind::Ldlt.to_string(), "LDLT");
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum FactorizationKind {
    /// LU factorization with partial pivoting.
    Lu,
    /// LDLT factorization without pivoting.
    Ldlt,
}

impl fmt::Display for FactorizationKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Lu => "LU",
            Self::Ldlt => "LDLT",
        })
    }
}

/// Reason a raw tolerance cannot become a [`crate::Tolerance`].
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// match LaError::invalid_tolerance(-1.0) {
///     LaError::InvalidTolerance {
///         reason: InvalidToleranceReason::Negative,
///         ..
///     } => {}
///     _ => unreachable!("a finite negative tolerance has a negative reason"),
/// }
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum InvalidToleranceReason {
    /// The tolerance is finite but negative.
    Negative,
    /// The tolerance is NaN or positive/negative infinity.
    NotFinite,
}

/// Location at which a non-finite value was observed.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// match LaError::non_finite_input_matrix(1, 2) {
///     LaError::NonFinite {
///         location: NonFiniteLocation::MatrixCell { row, col, .. },
///         ..
///     } => assert_eq!((row, col), (1, 2)),
///     _ => unreachable!("constructor returns a matrix-cell location"),
/// }
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum NonFiniteLocation {
    /// Cell `(row, col)` in matrix-shaped storage or computation.
    #[non_exhaustive]
    MatrixCell {
        /// Matrix row.
        row: usize,
        /// Matrix column.
        col: usize,
    },
    /// Entry in a vector input.
    #[non_exhaustive]
    VectorEntry {
        /// Vector index.
        index: usize,
    },
    /// Indexed step in a factorization, solve, or reduction.
    #[non_exhaustive]
    Step {
        /// Step index.
        index: usize,
    },
    /// Scalar value without a meaningful matrix or vector coordinate.
    Scalar,
}

/// Provenance of a non-finite value.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// let err = LaError::non_finite_computation_scalar(ArithmeticOperation::Determinant);
/// match err {
///     LaError::NonFinite {
///         origin: NonFiniteOrigin::Computation { operation, .. },
///         ..
///     } => assert_eq!(operation, ArithmeticOperation::Determinant),
///     _ => unreachable!("constructor preserves computation provenance"),
/// }
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum NonFiniteOrigin {
    /// The caller supplied a non-finite input.
    Input,
    /// Finite inputs produced a non-finite arithmetic result.
    #[non_exhaustive]
    Computation {
        /// Operation that produced the value.
        operation: ArithmeticOperation,
    },
}

/// Reason a symmetric matrix is outside the positive-semidefinite LDLT domain.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// match LaError::not_positive_semidefinite_negative(1, -3.0) {
///     LaError::NotPositiveSemidefinite {
///         violation: PositiveSemidefiniteViolation::NegativePivot { value, .. },
///         ..
///     } => assert_eq!(value, -3.0),
///     _ => unreachable!("constructor preserves the PSD violation"),
/// }
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum PositiveSemidefiniteViolation {
    /// LDLT produced a strictly negative diagonal pivot.
    #[non_exhaustive]
    NegativePivot {
        /// Observed negative pivot value.
        value: f64,
    },
    /// A zero diagonal pivot still has a non-zero coupling below it.
    #[non_exhaustive]
    ZeroPivotCoupling {
        /// Row containing the non-zero coupling.
        row: usize,
        /// Observed coupling value.
        value: f64,
    },
}

/// Mathematical or numerical reason a matrix was classified as singular.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// match LaError::singular_numerical(1, FactorizationKind::Lu, 0.0, 1e-12) {
///     LaError::Singular {
///         reason: SingularityReason::Numerical { factorization, .. }, ..
///     } => assert_eq!(factorization, FactorizationKind::Lu),
///     _ => unreachable!("constructor preserves the singularity reason"),
/// }
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum SingularityReason {
    /// The algorithm proved that the pivot is exactly zero.
    Exact,
    /// A finite pivot was rejected by the factorization tolerance.
    #[non_exhaustive]
    Numerical {
        /// Factorization applying the tolerance policy.
        factorization: FactorizationKind,
        /// Absolute magnitude of the rejected pivot.
        pivot_magnitude: f64,
        /// Tolerance against which the magnitude was compared.
        tolerance: f64,
    },
}

/// Reason an exact result cannot satisfy an exact-to-`f64` conversion contract.
///
/// `RequiresRounding` is recoverable when the caller is willing to opt into a
/// rounded exact-to-`f64` API. `NotFinite` means no finite `f64` can represent
/// the result even after rounding.
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
    /// No finite `f64` can represent the exact value after rounding.
    NotFinite,
}

/// Linear algebra errors.
///
/// This enum and each struct-style variant are `#[non_exhaustive]` so downstream
/// matches must retain a wildcard for future error categories and fields.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// match LaError::singular_exact(2) {
///     LaError::Singular {
///         pivot_col,
///         reason: SingularityReason::Exact,
///         ..
///     } => assert_eq!(pivot_col, 2),
///     _ => unreachable!("constructor returns an exact singularity"),
/// }
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum LaError {
    /// A matrix is exactly or numerically singular.
    #[non_exhaustive]
    Singular {
        /// Factorization column or step where a usable pivot was unavailable.
        pivot_col: usize,
        /// Typed reason for the singularity classification.
        reason: SingularityReason,
    },
    /// A caller input or arithmetic result is NaN or infinite.
    #[non_exhaustive]
    NonFinite {
        /// Typed location of the value.
        location: NonFiniteLocation,
        /// Whether the value came from input or a particular computation.
        origin: NonFiniteOrigin,
    },
    /// An exact result cannot satisfy the requested finite-`f64` conversion.
    #[non_exhaustive]
    Unrepresentable {
        /// Failed vector component, or `None` for a scalar result.
        index: Option<usize>,
        /// Reason the conversion contract cannot be satisfied.
        reason: UnrepresentableReason,
    },
    /// Exact determinant scaling overflowed the internal exponent representation.
    #[non_exhaustive]
    DeterminantScaleOverflow {
        /// Matrix dimension `D`.
        dim: usize,
        /// Minimum decomposed binary64 exponent among non-zero entries.
        min_exponent: i32,
    },
    /// A runtime matrix dimension has no stack-dispatch arm.
    #[non_exhaustive]
    UnsupportedDimension {
        /// Runtime dimension requested by the caller.
        requested: usize,
        /// Largest dimension supported by the dispatch helper.
        max: usize,
    },
    /// A matrix index is outside the `D×D` storage domain.
    #[non_exhaustive]
    IndexOutOfBounds {
        /// Requested row.
        row: usize,
        /// Requested column.
        col: usize,
        /// Matrix dimension `D`; valid indices are less than this value.
        dim: usize,
    },
    /// A raw tolerance is negative or non-finite.
    #[non_exhaustive]
    InvalidTolerance {
        /// Raw value supplied by the caller.
        value: f64,
        /// Typed reason the value violates the tolerance invariant.
        reason: InvalidToleranceReason,
    },
    /// A matrix required to be symmetric has an asymmetric off-diagonal pair.
    #[non_exhaustive]
    Asymmetric {
        /// Row of the upper-triangular entry.
        row: usize,
        /// Column of the upper-triangular entry.
        col: usize,
        /// Matrix dimension `D`.
        dim: usize,
        /// Observed entry at `(row, col)`.
        upper: f64,
        /// Observed entry at `(col, row)`.
        lower: f64,
        /// Maximum absolute difference allowed by the symmetry check.
        allowed_abs_diff: f64,
    },
    /// A symmetric matrix is outside the positive-semidefinite LDLT domain.
    #[non_exhaustive]
    NotPositiveSemidefinite {
        /// LDLT pivot column or step where the violation was detected.
        pivot_col: usize,
        /// Typed PSD-domain violation.
        violation: PositiveSemidefiniteViolation,
    },
}

impl LaError {
    /// Construct a [`LaError::Singular`] error proving that the pivot at
    /// `pivot_col` is exactly zero.
    #[inline]
    #[must_use]
    pub const fn singular_exact(pivot_col: usize) -> Self {
        Self::Singular {
            pivot_col,
            reason: SingularityReason::Exact,
        }
    }

    /// Construct a [`LaError::Singular`] error for a pivot rejected by a
    /// factorization tolerance, preserving the factorization, pivot magnitude,
    /// and tolerance in [`SingularityReason::Numerical`].
    #[inline]
    #[must_use]
    pub const fn singular_numerical(
        pivot_col: usize,
        factorization: FactorizationKind,
        pivot_magnitude: f64,
        tolerance: f64,
    ) -> Self {
        Self::Singular {
            pivot_col,
            reason: SingularityReason::Numerical {
                factorization,
                pivot_magnitude,
                tolerance,
            },
        }
    }

    /// Construct a [`LaError::NonFinite`] input error located at matrix cell
    /// `(row, col)`.
    #[inline]
    #[must_use]
    pub const fn non_finite_input_matrix(row: usize, col: usize) -> Self {
        Self::NonFinite {
            location: NonFiniteLocation::MatrixCell { row, col },
            origin: NonFiniteOrigin::Input,
        }
    }

    /// Construct a [`LaError::NonFinite`] input error located at vector entry
    /// `index`.
    #[inline]
    #[must_use]
    pub const fn non_finite_input_vector(index: usize) -> Self {
        Self::NonFinite {
            location: NonFiniteLocation::VectorEntry { index },
            origin: NonFiniteOrigin::Input,
        }
    }

    /// Construct a [`LaError::NonFinite`] input error for a scalar without an
    /// index or matrix coordinate.
    #[inline]
    #[must_use]
    pub const fn non_finite_input_scalar() -> Self {
        Self::NonFinite {
            location: NonFiniteLocation::Scalar,
            origin: NonFiniteOrigin::Input,
        }
    }

    /// Construct a [`LaError::NonFinite`] computation error at matrix cell
    /// `(row, col)`, retaining the originating `operation`.
    #[inline]
    #[must_use]
    pub const fn non_finite_computation_matrix(
        operation: ArithmeticOperation,
        row: usize,
        col: usize,
    ) -> Self {
        Self::NonFinite {
            location: NonFiniteLocation::MatrixCell { row, col },
            origin: NonFiniteOrigin::Computation { operation },
        }
    }

    /// Construct a [`LaError::NonFinite`] computation error at `index`,
    /// retaining the originating `operation`.
    #[inline]
    #[must_use]
    pub const fn non_finite_computation_step(operation: ArithmeticOperation, index: usize) -> Self {
        Self::NonFinite {
            location: NonFiniteLocation::Step { index },
            origin: NonFiniteOrigin::Computation { operation },
        }
    }

    /// Construct a scalar [`LaError::NonFinite`] computation error retaining
    /// the originating `operation`.
    #[inline]
    #[must_use]
    pub const fn non_finite_computation_scalar(operation: ArithmeticOperation) -> Self {
        Self::NonFinite {
            location: NonFiniteLocation::Scalar,
            origin: NonFiniteOrigin::Computation { operation },
        }
    }

    /// Construct a [`LaError::Unrepresentable`] conversion failure for a scalar
    /// (`index = None`) or vector component (`index = Some(_)`).
    #[inline]
    #[must_use]
    pub const fn unrepresentable(index: Option<usize>, reason: UnrepresentableReason) -> Self {
        Self::Unrepresentable { index, reason }
    }

    /// Return the typed exact-to-`f64` conversion reason, or `None` for every
    /// other error variant.
    #[inline]
    #[must_use]
    pub const fn unrepresentable_reason(&self) -> Option<UnrepresentableReason> {
        match self {
            Self::Unrepresentable { reason, .. } => Some(*reason),
            _ => None,
        }
    }

    /// Return whether this is a `RequiresRounding` conversion failure for which
    /// retrying with an explicit rounded API may succeed.
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

    /// Construct a [`LaError::DeterminantScaleOverflow`] retaining the matrix
    /// dimension and minimum decomposed entry exponent.
    #[inline]
    #[must_use]
    pub const fn determinant_scale_overflow(dim: usize, min_exponent: i32) -> Self {
        Self::DeterminantScaleOverflow { dim, min_exponent }
    }

    /// Construct a [`LaError::UnsupportedDimension`] retaining the requested
    /// and maximum supported dimensions.
    #[inline]
    #[must_use]
    pub const fn unsupported_dimension(requested: usize, max: usize) -> Self {
        Self::UnsupportedDimension { requested, max }
    }

    /// Construct a [`LaError::IndexOutOfBounds`] retaining the requested matrix
    /// coordinates and dimension.
    #[inline]
    #[must_use]
    pub const fn index_out_of_bounds(row: usize, col: usize, dim: usize) -> Self {
        Self::IndexOutOfBounds { row, col, dim }
    }

    /// Construct an invalid-tolerance error and classify its typed reason.
    ///
    /// This low-level constructor assumes `value` has already failed the
    /// tolerance invariant; raw caller input should normally be parsed through
    /// [`crate::Tolerance::try_new`].
    /// Non-finiteness takes precedence over negativity, so negative infinity is
    /// classified as [`InvalidToleranceReason::NotFinite`].
    #[inline]
    #[must_use]
    pub const fn invalid_tolerance(value: f64) -> Self {
        let reason = if value.is_finite() {
            InvalidToleranceReason::Negative
        } else {
            InvalidToleranceReason::NotFinite
        };
        Self::InvalidTolerance { value, reason }
    }

    /// Construct a [`LaError::Asymmetric`] error for the pair `(row, col)` and
    /// `(col, row)`, retaining both observed values and the effective absolute
    /// difference bound.
    #[inline]
    #[must_use]
    pub const fn asymmetric(
        row: usize,
        col: usize,
        dim: usize,
        upper: f64,
        lower: f64,
        allowed_abs_diff: f64,
    ) -> Self {
        Self::Asymmetric {
            row,
            col,
            dim,
            upper,
            lower,
            allowed_abs_diff,
        }
    }

    /// Construct a [`LaError::NotPositiveSemidefinite`] error for a negative
    /// LDLT diagonal pivot.
    #[inline]
    #[must_use]
    pub const fn not_positive_semidefinite_negative(pivot_col: usize, value: f64) -> Self {
        Self::NotPositiveSemidefinite {
            pivot_col,
            violation: PositiveSemidefiniteViolation::NegativePivot { value },
        }
    }

    /// Construct a [`LaError::NotPositiveSemidefinite`] error for a zero LDLT
    /// diagonal with a non-zero coupling at `row`, distinguishing the observed
    /// positive-semidefinite violation from an uncoupled singular pivot.
    #[inline]
    #[must_use]
    pub const fn not_positive_semidefinite_zero_coupling(
        pivot_col: usize,
        row: usize,
        value: f64,
    ) -> Self {
        Self::NotPositiveSemidefinite {
            pivot_col,
            violation: PositiveSemidefiniteViolation::ZeroPivotCoupling { row, value },
        }
    }
}

/// Write the structured location portion of [`LaError::NonFinite`]'s public
/// display contract without allocating an intermediate string.
///
fn write_non_finite_location(
    f: &mut fmt::Formatter<'_>,
    location: NonFiniteLocation,
) -> fmt::Result {
    match location {
        NonFiniteLocation::MatrixCell { row, col } => {
            write!(f, "matrix cell ({row}, {col})")
        }
        NonFiniteLocation::VectorEntry { index } => write!(f, "vector entry {index}"),
        NonFiniteLocation::Step { index } => write!(f, "step {index}"),
        NonFiniteLocation::Scalar => f.write_str("scalar value"),
    }
}

/// Write a [`LaError::NonFinite`] message from its structured location and
/// origin while distinguishing scalar input from a computed scalar result.
fn write_non_finite(
    f: &mut fmt::Formatter<'_>,
    location: NonFiniteLocation,
    origin: NonFiniteOrigin,
) -> fmt::Result {
    match (location, origin) {
        (NonFiniteLocation::Scalar, NonFiniteOrigin::Input) => {
            f.write_str("non-finite scalar input")
        }
        (NonFiniteLocation::Scalar, NonFiniteOrigin::Computation { operation }) => {
            write!(f, "non-finite scalar result computed during {operation}")
        }
        (location, NonFiniteOrigin::Input) => {
            f.write_str("non-finite input value at ")?;
            write_non_finite_location(f, location)
        }
        (location, NonFiniteOrigin::Computation { operation }) => {
            write!(f, "non-finite value computed during {operation} at ")?;
            write_non_finite_location(f, location)
        }
    }
}

impl fmt::Display for LaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::Singular {
                pivot_col,
                reason: SingularityReason::Exact,
            } => write!(f, "matrix is exactly singular at pivot column {pivot_col}"),
            Self::Singular {
                pivot_col,
                reason:
                    SingularityReason::Numerical {
                        factorization,
                        pivot_magnitude,
                        tolerance,
                    },
            } => write!(
                f,
                "matrix is numerically singular during {factorization} factorization at pivot column {pivot_col}: pivot magnitude {pivot_magnitude} <= tolerance {tolerance}"
            ),
            Self::NonFinite { location, origin } => write_non_finite(f, location, origin),
            Self::Unrepresentable {
                index: Some(index),
                reason: UnrepresentableReason::RequiresRounding,
            } => write!(
                f,
                "exact result requires rounding to fit finite f64 at index {index}"
            ),
            Self::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::RequiresRounding,
            } => f.write_str("exact result requires rounding to fit finite f64"),
            Self::Unrepresentable {
                index: Some(index),
                reason: UnrepresentableReason::NotFinite,
            } => write!(
                f,
                "exact result has no finite f64 representation after rounding at index {index}"
            ),
            Self::Unrepresentable {
                index: None,
                reason: UnrepresentableReason::NotFinite,
            } => f.write_str("exact result has no finite f64 representation after rounding"),
            Self::DeterminantScaleOverflow { dim, min_exponent } => write!(
                f,
                "exact determinant scale exponent overflows for dimension {dim} with minimum entry exponent {min_exponent}"
            ),
            Self::UnsupportedDimension { requested, max } => write!(
                f,
                "unsupported matrix dimension {requested}; maximum stack-dispatch dimension is {max}"
            ),
            Self::IndexOutOfBounds { row, col, dim } => write!(
                f,
                "matrix index ({row}, {col}) is out of bounds for dimension {dim}"
            ),
            Self::InvalidTolerance {
                value,
                reason: InvalidToleranceReason::Negative,
            } => write!(f, "invalid tolerance {value}; expected value >= 0"),
            Self::InvalidTolerance {
                value,
                reason: InvalidToleranceReason::NotFinite,
            } => write!(f, "invalid tolerance {value}; expected a finite value"),
            Self::Asymmetric {
                row,
                col,
                dim,
                upper,
                lower,
                allowed_abs_diff,
            } => write!(
                f,
                "matrix is not symmetric for dimension {dim}: entry ({row}, {col}) = {upper} and entry ({col}, {row}) = {lower} differ by more than allowed absolute difference {allowed_abs_diff}"
            ),
            Self::NotPositiveSemidefinite {
                pivot_col,
                violation: PositiveSemidefiniteViolation::NegativePivot { value },
            } => write!(
                f,
                "matrix is not positive semidefinite at LDLT pivot column {pivot_col}: diagonal value {value} < 0"
            ),
            Self::NotPositiveSemidefinite {
                pivot_col,
                violation: PositiveSemidefiniteViolation::ZeroPivotCoupling { row, value },
            } => write!(
                f,
                "matrix is not positive semidefinite at LDLT pivot column {pivot_col}: zero diagonal has non-zero coupling at row {row} with value {value}"
            ),
        }
    }
}

impl std::error::Error for LaError {}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use super::*;
    use crate::MAX_STACK_MATRIX_DISPATCH_DIM;

    #[test]
    fn category_displays_are_concise() {
        assert_eq!(FactorizationKind::Lu.to_string(), "LU");
        assert_eq!(FactorizationKind::Ldlt.to_string(), "LDLT");
        assert_eq!(
            ArithmeticOperation::MatrixInfinityNorm.to_string(),
            "matrix infinity norm"
        );
        assert_eq!(
            ArithmeticOperation::SymmetryCheck.to_string(),
            "symmetry check"
        );
        assert_eq!(
            ArithmeticOperation::LuFactorization.to_string(),
            "LU factorization"
        );
        assert_eq!(
            ArithmeticOperation::LdltFactorization.to_string(),
            "LDLT factorization"
        );
        assert_eq!(ArithmeticOperation::LuSolve.to_string(), "LU solve");
        assert_eq!(ArithmeticOperation::LdltSolve.to_string(), "LDLT solve");
        assert_eq!(ArithmeticOperation::Determinant.to_string(), "determinant");
        assert_eq!(
            ArithmeticOperation::DeterminantErrorBound.to_string(),
            "determinant error bound"
        );
        assert_eq!(
            ArithmeticOperation::VectorDotProduct.to_string(),
            "vector dot product"
        );
        assert_eq!(
            ArithmeticOperation::VectorSquaredNorm.to_string(),
            "vector squared norm"
        );
    }

    #[test]
    fn singular_constructors_and_displays_preserve_reason() {
        let exact = LaError::singular_exact(3);
        assert_eq!(
            exact,
            LaError::Singular {
                pivot_col: 3,
                reason: SingularityReason::Exact,
            }
        );
        assert_eq!(
            exact.to_string(),
            "matrix is exactly singular at pivot column 3"
        );

        let numerical = LaError::singular_numerical(2, FactorizationKind::Lu, 1e-14, 1e-12);
        assert_eq!(
            numerical,
            LaError::Singular {
                pivot_col: 2,
                reason: SingularityReason::Numerical {
                    factorization: FactorizationKind::Lu,
                    pivot_magnitude: 1e-14,
                    tolerance: 1e-12,
                },
            }
        );
        assert_eq!(
            numerical.to_string(),
            "matrix is numerically singular during LU factorization at pivot column 2: pivot magnitude 0.00000000000001 <= tolerance 0.000000000001"
        );
    }

    #[test]
    fn non_finite_constructors_preserve_location_and_origin() {
        assert_eq!(
            LaError::non_finite_input_matrix(1, 2),
            LaError::NonFinite {
                location: NonFiniteLocation::MatrixCell { row: 1, col: 2 },
                origin: NonFiniteOrigin::Input,
            }
        );
        assert_eq!(
            LaError::non_finite_input_vector(3).to_string(),
            "non-finite input value at vector entry 3"
        );
        assert_eq!(
            LaError::non_finite_input_scalar().to_string(),
            "non-finite scalar input"
        );
        assert_eq!(
            LaError::non_finite_computation_matrix(ArithmeticOperation::LuFactorization, 2, 1)
                .to_string(),
            "non-finite value computed during LU factorization at matrix cell (2, 1)"
        );
        assert_eq!(
            LaError::non_finite_computation_step(ArithmeticOperation::LuSolve, 1).to_string(),
            "non-finite value computed during LU solve at step 1"
        );
        assert_eq!(
            LaError::non_finite_computation_scalar(ArithmeticOperation::Determinant).to_string(),
            "non-finite scalar result computed during determinant"
        );
    }

    #[test]
    fn unrepresentable_helpers_preserve_recovery_reason() {
        let rounding = LaError::unrepresentable(Some(2), UnrepresentableReason::RequiresRounding);
        let not_finite = LaError::unrepresentable(None, UnrepresentableReason::NotFinite);
        assert_eq!(
            rounding.unrepresentable_reason(),
            Some(UnrepresentableReason::RequiresRounding)
        );
        assert!(rounding.requires_rounding());
        assert_eq!(
            not_finite.to_string(),
            "exact result has no finite f64 representation after rounding"
        );
        assert!(!not_finite.requires_rounding());
        assert_eq!(LaError::singular_exact(0).unrepresentable_reason(), None);
    }

    #[test]
    fn invalid_tolerance_classifies_non_finite_before_negative() {
        assert_eq!(
            LaError::invalid_tolerance(-1.0),
            LaError::InvalidTolerance {
                value: -1.0,
                reason: InvalidToleranceReason::Negative,
            }
        );
        assert_eq!(
            LaError::invalid_tolerance(f64::NEG_INFINITY),
            LaError::InvalidTolerance {
                value: f64::NEG_INFINITY,
                reason: InvalidToleranceReason::NotFinite,
            }
        );
        assert_eq!(
            LaError::invalid_tolerance(-1.0).to_string(),
            "invalid tolerance -1; expected value >= 0"
        );
    }

    #[test]
    fn asymmetric_error_retains_observed_values_and_bound() {
        let err = LaError::asymmetric(0, 2, 3, 1.0, 1.5, 1e-12);
        assert_eq!(
            err,
            LaError::Asymmetric {
                row: 0,
                col: 2,
                dim: 3,
                upper: 1.0,
                lower: 1.5,
                allowed_abs_diff: 1e-12,
            }
        );
        assert_eq!(
            err.to_string(),
            "matrix is not symmetric for dimension 3: entry (0, 2) = 1 and entry (2, 0) = 1.5 differ by more than allowed absolute difference 0.000000000001"
        );
    }

    #[test]
    fn positive_semidefinite_errors_preserve_distinct_violations() {
        assert_eq!(
            LaError::not_positive_semidefinite_negative(1, -3.0).to_string(),
            "matrix is not positive semidefinite at LDLT pivot column 1: diagonal value -3 < 0"
        );
        assert_eq!(
            LaError::not_positive_semidefinite_zero_coupling(0, 1, 2.0).to_string(),
            "matrix is not positive semidefinite at LDLT pivot column 0: zero diagonal has non-zero coupling at row 1 with value 2"
        );
    }

    #[test]
    fn remaining_helpers_and_displays_preserve_fields() {
        assert_eq!(
            LaError::determinant_scale_overflow(3, -1074).to_string(),
            "exact determinant scale exponent overflows for dimension 3 with minimum entry exponent -1074"
        );
        assert_eq!(
            LaError::unsupported_dimension(8, MAX_STACK_MATRIX_DISPATCH_DIM).to_string(),
            "unsupported matrix dimension 8; maximum stack-dispatch dimension is 7"
        );
        assert_eq!(
            LaError::index_out_of_bounds(3, 0, 3).to_string(),
            "matrix index (3, 0) is out of bounds for dimension 3"
        );
    }

    #[test]
    fn is_std_error_with_no_source() {
        let err = LaError::singular_exact(0);
        let error: &dyn Error = &err;
        assert!(error.source().is_none());
    }
}
