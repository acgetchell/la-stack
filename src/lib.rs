#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![doc = include_str!("../README.md")]

#[cfg(doc)]
mod readme_doctests {
    //! Executable versions of README examples.
    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // This system requires pivoting (a[0][0] = 0), so it's a good LU demo.
    /// let a = Matrix::<5>::try_from_rows([
    ///     [0.0, 1.0, 1.0, 1.0, 1.0],
    ///     [1.0, 0.0, 1.0, 1.0, 1.0],
    ///     [1.0, 1.0, 0.0, 1.0, 1.0],
    ///     [1.0, 1.0, 1.0, 0.0, 1.0],
    ///     [1.0, 1.0, 1.0, 1.0, 0.0],
    /// ])?;
    ///
    /// let b = Vector::<5>::try_new([14.0, 13.0, 12.0, 11.0, 10.0])?;
    ///
    /// let lu = a.lu(DEFAULT_SINGULAR_TOL)?;
    /// let x = lu.solve(b)?.into_array();
    ///
    /// // Floating-point rounding is expected; compare with a tolerance.
    /// let expected = [1.0, 2.0, 3.0, 4.0, 5.0];
    /// for (x_i, e_i) in x.iter().zip(expected.iter()) {
    ///     assert!((*x_i - *e_i).abs() <= 1e-12);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    fn solve_5x5_example() {}

    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // This matrix is symmetric positive-definite (A = L*L^T) so LDLT works without pivoting.
    /// let a = Matrix::<5>::try_from_rows([
    ///     [1.0, 1.0, 0.0, 0.0, 0.0],
    ///     [1.0, 2.0, 1.0, 0.0, 0.0],
    ///     [0.0, 1.0, 2.0, 1.0, 0.0],
    ///     [0.0, 0.0, 1.0, 2.0, 1.0],
    ///     [0.0, 0.0, 0.0, 1.0, 2.0],
    /// ])?;
    ///
    /// let ldlt = match a.ldlt(DEFAULT_SINGULAR_TOL) {
    ///     Ok(ldlt) => ldlt,
    ///     Err(err @ LaError::Asymmetric { row, col, .. }) => {
    ///         eprintln!("LDLT requires symmetry; first mismatch at ({row}, {col})");
    ///         return Err(err);
    ///     }
    ///     Err(err) => return Err(err),
    /// };
    ///
    /// let det = ldlt.det()?;
    /// assert!((det - 1.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    fn det_5x5_ldlt_example() {}

    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// // Evaluated entirely at compile time — no runtime cost.
    /// const DET: Result<Option<f64>, LaError> = match Matrix::<4>::try_from_rows([
    ///     [2.0, 0.0, 0.0, 0.0],
    ///     [0.0, 3.0, 0.0, 0.0],
    ///     [0.0, 0.0, 5.0, 0.0],
    ///     [0.0, 0.0, 0.0, 7.0],
    /// ]) {
    ///     Ok(matrix) => matrix.det_direct(),
    ///     Err(err) => Err(err),
    /// };
    ///
    /// # fn main() -> Result<(), LaError> {
    /// assert_eq!(DET?, Some(210.0));
    /// # Ok(())
    /// # }
    /// ```
    fn det_direct_4x4_const_example() {}

    #[cfg(feature = "exact")]
    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // Exact determinant
    /// let m = Matrix::<3>::try_from_rows([
    ///     [1.0, 2.0, 3.0],
    ///     [4.0, 5.0, 6.0],
    ///     [7.0, 8.0, 9.0],
    /// ])?;
    /// assert_eq!(m.det_sign_exact()?, 0); // exactly singular
    ///
    /// let det = m.det_exact()?;
    /// assert_eq!(det, BigRational::from_integer(0.into())); // exact zero
    /// let det_f64 = m.det_exact_f64()?;
    /// assert_eq!(det_f64, 0.0);
    ///
    /// // If the exact determinant cannot fit in f64, keep the BigRational value.
    /// let big = f64::MAX / 2.0;
    /// let huge = Matrix::<3>::try_from_rows([
    ///     [0.0, 0.0, 1.0],
    ///     [big, 0.0, 1.0],
    ///     [0.0, big, 1.0],
    /// ])?;
    /// let huge_det = huge.det_exact()?;
    /// assert_eq!(huge.det_exact_f64(), Err(LaError::Overflow { index: None }));
    /// println!("exact determinant = {huge_det}");
    ///
    /// // Exact linear system solve
    /// let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let b = Vector::<2>::try_new([5.0, 11.0])?;
    /// let x = a.solve_exact_f64(b)?.into_array();
    /// assert!((x[0] - 1.0).abs() <= f64::EPSILON);
    /// assert!((x[1] - 2.0).abs() <= f64::EPSILON);
    /// # Ok(())
    /// # }
    /// ```
    fn exact_arithmetic_example() {}
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
/// This constant is not a caller-tuned tolerance. It is the dimension-specific
/// multiplier that turns the matrix's absolute Leibniz sum into a conservative
/// bound on floating-point roundoff in the closed-form 2×2 determinant formula.
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
/// `Matrix::det_sign_exact` (under the `exact` feature) for the reference
/// adaptive-precision filter.
///
/// # Example
/// ```
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
/// let Some(det) = m.det_direct()? else {
///     return Ok(());
/// };
/// assert_eq!(det, -2.0);
/// // Compute the bound from the raw constant for illustration; most
/// // callers would match on `m.det_errbound()?` instead.
/// let p = (1.0_f64 * 4.0).abs() + (2.0_f64 * 3.0).abs();
/// let bound = ERR_COEFF_2 * p;
/// if det.abs() > bound {
///     // The f64 sign is provably correct without exact arithmetic.
/// }
/// # Ok(())
/// # }
/// ```
pub const ERR_COEFF_2: f64 = 3.0 * EPS + 16.0 * EPS * EPS;

/// Absolute error coefficient for [`Matrix::<3>::det_direct`](crate::Matrix::det_direct).
///
/// This constant is not a caller-tuned tolerance. It is the dimension-specific
/// multiplier that turns the matrix's absolute Leibniz sum into a conservative
/// bound on floating-point roundoff in the closed-form 3×3 determinant formula.
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
/// This constant is not a caller-tuned tolerance. It is the dimension-specific
/// multiplier that turns the matrix's absolute Leibniz sum into a conservative
/// bound on floating-point roundoff in the closed-form 4×4 determinant formula.
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

/// Largest dimension supported by [`try_with_stack_matrix!`].
///
/// The crate can represent `Matrix<D>` for any compile-time `D`, but runtime
/// dispatch must enumerate a finite set of concrete stack types.  Dimensions
/// `0..=7` cover downstream geometric predicate matrices while keeping the
/// dispatch surface explicit.
pub const MAX_STACK_MATRIX_DISPATCH_DIM: usize = 7;

/// Finite, non-negative tolerance used by numerical predicates and factorizations.
///
/// Construct with [`Tolerance::new`] when accepting raw caller input. Once
/// constructed, the stored value is guaranteed to be finite and `>= 0`, so
/// downstream algorithms do not need to revalidate the tolerance.
///
/// This is the crate-wide tolerance contract: raw negative, NaN, and infinite
/// values are rejected with [`LaError::InvalidTolerance`] at construction time.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tolerance {
    value: f64,
}

impl Tolerance {
    /// Construct a tolerance without checking the raw value.
    ///
    /// This crate-internal escape hatch is only for constants whose finite,
    /// non-negative value is visible at the call site. Public callers should
    /// use [`Tolerance::new`] so the returned value carries the validation
    /// proof.
    pub(crate) const fn new_unchecked(value: f64) -> Self {
        Self { value }
    }

    /// Construct a finite, non-negative tolerance.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let tol = Tolerance::new(1e-12)?;
    /// assert_eq!(tol.get(), 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::InvalidTolerance`] when `value` is NaN, infinite, or
    /// negative.
    #[inline]
    pub const fn new(value: f64) -> Result<Self, LaError> {
        if value >= 0.0 && value.is_finite() {
            Ok(Self::new_unchecked(value))
        } else {
            Err(LaError::invalid_tolerance(value))
        }
    }

    /// Return the raw finite, non-negative tolerance value.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let tol = Tolerance::new(0.0)?;
    /// assert_eq!(tol.get(), 0.0);
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    #[must_use]
    pub const fn get(self) -> f64 {
        self.value
    }
}

/// Default absolute threshold used for singularity/degeneracy detection.
///
/// This is intentionally conservative for geometric predicates and small systems.
///
/// Conceptually, this is an absolute bound for deciding when a scalar should be treated
/// as "numerically zero" (e.g. LU pivots, LDLT diagonal entries).
pub const DEFAULT_SINGULAR_TOL: Tolerance = Tolerance::new_unchecked(1e-12);

/// Relative tolerance used to validate matrices for LDLT factorization.
pub(crate) const LDLT_SYMMETRY_REL_TOL: Tolerance = Tolerance::new_unchecked(1e-12);

/// Linear algebra errors.
///
/// This enum is `#[non_exhaustive]` — downstream `match` arms must include a
/// wildcard (`_`) pattern to compile, allowing new variants to be added in
/// future minor releases without breaking existing code.
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
    ///   determinant product, tolerance-scale accumulator, solve accumulator, or
    ///   other scalar/intermediate that has no matrix row coordinate.
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
    /// Use this for non-finite values read from a stored `Matrix` entry or
    /// factorization cell, and for non-finite factorization updates that would
    /// be stored at `(row, col)` if accepted.  The resulting error has
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
    /// Use this for non-finite values in a `Vector` input, determinant scalar,
    /// tolerance-scale accumulator, or solve accumulator that overflowed during
    /// forward/back substitution.  The resulting error has `row: None, col`,
    /// matching the vector/scalar-intermediate convention documented on
    /// [`NonFinite`](Self::NonFinite).  For stored matrix cells or computed
    /// factorization updates tied to a matrix cell, use
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
            Self::Overflow { index: Some(i) } => {
                write!(
                    f,
                    "exact result overflows the target representation at index {i}"
                )
            }
            Self::Overflow { index: None } => {
                write!(f, "exact result overflows the target representation")
            }
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

pub use ldlt::Ldlt;
pub use lu::Lu;
pub use matrix::Matrix;
pub use vector::Vector;

/// Fallibly dispatch a runtime dimension to a concrete stack-allocated matrix.
///
/// The macro creates a zero matrix with type `Matrix<N>` for the selected
/// runtime dimension `N`, then evaluates the supplied closure body.  Supported
/// runtime dimensions run from `0` through [`MAX_STACK_MATRIX_DISPATCH_DIM`].
/// Unsupported dimensions return
/// `Err(LaError::UnsupportedDimension { requested, max })` converted with
/// `From<LaError>`, so downstream crates can use their own public error type.
///
/// # Errors
/// Returns [`LaError::UnsupportedDimension`] (converted through `From<LaError>`)
/// when the requested runtime dimension is greater than
/// [`MAX_STACK_MATRIX_DISPATCH_DIM`].  The closure body may return any other
/// error representable by its declared `Result` type.
///
/// # Examples
/// ```
/// use la_stack::prelude::*;
///
/// # fn main() -> Result<(), LaError> {
/// let requested = 2usize;
/// let det = try_with_stack_matrix!(requested, |mut m| -> Result<f64, LaError> {
///     m.set_checked(0, 0, 1.0)?;
///     m.set_checked(1, 1, 1.0)?;
///     m.det()
/// })?;
///
/// assert_eq!(det, 1.0);
/// # Ok(())
/// # }
/// ```
#[macro_export]
macro_rules! try_with_stack_matrix {
    ($dim:expr, |$matrix:ident| -> $ret:ty $body:block $(,)?) => {{
        let __la_stack_requested_dim: usize = $dim;
        match __la_stack_requested_dim {
            0 => $crate::try_with_stack_matrix!(@arm 0, $matrix, $ret, $body),
            1 => $crate::try_with_stack_matrix!(@arm 1, $matrix, $ret, $body),
            2 => $crate::try_with_stack_matrix!(@arm 2, $matrix, $ret, $body),
            3 => $crate::try_with_stack_matrix!(@arm 3, $matrix, $ret, $body),
            4 => $crate::try_with_stack_matrix!(@arm 4, $matrix, $ret, $body),
            5 => $crate::try_with_stack_matrix!(@arm 5, $matrix, $ret, $body),
            6 => $crate::try_with_stack_matrix!(@arm 6, $matrix, $ret, $body),
            7 => $crate::try_with_stack_matrix!(@arm 7, $matrix, $ret, $body),
            requested => Err(::core::convert::From::from(
                $crate::LaError::unsupported_dimension(
                    requested,
                    $crate::MAX_STACK_MATRIX_DISPATCH_DIM,
                ),
            )),
        }
    }};
    ($dim:expr, |mut $matrix:ident| -> $ret:ty $body:block $(,)?) => {{
        let __la_stack_requested_dim: usize = $dim;
        match __la_stack_requested_dim {
            0 => $crate::try_with_stack_matrix!(@arm_mut 0, $matrix, $ret, $body),
            1 => $crate::try_with_stack_matrix!(@arm_mut 1, $matrix, $ret, $body),
            2 => $crate::try_with_stack_matrix!(@arm_mut 2, $matrix, $ret, $body),
            3 => $crate::try_with_stack_matrix!(@arm_mut 3, $matrix, $ret, $body),
            4 => $crate::try_with_stack_matrix!(@arm_mut 4, $matrix, $ret, $body),
            5 => $crate::try_with_stack_matrix!(@arm_mut 5, $matrix, $ret, $body),
            6 => $crate::try_with_stack_matrix!(@arm_mut 6, $matrix, $ret, $body),
            7 => $crate::try_with_stack_matrix!(@arm_mut 7, $matrix, $ret, $body),
            requested => Err(::core::convert::From::from(
                $crate::LaError::unsupported_dimension(
                    requested,
                    $crate::MAX_STACK_MATRIX_DISPATCH_DIM,
                ),
            )),
        }
    }};
    (@arm $d:literal, $matrix:ident, $ret:ty, $body:block) => {{
        let __la_stack_body = |$matrix: $crate::Matrix<$d>| -> $ret { $body };
        __la_stack_body($crate::Matrix::<$d>::zero())
    }};
    (@arm_mut $d:literal, $matrix:ident, $ret:ty, $body:block) => {{
        let __la_stack_body = |mut $matrix: $crate::Matrix<$d>| -> $ret { $body };
        __la_stack_body($crate::Matrix::<$d>::zero())
    }};
}

/// Common imports for ergonomic usage.
///
/// This prelude re-exports the primary types and constants: [`Matrix`],
/// [`Vector`], [`Lu`], [`Ldlt`], [`Tolerance`], [`LaError`],
/// [`DEFAULT_SINGULAR_TOL`], and the determinant error bound coefficients
/// [`ERR_COEFF_2`], [`ERR_COEFF_3`], and [`ERR_COEFF_4`]. It also re-exports
/// [`MAX_STACK_MATRIX_DISPATCH_DIM`] and
/// [`try_with_stack_matrix!`] for runtime-to-const matrix dispatch.
///
/// When the `exact` feature is enabled, `BigInt` and `BigRational` are also
/// re-exported so callers can construct exact values (e.g. as the expected
/// result of `Matrix::det_exact`) without adding `num-bigint` / `num-rational`
/// to their own dependencies. The most commonly needed `num-traits` items are
/// re-exported alongside them: `FromPrimitive` for `BigRational::from_f64` /
/// `from_i64`, `ToPrimitive` for `BigRational::to_f64` / `to_i64`, and `Signed`
/// for `.is_positive()` / `.is_negative()` / `.abs()`.
pub mod prelude {
    pub use crate::{
        DEFAULT_SINGULAR_TOL, ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4, LaError, Ldlt, Lu,
        MAX_STACK_MATRIX_DISPATCH_DIM, Matrix, Tolerance, Vector, try_with_stack_matrix,
    };

    #[cfg(feature = "exact")]
    pub use crate::{BigInt, BigRational, FromPrimitive, Signed, ToPrimitive};
}

#[cfg(test)]
mod tests {
    use core::assert_matches;

    use super::*;

    use approx::assert_abs_diff_eq;

    #[test]
    fn default_singular_tol_is_expected() {
        assert_abs_diff_eq!(DEFAULT_SINGULAR_TOL.get(), 1e-12, epsilon = 0.0);
    }

    #[test]
    fn tolerance_new_accepts_finite_non_negative_values() {
        assert_eq!(
            Tolerance::new(0.0).unwrap().get().to_bits(),
            0.0f64.to_bits()
        );
        assert_eq!(
            Tolerance::new(1e-12).unwrap().get().to_bits(),
            1e-12f64.to_bits()
        );
        assert_eq!(
            Tolerance::new(f64::MAX).unwrap().get().to_bits(),
            f64::MAX.to_bits()
        );
    }

    #[test]
    fn tolerance_new_rejects_negative_nan_and_infinity() {
        assert_eq!(
            Tolerance::new(-1.0),
            Err(LaError::InvalidTolerance { value: -1.0 })
        );
        assert_matches!(
            Tolerance::new(f64::NAN),
            Err(LaError::InvalidTolerance { value }) if value.is_nan()
        );
        assert_eq!(
            Tolerance::new(f64::INFINITY),
            Err(LaError::InvalidTolerance {
                value: f64::INFINITY,
            })
        );
        assert_eq!(
            Tolerance::new(f64::NEG_INFINITY),
            Err(LaError::InvalidTolerance {
                value: f64::NEG_INFINITY,
            })
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
            max: MAX_STACK_MATRIX_DISPATCH_DIM,
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

    #[test]
    fn prelude_reexports_compile_and_work() -> Result<(), LaError> {
        use crate::prelude::*;

        // Use the items so we know they are in scope and usable.
        let m = Matrix::<2>::identity();
        let v = Vector::<2>::try_new([1.0, 2.0])?;
        assert_abs_diff_eq!(m.inf_norm()?, 1.0, epsilon = 0.0);
        assert_abs_diff_eq!(v.norm2_sq()?, 5.0, epsilon = 0.0);
        let _ = m.lu(DEFAULT_SINGULAR_TOL)?.solve(v)?;
        let _ = m.ldlt(DEFAULT_SINGULAR_TOL)?.solve(v)?;
        assert_eq!(MAX_STACK_MATRIX_DISPATCH_DIM, 7);
        Ok(())
    }

    macro_rules! gen_stack_matrix_dispatch_tests {
        ($d:literal) => {
            pastey::paste! {
                #[test]
                fn [<try_with_stack_matrix_dispatches_ $d d>]() {
                    let requested = $d;
                    let got = try_with_stack_matrix!(requested, |mut m| -> Result<usize, LaError> {
                        if $d > 0 {
                            m.set_checked($d - 1, $d - 1, f64::from($d))?;
                            assert_abs_diff_eq!(
                                m.get_checked($d - 1, $d - 1)?,
                                f64::from($d),
                                epsilon = 0.0
                            );
                        }
                        Ok($d)
                    });

                    assert_eq!(got, Ok($d));
                }
            }
        };
    }

    gen_stack_matrix_dispatch_tests!(2);
    gen_stack_matrix_dispatch_tests!(3);
    gen_stack_matrix_dispatch_tests!(4);
    gen_stack_matrix_dispatch_tests!(5);
    gen_stack_matrix_dispatch_tests!(6);
    gen_stack_matrix_dispatch_tests!(7);

    #[test]
    fn try_with_stack_matrix_supports_zero_dimension() {
        let got = try_with_stack_matrix!(0usize, |m| -> Result<Option<f64>, LaError> {
            m.det_direct()
        });

        assert_eq!(got, Ok(Some(1.0)));
    }

    #[test]
    fn try_with_stack_matrix_reports_unsupported_dimension() {
        let got = try_with_stack_matrix!(8usize, |m| -> Result<f64, LaError> { m.det() });

        assert_eq!(
            got,
            Err(LaError::UnsupportedDimension {
                requested: 8,
                max: MAX_STACK_MATRIX_DISPATCH_DIM,
            })
        );
    }

    #[derive(Debug, PartialEq)]
    struct DownstreamError(LaError);

    impl From<LaError> for DownstreamError {
        fn from(err: LaError) -> Self {
            Self(err)
        }
    }

    #[test]
    fn try_with_stack_matrix_converts_unsupported_dimension_error() {
        let got = try_with_stack_matrix!(9usize, |m| -> Result<usize, DownstreamError> {
            assert_abs_diff_eq!(m.inf_norm()?, 0.0, epsilon = 0.0);
            Ok(0)
        });

        assert_eq!(
            got,
            Err(DownstreamError(LaError::UnsupportedDimension {
                requested: 9,
                max: MAX_STACK_MATRIX_DISPATCH_DIM,
            }))
        );
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
        let half = BigRational::new(BigInt::from(1), BigInt::from(2));
        let two = BigRational::from_integer(BigInt::from(2));
        assert_eq!(BigRational::from_f64(0.5), Some(half.clone()));
        assert_eq!(BigRational::from_i64(2), Some(two.clone()));
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
        assert_eq!(half.to_f64(), Some(0.5));
        assert_eq!(two.to_i64(), Some(2));
    }
}
