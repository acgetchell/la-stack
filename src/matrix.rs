#![forbid(unsafe_code)]

//! Fixed-size, stack-allocated square matrices.

use core::hint::cold_path;

use crate::ldlt::Ldlt;
use crate::lu::Lu;
use crate::{ArithmeticOperation, ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4, LaError, Tolerance};

/// A closed-form determinant and its certified absolute error bound.
///
/// Values of this type are produced by
/// [`Matrix::det_direct_with_errbound`]. The paired result guarantees that the
/// determinant was evaluated once and that its matching bound was computed for
/// the same matrix in one call. The guarantee is unavailable when gradual
/// underflow could invalidate the relative-error analysis or when the matrix
/// dimension exceeds the closed-form D ≤ 4 scope.
#[must_use]
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeterminantWithErrorBound {
    determinant: f64,
    absolute_error_bound: f64,
}

impl DeterminantWithErrorBound {
    /// Return the closed-form determinant approximation.
    #[inline]
    #[must_use]
    pub const fn determinant(self) -> f64 {
        self.determinant
    }

    /// Return the certified absolute error bound.
    ///
    /// The exact determinant lies in
    /// `[determinant - bound, determinant + bound]`.
    #[inline]
    #[must_use]
    pub const fn absolute_error_bound(self) -> f64 {
        self.absolute_error_bound
    }
}

/// Finite fixed-size square matrix `D×D`, stored inline.
///
/// `Matrix` is designed for small, robustness-sensitive systems where stack
/// allocation and const-generic dimensions are useful. For large, dynamic, sparse,
/// or parallel workloads, prefer a broader linear-algebra crate such as
/// [`nalgebra`](https://crates.io/crates/nalgebra) or
/// [`faer`](https://crates.io/crates/faer).
///
/// Public construction and mutation reject NaN and infinity through
/// [`try_from_rows`](Self::try_from_rows) and [`set`](Self::set). The storage
/// field is private, so a
/// `Matrix` value carries the invariant that every stored entry is finite.
/// Algorithms therefore do not re-scan stored entries at every use; user-visible
/// non-finite errors come from construction/mutation boundaries or from values
/// computed during arithmetic, such as overflowed elimination or determinant
/// intermediates.
///
/// Direct field construction is intentionally unavailable to downstream callers:
///
/// ```compile_fail
/// use la_stack::Matrix;
///
/// let _ = Matrix::<2> {
///     rows: [[1.0, f64::NAN], [0.0, 1.0]],
/// };
/// ```
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Matrix<const D: usize> {
    rows: [[f64; D]; D],
}

/// A finite [`Matrix`] proven exactly symmetric for LDLT factorization.
///
/// Mirrored entries have equal numeric values; IEEE-754 signed zeros may have
/// different bit patterns because `+0.0 == -0.0`.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SymmetricMatrix<const D: usize> {
    matrix: Matrix<D>,
}

/// Rounded arithmetic result together with proof that gradual underflow could
/// not have changed that operation's result.
///
/// The determinant filter may only use its relative-error coefficients while
/// every rounded operation in both the determinant and absolute-Leibniz trees
/// stays in the normal range. Exact structural zeros are safe; cancellation to
/// zero is conservatively treated as inconclusive.
#[derive(Clone, Copy, Debug, PartialEq)]
struct FilterArithmetic<const TRACK_UNDERFLOW: bool> {
    value: f64,
    underflow_safe: bool,
}

impl<const TRACK_UNDERFLOW: bool> FilterArithmetic<TRACK_UNDERFLOW> {
    /// Return whether a rounded result is normal or non-finite.
    ///
    /// A single exponent-field test keeps the overwhelmingly common normal
    /// path cheap. Callers inspect operands only when the result is zero or
    /// subnormal so they can distinguish structural zero from range loss.
    #[expect(
        clippy::inline_always,
        reason = "determinant hot-path specialization must eliminate unused safety state"
    )]
    #[inline(always)]
    const fn has_nonzero_exponent(value: f64) -> bool {
        value.to_bits() & 0x7ff0_0000_0000_0000 != 0
    }

    /// Ordinary floating-point multiplication.
    #[expect(
        clippy::inline_always,
        reason = "determinant hot-path specialization must eliminate unused safety state"
    )]
    #[inline(always)]
    const fn multiply(lhs: f64, rhs: f64) -> Self {
        let value = lhs * rhs;
        Self {
            value,
            underflow_safe: !TRACK_UNDERFLOW
                || Self::has_nonzero_exponent(value)
                || lhs == 0.0
                || rhs == 0.0,
        }
    }

    /// Ordinary addition of the non-negative terms used by the error-bound tree.
    #[expect(
        clippy::inline_always,
        reason = "determinant hot-path specialization must eliminate unused safety state"
    )]
    #[inline(always)]
    const fn add_non_negative(lhs: f64, rhs: f64) -> Self {
        let value = lhs + rhs;
        Self {
            value,
            underflow_safe: !TRACK_UNDERFLOW
                || Self::has_nonzero_exponent(value)
                || (lhs == 0.0 && rhs == 0.0),
        }
    }

    /// Fused multiply-add.
    #[expect(
        clippy::inline_always,
        reason = "determinant hot-path specialization must eliminate unused safety state"
    )]
    #[inline(always)]
    const fn mul_add(lhs: f64, rhs: f64, addend: f64) -> Self {
        let value = lhs.mul_add(rhs, addend);
        Self {
            value,
            underflow_safe: !TRACK_UNDERFLOW
                || Self::has_nonzero_exponent(value)
                || ((lhs == 0.0 || rhs == 0.0) && addend == 0.0),
        }
    }
}

/// A finite D=4 matrix proven safe for shared-minor determinant and permanent
/// evaluation.
///
/// Construction proves both the fixed dimension and that every coefficient in
/// the first two rows is non-zero. The latter makes every shared 2×2 minor part
/// of an active Leibniz term, so the dense kernel cannot evaluate an overflowing
/// minor solely for a mathematically absent term.
#[repr(transparent)]
#[derive(Clone, Copy)]
struct Det4SharedMinorInput<'a, const D: usize> {
    matrix: &'a Matrix<D>,
}

impl<'a, const D: usize> Det4SharedMinorInput<'a, D> {
    /// Parse a matrix into the shared-minor D=4 domain.
    ///
    /// `None` selects the guarded determinant path; it does not represent an
    /// invalid public matrix.
    #[expect(
        clippy::inline_always,
        reason = "the D=4 determinant hot path must eliminate its proof wrapper"
    )]
    #[inline(always)]
    const fn try_new(matrix: &'a Matrix<D>) -> Option<Self> {
        if D != 4 {
            return None;
        }

        let r = &matrix.rows;
        let shared_minors_are_active = (r[0][0] != 0.0)
            && (r[0][1] != 0.0)
            && (r[0][2] != 0.0)
            && (r[0][3] != 0.0)
            && (r[1][0] != 0.0)
            && (r[1][1] != 0.0)
            && (r[1][2] != 0.0)
            && (r[1][3] != 0.0);

        if shared_minors_are_active {
            Some(Self { matrix })
        } else {
            None
        }
    }
}

impl<const D: usize> SymmetricMatrix<D> {
    /// Consume the wrapper and return the underlying matrix.
    #[inline]
    pub(crate) const fn into_matrix(self) -> Matrix<D> {
        self.matrix
    }

    /// Construct a symmetric matrix proof without checking the invariant.
    ///
    /// This constructor is only for paths that have already validated exact
    /// mirrored-entry equality with the same predicate as
    /// [`try_new`](Self::try_new). Finiteness is carried by [`Matrix`].
    #[inline]
    const fn new_unchecked(matrix: Matrix<D>) -> Self {
        Self { matrix }
    }

    /// Validate that every mirrored pair has exactly the same finite value.
    ///
    /// IEEE-754 signed zeros compare equal, so `+0.0` and `-0.0` satisfy this
    /// mathematical-symmetry proof even though their bit patterns differ.
    ///
    /// # Errors
    /// Returns [`LaError::Asymmetric`] with `allowed_abs_diff == 0.0` when the
    /// first off-diagonal pair is not exactly equal.
    #[inline]
    #[expect(
        clippy::float_cmp,
        reason = "LDLT requires exact mirrored-entry equality to factor the supplied operator"
    )]
    fn try_new(matrix: Matrix<D>) -> Result<Self, LaError> {
        for row in 0..D {
            for col in (row + 1)..D {
                let upper = matrix.rows[row][col];
                let lower = matrix.rows[col][row];
                if upper != lower {
                    cold_path();
                    return Err(LaError::asymmetric(row, col, D, upper, lower, 0.0));
                }
            }
        }

        Ok(Self::new_unchecked(matrix))
    }
}

impl<const D: usize> Matrix<D> {
    /// Try to create a finite matrix from row-major storage.
    ///
    /// This is the public raw-storage boundary for matrices. Successful
    /// construction makes the returned [`Matrix`] a finite-storage proof.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// assert_eq!(m.get(0, 1), Some(2.0));
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with matrix coordinates for the first
    /// offending entry in row-major order when `rows` contains NaN or infinity.
    #[inline]
    pub const fn try_from_rows(rows: [[f64; D]; D]) -> Result<Self, LaError> {
        if let Some((row, col)) = Self::first_non_finite_cell(&rows) {
            Err(LaError::non_finite_input_matrix(row, col))
        } else {
            Ok(Self::from_rows_unchecked(rows))
        }
    }

    /// Construct a matrix without checking that entries are finite.
    ///
    /// This module-private escape hatch is reserved for finite literals and
    /// algorithm outputs whose finite invariant is visible at the call site.
    /// Computed outputs must be validated before becoming observable API values.
    #[inline]
    const fn from_rows_unchecked(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }

    /// Borrow the finite row-major backing array.
    ///
    /// The returned view is tied to this [`Matrix`], so callers can inspect the
    /// canonical storage without copying it or bypassing the finite-value
    /// invariant.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let matrix = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// assert_eq!(matrix.as_rows(), &[[1.0, 2.0], [3.0, 4.0]]);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// A live view keeps the matrix immutably borrowed, so validated mutation
    /// cannot occur until the view is no longer used:
    ///
    /// ```compile_fail
    /// use la_stack::Matrix;
    ///
    /// let mut matrix = Matrix::<2>::identity();
    /// let rows = matrix.as_rows();
    /// assert!(matrix.set(0, 0, 5.0).is_ok());
    /// assert_eq!(rows[0][0], 1.0);
    /// ```
    #[inline]
    #[must_use]
    pub const fn as_rows(&self) -> &[[f64; D]; D] {
        &self.rows
    }

    /// Consume this matrix and return its finite row-major backing array.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let matrix = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// assert_eq!(matrix.into_rows(), [[1.0, 2.0], [3.0, 4.0]]);
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    #[must_use]
    pub const fn into_rows(self) -> [[f64; D]; D] {
        self.rows
    }

    /// All-zeros finite matrix.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let z = Matrix::<2>::zero();
    /// assert_eq!(z.get(1, 1), Some(0.0));
    /// ```
    #[inline]
    pub const fn zero() -> Self {
        Self::from_rows_unchecked([[0.0; D]; D])
    }

    /// Finite identity matrix.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// let i = Matrix::<3>::identity();
    /// assert_eq!(i.get(0, 0), Some(1.0));
    /// assert_eq!(i.get(0, 1), Some(0.0));
    /// assert_eq!(i.get(2, 2), Some(1.0));
    /// ```
    #[inline]
    pub const fn identity() -> Self {
        let mut m = Self::zero();

        let mut i = 0;
        while i < D {
            m.rows[i][i] = 1.0;
            i += 1;
        }

        m
    }

    /// Get a finite element with bounds checking.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// assert_eq!(m.get(1, 0), Some(3.0));
    /// assert_eq!(m.get(2, 0), None);
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    #[must_use]
    pub const fn get(&self, row: usize, col: usize) -> Option<f64> {
        if row < D && col < D {
            Some(self.rows[row][col])
        } else {
            None
        }
    }

    /// Get a finite element, preserving index context on failure.
    ///
    /// Prefer [`get`](Self::get) for const or hot paths that only need
    /// `Option`-style absence.  Use this method at public runtime boundaries
    /// where row, column, and dimension context should survive in a typed error.
    ///
    /// # Examples
    /// ```
    /// use core::assert_matches;
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// assert_eq!(m.try_get(1, 0)?, 3.0);
    /// assert_matches!(
    ///     m.try_get(2, 0),
    ///     Err(LaError::IndexOutOfBounds {
    ///         row: 2,
    ///         col: 0,
    ///         dim: 2,
    ///         ..
    ///     })
    /// );
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::IndexOutOfBounds`] when either index is not `< D`.
    #[inline]
    pub const fn try_get(&self, row: usize, col: usize) -> Result<f64, LaError> {
        if row < D && col < D {
            Ok(self.rows[row][col])
        } else {
            Err(LaError::index_out_of_bounds(row, col, D))
        }
    }

    /// Set a finite element with bounds checking.
    ///
    /// # Examples
    /// ```
    /// use core::assert_matches;
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let mut m = Matrix::<2>::zero();
    /// assert_eq!(m.set(0, 1, 2.5), Ok(()));
    /// assert_eq!(m.get(0, 1), Some(2.5));
    /// assert_matches!(
    ///     m.set(10, 0, 1.0),
    ///     Err(LaError::IndexOutOfBounds {
    ///         row: 10,
    ///         col: 0,
    ///         dim: 2,
    ///         ..
    ///     })
    /// );
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::IndexOutOfBounds`] when either index is not `< D`.
    /// Returns [`LaError::NonFinite`] when `value` is NaN or infinity.
    #[inline]
    pub const fn set(&mut self, row: usize, col: usize, value: f64) -> Result<(), LaError> {
        if row >= D || col >= D {
            return Err(LaError::index_out_of_bounds(row, col, D));
        }
        if !value.is_finite() {
            return Err(LaError::non_finite_input_matrix(row, col));
        }
        self.rows[row][col] = value;
        Ok(())
    }

    /// Infinity norm (maximum absolute row sum).
    ///
    /// # Non-finite handling
    /// Public constructors and setters reject raw non-finite entries, but
    /// `Matrix` values are finite by construction. `inf_norm` returns
    /// [`LaError::NonFinite`] with the matrix cell whose addition first makes a
    /// row sum non-finite.
    ///
    /// Row sums are accumulated in `f64` with ordinary addition.  This method
    /// checks for overflowed accumulators, but it does not provide a certified
    /// absolute rounding bound for the returned norm.
    ///
    /// # Examples
    /// ```
    /// use core::assert_matches;
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([[1.0, -2.0], [3.0, 4.0]])?;
    /// assert!((m.inf_norm()? - 7.0).abs() <= 1e-12);
    ///
    /// // Raw NaN entries are rejected with coordinates.
    /// assert_matches!(
    ///     Matrix::<2>::try_from_rows([[f64::NAN, 1.0], [2.0, 3.0]]),
    ///     Err(LaError::NonFinite {
    ///         location: NonFiniteLocation::MatrixCell { row: 0, col: 0, .. },
    ///         origin: NonFiniteOrigin::Input,
    ///         ..
    ///     })
    /// );
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with matrix coordinates when a row sum
    /// overflows to NaN or infinity.
    #[inline]
    pub const fn inf_norm(&self) -> Result<f64, LaError> {
        let mut max_row_sum: f64 = 0.0;

        let mut r = 0;
        while r < D {
            let row = &self.rows[r];
            let mut row_sum: f64 = 0.0;
            let mut c = 0;
            while c < D {
                row_sum += row[c].abs();
                c += 1;
            }
            if !row_sum.is_finite() {
                cold_path();
                return Err(Self::inf_norm_overflow_error(row, r));
            }
            if row_sum > max_row_sum {
                max_row_sum = row_sum;
            }
            r += 1;
        }

        Ok(max_row_sum)
    }

    /// Replay an overflowed infinity-norm row to locate the first non-finite sum.
    ///
    /// This runs only after the success-path traversal has found a non-finite
    /// completed row sum. Because stored entries are finite and their absolute
    /// values are non-negative, replaying the same additions must find the
    /// first column whose addition overflowed; if every earlier prefix is
    /// finite, the final column is that first failure.
    #[cold]
    const fn inf_norm_overflow_error(row: &[f64; D], row_index: usize) -> LaError {
        let mut row_sum = 0.0;
        let mut col = 0;
        let last_col = D.saturating_sub(1);
        while col < last_col {
            row_sum += row[col].abs();
            if !row_sum.is_finite() {
                return LaError::non_finite_computation_matrix(
                    ArithmeticOperation::MatrixInfinityNorm,
                    row_index,
                    col,
                );
            }
            col += 1;
        }

        LaError::non_finite_computation_matrix(
            ArithmeticOperation::MatrixInfinityNorm,
            row_index,
            last_col,
        )
    }

    /// Returns `true` if the matrix is approximately symmetric within a relative tolerance.
    ///
    /// Two entries `self[r][c]` and `self[c][r]` are considered equal (for the
    /// purposes of symmetry) when
    /// `|self[r][c] - self[c][r]| <= rel_tol * max(1.0, inf_norm(self))`.
    /// This is a diagnostic predicate for applications that have an
    /// approximation-specific symmetry threshold. It is not the precondition
    /// used by [`ldlt`](Self::ldlt), which requires exact mirrored-entry
    /// equality so the returned factors represent the original matrix.
    ///
    /// Use [`first_asymmetry`](Self::first_asymmetry) to locate the first
    /// offending pair when this returns `Ok(false)`.
    ///
    /// The `rel_tol` argument is a [`Tolerance`], so raw caller input must be
    /// finite and non-negative before it can reach this predicate. Use
    /// [`Tolerance::try_new`] when accepting a raw `f64`; negative, NaN, and
    /// infinite tolerances return
    /// [`LaError::InvalidTolerance`].
    ///
    /// # Overflow handling
    /// A finite matrix can return [`LaError::NonFinite`] with matrix coordinates
    /// if computing the scaled symmetry tolerance overflows to NaN or infinity.
    /// If both stored entries are finite but their difference overflows to ±∞,
    /// the pair is reported as asymmetric.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::try_from_rows([[4.0, 2.0], [2.0, 3.0]])?;
    /// let tol = Tolerance::try_new(1e-12)?;
    /// assert!(a.is_symmetric(tol)?);
    ///
    /// let b = Matrix::<2>::try_from_rows([[4.0, 2.0], [3.0, 3.0]])?;
    /// assert!(!b.is_symmetric(tol)?);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with matrix coordinates when computing the
    /// scaled symmetry tolerance overflows to NaN or infinity.
    #[inline]
    pub fn is_symmetric(&self, rel_tol: Tolerance) -> Result<bool, LaError> {
        Ok(self.first_asymmetry(rel_tol)?.is_none())
    }

    /// Returns the indices `(r, c)` (with `r < c`) of the first off-diagonal
    /// pair that violates approximate symmetry, or `None` if the matrix is
    /// symmetric within `rel_tol`.
    ///
    /// Iteration order is row-major over the strict upper triangle, so the
    /// returned indices are the lexicographically smallest such pair.  The
    /// predicate is the same as [`is_symmetric`](Self::is_symmetric):
    /// `|self[r][c] - self[c][r]| <= rel_tol * max(1.0, inf_norm(self))`.
    /// It is intentionally distinct from the exact equality required by
    /// [`ldlt`](Self::ldlt).
    ///
    /// A finite matrix can return [`LaError::NonFinite`] with matrix coordinates
    /// if computing the scaled symmetry tolerance overflows to NaN or infinity.
    /// If both stored entries are finite but their difference overflows to ±∞,
    /// the pair is reported as asymmetric.
    ///
    /// The `rel_tol` argument is a [`Tolerance`], so raw caller input must be
    /// finite and non-negative before it can reach this predicate. Use
    /// [`Tolerance::try_new`] when accepting a raw `f64`; negative, NaN, and
    /// infinite tolerances return
    /// [`LaError::InvalidTolerance`].
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<3>::try_from_rows([
    ///     [1.0, 2.0, 0.0],
    ///     [2.0, 4.0, 5.0],
    ///     [0.0, 6.0, 9.0], // 6.0 breaks symmetry with a[1][2] = 5.0
    /// ])?;
    /// let tol = Tolerance::try_new(1e-12)?;
    /// assert_eq!(a.first_asymmetry(tol)?, Some((1, 2)));
    /// assert_eq!(Matrix::<3>::identity().first_asymmetry(tol)?, None);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with matrix coordinates when computing the
    /// scaled symmetry tolerance overflows to NaN or infinity.
    #[inline]
    pub fn first_asymmetry(&self, rel_tol: Tolerance) -> Result<Option<(usize, usize)>, LaError> {
        let eps = self.symmetry_epsilon(rel_tol)?;
        for r in 0..D {
            for c in (r + 1)..D {
                let upper = self.rows[r][c];
                let lower = self.rows[c][r];

                let diff = (upper - lower).abs();
                if !diff.is_finite() || diff > eps {
                    cold_path();
                    return Ok(Some((r, c)));
                }
            }
        }
        Ok(None)
    }

    /// Compute an LU decomposition with partial pivoting.
    ///
    /// `D = 0` follows the empty-matrix convention: factorization succeeds,
    /// [`Lu::det`](crate::Lu::det) returns `1.0`, and solving a length-zero
    /// right-hand side returns a length-zero [`Vector`](crate::Vector).
    /// Partial pivoting is a practical finite-precision strategy, not a
    /// certified accuracy guarantee; see `REFERENCES.md` \[1-3, 11-12\].
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let lu = a.lu(DEFAULT_SINGULAR_TOL)?;
    ///
    /// let b = Vector::<2>::try_new([5.0, 11.0])?;
    /// let x = lu.solve(b)?.into_array();
    ///
    /// assert!((x[0] - 1.0).abs() <= 1e-12);
    /// assert!((x[1] - 2.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// Empty matrices use the standard empty-product convention:
    ///
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let lu = Matrix::<0>::zero().lu(DEFAULT_SINGULAR_TOL)?;
    ///
    /// assert_eq!(lu.det()?, 1.0);
    /// assert!(lu.solve(Vector::<0>::zero())?.into_array().is_empty());
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// The `tol` argument is a [`Tolerance`], so raw caller input must be
    /// finite and non-negative before it can reach factorization. Use
    /// [`Tolerance::try_new`] when accepting a raw `f64`; negative, NaN, and
    /// infinite tolerances return
    /// [`LaError::InvalidTolerance`].
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if, for some column `k`, the largest-magnitude candidate pivot
    /// in that column satisfies `|pivot| <= tol` (so no numerically usable pivot exists).
    /// Returns [`LaError::NonFinite`] if an elimination intermediate overflows
    /// to NaN/∞ before it can be stored in the returned [`Lu`].
    #[inline]
    pub fn lu(self, tol: Tolerance) -> Result<Lu<D>, LaError> {
        Lu::factor_finite(self, tol)
    }

    /// Compute an LDLT factorization (`A = L D Lᵀ`) without pivoting.
    ///
    /// `D = 0` follows the empty-matrix convention: factorization succeeds,
    /// [`Ldlt::det`](crate::Ldlt::det) returns `1.0`, and solving a length-zero
    /// right-hand side returns a length-zero [`Vector`](crate::Vector).
    ///
    /// This is intended for exactly symmetric positive-definite matrices such
    /// as nonsingular Gram matrices. Computed zero and tolerance-small positive
    /// pivots are diagnosed as singular rather than returned in a usable
    /// factorization. Because pivots are computed in binary64, success is not
    /// an exact proof that the stored matrix is positive definite.
    /// See `REFERENCES.md` \[4-6, 11-12\] for Cholesky/LDLT background and the
    /// pivoted symmetric-indefinite alternative.
    ///
    /// # Symmetry validation
    /// The input matrix `self` must be exactly symmetric: every mirrored pair
    /// must satisfy `self[i][j] == self[j][i]`. IEEE-754 signed zeros compare
    /// equal and are therefore accepted. Exact equality is a correctness
    /// invariant, not merely a performance hint: LDLT reads only the lower
    /// triangle, so accepting an approximate mismatch would factor a different
    /// operator than the matrix supplied by the caller. Asymmetric inputs return
    /// [`LaError::Asymmetric`] with an allowed absolute difference of `0.0`
    /// before factorization starts.
    ///
    /// [`is_symmetric`](Self::is_symmetric) remains available as a
    /// tolerance-based diagnostic, but `Ok(true)` from that method does not
    /// establish this exact LDLT precondition. If you need a general-purpose
    /// factorization for a non-symmetric matrix, use [`lu`](Self::lu) instead.
    ///
    /// The `tol` argument is a [`Tolerance`], so raw caller input must be
    /// finite and non-negative before it can reach factorization. Use
    /// [`Tolerance::try_new`] when accepting a raw `f64`; negative, NaN, and
    /// infinite tolerances return
    /// [`LaError::InvalidTolerance`].
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // Note the symmetric layout: a[0][1] == a[1][0] == 2.0.
    /// let a = Matrix::<2>::try_from_rows([[4.0, 2.0], [2.0, 3.0]])?;
    /// let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL)?;
    ///
    /// // det(A) = 8
    /// assert!((ldlt.det()? - 8.0).abs() <= 1e-12);
    ///
    /// // Solve A x = b
    /// let b = Vector::<2>::try_new([1.0, 2.0])?;
    /// let x = ldlt.solve(b)?.into_array();
    /// assert!((x[0] - (-0.125)).abs() <= 1e-12);
    /// assert!((x[1] - 0.75).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// Empty matrices use the standard empty-product convention:
    ///
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let ldlt = Matrix::<0>::zero().ldlt(DEFAULT_SINGULAR_TOL)?;
    ///
    /// assert_eq!(ldlt.det()?, 1.0);
    /// assert!(ldlt.solve(Vector::<0>::zero())?.into_array().is_empty());
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NotPositiveSemidefinite`] if a pivot is negative or a
    /// zero pivot retains a non-zero coupling below it.
    /// Returns [`LaError::Singular`] if a zero pivot has no remaining coupling,
    /// or if a positive pivot satisfies `d <= tol`, treating PSD degeneracy as
    /// singular.
    /// Returns [`LaError::NonFinite`] if factorization computes a non-finite
    /// intermediate.
    /// Returns [`LaError::Asymmetric`] if the input matrix is not symmetric.
    #[inline]
    pub fn ldlt(self, tol: Tolerance) -> Result<Ldlt<D>, LaError> {
        Ldlt::factor_symmetric(SymmetricMatrix::try_new(self)?, tol)
    }

    /// Return the first non-finite stored cell in row-major order.
    const fn first_non_finite_cell(rows: &[[f64; D]; D]) -> Option<(usize, usize)> {
        let mut r = 0;
        while r < D {
            let mut c = 0;
            while c < D {
                if !rows[r][c].is_finite() {
                    return Some((r, c));
                }
                c += 1;
            }
            r += 1;
        }
        None
    }

    /// Compute the approximate-symmetry tolerance scale for a finite matrix.
    ///
    /// This helper protects the public [`is_symmetric`](Self::is_symmetric) and
    /// [`first_asymmetry`](Self::first_asymmetry) diagnostic contracts: the
    /// documented norm-first formula is used whenever its intermediate is
    /// representable, while an overflow-safe termwise fallback reports the
    /// matrix cell that makes the scaled tolerance non-finite.
    fn symmetry_epsilon(&self, rel_tol: Tolerance) -> Result<f64, LaError> {
        let rel_tol = rel_tol.get();

        if rel_tol == 0.0 {
            return Ok(rel_tol);
        }

        if let Ok(norm) = self.inf_norm() {
            let scale = if norm > 1.0 { norm } else { 1.0 };
            let eps = rel_tol * scale;
            if eps.is_finite() {
                return Ok(eps);
            }
        }

        // If the unscaled row sum or the final multiplication overflows, apply
        // the tolerance to each non-negative contribution before summing. A row
        // can overflow only at magnitudes where multiplication by the smallest
        // positive tolerance is normal, so this fallback cannot introduce the
        // gradual-underflow discrepancy avoided by the direct path above.
        let mut eps = rel_tol;

        for r in 0..D {
            let mut row_eps = 0.0;
            for c in 0..D {
                row_eps = rel_tol.mul_add(self.rows[r][c].abs(), row_eps);
                if !row_eps.is_finite() {
                    cold_path();
                    return Err(LaError::non_finite_computation_matrix(
                        ArithmeticOperation::SymmetryCheck,
                        r,
                        c,
                    ));
                }
            }
            if row_eps > eps {
                eps = row_eps;
            }
        }

        Ok(eps)
    }

    /// Closed-form determinant for dimensions 0–4, bypassing LU factorization.
    ///
    /// Returns `Ok(Some(det))` for `D` ∈ {0, 1, 2, 3, 4}, `Ok(None)` for D ≥ 5.
    /// `D = 0` returns `Ok(Some(1.0))` (empty product).
    /// This is a `const fn` (Rust 1.94+) and uses fused multiply-add (`mul_add`)
    /// for improved accuracy and performance.
    ///
    /// For a determinant that works for any dimension (falling back to LU for D ≥ 5),
    /// use [`det`](Self::det).
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// assert_eq!(m.det_direct()?, Some(-2.0));
    ///
    /// // D = 0 is the empty product.
    /// assert_eq!(Matrix::<0>::zero().det_direct()?, Some(1.0));
    ///
    /// // D ≥ 5 returns None.
    /// assert!(Matrix::<5>::identity().det_direct()?.is_none());
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the closed-form determinant overflows
    /// to NaN or infinity.
    #[inline]
    pub const fn det_direct(&self) -> Result<Option<f64>, LaError> {
        let Some(det) = self.det_direct_arithmetic::<false>() else {
            cold_path();
            return Ok(None);
        };

        Self::computed_scalar_result(ArithmeticOperation::Determinant, det.value)
    }

    /// Evaluate the closed-form determinant while certifying every rounded
    /// operation against gradual underflow.
    #[expect(
        clippy::inline_always,
        reason = "det_direct callers must eliminate unused filter-safety bookkeeping"
    )]
    #[inline(always)]
    const fn det_direct_arithmetic<const TRACK_UNDERFLOW: bool>(
        &self,
    ) -> Option<FilterArithmetic<TRACK_UNDERFLOW>> {
        match D {
            0 => Some(FilterArithmetic {
                value: 1.0,
                underflow_safe: true,
            }),
            1 => Some(FilterArithmetic {
                value: self.rows[0][0],
                underflow_safe: true,
            }),
            2 => {
                let a = self.rows[0][0];
                let b = self.rows[0][1];
                let c = self.rows[1][0];
                let d = self.rows[1][1];
                let subtrahend = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(b, c);
                let mut det = FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(a, d, -subtrahend.value);
                det.underflow_safe &= subtrahend.underflow_safe;
                Some(det)
            }
            3 => Some(Self::det3_elements::<TRACK_UNDERFLOW>(
                [self.rows[0][0], self.rows[0][1], self.rows[0][2]],
                [self.rows[1][0], self.rows[1][1], self.rows[1][2]],
                [self.rows[2][0], self.rows[2][1], self.rows[2][2]],
            )),
            4 => {
                if !TRACK_UNDERFLOW && let Some(input) = Det4SharedMinorInput::try_new(self) {
                    return Some(FilterArithmetic {
                        value: Self::det4_dense_elements(input),
                        underflow_safe: true,
                    });
                }

                let r = &self.rows;
                let mut det = if r[0][3] == 0.0 {
                    FilterArithmetic {
                        value: 0.0,
                        underflow_safe: true,
                    }
                } else {
                    let c03 = Self::det3_elements::<TRACK_UNDERFLOW>(
                        [r[1][0], r[1][1], r[1][2]],
                        [r[2][0], r[2][1], r[2][2]],
                        [r[3][0], r[3][1], r[3][2]],
                    );
                    let mut term =
                        FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r[0][3], c03.value);
                    term.value = -term.value;
                    term.underflow_safe &= c03.underflow_safe;
                    term
                };
                if r[0][2] != 0.0 {
                    let c02 = Self::det3_elements::<TRACK_UNDERFLOW>(
                        [r[1][0], r[1][1], r[1][3]],
                        [r[2][0], r[2][1], r[2][3]],
                        [r[3][0], r[3][1], r[3][3]],
                    );
                    let prior_safe = det.underflow_safe && c02.underflow_safe;
                    det =
                        FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(r[0][2], c02.value, det.value);
                    det.underflow_safe &= prior_safe;
                }
                if r[0][1] != 0.0 {
                    let c01 = Self::det3_elements::<TRACK_UNDERFLOW>(
                        [r[1][0], r[1][2], r[1][3]],
                        [r[2][0], r[2][2], r[2][3]],
                        [r[3][0], r[3][2], r[3][3]],
                    );
                    let prior_safe = det.underflow_safe && c01.underflow_safe;
                    det = FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(
                        -r[0][1], c01.value, det.value,
                    );
                    det.underflow_safe &= prior_safe;
                }
                if r[0][0] != 0.0 {
                    let c00 = Self::det3_elements::<TRACK_UNDERFLOW>(
                        [r[1][1], r[1][2], r[1][3]],
                        [r[2][1], r[2][2], r[2][3]],
                        [r[3][1], r[3][2], r[3][3]],
                    );
                    let prior_safe = det.underflow_safe && c00.underflow_safe;
                    det =
                        FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(r[0][0], c00.value, det.value);
                    det.underflow_safe &= prior_safe;
                }

                Some(det)
            }
            _ => None,
        }
    }

    /// Evaluate the proof-bearing 4×4 cofactor expansion with shared 2×2 minors.
    ///
    /// When no intermediate undergoes gradual underflow, the rounding error is
    /// bounded by `ERR_COEFF_4 · p(|A|)`, where `p(|A|)` is the absolute Leibniz
    /// sum. This helper returns only the determinant; use
    /// [`Self::det_errbound`] or [`Self::det_direct_with_errbound`] to obtain the
    /// certified bound.
    #[expect(
        clippy::inline_always,
        reason = "the D=4 determinant hot path must inline its shared-minor expansion"
    )]
    #[inline(always)]
    const fn det4_dense_elements(input: Det4SharedMinorInput<'_, D>) -> f64 {
        let r = &input.matrix.rows;
        let s23 = r[2][2].mul_add(r[3][3], -(r[2][3] * r[3][2]));
        let s13 = r[2][1].mul_add(r[3][3], -(r[2][3] * r[3][1]));
        let s12 = r[2][1].mul_add(r[3][2], -(r[2][2] * r[3][1]));
        let s03 = r[2][0].mul_add(r[3][3], -(r[2][3] * r[3][0]));
        let s02 = r[2][0].mul_add(r[3][2], -(r[2][2] * r[3][0]));
        let s01 = r[2][0].mul_add(r[3][1], -(r[2][1] * r[3][0]));

        let c00 = r[1][1].mul_add(s23, (-r[1][2]).mul_add(s13, r[1][3] * s12));
        let c01 = r[1][0].mul_add(s23, (-r[1][2]).mul_add(s03, r[1][3] * s02));
        let c02 = r[1][0].mul_add(s13, (-r[1][1]).mul_add(s03, r[1][3] * s01));
        let c03 = r[1][0].mul_add(s12, (-r[1][1]).mul_add(s02, r[1][2] * s01));

        r[0][0].mul_add(
            c00,
            (-r[0][1]).mul_add(c01, r[0][2].mul_add(c02, -(r[0][3] * c03))),
        )
    }

    /// Evaluate the dense 4×4 absolute permanent with shared 2×2 minors.
    ///
    /// The proof carried by `input` makes every shared minor part of an active
    /// Leibniz term. The caller separately establishes a wide exponent margin,
    /// so this branch-free kernel cannot hide gradual underflow or evaluate an
    /// overflowing minor for a mathematically absent term.
    #[expect(
        clippy::inline_always,
        reason = "the D=4 determinant filter must inline its shared-minor permanent"
    )]
    #[inline(always)]
    const fn det4_dense_abs_permanent_elements(input: Det4SharedMinorInput<'_, D>) -> f64 {
        let r = &input.matrix.rows;
        let sp23 = (r[2][2] * r[3][3]).abs() + (r[2][3] * r[3][2]).abs();
        let sp13 = (r[2][1] * r[3][3]).abs() + (r[2][3] * r[3][1]).abs();
        let sp12 = (r[2][1] * r[3][2]).abs() + (r[2][2] * r[3][1]).abs();
        let sp03 = (r[2][0] * r[3][3]).abs() + (r[2][3] * r[3][0]).abs();
        let sp02 = (r[2][0] * r[3][2]).abs() + (r[2][2] * r[3][0]).abs();
        let sp01 = (r[2][0] * r[3][1]).abs() + (r[2][1] * r[3][0]).abs();

        let pc0 = r[1][3]
            .abs()
            .mul_add(sp12, r[1][2].abs().mul_add(sp13, r[1][1].abs() * sp23));
        let pc1 = r[1][3]
            .abs()
            .mul_add(sp02, r[1][2].abs().mul_add(sp03, r[1][0].abs() * sp23));
        let pc2 = r[1][3]
            .abs()
            .mul_add(sp01, r[1][1].abs().mul_add(sp03, r[1][0].abs() * sp13));
        let pc3 = r[1][2]
            .abs()
            .mul_add(sp01, r[1][1].abs().mul_add(sp02, r[1][0].abs() * sp12));

        r[0][3].abs().mul_add(
            pc3,
            r[0][2]
                .abs()
                .mul_add(pc2, r[0][1].abs().mul_add(pc1, r[0][0].abs() * pc0)),
        )
    }

    /// Floating-point determinant, using closed-form formulas for D ≤ 4 and
    /// LU decomposition for D ≥ 5.
    ///
    /// For D ∈ {1, 2, 3, 4}, this bypasses LU factorization entirely for a significant
    /// speedup (see [`det_direct`](Self::det_direct)).
    ///
    /// Because this method mixes closed-form paths from
    /// [`det_direct`](Self::det_direct) with an LU fallback, the returned value has
    /// no certified absolute error bound. Use
    /// [`det_errbound`](Self::det_errbound) for D ≤ 4 bounds, or the exact
    /// determinant APIs when exact singularity classification or certified values
    /// matter. For D ≥ 5, the zero-tolerance LU fallback surfaces
    /// [`LaError::Singular`] when elimination cannot produce a non-zero pivot.
    /// Floating-point elimination cannot in general distinguish an exactly
    /// singular matrix from a non-singular matrix whose intermediate pivot
    /// rounded to zero, so this method never converts that numerical failure into
    /// an exact `0.0` result.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let det = Matrix::<3>::identity().det()?;
    /// assert!((det - 1.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// The LU fallback accumulates its diagonal product with power-of-two
    /// scaling, so factor order cannot cause premature overflow or underflow in
    /// the final product. Elimination intermediates remain subject to binary64
    /// rounding and range limits.
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if the D ≥ 5 LU fallback cannot produce a
    /// non-zero pivot, including when a non-zero mathematical intermediate rounds
    /// to zero during elimination. Returns [`LaError::NonFinite`] if a D ≤ 4
    /// closed-form result is non-finite, if the LU fallback computes a
    /// non-finite factorization cell, or if its final scaled determinant cannot
    /// be represented as a finite `f64`.
    #[inline]
    pub fn det(self) -> Result<f64, LaError> {
        if let Some(d) = self.det_direct()? {
            return Ok(d);
        }
        self.lu(Tolerance::ZERO)?.det()
    }

    /// Evaluate `det_direct()` and its absolute error bound together.
    ///
    /// Returns `Ok(Some(result))` for D ≤ 4 when the relative-error analysis
    /// is valid. The result contains the closed-form determinant and a bound
    /// such that `|result.determinant() - det_exact| ≤
    /// result.absolute_error_bound()`. Returns `Ok(None)` when gradual
    /// underflow could invalidate that analysis or for D ≥ 5, where no
    /// closed-form bound is available.
    ///
    /// This is the preferred API when both values are needed: it evaluates the
    /// determinant arithmetic tree once, then computes the matching bound for
    /// the same matrix within that call.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let matrix = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// if let Some(estimate) = matrix.det_direct_with_errbound()? {
    ///     assert_eq!(estimate.determinant(), -2.0);
    ///     assert!(estimate.absolute_error_bound() >= 0.0);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the determinant or bound computation
    /// overflows to NaN or infinity. Underflow-sensitive finite computations
    /// return `Ok(None)` because they remain valid inputs for an exact fallback.
    #[inline]
    pub const fn det_direct_with_errbound(
        &self,
    ) -> Result<Option<DeterminantWithErrorBound>, LaError> {
        if self.det_bound_inputs_have_wide_exponent_margin() {
            let Some(det) = self.det_direct_arithmetic::<false>() else {
                cold_path();
                return Ok(None);
            };
            return self.det_direct_with_errbound_from_arithmetic(det);
        }

        let Some(det) = self.det_direct_arithmetic::<true>() else {
            cold_path();
            return Ok(None);
        };
        self.det_direct_with_errbound_from_arithmetic(det)
    }

    /// Conservative absolute error bound for `det_direct()`.
    ///
    /// Returns `Ok(Some(bound))` such that `|det_direct() - det_exact| ≤ bound`
    /// when every rounded intermediate used by the closed-form determinant and
    /// bound is normal (or an exact structural zero). Returns `Ok(None)` when
    /// gradual underflow could invalidate the relative-error analysis, or for
    /// D ≥ 5 where no fast bound is available.
    ///
    /// For D ≤ 4, the bound is derived from the absolute Leibniz sum using
    /// Shewchuk-style error analysis (see `REFERENCES.md` \[8\] and the
    /// per-constant docs on [`ERR_COEFF_2`], [`ERR_COEFF_3`], and
    /// [`ERR_COEFF_4`]). For D = 0 or 1, returns
    /// `Some(0.0)` since the determinant computation is exact (no
    /// arithmetic).
    ///
    /// This method does NOT require the `exact` feature — the bounds use
    /// pure f64 arithmetic and are useful for custom adaptive-precision logic.
    ///
    /// # When to use
    ///
    /// Use [`det_direct_with_errbound`](Self::det_direct_with_errbound) when the
    /// determinant and bound are both needed. This accessor is convenient when
    /// only the bound is needed.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<3>::try_from_rows([
    ///     [1.0, 2.0, 3.0],
    ///     [4.0, 5.0, 6.0],
    ///     [7.0, 8.0, 9.0],
    /// ])?;
    /// if let Some(bound) = m.det_errbound()? {
    ///     assert!(bound >= 0.0);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Adaptive precision pattern (requires `exact` feature)
    /// ```ignore
    /// use la_stack::prelude::*;
    ///
    /// fn adaptive_det_sign<const D: usize>(
    ///     matrix: &Matrix<D>,
    /// ) -> DeterminantSign {
    ///     if let Ok(Some(estimate)) = matrix.det_direct_with_errbound() {
    ///         if estimate.determinant().abs() > estimate.absolute_error_bound() {
    ///             return if estimate.determinant() > 0.0 {
    ///                 DeterminantSign::Positive
    ///             } else {
    ///                 DeterminantSign::Negative
    ///             };
    ///         }
    ///     }
    ///
    ///     matrix.det_sign_exact()
    /// }
    ///
    /// fn main() -> Result<(), LaError> {
    ///     assert_eq!(
    ///         adaptive_det_sign(&Matrix::<3>::identity()),
    ///         DeterminantSign::Positive
    ///     );
    ///
    ///     let big = f64::MAX / 2.0;
    ///     let overflowing = Matrix::<3>::try_from_rows([
    ///         [0.0, 0.0, 1.0],
    ///         [big, 0.0, 1.0],
    ///         [0.0, big, 1.0],
    ///     ])?;
    ///     assert_eq!(
    ///         adaptive_det_sign(&overflowing),
    ///         DeterminantSign::Positive
    ///     );
    ///     Ok(())
    /// }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the bound computation overflows to
    /// NaN or infinity. Underflow-sensitive finite computations return
    /// `Ok(None)` instead because they are valid inputs for an exact fallback.
    #[inline]
    pub const fn det_errbound(&self) -> Result<Option<f64>, LaError> {
        match self.det_direct_with_errbound() {
            Ok(Some(result)) => Ok(Some(result.absolute_error_bound)),
            Ok(None) => Ok(None),
            Err(error) => Err(error),
        }
    }

    /// Return whether every non-zero entry is large enough that the complete
    /// D≤4 determinant and permanent trees cannot gradually underflow.
    ///
    /// The `2^-16` threshold leaves hundreds of binary exponent bits of margin
    /// even after the D=4 tree's products, FMAs, and binary64 rounding steps.
    /// Overflow remains possible and is classified after evaluation. Inputs
    /// below this conservative threshold use per-operation tracking instead.
    const fn det_bound_inputs_have_wide_exponent_margin(&self) -> bool {
        const MIN_MAGNITUDE_BITS: u64 = 1007_u64 << 52; // 2^-16
        const MAGNITUDE_MASK: u64 = !(1_u64 << 63);

        if D > 4 {
            return false;
        }

        let mut row = 0;
        while row < D {
            let mut col = 0;
            while col < D {
                let magnitude_bits = self.rows[row][col].to_bits() & MAGNITUDE_MASK;
                if magnitude_bits != 0 && magnitude_bits < MIN_MAGNITUDE_BITS {
                    return false;
                }
                col += 1;
            }
            row += 1;
        }
        true
    }

    /// Classify a completed determinant tree and construct its matching bound.
    const fn det_direct_with_errbound_from_arithmetic<const TRACK_UNDERFLOW: bool>(
        &self,
        det: FilterArithmetic<TRACK_UNDERFLOW>,
    ) -> Result<Option<DeterminantWithErrorBound>, LaError> {
        let bound = match self.det_errbound_from_arithmetic(det) {
            Ok(Some(bound)) => bound,
            Ok(None) => return Ok(None),
            Err(error) => return Err(error),
        };
        if !det.value.is_finite() {
            cold_path();
            return Err(LaError::non_finite_computation_scalar(
                ArithmeticOperation::Determinant,
            ));
        }
        Ok(Some(DeterminantWithErrorBound {
            determinant: det.value,
            absolute_error_bound: bound,
        }))
    }

    /// Compute a bound after the matching determinant tree has been evaluated.
    const fn det_errbound_from_arithmetic<const TRACK_UNDERFLOW: bool>(
        &self,
        det: FilterArithmetic<TRACK_UNDERFLOW>,
    ) -> Result<Option<f64>, LaError> {
        if !det.underflow_safe {
            cold_path();
            return Ok(None);
        }

        match D {
            0 | 1 => Self::computed_scalar_result(ArithmeticOperation::DeterminantErrorBound, 0.0),
            2 => {
                let r = &self.rows;
                let product_0 = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r[0][0], r[1][1]);
                let product_1 = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r[0][1], r[1][0]);
                let mut permanent = FilterArithmetic::<TRACK_UNDERFLOW>::add_non_negative(
                    product_0.value.abs(),
                    product_1.value.abs(),
                );
                permanent.underflow_safe &= product_0.underflow_safe && product_1.underflow_safe;
                Self::certified_error_bound(ERR_COEFF_2, permanent)
            }
            3 => {
                let r = &self.rows;
                let permanent = Self::det3_abs_permanent_elements::<TRACK_UNDERFLOW>(
                    [r[0][0], r[0][1], r[0][2]],
                    [r[1][0], r[1][1], r[1][2]],
                    [r[2][0], r[2][1], r[2][2]],
                );
                Self::certified_error_bound(ERR_COEFF_3, permanent)
            }
            4 => self.det4_errbound::<TRACK_UNDERFLOW>(),
            _ => {
                cold_path();
                Ok(None)
            }
        }
    }

    /// Compute the D=4 determinant error bound after the dimension dispatch.
    const fn det4_errbound<const TRACK_UNDERFLOW: bool>(&self) -> Result<Option<f64>, LaError> {
        if !TRACK_UNDERFLOW && let Some(input) = Det4SharedMinorInput::try_new(self) {
            return Self::certified_error_bound(
                ERR_COEFF_4,
                FilterArithmetic::<TRACK_UNDERFLOW> {
                    value: Self::det4_dense_abs_permanent_elements(input),
                    underflow_safe: true,
                },
            );
        }

        let r = &self.rows;
        let mut permanent = if r[0][3] == 0.0 {
            FilterArithmetic {
                value: 0.0,
                underflow_safe: true,
            }
        } else {
            let pc3 = Self::det3_abs_permanent_elements::<TRACK_UNDERFLOW>(
                [r[1][0], r[1][1], r[1][2]],
                [r[2][0], r[2][1], r[2][2]],
                [r[3][0], r[3][1], r[3][2]],
            );
            let mut term = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r[0][3].abs(), pc3.value);
            term.underflow_safe &= pc3.underflow_safe;
            term
        };
        if r[0][2] != 0.0 {
            let pc2 = Self::det3_abs_permanent_elements::<TRACK_UNDERFLOW>(
                [r[1][0], r[1][1], r[1][3]],
                [r[2][0], r[2][1], r[2][3]],
                [r[3][0], r[3][1], r[3][3]],
            );
            let prior_safe = permanent.underflow_safe && pc2.underflow_safe;
            permanent = FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(
                r[0][2].abs(),
                pc2.value,
                permanent.value,
            );
            permanent.underflow_safe &= prior_safe;
        }
        if r[0][1] != 0.0 {
            let pc1 = Self::det3_abs_permanent_elements::<TRACK_UNDERFLOW>(
                [r[1][0], r[1][2], r[1][3]],
                [r[2][0], r[2][2], r[2][3]],
                [r[3][0], r[3][2], r[3][3]],
            );
            let prior_safe = permanent.underflow_safe && pc1.underflow_safe;
            permanent = FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(
                r[0][1].abs(),
                pc1.value,
                permanent.value,
            );
            permanent.underflow_safe &= prior_safe;
        }
        if r[0][0] != 0.0 {
            let pc0 = Self::det3_abs_permanent_elements::<TRACK_UNDERFLOW>(
                [r[1][1], r[1][2], r[1][3]],
                [r[2][1], r[2][2], r[2][3]],
                [r[3][1], r[3][2], r[3][3]],
            );
            let prior_safe = permanent.underflow_safe && pc0.underflow_safe;
            permanent = FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(
                r[0][0].abs(),
                pc0.value,
                permanent.value,
            );
            permanent.underflow_safe &= prior_safe;
        }
        Self::certified_error_bound(ERR_COEFF_4, permanent)
    }

    /// Evaluate a 3×3 determinant expansion with a guarded sparse fallback.
    ///
    /// When all three first-row coefficients are non-zero, one branch-free
    /// closed form is used. The sparse fallback protects the public
    /// [`det_direct`](Self::det_direct) contract: a mathematically absent term
    /// must not compute an overflowing minor and poison the determinant with
    /// `0.0 * inf == NaN`. Nonzero terms keep the same fused multiply-add
    /// ordering as the closed-form expansion.
    #[expect(
        clippy::inline_always,
        reason = "det_direct callers must eliminate unused filter-safety bookkeeping"
    )]
    #[inline(always)]
    const fn det3_elements<const TRACK_UNDERFLOW: bool>(
        r0: [f64; 3],
        r1: [f64; 3],
        r2: [f64; 3],
    ) -> FilterArithmetic<TRACK_UNDERFLOW> {
        let dense = (r0[0] != 0.0) && (r0[1] != 0.0) && (r0[2] != 0.0);
        if !TRACK_UNDERFLOW && dense {
            let m00 = r1[1].mul_add(r2[2], -(r1[2] * r2[1]));
            let m01 = r1[0].mul_add(r2[2], -(r1[2] * r2[0]));
            let m02 = r1[0].mul_add(r2[1], -(r1[1] * r2[0]));
            return FilterArithmetic {
                value: r0[0].mul_add(m00, (-r0[1]).mul_add(m01, r0[2] * m02)),
                underflow_safe: true,
            };
        }

        let mut det = if r0[2] == 0.0 {
            FilterArithmetic {
                value: 0.0,
                underflow_safe: true,
            }
        } else {
            let subtrahend = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r1[1], r2[0]);
            let mut m02 =
                FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(r1[0], r2[1], -subtrahend.value);
            m02.underflow_safe &= subtrahend.underflow_safe;
            let mut term = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r0[2], m02.value);
            term.underflow_safe &= m02.underflow_safe;
            term
        };
        if r0[1] != 0.0 {
            let subtrahend = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r1[2], r2[0]);
            let mut m01 =
                FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(r1[0], r2[2], -subtrahend.value);
            m01.underflow_safe &= subtrahend.underflow_safe;
            let prior_safe = det.underflow_safe && m01.underflow_safe;
            det = FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(-r0[1], m01.value, det.value);
            det.underflow_safe &= prior_safe;
        }
        if r0[0] != 0.0 {
            let subtrahend = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r1[2], r2[1]);
            let mut m00 =
                FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(r1[1], r2[2], -subtrahend.value);
            m00.underflow_safe &= subtrahend.underflow_safe;
            let prior_safe = det.underflow_safe && m00.underflow_safe;
            det = FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(r0[0], m00.value, det.value);
            det.underflow_safe &= prior_safe;
        }
        det
    }

    /// Evaluate a 3×3 absolute permanent while skipping zero coefficients.
    ///
    /// This mirrors [`det3_elements`](Self::det3_elements) for error-bound
    /// computation: absent determinant terms should not force evaluation of an
    /// overflowing absolute minor.
    #[expect(
        clippy::inline_always,
        reason = "error-bound call-site specialization avoids tracked-helper overhead"
    )]
    #[inline(always)]
    const fn det3_abs_permanent_elements<const TRACK_UNDERFLOW: bool>(
        r0: [f64; 3],
        r1: [f64; 3],
        r2: [f64; 3],
    ) -> FilterArithmetic<TRACK_UNDERFLOW> {
        let dense = (r0[0] != 0.0) && (r0[1] != 0.0) && (r0[2] != 0.0);
        if !TRACK_UNDERFLOW && dense {
            let pm00 = (r1[1] * r2[2]).abs() + (r1[2] * r2[1]).abs();
            let pm01 = (r1[0] * r2[2]).abs() + (r1[2] * r2[0]).abs();
            let pm02 = (r1[0] * r2[1]).abs() + (r1[1] * r2[0]).abs();
            return FilterArithmetic {
                value: r0[2]
                    .abs()
                    .mul_add(pm02, r0[1].abs().mul_add(pm01, r0[0].abs() * pm00)),
                underflow_safe: true,
            };
        }

        let mut permanent = if r0[2] == 0.0 {
            FilterArithmetic {
                value: 0.0,
                underflow_safe: true,
            }
        } else {
            let product_0 = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r1[0], r2[1]);
            let product_1 = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r1[1], r2[0]);
            let mut pm02 = FilterArithmetic::<TRACK_UNDERFLOW>::add_non_negative(
                product_0.value.abs(),
                product_1.value.abs(),
            );
            pm02.underflow_safe &= product_0.underflow_safe && product_1.underflow_safe;
            let mut term = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r0[2].abs(), pm02.value);
            term.underflow_safe &= pm02.underflow_safe;
            term
        };
        if r0[1] != 0.0 {
            let product_0 = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r1[0], r2[2]);
            let product_1 = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r1[2], r2[0]);
            let mut pm01 = FilterArithmetic::<TRACK_UNDERFLOW>::add_non_negative(
                product_0.value.abs(),
                product_1.value.abs(),
            );
            pm01.underflow_safe &= product_0.underflow_safe && product_1.underflow_safe;
            let prior_safe = permanent.underflow_safe && pm01.underflow_safe;
            permanent = FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(
                r0[1].abs(),
                pm01.value,
                permanent.value,
            );
            permanent.underflow_safe &= prior_safe;
        }
        if r0[0] != 0.0 {
            let product_0 = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r1[1], r2[2]);
            let product_1 = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(r1[2], r2[1]);
            let mut pm00 = FilterArithmetic::<TRACK_UNDERFLOW>::add_non_negative(
                product_0.value.abs(),
                product_1.value.abs(),
            );
            pm00.underflow_safe &= product_0.underflow_safe && product_1.underflow_safe;
            let prior_safe = permanent.underflow_safe && pm00.underflow_safe;
            permanent = FilterArithmetic::<TRACK_UNDERFLOW>::mul_add(
                r0[0].abs(),
                pm00.value,
                permanent.value,
            );
            permanent.underflow_safe &= prior_safe;
        }
        permanent
    }

    /// Finish a determinant error bound only when its full arithmetic tree is
    /// outside the gradual-underflow regime.
    const fn certified_error_bound<const TRACK_UNDERFLOW: bool>(
        coefficient: f64,
        permanent: FilterArithmetic<TRACK_UNDERFLOW>,
    ) -> Result<Option<f64>, LaError> {
        let mut bound = FilterArithmetic::<TRACK_UNDERFLOW>::multiply(coefficient, permanent.value);
        bound.underflow_safe &= permanent.underflow_safe;
        if bound.underflow_safe {
            Self::computed_scalar_result(ArithmeticOperation::DeterminantErrorBound, bound.value)
        } else {
            cold_path();
            Ok(None)
        }
    }

    /// Return a computed scalar result for a matrix with finite stored entries.
    const fn computed_scalar_result(
        operation: ArithmeticOperation,
        value: f64,
    ) -> Result<Option<f64>, LaError> {
        if value.is_finite() {
            Ok(Some(value))
        } else {
            Err(LaError::non_finite_computation_scalar(operation))
        }
    }
}

impl<const D: usize> Default for Matrix<D> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(all(doc, feature = "exact"))]
mod det_errbound_doctests {
    /// ```rust
    /// use la_stack::prelude::*;
    ///
    /// fn adaptive_det_sign<const D: usize>(
    ///     matrix: &Matrix<D>,
    /// ) -> DeterminantSign {
    ///     if let Ok(Some(estimate)) = matrix.det_direct_with_errbound() {
    ///         if estimate.determinant().abs() > estimate.absolute_error_bound() {
    ///             return if estimate.determinant() > 0.0 {
    ///                 DeterminantSign::Positive
    ///             } else {
    ///                 DeterminantSign::Negative
    ///             };
    ///         }
    ///     }
    ///
    ///     matrix.det_sign_exact()
    /// }
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let identity = Matrix::<3>::identity();
    /// assert_eq!(
    ///     adaptive_det_sign(&identity),
    ///     DeterminantSign::Positive
    /// );
    ///
    /// let singular = Matrix::<3>::try_from_rows([
    ///     [1.0, 2.0, 3.0],
    ///     [4.0, 5.0, 6.0],
    ///     [7.0, 8.0, 9.0],
    /// ])?;
    /// assert_eq!(adaptive_det_sign(&singular), DeterminantSign::Zero);
    ///
    /// let big = f64::MAX / 2.0;
    /// let overflowing = Matrix::<3>::try_from_rows([
    ///     [0.0, 0.0, 1.0],
    ///     [big, 0.0, 1.0],
    ///     [0.0, big, 1.0],
    /// ])?;
    /// assert_eq!(
    ///     adaptive_det_sign(&overflowing),
    ///     DeterminantSign::Positive
    /// );
    /// # Ok(())
    /// # }
    /// ```
    fn adaptive_precision_pattern() {}
}

#[cfg(test)]
mod tests {
    use core::hint::black_box;

    use approx::assert_abs_diff_eq;
    use pastey::paste;

    use super::*;
    use crate::{DEFAULT_SINGULAR_TOL, FactorizationKind, Vector};

    macro_rules! gen_matrix_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<matrix_try_from_rows_get_set_bounds_checked_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = 1.0;
                    rows[$d - 1][$d - 1] = -2.0;

                    let mut m = Matrix::<$d>::try_from_rows(rows).unwrap();

                    assert_eq!(m.get(0, 0), Some(1.0));
                    assert_eq!(m.get($d - 1, $d - 1), Some(-2.0));
                    assert_eq!(m.try_get(0, 0), Ok(1.0));
                    assert_eq!(m.try_get($d - 1, $d - 1), Ok(-2.0));

                    // Out-of-bounds is None.
                    assert_eq!(m.get($d, 0), None);
                    assert_eq!(
                        m.try_get($d, 0),
                        Err(LaError::IndexOutOfBounds {
                            row: $d,
                            col: 0,
                            dim: $d,
                        })
                    );

                    // Out-of-bounds set fails.
                    let before_failed_set = m;
                    assert_eq!(
                        m.set($d, 0, 3.0),
                        Err(LaError::IndexOutOfBounds {
                            row: $d,
                            col: 0,
                            dim: $d,
                        })
                    );
                    assert_eq!(m, before_failed_set);
                    assert_eq!(
                        m.set(0, $d, 3.0),
                        Err(LaError::IndexOutOfBounds {
                            row: 0,
                            col: $d,
                            dim: $d,
                        })
                    );
                    assert_eq!(m, before_failed_set);
                    assert_eq!(m.get(0, 0), Some(1.0));

                    // In-bounds set works.
                    assert_eq!(m.set(0, $d - 1, 3.0), Ok(()));
                    assert_eq!(m.get(0, $d - 1), Some(3.0));
                    assert_eq!(m.set($d - 1, 0, 4.0), Ok(()));
                    assert_eq!(m.try_get($d - 1, 0), Ok(4.0));
                }

                #[test]
                fn [<matrix_set_rejects_non_finite_and_preserves_storage_ $d d>]() {
                    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                        let mut m = Matrix::<$d>::identity();
                        let before = m;
                        assert_eq!(
                            m.set($d - 1, 0, value),
                            Err(LaError::non_finite_input_matrix($d - 1, 0))
                        );
                        assert_eq!(m, before);
                    }
                }

                #[test]
                fn [<matrix_try_from_rows_rejects_non_finite_ $d d>]() {
                    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                        let mut rows = [[0.0f64; $d]; $d];
                        rows[$d - 1][$d - 1] = value;
                        assert_eq!(
                            Matrix::<$d>::try_from_rows(rows),
                            Err(LaError::non_finite_input_matrix($d - 1, $d - 1))
                        );
                    }

                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][$d - 1] = f64::INFINITY;
                    rows[$d - 1][0] = f64::NAN;
                    assert_eq!(
                        Matrix::<$d>::try_from_rows(rows),
                        Err(LaError::non_finite_input_matrix(0, $d - 1))
                    );
                }

                #[test]
                fn [<matrix_zero_and_default_are_zero_ $d d>]() {
                    let z = Matrix::<$d>::zero();
                    assert_abs_diff_eq!(z.inf_norm().unwrap(), 0.0, epsilon = 0.0);

                    let d = Matrix::<$d>::default();
                    assert_abs_diff_eq!(d.inf_norm().unwrap(), 0.0, epsilon = 0.0);
                }

                #[test]
                fn [<matrix_inf_norm_max_row_sum_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];

                    // Row 0 has a smaller absolute row sum.
                    for c in 0..$d {
                        rows[0][c] = 0.5;
                    }

                    // The last row has absolute row sum = D.
                    for c in 0..$d {
                        rows[$d - 1][c] = -1.0;
                    }

                    let m = Matrix::<$d>::try_from_rows(rows).unwrap();
                    assert_abs_diff_eq!(m.inf_norm().unwrap(), f64::from($d), epsilon = 0.0);
                }

                #[test]
                fn [<matrix_inf_norm_reports_first_overflowing_column_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[$d - 1][0] = f64::MAX;
                    rows[$d - 1][1] = f64::MAX;

                    let m = Matrix::<$d>::try_from_rows(rows).unwrap();
                    assert_eq!(
                        m.inf_norm(),
                        Err(LaError::non_finite_computation_matrix(
                            ArithmeticOperation::MatrixInfinityNorm,
                            $d - 1,
                            1,
                        ))
                    );
                }

                #[test]
                fn [<matrix_inf_norm_reports_first_overflowing_row_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = f64::MAX;
                    rows[0][$d - 1] = f64::MAX;
                    rows[$d - 1][0] = f64::MAX;
                    rows[$d - 1][1] = f64::MAX;

                    let m = Matrix::<$d>::try_from_rows(rows).unwrap();
                    assert_eq!(
                        m.inf_norm(),
                        Err(LaError::non_finite_computation_matrix(
                            ArithmeticOperation::MatrixInfinityNorm,
                            0,
                            $d - 1,
                        ))
                    );
                }

                #[test]
                fn [<matrix_identity_lu_det_solve_ $d d>]() {
                    let m = Matrix::<$d>::identity();

                    // Identity has ones on diag and zeros off diag.
                    for r in 0..$d {
                        for c in 0..$d {
                            let expected = if r == c { 1.0 } else { 0.0 };
                            assert_abs_diff_eq!(m.get(r, c).unwrap(), expected, epsilon = 0.0);
                        }
                    }

                    // Determinant is 1.
                    let det = m.det().unwrap();
                    assert_abs_diff_eq!(det, 1.0, epsilon = 1e-12);

                    // LU solve on identity returns the RHS.
                    let lu = m.lu(DEFAULT_SINGULAR_TOL).unwrap();

                    let b_arr = {
                        let mut arr = [0.0f64; $d];
                        let values = [1.0f64, 2.0, 3.0, 4.0, 5.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = *src;
                        }
                        arr
                    };

                    let b = Vector::<$d>::new(b_arr);
                    let x = lu.solve(b).unwrap().into_array();

                    for (x_i, b_i) in x.iter().zip(b_arr.iter()) {
                        assert_abs_diff_eq!(*x_i, *b_i, epsilon = 1e-12);
                    }
                }

            }
        };
    }

    // Mirror delaunay-style multi-dimension tests.
    gen_matrix_tests!(2);
    gen_matrix_tests!(3);
    gen_matrix_tests!(4);
    gen_matrix_tests!(5);

    #[test]
    fn matrix_inf_norm_preserves_left_to_right_row_sum_order() {
        let large = 9_007_199_254_740_992.0;
        let matrix =
            Matrix::<4>::try_from_rows([[large, 1.0, 1.0, 1.0], [0.0; 4], [0.0; 4], [0.0; 4]])
                .unwrap();

        assert_eq!(matrix.inf_norm(), Ok(large));
    }

    // === det_direct tests ===

    #[test]
    fn det_direct_d0_is_one() {
        assert_eq!(Matrix::<0>::zero().det_direct(), Ok(Some(1.0)));
    }

    #[test]
    fn det_direct_d1_returns_element() {
        let m = Matrix::<1>::try_from_rows([[42.0]]).unwrap();
        assert_eq!(m.det_direct(), Ok(Some(42.0)));
    }

    #[test]
    fn det_direct_d2_known_value() {
        // [[1,2],[3,4]] → det = 1*4 - 2*3 = -2
        // black_box prevents compile-time constant folding of the const fn.
        let m = black_box(Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap());
        assert_abs_diff_eq!(m.det_direct().unwrap().unwrap(), -2.0, epsilon = 1e-15);
    }

    #[test]
    fn det_direct_d3_known_value() {
        // Classic 3×3: det = 0
        let m = black_box(
            Matrix::<3>::try_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
                .unwrap(),
        );
        assert_abs_diff_eq!(m.det_direct().unwrap().unwrap(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn det_direct_d3_dense_known_value() {
        // det = 1*(5*8 - 7*6) - 2*(4*8 - 7*2) + 3*(4*6 - 5*2) = 4
        let m = black_box(
            Matrix::<3>::try_from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 7.0], [2.0, 6.0, 8.0]])
                .unwrap(),
        );
        let direct = m.det_direct().unwrap().unwrap();
        let paired = m.det_direct_with_errbound().unwrap().unwrap();

        assert_abs_diff_eq!(direct, 4.0, epsilon = 1e-12);
        assert_eq!(paired.determinant().to_bits(), direct.to_bits());
    }

    #[test]
    fn det_direct_d3_dense_reports_legitimate_overflow() {
        // The unscaled matrix has determinant 54, so scaling every entry by
        // 1.6e102 gives a determinant of approximately 2.21e308.
        let scale = 1.6e102;
        let m = black_box(
            Matrix::<3>::try_from_rows([
                [4.0 * scale, scale, scale],
                [scale, 4.0 * scale, scale],
                [scale, scale, 4.0 * scale],
            ])
            .unwrap(),
        );
        let expected = LaError::non_finite_computation_scalar(ArithmeticOperation::Determinant);

        assert_eq!(m.det_direct(), Err(expected));
        assert_eq!(m.det(), Err(expected));
    }

    #[test]
    fn det_errbound_d3_dense_reports_legitimate_overflow() {
        // The unscaled matrix has determinant 54, so scaling every entry by
        // 1.6e102 gives a determinant of approximately 2.21e308.
        let scale = 1.6e102;
        let m = black_box(
            Matrix::<3>::try_from_rows([
                [4.0 * scale, scale, scale],
                [scale, 4.0 * scale, scale],
                [scale, scale, 4.0 * scale],
            ])
            .unwrap(),
        );
        let expected =
            LaError::non_finite_computation_scalar(ArithmeticOperation::DeterminantErrorBound);

        assert_eq!(m.det_errbound(), Err(expected));
        assert_eq!(m.det_direct_with_errbound(), Err(expected));
    }

    #[test]
    fn det_direct_d3_nonsingular() {
        // [[2,1,0],[0,3,1],[1,0,2]] → det = 2*(6-0) - 1*(0-1) + 0 = 13
        let m = black_box(
            Matrix::<3>::try_from_rows([[2.0, 1.0, 0.0], [0.0, 3.0, 1.0], [1.0, 0.0, 2.0]])
                .unwrap(),
        );
        assert_abs_diff_eq!(m.det_direct().unwrap().unwrap(), 13.0, epsilon = 1e-12);
    }

    #[test]
    fn det_direct_d3_skips_zero_coefficient_minor_that_would_overflow() {
        let m = black_box(
            Matrix::<3>::try_from_rows([
                [1.0, 0.0, 0.0],
                [1.0e300, 1.0, 1.0e300],
                [1.0e300, 0.0, 1.0e300],
            ])
            .unwrap(),
        );
        assert_eq!(m.det_direct(), Ok(Some(1.0e300)));
    }

    #[test]
    fn det_direct_d4_known_value() {
        // Diagonal matrix: det = product of diagonal entries.
        let mut rows = [[0.0f64; 4]; 4];
        rows[0][0] = 2.0;
        rows[1][1] = 3.0;
        rows[2][2] = 5.0;
        rows[3][3] = 7.0;
        let m = black_box(Matrix::<4>::try_from_rows(rows).unwrap());
        assert_abs_diff_eq!(m.det_direct().unwrap().unwrap(), 210.0, epsilon = 1e-12);
    }

    #[test]
    fn det_direct_d4_dense_known_value() {
        let m = black_box(
            Matrix::<4>::try_from_rows([
                [4.0, 1.0, 3.0, 2.0],
                [1.0, 5.0, 2.0, 1.0],
                [7.0, 2.0, 6.0, 3.0],
                [1.0, 8.0, 4.0, 9.0],
            ])
            .unwrap(),
        );
        let direct = m.det_direct().unwrap().unwrap();
        let paired = m.det_direct_with_errbound().unwrap().unwrap();

        assert_abs_diff_eq!(direct, 112.0, epsilon = 1e-12);
        assert_eq!(paired.determinant().to_bits(), direct.to_bits());
    }

    #[test]
    fn det_direct_d4_dense_reports_legitimate_overflow() {
        // The unscaled matrix has determinant 189, so scaling every entry by
        // 3.2e76 gives a determinant of approximately 1.98e308.
        let scale = 3.2e76;
        let m = black_box(
            Matrix::<4>::try_from_rows([
                [4.0 * scale, scale, scale, scale],
                [scale, 4.0 * scale, scale, scale],
                [scale, scale, 4.0 * scale, scale],
                [scale, scale, scale, 4.0 * scale],
            ])
            .unwrap(),
        );
        let expected = LaError::non_finite_computation_scalar(ArithmeticOperation::Determinant);

        assert_eq!(m.det_direct(), Err(expected));
        assert_eq!(m.det(), Err(expected));
    }

    #[test]
    fn det_errbound_d4_dense_reports_legitimate_overflow() {
        // The unscaled matrix has determinant 189, so scaling every entry by
        // 3.2e76 gives a determinant of approximately 1.98e308.
        let scale = 3.2e76;
        let m = black_box(
            Matrix::<4>::try_from_rows([
                [4.0 * scale, scale, scale, scale],
                [scale, 4.0 * scale, scale, scale],
                [scale, scale, 4.0 * scale, scale],
                [scale, scale, scale, 4.0 * scale],
            ])
            .unwrap(),
        );
        let expected =
            LaError::non_finite_computation_scalar(ArithmeticOperation::DeterminantErrorBound);

        assert_eq!(m.det_errbound(), Err(expected));
        assert_eq!(m.det_direct_with_errbound(), Err(expected));
    }

    #[test]
    fn det_direct_d4_skips_zero_coefficient_cofactors_that_would_overflow() {
        let m = black_box(
            Matrix::<4>::try_from_rows([
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [1.0e300, 0.0, 1.0e300, 1.0e300],
                [1.0e300, 0.0, 1.0e300, -1.0e300],
            ])
            .unwrap(),
        );
        assert_eq!(m.det_direct(), Ok(Some(0.0)));
    }

    #[test]
    fn det_direct_d4_sparse_second_row_skips_inactive_overflowing_minors() {
        let m = black_box(
            Matrix::<4>::try_from_rows([
                [1.0e-300, 1.0, 1.0, 1.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0e300, 1.0, 1.0e300],
                [0.0, 1.0e300, 0.0, 1.0e300],
            ])
            .unwrap(),
        );

        assert_eq!(m.det_direct(), Ok(Some(1.0)));
        assert_eq!(m.det(), Ok(1.0));
    }

    #[test]
    fn det_direct_d5_returns_none() {
        assert_eq!(Matrix::<5>::identity().det_direct(), Ok(None));
    }

    #[test]
    fn det_direct_d8_returns_none() {
        assert_eq!(Matrix::<8>::zero().det_direct(), Ok(None));
    }

    #[test]
    fn det_direct_rejects_computed_overflow() {
        let m = Matrix::<2>::try_from_rows([[1e300, 0.0], [0.0, 1e300]]).unwrap();
        assert_eq!(
            m.det_direct(),
            Err(LaError::non_finite_computation_scalar(
                ArithmeticOperation::Determinant
            ))
        );
    }

    #[test]
    fn det_d5_rejects_lu_product_overflow() {
        let m = Matrix::<5>::try_from_rows([
            [1.0e100, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0e100, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0e100, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0e100, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0e100],
        ])
        .unwrap();
        assert_eq!(
            m.det(),
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::Determinant,
                4
            ))
        );
    }

    #[test]
    fn det_d5_rejects_lu_trailing_update_overflow() {
        let m = Matrix::<5>::try_from_rows([
            [1.0, f64::MAX, 0.0, 0.0, 0.0],
            [-1.0, f64::MAX, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ])
        .unwrap();

        assert_eq!(
            m.det(),
            Err(LaError::non_finite_computation_matrix(
                ArithmeticOperation::LuFactorization,
                1,
                1
            ))
        );
    }

    macro_rules! gen_det_direct_agrees_with_lu {
        ($d:literal) => {
            paste! {
                #[test]
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "r, c, and D are tiny test integers exactly representable as f64"
                )]
                fn [<det_direct_agrees_with_lu_ $d d>]() {
                    // Well-conditioned matrix: diagonally dominant.
                    let mut rows = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            rows[r][c] = if r == c {
                                (r as f64) + f64::from($d) + 1.0
                            } else {
                                0.1 / ((r + c + 1) as f64)
                            };
                        }
                    }
                    let m = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let direct = m.det_direct().unwrap().unwrap();
                    let lu_det = m.lu(DEFAULT_SINGULAR_TOL).unwrap().det().unwrap();
                    let eps = lu_det.abs().mul_add(1e-12, 1e-12);
                    assert_abs_diff_eq!(direct, lu_det, epsilon = eps);
                }
            }
        };
    }

    gen_det_direct_agrees_with_lu!(1);
    gen_det_direct_agrees_with_lu!(2);
    gen_det_direct_agrees_with_lu!(3);
    gen_det_direct_agrees_with_lu!(4);

    #[test]
    fn det_direct_identity_all_dims() {
        assert_abs_diff_eq!(
            Matrix::<1>::identity().det_direct().unwrap().unwrap(),
            1.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<2>::identity().det_direct().unwrap().unwrap(),
            1.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<3>::identity().det_direct().unwrap().unwrap(),
            1.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<4>::identity().det_direct().unwrap().unwrap(),
            1.0,
            epsilon = 0.0
        );
    }

    #[test]
    fn det_direct_zero_matrix() {
        assert_abs_diff_eq!(
            Matrix::<2>::zero().det_direct().unwrap().unwrap(),
            0.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<3>::zero().det_direct().unwrap().unwrap(),
            0.0,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<4>::zero().det_direct().unwrap().unwrap(),
            0.0,
            epsilon = 0.0
        );
    }

    macro_rules! gen_det_singular_zero_matrix_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<det_singular_zero_matrix_returns_zero_ $d d>]() {
                    assert_abs_diff_eq!(
                        Matrix::<$d>::zero().det().unwrap(),
                        0.0,
                        epsilon = 0.0
                    );
                }
            }
        };
    }

    gen_det_singular_zero_matrix_tests!(2);
    gen_det_singular_zero_matrix_tests!(3);
    gen_det_singular_zero_matrix_tests!(4);

    #[test]
    fn det_singular_zero_matrix_d5_preserves_lu_error() {
        assert_eq!(
            Matrix::<5>::zero().det(),
            Err(LaError::singular_numerical(
                0,
                FactorizationKind::Lu,
                0.0,
                0.0
            ))
        );
    }

    #[test]
    fn det_d5_does_not_turn_elimination_underflow_into_exact_zero() {
        let min_subnormal = f64::from_bits(1);
        let two_pow_800 = f64::from_bits(1823_u64 << 52);
        let m = Matrix::<5>::try_from_rows([
            [2.0, min_subnormal, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, two_pow_800, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ])
        .unwrap();

        assert_eq!(
            m.det(),
            Err(LaError::singular_numerical(
                1,
                FactorizationKind::Lu,
                0.0,
                0.0
            ))
        );
    }

    #[test]
    fn det_d5_ignores_pivot_tolerance_for_tiny_nonsingular_matrix() {
        // A small nonzero determinant is still a determinant. `det` must not
        // flatten the value to zero merely because the default LU tolerance
        // would reject a pivot this small.
        let m = Matrix::<5>::try_from_rows([
            [1e-13, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ])
        .unwrap();

        assert_abs_diff_eq!(m.det().unwrap(), 1e-13, epsilon = 0.0);
        assert_eq!(
            m.lu(DEFAULT_SINGULAR_TOL),
            Err(LaError::singular_numerical(
                0,
                FactorizationKind::Lu,
                1e-13,
                DEFAULT_SINGULAR_TOL.get()
            ))
        );
    }

    #[test]
    fn det_returns_non_finite_error_for_overflow_with_finite_entries() {
        // det_direct produces an overflowing f64 (1e300 * 1e300 = ∞) even
        // though every matrix entry is finite. The entry scan in `det`
        // falls through and reports a computed determinant overflow rather
        // than a NaN/∞ input.
        let m = Matrix::<2>::try_from_rows([[1e300, 0.0], [0.0, 1e300]]).unwrap();
        assert_eq!(
            m.det(),
            Err(LaError::non_finite_computation_scalar(
                ArithmeticOperation::Determinant
            ))
        );
    }

    // === det_direct const-evaluability tests (D = 2..=5) ===
    //
    // Every dimension hits a distinct arm of the `match D { … }` body inside
    // `det_direct`, so exercising each at compile time is the tightest
    // const-fn proof available.

    macro_rules! gen_det_direct_const_eval_tests {
        ($d:literal) => {
            paste! {
                /// `Matrix::<D>::det_direct()` on the identity must const-evaluate
                /// to `Ok(Some(1.0))` for every closed-form dimension `D ∈ {1, 2, 3, 4}`.
                #[test]
                fn [<det_direct_const_eval_ $d d>]() {
                    const DET: Result<Option<f64>, LaError> = Matrix::<$d>::identity().det_direct();
                    assert_eq!(DET, Ok(Some(1.0)));
                }
            }
        };
    }

    gen_det_direct_const_eval_tests!(2);
    gen_det_direct_const_eval_tests!(3);
    gen_det_direct_const_eval_tests!(4);

    #[test]
    fn det_direct_const_eval_d5_is_none() {
        // D ≥ 5 has no closed-form arm; `det_direct` returns `Ok(None)`.  Verify
        // that the wildcard arm is reachable in a `const { … }` context.
        const DET: Result<Option<f64>, LaError> = Matrix::<5>::identity().det_direct();
        assert_eq!(DET, Ok(None));
    }

    // === det_errbound tests (no `exact` feature required) ===

    #[test]
    fn det_errbound_matches_documented_coefficient_scale() {
        let m2 = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        let expected_2 = ERR_COEFF_2 * ((1.0_f64 * 4.0).abs() + (2.0_f64 * 3.0).abs());
        assert_abs_diff_eq!(
            m2.det_errbound().unwrap().unwrap(),
            expected_2,
            epsilon = 0.0
        );

        assert_abs_diff_eq!(
            Matrix::<3>::identity().det_errbound().unwrap().unwrap(),
            ERR_COEFF_3,
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            Matrix::<4>::identity().det_errbound().unwrap().unwrap(),
            ERR_COEFF_4,
            epsilon = 0.0
        );
    }

    #[test]
    fn det_errbound_d3_skips_zero_coefficient_minor_that_would_overflow() {
        let m = Matrix::<3>::try_from_rows([
            [1.0, 0.0, 0.0],
            [1.0e300, 1.0, 1.0e300],
            [1.0e300, 0.0, 1.0e300],
        ])
        .unwrap();

        assert_eq!(m.det_errbound(), Ok(Some(ERR_COEFF_3 * 1.0e300)));
    }

    #[test]
    fn det_errbound_d4_skips_zero_coefficient_cofactors_that_would_overflow() {
        let m = Matrix::<4>::try_from_rows([
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [1.0e300, 0.0, 1.0e300, 1.0e300],
            [1.0e300, 0.0, 1.0e300, -1.0e300],
        ])
        .unwrap();

        assert_eq!(m.det_errbound(), Ok(Some(0.0)));
    }

    #[test]
    fn det_errbound_d5_returns_none() {
        // D=5 has no fast filter
        assert_eq!(Matrix::<5>::identity().det_errbound(), Ok(None));
    }

    #[test]
    fn combined_det_bound_wide_exponent_fast_path_matches_tracked_arithmetic() {
        let threshold = f64::from_bits(1007_u64 << 52); // 2^-16
        let at_threshold = Matrix::<2>::try_from_rows([[threshold, 0.0], [0.0, 2.0]]).unwrap();
        assert!(at_threshold.det_bound_inputs_have_wide_exponent_margin());

        let tracked = at_threshold
            .det_direct_with_errbound_from_arithmetic(
                at_threshold
                    .det_direct_arithmetic::<true>()
                    .expect("D=2 has direct arithmetic"),
            )
            .unwrap();
        assert_eq!(at_threshold.det_direct_with_errbound().unwrap(), tracked);

        let just_below = f64::from_bits(threshold.to_bits() - 1);
        let below_threshold = Matrix::<2>::try_from_rows([[just_below, 0.0], [0.0, 2.0]]).unwrap();
        assert!(!below_threshold.det_bound_inputs_have_wide_exponent_margin());
        assert!(!Matrix::<5>::identity().det_bound_inputs_have_wide_exponent_margin());
    }

    #[test]
    fn det_direct_with_errbound_covers_zero_and_one_dimensions() {
        let empty = Matrix::<0>::zero()
            .det_direct_with_errbound()
            .unwrap()
            .unwrap();
        assert_abs_diff_eq!(empty.determinant(), 1.0, epsilon = 0.0);
        assert_abs_diff_eq!(empty.absolute_error_bound(), 0.0, epsilon = 0.0);

        let scalar = Matrix::<1>::try_from_rows([[-7.0]])
            .unwrap()
            .det_direct_with_errbound()
            .unwrap()
            .unwrap();
        assert_abs_diff_eq!(scalar.determinant(), -7.0, epsilon = 0.0);
        assert_abs_diff_eq!(scalar.absolute_error_bound(), 0.0, epsilon = 0.0);
    }

    #[test]
    fn det_direct_with_errbound_pairs_the_closed_form_values() {
        let matrix = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]]).unwrap();
        let estimate = matrix.det_direct_with_errbound().unwrap().unwrap();

        assert_abs_diff_eq!(
            estimate.determinant(),
            matrix.det_direct().unwrap().unwrap(),
            epsilon = 0.0
        );
        assert_abs_diff_eq!(
            estimate.absolute_error_bound(),
            ERR_COEFF_2 * (4.0_f64 + 6.0_f64),
            epsilon = 0.0
        );
    }

    #[test]
    fn det_direct_with_errbound_d5_returns_none() {
        assert_eq!(Matrix::<5>::identity().det_direct_with_errbound(), Ok(None));
    }

    #[test]
    fn det_errbound_rejects_computed_overflow() {
        let m = Matrix::<2>::try_from_rows([[1e300, 0.0], [0.0, 1e300]]).unwrap();
        assert_eq!(
            m.det_errbound(),
            Err(LaError::non_finite_computation_scalar(
                ArithmeticOperation::DeterminantErrorBound
            ))
        );
    }

    // === det_errbound const-evaluability tests (D = 2..=5) ===

    macro_rules! gen_det_errbound_const_eval_tests {
        ($d:literal) => {
            paste! {
                /// `Matrix::<D>::det_errbound()` on the identity must const-evaluate
                /// to `Ok(Some(bound))` with `bound > 0` for every closed-form dimension
                /// `D ∈ {2, 3, 4}`.  Each dimension hits a distinct arm of
                /// `det_errbound` with a dimension-specific permanent computation.
                #[test]
                fn [<det_errbound_const_eval_ $d d>]() {
                    const BOUND: Result<Option<f64>, LaError> = Matrix::<$d>::identity().det_errbound();
                    assert!(BOUND.unwrap().unwrap() > 0.0);
                }
            }
        };
    }

    gen_det_errbound_const_eval_tests!(2);
    gen_det_errbound_const_eval_tests!(3);
    gen_det_errbound_const_eval_tests!(4);

    #[test]
    fn det_errbound_const_eval_d5_is_none() {
        // D ≥ 5 has no fast-filter bound; `det_errbound` returns `Ok(None)`.
        const BOUND: Result<Option<f64>, LaError> = Matrix::<5>::identity().det_errbound();
        assert_eq!(BOUND, Ok(None));
    }

    // === inf_norm const-evaluability tests (D = 2..=5) ===

    macro_rules! gen_inf_norm_const_eval_tests {
        ($d:literal) => {
            paste! {
                /// `Matrix::<D>::inf_norm()` on the identity must const-evaluate
                /// to `1.0` for every `D ≥ 1` — each row has a single `1.0`
                /// entry, so the max absolute row sum is exactly `1.0`.
                #[test]
                fn [<inf_norm_const_eval_ $d d>]() {
                    const NORM: Result<f64, LaError> = Matrix::<$d>::identity().inf_norm();
                    assert!((NORM.unwrap() - 1.0).abs() <= 1e-12);
                }
            }
        };
    }

    gen_inf_norm_const_eval_tests!(2);
    gen_inf_norm_const_eval_tests!(3);
    gen_inf_norm_const_eval_tests!(4);
    gen_inf_norm_const_eval_tests!(5);

    // === is_symmetric / first_asymmetry (public LDLT preconditions helpers) ===

    macro_rules! gen_is_symmetric_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<is_symmetric_true_for_identity_ $d d>]() {
                    let m = Matrix::<$d>::identity();
                    assert!(m.is_symmetric(Tolerance::try_new(1e-12).unwrap()).unwrap());
                    assert_eq!(m.first_asymmetry(Tolerance::try_new(1e-12).unwrap()).unwrap(), None);
                }

                #[test]
                fn [<is_symmetric_true_for_zero_ $d d>]() {
                    let m = Matrix::<$d>::zero();
                    assert!(m.is_symmetric(Tolerance::try_new(1e-12).unwrap()).unwrap());
                    assert_eq!(m.first_asymmetry(Tolerance::try_new(1e-12).unwrap()).unwrap(), None);
                }

                #[test]
                fn [<is_symmetric_true_for_constructed_symmetric_ $d d>]() {
                    // Construct A = M + Mᵀ so A is provably symmetric.
                    let mut m = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            #[expect(
                                clippy::cast_precision_loss,
                                reason = "matrix test indices are at most five and exactly representable as f64"
                            )]
                            {
                                m[r][c] = (r * $d + c) as f64;
                            }
                        }
                    }
                    let mut sym = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            sym[r][c] = m[r][c] + m[c][r];
                        }
                    }
                    let a = Matrix::<$d>::try_from_rows(sym).unwrap();
                    assert!(a.is_symmetric(Tolerance::try_new(1e-12).unwrap()).unwrap());
                    assert_eq!(a.first_asymmetry(Tolerance::try_new(1e-12).unwrap()).unwrap(), None);
                }

                #[test]
                fn [<is_symmetric_false_for_asymmetric_offdiagonal_ $d d>]() {
                    // Perturb a single off-diagonal entry so symmetry fails.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 1.0;
                    }
                    rows[0][$d - 1] = 1.0;
                    rows[$d - 1][0] = -1.0; // breaks symmetry
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();
                    assert!(!a.is_symmetric(Tolerance::try_new(1e-12).unwrap()).unwrap());
                    assert_eq!(
                        a.first_asymmetry(Tolerance::try_new(1e-12).unwrap()).unwrap(),
                        Some((0, $d - 1))
                    );
                }

            }
        };
    }

    gen_is_symmetric_tests!(2);
    gen_is_symmetric_tests!(3);
    gen_is_symmetric_tests!(4);
    gen_is_symmetric_tests!(5);

    macro_rules! gen_ldlt_symmetry_proof_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<matrix_ldlt_accepts_exact_symmetric_spd_ $d d>]() {
                    // This exactly mirrored, strictly diagonally dominant
                    // tridiagonal matrix is positive definite.
                    let mut rows = [[0.0_f64; $d]; $d];
                    for (index, row) in rows.iter_mut().enumerate() {
                        row[index] = 2.0;
                    }
                    for index in 1..$d {
                        rows[index - 1][index] = 0.5;
                        rows[index][index - 1] = 0.5;
                    }

                    let matrix = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let symmetric = SymmetricMatrix::try_new(matrix).unwrap();
                    assert_eq!(symmetric.into_matrix(), matrix);

                    let ldlt = matrix.ldlt(DEFAULT_SINGULAR_TOL).unwrap();
                    assert!(ldlt.det().unwrap() > 0.0);
                }

                #[test]
                fn [<symmetric_matrix_try_new_rejects_finite_asymmetric_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    for (i, row) in rows.iter_mut().enumerate() {
                        row[i] = 1.0;
                    }
                    rows[0][$d - 1] = 1.0;
                    rows[$d - 1][0] = -1.0;

                    assert_eq!(
                        Matrix::<$d>::try_from_rows(rows).and_then(SymmetricMatrix::try_new),
                        Err(LaError::asymmetric(0, $d - 1, $d, 1.0, -1.0, 0.0))
                    );
                }
            }
        };
    }

    gen_ldlt_symmetry_proof_tests!(2);
    gen_ldlt_symmetry_proof_tests!(3);
    gen_ldlt_symmetry_proof_tests!(4);
    gen_ldlt_symmetry_proof_tests!(5);

    #[test]
    fn symmetric_matrix_into_matrix_roundtrips_storage_internally() {
        let a = Matrix::<2>::try_from_rows([[2.0, 1.0], [1.0, 3.0]]).unwrap();
        let symmetric = SymmetricMatrix::try_new(a).unwrap();

        assert_eq!(symmetric.into_matrix(), a);
    }

    #[test]
    fn matrix_ldlt_accepts_opposite_signed_zero_mirrors() {
        let matrix = Matrix::<2>::try_from_rows([[2.0, 0.0], [-0.0, 2.0]]).unwrap();
        let ldlt = matrix.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        assert_eq!(ldlt.det(), Ok(4.0));
    }

    #[test]
    fn is_symmetric_tolerance_scales_with_inf_norm() {
        // Off-diagonal entries differ by 1e-6.  With inf_norm ≈ 2e6, the
        // relative tolerance 1e-12 yields eps ≈ 2e-6, which accepts the gap;
        // a stricter tol of 1e-15 rejects it.
        let a = Matrix::<2>::try_from_rows([[1.0e6, 1.0e6 + 1.0e-6], [1.0e6, 1.0e6]]).unwrap();
        assert!(a.is_symmetric(Tolerance::try_new(1e-12).unwrap()).unwrap());
        assert!(!a.is_symmetric(Tolerance::try_new(1e-15).unwrap()).unwrap());
    }

    #[test]
    fn symmetry_epsilon_multiplies_after_row_sum_near_subnormal_boundary() {
        let min_subnormal = f64::from_bits(1);
        let mut rows = [[0.0; 5]; 5];
        let mut col = 0;
        while col < 4 {
            rows[0][col] = 0.4;
            rows[col][0] = 0.4;
            col += 1;
        }
        rows[0][4] = 2.0 * min_subnormal;
        rows[4][0] = 0.0;

        let matrix = Matrix::<5>::try_from_rows(rows).unwrap();
        let tolerance = Tolerance::try_new(min_subnormal).unwrap();
        let expected_epsilon = tolerance.get() * matrix.inf_norm().unwrap().max(1.0);

        assert_eq!(expected_epsilon.to_bits(), 2);
        assert_eq!(matrix.first_asymmetry(tolerance), Ok(None));
        assert_eq!(matrix.is_symmetric(tolerance), Ok(true));
    }

    #[test]
    fn symmetry_epsilon_scales_terms_when_row_sum_overflows() {
        let matrix =
            Matrix::<2>::try_from_rows([[f64::MAX, f64::MAX], [f64::MAX / 2.0, f64::MAX]]).unwrap();

        assert_eq!(
            matrix.inf_norm(),
            Err(LaError::non_finite_computation_matrix(
                ArithmeticOperation::MatrixInfinityNorm,
                0,
                1
            ))
        );
        assert_eq!(
            matrix.first_asymmetry(Tolerance::try_new(0.25).unwrap()),
            Ok(None)
        );
        assert_eq!(
            matrix.first_asymmetry(Tolerance::try_new(0.125).unwrap()),
            Ok(Some((0, 1)))
        );
    }

    #[test]
    fn first_asymmetry_returns_lexicographically_first_pair() {
        // Two asymmetric pairs: (0, 2) and (1, 2).  We must get (0, 2) first.
        let a = Matrix::<3>::try_from_rows([[1.0, 0.0, 2.0], [0.0, 1.0, 3.0], [-2.0, -3.0, 1.0]])
            .unwrap();
        assert_eq!(
            a.first_asymmetry(Tolerance::try_new(1e-12).unwrap())
                .unwrap(),
            Some((0, 2))
        );
    }

    #[test]
    fn first_asymmetry_strict_tol_survives_row_sum_overflow() {
        let a = Matrix::<3>::try_from_rows([
            [1.0, 1.0, 0.0],
            [2.0, f64::MAX, f64::MAX],
            [0.0, 0.0, 1.0],
        ])
        .unwrap();

        assert_eq!(
            a.inf_norm(),
            Err(LaError::non_finite_computation_matrix(
                ArithmeticOperation::MatrixInfinityNorm,
                1,
                2
            ))
        );
        assert_eq!(
            a.first_asymmetry(Tolerance::try_new(0.0).unwrap()).unwrap(),
            Some((0, 1))
        );
        assert!(!a.is_symmetric(Tolerance::try_new(0.0).unwrap()).unwrap());
    }

    #[test]
    fn first_asymmetry_rejects_scaled_epsilon_overflow() {
        let a = Matrix::<2>::try_from_rows([[0.0, 0.0], [2.0, 1.0]]).unwrap();
        let tol = Tolerance::try_new(f64::MAX).unwrap();

        assert_eq!(
            a.first_asymmetry(tol),
            Err(LaError::non_finite_computation_matrix(
                ArithmeticOperation::SymmetryCheck,
                1,
                0
            ))
        );
        assert_eq!(
            a.is_symmetric(tol),
            Err(LaError::non_finite_computation_matrix(
                ArithmeticOperation::SymmetryCheck,
                1,
                0
            ))
        );
    }

    #[test]
    fn first_asymmetry_flags_overflowed_finite_difference() {
        let a = Matrix::<2>::try_from_rows([[1.0, f64::MAX], [-f64::MAX, 1.0]]).unwrap();
        assert_eq!(
            a.first_asymmetry(Tolerance::try_new(1e-12).unwrap())
                .unwrap(),
            Some((0, 1))
        );
        assert!(!a.is_symmetric(Tolerance::try_new(1e-12).unwrap()).unwrap());
    }
}
