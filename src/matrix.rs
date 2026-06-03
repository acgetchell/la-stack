//! Fixed-size, stack-allocated square matrices.

use core::hint::cold_path;

use crate::ldlt::Ldlt;
use crate::lu::Lu;
use crate::{ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4, LDLT_SYMMETRY_REL_TOL, LaError, Tolerance};

/// Fixed-size square matrix `D×D`, stored inline.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Matrix<const D: usize> {
    pub(crate) rows: [[f64; D]; D],
}

/// Fixed-size square matrix whose stored entries are all finite.
///
/// This proof-bearing wrapper makes the finite invariant explicit for internal
/// algorithms that should not repeatedly check stored entries for NaN or
/// infinity.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
#[allow(clippy::redundant_pub_crate)]
pub(crate) struct FiniteMatrix<const D: usize> {
    matrix: Matrix<D>,
}

impl<const D: usize> FiniteMatrix<D> {
    /// Construct a finite matrix without checking the invariant.
    ///
    /// This crate-internal escape hatch is only for paths with a local proof
    /// that every stored entry is finite.
    #[inline]
    pub(crate) const fn new_unchecked(matrix: Matrix<D>) -> Self {
        Self { matrix }
    }

    /// Validate a matrix and wrap it for algorithms that carry the finite
    /// invariant explicitly.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] with matrix coordinates for the first
    /// stored NaN or infinity in row-major order.
    #[inline]
    pub(crate) const fn new(matrix: Matrix<D>) -> Result<Self, LaError> {
        if let Some((row, col)) = Matrix::<D>::first_non_finite_cell_in(&matrix.rows) {
            Err(LaError::non_finite_cell(row, col))
        } else {
            Ok(Self::new_unchecked(matrix))
        }
    }

    /// Validate raw row-major storage and construct a finite matrix.
    /// # Errors
    /// Returns [`LaError::NonFinite`] with matrix coordinates for the first
    /// offending entry in row-major order when `rows` contains NaN or infinity.
    #[inline]
    pub(crate) const fn from_rows(rows: [[f64; D]; D]) -> Result<Self, LaError> {
        match Matrix::try_from_rows(rows) {
            Ok(matrix) => Ok(Self::new_unchecked(matrix)),
            Err(err) => Err(err),
        }
    }

    /// All-zeros finite matrix.
    #[inline]
    pub(crate) const fn zero() -> Self {
        Self::new_unchecked(Matrix::zero())
    }

    /// Consume the wrapper and return the underlying raw matrix.
    #[inline]
    pub(crate) const fn into_matrix(self) -> Matrix<D> {
        self.matrix
    }

    /// Borrow the underlying raw matrix without revalidating stored entries.
    #[cfg(feature = "exact")]
    #[inline]
    pub(crate) const fn as_matrix(&self) -> &Matrix<D> {
        &self.matrix
    }

    /// Infinity norm (maximum absolute row sum) for a finite matrix.
    ///
    /// Stored entries are known finite, so this path only checks whether a row
    /// sum overflows to NaN or infinity.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when a row sum overflows to NaN or infinity.
    #[inline]
    pub(crate) const fn inf_norm(&self) -> Result<f64, LaError> {
        let mut max_row_sum: f64 = 0.0;

        let mut r = 0;
        while r < D {
            let row = &self.matrix.rows[r];
            let mut row_sum: f64 = 0.0;
            let mut c = 0;
            while c < D {
                row_sum += row[c].abs();
                c += 1;
            }
            if !row_sum.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(r));
            }
            if row_sum > max_row_sum {
                max_row_sum = row_sum;
            }
            r += 1;
        }

        Ok(max_row_sum)
    }

    /// Returns `true` if the finite matrix is symmetric within a relative tolerance.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when computing the scaled symmetry tolerance
    /// overflows to NaN or infinity.
    #[inline]
    pub(crate) fn is_symmetric(&self, rel_tol: Tolerance) -> Result<bool, LaError> {
        Ok(self.first_asymmetry(rel_tol)?.is_none())
    }

    /// Returns the first asymmetric off-diagonal pair, if any.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when computing the scaled symmetry tolerance
    /// overflows to NaN or infinity.
    #[inline]
    pub(crate) fn first_asymmetry(
        &self,
        rel_tol: Tolerance,
    ) -> Result<Option<(usize, usize)>, LaError> {
        let eps = self.symmetry_epsilon(rel_tol)?;
        for r in 0..D {
            for c in (r + 1)..D {
                let upper = self.matrix.rows[r][c];
                let lower = self.matrix.rows[c][r];

                let diff = (upper - lower).abs();
                if !diff.is_finite() || diff > eps {
                    cold_path();
                    return Ok(Some((r, c)));
                }
            }
        }
        Ok(None)
    }

    /// Compute the symmetry tolerance scale for a finite matrix.
    fn symmetry_epsilon(&self, rel_tol: Tolerance) -> Result<f64, LaError> {
        let rel_tol = rel_tol.get();
        let mut eps = rel_tol;

        for r in 0..D {
            let mut row_eps = 0.0;
            for c in 0..D {
                row_eps = rel_tol.mul_add(self.matrix.rows[r][c].abs(), row_eps);
                if !row_eps.is_finite() {
                    cold_path();
                    return Err(LaError::non_finite_at(c));
                }
            }
            if row_eps > eps {
                eps = row_eps;
            }
        }

        Ok(eps)
    }

    /// Compute an LU decomposition with partial pivoting.
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] for an unusable pivot, or
    /// [`LaError::NonFinite`] if elimination computes a non-finite intermediate.
    #[inline]
    pub(crate) fn lu(self, tol: Tolerance) -> Result<Lu<D>, LaError> {
        Lu::factor_finite(self, tol)
    }

    /// Compute an LDLT factorization (`A = L D Lᵀ`) without pivoting.
    ///
    /// # Errors
    /// Returns [`LaError::Asymmetric`] if the matrix is not symmetric,
    /// [`LaError::NotPositiveSemidefinite`] for a negative LDLT diagonal,
    /// [`LaError::Singular`] for a zero or too-small non-negative diagonal, or
    /// [`LaError::NonFinite`] if factorization computes a non-finite intermediate.
    #[inline]
    pub(crate) fn ldlt(self, tol: Tolerance) -> Result<Ldlt<D>, LaError> {
        Ldlt::factor_symmetric(SymmetricMatrix::try_new(self)?, tol)
    }

    /// Closed-form determinant for dimensions 0–4, bypassing LU factorization.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the closed-form determinant overflows
    /// to NaN or infinity.
    #[inline]
    pub(crate) const fn det_direct(&self) -> Result<Option<f64>, LaError> {
        match D {
            0 => Ok(Some(1.0)),
            1 => Self::computed_scalar_result(Some(self.matrix.rows[0][0])),
            2 => Self::computed_scalar_result(Some(self.matrix.rows[0][0].mul_add(
                self.matrix.rows[1][1],
                -(self.matrix.rows[0][1] * self.matrix.rows[1][0]),
            ))),
            3 => {
                let m00 = self.matrix.rows[1][1].mul_add(
                    self.matrix.rows[2][2],
                    -(self.matrix.rows[1][2] * self.matrix.rows[2][1]),
                );
                let m01 = self.matrix.rows[1][0].mul_add(
                    self.matrix.rows[2][2],
                    -(self.matrix.rows[1][2] * self.matrix.rows[2][0]),
                );
                let m02 = self.matrix.rows[1][0].mul_add(
                    self.matrix.rows[2][1],
                    -(self.matrix.rows[1][1] * self.matrix.rows[2][0]),
                );
                Self::computed_scalar_result(Some(self.matrix.rows[0][0].mul_add(
                    m00,
                    (-self.matrix.rows[0][1]).mul_add(m01, self.matrix.rows[0][2] * m02),
                )))
            }
            4 => {
                let r = &self.matrix.rows;

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

                Self::computed_scalar_result(Some(r[0][0].mul_add(
                    c00,
                    (-r[0][1]).mul_add(c01, r[0][2].mul_add(c02, -(r[0][3] * c03))),
                )))
            }
            _ => {
                cold_path();
                Ok(None)
            }
        }
    }

    /// Floating-point determinant for a finite matrix.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if a computed determinant or LU fallback
    /// intermediate overflows to NaN or infinity.
    #[inline]
    pub(crate) fn det(self) -> Result<f64, LaError> {
        if let Some(d) = self.det_direct()? {
            return Ok(d);
        }
        match self.lu(Tolerance::new_unchecked(0.0)) {
            Ok(lu) => lu.det(),
            Err(LaError::Singular { .. }) => Ok(0.0),
            Err(err) => Err(err),
        }
    }

    /// Conservative absolute error bound for [`det_direct`](Self::det_direct).
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when the bound computation overflows to NaN
    /// or infinity.
    #[inline]
    pub(crate) const fn det_errbound(&self) -> Result<Option<f64>, LaError> {
        match D {
            0 | 1 => Self::computed_scalar_result(Some(0.0)),
            2 => {
                let r = &self.matrix.rows;
                let permanent = (r[0][0] * r[1][1]).abs() + (r[0][1] * r[1][0]).abs();
                Self::computed_scalar_result(Some(ERR_COEFF_2 * permanent))
            }
            3 => {
                let r = &self.matrix.rows;
                let pm00 = (r[1][1] * r[2][2]).abs() + (r[1][2] * r[2][1]).abs();
                let pm01 = (r[1][0] * r[2][2]).abs() + (r[1][2] * r[2][0]).abs();
                let pm02 = (r[1][0] * r[2][1]).abs() + (r[1][1] * r[2][0]).abs();
                let permanent = r[0][2]
                    .abs()
                    .mul_add(pm02, r[0][1].abs().mul_add(pm01, r[0][0].abs() * pm00));
                Self::computed_scalar_result(Some(ERR_COEFF_3 * permanent))
            }
            4 => {
                let r = &self.matrix.rows;
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
                let permanent = r[0][3].abs().mul_add(
                    pc3,
                    r[0][2]
                        .abs()
                        .mul_add(pc2, r[0][1].abs().mul_add(pc1, r[0][0].abs() * pc0)),
                );
                Self::computed_scalar_result(Some(ERR_COEFF_4 * permanent))
            }
            _ => {
                cold_path();
                Ok(None)
            }
        }
    }

    /// Return a computed scalar result for a matrix with finite stored entries.
    const fn computed_scalar_result(value: Option<f64>) -> Result<Option<f64>, LaError> {
        match value {
            Some(value) if value.is_finite() => Ok(Some(value)),
            Some(_) => Err(LaError::non_finite_at(0)),
            None => Ok(None),
        }
    }
}

impl<const D: usize> Default for FiniteMatrix<D> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl<const D: usize> From<FiniteMatrix<D>> for Matrix<D> {
    #[inline]
    fn from(value: FiniteMatrix<D>) -> Self {
        value.into_matrix()
    }
}

impl<const D: usize> TryFrom<Matrix<D>> for FiniteMatrix<D> {
    type Error = LaError;

    #[inline]
    fn try_from(value: Matrix<D>) -> Result<Self, Self::Error> {
        Self::new(value)
    }
}

impl<const D: usize> TryFrom<[[f64; D]; D]> for FiniteMatrix<D> {
    type Error = LaError;

    #[inline]
    fn try_from(value: [[f64; D]; D]) -> Result<Self, Self::Error> {
        Self::from_rows(value)
    }
}

/// Matrix proven finite and symmetric under the crate's LDLT symmetry tolerance.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
#[allow(clippy::redundant_pub_crate)]
pub(crate) struct SymmetricMatrix<const D: usize> {
    matrix: Matrix<D>,
}

impl<const D: usize> SymmetricMatrix<D> {
    /// Construct a symmetric matrix proof without checking the invariant.
    ///
    /// This constructor is only for paths that have already validated finite
    /// entries and LDLT symmetry with the same predicate as
    /// [`try_new`](Self::try_new).
    #[inline]
    pub(crate) const fn new_unchecked(matrix: Matrix<D>) -> Self {
        Self { matrix }
    }

    /// Validate that a finite matrix is symmetric under the LDLT symmetry tolerance.
    ///
    /// The predicate is the same one used by [`Matrix::ldlt`]:
    /// `|A[i][j] - A[j][i]| <= 1e-12 * max(1, inf_norm(A))`, with scaling that
    /// preserves strict tolerances when an unscaled row sum would overflow.
    ///
    /// # Errors
    /// Returns [`LaError::Asymmetric`] when the first off-diagonal pair violates
    /// the LDLT symmetry predicate.
    ///
    /// Returns [`LaError::NonFinite`] when computing the scaled symmetry
    /// tolerance overflows to NaN or infinity.
    #[inline]
    pub(crate) fn try_new(matrix: FiniteMatrix<D>) -> Result<Self, LaError> {
        if let Some((row, col)) = matrix.first_asymmetry(LDLT_SYMMETRY_REL_TOL)? {
            cold_path();
            Err(LaError::asymmetric(row, col, D))
        } else {
            Ok(Self::new_unchecked(matrix.into_matrix()))
        }
    }

    /// Consume the wrapper and return the underlying matrix.
    #[inline]
    pub(crate) const fn into_matrix(self) -> Matrix<D> {
        self.matrix
    }
}

impl<const D: usize> Matrix<D> {
    /// Test-only infallible constructor for finite literal fixtures.
    #[cfg(test)]
    #[inline]
    pub(crate) const fn from_rows(rows: [[f64; D]; D]) -> Self {
        match Self::try_from_rows(rows) {
            Ok(matrix) => matrix,
            Err(_) => panic!("Matrix::from_rows requires finite entries"),
        }
    }

    /// Try to create a finite matrix from row-major storage.
    ///
    /// This is the public raw-storage boundary for matrices. Public compute
    /// methods parse stored rows into crate-internal proof-bearing types before
    /// arithmetic, including when crate-internal unchecked storage exists.
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
        if let Some((row, col)) = Self::first_non_finite_cell_in(&rows) {
            Err(LaError::non_finite_cell(row, col))
        } else {
            Ok(Self::from_rows_unchecked(rows))
        }
    }

    /// Construct a matrix without checking that entries are finite.
    ///
    /// This crate-internal escape hatch is reserved for literals and algorithm
    /// outputs whose finite invariant is visible at the call site.
    #[inline]
    pub(crate) const fn from_rows_unchecked(rows: [[f64; D]; D]) -> Self {
        Self { rows }
    }

    /// All-zeros matrix.
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

    /// Identity matrix.
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

    /// Get an element with bounds checking.
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
    pub const fn get(&self, r: usize, c: usize) -> Option<f64> {
        if r < D && c < D {
            Some(self.rows[r][c])
        } else {
            None
        }
    }

    /// Get an element, preserving index context on failure.
    ///
    /// Prefer [`get`](Self::get) for const or hot paths that only need
    /// `Option`-style absence.  Use this method at public runtime boundaries
    /// where row, column, and dimension context should survive in a typed error.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// assert_eq!(m.get_checked(1, 0)?, 3.0);
    /// assert_eq!(
    ///     m.get_checked(2, 0),
    ///     Err(LaError::IndexOutOfBounds {
    ///         row: 2,
    ///         col: 0,
    ///         dim: 2,
    ///     })
    /// );
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::IndexOutOfBounds`] when either index is not `< D`.
    #[inline]
    pub const fn get_checked(&self, row: usize, col: usize) -> Result<f64, LaError> {
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
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let mut m = Matrix::<2>::zero();
    /// assert_eq!(m.set(0, 1, 2.5), Ok(()));
    /// assert_eq!(m.get(0, 1), Some(2.5));
    /// assert_eq!(
    ///     m.set(10, 0, 1.0),
    ///     Err(LaError::IndexOutOfBounds {
    ///         row: 10,
    ///         col: 0,
    ///         dim: 2,
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
        self.set_checked(row, col, value)
    }

    /// Set an element, preserving index context on failure.
    ///
    /// The matrix is mutated only when `(row, col)` is in bounds and `value` is
    /// finite.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let mut m = Matrix::<2>::zero();
    /// m.set_checked(0, 1, 2.5)?;
    /// assert_eq!(m.get_checked(0, 1)?, 2.5);
    ///
    /// assert_eq!(
    ///     m.set_checked(10, 0, 1.0),
    ///     Err(LaError::IndexOutOfBounds {
    ///         row: 10,
    ///         col: 0,
    ///         dim: 2,
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
    pub const fn set_checked(&mut self, row: usize, col: usize, value: f64) -> Result<(), LaError> {
        if row >= D || col >= D {
            return Err(LaError::index_out_of_bounds(row, col, D));
        }
        if !value.is_finite() {
            return Err(LaError::non_finite_cell(row, col));
        }
        self.rows[row][col] = value;
        Ok(())
    }

    /// Infinity norm (maximum absolute row sum).
    ///
    /// # Non-finite handling
    /// Raw non-finite entries are rejected by [`try_from_rows`](Self::try_from_rows)
    /// and [`set`](Self::set) before they can be stored.
    ///
    /// Row sums are accumulated in `f64` with ordinary addition.  This method
    /// checks for overflowed accumulators, but it does not provide a certified
    /// absolute rounding bound for the returned norm.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let m = Matrix::<2>::try_from_rows([[1.0, -2.0], [3.0, 4.0]])?;
    /// assert!((m.inf_norm()? - 7.0).abs() <= 1e-12);
    ///
    /// // Raw NaN entries are rejected with coordinates.
    /// assert_eq!(
    ///     Matrix::<2>::try_from_rows([[f64::NAN, 1.0], [2.0, 3.0]]),
    ///     Err(LaError::NonFinite {
    ///         row: Some(0),
    ///         col: 0,
    ///     })
    /// );
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when stored entries are NaN/infinity or a
    /// row sum overflows to NaN or infinity.
    #[inline]
    pub const fn inf_norm(&self) -> Result<f64, LaError> {
        match FiniteMatrix::new(*self) {
            Ok(matrix) => matrix.inf_norm(),
            Err(err) => Err(err),
        }
    }

    /// Returns `true` if the matrix is symmetric within a relative tolerance.
    ///
    /// Two entries `self[r][c]` and `self[c][r]` are considered equal (for the
    /// purposes of symmetry) when
    /// `|self[r][c] - self[c][r]| <= rel_tol * max(1.0, inf_norm(self))`.
    /// This mirrors the predicate used internally by [`ldlt`](Self::ldlt), so
    /// callers can pre-validate matrices that may come from untrusted sources.
    ///
    /// Use [`first_asymmetry`](Self::first_asymmetry) to locate the first
    /// offending pair when this returns `Ok(false)`.
    ///
    /// The `rel_tol` argument is a [`Tolerance`], so raw caller input must be
    /// finite and non-negative before it can reach this predicate. Use
    /// [`Tolerance::new`] or [`LaError::validate_tolerance`] when accepting a
    /// raw `f64`; negative, NaN, and infinite tolerances return
    /// [`LaError::InvalidTolerance`].
    ///
    /// # Overflow handling
    /// A finite matrix can return [`LaError::NonFinite`] if computing the scaled
    /// symmetry tolerance overflows to NaN or infinity.  If both stored entries
    /// are finite but their difference overflows to ±∞, the pair is reported as
    /// asymmetric.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::try_from_rows([[4.0, 2.0], [2.0, 3.0]])?;
    /// let tol = Tolerance::new(1e-12)?;
    /// assert!(a.is_symmetric(tol)?);
    ///
    /// let b = Matrix::<2>::try_from_rows([[4.0, 2.0], [3.0, 3.0]])?;
    /// assert!(!b.is_symmetric(tol)?);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when stored entries are NaN/infinity or
    /// computing the scaled symmetry tolerance overflows to NaN or infinity.
    #[inline]
    pub fn is_symmetric(&self, rel_tol: Tolerance) -> Result<bool, LaError> {
        FiniteMatrix::new(*self)?.is_symmetric(rel_tol)
    }

    /// Returns the indices `(r, c)` (with `r < c`) of the first off-diagonal
    /// pair that violates symmetry, or `None` if the matrix is symmetric
    /// within `rel_tol`.
    ///
    /// Iteration order is row-major over the strict upper triangle, so the
    /// returned indices are the lexicographically smallest such pair.  The
    /// predicate is the same as [`is_symmetric`](Self::is_symmetric):
    /// `|self[r][c] - self[c][r]| <= rel_tol * max(1.0, inf_norm(self))`.
    ///
    /// A finite matrix can return [`LaError::NonFinite`] if computing the scaled
    /// symmetry tolerance overflows to NaN or infinity. If both stored entries
    /// are finite but their difference overflows to ±∞, the pair is reported as
    /// asymmetric.
    ///
    /// The `rel_tol` argument is a [`Tolerance`], so raw caller input must be
    /// finite and non-negative before it can reach this predicate. Use
    /// [`Tolerance::new`] or [`LaError::validate_tolerance`] when accepting a
    /// raw `f64`; negative, NaN, and infinite tolerances return
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
    /// let tol = Tolerance::new(1e-12)?;
    /// assert_eq!(a.first_asymmetry(tol)?, Some((1, 2)));
    /// assert_eq!(Matrix::<3>::identity().first_asymmetry(tol)?, None);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when stored entries are NaN/infinity or
    /// computing the scaled symmetry tolerance overflows to NaN or infinity.
    #[inline]
    pub fn first_asymmetry(&self, rel_tol: Tolerance) -> Result<Option<(usize, usize)>, LaError> {
        FiniteMatrix::new(*self)?.first_asymmetry(rel_tol)
    }

    /// Compute an LU decomposition with partial pivoting.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let lu = a.lu(DEFAULT_PIVOT_TOL)?;
    ///
    /// let b = Vector::<2>::try_new([5.0, 11.0])?;
    /// let x = lu.solve_vec(b)?.into_array();
    ///
    /// assert!((x[0] - 1.0).abs() <= 1e-12);
    /// assert!((x[1] - 2.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// The `tol` argument is a [`Tolerance`], so raw caller input must be
    /// finite and non-negative before it can reach factorization. Use
    /// [`Tolerance::new`] or [`LaError::validate_tolerance`] when accepting a
    /// raw `f64`; negative, NaN, and infinite tolerances return
    /// [`LaError::InvalidTolerance`].
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if, for some column `k`, the largest-magnitude candidate pivot
    /// in that column satisfies `|pivot| <= tol` (so no numerically usable pivot exists).
    /// Returns [`LaError::NonFinite`] if stored entries are NaN/infinity or an
    /// elimination intermediate overflows to NaN/∞ before it can be stored in
    /// the returned [`Lu`].
    #[inline]
    pub fn lu(self, tol: Tolerance) -> Result<Lu<D>, LaError> {
        FiniteMatrix::new(self)?.lu(tol)
    }

    /// Compute an LDLT factorization (`A = L D Lᵀ`) without pivoting.
    ///
    /// This is intended for symmetric positive definite (SPD) and positive semi-definite (PSD)
    /// matrices such as Gram matrices.
    ///
    /// # Symmetry validation
    /// The input matrix `self` must be symmetric — that is,
    /// `self[i][j] == self[j][i]` within the crate's LDLT symmetry tolerance
    /// (`1e-12`, scaled like [`is_symmetric`](Self::is_symmetric)).  This is a
    /// correctness invariant, not merely a performance hint, so asymmetric inputs return
    /// [`LaError::Asymmetric`] before factorization starts.  If you need a
    /// general-purpose factorization that tolerates non-symmetric inputs, use
    /// [`lu`](Self::lu) instead.
    ///
    /// The `tol` argument is a [`Tolerance`], so raw caller input must be
    /// finite and non-negative before it can reach factorization. Use
    /// [`Tolerance::new`] or [`LaError::validate_tolerance`] when accepting a
    /// raw `f64`; negative, NaN, and infinite tolerances return
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
    /// let x = ldlt.solve_vec(b)?.into_array();
    /// assert!((x[0] - (-0.125)).abs() <= 1e-12);
    /// assert!((x[1] - 0.75).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NotPositiveSemidefinite`] if, for some step `k`, the required
    /// diagonal entry `d = D[k,k]` is negative.
    /// Returns [`LaError::Singular`] if `0 <= d <= tol`, treating PSD degeneracy
    /// as singular/degenerate.
    /// Returns [`LaError::NonFinite`] if stored entries are NaN/infinity or
    /// factorization computes a non-finite intermediate.
    /// Returns [`LaError::Asymmetric`] if the input matrix is not symmetric.
    #[inline]
    pub fn ldlt(self, tol: Tolerance) -> Result<Ldlt<D>, LaError> {
        FiniteMatrix::new(self)?.ldlt(tol)
    }

    /// Return the first non-finite stored cell in row-major order.
    const fn first_non_finite_cell_in(rows: &[[f64; D]; D]) -> Option<(usize, usize)> {
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
    /// Returns [`LaError::NonFinite`] when stored entries are NaN/infinity or
    /// the closed-form determinant overflows to NaN or infinity.
    #[inline]
    pub const fn det_direct(&self) -> Result<Option<f64>, LaError> {
        match FiniteMatrix::new(*self) {
            Ok(matrix) => matrix.det_direct(),
            Err(err) => Err(err),
        }
    }

    /// Floating-point determinant, using closed-form formulas for D ≤ 4 and
    /// LU decomposition for D ≥ 5.
    ///
    /// For D ∈ {1, 2, 3, 4}, this bypasses LU factorization entirely for a significant
    /// speedup (see [`det_direct`](Self::det_direct)).
    ///
    /// Finite inputs return a floating-point determinant estimate in every dimension;
    /// this method does not surface [`LaError::Singular`]. Because it mixes
    /// closed-form paths from [`det_direct`](Self::det_direct) with an LU fallback,
    /// the returned value has no certified absolute error bound. Use
    /// [`det_errbound`](Self::det_errbound) for D ≤ 4 bounds, or the exact
    /// determinant APIs when exact singularity classification or certified values
    /// matter. For D ≥ 5, the LU fallback only maps an exactly zero pivot to
    /// `Ok(0.0)`. Use [`lu`](Self::lu) directly when you need tolerance-aware
    /// singularity detection or the pivot column.
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
    /// # Errors
    /// Returns [`LaError::NonFinite`] if stored entries are NaN/infinity, the
    /// LU fallback computes a non-finite factorization cell, or the determinant
    /// product overflows to NaN or infinity.
    #[inline]
    pub fn det(self) -> Result<f64, LaError> {
        FiniteMatrix::new(self)?.det()
    }

    /// Conservative absolute error bound for `det_direct()`.
    ///
    /// Returns `Ok(Some(bound))` such that `|det_direct() - det_exact| ≤ bound`,
    /// or `Ok(None)` for D ≥ 5 where no fast bound is available.
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
    /// Use this to build adaptive-precision logic: if `|det_direct()| > bound`,
    /// the f64 sign is provably correct. Otherwise fall back to exact arithmetic.
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
    /// if let (Some(bound), Some(det_approx)) = (m.det_errbound()?, m.det_direct()?) {
    ///     // If |det_approx| > bound, the sign is guaranteed correct.
    ///     let sign_is_certified = det_approx.abs() > bound;
    ///     assert!(!sign_is_certified);
    /// }
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Adaptive precision pattern (requires `exact` feature)
    /// ```ignore
    /// use la_stack::prelude::*;
    ///
    /// let m = Matrix::<3>::identity();
    /// if let Some(bound) = m.det_errbound()? {
    ///     if let Some(det) = m.det_direct()? {
    ///         if det.abs() > bound {
    ///             // f64 sign is guaranteed correct
    ///             let sign = det.signum() as i8;
    ///         } else {
    ///             // Fall back to exact arithmetic (requires `exact` feature)
    ///             let sign = m.det_sign_exact()?;
    ///         }
    ///     }
    /// } else {
    ///     // D ≥ 5: no fast filter, use exact directly
    ///     let sign = m.det_sign_exact()?;
    /// }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] when stored entries are NaN/infinity or
    /// the bound computation overflows to NaN or infinity.
    #[inline]
    pub const fn det_errbound(&self) -> Result<Option<f64>, LaError> {
        match FiniteMatrix::new(*self) {
            Ok(matrix) => matrix.det_errbound(),
            Err(err) => Err(err),
        }
    }
}

impl<const D: usize> Default for Matrix<D> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_PIVOT_TOL;
    use crate::vector::{FiniteVector, Vector};

    use approx::assert_abs_diff_eq;
    use core::assert_matches;
    use pastey::paste;
    use std::hint::black_box;

    fn assert_matrix_abs_eq<const D: usize>(actual: &Matrix<D>, expected: &Matrix<D>) {
        for r in 0..D {
            for c in 0..D {
                assert_abs_diff_eq!(
                    actual.get(r, c).unwrap(),
                    expected.get(r, c).unwrap(),
                    epsilon = 0.0
                );
            }
        }
    }

    fn assert_array_abs_eq<const D: usize>(actual: &[f64; D], expected: &[f64; D]) {
        for i in 0..D {
            assert_abs_diff_eq!(actual[i], expected[i], epsilon = 0.0);
        }
    }

    macro_rules! gen_public_api_matrix_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_api_matrix_from_rows_get_set_bounds_checked_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = 1.0;
                    rows[$d - 1][$d - 1] = -2.0;

                    let mut m = Matrix::<$d>::from_rows(rows);

                    assert_eq!(m.get(0, 0), Some(1.0));
                    assert_eq!(m.get($d - 1, $d - 1), Some(-2.0));
                    assert_eq!(m.get_checked(0, 0), Ok(1.0));
                    assert_eq!(m.get_checked($d - 1, $d - 1), Ok(-2.0));

                    // Out-of-bounds is None.
                    assert_eq!(m.get($d, 0), None);
                    assert_eq!(
                        m.get_checked($d, 0),
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
                        m.set_checked($d, 0, 3.0),
                        Err(LaError::IndexOutOfBounds {
                            row: $d,
                            col: 0,
                            dim: $d,
                        })
                    );
                    assert_eq!(m, before_failed_set);
                    assert_eq!(
                        m.set_checked(0, $d, 3.0),
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
                    assert_eq!(m.set_checked($d - 1, 0, 4.0), Ok(()));
                    assert_eq!(m.get_checked($d - 1, 0), Ok(4.0));
                }

                #[test]
                fn [<public_api_matrix_zero_and_default_are_zero_ $d d>]() {
                    let z = Matrix::<$d>::zero();
                    assert_abs_diff_eq!(z.inf_norm().unwrap(), 0.0, epsilon = 0.0);

                    let d = Matrix::<$d>::default();
                    assert_abs_diff_eq!(d.inf_norm().unwrap(), 0.0, epsilon = 0.0);
                }

                #[test]
                fn [<public_api_matrix_inf_norm_max_row_sum_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];

                    // Row 0 has absolute row sum = D.
                    for c in 0..$d {
                        rows[0][c] = -1.0;
                    }

                    // Row 1 has smaller absolute row sum.
                    for c in 0..$d {
                        rows[1][c] = 0.5;
                    }

                    let m = Matrix::<$d>::from_rows(rows);
                    assert_abs_diff_eq!(m.inf_norm().unwrap(), f64::from($d), epsilon = 0.0);
                }

                #[test]
                fn [<public_api_matrix_identity_lu_det_solve_vec_ $d d>]() {
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
                    let lu = m.lu(DEFAULT_PIVOT_TOL).unwrap();

                    let b_arr = {
                        let mut arr = [0.0f64; $d];
                        let values = [1.0f64, 2.0, 3.0, 4.0, 5.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = *src;
                        }
                        arr
                    };

                    let b = Vector::<$d>::new(b_arr);
                    let x = lu.solve_vec(b).unwrap().into_array();

                    for (x_i, b_i) in x.iter().zip(b_arr.iter()) {
                        assert_abs_diff_eq!(*x_i, *b_i, epsilon = 1e-12);
                    }
                }

                #[test]
                fn [<finite_matrix_accepts_and_roundtrips_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        rows[r][r] = 1.0;
                    }

                    let finite = FiniteMatrix::<$d>::from_rows(rows).unwrap();

                    assert_matrix_abs_eq(&finite.into_matrix(), &Matrix::<$d>::from_rows(rows));
                    assert_matrix_abs_eq(&FiniteMatrix::<$d>::zero().into_matrix(), &Matrix::<$d>::zero());
                    assert_matrix_abs_eq(&FiniteMatrix::<$d>::default().into_matrix(), &Matrix::<$d>::zero());
                    assert_eq!(finite.into_matrix().get(0, 0), Some(1.0));
                    assert_eq!(finite.into_matrix().get($d, 0), None);
                    assert_eq!(
                        finite.into_matrix().get_checked($d, 0),
                        Err(LaError::IndexOutOfBounds {
                            row: $d,
                            col: 0,
                            dim: $d,
                        })
                    );
                    assert_matrix_abs_eq(&Matrix::from(finite), &Matrix::<$d>::from_rows(rows));
                    assert_matrix_abs_eq(
                        &FiniteMatrix::<$d>::try_from(Matrix::<$d>::from_rows(rows))
                            .unwrap()
                            .into_matrix(),
                        &finite.into_matrix()
                    );
                    assert_matrix_abs_eq(
                        &FiniteMatrix::<$d>::try_from(rows).unwrap().into_matrix(),
                        &finite.into_matrix()
                    );
                }

                #[test]
                fn [<finite_matrix_rejects_nonfinite_with_coordinates_ $d d>]() {
                    let mut rows = [[1.0f64; $d]; $d];
                    rows[$d - 1][0] = f64::NAN;

                    assert_eq!(
                        FiniteMatrix::<$d>::from_rows(rows),
                        Err(LaError::NonFinite {
                            row: Some($d - 1),
                            col: 0,
                        })
                    );
                }

                #[test]
                fn [<finite_matrix_try_from_raw_matrix_revalidates_entries_ $d d>]() {
                    let mut rows = [[1.0f64; $d]; $d];
                    rows[$d - 1][$d - 1] = f64::INFINITY;
                    let raw = Matrix::<$d>::from_rows_unchecked(rows);

                    assert_eq!(
                        FiniteMatrix::<$d>::try_from(raw),
                        Err(LaError::NonFinite {
                            row: Some($d - 1),
                            col: $d - 1,
                        })
                    );
                }

                #[test]
                fn [<finite_matrix_algorithms_match_raw_boundary_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    let diag_values = [2.0f64, 3.0, 5.0, 7.0, 11.0];
                    for i in 0..$d {
                        rows[i][i] = diag_values[i];
                    }

                    let raw = Matrix::<$d>::from_rows(rows);
                    let finite = FiniteMatrix::<$d>::from_rows(rows).unwrap();

                    assert_abs_diff_eq!(finite.inf_norm().unwrap(), raw.inf_norm().unwrap(), epsilon = 0.0);
                    assert_eq!(finite.det_direct(), raw.det_direct());
                    assert_abs_diff_eq!(finite.det().unwrap(), raw.det().unwrap(), epsilon = 0.0);
                    assert_eq!(finite.det_errbound(), raw.det_errbound());
                    assert_eq!(
                        finite.first_asymmetry(Tolerance::new(1e-12).unwrap()).unwrap(),
                        raw.first_asymmetry(Tolerance::new(1e-12).unwrap()).unwrap()
                    );
                    assert_eq!(
                        finite.is_symmetric(Tolerance::new(1e-12).unwrap()).unwrap(),
                        raw.is_symmetric(Tolerance::new(1e-12).unwrap()).unwrap()
                    );
                }

                #[test]
                fn [<finite_matrix_lu_and_ldlt_accept_finite_rhs_ $d d>]() {
                    let finite = FiniteMatrix::<$d>::new(Matrix::<$d>::identity()).unwrap();
                    let rhs = {
                        let mut arr = [0.0f64; $d];
                        let values = [1.0f64, 2.0, 3.0, 4.0, 5.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = *src;
                        }
                        FiniteVector::<$d>::from_array(arr).unwrap()
                    };

                    let lu = finite.lu(DEFAULT_PIVOT_TOL).unwrap();
                    assert_array_abs_eq(
                        &lu.solve_finite_vec(rhs).unwrap().into_array(),
                        &rhs.into_array()
                    );

                    let ldlt = finite.ldlt(DEFAULT_PIVOT_TOL).unwrap();
                    assert_array_abs_eq(
                        &ldlt.solve_finite_vec(rhs).unwrap().into_array(),
                        &rhs.into_array()
                    );
                }
            }
        };
    }

    // Mirror delaunay-style multi-dimension tests.
    gen_public_api_matrix_tests!(2);
    gen_public_api_matrix_tests!(3);
    gen_public_api_matrix_tests!(4);
    gen_public_api_matrix_tests!(5);

    macro_rules! gen_public_matrix_revalidation_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_matrix_compute_methods_revalidate_unchecked_storage_ $d d>]() {
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 1.0;
                    }
                    rows[0][$d - 1] = f64::NAN;
                    let raw = Matrix::<$d>::from_rows_unchecked(rows);
                    let expected = LaError::NonFinite {
                        row: Some(0),
                        col: $d - 1,
                    };

                    assert_eq!(raw.inf_norm(), Err(expected));
                    assert_eq!(
                        raw.is_symmetric(Tolerance::new(0.0).unwrap()),
                        Err(expected)
                    );
                    assert_eq!(
                        raw.first_asymmetry(Tolerance::new(0.0).unwrap()),
                        Err(expected)
                    );
                    assert_eq!(raw.lu(DEFAULT_PIVOT_TOL), Err(expected));
                    assert_eq!(raw.ldlt(DEFAULT_PIVOT_TOL), Err(expected));
                    assert_eq!(raw.det_direct(), Err(expected));
                    assert_eq!(raw.det(), Err(expected));
                    assert_eq!(raw.det_errbound(), Err(expected));
                }
            }
        };
    }

    gen_public_matrix_revalidation_tests!(2);
    gen_public_matrix_revalidation_tests!(3);
    gen_public_matrix_revalidation_tests!(4);
    gen_public_matrix_revalidation_tests!(5);

    #[test]
    fn public_matrix_fast_none_paths_revalidate_unchecked_storage() {
        let mut rows = [[1.0f64; 5]; 5];
        rows[4][1] = f64::INFINITY;
        let raw = Matrix::<5>::from_rows_unchecked(rows);

        let expected = Err(LaError::NonFinite {
            row: Some(4),
            col: 1,
        });

        assert_eq!(raw.det_direct(), expected);
        assert_eq!(raw.det_errbound(), expected);
    }

    // === det_direct tests ===

    #[test]
    fn det_direct_d0_is_one() {
        assert_eq!(Matrix::<0>::zero().det_direct(), Ok(Some(1.0)));
    }

    #[test]
    fn det_direct_d1_returns_element() {
        let m = Matrix::<1>::from_rows([[42.0]]);
        assert_eq!(m.det_direct(), Ok(Some(42.0)));
    }

    #[test]
    fn det_direct_d2_known_value() {
        // [[1,2],[3,4]] → det = 1*4 - 2*3 = -2
        // black_box prevents compile-time constant folding of the const fn.
        let m = black_box(Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]));
        assert_abs_diff_eq!(m.det_direct().unwrap().unwrap(), -2.0, epsilon = 1e-15);
    }

    #[test]
    fn det_direct_d3_known_value() {
        // Classic 3×3: det = 0
        let m = black_box(Matrix::<3>::from_rows([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]));
        assert_abs_diff_eq!(m.det_direct().unwrap().unwrap(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn det_direct_d3_nonsingular() {
        // [[2,1,0],[0,3,1],[1,0,2]] → det = 2*(6-0) - 1*(0-1) + 0 = 13
        let m = black_box(Matrix::<3>::from_rows([
            [2.0, 1.0, 0.0],
            [0.0, 3.0, 1.0],
            [1.0, 0.0, 2.0],
        ]));
        assert_abs_diff_eq!(m.det_direct().unwrap().unwrap(), 13.0, epsilon = 1e-12);
    }

    #[test]
    fn det_direct_d4_identity() {
        let m = black_box(Matrix::<4>::identity());
        assert_abs_diff_eq!(m.det_direct().unwrap().unwrap(), 1.0, epsilon = 1e-15);
    }

    #[test]
    fn det_direct_d4_known_value() {
        // Diagonal matrix: det = product of diagonal entries.
        let mut rows = [[0.0f64; 4]; 4];
        rows[0][0] = 2.0;
        rows[1][1] = 3.0;
        rows[2][2] = 5.0;
        rows[3][3] = 7.0;
        let m = black_box(Matrix::<4>::from_rows(rows));
        assert_abs_diff_eq!(m.det_direct().unwrap().unwrap(), 210.0, epsilon = 1e-12);
    }

    #[test]
    fn det_direct_d5_returns_none() {
        assert_eq!(Matrix::<5>::identity().det_direct(), Ok(None));
    }

    #[test]
    fn det_direct_d5_rejects_nonfinite_before_returning_none() {
        let mut m = Matrix::<5>::identity();
        assert_eq!(
            m.set(3, 4, f64::NAN),
            Err(LaError::NonFinite {
                row: Some(3),
                col: 4,
            })
        );
    }

    #[test]
    fn det_direct_d8_returns_none() {
        assert_eq!(Matrix::<8>::zero().det_direct(), Ok(None));
    }

    #[test]
    fn det_direct_rejects_nonfinite_entry_with_coordinates() {
        assert_eq!(
            Matrix::<3>::try_from_rows([[1.0, 0.0, 0.0], [0.0, f64::NAN, 0.0], [0.0, 0.0, 1.0]]),
            Err(LaError::NonFinite {
                row: Some(1),
                col: 1,
            })
        );
    }

    #[test]
    fn det_direct_rejects_computed_overflow() {
        let m = Matrix::<2>::from_rows([[1e300, 0.0], [0.0, 1e300]]);
        assert_eq!(
            m.det_direct(),
            Err(LaError::NonFinite { row: None, col: 0 })
        );
    }

    #[test]
    fn det_d5_rejects_lu_product_overflow() {
        let m = Matrix::<5>::from_rows([
            [1.0e100, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0e100, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0e100, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0e100, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0e100],
        ]);
        assert_eq!(m.det(), Err(LaError::NonFinite { row: None, col: 3 }));
    }

    #[test]
    fn det_d5_rejects_lu_trailing_update_overflow() {
        let m = Matrix::<5>::from_rows([
            [1.0, f64::MAX, 0.0, 0.0, 0.0],
            [-1.0, f64::MAX, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);

        assert_eq!(
            m.det(),
            Err(LaError::NonFinite {
                row: Some(1),
                col: 1,
            })
        );
    }

    macro_rules! gen_det_direct_agrees_with_lu {
        ($d:literal) => {
            paste! {
                #[test]
                #[allow(clippy::cast_precision_loss)] // r, c, D are tiny integers
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
                    let m = Matrix::<$d>::from_rows(rows);
                    let direct = m.det_direct().unwrap().unwrap();
                    let lu_det = m.lu(DEFAULT_PIVOT_TOL).unwrap().det().unwrap();
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
    gen_det_singular_zero_matrix_tests!(5);

    #[test]
    fn det_d5_ignores_pivot_tolerance_for_tiny_nonsingular_matrix() {
        // A small nonzero determinant is still a determinant. `det` must not
        // flatten the value to zero merely because the default LU tolerance
        // would reject a pivot this small.
        let m = Matrix::<5>::from_rows([
            [1e-13, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]);

        assert_abs_diff_eq!(m.det().unwrap(), 1e-13, epsilon = 0.0);
        assert_eq!(
            m.lu(DEFAULT_PIVOT_TOL),
            Err(LaError::Singular { pivot_col: 0 })
        );
    }

    #[test]
    fn det_returns_nonfinite_error_for_nan_d2() {
        assert_eq!(
            Matrix::<2>::try_from_rows([[f64::NAN, 1.0], [1.0, 1.0]]),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_returns_nonfinite_error_for_inf_d3() {
        assert_eq!(
            Matrix::<3>::try_from_rows([
                [f64::INFINITY, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0]
            ]),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0
            })
        );
    }

    #[test]
    fn det_returns_nonfinite_error_for_overflow_with_finite_entries() {
        // det_direct produces an overflowing f64 (1e300 * 1e300 = ∞) even
        // though every matrix entry is finite.  The entry scan in `det`
        // falls through and returns NonFinite { row: None, col: 0 } to signal
        // a computed overflow rather than a NaN/∞ input.
        let m = Matrix::<2>::from_rows([[1e300, 0.0], [0.0, 1e300]]);
        assert_eq!(m.det(), Err(LaError::NonFinite { row: None, col: 0 }));
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
    fn det_errbound_available_without_exact_feature() {
        // Verify det_errbound is accessible without exact feature
        let m = Matrix::<3>::identity();
        let bound = m.det_errbound().unwrap();
        assert!(bound.is_some());
        assert!(bound.unwrap() > 0.0);
    }

    #[test]
    fn det_errbound_matches_documented_coefficient_scale() {
        let m2 = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
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
    fn det_errbound_d5_returns_none() {
        // D=5 has no fast filter
        assert_eq!(Matrix::<5>::identity().det_errbound(), Ok(None));
    }

    #[test]
    fn det_errbound_d1_rejects_nonfinite_even_with_zero_bound() {
        assert_eq!(
            Matrix::<1>::try_from_rows([[f64::INFINITY]]),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0,
            })
        );
    }

    #[test]
    fn det_errbound_d5_rejects_nonfinite_before_returning_none() {
        let mut m = Matrix::<5>::identity();
        assert_eq!(
            m.set(4, 1, f64::NAN),
            Err(LaError::NonFinite {
                row: Some(4),
                col: 1,
            })
        );
    }

    #[test]
    fn det_errbound_rejects_nonfinite_entry_with_coordinates() {
        assert_eq!(
            Matrix::<2>::try_from_rows([[1.0, f64::INFINITY], [0.0, 1.0]]),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 1,
            })
        );
    }

    #[test]
    fn det_errbound_rejects_computed_overflow() {
        let m = Matrix::<2>::from_rows([[1e300, 0.0], [0.0, 1e300]]);
        assert_eq!(
            m.det_errbound(),
            Err(LaError::NonFinite { row: None, col: 0 })
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

    // === inf_norm NaN / Inf rejection (regression tests for #85) ===

    macro_rules! gen_inf_norm_nonfinite_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<inf_norm_all_nan_returns_nonfinite_error_ $d d>]() {
                    // Before the fix, `NaN > max_row_sum` was always false, so a
                    // matrix full of NaN silently produced inf_norm == 0.0.
                    assert_eq!(
                        Matrix::<$d>::try_from_rows([[f64::NAN; $d]; $d]),
                        Err(LaError::NonFinite {
                            row: Some(0),
                            col: 0,
                        })
                    );
                }

                #[test]
                fn [<inf_norm_single_nan_entry_returns_nonfinite_error_ $d d>]() {
                    // A single NaN entry must surface with its source coordinates.
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = f64::NAN;
                    rows[$d - 1][$d - 1] = 1.0;
                    assert_eq!(
                        Matrix::<$d>::try_from_rows(rows),
                        Err(LaError::NonFinite {
                            row: Some(0),
                            col: 0,
                        })
                    );
                }

                #[test]
                fn [<inf_norm_infinity_entry_returns_nonfinite_error_ $d d>]() {
                    // Infinity entries should be rejected with their source coordinates.
                    let mut rows = [[0.0f64; $d]; $d];
                    rows[0][0] = f64::INFINITY;
                    assert_eq!(
                        Matrix::<$d>::try_from_rows(rows),
                        Err(LaError::NonFinite {
                            row: Some(0),
                            col: 0,
                        })
                    );
                }
            }
        };
    }

    gen_inf_norm_nonfinite_tests!(2);
    gen_inf_norm_nonfinite_tests!(3);
    gen_inf_norm_nonfinite_tests!(4);
    gen_inf_norm_nonfinite_tests!(5);

    // === is_symmetric / first_asymmetry (public LDLT preconditions helpers) ===

    macro_rules! gen_is_symmetric_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<is_symmetric_true_for_identity_ $d d>]() {
                    let m = Matrix::<$d>::identity();
                    assert!(m.is_symmetric(Tolerance::new(1e-12).unwrap()).unwrap());
                    assert_eq!(m.first_asymmetry(Tolerance::new(1e-12).unwrap()).unwrap(), None);
                }

                #[test]
                fn [<is_symmetric_true_for_zero_ $d d>]() {
                    let m = Matrix::<$d>::zero();
                    assert!(m.is_symmetric(Tolerance::new(1e-12).unwrap()).unwrap());
                    assert_eq!(m.first_asymmetry(Tolerance::new(1e-12).unwrap()).unwrap(), None);
                }

                #[test]
                fn [<is_symmetric_true_for_constructed_symmetric_ $d d>]() {
                    // Construct A = M + Mᵀ so A is provably symmetric.
                    let mut m = [[0.0f64; $d]; $d];
                    for r in 0..$d {
                        for c in 0..$d {
                            #[allow(clippy::cast_precision_loss)]
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
                    let a = Matrix::<$d>::from_rows(sym);
                    assert!(a.is_symmetric(Tolerance::new(1e-12).unwrap()).unwrap());
                    assert_eq!(a.first_asymmetry(Tolerance::new(1e-12).unwrap()).unwrap(), None);
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
                    let a = Matrix::<$d>::from_rows(rows);
                    assert!(!a.is_symmetric(Tolerance::new(1e-12).unwrap()).unwrap());
                    assert_eq!(
                        a.first_asymmetry(Tolerance::new(1e-12).unwrap()).unwrap(),
                        Some((0, $d - 1))
                    );
                }

                #[test]
                fn [<is_symmetric_rejects_nan_offdiagonal_ $d d>]() {
                    // A NaN off-diagonal is a stored non-finite matrix value,
                    // not merely a symmetry mismatch.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 1.0;
                    }
                    rows[0][1] = f64::NAN;
                    rows[1][0] = f64::NAN;
                    assert_eq!(
                        Matrix::<$d>::try_from_rows(rows),
                        Err(LaError::NonFinite {
                            row: Some(0),
                            col: 1,
                        })
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
                fn [<matrix_ldlt_accepts_identity_ $d d>]() {
                    let ldlt = Matrix::<$d>::identity().ldlt(DEFAULT_PIVOT_TOL).unwrap();
                    assert_abs_diff_eq!(ldlt.det().unwrap(), 1.0, epsilon = 0.0);
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
                        Matrix::<$d>::try_from_rows(rows)
                            .and_then(FiniteMatrix::new)
                            .and_then(SymmetricMatrix::try_new),
                        Err(LaError::Asymmetric {
                            row: 0,
                            col: $d - 1,
                            dim: $d,
                        })
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
    fn matrix_ldlt_rejects_nonfinite_before_asymmetry() {
        assert_eq!(
            Matrix::<2>::try_from_rows([[1.0, f64::NAN], [0.0, 1.0]]),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 1,
            })
        );
    }

    #[test]
    fn symmetric_matrix_into_matrix_roundtrips_storage_internally() {
        let a = Matrix::<2>::from_rows([[2.0, 1.0], [1.0, 3.0]]);
        let symmetric = SymmetricMatrix::try_new(FiniteMatrix::new(a).unwrap()).unwrap();

        assert_eq!(symmetric.into_matrix(), a);
    }

    #[test]
    fn is_symmetric_tolerance_scales_with_inf_norm() {
        // Off-diagonal entries differ by 1e-6.  With inf_norm ≈ 2e6, the
        // relative tolerance 1e-12 yields eps ≈ 2e-6, which accepts the gap;
        // a stricter tol of 1e-15 rejects it.
        let a = Matrix::<2>::from_rows([[1.0e6, 1.0e6 + 1.0e-6], [1.0e6, 1.0e6]]);
        assert!(a.is_symmetric(Tolerance::new(1e-12).unwrap()).unwrap());
        assert!(!a.is_symmetric(Tolerance::new(1e-15).unwrap()).unwrap());
    }

    #[test]
    fn first_asymmetry_returns_lexicographically_first_pair() {
        // Two asymmetric pairs: (0, 2) and (1, 2).  We must get (0, 2) first.
        let a = Matrix::<3>::from_rows([[1.0, 0.0, 2.0], [0.0, 1.0, 3.0], [-2.0, -3.0, 1.0]]);
        assert_eq!(
            a.first_asymmetry(Tolerance::new(1e-12).unwrap()).unwrap(),
            Some((0, 2))
        );
    }

    #[test]
    fn first_asymmetry_rejects_infinite_offdiagonal() {
        assert_eq!(
            Matrix::<2>::try_from_rows([[1.0, f64::INFINITY], [0.0, 1.0]]),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 1,
            })
        );
    }

    #[test]
    fn first_asymmetry_rejects_nan_diagonal() {
        assert_eq!(
            Matrix::<2>::try_from_rows([[f64::NAN, 1.0], [1.0, 1.0]]),
            Err(LaError::NonFinite {
                row: Some(0),
                col: 0,
            })
        );
    }

    #[test]
    fn first_asymmetry_strict_tol_survives_row_sum_overflow() {
        let a = Matrix::<3>::from_rows([
            [f64::MAX, 1.0, f64::MAX],
            [2.0, f64::MAX, 0.0],
            [f64::MAX, 0.0, f64::MAX],
        ]);

        assert_eq!(a.inf_norm(), Err(LaError::NonFinite { row: None, col: 0 }));
        assert_eq!(
            a.first_asymmetry(Tolerance::new(0.0).unwrap()).unwrap(),
            Some((0, 1))
        );
        assert!(!a.is_symmetric(Tolerance::new(0.0).unwrap()).unwrap());
    }

    #[test]
    fn first_asymmetry_rejects_scaled_epsilon_overflow() {
        let a = Matrix::<2>::from_rows([[2.0, 0.0], [0.0, 1.0]]);
        let tol = Tolerance::new(f64::MAX).unwrap();

        assert_eq!(
            a.first_asymmetry(tol),
            Err(LaError::NonFinite { row: None, col: 0 })
        );
        assert_eq!(
            a.is_symmetric(tol),
            Err(LaError::NonFinite { row: None, col: 0 })
        );
    }

    #[test]
    fn first_asymmetry_flags_overflowed_finite_difference() {
        let a = Matrix::<2>::from_rows([[1.0, f64::MAX], [-f64::MAX, 1.0]]);
        assert_eq!(
            a.first_asymmetry(Tolerance::new(1e-12).unwrap()).unwrap(),
            Some((0, 1))
        );
        assert!(!a.is_symmetric(Tolerance::new(1e-12).unwrap()).unwrap());
    }

    #[test]
    fn is_symmetric_rejects_invalid_tol() {
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
    }

    #[test]
    fn first_asymmetry_rejects_negative_tol() {
        assert_eq!(
            Tolerance::new(-1.0),
            Err(LaError::InvalidTolerance { value: -1.0 })
        );
    }

    #[test]
    fn first_asymmetry_rejects_nonfinite_tol() {
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
    }
}
