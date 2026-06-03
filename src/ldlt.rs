//! LDLT factorization and solves.
//!
//! This module provides a stack-allocated LDLT factorization (`A = L D Lᵀ`) intended for
//! symmetric positive definite (SPD) and positive semi-definite (PSD) matrices (e.g. Gram
//! matrices) without pivoting.
//!
//! # Preconditions
//! The input matrix must be **symmetric**.  This is a correctness contract, not a hint:
//! the factorization algorithm reads only the lower triangle and implicitly assumes the
//! upper triangle mirrors it.  Asymmetric inputs return [`LaError::Asymmetric`]
//! before factorization starts.  Callers who know their matrices may not be
//! symmetric at all should use [`crate::Lu`] instead.

use core::hint::cold_path;

use crate::matrix::{Matrix, SymmetricMatrix};
use crate::vector::{FiniteVector, Vector};
use crate::{LaError, Tolerance};

/// LDLT factorization (`A = L D Lᵀ`) for symmetric positive (semi)definite matrices.
///
/// This factorization is **not** a general-purpose symmetric-indefinite LDLT (no pivoting).
/// It assumes the input matrix is symmetric and (numerically) SPD/PSD.
///
/// # Preconditions
/// The source matrix passed to [`Matrix::ldlt`](crate::Matrix::ldlt) must be
/// symmetric (`A[i][j] == A[j][i]` within rounding).  Asymmetric inputs return
/// [`LaError::Asymmetric`] before factorization starts; see
/// [`Matrix::ldlt`](crate::Matrix::ldlt) for details and alternatives.
///
/// # Storage
/// The factors are stored in a single [`Matrix`]:
/// - `D` is stored on the diagonal.
/// - The strict lower triangle stores the multipliers of `L`.
/// - The diagonal of `L` is implicit ones.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ldlt<const D: usize> {
    factors: LdltFactors<D>,
}

/// In-place LDLT factor storage whose diagonal entries are finite and usable.
///
/// Construction through [`Ldlt::factor_symmetric`] proves every stored entry is
/// finite and every diagonal satisfies the factorization tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
struct LdltFactors<const D: usize> {
    storage: Matrix<D>,
}

impl<const D: usize> LdltFactors<D> {
    /// Construct factors after LDLT factorization has proven the storage invariant.
    #[inline]
    const fn new_unchecked(storage: Matrix<D>) -> Self {
        Self { storage }
    }

    /// Return a copied factor row.
    #[inline]
    #[must_use]
    const fn row(&self, index: usize) -> [f64; D] {
        self.storage.rows[index]
    }

    /// Return a factor entry.
    #[inline]
    #[must_use]
    const fn entry(&self, row: usize, col: usize) -> f64 {
        self.storage.rows[row][col]
    }

    /// Return a diagonal entry of `D`.
    #[inline]
    #[must_use]
    const fn diag(&self, index: usize) -> f64 {
        self.storage.rows[index][index]
    }
}

impl<const D: usize> Ldlt<D> {
    /// Factor a matrix that has already passed LDLT symmetry validation.
    #[inline]
    pub(crate) fn factor_symmetric(a: SymmetricMatrix<D>, tol: Tolerance) -> Result<Self, LaError> {
        let mut f = a.into_matrix();

        // LDLT via symmetric rank-1 updates, using only the lower triangle.
        for j in 0..D {
            let d = f.rows[j][j];
            if !d.is_finite() {
                cold_path();
                return Err(LaError::non_finite_cell(j, j));
            }
            if d < 0.0 {
                cold_path();
                return Err(LaError::not_positive_semidefinite(j, d));
            }
            if d <= tol.get() {
                cold_path();
                return Err(LaError::Singular { pivot_col: j });
            }

            // Compute L multipliers below the diagonal in column j.
            for i in (j + 1)..D {
                let l = f.rows[i][j] / d;
                if !l.is_finite() {
                    cold_path();
                    return Err(LaError::non_finite_cell(i, j));
                }
                f.rows[i][j] = l;
            }

            // Update the trailing submatrix (lower triangle): A := A - (L_col * d) * L_col^T.
            for i in (j + 1)..D {
                let l_i = f.rows[i][j];
                let l_i_d = l_i * d;

                for k in (j + 1)..=i {
                    let l_k = f.rows[k][j];
                    let new_val = (-l_i_d).mul_add(l_k, f.rows[i][k]);
                    if !new_val.is_finite() {
                        cold_path();
                        return Err(LaError::non_finite_cell(i, k));
                    }
                    f.rows[i][k] = new_val;
                }
            }
        }

        Ok(Self {
            factors: LdltFactors::new_unchecked(f),
        })
    }

    /// Determinant of the original matrix.
    ///
    /// For SPD/PSD matrices, this is the product of the diagonal terms of `D`.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// // Symmetric SPD matrix.
    /// let a = Matrix::<2>::try_from_rows([[4.0, 2.0], [2.0, 3.0]])?;
    /// let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL)?;
    ///
    /// assert!((ldlt.det()? - 8.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if the determinant product overflows to
    /// NaN or infinity.
    #[inline]
    pub const fn det(&self) -> Result<f64, LaError> {
        let mut det = 1.0;
        let mut i = 0;
        while i < D {
            det *= self.factors.diag(i);
            if !det.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            i += 1;
        }
        Ok(det)
    }

    /// Solve `A x = b` using this LDLT factorization.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::try_from_rows([[4.0, 2.0], [2.0, 3.0]])?;
    /// let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL)?;
    ///
    /// let b = Vector::<2>::try_new([1.0, 2.0])?;
    /// let x = ldlt.solve_vec(b)?.into_array();
    ///
    /// assert!((x[0] - (-0.125)).abs() <= 1e-12);
    /// assert!((x[1] - 0.75).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if the right-hand side contains
    /// NaN/infinity or a computed substitution intermediate overflows to NaN or
    /// infinity. The right-hand side is revalidated at this boundary before
    /// entering substitution; public callers normally encounter the same
    /// rejection earlier through [`Vector::try_new`](crate::Vector::try_new).
    #[inline]
    pub const fn solve_vec(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        let b = match FiniteVector::new(b) {
            Ok(b) => b,
            Err(err) => return Err(err),
        };
        match self.solve_finite_vec(b) {
            Ok(x) => Ok(x.into_vector()),
            Err(err) => Err(err),
        }
    }

    /// Solve `A x = b` using this LDLT factorization and a finite right-hand side.
    ///
    /// The right-hand side entries and stored factors are known finite, so this
    /// path only checks computed substitution overflows.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if a computed substitution intermediate
    /// overflows to NaN or infinity.
    #[inline]
    pub(crate) const fn solve_finite_vec(
        &self,
        b: FiniteVector<D>,
    ) -> Result<FiniteVector<D>, LaError> {
        let mut x = b.into_array();

        // Forward substitution: L y = b (L has unit diagonal).
        let mut i = 0;
        while i < D {
            let mut sum = x[i];
            let row = self.factors.row(i);
            let mut j = 0;
            while j < i {
                sum = (-row[j]).mul_add(x[j], sum);
                j += 1;
            }
            if !sum.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            x[i] = sum;
            i += 1;
        }

        // Diagonal solve: D z = y.
        let mut i = 0;
        while i < D {
            let diag = self.factors.diag(i);

            let quotient = x[i] / diag;
            if !quotient.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            x[i] = quotient;
            i += 1;
        }

        // Back substitution: Lᵀ x = z.
        let mut ii = 0;
        while ii < D {
            let i = D - 1 - ii;
            let mut sum = x[i];
            let mut j = i + 1;
            while j < D {
                sum = (-self.factors.entry(j, i)).mul_add(x[j], sum);
                j += 1;
            }
            if !sum.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            x[i] = sum;
            ii += 1;
        }

        Ok(FiniteVector::new_unchecked(Vector::new_unchecked(x)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::DEFAULT_SINGULAR_TOL;
    use crate::matrix::FiniteMatrix;

    use core::assert_matches;
    use core::hint::black_box;

    use approx::assert_abs_diff_eq;
    use pastey::paste;

    macro_rules! gen_public_api_ldlt_identity_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_api_ldlt_det_and_solve_identity_ $d d>]() {
                    let a = Matrix::<$d>::identity();
                    let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

                    assert_abs_diff_eq!(ldlt.det().unwrap(), 1.0, epsilon = 1e-12);

                    let b_arr = {
                        let mut arr = [0.0f64; $d];
                        let values = [1.0f64, 2.0, 3.0, 4.0, 5.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = *src;
                        }
                        arr
                    };
                    let b = Vector::<$d>::new(black_box(b_arr));
                    let x = ldlt.solve_vec(b).unwrap().into_array();

                    for i in 0..$d {
                        assert_abs_diff_eq!(x[i], b_arr[i], epsilon = 1e-12);
                    }
                }
            }
        };
    }

    gen_public_api_ldlt_identity_tests!(2);
    gen_public_api_ldlt_identity_tests!(3);
    gen_public_api_ldlt_identity_tests!(4);
    gen_public_api_ldlt_identity_tests!(5);

    macro_rules! gen_public_api_ldlt_diagonal_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_api_ldlt_det_and_solve_diagonal_spd_ $d d>]() {
                    let diag = {
                        let mut arr = [0.0f64; $d];
                        let values = [1.0f64, 2.0, 3.0, 4.0, 5.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = *src;
                        }
                        arr
                    };

                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = diag[i];
                    }

                    let a = Matrix::<$d>::from_rows(black_box(rows));
                    let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

                    let expected_det = {
                        let mut acc = 1.0;
                        for i in 0..$d {
                            acc *= diag[i];
                        }
                        acc
                    };
                    assert_abs_diff_eq!(ldlt.det().unwrap(), expected_det, epsilon = 1e-12);

                    let b_arr = {
                        let mut arr = [0.0f64; $d];
                        let values = [5.0f64, 4.0, 3.0, 2.0, 1.0];
                        for (dst, src) in arr.iter_mut().zip(values.iter()) {
                            *dst = *src;
                        }
                        arr
                    };

                    let b = Vector::<$d>::new(black_box(b_arr));
                    let x = ldlt.solve_vec(b).unwrap().into_array();

                    for i in 0..$d {
                        assert_abs_diff_eq!(x[i], b_arr[i] / diag[i], epsilon = 1e-12);
                    }
                }
            }
        };
    }

    gen_public_api_ldlt_diagonal_tests!(2);
    gen_public_api_ldlt_diagonal_tests!(3);
    gen_public_api_ldlt_diagonal_tests!(4);
    gen_public_api_ldlt_diagonal_tests!(5);

    #[test]
    fn solve_2x2_known_spd() {
        let a = Matrix::<2>::from_rows(black_box([[4.0, 2.0], [2.0, 3.0]]));
        let ldlt = FiniteMatrix::new(a)
            .unwrap()
            .ldlt(DEFAULT_SINGULAR_TOL)
            .unwrap();

        let b = Vector::<2>::new(black_box([1.0, 2.0]));
        let x = ldlt.solve_vec(b).unwrap().into_array();

        assert_abs_diff_eq!(x[0], -0.125, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], 0.75, epsilon = 1e-12);
        assert_abs_diff_eq!(ldlt.det().unwrap(), 8.0, epsilon = 1e-12);
    }

    #[test]
    fn solve_3x3_spd_tridiagonal_smoke() {
        let a = Matrix::<3>::from_rows(black_box([
            [2.0, -1.0, 0.0],
            [-1.0, 2.0, -1.0],
            [0.0, -1.0, 2.0],
        ]));
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        // Choose x = 1 so b = A x is simple: [1, 0, 1].
        let b = Vector::<3>::new(black_box([1.0, 0.0, 1.0]));
        let x = ldlt.solve_vec(b).unwrap().into_array();

        for &x_i in &x {
            assert_abs_diff_eq!(x_i, 1.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn singular_detected_for_degenerate_psd() {
        // Rank-1 Gram-like matrix.
        let a = Matrix::<2>::from_rows(black_box([[1.0, 1.0], [1.0, 1.0]]));
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(err, LaError::Singular { pivot_col: 1 });
    }

    #[test]
    fn negative_initial_diagonal_reports_not_positive_semidefinite() {
        let a = Matrix::<2>::from_rows(black_box([[-1.0, 0.0], [0.0, 1.0]]));
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::NotPositiveSemidefinite {
                pivot_col: 0,
                value: -1.0,
            }
        );
    }

    #[test]
    fn negative_updated_diagonal_reports_not_positive_semidefinite() {
        let a = Matrix::<2>::from_rows(black_box([[1.0, 2.0], [2.0, 1.0]]));
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::NotPositiveSemidefinite {
                pivot_col: 1,
                value: -3.0,
            }
        );
    }

    #[test]
    fn matrix_constructor_rejects_nonfinite_diagonal() {
        let err = Matrix::<2>::try_from_rows([[f64::NAN, 0.0], [0.0, 1.0]]).unwrap_err();
        assert_eq!(
            err,
            LaError::NonFinite {
                row: Some(0),
                col: 0
            }
        );
    }

    #[test]
    fn matrix_constructor_rejects_nonfinite_offdiagonal_before_asymmetry() {
        let err = Matrix::<2>::try_from_rows([[1.0, f64::NAN], [0.0, 1.0]]).unwrap_err();
        assert_eq!(
            err,
            LaError::NonFinite {
                row: Some(0),
                col: 1,
            }
        );
    }

    #[test]
    fn nonfinite_l_multiplier_overflow() {
        // d = 1e-11 > tol, but l = 1e300 / 1e-11 = 1e311 overflows f64.
        let a = Matrix::<2>::from_rows([[1e-11, 1e300], [1e300, 1.0]]);
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::NonFinite {
                row: Some(1),
                col: 0
            }
        );
    }

    #[test]
    fn nonfinite_trailing_submatrix_overflow() {
        // L multiplier is finite (1e200), but the rank-1 update
        // (-1e200 * 1.0) * 1e200 + 1.0 overflows.
        let a = Matrix::<2>::from_rows([[1.0, 1e200], [1e200, 1.0]]);
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::NonFinite {
                row: Some(1),
                col: 1
            }
        );
    }

    #[test]
    fn nonfinite_solve_vec_forward_substitution_overflow() {
        // SPD matrix with large L multiplier: L[1,0] = 1e153.
        // Forward substitution overflows: y[1] = 0 - 1e153 * 1e156 = -inf.
        let a = Matrix::<3>::from_rows([
            [1.0, 1e153, 0.0],
            [1e153, 1e306 + 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<3>::new([1e156, 0.0, 0.0]);
        let err = ldlt.solve_vec(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { row: None, col: 1 });
    }

    #[test]
    fn nonfinite_solve_vec_back_substitution_overflow() {
        // SPD matrix: [[1,0,0],[0,1,2],[0,2,5]] has LDLT factors
        // D=[1,1,1], L[2,1]=2.  Forward sub and diagonal solve produce
        // z=[0,0,1e308].  Back-substitution: x[2]=1e308 then
        // x[1] = 0 - 2*1e308 = -inf (overflows f64).
        let a = Matrix::<3>::from_rows([[1.0, 0.0, 0.0], [0.0, 1.0, 2.0], [0.0, 2.0, 5.0]]);
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<3>::new([0.0, 0.0, 1e308]);
        let err = ldlt.solve_vec(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { row: None, col: 1 });
    }

    #[test]
    fn nonfinite_solve_vec_diagonal_solve_overflow() {
        // Diagonal SPD matrix with a tiny diagonal entry just above the
        // singularity tolerance.  Forward substitution passes through the
        // large RHS unchanged, then the diagonal solve z[1] = y[1] / D[1]
        // = 1e300 / 1e-11 = 1e311 overflows f64, exercising the
        // `!v.is_finite()` branch of the diagonal solve.
        let a = Matrix::<2>::from_rows([[1.0, 0.0], [0.0, 1.0e-11]]);
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<2>::new([0.0, 1.0e300]);
        let err = ldlt.solve_vec(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { row: None, col: 1 });
    }

    #[test]
    fn det_rejects_product_overflow() {
        let a = Matrix::<5>::from_rows([
            [1.0e100, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0e100, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0e100, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0e100, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0e100],
        ]);
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();
        assert_eq!(ldlt.det(), Err(LaError::NonFinite { row: None, col: 3 }));
    }

    #[test]
    fn asymmetric_input_returns_typed_error() {
        // a[0][1] = 2.0 but a[1][0] = -2.0 → clearly asymmetric.
        let a = Matrix::<3>::from_rows([[4.0, 2.0, 0.0], [-2.0, 5.0, 1.0], [0.0, 1.0, 3.0]]);
        assert_eq!(
            a.ldlt(DEFAULT_SINGULAR_TOL),
            Err(LaError::Asymmetric {
                row: 0,
                col: 1,
                dim: 3,
            })
        );
    }

    #[test]
    fn invalid_tolerance_rejected() {
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

    macro_rules! gen_solve_vec_boundary_tests {
        ($d:literal) => {
            paste! {
                /// Raw non-finite right-hand sides are rejected before a
                /// public caller can construct a `Vector`.
                #[test]
                fn [<solve_vec_rhs_boundary_rejects_non_finite_ $d d>]() {
                    let mut rhs = [1.0; $d];
                    rhs[$d - 1] = f64::NAN;

                    assert_eq!(
                        Vector::<$d>::try_new(rhs),
                        Err(LaError::NonFinite {
                            row: None,
                            col: $d - 1,
                        })
                    );
                }

                #[test]
                fn [<solve_vec_revalidates_unchecked_rhs_storage_ $d d>]() {
                    let ldlt = Matrix::<$d>::identity().ldlt(DEFAULT_SINGULAR_TOL).unwrap();
                    let mut rhs = [1.0; $d];
                    rhs[$d - 1] = f64::NAN;
                    let rhs = Vector::<$d>::new_unchecked(rhs);

                    assert_eq!(
                        ldlt.solve_vec(rhs),
                        Err(LaError::NonFinite {
                            row: None,
                            col: $d - 1,
                        })
                    );
                }
            }
        };
    }

    gen_solve_vec_boundary_tests!(2);
    gen_solve_vec_boundary_tests!(3);
    gen_solve_vec_boundary_tests!(4);
    gen_solve_vec_boundary_tests!(5);

    // -----------------------------------------------------------------------
    // Const-evaluability tests.
    //
    // These prove that `Ldlt::det` and `Ldlt::solve_vec` are truly `const fn`
    // by forcing the compiler to evaluate them inside a `const` initializer.
    // `Ldlt::factor` is not (yet) `const fn` because the rank-1 update loop
    // uses array indexing patterns that still require non-const helpers on
    // some toolchains; we therefore construct `Ldlt<D>` directly.
    // -----------------------------------------------------------------------

    macro_rules! gen_ldlt_const_eval_tests {
        ($d:literal) => {
            paste! {
                /// `Ldlt::det` must be fully const-evaluable. Setting
                /// `factors[0][0] = 2.0` and leaving the remaining identity
                /// diagonals at `1.0` gives `det = 2.0` for every `D ≥ 1`,
                /// exercising the multiply-accumulate loop at each dimension.
                #[test]
                fn [<ldlt_det_const_eval_ $d d>]() {
                    const DET: Result<f64, LaError> = {
                        let mut factors = Matrix::<$d>::identity();
                        factors.rows[0][0] = 2.0;
                        let ldlt = Ldlt::<$d> {
                            factors: LdltFactors::new_unchecked(factors),
                        };
                        ldlt.det()
                    };
                    assert_eq!(DET, Ok(2.0));
                }

                /// `Ldlt::solve_vec` must be fully const-evaluable. Identity
                /// factors with RHS `b = [1.0, 2.0, …, D]` round-trips `b`
                /// unchanged, exercising the full forward sub / diagonal solve
                /// / back sub pipeline inside a `const { … }` initializer.
                #[test]
                fn [<ldlt_solve_vec_const_eval_ $d d>]() {
                    #[allow(clippy::cast_precision_loss)]
                        const X: [f64; $d] = {
                            let ldlt = Ldlt::<$d> {
                                factors: LdltFactors::new_unchecked(Matrix::<$d>::identity()),
                            };
                        let mut b_arr = [0.0f64; $d];
                        let mut i = 0;
                        while i < $d {
                            b_arr[i] = i as f64 + 1.0;
                            i += 1;
                        }
                        let b = Vector::<$d>::new(b_arr);
                        match ldlt.solve_vec(b) {
                            Ok(v) => v.into_array(),
                            Err(_) => [0.0f64; $d],
                        }
                    };
                    #[allow(clippy::cast_precision_loss)]
                    for i in 0..$d {
                        let expected = i as f64 + 1.0;
                        assert!((X[i] - expected).abs() <= 1e-12);
                    }
                }
            }
        };
    }

    gen_ldlt_const_eval_tests!(2);
    gen_ldlt_const_eval_tests!(3);
    gen_ldlt_const_eval_tests!(4);
    gen_ldlt_const_eval_tests!(5);
}
