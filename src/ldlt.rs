//! LDLT factorization and solves.
//!
//! This module provides a stack-allocated LDLT factorization (`A = L D Lᵀ`) intended for
//! symmetric positive definite (SPD) and positive semi-definite (PSD) matrices (e.g. Gram
//! matrices) without pivoting.

use crate::LaError;
use crate::matrix::Matrix;
use crate::vector::Vector;

/// LDLT factorization (`A = L D Lᵀ`) for symmetric positive (semi)definite matrices.
///
/// This factorization is **not** a general-purpose symmetric-indefinite LDLT (no pivoting).
/// It assumes the input matrix is symmetric and (numerically) SPD/PSD.
///
/// # Storage
/// The factors are stored in a single [`Matrix`]:
/// - `D` is stored on the diagonal.
/// - The strict lower triangle stores the multipliers of `L`.
/// - The diagonal of `L` is implicit ones.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ldlt<const D: usize> {
    factors: Matrix<D>,
    tol: f64,
}

impl<const D: usize> Ldlt<D> {
    #[inline]
    pub(crate) fn factor(a: Matrix<D>, tol: f64) -> Result<Self, LaError> {
        debug_assert!(tol >= 0.0, "tol must be non-negative");

        #[cfg(debug_assertions)]
        debug_assert_symmetric(&a);

        let mut f = a;

        // LDLT via symmetric rank-1 updates, using only the lower triangle.
        for j in 0..D {
            let d = f.rows[j][j];
            if !d.is_finite() {
                return Err(LaError::NonFinite { pivot_col: j });
            }
            if d <= tol {
                return Err(LaError::Singular { pivot_col: j });
            }

            // Compute L multipliers below the diagonal in column j.
            for i in (j + 1)..D {
                let l = f.rows[i][j] / d;
                if !l.is_finite() {
                    return Err(LaError::NonFinite { pivot_col: j });
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
                        return Err(LaError::NonFinite { pivot_col: j });
                    }
                    f.rows[i][k] = new_val;
                }
            }
        }

        Ok(Self { factors: f, tol })
    }

    /// Determinant of the original matrix.
    ///
    /// For SPD/PSD matrices, this is the product of the diagonal terms of `D`.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// // Symmetric SPD matrix.
    /// let a = Matrix::<2>::from_rows([[4.0, 2.0], [2.0, 3.0]]);
    /// let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();
    ///
    /// assert!((ldlt.det() - 8.0).abs() <= 1e-12);
    /// ```
    #[inline]
    #[must_use]
    pub fn det(&self) -> f64 {
        let mut det = 1.0;
        for i in 0..D {
            det *= self.factors.rows[i][i];
        }
        det
    }

    /// Solve `A x = b` using this LDLT factorization.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::from_rows([[4.0, 2.0], [2.0, 3.0]]);
    /// let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL)?;
    ///
    /// let b = Vector::<2>::new([1.0, 2.0]);
    /// let x = ldlt.solve_vec(b)?.into_array();
    ///
    /// assert!((x[0] - (-0.125)).abs() <= 1e-12);
    /// assert!((x[1] - 0.75).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if a diagonal entry `d = D[i,i]` satisfies `d <= tol`
    /// (non-positive or too small), where `tol` is the tolerance that was used during factorization.
    /// Returns [`LaError::NonFinite`] if NaN/∞ is detected.
    #[inline]
    pub fn solve_vec(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        let mut x = b.data;

        // Forward substitution: L y = b (L has unit diagonal).
        for i in 0..D {
            let mut sum = x[i];
            let row = self.factors.rows[i];
            for (j, x_j) in x.iter().enumerate().take(i) {
                sum = (-row[j]).mul_add(*x_j, sum);
            }
            if !sum.is_finite() {
                return Err(LaError::NonFinite { pivot_col: i });
            }
            x[i] = sum;
        }

        // Diagonal solve: D z = y.
        for (i, x_i) in x.iter_mut().enumerate().take(D) {
            let diag = self.factors.rows[i][i];
            if !diag.is_finite() {
                return Err(LaError::NonFinite { pivot_col: i });
            }
            if diag <= self.tol {
                return Err(LaError::Singular { pivot_col: i });
            }

            let v = *x_i / diag;
            if !v.is_finite() {
                return Err(LaError::NonFinite { pivot_col: i });
            }
            *x_i = v;
        }

        // Back substitution: Lᵀ x = z.
        for ii in 0..D {
            let i = D - 1 - ii;
            let mut sum = x[i];
            for (j, x_j) in x.iter().enumerate().skip(i + 1) {
                sum = (-self.factors.rows[j][i]).mul_add(*x_j, sum);
            }
            if !sum.is_finite() {
                return Err(LaError::NonFinite { pivot_col: i });
            }
            x[i] = sum;
        }

        Ok(Vector::new(x))
    }
}

#[cfg(debug_assertions)]
fn debug_assert_symmetric<const D: usize>(a: &Matrix<D>) {
    let scale = a.inf_norm().max(1.0);
    let eps = 1e-12 * scale;

    for r in 0..D {
        for c in (r + 1)..D {
            let diff = (a.rows[r][c] - a.rows[c][r]).abs();
            debug_assert!(
                diff <= eps,
                "matrix must be symmetric (diff={diff}, eps={eps}) at ({r}, {c})"
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::DEFAULT_SINGULAR_TOL;

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

                    assert_abs_diff_eq!(ldlt.det(), 1.0, epsilon = 1e-12);

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
                    assert_abs_diff_eq!(ldlt.det(), expected_det, epsilon = 1e-12);

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
        let ldlt = (black_box(Ldlt::<2>::factor))(a, DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<2>::new(black_box([1.0, 2.0]));
        let x = ldlt.solve_vec(b).unwrap().into_array();

        assert_abs_diff_eq!(x[0], -0.125, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], 0.75, epsilon = 1e-12);
        assert_abs_diff_eq!(ldlt.det(), 8.0, epsilon = 1e-12);
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
    fn nonfinite_detected() {
        let a = Matrix::<2>::from_rows([[f64::NAN, 0.0], [0.0, 1.0]]);
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(err, LaError::NonFinite { pivot_col: 0 });
    }
}
