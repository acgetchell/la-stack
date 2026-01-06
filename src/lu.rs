//! LU decomposition and solves.

use crate::LaError;
use crate::matrix::Matrix;
use crate::vector::Vector;

/// LU decomposition (PA = LU) with partial pivoting.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lu<const D: usize> {
    factors: Matrix<D>,
    piv: [usize; D],
    piv_sign: f64,
    tol: f64,
}

impl<const D: usize> Lu<D> {
    #[inline]
    pub(crate) fn factor(a: Matrix<D>, tol: f64) -> Result<Self, LaError> {
        let mut lu = a;

        let mut piv = [0usize; D];
        for (i, p) in piv.iter_mut().enumerate() {
            *p = i;
        }

        let mut piv_sign = 1.0;

        for k in 0..D {
            // Choose pivot row.
            let mut pivot_row = k;
            let mut pivot_abs = lu.rows[k][k].abs();
            if !pivot_abs.is_finite() {
                return Err(LaError::NonFinite { pivot_col: k });
            }

            for r in (k + 1)..D {
                let v = lu.rows[r][k].abs();
                if !v.is_finite() {
                    return Err(LaError::NonFinite { pivot_col: k });
                }
                if v > pivot_abs {
                    pivot_abs = v;
                    pivot_row = r;
                }
            }

            if pivot_abs <= tol {
                return Err(LaError::Singular { pivot_col: k });
            }

            if pivot_row != k {
                lu.rows.swap(k, pivot_row);
                piv.swap(k, pivot_row);
                piv_sign = -piv_sign;
            }

            let pivot = lu.rows[k][k];
            if !pivot.is_finite() {
                return Err(LaError::NonFinite { pivot_col: k });
            }

            // Eliminate below pivot.
            for r in (k + 1)..D {
                let mult = lu.rows[r][k] / pivot;
                if !mult.is_finite() {
                    return Err(LaError::NonFinite { pivot_col: k });
                }
                lu.rows[r][k] = mult;

                for c in (k + 1)..D {
                    lu.rows[r][c] = (-mult).mul_add(lu.rows[k][c], lu.rows[r][c]);
                }
            }
        }

        Ok(Self {
            factors: lu,
            piv,
            piv_sign,
            tol,
        })
    }

    /// Solve `A x = b` using this LU factorization.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let lu = a.lu(DEFAULT_PIVOT_TOL)?;
    ///
    /// let b = Vector::<2>::new([5.0, 11.0]);
    /// let x = lu.solve_vec(b)?.into_array();
    ///
    /// assert!((x[0] - 1.0).abs() <= 1e-12);
    /// assert!((x[1] - 2.0).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if a diagonal entry of `U` satisfies `|u_ii| <= tol`, where
    /// `tol` is the tolerance that was used during factorization.
    /// Returns [`LaError::NonFinite`] if NaN/∞ is detected.
    #[inline]
    pub fn solve_vec(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        let mut x = [0.0; D];
        for (i, x_i) in x.iter_mut().enumerate() {
            *x_i = b.data[self.piv[i]];
        }

        // Forward substitution for L (unit diagonal).
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

        // Back substitution for U.
        for ii in 0..D {
            let i = D - 1 - ii;
            let mut sum = x[i];
            let row = self.factors.rows[i];
            for (j, x_j) in x.iter().enumerate().skip(i + 1) {
                sum = (-row[j]).mul_add(*x_j, sum);
            }

            let diag = row[i];
            if !diag.is_finite() || !sum.is_finite() {
                return Err(LaError::NonFinite { pivot_col: i });
            }
            if diag.abs() <= self.tol {
                return Err(LaError::Singular { pivot_col: i });
            }

            x[i] = sum / diag;
        }

        Ok(Vector::new(x))
    }

    /// Determinant of the original matrix.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
    /// let lu = a.lu(DEFAULT_PIVOT_TOL)?;
    ///
    /// let det = lu.det();
    /// assert!((det - (-2.0)).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    #[must_use]
    pub fn det(&self) -> f64 {
        let mut det = self.piv_sign;
        for i in 0..D {
            det *= self.factors.rows[i][i];
        }
        det
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_PIVOT_TOL;

    use core::hint::black_box;

    use approx::assert_abs_diff_eq;
    use pastey::paste;

    macro_rules! gen_public_api_pivoting_solve_vec_and_det_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_api_lu_solve_vec_pivoting_ $d d>]() {
                    // Public API path under test:
                    // Matrix::lu (pub) -> Lu::solve_vec (pub).

                    // Permutation matrix that swaps the first two basis vectors.
                    // This forces pivoting in column 0 for any D >= 2.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 1.0;
                    }
                    rows.swap(0, 1);

                    let a = Matrix::<$d>::from_rows(black_box(rows));
                    let lu_fn: fn(Matrix<$d>, f64) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_PIVOT_TOL).unwrap();

                    // Pick a simple RHS with unique entries, so the expected swap is obvious.
                    let b_arr = {
                        let mut arr = [0.0f64; $d];
                        let mut val = 1.0f64;
                        for dst in arr.iter_mut() {
                            *dst = val;
                            val += 1.0;
                        }
                        arr
                    };
                    let mut expected = b_arr;
                    expected.swap(0, 1);
                    let b = Vector::<$d>::new(black_box(b_arr));

                    let solve_fn: fn(&Lu<$d>, Vector<$d>) -> Result<Vector<$d>, LaError> =
                        black_box(Lu::<$d>::solve_vec);
                    let x = solve_fn(&lu, b).unwrap().into_array();

                    for i in 0..$d {
                        assert_abs_diff_eq!(x[i], expected[i], epsilon = 1e-12);
                    }
                }

                #[test]
                fn [<public_api_lu_det_pivoting_ $d d>]() {
                    // Public API path under test:
                    // Matrix::lu (pub) -> Lu::det (pub).

                    // Permutation matrix that swaps the first two basis vectors.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 1.0;
                    }
                    rows.swap(0, 1);

                    let a = Matrix::<$d>::from_rows(black_box(rows));
                    let lu_fn: fn(Matrix<$d>, f64) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_PIVOT_TOL).unwrap();

                    // Row swap ⇒ determinant sign flip.
                    let det_fn: fn(&Lu<$d>) -> f64 = black_box(Lu::<$d>::det);
                    assert_abs_diff_eq!(det_fn(&lu), -1.0, epsilon = 1e-12);
                }
            }
        };
    }

    gen_public_api_pivoting_solve_vec_and_det_tests!(2);
    gen_public_api_pivoting_solve_vec_and_det_tests!(3);
    gen_public_api_pivoting_solve_vec_and_det_tests!(4);
    gen_public_api_pivoting_solve_vec_and_det_tests!(5);

    macro_rules! gen_public_api_tridiagonal_smoke_solve_vec_and_det_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_api_lu_solve_vec_tridiagonal_smoke_ $d d>]() {
                    // Public API path under test:
                    // Matrix::lu (pub) -> Lu::solve_vec (pub).

                    // Classic SPD tridiagonal: 2 on diagonal, -1 on sub/super-diagonals.
                    #[allow(clippy::large_stack_arrays)]
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 2.0;
                        if i > 0 {
                            rows[i][i - 1] = -1.0;
                        }
                        if i + 1 < $d {
                            rows[i][i + 1] = -1.0;
                        }
                    }

                    let a = Matrix::<$d>::from_rows(black_box(rows));
                    let lu_fn: fn(Matrix<$d>, f64) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_PIVOT_TOL).unwrap();

                    // Choose x = 1, so b = A x is simple: [1, 0, 0, ..., 0, 1].
                    let mut b_arr = [0.0f64; $d];
                    b_arr[0] = 1.0;
                    b_arr[$d - 1] = 1.0;
                    let b = Vector::<$d>::new(black_box(b_arr));

                    let solve_fn: fn(&Lu<$d>, Vector<$d>) -> Result<Vector<$d>, LaError> =
                        black_box(Lu::<$d>::solve_vec);
                    let x = solve_fn(&lu, b).unwrap().into_array();

                    for &x_i in &x {
                        assert_abs_diff_eq!(x_i, 1.0, epsilon = 1e-9);
                    }
                }

                #[test]
                fn [<public_api_lu_det_tridiagonal_smoke_ $d d>]() {
                    // Public API path under test:
                    // Matrix::lu (pub) -> Lu::det (pub).

                    // Classic SPD tridiagonal: 2 on diagonal, -1 on sub/super-diagonals.
                    // Determinant is known exactly: det = D + 1.
                    #[allow(clippy::large_stack_arrays)]
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 2.0;
                        if i > 0 {
                            rows[i][i - 1] = -1.0;
                        }
                        if i + 1 < $d {
                            rows[i][i + 1] = -1.0;
                        }
                    }

                    let a = Matrix::<$d>::from_rows(black_box(rows));
                    let lu_fn: fn(Matrix<$d>, f64) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_PIVOT_TOL).unwrap();

                    let det_fn: fn(&Lu<$d>) -> f64 = black_box(Lu::<$d>::det);
                    assert_abs_diff_eq!(det_fn(&lu), f64::from($d) + 1.0, epsilon = 1e-8);
                }
            }
        };
    }

    gen_public_api_tridiagonal_smoke_solve_vec_and_det_tests!(16);
    gen_public_api_tridiagonal_smoke_solve_vec_and_det_tests!(32);
    gen_public_api_tridiagonal_smoke_solve_vec_and_det_tests!(64);

    #[test]
    fn solve_1x1() {
        let a = Matrix::<1>::from_rows(black_box([[2.0]]));
        let lu = (black_box(Lu::<1>::factor))(a, DEFAULT_PIVOT_TOL).unwrap();

        let b = Vector::<1>::new(black_box([6.0]));
        let solve_fn: fn(&Lu<1>, Vector<1>) -> Result<Vector<1>, LaError> =
            black_box(Lu::<1>::solve_vec);
        let x = solve_fn(&lu, b).unwrap().into_array();
        assert_abs_diff_eq!(x[0], 3.0, epsilon = 1e-12);

        let det_fn: fn(&Lu<1>) -> f64 = black_box(Lu::<1>::det);
        assert_abs_diff_eq!(det_fn(&lu), 2.0, epsilon = 0.0);
    }

    #[test]
    fn solve_2x2_basic() {
        let a = Matrix::<2>::from_rows(black_box([[1.0, 2.0], [3.0, 4.0]]));
        let lu = (black_box(Lu::<2>::factor))(a, DEFAULT_PIVOT_TOL).unwrap();
        let b = Vector::<2>::new(black_box([5.0, 11.0]));

        let solve_fn: fn(&Lu<2>, Vector<2>) -> Result<Vector<2>, LaError> =
            black_box(Lu::<2>::solve_vec);
        let x = solve_fn(&lu, b).unwrap().into_array();

        assert_abs_diff_eq!(x[0], 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], 2.0, epsilon = 1e-12);
    }

    #[test]
    fn det_2x2_basic() {
        let a = Matrix::<2>::from_rows(black_box([[1.0, 2.0], [3.0, 4.0]]));
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

        let det_fn: fn(&Lu<2>) -> f64 = black_box(Lu::<2>::det);
        assert_abs_diff_eq!(det_fn(&lu), -2.0, epsilon = 1e-12);
    }

    #[test]
    fn det_requires_pivot_sign() {
        // Row swap ⇒ determinant sign flip.
        let a = Matrix::<2>::from_rows(black_box([[0.0, 1.0], [1.0, 0.0]]));
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

        let det_fn: fn(&Lu<2>) -> f64 = black_box(Lu::<2>::det);
        assert_abs_diff_eq!(det_fn(&lu), -1.0, epsilon = 0.0);
    }

    #[test]
    fn solve_requires_pivoting() {
        let a = Matrix::<2>::from_rows(black_box([[0.0, 1.0], [1.0, 0.0]]));
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();
        let b = Vector::<2>::new(black_box([1.0, 2.0]));

        let solve_fn: fn(&Lu<2>, Vector<2>) -> Result<Vector<2>, LaError> =
            black_box(Lu::<2>::solve_vec);
        let x = solve_fn(&lu, b).unwrap().into_array();

        assert_abs_diff_eq!(x[0], 2.0, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn singular_detected() {
        let a = Matrix::<2>::from_rows(black_box([[1.0, 2.0], [2.0, 4.0]]));
        let err = a.lu(DEFAULT_PIVOT_TOL).unwrap_err();
        assert_eq!(err, LaError::Singular { pivot_col: 1 });
    }

    #[test]
    fn singular_due_to_tolerance_at_first_pivot() {
        // Not exactly singular, but below DEFAULT_PIVOT_TOL.
        let a = Matrix::<2>::from_rows(black_box([[1e-13, 0.0], [0.0, 1.0]]));
        let err = a.lu(DEFAULT_PIVOT_TOL).unwrap_err();
        assert_eq!(err, LaError::Singular { pivot_col: 0 });
    }

    #[test]
    fn nonfinite_detected_on_pivot_entry() {
        let a = Matrix::<2>::from_rows([[f64::NAN, 0.0], [0.0, 1.0]]);
        let err = a.lu(DEFAULT_PIVOT_TOL).unwrap_err();
        assert_eq!(err, LaError::NonFinite { pivot_col: 0 });
    }

    #[test]
    fn nonfinite_detected_in_pivot_column_scan() {
        let a = Matrix::<2>::from_rows([[1.0, 0.0], [f64::INFINITY, 1.0]]);
        let err = a.lu(DEFAULT_PIVOT_TOL).unwrap_err();
        assert_eq!(err, LaError::NonFinite { pivot_col: 0 });
    }

    #[test]
    fn solve_vec_nonfinite_forward_substitution_overflow() {
        // L has a -1 multiplier, and a large RHS makes forward substitution overflow.
        let a = Matrix::<3>::from_rows([[1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

        let b = Vector::<3>::new([1.0e308, 1.0e308, 0.0]);
        let err = lu.solve_vec(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { pivot_col: 1 });
    }

    #[test]
    fn solve_vec_nonfinite_back_substitution_overflow() {
        // Make x[1] overflow during back substitution, then ensure it is detected on the next row.
        let a = Matrix::<2>::from_rows([[1.0, 1.0], [0.0, 2.0e-12]]);
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

        let b = Vector::<2>::new([0.0, 1.0e300]);
        let err = lu.solve_vec(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { pivot_col: 0 });
    }
}
