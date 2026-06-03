//! LU decomposition and solves.

use core::hint::cold_path;

use crate::matrix::{FiniteMatrix, Matrix};
use crate::vector::{FiniteVector, Vector};
use crate::{LaError, Tolerance};

/// LU decomposition (PA = LU) with partial pivoting.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lu<const D: usize> {
    factors: Matrix<D>,
    piv: [usize; D],
    piv_sign: f64,
    tol: Tolerance,
}

impl<const D: usize> Lu<D> {
    /// Factor a finite square matrix into in-place LU storage for
    /// `FiniteMatrix::lu`.
    ///
    /// The input has already proven finite entries, so LU construction rejects
    /// numerically singular pivots and non-finite elimination intermediates
    /// before callers can observe a [`Lu`] value. Computed trailing updates are
    /// checked before storage so successful factors do not contain a non-finite
    /// value produced during elimination.
    #[inline]
    pub(crate) fn factor_finite(a: FiniteMatrix<D>, tol: Tolerance) -> Result<Self, LaError> {
        let mut lu = a.into_matrix();

        let mut piv = [0usize; D];
        for (i, p) in piv.iter_mut().enumerate() {
            *p = i;
        }

        let mut piv_sign = 1.0;

        for k in 0..D {
            // Choose pivot row.
            let mut pivot_row = k;
            let mut pivot_abs = lu.rows[k][k].abs();

            for r in (k + 1)..D {
                let v = lu.rows[r][k].abs();
                if v > pivot_abs {
                    pivot_abs = v;
                    pivot_row = r;
                }
            }

            if pivot_abs <= tol.get() {
                cold_path();
                return Err(LaError::Singular { pivot_col: k });
            }

            if pivot_row != k {
                lu.rows.swap(k, pivot_row);
                piv.swap(k, pivot_row);
                piv_sign = -piv_sign;
            }

            let pivot = lu.rows[k][k];

            // Eliminate below pivot.
            for r in (k + 1)..D {
                let mult = lu.rows[r][k] / pivot;
                if !mult.is_finite() {
                    cold_path();
                    return Err(LaError::non_finite_cell(r, k));
                }
                lu.rows[r][k] = mult;

                for c in (k + 1)..D {
                    let updated = (-mult).mul_add(lu.rows[k][c], lu.rows[r][c]);
                    if !updated.is_finite() {
                        cold_path();
                        return Err(LaError::non_finite_cell(r, c));
                    }
                    lu.rows[r][c] = updated;
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
    ///
    /// Returns [`LaError::NonFinite`] if NaN/∞ is detected. The `row`/`col` coordinates
    /// follow the convention documented on [`LaError::NonFinite`]:
    ///
    /// - `row: Some(i), col: i` — the stored `U` diagonal at `(i, i)` is non-finite
    ///   (only reachable via direct `Lu` construction; [`Matrix::lu`](crate::Matrix::lu)
    ///   rejects such factorizations).
    /// - `row: None, col: i` — a computed intermediate (forward/back-substitution
    ///   accumulator or the quotient `sum / diag`) overflowed to NaN/∞ at step `i`.
    #[inline]
    pub const fn solve_vec(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        match FiniteVector::try_new(b) {
            Ok(finite) => match self.solve_finite_vec(finite) {
                Ok(x) => Ok(x.into_vector()),
                Err(err) => Err(err),
            },
            Err(err) => Err(err),
        }
    }

    /// Solve `A x = b` using this LU factorization and a finite right-hand side.
    ///
    /// The right-hand side entries are known finite, so this path only checks
    /// factorization defensive invariants and computed substitution overflows.
    ///
    /// # Errors
    /// Returns [`LaError::Singular`] if a diagonal entry of `U` satisfies
    /// `|u_ii| <= tol`, where `tol` is the tolerance that was used during
    /// factorization.
    ///
    /// Returns [`LaError::NonFinite`] if a stored factorization diagonal is
    /// corrupt or if a computed substitution intermediate overflows to NaN or
    /// infinity.
    #[inline]
    pub(crate) const fn solve_finite_vec(
        &self,
        b: FiniteVector<D>,
    ) -> Result<FiniteVector<D>, LaError> {
        let mut x = [0.0; D];
        let mut i = 0;
        let b = b.into_array();
        while i < D {
            x[i] = b[self.piv[i]];
            i += 1;
        }

        // Forward substitution for L (unit diagonal).
        let mut i = 0;
        while i < D {
            let mut sum = x[i];
            let row = self.factors.rows[i];
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

        // Back substitution for U.
        let mut ii = 0;
        while ii < D {
            let i = D - 1 - ii;
            let mut sum = x[i];
            let row = self.factors.rows[i];
            let mut j = i + 1;
            while j < D {
                sum = (-row[j]).mul_add(x[j], sum);
                j += 1;
            }

            let diag = row[i];
            // Distinguish a corrupt stored pivot (row: Some(i), col: i) from
            // a computed intermediate overflow (row: None, col: i) so callers
            // can diagnose the failure source without inspecting internals.
            if !diag.is_finite() {
                cold_path();
                return Err(LaError::non_finite_cell(i, i));
            }
            if !sum.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            if diag.abs() <= self.tol.get() {
                cold_path();
                return Err(LaError::Singular { pivot_col: i });
            }

            let quotient = sum / diag;
            if !quotient.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            x[i] = quotient;
            ii += 1;
        }

        Ok(FiniteVector::new_unchecked(Vector::new(x)))
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
    /// let det = lu.det()?;
    /// assert!((det - (-2.0)).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if the determinant product overflows to
    /// NaN or infinity.
    #[inline]
    pub const fn det(&self) -> Result<f64, LaError> {
        let mut det = self.piv_sign;
        let mut i = 0;
        while i < D {
            det *= self.factors.rows[i][i];
            if !det.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            i += 1;
        }
        Ok(det)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DEFAULT_PIVOT_TOL;

    use core::assert_matches;
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
                    let lu_fn: fn(Matrix<$d>, Tolerance) -> Result<Lu<$d>, LaError> =
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
                    let lu_fn: fn(Matrix<$d>, Tolerance) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_PIVOT_TOL).unwrap();

                    // Row swap ⇒ determinant sign flip.
                    let det_fn: fn(&Lu<$d>) -> Result<f64, LaError> =
                        black_box(Lu::<$d>::det);
                    assert_abs_diff_eq!(det_fn(&lu).unwrap(), -1.0, epsilon = 1e-12);
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
                    let lu_fn: fn(Matrix<$d>, Tolerance) -> Result<Lu<$d>, LaError> =
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
                    let lu_fn: fn(Matrix<$d>, Tolerance) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_PIVOT_TOL).unwrap();

                    let det_fn: fn(&Lu<$d>) -> Result<f64, LaError> =
                        black_box(Lu::<$d>::det);
                    assert_abs_diff_eq!(det_fn(&lu).unwrap(), f64::from($d) + 1.0, epsilon = 1e-8);
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
        let lu = FiniteMatrix::try_new(a)
            .unwrap()
            .lu(DEFAULT_PIVOT_TOL)
            .unwrap();

        let b = Vector::<1>::new(black_box([6.0]));
        let solve_fn: fn(&Lu<1>, Vector<1>) -> Result<Vector<1>, LaError> =
            black_box(Lu::<1>::solve_vec);
        let x = solve_fn(&lu, b).unwrap().into_array();
        assert_abs_diff_eq!(x[0], 3.0, epsilon = 1e-12);

        let det_fn: fn(&Lu<1>) -> Result<f64, LaError> = black_box(Lu::<1>::det);
        assert_abs_diff_eq!(det_fn(&lu).unwrap(), 2.0, epsilon = 0.0);
    }

    #[test]
    fn solve_2x2_basic() {
        let a = Matrix::<2>::from_rows(black_box([[1.0, 2.0], [3.0, 4.0]]));
        let lu = FiniteMatrix::try_new(a)
            .unwrap()
            .lu(DEFAULT_PIVOT_TOL)
            .unwrap();
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

        let det_fn: fn(&Lu<2>) -> Result<f64, LaError> = black_box(Lu::<2>::det);
        assert_abs_diff_eq!(det_fn(&lu).unwrap(), -2.0, epsilon = 1e-12);
    }

    #[test]
    fn det_requires_pivot_sign() {
        // Row swap ⇒ determinant sign flip.
        let a = Matrix::<2>::from_rows(black_box([[0.0, 1.0], [1.0, 0.0]]));
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

        let det_fn: fn(&Lu<2>) -> Result<f64, LaError> = black_box(Lu::<2>::det);
        assert_abs_diff_eq!(det_fn(&lu).unwrap(), -1.0, epsilon = 0.0);
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

    #[test]
    fn nonfinite_detected_on_pivot_entry() {
        let a = Matrix::<2>::from_rows([[f64::NAN, 0.0], [0.0, 1.0]]);
        let err = a.lu(DEFAULT_PIVOT_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::NonFinite {
                row: Some(0),
                col: 0
            }
        );
    }

    #[test]
    fn nonfinite_detected_in_pivot_column_scan() {
        let a = Matrix::<2>::from_rows([[1.0, 0.0], [f64::INFINITY, 1.0]]);
        let err = a.lu(DEFAULT_PIVOT_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::NonFinite {
                row: Some(1),
                col: 0
            }
        );
    }

    #[test]
    fn nonfinite_detected_in_trailing_update() {
        let a =
            Matrix::<3>::from_rows([[1.0, f64::MAX, 0.0], [-1.0, f64::MAX, 0.0], [0.0, 0.0, 1.0]]);

        let err = a.lu(DEFAULT_PIVOT_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::NonFinite {
                row: Some(1),
                col: 1,
            }
        );
    }

    #[test]
    fn solve_vec_nonfinite_forward_substitution_overflow() {
        // L has a -1 multiplier, and a large RHS makes forward substitution overflow.
        let a = Matrix::<3>::from_rows([[1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

        let b = Vector::<3>::new([1.0e308, 1.0e308, 0.0]);
        let err = lu.solve_vec(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { row: None, col: 1 });
    }

    #[test]
    fn solve_vec_nonfinite_back_substitution_overflow() {
        // Make x[1] overflow during back substitution, then ensure it is detected on the next row.
        let a = Matrix::<2>::from_rows([[1.0, 1.0], [0.0, 2.0e-12]]);
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

        let b = Vector::<2>::new([0.0, 1.0e300]);
        let err = lu.solve_vec(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { row: None, col: 1 });
    }

    #[test]
    fn solve_vec_nonfinite_back_substitution_sum_overflow() {
        // Upper-triangular U with a very large off-diagonal in row 1 and a
        // very large x[2] produced by the RHS.  The back-substitution
        // accumulator `sum = (-row[j]).mul_add(x[j], sum)` overflows while
        // reducing row 1, so the failure is detected via the `!sum.is_finite()`
        // branch of the combined diag/sum check (distinct from the
        // `q = sum / diag` overflow path covered above).
        let a = Matrix::<3>::from_rows([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0e200], [0.0, 0.0, 1.0]]);
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

        let b = Vector::<3>::new([0.0, 0.0, 1.0e200]);
        let err = lu.solve_vec(b).unwrap_err();
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
        let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();
        assert_eq!(lu.det(), Err(LaError::NonFinite { row: None, col: 3 }));
    }

    // -----------------------------------------------------------------------
    // Defensive-path coverage for `solve_vec`.
    //
    // `Lu::factor` guarantees that every stored U diagonal is finite and
    // satisfies `|U[i,i]| > tol`.  `solve_vec` still re-checks during
    // back-substitution as a safety net (see the `!diag.is_finite()` and
    // `diag.abs() <= self.tol` guards).  Those branches are unreachable
    // through the public API, so the only way to exercise them is to
    // construct `Lu` directly with a corrupt U.  The tests below document
    // and verify that the safety nets return the documented error variants
    // with coordinates that locate the offending stored cell.
    // -----------------------------------------------------------------------

    macro_rules! gen_solve_vec_defensive_tests {
        ($d:literal) => {
            paste! {
                /// `solve_vec` must surface `Singular` when a stored U
                /// diagonal is at or below the recorded tolerance, even
                /// though `factor` cannot produce such a factorization.
                #[test]
                fn [<solve_vec_defensive_sub_tolerance_diagonal_ $d d>]() {
                    let mut factors = Matrix::<$d>::identity();
                    factors.rows[$d - 1][$d - 1] = 0.0;

                    let mut piv = [0usize; $d];
                    for (i, p) in piv.iter_mut().enumerate() {
                        *p = i;
                    }

                    let lu = Lu::<$d> {
                        factors,
                        piv,
                        piv_sign: 1.0,
                        tol: DEFAULT_PIVOT_TOL,
                    };
                    let b = Vector::<$d>::new([0.0; $d]);
                    let err = lu.solve_vec(b).unwrap_err();
                    assert_eq!(err, LaError::Singular { pivot_col: $d - 1 });
                }

                /// `solve_vec` must surface `NonFinite` with the corrupt
                /// cell's coordinates when a stored U diagonal is NaN,
                /// even though `factor` cannot produce such a
                /// factorization. The error must pinpoint `(D-1, D-1)`
                /// per the [`LaError::NonFinite`] convention.
                #[test]
                fn [<solve_vec_defensive_non_finite_diagonal_ $d d>]() {
                    let mut factors = Matrix::<$d>::identity();
                    factors.rows[$d - 1][$d - 1] = f64::NAN;

                    let mut piv = [0usize; $d];
                    for (i, p) in piv.iter_mut().enumerate() {
                        *p = i;
                    }

                    let lu = Lu::<$d> {
                        factors,
                        piv,
                        piv_sign: 1.0,
                        tol: DEFAULT_PIVOT_TOL,
                    };
                    let b = Vector::<$d>::new([1.0; $d]);
                    let err = lu.solve_vec(b).unwrap_err();
                    assert_eq!(
                        err,
                        LaError::NonFinite {
                            row: Some($d - 1),
                            col: $d - 1,
                        }
                    );
                }
            }
        };
    }

    gen_solve_vec_defensive_tests!(2);
    gen_solve_vec_defensive_tests!(3);
    gen_solve_vec_defensive_tests!(4);
    gen_solve_vec_defensive_tests!(5);

    // -----------------------------------------------------------------------
    // Const-evaluability tests.
    //
    // These prove that `Lu::det` and `Lu::solve_vec` are truly `const fn` by
    // forcing the compiler to evaluate them inside a `const` initializer.
    // `Lu::factor` is not (yet) `const fn` because it relies on `<[T]>::swap`,
    // which is not const-stable; we therefore construct `Lu<D>` directly.
    // -----------------------------------------------------------------------

    #[test]
    fn lu_det_const_eval_d2() {
        const DET: Result<f64, LaError> = {
            // Triangular factors with diag [2.0, 3.0] and no row swaps.
            let mut factors = Matrix::<2>::identity();
            factors.rows[0][0] = 2.0;
            factors.rows[1][1] = 3.0;
            let lu = Lu::<2> {
                factors,
                piv: [0, 1],
                piv_sign: 1.0,
                tol: DEFAULT_PIVOT_TOL,
            };
            lu.det()
        };
        assert_eq!(DET, Ok(6.0));
    }

    #[test]
    fn lu_det_const_eval_d3_row_swap() {
        const DET: Result<f64, LaError> = {
            // Identity factors but `piv_sign = -1.0` encoding a single row swap;
            // the determinant magnitude is 1 but the sign flips.
            let lu = Lu::<3> {
                factors: Matrix::<3>::identity(),
                piv: [1, 0, 2],
                piv_sign: -1.0,
                tol: DEFAULT_PIVOT_TOL,
            };
            lu.det()
        };
        assert_eq!(DET, Ok(-1.0));
    }

    #[test]
    fn lu_solve_vec_const_eval_d2() {
        // Identity LU ⇒ solve_vec returns the permuted RHS untouched.
        const X: [f64; 2] = {
            let lu = Lu::<2> {
                factors: Matrix::<2>::identity(),
                piv: [0, 1],
                piv_sign: 1.0,
                tol: DEFAULT_PIVOT_TOL,
            };
            let b = Vector::<2>::new([1.0, 2.0]);
            match lu.solve_vec(b) {
                Ok(v) => v.into_array(),
                Err(_) => [0.0, 0.0],
            }
        };
        assert!((X[0] - 1.0).abs() <= 1e-12);
        assert!((X[1] - 2.0).abs() <= 1e-12);
    }
}
