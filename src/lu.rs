//! LU decomposition and solves.

use core::hint::cold_path;

use crate::matrix::Matrix;
use crate::vector::Vector;
use crate::{LaError, Tolerance};

/// LU decomposition (PA = LU) with partial pivoting.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lu<const D: usize> {
    factors: LuFactors<D>,
    piv: [usize; D],
    piv_sign: f64,
}

/// In-place LU factor storage whose `U` diagonal is finite and usable.
///
/// Construction through [`Lu::factor_finite`] proves every stored entry is
/// finite and every `U[i,i]` satisfies the factorization tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
struct LuFactors<const D: usize> {
    storage: Matrix<D>,
}

impl<const D: usize> LuFactors<D> {
    /// Construct factors after LU factorization has proven the storage invariant.
    #[inline]
    const fn new_unchecked(storage: Matrix<D>) -> Self {
        Self { storage }
    }

    /// Borrow a factor row.
    #[inline]
    #[must_use]
    const fn row(&self, index: usize) -> &[f64; D] {
        &self.storage.rows()[index]
    }

    /// Return a diagonal entry of `U`.
    #[inline]
    #[must_use]
    const fn diag(&self, index: usize) -> f64 {
        self.storage.rows()[index][index]
    }
}

impl<const D: usize> Lu<D> {
    /// Factor a finite square matrix into in-place LU storage for
    /// [`Matrix::lu`].
    ///
    /// The input has already proven finite entries, so LU construction rejects
    /// numerically singular pivots and non-finite elimination intermediates
    /// before callers can observe a [`Lu`] value. Completed factor storage is
    /// checked before return so successful factors do not contain a non-finite
    /// value produced during elimination.
    #[inline]
    #[allow(clippy::needless_range_loop)]
    pub(crate) fn factor_finite(a: Matrix<D>, tol: Tolerance) -> Result<Self, LaError> {
        let mut lu = a;
        let tol = tol.get();

        let mut piv = [0usize; D];
        for (i, p) in piv.iter_mut().enumerate() {
            *p = i;
        }

        let mut piv_sign = 1.0;

        {
            let rows = lu.rows_mut_unchecked();

            for k in 0..D {
                // Choose pivot row.
                let mut pivot_row = k;
                let mut pivot_abs = rows[k][k].abs();

                for r in (k + 1)..D {
                    let v = rows[r][k].abs();
                    if v > pivot_abs {
                        pivot_abs = v;
                        pivot_row = r;
                    }
                }

                if pivot_abs <= tol {
                    cold_path();
                    return Err(LaError::Singular { pivot_col: k });
                }

                if pivot_row != k {
                    rows.swap(k, pivot_row);
                    piv.swap(k, pivot_row);
                    piv_sign = -piv_sign;
                }

                let pivot = rows[k][k];

                // Eliminate below pivot.
                for r in (k + 1)..D {
                    let mult = rows[r][k] / pivot;
                    rows[r][k] = mult;

                    for c in (k + 1)..D {
                        let updated = (-mult).mul_add(rows[k][c], rows[r][c]);
                        rows[r][c] = updated;
                    }
                }
            }
        }

        let lu = lu.validate_finite()?;

        Ok(Self {
            factors: LuFactors::new_unchecked(lu),
            piv,
            piv_sign,
        })
    }

    /// Solve `A x = b` using this LU factorization.
    ///
    /// [`Vector`] is finite by construction, so this method only checks computed
    /// substitution overflows. It performs floating-point forward/back
    /// substitution and does not provide a certified absolute rounding-error
    /// bound for the returned solution.
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
    /// # Errors
    /// Returns [`LaError::NonFinite`] if a computed substitution intermediate
    /// overflows to NaN or infinity.
    #[inline]
    pub const fn solve(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        self.solve_finite(b)
    }

    /// Solve `A x = b` using this LU factorization and a finite right-hand side.
    ///
    /// The right-hand side entries and stored factors are known finite, so this
    /// path only checks computed substitution overflows.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if a computed substitution intermediate
    /// overflows to NaN or infinity.
    #[inline]
    pub(crate) const fn solve_finite(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
        let mut x = [0.0; D];
        let b = b.as_array();
        let mut i = 0;

        if D <= 4 {
            while i < D {
                x[i] = b[self.piv[i]];
                i += 1;
            }

            // Tiny matrices benchmark better when pivoted RHS materialization
            // stays separate from forward substitution.
            i = 0;
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
        } else {
            // Larger fixed dimensions avoid an extra pass by reading the
            // pivoted right-hand side directly into forward substitution.
            while i < D {
                let mut sum = b[self.piv[i]];
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
        }

        // Back substitution for U.
        let mut ii = 0;
        while ii < D {
            let i = D - 1 - ii;
            let mut sum = x[i];
            let row = self.factors.row(i);
            let mut j = i + 1;
            while j < D {
                sum = (-row[j]).mul_add(x[j], sum);
                j += 1;
            }

            let diag = row[i];
            if !sum.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }

            let quotient = sum / diag;
            if !quotient.is_finite() {
                cold_path();
                return Err(LaError::non_finite_at(i));
            }
            x[i] = quotient;
            ii += 1;
        }

        Ok(Vector::new_unchecked(x))
    }

    /// Determinant of the original matrix.
    ///
    /// # Examples
    /// ```
    /// use la_stack::prelude::*;
    ///
    /// # fn main() -> Result<(), LaError> {
    /// let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    /// let lu = a.lu(DEFAULT_SINGULAR_TOL)?;
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
            det *= self.factors.diag(i);
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
    use crate::DEFAULT_SINGULAR_TOL;

    use core::hint::black_box;

    use approx::assert_abs_diff_eq;
    use pastey::paste;

    macro_rules! gen_public_api_pivoting_solve_and_det_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_api_lu_solve_pivoting_ $d d>]() {
                    // Public API path under test:
                    // Matrix::lu (pub) -> Lu::solve (pub).

                    // Permutation matrix that swaps the first two basis vectors.
                    // This forces pivoting in column 0 for any D >= 2.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = 1.0;
                    }
                    rows.swap(0, 1);

                    let a = Matrix::<$d>::try_from_rows(black_box(rows)).unwrap();
                    let lu_fn: fn(Matrix<$d>, Tolerance) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_SINGULAR_TOL).unwrap();

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
                        black_box(Lu::<$d>::solve);
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

                    let a = Matrix::<$d>::try_from_rows(black_box(rows)).unwrap();
                    let lu_fn: fn(Matrix<$d>, Tolerance) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_SINGULAR_TOL).unwrap();

                    // Row swap ⇒ determinant sign flip.
                    let det_fn: fn(&Lu<$d>) -> Result<f64, LaError> =
                        black_box(Lu::<$d>::det);
                    assert_abs_diff_eq!(det_fn(&lu).unwrap(), -1.0, epsilon = 1e-12);
                }
            }
        };
    }

    gen_public_api_pivoting_solve_and_det_tests!(2);
    gen_public_api_pivoting_solve_and_det_tests!(3);
    gen_public_api_pivoting_solve_and_det_tests!(4);
    gen_public_api_pivoting_solve_and_det_tests!(5);

    macro_rules! gen_public_api_tridiagonal_smoke_solve_and_det_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<public_api_lu_solve_tridiagonal_smoke_ $d d>]() {
                    // Public API path under test:
                    // Matrix::lu (pub) -> Lu::solve (pub).

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

                    let a = Matrix::<$d>::try_from_rows(black_box(rows)).unwrap();
                    let lu_fn: fn(Matrix<$d>, Tolerance) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_SINGULAR_TOL).unwrap();

                    // Choose x = 1, so b = A x is simple: [1, 0, 0, ..., 0, 1].
                    let mut b_arr = [0.0f64; $d];
                    b_arr[0] = 1.0;
                    b_arr[$d - 1] = 1.0;
                    let b = Vector::<$d>::new(black_box(b_arr));

                    let solve_fn: fn(&Lu<$d>, Vector<$d>) -> Result<Vector<$d>, LaError> =
                        black_box(Lu::<$d>::solve);
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

                    let a = Matrix::<$d>::try_from_rows(black_box(rows)).unwrap();
                    let lu_fn: fn(Matrix<$d>, Tolerance) -> Result<Lu<$d>, LaError> =
                        black_box(Matrix::<$d>::lu);
                    let lu = lu_fn(a, DEFAULT_SINGULAR_TOL).unwrap();

                    let det_fn: fn(&Lu<$d>) -> Result<f64, LaError> =
                        black_box(Lu::<$d>::det);
                    assert_abs_diff_eq!(det_fn(&lu).unwrap(), f64::from($d) + 1.0, epsilon = 1e-8);
                }
            }
        };
    }

    gen_public_api_tridiagonal_smoke_solve_and_det_tests!(16);
    gen_public_api_tridiagonal_smoke_solve_and_det_tests!(32);
    gen_public_api_tridiagonal_smoke_solve_and_det_tests!(64);

    #[test]
    fn solve_1x1() {
        let a = Matrix::<1>::try_from_rows(black_box([[2.0]])).unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<1>::new(black_box([6.0]));
        let solve_fn: fn(&Lu<1>, Vector<1>) -> Result<Vector<1>, LaError> =
            black_box(Lu::<1>::solve);
        let x = solve_fn(&lu, b).unwrap().into_array();
        assert_abs_diff_eq!(x[0], 3.0, epsilon = 1e-12);

        let det_fn: fn(&Lu<1>) -> Result<f64, LaError> = black_box(Lu::<1>::det);
        assert_abs_diff_eq!(det_fn(&lu).unwrap(), 2.0, epsilon = 0.0);
    }

    #[test]
    fn solve_2x2_basic() {
        let a = Matrix::<2>::try_from_rows(black_box([[1.0, 2.0], [3.0, 4.0]])).unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();
        let b = Vector::<2>::new(black_box([5.0, 11.0]));

        let solve_fn: fn(&Lu<2>, Vector<2>) -> Result<Vector<2>, LaError> =
            black_box(Lu::<2>::solve);
        let x = solve_fn(&lu, b).unwrap().into_array();

        assert_abs_diff_eq!(x[0], 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], 2.0, epsilon = 1e-12);
    }

    #[test]
    fn det_2x2_basic() {
        let a = Matrix::<2>::try_from_rows(black_box([[1.0, 2.0], [3.0, 4.0]])).unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        let det_fn: fn(&Lu<2>) -> Result<f64, LaError> = black_box(Lu::<2>::det);
        assert_abs_diff_eq!(det_fn(&lu).unwrap(), -2.0, epsilon = 1e-12);
    }

    #[test]
    fn det_requires_pivot_sign() {
        // Row swap ⇒ determinant sign flip.
        let a = Matrix::<2>::try_from_rows(black_box([[0.0, 1.0], [1.0, 0.0]])).unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        let det_fn: fn(&Lu<2>) -> Result<f64, LaError> = black_box(Lu::<2>::det);
        assert_abs_diff_eq!(det_fn(&lu).unwrap(), -1.0, epsilon = 0.0);
    }

    #[test]
    fn solve_requires_pivoting() {
        let a = Matrix::<2>::try_from_rows(black_box([[0.0, 1.0], [1.0, 0.0]])).unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();
        let b = Vector::<2>::new(black_box([1.0, 2.0]));

        let solve_fn: fn(&Lu<2>, Vector<2>) -> Result<Vector<2>, LaError> =
            black_box(Lu::<2>::solve);
        let x = solve_fn(&lu, b).unwrap().into_array();

        assert_abs_diff_eq!(x[0], 2.0, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn singular_detected() {
        let a = Matrix::<2>::try_from_rows(black_box([[1.0, 2.0], [2.0, 4.0]])).unwrap();
        let err = a.lu(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(err, LaError::Singular { pivot_col: 1 });
    }

    #[test]
    fn singular_due_to_tolerance_at_first_pivot() {
        // Not exactly singular, but below DEFAULT_SINGULAR_TOL.
        let a = Matrix::<2>::try_from_rows(black_box([[1e-13, 0.0], [0.0, 1.0]])).unwrap();
        let err = a.lu(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(err, LaError::Singular { pivot_col: 0 });
    }

    #[test]
    fn matrix_constructor_rejects_nonfinite_pivot_entry() {
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
    fn matrix_constructor_rejects_nonfinite_pivot_column_entry() {
        let err = Matrix::<2>::try_from_rows([[1.0, 0.0], [f64::INFINITY, 1.0]]).unwrap_err();
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
        let a = Matrix::<3>::try_from_rows([
            [1.0, f64::MAX, 0.0],
            [-1.0, f64::MAX, 0.0],
            [0.0, 0.0, 1.0],
        ])
        .unwrap();

        let err = a.lu(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::NonFinite {
                row: Some(1),
                col: 1,
            }
        );
    }

    #[test]
    fn solve_nonfinite_forward_substitution_overflow() {
        // L has a -1 multiplier, and a large RHS makes forward substitution overflow.
        let a = Matrix::<3>::try_from_rows([[1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<3>::new([1.0e308, 1.0e308, 0.0]);
        let err = lu.solve(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { row: None, col: 1 });
    }

    #[test]
    fn solve_nonfinite_forward_substitution_overflow_fused_branch_5d() {
        // Exercises the D >= 5 fused pivot/forward-substitution branch with the
        // same overflowing L multiplier as the D3 test.
        let a = Matrix::<5>::try_from_rows([
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [-1.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ])
        .unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<5>::new([1.0e308, 1.0e308, 0.0, 0.0, 0.0]);
        let err = lu.solve(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { row: None, col: 1 });
    }

    #[test]
    fn solve_nonfinite_back_substitution_overflow() {
        // Make x[1] overflow during back substitution, then ensure it is detected on the next row.
        let a = Matrix::<2>::try_from_rows([[1.0, 1.0], [0.0, 2.0e-12]]).unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<2>::new([0.0, 1.0e300]);
        let err = lu.solve(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { row: None, col: 1 });
    }

    #[test]
    fn solve_nonfinite_back_substitution_sum_overflow() {
        // Upper-triangular U with a very large off-diagonal in row 1 and a
        // very large x[2] produced by the RHS.  The back-substitution
        // accumulator `sum = (-row[j]).mul_add(x[j], sum)` overflows while
        // reducing row 1, so the failure is detected via the `!sum.is_finite()`
        // branch of the combined diag/sum check (distinct from the
        // `q = sum / diag` overflow path covered above).
        let a = Matrix::<3>::try_from_rows([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0e200], [0.0, 0.0, 1.0]])
            .unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<3>::new([0.0, 0.0, 1.0e200]);
        let err = lu.solve(b).unwrap_err();
        assert_eq!(err, LaError::NonFinite { row: None, col: 1 });
    }

    #[test]
    fn det_rejects_product_overflow() {
        let a = Matrix::<5>::try_from_rows([
            [1.0e100, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0e100, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0e100, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0e100, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0e100],
        ])
        .unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();
        assert_eq!(lu.det(), Err(LaError::NonFinite { row: None, col: 3 }));
    }

    macro_rules! gen_solve_boundary_tests {
        ($d:literal) => {
            paste! {
                /// Raw non-finite right-hand sides are rejected before a
                /// public caller can construct a `Vector`.
                #[test]
                fn [<solve_rhs_constructor_rejects_non_finite_ $d d>]() {
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
            }
        };
    }

    gen_solve_boundary_tests!(2);
    gen_solve_boundary_tests!(3);
    gen_solve_boundary_tests!(4);
    gen_solve_boundary_tests!(5);

    // -----------------------------------------------------------------------
    // Const-evaluability tests.
    //
    // These prove that `Lu::det` and `Lu::solve` are truly `const fn` by
    // forcing the compiler to evaluate them inside a `const` initializer.
    // `Lu::factor` is not (yet) `const fn` because it relies on `<[T]>::swap`,
    // which is not const-stable; we therefore construct `Lu<D>` directly.
    // -----------------------------------------------------------------------

    #[test]
    fn lu_det_const_eval_d2() {
        const DET: Result<f64, LaError> = {
            // Triangular factors with diag [2.0, 3.0] and no row swaps.
            let factors = Matrix::<2>::from_rows_unchecked([[2.0, 0.0], [0.0, 3.0]]);
            let lu = Lu::<2> {
                factors: LuFactors::new_unchecked(factors),
                piv: [0, 1],
                piv_sign: 1.0,
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
                factors: LuFactors::new_unchecked(Matrix::<3>::identity()),
                piv: [1, 0, 2],
                piv_sign: -1.0,
            };
            lu.det()
        };
        assert_eq!(DET, Ok(-1.0));
    }

    #[test]
    fn lu_solve_const_eval_d2() {
        // Identity LU ⇒ solve returns the permuted RHS untouched.
        const X: [f64; 2] = {
            let lu = Lu::<2> {
                factors: LuFactors::new_unchecked(Matrix::<2>::identity()),
                piv: [0, 1],
                piv_sign: 1.0,
            };
            let b = Vector::<2>::new([1.0, 2.0]);
            match lu.solve(b) {
                Ok(v) => v.into_array(),
                Err(_) => [0.0, 0.0],
            }
        };
        assert!((X[0] - 1.0).abs() <= 1e-12);
        assert!((X[1] - 2.0).abs() <= 1e-12);
    }
}
