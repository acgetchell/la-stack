#![forbid(unsafe_code)]

//! LU decomposition and solves.
//!
//! The implementation computes `P A = L U` with partial pivoting. Partial
//! pivoting is a practical finite-precision strategy rather than an
//! unconditional accuracy guarantee; see `REFERENCES.md` \[1-3, 11-12\] for
//! stability analysis and standard algorithmic background.

use core::hint::cold_path;

use crate::matrix::Matrix;
use crate::scaled_product::{RangeCheckedProduct, ScaledProduct, range_checked_product};
use crate::vector::Vector;
use crate::{ArithmeticOperation, FactorizationKind, LaError, Tolerance};

/// LU decomposition (PA = LU) with partial pivoting.
///
/// `Lu<0>` represents the empty factorization. Its determinant is the empty
/// product `1.0`, and solving against [`Vector<0>`] returns [`Vector<0>`].
/// Numerical solves and determinants remain subject to binary64 rounding and
/// matrix conditioning; this type does not provide a certified error bound.
#[must_use]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lu<const D: usize> {
    factors: LuFactors<D>,
    permutation: RowPermutation<D>,
}

/// Finite LU factor storage.
///
/// [`Lu::factor_finite`] separately proves that every `U[i,i]` satisfies the
/// factorization tolerance before this storage becomes part of a [`Lu`].
#[derive(Clone, Copy, Debug, PartialEq)]
struct LuFactors<const D: usize> {
    storage: [[f64; D]; D],
}

impl<const D: usize> LuFactors<D> {
    /// Validate and finalize raw factorization work storage as finite factors.
    #[inline]
    const fn try_from_computation(storage: [[f64; D]; D]) -> Result<Self, LaError> {
        let mut row = 0;
        while row < D {
            let mut col = 0;
            while col < D {
                if !storage[row][col].is_finite() {
                    return Err(LaError::non_finite_computation_matrix(
                        ArithmeticOperation::LuFactorization,
                        row,
                        col,
                    ));
                }
                col += 1;
            }
            row += 1;
        }

        Ok(Self { storage })
    }

    /// Borrow a factor row.
    #[inline]
    #[must_use]
    const fn row(&self, index: usize) -> &[f64; D] {
        &self.storage[index]
    }

    /// Return a diagonal entry of `U`.
    #[inline]
    #[must_use]
    const fn diag(&self, index: usize) -> f64 {
        self.storage[index][index]
    }
}

/// Source-row permutation and its determinant parity.
///
/// Starting from identity and permitting only synchronized swaps makes every
/// stored source row in-bounds and unique while keeping parity inseparable from
/// the index mapping.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct RowPermutation<const D: usize> {
    source_rows: [usize; D],
    odd: bool,
}

impl<const D: usize> RowPermutation<D> {
    /// Construct the identity permutation.
    const fn identity() -> Self {
        let mut source_rows = [0; D];
        let mut row = 0;
        while row < D {
            source_rows[row] = row;
            row += 1;
        }
        Self {
            source_rows,
            odd: false,
        }
    }

    /// Apply one row swap and update parity atomically.
    const fn swap(&mut self, left: usize, right: usize) {
        if left != right {
            let source_row = self.source_rows[left];
            self.source_rows[left] = self.source_rows[right];
            self.source_rows[right] = source_row;
            self.odd = !self.odd;
        }
    }

    /// Return the original source row now occupying `row`.
    const fn source_row(&self, row: usize) -> usize {
        self.source_rows[row]
    }

    /// Return whether the permutation contains an odd number of swaps.
    const fn is_odd(&self) -> bool {
        self.odd
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
    pub(crate) fn factor_finite(a: Matrix<D>, tol: Tolerance) -> Result<Self, LaError> {
        let mut rows = a.into_rows();
        let tolerance = tol.get();
        let mut permutation = RowPermutation::identity();

        {
            let rows = &mut rows;

            for k in 0..D {
                // Choose pivot row.
                let mut pivot_row = k;
                let mut pivot_abs = rows[k][k].abs();

                #[expect(
                    clippy::needless_range_loop,
                    reason = "the row index identifies the pivot later used for synchronized matrix and permutation swaps"
                )]
                for r in (k + 1)..D {
                    let v = rows[r][k].abs();
                    if v > pivot_abs {
                        pivot_abs = v;
                        pivot_row = r;
                    }
                }

                if pivot_abs <= tolerance {
                    cold_path();

                    // A non-finite value produced in an earlier update does not
                    // participate in `v > pivot_abs` comparisons. Scan only on
                    // this cold failure path so it cannot be masked as singular.
                    for (row, values) in rows.iter().enumerate() {
                        for (col, value) in values.iter().enumerate() {
                            if !value.is_finite() {
                                return Err(LaError::non_finite_computation_matrix(
                                    ArithmeticOperation::LuFactorization,
                                    row,
                                    col,
                                ));
                            }
                        }
                    }

                    return Err(LaError::singular_numerical(
                        k,
                        FactorizationKind::Lu,
                        pivot_abs,
                        tolerance,
                    ));
                }

                if pivot_row != k {
                    rows.swap(k, pivot_row);
                    permutation.swap(k, pivot_row);
                }

                let pivot = rows[k][k];

                // Eliminate below pivot.
                for r in (k + 1)..D {
                    let mult = rows[r][k] / pivot;
                    rows[r][k] = mult;

                    #[expect(
                        clippy::needless_range_loop,
                        reason = "the column index pairs pivot-row reads with eliminated-row writes in the in-place update"
                    )]
                    for c in (k + 1)..D {
                        let updated = (-mult).mul_add(rows[k][c], rows[r][c]);
                        rows[r][c] = updated;
                    }
                }
            }
        }

        let factors = LuFactors::try_from_computation(rows)?;

        Ok(Self {
            factors,
            permutation,
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
        let mut x = [0.0; D];
        let b = b.as_array();
        let mut i = 0;

        if D <= 4 {
            while i < D {
                x[i] = b[self.permutation.source_row(i)];
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
                    return Err(LaError::non_finite_computation_step(
                        ArithmeticOperation::LuSolve,
                        i,
                    ));
                }
                x[i] = sum;
                i += 1;
            }
        } else {
            // Larger fixed dimensions avoid an extra pass by reading the
            // pivoted right-hand side directly into forward substitution.
            while i < D {
                let mut sum = b[self.permutation.source_row(i)];
                let row = self.factors.row(i);
                let mut j = 0;
                while j < i {
                    sum = (-row[j]).mul_add(x[j], sum);
                    j += 1;
                }
                if !sum.is_finite() {
                    cold_path();
                    return Err(LaError::non_finite_computation_step(
                        ArithmeticOperation::LuSolve,
                        i,
                    ));
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
                return Err(LaError::non_finite_computation_step(
                    ArithmeticOperation::LuSolve,
                    i,
                ));
            }

            let quotient = sum / diag;
            if !quotient.is_finite() {
                cold_path();
                return Err(LaError::non_finite_computation_step(
                    ArithmeticOperation::LuSolve,
                    i,
                ));
            }
            x[i] = quotient;
            ii += 1;
        }

        Vector::from_computation(x, ArithmeticOperation::LuSolve)
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
    /// Diagonal pivots are multiplied directly while each non-zero running
    /// product remains finite and normal. If direct accumulation detects range
    /// loss, all pivots are recomputed with power-of-two scaling before a
    /// premature overflow or underflow can affect the returned determinant.
    /// The final product is rounded to `f64`; a non-zero magnitude below the
    /// binary64 range may round to zero. No certified absolute error bound is
    /// provided.
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if the final scaled determinant cannot be
    /// represented as a finite `f64`.
    #[inline]
    pub const fn det(&self) -> Result<f64, LaError> {
        let mut det = if self.permutation.is_odd() { -1.0 } else { 1.0 };
        let mut i = 0;
        while i < D {
            let factor = self.factors.diag(i);
            match range_checked_product(det, factor) {
                RangeCheckedProduct::Safe(next) => det = next,
                RangeCheckedProduct::NeedsScaling => {
                    cold_path();
                    return self.scaled_det();
                }
            }
            i += 1;
        }
        Ok(det)
    }

    /// Recompute the determinant with normalized mantissa/exponent scaling.
    #[cold]
    const fn scaled_det(&self) -> Result<f64, LaError> {
        let mut product = ScaledProduct::new(self.permutation.is_odd());
        let mut i = 0;
        while i < D {
            product.multiply(self.factors.diag(i));
            i += 1;
        }

        if let Some(det) = product.finish() {
            Ok(det)
        } else {
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::Determinant,
                D.saturating_sub(1),
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use core::hint::black_box;

    use approx::assert_abs_diff_eq;
    use pastey::paste;

    use super::*;
    use crate::DEFAULT_SINGULAR_TOL;

    const TWO_NEG_800: f64 = f64::from_bits(223_u64 << 52);
    const TWO_POS_800: f64 = f64::from_bits(1823_u64 << 52);

    #[test]
    fn row_permutation_keeps_mapping_and_parity_synchronized() {
        let mut permutation = RowPermutation::<4>::identity();
        assert_eq!(
            core::array::from_fn(|row| permutation.source_row(row)),
            [0, 1, 2, 3]
        );
        assert!(!permutation.is_odd());

        permutation.swap(0, 3);
        assert_eq!(
            core::array::from_fn(|row| permutation.source_row(row)),
            [3, 1, 2, 0]
        );
        assert!(permutation.is_odd());

        permutation.swap(1, 2);
        assert_eq!(
            core::array::from_fn(|row| permutation.source_row(row)),
            [3, 2, 1, 0]
        );
        assert!(!permutation.is_odd());
    }

    macro_rules! gen_pivoting_solve_and_det_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<lu_solve_pivoting_ $d d>]() {
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
                fn [<lu_det_pivoting_ $d d>]() {
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

    gen_pivoting_solve_and_det_tests!(2);
    gen_pivoting_solve_and_det_tests!(3);
    gen_pivoting_solve_and_det_tests!(4);
    gen_pivoting_solve_and_det_tests!(5);

    macro_rules! gen_tridiagonal_smoke_solve_and_det_tests {
        ($d:literal $(, #[$stack_array_expectation:meta])?) => {
            paste! {
                #[test]
                fn [<lu_solve_tridiagonal_smoke_ $d d>]() {
                    // Public API path under test:
                    // Matrix::lu (pub) -> Lu::solve (pub).

                    // Classic SPD tridiagonal: 2 on diagonal, -1 on sub/super-diagonals.
                    $(#[$stack_array_expectation])?
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
                fn [<lu_det_tridiagonal_smoke_ $d d>]() {
                    // Public API path under test:
                    // Matrix::lu (pub) -> Lu::det (pub).

                    // Classic SPD tridiagonal: 2 on diagonal, -1 on sub/super-diagonals.
                    // Determinant is known exactly: det = D + 1.
                    $(#[$stack_array_expectation])?
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

    gen_tridiagonal_smoke_solve_and_det_tests!(16);
    gen_tridiagonal_smoke_solve_and_det_tests!(32);
    gen_tridiagonal_smoke_solve_and_det_tests!(
        64,
        #[expect(
            clippy::large_stack_arrays,
            reason = "the test deliberately exercises the crate's stack-allocated matrix storage"
        )]
    );

    #[test]
    fn solve_0x0_returns_empty_vector_and_unit_det() {
        let a = Matrix::<0>::zero();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        assert_eq!(lu.det(), Ok(1.0));
        assert!(
            lu.solve(Vector::<0>::zero())
                .unwrap()
                .into_array()
                .is_empty()
        );
    }

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
    fn det_ordinary_factors_matches_direct_product_bits() {
        let diagonal = [1.5, -2.0, 0.25, 8.0];
        let mut rows = [[0.0; 4]; 4];
        let mut expected = 1.0;
        for (i, factor) in diagonal.into_iter().enumerate() {
            rows[i][i] = factor;
            expected *= factor;
        }

        let lu = Matrix::<4>::try_from_rows(rows)
            .unwrap()
            .lu(DEFAULT_SINGULAR_TOL)
            .unwrap();
        assert_eq!(lu.det().unwrap().to_bits(), expected.to_bits());
    }

    #[test]
    fn singular_detected() {
        let a = Matrix::<2>::try_from_rows(black_box([[1.0, 2.0], [2.0, 4.0]])).unwrap();
        let err = a.lu(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::singular_numerical(1, FactorizationKind::Lu, 0.0, DEFAULT_SINGULAR_TOL.get())
        );
    }

    #[test]
    fn singular_due_to_tolerance_at_first_pivot() {
        // Not exactly singular, but below DEFAULT_SINGULAR_TOL.
        let a = Matrix::<2>::try_from_rows(black_box([[1e-13, 0.0], [0.0, 1.0]])).unwrap();
        let err = a.lu(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::singular_numerical(
                0,
                FactorizationKind::Lu,
                1e-13,
                DEFAULT_SINGULAR_TOL.get()
            )
        );
    }

    #[test]
    fn non_finite_detected_in_trailing_update() {
        let a = Matrix::<3>::try_from_rows([
            [1.0, f64::MAX, 0.0],
            [-1.0, f64::MAX, 0.0],
            [0.0, 0.0, 1.0],
        ])
        .unwrap();

        let err = a.lu(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_matrix(ArithmeticOperation::LuFactorization, 1, 1)
        );
    }

    #[test]
    fn generated_non_finite_takes_precedence_over_later_singular_pivot() {
        // The first update generates infinities, and the next generates NaN.
        // NaN does not win a pivot comparison and must not be masked as singular.
        let a = Matrix::<4>::try_from_rows([
            [1.0, f64::MAX, 0.0, 0.0],
            [1.0, f64::MAX, 0.0, 0.0],
            [-1.0, f64::MAX, 0.0, 0.0],
            [-1.0, f64::MAX, 0.0, 0.0],
        ])
        .unwrap();

        let err = a.lu(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_matrix(ArithmeticOperation::LuFactorization, 1, 1,)
        );
    }

    #[test]
    fn solve_non_finite_forward_substitution_overflow() {
        // L has a -1 multiplier, and a large RHS makes forward substitution overflow.
        let a = Matrix::<3>::try_from_rows([[1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<3>::new([1.0e308, 1.0e308, 0.0]);
        let err = lu.solve(b).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_step(ArithmeticOperation::LuSolve, 1)
        );
    }

    #[test]
    fn solve_non_finite_forward_substitution_overflow_fused_branch_5d() {
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
        assert_eq!(
            err,
            LaError::non_finite_computation_step(ArithmeticOperation::LuSolve, 1)
        );
    }

    #[test]
    fn solve_non_finite_back_substitution_overflow() {
        // Make x[1] overflow during back substitution, then ensure it is detected on the next row.
        let a = Matrix::<2>::try_from_rows([[1.0, 1.0], [0.0, 2.0e-12]]).unwrap();
        let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<2>::new([0.0, 1.0e300]);
        let err = lu.solve(b).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_step(ArithmeticOperation::LuSolve, 1)
        );
    }

    #[test]
    fn solve_non_finite_back_substitution_sum_overflow() {
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
        assert_eq!(
            err,
            LaError::non_finite_computation_step(ArithmeticOperation::LuSolve, 1)
        );
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
        assert_eq!(
            lu.det(),
            Err(LaError::non_finite_computation_step(
                ArithmeticOperation::Determinant,
                4
            ))
        );
    }

    #[test]
    fn det_balances_extreme_diagonals_independently_of_storage_order() {
        let zero_tolerance = Tolerance::try_new(0.0).unwrap();
        for diagonal in [
            [TWO_NEG_800, TWO_NEG_800, TWO_POS_800, TWO_POS_800],
            [TWO_POS_800, TWO_POS_800, TWO_NEG_800, TWO_NEG_800],
        ] {
            let mut rows = [[0.0; 4]; 4];
            for (i, value) in diagonal.into_iter().enumerate() {
                rows[i][i] = value;
            }

            let lu = Matrix::<4>::try_from_rows(rows)
                .unwrap()
                .lu(zero_tolerance)
                .unwrap();
            assert_eq!(lu.det(), Ok(1.0));
        }
    }

    #[test]
    fn matrix_det_fallback_inherits_balanced_extreme_accumulation() {
        let zero_tolerance = Tolerance::try_new(0.0).unwrap();
        for diagonal in [
            [TWO_NEG_800, TWO_NEG_800, TWO_POS_800, TWO_POS_800, 1.0, 1.0],
            [TWO_POS_800, TWO_POS_800, TWO_NEG_800, TWO_NEG_800, 1.0, 1.0],
        ] {
            let mut rows = [[0.0; 6]; 6];
            for (i, value) in diagonal.into_iter().enumerate() {
                rows[i][i] = value;
            }

            let matrix = Matrix::<6>::try_from_rows(rows).unwrap();
            assert_eq!(matrix.det(), Ok(1.0));
            assert_eq!(matrix.lu(zero_tolerance).unwrap().det(), Ok(1.0));
        }
    }

    #[test]
    fn det_rounds_final_tiny_magnitude_to_zero() {
        let zero_tolerance = Tolerance::try_new(0.0).unwrap();
        let positive =
            Matrix::<2>::try_from_rows([[TWO_NEG_800, 0.0], [0.0, TWO_NEG_800]]).unwrap();
        let positive_det = positive.lu(zero_tolerance).unwrap().det().unwrap();
        assert_eq!(positive_det.to_bits(), 0.0f64.to_bits());

        let negative =
            Matrix::<2>::try_from_rows([[-TWO_NEG_800, 0.0], [0.0, TWO_NEG_800]]).unwrap();
        let negative_det = negative.lu(zero_tolerance).unwrap().det().unwrap();
        assert_eq!(negative_det.to_bits(), (-0.0f64).to_bits());
    }

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
            let Ok(factors) = LuFactors::try_from_computation([[2.0, 0.0], [0.0, 3.0]]) else {
                panic!("LU test factors must be finite");
            };
            let lu = Lu::<2> {
                factors,
                permutation: RowPermutation::identity(),
            };
            lu.det()
        };
        assert_eq!(DET, Ok(6.0));
    }

    #[test]
    fn lu_det_const_eval_d3_row_swap() {
        const DET: Result<f64, LaError> = {
            // Identity factors with odd row-swap parity;
            // the determinant magnitude is 1 but the sign flips.
            let Ok(factors) = LuFactors::try_from_computation(Matrix::<3>::identity().into_rows())
            else {
                panic!("LU test factors must be usable");
            };
            let mut permutation = RowPermutation::identity();
            permutation.swap(0, 1);
            let lu = Lu::<3> {
                factors,
                permutation,
            };
            lu.det()
        };
        assert_eq!(DET, Ok(-1.0));
    }

    #[test]
    fn lu_solve_const_eval_d2() {
        // Identity LU ⇒ solve returns the permuted RHS untouched.
        const X: Result<Vector<2>, LaError> = {
            let Ok(factors) = LuFactors::try_from_computation(Matrix::<2>::identity().into_rows())
            else {
                panic!("LU test factors must be usable");
            };
            let lu = Lu::<2> {
                factors,
                permutation: RowPermutation::identity(),
            };
            let b = Vector::<2>::new([1.0, 2.0]);
            lu.solve(b)
        };
        let x = X.unwrap().into_array();
        assert!((x[0] - 1.0).abs() <= 1e-12);
        assert!((x[1] - 2.0).abs() <= 1e-12);
    }
}
