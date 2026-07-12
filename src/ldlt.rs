#![forbid(unsafe_code)]

//! LDLT factorization and solves.
//!
//! This module provides a stack-allocated LDLT factorization (`A = L D Lᵀ`)
//! without pivoting. Successful factors require an exactly symmetric input and
//! every computed diagonal pivot to be positive and above the caller's
//! tolerance. Computed zero and tolerance-small positive pivots are diagnosed
//! rather than returned in a usable factor. See `REFERENCES.md` \[4-6, 11-12\] for
//! Cholesky/LDLT background and pivoted symmetric-indefinite alternatives.
//!
//! # Preconditions
//! The input matrix must be **symmetric**.  This is a correctness contract, not a hint:
//! the factorization algorithm reads only the lower triangle and implicitly assumes the
//! upper triangle mirrors it exactly. Asymmetric inputs return [`LaError::Asymmetric`]
//! with an allowed absolute difference of `0.0` before factorization starts. IEEE-754
//! signed zeros compare equal and are accepted. Callers who know their matrices may
//! not be symmetric at all should use [`crate::Lu`] instead.

use core::hint::cold_path;

use crate::matrix::SymmetricMatrix;
use crate::scaled_product::{RangeCheckedProduct, ScaledProduct, range_checked_product};
use crate::vector::Vector;
use crate::{ArithmeticOperation, FactorizationKind, LaError, Tolerance};

/// LDLT factorization (`A = L D Lᵀ`) for exactly symmetric positive-definite matrices.
///
/// `Ldlt<0>` represents the empty factorization. Its determinant is the empty
/// product `1.0`, and solving against [`Vector<0>`] returns [`Vector<0>`].
///
/// This factorization is **not** a general-purpose symmetric-indefinite LDLT (no pivoting).
/// It assumes the input matrix is exactly symmetric and numerically positive
/// definite under the caller's absolute pivot tolerance. An uncoupled computed
/// zero or a tolerance-small positive pivot returns [`LaError::Singular`]; a
/// computed zero with non-zero remaining coupling returns
/// [`LaError::NotPositiveSemidefinite`]. Because pivots are computed in
/// binary64, success is not an exact proof that the stored matrix is positive
/// definite.
///
/// # Preconditions
/// The source matrix passed to [`Matrix::ldlt`](crate::Matrix::ldlt) must be
/// exactly symmetric (`A[i][j] == A[j][i]` for every mirrored pair). Asymmetric
/// inputs return [`LaError::Asymmetric`] before factorization starts; see
/// [`Matrix::ldlt`](crate::Matrix::ldlt) for details and alternatives.
///
/// # Storage
/// The factors are stored in one inline row-major array:
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
    storage: [[f64; D]; D],
}

impl<const D: usize> LdltFactors<D> {
    /// Store rows after the factorization loop has proven all factor invariants.
    #[inline]
    const fn from_proven_rows(storage: [[f64; D]; D]) -> Self {
        Self { storage }
    }

    /// Borrow a factor row.
    #[inline]
    #[must_use]
    const fn row(&self, index: usize) -> &[f64; D] {
        &self.storage[index]
    }

    /// Return a diagonal entry of `D`.
    #[inline]
    #[must_use]
    const fn diag(&self, index: usize) -> f64 {
        self.storage[index][index]
    }
}

impl<const D: usize> Ldlt<D> {
    /// Factor a finite, symmetry-proven matrix for
    /// [`Matrix::ldlt`](crate::Matrix::ldlt).
    ///
    /// Consuming [`SymmetricMatrix`] lets the factorization read only the lower
    /// triangle without revalidating symmetry. A successful result contains
    /// only finite factor storage with diagonals above `tol`.
    ///
    /// # Errors
    /// Returns [`LaError::NotPositiveSemidefinite`] for a negative pivot or a
    /// zero pivot with non-zero coupling, [`LaError::Singular`] for an uncoupled
    /// zero pivot or a positive pivot at or below `tol`, and [`LaError::NonFinite`]
    /// when a pivot, multiplier, or update is not finite.
    #[inline]
    pub(crate) fn factor_symmetric(a: SymmetricMatrix<D>, tol: Tolerance) -> Result<Self, LaError> {
        let mut rows = a.into_matrix().into_rows();
        let tolerance = tol.get();

        {
            let rows = &mut rows;

            // LDLT via symmetric rank-1 updates, using only the lower triangle.
            for j in 0..D {
                let d = rows[j][j];
                if !d.is_finite() {
                    cold_path();
                    return Err(LaError::non_finite_computation_matrix(
                        ArithmeticOperation::LdltFactorization,
                        j,
                        j,
                    ));
                }
                if d < 0.0 {
                    cold_path();
                    if let Some(error) = Self::non_finite_factor_error(rows) {
                        return Err(error);
                    }
                    return Err(LaError::not_positive_semidefinite_negative(j, d));
                }
                if d == 0.0 {
                    cold_path();
                    return Err(Self::zero_pivot_failure(rows, j, tolerance));
                }
                if d <= tolerance {
                    cold_path();
                    if let Some(error) = Self::non_finite_factor_error(rows) {
                        return Err(error);
                    }
                    return Err(LaError::singular_numerical(
                        j,
                        FactorizationKind::Ldlt,
                        d,
                        tolerance,
                    ));
                }
                if D <= 5 {
                    // Tiny matrices benchmark better when column normalization stays
                    // separate from the trailing update.
                    #[expect(
                        clippy::needless_range_loop,
                        reason = "the row index identifies the lower-triangle entry and any reported non-finite coordinate"
                    )]
                    for i in (j + 1)..D {
                        let l = rows[i][j] / d;
                        if !l.is_finite() {
                            cold_path();
                            return Err(LaError::non_finite_computation_matrix(
                                ArithmeticOperation::LdltFactorization,
                                i,
                                j,
                            ));
                        }
                        rows[i][j] = l;
                    }

                    for i in (j + 1)..D {
                        let l_i = rows[i][j];
                        let l_i_d = l_i * d;

                        #[expect(
                            clippy::needless_range_loop,
                            reason = "the triangular column index coordinates multiplier reads with in-place trailing-row writes"
                        )]
                        for k in (j + 1)..=i {
                            let l_k = rows[k][j];
                            let new_val = (-l_i_d).mul_add(l_k, rows[i][k]);
                            rows[i][k] = new_val;
                        }
                    }
                } else {
                    // Larger fixed dimensions avoid an extra column walk by updating
                    // each lower-triangular row prefix as soon as its multiplier is finite.
                    for i in (j + 1)..D {
                        let l_i = rows[i][j] / d;
                        if !l_i.is_finite() {
                            cold_path();
                            return Err(LaError::non_finite_computation_matrix(
                                ArithmeticOperation::LdltFactorization,
                                i,
                                j,
                            ));
                        }
                        rows[i][j] = l_i;

                        let l_i_d = l_i * d;

                        #[expect(
                            clippy::needless_range_loop,
                            reason = "the triangular column index coordinates normalized-column reads with the fused in-place update"
                        )]
                        for k in (j + 1)..=i {
                            let l_k = rows[k][j];
                            let new_val = (-l_i_d).mul_add(l_k, rows[i][k]);
                            rows[i][k] = new_val;
                        }
                    }
                }
            }
        }

        // Every computed lower-triangular entry is checked when it becomes a
        // pivot or multiplier; the untouched upper triangle remains finite input.
        Ok(Self {
            factors: LdltFactors::from_proven_rows(rows),
        })
    }

    /// Return the first non-finite factor cell in row-major order.
    fn non_finite_factor_error(rows: &[[f64; D]; D]) -> Option<LaError> {
        for (row, values) in rows.iter().enumerate() {
            for (col, value) in values.iter().enumerate() {
                if !value.is_finite() {
                    return Some(LaError::non_finite_computation_matrix(
                        ArithmeticOperation::LdltFactorization,
                        row,
                        col,
                    ));
                }
            }
        }
        None
    }

    /// Classify a zero pivot after checking factor storage and every coupling.
    fn zero_pivot_failure(rows: &[[f64; D]; D], pivot_col: usize, tolerance: f64) -> LaError {
        if let Some(error) = Self::non_finite_factor_error(rows) {
            return error;
        }
        for (row, values) in rows.iter().enumerate().skip(pivot_col + 1) {
            let coupling = values[pivot_col];
            if coupling != 0.0 {
                return LaError::not_positive_semidefinite_zero_coupling(pivot_col, row, coupling);
            }
        }
        LaError::singular_numerical(pivot_col, FactorizationKind::Ldlt, 0.0, tolerance)
    }

    /// Determinant of the original matrix.
    ///
    /// For a successfully constructed factorization, this is the product of
    /// the diagonal terms of `D`.
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
        let mut det = 1.0;
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
        let mut product = ScaledProduct::new(false);
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

    /// Solve `A x = b` using this LDLT factorization.
    ///
    /// [`Vector`] is finite by construction, so this method only checks computed
    /// substitution overflows. It performs floating-point substitution and does
    /// not provide a certified absolute rounding-error bound for the returned
    /// solution.
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
    /// let x = ldlt.solve(b)?.into_array();
    ///
    /// assert!((x[0] - (-0.125)).abs() <= 1e-12);
    /// assert!((x[1] - 0.75).abs() <= 1e-12);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    /// Returns [`LaError::NonFinite`] if a computed substitution intermediate
    /// overflows to NaN or infinity.
    #[inline]
    pub const fn solve(&self, b: Vector<D>) -> Result<Vector<D>, LaError> {
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
                return Err(LaError::non_finite_computation_step(
                    ArithmeticOperation::LdltSolve,
                    i,
                ));
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
                return Err(LaError::non_finite_computation_step(
                    ArithmeticOperation::LdltSolve,
                    i,
                ));
            }
            x[i] = quotient;
            i += 1;
        }

        if D <= 4 {
            // Tiny matrices benchmark better with the direct textbook dot
            // product for each row of Lᵀ.
            let mut ii = 0;
            while ii < D {
                let i = D - 1 - ii;
                let mut sum = x[i];
                let mut j = i + 1;
                while j < D {
                    sum = (-self.factors.row(j)[i]).mul_add(x[j], sum);
                    j += 1;
                }
                if !sum.is_finite() {
                    cold_path();
                    return Err(LaError::non_finite_computation_step(
                        ArithmeticOperation::LdltSolve,
                        i,
                    ));
                }
                x[i] = sum;
                ii += 1;
            }
        } else {
            // Larger fixed dimensions benchmark better by walking finalized
            // rows downward and scattering contributions into the remaining
            // contiguous lower-triangular row prefix.
            let mut jj = D;
            while jj > 0 {
                jj -= 1;

                let x_j = x[jj];
                if !x_j.is_finite() {
                    cold_path();
                    return Err(LaError::non_finite_computation_step(
                        ArithmeticOperation::LdltSolve,
                        jj,
                    ));
                }

                let row = self.factors.row(jj);
                let mut i = 0;
                while i < jj {
                    x[i] = (-row[i]).mul_add(x_j, x[i]);
                    i += 1;
                }
            }
        }

        Vector::from_computation(x, ArithmeticOperation::LdltSolve)
    }
}

#[cfg(test)]
mod tests {
    use core::hint::black_box;

    use approx::assert_abs_diff_eq;
    use pastey::paste;

    use super::*;
    use crate::DEFAULT_SINGULAR_TOL;
    use crate::matrix::Matrix;

    const TWO_NEG_800: f64 = f64::from_bits(223_u64 << 52);
    const TWO_NEG_38: f64 = f64::from_bits(985_u64 << 52);
    const TWO_POS_43: f64 = f64::from_bits(1066_u64 << 52);
    const TWO_POS_800: f64 = f64::from_bits(1823_u64 << 52);

    macro_rules! gen_ldlt_identity_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<ldlt_det_and_solve_identity_ $d d>]() {
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
                    let x = ldlt.solve(b).unwrap().into_array();

                    for i in 0..$d {
                        assert_abs_diff_eq!(x[i], b_arr[i], epsilon = 1e-12);
                    }
                }
            }
        };
    }

    gen_ldlt_identity_tests!(2);
    gen_ldlt_identity_tests!(3);
    gen_ldlt_identity_tests!(4);
    gen_ldlt_identity_tests!(5);

    macro_rules! gen_ldlt_diagonal_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<ldlt_det_and_solve_diagonal_spd_ $d d>]() {
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

                    let a = Matrix::<$d>::try_from_rows(black_box(rows)).unwrap();
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
                    let x = ldlt.solve(b).unwrap().into_array();

                    for i in 0..$d {
                        assert_abs_diff_eq!(x[i], b_arr[i] / diag[i], epsilon = 1e-12);
                    }
                }
            }
        };
    }

    gen_ldlt_diagonal_tests!(2);
    gen_ldlt_diagonal_tests!(3);
    gen_ldlt_diagonal_tests!(4);
    gen_ldlt_diagonal_tests!(5);

    #[test]
    fn solve_0x0_returns_empty_vector_and_unit_det() {
        let a = Matrix::<0>::zero();
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        assert_eq!(ldlt.det(), Ok(1.0));
        assert!(
            ldlt.solve(Vector::<0>::zero())
                .unwrap()
                .into_array()
                .is_empty()
        );
    }

    #[test]
    fn solve_2x2_known_spd() {
        let a = Matrix::<2>::try_from_rows(black_box([[4.0, 2.0], [2.0, 3.0]])).unwrap();
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<2>::new(black_box([1.0, 2.0]));
        let x = ldlt.solve(b).unwrap().into_array();

        assert_abs_diff_eq!(x[0], -0.125, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], 0.75, epsilon = 1e-12);
        assert_abs_diff_eq!(ldlt.det().unwrap(), 8.0, epsilon = 1e-12);
    }

    #[test]
    fn det_ordinary_factors_matches_direct_product_bits() {
        let diagonal = [1.5, 2.0, 0.25, 8.0];
        let mut rows = [[0.0; 4]; 4];
        let mut expected = 1.0;
        for (i, factor) in diagonal.into_iter().enumerate() {
            rows[i][i] = factor;
            expected *= factor;
        }

        let ldlt = Matrix::<4>::try_from_rows(rows)
            .unwrap()
            .ldlt(DEFAULT_SINGULAR_TOL)
            .unwrap();
        assert_eq!(ldlt.det().unwrap().to_bits(), expected.to_bits());
    }

    #[test]
    fn solve_3x3_spd_tridiagonal_smoke() {
        let a = Matrix::<3>::try_from_rows(black_box([
            [2.0, -1.0, 0.0],
            [-1.0, 2.0, -1.0],
            [0.0, -1.0, 2.0],
        ]))
        .unwrap();
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        // Choose x = 1 so b = A x is simple: [1, 0, 1].
        let b = Vector::<3>::new(black_box([1.0, 0.0, 1.0]));
        let x = ldlt.solve(b).unwrap().into_array();

        for &x_i in &x {
            assert_abs_diff_eq!(x_i, 1.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn singular_detected_for_degenerate_psd() {
        // Rank-1 Gram-like matrix.
        let a = Matrix::<2>::try_from_rows(black_box([[1.0, 1.0], [1.0, 1.0]])).unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::singular_numerical(
                1,
                FactorizationKind::Ldlt,
                0.0,
                DEFAULT_SINGULAR_TOL.get()
            )
        );
    }

    #[test]
    fn zero_pivot_with_nonzero_coupling_is_not_reported_as_singular() {
        let a = Matrix::<2>::try_from_rows(black_box([[0.0, 1.0], [1.0, 0.0]])).unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::not_positive_semidefinite_zero_coupling(0, 1, 1.0)
        );
    }

    #[test]
    fn zero_pivot_reports_non_finite_coupling_before_domain_violation() {
        let a = Matrix::<3>::try_from_rows(black_box([
            [1.0, 1.0, f64::MAX],
            [1.0, 1.0, -f64::MAX],
            [f64::MAX, -f64::MAX, 1.0],
        ]))
        .unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_matrix(ArithmeticOperation::LdltFactorization, 2, 1,)
        );
    }

    #[test]
    fn small_positive_pivot_reports_numerical_singularity() {
        let a = Matrix::<2>::try_from_rows(black_box([[1e-13, 0.0], [0.0, 1.0]])).unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::singular_numerical(
                0,
                FactorizationKind::Ldlt,
                1e-13,
                DEFAULT_SINGULAR_TOL.get()
            )
        );
    }

    #[test]
    fn small_positive_pivot_does_not_mask_earlier_non_finite_update() {
        let a = Matrix::<3>::try_from_rows(black_box([
            [1.0, 1.0, f64::MAX],
            [1.0, 1.0 + f64::EPSILON, -f64::MAX],
            [f64::MAX, -f64::MAX, 1.0],
        ]))
        .unwrap();

        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_matrix(ArithmeticOperation::LdltFactorization, 2, 1)
        );
    }

    #[test]
    fn negative_initial_diagonal_reports_not_positive_semidefinite() {
        let a = Matrix::<2>::try_from_rows(black_box([[-1.0, 0.0], [0.0, 1.0]])).unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(err, LaError::not_positive_semidefinite_negative(0, -1.0));
    }

    #[test]
    fn negative_updated_diagonal_reports_not_positive_semidefinite() {
        let a = Matrix::<2>::try_from_rows(black_box([[1.0, 2.0], [2.0, 1.0]])).unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(err, LaError::not_positive_semidefinite_negative(1, -3.0));
    }

    #[test]
    fn negative_pivot_does_not_mask_earlier_non_finite_update() {
        let a = Matrix::<3>::try_from_rows(black_box([
            [1.0, 2.0, f64::MAX],
            [2.0, 1.0, 0.0],
            [f64::MAX, 0.0, 1.0],
        ]))
        .unwrap();

        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_matrix(ArithmeticOperation::LdltFactorization, 2, 1,)
        );
    }

    #[test]
    fn non_finite_l_multiplier_overflow() {
        // d = 1e-11 > tol, but l = 1e300 / 1e-11 = 1e311 overflows f64.
        let a = Matrix::<2>::try_from_rows([[1e-11, 1e300], [1e300, 1.0]]).unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_matrix(ArithmeticOperation::LdltFactorization, 1, 0)
        );
    }

    #[test]
    fn non_finite_l_multiplier_overflow_fused_branch_6d() {
        // D > 5 uses the fused LDLT update path. Keep the same overflow shape
        // as the 2D test while forcing that branch.
        let mut rows = [[0.0; 6]; 6];
        for (i, row) in rows.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        rows[0][0] = 1e-11;
        rows[0][5] = 1e300;
        rows[5][0] = 1e300;

        let a = Matrix::<6>::try_from_rows(rows).unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_matrix(ArithmeticOperation::LdltFactorization, 5, 0)
        );
    }

    #[test]
    fn non_finite_trailing_submatrix_overflow() {
        // L multiplier is finite (1e200), but the rank-1 update
        // (-1e200 * 1.0) * 1e200 + 1.0 overflows.
        let a = Matrix::<2>::try_from_rows([[1.0, 1e200], [1e200, 1.0]]).unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_matrix(ArithmeticOperation::LdltFactorization, 1, 1)
        );
    }

    #[test]
    fn non_finite_trailing_submatrix_overflow_fused_branch_6d() {
        // D > 5 uses the fused LDLT update path. The overflowing trailing
        // diagonal is detected when it later becomes a pivot.
        let mut rows = [[0.0; 6]; 6];
        for (i, row) in rows.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        rows[0][5] = 1e200;
        rows[5][0] = 1e200;

        let a = Matrix::<6>::try_from_rows(rows).unwrap();
        let err = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_matrix(ArithmeticOperation::LdltFactorization, 5, 5)
        );
    }

    #[test]
    fn non_finite_solve_forward_substitution_overflow() {
        // SPD matrix with large L multiplier: L[1,0] = 1e153.
        // Forward substitution overflows: y[1] = 0 - 1e153 * 1e156 = -inf.
        let a = Matrix::<3>::try_from_rows([
            [1.0, 1e153, 0.0],
            [1e153, 1e306 + 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
        .unwrap();
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<3>::new([1e156, 0.0, 0.0]);
        let err = ldlt.solve(b).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_step(ArithmeticOperation::LdltSolve, 1)
        );
    }

    #[test]
    fn non_finite_solve_back_substitution_overflow() {
        // SPD matrix: [[1,0,0],[0,1,2],[0,2,5]] has LDLT factors
        // D=[1,1,1], L[2,1]=2.  Forward sub and diagonal solve produce
        // z=[0,0,1e308].  Back-substitution: x[2]=1e308 then
        // x[1] = 0 - 2*1e308 = -inf (overflows f64).
        let a = Matrix::<3>::try_from_rows([[1.0, 0.0, 0.0], [0.0, 1.0, 2.0], [0.0, 2.0, 5.0]])
            .unwrap();
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<3>::new([0.0, 0.0, 1e308]);
        let err = ldlt.solve(b).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_step(ArithmeticOperation::LdltSolve, 1)
        );
    }

    #[test]
    fn non_finite_solve_back_substitution_overflow_scatter_branch_5d() {
        // Exercises the D >= 5 row-prefix scatter branch with the same
        // bottom-right 2x2 SPD block used by the D3 back-substitution test.
        let a = Matrix::<5>::try_from_rows([
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 2.0],
            [0.0, 0.0, 0.0, 2.0, 5.0],
        ])
        .unwrap();
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<5>::new([0.0, 0.0, 0.0, 0.0, 1e308]);
        let err = ldlt.solve(b).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_step(ArithmeticOperation::LdltSolve, 3)
        );
    }

    #[test]
    fn non_finite_solve_diagonal_solve_overflow() {
        // Diagonal SPD matrix with a tiny diagonal entry just above the
        // singularity tolerance.  Forward substitution passes through the
        // large RHS unchanged, then the diagonal solve z[1] = y[1] / D[1]
        // = 1e300 / 1e-11 = 1e311 overflows f64, exercising the
        // `!v.is_finite()` branch of the diagonal solve.
        let a = Matrix::<2>::try_from_rows([[1.0, 0.0], [0.0, 1.0e-11]]).unwrap();
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

        let b = Vector::<2>::new([0.0, 1.0e300]);
        let err = ldlt.solve(b).unwrap_err();
        assert_eq!(
            err,
            LaError::non_finite_computation_step(ArithmeticOperation::LdltSolve, 1)
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
        let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();
        assert_eq!(
            ldlt.det(),
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

            let ldlt = Matrix::<4>::try_from_rows(rows)
                .unwrap()
                .ldlt(zero_tolerance)
                .unwrap();
            assert_eq!(ldlt.det(), Ok(1.0));
        }
    }

    #[test]
    fn det_balances_extreme_diagonals_in_large_dimension() {
        let zero_tolerance = Tolerance::try_new(0.0).unwrap();
        for diagonal in [
            [TWO_NEG_800, TWO_NEG_800, TWO_POS_800, TWO_POS_800, 1.0, 1.0],
            [TWO_POS_800, TWO_POS_800, TWO_NEG_800, TWO_NEG_800, 1.0, 1.0],
        ] {
            let mut rows = [[0.0; 6]; 6];
            for (i, value) in diagonal.into_iter().enumerate() {
                rows[i][i] = value;
            }

            let ldlt = Matrix::<6>::try_from_rows(rows)
                .unwrap()
                .ldlt(zero_tolerance)
                .unwrap();
            assert_eq!(ldlt.det(), Ok(1.0));
        }
    }

    #[test]
    fn det_rounds_final_tiny_magnitude_to_zero() {
        let zero_tolerance = Tolerance::try_new(0.0).unwrap();
        let matrix = Matrix::<2>::try_from_rows([[TWO_NEG_800, 0.0], [0.0, TWO_NEG_800]]).unwrap();
        let det = matrix.ldlt(zero_tolerance).unwrap().det().unwrap();

        assert_eq!(det.to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn ldlt_d1_classifies_positive_zero_and_negative_inputs() {
        let positive = Matrix::<1>::try_from_rows([[2.0]])
            .unwrap()
            .ldlt(DEFAULT_SINGULAR_TOL)
            .unwrap();
        assert_eq!(positive.det(), Ok(2.0));
        assert_abs_diff_eq!(
            positive
                .solve(Vector::<1>::new([6.0]))
                .unwrap()
                .into_array()[0],
            3.0,
            epsilon = 0.0
        );

        let zero = Matrix::<1>::try_from_rows([[0.0]])
            .unwrap()
            .ldlt(DEFAULT_SINGULAR_TOL);
        assert_eq!(
            zero,
            Err(LaError::singular_numerical(
                0,
                FactorizationKind::Ldlt,
                0.0,
                DEFAULT_SINGULAR_TOL.get()
            ))
        );

        let negative = Matrix::<1>::try_from_rows([[-1.0]])
            .unwrap()
            .ldlt(DEFAULT_SINGULAR_TOL);
        assert_eq!(
            negative,
            Err(LaError::not_positive_semidefinite_negative(0, -1.0))
        );
    }

    /// Construct an exactly representable tridiagonal SPD system from a unit
    /// lower-bidiagonal `L` and positive integer diagonal `D`.
    fn nontrivial_spd_system<const D: usize>() -> (Matrix<D>, Vector<D>, [f64; D], f64) {
        let mut rows = [[0.0_f64; D]; D];
        let mut expected_det = 1.0_f64;
        let mut diagonal = 1.0_f64;
        let mut k = 0;
        while k < D {
            rows[k][k] += diagonal;
            expected_det *= diagonal;
            if k + 1 < D {
                let off_diagonal = 0.5 * diagonal;
                rows[k][k + 1] += off_diagonal;
                rows[k + 1][k] += off_diagonal;
                rows[k + 1][k + 1] = 0.25_f64.mul_add(diagonal, rows[k + 1][k + 1]);
            }
            diagonal += 1.0;
            k += 1;
        }

        let mut expected_x = [0.0_f64; D];
        let mut value = 1.0_f64;
        for entry in &mut expected_x {
            *entry = value;
            value += 1.0;
        }
        let rhs = core::array::from_fn(|row| {
            rows[row]
                .iter()
                .zip(expected_x.iter())
                .fold(0.0_f64, |sum, (&coefficient, &x)| {
                    coefficient.mul_add(x, sum)
                })
        });

        (
            Matrix::<D>::try_from_rows(rows).unwrap(),
            Vector::<D>::try_new(rhs).unwrap(),
            expected_x,
            expected_det,
        )
    }

    macro_rules! gen_nontrivial_large_ldlt_tests {
        ($d:literal) => {
            paste! {
                #[test]
                fn [<ldlt_nontrivial_success_solve_and_det_agree_ $d d>]() {
                    let (matrix, rhs, expected_x, expected_det) =
                        nontrivial_spd_system::<$d>();
                    let ldlt = matrix.ldlt(DEFAULT_SINGULAR_TOL).unwrap();
                    let solution = ldlt.solve(rhs).unwrap().into_array();

                    for (actual, expected) in solution.into_iter().zip(expected_x) {
                        assert_abs_diff_eq!(actual, expected, epsilon = 1e-12);
                    }
                    assert_abs_diff_eq!(ldlt.det().unwrap(), expected_det, epsilon = 1e-12);
                    assert_abs_diff_eq!(matrix.det().unwrap(), expected_det, epsilon = 1e-10);
                }
            }
        };
    }

    gen_nontrivial_large_ldlt_tests!(6);
    gen_nontrivial_large_ldlt_tests!(8);

    #[test]
    fn asymmetric_input_returns_typed_error() {
        // a[0][1] = 2.0 but a[1][0] = -2.0 → clearly asymmetric.
        let a = Matrix::<3>::try_from_rows([[4.0, 2.0, 0.0], [-2.0, 5.0, 1.0], [0.0, 1.0, 3.0]])
            .unwrap();
        assert_eq!(
            a.ldlt(DEFAULT_SINGULAR_TOL),
            Err(LaError::asymmetric(0, 1, 3, 2.0, -2.0, 0.0))
        );
    }

    #[test]
    fn approximately_symmetric_input_is_rejected_before_factoring_another_operator() {
        // The tolerance-based diagnostic accepts this exact power-of-two
        // counterexample because 4 <= 1e-12 * 2^43. Factoring only its lower
        // triangle would instead replace the upper zero with 4: the original
        // determinant is 32, while that projected matrix has determinant 16.
        let matrix = Matrix::<2>::try_from_rows([[TWO_POS_43, 0.0], [4.0, TWO_NEG_38]]).unwrap();
        let diagnostic_tolerance = Tolerance::try_new(1e-12).unwrap();

        assert_eq!(matrix.det(), Ok(32.0));
        assert_eq!(matrix.is_symmetric(diagnostic_tolerance), Ok(true));
        assert_eq!(
            matrix.ldlt(DEFAULT_SINGULAR_TOL),
            Err(LaError::asymmetric(0, 1, 2, 0.0, 4.0, 0.0))
        );
    }

    // -----------------------------------------------------------------------
    // Const-evaluability tests.
    //
    // These prove that `Ldlt::det` and `Ldlt::solve` are truly `const fn`
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
                        let mut rows = [[0.0f64; $d]; $d];
                        let mut i = 0;
                        while i < $d {
                            rows[i][i] = 1.0;
                            i += 1;
                        }
                        rows[0][0] = 2.0;
                        let factors = LdltFactors::from_proven_rows(rows);
                        let ldlt = Ldlt::<$d> { factors };
                        ldlt.det()
                    };
                    assert_eq!(DET, Ok(2.0));
                }

                /// `Ldlt::solve` must be fully const-evaluable. Identity
                /// factors with RHS `b = [1.0, 2.0, …, D]` round-trips `b`
                /// unchanged, exercising the full forward sub / diagonal solve
                /// / back sub pipeline inside a `const { … }` initializer.
                #[test]
                fn [<ldlt_solve_const_eval_ $d d>]() {
                    #[expect(
                        clippy::cast_precision_loss,
                        reason = "test indices are at most five and exactly representable as f64"
                    )]
                    const X: Result<Vector<$d>, LaError> = {
                        let factors = LdltFactors::from_proven_rows(
                            Matrix::<$d>::identity().into_rows()
                        );
                        let ldlt = Ldlt::<$d> { factors };
                        let mut b_arr = [0.0f64; $d];
                        let mut i = 0;
                        while i < $d {
                            b_arr[i] = i as f64 + 1.0;
                            i += 1;
                        }
                        let b = Vector::<$d>::new(b_arr);
                        ldlt.solve(b)
                    };
                    let x = X.unwrap().into_array();
                    #[expect(
                        clippy::cast_precision_loss,
                        reason = "test indices are at most five and exactly representable as f64"
                    )]
                    for i in 0..$d {
                        let expected = i as f64 + 1.0;
                        assert!((x[i] - expected).abs() <= 1e-12);
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
