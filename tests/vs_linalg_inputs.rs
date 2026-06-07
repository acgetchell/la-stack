//! Smoke tests for the deterministic inputs used by the `vs_linalg` benchmark.

#![cfg(feature = "bench")]

use faer::linalg::solvers::Solve;
use faer::{Mat, Side};
use nalgebra::{SMatrix, SVector};

use la_stack::{DEFAULT_SINGULAR_TOL, Matrix, Vector};

#[path = "../benches/common/vs_linalg.rs"]
mod vs_linalg;

use vs_linalg::{
    faer_det_from_ldlt, faer_det_from_partial_piv_lu, make_matrix_rows, make_vector_array,
    matrix_entry, nalgebra_inf_norm, vector_entry,
};

/// Assert scalar agreement with a tolerance that scales for larger magnitudes.
fn assert_close(label: &str, actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    let diff = (actual - expected).abs();
    assert!(
        diff <= 1.0e-9 * scale,
        "{label}: actual={actual:?}, expected={expected:?}, diff={diff:?}, scale={scale:?}",
    );
}

/// Assert componentwise agreement for fixed-size vector results from peer crates.
fn assert_vector_close<const D: usize>(label: &str, actual: [f64; D], expected: [f64; D]) {
    for i in 0..D {
        assert_close(&format!("{label}[{i}]"), actual[i], expected[i]);
    }
}

/// Copy a faer single-column solve result into an array for componentwise checks.
fn faer_column_to_array<const D: usize>(m: &Mat<f64>) -> [f64; D] {
    let mut data = [0.0; D];
    for i in 0..D {
        data[i] = m[(i, 0)];
    }
    data
}

/// Copy a nalgebra stack vector into an array for componentwise checks.
fn nalgebra_vector_to_array<const D: usize>(v: &SVector<f64, D>) -> [f64; D] {
    let mut data = [0.0; D];
    data.copy_from_slice(v.as_slice());
    data
}

/// Generate one cross-crate smoke test for a concrete benchmark dimension.
macro_rules! gen_smoke_test {
    ($name:ident, $d:literal) => {
        #[test]
        #[allow(clippy::too_many_lines)]
        fn $name() {
            let a = Matrix::<$d>::try_from_rows(make_matrix_rows::<$d>())
                .unwrap_or_else(|err| panic!("la_stack matrix construction failed: {err}"));
            let rhs = Vector::<$d>::try_new(make_vector_array::<$d>(0.0))
                .unwrap_or_else(|err| panic!("la_stack RHS vector construction failed: {err}"));
            let v1 = Vector::<$d>::try_new(make_vector_array::<$d>(0.0))
                .unwrap_or_else(|err| panic!("la_stack vector construction failed: {err}"));
            let v2 = Vector::<$d>::try_new(make_vector_array::<$d>(1.0))
                .unwrap_or_else(|err| panic!("la_stack vector construction failed: {err}"));

            let na = SMatrix::<f64, $d, $d>::from_fn(|r, c| matrix_entry::<$d>(r, c));
            let nrhs = SVector::<f64, $d>::from_fn(|i, _| vector_entry(i, 0.0));
            let nv1 = SVector::<f64, $d>::from_fn(|i, _| vector_entry(i, 0.0));
            let nv2 = SVector::<f64, $d>::from_fn(|i, _| vector_entry(i, 1.0));

            let fa = Mat::<f64>::from_fn($d, $d, |r, c| matrix_entry::<$d>(r, c));
            let frhs = Mat::<f64>::from_fn($d, 1, |i, _| vector_entry(i, 0.0));
            let fv1 = Mat::<f64>::from_fn($d, 1, |i, _| vector_entry(i, 0.0));
            let fv2 = Mat::<f64>::from_fn($d, 1, |i, _| vector_entry(i, 1.0));

            let la_lu = a
                .lu(DEFAULT_SINGULAR_TOL)
                .unwrap_or_else(|err| panic!("la_stack LU factorization failed: {err}"));
            let na_lu = na.lu();
            let fa_lu = fa.partial_piv_lu();

            let la_lu_det = la_lu
                .det()
                .unwrap_or_else(|err| panic!("la_stack LU determinant failed: {err}"));
            assert_close("nalgebra_det_from_lu", na_lu.determinant(), la_lu_det);
            assert_close(
                "faer_det_from_lu",
                faer_det_from_partial_piv_lu(&fa_lu),
                la_lu_det,
            );

            let la_lu_x = la_lu
                .solve(rhs)
                .unwrap_or_else(|err| panic!("la_stack LU solve failed: {err}"));
            let na_lu_x = na_lu
                .solve(&nrhs)
                .unwrap_or_else(|| panic!("nalgebra LU solve returned no result"));
            let fa_lu_x = fa_lu.solve(&frhs);
            let la_lu_x = la_lu_x.into_array();
            assert_vector_close(
                "nalgebra_solve_from_lu",
                nalgebra_vector_to_array(&na_lu_x),
                la_lu_x,
            );
            assert_vector_close(
                "faer_solve_from_lu",
                faer_column_to_array(&fa_lu_x),
                la_lu_x,
            );

            let la_ldlt = a
                .ldlt(DEFAULT_SINGULAR_TOL)
                .unwrap_or_else(|err| panic!("la_stack LDLT factorization failed: {err}"));
            let na_cholesky = na
                .cholesky()
                .unwrap_or_else(|| panic!("nalgebra Cholesky factorization returned no result"));
            let fa_ldlt = fa
                .ldlt(Side::Lower)
                .unwrap_or_else(|err| panic!("faer LDLT factorization failed: {err}"));

            let la_ldlt_det = la_ldlt
                .det()
                .unwrap_or_else(|err| panic!("la_stack LDLT determinant failed: {err}"));
            assert_close(
                "nalgebra_det_from_cholesky",
                na_cholesky.determinant(),
                la_ldlt_det,
            );
            assert_close(
                "faer_det_from_ldlt",
                faer_det_from_ldlt(&fa_ldlt),
                la_ldlt_det,
            );

            let la_ldlt_x = la_ldlt
                .solve(rhs)
                .unwrap_or_else(|err| panic!("la_stack LDLT solve failed: {err}"));
            let na_cholesky_x = na_cholesky.solve(&nrhs);
            let fa_ldlt_x = fa_ldlt.solve(&frhs);
            let la_ldlt_x = la_ldlt_x.into_array();
            assert_vector_close(
                "nalgebra_solve_from_cholesky",
                nalgebra_vector_to_array(&na_cholesky_x),
                la_ldlt_x,
            );
            assert_vector_close(
                "faer_solve_from_ldlt",
                faer_column_to_array(&fa_ldlt_x),
                la_ldlt_x,
            );

            let la_dot = v1
                .dot(v2)
                .unwrap_or_else(|err| panic!("la_stack dot failed: {err}"));
            assert_close("nalgebra_dot", nv1.dot(&nv2), la_dot);
            let mut fa_dot = 0.0;
            for i in 0..$d {
                fa_dot = fv1[(i, 0)].mul_add(fv2[(i, 0)], fa_dot);
            }
            assert_close("faer_dot", fa_dot, la_dot);

            let la_norm = a
                .inf_norm()
                .unwrap_or_else(|err| panic!("la_stack inf_norm failed: {err}"));
            assert_close("nalgebra_inf_norm", nalgebra_inf_norm(&na), la_norm);
            let mut fa_norm = 0.0;
            for r in 0..$d {
                let mut row_sum = 0.0;
                for c in 0..$d {
                    row_sum += fa[(r, c)].abs();
                }
                if row_sum > fa_norm {
                    fa_norm = row_sum;
                }
            }
            assert_close("faer_inf_norm", fa_norm, la_norm);
        }
    };
}

gen_smoke_test!(vs_linalg_shared_inputs_agree_2d, 2);
gen_smoke_test!(vs_linalg_shared_inputs_agree_3d, 3);
gen_smoke_test!(vs_linalg_shared_inputs_agree_4d, 4);
gen_smoke_test!(vs_linalg_shared_inputs_agree_5d, 5);
