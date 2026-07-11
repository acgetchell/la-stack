#![forbid(unsafe_code)]

//! Smoke tests for the deterministic inputs used by the `vs_linalg` benchmark.

#![cfg(feature = "bench")]

use faer::linalg::solvers::Solve;
use faer::mat::AsMatRef;
use faer::perm::PermRef;
use faer::{Mat, Side};
use nalgebra::{Const, DimMin, SMatrix, SVector};

use la_stack::{DEFAULT_SINGULAR_TOL, Matrix, Vector};

#[path = "../benches/common/vs_linalg.rs"]
pub mod vs_linalg_common;

use vs_linalg_common::{
    PreparedFaerLuDet, faer_det_from_ldlt, faer_perm_sign, la_stack_dot, la_stack_tolerance,
    make_balanced_dynamic_range_rows, make_ill_conditioned_matrix_rows, make_matrix_rows,
    make_pivoting_matrix_rows, make_vector_array, matrix_entry, nalgebra_inf_norm, vector_entry,
};

/// Assert scalar agreement with a tolerance that scales for larger magnitudes.
fn assert_close(label: &str, actual: f64, expected: f64) {
    assert!(
        actual.is_finite() && expected.is_finite(),
        "{label}: comparison requires finite values, actual={actual:?}, expected={expected:?}",
    );
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

/// Check LU determinant and solve agreement for one benchmark dimension.
fn assert_lu_agreement<const D: usize>()
where
    Const<D>: DimMin<Const<D>, Output = Const<D>>,
{
    let a = Matrix::<D>::try_from_rows(make_matrix_rows::<D>())
        .unwrap_or_else(|err| panic!("la_stack matrix construction failed: {err}"));
    let rhs = Vector::<D>::try_new(make_vector_array::<D>(0.0))
        .unwrap_or_else(|err| panic!("la_stack RHS vector construction failed: {err}"));
    let na = SMatrix::<f64, D, D>::from_fn(matrix_entry::<D>);
    let nrhs = SVector::<f64, D>::from_fn(|i, _| vector_entry(i, 0.0));
    let fa = Mat::<f64>::from_fn(D, D, matrix_entry::<D>);
    let frhs = Mat::<f64>::from_fn(D, 1, |i, _| vector_entry(i, 0.0));

    let la_lu = a
        .lu(DEFAULT_SINGULAR_TOL)
        .unwrap_or_else(|err| panic!("la_stack LU factorization failed: {err}"));
    let na_lu = na.lu();
    let fa_lu = fa.partial_piv_lu();

    let la_lu_det = la_lu
        .det()
        .unwrap_or_else(|err| panic!("la_stack LU determinant failed: {err}"));
    let la_matrix_det = a
        .det()
        .unwrap_or_else(|err| panic!("la_stack Matrix determinant failed: {err}"));
    assert_close("la_stack_det", la_matrix_det, la_lu_det);
    assert_close("nalgebra_det_from_lu", na_lu.determinant(), la_lu_det);
    assert_close(
        "faer_det_from_lu",
        PreparedFaerLuDet::new(&fa_lu).det(),
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
}

/// Check LDLT/Cholesky determinant and solve agreement for one dimension.
fn assert_ldlt_agreement<const D: usize>() {
    let a = Matrix::<D>::try_from_rows(make_matrix_rows::<D>())
        .unwrap_or_else(|err| panic!("la_stack matrix construction failed: {err}"));
    let rhs = Vector::<D>::try_new(make_vector_array::<D>(0.0))
        .unwrap_or_else(|err| panic!("la_stack RHS vector construction failed: {err}"));
    let na = SMatrix::<f64, D, D>::from_fn(matrix_entry::<D>);
    let nrhs = SVector::<f64, D>::from_fn(|i, _| vector_entry(i, 0.0));
    let fa = Mat::<f64>::from_fn(D, D, matrix_entry::<D>);
    let frhs = Mat::<f64>::from_fn(D, 1, |i, _| vector_entry(i, 0.0));

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
}

/// Check vector dot-product and squared-norm agreement for one benchmark dimension.
fn assert_vector_operation_agreement<const D: usize>() {
    let v1 = Vector::<D>::try_new(make_vector_array::<D>(0.0))
        .unwrap_or_else(|err| panic!("la_stack vector construction failed: {err}"));
    let v2 = Vector::<D>::try_new(make_vector_array::<D>(1.0))
        .unwrap_or_else(|err| panic!("la_stack vector construction failed: {err}"));
    let nv1 = SVector::<f64, D>::from_fn(|i, _| vector_entry(i, 0.0));
    let nv2 = SVector::<f64, D>::from_fn(|i, _| vector_entry(i, 1.0));
    let fv1 = Mat::<f64>::from_fn(D, 1, |i, _| vector_entry(i, 0.0));
    let fv2 = Mat::<f64>::from_fn(D, 1, |i, _| vector_entry(i, 1.0));

    let la_dot = la_stack_dot(&v1, &v2).unwrap_or_else(|err| panic!("la_stack dot failed: {err}"));
    assert_close("nalgebra_dot", nv1.dot(&nv2), la_dot);
    let mut fa_dot = 0.0;
    for i in 0..D {
        fa_dot = fv1[(i, 0)].mul_add(fv2[(i, 0)], fa_dot);
    }
    assert_close("faer_dot", fa_dot, la_dot);

    let la_norm2_sq = v1
        .norm2_sq()
        .unwrap_or_else(|err| panic!("la_stack norm2_sq failed: {err}"));
    assert_close("nalgebra_norm_squared", nv1.norm_squared(), la_norm2_sq);
    let fa_norm2_sq = fv1.as_mat_ref().squared_norm_l2();
    assert_close("faer_norm2_sq", fa_norm2_sq, la_norm2_sq);
}

#[test]
fn faer_lu_determinant_includes_odd_row_permutation_sign() {
    let matrix = Mat::<f64>::from_fn(2, 2, |row, col| [[0.0, 2.0], [3.0, 4.0]][row][col]);
    let lu = matrix.partial_piv_lu();

    assert_close(
        "faer odd row-permutation sign",
        faer_perm_sign(lu.P()),
        -1.0,
    );
    assert_close(
        "faer determinant with one pivot swap",
        PreparedFaerLuDet::new(&lu).det(),
        -6.0,
    );
}

#[test]
fn scalar_agreement_rejects_non_finite_values() {
    for (actual, expected) in [
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, -1.0),
        (f64::NAN, 1.0),
        (1.0, f64::NAN),
    ] {
        assert!(
            std::panic::catch_unwind(|| assert_close("non-finite regression", actual, expected))
                .is_err()
        );
    }
}

#[test]
fn faer_permutation_sign_handles_valid_cycle_parities() {
    let empty = PermRef::new_checked(&[], &[], 0);
    let identity = PermRef::new_checked(&[0, 1, 2], &[0, 1, 2], 3);
    let transposition = PermRef::new_checked(&[1, 0], &[1, 0], 2);
    let three_cycle = PermRef::new_checked(&[1, 2, 0], &[2, 0, 1], 3);

    assert_close("empty permutation sign", faer_perm_sign(empty), 1.0);
    assert_close("identity permutation sign", faer_perm_sign(identity), 1.0);
    assert_close(
        "transposition permutation sign",
        faer_perm_sign(transposition),
        -1.0,
    );
    assert_close(
        "three-cycle permutation sign",
        faer_perm_sign(three_cycle),
        1.0,
    );
}

#[test]
fn faer_permutation_sign_handles_large_permutations_without_allocation() {
    let mut forward: [usize; 129] = std::array::from_fn(|index| index);
    forward.swap(127, 128);
    let permutation = PermRef::new_checked(&forward, &forward, forward.len());

    assert_close(
        "large transposition permutation sign",
        faer_perm_sign(permutation),
        -1.0,
    );
}

#[test]
fn ill_conditioned_fixture_is_fixed_positive_definite_d8() {
    let rows = make_ill_conditioned_matrix_rows();

    for (row_index, row) in rows.iter().enumerate() {
        for (col_index, &value) in row.iter().enumerate() {
            if row_index == col_index {
                assert!(value.is_normal() && value.is_sign_positive());
            } else {
                assert_eq!(value.to_bits(), 0.0f64.to_bits());
            }
        }
    }
    assert_eq!(
        rows[7][7].to_bits(),
        f64::from_bits(911_u64 << 52).to_bits()
    );
}

#[test]
fn stress_inputs_exercise_pivoting_conditioning_and_scaled_products() {
    let zero_tolerance = la_stack_tolerance(0.0).unwrap();

    let pivoting_rows = make_pivoting_matrix_rows::<8>();
    assert!(pivoting_rows[1][0].abs() > pivoting_rows[0][0].abs());
    let pivoting_lu = Matrix::<8>::try_from_rows(pivoting_rows)
        .unwrap()
        .lu(zero_tolerance)
        .unwrap();
    assert!(pivoting_lu.det().unwrap().is_finite());

    let ill_conditioned = Matrix::<8>::try_from_rows(make_ill_conditioned_matrix_rows()).unwrap();
    let expected_ill_conditioned_det = f64::from_bits((1023_u64 - 448) << 52);
    assert_eq!(
        ill_conditioned
            .lu(zero_tolerance)
            .unwrap()
            .det()
            .unwrap()
            .to_bits(),
        expected_ill_conditioned_det.to_bits()
    );
    assert_eq!(
        ill_conditioned
            .ldlt(zero_tolerance)
            .unwrap()
            .det()
            .unwrap()
            .to_bits(),
        expected_ill_conditioned_det.to_bits()
    );

    #[cfg(not(la_stack_v0_4_3_api))]
    {
        let balanced = Matrix::<8>::try_from_rows(make_balanced_dynamic_range_rows()).unwrap();
        assert_eq!(balanced.lu(zero_tolerance).unwrap().det(), Ok(1.0));
        assert_eq!(balanced.ldlt(zero_tolerance).unwrap().det(), Ok(1.0));
    }
}

/// Check matrix infinity-norm agreement for one benchmark dimension.
fn assert_matrix_inf_norm_agreement<const D: usize>() {
    let a = Matrix::<D>::try_from_rows(make_matrix_rows::<D>())
        .unwrap_or_else(|err| panic!("la_stack matrix construction failed: {err}"));
    let na = SMatrix::<f64, D, D>::from_fn(matrix_entry::<D>);
    let fa = Mat::<f64>::from_fn(D, D, matrix_entry::<D>);

    let la_norm = a
        .inf_norm()
        .unwrap_or_else(|err| panic!("la_stack inf_norm failed: {err}"));
    assert_close("nalgebra_inf_norm", nalgebra_inf_norm(&na), la_norm);
    let mut fa_norm = 0.0;
    for r in 0..D {
        let mut row_sum = 0.0;
        for c in 0..D {
            row_sum += fa[(r, c)].abs();
        }
        if row_sum > fa_norm {
            fa_norm = row_sum;
        }
    }
    assert_close("faer_inf_norm", fa_norm, la_norm);
}

/// Generate focused cross-crate smoke tests for a benchmark dimension.
macro_rules! gen_smoke_tests {
    ($d:literal, $lu:ident, $ldlt:ident, $vector:ident, $norm:ident) => {
        #[test]
        fn $lu() {
            assert_lu_agreement::<$d>();
        }

        #[test]
        fn $ldlt() {
            assert_ldlt_agreement::<$d>();
        }

        #[test]
        fn $vector() {
            assert_vector_operation_agreement::<$d>();
        }

        #[test]
        fn $norm() {
            assert_matrix_inf_norm_agreement::<$d>();
        }
    };
}

gen_smoke_tests!(
    2,
    vs_linalg_lu_agrees_2d,
    vs_linalg_ldlt_agrees_2d,
    vs_linalg_vector_operations_agree_2d,
    vs_linalg_matrix_inf_norm_agrees_2d
);
gen_smoke_tests!(
    3,
    vs_linalg_lu_agrees_3d,
    vs_linalg_ldlt_agrees_3d,
    vs_linalg_vector_operations_agree_3d,
    vs_linalg_matrix_inf_norm_agrees_3d
);
gen_smoke_tests!(
    4,
    vs_linalg_lu_agrees_4d,
    vs_linalg_ldlt_agrees_4d,
    vs_linalg_vector_operations_agree_4d,
    vs_linalg_matrix_inf_norm_agrees_4d
);
gen_smoke_tests!(
    5,
    vs_linalg_lu_agrees_5d,
    vs_linalg_ldlt_agrees_5d,
    vs_linalg_vector_operations_agree_5d,
    vs_linalg_matrix_inf_norm_agrees_5d
);
gen_smoke_tests!(
    8,
    vs_linalg_lu_agrees_8d,
    vs_linalg_ldlt_agrees_8d,
    vs_linalg_vector_operations_agree_8d,
    vs_linalg_matrix_inf_norm_agrees_8d
);
gen_smoke_tests!(
    16,
    vs_linalg_lu_agrees_16d,
    vs_linalg_ldlt_agrees_16d,
    vs_linalg_vector_operations_agree_16d,
    vs_linalg_matrix_inf_norm_agrees_16d
);
gen_smoke_tests!(
    32,
    vs_linalg_lu_agrees_32d,
    vs_linalg_ldlt_agrees_32d,
    vs_linalg_vector_operations_agree_32d,
    vs_linalg_matrix_inf_norm_agrees_32d
);
gen_smoke_tests!(
    64,
    vs_linalg_lu_agrees_64d,
    vs_linalg_ldlt_agrees_64d,
    vs_linalg_vector_operations_agree_64d,
    vs_linalg_matrix_inf_norm_agrees_64d
);
