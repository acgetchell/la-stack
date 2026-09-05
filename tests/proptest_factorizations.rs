#![forbid(unsafe_code)]

//! Property-based tests for LU/LDLT factorization APIs.
//!
//! These tests construct matrices from known factors so we have a reliable oracle for
//! determinant and solve behavior.

use core::{array::from_fn, cmp::Ordering};

use approx::assert_abs_diff_eq;
use pastey::paste;
use proptest::{array, prelude::*};

use la_stack::prelude::*;

#[path = "common/proptest_config.rs"]
mod proptest_config;
use proptest_config::with_default_cases;

fn small_f64() -> impl Strategy<Value = f64> {
    (-1000i16..=1000i16).prop_map(|x| f64::from(x) / 10.0)
}

fn small_factor_entry() -> impl Strategy<Value = f64> {
    // Keep entries small so constructed matrices are reasonably conditioned.
    (-50i16..=50i16).prop_map(|x| f64::from(x) / 100.0)
}

fn positive_diag_entry() -> impl Strategy<Value = f64> {
    // Strictly positive diagonal, comfortably above DEFAULT_SINGULAR_TOL.
    (1i16..=20i16).prop_map(|x| f64::from(x) / 10.0)
}

fn nonzero_diag_entry() -> impl Strategy<Value = f64> {
    // Strictly non-zero diagonal with a margin from 0.
    prop_oneof![(-20i16..=-1i16), (1i16..=20i16)].prop_map(|x| f64::from(x) / 10.0)
}

fn unit_lower<const D: usize>(raw: &[[f64; D]; D]) -> [[f64; D]; D] {
    from_fn(|i| {
        from_fn(|j| match i.cmp(&j) {
            Ordering::Equal => 1.0,
            Ordering::Greater => raw[i][j],
            Ordering::Less => 0.0,
        })
    })
}

/// Construct A = L * U with unit-lower L and the supplied diagonal of U.
fn lu_product<const D: usize>(
    l_raw: &[[f64; D]; D],
    u_raw: &[[f64; D]; D],
    u_diag: &[f64; D],
) -> [[f64; D]; D] {
    let l = unit_lower(l_raw);
    let u: [[f64; D]; D] = from_fn(|i| {
        from_fn(|j| match i.cmp(&j) {
            Ordering::Equal => u_diag[i],
            Ordering::Less => u_raw[i][j],
            Ordering::Greater => 0.0,
        })
    });
    from_fn(|i| {
        from_fn(|j| {
            // L[i][k] is zero for k > i; U[k][j] is zero for k > j.
            (0..=i.min(j)).fold(0.0, |sum, k| l[i][k].mul_add(u[k][j], sum))
        })
    })
}

fn matvec<const D: usize>(rows: &[[f64; D]; D], vector: &[f64; D]) -> [f64; D] {
    from_fn(|i| {
        rows[i]
            .iter()
            .zip(vector)
            .fold(0.0, |sum, (&coefficient, &value)| {
                coefficient.mul_add(value, sum)
            })
    })
}

fn check_lu<const D: usize>(rows: &[[f64; D]; D], expected_det: f64, x_true: &[f64; D]) {
    let b_arr = matvec(rows, x_true);
    let a = Matrix::try_from_rows(*rows).unwrap();
    let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();
    assert_abs_diff_eq!(lu.det().unwrap(), expected_det, epsilon = 1e-8);

    let b = Vector::try_new(b_arr).unwrap();
    let x = lu.solve(b).unwrap().into_array();
    for (actual, expected) in x.iter().zip(x_true) {
        assert_abs_diff_eq!(actual, expected, epsilon = 1e-8);
    }
}

macro_rules! gen_factorization_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(with_default_cases(64))]

                #[test]
                fn [<ldlt_det_and_solve_match_constructed_factors_ $d d>](
                    l_raw in array::[<uniform $d>](
                        array::[<uniform $d>](small_factor_entry()),
                    ),
                    d_diag in array::[<uniform $d>](positive_diag_entry()),
                    x_true in array::[<uniform $d>](small_f64()),
                ) {
                    // Construct A = L * diag(D) * L^T, where L is unit-lower-triangular.
                    let l = unit_lower(&l_raw);

                    let mut a_rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        for j in 0..=i {
                            let mut sum = 0.0;
                            // L[j][k] is zero for k > j.
                            for k in 0..=j {
                                sum = (l[i][k] * d_diag[k]).mul_add(l[j][k], sum);
                            }
                            a_rows[i][j] = sum;
                            a_rows[j][i] = sum;
                        }
                    }

                    let expected_det = {
                        let mut acc = 1.0;
                        for i in 0..$d {
                            acc *= d_diag[i];
                        }
                        acc
                    };

                    let b_arr = matvec(&a_rows, &x_true);

                    let a = Matrix::<$d>::try_from_rows(a_rows).unwrap();
                    let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

                    let det_ldlt = ldlt.det().unwrap();
                    let det_lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap().det().unwrap();
                    assert_abs_diff_eq!(det_ldlt, expected_det, epsilon = 1e-8);
                    assert_abs_diff_eq!(det_ldlt, det_lu, epsilon = 1e-8);

                    let b = Vector::<$d>::try_new(b_arr).unwrap();
                    let x = ldlt.solve(b).unwrap().into_array();
                    for i in 0..$d {
                        assert_abs_diff_eq!(x[i], x_true[i], epsilon = 1e-8);
                    }
                }

                #[test]
                fn [<lu_det_and_solve_match_constructed_factors_no_perm_ $d d>](
                    l_raw in array::[<uniform $d>](
                        array::[<uniform $d>](small_factor_entry()),
                    ),
                    u_raw in array::[<uniform $d>](
                        array::[<uniform $d>](small_factor_entry()),
                    ),
                    u_diag in array::[<uniform $d>](nonzero_diag_entry()),
                    x_true in array::[<uniform $d>](small_f64()),
                ) {
                    let a_rows = lu_product(&l_raw, &u_raw, &u_diag);
                    let expected_det = u_diag.iter().fold(1.0, |acc, &value| acc * value);
                    check_lu(&a_rows, expected_det, &x_true);
                }

                #[test]
                fn [<lu_det_and_solve_match_constructed_factors_row_swap_ $d d>](
                    l_raw in array::[<uniform $d>](
                        array::[<uniform $d>](small_factor_entry()),
                    ),
                    u_raw in array::[<uniform $d>](
                        array::[<uniform $d>](small_factor_entry()),
                    ),
                    u_diag in array::[<uniform $d>](nonzero_diag_entry()),
                    x_true in array::[<uniform $d>](small_f64()),
                ) {
                    // A = P^{-1} * L * U: swapping rows 0 and 1 flips det(LU)'s sign.
                    let mut a_rows = lu_product(&l_raw, &u_raw, &u_diag);
                    a_rows.swap(0, 1);
                    let expected_det = -u_diag.iter().fold(1.0, |acc, &value| acc * value);
                    check_lu(&a_rows, expected_det, &x_true);
                }
            }
        }
    };
}

// Mirror delaunay-style multi-dimension tests.
gen_factorization_proptests!(2);
gen_factorization_proptests!(3);
gen_factorization_proptests!(4);
gen_factorization_proptests!(5);
// Exercise the D > 5 factorization and solve branches with randomized inputs.
gen_factorization_proptests!(8);
