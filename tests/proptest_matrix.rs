#![forbid(unsafe_code)]

//! Property-based tests for the `Matrix` public API.

use approx::assert_abs_diff_eq;
use pastey::paste;
use proptest::{array, prelude::*};

use la_stack::prelude::*;

fn small_f64() -> impl Strategy<Value = f64> {
    (-1000i16..=1000i16).prop_map(|x| f64::from(x) / 10.0)
}

fn small_nonzero_f64() -> impl Strategy<Value = f64> {
    prop_oneof![(-1000i16..=-1i16), (1i16..=1000i16)].prop_map(|x| f64::from(x) / 10.0)
}

macro_rules! gen_matrix_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<matrix_try_from_rows_views_and_get_roundtrip_ $d d>](
                    rows in array::[<uniform $d>](
                        array::[<uniform $d>](small_f64()),
                    ),
                ) {
                    let m = Matrix::<$d>::try_from_rows(rows).unwrap();
                    prop_assert_eq!(m.as_rows(), &rows);

                    for r in 0..$d {
                        for c in 0..$d {
                            assert_abs_diff_eq!(m.get(r, c).unwrap(), rows[r][c], epsilon = 0.0);
                            assert_abs_diff_eq!(m.try_get(r, c).unwrap(), rows[r][c], epsilon = 0.0);
                        }
                    }

                    // Out-of-bounds is None.
                    prop_assert_eq!(m.get($d, 0), None);
                    prop_assert_eq!(m.get(0, $d), None);
                    let row_out_of_bounds = matches!(
                        m.try_get($d, 0),
                        Err(LaError::IndexOutOfBounds {
                            row,
                            col,
                            dim,
                            ..
                        })
                            if row == $d && col == 0 && dim == $d
                    );
                    prop_assert!(row_out_of_bounds);
                    let col_out_of_bounds = matches!(
                        m.try_get(0, $d),
                        Err(LaError::IndexOutOfBounds {
                            row,
                            col,
                            dim,
                            ..
                        })
                            if row == 0 && col == $d && dim == $d
                    );
                    prop_assert!(col_out_of_bounds);

                    prop_assert_eq!(m.into_rows(), rows);
                }

                #[test]
                fn [<matrix_set_get_in_bounds_ $d d>](
                    r in 0usize..$d,
                    c in 0usize..$d,
                    v in small_f64(),
                ) {
                    let mut m = Matrix::<$d>::zero();
                    prop_assert_eq!(m.set(r, c, v), Ok(()));
                    assert_abs_diff_eq!(m.get(r, c).unwrap(), v, epsilon = 0.0);
                    prop_assert_eq!(m.set(r, c, -v), Ok(()));
                    assert_abs_diff_eq!(m.try_get(r, c).unwrap(), -v, epsilon = 0.0);
                }

                #[test]
                fn [<matrix_set_out_of_bounds_preserves_matrix_ $d d>](
                    rows in array::[<uniform $d>](
                        array::[<uniform $d>](small_f64()),
                    ),
                    v in small_f64(),
                ) {
                    let mut m = Matrix::<$d>::try_from_rows(rows).unwrap();
                    let original = m;

                    let row_out_of_bounds = matches!(
                        m.set($d, 0, v),
                        Err(LaError::IndexOutOfBounds {
                            row,
                            col,
                            dim,
                            ..
                        })
                            if row == $d && col == 0 && dim == $d
                    );
                    prop_assert!(row_out_of_bounds);
                    prop_assert_eq!(m, original);
                    let col_out_of_bounds = matches!(
                        m.set(0, $d, v),
                        Err(LaError::IndexOutOfBounds {
                            row,
                            col,
                            dim,
                            ..
                        })
                            if row == 0 && col == $d && dim == $d
                    );
                    prop_assert!(col_out_of_bounds);
                    prop_assert_eq!(m, original);
                }

                #[test]
                fn [<matrix_inf_norm_matches_max_abs_row_sum_ $d d>](
                    rows in array::[<uniform $d>](
                        array::[<uniform $d>](small_f64()),
                    ),
                ) {
                    let m = Matrix::<$d>::try_from_rows(rows).unwrap();

                    let expected = rows
                        .iter()
                        .map(|row| row.iter().map(|&x| x.abs()).sum::<f64>())
                        .fold(0.0f64, f64::max);

                    let actual = m.inf_norm().unwrap();
                    assert_abs_diff_eq!(actual, expected, epsilon = 0.0);
                    prop_assert!(actual >= 0.0);
                }

                #[test]
                fn [<matrix_det_and_solve_for_diagonal_ $d d>](
                    diag in array::[<uniform $d>](small_nonzero_f64()),
                    b_arr in array::[<uniform $d>](small_f64()),
                ) {
                    // Diagonal matrix: det is product of diagonal, and solve is element-wise division.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = diag[i];
                    }
                    let a = Matrix::<$d>::try_from_rows(rows).unwrap();

                    let det = a.det().unwrap();
                    let expected_det = {
                        let mut acc = 1.0;
                        for i in 0..$d {
                            acc *= diag[i];
                        }
                        acc
                    };
                    // The closed-form and LU paths evaluate the diagonal
                    // product in different orders, so we allow a few ULPs of
                    // relative error (floor 1e-12 for near-zero determinants).
                    let eps = expected_det.abs().mul_add(1e-12, 1e-12);
                    assert_abs_diff_eq!(det, expected_det, epsilon = eps);

                    let lu = a.lu(DEFAULT_SINGULAR_TOL).unwrap();
                    let b = Vector::<$d>::try_new(b_arr).unwrap();
                    let x = lu.solve(b).unwrap().into_array();

                    for i in 0..$d {
                        let expected_x = b_arr[i] / diag[i];
                        assert_abs_diff_eq!(x[i], expected_x, epsilon = 1e-12);
                    }
                }

            }
        }
    };
}

// Mirror delaunay-style multi-dimension tests.
gen_matrix_proptests!(2);
gen_matrix_proptests!(3);
gen_matrix_proptests!(4);
gen_matrix_proptests!(5);
