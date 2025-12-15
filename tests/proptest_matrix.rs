//! Property-based tests for the `Matrix` public API.

use approx::assert_abs_diff_eq;
use pastey::paste;
use proptest::prelude::*;

use la_stack::prelude::*;

fn small_f64() -> impl Strategy<Value = f64> {
    (-1000i16..=1000i16).prop_map(|x| f64::from(x) / 10.0)
}

fn small_nonzero_f64() -> impl Strategy<Value = f64> {
    prop_oneof![(-1000i16..=-1i16), (1i16..=1000i16)].prop_map(|x| f64::from(x) / 10.0)
}

macro_rules! gen_public_api_matrix_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<matrix_from_rows_get_roundtrip_ $d d>](
                    rows in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_f64()),
                    ),
                ) {
                    let m = Matrix::<$d>::from_rows(rows);

                    for r in 0..$d {
                        for c in 0..$d {
                            assert_abs_diff_eq!(m.get(r, c).unwrap(), rows[r][c], epsilon = 0.0);
                        }
                    }

                    // Out-of-bounds is None.
                    prop_assert_eq!(m.get($d, 0), None);
                    prop_assert_eq!(m.get(0, $d), None);
                }

                #[test]
                fn [<matrix_set_get_in_bounds_ $d d>](
                    r in 0usize..$d,
                    c in 0usize..$d,
                    v in small_f64(),
                ) {
                    let mut m = Matrix::<$d>::zero();
                    prop_assert!(m.set(r, c, v));
                    assert_abs_diff_eq!(m.get(r, c).unwrap(), v, epsilon = 0.0);
                }

                #[test]
                fn [<matrix_inf_norm_matches_max_abs_row_sum_ $d d>](
                    rows in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_f64()),
                    ),
                ) {
                    let m = Matrix::<$d>::from_rows(rows);

                    let expected = rows
                        .iter()
                        .map(|row| row.iter().map(|&x| x.abs()).sum::<f64>())
                        .fold(0.0f64, f64::max);

                    assert_abs_diff_eq!(m.inf_norm(), expected, epsilon = 0.0);
                    prop_assert!(m.inf_norm() >= 0.0);
                }

                #[test]
                fn [<matrix_det_and_solve_vec_for_diagonal_ $d d>](
                    diag in proptest::array::[<uniform $d>](small_nonzero_f64()),
                    b_arr in proptest::array::[<uniform $d>](small_f64()),
                ) {
                    // Diagonal matrix: det is product of diagonal, and solve is element-wise division.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = diag[i];
                    }
                    let a = Matrix::<$d>::from_rows(rows);

                    let det = a.det(DEFAULT_PIVOT_TOL).unwrap();
                    let expected_det = {
                        let mut acc = 1.0;
                        for i in 0..$d {
                            acc *= diag[i];
                        }
                        acc
                    };
                    assert_abs_diff_eq!(det, expected_det, epsilon = 1e-12);

                    let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();
                    let b = Vector::<$d>::new(b_arr);
                    let x = lu.solve_vec(b).unwrap().into_array();

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
gen_public_api_matrix_proptests!(2);
gen_public_api_matrix_proptests!(3);
gen_public_api_matrix_proptests!(4);
gen_public_api_matrix_proptests!(5);
