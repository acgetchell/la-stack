//! Property-based tests for LU/LDLT factorization APIs.
//!
//! These tests construct matrices from known factors so we have a reliable oracle for
//! determinant and solve behavior.

use approx::assert_abs_diff_eq;
use pastey::paste;
use proptest::prelude::*;

use la_stack::prelude::*;

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

macro_rules! gen_factorization_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<ldlt_det_and_solve_match_constructed_factors_ $d d>](
                    l_raw in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_factor_entry()),
                    ),
                    d_diag in proptest::array::[<uniform $d>](positive_diag_entry()),
                    x_true in proptest::array::[<uniform $d>](small_f64()),
                ) {
                    // Construct A = L * diag(D) * L^T, where L is unit-lower-triangular.
                    let mut l = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        for j in 0..$d {
                            l[i][j] = if i == j {
                                1.0
                            } else if i > j {
                                l_raw[i][j]
                            } else {
                                0.0
                            };
                        }
                    }

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

                    let mut b_arr = [0.0f64; $d];
                    for i in 0..$d {
                        let mut sum = 0.0;
                        for j in 0..$d {
                            sum = a_rows[i][j].mul_add(x_true[j], sum);
                        }
                        b_arr[i] = sum;
                    }

                    let a = Matrix::<$d>::from_rows(a_rows);
                    let ldlt = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap();

                    assert_abs_diff_eq!(ldlt.det(), expected_det, epsilon = 1e-8);

                    let b = Vector::<$d>::new(b_arr);
                    let x = ldlt.solve_vec(b).unwrap().into_array();
                    for i in 0..$d {
                        assert_abs_diff_eq!(x[i], x_true[i], epsilon = 1e-8);
                    }
                }

                #[test]
                fn [<lu_det_and_solve_match_constructed_factors_no_perm_ $d d>](
                    l_raw in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_factor_entry()),
                    ),
                    u_raw in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_factor_entry()),
                    ),
                    u_diag in proptest::array::[<uniform $d>](nonzero_diag_entry()),
                    x_true in proptest::array::[<uniform $d>](small_f64()),
                ) {
                    // Construct A = L * U, where L is unit-lower-triangular and U is upper-triangular.
                    let mut l = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        for j in 0..$d {
                            l[i][j] = if i == j {
                                1.0
                            } else if i > j {
                                l_raw[i][j]
                            } else {
                                0.0
                            };
                        }
                    }

                    let mut u = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        for j in 0..$d {
                            u[i][j] = if i == j {
                                u_diag[i]
                            } else if i < j {
                                u_raw[i][j]
                            } else {
                                0.0
                            };
                        }
                    }

                    let mut a_rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        for j in 0..$d {
                            let mut sum = 0.0;
                            // L[i][k] is zero for k > i; U[k][j] is zero for k > j.
                            let k_max = if i < j { i } else { j };
                            for k in 0..=k_max {
                                sum = l[i][k].mul_add(u[k][j], sum);
                            }
                            a_rows[i][j] = sum;
                        }
                    }

                    let expected_det = {
                        let mut acc = 1.0;
                        for i in 0..$d {
                            acc *= u_diag[i];
                        }
                        acc
                    };

                    let mut b_arr = [0.0f64; $d];
                    for i in 0..$d {
                        let mut sum = 0.0;
                        for j in 0..$d {
                            sum = a_rows[i][j].mul_add(x_true[j], sum);
                        }
                        b_arr[i] = sum;
                    }

                    let a = Matrix::<$d>::from_rows(a_rows);
                    let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

                    assert_abs_diff_eq!(lu.det(), expected_det, epsilon = 1e-8);

                    let b = Vector::<$d>::new(b_arr);
                    let x = lu.solve_vec(b).unwrap().into_array();
                    for i in 0..$d {
                        assert_abs_diff_eq!(x[i], x_true[i], epsilon = 1e-8);
                    }
                }

                #[test]
                fn [<lu_det_and_solve_match_constructed_factors_row_swap_ $d d>](
                    l_raw in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_factor_entry()),
                    ),
                    u_raw in proptest::array::[<uniform $d>](
                        proptest::array::[<uniform $d>](small_factor_entry()),
                    ),
                    u_diag in proptest::array::[<uniform $d>](nonzero_diag_entry()),
                    x_true in proptest::array::[<uniform $d>](small_f64()),
                ) {
                    // Construct A = P^{-1} * L * U, where P swaps the first two rows.
                    // This ensures det(A) has an extra sign flip vs det(LU).
                    prop_assume!($d >= 2);

                    let mut l = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        for j in 0..$d {
                            l[i][j] = if i == j {
                                1.0
                            } else if i > j {
                                l_raw[i][j]
                            } else {
                                0.0
                            };
                        }
                    }

                    let mut u = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        for j in 0..$d {
                            u[i][j] = if i == j {
                                u_diag[i]
                            } else if i < j {
                                u_raw[i][j]
                            } else {
                                0.0
                            };
                        }
                    }

                    let mut lu_rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        for j in 0..$d {
                            let mut sum = 0.0;
                            let k_max = if i < j { i } else { j };
                            for k in 0..=k_max {
                                sum = l[i][k].mul_add(u[k][j], sum);
                            }
                            lu_rows[i][j] = sum;
                        }
                    }

                    // Apply P^{-1}: swap rows 0 and 1.
                    let mut a_rows = lu_rows;
                    a_rows.swap(0, 1);

                    let expected_det = {
                        let mut acc = 1.0;
                        for i in 0..$d {
                            acc *= u_diag[i];
                        }
                        -acc
                    };

                    let mut b_arr = [0.0f64; $d];
                    for i in 0..$d {
                        let mut sum = 0.0;
                        for j in 0..$d {
                            sum = a_rows[i][j].mul_add(x_true[j], sum);
                        }
                        b_arr[i] = sum;
                    }

                    let a = Matrix::<$d>::from_rows(a_rows);
                    let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();

                    assert_abs_diff_eq!(lu.det(), expected_det, epsilon = 1e-8);

                    let b = Vector::<$d>::new(b_arr);
                    let x = lu.solve_vec(b).unwrap().into_array();
                    for i in 0..$d {
                        assert_abs_diff_eq!(x[i], x_true[i], epsilon = 1e-8);
                    }
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
