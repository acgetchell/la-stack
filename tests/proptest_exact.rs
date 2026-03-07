//! Property-based tests for the `det_sign_exact` API (requires `exact` feature).

#![cfg(feature = "exact")]

use pastey::paste;
use proptest::prelude::*;

use la_stack::prelude::*;

fn small_nonzero_f64() -> impl Strategy<Value = f64> {
    prop_oneof![(-1000i16..=-1i16), (1i16..=1000i16)].prop_map(|x| f64::from(x) / 10.0)
}

macro_rules! gen_det_sign_exact_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<det_sign_exact_agrees_with_det_for_diagonal_ $d d>](
                    diag in proptest::array::[<uniform $d>](small_nonzero_f64()),
                ) {
                    // Diagonal matrix: determinant sign = product of diagonal signs.
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = diag[i];
                    }
                    let m = Matrix::<$d>::from_rows(rows);

                    let exact_sign = m.det_sign_exact();

                    // Expected sign from the product of diagonal entries.
                    let neg_count = diag.iter().filter(|&&x| x < 0.0).count();
                    let expected_sign: i8 = if neg_count % 2 == 0 { 1 } else { -1 };

                    prop_assert_eq!(exact_sign, expected_sign);
                }

                #[test]
                fn [<det_sign_exact_agrees_with_det_signum_ $d d>](
                    diag in proptest::array::[<uniform $d>](small_nonzero_f64()),
                ) {
                    // For well-conditioned diagonal matrices, det().signum()
                    // should agree with det_sign_exact().
                    let mut rows = [[0.0f64; $d]; $d];
                    for i in 0..$d {
                        rows[i][i] = diag[i];
                    }
                    let m = Matrix::<$d>::from_rows(rows);

                    let exact_sign = m.det_sign_exact();
                    let fp_det = m.det(DEFAULT_PIVOT_TOL).unwrap();
                    let fp_sign: i8 = if fp_det > 0.0 {
                        1
                    } else if fp_det < 0.0 {
                        -1
                    } else {
                        0
                    };

                    prop_assert_eq!(exact_sign, fp_sign);
                }
            }
        }
    };
}

gen_det_sign_exact_proptests!(2);
gen_det_sign_exact_proptests!(3);
gen_det_sign_exact_proptests!(4);
gen_det_sign_exact_proptests!(5);
