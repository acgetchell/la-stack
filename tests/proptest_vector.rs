//! Property-based tests for the `Vector` public API.

use approx::assert_abs_diff_eq;
use pastey::paste;
use proptest::prelude::*;

use la_stack::prelude::*;

fn small_f64() -> impl Strategy<Value = f64> {
    (-1000i16..=1000i16).prop_map(|x| f64::from(x) / 10.0)
}

macro_rules! gen_public_api_vector_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<vector_new_as_array_into_array_roundtrip_ $d d>](
                    arr in proptest::array::[<uniform $d>](small_f64()),
                ) {
                    let v = Vector::<$d>::new(arr);

                    for i in 0..$d {
                        assert_abs_diff_eq!(v.as_array()[i], arr[i], epsilon = 0.0);
                    }

                    let out = v.into_array();
                    for i in 0..$d {
                        assert_abs_diff_eq!(out[i], arr[i], epsilon = 0.0);
                    }
                }

                #[test]
                fn [<vector_dot_commutes_and_norm2_sq_matches_dot_self_ $d d>](
                    a_arr in proptest::array::[<uniform $d>](small_f64()),
                    b_arr in proptest::array::[<uniform $d>](small_f64()),
                ) {
                    let a = Vector::<$d>::new(a_arr);
                    let b = Vector::<$d>::new(b_arr);

                    let dot_ab = a.dot(b);
                    let dot_ba = b.dot(a);
                    assert_abs_diff_eq!(dot_ab, dot_ba, epsilon = 0.0);

                    let dot_aa = a.dot(a);
                    assert_abs_diff_eq!(a.norm2_sq(), dot_aa, epsilon = 0.0);

                    // Squared norm is always non-negative for finite inputs.
                    prop_assert!(a.norm2_sq() >= 0.0);

                    // Dot with zero vector is zero.
                    let z = Vector::<$d>::zero();
                    assert_abs_diff_eq!(a.dot(z), 0.0, epsilon = 0.0);
                }
            }
        }
    };
}

// Mirror delaunay-style multi-dimension tests.
gen_public_api_vector_proptests!(2);
gen_public_api_vector_proptests!(3);
gen_public_api_vector_proptests!(4);
gen_public_api_vector_proptests!(5);
