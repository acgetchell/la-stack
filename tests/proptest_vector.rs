#![forbid(unsafe_code)]

//! Property-based tests for the `Vector` public API.

use approx::assert_abs_diff_eq;
use pastey::paste;
use proptest::{array, prelude::*};

use la_stack::prelude::*;

fn small_f64() -> impl Strategy<Value = f64> {
    (-1000i16..=1000i16).prop_map(|x| f64::from(x) / 10.0)
}

macro_rules! gen_vector_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(64))]

                #[test]
                fn [<vector_try_new_as_array_into_array_roundtrip_ $d d>](
                    arr in array::[<uniform $d>](small_f64()),
                ) {
                    let v = Vector::<$d>::try_new(arr).unwrap();

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
                    a_arr in array::[<uniform $d>](small_f64()),
                    b_arr in array::[<uniform $d>](small_f64()),
                ) {
                    let a = Vector::<$d>::try_new(a_arr).unwrap();
                    let b = Vector::<$d>::try_new(b_arr).unwrap();

                    let dot_ab = a.dot(&b).unwrap();
                    let dot_reversed = b.dot(&a).unwrap();
                    assert_abs_diff_eq!(dot_ab, dot_reversed, epsilon = 1e-14);

                    let dot_aa = a.dot(&a).unwrap();
                    let norm2_sq = a.norm2_sq().unwrap();
                    assert_abs_diff_eq!(norm2_sq, dot_aa, epsilon = 0.0);

                    // Squared norm is always non-negative for finite inputs.
                    prop_assert!(norm2_sq >= 0.0);

                    // Dot with zero vector is zero.
                    let z = Vector::<$d>::zero();
                    assert_abs_diff_eq!(a.dot(&z).unwrap(), 0.0, epsilon = 1e-14);
                }
            }
        }
    };
}

// Mirror delaunay-style multi-dimension tests.
gen_vector_proptests!(1);
gen_vector_proptests!(2);
gen_vector_proptests!(3);
gen_vector_proptests!(4);
gen_vector_proptests!(5);

#[test]
fn zero_dimension_vector_obeys_empty_sum_contracts() {
    let vector = Vector::<0>::try_new([]).unwrap();

    assert!(vector.as_array().is_empty());
    assert!(vector.into_array().is_empty());
    assert_eq!(vector.dot(&Vector::zero()), Ok(0.0));
    assert_eq!(vector.norm2_sq(), Ok(0.0));
}
