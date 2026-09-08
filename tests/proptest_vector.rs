#![forbid(unsafe_code)]

//! Property-based tests for the `Vector` public API.

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

macro_rules! gen_vector_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(with_default_cases(64))]

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
                fn [<vector_dot_and_squared_norm_match_integer_oracles_ $d d>](
                    left in array::[<uniform $d>](-1000i16..=1000),
                    right in array::[<uniform $d>](-1000i16..=1000),
                ) {
                    // At D<=8 these integer products and sums fit in i32 and
                    // are exact in binary64, independently of the FMA kernel.
                    let dot: i32 = left.iter().zip(&right)
                        .map(|(&a, &b)| i32::from(a) * i32::from(b))
                        .sum();
                    let squared_norm: i32 = left.iter()
                        .map(|&value| i32::from(value).pow(2))
                        .sum();
                    let a = Vector::<$d>::try_new(left.map(f64::from)).unwrap();
                    let b = Vector::<$d>::try_new(right.map(f64::from)).unwrap();

                    prop_assert_eq!(a.dot(&b), Ok(f64::from(dot)));
                    prop_assert_eq!(a.norm_squared(), Ok(f64::from(squared_norm)));
                }

                #[test]
                fn [<vector_dot_commutes_and_norm_squared_matches_dot_self_ $d d>](
                    a_arr in array::[<uniform $d>](small_f64()),
                    b_arr in array::[<uniform $d>](small_f64()),
                ) {
                    let a = Vector::<$d>::try_new(a_arr).unwrap();
                    let b = Vector::<$d>::try_new(b_arr).unwrap();

                    let dot_ab = a.dot(&b).unwrap();
                    let dot_reversed = b.dot(&a).unwrap();
                    assert_abs_diff_eq!(dot_ab, dot_reversed, epsilon = 1e-14);

                    let dot_aa = a.dot(&a).unwrap();
                    let norm_squared = a.norm_squared().unwrap();
                    assert_abs_diff_eq!(norm_squared, dot_aa, epsilon = 0.0);

                    // Squared norm is always non-negative for finite inputs.
                    prop_assert!(norm_squared >= 0.0);

                    // The safe norm agrees with an independently ordered chain
                    // of binary64 hypot operations on this moderate domain.
                    let hypot_norm = a_arr
                        .iter()
                        .fold(0.0f64, |accumulator, &value| accumulator.hypot(value));
                    let norm = a.norm().unwrap();
                    assert_abs_diff_eq!(norm, hypot_norm, epsilon = 1e-12);
                    prop_assert!(norm >= 0.0 && norm.is_finite());

                    let negated = Vector::<$d>::try_new(a_arr.map(|value| -value)).unwrap();
                    prop_assert_eq!(norm.to_bits(), negated.norm().unwrap().to_bits());

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
gen_vector_proptests!(6);
gen_vector_proptests!(7);
gen_vector_proptests!(8);

#[test]
fn vector_norm_zero_dimension_is_always_positive_zero() {
    let norm = Vector::<0>::zero().norm().unwrap();
    assert_eq!(norm.to_bits(), 0.0f64.to_bits());
}
