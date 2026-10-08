#![forbid(unsafe_code)]

//! Integer Gram-determinant oracles, symmetry, and exact rescaling properties.

use core::f64::consts::PI;

use la_stack::prelude::*;
use pastey::paste;
use proptest::{array, prelude::*};

#[path = "common/proptest_config.rs"]
mod proptest_config;
use proptest_config::with_default_cases;

/// A 2x2 Gram determinant gives the squared sine times squared norm products.
/// Integer arithmetic makes its cancellation exact. This does not enumerate
/// exterior minors, scale directions, or use compensated floating arithmetic.
#[expect(
    clippy::cast_precision_loss,
    reason = "D <= 6 and |coordinate| <= 1000 keep the Gram determinant below 2^53"
)]
fn reference<const D: usize>(left: &[i16; D], right: &[i16; D]) -> f64 {
    let mut dot = 0_i32;
    let mut left_squared = 0_i64;
    let mut right_squared = 0_i64;
    for i in 0..D {
        dot += i32::from(left[i]) * i32::from(right[i]);
        left_squared += i64::from(left[i]).pow(2);
        right_squared += i64::from(right[i]).pow(2);
    }
    let gram = left_squared * right_squared - i64::from(dot).pow(2);
    (gram as f64).sqrt().atan2(f64::from(dot))
}

macro_rules! gen_angle_properties {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(with_default_cases(256))]

                #[test]
                fn [<angle_integer_oracle_and_rescaling_ $d d>](
                    left in array::[<uniform $d>](-1000_i16..=1000),
                    right in array::[<uniform $d>](-1000_i16..=1000),
                    left_exponent in -1000_i32..=1000,
                    right_exponent in -1000_i32..=1000,
                ) {
                    prop_assume!(left.iter().any(|&x| x != 0));
                    prop_assume!(right.iter().any(|&x| x != 0));
                    let expected = reference(&left, &right);
                    let a = left.map(|x| f64::from(x) * 2.0_f64.powi(left_exponent));
                    let b = right.map(|x| f64::from(x) * 2.0_f64.powi(right_exponent));
                    let actual = angle_between(&a, &b).unwrap();
                    prop_assert!((actual - expected).abs() <= 16.0 * f64::EPSILON);
                    prop_assert_eq!(actual.to_bits(), angle_between(&b, &a).unwrap().to_bits());
                    let unscaled = angle_between(&left.map(f64::from), &right.map(f64::from)).unwrap();
                    prop_assert_eq!(actual.to_bits(), unscaled.to_bits());
                    prop_assert_eq!(angle_between(&a, &a).unwrap().to_bits(), 0);
                    prop_assert_eq!(angle_between(&a, &a.map(|x| -x)).unwrap().to_bits(), PI.to_bits());
                }

                #[test]
                fn [<angle_full_binary64_range_ $d d>](
                    left in array::[<uniform $d>](any::<u64>()),
                    right in array::[<uniform $d>](any::<u64>()),
                ) {
                    let a = left.map(f64::from_bits);
                    let b = right.map(f64::from_bits);
                    prop_assume!(a.iter().all(|x| x.is_finite()) && a.iter().any(|&x| x != 0.0));
                    prop_assume!(b.iter().all(|x| x.is_finite()) && b.iter().any(|&x| x != 0.0));
                    let angle = angle_between(&a, &b).unwrap();
                    prop_assert!(angle.is_finite() && (0.0..=PI).contains(&angle));
                    prop_assert_eq!(angle.to_bits(), angle_between(&b, &a).unwrap().to_bits());
                    prop_assert_eq!(angle_between(&a, &a).unwrap().to_bits(), 0);
                    prop_assert_eq!(angle_between(&a, &a.map(|x| -x)).unwrap().to_bits(), PI.to_bits());
                    prop_assert!((angle + angle_between(&a, &b.map(|x| -x)).unwrap() - PI).abs() <= 4.0 * f64::EPSILON);
                }
            }
        }
    };
}

gen_angle_properties!(2);
gen_angle_properties!(3);
gen_angle_properties!(4);
gen_angle_properties!(5);
gen_angle_properties!(6);
