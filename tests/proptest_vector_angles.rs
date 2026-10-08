#![forbid(unsafe_code)]

//! Exact Gram-determinant oracles, symmetry, and exact rescaling properties.

use core::f64::consts::PI;

use la_stack::prelude::*;
use pastey::paste;
use proptest::{array, prelude::*};

#[path = "common/proptest_config.rs"]
mod proptest_config;
use proptest_config::with_default_cases;

#[cfg(feature = "exact")]
mod rational_reference {
    use core::f64::consts::{FRAC_PI_2, PI};

    use num_bigint::BigInt;
    use num_rational::BigRational;
    use num_traits::{FromPrimitive, One, Signed, ToPrimitive, Zero};

    /// Approximate sqrt(q) for 0 < q <= 1 without rounding q to binary64.
    /// Integer square root gives a lower bound with relative error < 2^-127
    /// before the final binary64 conversion, even when q itself underflows.
    fn sqrt_fraction(value: &BigRational) -> f64 {
        assert!(value.is_positive() && value <= &BigRational::one());
        let numerator = value.numer().magnitude();
        let denominator = value.denom().magnitude();
        let exponent =
            i64::try_from(numerator.bits()).unwrap() - i64::try_from(denominator.bits()).unwrap();
        let shift = usize::try_from(128 - exponent.div_euclid(2)).unwrap();
        let square = (numerator << (2 * shift)) / denominator;
        let root = square.sqrt();
        BigRational::new(BigInt::from(root), BigInt::one() << shift)
            .to_f64()
            .unwrap()
    }

    /// Exact rational Gram products are independent of the production minors,
    /// FMA compensation, direction scaling, and hypot reduction. The square root
    /// uses the bounded approximation above; the final binary64 angle evaluation,
    /// including quadrant reconstruction, is rounded.
    pub(super) fn angle<const D: usize>(left: &[f64; D], right: &[f64; D]) -> f64 {
        let mut dot = BigRational::zero();
        let mut left_squared = BigRational::zero();
        let mut right_squared = BigRational::zero();
        for (&x, &y) in left.iter().zip(right) {
            let x = BigRational::from_f64(x).unwrap();
            let y = BigRational::from_f64(y).unwrap();
            dot += &x * &y;
            left_squared += &x * &x;
            right_squared += &y * &y;
        }
        assert!(left_squared.is_positive() && right_squared.is_positive());
        if dot.is_zero() {
            return FRAC_PI_2;
        }
        let dot_squared = &dot * &dot;
        let gram = left_squared * right_squared - &dot_squared;
        assert!(!gram.is_negative());
        if gram.is_zero() {
            return if dot.is_negative() { PI } else { 0.0 };
        }
        // Use the reciprocal for angles above pi/4 to keep the root in [0, 1].
        let reciprocal = gram > dot_squared;
        let squared_ratio = if reciprocal {
            dot_squared / gram
        } else {
            gram / dot_squared
        };
        let ratio = sqrt_fraction(&squared_ratio);
        // Below 2^-30, |atan(r)-r| <= r^3/3 < r*2^-61. This independent
        // reference also avoids a platform atan losing subnormal outputs.
        let acute = if ratio < 2.0_f64.powi(-30) {
            ratio
        } else {
            ratio.atan()
        };
        let acute = if reciprocal { FRAC_PI_2 - acute } else { acute };
        if dot.is_negative() { PI - acute } else { acute }
    }
}

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
                    let actual = a.as_slice().angle(&b).unwrap();
                    prop_assert!((actual - expected).abs() <= 16.0 * f64::EPSILON);
                    prop_assert_eq!(actual.to_bits(), b.as_slice().angle(&a).unwrap().to_bits());
                    let unscaled = left.map(f64::from).as_slice().angle(&right.map(f64::from)).unwrap();
                    prop_assert_eq!(actual.to_bits(), unscaled.to_bits());
                    prop_assert_eq!(a.as_slice().angle(&a).unwrap().to_bits(), 0);
                    prop_assert_eq!(a.as_slice().angle(&a.map(|x| -x)).unwrap().to_bits(), PI.to_bits());
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
                    let angle = a.as_slice().angle(&b).unwrap();
                    prop_assert!(angle.is_finite() && (0.0..=PI).contains(&angle));
                    prop_assert_eq!(angle.to_bits(), b.as_slice().angle(&a).unwrap().to_bits());
                    prop_assert_eq!(a.as_slice().angle(&a).unwrap().to_bits(), 0);
                    prop_assert_eq!(a.as_slice().angle(&a.map(|x| -x)).unwrap().to_bits(), PI.to_bits());
                    prop_assert!((angle + a.as_slice().angle(&b.map(|x| -x)).unwrap() - PI).abs() <= 4.0 * f64::EPSILON);
                    #[cfg(feature = "exact")]
                    {
                        let expected = rational_reference::angle(&a, &b);
                        let tolerance = expected.mul_add(32.0 * f64::EPSILON, f64::from_bits(1));
                        prop_assert!((angle - expected).abs() <= tolerance,
                            "{a:?}, {b:?}: {angle:e}, expected {expected:e}");
                    }
                }

                #[cfg(feature = "exact")]
                #[test]
                fn [<angle_nearly_parallel_rational_oracle_ $d d>](
                    left in array::[<uniform $d>](0.5_f64..2.0),
                    index in 0_usize..$d,
                    left_exponent in -1000_i32..=1000,
                    right_exponent in -1000_i32..=1000,
                ) {
                    let mut right = left;
                    right[index] = right[index].next_up();
                    let a = left.map(|x| x * 2.0_f64.powi(left_exponent));
                    let b = right.map(|x| x * 2.0_f64.powi(right_exponent));
                    let expected = rational_reference::angle(&a, &b);
                    let actual = a.as_slice().angle(&b).unwrap();
                    prop_assert!(expected > 0.0);
                    prop_assert!((actual - expected).abs() <= expected * (32.0 * f64::EPSILON),
                        "{a:?}, {b:?}: {actual:e}, expected {expected:e}");
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

#[cfg(feature = "exact")]
#[test]
fn rational_reference_and_api_match_analytical_extremes() {
    let tiny = f64::from_bits(1);
    let n = 2.0_f64.powi(51);
    for (left, right, expected) in [
        ([1.0, 0.0, 0.0], [1.0, tiny, 0.0], tiny),
        ([3.0, tiny, 0.0], [3.0, -tiny, 0.0], tiny),
        ([1.0, 1.0, 0.0], [1.0, 1.0, tiny], tiny),
        ([1.0, 0.0, 0.0], [1.0, 3.0 * tiny, 4.0 * tiny], 5.0 * tiny),
        ([n, n - 1.0, 0.0], [n + 1.0, n, 0.0], 2.0_f64.powi(-103)),
        ([f64::MAX, 0.0, 0.0], [f64::MAX, tiny, 0.0], 0.0),
        ([f64::MAX; 3], [tiny; 3], 0.0),
        ([f64::MAX; 3], [-tiny; 3], PI),
        ([1.0, 0.0, 0.0], [0.0, -f64::MAX, 0.0], PI / 2.0),
        ([1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], 3.0 * PI / 4.0),
    ] {
        for actual in [
            rational_reference::angle(&left, &right),
            left.as_slice().angle(&right).unwrap(),
        ] {
            assert!(
                (actual - expected).abs() <= expected * (4.0 * f64::EPSILON),
                "{left:?}, {right:?}: {actual:e}, expected {expected:e}"
            );
        }
    }
}
