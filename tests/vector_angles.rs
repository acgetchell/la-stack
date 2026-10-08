#![forbid(unsafe_code)]

//! Independent analytical references and downstream error/allocation contracts.

use core::assert_matches;
use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};
use std::hint::black_box;

use la_stack::prelude::*;
use pastey::paste;

fn check<const D: usize>(left: [f64; D], right: [f64; D], expected: f64, tolerance: f64) {
    let actual = angle_between(&left, &right).unwrap();
    assert!(actual.is_finite() && (0.0..=PI).contains(&actual));
    assert!(
        (actual - expected).abs() <= tolerance,
        "{left:?}, {right:?}: {actual:e}, expected {expected:e} ± {tolerance:e}"
    );
    assert_eq!(
        actual.to_bits(),
        angle_between(&right, &left).unwrap().to_bits()
    );
    let left = Vector::try_new(left).unwrap();
    let right = Vector::try_new(right).unwrap();
    assert_eq!(actual.to_bits(), left.angle(&right).unwrap().to_bits());
    if expected == 0.0 {
        assert_eq!(actual.to_bits(), 0);
    }
}

fn planar<const D: usize>(x: f64, y: f64) -> [f64; D] {
    let mut result = [0.0; D];
    result[0] = x;
    result[D - 1] = y;
    result
}

fn known_angles<const D: usize>() {
    let axis = planar::<D>(1.0, 0.0);
    let t = f64::from_bits(0x3e10_0000_0000_0000); // Exactly 2^-30.
    check(axis, axis, 0.0, 0.0);
    check([1.0; D], [1.0; D], 0.0, 0.0);
    check(planar::<D>(1.0, 2.0), planar(3.0, 6.0), 0.0, 0.0);
    check(axis, planar(-1.0, -0.0), PI, 0.0);
    check([1.0; D], [-1.0; D], PI, 0.0);
    check(axis, planar(-0.0, -1.0), FRAC_PI_2, 2.0 * f64::EPSILON);
    check(axis, planar(1.0, 1.0), FRAC_PI_4, 2.0 * f64::EPSILON);
    check(axis, planar(1.0, t), t.atan(), t * 4.0 * f64::EPSILON);
    check(axis, planar(-1.0, t), PI - t.atan(), 2.0 * f64::EPSILON);

    // Non-axis directions exercise cancellation in both weighted coordinates.
    // Their exact planar determinant is t and dot product is 2+t.
    let diagonal_angle = (t / (2.0 + t)).atan();
    check(
        planar::<D>(1.0, 1.0),
        planar(1.0, 1.0 + t),
        diagonal_angle,
        2.0 * f64::EPSILON,
    );
    check(
        planar::<D>(1.0, 1.0),
        planar(-1.0, -1.0 - t),
        PI - diagonal_angle,
        2.0 * f64::EPSILON,
    );
    let n = 2.0_f64.powi(51);
    let cancellation_angle = 2.0_f64.powi(-103).atan();
    check(
        planar::<D>(n, n - 1.0),
        planar(n + 1.0, n),
        cancellation_angle,
        cancellation_angle * 4.0 * f64::EPSILON,
    );

    // The integer cross and dot are formed exactly before conversion; no
    // normalization, norm weighting, or half-angle identity is used here.
    for (a, b, c, d) in [(3_i32, 4, -4, 3), (12, -5, 7, 24), (999, 1000, 1000, 1001)] {
        let expected = f64::from((a * d - b * c).abs()).atan2(f64::from(a * c + b * d));
        check(
            planar::<D>(f64::from(a), f64::from(b)),
            planar(f64::from(c), f64::from(d)),
            expected,
            8.0 * f64::EPSILON,
        );
    }
}

fn extreme_angles<const D: usize>() {
    let tiny = f64::from_bits(1);
    for scale in [tiny, f64::MIN_POSITIVE, 1e-300, 1.0, 1e300, f64::MAX] {
        check(
            planar::<D>(scale, 0.0),
            planar(scale, scale),
            FRAC_PI_4,
            2.0 * f64::EPSILON,
        );
        check([scale; D], [scale; D], 0.0, 0.0);
        check([scale; D], [-scale; D], PI, 0.0);
    }
    check(
        planar::<D>(f64::MAX, 0.0),
        planar(tiny, tiny),
        FRAC_PI_4,
        2.0 * f64::EPSILON,
    );
    check(
        planar::<D>(f64::MAX, tiny),
        planar(tiny, f64::MAX),
        FRAC_PI_2,
        2.0 * f64::EPSILON,
    );
    check(
        planar::<D>(f64::MAX, 0.0),
        planar(f64::MAX, f64::MAX * 2.0_f64.powi(-30)),
        2.0_f64.powi(-30).atan(),
        1e-24,
    );
    // The norm overflows, but the angle is still well-defined.
    check([f64::MAX; D], [tiny; D], 0.0, 0.0);
}

fn subnormal_angles<const D: usize>() {
    for bits in [1, 2, 3, 4, 5, 7, 31, 1023, (1_u64 << 51) - 1] {
        let tiny = f64::from_bits(bits);
        for sign in [-1.0, 1.0] {
            check(planar::<D>(1.0, -0.0), planar(1.0, sign * tiny), tiny, 0.0);
        }
        // Each normalized small component is only half the final angle.
        check(planar::<D>(2.0, tiny), planar(2.0, -tiny), tiny, 0.0);
        // Dividing each coordinate by 3 first also loses information. The
        // analytical angle is 2*atan(bits*2^-1074 / 3); its cubic correction
        // is far below half an ulp, so integer rounding gives the reference.
        let expected = f64::from_bits((2 * bits + 1) / 3);
        check(planar::<D>(3.0, tiny), planar(3.0, -tiny), expected, 0.0);
    }
    // Subnormal input coordinates can also describe a large ordinary angle.
    let tiny = f64::from_bits(1);
    check(
        planar::<D>(tiny, tiny),
        planar(tiny, -tiny),
        FRAC_PI_2,
        2.0 * f64::EPSILON,
    );
    // A truly unrepresentable angle rounds to zero.
    check(planar::<D>(f64::MAX, 0.0), planar(f64::MAX, tiny), 0.0, 0.0);
}

fn small_angle_transition<const D: usize>() {
    let boundary = 2.0_f64.powi(-27);
    for slope in [
        boundary.next_down(),
        boundary,
        boundary.next_up(),
        2.0 * boundary,
        2.0_f64.powi(-20),
    ] {
        // The alternating series bounds atan(slope) between this cubic
        // approximation and itself plus slope^5/5 (< 2e-31 here).
        let expected = slope - slope.powi(3) / 3.0;
        check(
            planar::<D>(1.0, 0.0),
            planar(1.0, slope),
            expected,
            expected * 4.0 * f64::EPSILON,
        );
    }
    for bits in [(1_u64 << 52) - 1, 1_u64 << 52, (1_u64 << 52) + 1] {
        let slope = f64::from_bits(bits);
        // The cubic correction cannot affect these adjacent floats on either
        // side of the normal/subnormal boundary.
        check(planar::<D>(1.0, 0.0), planar(1.0, slope), slope, 0.0);
    }
}

macro_rules! gen_angle_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<known_angles_ $d d>]() { known_angles::<$d>(); }
            #[test]
            fn [<extreme_angles_ $d d>]() { extreme_angles::<$d>(); }
            #[test]
            fn [<subnormal_angles_ $d d>]() { subnormal_angles::<$d>(); }
            #[test]
            fn [<small_angle_transition_ $d d>]() { small_angle_transition::<$d>(); }
        }
    };
}

gen_angle_tests!(2);
gen_angle_tests!(3);
gen_angle_tests!(4);
gen_angle_tests!(5);
gen_angle_tests!(6);
gen_angle_tests!(8);
gen_angle_tests!(64);

#[test]
fn subnormal_transverse_norm_is_not_squared_away() {
    let tiny = f64::from_bits(1);
    check(
        [1.0, 0.0, 0.0],
        [1.0, 3.0 * tiny, 4.0 * tiny],
        5.0 * tiny,
        0.0,
    );
    // atan(t/sqrt(2)) and sqrt(2)*t both round to t at the least subnormal.
    check([1.0, 1.0, 0.0], [1.0, 1.0, tiny], tiny, 0.0);
    check([1.0, 1.0, tiny], [1.0, 1.0, -tiny], tiny, 0.0);
    check([2.0, tiny, 0.0], [2.0, 0.0, tiny], tiny, 0.0);
}

#[test]
fn one_dimensional_angles_and_signed_zeros() {
    check([f64::MAX], [f64::from_bits(1)], 0.0, 0.0);
    check([-f64::MAX], [f64::from_bits(1)], PI, 0.0);
    check([1.0, 0.0, -0.0], [1.0, -0.0, 0.0], 0.0, 0.0);
    check([1.0, 0.0, -0.0], [-1.0, -0.0, 0.0], PI, 0.0);
}

#[test]
fn nearly_proportional_directions_keep_product_cancellation_residuals() {
    for exponent in 10..=51 {
        let n = 2.0_f64.powi(exponent);
        // det([n, n-1], [n+1, n]) = 1 and dot = 2*n² exactly.
        let expected = (0.5 / (n * n)).atan();
        check(
            [n, n - 1.0],
            [n + 1.0, n],
            expected,
            expected * 4.0 * f64::EPSILON,
        );
        for scale in [2.0_f64.powi(-1000), 2.0_f64.powi(900)] {
            check(
                [n * scale, (n - 1.0) * scale],
                [n + 1.0, n],
                expected,
                expected * 4.0 * f64::EPSILON,
            );
        }
    }
}

#[test]
fn typed_shape_and_zero_errors_have_deterministic_precedence() {
    for (left, right) in [(&[][..], &[1.0][..]), (&[f64::NAN][..], &[1.0, 2.0][..])] {
        assert_matches!(angle_between(left, right), Err(LaError::DimensionMismatch { left: n, right: m, .. }) if n == left.len() && m == right.len());
    }
    assert_matches!(angle_between(&[], &[]), Err(LaError::EmptyVector));
    assert_matches!(
        Vector::<0>::zero().angle(&Vector::zero()),
        Err(LaError::EmptyVector)
    );
    for (left, right, expected) in [
        ([0.0, -0.0], [0.0, 0.0], VectorOperand::Left),
        ([0.0, -0.0], [1.0, 0.0], VectorOperand::Left),
        ([1.0, 0.0], [-0.0, 0.0], VectorOperand::Right),
    ] {
        assert_matches!(angle_between(&left, &right), Err(LaError::ZeroVector { operand, .. }) if operand == expected);
        assert_matches!(Vector::try_new(left).unwrap().angle(&Vector::try_new(right).unwrap()), Err(LaError::ZeroVector { operand, .. }) if operand == expected);
    }
}

#[test]
fn non_finite_errors_identify_operand_and_coordinate_before_zero_errors() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for index in 0..4 {
            let mut bad = [1.0; 4];
            bad[index] = value;
            for (left, right, expected) in [
                (bad, [f64::NAN; 4], VectorOperand::Left),
                ([0.0; 4], bad, VectorOperand::Right),
            ] {
                assert_matches!(angle_between(&left, &right), Err(LaError::NonFinite {
                    location: NonFiniteLocation::VectorOperandEntry { operand, index: actual, .. },
                    origin: NonFiniteOrigin::Input, ..
                }) if operand == expected && actual == index);
            }
        }
    }
}

#[test]
fn angle_error_displays_preserve_context() {
    assert_eq!(
        angle_between(&[], &[1.0]).unwrap_err().to_string(),
        "vector dimension mismatch: left length 0, right length 1"
    );
    assert_eq!(
        angle_between(&[], &[]).unwrap_err().to_string(),
        "an empty vector has no direction"
    );
    assert_eq!(
        angle_between(&[1.0], &[-0.0]).unwrap_err().to_string(),
        "the right zero vector has no direction"
    );
    assert_eq!(
        angle_between(&[0.0], &[1.0]).unwrap_err().to_string(),
        "the left zero vector has no direction"
    );
    assert_eq!(
        angle_between(&[1.0, 1.0], &[1.0, f64::INFINITY])
            .unwrap_err()
            .to_string(),
        "non-finite input value at right vector entry 1"
    );
    assert_eq!(ArithmeticOperation::VectorAngle.to_string(), "vector angle");
}

#[test]
fn borrowed_and_fixed_angles_do_not_allocate() {
    let left = [1.0; 6];
    let right = [2.0; 6];
    let a = Vector::try_new(left).unwrap();
    let b = Vector::try_new(right).unwrap();
    let counts = allocation_counter::measure(|| {
        assert_eq!(angle_between(black_box(&left), black_box(&right)), Ok(0.0));
        assert_eq!(black_box(&a).angle(black_box(&b)), Ok(0.0));
    });
    assert_eq!(counts.count_total, 0);
}
