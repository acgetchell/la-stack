#![forbid(unsafe_code)]

//! Regression tests for bugs caught in public API behavior.

use la_stack::ERR_COEFF_3;
use la_stack::prelude::*;

#[cfg(feature = "exact")]
use proptest::prelude::*;

#[cfg(feature = "exact")]
#[path = "common/proptest_config.rs"]
mod proptest_config;

/// Pad the independently identified near-overflow pair to a fixed dimension.
fn norm_boundary_vector<const D: usize>(left: f64, right: f64) -> Vector<D> {
    let mut values = [0.0; D];
    values[0] = left;
    values[1] = right;
    Vector::try_new(values).expect("the boundary fixture is finite")
}

/// Preserve finite answers that the scaled recurrence previously overflowed.
fn assert_norm_boundary_regressions<const D: usize>() {
    // Expected bits are independently checked with exact squared midpoints in
    // the exact-feature tests below, rather than another floating norm kernel.
    for (left_bits, right_bits, expected_bits) in [
        (
            0x7f55_4c98_5f06_f693,
            0x7fef_fffe_3a58_0905,
            0x7fef_ffff_ffff_ffff,
        ),
        (
            0x7f61_3404_ea4a_8c14,
            0x7fef_fffb_6032_c601,
            0x7fef_ffff_ffff_fffe,
        ),
        (
            0x7f64_7ae1_47ae_147a,
            0x7fef_fff9_7246_996b,
            0x7fef_ffff_ffff_ffff,
        ),
        (
            0x7f69_652b_d3c3_6112,
            0x7fef_fff5_ec54_3df0,
            0x7fef_ffff_ffff_ffff,
        ),
    ] {
        let left = f64::from_bits(left_bits);
        let right = f64::from_bits(right_bits);
        for (left, right) in [(left, right), (right, left), (-left, right), (right, -left)] {
            let vector = norm_boundary_vector::<D>(left, right);
            #[cfg(feature = "exact")]
            assert_high_range_norm_rounding(&vector);
            assert_eq!(
                vector.norm2().unwrap().to_bits(),
                expected_bits,
                "boundary pair ({left_bits:#x}, {right_bits:#x})"
            );
        }
    }

    let two_997 = f64::from_bits(2020_u64 << 52);
    let two_998 = f64::from_bits(2021_u64 << 52);
    let finite = norm_boundary_vector::<D>(f64::MAX, two_997);
    assert_eq!(finite.norm2(), Ok(f64::MAX));
    #[cfg(feature = "exact")]
    assert_high_range_norm_rounding(&finite);

    // The old rounded sqrt could equal 1, hiding true overflow in this case.
    for (left, right) in [(f64::MAX, two_998), (f64::MAX, f64::MAX)] {
        let vector = norm_boundary_vector::<D>(left, right);
        assert_eq!(
            vector.norm2(),
            Err(LaError::non_finite_computation_scalar(
                ArithmeticOperation::VectorNorm
            )),
        );
        #[cfg(feature = "exact")]
        assert_high_range_norm_rounding(&vector);
    }
}

fn assert_norm_guard_transition<const D: usize>() {
    let dimension_bits = usize::BITS - D.leading_zeros();
    let boundary = f64::from_bits(u64::from(2046 - dimension_bits) << 52);

    // Single-coordinate norms are exact on either side of the conservative
    // dispatch boundary, regardless of coordinate position or sign.
    for magnitude in [boundary.next_down(), boundary, boundary.next_up()] {
        for index in 0..D {
            for sign in [-1.0, 1.0] {
                let mut values = [0.0; D];
                values[index] = sign * magnitude;
                let vector = Vector::try_new(values).unwrap();
                assert_eq!(vector.norm2(), Ok(magnitude));
                #[cfg(feature = "exact")]
                assert_high_range_norm_rounding(&vector);
            }
        }
    }

    // A 3-4-5 triple has its largest coordinate at the boundary but its norm
    // above it. Both the recurrence and the exact path must preserve it.
    let unit = boundary / 4.0;
    for (left, right) in [(3.0, 4.0), (4.0, 3.0), (-3.0, 4.0), (4.0, -3.0)] {
        let vector = norm_boundary_vector::<D>(left * unit, right * unit);
        assert_eq!(vector.norm2(), Ok(5.0 * unit));
        #[cfg(feature = "exact")]
        assert_high_range_norm_rounding(&vector);
    }
}

macro_rules! gen_norm_boundary_regressions {
    ($d:literal) => {
        pastey::paste! {
            #[test]
            fn [<norm2_preserves_finite_range_and_true_overflow_ $d d>]() {
                assert_norm_boundary_regressions::<$d>();
            }

            #[test]
            fn [<norm2_guard_transition_preserves_known_values_ $d d>]() {
                assert_norm_guard_transition::<$d>();
            }

            #[cfg(feature = "exact")]
            proptest! {
                #![proptest_config(proptest_config::with_default_cases(64))]

                #[test]
                fn [<norm2_high_range_matches_exact_squared_midpoints_ $d d>](
                    bits in any::<[u64; $d]>(),
                ) {
                    // Mix upper-range and arbitrary finite coordinates. The
                    // first coordinate always selects the exact fallback;
                    // other coordinates exercise widely separated squares,
                    // dense carry propagation, signs, and true overflow.
                    let values = std::array::from_fn(|index| {
                        let raw = bits[index];
                        let exponent = if index == 0 {
                            2046
                        } else if raw & 1 == 0 {
                            2043 + ((raw >> 52) % 4)
                        } else {
                            (raw >> 52) % 2047
                        };
                        f64::from_bits((raw & 0x800f_ffff_ffff_ffff) | (exponent << 52))
                    });
                    let vector = Vector::<$d>::try_new(values).unwrap();
                    assert_high_range_norm_rounding(&vector);
                }
            }
        }
    };
}

gen_norm_boundary_regressions!(2);
gen_norm_boundary_regressions!(3);
gen_norm_boundary_regressions!(4);
gen_norm_boundary_regressions!(5);
gen_norm_boundary_regressions!(6);
gen_norm_boundary_regressions!(7);
gen_norm_boundary_regressions!(8);

/// Validate a high-range result by squaring exact rational rounding midpoints.
#[cfg(feature = "exact")]
fn assert_high_range_norm_rounding<const D: usize>(vector: &Vector<D>) {
    let exact = |value| BigRational::from_f64(value).expect("finite oracle input");
    let square: BigRational = vector
        .as_array()
        .iter()
        .map(|&value| exact(value).pow(2))
        .sum();
    let overflow_midpoint = exact(f64::MAX) + BigRational::from_integer(BigInt::from(1) << 970);
    match vector.norm2() {
        Ok(norm) => {
            assert!(square < overflow_midpoint.pow(2));
            let lower_midpoint = (exact(norm.next_down()) + exact(norm)) / exact(2.0);
            let upper_midpoint = if norm.to_bits() == f64::MAX.to_bits() {
                overflow_midpoint
            } else {
                (exact(norm) + exact(norm.next_up())) / exact(2.0)
            };
            let lower_square = lower_midpoint.pow(2);
            let upper_square = upper_midpoint.pow(2);
            assert!(lower_square <= square && square <= upper_square);
            if square == lower_square || square == upper_square {
                assert_eq!(
                    norm.to_bits() & 1,
                    0,
                    "ties must select the even significand"
                );
            }
        }
        Err(error) => {
            assert_eq!(
                error,
                LaError::non_finite_computation_scalar(ArithmeticOperation::VectorNorm)
            );
            assert!(
                square >= overflow_midpoint.pow(2),
                "finite norm misclassified as overflow"
            );
        }
    }
}

#[test]
fn norm2_upper_range_ties_and_subnormal_tail() {
    let unit = f64::from_bits(1993_u64 << 52); // 2^970
    for offset in [1.0, 3.0] {
        let k = 2.0_f64.powi(51) + offset;
        // A scaled 3-4-5 triple has exact norm 5*k*2^970. Odd 54-bit 5*k
        // lies exactly halfway between two binary64 values; the two offsets
        // exercise both directions of ties-to-even.
        let vector = norm_boundary_vector::<3>((3.0 * k) * unit, (4.0 * k) * unit);
        let expected = (5.0 * k) * unit;
        assert_eq!(vector.norm2().unwrap().to_bits(), expected.to_bits());
        #[cfg(feature = "exact")]
        assert_high_range_norm_rounding(&vector);

        let mut values = vector.into_array();
        values[2] = f64::from_bits(1);
        let with_tail = Vector::try_new(values).unwrap();
        let expected_with_tail = if offset < 2.0 {
            expected.next_up()
        } else {
            expected
        };
        assert_eq!(
            with_tail.norm2().unwrap().to_bits(),
            expected_with_tail.to_bits()
        );
        #[cfg(feature = "exact")]
        assert_high_range_norm_rounding(&with_tail);
    }
}

#[test]
#[cfg(feature = "exact")]
fn det_exact_f64_preserves_min_positive_subnormal() -> Result<(), LaError> {
    let tiny = f64::from_bits(1);
    let m = Matrix::<1>::try_from_rows([[tiny]])?;

    assert_eq!(m.det_exact_f64().unwrap().to_bits(), tiny.to_bits());
    Ok(())
}

#[test]
#[cfg(feature = "exact")]
fn det_exact_f64_strict_vs_rounded_inexact_det() -> Result<(), LaError> {
    let m = Matrix::<2>::try_from_rows([[1.0 + f64::EPSILON, 0.0], [0.0, 1.0 - f64::EPSILON]])?;

    assert!(matches!(
        m.det_exact_f64(),
        Err(LaError::Unrepresentable {
            index: None,
            reason: UnrepresentableReason::RequiresRounding,
            ..
        })
    ));
    assert_eq!(
        m.det_exact_rounded_f64().unwrap().to_bits(),
        1.0f64.to_bits()
    );
    Ok(())
}

#[test]
#[cfg(feature = "exact")]
fn solve_exact_f64_strict_vs_rounded_non_dyadic() -> Result<(), LaError> {
    let a = Matrix::<1>::try_from_rows([[3.0]])?;
    let b = Vector::<1>::try_new([1.0])?;

    assert!(matches!(
        a.solve_exact_f64(b),
        Err(LaError::Unrepresentable {
            index: Some(0),
            reason: UnrepresentableReason::RequiresRounding,
            ..
        })
    ));
    assert_eq!(
        a.solve_exact_rounded_f64(b).unwrap().into_array()[0].to_bits(),
        (1.0f64 / 3.0).to_bits()
    );
    Ok(())
}

#[test]
#[cfg(feature = "exact")]
fn requires_rounding_error_can_fall_back_to_rounded_solve() -> Result<(), LaError> {
    let a = Matrix::<1>::try_from_rows([[3.0]])?;
    let b = Vector::<1>::try_new([1.0])?;

    let rounded = match a.solve_exact_f64(b) {
        Ok(x) => x,
        Err(err) if err.requires_rounding() => a.solve_exact_rounded_f64(b)?,
        Err(err) => return Err(err),
    };

    assert_eq!(rounded.into_array()[0].to_bits(), (1.0f64 / 3.0).to_bits());
    Ok(())
}

#[test]
#[cfg(feature = "exact")]
fn exact_determinant_overflow_midpoint_is_not_recoverable_by_rounding() -> Result<(), LaError> {
    let above_overflow_midpoint = 3.0 * 2.0_f64.powi(969);
    let m = Matrix::<2>::try_from_rows([[f64::MAX, -above_overflow_midpoint], [1.0, 1.0]])?;

    let strict = m.det_exact_f64().unwrap_err();
    assert!(matches!(
        strict,
        LaError::Unrepresentable {
            index: None,
            reason: UnrepresentableReason::NotFinite,
            ..
        }
    ));
    assert!(!strict.requires_rounding());
    assert!(matches!(
        m.det_exact_rounded_f64(),
        Err(LaError::Unrepresentable {
            reason: UnrepresentableReason::NotFinite,
            ..
        })
    ));

    Ok(())
}

#[test]
fn det_direct_skips_zero_coefficient_terms_that_would_overflow() -> Result<(), LaError> {
    let d3 = Matrix::<3>::try_from_rows([
        [1.0, 0.0, 0.0],
        [1.0e300, 1.0, 1.0e300],
        [1.0e300, 0.0, 1.0e300],
    ])?;
    assert_eq!(d3.det_direct(), Ok(Some(1.0e300)));

    let d4 = Matrix::<4>::try_from_rows([
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [1.0e300, 0.0, 1.0e300, 1.0e300],
        [1.0e300, 0.0, 1.0e300, -1.0e300],
    ])?;
    assert_eq!(d4.det_direct(), Ok(Some(0.0)));

    Ok(())
}

#[test]
fn det_errbound_skips_zero_coefficient_terms_that_would_overflow() -> Result<(), LaError> {
    let d3 = Matrix::<3>::try_from_rows([
        [1.0, 0.0, 0.0],
        [1.0e300, 1.0, 1.0e300],
        [1.0e300, 0.0, 1.0e300],
    ])?;
    assert_eq!(d3.det_errbound(), Ok(Some(ERR_COEFF_3 * 1.0e300)));

    let d4 = Matrix::<4>::try_from_rows([
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [1.0e300, 0.0, 1.0e300, 1.0e300],
        [1.0e300, 0.0, 1.0e300, -1.0e300],
    ])?;
    assert_eq!(d4.det_errbound(), Ok(Some(0.0)));

    Ok(())
}

#[test]
fn det_errbound_reports_overflow_even_when_another_term_underflows() -> Result<(), LaError> {
    let matrix = Matrix::<2>::try_from_rows([[1.0e308, 1.0e-308], [1.0e-308, 2.0]])?;
    let expected = LaError::non_finite_computation_scalar(ArithmeticOperation::Determinant);

    assert_eq!(matrix.det_direct(), Err(expected));
    assert_eq!(matrix.det_direct_with_errbound(), Err(expected));
    assert_eq!(matrix.det_errbound(), Err(expected));

    Ok(())
}
