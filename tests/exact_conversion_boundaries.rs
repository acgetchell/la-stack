//! Independent bit-pattern coverage for exact-to-binary64 conversion boundaries.

#![forbid(unsafe_code)]
#![cfg(feature = "exact")]

use core::array::from_fn;
use core::cmp::Ordering;

use pastey::paste;

use la_stack::prelude::*;

const POSITIVE_ZERO_BITS: u64 = 0;
const NEGATIVE_ZERO_BITS: u64 = 1_u64 << 63;
const BELOW_OVERFLOW_MIDPOINT_INCREMENT: f64 = f64::from_bits(1992_u64 << 52); // 2^969
const AT_OVERFLOW_MIDPOINT_INCREMENT: f64 = f64::from_bits(1993_u64 << 52); // 2^970

fn assert_unrepresentable<T>(
    result: &Result<T, LaError>,
    index: Option<usize>,
    reason: UnrepresentableReason,
) {
    assert!(matches!(
        result,
        Err(LaError::Unrepresentable {
            index: actual_index,
            reason: actual_reason,
            ..
        }) if *actual_index == index && *actual_reason == reason
    ));
}

fn diagonal<const D: usize>(values: [f64; D]) -> Matrix<D> {
    let mut rows = [[0.0; D]; D];
    for (index, value) in values.into_iter().enumerate() {
        rows[index][index] = value;
    }
    Matrix::try_from_rows(rows).expect("diagonal fixture is finite")
}

fn determinant_near_overflow(increment: f64, negative: bool) -> Matrix<2> {
    let rows = if negative {
        [[-f64::MAX, increment], [1.0, 1.0]]
    } else {
        [[f64::MAX, -increment], [1.0, 1.0]]
    };
    Matrix::try_from_rows(rows).expect("overflow-boundary fixture is finite")
}

fn raw_rational(numerator: i32, denominator: i32) -> BigRational {
    BigRational::new_raw(BigInt::from(numerator), BigInt::from(denominator))
}

fn canonical_vector_conversions_preserve_raw_contract<const D: usize>() {
    let success = [
        (raw_rational(0, -7), 0.0_f64.to_bits()),
        (raw_rational(3, -6), (-0.5_f64).to_bits()),
        (
            BigRational::new(1.into(), BigInt::from(1_u8) << 1074_u32),
            1,
        ),
        (
            BigRational::from_float(f64::MAX).unwrap(),
            f64::MAX.to_bits(),
        ),
    ];
    for (value, expected) in success {
        let raw = from_fn::<_, D, _>(|_| value.clone());
        let canonical = RationalVector::try_new(raw.clone()).unwrap();
        for actual in [
            canonical.try_to_f64(),
            canonical.to_rounded_f64(),
            raw.try_to_f64(),
            raw.to_rounded_f64(),
        ] {
            assert_eq!(
                actual.unwrap().into_array().map(f64::to_bits),
                [expected; D]
            );
        }
    }

    let failures = [
        (
            raw_rational(2, 6),
            UnrepresentableReason::RequiresRounding,
            Ok(0x3fd5_5555_5555_5555), // nearest binary64 to 1/3
        ),
        (
            BigRational::new((-1).into(), BigInt::from(1_u8) << 1075_u32),
            UnrepresentableReason::RequiresRounding,
            Ok(NEGATIVE_ZERO_BITS),
        ),
        (
            BigRational::from_integer(BigInt::from(1_u8) << 1024_u32),
            UnrepresentableReason::NotFinite,
            Err(UnrepresentableReason::NotFinite),
        ),
    ];
    for (value, reason, rounded) in failures {
        for index in 0..D {
            let raw = from_fn(|i| match i.cmp(&index) {
                Ordering::Equal => value.clone(),
                Ordering::Greater => BigRational::from_integer(BigInt::from(1_u8) << 1024_u32),
                Ordering::Less => raw_rational(3, -6),
            });
            let canonical = RationalVector::<D>::try_new(raw.clone()).unwrap();
            assert_unrepresentable(&canonical.try_to_f64(), Some(index), reason);
            assert_unrepresentable(&raw.try_to_f64(), Some(index), reason);
            // Rounding may recover this component, but must then report the
            // first later overflow rather than retaining the strict error.
            let expected_rounded = match rounded {
                Err(reason) => Err(LaError::unrepresentable(Some(index), reason)),
                Ok(_) if index + 1 < D => Err(LaError::unrepresentable(
                    Some(index + 1),
                    UnrepresentableReason::NotFinite,
                )),
                Ok(bits) => {
                    let mut expected = [(-0.5_f64).to_bits(); D];
                    expected[index] = bits;
                    Ok(expected)
                }
            };
            for (path, actual) in [
                ("canonical", canonical.to_rounded_f64()),
                ("raw", raw.to_rounded_f64()),
            ] {
                assert_eq!(
                    actual.map(|v| v.into_array().map(f64::to_bits)),
                    expected_rounded,
                    "{path}, D={D}, index={index}, value={value}",
                );
            }
        }
    }
}

macro_rules! gen_canonical_conversion_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<canonical_vector_conversion_contract_ $d d>]() {
                canonical_vector_conversions_preserve_raw_contract::<$d>();
            }
        }
    };
}

gen_canonical_conversion_tests!(2);
gen_canonical_conversion_tests!(3);
gen_canonical_conversion_tests!(4);
gen_canonical_conversion_tests!(5);

#[test]
fn raw_rational_conversion_uses_the_mathematical_quotient() {
    let cases = [
        (raw_rational(1, -2), (-0.5_f64).to_bits()),
        (raw_rational(-1, -2), 0.5_f64.to_bits()),
        (raw_rational(3, 6), 0.5_f64.to_bits()),
    ];

    for (exact, expected_bits) in &cases {
        assert_eq!(exact.try_to_f64().unwrap().to_bits(), *expected_bits);
        assert_eq!(exact.to_rounded_f64().unwrap().to_bits(), *expected_bits);
    }

    let exact = cases.map(|(value, _)| value);
    let strict = exact.try_to_f64().unwrap().into_array().map(f64::to_bits);
    let rounded = exact
        .to_rounded_f64()
        .unwrap()
        .into_array()
        .map(f64::to_bits);
    let expected = [(-0.5_f64).to_bits(), 0.5_f64.to_bits(), 0.5_f64.to_bits()];
    assert_eq!(strict, expected);
    assert_eq!(rounded, expected);
}

#[test]
fn raw_zero_denominator_is_not_finite_and_arrays_report_its_first_index() {
    for exact in [raw_rational(0, 0), raw_rational(1, 0)] {
        assert_unrepresentable(&exact.try_to_f64(), None, UnrepresentableReason::NotFinite);
        assert_unrepresentable(
            &exact.to_rounded_f64(),
            None,
            UnrepresentableReason::NotFinite,
        );
    }

    let exact = [
        raw_rational(1, -2),
        raw_rational(-1, -2),
        raw_rational(3, 6),
        raw_rational(0, 0),
        raw_rational(1, 0),
    ];
    assert_unrepresentable(
        &exact.try_to_f64(),
        Some(3),
        UnrepresentableReason::NotFinite,
    );
    assert_unrepresentable(
        &exact.to_rounded_f64(),
        Some(3),
        UnrepresentableReason::NotFinite,
    );
}

#[test]
fn d0_exact_strict_and_rounded_outputs_follow_empty_product_conventions() {
    let matrix = Matrix::<0>::zero();
    let rhs = Vector::<0>::zero();

    let determinant = matrix.det_exact().unwrap();
    assert_eq!(determinant.try_to_f64(), Ok(1.0));
    assert_eq!(determinant.to_rounded_f64(), Ok(1.0));
    assert_eq!(matrix.det_exact_f64(), Ok(1.0));
    assert_eq!(matrix.det_exact_rounded_f64(), Ok(1.0));

    let solution = matrix.solve_exact(rhs).unwrap();
    assert!(solution.try_to_f64().unwrap().as_array().is_empty());
    assert!(solution.to_rounded_f64().unwrap().as_array().is_empty());
    assert!(matrix.solve_exact_f64(rhs).unwrap().as_array().is_empty());
    assert!(
        matrix
            .solve_exact_rounded_f64(rhs)
            .unwrap()
            .as_array()
            .is_empty()
    );
}

#[test]
fn overflow_midpoint_classification_is_symmetric_and_bit_exact() {
    for negative in [false, true] {
        let below = determinant_near_overflow(BELOW_OVERFLOW_MIDPOINT_INCREMENT, negative);
        let exact_below = below.det_exact().unwrap();
        assert_unrepresentable(
            &exact_below.try_to_f64(),
            None,
            UnrepresentableReason::RequiresRounding,
        );
        let expected = if negative { -f64::MAX } else { f64::MAX };
        assert_eq!(
            exact_below.to_rounded_f64().unwrap().to_bits(),
            expected.to_bits()
        );
        assert_eq!(
            below.det_exact_rounded_f64().unwrap().to_bits(),
            expected.to_bits()
        );

        let midpoint = determinant_near_overflow(AT_OVERFLOW_MIDPOINT_INCREMENT, negative);
        let exact_midpoint = midpoint.det_exact().unwrap();
        assert_unrepresentable(
            &exact_midpoint.try_to_f64(),
            None,
            UnrepresentableReason::NotFinite,
        );
        assert_unrepresentable(
            &exact_midpoint.to_rounded_f64(),
            None,
            UnrepresentableReason::NotFinite,
        );
        assert_unrepresentable(
            &midpoint.det_exact_rounded_f64(),
            None,
            UnrepresentableReason::NotFinite,
        );
    }
}

#[test]
fn subnormal_halfway_cases_round_to_even_with_signed_zero() {
    let tiny = f64::from_bits(1);
    let cases = [
        ([tiny, 0.5], POSITIVE_ZERO_BITS),
        ([-tiny, 0.5], NEGATIVE_ZERO_BITS),
        ([f64::from_bits(3), 0.5], 2),
        ([-f64::from_bits(3), 0.5], NEGATIVE_ZERO_BITS | 2),
    ];

    for (diagonal_values, expected_bits) in cases {
        let matrix = diagonal(diagonal_values);
        let exact = matrix.det_exact().unwrap();
        assert_unrepresentable(
            &exact.try_to_f64(),
            None,
            UnrepresentableReason::RequiresRounding,
        );
        assert_eq!(exact.to_rounded_f64().unwrap().to_bits(), expected_bits);
        assert_eq!(
            matrix.det_exact_rounded_f64().unwrap().to_bits(),
            expected_bits
        );
    }
}

#[test]
fn solution_conversion_preserves_first_index_and_negative_underflow_zero() {
    let tiny = f64::from_bits(1);
    let matrix = diagonal([1.0, 2.0]);
    let rhs = Vector::<2>::try_new([0.0, -tiny]).unwrap();
    let exact = matrix.solve_exact(rhs).unwrap();

    assert_unrepresentable(
        &exact.try_to_f64(),
        Some(1),
        UnrepresentableReason::RequiresRounding,
    );
    let rounded = exact.to_rounded_f64().unwrap().into_array();
    assert_eq!(rounded[0].to_bits(), POSITIVE_ZERO_BITS);
    assert_eq!(rounded[1].to_bits(), NEGATIVE_ZERO_BITS);

    assert_unrepresentable(
        &matrix.solve_exact_f64(rhs),
        Some(1),
        UnrepresentableReason::RequiresRounding,
    );
    assert_eq!(
        matrix.solve_exact_rounded_f64(rhs).unwrap().as_array()[1].to_bits(),
        NEGATIVE_ZERO_BITS
    );
}

#[test]
fn d5_exact_sign_and_conversions_handle_final_underflow() {
    let tiny = f64::from_bits(1);
    let matrix = diagonal([tiny, tiny, 1.0, 1.0, 1.0]);

    assert_eq!(matrix.det_sign_exact(), DeterminantSign::Positive);
    let exact = matrix.det_exact().unwrap();
    assert_unrepresentable(
        &exact.try_to_f64(),
        None,
        UnrepresentableReason::RequiresRounding,
    );
    assert_eq!(
        exact.to_rounded_f64().unwrap().to_bits(),
        POSITIVE_ZERO_BITS
    );
}
