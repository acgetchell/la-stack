#![forbid(unsafe_code)]

//! Property tests for outward-rounded interval arithmetic and determinant signs.
//!
//! The oracle converts binary64 inputs independently to `BigRational` and uses
//! rational Gaussian elimination. Production interval determinants instead use
//! division-free subset dynamic programming.

#![cfg(feature = "exact")]

use std::array::from_fn;

use pastey::paste;
use proptest::prelude::*;

use la_stack::prelude::*;

#[path = "common/proptest_config.rs"]
mod proptest_config;
use proptest_config::with_default_cases;

fn exact_f64(value: f64) -> BigRational {
    BigRational::from_f64(value).expect("the generated binary64 value is finite")
}

fn interval_contains_exact(interval: Interval, value: &BigRational) -> bool {
    exact_f64(interval.lower()) <= *value && *value <= exact_f64(interval.upper())
}

fn finite_f64() -> impl Strategy<Value = f64> {
    any::<u64>()
        .prop_map(f64::from_bits)
        .prop_filter("binary64 value must be finite", |value| value.is_finite())
}

fn exact_fits_finite_interval(value: &BigRational) -> bool {
    let maximum = exact_f64(f64::MAX);
    -&maximum <= *value && *value <= maximum
}

fn assert_outward_result(
    result: Result<Interval, LaError>,
    exact: &BigRational,
    expected_operation: ArithmeticOperation,
) -> Result<(), TestCaseError> {
    match result {
        Ok(interval) => {
            prop_assert!(interval_contains_exact(interval, exact));
            prop_assert!(exact_fits_finite_interval(exact));
        }
        Err(LaError::IntervalRangeExhausted { operation, .. }) => {
            prop_assert_eq!(operation, expected_operation);
            prop_assert!(!exact_fits_finite_interval(exact));
        }
        Err(error) => {
            return Err(TestCaseError::fail(format!(
                "unexpected interval error: {error}"
            )));
        }
    }
    Ok(())
}

fn assert_conclusive_sign_matches(
    sign: IntervalDeterminantSign,
    exact: &BigRational,
) -> Result<(), TestCaseError> {
    match sign {
        IntervalDeterminantSign::Positive => prop_assert!(exact.is_positive()),
        IntervalDeterminantSign::Negative => prop_assert!(exact.is_negative()),
        IntervalDeterminantSign::Zero => {
            prop_assert_eq!(exact, &BigRational::from_integer(0.into()));
        }
        IntervalDeterminantSign::Inconclusive => {}
        _ => prop_assert!(false, "unknown interval determinant sign"),
    }
    Ok(())
}

/// Independent determinant oracle using rational Gaussian elimination.
fn rational_det<const D: usize>(rows: &[[f64; D]; D]) -> BigRational {
    let mut work: [[BigRational; D]; D] =
        from_fn(|row| from_fn(|column| exact_f64(rows[row][column])));
    let mut negative = false;

    for column in 0..D {
        let mut pivot_row = column;
        while pivot_row < D && work[pivot_row][column] == BigRational::from_integer(0.into()) {
            pivot_row += 1;
        }
        if pivot_row == D {
            return BigRational::from_integer(0.into());
        }
        if pivot_row != column {
            work.swap(pivot_row, column);
            negative = !negative;
        }

        let pivot = work[column][column].clone();
        let pivot_row = work[column].clone();
        for row in work.iter_mut().skip(column + 1) {
            let factor = &row[column] / &pivot;
            for (entry, pivot_entry) in row.iter_mut().zip(pivot_row.iter()).skip(column) {
                let reduction = &factor * pivot_entry;
                *entry -= reduction;
            }
        }
    }

    let mut determinant = BigRational::from_integer(1.into());
    for (index, row) in work.iter().enumerate() {
        determinant *= &row[index];
    }
    if negative { -determinant } else { determinant }
}

fn f64_rows<const D: usize>(rows: [[i8; D]; D]) -> [[f64; D]; D] {
    rows.map(|row| row.map(f64::from))
}

proptest! {
    #![proptest_config(with_default_cases(128))]

    #[test]
    fn scalar_interval_operations_contain_exact_rational_results(
        left in -10_000i16..=10_000,
        right in -10_000i16..=10_000,
    ) {
        let left = f64::from(left) / 10.0;
        let right = f64::from(right) / 10.0;
        let exact_left = exact_f64(left);
        let exact_right = exact_f64(right);
        let left_interval = Interval::point(left)?;
        let right_interval = Interval::point(right)?;

        let subtraction = Interval::try_from_subtraction(left, right)?;
        prop_assert!(interval_contains_exact(
            subtraction,
            &(&exact_left - &exact_right),
        ));

        let addition = left_interval.try_add(&right_interval)?;
        prop_assert!(interval_contains_exact(
            addition,
            &(&exact_left + &exact_right),
        ));

        let product = left_interval.try_mul(&right_interval)?;
        prop_assert!(interval_contains_exact(
            product,
            &(&exact_left * &exact_right),
        ));

        let square = left_interval.try_square()?;
        prop_assert!(interval_contains_exact(
            square,
            &(&exact_left * &exact_left),
        ));
    }

    #[test]
    fn wide_interval_operations_enclose_all_endpoint_extrema(
        left_center in -100i16..=100,
        left_radius in 0u8..=10,
        right_center in -100i16..=100,
        right_radius in 0u8..=10,
    ) {
        let left_lower = f64::from(left_center - i16::from(left_radius));
        let left_upper = f64::from(left_center + i16::from(left_radius));
        let right_lower = f64::from(right_center - i16::from(right_radius));
        let right_upper = f64::from(right_center + i16::from(right_radius));
        let left = Interval::try_new(left_lower, left_upper)?;
        let right = Interval::try_new(right_lower, right_upper)?;

        let sum = left.try_add(&right)?;
        prop_assert!(interval_contains_exact(
            sum,
            &(exact_f64(left_lower) + exact_f64(right_lower)),
        ));
        prop_assert!(interval_contains_exact(
            sum,
            &(exact_f64(left_upper) + exact_f64(right_upper)),
        ));

        let product = left.try_mul(&right)?;
        for (left_endpoint, right_endpoint) in [
            (left_lower, right_lower),
            (left_lower, right_upper),
            (left_upper, right_lower),
            (left_upper, right_upper),
        ] {
            prop_assert!(interval_contains_exact(
                product,
                &(exact_f64(left_endpoint) * exact_f64(right_endpoint)),
            ));
        }

        let square = left.try_square()?;
        prop_assert!(interval_contains_exact(
            square,
            &(exact_f64(left_lower) * exact_f64(left_lower)),
        ));
        prop_assert!(interval_contains_exact(
            square,
            &(exact_f64(left_upper) * exact_f64(left_upper)),
        ));
        if left.contains(0.0) {
            prop_assert!(square.contains(0.0));
        }
    }
}

proptest! {
    #![proptest_config(with_default_cases(512))]

    #[test]
    fn arbitrary_finite_point_operations_are_outward_or_report_true_range_loss(
        left in finite_f64(),
        right in finite_f64(),
    ) {
        let exact_left = exact_f64(left);
        let exact_right = exact_f64(right);
        let left_interval = Interval::point(left)?;
        let right_interval = Interval::point(right)?;

        assert_outward_result(
            Interval::try_from_subtraction(left, right),
            &(&exact_left - &exact_right),
            ArithmeticOperation::IntervalSubtraction,
        )?;
        assert_outward_result(
            left_interval.try_add(&right_interval),
            &(&exact_left + &exact_right),
            ArithmeticOperation::IntervalAddition,
        )?;
        assert_outward_result(
            left_interval.try_mul(&right_interval),
            &(&exact_left * &exact_right),
            ArithmeticOperation::IntervalMultiplication,
        )?;
        assert_outward_result(
            left_interval.try_square(),
            &(&exact_left * &exact_left),
            ArithmeticOperation::IntervalSquare,
        )?;
    }
}

macro_rules! gen_interval_determinant_proptests {
    ($d:literal) => {
        paste! {
            proptest! {
                #![proptest_config(with_default_cases(24))]

                #[test]
                fn [<point_interval_determinant_matches_rational_oracle_ $d d>](
                    raw_rows in any::<[[i8; $d]; $d]>(),
                ) {
                    let rows = f64_rows(raw_rows);
                    let matrix = IntervalMatrix::<$d>::try_from_point_rows(rows)?;
                    let determinant = matrix.det()?;
                    let exact = rational_det(&rows);
                    prop_assert!(interval_contains_exact(determinant, &exact));

                    match matrix.det_sign()? {
                        IntervalDeterminantSign::Positive => prop_assert!(exact.is_positive()),
                        IntervalDeterminantSign::Negative => prop_assert!(exact.is_negative()),
                        IntervalDeterminantSign::Zero => prop_assert_eq!(
                            exact,
                            BigRational::from_integer(0.into()),
                        ),
                        IntervalDeterminantSign::Inconclusive => {
                            prop_assert!(determinant.contains(0.0));
                        }
                        _ => prop_assert!(false, "unknown interval determinant sign"),
                    }
                }

                #[test]
                fn [<wide_interval_determinant_contains_selected_oracles_ $d d>](
                    raw_rows in any::<[[i8; $d]; $d]>(),
                    raw_radii in any::<[[u8; $d]; $d]>(),
                ) {
                    let center_rows = f64_rows(raw_rows);
                    let lower_rows = from_fn(|row| from_fn(|column| {
                        let center = i16::from(raw_rows[row][column]);
                        let radius = i16::from(raw_radii[row][column] % 3);
                        f64::from(center - radius)
                    }));
                    let upper_rows = from_fn(|row| from_fn(|column| {
                        let center = i16::from(raw_rows[row][column]);
                        let radius = i16::from(raw_radii[row][column] % 3);
                        f64::from(center + radius)
                    }));
                    let interval_rows = from_fn(|row| from_fn(|column| {
                        Interval::try_new(lower_rows[row][column], upper_rows[row][column])
                            .expect("integer endpoints are finite and ordered")
                    }));
                    let matrix = IntervalMatrix::<$d>::from_rows(interval_rows);
                    let determinant = matrix.det()?;
                    let sign = matrix.det_sign()?;

                    for selected_rows in [&lower_rows, &center_rows, &upper_rows] {
                        let exact = rational_det(selected_rows);
                        prop_assert!(interval_contains_exact(determinant, &exact));
                        assert_conclusive_sign_matches(sign, &exact)?;
                    }
                }
            }
        }
    };
}

gen_interval_determinant_proptests!(2);
gen_interval_determinant_proptests!(3);
gen_interval_determinant_proptests!(4);
gen_interval_determinant_proptests!(5);
gen_interval_determinant_proptests!(6);
gen_interval_determinant_proptests!(7);

proptest! {
    #![proptest_config(with_default_cases(48))]

    #[test]
    fn wide_2d_determinant_contains_every_endpoint_matrix(
        raw_rows in any::<[[i8; 2]; 2]>(),
        raw_radii in any::<[[u8; 2]; 2]>(),
    ) {
        let lower_rows: [[f64; 2]; 2] = from_fn(|row| from_fn(|column| {
            let center = i16::from(raw_rows[row][column]);
            let radius = i16::from(raw_radii[row][column] % 3);
            f64::from(center - radius)
        }));
        let upper_rows: [[f64; 2]; 2] = from_fn(|row| from_fn(|column| {
            let center = i16::from(raw_rows[row][column]);
            let radius = i16::from(raw_radii[row][column] % 3);
            f64::from(center + radius)
        }));
        let matrix = IntervalMatrix::<2>::from_rows(from_fn(|row| from_fn(|column| {
            Interval::try_new(lower_rows[row][column], upper_rows[row][column])
                .expect("integer endpoints are finite and ordered")
        })));
        let determinant = matrix.det()?;
        let sign = matrix.det_sign()?;

        for endpoint_mask in 0_u8..16 {
            let selected_rows: [[f64; 2]; 2] = from_fn(|row| from_fn(|column| {
                let bit = 1_u8 << (row * 2 + column);
                if endpoint_mask & bit == 0 {
                    lower_rows[row][column]
                } else {
                    upper_rows[row][column]
                }
            }));
            let exact = rational_det(&selected_rows);
            prop_assert!(interval_contains_exact(determinant, &exact));
            assert_conclusive_sign_matches(sign, &exact)?;
        }
    }
}
