//! Tests for exact benchmark input-generation configuration.

#![cfg(all(feature = "bench", feature = "exact"))]
#![forbid(unsafe_code)]

#[path = "../benches/common/exact.rs"]
pub mod exact_bench;

use core::array::from_fn;
use std::error::Error;

use exact_bench::{
    ExactBenchConfigError, ExactInput, I16Range, SplitMix64, ValidatedExactInput, hilbert_input,
    large_entries_3x3_input, make_matrix_rows, make_random_input_corpus, make_vector_array,
    near_singular_3x3_input, validate_exact_fixture,
};
use la_stack::{Matrix, Vector};
use pastey::paste;

fn baseline_input<const D: usize>() -> ExactInput<D> {
    let Ok(matrix) = Matrix::<D>::try_from_rows(make_matrix_rows::<D>()) else {
        panic!("baseline benchmark matrix should be finite");
    };
    let Ok(rhs) = Vector::<D>::try_new(make_vector_array::<D>()) else {
        panic!("baseline benchmark RHS should be finite");
    };
    ExactInput { matrix, rhs }
}

fn validate_baseline_and_random_corpus<const D: usize>() {
    let _ = validate_exact_fixture(baseline_input::<D>());
    for input in make_random_input_corpus::<D>() {
        let _ = validate_exact_fixture(input);
    }
}

#[test]
fn i16_range_rejects_unordered_bounds() {
    let Err(err) = I16Range::try_new(5, 4) else {
        panic!("unordered range should be rejected");
    };
    assert_eq!(
        err,
        ExactBenchConfigError::UnorderedRange { min: 5, max: 4 }
    );
}

#[test]
fn exact_bench_config_error_is_std_error() {
    let err = ExactBenchConfigError::UnorderedRange { min: 5, max: 4 };
    let as_error: &dyn Error = &err;

    assert_eq!(
        as_error.to_string(),
        "random integer range must be ordered: 5..=4"
    );
    assert!(as_error.source().is_none());
}

#[test]
fn i16_range_single_value_always_draws_that_value() {
    let Ok(range) = I16Range::try_new(7, 7) else {
        panic!("single-value range should be valid");
    };
    let mut rng = SplitMix64::new(0);

    for _ in 0..32 {
        assert_eq!(rng.next_i16(range), 7);
    }
}

#[test]
fn i16_range_draws_stay_inside_inclusive_bounds() {
    let Ok(range) = I16Range::try_new(-10, 10) else {
        panic!("ordered range should be valid");
    };
    let mut rng = SplitMix64::new(0xCAFE_F00D);

    for _ in 0..1024 {
        assert!((-10..=10).contains(&rng.next_i16(range)));
    }
}

#[test]
fn i16_range_full_domain_reaches_both_endpoints() {
    let Ok(range) = I16Range::try_new(i16::MIN, i16::MAX) else {
        panic!("the full i16 domain should be a valid range");
    };

    // These seeds make the first modulo offsets 0 and 65,535, respectively.
    let mut minimum_rng = SplitMix64::new(59_587);
    let mut maximum_rng = SplitMix64::new(16_165);

    assert_eq!(minimum_rng.next_i16(range), i16::MIN);
    assert_eq!(maximum_rng.next_i16(range), i16::MAX);
}

#[test]
fn splitmix64_sequence_is_stable_for_benchmark_seed() {
    let Ok(range) = I16Range::try_new(-10, 10) else {
        panic!("ordered range should be valid");
    };
    let mut rng = SplitMix64::new(0xCAFE_F00D);

    let draws = from_fn::<_, 16, _>(|_| rng.next_i16(range));

    assert_eq!(
        draws,
        [-8, 4, -3, -8, 8, -4, 4, -6, -9, -8, 6, 1, 0, 3, -8, -1]
    );
}

macro_rules! gen_exact_benchmark_fixture_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<exact_benchmark_fixtures_are_correct_ $d d>]() {
                validate_baseline_and_random_corpus::<$d>();
            }
        }
    };
}

gen_exact_benchmark_fixture_tests!(2);
gen_exact_benchmark_fixture_tests!(3);
gen_exact_benchmark_fixture_tests!(4);
gen_exact_benchmark_fixture_tests!(5);

#[test]
fn exact_adversarial_benchmark_fixtures_are_correct() {
    for input in [near_singular_3x3_input(), large_entries_3x3_input()] {
        let _ = validate_exact_fixture(input);
    }
    let _ = validate_exact_fixture(hilbert_input::<4>());
    let _ = validate_exact_fixture(hilbert_input::<5>());
}

#[test]
fn validated_fixture_exposes_only_checked_inputs() {
    let raw = baseline_input::<3>();
    let expected_matrix = raw.matrix;
    let expected_rhs = raw.rhs;

    let validated: ValidatedExactInput<3> = validate_exact_fixture(raw);

    assert_eq!(validated.matrix(), &expected_matrix);
    assert_eq!(validated.rhs(), expected_rhs);
}

#[test]
#[should_panic(expected = "exact solve oracle check failed")]
fn fixture_validation_rejects_singular_solve_input() {
    let matrix = Matrix::<2>::try_from_rows([[1.0, 2.0], [2.0, 4.0]])
        .unwrap_or_else(|error| panic!("singular fixture must still be finite: {error}"));
    let rhs = Vector::<2>::try_new([1.0, 2.0])
        .unwrap_or_else(|error| panic!("singular fixture RHS must be finite: {error}"));

    let _ = validate_exact_fixture(ExactInput { matrix, rhs });
}
