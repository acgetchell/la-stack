//! Tests for exact benchmark input-generation configuration.

#![cfg(all(feature = "bench", feature = "exact"))]
#![forbid(unsafe_code)]

#[path = "../benches/common/exact.rs"]
mod exact_bench;

use core::array::from_fn;

use exact_bench::{ExactBenchConfigError, I16Range, SplitMix64};

#[test]
fn i16_range_rejects_unordered_bounds() {
    let Err(err) = I16Range::new(5, 4) else {
        panic!("unordered range should be rejected");
    };
    assert_eq!(
        err,
        ExactBenchConfigError::UnorderedRange { min: 5, max: 4 }
    );
}

#[test]
fn empty_corpus_error_message_names_requirement() {
    assert_eq!(
        ExactBenchConfigError::EmptyCorpus.to_string(),
        "random input corpus must be nonempty"
    );
}

#[test]
fn exact_bench_config_error_is_std_error() {
    let err = ExactBenchConfigError::UnorderedRange { min: 5, max: 4 };
    let as_error: &dyn std::error::Error = &err;

    assert_eq!(
        as_error.to_string(),
        "random integer range must be ordered: 5..=4"
    );
    assert!(as_error.source().is_none());
}

#[test]
fn i16_range_single_value_always_draws_that_value() {
    let Ok(range) = I16Range::new(7, 7) else {
        panic!("single-value range should be valid");
    };
    let mut rng = SplitMix64::new(0);

    for _ in 0..32 {
        assert_eq!(rng.next_i16(range), 7);
    }
}

#[test]
fn i16_range_draws_stay_inside_inclusive_bounds() {
    let Ok(range) = I16Range::new(-10, 10) else {
        panic!("ordered range should be valid");
    };
    let mut rng = SplitMix64::new(0xCAFE_F00D);

    for _ in 0..1024 {
        assert!((-10..=10).contains(&rng.next_i16(range)));
    }
}

#[test]
fn splitmix64_sequence_is_stable_for_benchmark_seed() {
    let Ok(range) = I16Range::new(-10, 10) else {
        panic!("ordered range should be valid");
    };
    let mut rng = SplitMix64::new(0xCAFE_F00D);

    let draws = from_fn::<_, 16, _>(|_| rng.next_i16(range));

    assert_eq!(
        draws,
        [-8, 4, -3, -8, 8, -4, 4, -6, -9, -8, 6, 1, 0, 3, -8, -1]
    );
}
