//! Run all independent benchmark fixture gates without Criterion timing.

#![cfg(feature = "bench")]
#![forbid(unsafe_code)]

#[path = "../benches/common/angle.rs"]
mod angle;
#[path = "../benches/common/bench_utils.rs"]
mod bench_utils;

fn validate<const D: usize>() {
    let names = [
        "dense",
        "dense_varied",
        "dense_identical",
        "dense_antipodal",
        "sparse_adjacent",
        "sparse_alternating",
        "sparse_identical",
        "sparse_antipodal",
        "subnormal_axes",
        "maximum_small_angle",
        "near_parallel",
        "near_antipodal",
        "mixed_scale",
        "subnormal",
    ];
    let fixtures = angle::fixtures::<D>();
    assert_eq!(fixtures.len(), names.len());
    for (fixture, name) in fixtures.into_iter().zip(names) {
        assert_eq!(fixture.name(), name);
        assert!(fixture.left().angle(fixture.right()).unwrap().is_finite());
    }
}

#[test]
fn angle_benchmark_inputs_match_independent_references() {
    validate::<3>();
    validate::<4>();
    validate::<5>();
    validate::<6>();
}
