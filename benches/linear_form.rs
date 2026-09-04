#![forbid(unsafe_code)]

//! Criterion coverage for certified dot products and affine differences.

use std::hint::black_box;

use criterion::Criterion;

use la_stack::Vector;

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

fn main() {
    let axis =
        Vector::<4>::try_new([2.0, -1.0, 3.0, 4.0]).or_abort("well-separated axis construction");
    let left = Vector::<4>::try_new([4.0, 1.0, 2.0, 3.0])
        .or_abort("well-separated left vector construction");
    let right = Vector::<4>::try_new([1.0, 3.0, 0.0, 2.0])
        .or_abort("well-separated right vector construction");
    let cancellation_left =
        Vector::<4>::try_new([1.0, 1.0, 0.0, 0.0]).or_abort("cancellation left construction");
    let cancellation_right =
        Vector::<4>::try_new([1.0, -1.0, 7.0, -9.0]).or_abort("cancellation right construction");

    let plain_dot = axis
        .dot(&left)
        .or_abort("well-separated plain dot validation");
    assert_eq!(plain_dot.to_bits(), 25.0_f64.to_bits());
    let bounded_dot = axis
        .dot_with_errbound(&left)
        .or_abort("well-separated bounded dot validation")
        .or_abort("well-separated bounded dot certificate");
    assert_eq!(bounded_dot.estimate().to_bits(), plain_dot.to_bits());
    assert!(bounded_dot.lower_bound() > 0.0);

    let bounded_difference = axis
        .dot_difference_with_errbound(&left, &right)
        .or_abort("well-separated bounded difference validation")
        .or_abort("well-separated bounded difference certificate");
    assert_eq!(bounded_difference.estimate().to_bits(), 18.0_f64.to_bits());
    assert!(bounded_difference.lower_bound() > 1.0);

    let inconclusive_dot = cancellation_left
        .dot_with_errbound(&cancellation_right)
        .or_abort("inconclusive bounded dot validation")
        .or_abort("inconclusive bounded dot certificate");
    assert_eq!(inconclusive_dot.estimate().to_bits(), 0.0_f64.to_bits());
    assert!(inconclusive_dot.lower_bound() < 0.0 && inconclusive_dot.upper_bound() > 0.0);

    let inconclusive_difference = axis
        .dot_difference_with_errbound(&left, &left)
        .or_abort("inconclusive bounded difference validation")
        .or_abort("inconclusive bounded difference certificate");
    assert_eq!(
        inconclusive_difference.estimate().to_bits(),
        0.0_f64.to_bits()
    );
    assert!(
        inconclusive_difference.lower_bound() < 0.0 && inconclusive_difference.upper_bound() > 0.0
    );

    let mut criterion = Criterion::default().configure_from_args();
    {
        let mut group = criterion.benchmark_group("linear_form_d4");
        group.bench_function("dot_plain_well_separated", |bencher| {
            bencher.iter(|| {
                let result = black_box(&axis)
                    .dot(black_box(&left))
                    .or_abort("well-separated plain dot");
                let _ = black_box(result);
            });
        });
        group.bench_function("dot_bounded_well_separated", |bencher| {
            bencher.iter(|| {
                let result = black_box(&axis)
                    .dot_with_errbound(black_box(&left))
                    .or_abort("well-separated bounded dot");
                let _ = black_box(result);
            });
        });
        group.bench_function("dot_bounded_inconclusive", |bencher| {
            bencher.iter(|| {
                let result = black_box(&cancellation_left)
                    .dot_with_errbound(black_box(&cancellation_right))
                    .or_abort("inconclusive bounded dot");
                let _ = black_box(result);
            });
        });
        group.bench_function("dot_difference_bounded_well_separated", |bencher| {
            bencher.iter(|| {
                let result = black_box(&axis)
                    .dot_difference_with_errbound(black_box(&left), black_box(&right))
                    .or_abort("well-separated bounded difference");
                let _ = black_box(result);
            });
        });
        group.bench_function("dot_difference_bounded_inconclusive", |bencher| {
            bencher.iter(|| {
                let result = black_box(&axis)
                    .dot_difference_with_errbound(black_box(&left), black_box(&left))
                    .or_abort("inconclusive bounded difference");
                let _ = black_box(result);
            });
        });
        group.finish();
    }
    criterion.final_summary();
}
