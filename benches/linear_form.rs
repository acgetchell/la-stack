#![forbid(unsafe_code)]

//! Criterion coverage for certified dot products and affine differences.

use std::hint::black_box;

use criterion::Criterion;

use la_stack::{ArithmeticOperation, Interval, LaError, Vector};

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

fn projection_batch<const D: usize>(
    axis: &Vector<D>,
    shared: &Vector<D>,
    positive: &[Vector<D>; D],
    negative: &[Vector<D>; D],
) -> Result<bool, LaError> {
    let Some(common) = axis.dot_with_errbound(shared)? else {
        return Ok(false);
    };
    for vertex in positive {
        let Some(value) = axis.dot_with_errbound(vertex)? else {
            return Ok(false);
        };
        if Interval::try_from_subtraction(value.lower_bound(), common.upper_bound())?.lower() <= 1.0
        {
            return Ok(false);
        }
    }
    for vertex in negative {
        let Some(value) = axis.dot_with_errbound(vertex)? else {
            return Ok(false);
        };
        if Interval::try_from_subtraction(common.lower_bound(), value.upper_bound())?.lower() <= 1.0
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn prepared_benchmarks<const D: usize>(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group(format!("linear_form_d{D}_prepared"));
    for (name, axis, left, right, expected) in [
        ("dense", [2.0; D], [3.0; D], [1.0; D], Some(true)),
        ("dense_inexact", [0.1; D], [0.3; D], [0.2; D], Some(true)),
        (
            "sparse",
            [2.0; D],
            core::array::from_fn(|i| if i == 0 { 1.0 } else { 0.0 }),
            [0.0; D],
            Some(true),
        ),
        (
            "cancellation",
            core::array::from_fn(|i| {
                if i == 0 {
                    1.0
                } else if i == 1 {
                    -1.0
                } else {
                    0.0
                }
            }),
            [1.0; D],
            [1.0; D],
            Some(false),
        ),
        (
            "underflow",
            [f64::MIN_POSITIVE; D],
            [0.5; D],
            [0.0; D],
            None,
        ),
    ] {
        let axis = Vector::try_new(axis).or_abort("axis");
        let left = Vector::try_new(left).or_abort("left");
        let right = Vector::try_new(right).or_abort("right");
        for value in [
            axis.dot_with_errbound(&left),
            axis.dot_difference_with_errbound(&left, &right),
        ] {
            let value = value.or_abort("prepared arithmetic");
            assert_eq!(value.map(|v| v.lower_bound() > 0.0), expected);
        }
        group.bench_function(format!("{name}/dot"), |b| {
            b.iter(|| black_box(&axis).dot_with_errbound(black_box(&left)));
        });
        group.bench_function(format!("{name}/dot_endpoints"), |b| {
            b.iter(|| {
                black_box(&axis)
                    .dot_with_errbound(black_box(&left))
                    .map(|value| value.map(|bound| (bound.lower_bound(), bound.upper_bound())))
            });
        });
        group.bench_function(format!("{name}/difference"), |b| {
            b.iter(|| {
                black_box(&axis).dot_difference_with_errbound(black_box(&left), black_box(&right))
            });
        });
        group.bench_function(format!("{name}/difference_endpoints"), |b| {
            b.iter(|| {
                black_box(&axis)
                    .dot_difference_with_errbound(black_box(&left), black_box(&right))
                    .map(|value| value.map(|bound| (bound.lower_bound(), bound.upper_bound())))
            });
        });
    }
    let huge = Vector::try_new([f64::MAX; D]).or_abort("range fixture");
    let two = Vector::try_new([2.0; D]).or_abort("range multiplier");
    assert_eq!(
        huge.dot_with_errbound(&two),
        Err(LaError::non_finite_computation_step(
            ArithmeticOperation::VectorDotProduct,
            0
        ))
    );
    group.bench_function("overflow/dot", |b| {
        b.iter(|| black_box(&huge).dot_with_errbound(black_box(&two)));
    });
    group.finish();
}

fn batch_benchmarks<const D: usize>(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group(format!("linear_form_d{D}_batch"));
    for (name, shared_rows) in [
        ("origin", [0.0; D]),
        (
            "translated",
            core::array::from_fn(|i| {
                f64::from(u32::try_from(i + 1).or_abort("small dimension")) / 8.0
            }),
        ),
    ] {
        let axis_rows = [2048.0; D];
        let positive_rows: [[f64; D]; D] = core::array::from_fn(|i| {
            core::array::from_fn(|j| shared_rows[j] + if i == j { 1.0 } else { 0.0 })
        });
        let negative_rows: [[f64; D]; D] = core::array::from_fn(|i| {
            core::array::from_fn(|j| shared_rows[j] - if i == j { 1.0 } else { 0.0 })
        });
        let axis = Vector::try_new(axis_rows).or_abort("batch axis");
        let shared = Vector::try_new(shared_rows).or_abort("batch shared");
        let positive = positive_rows.map(|row| Vector::try_new(row).or_abort("positive vertex"));
        let negative = negative_rows.map(|row| Vector::try_new(row).or_abort("negative vertex"));
        assert!(
            projection_batch(&axis, &shared, &positive, &negative).or_abort("batch validation")
        );
        group.bench_function(format!("{name}/construction"), |b| {
            b.iter(|| Vector::try_new(black_box(shared_rows)));
        });
        group.bench_function(format!("{name}/prepared"), |b| {
            b.iter(|| {
                projection_batch(
                    black_box(&axis),
                    black_box(&shared),
                    black_box(&positive),
                    black_box(&negative),
                )
            });
        });
        group.bench_function(format!("{name}/with_construction"), |b| {
            b.iter(|| {
                let axis = Vector::try_new(black_box(axis_rows)).or_abort("axis construction");
                let shared =
                    Vector::try_new(black_box(shared_rows)).or_abort("shared construction");
                let positive = black_box(positive_rows)
                    .map(|row| Vector::try_new(row).or_abort("positive construction"));
                let negative = black_box(negative_rows)
                    .map(|row| Vector::try_new(row).or_abort("negative construction"));
                projection_batch(&axis, &shared, &positive, &negative)
            });
        });
    }
    group.finish();
}

fn dimension_benchmarks<const D: usize>(criterion: &mut Criterion) {
    prepared_benchmarks::<D>(criterion);
    batch_benchmarks::<D>(criterion);
}

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
    dimension_benchmarks::<2>(&mut criterion);
    dimension_benchmarks::<3>(&mut criterion);
    dimension_benchmarks::<4>(&mut criterion);
    dimension_benchmarks::<5>(&mut criterion);
    dimension_benchmarks::<6>(&mut criterion);
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
