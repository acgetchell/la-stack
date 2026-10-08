#![forbid(unsafe_code)]

//! Equivalent stable angle evaluations, with caller/adapter costs distinguished.

use std::hint::black_box;

use criterion::Criterion;
use la_stack::{Vector, VectorAngle};

#[path = "common/angle.rs"]
mod angle;
#[path = "common/bench_utils.rs"]
mod bench_utils;

fn register<const D: usize>(criterion: &mut Criterion) {
    for input in angle::fixtures::<D>() {
        let mut group = criterion.benchmark_group(format!("angle_d{D}/{}", input.name()));
        let left = input.left();
        let right = input.right();
        group.bench_function("stable_control", |b| {
            b.iter(|| {
                angle::stable_control(black_box(left.as_array()), black_box(right.as_array()))
            });
        });
        group.bench_function("prepared_vector", |b| {
            b.iter(|| black_box(left).angle(black_box(right)));
        });
        group.bench_function("borrowed_slice", |b| {
            b.iter(|| {
                black_box(left.as_array().as_slice()).angle(black_box(right.as_array().as_slice()))
            });
        });
        group.bench_function("construct_then_angle", |b| {
            b.iter(|| {
                let a = Vector::try_new(*black_box(left.as_array()))?;
                let b = Vector::try_new(*black_box(right.as_array()))?;
                a.angle(&b)
            });
        });
        group.finish();
    }
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    register::<3>(&mut criterion);
    register::<4>(&mut criterion);
    register::<5>(&mut criterion);
    register::<6>(&mut criterion);
    criterion.final_summary();
}
