#![forbid(unsafe_code)]

//! Exact dot products with analytically checked reciprocal and cancellation inputs.

use std::hint::black_box;

use criterion::Criterion;
use la_stack::{BigInt, BigRational, RationalVector, rational_dot};

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

fn register<const D: usize>(c: &mut Criterion) {
    for (name, bits, alternating) in [
        ("non_dyadic", 0, false),
        ("cancellation", 0, true),
        ("large", 512, false),
    ] {
        let left = RationalVector::<D>::try_from_fn(|i| {
            BigRational::new((BigInt::from(i + 1) << bits) + 1, BigInt::from(2 * i + 3))
        })
        .or_abort("rational dot input");
        let right = RationalVector::<D>::try_from_fn(|i| {
            let inverse = left.as_array()[i].recip();
            if alternating && i % 2 == 1 {
                -inverse
            } else {
                inverse
            }
        })
        .or_abort("rational dot reciprocal");
        let expected = BigRational::from_integer(BigInt::from(if alternating { D % 2 } else { D }));
        assert_eq!(left.dot(&right), expected);
        assert_eq!(
            rational_dot(left.as_array(), right.as_array()).or_abort("adapter validation"),
            expected
        );
        let mut group = c.benchmark_group(format!("rational_dot/{name}/d{D}"));
        group.bench_function("prepared", |b| {
            b.iter(|| black_box(black_box(&left).dot(black_box(&right))));
        });
        group.bench_function("slice", |b| {
            b.iter(|| {
                black_box(
                    rational_dot(black_box(left.as_array()), black_box(right.as_array()))
                        .or_abort("rational dot adapter"),
                )
            });
        });
        group.finish();
    }
}

fn main() {
    let mut c = Criterion::default().configure_from_args();
    register::<2>(&mut c);
    register::<3>(&mut c);
    register::<4>(&mut c);
    register::<5>(&mut c);
    register::<6>(&mut c);
    register::<7>(&mut c);
    register::<8>(&mut c);
    register::<14>(&mut c);
    register::<22>(&mut c);
    c.final_summary();
}
