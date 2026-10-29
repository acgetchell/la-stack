#![forbid(unsafe_code)]

//! Exact determinant and leading-term extraction for coincident lifted points.

use std::hint::black_box;

use criterion::Criterion;
use la_stack::polynomial_determinant;

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

#[path = "common/polynomial.rs"]
mod polynomial;
use polynomial::{BUDGET, coincident_lifted};

fn register<const N: usize>(criterion: &mut Criterion) {
    let input = coincident_lifted::<N>();
    let mut group = criterion.benchmark_group(format!("polynomial/coincident_lifted/order{N}"));
    group.bench_function("determinant", |bencher| {
        bencher.iter(|| {
            let determinant = polynomial_determinant(black_box(&input.matrix), black_box(BUDGET))
                .or_abort("validated fixture fits the consumer budget");
            black_box(determinant.leading_term());
            black_box(determinant)
        });
    });
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    register::<4>(&mut criterion);
    register::<5>(&mut criterion);
    register::<6>(&mut criterion);
    register::<7>(&mut criterion);
    register::<8>(&mut criterion);
    criterion.final_summary();
}
