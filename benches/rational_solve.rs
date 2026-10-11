#![forbid(unsafe_code)]

//! Focused structure-sensitive exact solve measurements, outside release reports.

use criterion::Criterion;
use std::hint::black_box;

#[path = "common/bench_utils.rs"]
mod bench_utils;
#[path = "common/rational_solve.rs"]
pub mod rational_solve;

use bench_utils::OrAbort;
use rational_solve::{Structure, solve_input};

fn bench_dimension<const D: usize>(criterion: &mut Criterion) {
    for structure in Structure::ALL {
        for dyadic in [false, true] {
            let input = solve_input::<D>(structure, dyadic);
            let coefficients = if dyadic { "dyadic" } else { "rational" };
            let mut group = criterion.benchmark_group(format!(
                "rational_solve/{}/{coefficients}_d{D}",
                structure.name()
            ));
            group.bench_function("legacy", |b| {
                b.iter(|| black_box(black_box(&input).legacy().or_abort("legacy solve")));
            });
            group.bench_function("construct", |b| {
                b.iter(|| black_box(black_box(&input).construct().or_abort("input construction")));
            });
            group.bench_function("prepared", |b| {
                b.iter(|| {
                    black_box(
                        black_box(input.matrix())
                            .solve(black_box(input.rhs()))
                            .or_abort("prepared solve"),
                    )
                });
            });
            group.bench_function("adapter", |b| {
                b.iter(|| black_box(black_box(&input).adapter().or_abort("adapter solve")));
            });
            group.finish();
        }
    }
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_dimension::<2>(&mut criterion);
    bench_dimension::<3>(&mut criterion);
    bench_dimension::<4>(&mut criterion);
    bench_dimension::<5>(&mut criterion);
    bench_dimension::<6>(&mut criterion);
    bench_dimension::<7>(&mut criterion);
    bench_dimension::<8>(&mut criterion);
    criterion.final_summary();
}
