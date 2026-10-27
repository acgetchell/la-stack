#![forbid(unsafe_code)]

//! Criterion diagnostics for interval reciprocal, division, and solution enclosures.

use criterion::Criterion;
use la_stack::Interval;
use std::hint::black_box;

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

#[path = "common/interval_solve.rs"]
mod interval_solve;

fn solve_benchmarks<const D: usize>(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group(format!("interval_solve_d{D}"));
    for family in interval_solve::FAMILIES {
        let inputs: Vec<_> = (0..interval_solve::CASES)
            .map(|index| interval_solve::input::<D>(family, index))
            .collect();
        for input in &inputs {
            interval_solve::validate(input);
        }
        let successes = inputs
            .iter()
            .filter(|input| input.expected == interval_solve::Outcome::Enclosed)
            .count();
        eprintln!(
            "interval_solve_d{D}/{family}: {successes}/{} certified",
            inputs.len()
        );
        group.bench_function(family, |b| {
            b.iter(|| {
                for input in &inputs {
                    let _ = black_box(black_box(&input.matrix).solve(black_box(&input.rhs)));
                }
            });
        });
    }
    group.finish();
}

fn division_benchmarks(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("interval_division");
    for (name, numerator, denominator, expected) in [
        ("positive", 1.0, 3.0, [1.0 / 3.0, (1.0_f64 / 3.0).next_up()]),
        (
            "negative",
            1.0,
            -3.0,
            [(-1.0_f64 / 3.0).next_down(), -1.0 / 3.0],
        ),
        (
            "subnormal",
            f64::from_bits(1),
            2.0,
            [0.0, f64::from_bits(1)],
        ),
        (
            "tiny_ratio",
            f64::from_bits(1),
            f64::from_bits(1),
            [1.0, 1.0],
        ),
    ] {
        let n = Interval::point(numerator).or_abort("division numerator");
        let d = Interval::point(denominator).or_abort("division denominator");
        assert_eq!(
            n.try_div(&d).or_abort("division fixture"),
            Interval::try_new(expected[0], expected[1]).or_abort("expected quotient")
        );
        group.bench_function(name, |b| b.iter(|| black_box(&n).try_div(black_box(&d))));
    }
    let divisor = Interval::try_new(-4.0, -2.0).or_abort("reciprocal divisor");
    assert_eq!(
        divisor.try_recip().or_abort("reciprocal fixture"),
        Interval::try_new(-0.5, -0.25).or_abort("expected reciprocal")
    );
    group.bench_function("reciprocal", |b| b.iter(|| black_box(&divisor).try_recip()));
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    division_benchmarks(&mut criterion);
    solve_benchmarks::<1>(&mut criterion);
    solve_benchmarks::<2>(&mut criterion);
    solve_benchmarks::<3>(&mut criterion);
    solve_benchmarks::<4>(&mut criterion);
    solve_benchmarks::<5>(&mut criterion);
    solve_benchmarks::<6>(&mut criterion);
    criterion.final_summary();
}
