#![forbid(unsafe_code)]

//! Criterion coverage for conclusive and inconclusive interval determinants.

use std::hint::black_box;

use criterion::Criterion;

use la_stack::{Interval, IntervalDeterminantSign, IntervalMatrix, LaError};

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

/// Assemble the relative-coordinate lifted matrix for a tetrahedral in-sphere
/// predicate whose interval determinant is conclusively negative.
fn conclusive_lifted_4x4() -> Result<IntervalMatrix<4>, LaError> {
    let x = Interval::try_from_subtraction(0.1, 0.0)?;
    let y = Interval::try_from_subtraction(0.1, 0.0)?;
    let z = Interval::try_from_subtraction(0.1, 0.0)?;
    let lifted = x
        .try_square()?
        .try_add(&y.try_square()?)?
        .try_add(&z.try_square()?)?;

    Ok(IntervalMatrix::from_rows([
        [Interval::ONE, Interval::ZERO, Interval::ZERO, Interval::ONE],
        [Interval::ZERO, Interval::ONE, Interval::ZERO, Interval::ONE],
        [Interval::ZERO, Interval::ZERO, Interval::ONE, Interval::ONE],
        [x, y, z, lifted],
    ]))
}

/// Assemble a lifted boundary case whose final coefficient retains one ULP of
/// expression uncertainty, forcing an inconclusive interval sign.
fn inconclusive_lifted_4x4() -> Result<IntervalMatrix<4>, LaError> {
    Ok(IntervalMatrix::from_rows([
        [Interval::ONE, Interval::ZERO, Interval::ZERO, Interval::ONE],
        [Interval::ZERO, Interval::ONE, Interval::ZERO, Interval::ONE],
        [Interval::ZERO, Interval::ZERO, Interval::ONE, Interval::ONE],
        [
            Interval::ONE,
            Interval::ONE,
            Interval::ONE,
            Interval::try_new(3.0_f64.next_down(), 3.0_f64.next_up())?,
        ],
    ]))
}

/// Assemble a six-coordinate lifted predicate matrix to exercise the maximum
/// supported D=7 subset-DP workload.
fn conclusive_lifted_7x7() -> Result<IntervalMatrix<7>, LaError> {
    let mut matrix = IntervalMatrix::zero();
    for index in 0..6 {
        matrix.set(index, index, Interval::ONE)?;
        matrix.set(index, 6, Interval::ONE)?;
    }

    let relative = Interval::try_from_subtraction(0.125, 0.0)?;
    let square = relative.try_square()?;
    let mut lifted = Interval::ZERO;
    for column in 0..6 {
        matrix.set(6, column, relative)?;
        lifted = lifted.try_add(&square)?;
    }
    matrix.set(6, 6, lifted)?;
    Ok(matrix)
}

fn main() {
    let conclusive_4 =
        conclusive_lifted_4x4().or_abort("conclusive D=4 interval fixture construction");
    let inconclusive_4 =
        inconclusive_lifted_4x4().or_abort("inconclusive D=4 interval fixture construction");
    let conclusive_7 =
        conclusive_lifted_7x7().or_abort("conclusive D=7 interval fixture construction");

    assert_eq!(
        conclusive_4
            .det_sign()
            .or_abort("conclusive D=4 interval fixture validation"),
        IntervalDeterminantSign::Negative,
    );
    assert_eq!(
        inconclusive_4
            .det_sign()
            .or_abort("inconclusive D=4 interval fixture validation"),
        IntervalDeterminantSign::Inconclusive,
    );
    assert_eq!(
        conclusive_7
            .det_sign()
            .or_abort("conclusive D=7 interval fixture validation"),
        IntervalDeterminantSign::Negative,
    );

    let mut criterion = Criterion::default().configure_from_args();
    {
        let mut group = criterion.benchmark_group("interval_det_sign");
        group.bench_function("d4_conclusive_lifted", |bencher| {
            bencher.iter(|| {
                let sign = black_box(&conclusive_4)
                    .det_sign()
                    .or_abort("D=4 conclusive interval determinant");
                let _ = black_box(sign);
            });
        });
        group.bench_function("d4_inconclusive_lifted", |bencher| {
            bencher.iter(|| {
                let sign = black_box(&inconclusive_4)
                    .det_sign()
                    .or_abort("D=4 inconclusive interval determinant");
                let _ = black_box(sign);
            });
        });
        group.bench_function("d7_conclusive_lifted", |bencher| {
            bencher.iter(|| {
                let sign = black_box(&conclusive_7)
                    .det_sign()
                    .or_abort("D=7 conclusive interval determinant");
                let _ = black_box(sign);
            });
        });
        group.finish();
    }
    criterion.final_summary();
}
