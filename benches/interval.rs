#![forbid(unsafe_code)]

//! Criterion coverage for conclusive and inconclusive interval determinants.

use std::hint::black_box;

use criterion::Criterion;

use la_stack::{Interval, IntervalDeterminantSign, IntervalMatrix, LaError};

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

fn scalar_benchmarks(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("interval_scalar");
    for (name, left, right) in [
        ("point", 0.1, 0.3),
        ("exact", 0.125, 2.0),
        ("zero", 0.0, 0.3),
        ("subnormal", f64::from_bits(1), 0.5),
        ("cancellation", 1.0, -1.0),
    ] {
        let x = Interval::point(left).or_abort("scalar left");
        let y = Interval::point(right).or_abort("scalar right");
        // Each result must enclose its rounded value; fixture checks stay out
        // of timing, including the underflow-to-zero square.
        for (result, rounded) in [
            (x.try_add(&y), left + right),
            (x.try_mul(&y), left * right),
            (x.try_square(), left * left),
        ] {
            let result = result.or_abort("scalar enclosure");
            assert!(result.lower() <= rounded && rounded <= result.upper());
        }
        group.bench_function(format!("{name}/add"), |b| {
            b.iter(|| black_box(&x).try_add(black_box(&y)));
        });
        group.bench_function(format!("{name}/multiply"), |b| {
            b.iter(|| black_box(&x).try_mul(black_box(&y)));
        });
        group.bench_function(format!("{name}/square"), |b| {
            b.iter(|| black_box(&x).try_square());
        });
    }
    group.finish();
}

/// Block matrix [I, 1; t, (N-1)t²], with determinant (N-1)(t²-t).
fn lifted<const N: usize>(coordinate: f64) -> Result<IntervalMatrix<N>, LaError> {
    let mut rows = [[Interval::ZERO; N]; N];
    for (index, row) in rows.iter_mut().enumerate().take(N - 1) {
        row[index] = Interval::ONE;
        row[N - 1] = Interval::ONE;
    }
    let mut norm = Interval::ZERO;
    for entry in &mut rows[N - 1][..N - 1] {
        *entry = Interval::try_from_subtraction(coordinate, 0.0)?;
        norm = norm.try_add(&entry.try_square()?)?;
    }
    rows[N - 1][N - 1] = norm;
    Ok(IntervalMatrix::from_rows(rows))
}

fn dimension_benchmarks<const N: usize>(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group(format!("interval_d{N}"));
    let coordinate = 0.1;
    let sparse = lifted::<N>(coordinate).or_abort("lifted fixture");
    assert_eq!(
        sparse.det_sign().or_abort("lifted sign"),
        IntervalDeterminantSign::Negative
    );
    // I + 11ᵀ has determinant N+1. Widen a repeated row to obtain an
    // inconclusive dense enclosure that still includes a singular matrix.
    let dense_rows = core::array::from_fn(|row| {
        core::array::from_fn(|column| {
            if row == column {
                Interval::point(2.0).or_abort("dense diagonal")
            } else {
                Interval::ONE
            }
        })
    });
    let dense = IntervalMatrix::<N>::from_rows(dense_rows);
    assert_eq!(
        dense.det_sign().or_abort("dense sign"),
        IntervalDeterminantSign::Positive
    );
    let mut uncertain_rows = dense_rows;
    uncertain_rows[N - 1] = dense_rows[0];
    uncertain_rows[N - 1][0] = Interval::try_new(1.99, 2.01).or_abort("uncertain entry");
    let uncertain = IntervalMatrix::from_rows(uncertain_rows);
    assert_eq!(
        uncertain.det_sign().or_abort("uncertain sign"),
        IntervalDeterminantSign::Inconclusive
    );
    group.bench_function("lifted_assembly", |b| {
        b.iter(|| lifted::<N>(black_box(coordinate)));
    });
    group.bench_function("lifted_assembly_and_sign", |b| {
        b.iter(|| lifted::<N>(black_box(coordinate)).and_then(|matrix| matrix.det_sign()));
    });
    for (name, matrix) in [
        ("sparse_sign", sparse),
        ("dense_sign", dense),
        ("dense_inconclusive", uncertain),
    ] {
        group.bench_function(name, |b| b.iter(|| black_box(&matrix).det_sign()));
    }
    group.finish();
}

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
    scalar_benchmarks(&mut criterion);
    dimension_benchmarks::<3>(&mut criterion);
    dimension_benchmarks::<4>(&mut criterion);
    dimension_benchmarks::<5>(&mut criterion);
    dimension_benchmarks::<6>(&mut criterion);
    dimension_benchmarks::<7>(&mut criterion);
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
