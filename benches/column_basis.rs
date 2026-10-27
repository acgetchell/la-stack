#![forbid(unsafe_code)]

//! Rectangular column selection with independently known ranks and pivot sets.

use std::{array::from_fn, hint::black_box};

use criterion::Criterion;
use la_stack::{Tolerance, column_basis};

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

fn register<const R: usize, const C: usize>(c: &mut Criterion) {
    for (name, scale, deficient) in [
        ("full_rank", 1.0, false),
        ("rank_one", 1.0, true),
        ("near_dependent", 1.0, false),
        ("small", 2.0_f64.powi(-900), false),
        ("large", 2.0_f64.powi(900), false),
    ] {
        // A leading identity (or all-ones rank-one matrix) gives known rank
        // without using another floating elimination as the oracle.
        let rows: [[f64; C]; R] = from_fn(|r| {
            from_fn(|col| {
                if name == "near_dependent" {
                    return if r == col {
                        1.0 + 2.0_f64.powi(-30)
                    } else {
                        1.0
                    };
                }
                if deficient || r == col || col >= R {
                    scale
                } else {
                    0.0
                }
            })
        });
        let tol = Tolerance::try_new(scale * 1e-12).or_abort("basis tolerance");
        let expected: Vec<_> = (0..if deficient { 1 } else { R.min(C) }).collect();
        assert_eq!(
            column_basis(&rows, tol)
                .or_abort("basis validation")
                .columns(),
            expected
        );
        c.bench_function(&format!("column_basis/{R}x{C}/{name}"), |b| {
            b.iter(|| {
                black_box(column_basis(black_box(&rows), black_box(tol)).or_abort("column basis"))
            });
        });
    }
}

fn affine<const R: usize, const C: usize>(c: &mut Criterion) {
    let rows: [[f64; C]; R] =
        from_fn(|r| from_fn(|col| f64::from(col == C - 1 || (r > 0 && col == r - 1))));
    let expected: Vec<_> = (0..R - 1).chain([C - 1]).collect();
    let tol = Tolerance::try_new(1e-12).or_abort("affine tolerance");
    assert_eq!(
        column_basis(&rows, tol)
            .or_abort("affine validation")
            .columns(),
        expected
    );
    c.bench_function(&format!("column_basis/{R}x{C}/shared_face"), |b| {
        b.iter(|| {
            black_box(column_basis(black_box(&rows), black_box(tol)).or_abort("affine basis"))
        });
    });
}

fn main() {
    let mut c = Criterion::default().configure_from_args();
    register::<2, 2>(&mut c);
    register::<3, 3>(&mut c);
    register::<4, 4>(&mut c);
    register::<5, 5>(&mut c);
    register::<5, 2>(&mut c);
    register::<4, 6>(&mut c);
    register::<5, 8>(&mut c);
    register::<6, 10>(&mut c);
    register::<7, 12>(&mut c);
    register::<8, 14>(&mut c);
    affine::<2, 3>(&mut c);
    affine::<3, 4>(&mut c);
    affine::<4, 5>(&mut c);
    affine::<5, 6>(&mut c);
    affine::<6, 7>(&mut c);
    c.final_summary();
}
