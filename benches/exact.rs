//! Benchmarks for exact arithmetic operations.
//!
//! These benchmarks measure the performance of the `exact` feature's
//! arbitrary-precision methods across dimensions D=2..5 (the primary
//! target for geometric predicates).

use criterion::Criterion;
use la_stack::{Matrix, Vector};
use pastey::paste;
use std::hint::black_box;

#[inline]
#[allow(clippy::cast_precision_loss)]
const fn matrix_entry<const D: usize>(r: usize, c: usize) -> f64 {
    if r == c {
        (r as f64).mul_add(1.0e-3, (D as f64) + 1.0)
    } else {
        0.1 / ((r + c + 1) as f64)
    }
}

#[inline]
const fn make_matrix_rows<const D: usize>() -> [[f64; D]; D] {
    let mut rows = [[0.0; D]; D];
    let mut r = 0;
    while r < D {
        let mut c = 0;
        while c < D {
            rows[r][c] = matrix_entry::<D>(r, c);
            c += 1;
        }
        r += 1;
    }
    rows
}

#[inline]
#[allow(clippy::cast_precision_loss)]
fn make_vector_array<const D: usize>() -> [f64; D] {
    let mut data = [0.0; D];
    let mut i = 0;
    while i < D {
        data[i] = (i as f64) + 1.0;
        i += 1;
    }
    data
}

/// Near-singular matrix: base singular matrix + tiny perturbation.
/// This forces the exact Bareiss fallback in `det_sign_exact` (the fast
/// f64 filter cannot resolve the sign).
#[inline]
fn near_singular_3x3() -> Matrix<3> {
    let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
    Matrix::<3>::from_rows([
        [1.0 + perturbation, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ])
}

macro_rules! gen_exact_benches_for_dim {
    ($c:expr, $d:literal) => {
        paste! {{
            let a = Matrix::<$d>::from_rows(make_matrix_rows::<$d>());
            let rhs = Vector::<$d>::new(make_vector_array::<$d>());

            let mut [<group_d $d>] = ($c).benchmark_group(concat!("exact_d", stringify!($d)));

            // === f64 baselines ===
            [<group_d $d>].bench_function("det", |bencher| {
                bencher.iter(|| {
                    let det = black_box(a)
                        .det(la_stack::DEFAULT_PIVOT_TOL)
                        .expect("should not fail");
                    black_box(det);
                });
            });

            [<group_d $d>].bench_function("det_direct", |bencher| {
                bencher.iter(|| {
                    let det = black_box(a).det_direct();
                    black_box(det);
                });
            });

            // === det_exact (BigRational result) ===
            [<group_d $d>].bench_function("det_exact", |bencher| {
                bencher.iter(|| {
                    let det = black_box(a).det_exact().expect("should not fail");
                    black_box(det);
                });
            });

            // === det_exact_f64 (exact → f64) ===
            [<group_d $d>].bench_function("det_exact_f64", |bencher| {
                bencher.iter(|| {
                    let det = black_box(a).det_exact_f64().expect("should not fail");
                    black_box(det);
                });
            });

            // === det_sign_exact (adaptive: fast filter + exact fallback) ===
            [<group_d $d>].bench_function("det_sign_exact", |bencher| {
                bencher.iter(|| {
                    let sign = black_box(a).det_sign_exact().expect("should not fail");
                    black_box(sign);
                });
            });

            // === solve_exact (BigRational result) ===
            [<group_d $d>].bench_function("solve_exact", |bencher| {
                bencher.iter(|| {
                    let x = black_box(a).solve_exact(black_box(rhs)).expect("should not fail");
                    black_box(x);
                });
            });

            // === solve_exact_f64 (exact → f64) ===
            [<group_d $d>].bench_function("solve_exact_f64", |bencher| {
                bencher.iter(|| {
                    let x = black_box(a).solve_exact_f64(black_box(rhs)).expect("should not fail");
                    black_box(x);
                });
            });

            [<group_d $d>].finish();
        }};
    };
}

fn main() {
    let mut c = Criterion::default().configure_from_args();

    #[allow(unused_must_use)]
    {
        gen_exact_benches_for_dim!(&mut c, 2);
        gen_exact_benches_for_dim!(&mut c, 3);
        gen_exact_benches_for_dim!(&mut c, 4);
        gen_exact_benches_for_dim!(&mut c, 5);
    }

    // Near-singular 3×3: forces Bareiss fallback in det_sign_exact.
    {
        let m = near_singular_3x3();
        let mut group = c.benchmark_group("exact_near_singular_3x3");

        group.bench_function("det_sign_exact", |bencher| {
            bencher.iter(|| {
                let sign = black_box(m).det_sign_exact().expect("should not fail");
                black_box(sign);
            });
        });

        group.bench_function("det_exact", |bencher| {
            bencher.iter(|| {
                let det = black_box(m).det_exact().expect("should not fail");
                black_box(det);
            });
        });

        group.finish();
    }

    c.final_summary();
}
