//! Benchmark comparison between la-stack and nalgebra.
//!
//! Goal: like-for-like comparisons of the operations la-stack supports across several
//! fixed dimensions.
//!
//! Notes:
//! - Determinant is benchmarked via LU on both sides (nalgebra uses closed-forms for 1×1/2×2/3×3).
//! - Matrix infinity norm is the maximum absolute row sum on both sides.

use criterion::Criterion;
use pastey::paste;
use std::hint::black_box;

#[inline]
#[allow(clippy::cast_precision_loss)] // D, r, c are small integers, precision loss is not an issue.
fn matrix_entry<const D: usize>(r: usize, c: usize) -> f64 {
    if r == c {
        // Strict diagonal dominance for stability.
        (r as f64).mul_add(1.0e-3, (D as f64) + 1.0)
    } else {
        // Small, varying off-diagonals.
        0.1 / ((r + c + 1) as f64)
    }
}

#[inline]
fn make_matrix_rows<const D: usize>() -> [[f64; D]; D] {
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
#[allow(clippy::cast_precision_loss)] // i is a small integer, precision loss is not an issue.
fn vector_entry(i: usize, offset: f64) -> f64 {
    (i as f64) + 1.0 + offset
}

#[inline]
fn make_vector_array<const D: usize>(offset: f64) -> [f64; D] {
    let mut data = [0.0; D];

    let mut i = 0;
    while i < D {
        data[i] = vector_entry(i, offset);
        i += 1;
    }

    data
}

#[inline]
fn nalgebra_inf_norm<const D: usize>(m: &nalgebra::SMatrix<f64, D, D>) -> f64 {
    // Infinity norm = max absolute row sum.
    let mut max_row_sum = 0.0;

    let mut r = 0;
    while r < D {
        let mut row_sum = 0.0;
        let mut c = 0;
        while c < D {
            row_sum += m[(r, c)].abs();
            c += 1;
        }
        if row_sum > max_row_sum {
            max_row_sum = row_sum;
        }
        r += 1;
    }

    max_row_sum
}

macro_rules! gen_vs_nalgebra_benches_for_dim {
    ($c:expr, $d:literal) => {
        paste! {{
            // Isolate each dimension's inputs to keep types and captures clean.
            {
                let a = la_stack::Matrix::<$d>::from_rows(make_matrix_rows::<$d>());
                let rhs = la_stack::Vector::<$d>::new(make_vector_array::<$d>(0.0));
                let v1 = la_stack::Vector::<$d>::new(make_vector_array::<$d>(0.0));
                let v2 = la_stack::Vector::<$d>::new(make_vector_array::<$d>(1.0));

                let na = nalgebra::SMatrix::<f64, $d, $d>::from_fn(|r, c| matrix_entry::<$d>(r, c));
                let nrhs = nalgebra::SVector::<f64, $d>::from_fn(|i, _| vector_entry(i, 0.0));
                let nv1 = nalgebra::SVector::<f64, $d>::from_fn(|i, _| vector_entry(i, 0.0));
                let nv2 = nalgebra::SVector::<f64, $d>::from_fn(|i, _| vector_entry(i, 1.0));

                // Precompute LU once for solve-only / det-only benchmarks.
                let a_lu = a
                    .lu(la_stack::DEFAULT_PIVOT_TOL)
                    .expect("matrix should be non-singular");
                let na_lu = na.clone().lu();

                let mut [<group_d $d>] = ($c).benchmark_group(concat!("d", stringify!($d)));

                // === Determinant via LU (factor + det) ===
                [<group_d $d>].bench_function("la_stack_det_via_lu", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(a)
                            .lu(la_stack::DEFAULT_PIVOT_TOL)
                            .expect("matrix should be non-singular");
                        let det = lu.det();
                        black_box(det);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_det_via_lu", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(na.clone()).lu();
                        let det = lu.determinant();
                        black_box(det);
                    });
                });

                // === LU factorization ===
                [<group_d $d>].bench_function("la_stack_lu", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(a)
                            .lu(la_stack::DEFAULT_PIVOT_TOL)
                            .expect("matrix should be non-singular");
                        let _ = black_box(lu);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_lu", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(na.clone()).lu();
                        black_box(lu);
                    });
                });

                // === LU solve (factor + solve) ===
                [<group_d $d>].bench_function("la_stack_lu_solve", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(a)
                            .lu(la_stack::DEFAULT_PIVOT_TOL)
                            .expect("matrix should be non-singular");
                        let x = lu
                            .solve_vec(black_box(rhs))
                            .expect("solve should succeed");
                        let _ = black_box(x);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_lu_solve", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(na.clone()).lu();
                        let x = lu
                            .solve(black_box(&nrhs))
                            .expect("solve should succeed");
                        black_box(x);
                    });
                });

                // === Solve using a precomputed LU ===
                [<group_d $d>].bench_function("la_stack_solve_from_lu", |bencher| {
                    bencher.iter(|| {
                        let x = a_lu
                            .solve_vec(black_box(rhs))
                            .expect("solve should succeed");
                        let _ = black_box(x);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_solve_from_lu", |bencher| {
                    bencher.iter(|| {
                        let x = na_lu
                            .solve(black_box(&nrhs))
                            .expect("solve should succeed");
                        black_box(x);
                    });
                });

                // === Determinant from a precomputed LU ===
                [<group_d $d>].bench_function("la_stack_det_from_lu", |bencher| {
                    bencher.iter(|| {
                        let det = a_lu.det();
                        black_box(det);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_det_from_lu", |bencher| {
                    bencher.iter(|| {
                        let det = na_lu.determinant();
                        black_box(det);
                    });
                });

                // === Vector dot product ===
                [<group_d $d>].bench_function("la_stack_dot", |bencher| {
                    bencher.iter(|| {
                        let result = black_box(v1).dot(black_box(v2));
                        black_box(result);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_dot", |bencher| {
                    bencher.iter(|| {
                        let result = black_box(&nv1).dot(black_box(&nv2));
                        black_box(result);
                    });
                });

                // === Vector norm squared ===
                [<group_d $d>].bench_function("la_stack_norm2_sq", |bencher| {
                    bencher.iter(|| {
                        let result = black_box(v1).norm2_sq();
                        black_box(result);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_norm_squared", |bencher| {
                    bencher.iter(|| {
                        let result = black_box(&nv1).norm_squared();
                        black_box(result);
                    });
                });

                // === Matrix infinity norm (max absolute row sum) ===
                [<group_d $d>].bench_function("la_stack_inf_norm", |bencher| {
                    bencher.iter(|| {
                        let result = black_box(a).inf_norm();
                        black_box(result);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_inf_norm", |bencher| {
                    bencher.iter(|| {
                        let result = nalgebra_inf_norm::<$d>(black_box(&na));
                        black_box(result);
                    });
                });

                [<group_d $d>].finish();
            }
        }}
    };
}

fn main() {
    let mut c = Criterion::default().configure_from_args();

    gen_vs_nalgebra_benches_for_dim!(&mut c, 2);
    gen_vs_nalgebra_benches_for_dim!(&mut c, 3);
    gen_vs_nalgebra_benches_for_dim!(&mut c, 4);
    gen_vs_nalgebra_benches_for_dim!(&mut c, 5);

    gen_vs_nalgebra_benches_for_dim!(&mut c, 8);
    gen_vs_nalgebra_benches_for_dim!(&mut c, 16);
    gen_vs_nalgebra_benches_for_dim!(&mut c, 32);
    gen_vs_nalgebra_benches_for_dim!(&mut c, 64);

    c.final_summary();
}
