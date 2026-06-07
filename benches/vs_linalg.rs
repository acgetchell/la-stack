//! Benchmark comparison between la-stack and other Rust linear algebra crates.
//!
//! Goal: like-for-like comparisons of the operations la-stack supports across several
//! fixed dimensions.
//!
//! Notes:
//! - Determinant is benchmarked via LU on all sides (nalgebra uses closed-forms for 1×1/2×2/3×3).
//! - Matrix infinity norm is the maximum absolute row sum on all sides.

use std::fmt::Display;
use std::hint::black_box;

use criterion::Criterion;
use faer::linalg::solvers::Solve;
use faer::{Mat, Side};
use nalgebra::{SMatrix, SVector};
use pastey::paste;

use la_stack::{DEFAULT_SINGULAR_TOL, Matrix, Vector};

mod common {
    pub mod vs_linalg;
}

use common::vs_linalg::{
    faer_det_from_ldlt, faer_det_from_partial_piv_lu, make_matrix_rows, make_vector_array,
    matrix_entry, nalgebra_inf_norm, vector_entry,
};

/// Return a successful benchmark operation result or panic with the named operation.
fn require_ok<T, E: Display>(result: Result<T, E>, operation: &str) -> T {
    match result {
        Ok(value) => value,
        Err(err) => panic!("{operation} failed: {err}"),
    }
}

/// Return a present third-party benchmark result or panic with the named operation.
fn require_some<T>(value: Option<T>, operation: &str) -> T {
    value.unwrap_or_else(|| panic!("{operation} returned no result"))
}

macro_rules! define_vs_linalg_benches_for_dim {
    ($fn_name:ident, $d:literal) => {
        paste! {
            #[allow(clippy::too_many_lines)]
            fn $fn_name(c: &mut Criterion) {
            // Isolate each dimension's inputs to keep types and captures clean.
            {
                let a = require_ok(
                    Matrix::<$d>::try_from_rows(make_matrix_rows::<$d>()),
                    "la_stack matrix construction",
                );
                let rhs = require_ok(
                    Vector::<$d>::try_new(make_vector_array::<$d>(0.0)),
                    "la_stack RHS vector construction",
                );
                let v1 = require_ok(
                    Vector::<$d>::try_new(make_vector_array::<$d>(0.0)),
                    "la_stack vector construction",
                );
                let v2 = require_ok(
                    Vector::<$d>::try_new(make_vector_array::<$d>(1.0)),
                    "la_stack vector construction",
                );
                let na = SMatrix::<f64, $d, $d>::from_fn(|r, c| matrix_entry::<$d>(r, c));
                let nrhs = SVector::<f64, $d>::from_fn(|i, _| vector_entry(i, 0.0));
                let nv1 = SVector::<f64, $d>::from_fn(|i, _| vector_entry(i, 0.0));
                let nv2 = SVector::<f64, $d>::from_fn(|i, _| vector_entry(i, 1.0));

                let fa = Mat::<f64>::from_fn($d, $d, |r, c| matrix_entry::<$d>(r, c));
                let frhs = Mat::<f64>::from_fn($d, 1, |i, _| vector_entry(i, 0.0));
                let fv1 = Mat::<f64>::from_fn($d, 1, |i, _| vector_entry(i, 0.0));
                let fv2 = Mat::<f64>::from_fn($d, 1, |i, _| vector_entry(i, 1.0));

                // Precompute LU once for solve-only / det-only benchmarks.
                let a_lu = require_ok(a.lu(DEFAULT_SINGULAR_TOL), "precomputed la_stack LU");
                let a_ldlt = require_ok(a.ldlt(DEFAULT_SINGULAR_TOL), "precomputed la_stack LDLT");
                let na_lu = na.clone().lu();
                let na_cholesky = require_some(na.clone().cholesky(), "precomputed nalgebra Cholesky");
                let fa_lu = fa.partial_piv_lu();
                let fa_ldlt = require_ok(fa.ldlt(Side::Lower), "precomputed faer LDLT");

                let mut [<group_d $d>] = c.benchmark_group(concat!("d", stringify!($d)));

                // === Determinant via LU (factor + det) ===
                [<group_d $d>].bench_function("la_stack_det_via_lu", |bencher| {
                    bencher.iter(|| {
                        let lu = require_ok(
                            black_box(a).lu(DEFAULT_SINGULAR_TOL),
                            "la_stack LU factorization",
                        );
                        let det = require_ok(lu.det(), "la_stack LU determinant");
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

                [<group_d $d>].bench_function("faer_det_via_lu", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(&fa).partial_piv_lu();
                        let det = faer_det_from_partial_piv_lu(&lu);
                        black_box(det);
                    });
                });

                // === Determinant via det() (closed-form for D≤4, LU for D≥5) ===
                [<group_d $d>].bench_function("la_stack_det", |bencher| {
                    bencher.iter(|| {
                        let det = require_ok(black_box(a).det(), "la_stack determinant");
                        black_box(det);
                    });
                });

                // === LU factorization ===
                [<group_d $d>].bench_function("la_stack_lu", |bencher| {
                    bencher.iter(|| {
                        let lu = require_ok(
                            black_box(a).lu(DEFAULT_SINGULAR_TOL),
                            "la_stack LU factorization",
                        );
                        let _ = black_box(lu);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_lu", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(na.clone()).lu();
                        black_box(lu);
                    });
                });

                [<group_d $d>].bench_function("faer_lu", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(&fa).partial_piv_lu();
                        black_box(lu);
                    });
                });

                // === SPD factorization (LDLT / Cholesky) ===
                [<group_d $d>].bench_function("la_stack_ldlt", |bencher| {
                    bencher.iter(|| {
                        let ldlt = require_ok(
                            black_box(a).ldlt(DEFAULT_SINGULAR_TOL),
                            "la_stack LDLT factorization",
                        );
                        let _ = black_box(ldlt);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_cholesky", |bencher| {
                    bencher.iter(|| {
                        let chol = require_some(
                            black_box(na.clone()).cholesky(),
                            "nalgebra Cholesky factorization",
                        );
                        black_box(chol);
                    });
                });

                [<group_d $d>].bench_function("faer_ldlt", |bencher| {
                    bencher.iter(|| {
                        let ldlt = require_ok(
                            black_box(&fa).ldlt(Side::Lower),
                            "faer LDLT factorization",
                        );
                        black_box(ldlt);
                    });
                });

                // === LU solve (factor + solve) ===
                [<group_d $d>].bench_function("la_stack_lu_solve", |bencher| {
                    bencher.iter(|| {
                        let lu = require_ok(
                            black_box(a).lu(DEFAULT_SINGULAR_TOL),
                            "la_stack LU factorization",
                        );
                        let x = require_ok(
                            lu.solve(black_box(rhs)),
                            "la_stack LU solve",
                        );
                        let _ = black_box(x);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_lu_solve", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(na.clone()).lu();
                        let x = require_some(lu.solve(black_box(&nrhs)), "nalgebra LU solve");
                        black_box(x);
                    });
                });

                [<group_d $d>].bench_function("faer_lu_solve", |bencher| {
                    bencher.iter(|| {
                        let lu = black_box(&fa).partial_piv_lu();
                        let x = lu.solve(black_box(&frhs));
                        black_box(x);
                    });
                });

                // === SPD solve (factor + solve) ===
                [<group_d $d>].bench_function("la_stack_ldlt_solve", |bencher| {
                    bencher.iter(|| {
                        let ldlt = require_ok(
                            black_box(a).ldlt(DEFAULT_SINGULAR_TOL),
                            "la_stack LDLT factorization",
                        );
                        let x = require_ok(
                            ldlt.solve(black_box(rhs)),
                            "la_stack LDLT solve",
                        );
                        let _ = black_box(x);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_cholesky_solve", |bencher| {
                    bencher.iter(|| {
                        let chol = require_some(
                            black_box(na.clone()).cholesky(),
                            "nalgebra Cholesky factorization",
                        );
                        let x = chol.solve(black_box(&nrhs));
                        black_box(x);
                    });
                });

                [<group_d $d>].bench_function("faer_ldlt_solve", |bencher| {
                    bencher.iter(|| {
                        let ldlt = require_ok(
                            black_box(&fa).ldlt(Side::Lower),
                            "faer LDLT factorization",
                        );
                        let x = ldlt.solve(black_box(&frhs));
                        black_box(x);
                    });
                });

                // === Solve using a precomputed LU ===
                [<group_d $d>].bench_function("la_stack_solve_from_lu", |bencher| {
                    bencher.iter(|| {
                        let x = require_ok(
                            a_lu.solve(black_box(rhs)),
                            "precomputed la_stack LU solve",
                        );
                        let _ = black_box(x);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_solve_from_lu", |bencher| {
                    bencher.iter(|| {
                        let x = require_some(
                            na_lu.solve(black_box(&nrhs)),
                            "precomputed nalgebra LU solve",
                        );
                        black_box(x);
                    });
                });

                [<group_d $d>].bench_function("faer_solve_from_lu", |bencher| {
                    bencher.iter(|| {
                        let x = fa_lu.solve(black_box(&frhs));
                        black_box(x);
                    });
                });

                // === Solve using a precomputed SPD factorization ===
                [<group_d $d>].bench_function("la_stack_solve_from_ldlt", |bencher| {
                    bencher.iter(|| {
                        let x = require_ok(
                            a_ldlt.solve(black_box(rhs)),
                            "precomputed la_stack LDLT solve",
                        );
                        let _ = black_box(x);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_solve_from_cholesky", |bencher| {
                    bencher.iter(|| {
                        let x = na_cholesky.solve(black_box(&nrhs));
                        black_box(x);
                    });
                });

                [<group_d $d>].bench_function("faer_solve_from_ldlt", |bencher| {
                    bencher.iter(|| {
                        let x = fa_ldlt.solve(black_box(&frhs));
                        black_box(x);
                    });
                });

                // === Determinant from a precomputed LU ===
                [<group_d $d>].bench_function("la_stack_det_from_lu", |bencher| {
                    bencher.iter(|| {
                        let det = require_ok(a_lu.det(), "precomputed la_stack LU determinant");
                        black_box(det);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_det_from_lu", |bencher| {
                    bencher.iter(|| {
                        let det = na_lu.determinant();
                        black_box(det);
                    });
                });

                [<group_d $d>].bench_function("faer_det_from_lu", |bencher| {
                    bencher.iter(|| {
                        let det = faer_det_from_partial_piv_lu(&fa_lu);
                        black_box(det);
                    });
                });

                // === Determinant from a precomputed SPD factorization ===
                [<group_d $d>].bench_function("la_stack_det_from_ldlt", |bencher| {
                    bencher.iter(|| {
                        let det = require_ok(a_ldlt.det(), "precomputed la_stack LDLT determinant");
                        black_box(det);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_det_from_cholesky", |bencher| {
                    bencher.iter(|| {
                        let det = na_cholesky.determinant();
                        black_box(det);
                    });
                });

                [<group_d $d>].bench_function("faer_det_from_ldlt", |bencher| {
                    bencher.iter(|| {
                        let det = faer_det_from_ldlt(&fa_ldlt);
                        black_box(det);
                    });
                });

                // === Vector dot product ===
                [<group_d $d>].bench_function("la_stack_dot", |bencher| {
                    bencher.iter(|| {
                        let result = require_ok(black_box(v1).dot(black_box(v2)), "la_stack dot");
                        black_box(result);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_dot", |bencher| {
                    bencher.iter(|| {
                        let result = black_box(&nv1).dot(black_box(&nv2));
                        black_box(result);
                    });
                });

                [<group_d $d>].bench_function("faer_dot", |bencher| {
                    bencher.iter(|| {
                        let mut sum = 0.0;
                        let a = black_box(&fv1);
                        let b = black_box(&fv2);
                        for i in 0..$d {
                            sum = a[(i, 0)].mul_add(b[(i, 0)], sum);
                        }
                        black_box(sum);
                    });
                });

                // === Vector norm squared ===
                [<group_d $d>].bench_function("la_stack_norm2_sq", |bencher| {
                    bencher.iter(|| {
                        let result = require_ok(black_box(v1).norm2_sq(), "la_stack norm2_sq");
                        black_box(result);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_norm_squared", |bencher| {
                    bencher.iter(|| {
                        let result = black_box(&nv1).norm_squared();
                        black_box(result);
                    });
                });

                [<group_d $d>].bench_function("faer_norm2_sq", |bencher| {
                    bencher.iter(|| {
                        let mut sum = 0.0;
                        let v = black_box(&fv1);
                        for i in 0..$d {
                            let x = v[(i, 0)];
                            sum = x.mul_add(x, sum);
                        }
                        black_box(sum);
                    });
                });

                // === Matrix infinity norm (max absolute row sum) ===
                [<group_d $d>].bench_function("la_stack_inf_norm", |bencher| {
                    bencher.iter(|| {
                        let result = require_ok(black_box(a).inf_norm(), "la_stack inf_norm");
                        black_box(result);
                    });
                });

                [<group_d $d>].bench_function("nalgebra_inf_norm", |bencher| {
                    bencher.iter(|| {
                        let result = nalgebra_inf_norm::<$d>(black_box(&na));
                        black_box(result);
                    });
                });

                [<group_d $d>].bench_function("faer_inf_norm", |bencher| {
                    bencher.iter(|| {
                        let m = black_box(&fa);
                        let mut max_row_sum = 0.0;

                        for r in 0..$d {
                            let mut row_sum = 0.0;
                            for c in 0..$d {
                                row_sum += m[(r, c)].abs();
                            }
                            if row_sum > max_row_sum {
                                max_row_sum = row_sum;
                            }
                        }

                        black_box(max_row_sum);
                    });
                });

                [<group_d $d>].finish();
            }
        }
    }
    };
}

define_vs_linalg_benches_for_dim!(bench_d2, 2);
define_vs_linalg_benches_for_dim!(bench_d3, 3);
define_vs_linalg_benches_for_dim!(bench_d4, 4);
define_vs_linalg_benches_for_dim!(bench_d5, 5);
define_vs_linalg_benches_for_dim!(bench_d8, 8);
define_vs_linalg_benches_for_dim!(bench_d16, 16);
define_vs_linalg_benches_for_dim!(bench_d32, 32);
define_vs_linalg_benches_for_dim!(bench_d64, 64);

fn main() {
    let mut c = Criterion::default().configure_from_args();

    bench_d2(&mut c);
    bench_d3(&mut c);
    bench_d4(&mut c);
    bench_d5(&mut c);

    bench_d8(&mut c);
    bench_d16(&mut c);
    bench_d32(&mut c);
    bench_d64(&mut c);

    c.final_summary();
}
