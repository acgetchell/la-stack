#![forbid(unsafe_code)]

//! Benchmark comparison between la-stack and other Rust linear algebra crates.
//!
//! Goal: like-for-like comparisons of the operations la-stack supports across several
//! fixed dimensions.
//!
//! Notes:
//! - Determinant groups distinguish factorization-inclusive LU, `Matrix::det`,
//!   precomputed LU determinant queries, and precomputed LDLT/Cholesky queries.
//! - Matrix infinity norm is the maximum absolute row sum on all sides.

use std::hint::black_box;

use criterion::measurement::WallTime;
use criterion::{BatchSize, BenchmarkGroup, Criterion};
use faer::linalg::solvers::Solve;
use faer::mat::AsMatRef;
use faer::{Mat, Side};
use nalgebra::{Const, DimMin, SMatrix, SVector};

use la_stack::{DEFAULT_SINGULAR_TOL, Matrix, Vector};

#[path = "common/bench_utils.rs"]
mod bench_utils;
#[path = "common/vs_linalg.rs"]
pub mod vs_linalg_common;

use bench_utils::OrAbort;
use vs_linalg_common::{
    PreparedFaerLuDet, faer_det_from_ldlt, la_stack_dot, la_stack_tolerance,
    make_balanced_dynamic_range_rows, make_ill_conditioned_matrix_rows, make_matrix_rows,
    make_pivoting_matrix_rows, make_vector_array, matrix_entry, nalgebra_inf_norm, vector_entry,
};

/// Build the deterministic la-stack matrix shared by a benchmark family.
fn la_matrix<const D: usize>() -> Matrix<D> {
    Matrix::try_from_rows(make_matrix_rows()).or_abort("la_stack matrix construction")
}

/// Build a deterministic la-stack vector with the requested offset.
fn la_vector<const D: usize>(offset: f64, operation: &str) -> Vector<D> {
    Vector::try_new(make_vector_array(offset)).or_abort(operation)
}

/// Build the deterministic nalgebra matrix shared by a benchmark family.
fn nalgebra_matrix<const D: usize>() -> SMatrix<f64, D, D> {
    SMatrix::from_fn(matrix_entry::<D>)
}

/// Build a deterministic nalgebra vector with the requested offset.
fn nalgebra_vector<const D: usize>(offset: f64) -> SVector<f64, D> {
    SVector::from_fn(|i, _| vector_entry(i, offset))
}

/// Build the deterministic faer matrix shared by a benchmark family.
fn faer_matrix<const D: usize>() -> Mat<f64> {
    Mat::from_fn(D, D, matrix_entry::<D>)
}

/// Build a deterministic faer column vector with the requested offset.
fn faer_vector<const D: usize>(offset: f64) -> Mat<f64> {
    Mat::from_fn(D, 1, |i, _| vector_entry(i, offset))
}

/// Register determinant benchmarks that include factorization work.
fn register_determinant_benchmarks<const D: usize>(group: &mut BenchmarkGroup<'_, WallTime>)
where
    Const<D>: DimMin<Const<D>, Output = Const<D>>,
{
    let a = la_matrix::<D>();
    let na = nalgebra_matrix::<D>();
    let fa = faer_matrix::<D>();

    group.bench_function("la_stack_det_via_lu", |bencher| {
        bencher.iter_batched(
            || a,
            |a| {
                let lu = black_box(a)
                    .lu(DEFAULT_SINGULAR_TOL)
                    .or_abort("la_stack LU factorization");
                let det = lu.det().or_abort("la_stack LU determinant");
                black_box(det);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("nalgebra_det_via_lu", |bencher| {
        bencher.iter_batched(
            || na,
            |na| {
                let lu = black_box(na).lu();
                let det = lu.determinant();
                black_box(det);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("faer_det_via_lu", |bencher| {
        bencher.iter_batched(
            || &fa,
            |fa| {
                let lu = black_box(fa).partial_piv_lu();
                let det = PreparedFaerLuDet::new(&lu).det();
                black_box(det);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("la_stack_det", |bencher| {
        bencher.iter_batched(
            || a,
            |a| {
                let det = black_box(a).det().or_abort("la_stack determinant");
                black_box(det);
            },
            BatchSize::SmallInput,
        );
    });
}

/// Register LU, LDLT, and Cholesky factorization benchmarks.
fn register_factorization_benchmarks<const D: usize>(group: &mut BenchmarkGroup<'_, WallTime>)
where
    Const<D>: DimMin<Const<D>, Output = Const<D>>,
{
    let a = la_matrix::<D>();
    let na = nalgebra_matrix::<D>();
    let fa = faer_matrix::<D>();

    group.bench_function("la_stack_lu", |bencher| {
        bencher.iter_batched(
            || a,
            |a| {
                let lu = black_box(a)
                    .lu(DEFAULT_SINGULAR_TOL)
                    .or_abort("la_stack LU factorization");
                let _ = black_box(lu);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("nalgebra_lu", |bencher| {
        bencher.iter_batched(
            || na,
            |na| {
                let lu = black_box(na).lu();
                let _ = black_box(lu);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("faer_lu", |bencher| {
        bencher.iter_batched(
            || &fa,
            |fa| {
                let lu = black_box(fa).partial_piv_lu();
                let _ = black_box(lu);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("la_stack_ldlt", |bencher| {
        bencher.iter_batched(
            || a,
            |a| {
                let ldlt = black_box(a)
                    .ldlt(DEFAULT_SINGULAR_TOL)
                    .or_abort("la_stack LDLT factorization");
                let _ = black_box(ldlt);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("nalgebra_cholesky", |bencher| {
        bencher.iter_batched(
            || na,
            |na| {
                let chol = black_box(na)
                    .cholesky()
                    .or_abort("nalgebra Cholesky factorization");
                black_box(chol);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("faer_ldlt", |bencher| {
        bencher.iter_batched(
            || &fa,
            |fa| {
                let ldlt = black_box(fa)
                    .ldlt(Side::Lower)
                    .or_abort("faer LDLT factorization");
                let _ = black_box(ldlt);
            },
            BatchSize::SmallInput,
        );
    });
}

/// Register solves that include LU factorization work.
fn register_lu_solve_benchmarks<const D: usize>(group: &mut BenchmarkGroup<'_, WallTime>)
where
    Const<D>: DimMin<Const<D>, Output = Const<D>>,
{
    let a = la_matrix::<D>();
    let rhs = la_vector::<D>(0.0, "la_stack RHS vector construction");
    let na = nalgebra_matrix::<D>();
    let nrhs = nalgebra_vector::<D>(0.0);
    let fa = faer_matrix::<D>();
    let frhs = faer_vector::<D>(0.0);

    group.bench_function("la_stack_lu_solve", |bencher| {
        bencher.iter(|| {
            let lu = black_box(a)
                .lu(DEFAULT_SINGULAR_TOL)
                .or_abort("la_stack LU factorization");
            let x = lu.solve(black_box(rhs)).or_abort("la_stack LU solve");
            let _ = black_box(x);
        });
    });

    group.bench_function("nalgebra_lu_solve", |bencher| {
        bencher.iter(|| {
            let lu = black_box(na).lu();
            let x = lu.solve(black_box(&nrhs)).or_abort("nalgebra LU solve");
            black_box(x);
        });
    });

    group.bench_function("faer_lu_solve", |bencher| {
        bencher.iter(|| {
            let lu = black_box(&fa).partial_piv_lu();
            let x = lu.solve(black_box(&frhs));
            black_box(x);
        });
    });
}

/// Register solves that include LDLT or Cholesky factorization work.
fn register_ldlt_solve_benchmarks<const D: usize>(group: &mut BenchmarkGroup<'_, WallTime>) {
    let a = la_matrix::<D>();
    let rhs = la_vector::<D>(0.0, "la_stack RHS vector construction");
    let na = nalgebra_matrix::<D>();
    let nrhs = nalgebra_vector::<D>(0.0);
    let fa = faer_matrix::<D>();
    let frhs = faer_vector::<D>(0.0);

    group.bench_function("la_stack_ldlt_solve", |bencher| {
        bencher.iter_batched(
            || (a, rhs),
            |(a, rhs)| {
                let ldlt = black_box(a)
                    .ldlt(DEFAULT_SINGULAR_TOL)
                    .or_abort("la_stack LDLT factorization");
                let x = ldlt.solve(black_box(rhs)).or_abort("la_stack LDLT solve");
                let _ = black_box(x);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("nalgebra_cholesky_solve", |bencher| {
        bencher.iter_batched(
            || (na, nrhs),
            |(na, nrhs)| {
                let chol = black_box(na)
                    .cholesky()
                    .or_abort("nalgebra Cholesky factorization");
                let x = chol.solve(black_box(&nrhs));
                black_box(x);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("faer_ldlt_solve", |bencher| {
        bencher.iter_batched(
            || (&fa, &frhs),
            |(fa, rhs)| {
                let ldlt = black_box(fa)
                    .ldlt(Side::Lower)
                    .or_abort("faer LDLT factorization");
                let x = ldlt.solve(black_box(rhs));
                black_box(x);
            },
            BatchSize::SmallInput,
        );
    });
}

/// Register solves using precomputed LU factorizations.
fn register_precomputed_lu_solve_benchmarks<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
) where
    Const<D>: DimMin<Const<D>, Output = Const<D>>,
{
    let a = la_matrix::<D>();
    let rhs = la_vector::<D>(0.0, "la_stack RHS vector construction");
    let na = nalgebra_matrix::<D>();
    let nrhs = nalgebra_vector::<D>(0.0);
    let fa = faer_matrix::<D>();
    let frhs = faer_vector::<D>(0.0);
    let a_lu = a
        .lu(DEFAULT_SINGULAR_TOL)
        .or_abort("precomputed la_stack LU");
    let na_lu = na.lu();
    let fa_lu = fa.partial_piv_lu();

    group.bench_function("la_stack_solve_from_lu", |bencher| {
        bencher.iter(|| {
            let x = black_box(&a_lu)
                .solve(black_box(rhs))
                .or_abort("precomputed la_stack LU solve");
            let _ = black_box(x);
        });
    });

    group.bench_function("nalgebra_solve_from_lu", |bencher| {
        bencher.iter(|| {
            let x = black_box(&na_lu)
                .solve(black_box(&nrhs))
                .or_abort("precomputed nalgebra LU solve");
            black_box(x);
        });
    });

    group.bench_function("faer_solve_from_lu", |bencher| {
        bencher.iter(|| {
            let x = black_box(&fa_lu).solve(black_box(&frhs));
            black_box(x);
        });
    });
}

/// Register solves using precomputed LDLT or Cholesky factorizations.
fn register_precomputed_ldlt_solve_benchmarks<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
) {
    let a = la_matrix::<D>();
    let rhs = la_vector::<D>(0.0, "la_stack RHS vector construction");
    let na = nalgebra_matrix::<D>();
    let nrhs = nalgebra_vector::<D>(0.0);
    let fa = faer_matrix::<D>();
    let frhs = faer_vector::<D>(0.0);
    let a_ldlt = a
        .ldlt(DEFAULT_SINGULAR_TOL)
        .or_abort("precomputed la_stack LDLT");
    let na_cholesky = na.cholesky().or_abort("precomputed nalgebra Cholesky");
    let fa_ldlt = fa.ldlt(Side::Lower).or_abort("precomputed faer LDLT");

    group.bench_function("la_stack_solve_from_ldlt", |bencher| {
        bencher.iter(|| {
            let x = black_box(&a_ldlt)
                .solve(black_box(rhs))
                .or_abort("precomputed la_stack LDLT solve");
            let _ = black_box(x);
        });
    });

    group.bench_function("nalgebra_solve_from_cholesky", |bencher| {
        bencher.iter(|| {
            let x = black_box(&na_cholesky).solve(black_box(&nrhs));
            black_box(x);
        });
    });

    group.bench_function("faer_solve_from_ldlt", |bencher| {
        bencher.iter(|| {
            let x = black_box(&fa_ldlt).solve(black_box(&frhs));
            black_box(x);
        });
    });
}

/// Register determinant queries using precomputed LU factorizations.
fn register_precomputed_lu_determinant_benchmarks<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
) where
    Const<D>: DimMin<Const<D>, Output = Const<D>>,
{
    let a = la_matrix::<D>();
    let na = nalgebra_matrix::<D>();
    let fa = faer_matrix::<D>();
    let a_lu = a
        .lu(DEFAULT_SINGULAR_TOL)
        .or_abort("precomputed la_stack LU");
    let na_lu = na.lu();
    let fa_lu = fa.partial_piv_lu();
    let fa_lu_det = PreparedFaerLuDet::new(&fa_lu);

    group.bench_function("la_stack_det_from_lu", |bencher| {
        bencher.iter(|| {
            let det = black_box(&a_lu)
                .det()
                .or_abort("precomputed la_stack LU determinant");
            black_box(det);
        });
    });

    group.bench_function("nalgebra_det_from_lu", |bencher| {
        bencher.iter(|| {
            let det = black_box(&na_lu).determinant();
            black_box(det);
        });
    });

    group.bench_function("faer_det_from_lu", |bencher| {
        bencher.iter(|| {
            let det = black_box(&fa_lu_det).det();
            black_box(det);
        });
    });
}

/// Register determinant queries using precomputed LDLT or Cholesky factorizations.
fn register_precomputed_ldlt_determinant_benchmarks<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
) {
    let a = la_matrix::<D>();
    let na = nalgebra_matrix::<D>();
    let fa = faer_matrix::<D>();
    let a_ldlt = a
        .ldlt(DEFAULT_SINGULAR_TOL)
        .or_abort("precomputed la_stack LDLT");
    let na_cholesky = na.cholesky().or_abort("precomputed nalgebra Cholesky");
    let fa_ldlt = fa.ldlt(Side::Lower).or_abort("precomputed faer LDLT");

    group.bench_function("la_stack_det_from_ldlt", |bencher| {
        bencher.iter(|| {
            let det = black_box(&a_ldlt)
                .det()
                .or_abort("precomputed la_stack LDLT determinant");
            black_box(det);
        });
    });

    group.bench_function("nalgebra_det_from_cholesky", |bencher| {
        bencher.iter(|| {
            let det = black_box(&na_cholesky).determinant();
            black_box(det);
        });
    });

    group.bench_function("faer_det_from_ldlt", |bencher| {
        bencher.iter(|| {
            let det = faer_det_from_ldlt(black_box(&fa_ldlt));
            black_box(det);
        });
    });
}

/// Register vector dot-product and squared-norm benchmarks.
fn register_vector_benchmarks<const D: usize>(group: &mut BenchmarkGroup<'_, WallTime>) {
    let v1 = la_vector::<D>(0.0, "la_stack vector construction");
    let v2 = la_vector::<D>(1.0, "la_stack vector construction");
    let nv1 = nalgebra_vector::<D>(0.0);
    let nv2 = nalgebra_vector::<D>(1.0);
    let fv1 = faer_vector::<D>(0.0);
    let fv2 = faer_vector::<D>(1.0);

    group.bench_function("la_stack_dot", |bencher| {
        bencher.iter(|| {
            let result = la_stack_dot(black_box(&v1), black_box(&v2)).or_abort("la_stack dot");
            black_box(result);
        });
    });

    group.bench_function("nalgebra_dot", |bencher| {
        bencher.iter(|| {
            let result = black_box(&nv1).dot(black_box(&nv2));
            black_box(result);
        });
    });

    group.bench_function("faer_dot", |bencher| {
        bencher.iter(|| {
            let mut sum = 0.0;
            let a = black_box(&fv1);
            let b = black_box(&fv2);
            for i in 0..D {
                sum = a[(i, 0)].mul_add(b[(i, 0)], sum);
            }
            black_box(sum);
        });
    });

    group.bench_function("la_stack_norm2_sq", |bencher| {
        bencher.iter(|| {
            let result = black_box(&v1).norm2_sq().or_abort("la_stack norm2_sq");
            black_box(result);
        });
    });

    group.bench_function("nalgebra_norm_squared", |bencher| {
        bencher.iter(|| {
            let result = black_box(&nv1).norm_squared();
            black_box(result);
        });
    });

    group.bench_function("faer_norm2_sq", |bencher| {
        bencher.iter(|| {
            let v = black_box(&fv1);
            let result = v.as_mat_ref().squared_norm_l2();
            black_box(result);
        });
    });
}

/// Register matrix infinity-norm benchmarks.
fn register_matrix_norm_benchmarks<const D: usize>(group: &mut BenchmarkGroup<'_, WallTime>) {
    let a = la_matrix::<D>();
    let na = nalgebra_matrix::<D>();
    let fa = faer_matrix::<D>();

    group.bench_function("la_stack_inf_norm", |bencher| {
        bencher.iter(|| {
            let result = black_box(&a).inf_norm().or_abort("la_stack inf_norm");
            black_box(result);
        });
    });

    group.bench_function("nalgebra_inf_norm", |bencher| {
        bencher.iter(|| {
            let result = nalgebra_inf_norm::<D>(black_box(&na));
            black_box(result);
        });
    });

    group.bench_function("faer_inf_norm", |bencher| {
        bencher.iter(|| {
            let m = black_box(&fa);
            let mut max_row_sum = 0.0;

            for r in 0..D {
                let mut row_sum = 0.0;
                for c in 0..D {
                    row_sum += m[(r, c)].abs();
                }
                if row_sum > max_row_sum {
                    max_row_sum = row_sum;
                }
            }

            black_box(max_row_sum);
        });
    });
}

/// Register D=8 stress cases that exercise pivoting, conditioning, and scaled products.
fn register_stress_benchmarks(group: &mut BenchmarkGroup<'_, WallTime>) {
    let zero_tolerance = la_stack_tolerance(0.0).or_abort("zero benchmark tolerance");
    let pivoting = Matrix::<8>::try_from_rows(make_pivoting_matrix_rows())
        .or_abort("pivoting benchmark matrix construction");
    let ill_conditioned = Matrix::<8>::try_from_rows(make_ill_conditioned_matrix_rows())
        .or_abort("ill-conditioned benchmark matrix construction");
    #[cfg(not(la_stack_v0_4_3_api))]
    let balanced = Matrix::<8>::try_from_rows(make_balanced_dynamic_range_rows())
        .or_abort("balanced-range benchmark matrix construction");

    group.bench_function("la_stack_lu_pivoting", |bencher| {
        bencher.iter_batched(
            || pivoting,
            |matrix| {
                let lu = black_box(matrix)
                    .lu(zero_tolerance)
                    .or_abort("pivoting LU factorization");
                let _ = black_box(lu);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("la_stack_lu_ill_conditioned", |bencher| {
        bencher.iter_batched(
            || ill_conditioned,
            |matrix| {
                let lu = black_box(matrix)
                    .lu(zero_tolerance)
                    .or_abort("ill-conditioned LU factorization");
                let _ = black_box(lu);
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("la_stack_ldlt_ill_conditioned", |bencher| {
        bencher.iter_batched(
            || ill_conditioned,
            |matrix| {
                let ldlt = black_box(matrix)
                    .ldlt(zero_tolerance)
                    .or_abort("ill-conditioned LDLT factorization");
                let _ = black_box(ldlt);
            },
            BatchSize::SmallInput,
        );
    });

    #[cfg(not(la_stack_v0_4_3_api))]
    let balanced_lu = balanced
        .lu(zero_tolerance)
        .or_abort("balanced-range LU factorization");
    #[cfg(not(la_stack_v0_4_3_api))]
    let balanced_ldlt = balanced
        .ldlt(zero_tolerance)
        .or_abort("balanced-range LDLT factorization");

    #[cfg(not(la_stack_v0_4_3_api))]
    group.bench_function("la_stack_det_from_lu_balanced_range", |bencher| {
        bencher.iter(|| {
            let det = black_box(&balanced_lu)
                .det()
                .or_abort("balanced-range LU determinant");
            black_box(det);
        });
    });

    #[cfg(not(la_stack_v0_4_3_api))]
    group.bench_function("la_stack_det_from_ldlt_balanced_range", |bencher| {
        bencher.iter(|| {
            let det = black_box(&balanced_ldlt)
                .det()
                .or_abort("balanced-range LDLT determinant");
            black_box(det);
        });
    });
}

macro_rules! define_vs_linalg_benches_for_dim {
    ($fn_name:ident, $d:literal $(, $register_stress:ident)?) => {
        fn $fn_name(c: &mut Criterion) {
            let mut group = c.benchmark_group(concat!("d", stringify!($d)));
            register_determinant_benchmarks::<$d>(&mut group);
            register_factorization_benchmarks::<$d>(&mut group);
            register_lu_solve_benchmarks::<$d>(&mut group);
            register_ldlt_solve_benchmarks::<$d>(&mut group);
            register_precomputed_lu_solve_benchmarks::<$d>(&mut group);
            register_precomputed_ldlt_solve_benchmarks::<$d>(&mut group);
            register_precomputed_lu_determinant_benchmarks::<$d>(&mut group);
            register_precomputed_ldlt_determinant_benchmarks::<$d>(&mut group);
            register_vector_benchmarks::<$d>(&mut group);
            register_matrix_norm_benchmarks::<$d>(&mut group);
            $(
                $register_stress(&mut group);
            )?
            group.finish();
        }
    };
}

define_vs_linalg_benches_for_dim!(bench_d2, 2);
define_vs_linalg_benches_for_dim!(bench_d3, 3);
define_vs_linalg_benches_for_dim!(bench_d4, 4);
define_vs_linalg_benches_for_dim!(bench_d5, 5);
define_vs_linalg_benches_for_dim!(bench_d8, 8, register_stress_benchmarks);
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
