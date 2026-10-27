#![forbid(unsafe_code)]

//! Public construction, access, mutation and symmetry contracts outside kernels.

use std::hint::black_box;

use criterion::measurement::WallTime;
use criterion::{BenchmarkGroup, Criterion};
use la_stack::{Interval, IntervalMatrix, Matrix, Tolerance, Vector};

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

fn register<const D: usize>(c: &mut Criterion) {
    let matrix = Matrix::<D>::identity();
    let rows = *matrix.as_rows();
    let vector = Vector::<D>::try_new([1.0; D]).or_abort("public vector fixture");
    let tolerance = Tolerance::try_new(0.0).or_abort("public tolerance fixture");
    assert_eq!(matrix.det().or_abort("identity determinant"), 1.0);
    assert!(matrix.is_symmetric(tolerance).or_abort("identity symmetry"));
    assert_eq!(
        matrix
            .first_asymmetry(tolerance)
            .or_abort("identity asymmetry"),
        None
    );

    let mut group = c.benchmark_group(format!("public_d{D}"));
    macro_rules! bench {
        ($name:literal, $expression:expr) => {
            group.bench_function($name, |b| b.iter(|| black_box($expression)));
        };
    }
    bench!(
        "matrix_from_rows",
        Matrix::try_from_rows(black_box(rows)).or_abort("matrix construction")
    );
    bench!("matrix_zero", Matrix::<D>::zero());
    bench!("matrix_identity", Matrix::<D>::identity());
    bench!(
        "det_direct",
        black_box(&matrix)
            .det_direct()
            .or_abort("direct determinant")
    );
    bench!(
        "det_bound",
        black_box(&matrix)
            .det_errbound()
            .or_abort("determinant bound")
    );
    bench!("det_with_bound", {
        black_box(&matrix)
            .det_direct_with_errbound()
            .or_abort("bounded determinant")
            .map(|value| (value.determinant(), value.absolute_error_bound()))
    });
    bench!(
        "matrix_access",
        (
            black_box(&matrix).get(0, 0),
            black_box(&matrix).try_get(0, 0),
            black_box(&matrix).as_rows(),
            black_box(matrix).into_rows()
        )
    );
    bench!("matrix_set", {
        let mut value = black_box(matrix);
        value.set(0, 0, black_box(2.0)).or_abort("matrix mutation");
        value
    });
    bench!(
        "is_symmetric",
        black_box(&matrix)
            .is_symmetric(black_box(tolerance))
            .or_abort("symmetry")
    );
    bench!(
        "first_asymmetry",
        black_box(&matrix)
            .first_asymmetry(black_box(tolerance))
            .or_abort("asymmetry")
    );
    bench!(
        "vector_new",
        Vector::try_new(black_box([1.0; D])).or_abort("vector construction")
    );
    bench!("vector_zero", Vector::<D>::zero());
    bench!(
        "vector_access",
        (
            black_box(&vector).as_array(),
            black_box(vector).into_array()
        )
    );
    register_intervals::<D>(&mut group, &matrix);
    group.finish();
}

fn register_intervals<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    matrix: &Matrix<D>,
) {
    let rows = *matrix.as_rows();
    let interval = Interval::point(1.0).or_abort("public interval fixture");
    let intervals = IntervalMatrix::from_matrix(matrix);
    assert!(
        intervals
            .det()
            .or_abort("identity interval determinant")
            .contains(1.0)
    );
    macro_rules! bench {
        ($name:literal, $expression:expr) => {
            group.bench_function($name, |b| b.iter(|| black_box($expression)));
        };
    }
    bench!(
        "interval_matrix_from_rows",
        IntervalMatrix::from_rows(black_box(*intervals.as_rows()))
    );
    bench!(
        "interval_matrix_from_points",
        IntervalMatrix::try_from_point_rows(black_box(rows))
            .or_abort("interval matrix construction")
    );
    bench!(
        "interval_matrix_from_matrix",
        IntervalMatrix::from_matrix(black_box(matrix))
    );
    bench!("interval_matrix_zero", IntervalMatrix::<D>::zero());
    bench!("interval_matrix_identity", IntervalMatrix::<D>::identity());
    bench!(
        "interval_matrix_access",
        (
            black_box(&intervals).get(0, 0),
            black_box(&intervals).try_get(0, 0),
            black_box(&intervals).as_rows(),
            black_box(intervals).into_rows()
        )
    );
    bench!("interval_matrix_set", {
        let mut value = black_box(intervals);
        value
            .set(0, 0, black_box(interval))
            .or_abort("interval matrix mutation");
        value
    });
    bench!(
        "interval_det",
        black_box(&intervals).det().or_abort("interval determinant")
    );
}

fn register_scalars(c: &mut Criterion) {
    let interval = Interval::point(1.0).or_abort("public scalar fixture");
    let mut group = c.benchmark_group("public_scalar");
    macro_rules! bench {
        ($name:literal, $expression:expr) => {
            group.bench_function($name, |b| b.iter(|| black_box($expression)));
        };
    }
    bench!(
        "tolerance",
        Tolerance::try_new(black_box(0.0))
            .or_abort("tolerance parsing")
            .get()
    );
    bench!(
        "interval_new",
        Interval::try_new(black_box(-1.0), black_box(1.0)).or_abort("interval construction")
    );
    bench!(
        "interval_point",
        Interval::point(black_box(1.0)).or_abort("interval point")
    );
    bench!(
        "interval_subtraction",
        Interval::try_from_subtraction(black_box(1.0), black_box(0.5))
            .or_abort("interval subtraction")
    );
    bench!(
        "interval_access",
        (
            black_box(interval).lower(),
            black_box(interval).upper(),
            black_box(interval).contains(black_box(1.0))
        )
    );
    bench!("interval_negate", black_box(&interval).negate());
    group.finish();
}

#[cfg(feature = "exact")]
fn register_exact<const D: usize>(c: &mut Criterion) {
    use la_stack::{BigRational, RationalMatrix, RationalVector};
    let entry = BigRational::from_integer(1.into());
    let matrix = RationalMatrix::<D>::try_from_fn(|i, j| {
        BigRational::from_integer(i32::from(i == j).into())
    })
    .or_abort("rational matrix fixture");
    let vector =
        RationalVector::<D>::try_from_fn(|_| entry.clone()).or_abort("rational vector fixture");
    assert_eq!(matrix.det(), entry);
    assert_eq!(
        matrix.solve(&vector).or_abort("rational identity solve"),
        vector
    );
    let mut group = c.benchmark_group(format!("public_exact_d{D}"));
    group.bench_function("matrix_from_rows", |b| {
        b.iter(|| {
            RationalMatrix::try_from_rows(black_box(matrix.as_rows().clone()))
                .or_abort("rational rows")
        });
    });
    group.bench_function("matrix_from_fn", |b| {
        b.iter(|| {
            RationalMatrix::<D>::try_from_fn(|i, j| {
                black_box(&matrix)
                    .get(i, j)
                    .or_abort("rational entry")
                    .clone()
            })
            .or_abort("rational function")
        });
    });
    group.bench_function("matrix_zero", |b| b.iter(RationalMatrix::<D>::zero));
    group.bench_function("matrix_access", |b| {
        b.iter(|| black_box((black_box(&matrix).as_rows(), black_box(&matrix).get(0, 0))));
    });
    group.bench_function("matrix_into_rows", |b| {
        b.iter(|| black_box(matrix.clone()).into_rows());
    });
    group.bench_function("vector_new", |b| {
        b.iter(|| {
            RationalVector::try_new(black_box(vector.as_array().clone()))
                .or_abort("rational vector")
        });
    });
    group.bench_function("vector_from_fn", |b| {
        b.iter(|| {
            RationalVector::<D>::try_from_fn(|i| {
                black_box(&vector)
                    .get(i)
                    .or_abort("rational vector entry")
                    .clone()
            })
            .or_abort("rational vector function")
        });
    });
    group.bench_function("vector_zero", |b| b.iter(RationalVector::<D>::zero));
    group.bench_function("vector_access", |b| {
        b.iter(|| black_box((black_box(&vector).as_array(), black_box(&vector).get(0))));
    });
    group.bench_function("vector_into_array", |b| {
        b.iter(|| black_box(vector.clone()).into_array());
    });
    group.finish();
}

#[cfg(feature = "exact")]
fn register_scalar(c: &mut Criterion) {
    use la_stack::{BigRational, ExactF64Conversion};
    let scalar = BigRational::new(1.into(), 2.into());
    assert_eq!(scalar.try_to_f64().or_abort("exact scalar fixture"), 0.5);
    assert_eq!(
        scalar.to_rounded_f64().or_abort("rounded scalar fixture"),
        0.5
    );
    let mut group = c.benchmark_group("public_exact_scalar");
    group.bench_function("strict", |b| {
        b.iter(|| black_box(&scalar).try_to_f64().or_abort("strict scalar"));
    });
    group.bench_function("rounded", |b| {
        b.iter(|| {
            black_box(&scalar)
                .to_rounded_f64()
                .or_abort("rounded scalar")
        });
    });
    group.finish();
}

fn main() {
    let mut c = Criterion::default().configure_from_args();
    register_scalars(&mut c);
    register::<2>(&mut c);
    register::<3>(&mut c);
    register::<4>(&mut c);
    register::<5>(&mut c);
    #[cfg(feature = "exact")]
    {
        register_scalar(&mut c);
        register_exact::<2>(&mut c);
        register_exact::<3>(&mut c);
        register_exact::<4>(&mut c);
        register_exact::<5>(&mut c);
    }
    c.final_summary();
}
