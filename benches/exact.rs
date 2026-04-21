//! Benchmarks for exact arithmetic operations.
//!
//! These benchmarks measure the performance of the `exact` feature's
//! arbitrary-precision methods.  They are organised into two classes:
//!
//! 1. **General-case benches** (`exact_d{2..5}`) — a single
//!    well-conditioned diagonally-dominant matrix per dimension.  These
//!    measure typical-case performance and track regressions against a
//!    reproducible input.
//! 2. **Adversarial / extreme-input benches** — matrices chosen to
//!    stress specific corners of the exact-arithmetic pipeline:
//!    near-singularity (forces the Bareiss fallback), large f64 entries
//!    (stresses intermediate `BigInt` growth), and Hilbert-style
//!    ill-conditioning (wide range of `(mantissa, exponent)` pairs in
//!    the `f64_decompose → BigInt` path).  These measure tail behaviour
//!    that fixed well-conditioned inputs miss and provide stronger
//!    empirical evidence for `docs/PERFORMANCE.md`.

use criterion::{BenchmarkGroup, Criterion, measurement::WallTime};
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
///
/// The base `[[1,2,3],[4,5,6],[7,8,9]]` is exactly singular; adding
/// `2^-50` to entry (0,0) makes `det = -3 × 2^-50 ≠ 0`.  The f64 filter
/// in `det_sign_exact` cannot resolve this sign, so Bareiss is forced;
/// `solve_exact` is the primary use case for near-degenerate inputs
/// (exact circumcenter etc.) and exercises the largest intermediate
/// `BigInt` values in the hybrid solve.
#[inline]
fn near_singular_3x3() -> Matrix<3> {
    let perturbation = f64::from_bits(0x3CD0_0000_0000_0000); // 2^-50
    Matrix::<3>::from_rows([
        [1.0 + perturbation, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ])
}

/// Large-entry 3×3: strictly diagonally-dominant matrix with diagonal
/// entries near `f64::MAX / 2` and ones elsewhere.
///
/// Each big entry decomposes into a 53-bit mantissa with exponent `~970`;
/// the unit off-diagonals have exponent `0`, so the shared `e_min = 0`
/// shift in `component_to_bigint` produces `BigInt`s of `~1023` bits for
/// the diagonal and small integers elsewhere.  Bareiss fraction-free
/// updates then multiply these together, stressing the big-integer
/// multiply and allocator along the full `O(D³)` elimination phase.  The
/// matrix is non-singular (det ≈ `big³`) so both `det_*` and `solve_*`
/// paths complete.
#[inline]
fn large_entries_3x3() -> Matrix<3> {
    let big = f64::MAX / 2.0;
    Matrix::<3>::from_rows([[big, 1.0, 1.0], [1.0, big, 1.0], [1.0, 1.0, big]])
}

/// Hilbert matrix `H[i][j] = 1 / (i + j + 1)`.
///
/// Most entries (`1/3`, `1/5`, `1/6`, `1/7`, …) are non-terminating in
/// binary, so every cell has a distinct 53-bit mantissa and a small
/// negative exponent.  `f64_decompose` therefore produces a wide mix of
/// `(mantissa, exponent)` pairs with no shared power-of-two factors,
/// and the scaling shift to the common `e_min` yields `BigInt` values
/// of varied bit-lengths — a different kind of adversarial input from
/// the large-entries case.  Hilbert matrices are also classically
/// ill-conditioned (condition number grows exponentially with D), so
/// they are a realistic stand-in for the near-degenerate geometric
/// predicate inputs that motivate exact arithmetic.
#[inline]
#[allow(clippy::cast_precision_loss)]
fn hilbert<const D: usize>() -> Matrix<D> {
    let mut rows = [[0.0; D]; D];
    let mut r = 0;
    while r < D {
        let mut c = 0;
        while c < D {
            rows[r][c] = 1.0 / ((r + c + 1) as f64);
            c += 1;
        }
        r += 1;
    }
    Matrix::<D>::from_rows(rows)
}

/// Populate a Criterion group with the four headline exact-arithmetic
/// benches on a single `(matrix, rhs)` pair: `det_sign_exact`,
/// `det_exact`, `solve_exact`, `solve_exact_f64`.
///
/// Used by every adversarial-input group so each one measures the same
/// operations, making the resulting tables directly comparable.
fn bench_extreme_group<const D: usize>(
    group: &mut BenchmarkGroup<'_, WallTime>,
    m: Matrix<D>,
    rhs: Vector<D>,
) {
    group.bench_function("det_sign_exact", |bencher| {
        bencher.iter(|| {
            let sign = black_box(m)
                .det_sign_exact()
                .expect("finite matrix entries");
            black_box(sign);
        });
    });

    group.bench_function("det_exact", |bencher| {
        bencher.iter(|| {
            let det = black_box(m).det_exact().expect("finite matrix entries");
            black_box(det);
        });
    });

    group.bench_function("solve_exact", |bencher| {
        bencher.iter(|| {
            let x = black_box(m)
                .solve_exact(black_box(rhs))
                .expect("non-singular matrix with finite entries");
            let _ = black_box(x);
        });
    });

    group.bench_function("solve_exact_f64", |bencher| {
        bencher.iter(|| {
            let x = black_box(m)
                .solve_exact_f64(black_box(rhs))
                .expect("solution representable in f64");
            let _ = black_box(x);
        });
    });
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
                        .expect("diagonally dominant matrix is non-singular");
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
                    let det = black_box(a).det_exact().expect("finite matrix entries");
                    black_box(det);
                });
            });

            // === det_exact_f64 (exact → f64) ===
            [<group_d $d>].bench_function("det_exact_f64", |bencher| {
                bencher.iter(|| {
                    let det = black_box(a)
                        .det_exact_f64()
                        .expect("det representable in f64");
                    black_box(det);
                });
            });

            // === det_sign_exact (adaptive: fast filter + exact fallback) ===
            [<group_d $d>].bench_function("det_sign_exact", |bencher| {
                bencher.iter(|| {
                    let sign = black_box(a)
                        .det_sign_exact()
                        .expect("finite matrix entries");
                    black_box(sign);
                });
            });

            // === solve_exact (BigRational result) ===
            [<group_d $d>].bench_function("solve_exact", |bencher| {
                bencher.iter(|| {
                    let x = black_box(a)
                        .solve_exact(black_box(rhs))
                        .expect("diagonally dominant matrix is non-singular");
                    black_box(x);
                });
            });

            // === solve_exact_f64 (exact → f64) ===
            [<group_d $d>].bench_function("solve_exact_f64", |bencher| {
                bencher.iter(|| {
                    let x = black_box(a)
                        .solve_exact_f64(black_box(rhs))
                        .expect("solution representable in f64");
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

    // === Adversarial / extreme-input groups ===
    //
    // Each group runs the same four exact-arithmetic benches
    // (`det_sign_exact`, `det_exact`, `solve_exact`, `solve_exact_f64`)
    // via `bench_extreme_group`, so the resulting tables are directly
    // comparable across input classes.

    // Near-singular 3×3: forces Bareiss fallback in det_sign_exact and
    // exercises the largest intermediate BigInt values in solve_exact
    // (the primary motivating use case for exact solve).
    {
        let mut group = c.benchmark_group("exact_near_singular_3x3");
        bench_extreme_group(
            &mut group,
            near_singular_3x3(),
            Vector::<3>::new([1.0, 2.0, 3.0]),
        );
        group.finish();
    }

    // Large-entry 3×3: diagonal entries near `f64::MAX / 2` stress
    // BigInt growth during Bareiss forward elimination.
    {
        let mut group = c.benchmark_group("exact_large_entries_3x3");
        bench_extreme_group(
            &mut group,
            large_entries_3x3(),
            Vector::<3>::new([1.0, 1.0, 1.0]),
        );
        group.finish();
    }

    // Hilbert 4×4 and 5×5: classically ill-conditioned matrices whose
    // entries span many orders of magnitude in `(mantissa, exponent)`
    // space, exercising the f64 → BigInt scaling path.
    {
        let mut group = c.benchmark_group("exact_hilbert_4x4");
        bench_extreme_group(&mut group, hilbert::<4>(), Vector::<4>::new([1.0; 4]));
        group.finish();
    }

    {
        let mut group = c.benchmark_group("exact_hilbert_5x5");
        bench_extreme_group(&mut group, hilbert::<5>(), Vector::<5>::new([1.0; 5]));
        group.finish();
    }

    c.final_summary();
}
