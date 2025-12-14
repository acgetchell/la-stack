//! Minimal benchmark harness.
//!
//! This exists so the `[[bench]]` target in `Cargo.toml` is valid and `cargo test`
//! can parse the manifest. It also provides a small comparison point against
//! `nalgebra` for determinant computation.

use criterion::Criterion;
use std::hint::black_box;

fn main() {
    let mut c = Criterion::default().configure_from_args();

    // Small, fixed-size system representative of the geometric use cases.
    let a = la_stack::Matrix::<3>::from_rows([[1.0, 2.0, 3.0], [0.5, -4.0, 0.25], [7.0, 0.0, 1.0]]);

    let na = nalgebra::Matrix3::new(
        1.0, 2.0, 3.0, //
        0.5, -4.0, 0.25, //
        7.0, 0.0, 1.0,
    );

    c.bench_function("la_stack_det_3x3", |b| {
        b.iter(|| {
            let det = black_box(a)
                .det(la_stack::DEFAULT_PIVOT_TOL)
                .expect("matrix should be non-singular");
            black_box(det);
        });
    });

    c.bench_function("nalgebra_det_3x3", |b| {
        b.iter(|| {
            let det = black_box(na).determinant();
            black_box(det);
        });
    });

    c.final_summary();
}
