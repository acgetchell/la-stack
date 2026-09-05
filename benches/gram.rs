#![forbid(unsafe_code)]

//! Gram construction versus checked hand-written assembly, excluding setup.

use core::array::from_fn;
use std::hint::black_box;

use criterion::Criterion;

use la_stack::{LaError, Matrix, Vector, gram_matrix};

#[path = "common/bench_utils.rs"]
mod bench_utils;
use bench_utils::OrAbort;

/// Fixture families whose labels and construction must agree.
enum Scenario {
    Orthogonal,
    Dependent,
    NearDependent,
    MixedScale,
}

impl Scenario {
    /// Stable label used in Criterion result paths.
    const fn name(&self) -> &'static str {
        match self {
            Self::Orthogonal => "orthogonal",
            Self::Dependent => "dependent",
            Self::NearDependent => "near_dependent",
            Self::MixedScale => "mixed_scale",
        }
    }

    /// Small integer entries keep the independent matrix-product oracle exact.
    fn entry(&self, row: usize, coordinate: usize) -> i16 {
        let diagonal = i16::from(row == coordinate);
        match self {
            Self::Orthogonal => diagonal,
            Self::Dependent => 1,
            Self::NearDependent => 256 + diagonal,
            Self::MixedScale => 257 * diagonal - 1,
        }
    }
}

fn hand_written<const M: usize, const N: usize>(
    vectors: &[Vector<N>; M],
) -> Result<Matrix<M>, LaError> {
    let mut matrix = Matrix::zero();
    for (i, left) in vectors.iter().enumerate() {
        for (j, right) in vectors.iter().enumerate().skip(i) {
            let value = left.dot(right)?;
            matrix.set(i, j, value)?;
            matrix.set(j, i, value)?;
        }
    }
    Ok(matrix)
}

fn register<const M: usize, const N: usize>(c: &mut Criterion) {
    for scenario in [
        Scenario::Orthogonal,
        Scenario::Dependent,
        Scenario::NearDependent,
        Scenario::MixedScale,
    ] {
        let integers: [[i16; N]; M] = from_fn(|i| from_fn(|k| scenario.entry(i, k)));
        let vectors = integers
            .map(|row| Vector::try_new(row.map(f64::from)).or_abort("Gram benchmark input"));
        let expected = Matrix::try_from_rows(from_fn(|i| {
            from_fn(|j| {
                let value: i32 = (0..N)
                    .map(|k| i32::from(integers[i][k]) * i32::from(integers[j][k]))
                    .sum();
                f64::from(value)
            })
        }))
        .or_abort("Gram oracle");
        assert_eq!(gram_matrix(&vectors).or_abort("Gram validation"), expected);
        assert_eq!(
            hand_written(&vectors).or_abort("hand-written validation"),
            expected
        );
        let mut group = c.benchmark_group(format!("gram/{M}x{N}/{}", scenario.name()));
        group.bench_function("la_stack", |b| {
            b.iter(|| black_box(gram_matrix(black_box(&vectors)).or_abort("Gram construction")));
        });
        group.bench_function("hand_written", |b| {
            b.iter(|| black_box(hand_written(black_box(&vectors)).or_abort("Gram assembly")));
        });
        group.finish();
    }
}

fn main() {
    let mut c = Criterion::default().configure_from_args();
    register::<2, 2>(&mut c);
    register::<1, 2>(&mut c);
    register::<3, 3>(&mut c);
    register::<2, 3>(&mut c);
    register::<4, 4>(&mut c);
    register::<3, 4>(&mut c);
    register::<5, 5>(&mut c);
    register::<4, 5>(&mut c);
    register::<6, 6>(&mut c);
    register::<5, 6>(&mut c);
    register::<7, 7>(&mut c);
    register::<6, 7>(&mut c);
    register::<8, 8>(&mut c);
    register::<7, 8>(&mut c);
    c.final_summary();
}
