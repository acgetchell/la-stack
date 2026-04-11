# la-stack

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18158926.svg)](https://doi.org/10.5281/zenodo.18158926)
[![Crates.io](https://img.shields.io/crates/v/la-stack.svg)](https://crates.io/crates/la-stack)
[![Downloads](https://img.shields.io/crates/d/la-stack.svg)](https://crates.io/crates/la-stack)
[![License](https://img.shields.io/crates/l/la-stack.svg)](./LICENSE)
[![Docs.rs](https://docs.rs/la-stack/badge.svg)](https://docs.rs/la-stack)
[![CI](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml)
[![rust-clippy analyze](https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml)
[![codecov](https://codecov.io/gh/acgetchell/la-stack/graph/badge.svg?token=4eKXa5QjuZ)](https://codecov.io/gh/acgetchell/la-stack)
[![Audit dependencies](https://github.com/acgetchell/la-stack/actions/workflows/audit.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/audit.yml)
[![Codacy Security Scan](https://github.com/acgetchell/la-stack/actions/workflows/codacy.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/codacy.yml)

Fast, stack-allocated linear algebra for fixed dimensions in Rust.

This crate grew from the need to support [`delaunay`](https://crates.io/crates/delaunay) with fast, stack-allocated linear algebra primitives and algorithms
while keeping the API intentionally small and explicit.

## 📐 Introduction

`la-stack` provides a handful of const-generic, stack-backed building blocks:

- `Vector<const D: usize>` for fixed-length vectors (`[f64; D]` today)
- `Matrix<const D: usize>` for fixed-size square matrices (`[[f64; D]; D]` today)
- `Lu<const D: usize>` for LU factorization with partial pivoting (solve + det)
- `Ldlt<const D: usize>` for LDLT factorization without pivoting (solve + det; symmetric SPD/PSD)

## ✨ Design goals

- ✅ `Copy` types where possible
- ✅ Const-generic dimensions (no dynamic sizes)
- ✅ `const fn` where possible (compile-time evaluation of determinants, dot products, etc.)
- ✅ Explicit algorithms (LU, solve, determinant)
- ✅ Robust geometric predicates via optional exact arithmetic (`det_sign_exact`, `det_errbound`)
- ✅ Exact linear system solve via optional arbitrary-precision arithmetic (`solve_exact`, `solve_exact_f64`)
- ✅ No runtime dependencies by default (optional features may add deps)
- ✅ Stack storage only (no heap allocation in core types)
- ✅ `unsafe` forbidden

See [CHANGELOG.md](CHANGELOG.md) for details.

## 🚫 Anti-goals

- Bare-metal performance: see [`blas-src`](https://crates.io/crates/blas-src), [`lapack-src`](https://crates.io/crates/lapack-src), [`openblas-src`](https://crates.io/crates/openblas-src)
- Comprehensive: use [`nalgebra`](https://crates.io/crates/nalgebra) if you need a full-featured library
- Large matrices/dimensions with parallelism: use [`faer`](https://crates.io/crates/faer) if you need this

## 🔢 Scalar types

The core types use `f64`. When f64 precision is insufficient (e.g. near-degenerate
geometric configurations), the optional `"exact"` feature provides arbitrary-precision
arithmetic via `BigRational` (see below).

## 🚀 Quickstart

Add this to your `Cargo.toml`:

```toml
[dependencies]
la-stack = "0.4.0"
```

Solve a 5×5 system via LU:

```rust
use la_stack::prelude::*;

// This system requires pivoting (a[0][0] = 0), so it's a good LU demo.
// A = J - I: zeros on diagonal, ones elsewhere.
let a = Matrix::<5>::from_rows([
    [0.0, 1.0, 1.0, 1.0, 1.0],
    [1.0, 0.0, 1.0, 1.0, 1.0],
    [1.0, 1.0, 0.0, 1.0, 1.0],
    [1.0, 1.0, 1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0, 1.0, 0.0],
]);

let b = Vector::<5>::new([14.0, 13.0, 12.0, 11.0, 10.0]);

let lu = a.lu(DEFAULT_PIVOT_TOL).unwrap();
let x = lu.solve_vec(b).unwrap().into_array();

// Floating-point rounding is expected; compare with a tolerance.
let expected = [1.0, 2.0, 3.0, 4.0, 5.0];
for (x_i, e_i) in x.iter().zip(expected.iter()) {
    assert!((*x_i - *e_i).abs() <= 1e-12);
}
```

Compute a determinant for a symmetric SPD matrix via LDLT (no pivoting).

For symmetric positive-definite matrices, `LDL^T` is essentially a square-root-free form of the Cholesky decomposition
(you can recover a Cholesky factor by absorbing `sqrt(D)` into `L`):

```rust
use la_stack::prelude::*;

// This matrix is symmetric positive-definite (A = L*L^T) so LDLT works without pivoting.
let a = Matrix::<5>::from_rows([
    [1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 2.0, 1.0, 0.0, 0.0],
    [0.0, 1.0, 2.0, 1.0, 0.0],
    [0.0, 0.0, 1.0, 2.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 2.0],
]);

let det = a.ldlt(DEFAULT_SINGULAR_TOL).unwrap().det();
assert!((det - 1.0).abs() <= 1e-12);
```

## ⚡ Compile-time determinants (D ≤ 4)

`det_direct()` is a `const fn` providing closed-form determinants for D=0–4,
using fused multiply-add where applicable. `Matrix::<0>::zero().det_direct()`
returns `Some(1.0)` (the empty-product convention). For D=1–4, cofactor
expansion bypasses LU factorization entirely. This enables compile-time
evaluation when inputs are known:

```rust
use la_stack::prelude::*;

// Evaluated entirely at compile time — no runtime cost.
const DET: Option<f64> = {
    let m = Matrix::<3>::from_rows([
        [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 5.0],
    ]);
    m.det_direct()
};
assert_eq!(DET, Some(30.0));
```

The public `det()` method automatically dispatches through the closed-form path
for D ≤ 4 and falls back to LU for D ≥ 5 — no API change needed.

## 🔬 Exact arithmetic (`"exact"` feature)

The default build has **zero runtime dependencies**.  Enable the optional
`exact` Cargo feature to add exact arithmetic methods using arbitrary-precision
rationals (this pulls in `num-bigint`, `num-rational`, and `num-traits` for
`BigRational`):

```toml
[dependencies]
la-stack = { version = "0.4.0", features = ["exact"] }
```

**Determinants:**

- **`det_exact()`** — returns the exact determinant as a `BigRational`
- **`det_exact_f64()`** — returns the exact determinant converted to the nearest `f64`
- **`det_sign_exact()`** — returns the provably correct sign (−1, 0, or +1)

**Linear system solve:**

- **`solve_exact(b)`** — solves `Ax = b` exactly, returning `[BigRational; D]`
- **`solve_exact_f64(b)`** — solves `Ax = b` exactly, converting the result to `Vector<D>` (f64)

```rust,ignore
use la_stack::prelude::*;

// Exact determinant
let m = Matrix::<3>::from_rows([
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
]);
assert_eq!(m.det_sign_exact().unwrap(), 0); // exactly singular

let det = m.det_exact().unwrap();
assert_eq!(det, BigRational::from_integer(0.into())); // exact zero

// Exact linear system solve
let a = Matrix::<2>::from_rows([[1.0, 2.0], [3.0, 4.0]]);
let b = Vector::<2>::new([5.0, 11.0]);
let x = a.solve_exact_f64(b).unwrap().into_array();
assert!((x[0] - 1.0).abs() <= f64::EPSILON);
assert!((x[1] - 2.0).abs() <= f64::EPSILON);
```

`BigRational` is re-exported from the crate root and prelude when the `exact`
feature is enabled, so consumers don't need to depend on `num-rational` directly.

For `det_sign_exact()`, D ≤ 4 matrices use a fast f64 filter (error-bounded
`det_direct()`) that resolves the sign without allocating.  Only near-degenerate
or large (D ≥ 5) matrices fall through to the exact Bareiss algorithm.

### Adaptive precision with `det_errbound()`

`det_errbound()` returns the conservative absolute error bound used by the fast
filter. This method does NOT require the `exact` feature — it uses pure f64
arithmetic and is available by default. This enables building custom
adaptive-precision logic for geometric predicates:

```rust,ignore
use la_stack::prelude::*;

let m = Matrix::<3>::identity();
if let Some(bound) = m.det_errbound() {
    let det = m.det_direct().unwrap();
    if det.abs() > bound {
        // f64 sign is guaranteed correct
        let sign = det.signum() as i8;
    } else {
        // Fall back to exact arithmetic (requires `exact` feature)
        let sign = m.det_sign_exact().unwrap();
    }
} else {
    // D ≥ 5: no fast filter, use exact directly (requires `exact` feature)
    let sign = m.det_sign_exact().unwrap();
}
```

The error coefficients (`ERR_COEFF_2`, `ERR_COEFF_3`, `ERR_COEFF_4`) are also
exposed for advanced use cases.

## 🧩 API at a glance

| Type | Storage | Purpose | Key methods |
|---|---|---|---|
| `Vector<D>` | `[f64; D]` | Fixed-length vector | `new`, `zero`, `dot`, `norm2_sq` |
| `Matrix<D>` | `[[f64; D]; D]` | Square matrix | See below |
| `Lu<D>` | `Matrix<D>` + pivot array | Factorization for solves/det | `solve_vec`, `det` |
| `Ldlt<D>` | `Matrix<D>` | Factorization for symmetric SPD/PSD solves/det | `solve_vec`, `det` |

Storage shown above reflects the `f64` implementation.

`Matrix<D>` key methods: `lu`, `ldlt`, `det`, `det_direct`, `det_errbound`,
`det_exact`¹, `det_exact_f64`¹, `det_sign_exact`¹, `solve_exact`¹, `solve_exact_f64`¹.

¹ Requires `features = ["exact"]`.

## 📋 Examples

The `examples/` directory contains small, runnable programs:

- **`solve_5x5`** — solve a 5×5 system via LU with partial pivoting
- **`det_5x5`** — determinant of a 5×5 matrix via LU
- **`ldlt_solve_3x3`** — solve a 3×3 symmetric positive definite system via LDLT
- **`const_det_4x4`** — compile-time 4×4 determinant via `det_direct()`
- **`exact_det_3x3`** — exact determinant value of a near-singular 3×3 matrix (requires `exact` feature)
- **`exact_sign_3x3`** — exact determinant sign of a near-singular 3×3 matrix (requires `exact` feature)
- **`exact_solve_3x3`** — exact solve of a near-singular 3×3 system vs f64 LU (requires `exact` feature)

```bash
just examples
# or individually:
cargo run --example solve_5x5
cargo run --example det_5x5
cargo run --example ldlt_solve_3x3
cargo run --example const_det_4x4
cargo run --features exact --example exact_det_3x3
cargo run --features exact --example exact_sign_3x3
cargo run --features exact --example exact_solve_3x3
```

## 🤝 Contributing

A short contributor workflow:

```bash
cargo install just
just setup        # install/verify dev tools + sync Python deps
just check        # lint/validate (non-mutating)
just fix          # apply auto-fixes (mutating)
just ci           # lint + tests + examples + bench compile
```

For the full set of developer commands, see `just --list` and `AGENTS.md`.

## 📝 Citation

If you use this library in academic work, please cite it using [CITATION.cff](CITATION.cff) (or GitHub's
"Cite this repository" feature). A Zenodo DOI will be added for tagged releases.

## 📚 References

For canonical references to the algorithms used by this crate, see [REFERENCES.md](REFERENCES.md).

## 📊 Benchmarks (vs nalgebra/faer)

![LU solve (factor + solve): median time vs dimension](docs/assets/bench/vs_linalg_lu_solve_median.svg)

Raw data: [docs/assets/bench/vs_linalg_lu_solve_median.csv](docs/assets/bench/vs_linalg_lu_solve_median.csv)

Summary (median time; lower is better). The “la-stack vs nalgebra/faer” columns show the % time reduction relative to each baseline (positive = la-stack faster):

<!-- BENCH_TABLE:lu_solve:median:new:BEGIN -->
| D | la-stack median (ns) | nalgebra median (ns) | faer median (ns) | la-stack vs nalgebra | la-stack vs faer |
|---:|--------------------:|--------------------:|----------------:|---------------------:|----------------:|
| 2 | 2.309 | 4.365 | 140.156 | +47.1% | +98.4% |
| 3 | 18.331 | 22.706 | 181.074 | +19.3% | +89.9% |
| 4 | 27.430 | 51.372 | 210.451 | +46.6% | +87.0% |
| 5 | 53.819 | 70.722 | 276.064 | +23.9% | +80.5% |
| 8 | 143.611 | 160.309 | 356.960 | +10.4% | +59.8% |
| 16 | 611.393 | 580.793 | 871.704 | -5.3% | +29.9% |
| 32 | 2,631.241 | 2,733.946 | 2,832.816 | +3.8% | +7.1% |
| 64 | 17,233.345 | 14,112.678 | 12,164.571 | -22.1% | -41.7% |
<!-- BENCH_TABLE:lu_solve:median:new:END -->

## 📄 License

BSD 3-Clause License. See [LICENSE](./LICENSE).
