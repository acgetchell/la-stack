# la-stack

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18158926.svg)](https://doi.org/10.5281/zenodo.18158926)
[![Crates.io](https://img.shields.io/crates/v/la-stack.svg)](https://crates.io/crates/la-stack)
[![Downloads](https://img.shields.io/crates/d/la-stack.svg)](https://crates.io/crates/la-stack)
[![License](https://img.shields.io/crates/l/la-stack.svg)](https://github.com/acgetchell/la-stack/blob/v0.4.3/LICENSE)
[![Docs.rs](https://docs.rs/la-stack/badge.svg)](https://docs.rs/la-stack)
[![CI](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml)
[![rust-clippy analyze][clippy-badge]][clippy-workflow]
[![codecov](https://codecov.io/gh/acgetchell/la-stack/graph/badge.svg?token=4eKXa5QjuZ)](https://codecov.io/gh/acgetchell/la-stack)
[![Audit dependencies][audit-badge]][audit-workflow]

![la-stack](https://raw.githubusercontent.com/acgetchell/la-stack/main/docs/assets/la-stack.jpg)

Fast, stack-allocated linear algebra for fixed dimensions in Rust.

This crate grew from the need to support [`delaunay`](https://crates.io/crates/delaunay) with fast, stack-allocated linear algebra primitives and algorithms
while keeping the API intentionally small and explicit.

## 📐 Introduction

`la-stack` provides a handful of const-generic, stack-backed building blocks:

- `Vector<const D: usize>` for fixed-length `f64` vectors backed by `[f64; D]`
- `Matrix<const D: usize>` for fixed-size square `f64` matrices backed by `[[f64; D]; D]`
- `Lu<const D: usize>` for LU factorization with partial pivoting (solve + det)
- `Ldlt<const D: usize>` for LDLT factorization without pivoting (solve + det; symmetric SPD/PSD)

## ✨ Design goals

- ✅ `Copy` types where possible
- ✅ Const-generic dimensions (no dynamic sizes)
- ✅ `const fn` where possible (compile-time evaluation of determinants, dot products, etc.)
- ✅ Explicit algorithms (LU, solve, determinant)
- ✅ Robust geometric predicates via optional exact arithmetic (`det_sign_exact`, `det_errbound`)
- ✅ Exact linear system solve via optional arbitrary-precision arithmetic (`solve_exact`, strict/rounded f64 conversions)
- ✅ No runtime dependencies by default (optional features may add deps)
- ✅ Stack storage only (no heap allocation in core types)
- ✅ `unsafe` forbidden

See [CHANGELOG.md](https://github.com/acgetchell/la-stack/blob/v0.4.3/CHANGELOG.md)
for release history and
[docs/roadmap.md](https://github.com/acgetchell/la-stack/blob/v0.4.3/docs/roadmap.md)
for current release planning.

## 🚫 Anti-goals

- Bare-metal performance: see [`blas-src`](https://crates.io/crates/blas-src),
  [`lapack-src`](https://crates.io/crates/lapack-src), or [`openblas-src`](https://crates.io/crates/openblas-src)
- Comprehensive: use [`nalgebra`](https://crates.io/crates/nalgebra) if you need a full-featured library
- Large matrices/dimensions with parallelism: use [`faer`](https://crates.io/crates/faer) if you need this
- Alternate floating-point scalar families: `la-stack` supports `f64` and optional exact arithmetic, not `f32` / `f16` APIs

## ✅ Use this crate when

- Your matrices and vectors have small, fixed dimensions known at compile time
- Stack allocation and `Copy` value semantics fit your data flow
- You want explicit LU / LDLT / determinant APIs rather than a broad algebra toolkit
- You need exact determinants, exact determinant signs, or exact linear solves
  for fixed-size systems
- Robust predicates matter for geometry-style workloads near degeneracy
- You prefer a default build with no runtime dependencies

## 🔢 Scalar types

The scalar model is intentionally limited to `f64` for floating-point work and
exact rationals behind the optional `"exact"` feature. This matches the crate's
focus on small, robustness-sensitive numerical and computational geometry
workloads. When `f64` precision is insufficient (e.g. near-degenerate geometric
configurations), the optional `"exact"` feature provides arbitrary-precision
arithmetic via `BigRational` (see below).

Lower-precision `f32` / `f16` throughput-oriented workloads are outside the
crate's scope; they usually indicate large-matrix or accelerator-oriented use
cases better served by broader linear-algebra libraries.

## 🚀 Quickstart

Add this to your `Cargo.toml`:

```toml
[dependencies]
la-stack = "0.4.3"
```

### Feature flags

- `default`: no runtime dependencies
- `exact`: `BigRational` exact determinant and solve APIs
- `bench`: cfg-only gate for benchmark fixtures and benchmark-input tests;
  benchmark libraries remain development dependencies

Solve a 5×5 system via LU:

```rust
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // This system requires pivoting (a[0][0] = 0), so it's a good LU demo.
    // A = J - I: zeros on diagonal, ones elsewhere.
    let a = Matrix::<5>::try_from_rows([
        [0.0, 1.0, 1.0, 1.0, 1.0],
        [1.0, 0.0, 1.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 0.0],
    ])?;

    let b = Vector::<5>::try_new([14.0, 13.0, 12.0, 11.0, 10.0])?;

    let lu = a.lu(DEFAULT_SINGULAR_TOL)?;
    let x = lu.solve(b)?.into_array();

    // Floating-point rounding is expected; compare with a tolerance.
    let expected = [1.0, 2.0, 3.0, 4.0, 5.0];
    for (x_i, e_i) in x.iter().zip(expected.iter()) {
        assert!((*x_i - *e_i).abs() <= 1e-12);
    }

    Ok(())
}
```

Compute a determinant for a symmetric SPD matrix via LDLT (no pivoting).

For symmetric positive-definite matrices, `LDL^T` is essentially a square-root-free form of the Cholesky decomposition
(you can recover a Cholesky factor by absorbing `sqrt(D)` into `L`):

```rust
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // This matrix is symmetric positive-definite (A = L*L^T) so LDLT works without pivoting.
    let a = Matrix::<5>::try_from_rows([
        [1.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 2.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 2.0, 1.0, 0.0],
        [0.0, 0.0, 1.0, 2.0, 1.0],
        [0.0, 0.0, 0.0, 1.0, 2.0],
    ])?;

    let ldlt = match a.ldlt(DEFAULT_SINGULAR_TOL) {
        Ok(ldlt) => ldlt,
        Err(err @ LaError::Asymmetric {
            row,
            col,
            upper,
            lower,
            allowed_abs_diff,
            ..
        }) => {
            eprintln!(
                "LDLT mismatch at ({row}, {col}): {upper} vs {lower} (allowed {allowed_abs_diff})"
            );
            return Err(err);
        }
        Err(err) => return Err(err),
    };

    let det = ldlt.det()?;
    assert!((det - 1.0).abs() <= 1e-12);

    Ok(())
}
```

> ⚠️ **LDLT invariant:** The input matrix must be **exactly symmetric**: every
> mirrored pair must compare equal (`+0.0 == -0.0` is accepted). Asymmetric
> inputs passed to
> [`Matrix::ldlt`](https://docs.rs/la-stack/latest/la_stack/struct.Matrix.html#method.ldlt)
> return a typed `LaError::Asymmetric` containing both observed values and the
> required allowed difference of zero. The tolerance-based
> [`Matrix::first_asymmetry`](https://docs.rs/la-stack/latest/la_stack/struct.Matrix.html#method.first_asymmetry)
> and `Matrix::is_symmetric` methods remain useful diagnostics, but do not prove
> the exact precondition required by LDLT. Fall back to `lu()` if your matrices
> may not be symmetric at all. A negative LDLT diagonal or a zero diagonal with nonzero
> remaining coupling returns `LaError::NotPositiveSemidefinite` with a typed
> `PositiveSemidefiniteViolation`. An uncoupled zero or other non-negative pivot
> at or below the caller's tolerance returns `LaError::Singular` with a
> numerical `SingularityReason`.

## ⚡ Compile-time determinants (D ≤ 4)

`det_direct()` is a `const fn` providing closed-form determinants for D=0–4,
using fused multiply-add where applicable. `Matrix::<0>::zero().det_direct()`
returns `Ok(Some(1.0))` (the empty-product convention). For D=1–4, cofactor
expansion bypasses LU factorization entirely. This enables compile-time
evaluation when inputs are known:

```rust
use la_stack::prelude::*;

// Evaluated entirely at compile time — no runtime cost.
const DET: Result<Option<f64>, LaError> = match Matrix::<4>::try_from_rows([
    [2.0, 0.0, 0.0, 0.0],
    [0.0, 3.0, 0.0, 0.0],
    [0.0, 0.0, 5.0, 0.0],
    [0.0, 0.0, 0.0, 7.0],
]) {
    Ok(matrix) => matrix.det_direct(),
    Err(err) => Err(err),
};

fn main() -> Result<(), LaError> {
    assert_eq!(DET?, Some(210.0));
    Ok(())
}
```

The public `det()` method automatically dispatches through the closed-form path
for D ≤ 4 and falls back to zero-tolerance LU for D ≥ 5. Tiny nonzero
determinants are not flattened by a configured pivot tolerance. The LU fallback
returns `LaError::Singular` when floating-point elimination cannot produce a
non-zero pivot; it does not misreport that numerical failure as an exact zero.
Use `lu()` directly when you need a different tolerance policy, and use the
exact determinant APIs when exact singularity classification matters.

## 🔬 Exact arithmetic (`"exact"` feature)

The default build has **zero runtime dependencies**. Enable the optional
`exact` Cargo feature to add exact arithmetic methods using arbitrary-precision
rationals (this pulls in `num-bigint`, `num-rational`, and `num-traits` for
`BigRational`):

```toml
[dependencies]
la-stack = { version = "0.4.3", features = ["exact"] }
```

**Determinants:**

- **`det_exact()`** — returns the exact determinant as a `BigRational`
- **`det_exact_f64()`** — returns the exact determinant as `f64` only when
  it is exactly representable (or `LaError::Unrepresentable` otherwise)
- **`det_exact_rounded_f64()`** — returns the exact determinant rounded to a
  finite `f64` using IEEE 754 round-to-nearest, ties-to-even
- **`det_sign_exact()`** — infallibly returns the provably correct
  `DeterminantSign` variant (`Negative`, `Zero`, or `Positive`)

**Linear system solve:**

- **`solve_exact(b)`** — solves `Ax = b` exactly, returning `[BigRational; D]`
- **`solve_exact_f64(b)`** — solves `Ax = b` exactly, returning `Vector<D>` only when
  every component is exactly representable as `f64`
- **`solve_exact_rounded_f64(b)`** — solves `Ax = b` exactly, returning each
  component rounded to finite `f64` using IEEE 754 round-to-nearest,
  ties-to-even
- **`ExactF64Conversion`** — converts an existing exact determinant or solution
  under the strict or rounded contract without repeating exact elimination

```rust,ignore
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // Exact determinant
    let m = Matrix::<3>::try_from_rows([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ])?;
    assert_eq!(m.det_sign_exact(), DeterminantSign::Zero); // exactly singular

    let det = m.det_exact()?;
    assert_eq!(det, BigRational::from_integer(0.into())); // exact zero
    let det_f64 = det.try_to_f64()?;
    assert_eq!(det_f64, 0.0);

    // If strict exact-to-f64 conversion would require rounding, opt in
    // explicitly with the rounded API.
    let inexact = Matrix::<2>::try_from_rows([
        [1.0 + f64::EPSILON, 0.0],
        [0.0, 1.0 - f64::EPSILON],
    ])?;
    let exact_det = inexact.det_exact()?;
    let rounded_det = match exact_det.try_to_f64() {
        Ok(det) => det,
        Err(err) if err.requires_rounding() => exact_det.to_rounded_f64()?,
        Err(err) => return Err(err),
    };
    assert_eq!(rounded_det.to_bits(), 1.0f64.to_bits());

    // If the exact determinant cannot fit in f64, keep the BigRational value.
    let big = f64::MAX / 2.0;
    let huge = Matrix::<3>::try_from_rows([
        [0.0, 0.0, 1.0],
        [big, 0.0, 1.0],
        [0.0, big, 1.0],
    ])?;
    let huge_det = huge.det_exact()?;
    assert_eq!(
        huge_det
            .try_to_f64()
            .err()
            .and_then(|err| err.unrepresentable_reason()),
        Some(UnrepresentableReason::NotFinite)
    );
    println!("exact determinant = {huge_det}");

    // Exact linear system solve
    let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    let b = Vector::<2>::try_new([5.0, 11.0])?;
    let exact_x = a.solve_exact(b)?;
    let x = exact_x.try_to_f64()?.into_array();
    assert!((x[0] - 1.0).abs() <= f64::EPSILON);
    assert!((x[1] - 2.0).abs() <= f64::EPSILON);

    Ok(())
}
```

With the `exact` feature enabled, `DeterminantSign`, `ExactF64Conversion`,
`BigInt`, and `BigRational` are re-exported from the crate root and prelude,
alongside the most commonly needed `num-traits` items (`FromPrimitive`,
`ToPrimitive`, `Signed`). This lets consumers construct exact values
(`BigRational::from_f64`, `from_i64`), query sign (`is_positive` /
`is_negative`), and convert back to `f64` (`try_to_f64`, `to_rounded_f64`, or
the raw `to_f64`) with a single
`use la_stack::prelude::*;` — no need to add `num-bigint`, `num-rational`,
or `num-traits` to their own `Cargo.toml`. Use
`DeterminantSign::as_i8()` only when numeric −1/0/+1 interoperability is
required.

For `det_sign_exact()`, D ≤ 4 matrices first use a fast f64 filter
(error-bounded `det_direct()`) when its rounded intermediates stay in the normal
range or are exact structural zeros. An inconclusive filter falls back to the
same direct determinant expansion in `BigInt`. D ≥ 5 skips the closed-form
filter and uses fraction-free Bareiss elimination in `BigInt`.
Because `Matrix` stores only finite entries, arithmetic range failures in the
filter are inconclusive rather than errors and the exact fallback is total.

### Adaptive precision with `det_direct_with_errbound()`

`det_direct_with_errbound()` returns a closed-form determinant together with
the conservative absolute error bound used by the fast filter, computed from
one shared traversal. It returns `None` when a D ≤ 4 computation may be
affected by gradual underflow, as well as for unsupported D ≥ 5 dimensions.
This method does NOT require the `exact` feature — it uses pure f64 arithmetic
and is available by default. Use `det_errbound()` when only the bound is needed.
The paired API enables custom adaptive-precision logic for geometric predicates:

```rust,ignore
use la_stack::prelude::*;

fn adaptive_det_sign<const D: usize>(
    matrix: &Matrix<D>,
) -> DeterminantSign {
    if let Ok(Some(estimate)) = matrix.det_direct_with_errbound() {
        if estimate.determinant().abs() > estimate.absolute_error_bound() {
            return if estimate.determinant() > 0.0 {
                DeterminantSign::Positive
            } else {
                DeterminantSign::Negative
            };
        }
    }

    matrix.det_sign_exact()
}

fn main() -> Result<(), LaError> {
    let identity = Matrix::<3>::identity();
    assert_eq!(
        adaptive_det_sign(&identity),
        DeterminantSign::Positive
    );

    // A zero determinant cannot pass the f64 sign filter, so this exercises
    // the exact fallback.
    let singular = Matrix::<3>::try_from_rows([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ])?;
    assert_eq!(adaptive_det_sign(&singular), DeterminantSign::Zero);

    // The f64 filter overflows for this finite matrix, but the exact fallback
    // still resolves its positive determinant sign.
    let big = f64::MAX / 2.0;
    let overflowing = Matrix::<3>::try_from_rows([
        [0.0, 0.0, 1.0],
        [big, 0.0, 1.0],
        [0.0, big, 1.0],
    ])?;
    assert_eq!(
        adaptive_det_sign(&overflowing),
        DeterminantSign::Positive
    );

    Ok(())
}
```

The error coefficients (`ERR_COEFF_2`, `ERR_COEFF_3`, `ERR_COEFF_4`) are the
dimension-specific constants behind that bound. In plain terms, they answer:
"how many machine-epsilon-sized rounding mistakes can this closed-form
determinant formula accumulate?" To get an absolute error bound, `det_errbound()`
multiplies the coefficient by a size measure of the matrix entries, the
**absolute Leibniz sum**, equivalently the permanent of `|A|`:

```text
p(|A|) = sum over determinant terms of product of absolute values
```

For a 2×2 matrix `[[a, b], [c, d]]`, that scale is `|a*d| + |b*c|`, so:

```text
|det_direct(A) - det_exact(A)| <= ERR_COEFF_2 * (|a*d| + |b*c|)
```

The coefficients are not tolerances and are not meant to be tuned by callers;
they are conservative constants derived from the fixed D ≤ 4 formulas and their
floating-point rounding chains when gradual underflow is absent. They are
explicit crate-root exports for
advanced users who want to compose the same bound themselves:
`use la_stack::{ERR_COEFF_2, ERR_COEFF_3, ERR_COEFF_4};`. They intentionally stay
out of the common prelude.

## 🧩 API at a glance

| Type | Storage | Purpose | Key methods |
|---|---|---|---|
| `Vector<D>` | `[f64; D]` | Finite fixed-length vector for input and computation | `try_new`, `as_array`, `into_array`, `dot`, `norm2_sq` |
| `Matrix<D>` | `[[f64; D]; D]` | Finite square matrix for input and computation | See below |
| `DeterminantWithErrorBound` | two private `f64` fields | Paired direct determinant and certified absolute bound | `determinant`, `absolute_error_bound` |
| `Lu<D>` | `Matrix<D>` + pivot array | Factorization for solves/det | `solve`, `det` |
| `Ldlt<D>` | `Matrix<D>` | Factorization for symmetric SPD/PSD solves/det | `solve`, `det` |
| `Tolerance` | finite non-negative `f64` | Validated numerical threshold | `try_new`, `get` |
| `LaError` | typed variants and reasons | Structured, actionable failure reporting | See error enums below |
| `DeterminantSign`¹ | enum | Exact determinant sign | `as_i8` |

Storage shown above reflects the intentional `f64` scalar model.

`Matrix<D>` key methods: `as_rows`, `into_rows`, `lu`, `ldlt`, `det`,
`det_direct`, `det_direct_with_errbound`, `det_errbound`,
`det_exact`¹, `det_exact_f64`¹, `det_exact_rounded_f64`¹, `det_sign_exact`¹,
`solve_exact`¹, `solve_exact_f64`¹, `solve_exact_rounded_f64`¹.
Matrix and vector constructors validate non-finite inputs at public API
boundaries. After construction, `Matrix<D>` and `Vector<D>` carry that
finite-storage invariant directly, so factorization kernels do not repeat an
O(D²) input scan. Computed factor matrices are still checked before they become
observable results.

`Matrix::as_rows` and `Vector::as_array` borrow their validated backing arrays;
`Matrix::into_rows` and `Vector::into_array` consume the value and return the
owned fixed-size arrays.

`Matrix::get` returns `Option` for bounds-only access; `Matrix::try_get`
preserves invalid coordinates in `LaError`. The single fallible `Matrix::set`
validates both bounds and finiteness before mutating the matrix.

`LaError` and its reason/location enums are non-exhaustive. Numerical
singularity records the [`FactorizationKind`](https://docs.rs/la-stack/latest/la_stack/enum.FactorizationKind.html),
observed pivot magnitude, and tolerance, while exact-arithmetic singularity is
identified separately. `LaError::NonFinite` retains the crate-wide non-finite
contract but uses `NonFiniteOrigin`, `NonFiniteLocation`, and
`ArithmeticOperation` to distinguish invalid inputs from computed overflow.
`InvalidToleranceReason` distinguishes negative from non-finite tolerances, and
`PositiveSemidefiniteViolation` distinguishes negative LDLT pivots from a zero
pivot with nonzero coupling. Match these public enums with a wildcard and use
`..` for struct-style variants so future error context can be added without
breaking callers.

¹ Requires `features = ["exact"]`.

## 📊 Benchmarks (vs nalgebra/faer)

![LU solve (factor + solve): median time vs dimension][lu-solve-benchmark]

Raw data:
[docs/assets/bench/vs_linalg_lu_solve_median.csv](https://github.com/acgetchell/la-stack/blob/v0.4.3/docs/assets/bench/vs_linalg_lu_solve_median.csv)
Historical provenance status:
[docs/assets/bench/vs_linalg_lu_solve_median.provenance.json][benchmark-provenance]

Representative benchmark: `lu_solve` factors the matrix and solves one
right-hand side. Median time is lower-is-better, and the “la-stack vs
nalgebra/faer” columns show the % time reduction relative to each baseline
(positive means the recorded la-stack median is lower). These are descriptive
point-estimate ratios, not statistical significance claims or an aggregate score
across operations.

Timings count only when the implementation preserves the documented
correctness guarantees and invariants. Performance claims require comparable
before-and-after evidence using the same inputs, configuration, and environment.
This v0.4.3 snapshot predates deterministic measurement-provenance capture, so
its CPU, operating system, Rust toolchain, exact measured source state,
dependency lock digest, and Criterion configuration are unavailable. The CSV
preserves confidence bounds, but without the missing configuration and
environment they do not make the result reproducible across environments. Treat
it as a historical snapshot, not reproducible cross-environment evidence. Future
`just plot-vs-linalg-readme` publications run the benchmark-input correctness
gate, require complete canonical-dimension coverage, and write deterministic
JSON provenance beside the CSV and SVG.

For the full per-kernel comparison methodology, input construction, and
release-comparison workflow details, see
[docs/BENCHMARKING.md](https://github.com/acgetchell/la-stack/blob/v0.4.3/docs/BENCHMARKING.md).
For the current release-to-release performance snapshot, see
[docs/PERFORMANCE.md](https://github.com/acgetchell/la-stack/blob/v0.4.3/docs/PERFORMANCE.md).

<!-- BENCH_TABLE:lu_solve:median:new:BEGIN -->

| D | la-stack median (ns) | nalgebra median (ns) | faer median (ns) | la-stack vs nalgebra | la-stack vs faer |
|---:|--------------------:|--------------------:|----------------:|---------------------:|----------------:|
| 2 | 2.044 | 4.542 | 143.958 | +55.0% | +98.6% |
| 3 | 9.596 | 23.599 | 185.466 | +59.3% | +94.8% |
| 4 | 23.338 | 50.717 | 210.976 | +54.0% | +88.9% |
| 5 | 45.368 | 69.065 | 277.564 | +34.3% | +83.7% |
| 8 | 127.861 | 164.412 | 364.864 | +22.2% | +65.0% |
| 16 | 631.997 | 663.822 | 882.674 | +4.8% | +28.4% |
| 32 | 2,745.604 | 2,424.540 | 2,867.431 | -13.2% | +4.2% |
| 64 | 17,543.034 | 14,747.731 | 12,266.271 | -19.0% | -43.0% |

<!-- BENCH_TABLE:lu_solve:median:new:END -->

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

Install Rust 1.97.0 through [rustup](https://rustup.rs/), Git, Python 3.14,
[`uv` 0.11.28](https://docs.astral.sh/uv/), and `jq`. Then install the pinned
`just` release from its locked dependency graph:

```bash
cargo install --locked just --version 1.56.0
just setup        # install/verify dev tools + sync Python deps
just check        # lint/validate (non-mutating)
just fix          # apply auto-fixes (mutating)
just ci           # lint + tests + examples + bench compile
```

The repository uses `cargo-nextest` for runnable Rust tests, `cargo-machete`
for unused-dependency checks, `rumdl` for Markdown, `dprint` plus `yamllint`
for YAML/CFF, `taplo` for TOML, and `typos` for spelling. Python 3.14 support
tooling is locked with `uv` and checked by Ruff, Ty, and Semgrep. GitHub Actions
references are SHA-pinned, restricted to an explicit allowlist, and kept with
readable version comments for review.

CI runs `just ci` on Ubuntu, macOS, and Windows to keep platform coverage
aligned with the local comprehensive validation path.

For coverage commands and report locations, see
[`docs/COVERAGE.md`](https://github.com/acgetchell/la-stack/blob/v0.4.3/docs/COVERAGE.md).
For the full contributor workflow, see
[CONTRIBUTING.md](https://github.com/acgetchell/la-stack/blob/v0.4.3/CONTRIBUTING.md).

## 📝 Citation

If you use this library in academic work, please cite it using
[CITATION.cff](https://github.com/acgetchell/la-stack/blob/v0.4.3/CITATION.cff)
(or GitHub's "Cite this repository" feature). Tagged releases are archived on
Zenodo.

## 📚 References

For canonical references to the algorithms used by this crate, see
[REFERENCES.md](https://github.com/acgetchell/la-stack/blob/v0.4.3/REFERENCES.md).

## 🤖 AI Agents

AI coding assistants should read
[AGENTS.md](https://github.com/acgetchell/la-stack/blob/v0.4.3/AGENTS.md)
before proposing or applying changes. See
[CONTRIBUTING.md](https://github.com/acgetchell/la-stack/blob/v0.4.3/CONTRIBUTING.md)
for the repository's AI-assisted development note.

## 📄 License

BSD 3-Clause License. See [LICENSE](https://github.com/acgetchell/la-stack/blob/v0.4.3/LICENSE).

[audit-badge]: https://github.com/acgetchell/la-stack/actions/workflows/audit.yml/badge.svg
[audit-workflow]: https://github.com/acgetchell/la-stack/actions/workflows/audit.yml
[benchmark-provenance]: https://github.com/acgetchell/la-stack/blob/v0.4.3/docs/assets/bench/vs_linalg_lu_solve_median.provenance.json
[clippy-badge]: https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml/badge.svg
[clippy-workflow]: https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml
[lu-solve-benchmark]: https://raw.githubusercontent.com/acgetchell/la-stack/v0.4.3/docs/assets/bench/vs_linalg_lu_solve_median.svg
