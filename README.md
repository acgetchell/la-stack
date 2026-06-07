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

See [CHANGELOG.md](CHANGELOG.md) for release history and
[docs/roadmap.md](docs/roadmap.md) for current release planning.

## 🚫 Anti-goals

- Bare-metal performance: see [`blas-src`](https://crates.io/crates/blas-src), [`lapack-src`](https://crates.io/crates/lapack-src), [`openblas-src`](https://crates.io/crates/openblas-src)
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
la-stack = "0.4.2"
```

### Feature flags

- `default`: no runtime dependencies
- `exact`: `BigRational` exact determinant and solve APIs
- `bench`: Criterion, nalgebra, and faer for internal benchmarks

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
        Err(err @ LaError::Asymmetric { row, col, .. }) => {
            eprintln!("LDLT requires symmetry; first mismatch at ({row}, {col})");
            return Err(err);
        }
        Err(err) => return Err(err),
    };

    let det = ldlt.det()?;
    assert!((det - 1.0).abs() <= 1e-12);

    Ok(())
}
```

> ⚠️ **LDLT invariant:** The input matrix must be **symmetric**. Asymmetric
> inputs passed to
> [`Matrix::ldlt`](https://docs.rs/la-stack/latest/la_stack/struct.Matrix.html#method.ldlt)
> return a typed `LaError::Asymmetric` before factorization starts. Use
> [`Matrix::first_asymmetry`](https://docs.rs/la-stack/latest/la_stack/struct.Matrix.html#method.first_asymmetry)
> to locate the offending pair, or fall back to `lu()` if your matrices may not
> be symmetric at all. Symmetric inputs with a negative LDLT diagonal return
> `LaError::NotPositiveSemidefinite`; zero or too-small non-negative diagonals
> return `LaError::Singular`.

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
for D ≤ 4 and falls back to LU for D ≥ 5. Finite inputs return a floating-point
determinant estimate in every dimension; `det()` does not surface
`LaError::Singular`. Tiny nonzero determinants are not flattened by a pivot
tolerance. Use `lu()` directly when you need tolerance-aware singularity
detection or the pivot-column diagnostic from the factorization, and use the
exact determinant APIs when exact singularity classification matters.

## 🔬 Exact arithmetic (`"exact"` feature)

The default build has **zero runtime dependencies**. Enable the optional
`exact` Cargo feature to add exact arithmetic methods using arbitrary-precision
rationals (this pulls in `num-bigint`, `num-rational`, and `num-traits` for
`BigRational`):

```toml
[dependencies]
la-stack = { version = "0.4.2", features = ["exact"] }
```

**Determinants:**

- **`det_exact()`** — returns the exact determinant as a `BigRational`
- **`det_exact_f64()`** — returns the exact determinant as `f64` only when
  it is exactly representable (or `LaError::Unrepresentable` otherwise)
- **`det_exact_rounded_f64()`** — returns the exact determinant rounded to a
  finite `f64`
- **`det_sign_exact()`** — returns the provably correct sign (−1, 0, or +1)

**Linear system solve:**

- **`solve_exact(b)`** — solves `Ax = b` exactly, returning `[BigRational; D]`
- **`solve_exact_f64(b)`** — solves `Ax = b` exactly, returning `Vector<D>` only when
  every component is exactly representable as `f64`
- **`solve_exact_rounded_f64(b)`** — solves `Ax = b` exactly, returning each
  component rounded to finite `f64`

```rust,ignore
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // Exact determinant
    let m = Matrix::<3>::try_from_rows([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0],
    ])?;
    assert_eq!(m.det_sign_exact()?, 0); // exactly singular

    let det = m.det_exact()?;
    assert_eq!(det, BigRational::from_integer(0.into())); // exact zero
    let det_f64 = m.det_exact_f64()?;
    assert_eq!(det_f64, 0.0);

    // If strict exact-to-f64 conversion would require rounding, opt in
    // explicitly with the rounded API.
    let inexact = Matrix::<2>::try_from_rows([
        [1.0 + f64::EPSILON, 0.0],
        [0.0, 1.0 - f64::EPSILON],
    ])?;
    let rounded_det = match inexact.det_exact_f64() {
        Ok(det) => det,
        Err(err) if err.requires_rounding() => inexact.det_exact_rounded_f64()?,
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
        huge.det_exact_f64()
            .err()
            .and_then(|err| err.unrepresentable_reason()),
        Some(UnrepresentableReason::NotFinite)
    );
    println!("exact determinant = {huge_det}");

    // Exact linear system solve
    let a = Matrix::<2>::try_from_rows([[1.0, 2.0], [3.0, 4.0]])?;
    let b = Vector::<2>::try_new([5.0, 11.0])?;
    let x = a.solve_exact_f64(b)?.into_array();
    assert!((x[0] - 1.0).abs() <= f64::EPSILON);
    assert!((x[1] - 2.0).abs() <= f64::EPSILON);

    Ok(())
}
```

With the `exact` feature enabled, `BigInt` and `BigRational` are re-exported
from the crate root and prelude, alongside the most commonly needed
`num-traits` items (`FromPrimitive`, `ToPrimitive`, `Signed`). This lets
consumers construct exact values (`BigRational::from_f64`, `from_i64`), query
sign (`is_positive` / `is_negative`), and convert back to `f64` (`to_f64`)
with a single `use la_stack::prelude::*;` — no need to add `num-bigint`,
`num-rational`, or `num-traits` to their own `Cargo.toml`.

For `det_sign_exact()`, D ≤ 4 matrices use a fast f64 filter (error-bounded
`det_direct()`) that resolves the sign without allocating. Only near-degenerate
or large (D ≥ 5) matrices fall through to the exact Bareiss algorithm.

### Adaptive precision with `det_errbound()`

`det_errbound()` returns the conservative absolute error bound used by the fast
filter. This method does NOT require the `exact` feature — it uses pure f64
arithmetic and is available by default. This enables building custom
adaptive-precision logic for geometric predicates:

```rust,ignore
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    let m = Matrix::<3>::identity();
    if let Some(bound) = m.det_errbound()? {
        if let Some(det) = m.det_direct()? {
            if det.abs() > bound {
                // f64 sign is guaranteed correct
                let sign = det.signum() as i8;
            } else {
                // Fall back to exact arithmetic (requires `exact` feature)
                let sign = m.det_sign_exact()?;
            }
        }
    } else {
        // D ≥ 5: no fast filter, use exact directly (requires `exact` feature)
        let sign = m.det_sign_exact()?;
    }

    Ok(())
}
```

The error coefficients (`ERR_COEFF_2`, `ERR_COEFF_3`, `ERR_COEFF_4`) are the
dimension-specific constants behind that bound. In plain terms, they answer:
"how many machine-epsilon-sized rounding mistakes can this closed-form
determinant formula accumulate?" To get an absolute error bound, `det_errbound()`
multiplies the coefficient by a size measure of the matrix entries, the
**absolute Leibniz sum**:

```text
p(|A|) = sum over determinant terms of product of absolute values
```

For a 2×2 matrix `[[a, b], [c, d]]`, that scale is `|a*d| + |b*c|`, so:

```text
|det_direct(A) - det_exact(A)| <= ERR_COEFF_2 * (|a*d| + |b*c|)
```

The coefficients are not tolerances and are not meant to be tuned by callers;
they are conservative constants derived from the fixed D ≤ 4 formulas and their
floating-point rounding chains. They are exposed for advanced users who want to
compose the same bound themselves.

## 🧩 API at a glance

| Type | Storage | Purpose | Key methods |
|---|---|---|---|
| `Vector<D>` | `[f64; D]` | Finite fixed-length vector for input and computation | `try_new`, `zero`, `dot`, `norm2_sq` |
| `Matrix<D>` | `[[f64; D]; D]` | Finite square matrix for input and computation | See below |
| `Lu<D>` | `Matrix<D>` + pivot array | Factorization for solves/det | `solve`, `det` |
| `Ldlt<D>` | `Matrix<D>` | Factorization for symmetric SPD/PSD solves/det | `solve`, `det` |

Storage shown above reflects the intentional `f64` scalar model.

`Matrix<D>` key methods: `lu`, `ldlt`, `det`, `det_direct`, `det_errbound`,
`det_exact`¹, `det_exact_f64`¹, `det_exact_rounded_f64`¹, `det_sign_exact`¹,
`solve_exact`¹, `solve_exact_f64`¹, `solve_exact_rounded_f64`¹.
Matrix and vector constructors validate non-finite inputs at public API
boundaries. After construction, `Matrix<D>` and `Vector<D>` carry that
finite-storage invariant directly, so kernels do not revalidate stored entries.

¹ Requires `features = ["exact"]`.

## 📊 Benchmarks (vs nalgebra/faer)

![LU solve (factor + solve): median time vs dimension](docs/assets/bench/vs_linalg_lu_solve_median.svg)

Raw data: [docs/assets/bench/vs_linalg_lu_solve_median.csv](docs/assets/bench/vs_linalg_lu_solve_median.csv)

Representative benchmark: `lu_solve` factors the matrix and solves one
right-hand side. Median time is lower-is-better, and the “la-stack vs
nalgebra/faer” columns show the % time reduction relative to each baseline
(positive = la-stack faster). This is not an aggregate score across all
operations.

For the full per-kernel comparison methodology, input construction, and
release-comparison workflow details, see [docs/BENCHMARKING.md](docs/BENCHMARKING.md).

<!-- BENCH_TABLE:lu_solve:median:new:BEGIN -->

| D | la-stack median (ns) | nalgebra median (ns) | faer median (ns) | la-stack vs nalgebra | la-stack vs faer |
|---:|--------------------:|--------------------:|----------------:|---------------------:|----------------:|
| 2 | 2.585 | 4.486 | 137.653 | +42.4% | +98.1% |
| 3 | 12.204 | 22.990 | 182.618 | +46.9% | +93.3% |
| 4 | 27.228 | 51.660 | 208.181 | +47.3% | +86.9% |
| 5 | 53.141 | 68.714 | 272.117 | +22.7% | +80.5% |
| 8 | 141.279 | 162.225 | 348.216 | +12.9% | +59.4% |
| 16 | 626.561 | 574.115 | 854.941 | -9.1% | +26.7% |
| 32 | 2,862.795 | 2,709.532 | 2,806.698 | -5.7% | -2.0% |
| 64 | 19,703.239 | 14,388.285 | 12,085.453 | -36.9% | -63.0% |

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

```bash
cargo install just
just setup        # install/verify dev tools + sync Python deps
just check        # lint/validate (non-mutating)
just fix          # apply auto-fixes (mutating)
just ci           # lint + tests + examples + bench compile
```

The repository uses Rust-native tooling for documentation and config checks:
`rumdl` for Markdown, `dprint` with `pretty_yaml` for YAML, `taplo` for TOML,
and `typos` for spelling. GitHub Actions references are SHA-pinned, restricted
to an explicit allowlist, and kept with readable version comments for review.

CI runs `just ci` on Ubuntu, macOS, and Windows to keep platform coverage
aligned with the local comprehensive validation path.

For coverage commands and report locations, see [`docs/COVERAGE.md`](docs/COVERAGE.md).
For the full contributor workflow, see [CONTRIBUTING.md](CONTRIBUTING.md).

## 📝 Citation

If you use this library in academic work, please cite it using [CITATION.cff](CITATION.cff) (or GitHub's
"Cite this repository" feature). Tagged releases are archived on Zenodo.

## 📚 References

For canonical references to the algorithms used by this crate, see [REFERENCES.md](REFERENCES.md).

## 🤖 AI Agents

AI coding assistants should read [AGENTS.md](AGENTS.md) before proposing or
applying changes. See [CONTRIBUTING.md](CONTRIBUTING.md) for the repository's
AI-assisted development note.

## 📄 License

BSD 3-Clause License. See [LICENSE](./LICENSE).
