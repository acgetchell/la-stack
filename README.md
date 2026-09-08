# la-stack

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18158926.svg)](https://doi.org/10.5281/zenodo.18158926)
[![Crates.io](https://badgen.net/crates/v/la-stack)](https://crates.io/crates/la-stack)
[![Downloads](https://badgen.net/crates/d/la-stack)](https://crates.io/crates/la-stack)
[![License](https://badgen.net/github/license/acgetchell/la-stack)](https://github.com/acgetchell/la-stack/blob/v0.4.6/LICENSE)
[![Docs.rs](https://docs.rs/la-stack/badge.svg)](https://docs.rs/la-stack)
[![CI](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml)
[![rust-clippy analyze][clippy-badge]][clippy-workflow]
[![codecov](https://codecov.io/gh/acgetchell/la-stack/graph/badge.svg?token=4eKXa5QjuZ)](https://codecov.io/gh/acgetchell/la-stack)
[![Audit dependencies][audit-badge]][audit-workflow]

![la-stack](https://raw.githubusercontent.com/acgetchell/la-stack/main/docs/assets/la-stack.jpg)

Fast, stack-allocated linear algebra for fixed dimensions in Rust.

This crate grew from the need to support [`delaunay`](https://crates.io/crates/delaunay) with fast, stack-allocated linear algebra primitives and algorithms
while keeping the API intentionally small and explicit.

## Contents

- [Introduction](#-introduction)
- [Use this crate when](#-use-this-crate-when)
- [Quickstart](#-quickstart)
- [Scalar and bounded-value types](#-scalar-and-bounded-value-types)
- [API at a glance](#-api-at-a-glance)
- [Features](#-features)
  - [Adaptive determinant filtering (D ≤ 4)](#adaptive-determinant-filtering-d--4)
  - [Certified dot products and affine differences](#certified-dot-products-and-affine-differences)
  - [Compile-time determinants (D ≤ 4)](#compile-time-determinants-d--4)
  - [Exact arithmetic](#exact-arithmetic-exact-feature)
  - [LDLT determinant](#ldlt-determinant)
  - [LU solve](#lu-solve)
  - [Outward-rounded interval determinants](#outward-rounded-interval-determinants)
  - [Overflow-safe Euclidean norms](#overflow-safe-euclidean-norms)
- [Mathematical basis](#-mathematical-basis)
- [Design goals](#-design-goals)
- [Anti-goals](#-anti-goals)
- [Documentation Map](#documentation-map)
- [Examples](#-examples)
- [Benchmarks](#-benchmarks-vs-nalgebrafaer)
- [Contributing](#-contributing)
- [Citation](#-citation)
- [References](#-references)
- [AI Agents](#-ai-agents)
- [License](#-license)

## 📐 Introduction

`la-stack` provides a handful of const-generic, stack-backed building blocks:

- `gram_matrix(&[Vector<N>; M])` for allocation-free `Matrix<M>` construction
  from pairwise vector inner products, with bit-for-bit symmetry. Gram matrices
  encode lengths and angles and support simplex/facet volume calculations; see
  [Gram matrices and geometric measures][refs-gram].
  Each independent dot product is checked once;
  rounding has no certified error bound, and positive definiteness or affine
  independence must still be established by factorization or the caller.
  Benchmark square simplex and rectangular facet inputs through dimension 8
  with `cargo bench --locked --features bench --bench gram`.
- `Interval` and `IntervalMatrix<const D: usize>` for outward-rounded,
  proof-bearing determinant filters through D=7
- `Ldlt<const D: usize>` for no-pivot factorization intended for exactly
  symmetric positive-definite matrices (solve + det; typed pivot diagnostics)
- `Lu<const D: usize>` for LU factorization with partial pivoting (solve + det)
- `Matrix<const D: usize>` for fixed-size square `f64` matrices backed by `[[f64; D]; D]`
- `RationalVector<const D: usize>` and `RationalMatrix<const D: usize>` for
  exact rational inputs behind the optional `"exact"` feature
- `ScalarWithErrorBound` for proof-bearing fixed-vector dot products and
  affine differences over finite `f64` inputs
- `Vector<const D: usize>` for fixed-length `f64` vectors backed by `[f64; D]`

## ✅ Use this crate when

- Robust predicates matter for geometry-style workloads near degeneracy
- Stack allocation and `Copy` value semantics fit your data flow
- You need a certified sign or threshold comparison for a fixed-vector dot
  product or `axis · (left - right)` expression
- You need a cheap, sound interval filter for determinant expressions assembled
  from rounded binary64 operations
- You need exact determinants, exact determinant signs, or exact linear solves
  for fixed-size systems
- You prefer a default build with no runtime dependencies
- You want explicit LU / LDLT / determinant APIs rather than a broad algebra toolkit
- Your matrices and vectors have small, fixed dimensions known at compile time

## 🚀 Quickstart

The minimum supported Rust version (MSRV) is 1.98.1.

Add this to your `Cargo.toml`:

```toml
[dependencies]
la-stack = "0.4.6"
```

### Solve a 5×5 system

This system has solution `[1, 2, 3, 4, 5]` and requires partial pivoting:

```rust
use la_stack::prelude::*;

fn main() -> Result<(), LaError> {
    // The zero leading entry requires LU pivoting.
    let a = Matrix::<5>::try_from_rows([
        [0.0, 2.0, -1.0, 1.0, 3.0],
        [4.0, -1.0, 2.0, 0.0, 1.0],
        [1.0, 3.0, 5.0, -2.0, 0.0],
        [2.0, 0.0, -1.0, 4.0, 1.0],
        [-1.0, 2.0, 0.0, 1.0, 6.0],
    ])?;
    let b = Vector::try_new([20.0, 13.0, 14.0, 20.0, 37.0])?;
    let lu = a.lu(DEFAULT_SINGULAR_TOL)?;
    let x = lu.solve(b)?;

    for (&actual, expected) in x.as_array().iter().zip([1.0, 2.0, 3.0, 4.0, 5.0]) {
        assert!((actual - expected).abs() <= 1e-12);
    }
    Ok(())
}
```

The assertion tolerance is suitable for this known example; LU does not
provide a certified solution error bound.

### Feature flags

- `bench`: repository-development gate used only by benchmark targets and
  benchmark-input tests; application crates should not enable it
- `default`: no runtime dependencies; includes outward-rounded `Interval` and
  `IntervalMatrix` APIs
- `exact`: exact determinant signs, determinant values, and solves over stored
  `f64` values or caller-supplied `BigRational` inputs

## 🔢 Scalar and bounded-value types

The public point-value scalar model deliberately has two input domains:

- arbitrary-precision `BigRational` through `RationalMatrix<D>` and
  `RationalVector<D>` behind the optional `"exact"` feature;
- finite `f64` through `Matrix<D>` and `Vector<D>` for floating-point work.

`Interval` is a separate bounded-value layer over finite `f64` endpoints. It
encloses exact-real values during a small set of outward-rounded operations and
feeds `IntervalMatrix<D>` determinant proofs; it does not make `Matrix` generic
over alternate scalars or provide a general interval package.

This is not a generic scalar-parameterized API. Exact support intentionally
covers the robustness-sensitive operations that require it: determinant sign,
determinant value, and linear solve, followed by explicit strict or rounded
conversion when an `f64` result is required. It does not promise a
`BigRational` counterpart for every floating-point helper or factorization.

Lower-precision `f32` / `f16` throughput-oriented workloads are outside the
crate's scope; they usually indicate large-matrix or accelerator-oriented use
cases better served by broader linear-algebra libraries.

## 🧩 API at a glance

Start with the capability you need; the [API reference][api-reference] lists
the complete public surface, and the [worked examples][api-guide] show how
to combine operations.

| Capability | Main entry points |
|---|---|
| Certified dot, affine-difference, and determinant estimates | [`ScalarWithErrorBound`][api-scalar-bound], [`DeterminantWithErrorBound`][api-det-bound] |
| Exact signs, determinants, solves, and output conversion¹ | [Exact arithmetic examples][api-exact] |
| Floating-point determinants and solves | [`Matrix<D>`][api-matrix], [`Lu<D>`][api-lu], [`Ldlt<D>`][api-ldlt] |
| Gram matrix construction | [`gram_matrix`][api-gram] |
| Interval expressions and determinant signs | [`Interval`][api-interval], [`IntervalMatrix<D>`][api-interval-matrix] |
| Runtime selection of a const-generic matrix dimension | [Dimension dispatch examples][api-dispatch] |
| Vector operations and norms | [`Vector<D>`][api-vector] |

[`Tolerance`][api-tolerance] validates numerical rejection thresholds.
[`LaError`][api-error] and its reason/location enums preserve structured
failure details; match non-exhaustive enums with a wildcard and struct-style
variants with `..`. See the [storage, access, and error guide][api-contracts]
for the full contracts.

¹ Requires `features = ["exact"]`.

[api-reference]: https://docs.rs/la-stack/latest/la_stack/
[api-guide]: https://docs.rs/la-stack/latest/la_stack/guide/index.html
[api-vector]: https://docs.rs/la-stack/latest/la_stack/struct.Vector.html
[api-matrix]: https://docs.rs/la-stack/latest/la_stack/struct.Matrix.html
[api-lu]: https://docs.rs/la-stack/latest/la_stack/struct.Lu.html
[api-ldlt]: https://docs.rs/la-stack/latest/la_stack/struct.Ldlt.html
[api-gram]: https://docs.rs/la-stack/latest/la_stack/fn.gram_matrix.html
[api-scalar-bound]: https://docs.rs/la-stack/latest/la_stack/struct.ScalarWithErrorBound.html
[api-det-bound]: https://docs.rs/la-stack/latest/la_stack/struct.DeterminantWithErrorBound.html
[api-interval]: https://docs.rs/la-stack/latest/la_stack/struct.Interval.html
[api-interval-matrix]: https://docs.rs/la-stack/latest/la_stack/struct.IntervalMatrix.html
[api-exact]: https://docs.rs/la-stack/latest/la_stack/guide/exact/index.html
[api-dispatch]: https://docs.rs/la-stack/latest/la_stack/guide/index.html#dimension-dispatch
[api-tolerance]: https://docs.rs/la-stack/latest/la_stack/struct.Tolerance.html
[api-error]: https://docs.rs/la-stack/latest/la_stack/enum.LaError.html
[api-contracts]: https://docs.rs/la-stack/latest/la_stack/guide/index.html#storage-access-and-errors

## ✨ Features

### Adaptive determinant filtering (D ≤ 4)

`det_direct_with_errbound()` pairs a determinant with its certified absolute
bound, without optional dependencies. Resolve the sign when `|det| > bound`;
otherwise an exact fallback is needed. With `exact`, `det_sign_exact()` handles
filtering and fallback automatically.
[Worked examples: the floating-point filter and exact fallback][guide-adaptive].

### Certified dot products and affine differences

`dot_with_errbound()` and `dot_difference_with_errbound()` return certified
bounds for dot products and `axis · (left - right)` over the original stored
coordinates. Their endpoints support sign and threshold proofs; a bound that
straddles the threshold or an unavailable certificate is inconclusive.
[Worked examples: dot-product signs and affine threshold tests][guide-certified].

### Compile-time determinants (D ≤ 4)

`det_direct()` evaluates closed-form determinants in `const` contexts through
D=4. `det()` selects those formulas automatically and uses zero-tolerance LU
for larger dimensions; a failed numerical pivot remains `LaError::Singular`.
[Compile-time example and dimension contracts][guide-compile-time].

### Exact arithmetic (`"exact"` feature)

Enable exact determinant signs, determinant values, and solves:

```toml
[dependencies]
la-stack = { version = "0.4.6", features = ["exact"] }
```

`Matrix` / `Vector` exact methods preserve stored `f64` values;
`RationalMatrix` / `RationalVector` also preserve rational expressions before
any `f64` rounding. Keep exact results or explicitly choose strict versus
rounded conversion with `ExactF64Conversion`.
[Worked examples: rational inputs, exact solves, and output conversion][api-exact].

### LDLT determinant

`Matrix::ldlt()` provides a square-root-free factorization for exactly symmetric
positive-definite matrices, supporting determinants and solves without pivoting.
Approximate symmetry is not sufficient, and floating-point success is not an
exact positive-definiteness certificate.
[Worked example and typed pivot diagnostics][guide-ldlt].

### LU solve

`Matrix::lu()` uses partial pivoting for general square systems. Reuse one
`Lu` factorization for multiple right-hand sides or a determinant; pivot
tolerances control rejection, not solution accuracy.
[Worked example: solving and reusing factors][guide-lu].

### Outward-rounded interval determinants

`Interval` preserves bounds while assembling differences, squares, and other
expressions. `IntervalMatrix::det_sign()` certifies determinant signs through
D=7: an enclosure separated from zero proves its sign, and `[0, 0]` proves
exact zero. Other overlaps with zero are inconclusive and may need exact fallback.
[Worked example: lifted coordinates, range errors, and fallback][guide-intervals].

### Overflow-safe Euclidean norms

`Vector::norm()` avoids unnecessary overflow and underflow from squaring
coordinates. `norm_squared()` computes the squared norm and can overflow even
when the norm is finite. Both results remain approximate, without a certified
error bound.
[Worked example and range contracts][guide-norms].

**v0.4.6 migration:** `Vector::norm2_sq()` is renamed to `Vector::norm_squared()`,
the unreleased `Vector::norm2()` API is named `Vector::norm()`, and
`Matrix::inf_norm()` is renamed to `Matrix::norm_inf()`. The old method names
are removed; their numerical behavior and error contracts are unchanged by
the renames. `Matrix::norm_inf()` remains the maximum absolute row sum.

[guide-lu]: https://docs.rs/la-stack/latest/la_stack/guide/index.html#solving-and-reusing-factors
[guide-ldlt]: https://docs.rs/la-stack/latest/la_stack/guide/ldlt/index.html
[guide-compile-time]: https://docs.rs/la-stack/latest/la_stack/guide/compile_time/index.html
[guide-intervals]: https://docs.rs/la-stack/latest/la_stack/guide/intervals/index.html
[guide-norms]: https://docs.rs/la-stack/latest/la_stack/guide/norms/index.html
[guide-certified]: https://docs.rs/la-stack/latest/la_stack/guide/certified/index.html
[guide-adaptive]: https://docs.rs/la-stack/latest/la_stack/guide/adaptive/index.html

## 🧮 Mathematical basis

`la-stack` operates on finite IEEE 754 binary64 values in small, fixed
dimensions. Its floating-point paths use [LU with partial pivoting][refs-lu],
[LDLT without pivoting][refs-ldlt] for exactly symmetric positive-definite matrices, and closed-form
determinants through D=4. These results remain subject to conditioning and
binary64 rounding;
factorization tolerances are rejection thresholds, not accuracy guarantees. For
D≤4, direct determinants can be paired with a
[conservative absolute roundoff bound][refs-det-bound] when its range
preconditions hold. [Fixed-vector dot products and direct affine differences][refs-reductions]
can likewise return a paired estimate and certified absolute
roundoff bound without enabling arbitrary-precision dependencies.

Derived binary64 expressions can instead be assembled with `Interval`
subtraction, addition, multiplication, negation, and square. The resulting
`IntervalMatrix<D>` [determinant sign][refs-interval] is certified through D=7 when its enclosure
separates zero; the singleton `[0, 0]` also certifies exact zero. Every other
overlap with zero is explicitly inconclusive. This default-feature surface is
distinct from arbitrary-precision exact arithmetic.

With `features = ["exact"]`, callers can either lift stored binary64 inputs
losslessly or supply already-exact rational inputs for
[exact determinant signs][refs-exact-sign], determinant values, and
[solves][refs-exact-solve]. Exactness over binary64 input starts at the
stored values and cannot recover information rounded away before construction.
See the
[mathematical basis](https://github.com/acgetchell/la-stack/blob/main/docs/mathematical_basis.md)
for the algorithms, validity boundaries, and supporting references.

## 🎯 Design goals

- ✅ `const fn` where possible (compile-time evaluation of determinants, dot products, etc.)
- ✅ Const-generic storage (no dynamically sized matrix or vector representation)
- ✅ `Copy` types where possible
- ✅ Defined binary64 arithmetic semantics: Rust's `f64::algebraic_*`
  operations are forbidden because their unspecified reassociation, precision,
  and special-value behavior is incompatible with the crate's error bounds,
  non-finite classification, exact fallbacks, and reproducibility contract;
  deliberate `f64::mul_add` remains allowed for its defined single-rounding
  semantics
- ✅ Error-bounded f64 dot, affine-difference, and determinant filtering plus
  optional exact signs (`dot_with_errbound`, `dot_difference_with_errbound`,
  `det_errbound`, `det_sign_exact`)
- ✅ Overflow- and underflow-safe Euclidean vector norms (`norm`)
- ✅ Outward-rounded interval expressions and division-free determinant signs
  through D=7, with explicit inconclusive evidence
- ✅ Exact determinant values and linear solves via optional arbitrary-precision
  arithmetic (`det_exact`, `solve_exact`, strict/rounded f64 conversions)
- ✅ Explicit algorithms (LU, solve, determinant)
- ✅ Inline, stack-backed storage for core types; optional arbitrary-precision
  exact values allocate as required
- ✅ No runtime dependencies by default (optional features may add deps)
- ✅ `unsafe` forbidden

See [CHANGELOG.md](https://github.com/acgetchell/la-stack/blob/v0.4.6/CHANGELOG.md)
for release history and
[docs/roadmap.md](https://github.com/acgetchell/la-stack/blob/v0.4.6/docs/roadmap.md)
for current release planning.

## 🚫 Anti-goals

- Alternate scalar families: `la-stack` deliberately supports finite `f64` and
  optional exact `BigRational` input domains, not `f32`, `f16`, complex, or
  generic scalar APIs
- Bare-metal performance: use [`blas`](https://crates.io/crates/blas) or
  [`lapack`](https://crates.io/crates/lapack) with a native backend selected
  through [`blas-src`](https://crates.io/crates/blas-src),
  [`lapack-src`](https://crates.io/crates/lapack-src), or
  [`openblas-src`](https://crates.io/crates/openblas-src)
- Broad general-purpose linear algebra: use [`nalgebra`](https://crates.io/crates/nalgebra)
- Large matrices/dimensions with parallelism: use [`faer`](https://crates.io/crates/faer)

<a name="documentation-map"></a>

## 🗺️ Documentation Map

- [API guide][api-guide] — worked examples, API selection, storage, and error contracts.
- [Benchmarking](https://github.com/acgetchell/la-stack/blob/v0.4.6/docs/BENCHMARKING.md) — benchmark suites, comparison workflows, and measurement methodology.
- [Coverage](https://github.com/acgetchell/la-stack/blob/main/docs/MEASURING_COVERAGE.md) — local and CI coverage commands and report locations.
- [Mathematical basis](https://github.com/acgetchell/la-stack/blob/main/docs/mathematical_basis.md) — algorithms, numerical guarantees, and limitations.
- [Performance reports](https://github.com/acgetchell/la-stack/blob/main/docs/performance.md) — release-to-release measurement results and provenance.
- [Releasing](https://github.com/acgetchell/la-stack/blob/v0.4.6/docs/RELEASING.md) — release preparation, validation, and publication.
- [Roadmap](https://github.com/acgetchell/la-stack/blob/v0.4.6/docs/roadmap.md) — release planning, future directions, and non-goals.

## 📋 Examples

The `examples/` directory contains small, runnable programs:

- **`const_det_4x4`** — compile-time 4×4 determinant via `det_direct()`
- **`det_5x5`** — determinant of a 5×5 matrix via LU
- **`exact_det_3x3`** — exact determinant value of a near-singular 3×3 matrix (requires `exact` feature)
- **`exact_sign_3x3`** — exact determinant sign of a near-singular 3×3 matrix (requires `exact` feature)
- **`exact_solve_3x3`** — exact solve of a near-singular 3×3 system vs f64 LU (requires `exact` feature)
- **`ldlt_solve_3x3`** — solve a 3×3 symmetric positive definite system via LDLT
- **`rational_input_5x5`** — exact rational solve of a 5×5 system that becomes singular as f64 (requires `exact` feature)
- **`solve_5x5`** — solve a 5×5 system via LU with partial pivoting

```bash
just examples
# or individually:
cargo run --example const_det_4x4
cargo run --example det_5x5
cargo run --features exact --example exact_det_3x3
cargo run --features exact --example exact_sign_3x3
cargo run --features exact --example exact_solve_3x3
cargo run --example ldlt_solve_3x3
cargo run --features exact --example rational_input_5x5
cargo run --example solve_5x5
```

## 📈 Benchmarks (vs nalgebra/faer)

![LU solve (factor + solve): median time vs dimension][lu-solve-benchmark]

Raw data:
[docs/assets/bench/vs_linalg_lu_solve_median.csv](https://github.com/acgetchell/la-stack/blob/v0.4.6/docs/assets/bench/vs_linalg_lu_solve_median.csv)
Measurement provenance:
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
This snapshot records the measured source state, available CPU model, operating system, Rust
toolchain, dependency lock and harness digests, Criterion command, and
correctness-gate result in the adjacent JSON sidecar. The publication workflow
requires complete canonical-dimension coverage and regenerates the CSV, SVG,
README table, and provenance together.

For the full per-kernel comparison methodology, algorithm citations, input
construction, and release-comparison workflow details, see
[docs/BENCHMARKING.md](https://github.com/acgetchell/la-stack/blob/v0.4.6/docs/BENCHMARKING.md).
For the current release-to-release performance snapshot, see
[docs/performance.md](https://github.com/acgetchell/la-stack/blob/main/docs/performance.md).
The exact release suite includes the already-exact rational-input groups for
D=2 through D=8. Those rows report `RationalMatrix::det_sign`, `det`, and
`solve` alongside straightforward `BigRational` Gaussian determinant and solve
references. Releases produced with the rational-input harness include Criterion
point estimates and confidence intervals for these rows; comparisons against a
pre-API baseline retain them as explicit current-only measurements.

The focused `interval` Criterion suite covers conclusive and inconclusive
relative-coordinate lifted determinant signs at D=4 and the maximum supported
D=7 workload. Run it with `just bench-interval`; fixture validation stays
outside the timed closures.

The focused `linear_form` Criterion suite compares plain and certified dot
products and covers both well-separated and inconclusive dot/affine-difference
filters at D=4. Run it with `just bench-linear-form`; exact small-integer fixture
expectations are validated outside the timed closures.

<!-- BENCH_TABLE:lu_solve:median:new:BEGIN -->

| D | la-stack median (ns) | nalgebra median (ns) | faer median (ns) | reduction vs nalgebra (point est.) | reduction vs faer (point est.) |
|---:|--------------------:|--------------------:|----------------:|---------------------:|----------------:|
| 2 | 2.021 | 4.479 | 182.532 | +54.9% | +98.9% |
| 3 | 9.977 | 22.795 | 217.567 | +56.2% | +95.4% |
| 4 | 22.101 | 51.739 | 241.525 | +57.3% | +90.8% |
| 5 | 46.070 | 69.003 | 326.848 | +33.2% | +85.9% |
| 8 | 128.364 | 165.933 | 401.400 | +22.6% | +68.0% |
| 16 | 643.532 | 569.446 | 903.140 | -13.0% | +28.7% |
| 32 | 2,755.773 | 2,710.235 | 2,919.290 | -1.7% | +5.6% |
| 64 | 17,518.634 | 14,313.979 | 12,108.783 | -22.4% | -44.7% |

<!-- BENCH_TABLE:lu_solve:median:new:END -->

## 🤝 Contributing

A short contributor workflow:

Install Rust 1.98.1 through [rustup](https://rustup.rs/), Git,
[GitHub CLI](https://cli.github.com/), Python 3.14,
[`uv` 0.12.5](https://docs.astral.sh/uv/), and `jq`. Then install the pinned
`just` release from its locked dependency graph:

```bash
cargo install --locked just --version 1.58.0
just setup        # install/verify dev tools + sync Python deps + build
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
[`docs/MEASURING_COVERAGE.md`](https://github.com/acgetchell/la-stack/blob/main/docs/MEASURING_COVERAGE.md).
For the full contributor workflow, see
[CONTRIBUTING.md](https://github.com/acgetchell/la-stack/blob/v0.4.6/CONTRIBUTING.md).

## 📚 Citation

If you use this library in academic work, please cite it using
[CITATION.cff](https://github.com/acgetchell/la-stack/blob/v0.4.6/CITATION.cff)
(or GitHub's "Cite this repository" feature). Tagged releases are archived on
Zenodo under the
[all-versions concept DOI](https://doi.org/10.5281/zenodo.18158926).

## 🔎 References

For canonical references to the algorithms used by this crate, see
[REFERENCES.md](https://github.com/acgetchell/la-stack/blob/v0.4.6/REFERENCES.md).

## 🤖 AI Agents

AI coding assistants should read
[AGENTS.md](https://github.com/acgetchell/la-stack/blob/v0.4.6/AGENTS.md)
before proposing or applying changes. See
[CONTRIBUTING.md](https://github.com/acgetchell/la-stack/blob/v0.4.6/CONTRIBUTING.md)
for the repository's AI-assisted development note.

## 📜 License

BSD 3-Clause License. See [LICENSE](https://github.com/acgetchell/la-stack/blob/v0.4.6/LICENSE).

[audit-badge]: https://github.com/acgetchell/la-stack/actions/workflows/audit.yml/badge.svg
[audit-workflow]: https://github.com/acgetchell/la-stack/actions/workflows/audit.yml
[benchmark-provenance]: https://github.com/acgetchell/la-stack/blob/v0.4.6/docs/assets/bench/vs_linalg_lu_solve_median.provenance.json
[clippy-badge]: https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml/badge.svg
[clippy-workflow]: https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml
[lu-solve-benchmark]: https://raw.githubusercontent.com/acgetchell/la-stack/v0.4.6/docs/assets/bench/vs_linalg_lu_solve_median.svg
[refs-det-bound]: https://github.com/acgetchell/la-stack/blob/main/REFERENCES.md#absolute-error-bound-for-closed-form-determinants
[refs-exact-sign]: https://github.com/acgetchell/la-stack/blob/main/REFERENCES.md#exact-determinant-sign-adaptive-precision-integer-arithmetic
[refs-exact-solve]: https://github.com/acgetchell/la-stack/blob/main/REFERENCES.md#exact-linear-system-solve-hybrid-bareiss--bigrational
[refs-gram]: https://github.com/acgetchell/la-stack/blob/main/REFERENCES.md#gram-matrices-and-geometric-measures
[refs-interval]: https://github.com/acgetchell/la-stack/blob/main/REFERENCES.md#outward-rounded-interval-determinant-sign
[refs-ldlt]: https://github.com/acgetchell/la-stack/blob/main/REFERENCES.md#ldlᵀ-factorization-exactly-symmetric-positive-definite-inputs
[refs-lu]: https://github.com/acgetchell/la-stack/blob/main/REFERENCES.md#lu-decomposition-gaussian-elimination-with-partial-pivoting
[refs-reductions]: https://github.com/acgetchell/la-stack/blob/main/REFERENCES.md#certified-fixed-vector-reductions
