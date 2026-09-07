# Code Organization

Module, feature, and file-placement guidance for the single Rust library crate.
Read [AGENTS.md](../AGENTS.md) for agent rules and use the
[contributor workflow](../CONTRIBUTING.md) for setup and validation commands.

## Contents

- [Library modules](#library-modules)
- [Feature boundaries](#feature-boundaries)
- [Tests and examples](#tests-and-examples)
- [Benchmarks and support tooling](#benchmarks-and-support-tooling)
- [Documentation owners](#documentation-owners)

## Library modules

`src/lib.rs` wires private implementation modules and public re-exports. It
also owns the documentation-only `guide` module, the private README doctest
mirrors, and the public prelude. There is no `src/main.rs`.

| Module | Owns |
|--------|------|
| [`src/error.rs`](../src/error.rs) | `LaError` and typed singularity, non-finite, positive-semidefinite, tolerance, factorization, arithmetic-operation, and exact-conversion categories |
| [`src/exact.rs`](../src/exact.rs) | Exact determinants and solves, determinant-sign filtering, and exact-to-`f64` conversion |
| [`src/gram.rs`](../src/gram.rs) | Fixed-size Gram construction from vector dot products |
| [`src/interval.rs`](../src/interval.rs) | Interval scalars, matrices, and determinant-sign certification |
| [`src/ldlt.rs`](../src/ldlt.rs) | Unpivoted `Ldlt<D>` factorization for exactly symmetric positive-definite matrices, solves, and determinants |
| [`src/lu.rs`](../src/lu.rs) | Partially pivoted `Lu<D>` factorization, solves, and determinants |
| [`src/matrix.rs`](../src/matrix.rs) | `Matrix<D>`, finite storage, accessors, matrix operations, norms, and direct determinant APIs |
| [`src/norm.rs`](../src/norm.rs) | Exact range-boundary fallback for the approximate Euclidean norm |
| [`src/rational.rs`](../src/rational.rs) | `RationalMatrix<D>` and `RationalVector<D>`, canonical rational inputs, and denominator clearing for exact arithmetic |
| [`src/rounding.rs`](../src/rounding.rs) | Shared binary64 rounding primitives for certified arithmetic |
| [`src/scaled_product.rs`](../src/scaled_product.rs) | Allocation-free scaled products for floating-point factor diagonals |
| [`src/tolerance.rs`](../src/tolerance.rs) | Validated singular-tolerance policy |
| [`src/vector.rs`](../src/vector.rs) | `Vector<D>`, finite storage, reductions, norms, and certified scalar results |

Keep shared arithmetic invariants in the lowest owning module and expose public
items through `src/lib.rs`. Mathematical derivations and algorithm references
belong in [Mathematical basis](mathematical_basis.md) and
[REFERENCES.md](../REFERENCES.md).

## Feature boundaries

The [Cargo manifest](../Cargo.toml) owns feature and dependency declarations.

- **`exact`** gates `src/exact.rs`, `src/rational.rs`, their additional tests,
  exact guides, and exact-arithmetic examples. It enables `det_exact()`,
  strict `det_exact_f64()`, rounded `det_exact_rounded_f64()`,
  `det_sign_exact()`, `solve_exact()`, strict `solve_exact_f64()`, and
  rounded `solve_exact_rounded_f64()` for finite binary64 inputs, plus the
  rational-input types.
- `det_sign_exact()` is infallible for finite-by-construction `Matrix` values.
  Exact-value, conversion, and solve APIs retain their genuine scale,
  representation, and singularity failures.
- `ExactF64Conversion` converts an already-computed determinant or solution
  under the strict or rounded contract without rerunning exact elimination.
  Feature-gated re-exports include `DeterminantSign`, `ExactF64Conversion`,
  `RationalMatrix`, `RationalVector`, `BigInt`, `BigRational`, `FromPrimitive`,
  `ToPrimitive`, and `Signed`.
- Typed error categories, including `UnrepresentableReason`, remain available
  without `exact`; matching errors must not require optional arithmetic
  dependencies.
- The scaled `BigInt` determinant core uses direct expansions for D≤4 and
  fraction-free Bareiss elimination for D≥5. Exact signs add a floating-point
  filter for D≤4. Exact solves use fraction-free forward elimination with
  first-non-zero pivoting and `BigRational` back-substitution. Rational inputs
  clear denominators before reusing the integer backend.
- **`bench`** is a cfg-only gate for benchmark targets and
  `tests/vs_linalg_inputs.rs`; benchmark libraries remain dev-dependencies.

## Tests and examples

- Rust unit tests live in inline `#[cfg(test)]` modules in `src/*.rs`.
- `tests/proptest_*.rs` covers matrix, vector, factorization, exact, rational,
  interval, and Gram properties. Exact and rational suites require `exact`.
- Other `tests/*.rs` suites cover regressions, API contracts, allocation
  behavior, conversion boundaries, and benchmark inputs. Shared property-test
  configuration lives in `tests/common/`.
- `tests/semgrep/` contains deliberate static-analysis fixtures, not ordinary
  runnable API examples.
- Python tests live under `scripts/tests/` and run with pytest through
  `just test-python` or the aggregate validation workflows.
- `examples/` contains complete runnable workflows. README mirrors and public
  guide doctests live in `src/lib.rs`.

Use [Testing guidance](dev/testing.md) for dimension coverage and test design,
and [Documentation guidance](dev/docs.md) for executable-example ownership.

## Benchmarks and support tooling

`benches/` contains Criterion suites for exact arithmetic, Gram construction,
intervals, linear forms, and nalgebra/faer comparisons. Helpers under
`benches/common/` own fixture and oracle validation. Exact benchmark helpers
accept only `ValidatedExactInput`, after independent validation outside timing.
Adversarial groups include near-singular, large-entry, and Hilbert inputs.

The [Benchmarking guide](BENCHMARKING.md) owns benchmark commands, methodology,
baselines, output locations, and report promotion. The [Scripts guide](../scripts/README.md)
owns the Python script inventory and entry points for comparisons, plotting,
release metadata, changelog generation/archiving, and tag preparation.
The [justfile](../justfile) owns executable development workflows.

`scripts/release_baseline.py` owns release-suite inventory and complete raw
Criterion validation. The release workflow packages only datasets that pass
that gate; its regression and archive tests live in
`scripts/tests/test_release_baseline.py`.

## Documentation owners

Use [Documentation guidance](dev/docs.md) for README, references, mathematical
background, rustdoc, and generated-file ownership. Focused operational rules
live in `docs/dev/`; human setup and validation guidance lives in
[CONTRIBUTING.md](../CONTRIBUTING.md). Release procedures belong in
[Releasing](RELEASING.md).

[Managing changes](dev/MANAGING_CHANGES.md) owns Git and GitHub procedures.
[Measuring coverage](MEASURING_COVERAGE.md) owns local and CI coverage execution.
The generated [performance report](performance.md) records measured results;
[Benchmarking](BENCHMARKING.md) owns the commands that produce and compare them.

When adding, removing, renaming, or moving files, update the applicable ownership
rows here. Prefer links to the detailed owner over copying its procedure into
this map or `AGENTS.md`.
