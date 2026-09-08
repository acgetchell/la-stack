# AGENTS.md

Essential guidance for AI assistants working in this repository. Keep this file
as the entry point; detailed rules belong in the focused documents below.

## Contents

- [Required Reading](#required-reading)
- [Priorities](#priorities)
- [Core Rules](#core-rules)
- [Scientific Invariants](#scientific-invariants)
- [Rust and API Design](#rust-and-api-design)
- [Testing and Performance](#testing-and-performance)
- [Documentation](#documentation)
- [Validation](#validation)
- [Python](#python)
- [Project Context](#project-context)
- [Agent Expectations](#agent-expectations)

## Required Reading

Read this file before making changes. Then read the focused guidance matching
the task; those rules remain binding when applicable. Load additional owners
if the scope expands. Unrelated guides need not be loaded for every task.

| Task | Read |
|------|------|
| GitHub issues, branches, or commit messages | [Git and GitHub guidance](docs/dev/MANAGING_CHANGES.md) |
| Tests, doctests, or benchmark design | [Testing guidance](docs/dev/testing.md) |
| README, rustdoc, citations, or documentation layout | [Documentation guidance](docs/dev/docs.md) |
| Module ownership, features, or file organization | [Code organization](docs/code_organization.md) |
| Setup, commands, or validation selection | [Contributor workflow](CONTRIBUTING.md#validation-workflow), [justfile](justfile), `just --list` |
| Benchmark execution or performance reports | [Benchmarking](docs/BENCHMARKING.md) |
| Python support scripts | [Scripts guide](scripts/README.md) |
| Maintainer-requested release work | [Releasing](docs/RELEASING.md) |

## Priorities

This is a scientific linear-algebra library. Design decisions prioritize:

1. Mathematical correctness and invariant preservation.
2. API stability.
3. Composability.
4. Idiomatic, well-tested Rust.
5. Performance within the documented scope.

Favor the invariant over a convenient edit or faster implementation.

## Core Rules

### Git and Editing

- **Never mutate version-control state.** Do not run `git commit`, `git push`,
  `git tag`, or other ref/index-mutating commands. Suggest those commands for
  the user to run manually.
- Use `git --no-pager` for read-only Git commands, including status, diff, log,
  show, and blame.
- Preserve user changes. The worktree may be dirty; work around overlapping
  edits and do not revert unrelated work.
- **Use the structured patch/edit tool for manual edits.** Never use `sed`,
  `awk`, Python, or `perl` to write files. Shell text tools may be used for
  read-only inspection. Keep patch rejects and backup files outside the
  repository and clean them immediately.
- Repository formatters and linters are allowed, including Cargo, Taplo,
  Ruff, rumdl, dprint, typos, and actionlint commands.

### Safety and Releases

- **Unsafe Rust is forbidden.** Preserve the manifest-level
  `unsafe_code = "forbid"` lint and crate/module `#![forbid(unsafe_code)]`.
- **Dead code is forbidden.** Remove unused items; never add
  `#[allow(dead_code)]` or `#![allow(dead_code)]`, including in tests and
  static-analysis fixtures.
- Version bumps are maintainer-driven release work. Do not update package
  versions, lockfile versions, or dependency snippets during ordinary work;
  change them only when explicitly asked to do release/version-bump work.
- Never hand-edit `CHANGELOG.md`; use its documented generation workflow.
  Release procedures and generated-artifact ownership live in
  [Releasing](docs/RELEASING.md).

## Scientific Invariants

- Arbitrary-precision paths (`det_exact`, `solve_exact`) never silently lose
  precision. Strict exact-to-`f64` methods return
  `LaError::Unrepresentable` instead of rounding:
  `UnrepresentableReason::RequiresRounding` means a finite result requires
  rounding; `NotFinite` means no finite `f64` can represent the result.
  Explicit `*_exact_rounded_f64` methods opt into rounding but still return
  `NotFinite` when rounding cannot produce a finite value.
- New or changed `f64` operations that accumulate rounding error document
  an absolute bound (`det_errbound`, `ERR_COEFF_*`) or explicitly state
  that no bound is provided.
- Non-finite inputs and intermediates surface as `LaError::NonFinite`
  with typed `NonFiniteOrigin` and `NonFiniteLocation` metadata. Computed
  failures identify their `ArithmeticOperation`; raw matrix/vector values
  retain input origins and exact source locations. Never silently propagate
  NaN or use `unwrap_or(f64::NAN)`.
- Keep exact singularity distinct from tolerance-based rejection through
  `SingularityReason`. Numerical failures retain the factorization,
  observed pivot magnitude, and tolerance.
- Parse raw tolerances through `Tolerance::try_new`; failures retain typed
  `InvalidToleranceReason`. Exact-to-`f64` failures remain
  `LaError::Unrepresentable`.
- `Matrix::det()` uses closed forms through D=4. Its D≥5 zero-tolerance LU
  fallback preserves `LaError::Singular` when elimination cannot produce a
  non-zero pivot. Do not relabel that numerical failure as an exact `0.0`;
  use exact determinant APIs for exact singularity classification.
- Algorithms cite their sources through `REFERENCES.md` and document
  conditioning behavior. Mathematical explanations belong in
  [Mathematical basis](docs/mathematical_basis.md).

## Rust and API Design

- Keep the MSRV and contributor/CI toolchain aligned across `Cargo.toml`,
  `rust-toolchain.toml`, and `clippy.toml`; the current baseline is 1.98.1.
- Rust's `f64::algebraic_*` operations are forbidden in all repository-owned
  Rust, including tests, examples, and benchmarks, except intentional
  Semgrep fixtures under `tests/semgrep/` whose `f64::algebraic_*` usage is
  validated by `just semgrep-test`. Their unspecified
  reassociation, precision, and special-value behavior can invalidate error
  bounds, exact fallbacks, non-finite classification, and reproducibility.
  Ordinary operators and deliberate `f64::mul_add` remain allowed. Any
  fast-math design requires a separate issue, opt-in contract, correctness
  analysis, and benchmark evidence.
- Use const-generic dimensions for core types; do not introduce runtime
  dimensions. Prefer stack storage and allocation-free paths. Heap use
  belongs behind a feature flag or where exact arithmetic requires it.
- Keep optional dependencies feature-isolated and default builds minimal.
- Prefer `const fn` wherever possible; compile-time evaluation reinforces
  pure operations.
- Public error enums and struct-style error variants are
  `#[non_exhaustive]`; downstream matches include wildcard arms and `..`.
  Public wrapper types are `#[must_use]`.
- New functionality is additive by default, with ergonomic prelude
  re-exports. Pre-1.0 breaks are allowed when they materially improve
  correctness, orthogonality, performance, or long-term clarity. Do not
  retain compatibility aliases that weaken the model; document intentional
  breaks in commit messages and release notes.
- Use `Result<_, LaError>` for fallible operations. Public library code
  must not panic on user input. Plain-value APIs must be infallible for all
  representable inputs; use `Result` or `Option` for observable failure
  instead of `panic!`, `assert!`, `unwrap`, or `expect`.
- Panics are reserved for unreachable internal invariant violations and
  must be documented when callers could observe them.
- Validate at the lowest owning type/module. Higher layers preserve and
  propagate typed errors rather than stringifying them.
- Borrow by default: accept `&T`, `&mut T`, or `&[T]` and return
  borrowed views where possible. Take ownership or allocate only when needed.
- Use textbook names such as `Matrix`, `Vector`, `Lu`, `Ldlt`,
  `solve`, `det`, and `norm_inf`.

## Testing and Performance

- Test known values, typed error paths, and dimension-generic behavior
  across D=2 through D=5 wherever possible.
- Match exact error variants, reasons, origins, locations, and structured
  fields. Do not substitute numeric sentinels or assert only `is_err()`.
- Property tests verify algebraic invariants; adversarial inputs accompany
  well-conditioned cases. Fast-filter/exact-fallback pairs must agree on the
  domain where both are defined.
- Performance stays within the crate's small, fixed-dimension scope:
  stack storage and closed forms where available. Large/dynamic dimensions,
  sparse matrices, and parallelism belong to other libraries.
- Within that scope, prefer allocation-free kernels and deliberate FMA where
  appropriate to the numerical contract.
- Performance-sensitive changes require comparable before-and-after
  measurements with the same command, inputs, features, and environment.
  Preserve provenance and distinguish descriptive ratios from supported
  performance claims. Invariant-violating runs are invalid evidence.
- Read [Testing guidance](docs/dev/testing.md) for dimension macros,
  Criterion measurement rules, and independent benchmark validation.

## Documentation

- Under `docs/`, use uppercase verb/gerund filenames for task guides with
  execution instructions (`RELEASING.md`, `BENCHMARKING.md`); prefer gerunds.
  Use lowercase descriptive filenames for invariants, principles, reference
  material, and reports. Classify by primary purpose; see the
  [filename conventions](docs/dev/docs.md#filename-conventions) for details.
- `README.md` owns concise orientation, quickstart, and API navigation;
  `REFERENCES.md` owns bibliographic provenance; `docs/mathematical_basis.md`
  owns mathematical explanations, assumptions, guarantees, and derivations.
  Link between owners instead of duplicating detail.
- Place "Use this crate when" immediately after Introduction. Maintain
  Contents navigation and put API selection and geometry scope near the top
  of the mathematical basis.
- Sort independent algorithm discussions and unordered capability bullets
  lexicographically within coherent groups; preserve prerequisite and
  procedural order. Keep the bibliography thematic and citation identifiers
  and deep links stable. Link claims to specific references.
- Read [Documentation guidance](docs/dev/docs.md) for README inclusion in
  rustdoc, guide placement, feature-gated doctests, link destinations,
  scientific notation, and generated-file rules.

## Validation

- Use `just check` during iterative review and fixes. Reserve `just ci` for
  final validation once those iterations are complete; core Rust or public
  behavior changes require that final comprehensive pass.
- Select additional focused validators proportionally to the changed surfaces.
  Compose each relevant focused validator once for mixed changes.
- Use [Contributor validation guidance](CONTRIBUTING.md#validation-workflow)
  for the surface-to-command mapping. The [justfile](justfile) and
  `just --list` own the full command catalog.
- Run `just spell-check` after editing. Add legitimate technical terms to
  `typos.toml` under `[default.extend-words]`.
- For Markdown changes, run `just markdown-fix` and `just markdown-ci`
  (Markdown and spelling checks). Changed guide docs also require
  `just doc-check`; changed executable examples require the default and
  matching feature doctests.

## Python

- Support tooling targets Python 3.14. Use `uv run --locked` for all
  Python scripts; never invoke `python` or `python3` directly.
- Use pytest, not unittest, and add type hints to new code.
- `just python-check` includes blocking type checking; all code must pass.
- Read the [Scripts guide](scripts/README.md) for maintenance rules and
  entry-point ownership.

## Project Context

The published Rust library is rooted at `src/lib.rs`. The unpublished
`benches/comparison` workspace member owns nalgebra/faer benchmark dependencies
and input tests. Core dimensions are compile-time constants. The `exact` feature
enables arbitrary-precision arithmetic; `bench` gates benchmark targets while
their dependencies remain dev-only. The [Code organization guide](docs/code_organization.md)
maps modules, features, tests, and tooling. Update that guide when file ownership
or layout changes.

## Agent Expectations

- Prefer small, focused patches and the simplest maintainable correct solution.
- Search existing documentation and nearby code before inventing conventions.
- Fix small, clearly related issues in touched areas when doing so improves
  correctness, clarity, tests, or maintainability.
- Avoid broad mechanical churn; separate repository-wide cleanup from focused
  work.
