# Roadmap

This roadmap records likely directions for the `la-stack` crate. It is not a
stability promise; release scope depends on mathematical correctness, API
maturity, downstream need, and validation quality. Concrete implementation work
is tracked as GitHub issues; keep this document focused on release direction,
ordering, and non-goals.

Until v1.0, API breaks are acceptable when they improve correctness,
performance, or orthogonality. The `v0.4.x` line should stay on stable Rust and
focus on documentation, downstream ergonomics, API-contract cleanup, validation,
invariant audits, benchmark coverage, and release/tooling hardening. The
`v0.5.0` line is reserved for work that depends on stabilized
`generic_const_exprs` or equivalent const-generic expressiveness.

## Scalar Scope

`la-stack` intentionally supports `f64` for floating-point APIs and exact
rationals behind the optional `"exact"` feature. This is not a temporary
implementation detail: the crate targets small, fixed-size,
robustness-sensitive numerical and computational geometry workloads where
`f64` plus exact arithmetic is the right tradeoff.

Other scalar families are non-goals. Users who need lower-precision `f32` /
`f16` support, throughput-oriented reduced precision, accelerator-heavy
computation, large dynamic matrices, sparse matrices, or broad decomposition
coverage should use larger linear-algebra ecosystems such as `nalgebra` or
`faer`.

## Current Release Sequence

### v0.4.2 Stable Rust Cleanup

The `v0.4.2` milestone collected work that could be done on stable Rust while
keeping the crate useful to downstream geometry crates:

Completed foundation work:

- [#100](https://github.com/acgetchell/la-stack/issues/100) clarified `f64`
  as the intended floating-point type in README and other documentation.
- [#109](https://github.com/acgetchell/la-stack/issues/109) added fallible
  runtime matrix dispatch and contextual index errors.
- [#83](https://github.com/acgetchell/la-stack/issues/83) changed
  `Matrix::set` to return `Option<()>` instead of `bool`.
- [#94](https://github.com/acgetchell/la-stack/issues/94) made tolerance
  parsing explicit and consistent across the crate.
- [#82](https://github.com/acgetchell/la-stack/issues/82) unified `det()`
  error behavior across all dimensions as far as stable Rust allows.
- [#120](https://github.com/acgetchell/la-stack/issues/120) completed the
  parse-don't-validate `NonZero*` audit.
- [#111](https://github.com/acgetchell/la-stack/issues/111),
  [#112](https://github.com/acgetchell/la-stack/issues/112),
  [#113](https://github.com/acgetchell/la-stack/issues/113), and
  [#117](https://github.com/acgetchell/la-stack/issues/117) cleaned up
  Markdown/YAML tooling, CI speed, and shared Rust workflow/security checks.

API-invariant cleanup:

- [#126](https://github.com/acgetchell/la-stack/issues/126) is resolved as an
  internal parse-don't-validate design rather than a public proof-bearing API.
  `Matrix<D>` and `Vector<D>` parse raw `f64` storage at construction and then
  carry the finite-entry proof directly. Crate-private finite wrappers remain
  only as implementation helpers at algorithm boundaries.
- The public prelude stays focused on downstream composition: raw boundary
  types, factorization handles, tolerances, crate errors, dispatch helpers, and
  documented constants. Proof-bearing wrappers remain crate-private, and
  exact-arithmetic integer/rational re-exports remain gated behind the `"exact"`
  feature.
- The public LDLT API remains `Matrix::ldlt`. Symmetry proof storage is kept
  internal, `SymmetricMatrix` is not exported, asymmetric inputs return
  `LaError::Asymmetric`, and negative LDLT pivots return
  `LaError::NotPositiveSemidefinite` rather than being folded into
  `LaError::Singular`.
- The determinant error-bound constants `ERR_COEFF_2`, `ERR_COEFF_3`, and
  `ERR_COEFF_4` are documented as dimension-specific roundoff multipliers over
  the absolute Leibniz sum, not caller-tuned tolerances.

Final release blockers:

- [#125](https://github.com/acgetchell/la-stack/issues/125) - Add a Semgrep
  guardrail against `unwrap` / `expect` in examples, benches, and doctests.
- [#98](https://github.com/acgetchell/la-stack/issues/98) - Add random-input
  percentile benchmarks to the exact arithmetic suite.

The broad shape is now: document scalar scope, add downstream dispatch
ergonomics, clean up small API contracts, tighten validation, encode reusable
invariants behind the public raw-boundary API, lock examples and benchmarks into
proper error handling, then finish with broader benchmark work.

### v0.4.3 Benchmark and Tooling Hardening

Before the const-generic API revision, tighten the benchmark and tooling story so
performance claims are auditable across releases and Python support scripts have
a modern typed baseline.

- [#137](https://github.com/acgetchell/la-stack/issues/137) - Investigate
  checked vector kernel performance for v0.4.3.
- [#138](https://github.com/acgetchell/la-stack/issues/138) - Add first-class
  cross-release performance comparison tooling.
- [#142](https://github.com/acgetchell/la-stack/issues/142) - Update Python
  tooling to 3.13 and parse scripts at boundaries.

### v0.5.0 Generic Const Expressions

`v0.5.0` is reserved for the post-stabilization const-generic API revision.
The current anchor issue is
[#123](https://github.com/acgetchell/la-stack/issues/123), which tracks
refactoring determinant APIs around stable `generic_const_exprs`.

That work should revisit determinant APIs after the language can express more
dimension-dependent support at the type level. In particular, it should look for
ways to reduce runtime `D` branching, make closed-form determinant support more
explicit, and unify direct, LU-backed, error-bound, and exact determinant paths
where doing so improves correctness and maintainability.

## Longer-Term Directions

- Keep the crate focused on fixed small dimensions and stack allocation.
- Preserve `const fn` surfaces where they make mathematical evaluation
  compile-time checkable.
- Grow exact-arithmetic support only where it strengthens robustness-sensitive
  workflows without making the default build dependency-heavy.
- Use benchmarks to validate performance claims against `nalgebra` and `faer`
  within the crate's intended dimensional scope.

## Non-Goals

- Dynamic or large matrix dimensions.
- Sparse matrices, parallel solvers, GPU acceleration, or out-of-core storage.
- Alternate floating-point scalar families such as `f32` or `f16`.
- Replacing full-featured linear-algebra libraries.
