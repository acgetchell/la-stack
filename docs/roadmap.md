# Roadmap

This roadmap records likely directions for the `la-stack` crate. It is not a
stability promise; release scope depends on mathematical correctness, API
maturity, downstream need, and validation quality. Concrete implementation work
is tracked as GitHub issues; keep this document focused on release direction,
ordering, and non-goals.

Pre-1.0 releases may still revise public APIs when that makes the crate more
correct or easier to use safely. The `v0.4.x` line should stay on stable Rust
and focus on documentation, downstream ergonomics, API-contract cleanup,
validation, invariant audits, and benchmark coverage. The `v0.5.0` line is
reserved for work that depends on stabilized `generic_const_exprs` or equivalent
const-generic expressiveness.

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

The `v0.4.2` milestone collects the work that can be done on today's stable
Rust while keeping the crate useful to downstream geometry crates. The GitHub
native Blocking / Is blocked by graph mirrors this order:

- [#100](https://github.com/acgetchell/la-stack/issues/100) - Clarify `f64`
  as the intended floating-point type in README and other documentation.
- [#109](https://github.com/acgetchell/la-stack/issues/109) - Add fallible
  runtime matrix dispatch and contextual index errors.
- [#83](https://github.com/acgetchell/la-stack/issues/83) - Change
  `Matrix::set` to return `Option<()>` instead of `bool`.
- [#94](https://github.com/acgetchell/la-stack/issues/94) - Validate
  tolerance arguments consistently across the crate.
- [#82](https://github.com/acgetchell/la-stack/issues/82) - Unify `det()`
  error behavior across all dimensions as far as stable Rust allows.
- [#120](https://github.com/acgetchell/la-stack/issues/120) - Run the
  parse-don't-validate `NonZero*` audit.
- [#126](https://github.com/acgetchell/la-stack/issues/126) - Add finite
  `Matrix` and `Vector` proof types.
- [#125](https://github.com/acgetchell/la-stack/issues/125) - Add a Semgrep
  guardrail against `unwrap` / `expect` in examples, benches, and doctests.
- [#98](https://github.com/acgetchell/la-stack/issues/98) - Add random-input
  percentile benchmarks to the exact arithmetic suite.

The broad shape is: document scope first, add downstream dispatch ergonomics,
clean up small API contracts, tighten validation, make reusable invariants
explicit in proof-carrying types, lock the public examples and benchmarks into
proper error handling, then finish with broader benchmark work.

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
