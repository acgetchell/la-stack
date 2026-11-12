# Roadmap

This roadmap records directions for `la-stack`, not release promises. Scope
depends on mathematical correctness, API maturity, downstream need, and
validation quality. Use [open issues](https://github.com/acgetchell/la-stack/issues)
and [milestones](https://github.com/acgetchell/la-stack/milestones) for scheduled
work. [CHANGELOG.md](../CHANGELOG.md) and
[GitHub releases](https://github.com/acgetchell/la-stack/releases) own release
history and migration notes.

Until v1.0, API breaks are acceptable when they materially improve correctness,
performance, or orthogonality. Releases remain on stable Rust; the
[Cargo manifest](../Cargo.toml) owns the minimum supported version. Version
changes and downstream migration belong to dedicated release preparation.

## Scalar scope

`la-stack` supports `f64` for floating-point APIs and exact rationals behind the
optional `exact` feature. The crate targets small, fixed-size,
robustness-sensitive numerical and computational geometry workloads where
these representations fit the required guarantees.

Other scalar families are non-goals. Users needing reduced precision,
accelerators, large dynamic matrices, sparse matrices, or broad decomposition
coverage should use larger linear-algebra ecosystems such as `nalgebra` or
`faer`.

## Longer-term directions

- Grow exact arithmetic where it strengthens robustness-sensitive workflows
  without adding dependencies to the default build.
- Keep core operations focused on small fixed dimensions and stack allocation.
- Preserve `const fn` surfaces where they make mathematical evaluation
  compile-time checkable.
- Validate performance claims against comparable workloads within the crate's
  intended dimensional scope, using the [benchmarking policy](BENCHMARKING.md).

### Expanded const generics

[#123](https://github.com/acgetchell/la-stack/issues/123) tracks determinant API
work in the [expanded const generics milestone](https://github.com/acgetchell/la-stack/milestone/11).
Stable const generics already support the crate's fixed dimensions. Richer
dimension-specific contracts depend on the particular language capabilities
needed to express them; no stabilization date or release is assumed here.

Before scheduling this work, confirm that the required capabilities are stable,
evaluate the API and MSRV impact, and establish downstream migration and
benchmark evidence. Unify determinant internals only where correctness,
maintainability, and caller usability justify it. A redesign that breaks an API
published in v1.x requires a major release.

## Non-goals

- Alternate floating-point scalar families such as `f32` or `f16`.
- Dynamic or large matrix dimensions.
- Replacing full-featured linear-algebra libraries.
- Sparse matrices, parallel solvers, GPU acceleration, or out-of-core storage.
