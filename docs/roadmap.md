# Roadmap

This roadmap records likely directions for the `la-stack` crate. It is not a
stability promise; release scope depends on mathematical correctness, API
maturity, downstream need, and validation quality. Concrete implementation work
is tracked as GitHub issues; keep this document focused on release direction,
ordering, and non-goals.

Until v1.0, API breaks are acceptable when they improve correctness,
performance, or orthogonality. The `v0.4.x` line stays on stable Rust and includes
focused robustness APIs, downstream ergonomics, API-contract cleanup,
documentation, validation, invariant audits, benchmark coverage, and
release/tooling hardening. The current release is
[`v0.5.0`](https://github.com/acgetchell/la-stack/releases/tag/v0.5.0), published
on 2026-10-17, consolidating the work previously targeted at `v0.4.7`.
It remains on stable Rust and makes the Rust 1.99.1 minimum an
explicit downstream upgrade boundary.
Expanded const-generic work is tracked separately without an assigned release;
it depends on the specific language capabilities required by the proposed API.

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

## Release History and Next Sequence

### v0.4.2 Stable Rust Cleanup (released)

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
  common defaults. Advanced determinant error-bound coefficients remain
  explicit crate-root exports. Proof-bearing wrappers remain crate-private,
  while the exact-arithmetic sign and integer/rational re-exports remain gated
  behind the `"exact"` feature.
- Fallible raw boundaries advertise parsing explicitly: `Tolerance::try_new`
  constructs validated tolerances, `Matrix::try_get` preserves index context,
  and `Matrix::set` validates atomically before mutation.
- Exact determinant signs use `DeterminantSign` rather than an invalidable raw
  integer. Because `Matrix` already carries the finite-entry proof and filter
  range failures fall back to exact integer arithmetic, `det_sign_exact` returns
  `DeterminantSign` infallibly; `as_i8` is reserved for numeric interoperability.
- The public LDLT API remains `Matrix::ldlt`. Symmetry proof storage is kept
  internal and `SymmetricMatrix` is not exported. Asymmetric errors retain both
  observed entries and their effective bound; negative pivots and zero pivots
  with nonzero coupling use distinct `PositiveSemidefiniteViolation` values.
- The public error model keeps semantic categories typed: exact and numerical
  singularity use distinct `SingularityReason` values, numerical rejection
  retains its factorization/pivot/tolerance, and `NonFiniteOrigin` plus
  `NonFiniteLocation` distinguish invalid inputs from arithmetic overflow.
- Strict exact-to-`f64` conversion reports `RequiresRounding` only when an
  explicit rounded fallback can return a finite value; values above the
  overflow-rounding midpoint report `NotFinite` consistently from strict and
  rounded APIs.
- The determinant error-bound constants `ERR_COEFF_2`, `ERR_COEFF_3`, and
  `ERR_COEFF_4` are documented as dimension-specific roundoff multipliers over
  the absolute Leibniz sum, not caller-tuned tolerances.

Completed release blockers:

- [#125](https://github.com/acgetchell/la-stack/issues/125) - Add a Semgrep
  guardrail against `unwrap` / `expect` in examples, benches, and doctests.
- [#98](https://github.com/acgetchell/la-stack/issues/98) - Add random-input
  percentile benchmarks to the exact arithmetic suite.

The broad shape is now: document scalar scope, add downstream dispatch
ergonomics, clean up small API contracts, tighten validation, encode reusable
invariants behind the public raw-boundary API, lock examples and benchmarks into
proper error handling, then finish with broader benchmark work.

### v0.4.3 Benchmark and Tooling Hardening (released)

This release tightened the benchmark and tooling story so performance claims are
auditable across releases and Python support scripts have a modern typed baseline.

- [#137](https://github.com/acgetchell/la-stack/issues/137) - Investigate
  checked vector kernel performance for v0.4.3.
- [#138](https://github.com/acgetchell/la-stack/issues/138) - Add first-class
  cross-release performance comparison tooling.
- [#142](https://github.com/acgetchell/la-stack/issues/142) - Update Python
  tooling to 3.13 and parse scripts at boundaries.

Release outcome:

- The release-signal comparison against `v0.4.2` showed broad
  improvement across LU, solve, determinant-via-LU, and vector helper rows.
- Treat the remaining `D=4` direct determinant regression as a tracked
  performance note rather than a release blocker because the LU-backed
  determinant and solve paths improved.
- Defer `Matrix::inf_norm` optimization to follow-up work after `v0.4.3`.
  Larger-dimension `vs_linalg` measurements suggest it is the most interesting
  leaf-kernel target, but it is not required for the release.

### v0.4.4 Focused Leaf-Kernel Performance (released)

This release used the improved benchmark workflow to investigate narrow
leaf-kernel performance gaps without broadening the crate's scope or weakening
the small fixed-dimension API model.

- [#154](https://github.com/acgetchell/la-stack/issues/154) - Investigate
  `Matrix::inf_norm` performance against `nalgebra` and `faer`.
- [#155](https://github.com/acgetchell/la-stack/issues/155) - Investigate
  `Vector::dot` and `Vector::norm2_sq` performance against `nalgebra` and
  `faer`.
- [#174](https://github.com/acgetchell/la-stack/issues/174) - Raise the stable
  Rust baseline to 1.97.0 and audit its new bit helpers. Retain the existing
  primitive and `num-bigint` operations where the new helpers do not simplify
  current hot paths or preserve benchmark performance.

Release outcome:

- `Matrix::inf_norm` moved finiteness checks off the ordinary per-cell success
  path while replaying only overflowed rows to preserve exact error locations.
- `Vector::dot` and `Vector::norm2_sq` adopted the same success-path reduction
  strategy while retaining left-to-right fused accumulation, `const fn`
  evaluation, and typed failure metadata.
- The MSRV moved to Rust 1.97.0 after auditing the new integer bit helpers; the
  existing exact-arithmetic operations remained where alternatives did not
  improve clarity or preserve measured performance.
- Direct and exact determinant hot paths were restored without weakening the
  numerical contracts or the fixed-size allocation model.

### v0.4.5 Rust 1.98 Numerical Policy (released)

This milestone completed stable-Rust maintenance without broadening the crate's
scalar or algorithm scope.

- [#208](https://github.com/acgetchell/la-stack/issues/208) raised the MSRV and
  pinned contributor/CI toolchain to Rust 1.98, audited the final stable release,
  and added a repository guard against the new algebraic floating-point
  operations in correctness-sensitive source, examples, and benchmarks.
- Local and release performance workflows were unified around retained,
  schema-versioned CSV/JSON inputs, validated rerendering, and transactional
  report promotion.
- Dependency, contributor-tool, and GitHub Action maintenance was refreshed
  while preserving explicit repository ownership of update scope.

The existing IEEE 754 operations, deterministic accumulation order, error
bounds, exact fallbacks, and typed non-finite behavior remain authoritative.
Any future approximate or fast-math API requires a separate opt-in design and
correctness analysis.

### v0.4.6 Exact Inputs and Certified Bounds (released)

This milestone extends robustness-sensitive workflows on stable Rust while
preserving fixed dimensions, `f64` floating-point storage, and feature-isolated
exact arithmetic.

Completed capability work:

- Certified dot-product and affine-difference bounds
  ([#220](https://github.com/acgetchell/la-stack/issues/220)) support sign and
  threshold decisions over stored coordinates, with explicit inconclusive
  outcomes when a certificate is unavailable or overlaps the threshold.
- Exact rational matrix and vector inputs
  ([#216](https://github.com/acgetchell/la-stack/issues/216)) preserve rational
  expressions before `f64` rounding and provide exact determinants and solves
  with explicit strict or rounded output conversion.
- Gram matrix construction
  ([#219](https://github.com/acgetchell/la-stack/issues/219)) accepts fixed
  vectors with independent const-generic input and output dimensions. The
  floating-point result is symmetric but is not a definiteness certificate.
- Interval expressions and determinant signs
  ([#218](https://github.com/acgetchell/la-stack/issues/218)) preserve outward
  bounds through determinant evaluation up to D=7, distinguish exact zero from
  inconclusive overlap, and report range exhaustion explicitly.
- Overflow-safe Euclidean norms
  ([#217](https://github.com/acgetchell/la-stack/issues/217)) avoid unnecessary
  overflow and underflow from squaring coordinates. Norm results remain
  approximate, without a certified error bound.

Performance and tooling outcomes:

- Benchmark release budgets and baseline validation
  ([#224](https://github.com/acgetchell/la-stack/issues/224)), attachment before
  immutable publication ([#225](https://github.com/acgetchell/la-stack/issues/225)),
  and cache isolation ([#226](https://github.com/acgetchell/la-stack/issues/226))
  now form one documented release workflow.
- Exact rational conversion and dense 4D determinant paths avoid redundant
  arithmetic while preserving exact results and typed conversion errors.
- Rational row-clearing measurements
  ([#233](https://github.com/acgetchell/la-stack/issues/233)) support retaining
  the existing borrowed component access; see the
  [allocation study](archive/performance/studies/rational-row-clearing.md). LU/LDLT solve
  finalization remains unchanged because its proposed replacement showed no
  repeatable speedup ([#234](https://github.com/acgetchell/la-stack/issues/234));
  see the [decision report](archive/performance/studies/solve-finalization.md).
- Tooling uses Rust 1.98.1 and Python 3.14, reconciles managed tool pins, rejects
  ambiguous updater output, and verifies the installed Cargo upgrade tool
  through its supported version command during setup.

Intentional compatibility changes include requiring Rust 1.98.1,
`Vector::norm2_sq()` becoming `Vector::norm_squared()`,
`Matrix::inf_norm()` becoming `Matrix::norm_inf()`,
and `Matrix::solve_exact()` returning `RationalVector<D>` instead of a raw
rational array. `Vector::norm()` names the new Euclidean norm; raw exact solution
storage remains available through `as_array()` and `into_array()`.

Further stable-Rust maintenance can remain in `v0.4.x`; it does not depend on
expanded const-generic support.

### v0.5.0 Stable Rust APIs and Performance (released)

The [v0.5.0 milestone](https://github.com/acgetchell/la-stack/milestone/3)
replaced the planned `v0.4.7` release. It combined the completed Rust 1.99
migration and stable vector-angle API with correctness-preserving performance
improvements and release-tooling adoption. Expanded const-generics stabilization
did not gate this release.

Completed work included in this release:

- Certified dot-product and affine-difference reductions
  ([#248](https://github.com/acgetchell/la-stack/issues/248)) reduce avoidable
  arithmetic and certificate overhead while preserving outward bounds, typed
  failures, and exact fallback.
- Interval arithmetic and determinant signs
  ([#247](https://github.com/acgetchell/la-stack/issues/247)) avoid redundant work
  while retaining containment guarantees and const evaluation.
- Rust 1.99.1 and LLVM 23
  ([#251](https://github.com/acgetchell/la-stack/issues/251)) establish the new
  MSRV and strict floating-point baseline, with retained correctness, codegen,
  and performance evidence.
- Shared performance orchestration
  ([#268](https://github.com/acgetchell/la-stack/issues/268)), complete-run reports
  ([#277](https://github.com/acgetchell/la-stack/pull/277)), and opt-in review
  recipes ([#253](https://github.com/acgetchell/la-stack/issues/253)) reduce local
  tooling duplication while preserving scientific policy and historical evidence.
- Stable unsigned vector angles
  ([#249](https://github.com/acgetchell/la-stack/issues/249)) support fixed vectors
  and allocation-free borrowed coordinate slices, including extreme magnitudes
  and representable subnormal angles. Results remain rounded without a certified
  error bound.

The bounded algebraic floating-point evaluation
([#250](https://github.com/acgetchell/la-stack/issues/250)) concluded with a
**retain-strict** decision on 2026-10-15. The selected matrix-norm reduction
improved local timing at D=4/5 but hid a required overflow error. The maintainer
rejected adoption and ended the remaining prototype investigation; no algebraic
API or feature is planned. The
[validation record](validation/algebraic-floating-point.md) documents the
decision, executed checks, and unperformed native Linux/Windows runs; the
[study](archive/performance/studies/algebraic-float.md) retains the raw evidence.

The release sequence is complete:

1. The triangular rational-solve optimization
   ([#246](https://github.com/acgetchell/la-stack/issues/246),
   [PR #278](https://github.com/acgetchell/la-stack/pull/278)) landed with the
   exact general-system fallback and independent regression evidence.
2. The published shared-tooling capabilities were adopted
   ([#276](https://github.com/acgetchell/la-stack/issues/276)), preserving the
   dimension plot, scientific eligibility rules, and historical artifacts.
3. Release validation and the `v0.4.6` versus `v0.5.0` performance record completed
   through the adopted workflow. The crates.io publishing workflow and publisher
   settings were verified ([#259](https://github.com/acgetchell/la-stack/issues/259)).
4. The reviewed `v0.5.0` release was published and registry availability verified;
   #259 closed after publication on 2026-10-17.

Consumers using a `0.4` Cargo requirement must explicitly adopt `0.5`.
Downstream adoption retains independent correctness and performance gates. Package versions
and generated release artifacts are updated during dedicated release preparation.

### v0.5.1 Follow-Up Work (planned)

The [v0.5.1 milestone](https://github.com/acgetchell/la-stack/milestone/10)
holds compatible follow-up work after `v0.5.0`. The negative #250 decision
does not transfer an algebraic floating-point adoption task to this milestone.
Changes incompatible with a published `v0.5` API require an appropriate later
minor release.

## Longer-Term Directions

- Keep the crate focused on fixed small dimensions and stack allocation.
- Preserve `const fn` surfaces where they make mathematical evaluation
  compile-time checkable.
- Grow exact-arithmetic support only where it strengthens robustness-sensitive
  workflows without making the default build dependency-heavy.
- Use benchmarks to validate performance claims against `nalgebra` and `faer`
  within the crate's intended dimensional scope.

### Expanded Const Generics (unscheduled)

[#123](https://github.com/acgetchell/la-stack/issues/123) tracks the determinant
API work in the [Future: expanded const generics milestone](https://github.com/acgetchell/la-stack/milestone/11).
Stable const generics already support the crate's fixed dimensions. The future
work requires richer expressions involving generic parameters and depends on
the exact capabilities needed to express dimension-specific support in types.

The Rust project's [October 2, 2026 update](https://blog.rust-lang.org/inside-rust/2026/10/02/generic-const-args-and-you/)
describes Generic Const Arguments (GCA) as the intended replacement for
`generic_const_exprs`. None of the discussed features were ready for stabilization
at that time, and the post gives no stabilization date. Evaluate the relevant
GCA features or equivalent support independently; do not assume that one
stabilized subset will enable every proposed determinant abstraction.

When sufficient support is stable, revisit determinant APIs to make closed-form
dimension support explicit and reduce runtime `D` branching where that improves
caller usability. Unify direct, LU-backed, error-bound, and exact determinant
internals where justified by correctness, maintainability, and benchmark evidence.
Assign a numbered release after confirming the required capabilities, API design,
MSRV, and downstream migration. This work is not a prerequisite for v0.5.0 or v1.0.
If a redesign breaks an API already published in v1.x, it requires a major release.

## Non-Goals

- Dynamic or large matrix dimensions.
- Sparse matrices, parallel solvers, GPU acceleration, or out-of-core storage.
- Alternate floating-point scalar families such as `f32` or `f16`.
- Replacing full-featured linear-algebra libraries.
