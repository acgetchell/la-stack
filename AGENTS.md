# AGENTS.md

Essential guidance for AI assistants working in this repository.

## Priorities

When making changes in this repo, prioritize (in order):

- Mathematical correctness and invariant preservation
- API stability and composability
- Idiomatic, well-tested Rust
- Performance within the documented scope

## Design Principles

This is a scientific linear-algebra library. Design decisions trade off in
roughly this priority: mathematical correctness → API stability →
composability → idiomatic Rust → performance within scope. The sections
below spell out what each means in practice; when in doubt, favour the
invariant over the convenient edit.

### Mathematical correctness as an invariant

- Arbitrary-precision paths (`det_exact`, `solve_exact`) never silently lose
  precision. Strict exact-to-`f64` methods (`det_exact_f64`,
  `solve_exact_f64`) return [`LaError::Unrepresentable`] rather than rounding:
  [`UnrepresentableReason::RequiresRounding`] means a finite `f64` is available
  only after rounding, while [`UnrepresentableReason::NotFinite`] means no
  finite `f64` can represent the result. The explicit `*_exact_rounded_f64`
  methods opt into rounding but still return `NotFinite` when rounding cannot
  produce a finite value.
- New or changed f64 operations that can accumulate rounding error document
  their absolute bound (`det_errbound`, `ERR_COEFF_*`) or explicitly state
  that no bound is provided.
- Non-finite matrix/vector inputs and arithmetic intermediates surface as
  `LaError::NonFinite` with typed `NonFiniteOrigin` and `NonFiniteLocation`
  metadata; do not silently propagate NaN or use `unwrap_or(f64::NAN)`.
  Computed failures name their `ArithmeticOperation`, while raw matrix/vector
  values use input origins and exact source locations.
- Exact singularity and tolerance-based rejection remain distinguishable through
  `SingularityReason`; numerical failures preserve the factorization, observed
  pivot magnitude, and tolerance.
- Parse raw tolerances through `Tolerance::try_new`; failures use typed
  `InvalidToleranceReason`. Exact-to-f64 output failures use
  `LaError::Unrepresentable`.
- Algorithms cite their source (Shewchuk, Bareiss, Goldberg, …) via
  `REFERENCES.md` and document their conditioning behaviour.
- `Matrix::det()` uses closed forms through D=4. Its D≥5 zero-tolerance LU
  fallback preserves `LaError::Singular` when elimination cannot produce a
  non-zero pivot; floating-point factorization must not relabel that numerical
  failure as an exact `0.0`. Use the exact determinant APIs when exact
  singularity classification is required.

### Public-API stability

- Public error enums and struct-style error variants are `#[non_exhaustive]`;
  downstream matches include wildcard arms and `..`. Public wrapper types are
  `#[must_use]`.
- New functionality is additive by default: use the prelude for ergonomic
  re-exports, and avoid churn for its own sake.
- Pre-1.0 semver: any `0.x.y` release may include breaking API changes when
  they materially improve correctness, orthogonality, performance, or
  long-term API clarity. Do not keep compatibility aliases that weaken the
  public model; document intentional breaks clearly in release notes and commit
  messages.
- Do **not** automatically update the library version in `Cargo.toml`,
  `Cargo.lock`, README dependency snippets, or related docs during ordinary
  feature, fix, review, or hygiene work. Version bumps are maintainer-driven
  release work performed manually as part of `docs/RELEASING.md`; only change
  them when the user explicitly asks for release/version-bump work.

### Composability

- Const-generic `D` for every core type (`Matrix<D>`, `Vector<D>`,
  `Lu<D>`, `Ldlt<D>`). No runtime dimension.
- Stack allocation by default; heap only behind a feature flag or where
  exact arithmetic inherently requires it (`BigInt` / `BigRational`).
- Feature flags isolate optional dependency weight; default builds stay
  dep-minimal.

### Idiomatic Rust as a proxy for mathematical clarity

- `const fn` wherever possible — not for micro-optimisation, but because
  compile-time evaluation forces a pure function of inputs.
- Use `Result<_, LaError>` for all fallible operations. Public library code
  must not panic on user input.
- Panics are reserved for truly unreachable internal invariant violations and
  must be documented when callers could observe them.
- Validation belongs to the lowest type or module that owns the invariant.
  Higher-level APIs preserve and propagate typed `LaError` values rather than
  stringifying them.
- Public APIs that return plain values must be genuinely infallible for all
  representable inputs. If callers can observe failure, return `Result` or
  `Option` instead of relying on `panic!`, `assert!`, `unwrap`, or `expect`.
- Borrow by default (`&T`, `&[T]`); return borrowed views when possible.
- Type and function names match textbook vocabulary (`Matrix`, `Vector`,
  `Lu`, `Ldlt`, `solve`, `det`, `inf_norm`). Avoid Rust-ecosystem
  abstractions that obscure the math.

### Scientific notation in docs

- Unicode math (×, ≤, ≥, ∈, Σ, ², `2^-50`, …) is welcome in doc
  comments — readability trumps ASCII-only preference.
- Reference literature via `REFERENCES.md` numbered citations (e.g.
  `\[8\]`, `\[9-10\]`).
- State invariants mathematically where possible
  (`|A[i][i]| > Σ_{j≠i} |A[i][j]|`) rather than prose-only.

### Performance within scope

- Performance is a design goal, but strictly subordinate to the
  principles above. Never trade correctness, stability, or clarity for
  speed; if the two conflict, re-scope the problem rather than
  compromise the invariant.
- The library earns its speed through *deliberate scope restriction*:
  fixed small dimensions via const generics, stack-allocated storage,
  and closed-form algorithms where available (D ≤ 4 for `det_direct` /
  `det_errbound`). Problems outside this scope — large or dynamic
  dimensions, sparse matrices, parallelism — belong to `nalgebra` or
  `faer` (see anti-goals in `README.md`).
- Within scope, prefer allocation-free paths, `const fn` wherever the inputs
  allow, and FMA where applicable.
- Performance-sensitive changes require comparable before-and-after evidence
  from the same representative benchmark command, inputs, features, and
  environment. Use `bench-vs-linalg` (vs nalgebra / faer) or `bench-exact`
  (exact arithmetic), as appropriate.
- Preserve benchmark provenance and distinguish descriptive point-estimate
  ratios from statistically supported performance claims. Marginal Criterion
  interval separation is not a paired confidence interval for the change.
- Measurements from runs that violate documented invariants are invalid
  performance evidence.

### Testing mirrors the principles

- Unit tests cover known values, error paths, and dimension-generic
  correctness across D=2..=5 (see **Dimension Coverage** below).
- Error-path tests match the exact variant, typed reason/origin/location, and
  structured fields; do not replace an unexpected error with a numeric sentinel
  or assert only `is_err()`.
- Proptests under `tests/proptest_*.rs` cover algebraic invariants
  (round-trip, residual, sign agreement) — not just "does it not panic".
- Adversarial inputs (near-singular, large-entry, Hilbert-style
  ill-conditioning) accompany well-conditioned inputs in both tests and
  benchmarks.
- When a public API has two paths for the same question (fast filter +
  exact fallback), a proptest verifies they agree on the domain where
  both are defined.

## Core Rules

### Git Operations

- **NEVER** run `git commit`, `git push`, `git tag`, or any git commands that modify version control state
- **ALLOWED**: Run read-only git commands (e.g. `git --no-pager status`, `git --no-pager diff`,
  `git --no-pager log`, `git --no-pager show`, `git --no-pager blame`) to inspect changes/history
- **ALWAYS** use `git --no-pager` when reading git output
- Suggest git commands that modify version control state for the user to run manually
- Do not revert user changes. The worktree may be dirty; preserve unrelated
  changes and work around overlapping edits.
- When suggesting branch names, prefer `{type}/{issue}-descriptor-or-two`, e.g. `fix/307-topology-validation`,
  `perf/315-bench-profile`, or `doc/329-branch-guidance`. If an environment requires an owner/tool prefix,
  keep this structure after the prefix, e.g. `codex/fix/307-topology-validation`.

### Commit Messages

When user requests commit message generation:

1. Run `git --no-pager diff --cached --stat`
2. Generate conventional commit format: `<type>: <brief summary>`
3. Types: `feat`, `fix`, `refactor`, `perf`, `docs`, `test`, `chore`, `style`, `ci`, `build`
4. Include body with organized bullet points and test results
5. Present in code block (no language) - user will commit manually

### Code Quality

- **Unsafe Rust is forbidden.** Keep the manifest-level `unsafe_code = "forbid"`
  lint and crate/module `#![forbid(unsafe_code)]` enforcement intact.
- **ALLOWED**: Run formatters/linters: `cargo fmt`, `cargo clippy`, `cargo doc`, `taplo fmt`, `taplo lint`,
  `uv run --locked ruff check --fix`, `uv run --locked ruff format`, `rumdl`, `dprint`,
  `typos`, `actionlint`
- **NEVER**: Use `sed`, `awk`, `perl` for code edits
- **ALWAYS**: Use the provided structured patch/edit tool for manual edits.
- **FALLBACK**: Direct patch rejects and backup files outside the repository
  and clean them immediately.
- **EXCEPTION**: Shell text tools OK for read-only analysis only

### Validation

- Select validators proportionally to the changed surfaces. Use focused recipes
  for documentation, configuration, Python, test-only, benchmark-only, or
  example-only changes; compose each relevant validator once when a patch spans
  multiple surfaces. Core Rust or public-behavior changes require final
  `just ci`.
- **JSON**: Validate with `jq empty <file>.json` after editing (or `just validate-json`)
- **TOML**: Lint/format with taplo: `just toml-lint`, `just toml-fmt-check`, `just toml-fmt`
- **GitHub Actions**: Validate workflows with `just action-lint` (uses `actionlint`)
- **Spell check**: Run `just spell-check` after editing; add legitimate technical terms to
  `typos.toml` under `[default.extend-words]`
- **Shell scripts**: Run `just shell-fix` and `just shell-check` after editing
- **YAML**: Use `just yaml-lint` and `just yaml-fix`
- **Markdown**: Use `just markdown-check` and `just markdown-fix`

### Rust

- The current MSRV and pinned contributor/CI toolchain are Rust 1.97.0. Keep
  `Cargo.toml`, `rust-toolchain.toml`, and `clippy.toml` aligned when that
  baseline changes deliberately.
- Prefer borrowed APIs by default:
  take references (`&T`, `&mut T`, `&[T]`) as arguments and return borrowed views (`&T`, `&[T]`) when possible.
  Only take ownership or return `Vec`/allocated data when required.

### Documentation

- `src/lib.rs` includes `README.md` with `#![doc = include_str!("../README.md")]`, so README examples are the
  docs.rs landing page examples.
- When changing Rust examples in `README.md`, mirror executable versions in the private `readme_doctests` module in
  `src/lib.rs`. Keep mirrors hidden/private so they do not duplicate the docs.rs landing page, but make them runnable
  by `cargo test --doc`.
- README examples that require optional features may remain `rust,ignore` in README for default-feature doctest
  compatibility, but must have a `#[cfg(feature = "...")]` hidden doctest mirror in `src/lib.rs` and be verified with
  the matching feature set (for example, `cargo test --features exact --doc`).
- When intentionally updating package versions or dependency snippets, keep README `la-stack` dependency examples in
  sync with the package `version` in `Cargo.toml`. Do not perform version bumps unless explicitly requested by the
  maintainer; see **Public-API stability** above.

### Dimension Coverage (2D–5D)

This library uses `const`-generic dimensions. Tests for dimension-generic code
**must cover D=2 through D=5** whenever possible.

#### Use macros for per-dimension test generation

Define a macro that accepts a dimension literal and generates the full set
of test functions for that dimension. Invoke it once per dimension:

```rust
macro_rules! gen_tests {
    ($d:literal) => {
        paste! {
            #[test]
            fn [<test_foo_ $d d>]() {
                // assertions …
            }
        }
    };
}

gen_tests!(2);
gen_tests!(3);
gen_tests!(4);
gen_tests!(5);
```

#### Keep core logic in generic helper functions

The macro body should be thin — primarily calling `const`-generic helpers and
asserting results. This keeps the macro readable and the helpers independently
testable.

#### Reference examples

- `src/matrix.rs` — `gen_matrix_tests!`
- `src/lu.rs` — `gen_pivoting_solve_and_det_tests!`, `gen_tridiagonal_smoke_solve_and_det_tests!`
- `src/ldlt.rs` — `gen_ldlt_identity_tests!`, `gen_ldlt_diagonal_tests!`
- `src/exact.rs` — `gen_det_exact_tests!`, `gen_det_exact_f64_tests!`, `gen_solve_exact_tests!`, `gen_solve_exact_f64_tests!`

#### When single-dimension tests are acceptable

Some tests are inherently dimension-specific (e.g. known values for a crafted
matrix, error-handling with a specific layout). These do not need
macro-ification.

### Python

- Python support tooling targets Python 3.14.
- Use `uv run --locked` for all Python scripts (never `python3` or `python` directly)
- Use pytest for tests (not unittest)
- **Type checking**: `just python-check` includes type checking (blocking - all code must pass type checks)
- Add type hints to new code

## Common Commands

```bash
just check            # Lint/validators (non-mutating)
just fix              # Apply formatters/auto-fixes (mutating)
just ci               # Full CI simulation (checks + tests + examples + bench compile)
just test             # Lib + doc tests (fast)
just test-all         # All tests (Rust, benchmark inputs, and Python)
just examples         # Run all examples
```

### Detailed Command Reference

- All tests (Rust, exact-feature doctests, benchmark-input smoke tests, and Python): `just test-all`
- Benchmark comparison (local report): `just bench-compare [baseline] [suite] [scope]`
- Benchmarks: `cargo bench --locked --features bench` (or `just bench`)
- Benchmarks (exact arithmetic): `just bench-exact`
- Benchmarks (la-stack vs nalgebra/faer): `just bench-vs-linalg [filter]` (full run) or `just bench-vs-linalg-quick [filter]` (reduced)
- Benchmarks (plot vs_linalg CSV/SVG/JSON provenance): `just plot-vs-linalg [metric] [stat] [sample] [log_y]`;
  publish a freshly gated full run to README with
  `just plot-vs-linalg-readme [metric] [stat] [sample] [log_y]`
- Benchmarks (save baseline): `just bench-save-baseline v0.4.1`
- Build (debug): `cargo build` (or `just build`)
- Build (release): `cargo build --release` (or `just build-release`)
- Changelog (generate full): `just changelog` (generates, post-processes, archives, and formats changelog files)
- Changelog (prepend unreleased): `just changelog-unreleased v0.4.1`
- Coverage (CI XML): `just coverage-ci`
- Coverage (HTML): `just coverage`
- Create release tag: `just tag v0.4.1` (creates annotated tag from CHANGELOG.md section) / `just tag-force v0.4.1` (recreate if the tag already exists)
- Fast compile check (no binary produced): `cargo check` (or `just check-fast`)
- Fast Rust tests (lib + doc): `just test`
- Format: `cargo fmt` (or `just fmt`)
- Integration tests: `just test-integration`
- Benchmark-input smoke tests: `just test-bench-inputs`
- Lint (Clippy, canonical default and all-feature passes): `just clippy`
- Lint (Clippy, focused exact-feature pass): `just clippy-exact`
- Lint/validate: `just check`
- Cargo manifest/lockfile synchronization: `just cargo-lock-check`
- Unused dependency check: `just unused-deps` (uses `cargo-machete`)
- Pre-commit validation / CI simulation: `just ci` (lint + tests + examples + bench compile)
- Python setup from the lockfile: `uv sync --locked --group dev` (or `just python-sync`)
- Python tests: `just test-python`
- Run one runnable test by substring: `cargo nextest run solve_2x2_basic`
  - For an exact full-path match, use `cargo nextest run -- --exact lu::tests::solve_2x2_basic`.
- Run exact-feature tests: `cargo nextest run --profile ci --features exact --verbose`
  (or `just test-exact`, which also runs exact-feature doctests)
- Run examples: `just examples` (or `cargo run --example det_5x5` / `cargo run --example solve_5x5` /
  `cargo run --example ldlt_solve_3x3` / `cargo run --example const_det_4x4` /
  `cargo run --features exact --example exact_det_3x3` /
  `cargo run --features exact --example exact_sign_3x3` /
  `cargo run --features exact --example exact_solve_3x3`)
- Spell check: `just spell-check` (uses `typos.toml` at repo root; add false positives to `[default.extend-words]`)

### Changelog

- Never edit `CHANGELOG.md` directly - it's auto-generated from git commits
- Use `just changelog` to regenerate
- Use `just changelog-unreleased <version>` to prepend unreleased changes

### GitHub CLI (`gh`)

When using `gh` to view issues, PRs, or other GitHub objects:

- **ALWAYS** use `--json` with `| cat` to avoid pager and scope errors:

  ```bash
  gh issue view 64 --repo acgetchell/la-stack --json title,body | cat
  ```

- To extract specific fields cleanly, combine `--json` with `--jq`:

  ```bash
  gh issue view 64 --repo acgetchell/la-stack --json title,body --jq '.title + "\n" + .body' | cat
  ```

- **AVOID** plain `gh issue view N` — it may fail with `read:project`
  scope errors or open a pager.

- For **arbitrary Markdown** (backticks, quotes, special characters) in
  comments, prefer `--body-file -` with a heredoc:

  ```bash
  gh issue comment 64 --repo acgetchell/la-stack --body-file - <<'EOF'
  ## Heading

  Body with `backticks`, **bold**, and apostrophes that's safe.
  EOF
  ```

### GitHub Issues

Use the `gh` CLI to read, create, and edit issues:

- **Read**: `gh issue view <number> --json title,body,labels,milestone | cat`
- **List**:
  `gh issue list --json number,title,labels --jq '.[] | "#\(.number) \(.title)"' | cat`
  (add `--label enhancement`, `--milestone v0.4.1`, etc. to filter)
- **Create**: `gh issue create --title "..." --body "..." --label enhancement --label rust`
- **Edit**: `gh issue edit <number> --add-label "..."`, `--milestone "..."`, `--title "..."`
- **Comment**: `gh issue comment <number> --body "..."`
- **Close**: `gh issue close <number>` (with optional `--reason completed` or `--reason "not planned"`)

When creating or updating issues:

- **Labels**: Use appropriate labels: `enhancement`, `bug`, `performance`, `documentation`, `rust`, `python`, etc.
- **Milestones**: Assign to the appropriate milestone (e.g., `v0.4.1`, `v0.5.0`)
- **Dependencies**: Document relationships in issue body and comments:
  - "Depends on: #XXX" - this issue cannot start until #XXX is complete
  - "Blocks: #YYY" - #YYY cannot start until this issue is complete
  - "Related: #ZZZ" - related work but not blocking
- **Relationships**: GitHub automatically parses blocking keywords in comments to create visual relationships:
  - Use `gh issue comment <number> --body "Blocked by #XXX"` to mark an issue as blocked
  - Use `gh issue comment <number> --body "Blocks #YYY"` to mark an issue as blocking another
  - GitHub will automatically create the relationship graph in the web UI
  - Example: `gh issue comment 217 --body "Blocked by #207"` creates a blocking dependency
- **Issue body format**: Include clear sections: Summary, Current State, Proposed Changes, Benefits, Implementation Notes
- **Cross-referencing**: Always reference related issues/PRs using #XXX notation for automatic linking

## Feature flags

- `exact` — enables exact arithmetic methods via `BigRational`:
  `det_exact()`, `det_exact_f64()`, `det_exact_rounded_f64()`,
  `det_sign_exact()`, `solve_exact()`, `solve_exact_f64()`, and
  `solve_exact_rounded_f64()`. `det_sign_exact()` is infallible for every
  finite-by-construction `Matrix`; the exact-value, conversion, and solve APIs
  remain fallible for their genuine scale, representation, and singularity
  failures. `ExactF64Conversion` converts an already-computed
  exact determinant or solution under the strict or rounded contract without
  rerunning exact elimination. Feature-gated re-exports include
  `DeterminantSign`, `ExactF64Conversion`, `BigInt`, `BigRational`, and the
  commonly needed `num-traits` items (`FromPrimitive`, `ToPrimitive`, and `Signed`).
  `UnrepresentableReason` and the other typed `LaError` category enums remain
  available without `exact`; callers should not need optional arithmetic
  dependencies merely to match errors.
  Gates `src/exact.rs`, additional tests, and the exact-arithmetic examples.
  Clippy, doc builds, and test commands have dedicated `--features exact`
  variants.
- `bench` — cfg-only gate required by the benchmark targets and
  `tests/vs_linalg_inputs.rs`. Benchmark libraries remain dev-dependencies.

## Code structure (big picture)

- This is a single Rust *library crate* (no `src/main.rs`). The crate root is `src/lib.rs`.
- The linear algebra implementation is split across:
  - `src/lib.rs`: crate root, public module wiring, and re-exports
  - `src/error.rs`: `LaError` plus typed singularity, non-finite,
    positive-semidefinite, tolerance, factorization, arithmetic-operation, and
    exact-conversion categories
  - `src/tolerance.rs`: validated singular-tolerance policy
  - `src/vector.rs`: `Vector<const D: usize>` (`[f64; D]`)
  - `src/matrix.rs`: `Matrix<const D: usize>` (`[[f64; D]; D]`) + helpers (`get`, `try_get`, `set`, `inf_norm`, `det`, `det_direct`)
  - `src/lu.rs`: `Lu<const D: usize>` factorization with partial pivoting (`solve`, `det`)
  - `src/ldlt.rs`: `Ldlt<const D: usize>` factorization without pivoting for exactly
    symmetric positive-definite matrices (`solve`, `det`)
  - `src/exact.rs`: exact arithmetic behind `features = ["exact"]`:
    - Determinants: `det_exact()`, strict `det_exact_f64()`, rounded
      `det_exact_rounded_f64()`, and `det_sign_exact()` via a scaled `BigInt`
      determinant core (`exact_det_int_finite`): direct expansions for D≤4 and
      fraction-free Bareiss elimination for D≥5. `det_sign_exact()` infallibly
      returns a `DeterminantSign` and adds a Shewchuk-style f64 filter for fast
      sign resolution in D≤4
    - Exact-to-`f64` conversion failures retain an `UnrepresentableReason` so
      callers can distinguish required rounding from non-finite output
    - Linear system solve: `solve_exact()`, strict `solve_exact_f64()`, and
      rounded `solve_exact_rounded_f64()` use fraction-free Bareiss forward
      elimination in `BigInt` with first-non-zero pivoting, followed by
      `BigRational` back-substitution
- Rust unit tests are inline `#[cfg(test)]` modules in each `src/*.rs` file.
- Property-based tests live under `tests/proptest_*.rs` (uses the `proptest`
  dev-dependency): `proptest_matrix.rs`, `proptest_vector.rs`,
  `proptest_factorizations.rs`, and `proptest_exact.rs` (the last gated on
  the `exact` feature). They run as integration tests via
  `just test-integration` or `just test-all`.
- Python tests live in `scripts/tests/` and run via `just test-python` (`uv run --locked pytest`).
- The public API re-exports these items from `src/lib.rs`.
- The `justfile` defines all dev workflows (see `just --list`).
- Dev-only benchmarks live in `benches/vs_linalg.rs` (Criterion + nalgebra/faer comparison)
  and `benches/exact.rs` (exact arithmetic across D=2–5, plus adversarial-input groups
  `exact_near_singular_3x3`, `exact_large_entries_3x3`, `exact_hilbert_4x4`, `exact_hilbert_5x5`).
  Exact Criterion helpers accept only `ValidatedExactInput`, so independent
  oracle validation is a type-checked prerequisite outside timed closures.
- Key Python scripts under `scripts/`:
  - `bench_compare.py`: exact and vs-linalg Criterion comparison reports under
    `target/bench-reports/`
  - `archive_performance.py`: promote and archive curated release performance reports
  - `criterion_dim_plot.py`: benchmark plotting and fail-closed README publication
    (CSV + SVG + JSON provenance + README table)
  - `tag_release.py`: annotated tag creation from CHANGELOG.md sections
  - `archive_changelog.py`: archive completed changelog minor series
  - `postprocess_changelog.py`: inject summaries, reflow and normalize Markdown,
    and strip trailing blank lines from git-cliff output
  - `subprocess_utils.py`: safe subprocess wrappers for git commands
- Release workflow is documented in `docs/RELEASING.md`.

## Agent Expectations

- Prefer small, focused patches and the simplest maintainable correct solution.
- Search existing documentation and nearby code before inventing conventions.
- Fix small, clearly related issues discovered in a touched area when doing so
  improves correctness, clarity, tests, or maintainability.
- Avoid broad mechanical churn; separate repository-wide cleanup from focused
  work.

## Publishing note

- If you publish this crate to crates.io, prefer updating documentation
  *before* publishing a new version (doc-only changes still require a version bump on crates.io).

## Editing tools policy

- Never use `sed`, `awk`, `python`, or `perl` to edit code or write file changes.
- These tools may be used for read-only inspection, parsing, or analysis, but never for writing.
