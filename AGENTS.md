# AGENTS.md

Essential guidance for AI assistants working in this repository.

## Priorities

When making changes in this repo, prioritize (in order):

- Correctness
- Speed
- Coverage (but keep the code idiomatic Rust)

## Design Principles

This is a scientific linear-algebra library.  Design decisions trade off in
roughly this priority: mathematical correctness → API stability →
composability → idiomatic Rust → performance within scope.  The sections
below spell out what each means in practice; when in doubt, favour the
invariant over the convenient edit.

### Mathematical correctness as an invariant

- Exact paths (`*_exact`) never silently lose precision.  When f64 output
  is required, a separate `*_exact_f64` method returns
  [`LaError::Overflow`] on unrepresentability — not a truncation.
- Any f64 operation that can accumulate rounding error either documents
  its absolute bound (`det_errbound`, `ERR_COEFF_*`) or explicitly states
  that no bound is provided.
- Non-finite values (NaN, ±∞) always surface as
  `LaError::NonFinite { row, col }` with source-location metadata.  No
  silent NaN propagation, no `unwrap_or(f64::NAN)`.
- Algorithms cite their source (Shewchuk, Bareiss, Goldberg, …) via
  `REFERENCES.md` and document their conditioning behaviour.

### Public-API stability

- Error enums are `#[non_exhaustive]`; public wrapper types are
  `#[must_use]`.
- New functionality is additive: use the prelude for ergonomic re-exports;
  never silently rename or remove a public item.
- Pre-1.0 semver: `0.x.Y` is a patch-level additive bump, `0.X.y` is a
  minor bump that may include breaking changes.  Conventional-commit
  types (`feat`, `fix`, `refactor`, …) mirror this convention.

### Composability

- Const-generic `D` for every core type (`Matrix<D>`, `Vector<D>`,
  `Lu<D>`, `Ldlt<D>`).  No runtime dimension.
- Stack allocation by default; heap only behind a feature flag or where
  exact arithmetic inherently requires it (`BigInt` / `BigRational`).
- Feature flags isolate optional dependency weight; default builds stay
  dep-minimal.

### Idiomatic Rust as a proxy for mathematical clarity

- `const fn` wherever possible — not for micro-optimisation, but because
  compile-time evaluation forces a pure function of inputs.
- `Result<_, LaError>` for all fallible operations.  Panics are reserved
  for debug-only precondition violations (e.g. LDLT symmetry check) and
  documented on the method.
- Borrow by default (`&T`, `&[T]`); return borrowed views when possible.
- Type and function names match textbook vocabulary (`Matrix`, `Vector`,
  `Lu`, `Ldlt`, `solve_vec`, `det`, `inf_norm`).  Avoid Rust-ecosystem
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
  principles above.  Never trade correctness, stability, or clarity for
  speed; if the two conflict, re-scope the problem rather than
  compromise the invariant.
- The library earns its speed through *deliberate scope restriction*:
  fixed small dimensions via const generics, stack-allocated storage,
  and closed-form algorithms where available (D ≤ 4 for `det_direct` /
  `det_errbound`).  Problems outside this scope — large or dynamic
  dimensions, sparse matrices, parallelism — belong to `nalgebra` or
  `faer` (see anti-goals in `README.md`).
- Within scope, prefer allocation-free paths, `const fn` wherever the
  inputs allow, and FMA where applicable.  Validate any performance
  claim against the `bench-vs-linalg` (vs nalgebra / faer) or
  `bench-exact` (exact-arithmetic) suites before relying on it.

### Testing mirrors the principles

- Unit tests cover known values, error paths, and dimension-generic
  correctness across D=2..=5 (see **Dimension Coverage** below).
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

### Commit Messages

When user requests commit message generation:

1. Run `git --no-pager diff --cached --stat`
2. Generate conventional commit format: `<type>: <brief summary>`
3. Types: `feat`, `fix`, `refactor`, `perf`, `docs`, `test`, `chore`, `style`, `ci`, `build`
4. Include body with organized bullet points and test results
5. Present in code block (no language) - user will commit manually

### Code Quality

- **ALLOWED**: Run formatters/linters: `cargo fmt`, `cargo clippy`, `cargo doc`, `taplo fmt`, `taplo lint`,
  `uv run ruff check --fix`, `uv run ruff format`, `shfmt -w`, `shellcheck -x`, `npx markdownlint --fix`,
  `typos`, `actionlint`
- **NEVER**: Use `sed`, `awk`, `perl` for code edits
- **ALWAYS**: Use `edit_files` tool for edits (and `create_file` for new files)
- **EXCEPTION**: Shell text tools OK for read-only analysis only

### Validation

- **JSON**: Validate with `jq empty <file>.json` after editing (or `just validate-json`)
- **TOML**: Lint/format with taplo: `just toml-lint`, `just toml-fmt-check`, `just toml-fmt`
- **GitHub Actions**: Validate workflows with `just action-lint` (uses `actionlint`)
- **Spell check**: Run `just spell-check` after editing; add legitimate technical terms to
  `typos.toml` under `[default.extend-words]`
- **Shell scripts**: Run `shfmt -w scripts/*.sh` and `shellcheck -x scripts/*.sh` after editing
- **YAML**: Use `just yaml-lint` and `just yaml-fix`
- **Markdown**: Use `just markdown-check` and `just markdown-fix`

### Rust

- Prefer borrowed APIs by default:
  take references (`&T`, `&mut T`, `&[T]`) as arguments and return borrowed views (`&T`, `&[T]`) when possible.
  Only take ownership or return `Vec`/allocated data when required.

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

- `src/matrix.rs` — `gen_public_api_matrix_tests!`
- `src/lu.rs` — `gen_public_api_pivoting_solve_vec_and_det_tests!`, `gen_public_api_tridiagonal_smoke_solve_vec_and_det_tests!`
- `src/ldlt.rs` — `gen_public_api_ldlt_identity_tests!`, `gen_public_api_ldlt_diagonal_tests!`
- `src/exact.rs` — `gen_det_exact_tests!`, `gen_det_exact_f64_tests!`, `gen_solve_exact_tests!`, `gen_solve_exact_f64_tests!`

#### When single-dimension tests are acceptable

Some tests are inherently dimension-specific (e.g. known values for a crafted
matrix, error-handling with a specific layout). These do not need
macro-ification.

### Python

- Use `uv run` for all Python scripts (never `python3` or `python` directly)
- Use pytest for tests (not unittest)
- **Type checking**: `just python-check` includes type checking (blocking - all code must pass type checks)
- Add type hints to new code

## Common Commands

```bash
just fix              # Apply formatters/auto-fixes (mutating)
just check            # Lint/validators (non-mutating)
just ci               # Full CI simulation (checks + tests + examples + bench compile)
just test             # Lib + doc tests (fast)
just test-all         # All tests (Rust + Python)
just examples         # Run all examples
```

### Detailed Command Reference

- All tests (Rust + Python): `just test-all`
- Benchmark comparison (generate `docs/PERFORMANCE.md`): `just bench-compare` (snapshot) or `just bench-compare v0.4.1` (vs baseline)
- Benchmarks: `cargo bench` (or `just bench`)
- Benchmarks (exact arithmetic): `just bench-exact`
- Benchmarks (la-stack vs nalgebra/faer): `just bench-vs-linalg [filter]` (full run) or `just bench-vs-linalg-quick [filter]` (reduced)
- Benchmarks (plot vs_linalg CSV/SVG): `just plot-vs-linalg [metric] [stat] [sample] [update_readme]` / `just plot-vs-linalg-readme [metric] [stat] [sample] [update_readme]`
- Benchmarks (save baseline): `just bench-save-baseline v0.4.1`
- Build (debug): `cargo build` (or `just build`)
- Build (release): `cargo build --release` (or `just build-release`)
- Changelog (generate full): `just changelog` (runs `git-cliff -o CHANGELOG.md` + post-processing)
- Changelog (prepend unreleased): `just changelog-unreleased v0.4.1`
- Coverage (CI XML): `just coverage-ci`
- Coverage (HTML): `just coverage`
- Create release tag: `just tag v0.4.1` (creates annotated tag from CHANGELOG.md section) / `just tag-force v0.4.1` (recreate if the tag already exists)
- Fast compile check (no binary produced): `cargo check` (or `just check-fast`)
- Fast Rust tests (lib + doc): `just test`
- Format: `cargo fmt` (or `just fmt`)
- Integration tests: `just test-integration`
- Lint (Clippy): `cargo clippy --all-targets --all-features -- -D warnings` (or `just clippy`)
- Lint (Clippy, exact feature): `cargo clippy --features exact --all-targets -- -D warnings` (or `just clippy-exact`)
- Lint/validate: `just check`
- Pre-commit validation / CI simulation: `just ci` (lint + tests + examples + bench compile)
- Python setup: `uv sync --group dev` (or `just python-sync`)
- Python tests: `just test-python`
- Run a single test (by name filter): `cargo test solve_2x2_basic` (or the full path: `cargo test lu::tests::solve_2x2_basic`)
- Run exact-feature tests: `cargo test --features exact --verbose` (or `just test-exact`)
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
- **List**: `gh issue list --json number,title,labels --jq '.[] | "#\(.number) \(.title)"' | cat` (add `--label enhancement`, `--milestone v0.4.1`, etc. to filter)
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
  `det_exact()`, `det_exact_f64()`, `det_sign_exact()`, `solve_exact()`, and `solve_exact_f64()`.
  Re-exports `BigInt`, `BigRational`, and the commonly needed `num-traits`
  items (`FromPrimitive`, `ToPrimitive`, and `Signed`) from the crate root and prelude
  (so consumers get usable `from_f64` / `to_f64` / `is_positive` etc. without adding
  `num-bigint` / `num-rational` / `num-traits` as their own deps).
  Gates `src/exact.rs`, additional tests, and the `exact_det_3x3`/`exact_sign_3x3`/`exact_solve_3x3` examples.
  Clippy, doc builds, and test commands have dedicated `--features exact` variants.

## Code structure (big picture)

- This is a single Rust *library crate* (no `src/main.rs`). The crate root is `src/lib.rs`.
- The linear algebra implementation is split across:
  - `src/lib.rs`: crate root + shared items (`LaError`, `DEFAULT_SINGULAR_TOL`, `DEFAULT_PIVOT_TOL`) + re-exports
  - `src/vector.rs`: `Vector<const D: usize>` (`[f64; D]`)
  - `src/matrix.rs`: `Matrix<const D: usize>` (`[[f64; D]; D]`) + helpers (`get`, `set`, `inf_norm`, `det`, `det_direct`)
  - `src/lu.rs`: `Lu<const D: usize>` factorization with partial pivoting (`solve_vec`, `det`)
  - `src/ldlt.rs`: `Ldlt<const D: usize>` factorization without pivoting for symmetric SPD/PSD matrices (`solve_vec`, `det`)
  - `src/exact.rs`: exact arithmetic behind `features = ["exact"]`:
    - Determinants: `det_exact()`, `det_exact_f64()`, `det_sign_exact()` via integer-only
      Bareiss in `BigInt` (`bareiss_det_int`); `det_sign_exact()` adds a Shewchuk-style
      f64 filter for fast sign resolution
    - Linear system solve: `solve_exact()`, `solve_exact_f64()` via Gaussian elimination
      with first-non-zero pivoting in `BigRational`
- Rust unit tests are inline `#[cfg(test)]` modules in each `src/*.rs` file.
- Property-based tests live under `tests/proptest_*.rs` (uses the `proptest`
  dev-dependency): `proptest_matrix.rs`, `proptest_vector.rs`,
  `proptest_factorizations.rs`, and `proptest_exact.rs` (the last gated on
  the `exact` feature). They run as integration tests via
  `just test-integration` or `just test-all`.
- Python tests live in `scripts/tests/` and run via `just test-python` (`uv run pytest`).
- The public API re-exports these items from `src/lib.rs`.
- The `justfile` defines all dev workflows (see `just --list`).
- Dev-only benchmarks live in `benches/vs_linalg.rs` (Criterion + nalgebra/faer comparison)
  and `benches/exact.rs` (exact arithmetic across D=2–5, plus adversarial-input groups
  `exact_near_singular_3x3`, `exact_large_entries_3x3`, `exact_hilbert_4x4`, `exact_hilbert_5x5`).
- Python scripts under `scripts/`:
  - `bench_compare.py`: exact-arithmetic benchmark comparison across releases (generates `docs/PERFORMANCE.md`)
  - `criterion_dim_plot.py`: benchmark plotting (CSV + SVG + README table update)
  - `tag_release.py`: annotated tag creation from CHANGELOG.md sections
  - `postprocess_changelog.py`: strips trailing blank lines from git-cliff output
  - `subprocess_utils.py`: safe subprocess wrappers for git commands
- Release workflow is documented in `docs/RELEASING.md`.

## Publishing note

- If you publish this crate to crates.io, prefer updating documentation *before* publishing a new version (doc-only changes still require a version bump on crates.io).

## Editing tools policy

- Never use `sed`, `awk`, `python`, or `perl` to edit code or write file changes.
- These tools may be used for read-only inspection, parsing, or analysis, but never for writing.
