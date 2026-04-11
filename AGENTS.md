# AGENTS.md

Essential guidance for AI assistants working in this repository.

## Priorities

When making changes in this repo, prioritize (in order):

- Correctness
- Speed
- Coverage (but keep the code idiomatic Rust)

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
- Benchmark comparison (generate `docs/PERFORMANCE.md`): `just bench-compare` (snapshot) or `just bench-compare v0.3.0` (vs baseline)
- Benchmarks: `cargo bench` (or `just bench`)
- Benchmarks (exact arithmetic): `just bench-exact`
- Benchmarks (save baseline): `just bench-save-baseline v0.3.0`
- Build (debug): `cargo build` (or `just build`)
- Build (release): `cargo build --release` (or `just build-release`)
- Changelog (generate full): `just changelog` (runs `git-cliff -o CHANGELOG.md` + post-processing)
- Changelog (prepend unreleased): `just changelog-unreleased v0.3.0`
- Coverage (CI XML): `just coverage-ci`
- Coverage (HTML): `just coverage`
- Create release tag: `just tag v0.3.0` (creates annotated tag from CHANGELOG.md section)
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

### GitHub Issues

Use the `gh` CLI to read, create, and edit issues:

- **Read**: `gh issue view <number>` (or `--json title,body,labels,milestone` for structured data)
- **List**: `gh issue list` (add `--label enhancement`, `--milestone v0.3.0`, etc. to filter)
- **Create**: `gh issue create --title "..." --body "..." --label enhancement --label rust`
- **Edit**: `gh issue edit <number> --add-label "..."`, `--milestone "..."`, `--title "..."`
- **Comment**: `gh issue comment <number> --body "..."`
- **Close**: `gh issue close <number>` (with optional `--reason completed` or `--reason "not planned"`)

When creating or updating issues:

- **Labels**: Use appropriate labels: `enhancement`, `bug`, `performance`, `documentation`, `rust`, `python`, etc.
- **Milestones**: Assign to the appropriate milestone (e.g., `v0.3.0`, `v0.4.0`)
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
  Also re-exports `BigRational` from the crate root and prelude.
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
    - Determinants: `det_exact()`, `det_exact_f64()`, `det_sign_exact()` via Bareiss in
      `BigRational`; `det_sign_exact()` adds a Shewchuk-style f64 filter for fast sign resolution
    - Linear system solve: `solve_exact()`, `solve_exact_f64()` via Gaussian elimination
      with first-non-zero pivoting in `BigRational`
- Rust tests are inline `#[cfg(test)]` modules in each `src/*.rs` file.
- Python tests live in `scripts/tests/` and run via `just test-python` (`uv run pytest`).
- The public API re-exports these items from `src/lib.rs`.
- The `justfile` defines all dev workflows (see `just --list`).
- Dev-only benchmarks live in `benches/vs_linalg.rs` (Criterion + nalgebra/faer comparison)
  and `benches/exact.rs` (exact arithmetic across D=2–5).
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
