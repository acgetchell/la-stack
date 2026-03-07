# AGENTS.md

This file provides guidance for automated agents (including Warp at warp.dev) when working with code in this repository.

## Priorities

When making changes in this repo, prioritize (in order):

- Correctness
- Speed
- Coverage (but keep the code idiomatic Rust)

## Common commands

- All tests (Rust + Python): `just test-all`
- Benchmarks: `cargo bench` (or `just bench`)
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
  `cargo run --example const_det_4x4` / `cargo run --features exact --example exact_sign_3x3`)
- Spell check: `just spell-check` (uses `typos.toml` at repo root; add false positives to `[default.extend-words]`)

## Feature flags

- `exact` — enables `det_sign_exact()` (adaptive-precision determinant sign via `BigRational`).
  Gates `src/exact.rs`, additional tests, and the `exact_sign_3x3` example.
  Clippy, doc builds, and test commands have dedicated `--features exact` variants.

## Code structure (big picture)

- This is a single Rust *library crate* (no `src/main.rs`). The crate root is `src/lib.rs`.
- The linear algebra implementation is split across:
  - `src/lib.rs`: crate root + shared items (`LaError`, `DEFAULT_SINGULAR_TOL`, `DEFAULT_PIVOT_TOL`) + re-exports
  - `src/vector.rs`: `Vector<const D: usize>` (`[f64; D]`)
  - `src/matrix.rs`: `Matrix<const D: usize>` (`[[f64; D]; D]`) + helpers (`get`, `set`, `inf_norm`, `det`, `det_direct`)
  - `src/lu.rs`: `Lu<const D: usize>` factorization with partial pivoting (`solve_vec`, `det`)
  - `src/ldlt.rs`: `Ldlt<const D: usize>` factorization without pivoting for symmetric SPD/PSD matrices (`solve_vec`, `det`)
  - `src/exact.rs`: `det_sign_exact()` — adaptive-precision determinant sign
    (Shewchuk-style f64 filter + Bareiss in `BigRational`); `features = ["exact"]`
- Rust tests are inline `#[cfg(test)]` modules in each `src/*.rs` file.
- Python tests live in `scripts/tests/` and run via `just test-python` (`uv run pytest`).
- The public API re-exports these items from `src/lib.rs`.
- The `justfile` defines all dev workflows (see `just --list`).
- Dev-only benchmarks live in `benches/vs_linalg.rs` (Criterion + nalgebra/faer comparison).
- Python scripts under `scripts/`:
  - `criterion_dim_plot.py`: benchmark plotting (CSV + SVG + README table update)
  - `tag_release.py`: annotated tag creation from CHANGELOG.md sections
  - `postprocess_changelog.py`: strips trailing blank lines from git-cliff output
  - `subprocess_utils.py`: safe subprocess wrappers for git commands
- Release workflow is documented in `docs/RELEASING.md`.

## Publishing note

- If you publish this crate to crates.io, prefer updating documentation *before* publishing a new version (doc-only changes still require a version bump on crates.io).
