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
- CI simulation (lint + tests + examples + bench compile): `just ci`
- Coverage (CI XML): `just coverage-ci`
- Coverage (HTML): `just coverage`
- Fast compile check (no binary produced): `cargo check` (or `just check-fast`)
- Fast Rust tests (lib + doc): `just test`
- Format: `cargo fmt` (or `just fmt`)
- Integration tests: `just test-integration`
- Lint (Clippy): `cargo clippy --all-targets --all-features -- -D warnings` (or `just clippy`)
- Lint/validate: `just check`
- Pre-commit validation: `just ci`
- Python tests: `just test-python`
- Run a single test (by name filter): `cargo test solve_2x2_basic` (or the full path: `cargo test lu::tests::solve_2x2_basic`)
- Run examples: `just examples` (or `cargo run --example det_5x5` / `cargo run --example solve_5x5`)
- Spell check: `just spell-check` (uses `cspell.json` at repo root; keep the `words` list sorted lexicographically)

## Code structure (big picture)

- This is a single Rust *library crate* (no `src/main.rs`). The crate root is `src/lib.rs`.
- The linear algebra implementation is split across:
  - `src/lib.rs`: crate root + shared items (`LaError`, `DEFAULT_SINGULAR_TOL`, `DEFAULT_PIVOT_TOL`) + re-exports
  - `src/vector.rs`: `Vector<const D: usize>` (`[f64; D]`)
  - `src/matrix.rs`: `Matrix<const D: usize>` (`[[f64; D]; D]`) + helpers (`get`, `set`, `inf_norm`, `det`)
  - `src/lu.rs`: `Lu<const D: usize>` factorization with partial pivoting (`solve_vec`, `det`)
  - `src/ldlt.rs`: `Ldlt<const D: usize>` factorization without pivoting for symmetric SPD/PSD matrices (`solve_vec`, `det`)
- A minimal `justfile` exists for common workflows (see `just --list`).
- The public API re-exports these items from `src/lib.rs`.
- Dev-only benchmarks live in `benches/vs_linalg.rs` (Criterion + nalgebra/faer comparison).

## Publishing note

- If you publish this crate to crates.io, prefer updating documentation *before* publishing a new version (doc-only changes still require a version bump on crates.io).
