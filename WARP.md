# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Common commands

- Build (debug): `cargo build` (or `just build`)
- Build (release): `cargo build --release` (or `just build-release`)
- Fast compile check (no binary produced): `cargo check` (or `just check`)
- Run tests: `cargo test` (or `just test`)
- Run a single test (by name filter): `cargo test solve_2x2_basic` (or the full path: `cargo test lu::tests::solve_2x2_basic`)
- Format: `cargo fmt` (or `just fmt`)
- Lint (Clippy): `cargo clippy --all-targets --all-features -- -D warnings` (or `just clippy`)
- Spell check: `just spell-check` (uses `cspell.json` at repo root; keep the `words` list sorted lexicographically)
- Run benchmarks: `cargo bench` (or `just bench`)
- Run examples: `just examples` (or `cargo run --example det_3x3`)
- CI simulation (lint + tests + bench compile): `just ci`
- Pre-commit validation: `just commit-check`

## Code structure (big picture)

- This is a single Rust *library crate* (no `src/main.rs`). The crate root is `src/lib.rs`.
- The linear algebra implementation is split across:
  - `src/lib.rs`: crate root + shared items (`LaError`, `DEFAULT_PIVOT_TOL`) + re-exports
  - `src/vector.rs`: `Vector<const D: usize>` (`[f64; D]`)
  - `src/matrix.rs`: `Matrix<const D: usize>` (`[[f64; D]; D]`) + helpers (`get`, `set`, `inf_norm`, `det`)
  - `src/lu.rs`: `Lu<const D: usize>` factorization with partial pivoting (`solve_vec`, `det`)
- A minimal `justfile` exists for common workflows (see `just --list`).
- The public API re-exports these items from `src/lib.rs`.
- Dev-only benchmarks live in `benches/vs_nalgebra.rs` (Criterion + nalgebra comparison).

## Publishing note

- If you publish this crate to crates.io, prefer updating documentation *before* publishing a new version (doc-only changes still require a version bump on crates.io).