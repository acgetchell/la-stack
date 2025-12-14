# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Common commands

- Build (debug): `cargo build`
- Build (release): `cargo build --release`
- Fast compile check (no binary produced): `cargo check`
- Run tests: `cargo test`
- Run a single test (by name filter): `cargo test solve_2x2_basic` (or the full path: `cargo test lu::tests::solve_2x2_basic`)
- Format: `cargo fmt`
- Lint (Clippy): `cargo clippy --all-targets --all-features -- -D warnings`
- Run benchmarks: `cargo bench`

## Code structure (big picture)

- This is a single Rust *library crate* (no `src/main.rs`). The crate root is `src/lib.rs`.
- The linear algebra implementation is split across:
  - `src/lib.rs`: crate root + shared items (`LaError`, `DEFAULT_PIVOT_TOL`) + re-exports
  - `src/vector.rs`: `Vector<const D: usize>` (`[f64; D]`)
  - `src/matrix.rs`: `Matrix<const D: usize>` (`[[f64; D]; D]`) + helpers (`get`, `set`, `inf_norm`, `det`)
  - `src/lu.rs`: `Lu<const D: usize>` factorization with partial pivoting (`solve_vec`, `det`)
- The public API re-exports these items from `src/lib.rs`.
- Dev-only benchmarks live in `benches/vs_nalgebra.rs` (Criterion + nalgebra comparison).

## Publishing note

- If you publish this crate to crates.io, prefer updating documentation *before* publishing a new version (doc-only changes still require a version bump on crates.io).