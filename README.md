# la-stack

[![Crates.io](https://img.shields.io/crates/v/la-stack.svg)](https://crates.io/crates/la-stack)
[![Downloads](https://img.shields.io/crates/d/la-stack.svg)](https://crates.io/crates/la-stack)
[![License](https://img.shields.io/crates/l/la-stack.svg)](LICENSE)
[![Docs.rs](https://docs.rs/la-stack/badge.svg)](https://docs.rs/la-stack)
[![CI](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/ci.yml)
[![rust-clippy analyze](https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/rust-clippy.yml)
[![Audit dependencies](https://github.com/acgetchell/la-stack/actions/workflows/audit.yml/badge.svg)](https://github.com/acgetchell/la-stack/actions/workflows/audit.yml)

Fast, stack-allocated linear algebra for fixed dimensions in Rust.

This crate exists primarily to support the [`delaunay`](https://crates.io/crates/delaunay) crate’s needs while keeping the API intentionally small and explicit.

## 📐 Introduction

`la-stack` provides a handful of const-generic, stack-backed building blocks:

- `Vector<const D: usize>` for fixed-length vectors (`[f64; D]` today)
- `Matrix<const D: usize>` for fixed-size square matrices (`[[f64; D]; D]` today)
- `Lu<const D: usize>` for LU factorization with partial pivoting (solve + det)

## ✨ Design goals

- ✅ Const-generic dimensions (no dynamic sizes)
- ✅ Stack storage only (no heap allocation in core types)
- ✅ `Copy` types where possible
- ✅ Explicit algorithms (LU, solve, determinant)
- ✅ `unsafe` forbidden
- ✅ No runtime dependencies (dev-dependencies are for contributors only)

## 🔢 Scalar types

Today, the core types are implemented for `f64`. The intent is to support `f32` and `f64` (and `f128` if/when Rust gains a stable primitive for it). Longer term, we may add optional arbitrary-precision support (e.g. via `rug`) depending on performance.

## 🚀 Quickstart

Add this to your `Cargo.toml`:

```toml
[dependencies]
la-stack = "0.1"
```

Solve a 3×3 system via LU:

```rust
use la_stack::{LaError, Matrix, Vector, DEFAULT_PIVOT_TOL};

fn main() -> Result<(), LaError> {
    // Requires pivoting (a[0][0] = 0), so it's a good LU demo.
    let a = Matrix::<3>::from_rows([[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]);
    let b = Vector::<3>::new([5.0, 4.0, 3.0]);

    let lu = a.lu(DEFAULT_PIVOT_TOL)?;
    let x = lu.solve_vec(b)?.into_array();

    println!("x = {x:?}");
    Ok(())
}
```

## 🧩 API at a glance

| Type | Storage | Purpose | Key methods |
|---|---|---|---|
| `Vector<D>` | `[f64; D]` | Fixed-length vector | `new`, `zero`, `dot`, `norm2_sq` |
| `Matrix<D>` | `[[f64; D]; D]` | Fixed-size square matrix | `from_rows`, `zero`, `identity`, `lu`, `det` |
| `Lu<D>` | `Matrix<D>` + pivot array | Factorization for solves/det | `solve_vec`, `det` |

Storage shown above reflects the current `f64` implementation.

## 📋 Examples

The `examples/` directory contains small, runnable programs:

```bash
just examples
# or:
cargo run --example solve_3x3
cargo run --example det_3x3
```

## 🤝 Contributing

A short contributor workflow:

```bash
cargo install just
just ci           # lint + fast tests + bench compile
just commit-check # lint + all tests + examples
```

For the full set of developer commands, see `just --list` and `WARP.md`.

## 📄 License

BSD 3-Clause License. See [LICENSE](LICENSE).
