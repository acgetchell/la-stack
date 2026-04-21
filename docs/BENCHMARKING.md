# Benchmarking

This guide covers how to run, compare, and track performance for la-stack.

## Benchmark suites

la-stack has two Criterion benchmark suites:

- **`vs_linalg`** (`benches/vs_linalg.rs`) — compares la-stack against
  nalgebra and faer across D=2–64 for LU, solve, det, dot, norm, etc.
  Use this to answer "why choose la-stack over other crates?"

- **`exact`** (`benches/exact.rs`) — measures exact-arithmetic methods
  (`det_exact`, `solve_exact`, `det_sign_exact`, etc.) alongside f64
  baselines (`det`, `det_direct`) across D=2–5. Use this to understand
  the cost of exact arithmetic and track optimization progress.
  In addition to the per-dimension groups (`exact_d{2..5}`), the suite
  includes four adversarial-input groups designed to stress specific
  corners of the pipeline:
  - `exact_near_singular_3x3` — a 2^-50 perturbation of a singular base
    matrix; forces the Bareiss fallback in `det_sign_exact` and
    exercises the largest intermediate `BigInt` values in `solve_exact`.
  - `exact_large_entries_3x3` — diagonal entries near `f64::MAX / 2`
    stress `BigInt` growth during Bareiss forward elimination.
  - `exact_hilbert_4x4` / `exact_hilbert_5x5` — classically
    ill-conditioned matrices whose non-terminating-in-binary entries
    stress the `f64_decompose → BigInt` scaling path.
  Each adversarial group runs the same four benches (`det_sign_exact`,
  `det_exact`, `solve_exact`, `solve_exact_f64`) so the resulting tables
  are directly comparable across input classes.

## Quick reference

```bash
# Run vs_linalg benchmarks
just bench-vs-linalg

# Run exact-arithmetic benchmarks
just bench-exact

# Save an exact baseline (e.g., before optimising)
just bench-save-baseline v0.4.1

# Compare current code against a saved baseline
just bench-compare v0.4.1

# Generate a snapshot without comparison
just bench-compare
```

## Comparing performance across releases

Criterion baselines are saved into `target/criterion/` and persist across
`git checkout` but **not** across `cargo clean`.

### Workflow

```bash
# 1. Check out the old release and save its baseline
git checkout v0.2.0
just bench-save-baseline v0.2.0

# 2. Switch to current code and run benchmarks
git checkout main   # or your feature branch
just bench-exact    # populates target/criterion/*/new/

# 3. Generate a comparison table in docs/PERFORMANCE.md
just bench-compare v0.2.0
```

You can save multiple baselines and compare against any of them.

### Output

`just bench-compare` writes `docs/PERFORMANCE.md` (gitignored — it contains
machine-specific timings). The file includes per-dimension tables showing
median times, percent change, and speedup for each benchmark.

## vs\_linalg plotting

The `criterion_dim_plot.py` script generates CSV/SVG plots and updates the
README benchmark table from vs\_linalg results:

```bash
# Run benchmarks + update README table and SVG plot
just bench-vs-linalg
just plot-vs-linalg-readme
```

See `scripts/criterion_dim_plot.py --help` for options.

## Release workflow

At release time, save a baseline so future work can compare against it:

```bash
just bench-save-baseline $TAG
```

See `docs/RELEASING.md` step 5 for where this fits in the release process.
