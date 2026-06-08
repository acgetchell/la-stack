# Exact Arithmetic Performance

**la-stack** v0.4.2 · `7e11f93` (HEAD) · 2026-06-08 20:39:03 UTC
**Statistic**: median

## Benchmark Results

Comparison against baseline **v0.4.1**:

Negative change = faster. Speedup > 1.00x = improvement.

### D=2

| Benchmark | v0.4.1 | Current | Change | Speedup |
|-----------|-------:|--------:|-------:|--------:|
| det | 0.6 ns | 0.9 ns | +61.1% | 0.62x |
| det_direct | 0.7 ns | 1.0 ns | +44.7% | 0.69x |
| det_exact | 315.5 ns | 318.4 ns | +0.9% | 0.99x |
| det_exact_f64 | 555.7 ns | 555.7 ns | -0.0% | 1.00x |
| det_sign_exact | 0.7 ns | 1.5 ns | +128.2% | 0.44x |
| solve_exact | 7.05 µs | 7.06 µs | +0.2% | 1.00x |
| solve_exact_f64 | 7.50 µs | 7.67 µs | +2.3% | 0.98x |

### D=3

| Benchmark | v0.4.1 | Current | Change | Speedup |
|-----------|-------:|--------:|-------:|--------:|
| det | 1.3 ns | 1.8 ns | +36.3% | 0.73x |
| det_direct | 4.7 ns | 2.2 ns | **-51.9%** | 2.08x |
| det_exact | 936.9 ns | 924.3 ns | **-1.3%** | 1.01x |
| det_exact_f64 | 1.18 µs | 1.19 µs | +1.1% | 0.99x |
| det_sign_exact | 2.4 ns | 4.2 ns | +78.1% | 0.56x |
| solve_exact | 27.02 µs | 27.41 µs | +1.5% | 0.99x |
| solve_exact_f64 | 28.06 µs | 27.98 µs | -0.3% | 1.00x |

### D=4

| Benchmark | v0.4.1 | Current | Change | Speedup |
|-----------|-------:|--------:|-------:|--------:|
| det | 2.4 ns | 3.3 ns | +36.8% | 0.73x |
| det_direct | 2.4 ns | 4.1 ns | +70.2% | 0.59x |
| det_exact | 2.33 µs | 2.33 µs | -0.0% | 1.00x |
| det_exact_f64 | 2.59 µs | 2.58 µs | -0.7% | 1.01x |
| det_sign_exact | 5.3 ns | 6.9 ns | +30.5% | 0.77x |
| solve_exact | 67.14 µs | 67.99 µs | +1.3% | 0.99x |
| solve_exact_f64 | 67.86 µs | 68.51 µs | +1.0% | 0.99x |

### D=5

| Benchmark | v0.4.1 | Current | Change | Speedup |
|-----------|-------:|--------:|-------:|--------:|
| det | 21.6 ns | 24.5 ns | +13.7% | 0.88x |
| det_direct | 2.3 ns | 4.7 ns | +104.8% | 0.49x |
| det_exact | 5.04 µs | 4.99 µs | -1.0% | 1.01x |
| det_exact_f64 | 5.32 µs | 5.31 µs | -0.1% | 1.00x |
| det_sign_exact | 4.97 µs | 4.99 µs | +0.3% | 1.00x |
| solve_exact | 134.99 µs | 136.04 µs | +0.8% | 0.99x |
| solve_exact_f64 | 137.11 µs | 138.97 µs | +1.4% | 0.99x |

### Near-singular 3x3

| Benchmark | v0.4.1 | Current | Change | Speedup |
|-----------|-------:|--------:|-------:|--------:|
| det_sign_exact | 871.8 ns | 877.6 ns | +0.7% | 0.99x |
| det_exact | 907.3 ns | 904.4 ns | -0.3% | 1.00x |
| solve_exact | 4.31 µs | 4.25 µs | **-1.5%** | 1.02x |
| solve_exact_f64 | 4.29 µs | 4.32 µs | +0.7% | 0.99x |

### Large entries 3x3

| Benchmark | v0.4.1 | Current | Change | Speedup |
|-----------|-------:|--------:|-------:|--------:|
| det_sign_exact | 3.14 µs | 3.09 µs | **-1.3%** | 1.01x |
| det_exact | 3.19 µs | 3.11 µs | **-2.3%** | 1.02x |
| solve_exact | 84.77 µs | 83.89 µs | **-1.0%** | 1.01x |
| solve_exact_f64 | 84.62 µs | 83.92 µs | -0.8% | 1.01x |

### Hilbert 4x4

| Benchmark | v0.4.1 | Current | Change | Speedup |
|-----------|-------:|--------:|-------:|--------:|
| det_sign_exact | 5.3 ns | 6.9 ns | +30.4% | 0.77x |
| det_exact | 2.39 µs | 2.31 µs | **-3.2%** | 1.03x |
| solve_exact | 51.69 µs | 52.27 µs | +1.1% | 0.99x |
| solve_exact_f64 | 52.90 µs | 53.26 µs | +0.7% | 0.99x |

### Hilbert 5x5

| Benchmark | v0.4.1 | Current | Change | Speedup |
|-----------|-------:|--------:|-------:|--------:|
| det_sign_exact | 5.03 µs | 4.88 µs | **-2.9%** | 1.03x |
| det_exact | 5.07 µs | 4.96 µs | **-2.1%** | 1.02x |
| solve_exact | 105.35 µs | 102.72 µs | **-2.5%** | 1.03x |
| solve_exact_f64 | 104.99 µs | 103.94 µs | -1.0% | 1.01x |

## How to Update

Release performance docs are generated in isolated temporary worktrees:

```bash
# Release PR: update docs/PERFORMANCE.md and archive the previous report
just performance-release <current-tag> <previous-tag>

# Historical published comparison
just performance-archive-published

# Explicit historical repair
just performance-archive-published <current-tag> <previous-tag>
```

For local scratch comparisons, use `just bench-latest` and `just bench-compare`.
Those write `target/bench-reports/performance.md`.

See `docs/BENCHMARKING.md` for the full comparison workflow.
