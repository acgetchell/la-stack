# Benchmark Performance

**la-stack** v0.4.3 · `45affa8` (HEAD) · 2026-06-09 08:41:32 UTC
**Statistic**: median
**Suite**: all
**Scope**: release-signal

## Benchmark Results

Comparison against baseline **v0.4.2**:

Negative change = faster. Speedup > 1.00x = improvement.

## Exact arithmetic

### D=2

| Benchmark | v0.4.2 | Latest | Change | Speedup |
|-----------|-------:|-------:|-------:|--------:|
| det | 0.9 ns | 0.7 ns | **-24.0%** | 1.32x |
| det_direct | 1.0 ns | 1.0 ns | +2.1% | 0.98x |
| det_exact | 248.9 ns | 195.6 ns | **-21.4%** | 1.27x |
| det_exact_f64_result (vs det_exact_f64) | 429.1 ns | 167.6 ns | **-60.9%** | 2.56x |
| det_exact_rounded_f64 (vs det_exact_f64) | 429.1 ns | 375.2 ns | **-12.6%** | 1.14x |
| det_sign_exact | 1.5 ns | 3.2 ns | +115.9% | 0.46x |
| solve_exact | 6.53 µs | 6.45 µs | **-1.1%** | 1.01x |
| solve_exact_f64_result (vs solve_exact_f64) | 6.90 µs | 6.60 µs | **-4.4%** | 1.05x |
| solve_exact_rounded_f64 (vs solve_exact_f64) | 6.90 µs | 7.02 µs | +1.7% | 0.98x |

### D=3

| Benchmark | v0.4.2 | Latest | Change | Speedup |
|-----------|-------:|-------:|-------:|--------:|
| det | 1.8 ns | 1.5 ns | **-19.4%** | 1.24x |
| det_direct | 2.0 ns | 2.0 ns | +2.4% | 0.98x |
| det_exact | 739.0 ns | 468.6 ns | **-36.6%** | 1.58x |
| det_exact_f64_result (vs det_exact_f64) | 913.1 ns | 435.6 ns | **-52.3%** | 2.10x |
| det_exact_rounded_f64 (vs det_exact_f64) | 913.1 ns | 648.1 ns | **-29.0%** | 1.41x |
| det_sign_exact | 4.2 ns | 5.5 ns | +30.9% | 0.76x |
| solve_exact | 25.69 µs | 25.16 µs | **-2.1%** | 1.02x |
| solve_exact_f64_result (vs solve_exact_f64) | 26.16 µs | 25.42 µs | **-2.8%** | 1.03x |
| solve_exact_rounded_f64 (vs solve_exact_f64) | 26.16 µs | 25.67 µs | **-1.9%** | 1.02x |

### D=4

| Benchmark | v0.4.2 | Latest | Change | Speedup |
|-----------|-------:|-------:|-------:|--------:|
| det | 3.3 ns | 4.5 ns | +38.1% | 0.72x |
| det_direct | 3.7 ns | 4.3 ns | +17.6% | 0.85x |
| det_exact | 1.87 µs | 1.47 µs | **-21.8%** | 1.28x |
| det_exact_f64_result (vs det_exact_f64) | 2.04 µs | 1.47 µs | **-27.9%** | 1.39x |
| det_exact_rounded_f64 (vs det_exact_f64) | 2.04 µs | 1.63 µs | **-19.8%** | 1.25x |
| det_sign_exact | 6.9 ns | 11.5 ns | +67.1% | 0.60x |
| solve_exact | 64.95 µs | 61.67 µs | **-5.1%** | 1.05x |
| solve_exact_f64_result (vs solve_exact_f64) | 66.35 µs | 62.37 µs | **-6.0%** | 1.06x |
| solve_exact_rounded_f64 (vs solve_exact_f64) | 66.35 µs | 63.59 µs | **-4.2%** | 1.04x |

### D=5

| Benchmark | v0.4.2 | Latest | Change | Speedup |
|-----------|-------:|-------:|-------:|--------:|
| det | 26.0 ns | 23.3 ns | **-10.6%** | 1.12x |
| det_direct | 4.5 ns | 2.5 ns | **-44.2%** | 1.79x |
| det_exact | 4.10 µs | 4.05 µs | **-1.3%** | 1.01x |
| det_exact_f64_result (vs det_exact_f64) | 4.21 µs | 4.02 µs | **-4.4%** | 1.05x |
| det_exact_rounded_f64 (vs det_exact_f64) | 4.21 µs | 4.33 µs | +2.8% | 0.97x |
| det_sign_exact | 3.94 µs | 3.96 µs | +0.6% | 0.99x |
| solve_exact | 130.82 µs | 126.75 µs | **-3.1%** | 1.03x |
| solve_exact_f64_result (vs solve_exact_f64) | 132.70 µs | 127.37 µs | **-4.0%** | 1.04x |
| solve_exact_rounded_f64 (vs solve_exact_f64) | 132.70 µs | 128.15 µs | **-3.4%** | 1.04x |

### Near-singular 3x3

| Benchmark | v0.4.2 | Latest | Change | Speedup |
|-----------|-------:|-------:|-------:|--------:|
| det_sign_exact | 705.2 ns | 444.2 ns | **-37.0%** | 1.59x |
| det_exact | 724.0 ns | 478.9 ns | **-33.9%** | 1.51x |
| solve_exact | 3.44 µs | 3.39 µs | **-1.6%** | 1.02x |
| solve_exact_f64_result (vs solve_exact_f64) | 3.47 µs | 3.36 µs | **-3.2%** | 1.03x |
| solve_exact_rounded_f64 (vs solve_exact_f64) | 3.47 µs | 3.39 µs | **-2.5%** | 1.03x |

### Large entries 3x3

| Benchmark | v0.4.2 | Latest | Change | Speedup |
|-----------|-------:|-------:|-------:|--------:|
| det_sign_exact | 2.91 µs | 402.4 ns | **-86.2%** | 7.23x |
| det_exact | 2.94 µs | 434.0 ns | **-85.2%** | 6.76x |
| solve_exact | 82.81 µs | 81.57 µs | **-1.5%** | 1.02x |
| solve_exact_f64_result (vs solve_exact_f64) | 84.32 µs | 81.66 µs | **-3.1%** | 1.03x |
| solve_exact_rounded_f64 (vs solve_exact_f64) | 84.32 µs | 82.04 µs | **-2.7%** | 1.03x |

### Hilbert 4x4

| Benchmark | v0.4.2 | Latest | Change | Speedup |
|-----------|-------:|-------:|-------:|--------:|
| det_sign_exact | 6.9 ns | 11.5 ns | +66.4% | 0.60x |
| det_exact | 1.91 µs | 1.50 µs | **-21.7%** | 1.28x |
| solve_exact | 49.42 µs | 47.77 µs | **-3.3%** | 1.03x |
| solve_exact_f64_result (vs solve_exact_f64) | 50.38 µs | 47.67 µs | **-5.4%** | 1.06x |
| solve_exact_rounded_f64 (vs solve_exact_f64) | 50.38 µs | 48.17 µs | **-4.4%** | 1.05x |

### Hilbert 5x5

| Benchmark | v0.4.2 | Latest | Change | Speedup |
|-----------|-------:|-------:|-------:|--------:|
| det_sign_exact | 4.09 µs | 3.91 µs | **-4.6%** | 1.05x |
| det_exact | 4.00 µs | 4.02 µs | +0.6% | 0.99x |
| solve_exact | 98.71 µs | 95.41 µs | **-3.4%** | 1.03x |
| solve_exact_f64_result (vs solve_exact_f64) | 99.88 µs | 98.14 µs | **-1.7%** | 1.02x |
| solve_exact_rounded_f64 (vs solve_exact_f64) | 99.88 µs | 97.50 µs | **-2.4%** | 1.02x |

## How to Update

Local performance reports are generated in isolated temporary worktrees:

```bash
# Local development: compare the current tree with the latest release
just performance-local

# Release PR: update docs/PERFORMANCE.md and archive the previous report
just performance-release

# GitHub Actions release assets
just performance-github-assets

# Explicit repair
just performance-release <current-tag> <previous-tag>
```

`just performance-local` writes `target/bench-reports/performance.md`.
`just performance-github-assets` writes `target/bench-reports/github-assets-performance.md`.

Older curated release-to-release reports are archived in `docs/archive/performance/`.

See `docs/BENCHMARKING.md` for the full comparison workflow.
