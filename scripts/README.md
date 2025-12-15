# scripts/

This directory contains Python utilities used during development of the `la-stack` crate.

## Setup

The Python tooling in this repo is managed with [`uv`](https://github.com/astral-sh/uv).

```bash
just python-sync
# or:
uv sync --group dev
```

## How to use it

### Plotting Criterion benchmarks (la-stack vs nalgebra)

The plotter reads Criterion output under:

- `target/criterion/d{D}/{benchmark}/{new|base}/estimates.json`

And writes:

- `docs/assets/bench/vs_nalgebra_{metric}_{stat}.csv`
- `docs/assets/bench/vs_nalgebra_{metric}_{stat}.svg`

To generate the single “time vs dimension” chart:

By default, the benchmark suite runs for dimensions 2–5, 8, 16, 32, and 64.

1. Run the benchmarks you want to plot (this produces `target/criterion/...`):

```bash
# full run (takes longer, better for README plots)
just bench-vs-nalgebra lu_solve

# or quick run (fast sanity check; still produces estimates.json)
just bench-vs-nalgebra-quick lu_solve
```

2. Generate the chart (median or mean):

```bash
# median (recommended)
just plot-vs-nalgebra lu_solve median new true

# median + update README's benchmark table (between BENCH_TABLE markers)
just plot-vs-nalgebra-readme lu_solve median new true

# or mean
just plot-vs-nalgebra lu_solve mean new true
```

This writes:

- `docs/assets/bench/vs_nalgebra_lu_solve_median.csv`
- `docs/assets/bench/vs_nalgebra_lu_solve_median.svg` (requires `gnuplot`)

(For `stat=mean`, the filenames end in `_mean` instead of `_median`.)

### More examples

Plot a different metric:

```bash
uv run criterion-dim-plot --metric dot --stat median --sample new
uv run criterion-dim-plot --metric inf_norm --stat median --sample new
```

Plot a different statistic:

```bash
uv run criterion-dim-plot --metric lu_solve --stat mean --sample new
```

Plot the previous (baseline) sample instead of the newest run:

```bash
uv run criterion-dim-plot --metric lu_solve --stat median --sample base
```

Use a log-scale y-axis:

```bash
uv run criterion-dim-plot --metric lu_solve --stat median --sample new --log-y
```

Write to custom output paths:

```bash
uv run criterion-dim-plot \
  --metric lu_solve --stat median --sample new \
  --csv docs/assets/bench/custom.csv \
  --out docs/assets/bench/custom.svg
```

CSV only (skip SVG/gnuplot):

```bash
uv run criterion-dim-plot --no-plot --metric lu_solve --stat median --sample new
```

### gnuplot

SVG rendering requires `gnuplot` to be installed and available on `PATH`.

Install (macOS/Homebrew):

```bash
brew install gnuplot
```

Verify the installed version:

```bash
gnuplot --version
```

This repo has been tested with `gnuplot 6.0 patchlevel 3` (Homebrew `gnuplot 6.0.3`).
