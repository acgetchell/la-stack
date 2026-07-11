# scripts/

This directory contains Python utilities used during development of the `la-stack` crate.

## Setup

The Python 3.14 support tooling in this repo is managed with
[`uv`](https://github.com/astral-sh/uv) and resolved from `uv.lock`.

```bash
just python-sync
# or:
uv sync --locked --group dev
```

## Python maintenance rules

- Keep Python 3.14 code precisely typed and add focused pytest coverage for
  changed behavior and error paths.
- Mock subprocess results as `subprocess.CompletedProcess[str]`, matching the
  production boundary.
- Catch only specific, recoverable exceptions; do not use broad
  `except Exception` handlers.
- Give writer/parser pairs round-trip tests and explicit malformed-input tests.
- Update this README whenever Python entry points in `pyproject.toml` change.

## How to use it

### Comparing performance

The comparison script reads Criterion output and writes a local report to
`target/bench-reports/performance.md` by default:

```bash
# Re-render from existing Criterion output
just bench-compare
```

Use `uv run --locked bench-compare --snapshot` for a no-baseline snapshot, or
`uv run --locked bench-compare <baseline>` to compare against a named saved
baseline.

Use the top-level `just` workflows for routine release and local comparisons:

```bash
# Local development: compare the current tree with the latest release
just performance-local

# Release PR: update docs/PERFORMANCE.md and archive the previous report
just performance-release

# GitHub Actions release assets, without local cargo benchmark runs
just performance-github-assets
```

The local release workflows run the independent benchmark-input correctness gate
and then measure both library revisions with one hashed current benchmark
harness. Reports record source-state, environment, toolchain, dependency,
Criterion, harness, and validation provenance and fail on incomplete selected
coverage. Direct comparisons of separately published artifacts retain their
original per-release harnesses and label unavailable historical measurement
metadata explicitly.

See `docs/BENCHMARKING.md` for the current command matrix, local saved-baseline
workflow, explicit tag arguments, output locations, and release-artifact
comparison details.

### Plotting Criterion benchmarks (la-stack vs nalgebra/faer)

The plotter reads Criterion output under:

- `target/criterion/d{D}/{benchmark}/{new|base}/estimates.json`

And writes:

- `docs/assets/bench/vs_linalg_{metric}_{stat}.csv`
- `docs/assets/bench/vs_linalg_{metric}_{stat}.svg`
- `docs/assets/bench/vs_linalg_{metric}_{stat}.provenance.json`

To generate the single “time vs dimension” chart:

By default, the benchmark suite runs for dimensions 2–5, 8, 16, 32, and 64.

1. For exploratory plots, run the benchmarks you want to plot (this produces
   `target/criterion/...`):

```bash
# full run (takes longer, better for README plots)
just bench-vs-linalg lu_solve

# or quick run (fast sanity check; still produces estimates.json)
just bench-vs-linalg-quick lu_solve
```

2. Generate an exploratory chart (median or mean):

```bash
# median (recommended)
just plot-vs-linalg lu_solve median new true

# or mean
just plot-vs-linalg lu_solve mean new true
```

Use the dedicated publication path to update README's benchmark table (between
`BENCH_TABLE` markers):

```bash
just plot-vs-linalg-readme lu_solve median new true
```

That recipe runs the benchmark-input gate and a fresh full `vs_linalg` benchmark,
requires la-stack/nalgebra/faer results at every canonical dimension, and then
publishes CSV, SVG, JSON provenance, and README together. Partial dimensions are
available only through the plotter's explicit `--allow-partial` exploratory
option and cannot update README.

This writes:

- `docs/assets/bench/vs_linalg_lu_solve_median.csv`
- `docs/assets/bench/vs_linalg_lu_solve_median.svg` (requires `gnuplot`)

(For `stat=mean`, the filenames end in `_mean` instead of `_median`.)

### More examples

Plot a different metric:

```bash
uv run --locked criterion-dim-plot --metric dot --stat median --sample new
uv run --locked criterion-dim-plot --metric inf_norm --stat median --sample new
```

Plot a different statistic:

```bash
uv run --locked criterion-dim-plot --metric lu_solve --stat mean --sample new
```

Plot the previous (baseline) sample instead of the newest run:

```bash
uv run --locked criterion-dim-plot --metric lu_solve --stat median --sample base
```

Use a log-scale y-axis:

```bash
uv run --locked criterion-dim-plot --metric lu_solve --stat median --sample new --log-y
```

Write to custom output paths:

```bash
uv run --locked criterion-dim-plot \
  --metric lu_solve --stat median --sample new \
  --csv docs/assets/bench/custom.csv \
  --out docs/assets/bench/custom.svg
```

CSV only (skip SVG/gnuplot):

```bash
uv run --locked criterion-dim-plot --no-plot --metric lu_solve --stat median --sample new
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

## Changelog and release tooling

### Generating the changelog

```bash
# Full regeneration from all history
just changelog

# Prepend only unreleased changes for a new version
just changelog-unreleased v0.3.0
```

`just changelog` runs `git-cliff -o CHANGELOG.md`, strips trailing blank
lines, archives completed changelog series, and formats the generated Markdown.
Configuration lives in `cliff.toml` at the repo root.

### Creating a release tag

```bash
just tag v0.3.0          # create annotated tag from CHANGELOG.md section
just tag-force v0.3.0    # recreate tag if it already exists
```

The `tag-release` CLI (in `tag_release.py`) extracts the matching version
section from `CHANGELOG.md`, validates semver, and handles GitHub's 125KB
tag-annotation size limit.

### Scripts overview

| Script | Purpose |
|---|---|
| `archive_performance.py` | Promote release performance docs and archive older comparisons |
| `bench_compare.py` | Compare Criterion benchmark baselines and render Markdown reports |
| `criterion_dim_plot.py` | Plot Criterion benchmark results (CSV + SVG + README table) |
| `tag_release.py` | Create annotated git tags from CHANGELOG.md sections |
| `postprocess_changelog.py` | Strip trailing blank lines from git-cliff output |
| `subprocess_utils.py` | Safe subprocess wrappers for git commands |

See `docs/RELEASING.md` for the full release workflow.
