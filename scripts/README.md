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

`performance-github-assets` always requires the GitHub CLI (`gh`) authenticated
for the repository because it downloads release assets. Local-generation
recipes require `gh` only when discovering published release tags.

```bash
# Local development: compare the current tree with the latest release
just performance-local

# Release PR: update docs/PERFORMANCE.md and archive the previous report
just performance-release

# Build release docs from retained CSV/JSON inputs
just performance-doc

# GitHub Actions release assets, without local cargo benchmark runs
just performance-github-assets
```

Local benchmark generation streams Cargo and Criterion progress while retaining
the existing fail-closed report and provenance checks. Lines prefixed with
`[performance]` identify the active validation or timing phase. Staged and
unstaged changes to tracked files participate. Untracked files are excluded;
stage a new file before running the command if it must participate. A local
current-vs-latest report may use the same package and release identifier because
its commit/ref and source-state provenance still distinguish the revisions.

The local release workflows run the independent benchmark-input correctness gate
and then measure both library revisions with one hashed current benchmark
harness. Reports record source-state, environment, toolchain, dependency,
Criterion, harness, and validation provenance and fail on incomplete selected
coverage. `performance-local` writes Markdown plus schema-versioned
`performance.csv` and `performance.provenance.json` inputs under
`target/bench-reports/` without promoting documentation. `performance-release`
does the same measurement and retention work, requires distinct releases, and
promotes the validated result. `performance-doc` consumes the retained pair
from either workflow without Cargo or temporary worktrees, then promotes the
result into `docs/PERFORMANCE.md` and the archive. Same-version local artifacts
remain valid comparison evidence but cannot be promoted as a release report.
These files are reproducible scratch and may be removed with `target/`; native
Criterion release archives remain the durable raw baselines. Direct comparisons
of separately published artifacts retain their original per-release harnesses
and label unavailable historical measurement metadata explicitly.

Operationally, `performance-release` is the atomic composition of
`performance-local` and `performance-doc`: measure, retain the common
comparison inputs, render, and promote. The narrowed
`performance-local-non-exact` view uses the same metrics and renderer with
nalgebra/faer context, but writes `performance-non-exact.*` so it does not
overwrite the canonical full comparison bundle.

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
just changelog-unreleased vX.Y.Z
```

`just changelog` runs `git-cliff -o CHANGELOG.md`, strips trailing blank
lines, archives completed changelog series, and formats the generated Markdown.
Configuration lives in `cliff.toml` at the repo root.

### Creating a release tag

```bash
just tag vX.Y.Z          # create an annotated tag matching Cargo.toml
just tag-force vX.Y.Z    # replace that tag only when explicitly repairing it
```

The `tag-release` CLI (in `tag_release.py`) extracts the matching version
section from `CHANGELOG.md`, requires the tag to match the Cargo package version,
validates SemVer, and handles GitHub's 125KB tag-annotation size limit.

### Scripts overview

| Script | Purpose |
|---|---|
| `archive_changelog.py` | Split completed changelog minor series into archives |
| `archive_performance.py` | Promote release performance docs and archive older comparisons |
| `performance_artifacts.py` | Validate and publish schema-versioned performance-comparison CSV/JSON inputs |
| `bench_compare.py` | Compare Criterion benchmark baselines and render Markdown reports |
| `check_docs_version_sync.py` | Verify versioned documentation links and snippets stay synchronized |
| `criterion_dim_plot.py` | Plot Criterion benchmark results (CSV + SVG + README table) |
| `tag_release.py` | Create annotated git tags from CHANGELOG.md sections |
| `postprocess_changelog.py` | Normalize and reflow generated git-cliff Markdown safely |
| `subprocess_utils.py` | Safe subprocess wrappers for git commands |

See `docs/RELEASING.md` for the full release workflow.
