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

### Updating dependencies and repository-owned tools

Run `just update` for the deliberate maintenance workflow. It updates Cargo and
exact Python development-tool declarations and their locks, upgrades only the
Cargo CLI packages owned by `setup-tools`, and then reconciles their installed
versions plus the active uv version with the root `justfile` atomically. All
required update tools are checked before the first dependency write, and the
maintenance workflow accepts a newer active uv so it can become the new pin.
The coupled `num-bigint` and `num-rational` requirements are excluded from
independent incompatible upgrades and must be advanced together. The Python
updater asks uv to resolve one cross-platform tool set before applying all
changed exact pins together; it does not change runtime or build-system
requirements.

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

New release comparisons require v0.4.4 or newer on both sides. v0.4.3 and older
are excluded; historical Markdown reports remain available in the archive.

```bash
# Local development: compare the current tree with the latest release
just performance-local

# Release PR: update docs/performance.md and archive the previous report
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
`target/bench-reports/`, plus `performance.full.csv` and its provenance JSON
containing every recorded case from both measurement phases.
`performance-release` does the same measurement work, requires distinct releases,
and preserves all four files under `docs/performance/<release-pair>/<run-digest>/`
when promoting the validated result. `performance-doc` consumes the retained pair
from either workflow without Cargo or temporary worktrees, then promotes the
result into `docs/performance.md` and the archive. Same-version local artifacts
remain valid comparison evidence but cannot be promoted as a release report.
Scratch files may be removed with `target/`. After cleanup, `performance-doc` and
`performance-readme` resolve `docs/performance/latest.json` to the most recently
promoted complete snapshot. Existing partial scratch inputs fail validation.
The [local summary index](../docs/performance/README.md) owns the saved-data
schema and retention contract. Native Criterion release archives remain the
durable raw baselines from the GitHub runner. Direct comparisons
of separately published artifacts retain their original per-release harnesses
and label unavailable historical measurement metadata explicitly.

The comparison suite lives in the unpublished `benches/comparison` package, so
exact benchmark builds do not compile nalgebra/faer. Local release comparisons
measure la-stack on each revision and reuse peer measurements from the baseline
phase. The filtered current run skips Criterion HTML generation to avoid reading
missing `new` samples for peers that were not rerun.

Operationally, `performance-release` is the atomic composition of
`performance-local` and `performance-doc`: measure, retain the common
comparison inputs, render, and promote. The narrowed
`performance-local-non-exact` view uses the same metrics and renderer with
nalgebra/faer context, but writes `performance-non-exact.*` so it does not
overwrite the canonical full comparison bundle.

See `docs/BENCHMARKING.md` for the current command matrix, local saved-baseline
workflow, explicit tag arguments, output locations, and release-artifact
comparison details.

### Plotting la-stack vs nalgebra/faer benchmarks

The exploratory plotter reads Criterion output under:

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
# full run (takes longer, better for exploratory plots)
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

For release publication, first retain the validated release comparison, then
update README's benchmark table (between `BENCH_TABLE` markers):

```bash
just performance-release
just performance-readme lu_solve median new true
```

The README recipe consumes `target/bench-reports/performance.csv` and its
adjacent provenance JSON, or the latest committed snapshot after cleanup;
it does not run the benchmark-input gate or Criterion
again. It uses the current la-stack timing and retained same-current-harness
nalgebra/faer timings, requires all three at every canonical dimension, and then
publishes CSV, SVG, derived JSON provenance, the README table, and its tag-pinned
benchmark artifact links together. Until publication succeeds, those links keep
referencing the previous published artifacts. Partial dimensions are available
only through the raw-Criterion plotter's explicit `--allow-partial` exploratory
option and cannot update README.

This writes:

- `docs/assets/bench/vs_linalg_lu_solve_median.csv`
- `docs/assets/bench/vs_linalg_lu_solve_median.svg` (requires `gnuplot`)
- `docs/assets/bench/vs_linalg_lu_solve_median.provenance.json`
- the benchmark table and tag-pinned asset links in `README.md`

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

### Updating release metadata

```bash
just update-version vX.Y.Z
```

The updater infers the previous stable release from published GitHub releases,
updates package, lockfile, citation, non-artifact README, and active
benchmark workflow version references transactionally, and records the current
UTC date in `CITATION.cff`. It leaves README benchmark artifact links for
`performance-readme`, and it does not upgrade dependencies. If the target
changelog heading already exists, the updater advances its date atomically with
the citation date.

### Generating the changelog

```bash
# Full regeneration from all history
just changelog

# Preview without replacing the root or archive files
just changelog-preview

# Generate a prospective release with an explicit ISO date
just changelog-unreleased vX.Y.Z YYYY-MM-DD
# Equivalent name:
just changelog-release vX.Y.Z YYYY-MM-DD

# Rotate existing notes without regenerating Git history
just changelog-archive
just changelog-check
just release-notes vX.Y.Z
```

The exact published `research-repo-tools==0.1.7` dependency is in `tooling`,
included by `dev`, and resolved from PyPI in `uv.lock`. Thin Just recipes use
the documented `research-repo-tools changelog` CLI. Its implementation and
common regressions belong to the shared package; this repository owns the
configuration, recipes, and `scripts/tests/test_changelog_integration.py`.
Internal shared modules are not a supported consumer API.

The consumer-owned `update-python` helper advances only exact requirements
declared directly in `dev`. Included groups retain their own upgrade policy;
the shared `tooling` pin and its lockfile change together through an intentional
dependency upgrade.

`just changelog` generates, normalizes, formats, and rotates completed minor
series in one operation. It retains Unreleased and the newest minor series in
`CHANGELOG.md`, with older releases in `docs/archives/changelog/MAJOR.MINOR.md`.
Existing dated headings remain authoritative. Prospective generation requires
a `v`-prefixed SemVer tag and an explicit ISO date, and leaves package and
citation metadata untouched. Prepare metadata separately with
`just update-version`; choose the same date recorded in `CITATION.cff`.
Preview validates root and archive candidates without publishing them.
Malformed versions/dates, duplicate releases, conflicting retained notes,
and formatter failures stop publication with diagnostics.

`cliff.toml` is materialized from the installed package's common 0.1.7 template.
It adopts the shared categories, SemVer tag grammar, compare links, and
Markdown-aware code handling. Its only semantic exception retains bodies for
`deps-dev` commits when grouped as Dependencies, preserving historical Ruff,
Ty, and setuptools release-note, changelog, and comparison links. A focused
integration check compares the parsed configuration with the installed template
plus that single condition, preventing unrelated local policy drift.
The exception can retire after adopting a published shared solution to
[research-repo-tools#62](https://github.com/acgetchell/research-repo-tools/issues/62).
The shared normalizer preserves complete breaking-change descriptions, Markdown
links, and literal code in the input while adding merged-PR summaries.
It retains embedded conventional headings that the old normalizer deduplicated.

`changelog-rumdl.toml` inherits the repository Markdown policy and disables
MD057 only while formatting unpublished candidates, whose future archive paths
do not yet exist. Repository Markdown checks and focused archive-link checks
still validate published destinations. The archive move preserves all release
dates and links; generated whitespace changes align existing series with the
shared formatter so later generation remains conflict-free.

The current local setup, release-metadata, dependency-update, scientific,
benchmark, and performance tooling remains consumer-owned. Shared toolchain
setup, updates, and other maintenance adoption belong in later PRs. Notebook
tooling remains outside this repository's current scope.

The same pinned release owns opt-in CodeRabbit review orchestration through
`research-repo-tools review branch --base=REF` and `review uncommitted`.
Thin Just wrappers retain the common implementation upstream; consumer checks
in `tests/test_review_integration.py` exercise recipe forwarding, instruction
discovery, freshness diagnostics, and failure propagation with local stubs.
CodeRabbit remains externally installed and authenticated. See the
[contributor review workflow](../CONTRIBUTING.md#coderabbit-review) for scopes
and invocation policy.

### Creating a release tag

```bash
just tag vX.Y.Z          # create an annotated tag matching Cargo.toml
just tag-force vX.Y.Z    # replace that tag only when explicitly repairing it
```

These recipes call `research-repo-tools changelog tag`. It searches the root
and canonical archives, validates the entire history, and requires the
`v`-prefixed tag to match the Cargo package version. The consumer's `declared`
date policy requires the dated release heading to match any `CITATION.cff` date,
so a release prepared before merge can be tagged on a later UTC day. Oversized annotations
link to the full notes to respect GitHub's 125KB limit. `--force` replaces an
existing local ref only after validation; it does not first delete the old tag.
Tagging never pushes or publishes a release. Use the CLI's `--dry-run` to
preview an annotation without creating a tag.

### Scripts overview

| Script | Purpose |
|---|---|
| `archive_performance.py` | Promote release performance docs and archive older comparisons |
| `performance_artifacts.py` | Validate and publish schema-versioned performance-comparison CSV/JSON inputs |
| `bench_compare.py` | Compare Criterion benchmark baselines and render Markdown reports |
| `benchmark_summaries.py` | Preserve every local case summary, bind provenance, and resolve saved report inputs after cleanup |
| `check_docs_version_sync.py` | Verify versioned documentation links and snippets stay synchronized |
| `criterion_dim_plot.py` | Plot Criterion benchmark results (CSV + SVG + README table) |
| `criterion_measurements.py` | Validate full Criterion sampling and estimates for local summaries and hosted archives |
| `release_baseline.py` | Inventory full Criterion suites and validate complete raw release baselines before packaging |
| `subprocess_utils.py` | Safe subprocess wrappers for git commands |
| `update_cargo_tool_pins.py` | Reconcile repository-owned Cargo and active uv tool pins with installed versions |
| `update_python_dev_pins.py` | Resolve and advance exact Python development-tool pins through uv |
| `update_release_version.py` | Transactionally update deterministic release-version metadata |

See `docs/RELEASING.md` for the full release workflow.
