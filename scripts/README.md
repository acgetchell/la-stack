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

Run `just update` for the deliberate shared maintenance workflow. Tool updates
run first: uv upgrades through its installation owner, declared managed Cargo
tools upgrade with verified TOML pins, and setup synchronizes the environment.
Just follows the shared package's `rust-just` dependency. Dependency updates then
advance Cargo requirements and exact direct Python `dev` pins, refresh both
locks, and explicitly synchronize `dev` with managed Rust available. A failed
step stops subsequent work without rolling back earlier package-manager steps.
The coupled `num-bigint` and `num-rational` requirements are excluded from
independent incompatible upgrades and must be advanced together. The Python
shared updater asks uv to resolve one cross-platform tool set before applying all
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

# Release PR: publish docs/performance.md and retain the complete run
just performance-release

# Build release docs from retained complete-run evidence (or historical inputs)
just performance-doc

# GitHub Actions release assets, without local cargo benchmark runs
just performance-github-assets
```

Local benchmark generation streams native Just, Cargo, and Criterion progress
while retaining the existing fail-closed report and provenance checks. Staged and
unstaged changes to tracked files participate. Untracked files are excluded;
stage a new file before running the command if it must participate. A local
current-vs-latest report may use the same package and release identifier because
its commit/ref and source-state provenance still distinguish the revisions.

The local release workflows run the independent benchmark-input correctness gate
and then measure both library revisions with one hashed current benchmark
harness. Retained evidence records source-state, environment, toolchain, dependency,
Criterion, harness, and validation provenance and fails on incomplete selected
coverage. `performance-local` writes `performance.md`, `performance.run.json`,
and `performance.evidence.json` under `target/bench-reports/`. The shared
complete-run payload preserves every semantic Criterion case, 100 raw samples,
mean and median estimates, and 95% intervals from both phases.
`performance-release` requires distinct releases and publishes immutable
`run.json`, `evidence.json`, and `report.md` files under
`docs/performance-runs/runs/<content-id>/`, with a validated index and latest
pointer. The shared full report is published to `docs/performance.md` in the
same transaction, using the path and title in `tooling/performance-report.toml`.
Its mean and median tables show named series, marginal intervals, and unavailable
measurements, with each series attributed to its original phase.
Full provenance remains in the evidence envelope. Repeated runs for one pair coexist; publication failures
preserve the previous reports and selection. `performance-doc` replays this
evidence without Cargo. Same-version local comparisons require distinct source
identities and cannot be promoted as release reports.

After `target/` cleanup, report and README commands select validated shared
history. Any partial scratch bundle fails instead of selecting older evidence.
Legacy CSV/JSON bytes, hashes, and rendered reports remain readable through
their original schema under `docs/performance/`; new runs never overwrite or
reinterpret that history. Native Criterion release archives remain the
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

The README recipe consumes the complete-run payload and evidence envelope,
or the validated latest committed run after cleanup;
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

The shared `research-repo-tools release update` CLI infers the previous stable
release from published GitHub releases and updates package, lockfile, citation,
and release metadata transactionally. Active README navigation remains on `main`;
measured artifact links retain their recorded revision. It records the current UTC
date in `CITATION.cff`. Pass `--previous-release vA.B.C` to avoid release discovery,
`--date YYYY-MM-DD` to declare a date, or `--dry-run` to preview validated changes:

```bash
just update-version vX.Y.Z --previous-release vA.B.C --date YYYY-MM-DD --dry-run
```

Consumer policy in `pyproject.toml` checks the concept DOI and the expected README
reference count. Measured benchmark links, retained performance reports, archived
documentation, and dependency versions stay unchanged. Benchmark command examples
use selected tag variables rather than mutable release literals. If the target
changelog heading already exists, its date advances with the citation date.
After metadata preparation, `just docs-version-check` requires generated notes
for the new version; generate them with the same declared date before validation.

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

The exact published `research-repo-tools==0.1.8` dependency is in `tooling`,
included by `dev`, and resolved from PyPI in `uv.lock`. Thin Just recipes use
the documented `research-repo-tools changelog` CLI. Its implementation and
common regressions belong to the shared package; this repository owns the
configuration, recipes, and `scripts/tests/test_changelog_integration.py`.
Internal shared modules are not a supported consumer API.

The shared `deps update-python` command advances only exact requirements
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

The shared changelog table in `pyproject.toml` declares the repository owner,
repository name, formatter, and `dependency-bodies = "preserve"`. The installed
v0.1.8 common template owns grouping, tag grammar, links, and Markdown handling;
there is no local `cliff.toml` or template exception. Installed-package tests
exercise Ruff, Ty, and setuptools authored text, links, and fenced code across
repeat generation, and check retained historical notes.
The shared normalizer preserves complete breaking-change descriptions, Markdown
links, and literal code in the input while adding merged-PR summaries.
It retains embedded conventional headings that the old normalizer deduplicated.

`changelog-rumdl.toml` inherits the repository Markdown policy and disables
MD057 only while formatting unpublished candidates, whose future archive paths
do not yet exist. Repository Markdown checks and focused archive-link checks
still validate published destinations. The archive move preserves all release
dates and links; generated whitespace changes align existing series with the
shared formatter so later generation remains conflict-free.

Setup, tool verification, managed execution, and updates now use the published
shared CLI. Consumer declarations stay in `.python-version`, `rust-toolchain.toml`,
and `pyproject.toml`; update recipes retain coupled Cargo exclusions and explicit
dev synchronization. `tests/test_toolchain_integration.py` exercises the actual
recipes with recording update boundaries and verifies native managed execution.
The shared CLI also owns release metadata and version checks, Markdown line
checks, and Semgrep fixture validation. `tests/test_maintenance_integration.py`
checks the actual release selectors, preserved scientific evidence, DOI policy,
and recipe forwarding. `tests/test_cargo_update_integration.py` executes native
Cargo upgrades against a disposable local registry; Python updates have a
matching real-uv fixture in the toolchain tests. Common parser, transaction,
Markdown, and fixture regressions belong to the shared package. Scientific
eligibility, Cargo commands/features, release adapters, selected rows, historical
report layouts, and dimension plots remain consumer-owned. The shared package
owns complete-run reports, common-harness phases, completeness, worktrees, run
identities, retention, and publication.
Notebook tooling remains outside scope.

The same pinned release owns opt-in CodeRabbit review orchestration through
`research-repo-tools review branch --base=REF` and `review uncommitted`.
Thin Just wrappers retain the common implementation upstream; consumer checks
in `tests/test_review_integration.py` exercise recipe forwarding, argument
quoting, and failure propagation with local stubs.
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
| `archive_performance.py` | Select eligible releases and compose shared measurement/publication APIs |
| `bench_compare.py` | Compare Criterion benchmark baselines and render Markdown reports |
| `benchmark_contract.py` | Select the Cargo checkout and hash the historical benchmark inventory contract |
| `benchmark_summaries.py` | Read historical complete CSV/JSON summaries without changing framing |
| `criterion_dim_plot.py` | Plot Criterion benchmark results (CSV + SVG + README table) |
| `performance_runs.py` | Declare scientific run policy and adapt retained series for dimension plots |
| `performance_artifacts.py` | Validate and publish schema-versioned performance-comparison CSV/JSON inputs |
| `release_baseline.py` | Inventory full Criterion suites and validate complete raw release baselines before packaging |

Shared process discovery, execution, byte transport, CPU detection, diagnostics,
and zizmor authentication belong to research-repo-tools. Callers use its public
process API directly. Shared Criterion estimates and comparisons also replace
local timing wrappers and numerical validators. Consumer tests use the shared
Just inspection API and retain benchmark contracts, caller file-policy coverage,
and recipe forwarding; common parser and review regressions belong upstream.

The `tooling/` directory owns declarative performance configuration:

- `performance.toml`: shared measurement inputs, harness files, tool/dependency
  probes, timeout, and provenance compatibility.
- `performance.just`: native Cargo inventory, independent gates, and timing
  recipes, including baseline-only peer measurements.
- `performance-report.toml`: shared report and immutable-history paths.

Both local measurement and release inventory use the native recipes. The shared
runner streams them in the measured checkout, with the recipes and measurement
configuration included in the harness fingerprint. Scientific row selection,
sampling requirements, historical API adapters, and figure layout remain Python.
Old evidence naming the retired Python phase driver can still be validated and
rendered offline; its recorded commands are never executed.

The superseded `criterion_measurements.py`, archive/worktree measurement engine,
summary retention writer, and their duplicate generic tests are removed.
`test_performance_workflow.py` exercises the installed package through consumer
policies, failure preservation, historical golden bytes, and offline replay.
The shared multi-series plotting extension is deferred, so the existing local
CSV/SVG/README adapter preserves figure format and baseline-phase peer labels.
Complete-run Markdown uses the shared renderer directly. Upstream follow-ups
cover [coordinate plots](https://github.com/acgetchell/research-repo-tools/issues/95),
[provenance summaries](https://github.com/acgetchell/research-repo-tools/issues/96), and
[complete-run document publication](https://github.com/acgetchell/research-repo-tools/issues/97).

Performance scripts also use the published shared Criterion parser and estimate
validation, comparison arithmetic, exact-byte digest verification, safe archive
extraction, document marker replacement, and multi-file publication transaction.
Raw Criterion numeric strings are rejected by the shared parser. Plotting and
release publication still require complete confidence intervals; hosted baselines
and full local summaries additionally require 100 samples and 95% intervals.
README marker replacement preserves bytes outside the selected section.
Publication candidates are fully rendered and validated before one transaction
replaces the report, evidence, archive index, and retained summaries. Shared
recovery errors identify preserved backups if rollback fails.

The remaining modules own la-stack's historical CSV/JSON schemas and fingerprints,
scientific eligibility, benchmark inventories, Cargo commands, release adapters,
release selection policy, and multi-library plot layout.
Existing evidence keeps its schema and digest framing. Generic parser and
transaction regressions belong upstream; consumer tests check retained-schema
round trips, complete output groups, failure preservation, and scientific policy.
The published v0.1.8 contract owns common-harness installation, independent
gates, measurement freshness, configurable completeness, and immutable retention.

See `docs/RELEASING.md` for the full release workflow.
