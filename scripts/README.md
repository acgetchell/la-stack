# scripts/

Only two production Python adapters remain. Generic tools run from the installed
`research-repo-tools==0.1.9` base package; their declarations belong in
`tooling/` and `pyproject.toml`. The project is a dependency-only environment,
with no local Python distribution, build backend or console entry points.

## Setup

Use `just setup`, then `uv run --locked` for every Python command. The support
environment targets Python 3.15. `just python-check` checks formatting, lint and
types; `just test-python` runs focused installed-package consumer tests.

Semgrep 1.180.0 has a separate locked Python 3.14 environment under
`tooling/semgrep/`, because its published wheels do not support Python 3.15.
Both local scans and hosted SARIF use that environment. Review its manifest and
refresh its lock deliberately when updating the scanner.

## Performance tooling

The [benchmark guide](../docs/BENCHMARKING.md) owns command selection, scientific
interpretation, output paths and release procedures.

| Owner | Responsibility |
| --- | --- |
| `archive_performance.py` | Release eligibility, tracked-only snapshots, historical raw-archive composition, and la-stack README markers/pointers |
| `performance_runs.py` | Historical API cfg/environment, scientific coverage, sampling and original reference-phase assertions |
| `tooling/performance-policy.toml` | Exact group selections, metric names, peers, dimensions and historical API availability |
| `tooling/performance.toml` | Shared sources, harness, probes, dependencies, timeouts and compatibility |
| `tooling/performance.just` | Native Cargo discovery, independent correctness gates and timing commands |
| `tooling/performance-report.toml` | Shared current report and immutable complete-run history |
| `tooling/performance-plot.toml` | Scientific coordinate, unit, series and axis declarations |
| `tooling/criterion-inventory.toml` | Shared native discovery and minimum coverage of all report/plot consumers |
| `tooling/criterion-saved.toml` | Shared complete saved-sample policy, including exact `tukey.json` bytes |

The adapters call public shared APIs for subprocesses, worktrees, Criterion
discovery, completeness, estimates, comparison arithmetic, evidence, rendering,
retention, publication, archive extraction and authenticated downloads. They
contain no generic renderer, CSV parser, legacy writer or measurement engine.

The remaining composition is required by la-stack's contracts: versions before
v0.4.4 are ineligible; v0.4.4/v0.4.5 need the captured old-API cfg and root-test
layout adapter; local measurements exclude untracked files; nalgebra/faer
reference views retain their baseline origin. The shared CLI's default worktree
snapshot includes untracked files and its declarations cannot select an adapter
from a release pair. Historical release archives lack the metadata required by
the shared asset-evidence CLI, so their small adapter renders descriptive
comparisons without inventing source/host/gate provenance.

`bench-compare` now calls the shared Criterion CLI directly. It compares every
available row, including additions and omissions, as an exploratory local report.
Its former suite/scope flags and snapshot CLI are retired. Use Criterion's native
reports for a single snapshot and `performance-local` for checked scientific
release comparisons.

`performance-local` and `performance-release` retain complete native measurements.
`performance-doc` renders those inputs or validated shared history offline.
`performance-readme` and `plot-vs-linalg` use the shared complete-run SVG, CSV and
table renderers. The former preserves current source/harness membership and
digests and advances the three existing measured-artifact pointers transactionally.
README candidates pass the shared Markdown formatter with repository policy
before any publication files are replaced.
The README adapter labels dimensions, keeps each confidence interval below its
timing, states the confidence level once, and selects the shared speedup columns.
Timing estimates and comparison arithmetic remain owned by the shared renderer.

```bash
just performance-release
just performance-readme lu_solve median new true
just plot-vs-linalg dot mean new true
```

The raw-Criterion/gnuplot publisher, historical summary writer and custom
report renderer are retired after installed-package replacement checks.
Original CSV/JSON, reports, fingerprints and archives remain unchanged.
They are read directly as frozen records; no regeneration requirement keeps
their writers alive. For verified historical conversion into shared comparison
evidence, use `research-repo-tools performance convert` with
`tooling/performance-legacy.toml`; this never reconstructs complete runs from
summaries. Do not promote converted historical summaries as new measurements.
The conversion declaration remains until every retained original has a verified
shared comparison companion and its original bytes remain independently available.

Generic parser, renderer, transaction and provider failure matrices belong
upstream. Consumer tests retain scientific selection/phase/gate checks, actual
native recipe wiring, complete release archive round trips, captured-target
handoff, offline publication and original artifact hashes.

## Changelog and release tooling

Use the shared commands exposed by `just --list`. The [release guide](../docs/RELEASING.md)
owns release metadata, changelog generation, tagging and registry publication.
Repository declarations in `pyproject.toml` retain the concept DOI, active link
policy, declared date policy, changelog formatter and dependency-body handling.

The publishing declaration selects only `la-stack`, requires the native platform
and scientific check results on the exact commit, and requires the durable
Criterion release asset. `release-gate`, `release-publish` and `release-verify`
use the shared release contract; `release-package` owns native default/exact
Cargo dry-runs. The local notes check binds the GitHub title/body to the generated
changelog. `scripts/tests/test_release_publishing.py` checks this consumer policy,
installed CLI failures and the actual workflow's upload barriers without network
access or registry credentials. Account registration and live release evidence
remain separate from fixture results.

The exact shared package pin lives in the `tooling` dependency group, included
by `dev`, and is resolved from PyPI in `uv.lock`. No sibling checkout supplies
runtime code. Shared `deps update-python` advances direct exact development
requirements; the shared package pin advances deliberately with its lockfile.

Setup, tool verification, managed execution, dependency updates, Markdown line
checks, Semgrep fixture validation, security scanning and opt-in CodeRabbit
review also use shared commands. Consumer tests check configuration and actual
recipe forwarding; shared implementation regressions belong upstream.

See the [adoption record](../docs/dev/performance-tooling-adoption.md) for deletion
dispositions, Python 3.15 review, size measurements and validation evidence.
