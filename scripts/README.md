# scripts/

Only two production Python adapters remain. Generic tools run from the installed
`research-repo-tools==0.2.0` base package; their declarations belong in
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
| `archive_performance.py` | Release eligibility, tracked-only snapshots, and la-stack README markers/pointers |
| `performance_runs.py` | Measurement environment, scientific coverage, sampling and original reference-phase assertions |
| `tooling/performance-policy.toml` | Exact group selections, metric names, peers, dimensions and minimum measurement release |
| `tooling/performance.toml` | Shared sources, harness, probes, dependencies, timeouts and compatibility |
| `tooling/performance.just` | Native Cargo discovery, independent correctness gates and timing commands |
| `tooling/performance-report.toml` | Shared current report and immutable complete-run history |
| `tooling/performance-plot.toml` | Scientific coordinate, unit, series and axis declarations |
| `tooling/criterion-inventory.toml` | Shared native discovery and minimum coverage of all report/plot consumers |
| `tooling/criterion-saved.toml` | Shared complete saved-sample policy, including exact `tukey.json` bytes |

The adapters call public shared APIs for subprocesses, worktrees, Criterion
discovery, completeness, estimates, comparison arithmetic, evidence, rendering,
retention and publication. They
contain no generic renderer, CSV parser, legacy writer or measurement engine.

New measurements require v0.5.0 or newer on both sides. Historical API cfg,
missing-row exemptions, root-test layout adapters and the metadata-free
`performance-github-assets` command are retired. Archived artifacts remain
unchanged; retained complete runs satisfying the modern policy still render
offline regardless of the measurement floor.

The remaining composition preserves tracked-only local snapshots, independent
scientific eligibility checks and baseline-origin nalgebra/faer reference views.
The shared CLI's default worktree snapshot includes untracked files.
[Upstream #115](https://github.com/acgetchell/research-repo-tools/issues/115)
tracks the declarative interfaces needed to retire these final adapters.

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
README publication configures the shared `RunTableStyle` API to label dimensions,
keep both interval bounds below each timing, state the recorded confidence level
once, and select ratio-only speedup columns with explicit numerator/denominator
semantics. Series labels supply the timing headings and measured versions.
The shared renderer owns formatting, timing estimates and comparison arithmetic.

```bash
just performance-release
just performance-readme lu_solve median new true
just plot-vs-linalg dot mean new true
```

The raw-Criterion/gnuplot publisher, historical summary writer and custom
report renderer are retired after installed-package replacement checks.
Original CSV/JSON, reports, fingerprints and archives remain unchanged.
They are read directly as frozen records; no regeneration requirement keeps
their writers alive. The obsolete CSV conversion declaration and pre-Just Python
phase format are retired. Report and README publication require complete-run
evidence with the declared Just driver. Historical summaries lack raw samples
and cannot supply that evidence; their original bytes remain available directly.

Generic parser, renderer, transaction and provider failure matrices belong
upstream. Consumer tests retain scientific selection/phase/gate checks, actual
native recipe wiring, complete release archive round trips, captured-target
handoff, offline publication and original artifact hashes.
`scripts/tests/test_toolchain_integration.py` owns the single published-package
contract: manifest pin, registry lock, installed version and wheel/sdist hashes.
Dependency-update and release-workflow tests derive that requirement from the
manifest and verify it survives their actual recipe execution.

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
