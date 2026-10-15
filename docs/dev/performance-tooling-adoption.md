# Performance tooling adoption

This records [#276](https://github.com/acgetchell/la-stack/issues/276) against
the published `research-repo-tools==0.1.9` base package. No extras are required.
Publication prerequisite [upstream #99](https://github.com/acgetchell/research-repo-tools/issues/99)
is closed. PyPI's unyanked wheel and sdist hashes match the exact lockfile entries;
runtime code comes from the installed package.

## Capability disposition

| Upstream capability | Consumer adoption | Local responsibility |
| --- | --- | --- |
| #94 captured targets | Installed release-target CLI passes versioned JSON from trusted preflight to writer | Dispatch/ref equality, current attempt, title, privilege separation, asset filename and explicit publication |
| #95 coordinate plots | Shared coordinate SVG, table and CSV renderers | Scientific metric mappings, dimensions, nanoseconds and original phases |
| #96 provenance reports | Shared renderer 2 and publication-relative evidence links | Scientific eligibility and frozen historical records |
| #97 complete-run sections | Shared presentation and publication transaction | README markers, measured-asset pointers and current source/harness checks |
| #98 Criterion verification | Declarative shared discovery, versioned inventory and saved-sample APIs | Full Cargo suites/features, minimum scientific coverage and correctness gates |

The full-report standardization from #277 remains. New reports expose recorded
provenance and evidence links. Shared mixed-renderer history preserves immutable
older complete-run reports while upgrading mutable indexes.

## Python deletion audit

The audit compared the dependency-only MCMC environment and declarative
`tooling/` commands with the installed v0.1.9 public APIs.

| Former helper | Disposition and supported replacement |
| --- | --- |
| `bench_compare.py` | Deleted; saved-baseline reports use the installed shared Criterion comparison CLI |
| `benchmark_contract.py` | Deleted; metric/group/peer/API declarations moved to `tooling/performance-policy.toml`; historical digests stay frozen |
| `benchmark_summaries.py` | Deleted; frozen historical summaries are records rather than regeneration inputs |
| `criterion_dim_plot.py` | Deleted; shared complete-run SVG/CSV/table presentation replaces raw-Criterion/gnuplot rendering |
| `performance_artifacts.py` | Deleted; shared evidence/publication and declarative `performance-legacy.toml` replace CSV/provenance writers and readers |
| `release_baseline.py` | Deleted; native discovery uses the shared CLI and saved-copy verification uses the shared configured API |
| Local wheel/sdist and console-entry-point packaging | Removed; `uv` manages a dependency-only environment |
| Duplicate parser, renderer, legacy-writer and packaging tests | Removed; actual installed-package/configuration and scientific integration tests remain |

Retirement conditions are satisfied by installed replacement checks, preservation
of original bytes/hashes and a documented supported workflow. The old report
layout need not be regenerated indefinitely. New reports use the shared format;
historical reports, CSV/JSON, indexes, plots and archives remain untouched.
Shared historical conversion verifies the original hash and retains opaque
metadata and original digest names. It produces comparison evidence, never
fabricated complete-run samples or newly captured host provenance. Retire its
conversion declaration only after verified companions exist for every original
and those original bytes remain independently available.

Only two production modules remain:

- `archive_performance.py`: la-stack release eligibility, tracked-only snapshots,
  historical raw archives without metadata, and README publication composition.
- `performance_runs.py`: captured old-API cfg/environment, scientific selection,
  sampling, originating reference phases and source/harness eligibility.

These are consumer policies rather than shared implementations. v0.4.4/v0.4.5
need the historical cfg and root-test layout adapter. Shared CLI snapshots include
untracked files; la-stack's established contract excludes them. The shared plan
cannot select cfg and API availability from a release pair. Raw historical release
archives lack the metadata required by the shared asset-evidence CLI. Their
small adapter uses shared parsing/rendering and explicitly identifies unrecorded
source, host and gates without constructing an evidence envelope.

Generic subprocesses, discovery parsing, worktrees, measurement, completeness,
estimates, comparison arithmetic, rendering, run identity, history, downloads,
extraction and publication transactions all remain upstream. The retained Python
calls public APIs. Scientific names and mappings live in `tooling/`.

## Python 3.15 review

Metadata, `.python-version`, the locked environment and CI select Python 3.15.
The remaining code uses deferred annotations, typed boundaries and explicit
byte/UTF-8 transport. Scientific catalogs parse into immutable typed records before
discovery; malformed dimensions, names, peer pairs and metric rows retain field
diagnostics. Existing evidence records are supplied by the shared
package. Deleting the historical serializer also removes its custom recursive
freezing/thawing machinery. No artificial syntax modernization replaces a
supported shared implementation.

Semgrep 1.180.0 wheels stop at Python 3.14. `tooling/semgrep/` pins the current
scanner and Python 3.14. The obsolete PyJWT override is retired because the pinned
scanner accepts patched PyJWT 2.15.1 with its crypto extra. Local scans,
fixture validation and hosted SARIF use that lock. Retire this exception after
current published wheels support Python 3.15 and native gates pass.

## Size and validation evidence

Counts use physical lines from `wc -l scripts/*.py scripts/tests/*.py` and
`uv run --locked pytest --collect-only -q`. Baseline is `905aba7`; test counts
describe consumer checks, not measurements.

| Measure | Before | After |
| --- | ---: | ---: |
| Production modules | 8 | 2 |
| Production Python lines | 5,323 | 788 |
| Test Python lines | 5,096 | 1,952 |
| Collected tests | 371 | 130 |

Production Python is reduced by about 85%. Focused installed-package tests cover
complete-run gates/coverage, original reference phases, offline report/plot replay,
source membership and unindexed-history freshness, captured release targets,
native archive round trips and unchanged historical hashes. The actual saved
comparison recipe and declarative historical conversion passed replacement checks.
Native discovery compiled/listed 533 comparison and 264 exact cases with all
417 required report/plot IDs, without collecting timings.

Local executor: macOS arm64, CPython 3.15.0; Semgrep uses CPython 3.14.7.
The final repository validation includes Python format/lint/types, consumer tests,
Rust tests/doctests, docs, examples, benchmark compilation and native static checks.
No new scientific measurements or regenerated historical artifacts are claimed.

Issue closure still requires an adoption commit, native Linux/macOS/Windows CI
for that exact commit, linked run evidence and merge before v0.5.0. Local macOS
validation does not substitute for native Linux or Windows results.
