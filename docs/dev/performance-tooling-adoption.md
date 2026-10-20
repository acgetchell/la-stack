# Performance tooling adoption

This records [#276](https://github.com/acgetchell/la-stack/issues/276) against
the published `research-repo-tools==0.1.9` base package. No extras are required.
Publication prerequisite [upstream #99](https://github.com/acgetchell/research-repo-tools/issues/99)
is closed. PyPI's unyanked wheel and sdist hashes match the exact lockfile entries;
runtime code comes from the installed package.

## Capability disposition

Follow-up [#285](https://github.com/acgetchell/la-stack/issues/285) adopts the
published `research-repo-tools==0.2.0` wheel and sdist, including
[upstream #110](https://github.com/acgetchell/research-repo-tools/issues/110).
The exact PyPI pin and lock replace v0.1.9; the package also advances its required
Just dependency to v1.59.0. `RunTableStyle` now supplies the six-column README
layout directly, so `_readme_table()` and its Markdown string/column rewrites
are removed. The shared renderer retains recorded intervals, the selected
statistic, confidence level and explicit ratio direction. Scientific selection,
original reference phases, README markers, evidence freshness and transactional
asset pointers remain consumer responsibilities.

Offline consumer checks exercise the actual `performance-readme` recipe with
mean and median estimates, original peer phases, retained history and the README
figure. No new measurements or regenerated immutable reports are required.
Offline rendering of the retained published run preserves its SVG and CSV bytes.
Current README republication still rejects the stale source digest and leaves
the figure and measured-asset pointers untouched.
The obsolete CSV conversion declaration and pre-Just Python phase replay are
removed. Retained complete runs use the declared Just driver; original legacy
CSV/JSON and Markdown records remain available without regeneration or conversion.
Closure of #285 requires native Linux/macOS/Windows CI on the adoption commit;
local validation alone does not establish those results.

| Upstream capability | Consumer adoption | Local responsibility |
| --- | --- | --- |
| #94 captured targets | Installed release-target CLI passes versioned JSON from trusted preflight to writer | Dispatch/ref equality, current attempt, title, privilege separation, asset filename and explicit publication |
| #95 coordinate plots | Shared coordinate SVG, table and CSV renderers | Scientific metric mappings, dimensions, nanoseconds and original phases |
| #96 provenance reports | Shared renderer 2 and publication-relative evidence links | Scientific eligibility and frozen historical records |
| #97 complete-run sections | Shared presentation and publication transaction | README markers, measured-asset pointers and current source/harness checks |
| #98 Criterion verification | Declarative shared discovery, versioned inventory and saved-sample APIs | Full Cargo suites/features, minimum scientific coverage and correctness gates |
| #110 compact table styles | Installed `RunTableStyle` selects headings, intervals, confidence caption and ratio-only columns | Scientific selections, ratio pairs and README publication policy |

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
| `performance_artifacts.py` | Deleted; shared evidence/publication replace writers and readers; historical CSV conversion retired in #285 |
| `release_baseline.py` | Deleted; native discovery uses the shared CLI and saved-copy verification uses the shared configured API |
| Local wheel/sdist and console-entry-point packaging | Removed; `uv` manages a dependency-only environment |
| Duplicate parser, renderer, legacy-writer and packaging tests | Removed; actual installed-package/configuration and scientific integration tests remain |

Retirement conditions are satisfied by installed replacement checks, preservation
of original bytes/hashes and a documented supported workflow. The old report
layout need not be regenerated indefinitely. New reports use the shared format;
historical reports, CSV/JSON, indexes, plots and archives remain untouched.
Historical conversion was verified during #276 and retired in #285. Original
summaries remain independently readable with their original hashes and opaque
metadata. They cannot provide complete-run samples or newly captured host provenance.

Only two production modules remain:

- `archive_performance.py`: la-stack release eligibility, tracked-only snapshots,
  and README publication composition.
- `performance_runs.py`: captured measurement environment, scientific selection,
  sampling, originating reference phases and source/harness eligibility.

Before v0.5.1, new measurements are restricted to v0.5.0 or newer on both sides.
The historical API cfg, missing-row exemptions, root-test layout adapter and
metadata-free `performance-github-assets` command are retired. Their Rust
benchmark branches are removed while the independent scientific gates remain.
Original reports, CSV/JSON, plots and archives retain their bytes and hashes.
The measurement floor does not block offline use of retained complete runs
satisfying the modern scientific policy.

Shared CLI snapshots include untracked files; la-stack's established contract
excludes them. Independent scientific eligibility, saved-sample storage checks
and publication composition also remain local.
[Upstream #115](https://github.com/acgetchell/research-repo-tools/issues/115)
targets declarative replacements in v0.2.1. Adoption and deletion of the final
production Python modules follow publication of that package.

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

The #285 consumer-test cleanup was checked against tests bundled in the published
v0.2.0 sdist. Shared implementation matrices are removed only where that release
already covers them; consumer policy and installed-package integration remain.

| Coverage | Published upstream owner | Retained consumer check |
| --- | --- | --- |
| Changelog rendering and repeatability | `tests/changelog/test_cliff_template.py` and `test_workflow.py` | Actual la-stack formatter/preservation policy, one authored-note generation and historical notes |
| Measurement completeness and failed gates | `tests/performance/public_complete_consumer.py` | Actual sampling policy, required scientific inventory, native gate commands and original peer phases |
| Package identity | Consumer-owned | One test in `test_toolchain_integration.py` checks pin, registry, installed version and wheel/sdist hashes; recipe tests verify preservation and forwarding |
| Release event identity and check-status matrices | `tests/releases/test_publishing.py` | Every configured native check remains required; actual workflow approval, immutability, notes and upload barriers remain tested |

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
