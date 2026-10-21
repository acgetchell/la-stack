# Tooling

The installed `research-repo-tools` package owns common development and
performance tooling. The exact released pin and registry hashes live in
`pyproject.toml` and `uv.lock`. This dependency-only environment has no local
Python distribution, console entry points, or production Python adapters.

## Setup and validation

Use `just setup`, then `uv run --locked` for Python commands. Python 3.15
runs the shared package and the consumer tests in `tooling/tests/`.
`just python-check` checks formatting, lint and types; `just test-python`
runs the installed-package integration tests. `just ci` is the final
repository validation. Native Linux, macOS and Windows CI is required for
adoption; local macOS results do not establish the other platforms.

Semgrep 1.180.0 retains its separate locked Python 3.14 environment in
`tooling/semgrep/` while its wheels exclude Python 3.15. Ordinary Rust and
Python tests remain scan targets; deliberate `tests/semgrep/` fixtures are
checked separately by `just semgrep-test`.

## Performance ownership

The [benchmark guide](../docs/BENCHMARKING.md) owns command selection,
interpretation and output paths. Every performance recipe calls the installed
CLI or a native Cargo command. There is no runtime configuration generator.

| Declaration | Responsibility |
| --- | --- |
| `measurements/full.toml`, `measurements/non-exact.toml` | Release floor, required cases, sampling, original peer phases, source/harness inventories, probes and environments |
| `performance-policy.toml` | Independent scientific catalog used by consumer acceptance tests |
| `performance.just` | Native discovery, correctness gates and measurement commands |
| `criterion-inventory.toml` | Release discovery and independent required scientific IDs |
| `criterion-saved.toml` | Required IDs, canonical paths, 100 samples and exact saved raw JSON bytes |
| `performance-report.toml` | Independent release eligibility, current report and immutable run history |
| `performance-readme.toml` | Canonical README selection, formatter, source freshness and tagged links |
| `plots/{non-exact,retained,scratch}/*.toml` | Named metric/statistic/axis selections with explicit input ownership |

Both labels of a new measurement must be v0.5.0 or newer. The full and
non-exact declarations require 100 samples, mean and median estimates,
95% intervals and nanoseconds. Both native correctness gates run before
timing. Reference libraries originate in the baseline phase.

Working-tree measurements include tracked changes and nonignored untracked
regular files. Scratch outputs belong under ignored `target/`.
Tag preparation is explicit: `just performance-prepare` fetches tags and
requires authorization to mutate Git. Measurement recipes create temporary
worktrees. Their declared environment fixes the Rust toolchain and
code-generation flags; native recipes select the host target. Change and
review the declarations for a different measurement environment instead of
relying on inherited Cargo flag overrides.

Input selection is explicit. `performance-promote` consumes the full scratch
pair; `performance-doc` selects retained history. A missing scratch file is
an error. The plot recipe selects a full scratch, retained, or non-exact
scratch configuration. Non-exact plots use separate asset paths; official
release promotion and README publication require the full measurement policy.
README publication requires archived release-eligible evidence, independently
supplied release/commit assertions, matching current source/harness membership
and digests, and a successful formatter before replacing any output.

Historical CSV/JSON, plots, reports and complete runs remain frozen with their
original bytes and links. This includes the v0.5.0/v0.4.6 comparison: it lacks
the newer acquisition declaration and is no longer a regeneration input.
Do not invent missing provenance or relax publication checks to replay it.
Future complete runs acquired with the declarations above support offline
report and plot rendering after scratch cleanup.

Shared implementation matrices live upstream, including
`tests/performance/public_complete_consumer.py` and
`public_criterion_consumer.py`, exercised from wheel and sdist installations.
Consumer tests retain independent scientific selections, actual recipe wiring,
saved archive round trips, publication boundaries and historical hashes.

## Changelog and release tooling

Use the shared commands exposed by `just --list`. The
[release guide](../docs/RELEASING.md) owns metadata, changelog generation,
tagging and registry publication. Declarations in `pyproject.toml` retain
the concept DOI, active link policy, declared date policy, formatter and
authored dependency notes.

The publishing declaration selects only `la-stack`, requires native platform
and scientific checks on the exact commit, and requires the durable Criterion
release asset. `release-gate`, `release-publish` and `release-verify`
use the shared release contract; `release-package` owns native default/exact
Cargo dry-runs. `tests/test_release_publishing.py` checks consumer policy and
upload barriers without registry credentials.

`tests/test_toolchain_integration.py` owns the package-identity check:
manifest pin, registry lock, installed version and distribution hashes.
Dependency updates and release tests derive their requirement from the
manifest. The shared package pin advances deliberately with its lockfile;
no sibling checkout supplies runtime code.

Setup, managed execution, dependency updates, Markdown checks, Semgrep,
security scans and opt-in CodeRabbit review also use shared commands.
See the [adoption record](../docs/dev/performance-tooling-adoption.md) for
deletion dispositions and release validation evidence.
