# Performance tooling adoption

[#296](https://github.com/acgetchell/la-stack/issues/296) completes the shared
performance-tooling migration begun in #276 and #285. The consumer retains
scientific declarations, native correctness gates, and installed-package
integration tests. Generic execution, validation, rendering, and publication
belong to research-repo-tools.

## Publication boundary

The dependency-only environment pins the published `research-repo-tools==0.2.4`
package. Release target capture and publication use the same exact version.
The suite verification, packaging and assembly commands adopted for
[#316](https://github.com/acgetchell/la-stack/issues/316), targeting v0.5.2,
use that release as well. Its verified tool updater is adopted without local
updater logic. A temporary scanner-selection override preserves the standard
Python 3.15 Semgrep environment until the shared fix ships in v0.2.5;
[#317](https://github.com/acgetchell/la-stack/issues/317) tracks its removal.
Upstream [#123](https://github.com/acgetchell/research-repo-tools/pull/123)
adds trusted main-validation artifact reuse for release promotion; v0.2.3
was published from commit `4085125526e21b2855eb7c61a52ee312782ae31a`.
It also refreshes explicitly declared independent Python projects while
preserving their selectors, exact pins, overrides, and separate environments.
The minimum-release acquisition contract introduced by v0.2.2 remains in place:
v0.2.1 could not express the consumer's v0.5.0 measurement floor.

[PyPI v0.2.4](https://pypi.org/project/research-repo-tools/0.2.4/) supplies both
distributions. Their SHA-256 hashes match the published
[GitHub release assets](https://github.com/acgetchell/research-repo-tools/releases/tag/v0.2.4)
and the exact registry lock entries:

| Distribution | SHA-256 |
| --- | --- |
| Wheel | `00aaab2b7efd680d79324c6debb64db4d8d1ff076a460ff46e24e0ef1c24bd7d` |
| Source distribution | `30b6d4d2f8d0cb98d5311656928dd65800088b651dd6c0c29f9a02526e9e902a` |

Runtime imports and acceptance checks use this installed package. No sibling
checkout or release-candidate installation supplies final evidence. The temporary
Semgrep adapter overrides only the installed scanner-selection function; all
benchmark execution, validation and publication use the public shared APIs.

The maintainer approved freezing earlier evidence instead of extending the
shared policy gate to accept missing acquisition declarations. This deliberately
supersedes #296's original v0.5.0/v0.4.6 replay requirement. The minimum release
applies to new acquisition; it does not redefine eligibility of already
declared complete runs.

## Entry-point replacements

The root Justfile calls the installed CLI directly. Static TOML describes
measurement and publication; the native Justfile contains Cargo commands.
There are no replacement Python wrappers, runtime configuration generators,
or shell implementations of the shared workflow.

| Consumer command | Installed owner and declaration |
| --- | --- |
| `performance-local` | `performance measure`, `measurements/full.toml`, current versus latest published release |
| `performance-local-non-exact` | `performance measure`, `measurements/non-exact.toml`, current versus latest published release |
| `performance-measure` | Explicit current/baseline labels and a named full/non-exact declaration |
| `performance-prepare` | `performance prepare` with `performance-release.toml`; independent saved assertions |
| `performance-release` | `performance release` with `performance-release.toml`; measure the prepared pair and promote |
| `performance-promote` | Explicit scratch payload and envelope, checked against `performance-report.toml` |
| `performance-doc` | Offline retained history through `performance promote` and the independent full policy |
| `performance-readme` | `performance readme` with `performance-release.toml`, prepared or explicit independent identities and transactional formatting |
| `plot-vs-linalg` | `performance publish --artifacts-only` with named full scratch, retained, or non-exact scratch views |
| `bench-release-check` | `performance verify-saved` with saved name, Criterion directory and `criterion-saved.toml` |
| `bench-release-check-suite` | `performance verify-saved --suite` with the complete independent policy and selected feature suite |
| Hosted suite packaging and assembly | `performance host`, `pack-suite` and `assemble-suites`; la-stack retains the native commands, matrix, artifact selection and budgets |
| `bench-compare` | Existing installed saved-Criterion comparison CLI |

Paths in the table are relative to `tooling/`.

### Release adapter migration

Published `research-repo-tools==0.2.4` now owns the former
`release_changelog.py`, `release_evidence.py` and `release_performance.py`
adapters. Citation dates are declarative in `pyproject.toml`; evidence, draft
assets and performance phases have dedicated TOML policies. The native Just
launchers observe checkout identity and forward arguments to the installed
commands. Evidence and asset verification use an isolated runtime outside the
checkout. Privileged benchmark jobs retain their existing installed asset
commands and captured target without checking out or executing repository code.

The upstream installed-consumer suites
[`public_evidence_consumer.py`](https://github.com/acgetchell/research-repo-tools/blob/v0.2.4/tests/releases/public_evidence_consumer.py)
and
[`public_dates_consumer.py`](https://github.com/acgetchell/research-repo-tools/blob/v0.2.4/tests/changelog/public_dates_consumer.py)
own generic output rollback, immutable handoffs, malformed metadata, signature
and date parsing regressions. Local tests retain the actual native/Cargo policy,
workflow approval barriers, real Just output transport, installed runtime
identity, stale preparation rejection and scientific publication fixtures.
The shared downloader's archive bound is 32 MiB; signed manifests remain
bounded at 1 MiB. Recovery uses the shared selection schema; legacy preparation
files must be replaced through a new independent preparation.
[Benchmarking](../BENCHMARKING.md) owns arguments and examples.

## Consumer contracts

- Both acquisition labels must be v0.5.0 or newer. The published parser rejects
  older labels before Git mutations or measurement output.
- Independent scientific catalogs require 417 full-report IDs. The non-exact
  declaration selects the comparison suite. Native discovery may capture
  additional diagnostic cases without weakening required coverage.
- Every phase retains 100 samples, mean and median estimates, 95% intervals,
  nanoseconds, and the exact discovered inventory. Both native scientific gates
  run before collection.
- nalgebra and faer references always originate in the baseline phase.
  Compatibility includes host, CPU, architecture, harness, Rust/Cargo/Just, and
  dependency identities.
- Measurements use the captured current harness with explicit Rust toolchain,
  native host target, compiler/wrapper settings, flags, and isolated output
  directories. Ambient code-generation overrides no longer silently change
  the measured configuration.
- Snapshots include tracked changes and nonignored untracked regular files.
  This intentionally replaces tracked-only snapshots. Source and harness
  membership and bytes are checked again before README publication.
- Saved samples require the independent inventory, canonical storage paths,
  matching fresh/saved raw files and complete statistical data. The configured
  raw inventory remains `benchmark.json`, `estimates.json`, `sample.json`, and
  `tukey.json`.
- Scratch and retained inputs are explicit choices. Incomplete scratch never
  falls back to an older run. Promotion and document publication validate
  independent declarations, not only a payload's self-consistency.
- All eight scientific plot metrics retain mean/median and linear/log choices,
  giving 32 named views for each input source. Non-exact scratch views use
  separate asset paths; official release promotion requires the full policy.
  The README has one canonical
  LU-solve median/log view and one shared publication section.
- README publication requires independently supplied releases and full Git
  revisions, a retained run, current source/harness freshness, package/report
  reference checks, tag checks, and successful formatting before any output
  changes. The shared publisher emits full immutable links, so README uses
  rumdl's Markdown-aware line policy; the extra raw-character gate continues
  to cover other active documents.

The minimum release, scientific rows, native commands, and publication choices
remain reviewable static consumer data. The shared engine owns their runtime
validation. `tooling/performance-policy.toml` is the independent scientific
catalog used by acceptance tests, not a local production parser.

## Historical evidence

Historical reports, plots, CSV/JSON, archives, and complete-run files retain
their original bytes and links. The existing v0.5.0/v0.4.6 complete run
`634bc9527c2caf0a431a9bdfac082c190e1094f0b29eb6e3f28c4b96961cdec5`
is a frozen readable record. It predates acquisition declarations and is no
longer an input for regeneration or current publication. Do not backfill its
provenance or relax the shared eligibility gate.

New declared complete runs support offline report and plot rendering after
scratch deletion. README publication additionally requires current-source
freshness. New axis-specific plot paths preserve existing historical assets;
release metadata updates exclude immutable archive contents and their README
links. Regression hashes protect both the legacy CSV/JSON evidence and the
frozen complete run.

## Deletion and test ownership

The retirement removes `archive_performance.py`, `performance_runs.py`,
their adapter tests, and their obsolete measurement/plot declarations.
Useful non-performance integration tests move from `scripts/tests/` to
`tests/tooling/`. The root `scripts/` directory then has no production or
test owner. Semgrep fixtures under `tests/semgrep/scripts/` and historical
security-exception paths keep their existing identities.

| Generic behavior | Upstream owner | Consumer evidence retained |
| --- | --- | --- |
| Floors, captured phases, complete samples and independent policy validation | `tests/performance/public_complete_consumer.py` | Actual full/non-exact declarations, release floor, scientific coverage and native commands |
| Transactional publication, formatting and provenance | `tests/publication/public_publication_consumer.py` and complete-run consumer suite | Actual recipes, stale source/harness, partial scratch, formatter failure, immutable links and frozen hashes |
| Saved-sample validation and archive transport | `tests/performance/public_performance_consumer.py` and Criterion baseline tests | Actual workflow archive, required IDs, canonical storage and saved/fresh bytes |
| Changelog, release target capture, toolchains and maintenance | Capability-specific installed-consumer suites | Existing la-stack workflow, package identity, release policy and configuration checks |

The shared package runs public consumer suites from isolated wheel and sdist
installations. la-stack acceptance imports public APIs only and invokes the
actual Just recipes against bounded synthetic samples. No historic benchmark
reruns or new performance claims are needed.

## Runtime and validation

The dependency-only environment remains Python 3.15. Python supports shared
tools and consumer tests; deleting production scripts does not remove it.
Semgrep initially used a separately locked Python 3.14 environment. Following
[MCMC #214](https://github.com/acgetchell/markov-chain-monte-carlo/pull/214),
Semgrep 1.182.1 scans and fixture checks share the root Python 3.15 environment and lockfile.
The root manifest retains its native-wheel requirement, compatibility overrides
and Python 3.15 Windows wheels. Local and hosted scans retain the Linux
stack-limit workaround. Consumer tests now live under `tests/tooling/`.

Initial migration discovery listed 533 comparison and 354 exact cases while verifying
the independent 417 required IDs. The scientific gates passed 26 exact-input
tests and 42 comparison-input tests. These checks compiled/listed the native
benchmarks without collecting timing measurements.
The expanded v0.5.2 selection and its current correctness gates are documented
in [Benchmarking](../BENCHMARKING.md#release-acquisition-policy).

Final acceptance also requires focused installed-package tests, `just check`,
`just ci`, Markdown/spelling/docs checks, and Linux/macOS/Windows CI for the
adoption commit. Local macOS validation does not establish other native results.
Release and adoption PRs retain the source-bound validation evidence.
