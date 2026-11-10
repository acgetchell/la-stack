# Release and validation integration

Historical implementation record for [#326](https://github.com/acgetchell/la-stack/issues/326),
[#330](https://github.com/acgetchell/la-stack/issues/330)
and [#338](https://github.com/acgetchell/la-stack/issues/338). The starting source
is `7c3a3625220773369c089e64e568a979abc0bb7f`. The authorized supplying revision is
research-repo-tools `cd8594f658b5bf3ade063a685834565cc62e0cb3` (PRs #207 and #209).
Its installed metadata remains 0.2.6; the Git commit, manifest, lock and isolated
runtime assertions distinguish this development rehearsal from publication.
The [hosted CI stages](#hosted-ci-stages) section records the subsequent published
v0.2.8 integration and workflow split in PR #343; earlier timings retain their
original source and supplier identities.

## Workstream readiness

[PR #335](https://github.com/acgetchell/la-stack/pull/335) merged at the starting
source above. Its Linux, macOS, Windows, Cargo packaging and static checks passed;
all six prerequisites of #326 (#320 through #325) are closed. This integration
completes the remaining #330 command/guide and automatic-draft work, plus #338's
validation graph and cache wiring. Required PR checks and review still gate merge.
Final release validation and published-tooling adoption belong to #333.

## Release contract

The [release command comparison](../RELEASING.md#release-command-alignment)
records the v0.2.6 reference revision, common operator commands and Cargo/benchmark
differences. The eligible final-main attempt owns native tests and default/exact
packaging. Step 2 inspects retained evidence before tagging. Tag-triggered draft
preparation uses the supplier's signed-selection and durable intent API; the
read-only job passes inert notes and evidence to an isolated writer through an
authenticated immutable artifact. Initial staging requires the authenticated
preparation target and never discovers a replacement draft. The supplier's
tag-aware receipt retrieval freezes the original provider run/attempt and artifact
before using its target. Step 1 registers one immutable external owner selection
per checkout and tag. Routine commands recover that owner and derive clean HEAD;
explicit owner/target/transfer inputs remain manual recovery options.

Generic draft lifecycle and retry logic remain upstream. Existing owner records,
scientific workloads, raw measurements, generated reports, archive identity,
publication approval and nonbuilding Cargo upload remain separate responsibilities.
The development-pin gate still prevents acquisition for release, signed evidence
creation and publication until a published supplying package is adopted.

## Validation graph

Previously `ci` listed all validators sequentially. The new graph runs shared
`python-sync` and `tools-check`, then `_ci-policy`, then the two dependencies of
`[parallel] _ci-execution`. One Just invocation coalesces shared prerequisites;
there are no nested Just processes for lane orchestration. `just --jobs 1 ci`
serializes the same graph.
This uses [Just's parallel dependency scheduling](https://just.systems/man/en/parallelism.html).

| Owner | Existing validators retained |
| --- | --- |
| Shared setup | `python-sync`, `tools-check`, their prerequisite checks |
| `_ci-policy` | `action-lint`, `workflow-allowlist-check`, `zizmor`, Markdown/spelling, private nonfinal version check, changelog, TOML, YAML, citation, JSON, Just formatting, Python static checks, Cargo lock/formatting, Semgrep/fixtures, unused dependencies, shell checks |
| `_ci-rust` | `clippy-all-targets`, `doc-check`, `test-rust-ci`, `test-doc`, `test-doc-exact`, `bench-compile`, `test-bench-exact`, `examples` |
| Parallel Python lane | `test-python`, with explicit `--no-sync` after shared setup |
| Separate main-CI packaging job | `publish-check` (formerly recipe `release-package`); hosted job and step names remain unchanged |
| Explicit acquisition | No timing recipe joins validation; benchmark collection remains opt-in and isolated |

All root Cargo output users stay sequential. Python consumer fixture workspaces
remain isolated. No source, installed-distribution or fresh-setup matrix is
invented for this private dependency-only Python environment. There was no
consumer fresh-install smoke job to remove. Existing native Linux/macOS/Windows
`build` jobs, `Run CI checks`, `release-package`, its packaging step, signer
dependencies and signer steps remain unchanged, so evidence policy and required
check names require no migration.

`tests/tooling/test_validation_graph.py` executes the real graph with inert
validator bodies. It requires identical complete coverage and one setup in
parallel and serial modes, proves lane overlap, and injects failures before and
inside each lane. The native scientific tests and package checks still run
through the owning recipes; graph fixtures do not replace them.

## Managed compilation cache

The supplying package owns the local Cargo build-cache implementation and exact
tool installation verification. Ordinary CI's existing setup action now caches
only `RUNNER_TEMP/rrt-cargo-build`, supplied through the documented absolute
`RESEARCH_REPO_TOOLS_CARGO_BUILD_HOME` override. Exact-version installation homes
under `RUNNER_TEMP/rrt` are separate and freshly verified. Cache keys include OS,
architecture, Rust toolchain, Python/tool declarations and both locks, without
broad restore prefixes. GitHub's branch/PR cache scopes remain in force; PR
entries cannot replace protected main entries. No remote cache crosses those
provider boundaries. This is compilation reuse, not validation reuse.

Release packaging, preparation, evidence signing, benchmark attachment and Cargo
publication keep caches disabled. No compiler flags, scientific inputs or
benchmark directories change.

## Hosted CI stages

PR #343 separates hosted validation while retaining the complete local `just ci`
graph. Each preparation job installs the declared tools through RRT and completes
the same policy validators. Its successful completion saves managed compilation
outputs. After all three preparation jobs succeed, native jobs start independently
and install verified tools using that compilation cache. Installation homes,
project outputs and validation results are never shared between these jobs.
A cache miss may rebuild tools; it cannot skip checks or satisfy provenance.

| Hosted owner, on each native OS | Canonical command and complete responsibility |
| --- | --- |
| `prepare (OS)` | RRT setup, then `just ci-policy`: every `_ci-policy` validator and shared setup check listed above |
| `native-rust (OS)` | `just ci-rust`: the complete sequential `_ci-rust` lane, including compilation, tests, docs, examples and benchmark validation |
| `native-python (OS)` | `just ci-python`: every Python test without the `release_candidate` marker, using four workers |
| `native-candidate (OS)` | `just ci-candidate`: the complete candidate capture/acceptance regression, using one worker on a dedicated runner |
| `native-dependencies (OS)` | `just ci-dependencies`: locked Cargo fetch, native build receipts and dependency-provenance verification |
| `build (OS)` | Required-check name preserved on a lightweight Linux gate; rejects any failed, skipped, cancelled or absent preparation/native result |
| `release-package` | Independent uncached dependency verification and default/exact Cargo package validation |
| `release-evidence` | Protected main-push signer; requires all aggregate gates and packaging, then verifies every actual native job, required step and OS label |

`tests/tooling/test_validation_graph.py` collects the real pytest suite under both
hosted selections and requires a disjoint, complete partition. Local `just ci`
still selects the full Python suite once. Gate regressions execute the workflow's
actual shell for success, failure, skipped, cancelled and missing results.
Release-evidence regressions reject incomplete steps in every native stage.
Scientific fixtures, assertions and benchmark sampling policy are unchanged.

The preceding monolithic run at `1548f59`,
[42858660611](https://github.com/acgetchell/la-stack/actions/runs/42858660611),
passed Linux, macOS and uncached packaging. Windows spent about 51 minutes in
setup and 62 minutes in `just ci`, including 465 passing Python tests, before
its two-hour limit interrupted Cargo provenance verification. Python tests took
3,632.39 seconds; candidate acceptance alone took 1,783.31 seconds plus its fixture
setup. These observations motivate resource separation; they do not establish a
speedup. Native timing and cache-hit evidence for the split workflow must be
reconciled on its final source revision.

## Timing evidence and limits

Baseline native data comes from successful main-push
[CI run 42443172694, attempt 1](https://github.com/acgetchell/la-stack/actions/runs/42443172694/attempts/1)
on the starting source. Durations below use API job/step timestamps and are
descriptive wall time, not benchmark results or a cold/warm comparison.

| Job | Queue from run start | Tool setup | CI or packaging step | Total job wall time |
| --- | --- | --- | --- | --- |
| Linux build | 3 s | 23 m 51 s | 26 m 27 s | 50 m 40 s |
| macOS build | 4 s | 41 m 59 s | 36 m 27 s | 80 m 35 s |
| Windows build | 3 s | 1 m 49 s | 62 m 18 s | 66 m 57 s |
| Release packaging | 3 s | 21 m 23 s | 4 s | 21 m 35 s |
| Evidence job | Waited for native jobs | — | Development-pin signing/upload skipped | 14 s |

The complete workflow occupied 80 m 56 s from creation to completion. The four
build/package job intervals total 219 m 47 s of runner wall time; this is distinct
from critical-path duration and excludes provider billing/rounding. Setup cache
hit state was not established, so none of these observations is labelled cold
or warm. The evidence job's success does not establish signed release evidence:
its development-pin signing and artifact upload steps were skipped.

Locally, the original `just check` graph at the new exact supplier pin completed
in 52.67 s elapsed (27.52 s user, 13.11 s system) on the existing macOS environment
before recipe restructuring. This is a warm development environment, not a clean
managed-install experiment. It is not directly comparable to full hosted CI.

Final local validation at the supplying revision above passed on macOS 27.2,
arm64, Python 3.15 and Rust 1.99.1, using the existing development environment
and four Python test workers:

| Command | Elapsed | User CPU | System CPU | Python suite |
| --- | --- | --- | --- | --- |
| `just --jobs 1 ci` | 905.35 s | 2,152.17 s | 365.49 s | 438 passed in 676.01 s |
| `just ci` | 736.74 s | 1,640.82 s | 340.51 s | 438 passed in 691.83 s |

Both runs completed the same validators, including 1,143 Rust tests, default/exact
doctests, benchmark compilation, the exact harness in test mode and examples.
The serial run compiled missing targets; the later parallel run reused them.
These are descriptive observations with different compilation cache states,
not an isolated measurement of parallelism or managed-cache speedup. The long
synthetic acquisition/acceptance consumer is a remaining bottleneck: its call
alone took 525.65 s in the serial run. No scientific measurements were acquired.

The focused owner/receipt/dispatch and graph regressions passed, including real
Just rejection of publication without prior review. Windows-targeted static
type checking and the exact lock check passed. Astra's source review found no
remaining trust-boundary blocker; it did not substitute for native validation.

The implementation PR closes #326, #330 and #338 on merge.
[#333](https://github.com/acgetchell/la-stack/issues/333) owns final native
validation, authenticated hosted retained-evidence handoff, matched cold/warm
managed installation timings and published-tooling adoption. No speedup is
claimed from configuration, graph-fixture timing or these unmatched cache
observations. No release, tag or benchmark campaign is run as part of this
integration.
