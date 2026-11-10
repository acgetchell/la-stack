# Tooling

The installed `research-repo-tools` package owns common development and
performance tooling. The exact supplier pin lives in `pyproject.toml` and
`uv.lock`. This dependency-only environment has no local
Python distribution or console entry points. Declarative release policies select
the native Cargo gates; shared tooling is not copied into this repository.

## Setup and validation

Use `just setup`, then `uv run --locked` for Python commands. Python 3.15
runs the shared package, Semgrep and the consumer tests in `tests/tooling/`.
`just python-check` delegates formatting, lint and strict types to the shared
complete tracked/nonignored Python inventory, including static-analysis fixtures.
The `research-repo-tools[python-tools]` extra owns the pytest, Ruff and Ty
versions. `toolchain.inherit-python-tools` verifies declarations, the locked
development graph and executable versions. Local lint rules remain in this
repository; the `dev` group does not carry competing tool pins.
`just test-python`
runs every installed-package integration test with four isolated pytest workers,
work-stealing scheduling and no worker restarts. Set `LA_STACK_TEST_WORKERS=0`
for serial diagnosis or another fixed count to bound local resource use. The
slowest test durations are reported; assertions and test selection are unchanged.
`just ci` is the final
repository validation. Native Linux, macOS and Windows CI is required for
adoption; local macOS results do not establish the other platforms.

After shared setup and fast policy checks, `just ci` overlaps Python consumer
tests with a single sequential Cargo lane. `just --jobs 1 ci` preserves serial
execution. No measurement recipes join either lane. The shared managed Cargo
build cache retains compilation separately from exact-version installations;
ordinary CI caches only that compilation directory, with exact OS/architecture
and tool/lock keys. Release jobs disable caches. See
[validation integration](../docs/dev/release-flow-integration.md).

For subprocess fixtures, run Bash and Python source from files, and pass paths
and payloads as arguments or environment data. Nested `-c` commands add Windows
argument parsing boundaries that POSIX shell quoting does not protect. Supply an
explicit child environment when testing empty values: Windows `putenv` removes
them from the inherited native environment. Exercise Windows path data and
quoting characters as well as native Windows execution. Negative tests must
prove the received arguments, intended rejection and unchanged output; a nonzero
exit alone also accepts a broken import or shell. When quiet mode hides
diagnostics, retain the verifier's arguments, exit status and diagnostic in a
fixture-owned result file.

The candidate and native evidence recipes use Just's isolated uvx/Python script interpreter
and the shared `isolated_release_runtime` context manager. Asset operations use
the shared `release --isolated-runtime` CLI. The supplier owns managed-tool
selection, checkout-path exclusion, and environment restoration; local recipes
only pass release arguments or invoke the candidate controller.
MSYS applies shell-like argument expansion when started by a native Windows
program, which can corrupt JSON passed as Bash positional arguments; see
[MSYS2 #522](https://github.com/msys2/MSYS2-packages/issues/522). The integration
fixture keeps the real Just, uv and Python launch and substitutes only the
entry point's provider responses. It checks the exact received identity and
rejects an explicitly empty pin rather than selecting fresh evidence.

Semgrep scans and fixture checks share the root Python 3.15 environment and lockfile. Preserve its
native-wheel requirement and compatibility overrides in `pyproject.toml`.
Local scans, fixture checks and SARIF retain the Linux stack-limit workaround
for [Semgrep #11979](https://github.com/semgrep/semgrep/issues/11979).
Ordinary Rust and Python tests remain scan targets; deliberate `tests/semgrep/`
fixtures are checked separately by `just semgrep-test`.
The shared package's `semgrep` extra declares the compatible scanner range.
Setup and `semgrep check-fixtures` verify the installed executable against
`uv.lock` and the root Python selector.
Intel macOS uses Semgrep 1.182.1, the last release with a wheel for that platform.
Other platforms resolve Semgrep 1.184.0 or newer within the shared range.

## Shared capability adoption

The tooling group and isolated `uvx --python 3.15` runtimes use the published
`research-repo-tools==0.2.8` package. Integration checks verify the manifest,
registry lock artifacts and installed version, including the isolated runtimes.
For development rehearsal, an exact full-commit uv Git source override remains
supported; keep the isolated runtime pins aligned and do not float on `main`.
Consumer identity tests derive the selected requirement from the manifest and
use RRT to verify Git declarations against the lock and installed provenance.
Both registry and verified Git selections can run full `just ci`; only the
published selection is eligible for release evidence and publication.
Before an actual release, replace any override and isolated Git pins with the
supplying published version, resolve/sync the lock and revalidate affected
evidence. The configured `release.published-dependencies` policy uses
the shared implementation from
[research-repo-tools#193](https://github.com/acgetchell/research-repo-tools/issues/193).
`just release-dependencies` checks that transition before release work. It
rejects development pins; ordinary checks and bounded acquisition remain
available during rehearsal. Release evidence and publication enforce the
same policy again, including the actual isolated runtime.

The source-built dependencies `actionlint-py`, `pyyaml`, and `ruamel-yaml-clib`
require native preparation with `just release-dependency-builds` before
`just release-dependencies`. The separate locked `build` group supplies their
build tools without adding those tools to `dev`. `dependency-sources.toml`
also pins actionlint's downloaded executable archives for each supported
platform; update those hashes and wheel member paths with `actionlint-py`.
The shared builder verifies the published sources and build tools, builds and
installs native wheels, then verifies the complete installed environment.

Receipts and wheels remain outside the checkout in
`~/.research-repo-tools/source-builds`. Each machine prepares its own receipts;
never copy their signing keys into CI artifacts or caches. Changed manifests,
locks, source policy, interpreter, or builder invalidate them. Rebuild before
the next release check after such changes. These receipts establish observed
source/build provenance, not reproducible builds or publisher-signed wheels.
See the supplier's [source-build contract](https://github.com/acgetchell/research-repo-tools/blob/v0.2.8/README.md#source-built-python-dependencies).

Native CI runs `just ci-dependencies` in the `dependencies` stage of the `native`
job alongside the Rust, Python, and candidate-acceptance stages after shared policy checks. Signed release evidence
requires every native stage. A verified exact development
tooling pin skips this release-only preparation; invalid pins and failed builds
or verification still fail. Release jobs with fresh project environments prepare
their own receipts before acquisition/acceptance verification. Local `just ci`
retains its ordinary validation scope.

Normal `just update` uses the supplier's exact Git-source support from
[research-repo-tools#201](https://github.com/acgetchell/research-repo-tools/pull/201),
which fixes [#200](https://github.com/acgetchell/research-repo-tools/issues/200).
The updater preserves the included tooling group's exact supplier pin, any
immutable source commit and inherited tool versions while updating unrelated
dependencies. Upgrading a published supplier uses the target release's adopter
outside the old project environment. For the current release:

```bash
uvx --no-config --isolated --managed-python --python 3.15 --from research-repo-tools==0.2.8 research-repo-tools toolchain adopt --dry-run
uvx --no-config --isolated --managed-python --python 3.15 --from research-repo-tools==0.2.8 research-repo-tools toolchain adopt --apply
```

The adopter reconciles requirements, lock and environment. Align the declared
source in `dependency-sources.toml`, isolated recipe/workflow pins and consumer
identity tests with that target, then run the normal checks. `just update` does
not perform this supplier migration, including when run by v0.2.8.
For an exact Git development pin, v0.2.8 also supports explicit attested
[promotion](https://github.com/acgetchell/research-repo-tools/blob/v0.2.8/docs/UPDATING_DEPENDENCIES.md#promoting-an-exact-development-pin);
ordinary updates continue to preserve the pin.

For a deliberate source transition, edit `[sources.tooling]` in
`dependency-sources.toml`, then run `just dependency-source tooling` to preview
the shared resolver's manifest/lock changes. Add `--apply` only to apply that
selection. A published source uses `version`; a development source uses
`github`, a full `commit`, and `directory` instead. The helper supports both
Python and Cargo selections. After applying, sync the environment and reconcile
isolated runtime pins and identity tests; the helper does not rewrite workflows
or automatically promote a development pin when a release appears.

The root `dependency-sources.toml` declares both dependency graphs for the
shared source-policy validation introduced by
[research-repo-tools#198](https://github.com/acgetchell/research-repo-tools/pull/198).
Cargo release ownership includes `la-stack` and its private comparison member;
the Python project is the private tooling environment. Other dependencies must
come from their registries for release. The file is retained in source exports
and both benchmark harness declarations, so policy changes invalidate evidence.

The capability inventory below uses the pinned CLI and
its [v0.2.8 contracts](https://github.com/acgetchell/research-repo-tools/blob/v0.2.8/README.md#common-workflows).
The shared updater installs and verifies candidate Cargo tools before publishing
new pins; no local updater implementation is needed.
Native Rust algorithms, assertions, benchmark selection and release decisions
remain owned here.

| Capability | Consumer use |
| --- | --- |
| Actions | `workflow-allowlist-check` gates external workflow dependencies; `update-actions` updates explicitly selected workflow pins |
| Changelog and release | Shared generation, archives, metadata rules, registry gates and verification; the Cargo adapter reuses signed native/package evidence |
| CI environment and toolchain | Shared setup/check/run/export and verified cached installations; `tools-clean` previews obsolete managed tools separately from Cargo `clean` |
| Coverage | Native `cargo llvm-cov` acquisition; `coverage-report` summarizes its existing Cobertura output |
| Dependencies | Shared uv/tool upgrades and `deps update-python --all-projects`; consumer Python dependencies share the root environment and lock |
| Documentation and files | Shared tracked/nonignored selection and bounded argument batches for formatting, lint, spelling and syntax checks; one active Markdown policy preserves archives |
| Performance | Declarative discovery, measurement, suite verification/packaging/assembly, immutable evidence, offline rendering, report/README publication and release assets |
| Python | Shared complete-inventory check/fix/typecheck with repository Ruff exceptions and strict ty diagnostics |
| Review | Shared opt-in CodeRabbit branch and uncommitted review |
| Security and workflow audit | Shared OSV, full-history Gitleaks and zizmor policy; consumer advisory exceptions and workflow permissions remain local |
| Semgrep | Shared fixture validation; native production scans retain `.semgrepignore`, explicit status gates and the reviewed in-source suppression policy |
| Validation | `tooling/examples.toml` selects the published Cargo package and exact-feature overrides; shared Cargo metadata discovery runs every binary example |

The example runner uses compiler-reported executable paths, including configured
native targets and output directories, with per-build and per-execution budgets.
Ordinary examples run with default features; exact examples opt in individually.
Their numerical assertions remain inside the executables. New examples join
automatically; a new required feature needs a matching configuration entry.

Notebook, paper/PDF and Tectonic capabilities do not apply to this repository.
SARIF splitting is unnecessary for the existing individual scanner uploads;
the native suppression filter retains the repository's specific GitHub policy.
Templates are reference material when changing declarations, not generated copies
to maintain here. Raw `files check-lines` intentionally does not replace
`docs check-lines`: the latter preserves the existing Markdown table exemption.

Use `just --list` for arguments. Actions updates and managed-tool removal are
explicit maintenance operations; checking their wiring does not run them.

`.github/actions-allowlist.json` owns exact action identities accepted by local
validation and the updater. Workflows retain full commit pins. This source policy
is separate from the historical GitHub settings payloads under `.github/settings/`,
which also constrain the publishing action to a specific commit. Reconcile live
settings before deploying a pin update; the local updater does not change them.
The allowlist recipe remains scoped to workflows. Semgrep rules check
composite actions under `.github/actions/`, including their
external dependencies and commit pins.

## Performance ownership

The [benchmark guide](../docs/BENCHMARKING.md) owns command selection,
interpretation and output paths. Measurement and publication use the installed
APIs or native Cargo commands. `performance-release.toml` declares preparation,
release measurement and README publication. `release_candidate.py` supplies
repository orchestration for independent owner records and hosted transport. Preparation retains
independent commit assertions and declaration identity; measurement/report
promotion and README publication remain separate retry boundaries.
There is no runtime configuration generator.

| Declaration | Responsibility |
| --- | --- |
| `measurements/full.toml`, `measurements/non-exact.toml` | Release floor, required cases, sampling, original peer phases, source/harness inventories, probes and environments |
| `performance-policy.toml` | Independent scientific catalog used by consumer acceptance tests |
| `performance.just` | Native discovery, correctness gates and measurement commands |
| `criterion-inventory.toml` | Optional saved-baseline discovery and independent required scientific IDs |
| `criterion-saved.toml` | Required IDs, canonical paths, 100 samples and exact saved raw JSON bytes |
| `performance-report.toml` | Independent release eligibility, current report and immutable run history |
| `performance-release.toml` | Preparation storage, publication chronology and separate phase declarations |
| `release_candidate.py` | Owner-controlled capture, final acceptance and exact hosted handoffs using public shared APIs |
| `performance-readme.toml` | Canonical README selection, formatter, source freshness and tagged links |
| `plots/{non-exact,retained,scratch}/*.toml` | Named metric/statistic/axis selections with explicit input ownership |

Both labels of a new measurement must be v0.5.0 or newer. The full and
non-exact declarations require 100 samples, mean and median estimates,
95% intervals and nanoseconds. Each phase runs its native correctness gate before
timing. Schema-6 availability resolves each suite before compilation, inside one
sequential trial. README and LU plots select trial 1. New-API
suites run their column-basis, rational-dot, interval-solve or polynomial gate
on every supported revision. The non-exact diagnostic retains its original gate
and baseline-phase reference libraries.

The expanded #301 policy measures all public la-stack workflows at the current
revision and on supported baselines. Older baselines retain 1,254 shared library
cases; the 149 cases introduced in v0.5.2 become required automatically at that
boundary. The paired LU suite acquires 16 peer cases on each phase alongside
the matching library cases; README ratios use current-phase observations. Hosted
release jobs reuse the accepted local acquisition without collecting timings. The full scratch/retained plot selections
therefore expose only LU solve; broad peer plots remain in the explicitly
requested non-exact diagnostic. The shared reporter prints the resolved suite
inventory and nominal budget before measurement, and `performance-release` records the recipe
body's elapsed seconds and failure status under ignored scratch outputs. The
[public workflow map](../docs/BENCHMARKING.md#release-acquisition-policy) owns
coverage decisions, feature configurations and runtime estimates.

`performance-release` forwards its `progress` and `cadence` arguments directly
to the installed reporter, defaulting to `human` and `300` seconds. JSON status
and `none` are supported without a local progress engine. Consumer fixtures
exercise validated totals, unknown ETA, captured logs and early failure; report
and README recipes retain their measurement-free retry boundaries. The
[workload review](../docs/dev/benchmark-workload-review.md) records the availability
decision, preserved coverage and supplier's corrected gate-failure event.

General workload availability across API revisions is tracked in
[research-repo-tools#166](https://github.com/acgetchell/research-repo-tools/issues/166)
for v0.2.5, with published-package adoption in
[la-stack#320](https://github.com/acgetchell/la-stack/issues/320) for v0.5.3.
The full declaration now uses inclusive `since = "v0.5.2"` suite bounds.
The normally installed 0.2.6 interface was checked against the consumer's peer
comparison contract. Astra recommends retaining its cross-phase different-ID
restriction and migrating to paired LU suites with current-phase peer ratios.
The published v0.2.7 supplier fixes the gate-status defect tracked in
[research-repo-tools#191](https://github.com/acgetchell/research-repo-tools/issues/191);
the preserved schema-2 and schema-5 regressions pass against it.
No permanent baseline exclusion remains; supported missing cases or failed
compilation remain errors.

The current pin also supplies schema-6 interrupted-trial recovery and
phase-specific sources from [PR #205](https://github.com/acgetchell/research-repo-tools/pull/205).
Baseline sources declare historical README compiler input; current gates prove
README is excluded from ordinary debug/release compilation. Current README is
one of the precise generated outputs. Historical current revisions that still
compile it are rejected by the full candidate plan.

Earlier adopted fixes for
[research-repo-tools#196](https://github.com/acgetchell/research-repo-tools/issues/196)
and [#197](https://github.com/acgetchell/research-repo-tools/issues/197) arrived through
[PR #199](https://github.com/acgetchell/research-repo-tools/pull/199).
The paired-suite publisher resolves dependency-version labels from the owning suite
observation; preparation checks each available suite's minimum release before
saving a source pair. The preserved public-API reproducers pass without changing
their assertions. The [workload review](../docs/dev/benchmark-workload-review.md#validation-evidence)
records downstream validation and the remaining publication/platform boundary.

Working-tree measurements include tracked changes and nonignored untracked
regular files. Scratch outputs belong under ignored `target/`.
`just performance-release` fetches tags and saves independent source assertions
before measurement and report promotion. It stops if preparation fails;
`performance-readme` remains a separate publication step. `performance-resume`
reuses the same preparation and checkpoint; it preserves completed suite pairs
inside an interrupted trial. A fresh acquisition refuses an existing checkpoint.
The owner record lives outside the checkout and freezes original bytes/digests
before promotion. Explicit capture retries never start measurements.
For preparation alone, `just performance-prepare` uses local tags;
`--fetch-tags --allow-git-mutations` explicitly enables fetching.
Measurement recipes create temporary
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

The completed v0.5.1/v0.5.0 run predates the #301 acquisition policy and remains
frozen. It is not reclassified as a new-policy run or rerun by these commands.
Current offline rendering requires a complete run acquired under the matching
declaration.

Historical CSV/JSON, plots, reports and complete runs remain frozen with their
original bytes and links. This includes the v0.5.0/v0.4.6 comparison: it lacks
the newer acquisition declaration and is no longer a regeneration input.
Do not invent missing provenance or relax publication checks to replay it.
Future complete runs acquired with the declarations above support offline
report and plot rendering after scratch cleanup.

Shared implementation matrices live upstream, including
`tests/performance/public_complete_consumer.py` and
`public_criterion_consumer.py`, exercised from wheel and sdist installations.
The shared package also owns `public_suites_consumer.py`, which covers suite selection, raw-file
digests, archive paths, original host retention and complete assembly failures.
Legacy suite archive consumer fixtures retain public `pack-suite` and
`assemble-suites` coverage; `bench-release-check-suite` remains diagnostic.
New releases derive archives from the original candidate through shared APIs.
There is no local suite verifier, comparison engine or archive merger.
Consumer tests retain independent scientific selections, actual recipe wiring,
saved archive round trips, publication boundaries and historical hashes.

## Changelog and release tooling

Use the shared commands exposed by `just --list`. The
[release guide](../docs/RELEASING.md) owns metadata, changelog generation,
tagging and registry publication. Declarations in `pyproject.toml` retain
the concept DOI, active link policy, declared date policy, formatter and
authored dependency notes.

The shared changelog command uses the recorded `CITATION.cff` date when
`changelog-release` or a tagged preview
omits the date. `changelog.date-source = "citation"` owns that policy;
explicit dates remain supported. Generation never advances citation metadata.

The publishing declaration selects only `la-stack`, requires native platform
and scientific checks on the exact commit, and requires the durable Criterion
release asset. `release-gate`, `release-publish` and `release-verify`
use the shared release contract; `publish-check` owns native default/exact
Cargo dry-runs. Main CI signs those packaging and native validation identities;
`release-evidence <owner>` binds the original acquisition, precise publisher
receipt, final clean commit and signed native acceptance. Named review records
are immutable; selecting a later native attempt requires an explicit new name.
`release-evidence.toml` owns the package, workflow, native job and signer policy;
`release-assets.toml` owns the draft-target policy. Evidence recipes use the pinned
installed runtime outside the observed checkout. The [release guide](../docs/RELEASING.md#reusable-ci-evidence)
owns the independent owner selection and frozen approval contract.
Tag-triggered `prepare-release.yml` reuses signed native evidence and generated
notes through the installed draft API. Its isolated writer retains the original
journal and captured target. `release-stage` retrieves that original target through
the installed tag-aware receipt verifier and freezes its ID before upload;
`release-publication` adds
canonical archive ID/hash assertions. Two maintainer dispatches use those
independent records. Write-scope hosted jobs execute no checkout code;
read-only verification authenticates full acceptance with an owner-reviewed
SHA-256 and freezes its bytes in an immutable Actions artifact across approval.
Explicit native-evidence renewal changes that reviewed receipt while retaining
the original immutable release archive and asset selections.
`tests/tooling/test_candidate_transfer.py` exercises the actual handoff blocks.
`tests/tooling/test_release_evidence.py` checks installed policy boundaries, while
`tests/tooling/test_release_publishing.py` checks consumer policy and upload barriers
without registry credentials. The release guide owns retention, complete-rerun
recovery and the distinction between reused checks and Cargo repackaging.
`tests/tooling/test_release_flow.py` exercises the actual draft writer with the
installed verifier and draft API, plus canonical `release-benchmarks TAG`
forwarding with automatic original-owner lookup or explicit `CANDIDATE_TRANSFER_FILE`.
Step 1 retains the immutable external owner selection; routine Step 2 recipes
need no `OWNER` or `RELEASE_COMMIT` shell variables. `release-upload` dispatches
the existing reviewed publication selection without approving deployment.
`release-version-check`
requires final metadata; ordinary CI retains a private nonfinal metadata check.
There are no `docs-version-check` or `release-package` recipe aliases.

`tests/tooling/test_toolchain_integration.py` owns the package-identity check:
manifest requirement, registry lock and installed published version, including
the isolated runtime. Negative fixtures retain the development-source rejection
at release boundaries. The shared package pin advances deliberately with its
lockfile and isolated requirements; no sibling checkout supplies runtime code.

Setup, managed execution, dependency updates, Markdown checks, Semgrep,
security scans and opt-in CodeRabbit review also use shared commands.
See the [adoption record](../docs/dev/performance-tooling-adoption.md) for
deletion dispositions and release validation evidence.
