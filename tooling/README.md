# Tooling

The installed `research-repo-tools` package owns common development and
performance tooling. The exact released pin and registry hashes live in
`pyproject.toml` and `uv.lock`. This dependency-only environment has no local
Python distribution or console entry points. `release_evidence.py` is the
repository-specific Cargo evidence adapter around the installed validation
helpers; shared tooling is not copied into this repository.

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
runs the installed-package integration tests. `just ci` is the final
repository validation. Native Linux, macOS and Windows CI is required for
adoption; local macOS results do not establish the other platforms.

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

The `release-evidence` recipe uses Just's native uv/Python script interpreter.
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
Published v0.2.4 hard-codes Python 3.14 for those fixture checks. The temporary
`shared_tooling.py` adapter selects the scanner installed beside the active
interpreter and delegates all fixture validation to the shared CLI. It also puts
that environment's scripts first on `PATH` for the native scanner launcher.
Setup uses the same adapter so it does not provision a separate scanner runtime.
The shared Python baseline must apply to Semgrep too; the fix is tracked in
[research-repo-tools#165](https://github.com/acgetchell/research-repo-tools/issues/165)
for v0.2.5. Remove the version-guarded override through
[#317](https://github.com/acgetchell/la-stack/issues/317) after adopting that fix.

## Shared capability adoption

The capability inventory below uses the installed v0.2.4 CLI and
its [versioned contracts](https://github.com/acgetchell/research-repo-tools/blob/v0.2.4/README.md#common-workflows).
Trusted release jobs use the same exact package pin through `uvx --python 3.15`.
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
CLI or native Cargo commands. `release_performance.py` temporarily captures the
shared release selector's independent commit assertions during
`performance-prepare` and reuses them during the separate `performance-readme`
step. `performance-release` remains measurement and report promotion only.
The adapter will be removed when the shared package publishes equivalent commands.
That handoff is tracked in
[research-repo-tools#129](https://github.com/acgetchell/research-repo-tools/issues/129).
There is no runtime configuration generator.

| Declaration | Responsibility |
| --- | --- |
| `measurements/full.toml`, `measurements/non-exact.toml` | Release floor, required cases, sampling, original peer phases, source/harness inventories, probes and environments |
| `performance-policy.toml` | Independent scientific catalog used by consumer acceptance tests |
| `performance.just` | Native discovery, correctness gates and measurement commands |
| `criterion-inventory.toml` | Hosted la-stack-only discovery and independent required scientific IDs |
| `criterion-saved.toml` | Required IDs, canonical paths, 100 samples and exact saved raw JSON bytes |
| `performance-report.toml` | Independent release eligibility, current report and immutable run history |
| `performance-readme.toml` | Canonical README selection, formatter, source freshness and tagged links |
| `plots/{non-exact,retained,scratch}/*.toml` | Named metric/statistic/axis selections with explicit input ownership |

Both labels of a new measurement must be v0.5.0 or newer. The full and
non-exact declarations require 100 samples, mean and median estimates,
95% intervals and nanoseconds. Each phase runs its native correctness gate before
timing. The full `gate-public` adds column bases, rational dot products, interval
solve oracles and polynomial coefficient checks to the diagnostic `gate`.
Reference libraries originate in the baseline phase.

The expanded #301 policy measures all public la-stack workflows at both revisions,
plus only 16 baseline-phase peer LU-factorization-and-solve cases. Hosted
acquisition omits peers entirely. The full scratch/retained plot selections
therefore expose only LU solve; broad peer plots remain in the explicitly
requested non-exact diagnostic. `performance.just` prints the independent local
case budget before measurement, and `performance-release` records the recipe
body's elapsed seconds and failure status under ignored scratch outputs. The
[public workflow map](../docs/BENCHMARKING.md#release-acquisition-policy) owns
coverage decisions, feature configurations and runtime estimates.

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
Published v0.2.4 also owns `public_suites_consumer.py`: suite selection, raw-file
digests, archive paths, original host retention and complete assembly failures.
The release workflow calls `performance pack-suite` and `assemble-suites`
directly; `bench-release-check-suite` is a diagnostic wrapper around
`performance verify-saved --suite`. No local suite verifier or merge engine remains.
Consumer tests retain independent scientific selections, actual recipe wiring,
saved archive round trips, publication boundaries and historical hashes.

## Changelog and release tooling

Use the shared commands exposed by `just --list`. The
[release guide](../docs/RELEASING.md) owns metadata, changelog generation,
tagging and registry publication. Declarations in `pyproject.toml` retain
the concept DOI, active link policy, declared date policy, formatter and
authored dependency notes.

`release_changelog.py` temporarily supplies the recorded `CITATION.cff`
`date-released` when `changelog-release` (also named `changelog-unreleased`) or
a tagged `changelog-preview` omits the date. Explicit dates remain supported.
The adapter never chooses a new date or updates citation metadata; the shared
CLI still owns generation, formatting, archives and publication.
The shared replacement is tracked in
[research-repo-tools#130](https://github.com/acgetchell/research-repo-tools/issues/130).

The publishing declaration selects only `la-stack`, requires native platform
and scientific checks on the exact commit, and requires the durable Criterion
release asset. `release-gate`, `release-publish` and `release-verify`
use the shared release contract; `release-package` owns native default/exact
Cargo dry-runs. Main CI signs those packaging and native validation identities;
`release-evidence` verifies the original main-push attempt and immutable signed
artifact before publication and rechecks its pinned identity after approval.
The adapter uses the installed `ci.export_environment` helper for explicit
`--github-output` export; the [release guide](../docs/RELEASING.md#reusable-ci-evidence)
owns the transport and approval contract. Cargo evidence policy stays local.
`tests/tooling/test_release_evidence.py` checks the adapter's trust boundaries, while
`tests/tooling/test_release_publishing.py` checks consumer policy and upload barriers
without registry credentials. The release guide owns retention, complete-rerun
recovery and the distinction between reused checks and Cargo repackaging.

`tests/tooling/test_toolchain_integration.py` owns the package-identity check:
manifest pin, registry lock, installed version and distribution hashes.
Dependency updates and release tests derive their requirement from the
manifest. The shared package pin advances deliberately with its lockfile;
no sibling checkout supplies runtime code.

Setup, managed execution, dependency updates, Markdown checks, Semgrep,
security scans and opt-in CodeRabbit review also use shared commands.
See the [adoption record](../docs/dev/performance-tooling-adoption.md) for
deletion dispositions and release validation evidence.
