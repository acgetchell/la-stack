# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2026-10-17

### ⚠️ Breaking Changes

- changelog-unreleased now requires TAG DATE and leaves
  release metadata preparation separate. Replace the removed standalone
  archive-changelog, postprocess-changelog and tag-release entry points with
  research-repo-tools changelog commands. Tags must use v-prefixed SemVer;
  tagging requires today's UTC release date matching CITATION.cff.
- contributor recipes use audit, shell-fix, toml-fix,
  and test-unit. Redundant aliases and help-workflows are removed; use
  bare just to discover the canonical commands.
- la-stack now requires Rust 1.99.0 instead of Rust 1.98.1.

### Merged Pull Requests

- Publish reviewed GitHub releases to crates.io [#282](https://github.com/acgetchell/la-stack/pull/282)
- Record algebraic floating-point rejection [#281](https://github.com/acgetchell/la-stack/pull/281)
- Bump the github-actions group with 6 updates [#280](https://github.com/acgetchell/la-stack/pull/280)
- Bump acgetchell/research-repo-tools/.github/workflows/dependabot-approve.yml [#270](https://github.com/acgetchell/la-stack/pull/270)
- Bump the github-actions group across 1 directory with 5 updates [#267](https://github.com/acgetchell/la-stack/pull/267)

### Added

- Adopt shared CodeRabbit review commands [`c4fa23c`](https://github.com/acgetchell/la-stack/commit/c4fa23c424849a0c3ad2bf2c42e3c3e842ddaacc)

  - Add opt-in branch and uncommitted review recipes using the pinned
    research-repo-tools CLI and repository instructions.
  - Verify the default origin/main base against the live remote, preserve
    local-base overrides, and propagate review output and failures.
  - Document external CLI prerequisites, review scopes, and explicit agent
    invocation outside routine validation gates.
  - Advance the uv pin to 0.12.21 and refresh locked glam and python-dotenv
    dependencies.
- Add stable unsigned vector angles [`09020aa`](https://github.com/acgetchell/la-stack/commit/09020aa18a2b233831a8975b343095d2de1670eb)

  - Add allocation-free borrowed-slice and fixed-vector angle APIs with typed
    shape, zero-vector, and non-finite operand errors.
  - Preserve extreme scales and tiny separations with power-of-two scaling,
    compensated exterior products, and full-angle atan2 evaluation.
  - Add independent numerical regressions and properties, allocation checks,
    stable benchmark comparisons, and retained measurement provenance.
  - Document rounded semantics and the quadratic-time, constant-space tradeoff.
  - Permit Git mutations in agent guidance when explicitly approved by the user.
  - Validate with full just ci: 955 Rust tests, 502 Python tests, default/exact
    doctests, lint checks, benchmark compilation, and examples.

### Changed

- [**breaking**] Adopt pinned shared changelog commands [`ce72843`](https://github.com/acgetchell/la-stack/commit/ce72843d3c7e5c2ed6d97955ed6ce3707acb784f)

  - Delegate generation, normalization, minor archiving, release notes and
    local tag preparation to the published research-repo-tools 0.1.7 CLI.
  - Preserve retained release history and dependency-note links while moving
    completed series to docs/archives/changelog.
  - Keep the historical git-cliff policy as a documented consumer exception
    and validate published changelog links alongside the shared checks.
- Stabilize disposable consumer fixtures [`2372f81`](https://github.com/acgetchell/la-stack/commit/2372f81639838b4f3a767ff021c3ca766aaba63e)

  - Pin the synthetic Cargo identity independently of the checkout's release
    version while retaining the real consumer configuration and policy.
  - Disable inherited commit and tag signing only inside disposable repositories
    so fixture operations do not require user signing keys.
- Adopt common template and draft reviews [`1572957`](https://github.com/acgetchell/la-stack/commit/1572957a4ada1f8fba900ac18f9378e9e0366ba0)

  - Use the shared 0.1.7 changelog format with a bounded exception for historical deps-dev release-note bodies.
  - Preserve literal Rust code during generation and retain complete release history, dates and links.
  - Enable automatic CodeRabbit draft reviews and align zizmor-action with the existing scanner pin.
- Complete shared maintenance adoption [`8551e70`](https://github.com/acgetchell/la-stack/commit/8551e7046c725eca7e415031f09ad61c1859460a)

  - Route setup, managed execution, dependency updates, release metadata,
    and documentation and fixture checks through pinned research-repo-tools.
  - Declare managed Cargo and security tools centrally and reuse shared setup
    in hosted CI while preserving uncached release benchmark jobs.
  - Remove duplicate maintenance helpers and their dedicated tests, retaining
    scientific scripts, benchmark evidence, and consumer policy adapters.
  - Add shared OSV and Gitleaks recipes, upgrade PyJWT with crypto support,
    and narrowly allow a historical compatibility identifier false positive.
  - Document the current shared setup contract and deliberate update policies.
- Adopt shared benchmark primitives [`1c9c66d`](https://github.com/acgetchell/la-stack/commit/1c9c66d75e9c8e46c520372ca35db5435b225600)

  - Delegate Criterion parsing and comparisons, digest verification, archive
    extraction, and README section replacement to research-repo-tools.
  - Publish complete report, evidence, plot, and retained-summary output groups
    through shared transactions while preserving historical schemas and hashes.
  - Report invalid timings and publication failures with context, preserving
    existing artifacts when an operation fails.
  - Remove duplicate tooling and Semgrep rules, refresh subprocess policy,
    and update the Ruff and Semgrep pins.
  - Accept the benchmark-only paste maintenance advisory until January 1, 2027,
    while retaining checks for all other advisories.
- Cover numerical edge cases and preserve benchmark evidence [`1858b5e`](https://github.com/acgetchell/la-stack/commit/1858b5ed11d5c6c1c3c0a3779bc3ada0559016bd)

  - Cover exact zero subtraction, tight square bounds, smallest interval
    determinants, and certified proof loss after underflow or range exhaustion.
  - Check independently derived scalar benchmark endpoint bits before timing.
  - Isolate each benchmark reproduction run's raw Criterion data and logs,
    resolve relative input paths, and document the export workflow.
  - Refresh development-tool pins and synchronize Cargo and Python locks.
- Cover retained report promotion paths [`73dd02f`](https://github.com/acgetchell/la-stack/commit/73dd02f64773f9727236d7a0b095db3fe6e0f186)

  Verify that absent and empty retained-output mappings leave reports
  unchanged, and reject retained summaries that alias the current report
  with the complete path-collision diagnostic before mutation.
- Expose borrowed vector angles through an extension trait [`9b18f05`](https://github.com/acgetchell/la-stack/commit/9b18f05c8ac317edef4df13d94b7e7a1858c3061)

  - Replace the unreleased angle_between function with VectorAngle::angle on
    coordinate slices, retaining Vector::angle and root/prelude exports.
  - Require validated Vector inputs at the internal finite-angle boundary.
  - Expand independent rational-reference and typed error-path coverage.
  - Replace the angle benchmark study with fresh measurements and provenance
    for the extension API, and update API documentation and examples.
  - Refresh locked Rust/Python dependencies and the cargo-nextest tool pin.
- Adopt shared performance and changelog workflows [`279dbdf`](https://github.com/acgetchell/la-stack/commit/279dbdfd073bdbf677dbbc37fa56b47a023b1f25)

  Pin published research-repo-tools 0.1.8 and replace the copied changelog template with the supported dependency-body policy. Delegate common-harness
  measurement, complete run retention, worktrees and transactional publication while retaining la-stack scientific selection and legacy evidence readers.
- Replace Python adapters with shared tooling [`6decbf8`](https://github.com/acgetchell/la-stack/commit/6decbf8e30c5e151738ea148caa8d10797a3a285)

  - Move benchmark commands, measurement inputs, provenance settings, and
    report paths into tooling configuration.
  - Use shared process, Criterion, and Just inspection APIs and remove
    duplicate local adapters and generic tests.
  - Preserve scientific report contracts, historical evidence replay, and
    configurable retained-history selection.
  - Refresh the uv pin and locked tooling dependencies with just update.
- Compare archive symlink targets by filesystem identity [`0d3532d`](https://github.com/acgetchell/la-stack/commit/0d3532d16bc308fd32cc64619bfd0455b448e743)

  - Preserve the observed link target across rejected report publication.
  - Cover absolute and relative symlinks without assuming Windows readlink
    returns the same path representation used when creating the link.
- Standardize complete performance reports [`1f33f14`](https://github.com/acgetchell/la-stack/commit/1f33f14631d95d525a36563f2d1813481a241b9a)

  - Use shared rendering and configured publication for complete-run reports.
  - Preserve dimension plots, historical report bytes, and immutable evidence.
  - Align report guidance with retained provenance and make the Cargo shell
    fixture use LF on Windows.
- Adopt shared performance tooling and Python 3.15 [`e64ff30`](https://github.com/acgetchell/la-stack/commit/e64ff30fdacfef0a3f5ebee8ffe598c928cb72a5)

  - Replace six Python helpers and local packaging with research-repo-tools
    0.1.9 and declarative tooling configuration.
  - Retain only scientific policy and release adapters, with typed catalog
    parsing and source/harness publication checks.
  - Share benchmark discovery, baseline verification, reporting, plotting,
    and release publication while preserving historical artifacts.
  - Retire local console entry points and partial raw-Criterion plotting
    in favor of documented just recipes and complete-run evidence.
  - Target Python 3.15 while isolating Semgrep on locked Python 3.14 for
    local scans, fixtures, and SARIF.
  - Refresh development dependencies and document the migration.

### Dependencies

- Bump development tools and refresh dependency locks [`8a06208`](https://github.com/acgetchell/la-stack/commit/8a06208337a1e5f9d0ad638c45be478ae17e8200)

  - Advance Cargo CLI pins and align the uv pin with version 0.12.19.
  - Upgrade actionlint-py, Ruff, Semgrep, and Ty development requirements.
  - Refresh Cargo and Python locked dependencies while retaining the research-repo-tools 0.1.7 pin.
- Bump the github-actions group across 1 directory with 5 updates [#267](https://github.com/acgetchell/la-stack/pull/267)
  [`6e5a1f0`](https://github.com/acgetchell/la-stack/commit/6e5a1f02297a683f837081464f9c6be881bdcb8d)

  Bumps the github-actions group with 5 updates in the / directory:

  | Package | From | To |
  | --- | --- | --- |
  | [actions-rust-lang/setup-rust-toolchain](https://github.com/actions-rust-lang/setup-rust-toolchain) | `1.17.0` | `2.0.0` |
  | [github/codeql-action/upload-sarif](https://github.com/github/codeql-action) | `4.37.9` | `4.38.2` |
  | [codecov/codecov-action](https://github.com/codecov/codecov-action) | `7.0.0` | `7.1.1` |
  | [github/codeql-action/init](https://github.com/github/codeql-action) | `4.37.9` | `4.38.2` |
  | [github/codeql-action/analyze](https://github.com/github/codeql-action) | `4.37.9` | `4.38.2` |

  Updates `actions-rust-lang/setup-rust-toolchain` from 1.17.0 to 2.0.0

  - [Release notes](https://github.com/actions-rust-lang/setup-rust-toolchain/releases)
  - [Changelog](https://github.com/actions-rust-lang/setup-rust-toolchain/blob/main/CHANGELOG.md)
  - [Commits](https://github.com/actions-rust-lang/setup-rust-toolchain/compare/166cdcfd11aee3cb47222f9ddb555ce30ddb9659...ecabd13d1c56bd1345c230e542e9144811ad706f)

  Updates `github/codeql-action/upload-sarif` from 4.37.9 to 4.38.2
  - [Release notes](https://github.com/github/codeql-action/releases)
  - [Changelog](https://github.com/github/codeql-action/blob/main/CHANGELOG.md)
  - [Commits](https://github.com/github/codeql-action/compare/cdf488f595d80d6e07e03d4674febd5ab45fa938...2892aa5e19bbd11bc0cff5427e3b750a04d9e3c2)

  Updates `codecov/codecov-action` from 7.0.0 to 7.1.1
  - [Release notes](https://github.com/codecov/codecov-action/releases)
  - [Changelog](https://github.com/codecov/codecov-action/blob/main/CHANGELOG.md)
  - [Commits](https://github.com/codecov/codecov-action/compare/fb8b3582c8e4def4969c97caa2f19720cb33a72f...303a32d7a59b442fa8d48b6a1cc6825c09c847a5)

  Updates `github/codeql-action/init` from 4.37.9 to 4.38.2
  - [Release notes](https://github.com/github/codeql-action/releases)
  - [Changelog](https://github.com/github/codeql-action/blob/main/CHANGELOG.md)
  - [Commits](https://github.com/github/codeql-action/compare/cdf488f595d80d6e07e03d4674febd5ab45fa938...2892aa5e19bbd11bc0cff5427e3b750a04d9e3c2)

  Updates `github/codeql-action/analyze` from 4.37.9 to 4.38.2
  - [Release notes](https://github.com/github/codeql-action/releases)
  - [Changelog](https://github.com/github/codeql-action/blob/main/CHANGELOG.md)
  - [Commits](https://github.com/github/codeql-action/compare/cdf488f595d80d6e07e03d4674febd5ab45fa938...2892aa5e19bbd11bc0cff5427e3b750a04d9e3c2)
- Bump acgetchell/research-repo-tools/.github/workflows/dependabot-approve.yml [#270](https://github.com/acgetchell/la-stack/pull/270)
  [`1615b1b`](https://github.com/acgetchell/la-stack/commit/1615b1bcad80a209ced1116e81ffd545e26d85b3)

  Bumps the github-actions group with 1 update:
  [acgetchell/research-repo-tools/.github/workflows/dependabot-approve.yml](https://github.com/acgetchell/research-repo-tools).

  Updates `acgetchell/research-repo-tools/.github/workflows/dependabot-approve.yml` from cbb2ea6dee8866b3f0547bca935aef48fdd71707 to
  b972cbfd717319abaa68d9b1373711ad6f163294

  - [Release notes](https://github.com/acgetchell/research-repo-tools/releases)
  - [Changelog](https://github.com/acgetchell/research-repo-tools/blob/main/CHANGELOG.md)
  - [Commits](https://github.com/acgetchell/research-repo-tools/compare/cbb2ea6dee8866b3f0547bca935aef48fdd71707...b972cbfd717319abaa68d9b1373711ad6f163294)
- Bump the github-actions group with 6 updates [#280](https://github.com/acgetchell/la-stack/pull/280)
  [`d1f85bb`](https://github.com/acgetchell/la-stack/commit/d1f85bb606735d5213a2fa5bd2692e5685a9adc4)

  Bumps the github-actions group with 6 updates:

  | Package | From | To |
  | --- | --- | --- |
  | [github/codeql-action/upload-sarif](https://github.com/github/codeql-action) | `4.38.2` | `4.38.3` |
  | [actions/upload-artifact](https://github.com/actions/upload-artifact) | `7.0.1` | `7.0.2` |
  | [actions/download-artifact](https://github.com/actions/download-artifact) | `8.0.1` | `8.0.2` |
  | [github/codeql-action/init](https://github.com/github/codeql-action) | `4.38.2` | `4.38.3` |
  | [github/codeql-action/analyze](https://github.com/github/codeql-action) | `4.38.2` | `4.38.3` |
  | [astral-sh/setup-uv](https://github.com/astral-sh/setup-uv) | `10.0.1` | `10.2.0` |

  Updates `github/codeql-action/upload-sarif` from 4.38.2 to 4.38.3

  - [Release notes](https://github.com/github/codeql-action/releases)
  - [Changelog](https://github.com/github/codeql-action/blob/main/CHANGELOG.md)
  - [Commits](https://github.com/github/codeql-action/compare/2892aa5e19bbd11bc0cff5427e3b750a04d9e3c2...24c54180a607b1449ed407dd24f251e4e9147c8d)

  Updates `actions/upload-artifact` from 7.0.1 to 7.0.2
  - [Release notes](https://github.com/actions/upload-artifact/releases)
  - [Commits](https://github.com/actions/upload-artifact/compare/043fb46d1a93c77aae656e7c1c64a875d1fc6a0a...cf430e030ddbb5b0abf93d22962f4752f3646cd9)

  Updates `actions/download-artifact` from 8.0.1 to 8.0.2
  - [Release notes](https://github.com/actions/download-artifact/releases)
  - [Commits](https://github.com/actions/download-artifact/compare/3e5f45b2cfb9172054b4087a40e8e0b5a5461e7c...9000827ccba6bdab643e8b6fd33ac0654aef8333)

  Updates `github/codeql-action/init` from 4.38.2 to 4.38.3
  - [Release notes](https://github.com/github/codeql-action/releases)
  - [Changelog](https://github.com/github/codeql-action/blob/main/CHANGELOG.md)
  - [Commits](https://github.com/github/codeql-action/compare/2892aa5e19bbd11bc0cff5427e3b750a04d9e3c2...24c54180a607b1449ed407dd24f251e4e9147c8d)

  Updates `github/codeql-action/analyze` from 4.38.2 to 4.38.3
  - [Release notes](https://github.com/github/codeql-action/releases)
  - [Changelog](https://github.com/github/codeql-action/blob/main/CHANGELOG.md)
  - [Commits](https://github.com/github/codeql-action/compare/2892aa5e19bbd11bc0cff5427e3b750a04d9e3c2...24c54180a607b1449ed407dd24f251e4e9147c8d)

  Updates `astral-sh/setup-uv` from 10.0.1 to 10.2.0
  - [Release notes](https://github.com/astral-sh/setup-uv/releases)
  - [Commits](https://github.com/astral-sh/setup-uv/compare/20cfd1bf945f4377ade1205e4dbc17946fc9a30d...c18668ad3cf93ea998bef934396af7bb5c839dc7)

### Documentation

- Expose all configured workflow badges [`6933c15`](https://github.com/acgetchell/la-stack/commit/6933c15bb2abb1a53b994b6624128ff0b559e1cb)

  - Show CodeQL, zizmor, repository-rule analysis, coverage execution,
    benchmark, release-benchmark, and Dependabot approval status in the README.
  - Keep package metadata and workflow status in separate badge groups.
- Align release roadmap and benchmark descriptions [`fb25684`](https://github.com/acgetchell/la-stack/commit/fb25684ff06612977f1bb2ad86241bb5bfe4ee01)

  - Plan the stable-Rust v0.5.0 release sequence and conditional algebraic
    floating-point adoption, with compatible follow-up work in v0.5.1.
  - Defer expanded const generics until the required Generic Const Arguments
    capabilities are stable.
  - Clarify that research-repo-tools is a contributor development dependency.
  - Describe dense benchmark controls as row-permuted while preserving
    measured inputs and retained performance evidence.
- Record algebraic floating-point rejection [#281](https://github.com/acgetchell/la-stack/pull/281)
  [`a275ddf`](https://github.com/acgetchell/la-stack/commit/a275ddffcd1669eaadc64756ed2bfae60ad879ec)

  - Retain strict arithmetic after the norm reduction prototype hid required overflow errors despite local timing gains.
  - Preserve the numerical counterexample, measurements, provenance, and reproducible experimental evidence.
  - Map floating-point, certified, interval, and exact arithmetic contracts to validation evidence and record closure of the remaining prototype work.
  - Update documentation navigation, release planning, and the reviewed README link inventory.

### Fixed

- Preserve consumer command and release contracts [`cb8c471`](https://github.com/acgetchell/la-stack/commit/cb8c471511bf4f7f5f50cdfe3d7a0296e3357f14)

  - Forward preview options to the shared CLI and keep generated archives under
    the dedicated changelog formatter and final link check.
  - Accept valid included dependency groups while advancing only direct dev
    pins, leaving the shared tooling dependency under its own upgrade policy.
  - Keep declared release dates valid after a later merge and include canonical
    changelog archives in the release staging guidance.
- Preserve zizmor failures and discard invalid SARIF [`477acf1`](https://github.com/acgetchell/la-stack/commit/477acf1882cbf49782c81f7a4ca732c171d5f9df)

  - Preserve the scanner exit status while checking its report before upload.
  - Remove empty or malformed SARIF output so failed scans cannot publish
    an unusable report; retain valid reports containing findings.
- [**breaking**] Adopt shared approvals and generated command help
  [`fec147b`](https://github.com/acgetchell/la-stack/commit/fec147b5d9a385a609ff51720d790378e95fc786)

  - Replace owner-issued CodeRabbit requests with the pinned MCMC
    Dependabot approval workflow and document the repository settings.
  - Fail workflow security checks on findings before exporting SARIF.
  - Generate sorted Just help and keep usage, contributor guidance, and
    active documentation navigation with their respective owners.
  - Replace duplicate process and zizmor implementations with the
    published research-repo-tools APIs.
- Honor inline Semgrep suppressions in uploaded SARIF [`c314554`](https://github.com/acgetchell/la-stack/commit/c314554327bf1521473f2535af4f001f0c27fb24)

  - Remove in-source-suppressed results before GitHub code scanning
    ingests the report so reviewed inline exceptions stay effective.
  - Preserve unsuppressed findings, their metadata, and the scanner's
    original exit status for the blocking security gate.
- Preserve fixture line endings on Windows [`f748db0`](https://github.com/acgetchell/la-stack/commit/f748db0dbc1683184d5e16ee5bd3c00ad07f87e3)

  - Declare LF and CRLF benchmark snapshot inputs explicitly instead of
    inheriting the host's text-write behavior.
  - Check exact worktree bytes and staged contents against the selected
    line ending while preserving tracked-file and index guarantees.
- Require a working Semgrep wheel on Windows [`0e74cc3`](https://github.com/acgetchell/la-stack/commit/0e74cc3d267a3fc48a3c3893911deccc3afddd7c)

  - Restore Semgrep 1.178.0, whose Windows wheel includes semgrep-core.exe.
  - Disable Semgrep source builds so missing platform wheels fail during
    installation instead of leaving an unusable scanner in the environment.
- Restore Python checks with ty 0.0.85 [`0b97b4b`](https://github.com/acgetchell/la-stack/commit/0b97b4be8df4f87ab230be4f836b897e5499e9e8)

  Use Generator annotations for context managers and explicitly distinguish
  absent optional values in performance tooling. Preserve empty baseline
  names and retained-output mappings without suppressing diagnostics.
- Preserve subnormal vector angles on Windows [`7ee4bde`](https://github.com/acgetchell/la-stack/commit/7ee4bdee0855aa7381ba6ace17f847feb37fa3ce)

  - Use the direct exterior/dot ratio for sufficiently small positive angles,
    with a documented arctangent approximation bound, avoiding platform atan2
    underflow observed by the Windows CI regression tests.
  - Test the small-angle transition and normal/subnormal boundary across
    dimensions 2 through 6, 8, and 64 without relaxing existing tiny-angle tests.
  - Verify that exact-conversion error formatting propagates writer failures.
  - Retain the defensive non-finite-result guard without artificial coverage
    tests, exclusions, or coverage-threshold changes.
  - Refresh all 80 benchmark records and retain the original measurements by
    commit for a documented before/after comparison.
  - Validate with just check, just coverage-ci (876 tests), and isolated just ci
    (963 Rust tests, 502 Python tests, doctests, lints, benchmarks, and examples).
- Preserve explicit same-version local comparisons [`25cd362`](https://github.com/acgetchell/la-stack/commit/25cd362ea19e47e711d3e96c4997f64dc894bc8c)

  Keep the consumer-only local comparison exception while release promotion still requires distinct labels and sources remain independently identified. Cover
  explicit selection and current-phase gate failure before either timing phase.
- Retain report test collection and release guidance [`deaecbc`](https://github.com/acgetchell/la-stack/commit/deaecbc81ba0acb575f73b4b93a0bdcdb3aba44b)

  Restore the peer-context test prefix and describe complete-run README inputs, staged retention files, and the supported changelog dependency-body policy.
- Harden shared performance publication and Windows fixtures [`0ef7b78`](https://github.com/acgetchell/la-stack/commit/0ef7b782b9f8fa4f65efc8a4e51797beab93b301)

  - Exclude unsupported norm cases from historical baseline inventories.
  - Preserve prior reports on archive conflicts and reject ambiguous identities.
  - Protect complete-run inputs during README publication and honor configured
    archive locations and equivalent input paths.
  - Keep immutable generated reports outside Markdown formatting and render
    curated provenance with valid list structure.
  - Compare command fixture directories by filesystem identity across Windows
    and POSIX shells, and align the managed rumdl pin to 0.2.79.
- Preserve benchmark publication identity and diagnostics [`86a2cf8`](https://github.com/acgetchell/la-stack/commit/86a2cf82fec9f29a398cb6431f7f0d4aec23882e)

  - Keep README asset provenance tied to validated retained runs across
    scratch cleanup while preserving input freshness checks.
  - Reject non-canonical storage paths for release-report benchmarks.
  - Handle Git execution failures through the plotting CLI diagnostic.
  - Derive the minimum supported release diagnostic from the parsed policy.

### Maintenance

- [**breaking**] Require Rust 1.99.0 [`78c289a`](https://github.com/acgetchell/la-stack/commit/78c289a913360422a6665219f6d5a3fcbebde0a2)

  - Align the library, comparison workspace, contributor toolchain, and Clippy
    baseline with Rust 1.99.0 and LLVM 23.
  - Resolve new lint diagnostics without changing strict floating-point,
    exact-arithmetic, interval, or certified-error contracts.
  - Record paired compiler measurements, code-generation findings, and retained
    evidence as the strict-f64 baseline for future optimization work.
  - Update uv, Semgrep, dprint, and locked Rust/Python dependencies.
- Publish reviewed GitHub releases to crates.io [#282](https://github.com/acgetchell/la-stack/pull/282)
  [`51868ea`](https://github.com/acgetchell/la-stack/commit/51868ea4725c1a4790bba4c59341d2113b18b2bb)

#### ci: publish reviewed GitHub releases to crates.io

- Gate root-crate publication on immutable release identity, exact-commit
  native checks, release notes, benchmark evidence, and deployment approval.
- Use short-lived OIDC credentials, serialize uploads, and verify registry
  visibility without blindly retrying a potentially successful upload.
- Leave benchmark-equipped release drafts for maintainer review and document
  the publishing setup, release sequence, and recovery procedure.
- Keep publishing action updates under manual review so their pins stay
  aligned with the selected-actions policy.

### Performance

- Reduce interval and certified reduction overhead [`dd909a7`](https://github.com/acgetchell/la-stack/commit/dd909a7ee9e8c7471211c88a2139c98895486305)

  - Reuse singleton interval arithmetic, select only needed sum endpoints,
    and size determinant workspaces to the compile-time dimension.
  - Use exact FMA residuals above the proved underflow threshold while
    preserving integer product comparison for smaller products.
  - Compact scalar certificates and reduce magnitude-accumulation overhead
    while preserving ordered estimates, finite endpoints, and signed zero.
  - Certify single normal unit products with zero error, including MAX.
  - Retain reproducible upstream and downstream performance evidence,
    adoption patches, rejected candidates, and measured tradeoffs.
- Accelerate triangular rational solves [`5921577`](https://github.com/acgetchell/la-stack/commit/5921577b94410c24a7565d530755c284cb2f6f89)

  - Recognize upper and lower triangular systems through row permutations
    and solve by exact substitution with one reduction per component.
  - Preserve the general Bareiss path, canonical rational inputs, and typed
    singularity diagnostics without changing the public API.
  - Add permanent structured benchmarks and retain Rust 1.99 measurements:
    11.15x faster prepared solves and 4.26x faster full adapters at D=8
    on the reproduced triangular family.
  - Refresh the staged cc, syn, uv, Ruff, and Ty dependency/tooling pins.
  - Validation: just ci passes, including 1,002 Rust tests, 371 Python
    tests, doctests, examples, and benchmark compilation.

## Archives

Older releases are archived by minor series:

- [0.4.x](docs/archives/changelog/0.4.md)
- [0.3.x](docs/archives/changelog/0.3.md)
- [0.2.x](docs/archives/changelog/0.2.md)
- [0.1.x](docs/archives/changelog/0.1.md)

[0.5.0]: https://github.com/acgetchell/la-stack/compare/v0.4.6...v0.5.0
