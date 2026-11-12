# Release and validation integration

This reference maps validation, cache, and release-handoff ownership.
[CONTRIBUTING.md](../../CONTRIBUTING.md#validation-workflow) owns validator
selection; [Releasing](../RELEASING.md) owns the operator procedure.
The [justfile](../../justfile) and [CI workflow](../../.github/workflows/ci.yml)
are the executable authorities.

## Local validation graph

`just ci` runs shared `python-sync` and `tools-check`, then `_ci-policy`,
then the parallel `_ci-rust` and `test-python` lanes. One Just invocation
coalesces shared prerequisites. `just --jobs 1 ci` serializes the same graph.

| Owner | Responsibility |
| --- | --- |
| Shared setup | Locked Python synchronization and declared tool verification |
| `_ci-policy` | Workflow, Markdown/spelling, version/changelog, TOML/YAML/CFF/JSON, Just, Python static, Cargo lock/format, Semgrep/fixtures, unused-dependency and shell checks |
| `_ci-rust` | Clippy, default/exact rustdoc and doctests, native Rust tests, benchmark compilation, exact benchmark test mode and examples |
| `test-python` | Complete Python consumer suite, isolated from root Cargo output |
| `publish-check` | Separate default/exact Cargo package validation |
| Acquisition | Explicit benchmark collection; never a dependency of validation |

Root Cargo output users remain sequential. Python fixture workspaces use
isolated storage. `tests/tooling/test_validation_graph.py` checks complete
validator coverage, shared setup, lane overlap, and failure propagation in
parallel and serial modes. These graph fixtures do not replace native checks.

## Hosted CI stages

Preparation jobs complete policy checks on Linux, macOS, and Windows. After
all preparation jobs succeed, independent native jobs run on each OS.

| Hosted owner | Canonical command and responsibility |
| --- | --- |
| `prepare (OS)` | Tool setup, then `just ci-policy` |
| `native-rust (OS)` | `just ci-rust`: the complete sequential Rust lane |
| `native-python (OS)` | `just ci-python`: Python tests without the `release_candidate` marker, using four workers |
| `native-candidate (OS)` | `just ci-candidate`: candidate capture/acceptance tests, using one worker |
| `native-dependencies (OS)` | `just ci-dependencies`: locked Cargo fetch, native build receipts and dependency-provenance verification |
| `build (OS)` | Aggregate required-check gate; rejects failed, skipped, cancelled or missing preparation/native results |
| `release-package` | Uncached `just ci-dependencies` and `just publish-check` |
| `release-evidence` | Protected main-push signer, requiring aggregate gates and packaging and verifying actual native jobs, steps and OS labels |

The hosted Python selections form a disjoint, complete partition of the local
suite. Required `build (OS)` names are Linux-hosted aggregate gates; native
evidence comes from the actual OS-specific jobs. Manual dispatch can diagnose a
commit but does not satisfy the main-push event requirement for release signing.

## Managed compilation cache

Ordinary CI caches managed Cargo compilation outputs through the absolute
`RESEARCH_REPO_TOOLS_CARGO_BUILD_HOME` override. Cache keys include OS,
architecture, Rust toolchain, tool declarations, and both locks. Exact-version
installation homes are separate and freshly verified. GitHub's cache scopes
separate PR entries from protected main entries.

Preparation populates this cache before native jobs start. Installation homes,
project outputs, and validation results are not shared between jobs. A cache
miss may rebuild tools; a hit cannot skip checks or establish provenance.
Release packaging, evidence signing, draft preparation, benchmark attachment,
and Cargo publication disable caches. Cache configuration alone is not evidence
of a measured speedup.

## Release handoffs

Step 1 retains the benchmark acquisition and one immutable external owner
selection per checkout and tag. Final acceptance binds that original evidence
to the reviewed candidate. Routine commands recover the registered owner;
explicit owner, target, and transfer inputs support recovery.

Tag-triggered preparation verifies native evidence before passing inert notes
and evidence through an authenticated artifact to an isolated draft writer.
Durable receipts identify the original provider run, attempt, artifact, and
draft target. Initial staging uses that target and never discovers a replacement
draft. Attachment consumes retained measurements without collecting new ones.

Publication is a separate maintainer-approved operation. Its workflow rechecks
the candidate and canonical assets, retains an approval receipt, and performs
nonbuilding Cargo upload using authenticated native validation. Published
supplier tooling is required at release boundaries.

Generic verification, draft lifecycle, and retry primitives belong to
research-repo-tools. Scientific selection, repository orchestration, and
consumer acceptance remain local. See the [tooling guide](../../tooling/README.md)
for declarations and tests, and the [development review archive](../archive/development/README.md)
for historical design evidence.
