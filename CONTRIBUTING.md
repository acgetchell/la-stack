# Contributing

Thanks for helping improve `la-stack`. This crate is intentionally small and
invariant-heavy, so changes should preserve mathematical correctness, API
clarity, and the fixed-dimension stack-allocation model.

## Getting Started

Install Rust 1.98.0 through [rustup](https://rustup.rs/), Git, Python 3.14,
[`uv` 0.12.5](https://docs.astral.sh/uv/), and `jq`. Install the repository's
pinned `just` version from its locked dependency graph:

```bash
cargo install --locked just --version 1.58.0
```

Set up the remaining development tools and validate the checkout:

```bash
just setup        # install or verify dev tools and sync Python dependencies
just check        # lint and validate without changing files
just ci           # run the comprehensive local CI path
```

Use `just fix` when you intentionally want formatters and automatic fixes to
change files. Run `just --list` for the full command surface.

Use `just update` for deliberate dependency and tool maintenance. It composes
`just update-dependencies`, which advances Cargo dependency requirements and
the Cargo/uv locks, with `just update-cargo-tools`, which upgrades only the
Cargo CLI packages owned by `setup-tools` and atomically reconciles their root
`justfile` pins. The tool updater requires `cargo-install-update` from the
`cargo-update` package and does not touch unrelated Cargo executables or uv's
user-global tool environments.

The repository uses `cargo-nextest` for runnable Rust tests, `cargo-machete`
for unused-dependency checks, and `just cargo-lock-check` to verify that the
committed Cargo lockfile matches the manifest. `rumdl` checks Markdown,
`dprint` plus `yamllint` check YAML and CFF, `taplo` checks TOML, and `typos`
checks spelling. Python support tooling is locked with `uv` and checked by
Ruff, Ty, and Semgrep. GitHub Actions references are SHA-pinned, restricted to
an explicit allowlist, and kept with readable version comments for review.

CI runs `just ci` on Ubuntu, macOS, and Windows to keep platform coverage
aligned with the local comprehensive validation path.

### Rust Toolchain Policy

The `bench-compile` recipe uses Cargo's `CARGO_BUILD_WARNINGS=deny` policy for
both the normal benchmark feature and the `exact` benchmark target. Unlike a
global `RUSTFLAGS='-D warnings'`, this denies lint warnings only for local
packages without changing rustc arguments and creating separate cache
artifacts. Rustdoc has a separate warning contract, so `doc-check` continues to
set `RUSTDOCFLAGS='-D warnings'` for the default and `exact` documentation
surfaces. See the [Rust 1.98 release notes](https://doc.rust-lang.org/stable/releases.html#version-1980-2026-08-20)
and [Cargo configuration reference](https://doc.rust-lang.org/cargo/reference/config.html#buildwarnings).

Rust 1.98's `f32`/`f64::algebraic_{add,sub,mul,div,rem}` operations are
forbidden in `src/`, examples, and benchmark implementations. They permit
real-number transformations such as reassociation with unspecified precision,
which can invalidate the crate's documented error bounds and deterministic
operation order. Their relaxed treatment of non-finite values and signed zero
also conflicts with typed non-finite classification, while varying results
would weaken reproducibility and benchmark comparability. Repository Semgrep
rules enforce this policy; retain ordinary IEEE-754 operations and `mul_add` in
the existing numerical paths.

Any future approximate or fast-math API requires a separate design issue, an
explicit public contract, independent correctness analysis, and representative
benchmarks. It must remain opt-in and must not change existing APIs.

The remaining Rust 1.98 changes require no repository configuration:

- `core::fmt::NumBuffer`, integer `format_into`, `NonZero::from_str_radix`,
  circumfix stripping, UTF-16 decoding, and mutable atomic-slice APIs do not
  simplify an existing path or address a demonstrated bottleneck.
- The compatibility changes for runtime symbols, trait-object lifetimes,
  ambiguous imports, attributes, structural equality, `assert_eq!` temporaries,
  and derived ordering require no source change after comprehensive validation.
- New or promoted targets do not alter the Linux, macOS, and Windows MSVC CI
  matrix. Cargo 1.98's stable changes are fixes for Windows credential-provider
  line endings and diagnostic capitalization; no configuration change is needed.

## Contributor Workflow

Before starting work, check existing GitHub issues for related bug reports,
feature requests, or design discussions. Open an issue before substantial API,
algorithm, invariant, dependency, or performance changes so the expected
behavior and mathematical context can be agreed on first. A focused typo or
similarly mechanical correction does not require advance discussion.

Human contributors should create focused branches. Prefer
`{type}/{issue}-descriptor-or-two`, using the issue number when one exists and
a concise type such as `fix`, `feat`, `perf`, `doc`, `test`, `refactor`, `ci`,
`build`, `chore`, or `style`:

```bash
git switch -c fix/307-exact-conversion
git switch -c perf/315-lu-solve
git switch -c doc/329-branch-guidance
```

Keep each change scoped to one coherent purpose. Update tests and documentation
with the behavior they support, avoid unrelated formatting churn, and cite
relevant literature for numerical or algorithmic work. Automation and AI
assistants must stop before version-control mutations and release operations; a
human contributor performs and reviews commits, pushes, tags, and releases.

### Validation Workflow

Start with the smallest test selection that exercises the changed behavior. Use
`cargo nextest run <test-name>` for a focused runnable test, `cargo test --doc
<filter>` for a focused doctest, or `cargo nextest run --test <crate>` for one
integration-test crate.

For final validation of a non-core change, compose each affected surface once:

- Documentation: `just markdown-ci`
- Configuration and workflows: `just lint-config`
- Python support tooling: `just python-ci`
- Core Rust checks: `just rust-core-check`
- Rust unit tests: `just test-unit`
- Doctests: `just test-doc` and, when relevant, `just test-doc-exact`
- Integration tests: `just test-integration`
- Benchmark inputs or harnesses: `just test-bench-inputs` or `just bench-compile`
- Examples: `just examples`

Run `just ci` for core Rust, public behavior, or GitHub-equivalent validation.
It composes leaf validators directly and runs `clippy-all-targets` to match the
GitHub Clippy SARIF workflow. Unit and integration tests still run together once
through the release-profile `test-rust-ci` bucket, and doctests remain separate
because nextest does not execute them. The test, example, and benchmark buckets
retain their execution or compile-contract roles because ordinary compilation
does not execute Clippy lints.

`la-stack` intentionally has no notebook validation bucket: the repository has
no notebooks or supported Python binding surface. Add notebook tooling only
alongside an actual maintained notebook interface rather than speculatively.

## Project References

Use the existing canonical documents instead of duplicating their guidance:

| Topic | Canonical reference |
|-------|---------------------|
| Agent rules and repository invariants | [`AGENTS.md`](AGENTS.md) |
| User-facing API, examples, and project scope | [`README.md`](README.md) |
| Mathematical basis and numerical validity | [`docs/mathematical_basis.md`](docs/mathematical_basis.md) |
| Package metadata, features, and dependencies | [`Cargo.toml`](Cargo.toml) |
| Commands and validation workflow | [`justfile`](justfile), `just --list` |
| Python support tooling | [`scripts/README.md`](scripts/README.md) |
| Benchmark methodology and baselines | [`docs/BENCHMARKING.md`](docs/BENCHMARKING.md) |
| Coverage workflow and reports | [`docs/COVERAGE.md`](docs/COVERAGE.md) |
| Citations and bibliography | [`CITATION.cff`](CITATION.cff), [`REFERENCES.md`](REFERENCES.md) |
| Security reporting and support | [`SECURITY.md`](SECURITY.md) |
| Releases and changelog generation | [`docs/RELEASING.md`](docs/RELEASING.md), [`CHANGELOG.md`](CHANGELOG.md) |

## Commit Message Format

Use conventional commits so the release tooling can generate useful changelog
entries:

```text
type(scope): short description

- Explain the important behavior or maintenance change.
- Include issue or pull-request references when useful.
```

Common types are `feat`, `fix`, `perf`, `refactor`, `build`, `ci`, `docs`,
`test`, `style`, and `chore`. Mark incompatible public API or behavior changes
explicitly:

```text
feat!: redesign exact conversion API

BREAKING CHANGE: strict conversion now returns a typed unrepresentable error.
```

Pull-request titles should use the same conventional format because merge
commits feed the generated changelog.

## Submitting Changes

Open a pull request with a descriptive conventional title and a concise
summary covering:

- **Problem:** the issue, behavior, or invariant the change addresses.
- **Solution:** how the implementation addresses it and the important design
  choices.
- **Testing:** the validators, tests, and feature combinations that were run.
- **Performance:** comparable before-and-after measurements for
  performance-sensitive work, or why measurement is not applicable.

Use the same machine, toolchain, features, inputs, and benchmark configuration
for before-and-after measurements. The local release comparison is:

```bash
just performance-local
```

It writes `target/bench-reports/performance.md` without changing committed
release documentation or requiring a version bump. Staged and unstaged changes
to tracked files participate; untracked files are excluded, so stage any new
benchmark-relevant file first. Treat regressions as design feedback. If a
slowdown is intentional, explain the correctness, API clarity, or composability
benefit that justifies it.

Core Rust, Cargo, or public-behavior changes must pass `just ci` before a pull
request is ready. Documentation, configuration, Python, test-only,
benchmark-only, and example-only changes use the matching focused validators
documented in [`AGENTS.md`](AGENTS.md). Pull requests are reviewed for correctness,
mathematical accuracy, tests, documentation, style, dependency impact, and
performance. Non-substantive whitespace or formatting churn may be declined
unless it is part of an intentional tooling cleanup.

## Types of Contributions

Bug fixes, new features, documentation, tests, benchmarks, performance work,
and infrastructure improvements are welcome. For algorithmic or numerical
work, update [`REFERENCES.md`](REFERENCES.md) as needed and document the
assumptions, invariants, conditioning behavior, and known limitations.

## AI-Assisted Development

This repository contains an [`AGENTS.md`](AGENTS.md) file, which defines the
canonical rules and invariants for AI coding assistants and autonomous agents
working on this codebase.

AI tools, including ChatGPT, Claude, CodeRabbit, Codex, KiloCode, and WARP, are
expected to read and follow `AGENTS.md` when proposing or applying changes.

Portions of this library were developed with the assistance of these tools:

- [ChatGPT](https://openai.com/chatgpt)
- [Claude](https://www.anthropic.com/claude)
- [CodeRabbit](https://coderabbit.ai/)
- [Codex](https://openai.com/codex/)
- [KiloCode](https://kilocode.ai/)
- [WARP](https://www.warp.dev)

All AI-assisted work must be reviewed and validated by a human maintainer
before it is merged.

For full tool citation metadata, see the
[AI-Assisted Development Tools](REFERENCES.md#ai-assisted-development-tools)
section of [`REFERENCES.md`](REFERENCES.md).

## Release Process

Releases are deliberate, maintainer-driven work. Ordinary feature, fix,
review, and hygiene changes must not update package versions, version-pinned
dependency snippets, citation release dates, generated changelogs, checked-in
release benchmark reports, tags, or other release artifacts. The maintainer
performs every version bump and release manually by
following [`docs/RELEASING.md`](docs/RELEASING.md). Do not substitute an
automated or abbreviated release path.

## Getting Help

Use GitHub Issues for bug reports, feature requests, design questions, and
general project help. Search existing issues before opening a new one. For a
bug, include:

- The crate version or commit and enabled features.
- Rust version, operating system, and relevant development-tool versions.
- A minimal matrix, vector, or code reproduction when possible.
- Expected and actual behavior, including the complete error or panic output.
- The validation commands already run.
- Performance measurements and benchmark configuration when relevant.

Report suspected vulnerabilities privately through the process in
[`SECURITY.md`](SECURITY.md), not in a public issue.
