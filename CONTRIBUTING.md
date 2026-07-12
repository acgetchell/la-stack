# Contributing

Thanks for helping improve `la-stack`. This crate is intentionally small and
invariant-heavy, so changes should preserve mathematical correctness, API
clarity, and the fixed-dimension stack-allocation model.

## Getting Started

Install Rust 1.97.0 through [rustup](https://rustup.rs/), Git, Python 3.14,
[`uv` 0.11.28](https://docs.astral.sh/uv/), and `jq`. Install the repository's
pinned `just` version from its locked dependency graph:

```bash
cargo install --locked just --version 1.56.0
```

Set up the remaining development tools and validate the checkout:

```bash
just setup        # install or verify dev tools and sync Python dependencies
just check        # lint and validate without changing files
just ci           # run the comprehensive local CI path
```

Use `just fix` when you intentionally want formatters and automatic fixes to
change files. Run `just --list` for the full command surface.

The repository uses `cargo-nextest` for runnable Rust tests, `cargo-machete`
for unused-dependency checks, and `just cargo-lock-check` to verify that the
committed Cargo lockfile matches the manifest. `rumdl` checks Markdown,
`dprint` plus `yamllint` check YAML and CFF, `taplo` checks TOML, and `typos`
checks spelling. Python support tooling is locked with `uv` and checked by
Ruff, Ty, and Semgrep. GitHub Actions references are SHA-pinned, restricted to
an explicit allowlist, and kept with readable version comments for review.

CI runs `just ci` on Ubuntu, macOS, and Windows to keep platform coverage
aligned with the local comprehensive validation path.

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

## Project References

Use the existing canonical documents instead of duplicating their guidance:

| Topic | Canonical reference |
|-------|---------------------|
| Agent rules and repository invariants | [`AGENTS.md`](AGENTS.md) |
| User-facing API, examples, and project scope | [`README.md`](README.md) |
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
release documentation. Treat regressions as design feedback. If a slowdown is
intentional, explain the correctness, API clarity, or composability benefit
that justifies it.

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
