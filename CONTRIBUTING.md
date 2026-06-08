# Contributing

Thanks for helping improve `la-stack`. This crate is intentionally small and
invariant-heavy, so changes should preserve mathematical correctness, API
clarity, and the fixed-dimension stack-allocation model.

## Workflow

```bash
cargo install just
just setup        # install/verify dev tools + sync Python deps
just check        # lint/validate (non-mutating)
just fix          # apply auto-fixes (mutating)
just ci           # lint + tests + examples + bench compile
```

The repository uses Rust-native tooling for documentation and config checks:
`rumdl` for Markdown, `dprint` with `pretty_yaml` for YAML, `taplo` for TOML,
and `typos` for spelling. GitHub Actions references are SHA-pinned, restricted
to an explicit allowlist, and kept with readable version comments for review.

CI runs `just ci` on Ubuntu, macOS, and Windows to keep platform coverage
aligned with the local comprehensive validation path.

## Performance checks

Performance-sensitive changes should compare the current tree against the
latest published release:

```bash
just performance-local
```

This writes `target/bench-reports/performance.md` without changing committed
release docs. Regressions are worth treating as design feedback: if a slowdown
is intentional, document the correctness, API clarity, or composability benefit
that justifies it.

For coverage commands and report locations, see [`docs/COVERAGE.md`](docs/COVERAGE.md).
For benchmark methodology, see [`docs/BENCHMARKING.md`](docs/BENCHMARKING.md).
For the full set of developer commands, run `just --list`.

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

All code was written and/or reviewed and validated by the author.

For full tool citation metadata, see the
[AI-Assisted Development Tools](REFERENCES.md#ai-assisted-development-tools)
section of [`REFERENCES.md`](REFERENCES.md).
