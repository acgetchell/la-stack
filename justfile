# shellcheck disable=SC2148
# Justfile for la-stack development workflow
# Install just: https://github.com/casey/just
# Usage: just <command> or just --list

# Use bash with strict error handling for all recipes
set shell := ["bash", "-euo", "pipefail", "-c"]

# Default recipe shows available commands
default:
    @just --list

# Build commands
build:
    cargo build

build-release:
    cargo build --release

check:
    cargo check

# Code quality and formatting
fmt:
    cargo fmt

fmt-check:
    cargo fmt --check

clippy:
    cargo clippy --all-targets --all-features -- -D warnings

doc-check:
    RUSTDOCFLAGS='-D warnings' cargo doc --no-deps

# Lint groups (delaunay-style)
lint: lint-code lint-docs lint-config

lint-code: fmt-check clippy doc-check

lint-docs: spell-check

lint-config: action-lint

# Testing (delaunay-style split)
# - test: lib + doc tests (fast)
# - test-integration: tests/ (if present)
# - test-all: everything in Rust

test:
    cargo test --lib --verbose
    cargo test --doc --verbose

test-integration:
    cargo test --tests --verbose

test-all: test test-integration
    @echo "✅ All tests passed"

# Benchmarks
bench:
    cargo bench

bench-compile:
    cargo bench --no-run

# Examples
examples:
    cargo run --quiet --example det_3x3
    cargo run --quiet --example solve_3x3

# CI simulation (matches delaunay's `just ci` shape)
ci: lint test test-integration bench-compile
    @echo "🎯 CI simulation complete!"

# Pre-commit workflow: comprehensive validation before committing
# Runs: linting + all Rust tests (lib + doc + integration) + examples
commit-check: lint test-all examples
    @echo "🚀 Ready to commit! All checks passed!"

# GitHub Actions workflow validation (optional)
action-lint:
    #!/usr/bin/env bash
    set -euo pipefail
    if ! command -v actionlint >/dev/null; then
        echo "⚠️ 'actionlint' not found. Install: https://github.com/rhysd/actionlint"
        exit 0
    fi
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '.github/workflows/*.yml' '.github/workflows/*.yaml')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 actionlint
    else
        echo "No workflow files found to lint."
    fi

# Spell check (cspell)
#
# Requires either:
# - `cspell` on PATH (recommended: `npm i -g cspell`), or
# - `npx` (will run cspell without a global install)
spell-check:
    #!/usr/bin/env bash
    set -euo pipefail
    if command -v cspell >/dev/null; then
        cspell lint --config cspell.json --no-progress --gitignore --cache --exclude cspell.json .
    elif command -v npx >/dev/null; then
        npx cspell lint --config cspell.json --no-progress --gitignore --cache --exclude cspell.json .
    else
        echo "❌ cspell not found. Install via npm (recommended): npm i -g cspell"
        echo "   Or ensure npx is available (Node.js)."
        exit 1
    fi

clean:
    cargo clean
