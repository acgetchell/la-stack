# shellcheck disable=SC2148
# Justfile for la-stack development workflow
# Install just: https://github.com/casey/just
# Usage: just <command> or just --list

# Use bash with strict error handling for all recipes
set shell := ["bash", "-euo", "pipefail", "-c"]

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
    done < <(git ls-files -z '.github/workflows/*.yaml' '.github/workflows/*.yml')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 actionlint
    else
        echo "No workflow files found to lint."
    fi

# Benchmarks
bench:
    cargo bench

bench-compile:
    cargo bench --no-run

# Build commands
build:
    cargo build

build-release:
    cargo build --release

check:
    cargo check

# CI simulation (matches delaunay's `just ci` shape)
ci: lint test test-integration bench-compile
    @echo "🎯 CI simulation complete!"

clean:
    cargo clean

# Code quality and formatting
clippy:
    cargo clippy --all-targets --all-features -- -D warnings

# Pre-commit workflow: comprehensive validation before committing
# Runs: linting + all Rust tests (lib + doc + integration) + examples
commit-check: lint test-all examples
    @echo "🚀 Ready to commit! All checks passed!"

# Coverage (cargo-tarpaulin)
#
# Common tarpaulin arguments for all coverage runs
# Note: -t 300 sets per-test timeout to 5 minutes (needed for slow CI environments)
_coverage_base_args := '''--exclude-files 'benches/*' --exclude-files 'examples/*' \
  --workspace --lib --tests --all-features \
  -t 300 --verbose --implicit-test-threads'''

# Coverage analysis for local development (HTML output)
coverage:
    #!/usr/bin/env bash
    set -euo pipefail

    if ! command -v cargo-tarpaulin >/dev/null 2>&1; then
        echo "cargo-tarpaulin not found. Install with: cargo install cargo-tarpaulin"
        exit 1
    fi

    mkdir -p target/tarpaulin
    cargo tarpaulin {{_coverage_base_args}} --out Html --output-dir target/tarpaulin
    echo "Coverage report generated: target/tarpaulin/tarpaulin-report.html"

# Coverage analysis for CI (XML output for codecov/codacy)
coverage-ci:
    #!/usr/bin/env bash
    set -euo pipefail

    if ! command -v cargo-tarpaulin >/dev/null 2>&1; then
        echo "cargo-tarpaulin not found. Install with: cargo install cargo-tarpaulin"
        exit 1
    fi

    mkdir -p coverage
    cargo tarpaulin {{_coverage_base_args}} --out Xml --output-dir coverage

# Default recipe shows available commands
default:
    @just --list

doc-check:
    RUSTDOCFLAGS='-D warnings' cargo doc --no-deps

# Examples
examples:
    cargo run --quiet --example det_3x3
    cargo run --quiet --example solve_3x3

fmt:
    cargo fmt --all

fmt-check:
    cargo fmt --check

# Lint groups (delaunay-style)
lint: lint-code lint-docs lint-config

lint-code: fmt-check clippy doc-check

lint-config: validate-json action-lint

lint-docs: markdown-lint spell-check

markdown-lint:
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.md')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 -n100 npx markdownlint --config .markdownlint.json --fix
    else
        echo "No markdown files found to lint."
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

# Testing (delaunay-style split)
# - test: lib + doc tests (fast)
# - test-all: everything in Rust
# - test-integration: tests/ (if present)

test:
    cargo test --lib --verbose
    cargo test --doc --verbose

test-all: test test-integration
    @echo "✅ All tests passed"

test-integration:
    cargo test --tests --verbose

# File validation
validate-json:
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.json')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 -n1 jq empty
    else
        echo "No JSON files found to validate."
    fi
