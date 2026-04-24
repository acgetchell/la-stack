# shellcheck disable=SC2148
# Justfile for la-stack development workflow
# Install just: https://github.com/casey/just
# Usage: just <command> or just --list

# Use bash with strict error handling for all recipes
set shell := ["bash", "-euo", "pipefail", "-c"]

cargo_llvm_cov_version := "0.8.5"

# Internal helpers: ensure external tooling is installed
_ensure-actionlint:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v actionlint >/dev/null || { echo "❌ 'actionlint' not found. See 'just setup' or https://github.com/rhysd/actionlint"; exit 1; }

_ensure-cargo-llvm-cov:
    #!/usr/bin/env bash
    set -euo pipefail
    if ! command -v cargo-llvm-cov >/dev/null; then
        echo "❌ 'cargo-llvm-cov' not found. See 'just setup-tools' or install:"
        echo "   cargo install --locked cargo-llvm-cov --version {{cargo_llvm_cov_version}}"
        exit 1
    fi

_ensure-git-cliff:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v git-cliff >/dev/null || { echo "❌ 'git-cliff' not found. See 'just setup-tools' or install: cargo install git-cliff"; exit 1; }

_ensure-jq:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v jq >/dev/null || { echo "❌ 'jq' not found. See 'just setup' or install: brew install jq"; exit 1; }

_ensure-npx:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v npx >/dev/null || { echo "❌ 'npx' not found. See 'just setup' or install Node.js (for npx tools): https://nodejs.org"; exit 1; }

_ensure-prettier-or-npx:
    #!/usr/bin/env bash
    set -euo pipefail
    if command -v prettier >/dev/null; then
        exit 0
    fi
    command -v npx >/dev/null || {
        echo "❌ Neither 'prettier' nor 'npx' found. Install via npm (recommended): npm i -g prettier"
        echo "   Or install Node.js (for npx): https://nodejs.org"
        exit 1
    }

_ensure-shellcheck:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v shellcheck >/dev/null || { echo "❌ 'shellcheck' not found. See 'just setup' or https://www.shellcheck.net"; exit 1; }

_ensure-shfmt:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v shfmt >/dev/null || { echo "❌ 'shfmt' not found. See 'just setup' or install: brew install shfmt"; exit 1; }

_ensure-taplo:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v taplo >/dev/null || { echo "❌ 'taplo' not found. See 'just setup' or install: brew install taplo (or: cargo install taplo-cli)"; exit 1; }

# Internal helper: ensure typos-cli is installed
_ensure-typos:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v typos >/dev/null || { echo "❌ 'typos' not found. See 'just setup-tools' or install: cargo install typos-cli"; exit 1; }

_ensure-uv:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v uv >/dev/null || { echo "❌ 'uv' not found. See 'just setup' or https://github.com/astral-sh/uv"; exit 1; }

_ensure-yamllint:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v yamllint >/dev/null || { echo "❌ 'yamllint' not found. See 'just setup' or install: brew install yamllint"; exit 1; }

# GitHub Actions workflow validation
action-lint: _ensure-actionlint
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '.github/workflows/*.yml' '.github/workflows/*.yaml')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 actionlint
    else
        echo "No workflow files found to lint."
    fi

# Benchmarks
bench:
    cargo bench --features bench

# Compile benchmarks without running them, treating warnings as errors.
# This catches bench/release-profile-only warnings that won't show up in normal debug-profile runs.
bench-compile:
    RUSTFLAGS='-D warnings' cargo bench --no-run --features bench

# Compare exact-arithmetic benchmarks against a saved baseline.
# Omit the baseline arg to generate a snapshot without comparison.
bench-compare baseline="": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    baseline="{{baseline}}"
    if [ -n "$baseline" ]; then
        uv run bench-compare "$baseline"
    else
        uv run bench-compare
    fi

# Run the exact-arithmetic benchmark suite.
bench-exact:
    cargo bench --features bench,exact --bench exact

# Save a Criterion baseline for the exact-arithmetic benchmarks.
bench-save-baseline tag:
    cargo bench --features bench,exact --bench exact -- --save-baseline {{tag}}

# Bench the la-stack vs nalgebra/faer comparison suite.
bench-vs-linalg filter="":
    #!/usr/bin/env bash
    set -euo pipefail
    filter="{{filter}}"
    if [ -n "$filter" ]; then
        cargo bench --features bench --bench vs_linalg -- "$filter"
    else
        cargo bench --features bench --bench vs_linalg
    fi

# Quick iteration (reduced runtime, no Criterion HTML).
bench-vs-linalg-quick filter="":
    #!/usr/bin/env bash
    set -euo pipefail
    filter="{{filter}}"
    if [ -n "$filter" ]; then
        cargo bench --features bench --bench vs_linalg -- "$filter" --quick --noplot
    else
        cargo bench --features bench --bench vs_linalg -- --quick --noplot
    fi

# Build commands
build:
    cargo build

build-release:
    cargo build --release

# Changelog generation (git-cliff + post-processing)
changelog: _ensure-git-cliff python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    git-cliff -o CHANGELOG.md
    uv run postprocess-changelog

# Prepend unreleased changes to CHANGELOG.md for the given version
changelog-unreleased version: _ensure-git-cliff python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    git-cliff --unreleased --tag {{version}} --prepend CHANGELOG.md
    uv run postprocess-changelog

# Check (non-mutating): run all linters/validators
check: lint
    @echo "✅ Checks complete!"

# Fast compile check (no binary produced)
check-fast:
    cargo check

# CI simulation: comprehensive validation (matches CI expectations)
# Runs: checks + all tests (Rust + Python) + examples + bench compile
ci: check bench-compile test-all examples
    @echo "🎯 CI checks complete!"

# Clean build artifacts
clean:
    cargo clean
    rm -rf target/llvm-cov
    rm -rf coverage

# Code quality and formatting
clippy:
    cargo clippy --workspace --all-targets --all-features -- -D warnings -W clippy::pedantic -W clippy::nursery -W clippy::cargo

# Clippy for the "exact" feature (catches feature-gated lint issues)
clippy-exact:
    cargo clippy --features exact --all-targets -- -D warnings -W clippy::pedantic

# Coverage (cargo-llvm-cov)
#
# Common cargo-llvm-cov arguments for all coverage runs.
_coverage_base_args := '''--features exact \
  --workspace --lib --tests \
  --verbose'''

# Coverage analysis for local development (HTML output)
coverage: _ensure-cargo-llvm-cov
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p target/llvm-cov
    cargo llvm-cov {{_coverage_base_args}} --open --output-dir target/llvm-cov
    echo "Coverage report generated: target/llvm-cov/html/index.html"

# Coverage analysis for CI (XML output for codecov/codacy)
coverage-ci: _ensure-cargo-llvm-cov
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p coverage
    cargo llvm-cov {{_coverage_base_args}} --cobertura --output-path coverage/cobertura.xml

# Default recipe shows available commands
default:
    @just --list

# Documentation build check (includes exact feature for full API coverage)
doc-check:
    RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --features exact

# Examples
examples:
    cargo run --quiet --example det_5x5
    cargo run --quiet --example solve_5x5
    cargo run --quiet --example ldlt_solve_3x3
    cargo run --quiet --example const_det_4x4
    cargo run --quiet --features exact --example exact_det_3x3
    cargo run --quiet --features exact --example exact_sign_3x3
    cargo run --quiet --features exact --example exact_solve_3x3

# Fix (mutating): apply formatters/auto-fixes
fix: toml-fmt fmt python-fix shell-fmt markdown-fix yaml-fix
    @echo "✅ Fixes applied!"

# Rust formatting
fmt:
    cargo fmt --all

fmt-check:
    cargo fmt --all -- --check

help-workflows:
    @echo "Common Just workflows:"
    @echo "  just check             # Run lint/validators (non-mutating)"
    @echo "  just check-fast        # Fast compile check (cargo check)"
    @echo "  just ci                # Full CI simulation (check + tests + examples + bench compile)"
    @echo "  just fix               # Apply formatters/auto-fixes (mutating)"
    @echo "  just setup             # Install/verify dev tools + sync Python deps"
    @echo ""
    @echo "Benchmarks:"
    @echo "  just bench                 # Run benchmarks"
    @echo "  just bench-compile          # Compile benches with warnings-as-errors"
    @echo "  just bench-vs-linalg        # Run vs_linalg bench (optional filter)"
    @echo "  just bench-vs-linalg-quick  # Quick vs_linalg bench (reduced samples)"
    @echo ""
    @echo "Benchmark plotting:"
    @echo "  just plot-vs-linalg         # Plot Criterion results (CSV + SVG)"
    @echo "  just plot-vs-linalg-readme  # Plot + update README benchmark table"
    @echo ""
    @echo "Changelog & releases:"
    @echo "  just changelog              # Regenerate CHANGELOG.md from full history"
    @echo "  just changelog-unreleased <ver>  # Prepend unreleased changes for a version"
    @echo "  just tag <ver>              # Create annotated tag from CHANGELOG.md"
    @echo "  just tag-force <ver>        # Recreate an existing tag"
    @echo ""
    @echo "Setup:"
    @echo "  just setup             # Setup project environment (depends on setup-tools)"
    @echo "  just setup-tools       # Install/verify external tooling (best-effort)"
    @echo ""
    @echo "Testing:"
    @echo "  just coverage          # Generate coverage report (HTML)"
    @echo "  just coverage-ci       # Generate coverage for CI (XML)"
    @echo "  just examples          # Run examples"
    @echo "  just test              # Lib + doc tests (fast)"
    @echo "  just test-all          # All tests (Rust + Python)"
    @echo "  just test-integration  # Integration tests"
    @echo "  just test-python       # Python tests only (pytest)"
    @echo ""
    @echo "Note: Some recipes require external tools. Run 'just setup-tools' (tooling) or 'just setup' (full env) first."

# Lint groups (delaunay-style)
lint: lint-code lint-docs lint-config

lint-code: fmt-check clippy doc-check python-check shell-check

lint-config: validate-json toml-lint toml-fmt-check yaml-lint action-lint

lint-docs: markdown-check spell-check

# Markdown
markdown-check: _ensure-npx
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.md')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 -n100 npx markdownlint --config .markdownlint.json
    else
        echo "No markdown files found to check."
    fi

markdown-fix: _ensure-npx
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.md')
    if [ "${#files[@]}" -gt 0 ]; then
        echo "📝 markdownlint --fix (${#files[@]} files)"
        printf '%s\0' "${files[@]}" | xargs -0 -n100 npx markdownlint --config .markdownlint.json --fix
    else
        echo "No markdown files found to format."
    fi

markdown-lint: markdown-check

# Plot: generate a single time-vs-dimension SVG from Criterion results.
plot-vs-linalg metric="lu_solve" stat="median" sample="new" log_y="false": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--metric "{{metric}}" --stat "{{stat}}" --sample "{{sample}}")
    if [ "{{log_y}}" = "true" ]; then
        args+=(--log-y)
    fi
    uv run criterion-dim-plot "${args[@]}"

# Plot + update the README benchmark table between BENCH_TABLE markers.
plot-vs-linalg-readme metric="lu_solve" stat="median" sample="new" log_y="true": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--metric "{{metric}}" --stat "{{stat}}" --sample "{{sample}}" --update-readme)
    if [ "{{log_y}}" = "true" ]; then
        args+=(--log-y)
    fi
    uv run criterion-dim-plot "${args[@]}"

# Python tooling (uv)
python-check: python-typecheck
    uv run ruff format --check scripts/
    uv run ruff check scripts/

python-fix: python-sync
    uv run ruff check scripts/ --fix
    uv run ruff format scripts/

python-lint: python-check

python-sync: _ensure-uv
    uv sync --group dev

python-typecheck: python-sync
    uv run ty check scripts/
    uv run mypy scripts/bench_compare.py scripts/criterion_dim_plot.py scripts/tag_release.py scripts/postprocess_changelog.py scripts/subprocess_utils.py

# Setup
setup: setup-tools
    #!/usr/bin/env bash
    set -euo pipefail
    echo "Setting up la-stack development environment..."
    echo "Note: Rust toolchain and components managed by rust-toolchain.toml (if present)"
    echo ""

    echo "Installing Python tooling..."
    uv sync --group dev
    echo ""

    echo "Building project..."
    cargo build
    echo "✅ Setup complete! Run 'just help-workflows' to see available commands."

# Development tooling installation (best-effort)
setup-tools:
    #!/usr/bin/env bash
    set -euo pipefail

    echo "🔧 Ensuring tooling required by just recipes is installed..."
    echo ""

    os="$(uname -s || true)"

    have() { command -v "$1" >/dev/null 2>&1; }

    install_with_brew() {
        local formula="$1"
        if brew list --versions "$formula" >/dev/null 2>&1; then
            echo "  ✓ $formula (brew)"
        else
            echo "  ⏳ Installing $formula (brew)..."
            HOMEBREW_NO_AUTO_UPDATE=1 brew install "$formula"
        fi
    }

    if have brew; then
        echo "Using Homebrew (brew) to install missing tools..."
        install_with_brew uv
        install_with_brew jq
        install_with_brew taplo
        install_with_brew yamllint
        install_with_brew shfmt
        install_with_brew shellcheck
        install_with_brew actionlint
        install_with_brew node
        echo ""
    else
        echo "⚠️  'brew' not found; skipping automatic tool installation."
        if [[ "$os" == "Darwin" ]]; then
            echo "Install Homebrew from https://brew.sh (recommended), or install the following tools manually:"
        else
            echo "Install the following tools via your system package manager:"
        fi
        echo "  uv, jq, taplo, yamllint, shfmt, shellcheck, actionlint, node+npx"
        echo ""
    fi

    echo "Ensuring Rust toolchain + components..."
    if ! have rustup; then
        echo "❌ 'rustup' not found. Install Rust via https://rustup.rs and re-run: just setup-tools"
        exit 1
    fi
    rustup component add clippy rustfmt rust-docs rust-src llvm-tools-preview
    echo ""

    echo "Ensuring cargo tools..."
    if ! have samply; then
        echo "  ⏳ Installing samply (cargo)..."
        cargo install --locked samply
    else
        echo "  ✓ samply"
    fi

    if ! have git-cliff; then
        echo "  ⏳ Installing git-cliff (cargo)..."
        cargo install --locked git-cliff
    else
        echo "  ✓ git-cliff"
    fi

    if ! have typos; then
        echo "  ⏳ Installing typos-cli (cargo)..."
        cargo install --locked typos-cli
    else
        echo "  ✓ typos"
    fi

    if ! have cargo-llvm-cov; then
        echo "  ⏳ Installing cargo-llvm-cov {{cargo_llvm_cov_version}} (cargo)..."
        cargo install --locked cargo-llvm-cov --version {{cargo_llvm_cov_version}}
    else
        echo "  ✓ cargo-llvm-cov"
    fi

    echo ""
    echo "Verifying required commands are available..."
    missing=0
    for cmd in uv jq taplo yamllint shfmt shellcheck actionlint node npx typos git-cliff cargo-llvm-cov; do
        if have "$cmd"; then
            echo "  ✓ $cmd"
        else
            echo "  ✗ $cmd"
            missing=1
        fi
    done
    if [ "$missing" -ne 0 ]; then
        echo ""
        echo "❌ Some required tools are still missing."
        echo "Fix the installs above and re-run: just setup-tools"
        exit 1
    fi

    echo ""
    echo "✅ Tooling setup complete."

# Shell scripts
shell-check: _ensure-shellcheck _ensure-shfmt
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.sh')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 -n4 shellcheck -x
        printf '%s\0' "${files[@]}" | xargs -0 shfmt -d
    else
        echo "No shell files found to check."
    fi

shell-fmt: _ensure-shfmt
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.sh')
    if [ "${#files[@]}" -gt 0 ]; then
        echo "🧹 shfmt -w (${#files[@]} files)"
        printf '%s\0' "${files[@]}" | xargs -0 -n1 shfmt -w
    else
        echo "No shell files found to format."
    fi

shell-lint: shell-check

# Spell check (typos)
spell-check: _ensure-typos
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    # Use -z for NUL-delimited output to handle filenames with spaces.
    #
    # Note: For renames/copies, `git status --porcelain -z` emits *two* NUL-separated paths.
    # The ordering can differ depending on the porcelain output, so we read both and
    # spell-check whichever one exists on disk.
    while IFS= read -r -d '' status_line; do
        status="${status_line:0:2}"
        filename="${status_line:3}"

        # For renames/copies, consume the second path token to keep parsing in sync.
        # Prefer the path that exists on disk to avoid passing stale paths to typos.
        if [[ "$status" == *"R"* || "$status" == *"C"* ]]; then
            if IFS= read -r -d '' other_path; then
                if [ ! -e "$filename" ] && [ -e "$other_path" ]; then
                    filename="$other_path"
                fi
            fi
        fi

        # Skip deletions (file may no longer exist).
        if [[ "$status" == *"D"* ]]; then
            continue
        fi

        files+=("$filename")
    done < <(git status --porcelain -z --ignored=no)
    if [ "${#files[@]}" -gt 0 ]; then
        # Exclude typos.toml itself: it intentionally contains allowlisted fragments.
        printf '%s\0' "${files[@]}" | xargs -0 -n100 typos --config typos.toml --force-exclude --exclude typos.toml --
    else
        echo "No modified files to spell-check."
    fi

# Create an annotated git tag from the CHANGELOG.md section for the given version
tag version: python-sync
    uv run tag-release {{version}}

# Recreate an existing tag (delete + recreate)
tag-force version: python-sync
    uv run tag-release {{version}} --force

# Testing
# test: runs only lib and doc tests (fast)
test:
    cargo test --lib --verbose
    cargo test --doc --verbose

test-all: test test-integration test-exact test-python
    @echo "✅ All tests passed"

# Tests for the "exact" feature (det_sign_exact + BigRational Bareiss)
test-exact:
    cargo test --features exact --verbose

test-integration:
    cargo test --tests --verbose

test-python: python-sync
    uv run pytest -q

# TOML
toml-fmt: _ensure-taplo
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.toml')
    if [ "${#files[@]}" -gt 0 ]; then
        taplo fmt "${files[@]}"
    else
        echo "No TOML files found to format."
    fi

toml-fmt-check: _ensure-taplo
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.toml')
    if [ "${#files[@]}" -gt 0 ]; then
        taplo fmt --check "${files[@]}"
    else
        echo "No TOML files found to check."
    fi

toml-lint: _ensure-taplo
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.toml')
    if [ "${#files[@]}" -gt 0 ]; then
        taplo lint "${files[@]}"
    else
        echo "No TOML files found to lint."
    fi

# File validation
validate-json: _ensure-jq
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

# YAML
yaml-fix: _ensure-prettier-or-npx
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.yml' '*.yaml')
    if [ "${#files[@]}" -gt 0 ]; then
        echo "📝 prettier --write (YAML, ${#files[@]} files)"

        cmd=()
        if command -v prettier >/dev/null; then
            cmd=(prettier --write --print-width 120)
        elif command -v npx >/dev/null; then
            # Prefer non-interactive installs when supported (newer npm/npx).
            # NOTE: With `set -u`, expanding an empty array like "${arr[@]}" can error on older bash.
            cmd=(npx)
            if npx --help 2>&1 | grep -q -- '--yes'; then
                cmd+=(--yes)
            fi
            cmd+=(prettier --write --print-width 120)
        else
            echo "❌ 'prettier' not found. Install via npm (recommended): npm i -g prettier"
            echo "   Or install Node.js (for npx): https://nodejs.org"
            exit 1
        fi

        # Use CLI flags instead of a repo-wide prettier config: keeps the scope to YAML only.
        printf '%s\0' "${files[@]}" | xargs -0 -n100 "${cmd[@]}"
    else
        echo "No YAML files found to format."
    fi

yaml-lint: _ensure-yamllint
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.yml' '*.yaml')
    if [ "${#files[@]}" -gt 0 ]; then
        echo "🔍 yamllint (${#files[@]} files)"
        yamllint --strict -c .yamllint "${files[@]}"
    else
        echo "No YAML files found to lint."
    fi
