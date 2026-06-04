# shellcheck disable=SC2148
# Justfile for la-stack development workflow
# Install just: https://github.com/casey/just
# Usage: just <command> or just --list

# Use bash with strict error handling for all recipes
set shell := ["bash", "-euo", "pipefail", "-c"]

cargo_nextest_version := "0.9.137"
cargo_llvm_cov_version := "0.8.7"
dprint_version := "0.54.0"
git_cliff_version := "2.13.1"
rumdl_version := "0.2.6"
taplo_version := "0.10.0"
typos_version := "1.47.1"
zizmor_version := "1.25.2"

# Internal helpers: ensure external tooling is installed
_ensure-actionlint:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v uv >/dev/null || { echo "❌ 'uv' not found. Install with the official installer: https://docs.astral.sh/uv/getting-started/installation/"; exit 1; }
    uv run actionlint -version >/dev/null

_ensure-cargo-llvm-cov:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v cargo-llvm-cov >/dev/null; then
        installed_version="$(cargo llvm-cov --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{cargo_llvm_cov_version}}" ]]; then
        echo "❌ 'cargo-llvm-cov' {{cargo_llvm_cov_version}} not found. Install with:"
        echo "   cargo install --locked cargo-llvm-cov --version {{cargo_llvm_cov_version}}"
        exit 1
    fi

_ensure-cargo-nextest:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if cargo nextest --version >/dev/null 2>&1; then
        installed_version="$(cargo nextest --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{cargo_nextest_version}}" ]]; then
        echo "❌ 'cargo-nextest' {{cargo_nextest_version}} not found. Install with:"
        echo "   cargo install --locked cargo-nextest --version {{cargo_nextest_version}}"
        exit 1
    fi

_ensure-dprint:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v dprint >/dev/null; then
        installed_version="$(dprint --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{dprint_version}}" ]]; then
        echo "❌ 'dprint' {{dprint_version}} not found. Install with:"
        echo "   cargo install --locked dprint --version {{dprint_version}}"
        exit 1
    fi

_ensure-git-cliff:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v git-cliff >/dev/null; then
        installed_version="$(git-cliff --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{git_cliff_version}}" ]]; then
        echo "❌ 'git-cliff' {{git_cliff_version}} not found. Install with:"
        echo "   cargo install --locked git-cliff --version {{git_cliff_version}}"
        exit 1
    fi

_ensure-jq:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v jq >/dev/null || { echo "❌ 'jq' not found. See 'just setup' or install: brew install jq"; exit 1; }

_ensure-rumdl:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v rumdl >/dev/null; then
        installed_version="$(rumdl --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{rumdl_version}}" ]]; then
        echo "❌ 'rumdl' {{rumdl_version}} not found. Install with:"
        echo "   cargo install --locked rumdl --version {{rumdl_version}}"
        exit 1
    fi

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
    installed_version=""
    if command -v taplo >/dev/null; then
        installed_version="$(taplo --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{taplo_version}}" ]]; then
        echo "❌ 'taplo' {{taplo_version}} not found. Install with:"
        echo "   cargo install --locked taplo-cli --version {{taplo_version}}"
        exit 1
    fi

# Internal helper: ensure typos-cli is installed
_ensure-typos:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v typos >/dev/null; then
        installed_version="$(typos --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{typos_version}}" ]]; then
        echo "❌ 'typos' {{typos_version}} not found. Install with:"
        echo "   cargo install --locked typos-cli --version {{typos_version}}"
        exit 1
    fi

_ensure-uv:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v uv >/dev/null || { echo "❌ 'uv' not found. See 'just setup' or https://github.com/astral-sh/uv"; exit 1; }

_ensure-zizmor:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v zizmor >/dev/null; then
        installed_version="$(zizmor --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{zizmor_version}}" ]]; then
        echo "❌ 'zizmor' {{zizmor_version}} not found. Install with:"
        echo "   cargo install --locked zizmor --version {{zizmor_version}}"
        exit 1
    fi

# GitHub Actions workflow validation
action-lint: _ensure-actionlint
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '.github/workflows/*.yml' '.github/workflows/*.yaml')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 uv run actionlint
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

# Changelog generation (git-cliff + post-processing + archiving)
changelog: _ensure-git-cliff _ensure-rumdl python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    GIT_CLIFF_OFFLINE=true git-cliff -o CHANGELOG.md
    uv run postprocess-changelog
    uv run archive-changelog
    archive_files=()
    if [ -d docs/archive/changelog ]; then
        while IFS= read -r -d '' file; do
            archive_files+=("$file")
        done < <(find docs/archive/changelog -name '*.md' -print0)
    fi
    if ((${#archive_files[@]})); then
        rumdl fmt --silent CHANGELOG.md "${archive_files[@]}"
    else
        rumdl fmt --silent CHANGELOG.md
    fi

# Prepend unreleased changes to CHANGELOG.md for the given version
changelog-unreleased version: _ensure-git-cliff _ensure-rumdl python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    GIT_CLIFF_OFFLINE=true git-cliff --tag {{version}} -o CHANGELOG.md
    uv run postprocess-changelog
    uv run archive-changelog
    archive_files=()
    if [ -d docs/archive/changelog ]; then
        while IFS= read -r -d '' file; do
            archive_files+=("$file")
        done < <(find docs/archive/changelog -name '*.md' -print0)
    fi
    if ((${#archive_files[@]})); then
        rumdl fmt --silent CHANGELOG.md "${archive_files[@]}"
    else
        rumdl fmt --silent CHANGELOG.md
    fi

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
    #!/usr/bin/env bash
    set -euo pipefail
    cargo build --features exact --examples

    exe_suffix=""
    if [[ "${OS:-}" == "Windows_NT" ]]; then
        exe_suffix=".exe"
    fi

    shopt -s nullglob
    for example_path in examples/*.rs; do
        [[ -f "${example_path}" ]] || continue
        example="${example_path##*/}"
        example="${example%.rs}"
        "target/debug/examples/${example}${exe_suffix}"
    done


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

lint-code: fmt-check clippy doc-check python-check shell-check semgrep semgrep-test

lint-config: validate-json toml-lint toml-fmt-check yaml-check citation-check action-lint zizmor

lint-docs: markdown-check spell-check

# Markdown
markdown-check: _ensure-rumdl
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.md')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 -n100 rumdl check
    else
        echo "No markdown files found to check."
    fi

markdown-fix: _ensure-rumdl
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.md')
    if [ "${#files[@]}" -gt 0 ]; then
        echo "📝 rumdl check --fix (${#files[@]} files)"
        printf '%s\0' "${files[@]}" | xargs -0 -n100 rumdl check --fix
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
    uv run mypy scripts/archive_changelog.py scripts/bench_compare.py scripts/check_semgrep_fixtures.py scripts/criterion_dim_plot.py scripts/tag_release.py scripts/postprocess_changelog.py scripts/subprocess_utils.py

# Repository-owned Semgrep rules for project-specific diagnostics.
semgrep: _ensure-uv
    uv run semgrep --metrics off --error --strict --timeout 30 --config semgrep.yaml .

# Fixture tests for repository-owned Semgrep rules.
semgrep-test: _ensure-uv
    #!/usr/bin/env bash
    set -euo pipefail

    check_semgrep_fixture() {
        target="$1"
        json="$(uv run semgrep scan --metrics off --json --quiet --strict --config semgrep.yaml "$target")"
        SEMGREP_JSON="$json" uv run scripts/check_semgrep_fixtures.py "$target"
    }

    while IFS= read -r -d '' fixture; do
        check_semgrep_fixture "$fixture"
    done < <(find tests/semgrep -type f ! -name '*.fixed' -print0)

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

    have() { command -v "$1" >/dev/null 2>&1; }

    echo "🔧 Ensuring tooling required by just recipes is installed..."
    echo ""
    echo "Ensuring Rust components..."
    if ! have rustup; then
        echo "❌ 'rustup' not found. Install Rust via https://rustup.rs and re-run: just setup-tools"
        exit 1
    fi
    rustup component add clippy rustfmt rust-docs rust-src llvm-tools-preview
    echo ""

    echo "Ensuring cargo tools..."
    cargo_llvm_cov_version="{{cargo_llvm_cov_version}}"
    if ! have cargo-llvm-cov || [[ "$(cargo llvm-cov --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$cargo_llvm_cov_version" ]]; then
        cargo install --locked cargo-llvm-cov --version "$cargo_llvm_cov_version"
    fi
    cargo_nextest_version="{{cargo_nextest_version}}"
    if ! cargo nextest --version >/dev/null 2>&1 || [[ "$(cargo nextest --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$cargo_nextest_version" ]]; then
        cargo install --locked cargo-nextest --version "$cargo_nextest_version"
    fi
    dprint_version="{{dprint_version}}"
    if ! have dprint || [[ "$(dprint --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$dprint_version" ]]; then
        cargo install --locked dprint --version "$dprint_version"
    fi
    git_cliff_version="{{git_cliff_version}}"
    if ! have git-cliff || [[ "$(git-cliff --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$git_cliff_version" ]]; then
        cargo install --locked git-cliff --version "$git_cliff_version"
    fi
    rumdl_version="{{rumdl_version}}"
    if ! have rumdl || [[ "$(rumdl --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$rumdl_version" ]]; then
        cargo install --locked rumdl --version "$rumdl_version"
    fi
    taplo_version="{{taplo_version}}"
    if ! have taplo || [[ "$(taplo --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$taplo_version" ]]; then
        cargo install --locked taplo-cli --version "$taplo_version"
    fi
    typos_version="{{typos_version}}"
    if ! have typos || [[ "$(typos --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$typos_version" ]]; then
        cargo install --locked typos-cli --version "$typos_version"
    fi
    zizmor_version="{{zizmor_version}}"
    if ! have zizmor || [[ "$(zizmor --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$zizmor_version" ]]; then
        cargo install --locked zizmor --version "$zizmor_version"
    fi
    echo ""

    if have uv; then
        echo "Ensuring uv-managed Python tools..."
        uv sync --group dev
        echo ""
    else
        echo "❌ uv missing; cannot install project-managed Python tools."
        echo "Install uv and re-run: just setup-tools"
        exit 1
    fi

    if ! have jq; then
        echo "❌ 'jq' not found. See 'just setup' or install: brew install jq"
        echo ""
    fi

    echo ""
    echo "Verifying required commands are available..."
    missing=0
    for cmd in cargo-llvm-cov dprint git-cliff jq rumdl taplo typos uv zizmor; do
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
    if cargo nextest --version >/dev/null 2>&1; then
        echo "  ✓ cargo nextest"
    else
        echo "  ✗ cargo nextest"
        exit 1
    fi
    uv run actionlint -version >/dev/null
    echo "  ✓ actionlint (uv)"
    uv run semgrep --version >/dev/null
    echo "  ✓ semgrep (uv)"

    echo ""
    echo "✅ Tooling setup complete."

# Shell scripts
shell-check:
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.sh')
    if [ "${#files[@]}" -gt 0 ]; then
        command -v shellcheck >/dev/null || { echo "❌ 'shellcheck' not found."; exit 1; }
        command -v shfmt >/dev/null || { echo "❌ 'shfmt' not found."; exit 1; }
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

shell-fix: shell-fmt

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

# Testing: runnable Rust tests use nextest; rustdoc doctests remain on cargo test.
test: test-lib test-doc

test-all: test test-integration test-exact test-python
    @echo "✅ All tests passed"

test-doc:
    cargo test --doc --verbose

# Tests for the "exact" feature (det_sign_exact + BigRational Bareiss)
test-exact: _ensure-cargo-nextest
    cargo nextest run --features exact --verbose

test-integration: _ensure-cargo-nextest
    cargo nextest run --test '*' --verbose

test-lib: _ensure-cargo-nextest
    cargo nextest run --lib --verbose

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

toml-check: toml-fmt-check toml-lint

toml-fix: toml-fmt

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
yaml-check: _ensure-dprint
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.yml' '*.yaml')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 dprint check
    else
        echo "No YAML files found to check."
    fi

yaml-fix: _ensure-dprint
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        files+=("$file")
    done < <(git ls-files -z '*.yml' '*.yaml')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 dprint fmt
    else
        echo "No YAML files found to format."
    fi

yaml-lint: yaml-check

# Validate CITATION.cff against the Citation File Format schema.
citation-check: _ensure-uv
    uvx --from cffconvert==2.0.0 cffconvert --validate -i CITATION.cff

# GitHub Actions security analysis
zizmor: _ensure-zizmor
    zizmor .github
