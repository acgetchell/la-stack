# shellcheck disable=SC2148
# Justfile for la-stack development workflow
# Install just: https://github.com/casey/just
# Usage: just <command> or just --list

# Use bash with strict error handling for all recipes
set shell := ["bash", "-euo", "pipefail", "-c"]

home_dir := env_var_or_default("HOME", env_var_or_default("USERPROFILE", ""))
cargo_home := env_var_or_default("CARGO_HOME", home_dir + "/.cargo")
path_separator := if os_family() == "windows" { ";" } else { ":" }
export PATH := cargo_home + "/bin" + path_separator + env_var("PATH")

cargo_machete_version := "0.9.2"
cargo_nextest_version := "0.9.140"
cargo_llvm_cov_version := "0.8.7"
dprint_version := "0.55.1"
git_cliff_version := "2.13.1"
just_version := "1.56.0"
rumdl_version := "0.2.30"
taplo_version := "0.10.0"
typos_version := "1.48.0"
uv_version := "0.11.28"
zizmor_version := "1.26.1"

# Internal helpers: ensure external tooling is installed
_ensure-actionlint:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v uv >/dev/null || { echo "❌ 'uv' not found. Install with the official installer: https://docs.astral.sh/uv/getting-started/installation/"; exit 1; }
    uv run --locked actionlint -version >/dev/null

_ensure-cargo-llvm-cov:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v cargo-llvm-cov >/dev/null; then
        installed_version="$(cargo llvm-cov --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{ cargo_llvm_cov_version }}" ]]; then
        echo "❌ 'cargo-llvm-cov' {{ cargo_llvm_cov_version }} not found. Install with:"
        echo "   cargo install --locked cargo-llvm-cov --version {{ cargo_llvm_cov_version }}"
        exit 1
    fi

_ensure-cargo-machete:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if cargo machete --version >/dev/null 2>&1; then
        installed_version="$(cargo machete --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{ cargo_machete_version }}" ]]; then
        echo "❌ 'cargo-machete' {{ cargo_machete_version }} not found. Install with:"
        echo "   cargo install --locked cargo-machete --version {{ cargo_machete_version }}"
        exit 1
    fi

_ensure-cargo-nextest:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if cargo nextest --version >/dev/null 2>&1; then
        installed_version="$(cargo nextest --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{ cargo_nextest_version }}" ]]; then
        echo "❌ 'cargo-nextest' {{ cargo_nextest_version }} not found. Install with:"
        echo "   cargo install --locked cargo-nextest --version {{ cargo_nextest_version }}"
        exit 1
    fi

_ensure-dprint:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v dprint >/dev/null; then
        installed_version="$(dprint --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{ dprint_version }}" ]]; then
        echo "❌ 'dprint' {{ dprint_version }} not found. Install with:"
        echo "   cargo install --locked dprint --version {{ dprint_version }}"
        exit 1
    fi

_ensure-git-cliff:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v git-cliff >/dev/null; then
        installed_version="$(git-cliff --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{ git_cliff_version }}" ]]; then
        echo "❌ 'git-cliff' {{ git_cliff_version }} not found. Install with:"
        echo "   cargo install --locked git-cliff --version {{ git_cliff_version }}"
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
    if [[ "$installed_version" != "{{ rumdl_version }}" ]]; then
        echo "❌ 'rumdl' {{ rumdl_version }} not found. Install with:"
        echo "   cargo install --locked rumdl --version {{ rumdl_version }}"
        exit 1
    fi

_ensure-shellcheck:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v uv >/dev/null || { echo "❌ 'uv' not found. See 'just setup' or https://docs.astral.sh/uv/"; exit 1; }
    uv run --locked shellcheck --version >/dev/null

_ensure-shfmt:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v uv >/dev/null || { echo "❌ 'uv' not found. See 'just setup' or https://docs.astral.sh/uv/"; exit 1; }
    uv run --locked shfmt --version >/dev/null

_ensure-taplo:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v taplo >/dev/null; then
        installed_version="$(taplo --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{ taplo_version }}" ]]; then
        echo "❌ 'taplo' {{ taplo_version }} not found. Install with:"
        echo "   cargo install --locked taplo-cli --version {{ taplo_version }}"
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
    if [[ "$installed_version" != "{{ typos_version }}" ]]; then
        echo "❌ 'typos' {{ typos_version }} not found. Install with:"
        echo "   cargo install --locked typos-cli --version {{ typos_version }}"
        exit 1
    fi

_ensure-uv:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v uv >/dev/null || { echo "❌ 'uv' not found. See 'just setup' or https://github.com/astral-sh/uv"; exit 1; }

_ensure-yamllint:
    #!/usr/bin/env bash
    set -euo pipefail
    command -v uv >/dev/null || { echo "❌ 'uv' not found. See 'just setup' or https://docs.astral.sh/uv/"; exit 1; }
    uv run --locked yamllint --version >/dev/null

_ensure-zizmor:
    #!/usr/bin/env bash
    set -euo pipefail
    installed_version=""
    if command -v zizmor >/dev/null; then
        installed_version="$(zizmor --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
    fi
    if [[ "$installed_version" != "{{ zizmor_version }}" ]]; then
        echo "❌ 'zizmor' {{ zizmor_version }} not found. Install with:"
        echo "   cargo install --locked zizmor --version {{ zizmor_version }}"
        exit 1
    fi

# GitHub Actions workflow validation
action-lint: _ensure-actionlint
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '.github/workflows/*.yml' '.github/workflows/*.yaml')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 uv run --locked actionlint
    else
        echo "No workflow files found to lint."
    fi

# Benchmarks
bench:
    cargo bench --locked --features bench

# Compare latest measurements against a saved baseline.
# Defaults to the `last` full-release baseline.
bench-compare baseline="last" suite="all" scope="release-signal": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    baseline="{{ baseline }}"
    uv run --locked bench-compare "$baseline" --suite "{{ suite }}" --scope "{{ scope }}"

# Compile benchmarks without running them, treating warnings as errors.
# This catches bench/release-profile-only warnings that won't show up in normal debug-profile runs.
bench-compile:
    RUSTFLAGS='-D warnings' cargo bench --locked --no-run --features bench
    RUSTFLAGS='-D warnings' cargo bench --locked --no-run --features bench,exact --bench exact

# Run the exact-arithmetic benchmark suite.
bench-exact:
    cargo bench --locked --features bench,exact --bench exact

# Run the cheaper latest measurements used for latest-vs-last reports.
bench-latest: bench-vs-linalg-la-stack bench-exact

# Run latest measurements and render the latest-vs-last performance report.
bench-latest-vs-last baseline="last": bench-latest python-sync
    uv run --locked bench-compare {{ baseline }}

# Run only la-stack vs_linalg measurements and render a non-exact performance report.
bench-vs-linalg-latest-vs baseline="last": bench-vs-linalg-la-stack python-sync
    uv run --locked bench-compare {{ baseline }} --suite vs_linalg --scope release-signal

# Save a Criterion baseline. Defaults to all release-signal benchmark suites.
bench-save-baseline tag suite="all":
    #!/usr/bin/env bash
    set -euo pipefail
    suite="{{ suite }}"
    case "$suite" in
        all)
            cargo bench --locked --features bench --bench vs_linalg -- --save-baseline {{ tag }}
            cargo bench --locked --features bench,exact --bench exact -- --save-baseline {{ tag }}
            ;;
        exact)
            cargo bench --locked --features bench,exact --bench exact -- --save-baseline {{ tag }}
            ;;
        vs_linalg)
            cargo bench --locked --features bench --bench vs_linalg -- --save-baseline {{ tag }}
            ;;
        *)
            echo "unknown benchmark suite: $suite" >&2
            exit 2
            ;;
    esac

# Save a full Criterion baseline for the previous release signal.
bench-save-last:
    just bench-save-baseline last

# Bench the la-stack vs nalgebra/faer comparison suite.
bench-vs-linalg filter="":
    #!/usr/bin/env bash
    set -euo pipefail
    filter="{{ filter }}"
    if [ -n "$filter" ]; then
        cargo bench --locked --features bench --bench vs_linalg -- "$filter"
    else
        cargo bench --locked --features bench --bench vs_linalg
    fi

# Bench only la-stack rows from the vs_linalg suite for cheap latest-vs-last comparisons.
bench-vs-linalg-la-stack:
    cargo bench --locked --features bench --bench vs_linalg -- la_stack

# Quick iteration (reduced runtime, no Criterion HTML).
bench-vs-linalg-quick filter="":
    #!/usr/bin/env bash
    set -euo pipefail
    filter="{{ filter }}"
    if [ -n "$filter" ]; then
        cargo bench --locked --features bench --bench vs_linalg -- "$filter" --quick --noplot
    else
        cargo bench --locked --features bench --bench vs_linalg -- --quick --noplot
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
    uv run --locked postprocess-changelog
    uv run --locked archive-changelog
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
    GIT_CLIFF_OFFLINE=true git-cliff --tag {{ version }} -o CHANGELOG.md
    uv run --locked postprocess-changelog
    uv run --locked archive-changelog
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

# Verify Cargo.toml and the committed Cargo.lock are synchronized.
cargo-lock-check:
    cargo metadata --locked --format-version 1 --no-deps > /dev/null

# CI simulation: comprehensive validation (matches CI expectations)
# Runs: checks + all tests (Rust + Python) + examples + bench compile
ci: check bench-compile test-all examples
    @echo "🎯 CI checks complete!"

# Validate CITATION.cff against the Citation File Format schema.
citation-check: _ensure-uv
    uvx --from cffconvert==2.0.0 cffconvert --validate -i CITATION.cff

# Clean build artifacts
clean:
    cargo clean
    rm -rf target/llvm-cov
    rm -rf coverage

# Code quality and formatting
clippy: clippy-all-targets

clippy-all-targets:
    cargo clippy --workspace --all-targets -- -D warnings -W clippy::pedantic -W clippy::nursery -W clippy::cargo
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
coverage: _ensure-cargo-llvm-cov _ensure-cargo-nextest
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p target/llvm-cov
    cargo llvm-cov nextest {{ _coverage_base_args }} --open --output-dir target/llvm-cov -P coverage
    echo "Coverage report generated: target/llvm-cov/html/index.html"

# Coverage analysis for CI (XML output for Codecov)
coverage-ci: _ensure-cargo-llvm-cov _ensure-cargo-nextest
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p coverage
    cargo llvm-cov nextest {{ _coverage_base_args }} --cobertura --output-path coverage/cobertura.xml -P coverage

# Default recipe shows available commands
default:
    @just --list

# Documentation build checks for the default and exact-feature public APIs.
doc-check:
    RUSTDOCFLAGS='-D warnings' cargo doc --no-deps
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
    @echo "  just bench-latest           # Run cheap latest measurements"
    @echo "  just bench-latest-vs-last   # Run latest and compare against last"
    @echo "  just bench-vs-linalg-latest-vs # Run non-exact latest and compare against last"
    @echo "  just performance-github-assets # Compare stored GitHub Actions release assets"
    @echo "  just performance-local      # Compare current tree against latest release locally"
    @echo "  just performance-local-vs-linalg # Compare current non-exact kernels locally"
    @echo "  just performance-release    # Promote local release performance docs"
    @echo "  just bench-save-last        # Save full baseline as 'last'"
    @echo "  just bench-vs-linalg        # Run vs_linalg bench (optional filter)"
    @echo "  just bench-vs-linalg-la-stack # Run la-stack rows from vs_linalg"
    @echo "  just bench-vs-linalg-quick  # Quick vs_linalg bench (reduced samples)"
    @echo ""
    @echo "Benchmark plotting:"
    @echo "  just plot-vs-linalg         # Plot Criterion results (CSV + SVG + provenance)"
    @echo "  just plot-vs-linalg-readme  # Gate, rerun, and publish canonical README assets/table"
    @echo ""
    @echo "Changelog & releases:"
    @echo "  just changelog              # Regenerate CHANGELOG.md from full history"
    @echo "  just changelog-unreleased <ver>  # Prepend unreleased changes for a version"
    @echo "  just tag <ver>              # Create annotated tag from CHANGELOG.md"
    @echo "  just tag-force <ver>        # Recreate an existing tag"
    @echo ""
    @echo "Setup:"
    @echo "  just setup             # Setup project environment (depends on setup-tools)"
    @echo "  just setup-tools       # Install/verify external tooling"
    @echo ""
    @echo "Testing:"
    @echo "  just coverage          # Generate coverage report (HTML)"
    @echo "  just coverage-ci       # Generate coverage for CI (XML)"
    @echo "  just examples          # Run examples"
    @echo "  just test              # Lib + doc tests (fast)"
    @echo "  just test-all          # All tests (Rust + Python)"
    @echo "  just test-bench-inputs # Benchmark input smoke tests"
    @echo "  just test-exact        # Exact-feature tests and doctests"
    @echo "  just test-integration  # Integration tests"
    @echo "  just test-python       # Python tests only (pytest)"
    @echo ""
    @echo "Note: Some recipes require external tools. Run 'just setup-tools' (tooling) or 'just setup' (full env) first."

# Lint groups (delaunay-style)
lint: lint-code lint-docs lint-config

lint-code: rust-core-check python-check shell-check

lint-config: json-check toml-ci yaml-ci github-actions-check justfile-fmt-check

lint-docs: markdown-ci

github-actions-check: action-lint zizmor
    @echo "✅ GitHub Actions checks complete!"

markdown-ci: markdown-check spell-check
    @echo "✅ Markdown checks complete!"

python-ci: python-check test-python
    @echo "✅ Python checks complete!"

rust-core-check: cargo-lock-check fmt-check clippy-all-targets doc-check semgrep semgrep-test unused-deps
    @echo "✅ Rust core checks complete!"

toml-ci: toml-check
    @echo "✅ TOML checks complete!"

yaml-ci: yaml-check citation-check
    @echo "✅ YAML/CFF checks complete!"

# Markdown
markdown-check: _ensure-rumdl
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        case "$file" in
            CHANGELOG.md|docs/archive/*) continue ;;
        esac
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.md')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 -n100 rumdl check
        violations=0
        for file in "${files[@]}"; do
            line_number=0
            while IFS= read -r line || [[ -n "$line" ]]; do
                line_number=$((line_number + 1))
                if [ "${#line}" -gt 160 ]; then
                    printf '%s:%d: line length %d exceeds 160\n' "$file" "$line_number" "${#line}" >&2
                    violations=$((violations + 1))
                fi
            done < "$file"
        done
        if [ "$violations" -gt 0 ]; then
            echo "Markdown raw line-length check failed." >&2
            exit 1
        fi
    else
        echo "No markdown files found to check."
    fi

markdown-fix: _ensure-rumdl
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        case "$file" in
            CHANGELOG.md|docs/archive/*) continue ;;
        esac
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.md')
    if [ "${#files[@]}" -gt 0 ]; then
        echo "📝 rumdl check --fix (${#files[@]} files)"
        printf '%s\0' "${files[@]}" | xargs -0 -n100 rumdl check --fix
    else
        echo "No markdown files found to format."
    fi

markdown-lint: markdown-check

# Backward-compatible alias for the GitHub Actions release-asset comparison.
performance-archive-published current_tag="" baseline_tag="":
    just performance-github-assets "{{ current_tag }}" "{{ baseline_tag }}"

# Compare stored GitHub Actions release benchmark assets without local cargo runs.
performance-github-assets current_tag="" baseline_tag="": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    current_tag="{{ current_tag }}"
    baseline_tag="{{ baseline_tag }}"
    if [[ -n "$current_tag" || -n "$baseline_tag" ]]; then
        if [[ -z "$current_tag" || -z "$baseline_tag" ]]; then
            echo "current_tag and baseline_tag must be provided together" >&2
            exit 2
        fi
        uv run --locked archive-performance "$current_tag" "$baseline_tag" --github-assets --generate-in-temp-worktree --worktree-ref "$current_tag" --output-only --output target/bench-reports/github-assets-performance.md
    else
        uv run --locked archive-performance --published-latest --github-assets --generate-in-temp-worktree --output-only --output target/bench-reports/github-assets-performance.md
    fi

# Compare the current tree against the latest published release locally.
performance-local: python-sync
    uv run --locked archive-performance --current-vs-latest --generate-in-temp-worktree --output-only --output target/bench-reports/performance.md

# Compare current non-exact kernels locally without rerunning current peer crates.
performance-local-vs-linalg current_tag="" baseline_tag="": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    current_tag="{{ current_tag }}"
    baseline_tag="{{ baseline_tag }}"
    if [[ -n "$current_tag" || -n "$baseline_tag" ]]; then
        if [[ -z "$current_tag" || -z "$baseline_tag" ]]; then
            echo "current_tag and baseline_tag must be provided together" >&2
            exit 2
        fi
        uv run --locked archive-performance "$current_tag" "$baseline_tag" --suite vs_linalg --generate-in-temp-worktree --worktree-ref HEAD --output-only --output target/bench-reports/performance.md
    else
        uv run --locked archive-performance --current-vs-latest --suite vs_linalg --generate-in-temp-worktree --output-only --output target/bench-reports/performance.md
    fi

# Generate local release-signal measurements in a temp worktree, then promote/archive docs.
performance-release current_tag="" baseline_tag="": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    current_tag="{{ current_tag }}"
    baseline_tag="{{ baseline_tag }}"
    if [[ -n "$current_tag" || -n "$baseline_tag" ]]; then
        if [[ -z "$current_tag" || -z "$baseline_tag" ]]; then
            echo "current_tag and baseline_tag must be provided together" >&2
            exit 2
        fi
        uv run --locked archive-performance "$current_tag" "$baseline_tag" --generate-in-temp-worktree --worktree-ref HEAD
    else
        uv run --locked archive-performance --infer-release --generate-in-temp-worktree --worktree-ref HEAD
    fi

# Plot: generate a single time-vs-dimension SVG from Criterion results.
plot-vs-linalg metric="lu_solve" stat="median" sample="new" log_y="false" allow_partial="false": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--metric "{{ metric }}" --stat "{{ stat }}" --sample "{{ sample }}")
    if [ "{{ log_y }}" = "true" ]; then
        args+=(--log-y)
    fi
    if [ "{{ allow_partial }}" = "true" ]; then
        args+=(--allow-partial)
    fi
    uv run --locked criterion-dim-plot "${args[@]}"

# Validate fixtures, rerun the full comparison, and atomically publish the canonical README assets/table.
plot-vs-linalg-readme metric="lu_solve" stat="median" sample="new" log_y="true": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--metric "{{ metric }}" --stat "{{ stat }}" --sample "{{ sample }}" --update-readme)
    if [ "{{ log_y }}" = "true" ]; then
        args+=(--log-y)
    fi
    uv run --locked criterion-dim-plot "${args[@]}"

# Python tooling (uv)
python-check: python-typecheck
    uv run --locked ruff format --check scripts/
    uv run --locked ruff check scripts/

python-fix: python-sync
    uv run --locked ruff check scripts/ --fix
    uv run --locked ruff format scripts/

python-lint: python-check

python-sync: _ensure-uv
    uv sync --locked --group dev

python-typecheck: python-sync
    uv run --locked ty check scripts/ --error all

# Repository-owned Semgrep rules for project-specific diagnostics.
semgrep: _ensure-uv
    uv run --locked semgrep --metrics off --error --strict --timeout 30 --config semgrep.yaml .
    uv run --locked check-docs-version-sync

# Fixture tests for repository-owned Semgrep rules.
semgrep-test: _ensure-uv
    #!/usr/bin/env bash
    set -euo pipefail

    check_semgrep_fixture() {
        target="$1"
        json="$(uv run --locked semgrep scan --metrics off --json --quiet --strict --config semgrep.yaml "$target")"
        SEMGREP_JSON="$json" uv run --locked scripts/check_semgrep_fixtures.py "$target"
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

    echo "Building project..."
    cargo build
    echo "✅ Setup complete! Run 'just help-workflows' to see available commands."

# Development tooling installation and verification
setup-tools:
    #!/usr/bin/env bash
    set -euo pipefail

    have() { command -v "$1" >/dev/null 2>&1; }

    installed_tool_version() {
        case "$1" in
            cargo-llvm-cov)
                cargo llvm-cov --version 2>/dev/null
                ;;
            cargo-machete)
                cargo machete --version 2>/dev/null
                ;;
            cargo-nextest)
                cargo nextest --version 2>/dev/null
                ;;
            *)
                "$1" --version 2>/dev/null
                ;;
        esac | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true
    }

    verify_tool_version() {
        local cmd="$1"
        local expected="$2"
        local actual=""
        local resolved=""

        actual="$(installed_tool_version "$cmd")"
        resolved="$(command -v "$cmd" 2>/dev/null || true)"
        if [[ "$actual" != "$expected" ]]; then
            echo "❌ '$cmd' resolves to '${resolved:-missing}' at version '${actual:-missing}', expected '$expected'." >&2
            return 1
        fi
        echo "  ✓ $cmd $actual"
    }

    echo "🔧 Ensuring tooling required by just recipes is installed..."
    echo ""
    echo "Ensuring Rust components..."
    if ! have rustup; then
        echo "❌ 'rustup' not found. Install Rust via https://rustup.rs and re-run: just setup-tools"
        exit 1
    fi
    rustup component add clippy rustfmt rust-src llvm-tools-preview
    echo ""

    echo "Ensuring cargo tools..."
    just_version="{{ just_version }}"
    if ! have just || [[ "$(just --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$just_version" ]]; then
        cargo install --locked just --version "$just_version"
    fi

    cargo_llvm_cov_version="{{ cargo_llvm_cov_version }}"
    if ! have cargo-llvm-cov || [[ "$(cargo llvm-cov --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$cargo_llvm_cov_version" ]]; then
        cargo install --locked cargo-llvm-cov --version "$cargo_llvm_cov_version"
    fi
    cargo_machete_version="{{ cargo_machete_version }}"
    if ! cargo machete --version >/dev/null 2>&1 || [[ "$(cargo machete --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$cargo_machete_version" ]]; then
        cargo install --locked cargo-machete --version "$cargo_machete_version"
    fi

    cargo_nextest_version="{{ cargo_nextest_version }}"
    if ! cargo nextest --version >/dev/null 2>&1 || [[ "$(cargo nextest --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$cargo_nextest_version" ]]; then
        cargo install --locked cargo-nextest --version "$cargo_nextest_version"
    fi
    dprint_version="{{ dprint_version }}"
    if ! have dprint || [[ "$(dprint --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$dprint_version" ]]; then
        cargo install --locked dprint --version "$dprint_version"
    fi
    git_cliff_version="{{ git_cliff_version }}"
    if ! have git-cliff || [[ "$(git-cliff --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$git_cliff_version" ]]; then
        cargo install --locked git-cliff --version "$git_cliff_version"
    fi
    rumdl_version="{{ rumdl_version }}"
    if ! have rumdl || [[ "$(rumdl --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$rumdl_version" ]]; then
        cargo install --locked rumdl --version "$rumdl_version"
    fi
    taplo_version="{{ taplo_version }}"
    if ! have taplo || [[ "$(taplo --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$taplo_version" ]]; then
        cargo install --locked taplo-cli --version "$taplo_version"
    fi
    typos_version="{{ typos_version }}"
    if ! have typos || [[ "$(typos --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$typos_version" ]]; then
        cargo install --locked typos-cli --version "$typos_version"
    fi
    zizmor_version="{{ zizmor_version }}"
    if ! have zizmor || [[ "$(zizmor --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)" != "$zizmor_version" ]]; then
        cargo install --locked zizmor --version "$zizmor_version"
    fi
    echo ""

    uv_version="{{ uv_version }}"
    if have uv; then
        echo "Ensuring uv-managed Python tools..."
        uv sync --locked --group dev
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
    echo "Verifying required commands and versions..."
    have jq || { echo "❌ 'jq' is still missing."; exit 1; }
    echo "  ✓ jq"
    verify_tool_version just "$just_version"
    verify_tool_version cargo-llvm-cov "$cargo_llvm_cov_version"
    verify_tool_version cargo-machete "$cargo_machete_version"
    verify_tool_version cargo-nextest "$cargo_nextest_version"
    verify_tool_version dprint "$dprint_version"
    verify_tool_version git-cliff "$git_cliff_version"
    verify_tool_version rumdl "$rumdl_version"
    verify_tool_version taplo "$taplo_version"
    verify_tool_version typos "$typos_version"
    verify_tool_version uv "$uv_version"
    verify_tool_version zizmor "$zizmor_version"
    uv run --locked actionlint -version >/dev/null
    echo "  ✓ actionlint (uv)"
    for cmd in pytest ruff semgrep shellcheck shfmt ty yamllint; do
        uv run --locked "$cmd" --version >/dev/null
        echo "  ✓ $cmd (uv)"
    done

    echo ""
    echo "✅ Tooling setup complete."

# Shell scripts
shell-check: _ensure-shellcheck _ensure-shfmt
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.sh')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 -n4 uv run --locked shellcheck -x
        printf '%s\0' "${files[@]}" | xargs -0 uv run --locked shfmt -d
    else
        echo "No shell files found to check."
    fi

shell-fix: shell-fmt

shell-fmt: _ensure-shfmt
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.sh')
    if [ "${#files[@]}" -gt 0 ]; then
        echo "🧹 shfmt -w (${#files[@]} files)"
        printf '%s\0' "${files[@]}" | xargs -0 -n1 uv run --locked shfmt -w
    else
        echo "No shell files found to format."
    fi

shell-lint: shell-check

# Spell check (typos)
spell-check: _ensure-typos
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    # Check every tracked file plus untracked, non-ignored additions. This keeps
    # clean CI checkouts covered while validating new files before they are staged.
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z)
    if [ "${#files[@]}" -gt 0 ]; then
        # Exclude typos.toml itself: it intentionally contains allowlisted fragments.
        printf '%s\0' "${files[@]}" | xargs -0 -n100 typos --config typos.toml --force-exclude --exclude typos.toml --
    else
        echo "No files found to spell-check."
    fi

# Create an annotated git tag from the CHANGELOG.md section for the given version
tag version: python-sync
    uv run --locked tag-release {{ version }}

# Recreate an existing tag (delete + recreate)
tag-force version: python-sync
    uv run --locked tag-release {{ version }} --force

# Testing: runnable Rust tests use nextest; rustdoc doctests remain on cargo test.
test: test-lib test-doc

test-all: test-rust test-python
    @echo "✅ All tests passed"

test-doc:
    cargo test --doc --verbose

test-doc-exact:
    cargo test --features exact --doc --verbose

# Tests for the "exact" feature (exact determinants, conversions, and Bareiss solves)
test-exact: _ensure-cargo-nextest test-doc-exact
    cargo nextest run --profile ci --features exact --verbose

# Smoke-test deterministic inputs and configuration shared with benchmark suites.
test-bench-inputs: _ensure-cargo-nextest
    cargo nextest run --profile ci --features bench,exact --test vs_linalg_inputs --test exact_bench_config --verbose

# Compile all integration-test targets without running them.
test-integration-compile: _ensure-cargo-nextest
    cargo nextest run --all-features --tests --no-run

test-integration: _ensure-cargo-nextest
    cargo nextest run --profile ci --tests --verbose

test-lib: _ensure-cargo-nextest
    cargo nextest run --profile ci --lib --verbose

# CI Rust bucket: all runnable unit/integration targets in one nextest pass.
test-rust-ci: _ensure-cargo-nextest
    cargo nextest run --profile ci --all-features --lib --tests --verbose

test-rust: test-rust-ci test-doc test-doc-exact
    @echo "✅ Rust tests passed"

test-unit: test-lib

test-python: python-sync
    uv run --locked pytest -q

# TOML
toml-check: toml-fmt-check toml-lint

toml-fix: toml-fmt

toml-fmt: _ensure-taplo
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.toml')
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
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.toml')
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
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.toml')
    if [ "${#files[@]}" -gt 0 ]; then
        taplo lint "${files[@]}"
    else
        echo "No TOML files found to lint."
    fi

# Check for unused direct Cargo dependencies.
unused-deps: _ensure-cargo-machete
    cargo machete

# File validation
json-check: validate-json

validate-json: _ensure-jq
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.json')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 -n1 jq empty
    else
        echo "No JSON files found to validate."
    fi

# YAML
yaml-check: yaml-fmt-check yaml-lint

yaml-fmt-check: _ensure-dprint
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.yml' '*.yaml' 'CITATION.cff')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 dprint check --incremental=false
    else
        echo "No YAML files found to check."
    fi

yaml-fix: _ensure-dprint
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.yml' '*.yaml' 'CITATION.cff')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 dprint fmt --incremental=false
    else
        echo "No YAML files found to format."
    fi

yaml-lint: _ensure-yamllint
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.yml' '*.yaml' 'CITATION.cff')
    if [ "${#files[@]}" -gt 0 ]; then
        echo "🔍 yamllint (${#files[@]} YAML/CFF files)"
        uv run --locked yamllint --strict -c .yamllint "${files[@]}"
    else
        echo "No YAML files found to lint."
    fi

# Keep the command-memory layer itself canonically formatted.
justfile-fmt-check:
    just --fmt --check

# GitHub Actions security analysis
zizmor: _ensure-zizmor
    zizmor .github
