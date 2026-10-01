# shellcheck disable=SC2148
# Justfile for la-stack development workflow
# Install tools: uv run --locked --managed-python --only-group tooling research-repo-tools setup
# Usage: just <command> or just --list

# Use bash with strict error handling for all recipes
set shell := ["bash", "-euo", "pipefail", "-c"]

_run := "uv run --locked --no-sync --no-python-downloads research-repo-tools toolchain run --"

# Coverage (cargo-llvm-cov)
#
# Common cargo-llvm-cov arguments for all coverage runs.
_coverage_base_args := '''--features exact \
  --workspace --lib --tests \
  --verbose'''
# System prerequisites remain consumer-owned.
_ensure-gh:
    @command -v gh >/dev/null || { echo "GitHub CLI is required on PATH." >&2; exit 1; }

_ensure-jq:
    @command -v jq >/dev/null || { echo "jq is required on PATH." >&2; exit 1; }

# uv enforces its exact declaration before running the locked CLI.
_ensure-uv:
    uv run --locked --no-sync --no-python-downloads research-repo-tools deps check-uv

_ensure-actionlint: _ensure-uv
    uv run --locked actionlint -version >/dev/null

_ensure-shellcheck: _ensure-uv
    uv run --locked shellcheck --version >/dev/null

_ensure-shfmt: _ensure-uv
    uv run --locked shfmt --version >/dev/null

_ensure-yamllint: _ensure-uv
    uv run --locked yamllint --version >/dev/null

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
    {{ _run }} cargo bench --locked --workspace --features bench

# Compare latest measurements against a saved baseline.
# Defaults to the `last` full-release baseline.
bench-compare baseline="last" suite="all" scope="release-signal": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    baseline={{ quote(baseline) }}
    {{ _run }} uv run --locked bench-compare "$baseline" --suite {{ quote(suite) }} --scope {{ quote(scope) }}

# Compile benchmarks without running them, treating warnings as errors through
# Cargo so warning policy does not create separate rustc cache artifacts.
# This catches bench/release-profile-only warnings that won't show up in normal debug-profile runs.
bench-compile:
    CARGO_BUILD_WARNINGS=deny {{ _run }} cargo bench --locked --workspace --no-run --features bench
    CARGO_BUILD_WARNINGS=deny {{ _run }} cargo bench --locked --no-run --features bench,exact --bench exact

# Run the exact-arithmetic benchmark suite.
bench-exact:
    {{ _run }} cargo bench --locked --features bench,exact --bench exact

# Run the outward-rounded interval determinant benchmark suite.
bench-interval:
    {{ _run }} cargo bench --locked --features bench --bench interval

# Run the certified dot-product and affine-difference benchmark suite.
bench-linear-form:
    {{ _run }} cargo bench --locked --features bench --bench linear_form

# Run the cheaper latest measurements used for latest-vs-last reports.
bench-latest: bench-vs-linalg-la-stack bench-exact

# Run latest measurements and render the latest-vs-last performance report.
bench-latest-vs-last baseline="last": bench-latest python-sync
    {{ _run }} uv run --locked bench-compare {{ quote(baseline) }}

# Discover all release benchmarks and report their expected measurement budget.
bench-release-inventory: _ensure-uv
    {{ _run }} uv run --locked scripts/release_baseline.py inventory

# Check every discovered benchmark before the release workflow packages it.
bench-release-check tag: _ensure-uv
    {{ _run }} uv run --locked scripts/release_baseline.py validate --baseline {{ quote(tag) }}

# Save a Criterion baseline. Defaults to all release-signal benchmark suites.
bench-save-baseline tag suite="all":
    #!/usr/bin/env bash
    set -euo pipefail
    suite={{ quote(suite) }}
    case "$suite" in
        all)
            {{ _run }} cargo bench --locked -p la-stack-comparison --features bench --bench vs_linalg -- --noplot --save-baseline {{ quote(tag) }}
            {{ _run }} cargo bench --locked --features bench,exact --bench exact -- --noplot --save-baseline {{ quote(tag) }}
            ;;
        exact)
            {{ _run }} cargo bench --locked --features bench,exact --bench exact -- --noplot --save-baseline {{ quote(tag) }}
            ;;
        vs_linalg)
            {{ _run }} cargo bench --locked -p la-stack-comparison --features bench --bench vs_linalg -- --noplot --save-baseline {{ quote(tag) }}
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
    filter={{ quote(filter) }}
    if [ -n "$filter" ]; then
        {{ _run }} cargo bench --locked -p la-stack-comparison --features bench --bench vs_linalg -- "$filter"
    else
        {{ _run }} cargo bench --locked -p la-stack-comparison --features bench --bench vs_linalg
    fi

# Bench only la-stack rows from the vs_linalg suite for cheap latest-vs-last comparisons.
# Filtered runs omit peer samples, so disable Criterion's complete-group HTML reports.
bench-vs-linalg-la-stack:
    {{ _run }} cargo bench --locked -p la-stack-comparison --features bench --bench vs_linalg -- la_stack --noplot

# Run only la-stack vs_linalg measurements and render a non-exact performance report.
bench-vs-linalg-latest-vs baseline="last": bench-vs-linalg-la-stack python-sync
    {{ _run }} uv run --locked bench-compare {{ quote(baseline) }} --suite vs_linalg --scope release-signal

# Quick iteration (reduced runtime, no Criterion HTML).
bench-vs-linalg-quick filter="":
    #!/usr/bin/env bash
    set -euo pipefail
    filter={{ quote(filter) }}
    if [ -n "$filter" ]; then
        {{ _run }} cargo bench --locked -p la-stack-comparison --features bench --bench vs_linalg -- "$filter" --quick --noplot
    else
        {{ _run }} cargo bench --locked -p la-stack-comparison --features bench --bench vs_linalg -- --quick --noplot
    fi

# Build commands
build:
    {{ _run }} cargo build

build-release:
    {{ _run }} cargo build --release

# Verify Cargo.toml and the committed Cargo.lock are synchronized.
cargo-lock-check:
    {{ _run }} cargo metadata --locked --format-version 1 --no-deps > /dev/null

# Generate, normalize, and rotate completed minor series with the pinned shared CLI.
changelog: tools-check python-sync
    {{ _run }} research-repo-tools changelog generate

# Rotate existing history without regenerating release notes.
changelog-archive: python-sync
    uv run --locked --group dev research-repo-tools changelog archive

# Check release headings and all archives without writing files.
changelog-check: tools-check python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob
    uv run --locked --group dev research-repo-tools changelog check
    {{ _run }} rumdl check --no-cache --config pyproject.toml CHANGELOG.md docs/archives/changelog/*.md

# Validate generation and print the root changelog without publishing candidates.
[positional-arguments]
changelog-preview *args: tools-check python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    {{ _run }} research-repo-tools changelog generate --dry-run "$@"

# Generate a prospective release using the explicit ISO date, without updating metadata.
changelog-release tag date: tools-check python-sync
    {{ _run }} research-repo-tools changelog generate --tag {{ quote(tag) }} --date {{ quote(date) }}

alias changelog-unreleased := changelog-release

# Check (non-mutating): run all linters/validators
check: lint
    @echo "✅ Checks complete!"

# Fast compile check (no binary produced)
check-fast:
    {{ _run }} cargo check

# CI simulation: flat GitHub-equivalent union of leaf validators.
# Keep this dependency list explicit so each validation surface runs once without
# re-entering broad check/test bundles. All Cargo targets match the SARIF lint scope.
ci: action-lint zizmor markdown-check spell-check docs-version-check changelog-check toml-parse-check toml-fmt-check toml-lint yaml-fmt-check yaml-lint citation-check validate-json justfile-fmt-check python-format-check python-lint python-fixture-lint python-typecheck test-python cargo-lock-check fmt-check clippy-all-targets doc-check semgrep semgrep-test unused-deps shell-check test-rust-ci test-doc test-doc-exact bench-compile examples
    @echo "🎯 CI checks complete!"

# Validate CITATION.cff against the Citation File Format schema.
citation-check: _ensure-uv
    uvx --from cffconvert==2.0.0 cffconvert --validate -i CITATION.cff

# Clean build artifacts
clean:
    {{ _run }} cargo clean
    rm -rf target/llvm-cov
    rm -rf coverage

# Full Cargo-target Clippy sweep used by `just ci` and the GitHub SARIF workflow.
clippy: clippy-all-targets

clippy-all-targets:
    {{ _run }} cargo clippy --workspace --all-targets
    {{ _run }} cargo clippy --workspace --all-targets --all-features

# Core library Clippy checks used by the orthogonal CI graph.
clippy-core:
    {{ _run }} cargo clippy --workspace --lib
    {{ _run }} cargo clippy --workspace --lib --all-features

# Clippy for the "exact" feature (catches feature-gated lint issues)
clippy-exact:
    {{ _run }} cargo clippy --features exact --all-targets

# Coverage analysis for local development (HTML output)
coverage: tools-check
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p target/llvm-cov
    {{ _run }} cargo llvm-cov nextest {{ _coverage_base_args }} --open --output-dir target/llvm-cov -P coverage
    echo "Coverage report generated: target/llvm-cov/html/index.html"

# Coverage analysis for CI (XML output for Codecov)
coverage-ci: tools-check
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p coverage
    {{ _run }} cargo llvm-cov nextest {{ _coverage_base_args }} --cobertura --output-path coverage/cobertura.xml -P coverage

# Default recipe shows available commands
default:
    @just --list

# Documentation build checks for the default and exact-feature public APIs.
doc-check:
    RUSTDOCFLAGS='-D warnings' {{ _run }} cargo doc --no-deps
    RUSTDOCFLAGS='-D warnings' {{ _run }} cargo doc --no-deps --features exact

docs-version-check: python-sync
    uv run --locked --group dev research-repo-tools release check

# Examples
examples:
    #!/usr/bin/env bash
    set -euo pipefail
    {{ _run }} cargo build --features exact --examples

    exe_suffix=""
    if [[ "${OS:-}" == "Windows_NT" ]]; then
        exe_suffix=".exe"
    fi

    target_dir="${CARGO_TARGET_DIR:-target}"

    shopt -s nullglob
    for example_path in examples/*.rs; do
        [[ -f "${example_path}" ]] || continue
        example="${example_path##*/}"
        example="${example%.rs}"
        "${target_dir}/debug/examples/${example}${exe_suffix}"
    done

# Fix (mutating): apply formatters/auto-fixes
fix: toml-fmt fmt python-fix shell-fmt markdown-fix yaml-fix
    @echo "✅ Fixes applied!"

# Rust formatting
fmt:
    {{ _run }} cargo fmt --all

fmt-check:
    {{ _run }} cargo fmt --all -- --check

github-actions-check: action-lint zizmor
    @echo "✅ GitHub Actions checks complete!"

help-workflows:
    @echo "Common Just workflows:"
    @echo "  just check             # Run lint/validators (non-mutating)"
    @echo "  just check-fast        # Fast compile check (cargo check)"
    @echo "  just ci                # Full CI simulation (check + tests + examples + bench compile)"
    @echo "  just fix               # Apply formatters/auto-fixes (mutating)"
    @echo "  just setup             # Install/verify dev tools + sync Python deps"
    @echo ""
    @echo "CodeRabbit review (opt-in):"
    @echo "  just review [base]     # Review branch and local changes; verify origin/main by default"
    @echo "  just review-uncommitted # Review only local changes, including untracked files"
    @echo ""
    @echo "Benchmarks:"
    @echo "  just bench                 # Run benchmarks"
    @echo "  just bench-compile          # Compile benches with warnings-as-errors"
    @echo "  just bench-latest           # Run cheap latest measurements"
    @echo "  just bench-latest-vs-last   # Run latest and compare against last"
    @echo "  just bench-exact            # Run exact-arithmetic benchmarks"
    @echo "  just bench-interval         # Run interval determinant benchmarks"
    @echo "  just bench-linear-form      # Run certified linear-form benchmarks"
    @echo "  just bench-save-last        # Save full baseline as 'last'"
    @echo "  just bench-vs-linalg        # Run vs_linalg bench (optional filter)"
    @echo "  just bench-vs-linalg-la-stack # Run la-stack rows from vs_linalg"
    @echo "  just bench-vs-linalg-latest-vs # Run non-exact latest and compare against last"
    @echo "  just bench-vs-linalg-quick  # Quick vs_linalg bench (reduced samples)"
    @echo "  just performance-doc        # Build release docs from retained CSV/JSON"
    @echo "  just performance-github-assets # Compare stored GitHub Actions release assets"
    @echo "  just performance-local      # Compare current tree against latest release locally"
    @echo "  just performance-local-non-exact # Compare current non-exact kernels locally"
    @echo "  just performance-readme     # Publish retained release data to README assets/table"
    @echo "  just performance-release    # Measure, retain, and publish release docs"
    @echo ""
    @echo "Benchmark plotting:"
    @echo "  just plot-vs-linalg         # Plot Criterion results (CSV + SVG + provenance)"
    @echo ""
    @echo "Changelog & releases:"
    @echo "  just changelog              # Regenerate CHANGELOG.md from full history"
    @echo "  just changelog-preview       # Preview generation without writing"
    @echo "  just changelog-release <tag> <date>  # Generate a release with an explicit date"
    @echo "  just changelog-unreleased <tag> <date>  # Alias for changelog-release"
    @echo "  just changelog-archive       # Rotate existing changelog history"
    @echo "  just changelog-check         # Validate root and archived release notes"
    @echo "  just release-notes <tag>     # Print root or archived release notes"
    @echo "  just tag <ver>              # Create annotated tag from CHANGELOG.md"
    @echo "  just tag-force <ver>        # Recreate an existing tag"
    @echo "  just update-version <tag>       # Update release metadata and infer the previous tag"
    @echo ""
    @echo "Setup:"
    @echo "  just setup             # Install declared tools, sync dev, and build"
    @echo "  just setup-tools       # Install/verify external tooling"
    @echo "  just update            # Update dependencies and repository-owned tool pins"
    @echo ""
    @echo "Testing:"
    @echo "  just coverage          # Generate coverage report (HTML)"
    @echo "  just coverage-ci       # Generate coverage for CI (XML)"
    @echo "  just examples          # Run examples"
    @echo "  just test              # Lib + doc tests (fast)"
    @echo "  just test-all          # All tests (Rust + Python)"
    @echo "  just test-bench-inputs # Benchmark input smoke tests"
    @echo "  just test-doc          # Default-feature doctests"
    @echo "  just test-exact        # Exact-feature tests and doctests"
    @echo "  just test-integration  # Integration tests"
    @echo "  just test-python       # Python tests only (pytest)"
    @echo "  just test-rust-ci      # Release-profile unit + integration CI bucket"
    @echo "  just test-unit         # Library unit tests"
    @echo ""
    @echo "Note: Some recipes require external tools. Run 'just setup-tools' (tooling) or 'just setup' (full env) first."

# File validation
json-check: validate-json

# Keep the command-memory layer itself canonically formatted.
justfile-fmt-check:
    just --fmt --check

# Lint groups (delaunay-style)
lint: lint-code lint-docs lint-config

lint-code: rust-core-check python-check shell-check

lint-config: json-check toml-ci yaml-ci github-actions-check justfile-fmt-check

lint-docs: markdown-ci docs-version-check changelog-check

# Markdown
markdown-check: tools-check _ensure-uv
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        case "$file" in
            CHANGELOG.md|docs/archive/*|docs/archives/changelog/*) continue ;;
        esac
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.md')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 -n100 {{ _run }} rumdl check
        uv run --locked --group dev research-repo-tools docs check-lines "${files[@]}"
    else
        echo "No markdown files found to check."
    fi

markdown-ci: markdown-check spell-check
    @echo "✅ Markdown checks complete!"

markdown-fix: tools-check
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        case "$file" in
            CHANGELOG.md|docs/archive/*|docs/archives/changelog/*) continue ;;
        esac
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.md')
    if [ "${#files[@]}" -gt 0 ]; then
        echo "📝 rumdl check --fix (${#files[@]} files)"
        printf '%s\0' "${files[@]}" | xargs -0 -n100 {{ _run }} rumdl check --fix
    else
        echo "No markdown files found to format."
    fi

markdown-lint: markdown-check

# Build release docs from retained scratch inputs or the latest docs/performance snapshot.
performance-doc: python-sync
    {{ _run }} uv run --locked archive-performance --promote-artifacts

# Compare stored GitHub Actions release benchmark assets without local cargo runs.
performance-github-assets current_tag="" baseline_tag="": _ensure-gh python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    current_tag={{ quote(current_tag) }}
    baseline_tag={{ quote(baseline_tag) }}
    if [[ -n "$current_tag" || -n "$baseline_tag" ]]; then
        if [[ -z "$current_tag" || -z "$baseline_tag" ]]; then
            echo "current_tag and baseline_tag must be provided together" >&2
            exit 2
        fi
        {{ _run }} uv run --locked archive-performance "$current_tag" "$baseline_tag" --github-assets --generate-in-temp-worktree --worktree-ref "$current_tag" --output-only --output target/bench-reports/github-assets-performance.md --artifact-csv target/bench-reports/github-assets-performance.csv --artifact-provenance target/bench-reports/github-assets-performance.provenance.json
    else
        {{ _run }} uv run --locked archive-performance --published-latest --github-assets --generate-in-temp-worktree --output-only --output target/bench-reports/github-assets-performance.md --artifact-csv target/bench-reports/github-assets-performance.csv --artifact-provenance target/bench-reports/github-assets-performance.provenance.json
    fi

# Compare the current tree against the latest release; untracked files are excluded.
performance-local: _ensure-gh python-sync
    {{ _run }} uv run --locked archive-performance --current-vs-latest --generate-in-temp-worktree --output-only --local-report --output target/bench-reports/performance.md

# Compare current non-exact kernels locally without rerunning current peer crates.
performance-local-non-exact current_tag="" baseline_tag="": _ensure-gh python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    current_tag={{ quote(current_tag) }}
    baseline_tag={{ quote(baseline_tag) }}
    if [[ -n "$current_tag" || -n "$baseline_tag" ]]; then
        if [[ -z "$current_tag" || -z "$baseline_tag" ]]; then
            echo "current_tag and baseline_tag must be provided together" >&2
            exit 2
        fi
        {{ _run }} uv run --locked archive-performance "$current_tag" "$baseline_tag" --suite vs_linalg --generate-in-temp-worktree --worktree-ref HEAD --output-only --local-report --output target/bench-reports/performance-non-exact.md --artifact-csv target/bench-reports/performance-non-exact.csv --artifact-provenance target/bench-reports/performance-non-exact.provenance.json
    else
        {{ _run }} uv run --locked archive-performance --current-vs-latest --suite vs_linalg --generate-in-temp-worktree --output-only --local-report --output target/bench-reports/performance-non-exact.md --artifact-csv target/bench-reports/performance-non-exact.csv --artifact-provenance target/bench-reports/performance-non-exact.provenance.json
    fi

# Publish README assets/table from retained measurements, including after target cleanup.
performance-readme metric="lu_solve" stat="median" sample="new" log_y="true": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--metric {{ quote(metric) }} --stat {{ quote(stat) }} --sample {{ quote(sample) }} --update-readme)
    if [ {{ quote(log_y) }} = "true" ]; then
        args+=(--log-y)
    fi
    {{ _run }} uv run --locked criterion-dim-plot "${args[@]}"

# Measure locally, preserve complete summaries in docs/performance, and promote/archive docs.
performance-release current_tag="" baseline_tag="": _ensure-gh python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    current_tag={{ quote(current_tag) }}
    baseline_tag={{ quote(baseline_tag) }}
    if [[ -n "$current_tag" || -n "$baseline_tag" ]]; then
        if [[ -z "$current_tag" || -z "$baseline_tag" ]]; then
            echo "current_tag and baseline_tag must be provided together" >&2
            exit 2
        fi
        {{ _run }} uv run --locked archive-performance "$current_tag" "$baseline_tag" --generate-in-temp-worktree --worktree-ref HEAD
    else
        {{ _run }} uv run --locked archive-performance --infer-release --generate-in-temp-worktree --worktree-ref HEAD
    fi

# Plot: generate a single time-vs-dimension SVG from Criterion results.
plot-vs-linalg metric="lu_solve" stat="median" sample="new" log_y="false" allow_partial="false": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--metric {{ quote(metric) }} --stat {{ quote(stat) }} --sample {{ quote(sample) }})
    if [ {{ quote(log_y) }} = "true" ]; then
        args+=(--log-y)
    fi
    if [ {{ quote(allow_partial) }} = "true" ]; then
        args+=(--allow-partial)
    fi
    {{ _run }} uv run --locked criterion-dim-plot "${args[@]}"

# Python tooling (uv)
python-check: python-format-check python-lint python-fixture-lint python-typecheck

python-ci: python-format-check python-lint python-fixture-lint python-typecheck test-python
    @echo "✅ Python checks complete!"

python-fix: python-sync
    uv run --locked ruff check scripts/ --fix
    uv run --locked ruff format scripts/ tests/semgrep/scripts/

python-format-check: python-sync
    uv run --locked ruff format --check scripts/ tests/semgrep/scripts/

python-lint: python-sync
    uv run --locked ruff check scripts/

python-fixture-lint: python-sync
    uv run --locked ruff check tests/semgrep/scripts/

python-sync: _ensure-uv
    uv sync --locked --group dev

python-typecheck: python-sync
    uv run --locked ty check scripts/ tests/semgrep/scripts/ --error all

# Print release notes from the root changelog or a completed minor archive.
release-notes tag: python-sync
    uv run --locked --group dev research-repo-tools changelog notes {{ quote(tag) }}

# Review branch and local changes against a verified origin/main, or an explicit local base.
review base="origin/main":
    uv run --locked --group dev research-repo-tools review branch --base={{ quote(base) }}

# Review staged, unstaged, and non-ignored untracked changes without a remote lookup.
review-uncommitted:
    uv run --locked --group dev research-repo-tools review uncommitted

# Run the shared dependency and full-history secret scans.
security: security-osv security-secrets

# Audit the repository's Python and Rust lockfiles with the managed OSV scanner.
security-osv:
    uv run --locked --group dev research-repo-tools security osv uv.lock Cargo.lock

# Scan reachable Git history and current tracked/nonignored files with redacted reports.
security-secrets:
    uv run --locked --group dev research-repo-tools security secrets

rust-core-check: cargo-lock-check fmt-check clippy-core doc-check semgrep semgrep-test unused-deps
    @echo "✅ Rust core checks complete!"

# Repository-owned Semgrep rules for project-specific diagnostics.
semgrep: _ensure-uv
    uv run --locked semgrep --metrics off --error --strict --timeout 30 --exclude tests/semgrep/src/project_rules/algebraic_float.rs --config semgrep.yaml .

# Fixture tests for repository-owned Semgrep rules.
semgrep-test: python-sync
    uv run --locked --group dev research-repo-tools semgrep check-fixtures

# Install declared managed tools and explicitly synchronize dev.
setup: _ensure-gh _ensure-jq
    uv run --locked --managed-python --only-group tooling research-repo-tools setup
    {{ _run }} cargo build

# Tool-only setup; compilation belongs to setup.
setup-tools: _ensure-gh _ensure-jq
    uv run --locked --managed-python --only-group tooling research-repo-tools setup

# Check installed tools without synchronization, downloads, or installation.
tools-check:
    uv run --locked --no-sync --no-python-downloads research-repo-tools toolchain check

# Export verified managed paths to GITHUB_ENV in hosted workflows.
tools-export:
    uv run --locked --no-sync --no-python-downloads research-repo-tools toolchain export

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
spell-check: tools-check
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
        printf '%s\0' "${files[@]}" | xargs -0 -n100 {{ _run }} typos --config typos.toml --force-exclude --exclude typos.toml --
    else
        echo "No files found to spell-check."
    fi

# Create an annotated git tag from the CHANGELOG.md section for the given version
tag version: python-sync
    uv run --locked --group dev research-repo-tools changelog tag {{ quote(version) }}

# Replace an existing local tag only after validating the release notes.
tag-force version: python-sync
    uv run --locked --group dev research-repo-tools changelog tag {{ quote(version) }} --force

# Testing: runnable Rust tests use nextest; rustdoc doctests remain on cargo test.
test: test-lib test-doc

test-all: test-rust test-python
    @echo "✅ All tests passed"

# Smoke-test deterministic inputs and configuration shared with benchmark suites.
test-bench-inputs: tools-check
    {{ _run }} cargo nextest run --workspace --profile ci --features bench,exact --test vs_linalg_inputs --test exact_bench_config --verbose

test-doc:
    {{ _run }} cargo test --doc --verbose

test-doc-exact:
    {{ _run }} cargo test --features exact --doc --verbose

# Tests for the "exact" feature (exact determinants, conversions, and Bareiss solves)
test-exact: tools-check test-doc-exact
    {{ _run }} cargo nextest run --profile ci --features exact --verbose

test-integration: tools-check
    {{ _run }} cargo nextest run --profile ci --tests --verbose

# Compile all integration-test targets without running them.
test-integration-compile: tools-check
    {{ _run }} cargo nextest run --all-features --tests --no-run

test-lib: tools-check
    {{ _run }} cargo nextest run --profile ci --lib --verbose

test-python: tools-check python-sync
    {{ _run }} uv run --locked pytest -q

test-rust: test-rust-ci test-doc test-doc-exact
    @echo "✅ Rust tests passed"

# CI Rust bucket: all runnable unit/integration targets in one nextest pass.
test-rust-ci: tools-check
    {{ _run }} cargo nextest run --workspace --release --profile ci --all-features --lib --tests --verbose

test-unit: test-lib

# TOML
toml-check: toml-parse-check toml-fmt-check toml-lint

toml-ci: toml-check
    @echo "✅ TOML checks complete!"

toml-fix: toml-fmt

toml-fmt: tools-check
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.toml')
    if [ "${#files[@]}" -gt 0 ]; then
        {{ _run }} taplo fmt "${files[@]}"
    else
        echo "No TOML files found to format."
    fi

toml-fmt-check: tools-check
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.toml')
    if [ "${#files[@]}" -gt 0 ]; then
        {{ _run }} taplo fmt --check "${files[@]}"
    else
        echo "No TOML files found to check."
    fi

toml-lint: tools-check
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.toml')
    if [ "${#files[@]}" -gt 0 ]; then
        {{ _run }} taplo lint "${files[@]}"
    else
        echo "No TOML files found to lint."
    fi

# Validate TOML syntax independently with Python's standard-library parser.
toml-parse-check: python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.toml')
    if [ "${#files[@]}" -gt 0 ]; then
        uv run --locked python -c 'import pathlib, sys, tomllib; [tomllib.loads(pathlib.Path(path).read_text(encoding="utf-8")) for path in sys.argv[1:]]' "${files[@]}"
    else
        echo "No TOML files found to parse."
    fi

# Check for unused direct Cargo dependencies.
unused-deps: tools-check
    {{ _run }} cargo machete

# Upgrade tools before updating dependency requirements and lock resolutions.
update: update-tools update-dependencies
    @echo "✅ Repository dependencies and tools updated."

# Upgrade uv through its owner, then declared Cargo tools, then synchronize setup.
update-tools: update-uv update-cargo-tools setup-tools

# Bootstrap outside the project's old uv-version requirement; do not sync here.
update-uv:
    uv run --no-config --no-sync --no-python-downloads research-repo-tools deps update-uv

# Upgrade only declared managed Cargo tools and publish verified TOML pins.
update-cargo-tools:
    uv run --locked --only-group tooling --inexact research-repo-tools toolchain upgrade

# Dependency-only updates leave uv and managed Cargo tool pins unchanged.
update-dependencies: update-cargo-dependencies update-python-dependencies

# num-bigint and num-rational share public types and must advance together.
update-cargo-dependencies:
    uv run --locked --only-group tooling --inexact research-repo-tools toolchain run -- cargo upgrade --incompatible allow --exclude num-bigint --exclude num-rational
    uv run --locked --only-group tooling --inexact research-repo-tools toolchain run -- cargo update

# Advance direct dev pins, refresh the whole lock, and explicitly synchronize dev.
update-python-dependencies:
    uv run --locked --only-group tooling --inexact research-repo-tools deps update-python
    uv lock --upgrade
    {{ _run }} uv sync --locked --managed-python --group dev

alias update-python-deps := update-python-dependencies

# Update deterministic release metadata, inferring the previous stable published GitHub release.
[doc('Update package, citation, lockfile, and non-artifact documentation release versions.')]
[positional-arguments]
update-version tag *args: _ensure-gh python-sync
    {{ _run }} research-repo-tools release update "$@"

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

yaml-ci: yaml-check citation-check
    @echo "✅ YAML/CFF checks complete!"

yaml-fix: tools-check
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.yml' '*.yaml' 'CITATION.cff')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 {{ _run }} dprint fmt --incremental=false
    else
        echo "No YAML files found to format."
    fi

yaml-fmt-check: tools-check
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        if [ -f "$file" ]; then
            files+=("$file")
        fi
    done < <(git ls-files -co --exclude-standard -z -- '*.yml' '*.yaml' 'CITATION.cff')
    if [ "${#files[@]}" -gt 0 ]; then
        printf '%s\0' "${files[@]}" | xargs -0 {{ _run }} dprint check --incremental=false
    else
        echo "No YAML files found to check."
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

# GitHub Actions security analysis
zizmor: tools-check
    @{{ _run }} bash scripts/run_zizmor.sh
