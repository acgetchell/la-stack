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

# List all public recipes with their arguments and descriptions.
[default]
[private]
_default:
    @just --list

_ensure-actionlint: _ensure-uv
    uv run --locked actionlint -version >/dev/null

# System prerequisites remain consumer-owned.
_ensure-gh:
    @command -v gh >/dev/null || { echo "GitHub CLI is required on PATH." >&2; exit 1; }

_ensure-jq:
    @command -v jq >/dev/null || { echo "jq is required on PATH." >&2; exit 1; }

_ensure-shellcheck: _ensure-uv
    uv run --locked shellcheck --version >/dev/null

_ensure-shfmt: _ensure-uv
    uv run --locked shfmt --version >/dev/null

# uv enforces its exact declaration before running the locked CLI.
_ensure-uv:
    uv run --locked --no-sync --no-python-downloads research-repo-tools deps check-uv

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

# Audit the repository's Python and Rust lockfiles with the managed OSV scanner.
audit:
    uv run --locked --group dev research-repo-tools security osv uv.lock Cargo.lock

# Run the fixed-dimension benchmark suites.
bench:
    {{ _run }} cargo bench --locked --workspace --features bench

# Compare latest measurements against a saved baseline.
# Defaults to the `last` full-release baseline.
bench-compare baseline="last" suite="all" scope="release-signal": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    baseline={{ quote(baseline) }}
    {{ _run }} uv run --locked bench-compare "$baseline" --suite {{ quote(suite) }} --scope {{ quote(scope) }}

# Compile benchmark targets with warnings denied; do not measure timings.
bench-compile:
    CARGO_BUILD_WARNINGS=deny {{ _run }} cargo bench --locked --workspace --no-run --features bench
    CARGO_BUILD_WARNINGS=deny {{ _run }} cargo bench --locked --no-run --features bench,exact --bench exact

# Run the exact-arithmetic benchmark suite.
bench-exact:
    {{ _run }} cargo bench --locked --features bench,exact --bench exact

# Run the outward-rounded interval determinant benchmark suite.
bench-interval:
    {{ _run }} cargo bench --locked --features bench --bench interval

# Run the cheaper latest measurements used for latest-vs-last reports.
bench-latest: bench-vs-linalg-la-stack bench-exact

# Run latest measurements and render the latest-vs-last performance report.
bench-latest-vs-last baseline="last": bench-latest python-sync
    {{ _run }} uv run --locked bench-compare {{ quote(baseline) }}

# Run the certified dot-product and affine-difference benchmark suite.
bench-linear-form:
    {{ _run }} cargo bench --locked --features bench --bench linear_form

# Check every discovered benchmark before the release workflow packages it.
bench-release-check tag: _ensure-uv
    {{ _run }} uv run --locked scripts/release_baseline.py validate --baseline {{ quote(tag) }}

# Discover all release benchmarks and report their expected measurement budget.
bench-release-inventory: _ensure-uv
    {{ _run }} uv run --locked scripts/release_baseline.py inventory

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

# Measure only la-stack rows without generating incomplete peer HTML reports.
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

# Build the library in the debug profile.
build:
    {{ _run }} cargo build

# Build the library in the release profile.
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

# Run the complete local CI validation, tests, examples, and benchmark compilation.
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

# Check every workspace target with default and all features.
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

# Documentation build checks for the default and exact-feature public APIs.
doc-check:
    RUSTDOCFLAGS='-D warnings' {{ _run }} cargo doc --no-deps
    RUSTDOCFLAGS='-D warnings' {{ _run }} cargo doc --no-deps --features exact

# Check synchronized release metadata without changing files.
docs-version-check: python-sync
    uv run --locked --group dev research-repo-tools release check

# Build and run all library examples, including exact arithmetic.
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
fix: toml-fix fmt python-fix shell-fix markdown-fix yaml-fix
    @echo "✅ Fixes applied!"

# Rust formatting
fmt:
    {{ _run }} cargo fmt --all

# Check Rust formatting without changing files.
fmt-check:
    {{ _run }} cargo fmt --all -- --check

# Check workflow syntax and security findings.
github-actions-check: action-lint zizmor
    @echo "✅ GitHub Actions checks complete!"

# File validation
# Keep the command-memory layer itself canonically formatted.
justfile-fmt-check:
    just --fmt --check
    just --justfile tooling/performance.just --fmt --check

# Check code, documentation, and configuration without changing sources.
lint: lint-code lint-docs lint-config

# Check Rust, Python, and shell code without changing sources.
lint-code: rust-core-check python-check shell-check

# Check JSON, TOML, YAML, workflow security, and Just formatting.
lint-config: validate-json toml-check yaml-check github-actions-check justfile-fmt-check

# Check Markdown, spelling, release metadata, and generated changelogs.
lint-docs: markdown-ci docs-version-check changelog-check

# Check active Markdown formatting, local links, and line lengths.
markdown-check: tools-check _ensure-uv
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        case "$file" in
            CHANGELOG.md|docs/archive/*|docs/archives/changelog/*|docs/performance-runs/*) continue ;;
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

# Check active Markdown and spelling without changing files.
markdown-ci: markdown-check spell-check
    @echo "✅ Markdown checks complete!"

# Format active Markdown files; preserve generated and archived records.
markdown-fix: tools-check
    #!/usr/bin/env bash
    set -euo pipefail
    files=()
    while IFS= read -r -d '' file; do
        case "$file" in
            CHANGELOG.md|docs/archive/*|docs/archives/changelog/*|docs/performance-runs/*) continue ;;
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

# Render release docs from complete runs, or read historical CSV/JSON evidence.
performance-doc: python-sync
    {{ _run }} uv run --locked archive-performance report

# Compare published release archives with their original per-release harnesses.
performance-github-assets current_tag="" baseline_tag="": _ensure-gh python-sync
    {{ _run }} uv run --locked archive-performance assets {{ quote(current_tag) }} {{ quote(baseline_tag) }}

# Measure tracked current inputs and the latest release with the captured current harness.
performance-local: _ensure-gh python-sync
    {{ _run }} uv run --locked archive-performance local

# Measure the comparison suite while retaining baseline-phase nalgebra/faer context.
performance-local-non-exact current_tag="" baseline_tag="": _ensure-gh python-sync
    {{ _run }} uv run --locked archive-performance local {{ quote(current_tag) }} {{ quote(baseline_tag) }} --suite vs_linalg --stem target/bench-reports/performance-non-exact

# Publish the existing README plot format from retained measurements after cleanup.
performance-readme metric="lu_solve" stat="median" sample="new" log_y="true": python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--metric {{ quote(metric) }} --stat {{ quote(stat) }} --sample {{ quote(sample) }} --update-readme)
    if [ {{ quote(log_y) }} = "true" ]; then
        args+=(--log-y)
    fi
    {{ _run }} uv run --locked criterion-dim-plot "${args[@]}"

# Measure, retain an immutable complete run, and promote distinct-release documentation.
performance-release current_tag="" baseline_tag="": _ensure-gh python-sync
    {{ _run }} uv run --locked archive-performance release {{ quote(current_tag) }} {{ quote(baseline_tag) }}

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

# Check Python formatting, lint, fixtures, and types.
python-check: python-format-check python-lint python-fixture-lint python-typecheck

# Check Python tooling and run its consumer integration tests.
python-ci: python-format-check python-lint python-fixture-lint python-typecheck test-python
    @echo "✅ Python checks complete!"

# Apply Python lint fixes and formatting.
python-fix: python-sync
    uv run --locked ruff check scripts/ --fix
    uv run --locked ruff format scripts/ tests/semgrep/scripts/

# Lint deliberate Python fixtures with their declared exceptions.
python-fixture-lint: python-sync
    uv run --locked ruff check tests/semgrep/scripts/

# Check Python formatting, including static-analysis fixtures.
python-format-check: python-sync
    uv run --locked ruff format --check scripts/ tests/semgrep/scripts/

# Lint Python support tools and their tests.
python-lint: python-sync
    uv run --locked ruff check scripts/

# Synchronize the locked Python development environment.
python-sync: _ensure-uv
    uv sync --locked --group dev

# Check all Python support and fixture types; fail on every diagnostic.
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

# Check library contracts, docs, static analysis, and dependencies.
rust-core-check: cargo-lock-check fmt-check clippy-core doc-check semgrep semgrep-test unused-deps
    @echo "✅ Rust core checks complete!"

# Run the shared dependency and full-history secret scans.
security: audit security-secrets

# Scan reachable Git history and current tracked/nonignored files with redacted reports.
security-secrets:
    uv run --locked --group dev research-repo-tools security secrets

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

# Check shell script lint and formatting.
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

# Format maintained shell scripts.
shell-fix: _ensure-shfmt
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

# Run library unit tests and default-feature doctests.
test: test-unit test-doc

# Run all Rust tests, doctests, and Python integration tests.
test-all: test-rust test-python
    @echo "✅ All tests passed"

# Smoke-test deterministic inputs and configuration shared with benchmark suites.
test-bench-inputs: tools-check
    {{ _run }} cargo nextest run --workspace --profile ci --features bench,exact --test vs_linalg_inputs --test exact_bench_config --verbose

# Run default-feature Rust doctests.
test-doc:
    {{ _run }} cargo test --doc --verbose

# Run Rust doctests with exact arithmetic enabled.
test-doc-exact:
    {{ _run }} cargo test --features exact --doc --verbose

# Tests for the "exact" feature (exact determinants, conversions, and Bareiss solves)
test-exact: tools-check test-doc-exact
    {{ _run }} cargo nextest run --profile ci --features exact --verbose

# Run Rust integration tests with nextest.
test-integration: tools-check
    {{ _run }} cargo nextest run --profile ci --tests --verbose

# Compile all integration-test targets without running them.
test-integration-compile: tools-check
    {{ _run }} cargo nextest run --all-features --tests --no-run

# Run Python consumer integration and benchmark-policy tests.
test-python: tools-check python-sync
    {{ _run }} uv run --locked pytest -q

# Run all runnable Rust tests and both doctest configurations.
test-rust: test-rust-ci test-doc test-doc-exact
    @echo "✅ Rust tests passed"

# CI Rust bucket: all runnable unit/integration targets in one nextest pass.
test-rust-ci: tools-check
    {{ _run }} cargo nextest run --workspace --release --profile ci --all-features --lib --tests --verbose

# Run library unit tests with nextest.
test-unit: tools-check
    {{ _run }} cargo nextest run --profile ci --lib --verbose

# Check TOML parsing, formatting, and lint without changing files.
toml-check: toml-parse-check toml-fmt-check toml-lint

# Format maintained TOML files.
toml-fix: tools-check
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

# Check TOML formatting without changing files.
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

# Lint maintained TOML files.
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

# Check installed tools without synchronization, downloads, or installation.
tools-check:
    uv run --locked --no-sync --no-python-downloads research-repo-tools toolchain check

# Export verified managed paths to GITHUB_ENV in hosted workflows.
tools-export:
    uv run --locked --no-sync --no-python-downloads research-repo-tools toolchain export

# Check for unused direct Cargo dependencies.
unused-deps: tools-check
    {{ _run }} cargo machete

# Upgrade tools before updating dependency requirements and lock resolutions.
update: update-tools update-dependencies
    @echo "✅ Repository dependencies and tools updated."

# num-bigint and num-rational share public types and must advance together.
update-cargo-dependencies:
    uv run --locked --only-group tooling --inexact research-repo-tools toolchain run -- cargo upgrade --incompatible allow --exclude num-bigint --exclude num-rational
    uv run --locked --only-group tooling --inexact research-repo-tools toolchain run -- cargo update

# Upgrade only declared managed Cargo tools and publish verified TOML pins.
update-cargo-tools:
    uv run --locked --only-group tooling --inexact research-repo-tools toolchain upgrade

# Dependency-only updates leave uv and managed Cargo tool pins unchanged.
update-dependencies: update-cargo-dependencies update-python-dependencies

# Advance direct dev pins, refresh the whole lock, and explicitly synchronize dev.
update-python-dependencies:
    uv run --locked --only-group tooling --inexact research-repo-tools deps update-python
    uv lock --upgrade
    {{ _run }} uv sync --locked --managed-python --group dev

alias update-python-deps := update-python-dependencies

# Upgrade uv through its owner, then declared Cargo tools, then synchronize setup.
update-tools: update-uv update-cargo-tools setup-tools

# Bootstrap outside the project's old uv-version requirement; do not sync here.
update-uv:
    uv run --no-config --no-sync --no-python-downloads research-repo-tools deps update-uv

# Update deterministic release metadata, inferring the previous stable published GitHub release.
[doc('Update package, citation, lockfile, and non-artifact documentation release versions.')]
[positional-arguments]
update-version tag *args: _ensure-gh python-sync
    {{ _run }} research-repo-tools release update "$@"

# Check maintained JSON syntax.
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

# Check YAML and CFF formatting and lint.
yaml-check: yaml-fmt-check yaml-lint

# Format maintained YAML and CFF files.
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

# Check YAML and CFF formatting without changing files.
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

# Lint maintained YAML and CFF files.
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

# Audit workflows with the declared scanner and shared authentication policy.
[positional-arguments]
zizmor *args:
    uv run --locked --group dev research-repo-tools zizmor check "$@"
