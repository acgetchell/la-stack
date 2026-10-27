# shellcheck disable=SC2148
# Justfile for la-stack development workflow
# Install tools: uv run --locked --managed-python --only-group tooling research-repo-tools setup
# Usage: just <command> or just --list

# Use bash with strict error handling for all recipes
set shell := ["bash", "-euo", "pipefail", "-c"]

_run := "uv run --locked --no-sync --no-python-downloads research-repo-tools toolchain run --"

# Semgrep misreads Python 3.15's positive RLIM_INFINITY on Linux.
# Keep its inherited stack finite until semgrep/semgrep#11979 is fixed.
semgrep_stack_limit := if os() == "linux" { 'if [[ "$(ulimit -Hs)" == unlimited ]]; then ulimit -s 1048576; fi;' } else { "" }

# One active-document selection for formatting and checks; retain measured history.
_markdown_inputs := "--include '*.md' --exclude CHANGELOG.md --exclude 'docs/archive/**' --exclude 'docs/archives/changelog/**' --exclude docs/performance.md --exclude 'docs/performance-runs/**'"

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
# Send prerequisite diagnostics to stderr so generated stdout stays clean.
_ensure-uv:
    uv run --locked --no-sync --no-python-downloads research-repo-tools deps check-uv >&2

_ensure-yamllint: _ensure-uv
    uv run --locked yamllint --version >/dev/null

# Require the GitHub release title and body to match the generated notes.
[private]
_release-notes-check tag:
    #!/usr/bin/env bash
    set -euo pipefail
    tag={{ quote(tag) }}
    notes="$(mktemp)"
    trap 'rm -f "$notes"' EXIT
    uv run --locked --no-sync research-repo-tools changelog notes "$tag" > "$notes"
    gh release view "$tag" --repo acgetchell/la-stack --json tagName,name,body |
        jq -e --arg tag "$tag" --rawfile notes "$notes" '.tagName == $tag and .name == $tag and ((.body | gsub("\r\n"; "\n") | rtrimstr("\n")) == ($notes | rtrimstr("\n")))' >/dev/null

# GitHub Actions workflow validation
action-lint: _ensure-actionlint
    uv run --locked research-repo-tools files run --include '.github/workflows/*.yml' --include '.github/workflows/*.yaml' -- actionlint

# Audit the repository's Python and Rust lockfiles with the managed OSV scanner.
audit:
    uv run --locked --group dev research-repo-tools security osv uv.lock Cargo.lock

# Run the fixed-dimension benchmark suites.
bench:
    {{ _run }} cargo bench --locked --workspace --features bench

# Compare every available Criterion row with a saved baseline; exploratory only.
bench-compare baseline="last": python-sync
    {{ _run }} uv run --locked research-repo-tools performance compare target/criterion target/criterion --baseline-sample {{ quote(baseline) }} --format markdown --output target/bench-reports/performance.md

# Compile benchmark targets with warnings denied; do not measure timings.
bench-compile:
    CARGO_BUILD_WARNINGS=deny {{ _run }} cargo bench --locked --workspace --no-run --features bench
    CARGO_BUILD_WARNINGS=deny {{ _run }} cargo bench --locked --no-run --features bench,exact --bench exact
    CARGO_BUILD_WARNINGS=deny {{ _run }} cargo bench --locked --no-run --features bench,exact --bench rational_solve

# Run the exact-arithmetic benchmark suite.
bench-exact:
    {{ _run }} cargo bench --locked --features bench,exact --bench exact

# Run the outward-rounded interval determinant benchmark suite.
bench-interval:
    {{ _run }} cargo bench --locked --features bench --bench interval

# Run the cheaper latest measurements used for latest-vs-last reports.
bench-latest: bench-vs-linalg-la-stack bench-exact

# Run latest measurements and render the latest-vs-last performance report.
bench-latest-vs-last baseline="last": bench-latest
    just bench-compare {{ quote(baseline) }}

# Run the certified dot-product and affine-difference benchmark suite.
bench-linear-form:
    {{ _run }} cargo bench --locked --features bench --bench linear_form

# Check every discovered benchmark before the release workflow packages it.
bench-release-check tag criterion_dir="target/criterion": _ensure-uv
    {{ _run }} uv run --locked research-repo-tools performance verify-saved tooling/criterion-saved.toml --saved-sample {{ quote(tag) }} --criterion-dir {{ quote(criterion_dir) }}

# Discover all release benchmarks and report their expected measurement budget.
bench-release-inventory: _ensure-uv
    #!/usr/bin/env bash
    set -euo pipefail
    [[ ! -e "${CRITERION_HOME:-target/criterion}" ]] || { echo "Release discovery requires a fresh Criterion directory" >&2; exit 1; }
    {{ _run }} uv run --locked research-repo-tools performance discover tooling/criterion-inventory.toml --output target/release-benchmark-inventory.json
    jq -r '.suites | to_entries[] | "[release-baseline] \(.key): \(.value | length) benchmarks; nominal warmup + measurement \((.value | length) * 8 / 60) min; planning estimate at 12 s/benchmark \((.value | length) * 12 / 60) min (not an upper bound)"' target/release-benchmark-inventory.json

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

# Run only la-stack vs_linalg measurements, then compare available Criterion rows.
bench-vs-linalg-latest-vs baseline="last": bench-vs-linalg-la-stack
    just bench-compare {{ quote(baseline) }}

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

# Preview changelog generation; tagged previews default to the prepared citation date.
[positional-arguments]
changelog-preview *args: tools-check python-sync
    #!/usr/bin/env bash
    set -euo pipefail
    {{ _run }} uv run --locked --no-sync python -m tooling.release_changelog --dry-run "$@"

# Generate a prospective release using the citation date unless an ISO date is supplied.
changelog-release tag date="": tools-check python-sync
    {{ _run }} uv run --locked --no-sync python -m tooling.release_changelog --tag {{ quote(tag) }} --date {{ quote(date) }}

alias changelog-unreleased := changelog-release

# Check (non-mutating): run all linters/validators
check: lint
    @echo "✅ Checks complete!"

# Fast compile check (no binary produced)
check-fast:
    {{ _run }} cargo check

# Run the complete local CI validation, tests, examples, and benchmark compilation.
ci: action-lint workflow-allowlist-check zizmor markdown-check spell-check docs-version-check changelog-check toml-parse-check toml-fmt-check toml-lint yaml-fmt-check yaml-lint citation-check validate-json justfile-fmt-check python-check test-python cargo-lock-check fmt-check clippy-all-targets doc-check semgrep semgrep-test unused-deps shell-check test-rust-ci test-doc test-doc-exact bench-compile examples
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

# Summarize an existing Cobertura report without rerunning coverage.
coverage-report report="coverage/cobertura.xml": python-sync
    uv run --locked research-repo-tools coverage report --report {{ quote(report) }}

# Documentation build checks for the default and exact-feature public APIs.
doc-check:
    RUSTDOCFLAGS='-D warnings' {{ _run }} cargo doc --no-deps
    RUSTDOCFLAGS='-D warnings' {{ _run }} cargo doc --no-deps --features exact

# Check synchronized release metadata without changing files.
docs-version-check: python-sync
    uv run --locked --group dev research-repo-tools release check

# Build and run all library examples, including exact arithmetic.
examples:
    {{ _run }} research-repo-tools validation cargo-examples tooling/examples.toml

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
github-actions-check: action-lint workflow-allowlist-check zizmor
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
    {{ _run }} research-repo-tools files run {{ _markdown_inputs }} -- rumdl check
    uv run --locked research-repo-tools files run {{ _markdown_inputs }} --exclude README.md -- research-repo-tools docs check-lines

# Check active Markdown and spelling without changing files.
markdown-ci: markdown-check spell-check
    @echo "✅ Markdown checks complete!"

# Format active Markdown files; preserve generated and archived records.
markdown-fix: tools-check
    {{ _run }} research-repo-tools files run {{ _markdown_inputs }} -- rumdl check --fix

# Render the retained release report offline, with independent eligibility checks.
performance-doc: python-sync
    {{ _run }} uv run --locked research-repo-tools performance promote tooling/performance-report.toml

# Measure tracked and nonignored untracked current inputs against the latest release.
performance-local: _ensure-gh python-sync
    {{ _run }} uv run --locked research-repo-tools performance measure tooling/measurements/full.toml --mode current-vs-latest --repository acgetchell/la-stack --allow-git-mutations --payload target/bench-reports/full/run.json --manifest target/bench-reports/full/evidence.json --report target/bench-reports/full/report.md

# Measure the comparison suite while retaining baseline-phase nalgebra/faer context.
performance-local-non-exact: _ensure-gh python-sync
    {{ _run }} uv run --locked research-repo-tools performance measure tooling/measurements/non-exact.toml --mode current-vs-latest --repository acgetchell/la-stack --allow-git-mutations --payload target/bench-reports/non-exact/run.json --manifest target/bench-reports/non-exact/evidence.json --report target/bench-reports/non-exact/report.md

# Measure two existing tags with the full or non-exact declaration.
performance-measure current_tag baseline_tag view="full": python-sync
    {{ _run }} uv run --locked research-repo-tools performance measure {{ quote("tooling/measurements/" + view + ".toml") }} {{ quote(current_tag) }} {{ quote(baseline_tag) }} --allow-git-mutations --payload {{ quote("target/bench-reports/" + view + "/run.json") }} --manifest {{ quote("target/bench-reports/" + view + "/evidence.json") }} --report {{ quote("target/bench-reports/" + view + "/report.md") }}

# Fetch tags and save independent README source assertions before measurement.
performance-prepare current_tag="" baseline_tag="": _ensure-gh python-sync
    git fetch origin --tags
    uv run --locked --no-sync python -m tooling.release_performance capture {{ quote(current_tag) }} {{ quote(baseline_tag) }}

# Promote the explicitly selected full scratch run; missing inputs never fall back.
performance-promote: python-sync
    {{ _run }} uv run --locked research-repo-tools performance promote tooling/performance-report.toml --payload target/bench-reports/full/run.json --manifest target/bench-reports/full/evidence.json

# Publish the README using prepared source assertions, or four explicit assertions.
performance-readme current_tag="" current_revision="" baseline_tag="" baseline_revision="": python-sync
    {{ _run }} uv run --locked --no-sync python -m tooling.release_performance publish {{ quote(current_tag) }} {{ quote(current_revision) }} {{ quote(baseline_tag) }} {{ quote(baseline_revision) }}

# Measure, retain an immutable complete run, and promote distinct-release documentation.
performance-release current_tag="" baseline_tag="": _ensure-gh python-sync
    {{ _run }} uv run --locked research-repo-tools performance measure tooling/measurements/full.toml {{ quote(current_tag) }} {{ quote(baseline_tag) }} --repository acgetchell/la-stack --allow-git-mutations --payload target/bench-reports/full/run.json --manifest target/bench-reports/full/evidence.json --report target/bench-reports/full/report.md
    {{ _run }} uv run --locked research-repo-tools performance promote tooling/performance-report.toml --payload target/bench-reports/full/run.json --manifest target/bench-reports/full/evidence.json

# Render a named metric/statistic/axis view from explicitly retained or scratch inputs.
plot-vs-linalg current_tag current_revision baseline_tag baseline_revision view="lu_solve-median-linear" input="retained": python-sync
    {{ _run }} uv run --locked research-repo-tools performance publish {{ quote("tooling/plots/" + input + "/" + view + ".toml") }} --artifacts-only --current-release {{ quote(current_tag) }} --current-revision {{ quote(current_revision) }} --baseline-release {{ quote(baseline_tag) }} --baseline-revision {{ quote(baseline_revision) }}

# Check Python formatting, lint, fixtures, and types.
python-check: python-sync
    uv run --locked research-repo-tools python check

# Check Python tooling and run its consumer integration tests.
python-ci: python-check test-python
    @echo "✅ Python checks complete!"

# Apply Python lint fixes and formatting.
python-fix: python-sync
    uv run --locked research-repo-tools python fix

# Lint deliberate Python fixtures with their declared exceptions.
python-fixture-lint: python-sync
    uv run --locked ruff check tests/semgrep/scripts/

# Check Python formatting, including static-analysis fixtures.
python-format-check: python-sync
    uv run --locked research-repo-tools files run --include '*.py' --include '*.pyi' -- ruff format --check --no-force-exclude --no-respect-gitignore

# Lint Python support tools and their tests.
python-lint: python-sync
    uv run --locked research-repo-tools files run --include '*.py' --include '*.pyi' -- ruff check --no-fix --no-fix-only --no-force-exclude --no-respect-gitignore

# Synchronize the locked Python development environment.
python-sync: _ensure-uv
    uv sync --locked --group dev

# Check all Python support and fixture types; fail on every diagnostic.
python-typecheck: python-sync
    uv run --locked research-repo-tools python typecheck

# Verify signed Cargo/native CI evidence; accept --identity JSON and --github-output PATH.
[positional-arguments]
[script("uv", "run", "--locked", "--no-sync", "python")]
release-evidence *args:
    import runpy
    import sys
    sys.argv[1:1] = ["verify"]
    runpy.run_path("tooling/release_evidence.py", run_name="__main__")

# Require the release event, reviewed commit, native checks, assets, and notes.
release-gate tag:
    uv run --locked --no-sync research-repo-tools release gate {{ quote(tag) }}
    gh api "repos/acgetchell/la-stack/releases/$(jq -er '.release.id' "$GITHUB_EVENT_PATH")" | jq -e '.immutable == true' >/dev/null
    just _release-notes-check {{ quote(tag) }}

# Print release notes from the root changelog or a completed minor archive.
release-notes tag: python-sync
    uv run --locked --group dev research-repo-tools changelog notes {{ quote(tag) }}

# Verify native Cargo packaging for the root crate, with default and exact features.
release-package:
    {{ _run }} cargo publish --locked --package la-stack --registry crates-io --dry-run
    {{ _run }} cargo publish --locked --package la-stack --registry crates-io --dry-run --features exact

# Approve the benchmark-equipped draft as the authenticated maintainer.
release-publish tag: (_release-notes-check tag) release-evidence
    uv run --locked --no-sync research-repo-tools release publish {{ quote(tag) }} --approve

# Verify the published GitHub release, required assets, and exact crates.io version.
release-verify tag:
    uv run --locked --no-sync research-repo-tools release verify {{ quote(tag) }} --attempts 7 --interval 10

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
    {{ semgrep_stack_limit }} uv run --locked --group dev semgrep --metrics off --error --strict --timeout 30 --config semgrep.yaml .

# Fixture tests for repository-owned Semgrep rules.
semgrep-test: python-sync
    {{ semgrep_stack_limit }} uv run --locked --group dev research-repo-tools semgrep check-fixtures

# Install declared managed tools and explicitly synchronize dev.
setup: _ensure-gh _ensure-jq
    uv run --locked --managed-python --only-group tooling research-repo-tools setup
    {{ _run }} cargo build

# Tool-only setup; compilation belongs to setup.
setup-tools: _ensure-gh _ensure-jq
    uv run --locked --managed-python --only-group tooling research-repo-tools setup

# Check shell script lint and formatting.
shell-check: _ensure-shellcheck _ensure-shfmt
    uv run --locked research-repo-tools files run --include '*.sh' --batch-size 4 -- shellcheck -x
    uv run --locked research-repo-tools files run --include '*.sh' -- shfmt -d

# Format maintained shell scripts.
shell-fix: _ensure-shfmt
    uv run --locked research-repo-tools files run --include '*.sh' --batch-size 1 -- shfmt -w

# Spell check (typos)
spell-check: tools-check
    {{ _run }} research-repo-tools files run --exclude typos.toml -- typos --config typos.toml --force-exclude --

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
    {{ semgrep_stack_limit }} {{ _run }} uv run --locked pytest -q

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
    {{ _run }} research-repo-tools files run --include '*.toml' -- taplo fmt

# Check TOML formatting without changing files.
toml-fmt-check: tools-check
    {{ _run }} research-repo-tools files run --include '*.toml' -- taplo fmt --check

# Lint maintained TOML files.
toml-lint: tools-check
    {{ _run }} research-repo-tools files run --include '*.toml' -- taplo lint

# Validate TOML syntax independently with Python's standard-library parser.
toml-parse-check: python-sync
    uv run --locked research-repo-tools files run --include '*.toml' -- uv run --locked python -c 'import pathlib, sys, tomllib; [tomllib.loads(pathlib.Path(path).read_text(encoding="utf-8")) for path in sys.argv[1:]]'

# Check installed tools without synchronization, downloads, or installation.
tools-check:
    uv run --locked --no-sync --no-python-downloads research-repo-tools toolchain check

# Preview obsolete managed installations; pass --apply for explicit removal.
[positional-arguments]
tools-clean *args:
    uv run --locked --no-sync research-repo-tools toolchain clean "$@"

# Export verified managed paths to GITHUB_ENV in hosted workflows.
tools-export:
    uv run --locked --no-sync --no-python-downloads research-repo-tools toolchain export

# Check for unused direct Cargo dependencies.
unused-deps: tools-check
    {{ _run }} cargo machete

# Upgrade tools before updating dependency requirements and lock resolutions.
update: update-tools update-dependencies
    @echo "✅ Repository dependencies and tools updated."

# Update explicitly selected Actions pins; --dry-run and --check only preview.
[positional-arguments]
update-actions +args:
    uv run --locked research-repo-tools actions update --policy .github/actions-allowlist.json "$@"

# num-bigint and num-rational share public types and must advance together.
update-cargo-dependencies:
    uv run --locked --only-group tooling --inexact research-repo-tools toolchain run -- cargo upgrade --incompatible allow --exclude num-bigint --exclude num-rational
    uv run --locked --only-group tooling --inexact research-repo-tools toolchain run -- cargo update

# Upgrade only declared managed Cargo tools and publish verified TOML pins.
update-cargo-tools:
    uv run --locked --only-group tooling --inexact research-repo-tools toolchain upgrade

# Dependency-only updates leave uv and managed Cargo tool pins unchanged.
update-dependencies: update-cargo-dependencies update-python-dependencies

# Advance dev pins, then refresh and sync the root Python environment.
update-python-dependencies:
    uv run --locked --only-group tooling --inexact research-repo-tools deps update-python --all-projects

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
    uv run --locked research-repo-tools files run --include '*.json' --batch-size 1 -- jq empty

# Validate external workflow and composite-action dependencies against the repository allowlist.
workflow-allowlist-check: python-sync
    uv run --locked research-repo-tools actions allowlist --policy .github/actions-allowlist.json .github/workflows
    # The v0.2.3 parser accepts workflows only; Semgrep checks composite actions.
    {{ semgrep_stack_limit }} uv run --locked --group dev semgrep --metrics off --error --strict --timeout 30 --config semgrep.yaml .github/actions

# Check YAML and CFF formatting and lint.
yaml-check: yaml-fmt-check yaml-lint

# Format maintained YAML and CFF files.
yaml-fix: tools-check
    {{ _run }} research-repo-tools files run --include '*.yml' --include '*.yaml' --include CITATION.cff -- dprint fmt --incremental=false

# Check YAML and CFF formatting without changing files.
yaml-fmt-check: tools-check
    {{ _run }} research-repo-tools files run --include '*.yml' --include '*.yaml' --include CITATION.cff -- dprint check --incremental=false

# Lint maintained YAML and CFF files.
yaml-lint: _ensure-yamllint
    uv run --locked research-repo-tools files run --include '*.yml' --include '*.yaml' --include CITATION.cff -- yamllint --strict -c .yamllint

# Audit workflows with the declared scanner and shared authentication policy.
[positional-arguments]
zizmor *args:
    uv run --locked --group dev research-repo-tools zizmor check "$@"
