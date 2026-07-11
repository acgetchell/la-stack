# Coverage

la-stack uses `cargo-llvm-cov` with `cargo-nextest` for local and CI coverage.
Both coverage recipes use Rust's LLVM source-based instrumentation, run the
same unit and integration test surface, and select nextest's `coverage`
profile:

```bash
cargo llvm-cov nextest --features exact --workspace --lib --tests -P coverage
```

## Local HTML

Generate the local developer report with:

```bash
just coverage
```

The HTML report is written to:

```text
target/llvm-cov/html/index.html
```

The report opens automatically after generation.

## CI reports

Generate the CI-compatible reports with:

```bash
just coverage-ci
```

The Cobertura coverage report is written to:

```text
coverage/cobertura.xml
```

Nextest also writes JUnit test results to:

```text
target/nextest/coverage/test-results/junit.xml
```

The Codecov workflow reads the `cargo-llvm-cov` and `cargo-nextest` versions
from the `justfile`, installs those exact versions plus Rust's
`llvm-tools-preview` component, and runs `just coverage-ci`. It verifies both
reports, uploads the Cobertura and JUnit files to Codecov, and archives each
report directory as a workflow artifact. Local setup via `just setup-tools`
installs the same Rust component and cargo subcommands.

## Migration Notes

- Keep `just coverage-ci` as the single source of truth for CI coverage
  arguments; workflows should install tools and upload artifacts, not duplicate
  the coverage command.
- Keep nextest's `coverage` profile deterministic: retries are disabled,
  timeouts allow for LLVM instrumentation overhead, and JUnit output is
  configured in `.config/nextest.toml`.
- Use `--cobertura --output-path coverage/cobertura.xml` for services that
  consume Cobertura XML.
- Use `--open --output-dir target/llvm-cov` for local reports.
- Preserve the crate's full coverage surface with the `exact` feature and
  `--workspace --lib --tests`.
- `cargo-llvm-cov` excludes workspace `tests/`, `examples/`, and `benches/`
  source files from reports by default, while still allowing integration tests
  to exercise library code. This matches the intended reporting surface here:
  library implementation coverage, not test harness coverage.
- Doc-test coverage remains intentionally disabled because `cargo-llvm-cov`
  marks that path as unstable.
