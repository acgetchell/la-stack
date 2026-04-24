# Coverage

la-stack uses `cargo-llvm-cov` for local and CI coverage. Coverage runs use
Rust's LLVM source-based instrumentation with the same core test selection in
both environments:

```bash
cargo llvm-cov --features exact --workspace --lib --tests
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

## CI XML

Generate the CI-compatible Cobertura report with:

```bash
just coverage-ci
```

The XML report is written to:

```text
coverage/cobertura.xml
```

The Codecov workflow installs Rust's `llvm-tools-preview` component, installs
`cargo-llvm-cov`, caches the installed cargo binary by version, runs
`just coverage-ci`, verifies `coverage/cobertura.xml`, uploads that file to
Codecov, and archives the full `coverage/` directory. Local setup via
`just setup-tools` installs the same Rust component and cargo subcommand.

## Migration Notes

- Keep `just coverage-ci` as the single source of truth for CI coverage
  arguments; workflows should install tools and upload artifacts, not duplicate
  the coverage command.
- Use `--cobertura --output-path coverage/cobertura.xml` for services that
  consume Cobertura XML.
- Use `--open --output-dir target/llvm-cov` for local reports.
- Preserve the crate's full coverage surface with `--features exact
  --workspace --lib --tests`.
- `cargo-llvm-cov` excludes workspace `tests/`, `examples/`, and `benches/`
  source files from reports by default, while still allowing integration tests
  to exercise library code. This matches the intended reporting surface here:
  library implementation coverage, not test harness coverage.
- Doc-test coverage remains intentionally disabled because `cargo-llvm-cov`
  marks that path as unstable.
