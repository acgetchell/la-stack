# Documentation Guidance

Detailed documentation rules for [AGENTS.md](../../AGENTS.md).

## Contents

- [Document ownership and ordering](#document-ownership-and-ordering)
- [Filename conventions](#filename-conventions)
- [README and rustdoc guides](#readme-and-rustdoc-guides)
- [Links and published destinations](#links-and-published-destinations)
- [Executable examples](#executable-examples)
- [Scientific notation and references](#scientific-notation-and-references)
- [Release metadata and generated files](#release-metadata-and-generated-files)

## Document ownership and ordering

- `README.md` owns concise orientation, quickstart, capabilities, and API
  navigation. Place "Use this crate when" immediately after Introduction.
- `REFERENCES.md` owns bibliographic provenance. Keep the bibliography thematic
  and preserve citation identifiers and deep links.
- `docs/mathematical_basis.md` owns mathematical explanations, assumptions,
  guarantees, and derivations. Put API selection and geometry scope near the
  top, before detailed algorithm discussions.
- Maintain Contents navigation. Sort independent algorithm discussions and
  unordered capability bullets lexicographically within coherent groups;
  preserve prerequisite and procedural order.
- Link between content owners instead of duplicating detail. Keep useful
  explanations and numerical limitations when moving or shortening material.
- `AGENTS.md` is the short entry point for agent rules. Keep detailed Git,
  testing, and documentation policy in the matching `docs/dev/` guide and route
  tasks there from the Required Reading table.
- [Code organization](../code_organization.md) owns module and file placement.
  Update its map when ownership or layout changes; keep command procedures in
  their existing contributor, benchmark, script, and release guides.

## Filename conventions

Under `docs/`, including its subdirectories, name documents by their primary
purpose:

| Purpose | Filename convention | Examples |
|---------|---------------------|----------|
| Task guide with commands and instructions for execution | Uppercase verb or verb phrase, preferably a gerund | `RELEASING.md`, `BENCHMARKING.md`, `MEASURING_COVERAGE.md`, `TUNING_PERFORMANCE.md` |
| Discussion of invariants, principles, policy, or architecture | Lowercase descriptive name | `mathematical_basis.md`, `code_organization.md`, `dev/testing.md` |
| Reference material, analysis, or results report | Lowercase descriptive name | `roadmap.md`, `performance.md` |

A few command examples do not turn a discussion into a task guide. For example,
`dev/testing.md` explains test-design policy; `TESTING.md` would describe how to
execute a testing workflow. Keep conventional directory indexes named
`README.md`. This convention does not change standard filenames at the repository
root, such as `AGENTS.md`, `CONTRIBUTING.md`, and `REFERENCES.md`.

Apply the convention to new documents and deliberate documentation renames.
`MEASURING_COVERAGE.md` describes coverage execution, and `performance.md`
reports measured results. Git and GitHub procedures live in
`dev/MANAGING_CHANGES.md`. Reserve `TUNING_PERFORMANCE.md` for a guide
to optimizing code; `BENCHMARKING.md` owns measurement and comparison workflows.
When renaming a document, update links, ownership maps, validation, and any
generators or configured output paths together. Preserve historical archive
paths; never rename a generated report without updating its source workflow.

## README and rustdoc guides

`src/lib.rs` includes the README with
`#![doc = include_str!("../README.md")]`, so README examples are also the
docs.rs landing-page examples.

Keep the README quickstart and brief capability descriptions discoverable. Put
fuller worked API examples and caller contracts in the documentation-only
`guide` module in `src/lib.rs`. Preserve numerical limitations and error
semantics in the linked guide when shortening a README section.

## Links and published destinations

- Link API references and worked API guides to docs.rs. Keep repository-owned
  mathematical background, benchmark reports, roadmap, contributing, and
  release instructions on GitHub, using absolute URLs to the intended revision.
- README links and Contents anchors must work on both GitHub and generated
  rustdoc. Repository-relative file links may resolve incorrectly from rustdoc;
  verify destinations in both contexts.
- Within Rust docs, use intra-doc links where the referenced item is available
  under the selected features.
- Choose docs.rs `latest` links for intentionally current guidance and explicit
  versions for release-specific contracts. docs.rs builds published crates;
  merging changes does not publish new guide pages. Local rendering does not
  establish availability on the published site.

## Executable examples

- When changing Rust examples in `README.md`, mirror executable versions in
  the private `readme_doctests` module in `src/lib.rs`. Keep mirrors hidden so
  they do not duplicate the landing page, but runnable by `cargo test --doc`.
- README examples requiring optional features may remain `rust,ignore` for
  default-feature doctest compatibility. Give them hidden mirrors gated by the
  matching `#[cfg(feature = "...")]`, and verify the matching feature set,
  for example `cargo test --features exact --doc`.
- Guide examples run directly as doctests. Gate feature-dependent guides with
  the matching feature; remove obsolete private mirrors when examples move out
  of README into guides.
- Run `just doc-check` for changed guide docs and inspect generated pages and
  anchors for explicit links, which rustdoc does not validate. Test changed
  executable examples with the default and matching feature doctest recipes.
- Ordinary Markdown changes use `just markdown-fix` and `just markdown-ci`;
  the latter includes Markdown and spelling checks.

## Scientific notation and references

Unicode mathematics is welcome when it improves readability, including
`×`, `≤`, `≥`, `∈`, `Σ`, `²`, and `2^-50`. State invariants mathematically
where possible, for example `|A[i][i]| > Σ_{j≠i} |A[i][j]|`.

Use numbered `REFERENCES.md` citations such as `\[8\]` and `\[9-10\]` in
algorithm documentation. Link scientific claims to specific reference sections
or entries and describe conditioning behavior. Keep source attribution with the
owning explanation when material moves.

## Release metadata and generated files

Do not bump package versions during ordinary documentation work. When the
maintainer explicitly requests release/version changes, keep README `la-stack`
dependency snippets synchronized with `Cargo.toml` and follow
[Releasing](../RELEASING.md). Update documentation before publishing a new crate
version: crates.io documentation changes also require a new version.

Never edit `CHANGELOG.md` directly. `just changelog` generates, post-processes,
archives, and formats the changelog; `just changelog-unreleased <version>`
prepends unreleased changes. Commit-message guidance lives in
[Git and GitHub guidance](MANAGING_CHANGES.md). Benchmark report and generated asset ownership
belongs in [Benchmarking](../BENCHMARKING.md).
