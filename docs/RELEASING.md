# Releasing la-stack

This guide documents the release flow for `vX.Y.Z`: prepare a dedicated
release PR, merge it, create the final annotated tag from the generated
changelog, publish to crates.io, and create the GitHub release.

The release changelog is generated with `git-cliff --tag` through
`just changelog-unreleased`, so no temporary local tag is needed.

Applies to versions vX.Y.Z. Prefer updating documentation before publishing
to crates.io.

---

## Conventions and environment

Set these variables to avoid repeating the version string:

```bash
# tag has the leading v, version does not
TAG=vX.Y.Z
VERSION=${TAG#v}
PREVIOUS_TAG=vA.B.C
```

Verify your git remotes:

```bash
git remote -v
```

Ensure your local `main` is up to date before beginning:

```bash
git switch main
git pull --ff-only
```

---

## Step 1: Create a clean release PR

This PR should primarily include version bumps, changelog updates, benchmark
comparison updates, and documentation updates. All major code changes should
already be on `main`.

Finalize release-facing metadata and documentation in this dedicated release
PR. Ordinary feature, fix, review, and hygiene work should not preemptively
bump versions or prepare release artifacts.

Small, critical fixes discovered during the release process may be included,
but keep them minimal and release-critical.

1. Create the release branch

```bash
git switch -c "release/$TAG"
```

2. Bump versions

Preferred, if `cargo-edit` is installed:

```bash
cargo set-version "$VERSION"
```

Alternative: edit `Cargo.toml` manually and update `version = "..."` under
`[package]`.

Update release metadata to match the crate version:

- `CITATION.cff`: update `version` and `date-released`
- `pyproject.toml`: update `[project] version` for the Python utility package

Review the citation identity fields at the same time: author name and contact,
ORCID, repository URL, and license. Preserve la-stack's Zenodo concept DOI
(`all versions`) unless the archival policy is deliberately changed; do not
replace it with a release-specific DOI.

Refresh both committed lockfiles after those manual metadata edits:

```bash
cargo metadata --format-version 1 --no-deps > /dev/null
uv lock
```

Review version references in documentation:

```bash
uv run --locked check-docs-version-sync
```

The automated check covers package metadata, lockfiles, README dependency
snippets, and release-pinned README links. Then review historical references
that intentionally remain on older versions:

```bash
rg -n "\bv?[0-9]+\.[0-9]+\.[0-9]+\b" README.md docs/ CITATION.cff pyproject.toml || true
```

3. Generate the release changelog

```bash
# Generates CHANGELOG.md as though TAG already exists, then applies
# markdown hygiene and archives completed minor release series.
just changelog-unreleased "$TAG"
```

`just changelog-unreleased` runs
`GIT_CLIFF_OFFLINE=true git-cliff --tag "$TAG" -o CHANGELOG.md`, then
`postprocess-changelog`, then `archive-changelog`. The root changelog keeps
Unreleased plus the active minor series; older completed minor series live
under `docs/archive/changelog/`.

4. Run benchmarks and update the README comparison table

```bash
# Validate inputs, run a fresh complete vs_linalg benchmark, and atomically
# update the README table plus CSV/SVG/JSON-provenance assets
just plot-vs-linalg-readme
```

Review the updated table in `README.md`, the plot and CSV in `docs/assets/`, and
the adjacent provenance JSON. The publication command fails if the independent
correctness gate, canonical-dimension coverage, or required provenance is
incomplete.

5. Update the release performance comparison

```bash
# Infers TAG from Cargo.toml, compares it against the previous stable published
# release, writes docs/PERFORMANCE.md, and archives the previous docs/PERFORMANCE.md
# under docs/archive/performance/.
just performance-release
```

Review `docs/PERFORMANCE.md` for the latest release-to-release comparison. Older
committed comparisons are archived under `docs/archive/performance/` with
lexicographically sorted filenames such as `v0.4.2-vs-v0.4.1.md`. Iterative
local reports still live under `target/bench-reports/`. For an explicit release
repair, run `just performance-release <current-tag> <previous-tag>`. To compare
the stored GitHub Actions release assets instead of running cargo locally, use
`just performance-github-assets`. The local release workflow validates and then
compiles both library revisions with the current checkout's hashed benchmark
harness, recording source-state, environment, toolchain, dependency, Criterion,
and validation provenance. Stored release assets retain their original
per-release harnesses; unavailable historical measurement metadata is labelled
explicitly rather than treated as an isolated library-code comparison.

After the GitHub Release is published, the `Release Benchmarks` workflow checks
out the release tag, runs the independent benchmark-input tests, saves a full
Criterion baseline, and attaches
`la-stack-$TAG-criterion-baseline.tar.gz` to the release. That release asset is
the durable archive for historical baseline comparisons; the workflow also
uploads a short-lived Actions artifact for debugging the run.

See `docs/BENCHMARKING.md` for local saved-baseline workflows and the full
comparison command reference.

6. Validate the release branch

```bash
just ci
just cargo-lock-check
just citation-check
cargo publish --locked --allow-dirty --dry-run
```

7. Stage and commit release artifacts

```bash
git add Cargo.toml Cargo.lock CITATION.cff pyproject.toml uv.lock CHANGELOG.md README.md docs/

git commit -m "chore(release): release $TAG

- Bump version to $TAG
- Update citation and utility package metadata
- Update changelog with latest changes
- Update benchmark comparison table and release performance report
- Update documentation for release"
```

8. Push the branch and open a PR

```bash
git push -u origin "release/$TAG"
```

PR metadata:

- Title: chore(release): release $TAG
- Description: Clean release PR with version bump, changelog, and
  documentation updates. No feature work.

### Handling fixes discovered during the release process

If you discover issues after generating the changelog:

1. For critical fixes that must be in this release, make and commit the fix,
   then regenerate the release changelog:

   ```bash
   just changelog-unreleased "$TAG"
   git add CHANGELOG.md docs/archive/changelog/
   git commit -m "docs: update changelog with release fixes"
   ```

2. For non-critical fixes, document them as known issues in the release notes
   or include them in the next release.

---

## Step 2: After the PR is merged into main

1. Sync your local `main` to the merge commit

```bash
git switch main
git pull --ff-only
```

2. Create the final annotated tag using the changelog content

```bash
# Creates the annotated tag from the matching CHANGELOG.md section.
# Archived versions are read from docs/archive/changelog/ automatically.
# For large changelogs (>125KB), the tag message points to the changelog
# section instead of embedding the full content.
just tag "$TAG"
```

3. Optional: verify the tag message content

```bash
git tag -l --format='%(contents)' "$TAG"
```

4. Push the tag

```bash
git push origin "$TAG"
```

5. Publish to crates.io

```bash
# Publish the crate (ensure docs are already updated on main via the PR)
cargo publish --locked
```

6. Create the GitHub release with notes from the tag annotation

```bash
# Requires GitHub CLI (gh) and authenticated session
gh release create "$TAG" --title "$TAG" --notes-from-tag
```

Always set the GitHub release title to the exact tag string, including the
leading `v`.

7. Verify the durable Criterion baseline asset

After the `Release Benchmarks` workflow completes, verify that the GitHub
release contains the expected long-lived baseline archive:

```bash
gh release view "$TAG" --json assets \
  --jq ".assets[] | select(.name == \"la-stack-$TAG-criterion-baseline.tar.gz\") | .name" | cat
```

The command must print `la-stack-$TAG-criterion-baseline.tar.gz`. An Actions
artifact alone is not a durable release baseline.

8. Clean up the merged release branch

After publishing and asset verification succeed, remove the release branch
locally and on the remote:

```bash
git branch -d "release/$TAG"
git push origin --delete "release/$TAG"
```

---

## Notes and tips

- Do not create a temporary local release tag for changelog generation; use
  `just changelog-unreleased "$TAG"`.
- Keep the release PR scoped to version, changelog, archive, benchmark
  comparison, and documentation changes.
- `just changelog` regenerates the current changelog from existing tags and may
  update `docs/archive/changelog/`.
- `just changelog-unreleased "$TAG"` is for release PR preparation before the
  final tag exists.
- `just tag "$TAG"` is for the final post-merge annotated tag.
- If multiple crates or files reference the version, confirm all of them are
  updated consistently.
