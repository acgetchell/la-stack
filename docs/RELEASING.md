# Releasing la-stack

Prepare each `vX.Y.Z` release in a dedicated PR. After that PR is merged,
create the annotated tag, publish to crates.io, and create the GitHub release.

The changelog is generated for the target tag before the tag exists, so the
release process does not require a temporary local tag.

Release recipes are content-idempotent where their inputs permit it. Repeating
`update-version` on the same UTC date or repeating `performance-readme` for the
same retained bundle produces no file changes. Crossing UTC midnight
intentionally advances the release date, and repeated `performance-release`
runs can produce different measurements.

## Prepare the environment

Set the target tag once:

```bash
TAG=vX.Y.Z
```

Verify the repository remotes and synchronize `main`:

```bash
gh auth status
git remote -v
git switch main
git pull --ff-only
```

Refresh dependency requirements, lockfiles, and repository-owned Cargo tool
pins before creating the release branch:

```bash
just update
```

Review any tracked changes and land them separately before continuing with the
release PR, then synchronize `main` again. This keeps dependency and tool
upgrades independently reviewable; `just update-version` deliberately does not
run `just update`.

`update-version` and the release performance recipes use GitHub's published
stable releases as the authoritative release history.

## Step 1: Prepare the release PR

Keep the release PR focused on version metadata, the generated changelog,
benchmark artifacts, and release documentation. Major code changes should
already be on `main`.

### 1. Create the release branch

```bash
git switch -c "release/$TAG"
```

### 2. Update release metadata

```bash
just update-version "$TAG"
```

The recipe requires a stable `vX.Y.Z` target that is not older than any
published stable GitHub release and that has at least one earlier published
stable release. It infers the previous release from GitHub and updates the Rust
and Python package metadata and lockfiles, `CITATION.cff`, README dependency and
non-artifact links, and active benchmark examples. `date-released` uses the
current UTC date; if the target changelog section already exists, its date is
updated in the same transaction. `CITATION.cff` retains the Zenodo all-versions
concept DOI.

The recipe validates the synchronized references but does not update dependency
requirements, generate the changelog, or run benchmarks. Review its diff before
continuing.

### 3. Generate the release changelog

```bash
just changelog-unreleased "$TAG"
```

This generates `CHANGELOG.md` as though the target tag already existed, archives
completed minor series under `docs/archive/changelog/`, and synchronizes the
changelog heading with the UTC preparation date recorded in `CITATION.cff`.
Review the generated changelog and any archive changes.

### 4. Generate the release performance comparison

Run this after the package version has been updated:

```bash
just performance-release
```

The no-argument form compares the current package version with the previous
stable published release. Review `docs/PERFORMANCE.md`, any archived comparison
under `docs/archive/performance/`, and the retained CSV and provenance JSON under
`target/bench-reports/`.

The temporary current worktree includes staged and unstaged changes to tracked
files, but excludes untracked files. Stage any new benchmark-relevant file before
running the comparison. Do not run `just clean` or `cargo clean` until the
retained report inputs have been reviewed.

### 5. Refresh the README benchmark comparison

```bash
just performance-readme
```

This consumes the validated CSV and provenance JSON retained by
`just performance-release`; it does not run benchmarks again. It atomically
updates the table and tag-pinned benchmark links in `README.md` with the CSV,
SVG, and provenance JSON under `docs/assets/bench/`. Until it succeeds, those
README links continue to reference the previous published artifacts.

See `docs/BENCHMARKING.md` for repair commands, local comparison modes, artifact
ownership, and provenance details.

### 6. Validate the release branch

```bash
just ci
cargo publish --locked --allow-dirty --dry-run
```

`just ci` includes the lockfile, citation, documentation, test, and benchmark
compile checks.

### 7. Review, stage, and commit the release artifacts

Inspect all changes before staging:

```bash
git status --short
git --no-pager diff
```

Expected release artifacts include package metadata and lockfiles,
`CITATION.cff`, `CHANGELOG.md`, `README.md`, `docs/PERFORMANCE.md`, and generated
files under `docs/archive/` and `docs/assets/bench/`. Stage only the reviewed
paths that were intentionally changed; do not stage the entire `docs/` tree.
Then inspect the staged diff and commit it:

```bash
git --no-pager diff --cached

git commit -m "chore(release): release $TAG

- Bump version to $TAG
- Update citation and utility package metadata
- Generate the release changelog
- Update benchmark and performance artifacts
- Update release documentation"
```

### 8. Push the branch and open the PR

```bash
git push -u origin "release/$TAG"
```

Use `chore(release): release $TAG` as the PR title and describe the PR as a
focused release preparation without feature work.

### Handling fixes found during preparation

For a critical fix that must be included, make and commit the fix, rerun
`just changelog-unreleased "$TAG"`, review and stage only the resulting changelog
files, and commit that generated update separately.

For a non-critical fix, file an issue and defer it to a later release. Do not
hand-edit the generated changelog to add a known-issue note.

## Step 2: Publish after the PR is merged

### 1. Synchronize `main`

```bash
git switch main
git pull --ff-only
```

### 2. Create and verify the annotated tag

```bash
just tag "$TAG"
git --no-pager tag -l --format='%(contents)' "$TAG"
```

`just tag` builds the annotation from the matching active or archived changelog
section. For a changelog larger than 125 kB, the annotation points to that
section instead of embedding it.

### 3. Push the tag

```bash
git push origin "$TAG"
```

### 4. Publish to crates.io

```bash
cargo publish --locked
```

### 5. Create the GitHub release

```bash
gh release create "$TAG" --title "$TAG" --notes-from-tag
```

Keep the release title identical to the tag, including its leading `v`.

### 6. Verify the durable Criterion baseline

After the `Release Benchmarks` workflow completes, verify that the release
contains its long-lived baseline archive:

```bash
gh release view "$TAG" --json assets \
  --jq ".assets[] | select(.name == \"la-stack-$TAG-criterion-baseline.tar.gz\") | .name" | cat
```

The command must print `la-stack-$TAG-criterion-baseline.tar.gz`. A short-lived
Actions artifact is not a substitute for this release asset.

### 7. Remove the merged release branch

After publication and baseline verification succeed:

```bash
git branch -d "release/$TAG"
git push origin --delete "release/$TAG"
```
