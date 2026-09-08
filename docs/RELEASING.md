# Releasing la-stack

Prepare each `vX.Y.Z` release in a dedicated PR. After that PR is merged,
create the annotated tag, publish to crates.io, and create a draft GitHub release.
You then start the `Release Benchmarks` workflow, which attaches the durable
baseline and **automatically publishes that same draft as the public release**.

The local benchmarks used to prepare the release PR do not require a draft.
The draft is needed later as the upload destination for the benchmark archive
produced on GitHub Actions.

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

Install or verify the pinned development tools before running maintenance
recipes. This includes the `cargo-update` package that provides
`cargo-install-update` for `just update`:

```bash
just setup
```

Refresh Cargo dependency requirements, exact Python development-tool pins,
lockfiles, repository-owned Cargo tool pins, and the active uv pin before
creating the release branch:

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
stable published release. Review `docs/performance.md`, any archived comparison
under `docs/archive/performance/`, and the complete versioned local summary
snapshot and `latest.json` under `docs/performance/`. Include those CSV and JSON
files in the release commit; they preserve every recorded case, including the
measurements outside the readable report's selection.

The temporary current worktree includes staged and unstaged changes to tracked
files, but excludes untracked files. Stage any new benchmark-relevant file before
running the comparison. Successful release promotion preserves the selected
report inputs and complete summaries under `docs/performance/`, so cleanup no
longer removes the data needed to regenerate the reports. Local experiments
from `performance-local` remain scratch until explicitly promoted.

### 5. Refresh the README benchmark comparison

```bash
just performance-readme
```

This consumes the validated CSV and provenance JSON retained by
`just performance-release`; it does not run benchmarks again. It atomically
updates the table and tag-pinned benchmark links in `README.md` with the CSV,
SVG, and provenance JSON under `docs/assets/bench/`. Until it succeeds, those
README links continue to reference the previous published artifacts.
If scratch inputs were cleaned, it uses the committed local snapshot. Existing
partial or corrupt scratch inputs remain errors and never trigger a fallback.

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
`CITATION.cff`, `CHANGELOG.md`, `README.md`, `docs/performance.md`, and generated
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

### 5. Create the draft GitHub release

Create an unpublished release where the workflow can attach its benchmark archive:

```bash
gh release create "$TAG" --verify-tag --draft --title "$TAG" --notes-from-tag
```

Keep the release title identical to the tag, including its leading `v`.
Leave it as a stable draft: publishing now would make the release immutable
before the benchmark archive can be attached. GitHub recommends this
[draft-first sequence for immutable releases](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases).

### 6. Run benchmarks and publish the draft

**Starting this workflow also authorizes automatic publication.** You do not
need a separate command to convert the draft into a published release.

Start the workflow from the release tag, passing the same tag as its input:

```bash
gh workflow run release-benchmarks.yml --ref "$TAG" -f tag="$TAG"
```

Once GitHub accepts the command, you can close the terminal and leave the
workflow running. It runs on GitHub's runners and needs no further input to
publish the release. Check its status later, when convenient, on the
[Release Benchmarks Actions page](https://github.com/acgetchell/la-stack/actions/workflows/release-benchmarks.yml).

The workflow runs from the tagged version of its definition, so the tag must
include this release workflow. Using the tag as the workflow ref keeps execution
in that tag's cache scope. Dispatching from `main` or a different tag is rejected.

The workflow automatically:

1. Runs the full benchmark suites for the tagged commit.
2. Uploads the archive to the draft and verifies its name, state, size, and
   SHA-256 digest.
3. Removes the draft status, making the same release public with the archive
   already attached.

A benchmark, upload, or verification failure leaves the release as a draft.
Follow [Recovering a failed run](#recovering-a-failed-run) to retry.

Before benchmarking, the workflow requires
exactly one mutable, non-prerelease draft whose tag and title match the stable
`vX.Y.Z` input, and resolves the existing tag to a commit. It benchmarks that
commit only if it matches the workflow's own commit, then rechecks the captured
release ID and tag commit before attaching
the archive and again before publication. Do not edit the draft, move the tag,
or publish manually while the workflow is running.

Runs for the same tag are serialized without cancelling an active run.

### 7. Verify publication and the durable Criterion baseline

After the workflow completes, open the repository's
[Releases page](https://github.com/acgetchell/la-stack/releases). Confirm that
`$TAG` is published and its assets include
`la-stack-$TAG-criterion-baseline.tar.gz`. This is a check of the completed
publication; the workflow does not wait for you to perform it.

For an optional command-line check:

```bash
gh release view "$TAG" --json isDraft,assets \
  --jq "select(.isDraft == false) | .assets[] | select(.name == \"la-stack-$TAG-criterion-baseline.tar.gz\") | .name" | cat
```

The command must print `la-stack-$TAG-criterion-baseline.tar.gz`. A short-lived
Actions artifact is not a substitute for this release asset.

The benchmark producer has read-only repository permissions and restores no
dependency caches, including tool binaries. It disables Rust toolchain and
`setup-just` caching and installs the pinned just and cargo-nextest versions
with `cargo install --locked`. The short preflight job receives `contents: write`
for draft visibility; the separate publisher receives it to attach the archive
and publish. Neither privileged job checks out or executes repository code.

The producer budgets 150 minutes for `vs_linalg` and 90 minutes for `exact`
within a 285-minute job. Both suites retain full release sampling. Inventory
and raw-data validation must succeed before the single complete archive is
packaged and uploaded; inspect the suite timing summary when a run fails.
See the [hosted runtime budget](BENCHMARKING.md#hosted-release-runtime-budget)
for measured history, capacity estimates, and headroom, and
[workflow validation](BENCHMARKING.md#validate-the-release-workflow) for local
regression checks and hosted evidence requirements.

### Recovering a failed run

If preflight or the producer failed, rerun all jobs with `gh run rerun <run-id>`
or dispatch again. The producer requires a preflight from the current attempt;
rerunning failed benchmark jobs alone stops before tool setup or measurement
because the earlier draft check is stale.

If the producer succeeded and only the publisher failed, preserve the existing
30-day Actions artifact and rerun failed jobs:

```bash
gh run rerun <run-id> --failed
```

The publisher reuses an already uploaded draft asset only if its state, size,
and digest match the downloaded archive. It never overwrites or deletes assets.
A conflicting or incomplete upload stops the run; inspect the draft and remove
only that conflicting asset manually before retrying. An expired temporary
artifact requires a new full run.

A full rerun or new dispatch creates a distinct Actions artifact named
`bench-baseline-$TAG-<run-id>-<producer-attempt>`. Fresh measurements may differ
from an earlier uploaded draft asset, so resolve any conflict before retrying
publication. If the release is already published, every rerun fails closed
without changing it, including when an earlier publication succeeded but its
response was lost. Verify the public release in that case. Existing immutable
releases with missing assets cannot be repaired by this workflow.

### 8. Remove the merged release branch

After publication and baseline verification succeed:

```bash
git branch -d "release/$TAG"
git push origin --delete "release/$TAG"
```
