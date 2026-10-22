# Releasing la-stack

Prepare each `vX.Y.Z` release in a dedicated PR. After that PR is merged,
create the annotated tag and a draft GitHub release. Start `Release Benchmarks`
to attach the durable baseline, then review and publish the draft with
`just release-publish`. The published release starts `Publish to crates.io`;
approve its `crates-io` deployment to permit the Cargo upload.

The local benchmarks used to prepare the release PR do not require a draft.
The draft is needed later as the upload destination for the benchmark archive
produced on GitHub Actions.

The changelog is generated for the target tag before the tag exists, so the
release process does not require a temporary local tag.

Release recipes are content-idempotent where their inputs permit it. Repeating
`update-version` on the same UTC date or repeating `performance-readme` for the
same retained bundle produces no file changes. Rerunning `update-version` on a
later UTC date intentionally advances the declared release date. Crossing
midnight alone requires no update: the `declared` date policy allows preparation
and tagging on different UTC days when citation and changelog dates agree.
Repeated `performance-release` runs can produce different measurements.

## Contents

- [Reusable CI evidence](#reusable-ci-evidence)
- [Prepare the environment](#prepare-the-environment)
- [Step 1: Prepare the release PR](#step-1-prepare-the-release-pr)
- [Step 2: Publish after the PR is merged](#step-2-publish-after-the-pr-is-merged)
- [Configure or verify Trusted Publishing](#configure-or-verify-trusted-publishing)

## Reusable CI evidence

The publishing process uses the installed `research-repo-tools==0.2.3`
validation helpers. `tooling/release_evidence.py` supplies the Cargo-specific
policy; it does not use the shared package's wheel/sdist release commands.
Every successful protected `main` push produces signed evidence after all three
native `just ci` jobs and the default/exact Cargo packaging job succeed.
The `release-evidence` job signs their exact commit, run, attempt, job IDs and
package identity, then retains the manifest and Sigstore bundle in an immutable
Actions artifact for 30 days.

`just release-evidence` verifies that evidence for the current checkout. It
requires one unambiguous successful `push` run on `main` from this repository's
active CI workflow, every required job and step from that same attempt, and a
nonexpired artifact retrieved by immutable ID. The shared downloader checks its
authoritative ZIP size and SHA-256 digest and exact file inventory. GitHub CLI
then verifies the signed manifest's workflow, source commit, `main` ref, hosted
runner, push event, run and attempt. The manifest must match the live job and
package identities. Pending matching runs receive a bounded wait; missing,
failed, ambiguous or expired evidence stops publication. PR, manual dispatch,
fork, local and partial-rerun results cannot replace it.

The publishing workflow reuses these completed checks, captures the run/attempt
and artifact ID/digest, and revalidates that same identity after deployment
approval. It also retains the independent release metadata, exact-commit SARIF,
durable Criterion asset, immutable GitHub release and registry-absence gates.
Release benchmarks remain a separate scientific measurement workflow; CI
evidence does not replace their full sampling or durable archive.

This reuses validation of the source. Cargo still packages the clean release
commit when uploading, so it does not promise byte-for-byte promotion of the
dry-run `.crate` archive. The shared tooling's Python distribution promotion
process has a different artifact contract.

Use a current GitHub CLI with `gh attestation verify` and the signer/source
constraints used by the adapter. See GitHub's
[attestation verification options](https://cli.github.com/manual/gh_attestation_verify).

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

Install or verify the declared development tools before running maintenance
recipes. Shared setup installs managed tools and explicitly synchronizes `dev`:

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
package metadata and `Cargo.lock`, `CITATION.cff`, and README dependency metadata.
The dependency-only Python environment retains its own project and lockfile
version independently of the Rust release.
Active navigation stays on `main`; measured artifact links retain their
recorded revision. `date-released` uses the
current UTC date; if the target changelog section already exists, its date is
updated in the same transaction. `CITATION.cff` retains the Zenodo all-versions
concept DOI.

The recipe validates the synchronized references but does not update dependency
requirements, generate the changelog, or run benchmarks. Review its diff before
continuing.

### 3. Generate the release changelog

```bash
RELEASE_DATE="$(date -u +%F)"
just changelog-unreleased "$TAG" "$RELEASE_DATE"
```

This generates `CHANGELOG.md` as though the target tag already existed, archives
completed minor series under `docs/archives/changelog/`, and uses the explicit
ISO release date. Verify that `RELEASE_DATE` equals the date recorded in
`CITATION.cff` by metadata preparation. If preparation crossed UTC midnight,
either reuse that declared date or deliberately advance it by rerunning
`just update-version "$TAG"` first. Generation preserves existing release dates
and does not update package or citation metadata.
Use `just changelog-preview --tag "$TAG" --date "$RELEASE_DATE"` for a read-only
preview, and `just changelog-check` to validate the root and every archive.
Review the generated changelog and any archive changes. The
[Tooling guide](../tooling/README.md#changelog-and-release-tooling) explains
the supported `dependency-bodies = "preserve"` policy in `pyproject.toml`, which
preserves authored dependency-note links, text, and code without a local template.

If generation reports `conflicting retained release`, compare the regenerated
section with the retained archive and its original tagged commit history. Both
versions already pass through the same formatter before comparison, so this
error can indicate a content-policy change rather than whitespace differences.
For example, adopting `dependency-bodies = "preserve"` can restore original
Dependabot release-note and comparison links omitted by an older template.
Generate an archive candidate in a disposable directory with the published CLI,
review the policy migration, and replace only the affected generated archive.
Preserve release dates, commit identities, and authored notes, then rerun the
original command and `just changelog-check`.

### 4. Generate the release performance comparison

After updating the package version, fetch the tags and select the previous
stable release shown by GitHub's publication history:

```bash
just performance-prepare
gh release list --exclude-drafts --exclude-pre-releases --json tagName,publishedAt
```

Set `PREVIOUS_TAG` to that reviewed baseline. Capture the source commit
assertions independently before measurement:

```bash
CURRENT_COMMIT=$(git --no-pager rev-parse HEAD)
PREVIOUS_COMMIT=$(git --no-pager rev-parse "$PREVIOUS_TAG^{commit}")
just performance-release
```

The no-argument measurement compares the prospective package version with the
latest stable published release. Review `docs/performance.md` and the complete
new run and index changes under `docs/performance-runs/`. Its current/baseline
tags and commits must match the selected assertions. Both labels must be v0.5.0
or newer.

Evidence retains both phases, mean and median estimates, 100 raw samples per
case and 95% intervals. Current declarations record acquisition policy as well
as source/harness provenance. Future runs support offline rendering after
scratch cleanup; `just performance-doc` explicitly selects retained history.
`just performance-promote` explicitly selects the full scratch pair.

Historical artifacts, including the v0.5.0/v0.4.6 comparison, remain frozen with
their original bytes and links. That older run lacks the new acquisition
declaration and is no longer a regeneration input. Never backfill provenance
or relax current eligibility checks to replay it.

### 5. Refresh the README benchmark comparison

```bash
just performance-readme "$TAG" "$CURRENT_COMMIT" "$PREVIOUS_TAG" "$PREVIOUS_COMMIT"
```

This consumes archived, release-eligible evidence without benchmarking again.
It verifies the independent identities, current source/harness membership and
hashes, package/report references and existing-tag artifact identity.
The canonical LU-solve median/log view is declared in
`tooling/performance-readme.toml`.

The installed publisher formats and checks the candidate README before
atomically replacing its benchmark section and CSV/SVG outputs. Links to the
complete-run payload and evidence envelope replace new plot-provenance sidecars;
original sidecars remain unchanged. Tag-pinned links become available when the
release tag is pushed. Until then, preview the SVG directly under
`docs/assets/bench/`.

See [Benchmarking](BENCHMARKING.md) for explicit scratch versus retained plot
inputs, named metric/statistic/axis views, output ownership and methodology.

### 6. Validate the release branch

```bash
just ci
cargo publish --locked --package la-stack --registry crates-io --allow-dirty --dry-run
cargo publish --locked --package la-stack --registry crates-io --allow-dirty --dry-run --features exact
```

`just ci` includes the lockfile, citation, documentation, test, and benchmark
compile checks.
The preparation dry-runs permit the uncommitted release artifacts. After review
and merge, main CI runs `just release-package` on the clean release source with
both default and `exact` features and explicitly selects only `la-stack`.
Those hosted results become the reusable packaging evidence.

### 7. Review, stage, and commit the release artifacts

Inspect all changes before staging:

```bash
git status --short
git --no-pager diff
```

Expected release artifacts include package metadata and lockfiles,
`CITATION.cff`, `CHANGELOG.md`, `README.md`, `docs/performance.md`, and generated
files under `docs/archive/`, `docs/archives/changelog/`, `docs/assets/bench/`, and
`docs/performance-runs/`. Stage only the reviewed
paths that were intentionally changed; do not stage the entire `docs/` tree.
Then inspect the staged diff and commit it:

```bash
git --no-pager diff --cached

git commit -m "chore(release): release $TAG

- Bump version to $TAG
- Update citation and README release metadata
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
`just changelog-unreleased "$TAG" "$RELEASE_DATE"`, review and stage only the resulting changelog
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

### 4. Verify the release commit

```bash
just release-evidence
```

Require successful native `CI` builds on Linux, macOS and Windows, plus Clippy
and repository-rule SARIF checks, for this exact commit before proceeding.
These builds include the declared Rust 1.99.1 baseline, default/exact doctests,
all-feature tests, interval and linear-form regressions, examples and benchmark
compilation. The publishing policy in `pyproject.toml` binds all five check names
to the GitHub Actions App. CodeRabbit and human review remain protected-PR
requirements; CodeRabbit does not run on the merged `main` commit. Missing,
pending, cancelled or failed native checks block publication. The evidence
command also checks the original main CI packaging and signing jobs. If no
eligible push run exists, stop and correct the event-delivery problem before
selecting the release commit. A manually dispatched CI run cannot produce
release evidence. A failed or expired main-push run requires a complete rerun
of that original run, with fresh evidence from the new attempt.

### 5. Create the draft GitHub release

Create an unpublished release where the workflow can attach its benchmark archive:

```bash
NOTES="$(mktemp)"
just release-notes "$TAG" > "$NOTES"
gh release create "$TAG" --verify-tag --draft --title "$TAG" --notes-file "$NOTES"
rm "$NOTES"
```

Keep the release title identical to the tag, including its leading `v`.
Leave it as a stable draft: publishing now would make the release immutable
before the benchmark archive can be attached. GitHub recommends this
[draft-first sequence for immutable releases](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases).

### 6. Run benchmarks and review the draft

Start the workflow from the release tag, passing the same tag as its input:

```bash
gh workflow run release-benchmarks.yml --ref "$TAG" -f tag="$TAG"
```

Once GitHub accepts the command, you can close the terminal and leave the
workflow running. Check its status later on the
[Release Benchmarks Actions page](https://github.com/acgetchell/la-stack/actions/workflows/release-benchmarks.yml).

The workflow runs from the tagged version of its definition, so the tag must
include this release workflow. Using the tag as the workflow ref keeps execution
in that tag's cache scope. Dispatching from `main` or a different tag is rejected.

The workflow automatically:

1. Runs the full benchmark suites for the tagged commit.
2. Uploads the archive to the draft and verifies its name, state, size, and
   SHA-256 digest.
3. Leaves the draft and verified archive ready for maintainer review.

A benchmark, upload, or verification failure leaves the release as a draft.
Follow [Recovering a failed run](#recovering-a-failed-run) to retry.

Before benchmarking, the workflow requires
exactly one mutable, non-prerelease draft whose tag and title match the stable
`vX.Y.Z` input, and resolves the existing tag to a commit. It benchmarks that
commit only if it matches the workflow's own commit, then rechecks the captured
release ID and tag commit before attaching
the archive. Do not edit the draft, move the tag,
or publish manually while the workflow is running.

The privileged jobs install the exact published `research-repo-tools==0.2.3`
package on Python 3.15 without repository checkout. Preflight emits versioned
target JSON with the captured ID, commit and title assertion. Its trusted job
output passes directly to the publisher, separately from the producer's inert
archive. The shared writer verifies downloaded asset bytes as well as provider
metadata. It leaves publication to the authenticated maintainer.

Runs for the same tag are serialized without cancelling an active run.

### 7. Publish the reviewed draft and approve Cargo publication

After the workflow completes, open the repository's
[Releases page](https://github.com/acgetchell/la-stack/releases). Review the
draft, generated notes and `la-stack-$TAG-criterion-baseline.tar.gz` asset.
From the clean tagged checkout, using your human GitHub CLI session, run:

```bash
just release-publish "$TAG"
gh run list --workflow publish-crates.yml --limit 5
```

This command authorizes GitHub Release publication. The shared gate checks the
clean checkout, stable tag/package metadata, remote tag commit, protected default
branch ancestry, successful exact-commit checks and required baseline asset.
The local recipe also verifies signed main CI evidence and requires the release
title and body to match the tag and generated notes. Human authentication is essential: publishing with a workflow's
`GITHUB_TOKEN` does not trigger the separate `release: published` workflow.

The publishing workflow checks out `github.sha`, rechecks the published event's
repository, release ID, tag/ref/commit and evidence, then verifies the original
signed main CI evidence without repeating the native checks or packaging builds.
It also reads the captured release ID to require GitHub's immutable state.
After it succeeds, approve the `crates-io` environment deployment in GitHub.
The upload job rechecks the same selected run, attempt and artifact ID/digest,
the release gate, and proven absence of the exact registry version after approval,
then obtains a short-lived OIDC token. Within the publishing workflow, only that
job has `id-token: write`; it has no repository write permission or dependency caches.
Uploads for `la-stack` are serialized without cancellation.

The original main CI job builds the default and `exact` packages without upload
credentials. After approval, Cargo repackages the same clean, revalidated commit
with `--locked --no-verify`, so build scripts do not run with the crates.io token.
The dry-run archive is not transferred between jobs, and byte-identical archives
are not asserted. This follows [Cargo's native publication model](https://doc.rust-lang.org/cargo/commands/cargo-publish.html).
No benchmark workspace member is published.

### 8. Verify registry availability and downstream resolution

After the publishing workflow completes:

```bash
just release-verify "$TAG"
```

This requires the published stable GitHub release, nonempty durable baseline
asset and exact non-yanked crates.io version with a checksum. It retries only
proven registry absence, never an upload. A short-lived Actions artifact or a
published GitHub release alone is insufficient.

Check a clean consumer with a fresh Cargo home before dependent Delaunay releases:

```bash
CONSUMER="$(mktemp -d)"
CARGO_HOME="$(mktemp -d)"; export CARGO_HOME
RELEASE_CRATE="la-stack@=${TAG#v}"
cargo init --vcs none --name la_stack_release_check "$CONSUMER"
cargo add --manifest-path "$CONSUMER/Cargo.toml" "$RELEASE_CRATE"
cargo check --locked --manifest-path "$CONSUMER/Cargo.toml"
cargo add --manifest-path "$CONSUMER/Cargo.toml" "$RELEASE_CRATE" --features exact
cargo check --locked --manifest-path "$CONSUMER/Cargo.toml"
```

Use a disposable shell for this check. Preserve the output and record the tag,
commit, benchmark run, publishing run, registry URL/checksum and consumer results
in the current release's tracking issue or merged release PR. Keep any release
tracking issue open until publication and consumer verification succeed.

The benchmark producer has read-only repository permissions and restores no
dependency caches, including tool binaries. It invokes `setup-tools` with
`cache: false`, installs declared managed Rust and Cargo tools through the locked
shared package, and exports verified paths. The short preflight job receives `contents: write`
for draft visibility; the separate asset writer receives it to attach the archive.
Neither privileged benchmark job checks out or executes repository code.

The producer budgets 150 minutes for `vs_linalg` and 90 minutes for `exact`
within a 350-minute job, including a 90-minute preparation budget for uncached
tool installation, compilation, inventory and input validation.
Both suites retain full release sampling. Inventory
and raw-data validation must succeed before the single complete archive is
packaged and uploaded; inspect the suite timing summary when a run fails.
See the [hosted runtime budget](BENCHMARKING.md#hosted-release-runtime-budget)
for measured history, capacity estimates, and headroom, and
[workflow validation](BENCHMARKING.md#validate-the-release-workflow) for local
regression checks and hosted evidence requirements.

### Recovering a failed run

For reusable CI evidence, rerun **all** jobs in the original main-push CI run.
Failed-jobs-only reruns omit producers or validators from the new attempt and
are rejected. A complete rerun creates a new signed artifact; it never changes
the earlier artifact. If CI was rerun while deployment approval was pending,
the pinned old attempt is stale: restart the publishing workflow's verification
and obtain approval for the newly selected evidence. Do not bypass expiration
or substitute manual-dispatch results. Local builds are useful diagnostics but
cannot renew hosted attestations.

The following recovery rules apply separately to `Release Benchmarks`.
If preflight or the producer failed, rerun all jobs with `gh run rerun <run-id>`
or dispatch again. The producer requires a preflight from the current attempt;
rerunning failed benchmark jobs alone stops before tool setup or measurement
because the earlier draft check is stale.

A rerun uses the original workflow at the release tag; it does not pick up a
workflow fix merged later. If the tagged workflow needs repair, keep the release
as a draft and have the maintainer review the corrected release commit and tag
before dispatching again. Never move a published tag.

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

If the crates.io upload fails or Cargo times out, inspect the original run and
exact version before retrying. The workflow checks visibility even after an
upload failure but preserves the failed job result. A visible version makes a
rerun stop before OIDC/upload; use `just release-verify` and the clean consumer
check to finish verification. A registry transport/authentication error does not
prove absence. After a visibility timeout, investigate the accepted-upload state
before authorizing a rerun for the unchanged tag; never blindly retry an upload.
Corrections to published source require a new version. Never move a published tag.

### 9. Remove the merged release branch

After publication and baseline verification succeed:

```bash
git branch -d "release/$TAG"
git push origin --delete "release/$TAG"
```

## Configure or verify Trusted Publishing

Follow the current [crates.io setup requirements](https://crates.io/docs/trusted-publishing).
Trusted Publishing requires an already-published crate and crate ownership;
initial publication uses an API token. `la-stack` already exists, so no throwaway
version or initial upload is needed. As its owner, verify this exact identity
under the crate's Settings → Trusted Publishing and register it if missing:

| Field | Value |
| --- | --- |
| Repository owner | `acgetchell` |
| Repository name | `la-stack` |
| Workflow filename | `publish-crates.yml` |
| Environment | `crates-io` |

For initial setup, create GitHub's `crates-io` environment with required maintainer reviewers,
administrator bypass disabled, and selected deployment tags matching `v*` only
(no branches). The `.github/settings/crates-io-environment.json` payload names
`acgetchell` (user ID 434849) as the required reviewer and permits that maintainer
to approve their own deployment; self-review prevention needs another eligible
reviewer. `.github/settings/crates-io-tags.json` declares the tag restriction.
The stable-tag and reviewed-commit
gates still apply after deployment approval. Add only the pinned
`rust-lang/crates-io-auth-action` entry from `.github/settings/actions-selected.json`
to the live selected-action allowlist, preserving newer entries. This reviewed
v1.0.5 pin comes from the shared consumer template; update it through a reviewed
dependency PR and reconcile the allowlist together. Both the publishing workflow
and the CI workflow that signs reusable evidence are excluded from Dependabot
automatic approval so these changes require maintainer
review. Do not add a long-lived
`CARGO_REGISTRY_TOKEN` secret.

Workflow YAML and settings payloads do not apply account configuration. The initial
Trusted Publishing rollout and verified v0.5.0 publication are recorded in
[#259](https://github.com/acgetchell/la-stack/issues/259). Before each release,
verify the actual crates.io publisher identity and GitHub environment reviewers/tag
policy, allowlist and permissions. Record that evidence in the current release's
tracking issue or release PR. Re-read live settings before changing them.
Validate the source with `just lint-config`, `just python-ci`, `just markdown-ci`
and `just doc-check`; PR tests use simulated providers and never upload crates.
Native CI and successful approved publication remain separate rollout evidence.
