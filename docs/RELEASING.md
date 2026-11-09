# Releasing la-stack

Prepare each `vX.Y.Z` release in a dedicated PR. Acquire the comparison once
locally during Step 1, keeping original assertions outside the checkout.
After merge and native CI, verify that acquisition against the final clean
commit before tagging. The tag push prepares a draft automatically. Stage accepted records and dispatch
`Release Benchmarks` to verify and attach the original archive. Publish the
reviewed draft, explicitly dispatch `Publish to crates.io` with the independent
publication selection, then approve its `crates-io` deployment.

Benchmark attachment and Cargo publication require the maintainer's authenticated
dispatch against the exact release tag. Preparation and both dispatches reuse evidence without measurements.
Development supplier pins support
rehearsal and ordinary checks; actual release work requires published tooling.

The changelog is generated for the target tag before the tag exists, so the
release process does not require a temporary local tag.

Release recipes are content-idempotent where their inputs permit it. Repeating
`update-version` on the same UTC date or repeating `performance-readme` for the
same retained bundle produces no file changes. Rerunning `update-version` on a
later UTC date intentionally advances the declared release date. Crossing
midnight alone requires no update: the `declared` date policy allows preparation
and tagging on different UTC days when citation and changelog dates agree.
Fresh acquisition requires a new checkpoint. Use `performance-resume` for an
interrupted campaign; completed suite pairs retain their original bytes.
With unchanged history and configuration, repeating `changelog-release`
reuses the declared citation date and preserves the generated contents,
including when a rerun happens after UTC midnight.

## Contents

- [Reusable CI evidence](#reusable-ci-evidence)
- [Release command alignment](#release-command-alignment)
- [Prepare the environment](#prepare-the-environment)
- [Step 1: Prepare the release PR](#step-1-prepare-the-release-pr)
- [Step 2: Publish after the PR is merged](#step-2-publish-after-the-pr-is-merged)
- [Configure or verify Trusted Publishing](#configure-or-verify-trusted-publishing)

## Reusable CI evidence

The publishing process uses the pinned shared tooling's release evidence
commands. `tooling/release-evidence.toml` declares the repository,
package, workflow, native jobs/steps and signer policy. The guarantee is
`signed-source`; it does not promise promotion of an identical Cargo archive.
Protected `main` pushes produce signed evidence only with published tooling,
after all three native `just ci` jobs and the default/exact Cargo packaging job
succeed. A development Git pin passes ordinary validation but cannot produce
release evidence. CI verifies its development identity and skips signing and
uploading the release-evidence artifact; ordinary validation still runs.
The `release-evidence` job signs their exact commit, run, attempt, job IDs and
package identity, then retains the manifest and Sigstore bundle in an immutable
Actions artifact for 30 days.

`just release-evidence` automatically resolves the original Step 1 owner and verifies original acquisition,
the precise publisher receipt, final clean tree and signed native evidence.
It requires one successful `push` run on `main`, every required job and step
from that attempt, and a nonexpired immutable artifact. The shared verifier
checks the ZIP digest, exact inventory, signature, workflow, source commit,
hosted runner and package identity. PR, manual-dispatch, fork and partial-rerun
results cannot replace main-push evidence.

Retain the owner directory printed by `performance-release` and an independent
backup under exclusive maintainer control. Its default location is
`~/.local/state/la-stack/acquisitions/<tag>/<checkpoint-path-sha256>`. Step 1
also freezes `owner.json` under `~/.local/state/la-stack/releases/<checkout-sha256>/<tag>/`,
binding the original owner and acquisition/publication digests. Routine Step 2
commands read that record and derive the full commit from the clean checkout;
they never search for the newest owner or reconstruct ownership from candidate files.
Retain both directories and an independent backup. An ignored
checkout directory is not an authentication boundary. Original raw evidence,
preparation, acquisition assertions, publication receipt and their digests
are immutable.

Final acceptance lives under `<owner-directory>/final/<commit>/<selection>/`.
The default selection is `reviewed`. Repeating verification revalidates the
full acceptance, including the original native run, attempt and artifact.
To explicitly review a newer complete CI attempt, use a new name such as
`just release-evidence "$OWNER" retry-2` and pass that same selection to
archive, staging and publication recipes before GitHub publication. The command
prints the full acceptance file's SHA-256. After GitHub publication, renew only
that acceptance and supply its reviewed digest to a fresh Cargo dispatch;
the original release assets and publication selection remain unchanged.
Existing records are never replaced.

These recipes launch the pinned supplier through `uvx` outside the observed
checkout. Step 2 does not measure, discover, rebuild or mutate Git. It writes
only requested receipts outside the checkout. `release-native-evidence`
remains a focused native-CI diagnostic; it does not establish candidate reuse.

The independently retained `transfer.json` selects the exact draft release
and four inert asset IDs/hashes. After attachment, `publication-transfer.json`
also selects the canonical archive ID and original hash. These operator-supplied
dispatch inputs authenticate the assertions; hashes from received candidate
files cannot authenticate those same files. The Cargo dispatch also requires
the independently reviewed full acceptance digest. Verification freezes those
exact bytes in an immutable Actions artifact before approval; the upload job
revalidates that receipt after approval. Both jobs check the original asset
selection and immutable release state; generated notes and registry absence
remain publication gates.

This reuses validation of the source. Cargo still packages the clean release
commit when uploading, so it does not promise byte-for-byte promotion of the
dry-run `.crate` archive. The shared tooling's Python distribution promotion
process has a different artifact contract.

Use a current GitHub CLI with `gh attestation verify` and the signer/source
constraints used by the installed verifier. See GitHub's
[attestation verification options](https://cli.github.com/manual/gh_attestation_verify).

## Release command alignment

The common Step 1/Step 2 flow follows the
[research-repo-tools v0.2.6 guide](https://github.com/acgetchell/research-repo-tools/blob/eb7e3eaf158024a045a6711898a8af199729131b/docs/RELEASING.md).
Draft preparation uses the public capability from upstream
[PR #207](https://github.com/acgetchell/research-repo-tools/pull/207), rehearsed here at
`cd8594f658b5bf3ade063a685834565cc62e0cb3`, including receipt retrieval from
[PR #209](https://github.com/acgetchell/research-repo-tools/pull/209). Published-package and final native adoption
remain required by [#333](https://github.com/acgetchell/la-stack/issues/333).

| Command | Common behavior and la-stack declaration |
| --- | --- |
| `update-version TAG [ARGS...]` | Prepare package/citation metadata; Cargo owns the Rust lock, the private Python environment retains its version |
| `changelog-release TAG [DATE]` | Generate prospective notes using the prepared citation date |
| `release-version-check` | Require synchronized final metadata, changelog and citation dates |
| `publish-check` | Verify clean Cargo source with default/exact dry-runs; Cargo requires builds here, unlike the reference's existing-distribution Twine check |
| `release-evidence [OWNER] [SELECTION]` | Inspect retained signed main CI and original benchmark acceptance without builds; defaults recover the original owner and `reviewed` selection |
| `tag TAG` | Run the final version check and create the annotated tag from generated notes |
| `release-notes TAG` | Print generated notes; tag-triggered preparation supplies these to the draft automatically |
| `release-benchmarks TAG` | Dispatch retained archive attachment from the original owner record; explicit `CANDIDATE_TRANSFER_FILE` remains a recovery override |
| `release-publish TAG [OWNER] [SELECTION]` | Revalidate acceptance, notes and canonical assets, then explicitly approve GitHub publication |
| `release-upload TAG [OWNER] [SELECTION]` | Dispatch Cargo publication using the existing reviewed acceptance and original publication transfer |
| `release-verify TAG` | Verify stable GitHub assets and the exact crates.io version |

This repository verifies evidence **before** creating or pushing the tag. The
reference's Python distribution inspection requires a pushed tag; that ordering
difference does not weaken la-stack's pre-tag gate. Preparation independently
revalidates its frozen native selection before creating an unpublished draft.
Benchmark acquisition and README generation belong only to Step 1. Cargo's
nonbuilding repackaging/upload, explicit transfer selection and `crates-io`
approval are the remaining registry and scientific-asset differences.

## Prepare the environment

Replace any development tooling pin with its supplying published release and
reconcile the manifest, lock and isolated runtime pins as described in
[Shared capability adoption](../tooling/README.md#shared-capability-adoption).
Before release acquisition, run `just release-dependencies`. This shared gate
rejects Git sources even when their version label matches a published package;
release evidence and publication repeat the check. Retain development results
as rehearsal evidence and revalidate the final published-pin inputs.
The root `dependency-sources.toml` declares the Python tooling project and the
Cargo workspace's release-owned packages. Keep the private comparison member
explicit; other dependency packages require verified registry artifacts. Include
this policy with manifests and locks in source and benchmark evidence.

Set the target tag once:

```bash
TAG=vX.Y.Z
```

Verify the repository remotes and synchronize `main`:

```bash
gh auth status
git --no-pager remote -v
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
just changelog-release "$TAG"
```

This generates `CHANGELOG.md` as though the target tag already existed, archives
completed minor series under `docs/archives/changelog/`, and defaults to the
`date-released` already recorded in `CITATION.cff` by `update-version`.
No `RELEASE_DATE` shell variable is needed, and generation does not read today's
date or update package or citation metadata. The citation version must match
the target tag, and its date must be a valid ISO `YYYY-MM-DD` date; missing or
invalid metadata stops generation with preparation guidance. Reruns reuse the
recorded date even after UTC midnight and preserve existing release dates.

The optional explicit date remains supported. For troubleshooting, populate
`RELEASE_DATE` manually and pass it as the second argument:

```bash
RELEASE_DATE="$(date -u +%F)"
just changelog-release "$TAG" "$RELEASE_DATE"
```

This override explicitly chooses today's UTC date; keep it consistent with
`CITATION.cff`. To intentionally change the declared release date, update the
metadata first with `just update-version "$TAG" --date "$RELEASE_DATE"`.
You can also set `RELEASE_DATE` to a reviewed `YYYY-MM-DD` value instead of today.

Use `just changelog-preview --tag "$TAG"` for a read-only
preview, and `just changelog-check` to validate the root and every archive.
Tagged previews use the same citation-date default; `--date YYYY-MM-DD` remains
available as an explicit override.
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

After updating the package version, run:

```bash
just performance-release
```

The recipe fetches tags, prepares the independent source assertions, then
measures and promotes the release report. Preparation selects the release pair
using GitHub's publication history, resolves the current and baseline commits,
and saves all four assertions to
`target/bench-reports/release-inputs.json`, relative to the repository root.
This Git-ignored file carries the tags and commits forward to
`just performance-readme`; it survives closing and reopening your terminal.
The selected tags and commits are printed before measurement. A failed fetch,
release lookup or preparation stops the recipe before benchmarking.
README publication remains a separate, retryable step below.

The optional progress arguments are `just performance-release human 300`
(the default) or, for captured JSON status, `just performance-release json 60`.
These commands acquire measurements. Shared status reports the active native
stage, validated counts, elapsed time and an ETA only when its observations
support one. Five-minute heartbeats do not limit native command runtime.
The [benchmark guide](BENCHMARKING.md#release-acquisition-policy) explains
unknown budgets/estimates and the pinned supplier's gate-failure event limitation.
The shell runtime record still retains the final duration and exit status.
After acquisition succeeds, retry report promotion with `just performance-promote`
and README publication separately; neither retry performs measurements.

For preparation alone, `just performance-prepare` uses local tags without
changing refs or running benchmarks. Add `--fetch-tags --allow-git-mutations`
to refresh tags explicitly. This separate command remains useful for inspecting
inputs before acquisition. Restore a lost original preparation from its
independent owner record; never derive replacement assertions from a candidate.
The shared preparation record also binds the package version and measurement
declaration. If you have an old adapter-format `release-inputs.json`, remove
that ignored scratch file once and prepare again. It is deliberately not
converted from measurement evidence. Moved refs or changed package/declaration
inputs require a new preparation before measurement or README publication.
Fresh acquisition refuses an existing checkpoint before preparation.
Resume explicitly reuses the original preparation and checkpoint:

```bash
just performance-resume
```

The default checkpoint is `target/bench-reports/full/checkpoint`; pass another
path as the first resume argument when needed. One sequential trial retains
completed suite pairs. An interrupted pair is recollected in full; completed
pairs, original sessions and the append-only journal remain intact. Changed
inputs, incompatible declarations or ambiguous work stop recovery.

Before report promotion, the adapter freezes the independent owner record.
After a capture failure, preserve scratch and retry
`just performance-capture "$OWNER"`. If raw bytes did not reach that record,
explicitly use `just performance-capture "$OWNER" --retained-scratch` while
you still control the original acquisition. Then retry promotion and README
rendering directly. Never treat incoming candidate files as owner assertions.

No shell variables or exports are needed. `performance-prepare` does not set
`PREVIOUS_TAG`, `CURRENT_COMMIT`, or `PREVIOUS_COMMIT` in your shell; those
variables may remain unset. Do not run the old manual `git rev-parse`
assignments or pass those variables to the ordinary README step below.

Review `docs/performance.md` and the complete
new run and index changes under `docs/performance-runs/`. Its current/baseline
tags and commits must match the selected assertions. Both labels must be v0.5.0
or newer.

Evidence retains every available phase, mean and median estimates, 100 raw
samples per case and 95% intervals. Availability uses the original revision:
APIs introduced in v0.5.2 are unavailable earlier and paired on v0.5.2-or-later
baselines. No supported missing case is treated as unavailable. All 1,403
library cases and 16 LU peer cases on each phase remain required.

The full schema-6 plan declares one sequential trial. Baseline inputs include
README because historical releases compile it. Current source includes README
only during rustdoc; current gates inspect debug and release dependency files
and reject a compiler dependency on README before timing. A historical current
revision that still compiles README cannot use this full candidate plan.
Use explicit non-exact diagnostics when appropriate.

Complete dependency pins, source, harness, declarations and release metadata
before acquisition. Later changes outside the precise output contract invalidate
candidate reuse. Reports and README select retained trial 1; no extra trial or
counterbalancing claim is implied.

Historical artifacts, including the v0.5.0/v0.4.6 comparison, remain frozen with
their original bytes and links. That older run lacks the new acquisition
declaration and is no longer a regeneration input. Never backfill provenance
or relax current eligibility checks to replay it.

### 5. Refresh the README benchmark comparison

```bash
just performance-readme
```

With no arguments, this reads the independent tags and commits from
`target/bench-reports/release-inputs.json`, saved during `performance-release` preparation.
It consumes archived, release-eligible evidence without benchmarking again or
requiring shell variables.
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

The four explicit tag/commit arguments remain supported as described in the
[benchmark guide](BENCHMARKING.md#update-the-readme-nalgebrafaer-table).
If saved inputs are missing or invalid, restore the original preparation from
the independent owner record or supply the four independently reviewed
assertions, then retry `performance-readme`. Do not start a new acquisition
solely to retry publication. A README failure leaves retained measurements
available for another publication attempt.

See [Benchmarking](BENCHMARKING.md) for explicit scratch versus retained plot
inputs, named metric/statistic/axis views, output ownership and methodology.

### 6. Validate the release branch

```bash
just ci
just release-version-check
```

`just ci` includes the lockfile, citation, documentation, test, and benchmark
compile checks.
The final version check requires the current changelog heading and citation date
to agree. After committing the reviewed artifacts below, run `just publish-check`
on that clean release branch. It performs native default/exact Cargo dry-runs,
selecting only `la-stack`. Main CI repeats this on the final merged commit;
those hosted results become the reusable packaging evidence. Step 2 reuses them
and never repeats package-verification builds.

### 7. Review, stage, and commit the release artifacts

Inspect all changes before staging:

```bash
git --no-pager status --short
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

just publish-check
```

### 8. Push the branch and open the PR

```bash
git push -u origin "release/$TAG"
```

Use `chore(release): release $TAG` as the PR title and describe the PR as a
focused release preparation without feature work.

### Handling fixes found during preparation

For a critical fix that must be included, make and commit the fix, rerun
`just changelog-release "$TAG"`, review and stage only the resulting changelog
files, and commit that generated update separately.

For a non-critical fix, file an issue and defer it to a later release. Do not
hand-edit the generated changelog to add a known-issue note.

## Step 2: Publish after the PR is merged

### 1. Synchronize `main` and wait for CI

```bash
git switch main
git pull --ff-only
test -z "$(git --no-pager status --porcelain)"
```

Use a clean checkout of the reviewed release commit. If `main` advanced, review
the intervening changes and refresh generated notes through a PR as needed.
Wait for that commit's complete main-push CI run. The next step derives the commit
and checks eligible evidence automatically; pending results stop it with a retry
diagnostic. Use GitHub's Actions page or the
[manual inspection commands](#manual-context-and-inspection) when troubleshooting.

Require the complete run: all three native `build` jobs, default/exact Cargo
packaging in `release-package`, and signed `release-evidence`. A development
supplier pin skips signing and cannot provide eligible release evidence.
Inspect the independent SARIF checks for this exact commit too:

```bash
gh run list --workflow rust-clippy.yml --event push --branch main --commit "$(git --no-pager rev-parse HEAD)" --limit 20
gh run list --workflow semgrep-sarif.yml --event push --branch main --commit "$(git --no-pager rev-parse HEAD)" --limit 20
```

Wait for successful `Clippy SARIF Analysis` and `Repository Rule SARIF Analysis`
from GitHub Actions (app ID `15368`). They are separate publication gates;
`release-evidence` verifies only the jobs and signer declared in its evidence
policy, not every independent workflow. PR review remains a protected-branch
requirement and cannot replace merged main-push evidence.

### 2. Verify the final release commit

The normal commands need no `OWNER` or `RELEASE_COMMIT` shell variables:

```bash
just release-evidence
just release-archive
```

Require native Linux, macOS and Windows CI, default/exact Cargo packaging and
signed main-push evidence, plus the independently checked SARIF results. The candidate verifier
also checks that the final tree differs from the measured candidate only through
the accepted publication contract. It preserves the original measured revision
separately from the final commit. Do not tag until verification succeeds.

If main-push evidence needs renewal, rerun every job in the original run and
explicitly review a new selection. Manual CI dispatch cannot supply release
evidence. Pending, failed, cancelled or expired evidence blocks publication.

### 3. Create and verify the annotated tag

```bash
just tag "$TAG"
git --no-pager tag -l --format='%(contents)' "$TAG"
test "$(git --no-pager rev-parse "$TAG^{commit}")" = "$(git --no-pager rev-parse HEAD)"
```

`just tag` first runs `release-version-check` and revalidates the existing reviewed
acceptance for clean HEAD. A different commit or missing receipt requires an
explicit `release-evidence` review before tagging. It then builds the annotation from the matching active or archived changelog
section. For a changelog larger than 125 kB, it links that section.

### 4. Push the tag

Only after evidence verification and tag inspection succeed, push the tag.
This starts `Prepare GitHub release`; pushing does not wait for or verify CI.

```bash
git push origin "$TAG"
```

### 5. Wait for the draft GitHub release

The tag push starts `prepare-release.yml`. Let its original attempt finish;
the next `release-stage` command automatically retrieves and verifies that exact
preparation receipt. Missing or pending preparation stops staging; wait for the
original run and retry. The manual inspection commands below can show that run.
Let preparation create the draft; manual creation is unnecessary.
The read-only job checks the exact tagged source, main ancestry, final metadata
and retained signed native evidence, then freezes that selection and generated
notes. An isolated pinned writer with no checkout authenticates the handoff by
immutable artifact ID, digest, run and attempt, revalidates the same native
selection and creates the unpublished draft through the shared API.

Retain the handoff and separate draft-receipt artifact IDs, ZIP digests and
producer run/attempts in independent release review notes. Preserve the captured
release ID in `target.json`. An existing exact matching draft is unchanged;
conflicting notes, title, lifecycle, moved tags or missing evidence stop work.
No tests, package builds, signing or benchmark acquisition run during preparation.

Keep the release title identical to the tag, including its leading `v`.
Leave it as a stable draft: publishing now would make the release immutable
before the benchmark archive can be attached. GitHub recommends this
[draft-first sequence for immutable releases](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases).

### 6. Attach accepted benchmarks and review the draft

Stage the four inert records:

```bash
just release-stage
```

Staging uses the installed `release assets retrieve` operation and
`tooling/draft-receipt.toml`. It authenticates the exact tag-push workflow source,
original successful attempt, required jobs/steps, immutable artifact ID and ZIP
digest before retaining `selection.json`, `target.json` and `draft-receipt.json`
under the owner's final record. Retries revalidate that frozen provider selection
instead of searching again. A replacement draft is rejected before any upload.

This revalidates acceptance, freezes the prepared draft ID before any upload, and prints
the independent `transfer.json` path. A retry uses that same draft and identical
bytes. Dispatch from the exact tag using the printed path:

```bash
just release-benchmarks "$TAG"
```

The recipe resolves the existing transfer from the original owner record and
clean HEAD, validates the explicitly supplied tag, and returns after dispatch.
It never substitutes another owner or tag, polls, creates a draft,
acquires measurements or runs `performance-readme`. Inspect the matching
`Release Benchmarks` run and require success before publishing.

The no-checkout retrieval job checks original release/asset identities and
downloads the selected bytes. A separate read-only checkout job verifies
acquisition, final binding and full frozen native acceptance. The no-checkout
attachment job receives an immutable artifact by exact ID, digest and current
attempt, rechecks the draft and original assertions, then attaches
`la-stack-$TAG-criterion-baseline.tar.gz`. No job measures or builds.

The archive retains original evidence, current Criterion files and a separate
final-release binding. Local ownership, measured revisions, hosts, sessions and
timestamps remain unchanged; acceptance does not relabel their source.

Review the draft title, generated notes and archive. The workflow leaves it
unpublished. Do not edit, replace or publish it while a job is running. Runs for
the same tag serialize without cancellation. Privileged jobs use isolated
pinned tooling and execute no checkout code.

### 7. Publish the reviewed draft and approve Cargo publication

From the clean tagged checkout and your human GitHub CLI session:

```bash
just release-publish "$TAG"
```

This verifies acceptance and remote notes, freezes the canonical attached archive
ID/hash in `publication-transfer.json`, then publishes the reviewed draft
through the shared gate. It authorizes GitHub Release publication. The release
must become immutable before Cargo publication.

Dispatch Cargo publication from that independently retained selection:

```bash
just release-upload "$TAG"
gh run list --workflow publish-crates.yml --limit 5
```

`release-upload` revalidates the already existing acceptance and original
publication transfer, checks the immutable release, and sends a snapshot of the
same assertions and acceptance digest. Missing review or changed evidence stops
dispatch; this command never creates a review selection or chooses a newer run.
It returns after dispatch and does not approve the deployment.

The first hosted job checks immutable release state, canonical and transport asset
IDs/hashes, original archive acceptance, final clean commit and generated notes.
It verifies current native evidence against the separately reviewed full
acceptance digest and freezes the exact receipt in an immutable Actions artifact.
After it succeeds, approve the `crates-io` deployment in GitHub. The upload job
retrieves that artifact by exact ID, digest and current run/attempt, repeats
verification against its frozen receipt, proves the registry version
absent, then requests an OIDC token. Only that job has `id-token: write`; it
has no repository write permission or dependency caches. Uploads serialize
without cancellation.

Original main CI built default and `exact` packages without upload credentials.
Cargo repackages the clean commit with `--locked --no-verify`, so build scripts
do not run with the crates.io token. This is signed source validation, not
byte-identical Cargo archive promotion. No benchmark workspace member is published.

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

Hosted jobs restore no dependency caches. The checkout verifier has read-only
repository permissions. Draft retrieval and attachment require write access
but execute no checkout code. See
[hosted verification](BENCHMARKING.md#hosted-release-verification) and
[workflow validation](BENCHMARKING.md#validate-the-release-workflow) for the
exact handoff, local regressions and separate hosted evidence requirements.

### Manual context and inspection

The normal commands recover context after restarting the shell. Explicit paths,
selections and GitHub commands remain available for older acquisitions, a moved
checkout or deliberate recovery. They do not replace or rewrite the normal
external owner selection. Do not infer missing assertions from received assets.

```bash
RELEASE_COMMIT=$(git --no-pager rev-parse HEAD)
OWNER="/path/to/original/independent/acquisition"
just release-evidence "$OWNER"
just release-archive "$OWNER"
```

The same explicit owner and optional selection arguments are supported by
`tag`, `tag-force`, `release-publish` and `release-upload` during recovery.

Use the recorded workflow/run or inspect the exact current commit:

```bash
gh run list --workflow ci.yml --event push --branch main --commit "$RELEASE_COMMIT" --limit 20
gh run list --workflow prepare-release.yml --event push --branch "$TAG" --commit "$RELEASE_COMMIT" --limit 20
gh run watch ORIGINAL_RUN_ID --exit-status
```

For manually reviewed preparation retrieval, the shared CLI preserves the same
source, original-attempt, job and artifact checks:

```bash
just release-assets retrieve tooling/draft-receipt.toml --tag "$TAG" --commit "$RELEASE_COMMIT" \
  --run-id ORIGINAL_PREPARATION_RUN_ID --output /path/to/independent/preparation
just release-stage "$OWNER" /path/to/independent/preparation/target.json
export CANDIDATE_TRANSFER_FILE="$OWNER/final/$RELEASE_COMMIT/reviewed/transfer.json"
just release-benchmarks "$TAG"
just release-publish "$TAG" "$OWNER"
just release-upload "$TAG" "$OWNER"
```

For direct dispatch troubleshooting, `publish-crates.yml` still accepts
`--ref "$TAG" -f tag="$TAG" -F candidate=@/original/publication-transfer.json`
and `-f acceptance=REVIEWED_FULL_ACCEPTANCE_SHA256`. The explicit selection must
be the same original owner record; do not copy a digest from incoming assets.

### Recovering a failed run

If HEAD changes, return to Step 2's synchronization and exact-commit verification.
If a native run/attempt changes, explicitly review a new named selection as
described below. Missing, expired, PR, dispatch or partial-rerun evidence never
authorizes fallback to a replacement run or a local rebuild.

#### Draft preparation recovery

Automatic preparation accepts only its original attempt and rejects prior runs
for the same tag. Do not rerun a failed preparation or recreate its tag to obtain
a fresh receipt. Preserve the original authenticated handoff and journal on
success and failure. A pending journal means creation may have happened; a
captured journal binds the existing release ID. A missing or edited draft is not
permission to create a replacement.

Recover the exact original artifacts by independently recorded immutable IDs,
verify their ZIP digests and producer run/attempts, and retain the original policy,
selection, notes and receipt. Follow the supplying package's
[recovery](https://github.com/acgetchell/research-repo-tools/blob/cd8594f658b5bf3ade063a685834565cc62e0cb3/docs/PREPARING_DRAFTS.md#retry-and-recovery).
After independent provider/audit review, its `release assets prepare` CLI accepts
`--recover-release-id REVIEWED_ID` with the original pending receipt and frozen
inputs; recovery performs no POST. A captured ID cannot be replaced. This is an
explicit maintainer recovery, not another automatic preparation attempt.

The workflow retains records for 30 days; hard runner termination or upload
failure may prevent journal retention. GitHub history guards are not a permanent
intent ledger. Keep workflow identity/history intact for retriggerable tags. If
history or the journal is unavailable, stop automation and reconcile intent and
provider state manually before authorizing any further creation. Recovery writers
must use the same repository/tag serialization group as preparation and attachment.

#### Benchmark and publication recovery

Resume local acquisition with `just performance-resume`; keep the original
checkpoint, preparation and owner record. Retry capture, report promotion and
README publication separately after acquisition completes. Source/harness drift
requires a new reviewed acquisition, not repaired receipts.

For native evidence, rerun **all** jobs in the original main-push run. Partial
attempts and manual CI dispatch remain ineligible. Accept a new complete native
attempt through an explicit new review selection:

```bash
just release-evidence "$OWNER" renewed-ci
```

Before GitHub publication, use the new selection for archive, staging and
attachment. After GitHub publication, keep the original immutable archive,
transport assets and `publication-transfer.json`. Review the new acceptance and
dispatch `publish-crates.yml` again with its full SHA-256 as `acceptance` and
the original publication selection as `candidate`. This renews CI evidence
without reacquiring timings or replacing release assets. A changed native
attempt after deployment approval requires a fresh reviewed dispatch.

For hosted recovery, rerun all jobs or dispatch again with the same independent
selection. Inert artifacts are selected by exact ID and attempt and retained for
30 days. Failed-jobs-only attachment reruns cannot reuse another attempt's
verification. These retries perform no timings.

Upload retries reuse the original draft ID and identical bytes. Conflicting
or incomplete assets stop without overwriting or deletion. Inspect failures
before manual correction. A replacement draft requires explicitly authorized
preparation recovery and an authenticated new target; a new review selection
alone does not authorize replacement. Never move a published tag or replace independent digests with values
from received files.

A rerun uses the tagged workflow. A later workflow fix requires a reviewed new
candidate. If Cargo reports an ambiguous timeout, inspect exact registry
visibility before retrying; a present version must not be uploaded again.
`just release-verify "$TAG"` checks the release and registry with bounded
retries. Keep protected deployment approval and all pre-upload checks enabled.

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
