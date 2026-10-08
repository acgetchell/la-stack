# Managing Git and GitHub Changes

Operational details for the Git rules in [AGENTS.md](../../AGENTS.md).

## Contents

- [Git operations and branch names](#git-operations-and-branch-names)
- [Commit messages](#commit-messages)
- [GitHub CLI](#github-cli)
- [Issue planning](#issue-planning)
- [Issue dependencies](#issue-dependencies)
- [Dependabot approval rollout](#dependabot-approval-rollout)

## Git operations and branch names

Agents use `git --no-pager` for read-only Git commands. Mutating refs or the
index requires explicit user approval covering the intended operations, such
as branch creation, staging, commits, pushes, or tags. Approval remains valid
for the authorized task; do not request it again unless the scope changes.
Without approval, suggest commands for the maintainer to run manually.
Preserve unrelated user changes.

Prefer branch names of the form `{type}/{issue}-descriptor-or-two`, for example
`fix/307-topology-validation`, `perf/315-bench-profile`, or
`doc/329-branch-guidance`. If an environment requires an owner/tool prefix,
retain that structure after it, for example
`codex/fix/307-topology-validation`.

## Commit messages

When asked to generate a commit message:

1. Run `git --no-pager diff --cached --stat` and inspect the staged diff.
2. Use `<type>: <brief summary>` with `feat`, `fix`, `refactor`, `perf`,
   `docs`, `test`, `chore`, `style`, `ci`, or `build`.
3. Include organized body bullets describing the changes and test results.
4. Unless committing is explicitly authorized, present the message in a code
   block with no language so the user can commit manually.

Document intentional API breaks explicitly. See the
[contributor commit-message guide](../../CONTRIBUTING.md#commit-message-format)
for Conventional Commit examples and breaking-change markers. Changelog
generation and release procedures belong in [Releasing](../RELEASING.md).

## GitHub CLI

Use structured output and `| cat` when reading GitHub objects to avoid pagers
and scope errors:

```bash
gh issue view 64 --repo acgetchell/la-stack --json title,body | cat
gh issue view 64 --repo acgetchell/la-stack --json title,body \
  --jq '.title + "\n" + .body' | cat
gh issue list --repo acgetchell/la-stack --json number,title,labels \
  --jq '.[] | "#\(.number) \(.title)"' | cat
```

Avoid plain `gh issue view N`, which may open a pager or fail with
`read:project` scope errors. Include `labels` and `milestone` in `--json`
when inspecting issue placement; use `--label` and `--milestone` to filter lists.

For arbitrary Markdown in issue bodies or comments, use `--body-file` with a
file or a quoted heredoc. For example:

```bash
gh issue comment 64 --repo acgetchell/la-stack --body-file - <<'EOF'
## Summary

Body with `backticks`, **bold**, and apostrophes that's safe.
EOF
```

Use `gh issue create`, `gh issue edit`, `gh issue comment`, and
`gh issue close` for the requested operation. Creation supports `--title`,
`--body-file`, `--label`, and `--milestone`; edits support those metadata
changes and `--add-label`. Close with the appropriate reason, `completed` or
`not planned`.

## Issue planning

- Use appropriate existing labels, such as `enhancement`, `bug`,
  `performance`, `documentation`, `rust`, or `python`.
- Assign the appropriate milestone and preserve the maintainer's requested
  release placement.
- Structure new issue bodies around Summary, Current State, Proposed Changes,
  Benefits, and Implementation Notes.
- Cross-reference related work with `#XXX`, or `owner/repo#XXX` across
  repositories. Distinguish actual prerequisites from related work.
- Record dependency intent clearly in prose, such as `Depends on: #XXX`,
  `Blocks: #YYY`, or `Related: #ZZZ`, and use native metadata for blocking
  relationships.
- Verify issue contents and metadata after creating or updating them.

## Issue dependencies

Create blocking relationships through GitHub's native dependency metadata.
Prose cross-references document intent; verify the relationship itself through
the [issue-dependency API](https://docs.github.com/en/rest/issues/issue-dependencies).
The API takes the blocking issue's internal integer ID, not its issue number.

For example, to make issue 217 blocked by issue 207:

```bash
# Inspect the blocker and existing relationships first.
gh api repos/acgetchell/la-stack/issues/207 --jq '.id' | cat
gh api repos/acgetchell/la-stack/issues/217/dependencies/blocked_by \
  --jq '[.[].number]' | cat

# Replace BLOCKING_ISSUE_ID with the returned integer before running.
gh api repos/acgetchell/la-stack/issues/217/dependencies/blocked_by \
  -X POST -F issue_id=BLOCKING_ISSUE_ID | cat

# Verify the relationship after adding it.
gh api repos/acgetchell/la-stack/issues/217/dependencies/blocked_by \
  --jq '[.[].number]' | cat
```

Use `-F` so the resolved numeric ID is encoded as an integer. Keep independent
tasks related without inventing a blocking dependency.

## Dependabot approval rollout

The caller in `.github/workflows/dependabot-auto-merge.yml` adopts the same
reviewed research-repo-tools revision as markov-chain-monte-carlo:
[`cbb2ea6`](https://github.com/acgetchell/research-repo-tools/commit/cbb2ea6dee8866b3f0547bca935aef48fdd71707).
It replaces the owner-issued CodeRabbit request that left bot PRs without an
approving review (#258). A successful CodeRabbit status on a skipped review does
not supply approval; the provider's underlying failure to process owner commands
remains unconfirmed. The new workflow removes that approval dependency.

The trusted base branch supplies exact ecosystem/file allowlists, including the
comparison workspace manifest. Signed Dependabot metadata, verified commits and
base merges, complete file inventories, current-head checks, and approval
deduplication belong to the shared workflow. It never checks out PR code or uses
personal tokens. `.github/dependabot.yml` continues to own dependency eligibility.
See the [shared approval contract][shared-approval].

[shared-approval]: https://github.com/acgetchell/research-repo-tools/blob/cbb2ea6dee8866b3f0547bca935aef48fdd71707/docs/AUTOMATING_DEPENDABOT.md

The settings payloads under `.github/settings/` capture the October 1, 2026
rollout. They enable Actions approvals with read-only default permissions, add
only the metadata action and reusable workflow to the existing selected-action
allowlist, and enable stale-review dismissal in ruleset `11113227`. All six
required checks, one required approval, resolved threads, strict checks, merge
methods, and the existing bypass actor are preserved. Auto-merge and squash
merging are already enabled. Actions approval is an intentional change from the
old #258 criterion, matching the maintainer's requested MCMC approach.

Before applying payloads later, reread current settings and reconcile intervening
changes; do not overwrite newer checks or allowlist entries with this snapshot.
The settings can then be applied with:

```bash
gh api repos/acgetchell/la-stack/actions/permissions/workflow \
  --method PUT --input .github/settings/actions-workflow.json
gh api repos/acgetchell/la-stack/actions/permissions/selected-actions \
  --method PUT --input .github/settings/actions-selected.json
gh api repos/acgetchell/la-stack/rulesets/11113227 \
  --method PUT --input .github/settings/main-ruleset.json
```

Merge the caller to `main` before expecting `pull_request_target` to use it.
An old workflow rerun retains its original definition; exercise a new eligible
Dependabot event after deployment. Verify an actual `github-actions[bot]`
approval whose commit ID is the current PR head, native squash merge only after
required checks and threads pass, and continued manual review for ineligible
updates. Record PR, head, approval, workflow-run, and merge links. Installing the
Python package does not deploy or configure this GitHub workflow.

Merges using `GITHUB_TOKEN` may not trigger push workflows. Dispatch `ci.yml`
against `main` and verify the run's `headSha` equals the merged commit before
claiming default-branch validation. Remove `CODERABBIT_REVIEW_TOKEN` from Actions
and Dependabot secrets only after the old caller has been replaced and no other
workflow consumes it. Local opt-in CodeRabbit CLI review remains separate.
