# Managing Git and GitHub Changes

Operational details for the Git rules in [AGENTS.md](../../AGENTS.md).

## Contents

- [Git operations and branch names](#git-operations-and-branch-names)
- [Commit messages](#commit-messages)
- [GitHub CLI](#github-cli)
- [Issue planning](#issue-planning)
- [Issue dependencies](#issue-dependencies)
- [Dependabot approval configuration](#dependabot-approval-configuration)

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

## Dependabot approval configuration

The [caller workflow](../../.github/workflows/dependabot-auto-merge.yml) owns the
reviewed research-repo-tools workflow pin and exact ecosystem/file allowlists.
The shared workflow verifies signed Dependabot metadata, commits, base merges,
complete file inventories, and the current head before approving. It never
checks out PR code or uses personal tokens.
[Dependency configuration](../../.github/dependabot.yml) owns update eligibility.

GitHub must allow Actions to approve PRs and permit the caller's pinned action
and reusable workflow. Required reviews, status checks, resolved threads, and
stale-review dismissal still gate native auto-merge. A successful CodeRabbit
status on a skipped review does not itself supply an approving review.
Installing the Python package does not deploy or configure this workflow.

The payloads under `.github/settings/` are historical snapshots, not live state.
Before applying settings, reread the repository's workflow permissions, selected
actions, and active branch rulesets. Reconcile the intended change against those
live responses; preserve newer required checks, allowlist entries, and review
policy. Do not replay a snapshot against a hardcoded ruleset ID.

After a workflow change reaches `main`, exercise a new eligible Dependabot
event. Rerunning an old job retains its original workflow definition. Verify an
actual `github-actions[bot]` approval at the current PR head and a merge only
after required checks and threads pass. Ineligible updates retain manual review.
Record the PR, head, approval, workflow-run, and merge links.

Merges using `GITHUB_TOKEN` may not trigger push workflows. A manual `ci.yml`
dispatch against `main` can diagnose the merged commit, but cannot produce
release-eligible signed evidence. Before selecting a release commit, require a
successful main-push run as described in [Releasing](../RELEASING.md#reusable-ci-evidence).
The CI signer and crates.io publishing workflows require manual dependency review.
Local opt-in CodeRabbit CLI review remains separate.
