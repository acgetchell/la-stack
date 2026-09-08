# Local Benchmark Summaries

This directory preserves complete local measurement summaries across releases
and `just clean`. The generated [performance report](../performance.md) is a
selected view of those measurements; the README chart uses its LU-solve subset.

## Contents

- [Saved runs](#saved-runs)
- [Regenerating reports](#regenerating-reports)
- [Comparing measurements](#comparing-measurements)
- [Historical availability](#historical-availability)

## Saved runs

`just performance-release` saves each successful local comparison under
`<current>-vs-<baseline>/<run-digest>/`. The digest identifies the complete
contents, so another run or machine creates a separate snapshot. Identical
promotion is idempotent. Earlier snapshots remain available.

Each snapshot contains:

| File | Contents |
| --- | --- |
| `performance.full.csv` | Every recorded baseline and current benchmark, including diagnostic and peer cases outside the report selection |
| `performance.full.provenance.json` | Full-summary schema, CSV digest, phase counts, and the source, harness, dependency, CPU, OS, toolchain, command, and validation provenance |
| `performance.csv` | Selected report inputs, including peer context |
| `performance.provenance.json` | Validated metadata and digest for the selected report inputs |

The full CSV stores the Criterion benchmark ID, measurement phase, mean and
median in nanoseconds, their 95% confidence intervals, and sample count. Numbers
retain floating-point precision rather than the rounded values displayed in
Markdown. Each recorded case must have valid metadata, estimates, and 100 raw
samples before publication. Baseline-only peers remain explicitly baseline
measurements; missing current measurements are never filled from an older run.

`latest.json` points to the most recently promoted complete snapshot. It is
updated with the report in the same rollback-protected publication operation.
Review and commit the snapshot and pointer along with the release report.

## Regenerating reports

After a successful release promotion, these commands can use the committed
snapshot when the scratch report inputs are absent:

```bash
just performance-doc
just performance-readme
```

Neither command reruns benchmarks. Existing scratch inputs take precedence;
partial or corrupt inputs fail validation rather than silently selecting an
older snapshot. README publication still checks the measured source and
benchmark contract against the current checkout, along with the lockfile and
release version. Committing only the saved artifacts does not invalidate a
snapshot: its original measurement commit remains in provenance. Legacy inputs
without a benchmark-contract digest also require the original commit.

`performance-local` initially keeps complete summaries beside its selected
inputs under `target/bench-reports/`. Promote a distinct-release comparison with
`performance-doc` before cleaning if it should become committed release history.
Same-version local experiments remain scratch output.

## Comparing measurements

Compare rows by benchmark identity and measurement phase, with the recorded
environment and harness held comparable. Baseline and current rows come from
separate executions. New cases can have only a current measurement, and peer
measurements reused by the report can have only a baseline measurement.

These CSVs support comparisons of recorded point estimates and intervals.
Recomputing Criterion statistics requires the raw samples, which are not stored
here. The GitHub release workflow separately publishes full raw Criterion
archives from its hosted runner; those are a different measurement environment
and do not substitute for local results.

## Historical availability

This retention workflow cannot recover measurements previously discarded with
temporary worktrees. If `latest.json` is absent, no complete local snapshot has
been promoted yet. Existing [archived Markdown reports](../archive/performance/README.md)
retain their selected summaries, but do not establish complete historical
coverage. A new authorized measurement run is required to fill that gap.

Completed optimization investigations and decisions are in
[archived studies](../archive/performance/studies/README.md).
The [Benchmarking guide](../BENCHMARKING.md) owns execution procedures.
