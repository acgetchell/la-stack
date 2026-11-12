# Release benchmark workload review

Historical review and development-pin validation evidence. Pending work and
supplier limitations below describe the reviewed revisions. For current commands
and policy, use [Benchmarking](../../BENCHMARKING.md) and
[Releasing](../../RELEASING.md).

Consumer review for [#322](https://github.com/acgetchell/la-stack/issues/322),
with availability and progress integration findings for
[#320](https://github.com/acgetchell/la-stack/issues/320) and
[#324](https://github.com/acgetchell/la-stack/issues/324).
Initial review on 2026-11-05 at source
`f520e2f59574774454a739c7dfa810f3d1845420` used the
`public-release-and-lu-peers-v3` policy. The 2026-11-06 integration uses
`public-release-and-paired-lu-peers-v4` without removing library workloads.

## Contents

- [Selection decision](#selection-decision)
- [Required coverage](#required-coverage)
- [High-cost combinations](#high-cost-combinations)
- [Optional and excluded measurements](#optional-and-excluded-measurements)
- [Availability adoption boundary](#availability-adoption-boundary)
- [Budget and progress](#budget-and-progress)
- [Validation evidence](#validation-evidence)
- [Remaining acceptance](#remaining-acceptance)

## Selection decision

Retain all 1,403 current library cases: 736 default and 667 exact. No workload
removal or sampling reduction is proposed. The exact required IDs remain in
[`performance-policy.toml`](../../../tooling/performance-policy.toml); this review
explains their coverage rather than introducing another selection mechanism.
The full native catalog remains available for explicit diagnostics.

The two largest repeated families measure different caller costs and numerical
regimes. Their counts alone do not establish equivalence. No comparable
end-to-end timing campaign was performed for this review, and no speedup is
claimed. Independent Astra review on 2026-11-05 accepted this no-removal
conclusion, subject to the required native discovery and correctness checks.
That review accepts scientific coverage, not minimality or measured runtime.

## Required coverage

Each row below is required in the current release selection. Construction and
conversion stay timed when they are part of the named public operation. Fixture
generation and independent oracle evaluation stay outside timing. The native
targets and filters are owned by [`performance.just`](../../../tooling/performance.just).

| Native family / feature | Cases | Calling forms, algorithm regimes and dimensions | Independent correctness evidence |
| --- | ---: | --- | --- |
| `angle` / default | 168 | Prepared vectors, checked borrowed slices, construction plus angle; 14 dense/sparse, endpoint, mixed-scale and subnormal families at D=3–6 | Analytical references and integer Gram/exterior-product oracles in `benches/common/angle.rs`; `tests/angle_bench_inputs.rs` and `tests/vector_angles.rs` |
| `column_basis` / default | 55 | Square/tall/wide, deficient, near-dependent, rescaled and affine constraints; shapes through 8 × 14 and 6 × 7 | Expected rank and pivot sets in the native harness; independent exact minor checks in `tests/column_basis.rs` |
| `exact` binary64 inputs / exact | 81 | Determinant, sign, solve, strict/rounded results; closed forms through D=4 and elimination at D=5; fixed, deterministic corpus, near-singular, large-entry and Hilbert inputs | Validated exact fixtures and independent rational Gaussian elimination; `tests/exact_bench_config.rs` |
| `exact` canonical conversion / exact | 128 | Four scalar/vector result forms, eight representability/size families at D=2–5; dyadic, non-dyadic, subnormal, underflow, overflow boundary and wide integers | Known exact conversion values and typed error expectations in the native fixture gates and exact tests |
| `exact` rational inputs / exact | 63 | Determinant, sign and solve; small, 256-bit and 1,024-bit inputs at D=2–8, preserving closed-form/elimination and integer-size regimes | Independent rational Gaussian determinant/solve and manufactured solutions |
| `gram` / default | 56 | Public Gram construction for square and embedded vectors through coordinate dimension 8; orthogonal, dependent, near-dependent and mixed-scale | Independent integer matrix-product oracle in `benches/gram.rs` |
| `interval` / default | 43 | Scalar arithmetic, lifted assembly and prepared determinant/sign; exact, cancellation, zero, underflow, sparse lifted and dense conclusive/inconclusive regimes at D=3–7 | Fixture-specific expected enclosures/signs in the native harness; interval unit/property tests |
| `interval_solve` / default | 35 | Five division/reciprocal cases and five 32-system solve corpora at D=1–6; dense, uncertain, near-singular, singular and range-exhaustion outcomes | Manufactured solutions plus independent rational Gauss-Jordan checks in `tests/interval_solve.rs` |
| `linear_form` / default | 140 | Dot/affine-difference certificates and endpoints at D=2–6; dense, sparse, cancellation, underflow/overflow; construction, prepared and construction-inclusive projection batches | Independent fixture expectations and certified reduction unit/property tests |
| `polynomial` / exact | 5 | Full sparse determinant polynomial plus borrowed leading term at orders 4–8; allocation/destruction included | Sparse permutation enumeration of every coefficient in `tests/polynomial_bench_inputs.rs` |
| `public_api` / default | 90 | Construction, access, reset-copy mutation, symmetry, tolerances, direct determinant/error-bound and interval queries at D=2–5; dimensionless scalar operations once | Known values and typed failure contracts in library/integration tests, with harness assertions |
| `public_api` / exact | 42 | Rational storage/construction/mutation and scalar strict/rounded conversion; D=2–5 where dimensional | Exact-value and conversion/error tests; native fixture assertions |
| `rational_dot` / exact | 54 | Prepared vector and checked slice; non-dyadic, cancellation and 512-bit integer families at lengths 2–8, 14 and 22 | Analytical reciprocal-product sums; `tests/rational_dot.rs` |
| `rational_solve` / exact | 294 | Construction, prepared solve and runtime-dispatch adapter; seven matrix structures, two coefficient families, D=2–8 | Manufactured exact solutions, independent Gaussian control and residuals; `tests/exact_bench_config.rs` and `tests/rational_structured_solve.rs` |
| `vs_linalg` library rows / default | 149 | LU/LDLT; factorization plus one solve/determinant and prepared operations; dot/norm regimes at D=2–64; five additional D=8 pivoting, conditioning and range cases | `benches/comparison/tests/vs_linalg_inputs.rs`, factorization unit/property tests |

The original README policy additionally required 16 baseline-origin peer cases:
nalgebra and faer factorization plus one solve at D=2, 3, 4, 5, 8, 16, 32 and 64. The matching
eight la-stack cases are already included above. Same-environment attribution
and the timed boundary must survive an orchestration change: LU factorization
plus one solve from prepared inputs, including required workspace and result
destruction. Matrix/vector input representations are constructed outside timing.
Hosted library-only evidence cannot replace the local denominator. The
[availability decision below](#availability-adoption-boundary) replaces the
baseline-only peer selection in new acquisitions with a paired suite.

Default constructors share their explicit construction workload. Formatting,
error factories, enum queries, indexing/dispatch macros and re-exported traits
are tested directly or exercised by their containing public workflow; they do
not require separate timing rows. This is a public-workflow coverage argument,
not a claim that every possible parameter combination has been measured.

## High-cost combinations

`rational_solve` has 7 structures × 2 coefficient families × 7 dimensions × 3
calling forms = 294 cases per revision. Diagonal, upper/lower triangular and
row-permuted triangular inputs exercise structure recognition and substitution;
sparse and dense inputs retain elimination behavior. Dyadic versus non-dyadic
coefficients change denominator work. D=2–8 retain separate const-generic
instantiations and the library's dimension/dispatch boundaries. The adapter
includes construction, dimension dispatch and output collection; the prepared
solve does not. Construction alone isolates a caller-visible cost. No equivalent
representative has been established that justifies removing a combination.

`angle` has 14 fixtures × 4 dimensions × 3 calling forms = 168 cases per revision.
Identical/antipodal endpoints, nearly parallel/antipodal inputs, sparse support
patterns, independent scaling and subnormal behavior are different numerical
questions. Prepared finite vectors avoid checked-slice validation; construction
plus angle includes two constructors. A single generic kernel does not prove
that distinct const-generic compiled dimensions or caller forms have equal cost.

Together these families contribute 924 paired cases, or 123.2 nominal minutes
at 8 seconds per case. That identifies review cost; it is not evidence of
redundancy or observed runtime. Correctness/property tests remain comprehensive
regardless of any future timing-selection decision.

## Optional and excluded measurements

| Classification | Cases / selectors | Reason |
| --- | --- | --- |
| Optional diagnostics, excluded from mandatory release timing | `angle_d*/stable_control`; hand-written Gram controls | Reference kernels omit or rearrange public validation/work and answer implementation questions |
| Optional diagnostics, excluded from mandatory release timing | `rational_solve/*/legacy`; rational Gaussian controls and exponent-spread studies in `exact` | Independent reference/implementation investigation; correctness oracles remain required outside timing |
| Optional diagnostics, excluded from mandatory release timing | Broad nalgebra/faer operations outside the 16 LU-solve peers | Existing README question is factorization plus one solve through D=64; full comparison remains explicitly runnable |
| Excluded only on historically unavailable revisions | 55 column-basis, 35 interval division/solve, 54 rational-dot and 5 polynomial cases before v0.5.2 | Missing public APIs; all 149 stay required on supported revisions, never optional because compilation failed |
| New exclusions proposed by this review | None | No scientifically justified reduction established |

## Availability adoption boundary

The initial review used **research-repo-tools 0.2.6**, resolved from PyPI.
The development branch now tests exact upstream commit
`055353c630604b1ae3438c1eaf534ba85e067a7e`; its unchanged 0.2.6 version label is
not a published-artifact identity. Follow the
[development pin policy](../../../tooling/README.md#shared-capability-adoption)
before release. On 2026-11-05 the
[exact-version registry metadata](https://pypi.org/pypi/research-repo-tools/0.2.6/json)
reported non-yanked wheel and source artifacts whose SHA-256 hashes match
the original `uv.lock`. The package-identity integration test now checks the
exact manifest, lock and installed commit, including isolated release runtimes.

Its [published schema-5 contract](https://github.com/acgetchell/research-repo-tools/blob/v0.2.6/docs/DECLARING_PERFORMANCE.md)
provides inclusive `since` bounds resolved from captured release/source receipts
before overlay or native commands. This is the intended replacement for the
permanent `[baseline].current-only` exclusions. The required resolved inventory
is 1,254 library cases on v0.5.0/v0.5.1 and 1,403 on v0.5.2 and later, plus the
16 peers on each phase. Reading an overlaid `Cargo.toml` cannot establish the baseline
release. Supported compilation/discovery failures must still fail.

**The original phase selection was incompatible with schema 5.** A bounded synthetic
probe of the installed public `parse_suite_plan`, `measure_suites`,
`validate_run_policy` and `comparison_incompatibilities` interfaces used one
accepted paired suite, identical commands, full raw samples and known local
host observations. Library-against-library cross-release comparison and
same-phase peer comparison were eligible. Comparing baseline
`d2/nalgebra_lu_solve` with current `d2/la_stack_lu_solve` returned:

```text
workload IDs differ within the semantic contract
```

This fails even before accounting for the different commands needed to acquire
peers only once. The interface also requires matching cross-phase commands and
rejects cross-unit ratios. Splitting peers into a separate unit, relabeling their
phase, reusing one ID for different implementations, or suppressing ratio checks
does not preserve the consumer contract. Measuring peers on both revisions
would change acquisition ownership and budget, and still would not make the
baseline-peer/current-library comparison eligible.

Independent Astra review recommended preserving this deliberate supplier
restriction. The new declaration puts the eight la-stack LU-solve cases and
16 peer cases in one paired schema-5 suite with identical native commands on
both revisions. Those eight library cases are excluded from the other suites to
preserve single ownership. The README selection compares current-phase peers
with current-phase la-stack; the suite also retains baseline/current library ratios.
Keep historical baseline-origin evidence unchanged and attribute each new
observation to its actual phase.

This adds 16 peer measurements, or 128 nominal seconds, about 0.6% of the full
planned campaign. The repeated peers serve distinct paired evidence. Closer
collection reduces exposure to phase/order drift, but ratios remain descriptive.
No widening of supplier ratio eligibility is needed. A future insistence on
baseline-only peers with cross-phase ratios would instead require a separately
reviewed contract for operation/input equivalence, implementation identities and
allowed selection differences, with all provenance safeguards preserved.

The full declaration now uses seven schema-5 suites: shared default, shared
exact, paired LU solve, and four suites introduced at v0.5.2. Native discovery,
correctness gates, independent inventory tests, series selections and progress
budgets have been reconciled. At commit `a4640b0`, end-to-end acceptance was blocked by
[research-repo-tools#196](https://github.com/acgetchell/research-repo-tools/issues/196):
the publisher read dependency-version labels from aggregate phase provenance,
but schema 5 records them in the owning suite observation. Both scratch and
retained plot publication failed despite valid captured dependency versions.
The publication regressions and version labels were preserved.

Preparation had a separate early-validation gap tracked in
[research-repo-tools#197](https://github.com/acgetchell/research-repo-tools/issues/197).
The schema-5 selector skipped active suites' minimum-release checks and could save
a successful preparation for v0.4.6 despite the v0.5.0 floor. A correctly
versioned temporary Git fixture and an independent public-API reproducer both
confirmed the gap. Acquisition still rejected that pair before native commands;
the finding did not establish invalid measurement acceptance.

Commit `a8dc578a7583f1e4d811a5bad1ceb47ae28949e4`
([PR #199](https://github.com/acgetchell/research-repo-tools/pull/199)) supplies
both fixes. The unchanged public-API reproducers pass against its installed
package. The current development pin also supplies the root dependency-source
policy from PR #198; that declaration is retained with the harness inputs.

This report does not claim #320 complete or authorize a long measurement
campaign. Upstream contract history is
[research-repo-tools#175](https://github.com/acgetchell/research-repo-tools/pull/175).
PR [#190](https://github.com/acgetchell/research-repo-tools/pull/190), inspected at
`19476650960c4a441121a0a8b98ee2cc0e75ebe4`, addresses committed publisher evidence
and release summaries. Its `suite_runs.py` is byte-identical to installed 0.2.6,
so it does not change this intentional comparison policy.

## Budget and progress

| Inventory | Measurements | Nominal warmup plus sampling |
| --- | ---: | ---: |
| Existing v0.5.1 baseline policy: 1,254 baseline + 16 peers + 1,403 current | 2,673 | 356.4 min |
| Reviewed paired-peer plan with v0.5.0/v0.5.1 baseline: 1,254 baseline + 1,403 current + 32 peers | 2,689 | 358.5 min |
| Reviewed paired-peer plan with v0.5.2-or-later baseline: 1,403 baseline + 1,403 current + 32 peers | 2,838 | 378.4 min |
| Hosted current library only | 1,403 | 187.1 min |

These totals use the existing 3-second warmup and 5-second measurement target;
they exclude preparation, compilation, gates, analysis and publication and are
not deadlines. No new observed acquisition cost is available.

The canonical release recipe exposes the installed reporter's human/JSON/none
mode and heartbeat cadence. Native diagnostics remain on their original streams;
the reporter reads no native logs and counts only validated phase results.
The full schema-5 plan reports resolved per-suite inventory, nominal budget,
completed/remaining cases, stage, elapsed seconds and unknown ETA for unseen
work. The former hardcoded native budget line has been removed. Synthetic
consumer checks verify 2,689 validated observations for a v0.5.1 baseline and
the suite/phase counts at the v0.5.2 availability boundary.
The shell runtime record remains the observed recipe-body duration and exit code.

The initial consumer gate-failure regression exposed a supplier defect: schema-2
acquisition raises `ValueError("baseline preflight gate failed")` correctly,
but its reporter closes the gate stage with `finished` before checking the
nonzero exit status. It emits no `failed` event through `measure_prepared_pair`.
The regression requires the failure event and zero accepted cases.
The defect is tracked in
[research-repo-tools#191](https://github.com/acgetchell/research-repo-tools/issues/191).
The inspected #190 head has exactly the same `common_measurement.py` blob as
installed 0.2.6 (`27f1692981cb590ab63482ae538239b324f1a9d2`), so that PR does not
fix this failure. Commit `ef783cac608bc0d2d4657052c20e663f2d45b232` fixes it;
the preserved prepared-pair and baseline/current schema-5 regressions pass.
Keep those regressions and adopt the supplying published release before final
publication; do not add consumer event rewriting or relax gate acceptance.

See [the benchmark guide](../../BENCHMARKING.md#release-acquisition-policy) for
commands and retry boundaries. Report/README retries consume retained evidence
without measurement. Heartbeats, nominal budgets and estimates never replace
raw-file, sampling, correctness or freshness validation.

## Validation evidence

The initial local macOS/Python 3.15 checks on this draft were:

- `just check`: passed, including documentation, spelling, configuration,
  Python types/lint, shell, Rust checks and repository static analysis.
- Both affected consumer modules: 45 passed, 1 failed in 279.73 seconds.
  The failure is the preserved gate-status regression described above;
  report/README retry and source/declaration freshness checks passed.
- Exact manifest/lock/installed-package identity test: passed.
- Windows-targeted strict type checking of both changed Python modules: passed.
  This is static evidence, not native Windows execution.
- `uv lock --check` and `git diff --check`: passed.

The 2026-11-06 rehearsal against installed commit `a4640b0` additionally verifies:

- 249 supplier public-consumer tests and 734 subtests passed outside the supplier
  checkout, covering acquisition, resume, rendering, release handoffs and the
  development-pin policy with synthetic samples.
- Eight focused consumer availability checks passed: v0.5.0, v0.5.1, exact
  v0.5.2 and later boundaries, captured source receipts, unsupported releases,
  and rejection of incompatible historical evidence without rewriting it.
- Native discovery and correctness gates passed on exported real v0.5.1 and
  current sources with separate compiler storage. Required inventories were
  exactly 1,270 baseline and 1,419 current cases, including 16 peers per phase;
  no benchmark timing samples were collected.
- Ordinary toolchain verification accepts the exact development commit.
  Release-only dependency validation rejects it before acquisition or fetching,
  as required by upstream #193. CI omits release signing for development pins.
- The public publication regression for upstream #196 fails in isolation with
  a captured dependency version and two synthetic samples. Downstream scratch
  and retained publication checks reproduce the same failure.
- The public preparation regression for upstream #197 fails without any native
  command or provider access. The full CLI separately saves an unsupported pair
  in a temporary Git fixture. Both upstream issues target v0.2.7.
- `just check` passed, including Markdown/spelling, guide documentation builds,
  Python lint/types, configuration/workflow checks and static analysis.
  Windows-targeted type checking passed; this is not native Windows execution.
- The full downstream Python run had 316 passes and seven failures. Focused
  reruns after local assertion/configuration fixes passed 69 preparation/release
  checks and both source/harness freshness checks. Four regression failures
  remained at that pin: three publication paths blocked by #196 and unsupported
  preparation accepted by #197. That downstream suite was not green.

The 2026-11-06 rehearsal against installed commit `a8dc578` verifies:

- Both unchanged public-API reproducers for #196 and #197 pass.
- 256 supplier public-consumer tests and 756 subtests pass outside the supplier
  checkout, covering synthetic acquisition, resume, rendering, release handoffs,
  dependency-source policy and the preserved regressions.
- All 125 focused downstream workflow, preparation, toolchain and publication
  checks pass. Scratch and retained publication preserve the captured dependency
  versions; unsupported preparation fails before native work.
- A separate consumer regression at the exact v0.5.2 availability boundary
  verifies that failed column-basis discovery remains a failure with zero
  accepted cases. It cannot silently exclude a supported workload or proceed
  to gates and acquisition.
- The root dependency-source declaration covers Python tooling and the Cargo
  workspace, including the private comparison member. Exact-pin checks accept
  the development identity; release-entry checks reject it without creating
  evidence or fetching source revisions.
- All 89 native Rust, manifest, lock, toolchain and driver inputs match the
  preserved native-discovery receipt above. The changed shared-tooling pin and
  dependency policy have new evidence identities; the old receipt is not
  relabeled as a new measurement.
- `just check`, Python lint/types, Windows-targeted strict type checking,
  `uv lock --check` and patch whitespace checks pass. Windows type checking is
  static evidence, not native Windows execution.
- Final `just ci` passes: 325 Python tests and 1,143 Rust tests, plus default/exact
  doctests, benchmark compilation, exact-benchmark test mode, examples and the
  repository's static checks. No timing campaign is part of this validation.

A subsequent maintainer-requested `just update` completed tool setup and Cargo
maintenance but exposed a separate shared-updater gap:
[research-repo-tools#200](https://github.com/acgetchell/research-repo-tools/issues/200).
Python pin discovery rejected the exact Git source in the preserved tooling
group before resolution. An installed public-API reproducer confirmed that
failure without modifying its manifest or creating a lock. Explicit uv lock
refresh and synchronization succeeded with unchanged pins and lock contents;
the normal update command at `a8dc578` required the upstream fix. This maintenance
path was not part of the preceding acquisition/publication acceptance suite.

On 2026-11-07, commit `13af6b60c9bc1712b329bfe9504c7a54840e043c` supplies
[the #200 fix](https://github.com/acgetchell/research-repo-tools/pull/201).
The unchanged reproducer passes, and the complete normal `just update` succeeds
while preserving the exact supplier commit and inherited tool pins. Its full
lock refresh advances only `sse-starlette` from 3.6.0 to 3.7.0. The manifest,
lock, isolated release runtimes and consumer identity assertions use the new
commit; prior rehearsal receipts retain their original identities.
The 72 focused downstream toolchain/publication tests and final `just check`
also pass on this pin. An initial online workflow audit hit a GitHub timeout;
the required online retry and subsequent complete check both passed. This
does not relabel the preceding full `just ci` receipt as a run on the new pin.

The 2026-11-07 continuation uses
`055353c630604b1ae3438c1eaf534ba85e067a7e` from upstream PR #205.
`just update`, locked synchronization and the three preserved #203/#204
regressions pass. The full declaration now uses schema 6 and one sequential
trial; an interrupted trial retains completed suite pairs. The consumer resume
fixture preserves completed LU observations, reacquires only the interrupted
column pair and rejects changed input before invoking commands.

Native discovery and correctness gates pass for all seven suites at the
captured v0.5.2 baseline and working-tree current source: 1,419 cases per phase,
including all 1,403 library cases and 16 LU peers. Dependency files confirm
README as a baseline compiler input and its absence from both current build
profiles. These receipts precede final hosted-renewal orchestration edits;
those edits do not alter the native workload, sampling or gate commands.
The independent current-input guard has six passing failure/success fixtures.

The focused hosted-transfer and publication suite has 90 passing tests, and
all nine command-discovery tests pass. Host and Windows-targeted Python type
checks, `just check`, Markdown validation and rustdoc checks pass. The final
comprehensive validation and synthetic renewal receipt are recorded in the
[candidate review](benchmark-candidate-review.md#current-adoption).

These results establish discovery and orchestration behavior, not performance.
No real timing campaign, native Windows/Linux execution, release publication or
repository Git mutation was performed. Temporary synthetic Git repositories
are test fixtures only.

## Remaining acceptance

Group 2's development integration now passes local bounded acceptance. The
remaining release requirements are:

- Before a real release, replace development pins with the normally published
  supplier version, resolve/sync the lock, and revalidate the changed inputs.
  Existing rehearsal evidence cannot be relabeled for the published pin.
- Complete native platform checks and any separately authorized real campaign.
  Preserve historical report/raw bytes under their original plan identities.

The bounded rehearsal does not establish release readiness. No issue is closed
by this report.

The subsequent [candidate/resume architecture checkpoint](benchmark-candidate-review.md)
records the separate prerequisites for #325 and #323. Those findings do not
change group 2's scientific selection or establish completion of later groups.
