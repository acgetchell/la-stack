# Benchmark candidate and resume review

Historical architecture and validation evidence. Adoption status and remaining
work below describe the reviewed revisions. The current operational procedure
belongs to [Releasing](../../RELEASING.md).

Independent Astra architecture checkpoint on 2026-11-07 for
[#325](https://github.com/acgetchell/la-stack/issues/325) and
[#323](https://github.com/acgetchell/la-stack/issues/323), requested by the
maintainer before workflow changes. Reviewed the working tree at HEAD
`f520e2f59574774454a739c7dfa810f3d1845420` and installed research-repo-tools
`13af6b60c9bc1712b329bfe9504c7a54840e043c`. This records a conditional design
decision and bounded experiments, not completed adoption or release readiness.

The follow-up below records adoption at the new pin. Earlier findings remain
as dated evidence and are not current blockers.

## Contents

- [Current adoption](#current-adoption)
- [Decision](#decision)
- [Resume contract](#resume-contract)
- [README input ownership](#readme-input-ownership)
- [Release ordering and trust](#release-ordering-and-trust)
- [Consumer acceptance](#consumer-acceptance)

## Current adoption

On 2026-11-07, upstream commit
`055353c630604b1ae3438c1eaf534ba85e067a7e`
([PR #205](https://github.com/acgetchell/research-repo-tools/pull/205)) supplies
the accepted #203 and #204 contracts. The preserved installed-package probes
pass unchanged. Schema 6 now retains completed suite pairs within an interrupted
trial; phase-specific sources truthfully declare baseline README input while
allowing current README output.

The consumer uses one sequential trial with all 1,403 library cases and 16
paired LU peers per phase. Current native gates inspect debug/release dependency
files and reject README compiler input before timing. Historical current tags
that still compile README are intentionally ineligible for this full plan;
non-exact diagnostics retain conservative input declarations.

The local adapter captures original raw bytes, preparation, candidate assertions
and publisher receipt before promotion in a separate owner-controlled directory.
Final acceptance uses an isolated installed supplier and a clean final commit.
Named immutable review selections permit explicit review of a newer native
attempt without overwriting the original acceptance. Capture recovery from
scratch is explicit; normal retries never adopt untrusted incoming assertions.

Astra accepted two maintainer dispatches as the authenticated handoff. The first
uses independent original asset IDs/hashes to retrieve, verify and attach inert
evidence; the second additionally pins the canonical attached archive and full
acceptance across deployment approval. Write-scope jobs execute no checkout code.
The native signer remains main-push-only. No new upstream capability was needed
for this consumer orchestration.

Follow-up review found and corrected local immutable-record bootstrap, runtime
isolation, capture recovery, final review selection, draft retry identity,
canonical asset verification, title/notes retention and Just variadic forwarding.
Review also identified a recovery gap after immutable GitHub publication:
native-evidence renewal could not replace the original acceptance asset.
The Cargo dispatch now takes an independently reviewed full acceptance digest,
freezes that exact receipt before approval and revalidates its bytes afterward.
The original archive, raw observations and all release assets stay unchanged.
Astra's final scoped review accepted this renewal mechanism without another
supplier change. Actual shell regression tests also require explicit failure
exits for the artifact guards, including on macOS's older Bash.
The full synthetic candidate path passes through capture, publication, final
binding, retry, archive derivation, hosted receipt verification and explicit
native-evidence renewal before and after approval. It uses the installed shared
candidate/archive verifiers and checks unchanged original bytes throughout;
provider responses and timings are synthetic. Both scratch and retained README
publication regressions pass. Final macOS `just ci` passes with 367 Python tests,
1,143 Rust tests, default/exact doctests, benchmark compilation and correctness
checks, examples and the required static validators. Windows-targeted Python
typing and explicit rule scans of the two new Python files also pass.
No additional supplier defect was found. Native consumer Linux/Windows CI,
a published supplier pin, real timing evidence and actual hosted release
handoffs remain separate release requirements. No dispatch or upload has run.
The release guide owns the operational commands and remaining publication gates.

For the v0.5.2-to-v0.5.3 availability case, acquisition falls from 2,838 local
observations plus 1,403 hosted repeats to the same 2,838 local observations.
At eight nominal warmup/measurement seconds per case, this removes 187.1
aggregate minutes (565.5 to 378.4). The old hosted feature suites ran in parallel,
so their nominal critical path was the slower 98.1-minute default suite plus
setup and assembly, not their summed sampling time. New verification/transfer
work remains. These are planning calculations; no comparable real campaign
establishes elapsed-time savings.

## Decision

Retain local maintainer-owned Step 1 acquisition. Preserve all 1,403 library
cases, availability, paired LU peers, sampling settings and original raw bytes.
Use one explicit sequential trial for an initial schema-6 migration; this
preserves the existing measurement budget and makes no counterbalancing claim.
README views must explicitly select trial 1. More trials require a scientific
reason and resource decision; they do not solve interruptions within a trial.

The shared candidate API supports the intended local acquisition, publication
receipt and final-commit binding. Accept that architecture subject to the input
and trusted-handoff conditions below. Joint resume/reuse completion requires
the supplier decisions in
[research-repo-tools#203](https://github.com/acgetchell/research-repo-tools/issues/203)
and [#204](https://github.com/acgetchell/research-repo-tools/issues/204), both
assigned to v0.2.7. Keep the existing release procedure until its replacement
passes acceptance. A dependency update or safe fallback does not close #323.

## Resume contract

The upstream [paired-suite decision](https://github.com/acgetchell/research-repo-tools/pull/175)
accepts a whole suite with both available phases as a durable schema-5 unit.
An incomplete pair is recollected in full after retained-unit validation;
completed raw bytes and original host/session attribution remain unchanged.

Candidate reuse requires schema 6 and complete-run v5. At the reviewed pin,
its journal retains only complete trials. A started trial without its complete
receipt rejects resume because it may contain unretained timings. The scheduled
engine does not pass a nested checkpoint to suite acquisition. With one trial,
no incomplete campaign can recover its already completed suites.

A bounded installed-package pytest probe uses two synthetic paired suites,
`shared` and `later`. It interrupts before `later/current`, after verifying
the command log contains `shared/baseline`, `shared/current` and
`later/baseline` measurements. The same success contract runs against both
schemas: schema 5 passes without repeating `shared`; schema 6 fails with
`interrupted trial may contain unretained timings`. The failing success
assertion is preserved, not changed to accept a fresh acquisition.

The probe and unchanged supplier fixture are retained under
`target/research-tools-rehearsal/13af6b60c9bc1712b329bfe9504c7a54840e043c/consumer/`.
They are ignored diagnostic artifacts, not ordinary test-suite inputs or
performance evidence. #203 contains the reproduction and required contract.

The supplier must independently review recovery granularity, original session
and timestamp retention, interruption history, phase order and incomplete-pair
recollection. Retain rejection of ambiguous or corrupt work; deleting the
existing guard would not establish safe recovery. Final resume wiring waits
for this accepted public contract.

## README input ownership

Both current source and v0.5.2 use an unconditional
`#![doc = include_str!("../README.md")]`. README is therefore a compiler input.
The candidate contract requires truthful input declarations and prohibits
input/output overlap, including documents read by builds. Its final comparison
checks inventories from every observation, including baseline phases.

Existing Rust dependency files and a fresh Rust 1.99.1 compiler probe confirm
the unconditional include reads README during ordinary compilation. A tiny
`cfg_attr(doc, doc = include_str!(...))` control excludes it from ordinary
compiler dependencies. This control does not establish native benchmark,
rustdoc or package equivalence for la-stack.

A separate installed-API probe loads the full la-stack plan, adds README to
the source inputs, and selects a sequential schema-6 schedule. Declaring README
as an output then correctly rejects with `declared output overlaps measured
inputs`. Omitting the input to bypass this check would violate the contract.

A current-source `cfg(doc)` change is an option to investigate, not an accepted
complete fix: the immutable baseline still includes README, and source patterns
are shared across phases. #204 requests a supported phase/source-aware contract
or demonstrated consumer architecture. Do not overlay changed library source
as harness, weaken input ownership, or remove crate documentation. Any accepted
consumer change needs both-phase compiler-input proof, default/exact rustdoc
and example validation, native gates and final package evidence.

## Release ordering and trust

The following is the reviewed target design; it is not yet executable guidance.
The [release guide](../../RELEASING.md) remains the operational owner.

1. Finish versions, dependency pins, lockfiles, source, harness, commands and
   output declarations before Step 1. Any subsequent non-output tree change
   invalidates reuse, including a new supplier pin.
2. Acquire once through the public release API. Its `before_promotion` hook
   freezes `CandidatePolicy(commit=None, report=...)`, the independently retained
   acquisition digest, publication receipt and separate receipt digest. The
   hook does not publish. A failed hook leaves raw evidence for a direct retry.
3. Keep the owner record outside candidate-controlled files with an immutable
   backup or authenticated transport. An ignored `target/` directory alone is
   not an authentication boundary. A checksum stored beside the candidate does
   not authenticate that candidate. Preserve local ownership throughout.
4. Render report/history through the exact shared-publisher receipt and README
   from the selected retained trial. Retry rendering or promotion directly;
   do not call acquisition again. Do not exempt arbitrary JSON or directories.
5. After merge and final native CI, run candidate verification before tagging.
   Independently select the final clean commit and supply the original trusted
   digests. Retain acceptance outside committed source. Step 2 must not measure,
   discover, preflight, rebuild, tag or publish; it may write its requested
   verification receipt.
6. Derive and verify the archive through public APIs. Replace the hosted second
   acquisition only after exact artifact selection and attachment pass tests.
   Keep original measured revisions, host, sessions and timestamps distinct
   from final-release acceptance.
7. Freeze the full acceptance and exact native CI/artifact selection before
   approval; revalidate that same selection after approval and before release
   credentials. Preserve protected environments, draft/tag checks, package
   guarantees, no-overwrite checks and credential separation.

The existing release publication event transfers native CI identity, not
independently authenticated candidate assertions. The recommended transport is
a maintainer-authorized workflow dispatch supplying acquisition, publication,
acceptance and archive digests from the separate owner record. Authenticate the
allowed workflow, dispatcher, event, ref, commit and exact run/attempt/artifact;
compare separately transported files against those assertions. Later jobs must
retrieve that frozen selection, never a latest-by-name artifact or digests
computed from the candidate being checked.

Actual handoff ownership and authentication require consumer implementation and
demonstration. No new service, signer or protected-environment policy is assumed
to exist. Keep generic verification upstream and use thin consumer orchestration.
Local capture/receipt prototypes can proceed independently, but they cannot
substitute for authenticated hosted publication.

## Consumer acceptance

- Preserve independent inventories, availability and same-suite peer attribution
  in bounded fixtures. Synthetic samples establish workflow behavior only.
- Exercise one Step 1 acquisition through publication, final clean merge,
  signed native evidence, candidate acceptance, archive and attachment retries.
  Compare retained raw bytes and original provenance throughout.
- Put failure sentinels on acquisition, discovery and preflight in Step 2,
  report/README retries, archive derivation and approval-time revalidation.
- Reject input, declaration, dependency, source, harness, compiler-relevant
  document and frozen generated-output drift before publication or new timing.
- Reject replaced owner assertions, receipts, acceptance files, run attempts
  and immutable artifacts, including replacement after initial approval checks.
- Prove eligible missing-suite recovery under the accepted supplier contract;
  retain partial-output, corruption, ordering and writer-ownership rejection.
- Complete focused consumer checks and native platform evidence before replacing
  the release workflow. Development pins remain forbidden at actual release
  boundaries until a published supplier is adopted and revalidated.

The original architecture checkpoint performed no real timing campaign, release
upload, publication, repository Git mutation or production resume/reuse wiring.
The current adoption section records the later consumer implementation. Group 2's
separate [validation record](benchmark-workload-review.md#validation-evidence)
retains its own source and supplier identities.
