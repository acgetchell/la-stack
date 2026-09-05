# Rational row-clearing allocation study (#233)

## Decision and scope

The component borrowing proposed in #233 was already implemented in
`a62e72a` (#219). This study isolates that change, supplies missing allocation
evidence, and adds reproducible wide-component fixtures. It makes no further
production-code or public-API change.

The baseline is current source with only `integer_at_scale` replaced by the
former positive-canonical cloning expression:

```rust
let denominator = value.denom().clone();
let multiplier = scale / denominator;
let numerator = value.numer().clone();
numerator * multiplier
```

The comparison uses the existing production expression:

```rust
value.numer() * (scale / value.denom())
```

Both use canonical positive denominators and the same exactly divisible row
scales. The GCD/LCM calculation, elimination, pivoting, reconstruction, and error
contracts are identical. This is an isolated component-cloning experiment, not
a comparison against the entire revision preceding #219 (which also changed
denominator handling in the LCM calculation).

## Fixtures and independent validation

Every family covers D=2 through D=8:

- `small`: the unchanged `rational_input_d{2..8}` diagonally-dominant fixture.
- `wide256`: multiply row `r` by `(2^k + 1) / (2^(k+1) - 1)`, with
  `k = 256 + 16r`.
- `wide1024`: the same construction with `k = 1024 + 16r`.

The positive row factors preserve strict diagonal dominance and existing zeros.
Non-zero entries stress large canonical numerators and denominators; reduction
can remove small factors. Right-hand sides are rebuilt exactly from the known
solution `x[i] = (i+1)/(i+2)`. The shared private-field
`ValidatedRationalInput` is constructed only after independent rational Gaussian
determinant/sign and solve checks agree with the production results and known
solution. Construction and validation run outside measured closures.

These deterministic families isolate component size. They are not a random
rational corpus or a survey of unrelated large coprime denominators, and they
do not establish performance for all exact-input workloads.

## Allocation evidence

The [complete allocation CSV](../assets/rational-row-clearing-allocations.csv)
retains total allocation counts and requested bytes for both implementations.
`det`, `det_sign`, and `solve` rows measure complete public calls, including
dropping their results. The remaining rows isolate component cloning and
integer scaling, with matrix-only or augmented matrix/RHS row scales computed
before measurement. Those scaling rows model the private expressions above;
they are not public row-clearing APIs and exclude LCM construction.

The small family allocates nothing for either component cloning or isolated
scaling. The pinned num-bigint 0.4.8 supports inline single-digit storage, which
these fixtures use, so counting `.clone()` calls would overstate their heap
cost. Larger small-family public calls still allocate during arithmetic.

For both wide families, borrowing saves the following allocations per public
call. The public-call differences exactly match the isolated scaling
differences, separating component copies from the remaining arithmetic cost.

| D | Saved per determinant/sign | Saved per solve |
|---|----------------------------|-----------------|
| 2 | 4 | 6 |
| 3 | 9 | 12 |
| 4 | 16 | 20 |
| 5 | 21 | 26 |
| 6 | 31 | 37 |
| 7 | 43 | 50 |
| 8 | 57 | 65 |

The net saving is smaller than the two explicit component clones per non-zero
wide entry: changing multiplication operand ownership also changes its storage
reuse. The measured difference, rather than a source-level clone count, is the
relevant result. At D=8, for example, wide1024 determinant allocations decrease
from 2,959 to 2,902, while solve allocations decrease from 4,203 to 4,138. Most
allocation traffic remains in exact arithmetic.

The allocation-counter 0.8.1 dependency is used only by the separate
`rational_allocations` test executable. It does not replace the allocator in
Criterion timing executables or the library. Each measured closure is checked
three times for identical counts and zero retained allocations/bytes. Counts
are observations for the pinned dependencies and platform, not public API
guarantees or fixed cross-version test expectations.

## Timing evidence

The [complete timing CSV](../assets/rational-row-clearing-timings.csv) retains
all 63 initial comparisons and two repeats, including Criterion mean estimates
and their marginal 95% confidence intervals in nanoseconds. Percentage change
is `100 × (borrowed mean / cloned mean - 1)`; negative values mean less time.

| Family | Cases | Median change | Minimum | Maximum |
|--------|-------|---------------|---------|---------|
| small | 21 | -0.89% | -6.00% | +3.60% |
| wide256 | 21 | -1.05% | -7.00% | +3.35% |
| wide1024 | 21 | -0.36% | -2.97% | +1.01% |

These are descriptive summaries of per-case mean ratios, not a pooled workload
speedup or a paired confidence interval for the implementation change. Small
movements on a working desktop can reflect run conditions. The two cases whose
initial means increased by more than 2% were repeated in reverse implementation
order (borrowed, then cloned), with identical Criterion settings:

| Case | Initial change | Repeat change |
|------|----------------|---------------|
| small D=8 solve | +3.60% | -0.51% |
| wide256 D=2 determinant | +3.35% | +1.33% |

The initial slowdown magnitudes did not persist. This does not prove identical
timing, but supplies no stable material regression supporting a production
change. The allocation reduction is repeatable, while the overall timing
benefit is modest and workload-dependent. Retain the existing borrowed
expression; these results do not justify another row-clearing ownership design
or a general speedup claim.

## Reproduction

Record allocation evidence with:

```bash
cargo test --locked --release --features bench,exact \
  --test rational_allocations -- --nocapture --test-threads=1
```

Each `allocation,...` record contains family, dimension, operation, allocation
count, and requested bytes. The isolated cloned and borrowed rows are available
in the same run. For complete public-call comparisons, run once with the
baseline expression above and once with the existing borrowed expression.

Record focused timings with the counting executable stopped:

```bash
CRITERION_HOME=target/criterion-issue-233 \
  cargo bench --locked --features bench,exact --bench exact -- \
  'rational_input.*/(det_sign|det|solve)_row_cleared_bareiss' \
  --sample-size 30 --warm-up-time 1 --measurement-time 2 --noplot \
  --save-baseline issue233-cloned
```

Repeat the identical command with the borrowed expression and baseline name
`issue233-borrowed`. Only the expression and output baseline name change.
`bencher.iter` measures each complete public operation; fixture creation and
oracle validation remain outside timing. The Gaussian reference timings are
excluded from this focused comparison. The separate Criterion directory avoids
overwriting release measurements.

For the two focused repeats, use the filter below with the same settings and
baseline names `issue233-borrowed-repeat` and `issue233-cloned-repeat`, in that
order:

```text
^(rational_input_d8/solve_row_cleared_bareiss|rational_input_wide256_d2/det_row_cleared_bareiss)$
```

## Environment

The study ran on 2026-09-05 UTC (2026-09-04 PDT), starting from
`25aa4c02beda16c5efa66f23d1ebfb8b54532ab1`, with the shared harness and allocation
probe in this patch. Hardware was an Apple M4 Max MacBook Pro, 16 CPU cores
(12 performance, 4 efficiency), 64 GB memory, running macOS 26.6.2 (25G83).
Rust was 1.98.0 (`88d9e12ae`, LLVM 22.1.8), targeting aarch64-apple-darwin.
Cargo's release profile used fat LTO and one codegen unit; no custom Rust flags
or target override were set. Criterion was 0.8.2, num-bigint 0.4.8, and
num-rational 0.4.2. Timing runs were sequential, without concurrent builds or
tests; the desktop environment was not otherwise isolated.

SHA-256 identities of the measured source and executables:

```text
Cargo.toml
353e511917fbde3e73f6a10bdd05ca4a6a44cf94a91995689868a1a9714913db
Cargo.lock
f07e3a6a4366712f5a96d1b806e5c1562542e9af2087b58b00b800d8f02080ee
benches/exact.rs
8302e233df39f4ed5931afd94f5e774f30e7224912619bc7e426b4624f9c3c70
benches/common/rational.rs
fb835a8c859f48d16f612165ff5ab18ca5313b968c36dcc0f8399109ca6c548c
tests/rational_allocations.rs
6bb87fd3b7e90ed3c42a7ebf0354d88f2f76d2b232314a85cc53aafd1d867239
src/rational.rs (cloned experiment)
af0aeb30059a61922b1ef78662a9ad45e12447f6eedf40fa8f84ea25dedff93c
src/rational.rs (existing borrowed implementation)
c707a22ab9abae32d92a9a781563ddfd22d9ec6b2cba77666c0bc7c85160df2d
exact benchmark executable (cloned experiment)
59b6f1bd13398276c57a744450e42b22be95163958c4fe89e6860ada45b0bb51
exact benchmark executable (borrowed implementation)
41f7337ba4ee2245d226901b690a13f7401fd2b40159fb2b40ff8527fcebccfd
```

## Validation

- `just test-exact`: 753 tests passed, both before the experiment and on the
  restored production implementation.
- The standalone allocation command: 10 tests passed for each implementation;
  the final borrowed run reproduced all 189 recorded allocation rows exactly.
- Focused Criterion runs: 63 cases per implementation, plus the two reverse-order
  repeats. All measured fixtures passed their independent checks.
- `just ci`: passed, including 824 Rust tests, 612 Python tests, default and
  exact-feature doctests, static checks, example execution, and benchmark
  compilation.
- Both CSVs were checked for complete row counts, valid timing intervals and
  ratios, and equality between public-call allocation reductions and the
  isolated scaling reductions. The production `src/` tree has no final diff.
