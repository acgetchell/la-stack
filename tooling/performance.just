# Captured native commands run in the measured checkout with its phase environment.
# Pass "list" to discover Criterion IDs without collecting samples.

# Cargo runs benches from each package root, including benches/comparison.
# Resolve the shared output tree before those working directories change.
export CRITERION_HOME := absolute_path(env("CRITERION_HOME", "target/criterion"))

# Compare native targets on the measured host regardless of a caller's target override.
host_target := `rustc --print host-tuple`

# Each public workload has one owner per phase. Peers are acquired once locally.
baseline-all mode="measure": budget (baseline-public-default mode) (baseline-public-exact mode) (baseline-peers mode)

baseline-peers mode="measure": (peers mode)

baseline-public-default mode="measure": (public-default mode)

baseline-public-exact mode="measure": (public-exact mode)

baseline-vs_linalg mode="measure":
    cargo bench --locked --target {{ host_target }} -p la-stack-comparison --features bench --bench vs_linalg -- --noplot {{ if mode == "list" { "--list" } else { "" } }}

# The independent case catalog describes the complete local acquisition budget.
budget:
    uv run --locked --no-sync python -c 'import pathlib, tomllib; policy = tomllib.loads(pathlib.Path("tooling/performance-policy.toml").read_text()); default, exact = (len(policy["release"][name]) for name in ("default", "exact")); total = 2 * (default + exact) + 16; print(f"Local release: {default} default + {exact} exact cases per revision, 16 peer cases once; {total} measurements, {total * 8 / 60:.1f} nominal minutes at 3 s warmup + 5 s measurement (excludes compilation, gates and analysis).")'

current-all mode="measure": (current-public-default mode) (current-public-exact mode)

current-public-default mode="measure": (public-default mode)

current-public-exact mode="measure": (public-exact mode)

current-vs_linalg mode="measure":
    cargo bench --locked --target {{ host_target }} -p la-stack-comparison --features bench --bench vs_linalg -- la_stack --noplot {{ if mode == "list" { "--list" } else { "" } }}

gate:
    cargo test --locked --target {{ host_target }} --workspace --features bench,exact --test vs_linalg_inputs --test exact_bench_config

la-stack mode="measure" sample="": (public-default mode sample) (public-exact mode sample)

# Shared native argument construction keeps discovery and measurement identical.
native suite filter features mode="measure" sample="":
    #!/usr/bin/env bash
    set -euo pipefail
    args=(--noplot)
    case {{ quote(mode) }} in
        list) args+=(--list) ;;
        test) args+=(--test) ;;
        measure) ;;
        *) echo "unknown acquisition mode" >&2; exit 2 ;;
    esac
    if [[ -n {{ quote(sample) }} ]]; then args+=(--save-baseline {{ quote(sample) }}); fi
    package=la-stack
    if [[ {{ quote(suite) }} == vs_linalg ]]; then package=la-stack-comparison; fi
    cargo bench --locked --target {{ host_target }} -p "$package" --features {{ quote(features) }} --bench {{ quote(suite) }} -- {{ quote(filter) }} "${args[@]}"

peers mode="measure" sample="": (native "vs_linalg" '^d(2|3|4|5|8|16|32|64)/(nalgebra|faer)_lu_solve$' "bench" mode sample)

public-default mode="measure" sample="": (native "vs_linalg" '^d[0-9]+/la_stack_' "bench" mode sample) (native "angle" '/(prepared_vector|borrowed_slice|construct_then_angle)$' "bench" mode sample) (native "gram" '/la_stack$' "bench" mode sample) (native "interval" "" "bench" mode sample) (native "linear_form" "" "bench" mode sample) (native "public_api" "" "bench" mode sample)

public-exact mode="measure" sample="": (native "exact" '^(canonical_conversion_|exact_(d[2-5]|random_corpus_d[2-5]|near_singular_3x3|large_entries_3x3|hilbert_[45]x[45])/|rational_input_.*d[2-8]/(det_sign_row_cleared_bareiss|det_row_cleared_bareiss|solve_row_cleared_bareiss)$)' "bench,exact" mode sample) (native "rational_solve" '/(construct|prepared|adapter)$' "bench,exact" mode sample) (native "public_api" '^public_exact_' "bench,exact" mode sample)
