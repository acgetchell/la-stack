# Captured native commands run in the measured checkout with its phase environment.
# Pass "list" to discover Criterion IDs without collecting samples.

# Cargo runs benches from each package root, including benches/comparison.
# Resolve the shared output tree before those working directories change.
export CRITERION_HOME := absolute_path(env("CRITERION_HOME", "target/criterion"))

# Compare native targets on the measured host regardless of a caller's target override.
host_target := `rustc --print host-tuple`

# Discovery and integration gates must both compile without the generated README.
[private]
[script("uv", "run", "--locked", "--no-sync", "python", "-I")]
_readme-input-check:
    import re
    from pathlib import Path

    root = Path("target/build").resolve(strict=True)
    profiles = set()
    files = sorted(root.rglob("la_stack-*.d"))
    if not files:
        raise SystemExit("missing la-stack compiler dependency files")
    for path in files:
        if not re.fullmatch(r"la_stack-[0-9a-f]+\.d", path.name):
            raise SystemExit(f"unexpected la-stack dependency filename: {path}")
        if path.is_symlink() or path.resolve(strict=True) != path:
            raise SystemExit(f"aliased compiler dependency file: {path}")
        data = path.read_bytes().lower().replace(b"\\", b"/")
        if b"src/lib.rs" not in data:
            raise SystemExit(f"missing library input in compiler dependency file: {path}")
        if b"readme.md" in data:
            raise SystemExit("current compilation reads README.md; full acquisition requires documentation-only inclusion; use non-exact for historical current tags")
        profiles.update(set(path.relative_to(root).parts) & {"debug", "release"})
    if profiles != {"debug", "release"}:
        raise SystemExit("missing benchmark-discovery or correctness-gate compiler dependencies")

baseline-vs_linalg mode="measure":
    cargo bench --locked --target {{ host_target }} -p la-stack-comparison --features bench --bench vs_linalg -- --noplot {{ if mode == "list" { "--list" } else { "" } }}

column-basis mode="measure": (native "column_basis" "" "bench" mode)

common-public-default mode="measure" sample="": (native "vs_linalg" '^d[0-9]+/la_stack_' "bench" mode sample) (default-operations mode sample)

common-public-exact mode="measure" sample="": (native "exact" '^(canonical_conversion_|exact_(d[2-5]|random_corpus_d[2-5]|near_singular_3x3|large_entries_3x3|hilbert_[45]x[45])/|rational_input_.*d[2-8]/(det_sign_row_cleared_bareiss|det_row_cleared_bareiss|solve_row_cleared_bareiss)$)' "bench,exact" mode sample) (native "rational_solve" '/(construct|prepared|adapter)$' "bench,exact" mode sample) (native "public_api" '^public_exact_' "bench,exact" mode sample)

# LU factor-and-solve and its peers have their own paired acquisition unit.
core-default mode="measure": (native "vs_linalg" '^d[0-9]+/la_stack_(det.*|dot|inf_norm.*|ldlt.*|lu|lu_pivoting|lu_ill_conditioned|norm2.*|solve.*)$' "bench" mode) (default-operations mode)

current-vs_linalg mode="measure":
    cargo bench --locked --target {{ host_target }} -p la-stack-comparison --features bench --bench vs_linalg -- la_stack --noplot {{ if mode == "list" { "--list" } else { "" } }}

default-operations mode="measure" sample="": (native "angle" '/(prepared_vector|borrowed_slice|construct_then_angle)$' "bench" mode sample) (native "gram" '/la_stack$' "bench" mode sample) (native "interval" "" "bench" mode sample) (native "linear_form" "" "bench" mode sample) (native "public_api" "" "bench" mode sample)

gate:
    cargo test --locked --target {{ host_target }} --workspace --features bench,exact --test vs_linalg_inputs --test exact_bench_config

gate-column-basis:
    cargo test --locked --target {{ host_target }} --features bench --test column_basis

gate-column-basis-current: gate-column-basis _readme-input-check

gate-current: gate _readme-input-check

gate-interval-solve:
    cargo test --locked --target {{ host_target }} --features bench,exact --test interval_solve

gate-interval-solve-current: gate-interval-solve _readme-input-check

gate-polynomial:
    cargo test --locked --target {{ host_target }} --features bench,exact --test polynomial_bench_inputs

gate-polynomial-current: gate-polynomial _readme-input-check

gate-public: gate gate-column-basis gate-interval-solve gate-polynomial gate-rational-dot

gate-rational-dot:
    cargo test --locked --target {{ host_target }} --features bench,exact --test rational_dot

gate-rational-dot-current: gate-rational-dot _readme-input-check
interval-solve mode="measure": (native "interval_solve" "" "bench" mode)

la-stack mode="measure" sample="": (public-default mode sample) (public-exact mode sample)

lu-solve mode="measure": (native "vs_linalg" '^d(2|3|4|5|8|16|32|64)/(la_stack|nalgebra|faer)_lu_solve$' "bench" mode)

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

polynomial mode="measure": (native "polynomial" "" "bench,exact" mode)

public-default mode="measure" sample="": (common-public-default mode sample) (native "column_basis" "" "bench" mode sample) (native "interval_solve" "" "bench" mode sample)

public-exact mode="measure" sample="": (common-public-exact mode sample) (native "polynomial" "" "bench,exact" mode sample) (native "rational_dot" "" "bench,exact" mode sample)

rational-dot mode="measure": (native "rational_dot" "" "bench,exact" mode)
