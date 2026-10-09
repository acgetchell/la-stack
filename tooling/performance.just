# Captured native commands run in the measured checkout with its phase environment.
# Pass "list" to discover Criterion IDs without collecting samples.

baseline-all mode="measure": (baseline-vs_linalg mode) (baseline-exact mode)

baseline-exact mode="measure":
    cargo bench --locked --features bench,exact --bench exact -- --noplot {{ if mode == "list" { "--list" } else { "" } }}

baseline-vs_linalg mode="measure":
    cargo bench --locked -p la-stack-comparison --features bench --bench vs_linalg -- --noplot {{ if mode == "list" { "--list" } else { "" } }}

current-all mode="measure": (current-vs_linalg mode) (current-exact mode)

current-exact mode="measure": (baseline-exact mode)

current-vs_linalg mode="measure":
    cargo bench --locked -p la-stack-comparison --features bench --bench vs_linalg -- la_stack --noplot {{ if mode == "list" { "--list" } else { "" } }}

gate:
    cargo test --locked --workspace --features bench,exact --test vs_linalg_inputs --test exact_bench_config
