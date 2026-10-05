#!/usr/bin/env bash
set -euo pipefail
export CRITERION_HOME="$PWD/target/toolchain-audit/repeat"
for toolchain in 1.99.0 1.98.1; do
	if [[ "$toolchain" == 1.99.0 ]]; then
		comparison=(--save-baseline rust-1.99.0)
	else
		comparison=(--baseline rust-1.99.0)
	fi
	for suite in vs_linalg exact interval; do
		cargo +"$toolchain" bench --locked --ignore-rust-version \
			--target-dir "target/toolchain-$toolchain" --workspace \
			--features bench,exact --bench "$suite" -- '^(d3/la_stack_dot|d3/la_stack_ldlt|d3/la_stack_ldlt_solve|d3/la_stack_lu|d3/la_stack_norm2|d3/la_stack_norm2_sq|d4/la_stack_dot|d4/la_stack_ldlt|d4/la_stack_ldlt_solve|d4/la_stack_norm2|d4/la_stack_norm2_scenario_sparse|d5/la_stack_ldlt|d5/la_stack_ldlt_solve|d5/la_stack_lu|d5/la_stack_lu_solve|interval_det_sign/d4_inconclusive_lifted|rational_input_d5/solve_row_cleared_bareiss)$' \
			--noplot --sample-size 100 --warm-up-time 2 --measurement-time 5 \
			"${comparison[@]}"
	done
done
