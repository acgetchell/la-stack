#!/usr/bin/env bash
set -euo pipefail
export CRITERION_HOME="$PWD/target/toolchain-audit/criterion"
for toolchain in 1.98.1 1.99.0; do
	if [[ "$toolchain" == 1.98.1 ]]; then
		comparison=(--save-baseline rust-1.98.1)
	else
		comparison=(--baseline rust-1.98.1)
	fi
	for suite in vs_linalg exact interval linear_form gram; do
		case "$suite" in
		vs_linalg) filter='^(d[2345]/la_stack_(dot|norm2|norm2_sq|det|lu|ldlt|lu_solve|ldlt_solve)|d4/la_stack_norm2_scenario_.*)$' ;;
		exact) filter='^(exact_d[45]/(det_exact|solve_exact)|rational_input_d5/(det_row_cleared_bareiss|solve_row_cleared_bareiss))$' ;;
		gram) filter='^gram/4x4/.*/la_stack$' ;;
		*) filter='' ;;
		esac
		cargo +"$toolchain" bench --locked --ignore-rust-version \
			--target-dir "target/toolchain-$toolchain" --workspace \
			--features bench,exact --bench "$suite" -- "$filter" \
			--noplot --sample-size 100 --warm-up-time 1 --measurement-time 3 \
			"${comparison[@]}"
	done
done
