#!/usr/bin/env bash
set -euo pipefail

# Run from the la-stack root, with an isolated Delaunay manifest as argument one.
manifest=${1:?Pass the isolated Delaunay Cargo.toml path}
mode=${2:-primary}
case "$mode" in
primary)
	export CRITERION_HOME="$PWD/target/toolchain-audit/downstream"
	toolchains=(1.98.1 1.99.0)
	baseline=rust-1.98.1
	warmup=1
	measurement=3
	;;
repeat)
	export CRITERION_HOME="$PWD/target/toolchain-audit/downstream-repeat"
	toolchains=(1.99.0 1.98.1)
	baseline=rust-1.99.0
	warmup=2
	measurement=5
	;;
*)
	printf 'Expected primary or repeat mode, got: %s\n' "$mode" >&2
	exit 2
	;;
esac

for toolchain in "${toolchains[@]}"; do
	if [[ "$toolchain" == "${toolchains[0]}" ]]; then
		comparison=(--save-baseline "$baseline")
	else
		comparison=(--baseline "$baseline")
	fi
	cargo +"$toolchain" bench --offline --locked --ignore-rust-version \
		--manifest-path "$manifest" \
		--config "patch.crates-io.la-stack.path=\"$PWD\"" \
		--target-dir "target/delaunay-$toolchain" --profile perf --bench math_kernels -- \
		'^math/(orientation/(well_conditioned|exact_nonzero|exact_zero)/[34]d|geometry/axis_simplex/(volume|circumradius)/3d)$' \
		--noplot --sample-size 100 --warm-up-time "$warmup" --measurement-time "$measurement" \
		"${comparison[@]}"
done
