#!/usr/bin/env bash
# Run prebuilt, correctness-checked binaries serially. See the parent report.
set -euo pipefail

study=${1:?directory containing saved binaries and receiving logs}
upstream=${2:?la-stack working directory}
downstream=${3:?isolated Delaunay working directory}
phase=${4:-measure}

predicates='predicates/(hot|exact_fallback)/insphere_lifted_[2-5]d'
realization='single_shared_vertex|realization_validation/[2-5]d'

if [[ "$phase" == measure || "$phase" == downstream ]]; then
	if [[ "$phase" == measure ]]; then
		cd "$upstream"
		for suite in interval linear; do
			"$study/$suite-before" --bench --warm-up-time 0.2 --measurement-time 1 \
				--sample-size 100 --nresamples 10000 --save-baseline issue-before --noplot \
				>"$study/final-$suite-before.log" 2>&1
		done
		for suite in interval linear; do
			"$study/$suite-after" --bench --warm-up-time 0.2 --measurement-time 1 \
				--sample-size 100 --nresamples 10000 --baseline issue-before --noplot \
				>"$study/final-$suite-after.log" 2>&1
		done
	fi
	cd "$downstream"
	for suite in predicates realization; do
		filter=$predicates
		if [[ "$suite" == realization ]]; then filter=$realization; fi
		"$study/$suite-private" --bench "$filter" --warm-up-time 1 --measurement-time 3 \
			--sample-size 100 --nresamples 10000 --save-baseline private-final --noplot \
			>"$study/final-$suite-private.log" 2>&1
		"$study/$suite-after" --bench "$filter" --warm-up-time 1 --measurement-time 3 \
			--sample-size 100 --nresamples 10000 --baseline private-final --noplot \
			>"$study/final-$suite-after.log" 2>&1
	done
elif [[ "$phase" == repeat ]]; then
	cd "$downstream"
	for suite in predicates realization; do
		filter=$predicates
		if [[ "$suite" == realization ]]; then filter=$realization; fi
		"$study/$suite-after" --bench "$filter" --warm-up-time 1 --measurement-time 3 \
			--sample-size 100 --nresamples 10000 --save-baseline adopted-repeat --noplot \
			>"$study/repeat-$suite-after.log" 2>&1
		"$study/$suite-private" --bench "$filter" --warm-up-time 1 --measurement-time 3 \
			--sample-size 100 --nresamples 10000 --baseline adopted-repeat --noplot \
			>"$study/repeat-$suite-private.log" 2>&1
	done
elif [[ "$phase" == confirm ]]; then
	cd "$upstream"
	for suite in interval linear; do
		filter='interval_scalar/(zero/(add|multiply)|cancellation/multiply)'
		if [[ "$suite" == linear ]]; then
			filter='linear_form_d2_batch/(origin|translated)/construction'
		fi
		"$study/$suite-after" --bench "$filter" --warm-up-time 0.2 --measurement-time 1 \
			--sample-size 100 --nresamples 10000 --save-baseline adopted-confirm --noplot \
			>"$study/confirm-$suite-after.log" 2>&1
		"$study/$suite-before" --bench "$filter" --warm-up-time 0.2 --measurement-time 1 \
			--sample-size 100 --nresamples 10000 --baseline adopted-confirm --noplot \
			>"$study/confirm-$suite-before.log" 2>&1
	done
	cd "$downstream"
	filter='single_shared_vertex/2d|realization_validation/(2d|5d)'
	for order in forward reverse; do
		first=private
		second=after
		if [[ "$order" == reverse ]]; then
			first=after
			second=private
		fi
		"$study/realization-$first" --bench "$filter" --warm-up-time 3 --measurement-time 5 \
			--sample-size 100 --nresamples 10000 --save-baseline "confirm-$order" --noplot \
			>"$study/confirm-$order-$first.log" 2>&1
		"$study/realization-$second" --bench "$filter" --warm-up-time 3 --measurement-time 5 \
			--sample-size 100 --nresamples 10000 --baseline "confirm-$order" --noplot \
			>"$study/confirm-$order-$second.log" 2>&1
		mkdir -p "$study/confirm-$order-data"
		cp -R "$downstream/target/criterion/." "$study/confirm-$order-data/"
	done
elif [[ "$phase" == control ]]; then
	cd "$downstream"
	"$study/realization-after" --bench 'realization_validation/3d/20v$' \
		--warm-up-time 3 --measurement-time 20 --sample-size 200 --nresamples 10000 \
		--save-baseline control-3d --noplot >"$study/control-3d-after.log" 2>&1
	"$study/realization-private" --bench 'realization_validation/3d/20v$' \
		--warm-up-time 3 --measurement-time 20 --sample-size 200 --nresamples 10000 \
		--baseline control-3d --noplot >"$study/control-3d-private.log" 2>&1
	mkdir -p "$study/control-3d-data"
	cp -R "$downstream/target/criterion/." "$study/control-3d-data/"
elif [[ "$phase" == flat ]]; then
	cd "$downstream"
	for order in forward reverse; do
		first=private
		second=after
		if [[ "$order" == reverse ]]; then
			first=after
			second=private
		fi
		"$study/realization-flat-$first" --bench 'single_shared_vertex/2d$' \
			--warm-up-time 3 --measurement-time 5 --sample-size 100 --nresamples 10000 \
			--save-baseline "flat-$order" --noplot >"$study/flat-$order-$first.log" 2>&1
		"$study/realization-flat-$second" --bench 'single_shared_vertex/2d$' \
			--warm-up-time 3 --measurement-time 5 --sample-size 100 --nresamples 10000 \
			--baseline "flat-$order" --noplot >"$study/flat-$order-$second.log" 2>&1
		mkdir -p "$study/flat-$order-data"
		cp -R "$downstream/target/criterion/." "$study/flat-$order-data/"
	done
else
	echo "unknown phase: $phase" >&2
	exit 2
fi
