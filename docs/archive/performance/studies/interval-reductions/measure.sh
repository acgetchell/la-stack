#!/usr/bin/env bash
# Run prebuilt, correctness-checked binaries serially. See the parent report.
set -euo pipefail

study=${1:?directory containing saved binaries and receiving logs}
upstream=${2:?la-stack working directory}
downstream=${3:?isolated Delaunay working directory}
phase=${4:-measure}

# Resolve every input against the invocation directory before any cd.
study=$(cd -- "$study" && pwd -P)
upstream=$(cd -- "$upstream" && pwd -P)
downstream=$(cd -- "$downstream" && pwd -P)
case "$phase" in
measure | downstream | repeat | confirm | control | flat) ;;
*)
	echo "unknown phase: $phase" >&2
	exit 2
	;;
esac

# Each invocation owns its logs and raw Criterion output. Later phases and
# reruns must not replace earlier new/change measurements or confidence bounds.
results=$(mktemp -d "$study/$phase.XXXXXX")
printf 'Measurement output: %s\n' "$results"

predicates='predicates/(hot|exact_fallback)/insphere_lifted_[2-5]d'
realization='single_shared_vertex|realization_validation/[2-5]d'

if [[ "$phase" == measure || "$phase" == downstream ]]; then
	if [[ "$phase" == measure ]]; then
		cd "$upstream"
		export CRITERION_HOME="$results/upstream-data"
		for suite in interval linear; do
			"$study/$suite-before" --bench --warm-up-time 0.2 --measurement-time 1 \
				--sample-size 100 --nresamples 10000 --save-baseline issue-before --noplot \
				>"$results/final-$suite-before.log" 2>&1
		done
		for suite in interval linear; do
			"$study/$suite-after" --bench --warm-up-time 0.2 --measurement-time 1 \
				--sample-size 100 --nresamples 10000 --baseline issue-before --noplot \
				>"$results/final-$suite-after.log" 2>&1
		done
	fi
	cd "$downstream"
	export CRITERION_HOME="$results/downstream-data"
	for suite in predicates realization; do
		filter=$predicates
		if [[ "$suite" == realization ]]; then filter=$realization; fi
		"$study/$suite-private" --bench "$filter" --warm-up-time 1 --measurement-time 3 \
			--sample-size 100 --nresamples 10000 --save-baseline private-final --noplot \
			>"$results/final-$suite-private.log" 2>&1
		"$study/$suite-after" --bench "$filter" --warm-up-time 1 --measurement-time 3 \
			--sample-size 100 --nresamples 10000 --baseline private-final --noplot \
			>"$results/final-$suite-after.log" 2>&1
	done
elif [[ "$phase" == repeat ]]; then
	cd "$downstream"
	export CRITERION_HOME="$results/downstream-repeat-data"
	for suite in predicates realization; do
		filter=$predicates
		if [[ "$suite" == realization ]]; then filter=$realization; fi
		"$study/$suite-after" --bench "$filter" --warm-up-time 1 --measurement-time 3 \
			--sample-size 100 --nresamples 10000 --save-baseline adopted-repeat --noplot \
			>"$results/repeat-$suite-after.log" 2>&1
		"$study/$suite-private" --bench "$filter" --warm-up-time 1 --measurement-time 3 \
			--sample-size 100 --nresamples 10000 --baseline adopted-repeat --noplot \
			>"$results/repeat-$suite-private.log" 2>&1
	done
elif [[ "$phase" == confirm ]]; then
	cd "$upstream"
	export CRITERION_HOME="$results/upstream-confirm-data"
	for suite in interval linear; do
		filter='interval_scalar/(zero/(add|multiply)|cancellation/multiply)'
		if [[ "$suite" == linear ]]; then
			filter='linear_form_d2_batch/(origin|translated)/construction'
		fi
		"$study/$suite-after" --bench "$filter" --warm-up-time 0.2 --measurement-time 1 \
			--sample-size 100 --nresamples 10000 --save-baseline adopted-confirm --noplot \
			>"$results/confirm-$suite-after.log" 2>&1
		"$study/$suite-before" --bench "$filter" --warm-up-time 0.2 --measurement-time 1 \
			--sample-size 100 --nresamples 10000 --baseline adopted-confirm --noplot \
			>"$results/confirm-$suite-before.log" 2>&1
	done
	cd "$downstream"
	filter='single_shared_vertex/2d|realization_validation/(2d|5d)'
	for order in forward reverse; do
		export CRITERION_HOME="$results/confirm-$order-data"
		first=private
		second=after
		if [[ "$order" == reverse ]]; then
			first=after
			second=private
		fi
		"$study/realization-$first" --bench "$filter" --warm-up-time 3 --measurement-time 5 \
			--sample-size 100 --nresamples 10000 --save-baseline "confirm-$order" --noplot \
			>"$results/confirm-$order-$first.log" 2>&1
		"$study/realization-$second" --bench "$filter" --warm-up-time 3 --measurement-time 5 \
			--sample-size 100 --nresamples 10000 --baseline "confirm-$order" --noplot \
			>"$results/confirm-$order-$second.log" 2>&1
	done
elif [[ "$phase" == control ]]; then
	cd "$downstream"
	export CRITERION_HOME="$results/control-3d-data"
	"$study/realization-after" --bench 'realization_validation/3d/20v$' \
		--warm-up-time 3 --measurement-time 20 --sample-size 200 --nresamples 10000 \
		--save-baseline control-3d --noplot >"$results/control-3d-after.log" 2>&1
	"$study/realization-private" --bench 'realization_validation/3d/20v$' \
		--warm-up-time 3 --measurement-time 20 --sample-size 200 --nresamples 10000 \
		--baseline control-3d --noplot >"$results/control-3d-private.log" 2>&1
elif [[ "$phase" == flat ]]; then
	cd "$downstream"
	for order in forward reverse; do
		export CRITERION_HOME="$results/flat-$order-data"
		first=private
		second=after
		if [[ "$order" == reverse ]]; then
			first=after
			second=private
		fi
		"$study/realization-flat-$first" --bench 'single_shared_vertex/2d$' \
			--warm-up-time 3 --measurement-time 5 --sample-size 100 --nresamples 10000 \
			--save-baseline "flat-$order" --noplot >"$results/flat-$order-$first.log" 2>&1
		"$study/realization-flat-$second" --bench 'single_shared_vertex/2d$' \
			--warm-up-time 3 --measurement-time 5 --sample-size 100 --nresamples 10000 \
			--baseline "flat-$order" --noplot >"$results/flat-$order-$second.log" 2>&1
	done
fi
