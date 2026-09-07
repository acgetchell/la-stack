#!/usr/bin/env bash
# Keep authentication out of command arguments and shell traces.
set +x
set -euo pipefail

if [[ "${ZIZMOR_OFFLINE:-false}" == true || "${ZIZMOR_NO_ONLINE_AUDITS:-false}" == true ]]; then
	echo "zizmor: offline audits requested; online audits are disabled."
	exec zizmor --offline --persona regular .github
fi

zizmor_token="${ZIZMOR_GITHUB_TOKEN:-${GH_TOKEN:-${GITHUB_TOKEN:-}}}"
if [[ -z "$zizmor_token" ]] && command -v gh >/dev/null; then
	if resolved_token="$(gh auth token 2>/dev/null)"; then
		zizmor_token="$resolved_token"
	fi
fi

if [[ -n "$zizmor_token" ]]; then
	echo "zizmor: running authenticated online audits (persona: regular)."
	# Normalize GH_TOKEN too: zizmor checks it before ZIZMOR_GITHUB_TOKEN.
	export GH_TOKEN="$zizmor_token" ZIZMOR_GITHUB_TOKEN="$zizmor_token"
	exec zizmor --persona regular .github
fi

echo "zizmor: no GitHub token available; using offline audits. Online findings are not checked."
exec zizmor --offline --persona regular .github
