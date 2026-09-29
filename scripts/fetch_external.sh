#!/usr/bin/env bash
# Re-clone every third-party repository surveyed in docs/github_survey.md at the
# exact commit pinned in research/external_pins.csv, into research/external/<name>.
#
# research/external/ is gitignored: third-party code is NOT redistributed with
# this repository (several of the repositories carry no licence at all, which
# means all rights reserved). This script only reproduces the local checkout
# needed by research/adapters/*.
#
# Usage:  scripts/fetch_external.sh [name ...]     (default: all pins)
#         FORCE=1 scripts/fetch_external.sh         (re-fetch even if present)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PINS="$ROOT/research/external_pins.csv"
DEST="$ROOT/research/external"
mkdir -p "$DEST"

want=("$@")
fail=0
while IFS=, read -r name url commit; do
  [[ "$name" == "name" || -z "$name" ]] && continue
  if (( ${#want[@]} )); then
    skip=1; for w in "${want[@]}"; do [[ "$w" == "$name" ]] && skip=0; done
    (( skip )) && continue
  fi
  dir="$DEST/$name"
  if [[ -d "$dir/.git" && -z "${FORCE:-}" ]]; then
    have="$(git -C "$dir" rev-parse HEAD 2>/dev/null || true)"
    if [[ "$have" == "$commit" ]]; then
      echo "ok      $name @ ${commit:0:12}"; continue
    fi
  fi
  rm -rf "$dir"; mkdir -p "$dir"
  (
    cd "$dir"
    git init -q
    git remote add origin "$url.git"
    # GitHub serves reachable commits by SHA; shallow fetch of exactly the pin.
    if ! git fetch -q --depth 1 origin "$commit" 2>/dev/null; then
      git fetch -q origin            # fallback: full history, then check out pin
    fi
    git -c advice.detachedHead=false checkout -q "$commit"
  ) && got="$(git -C "$dir" rev-parse HEAD)" || got="FAILED"
  if [[ "$got" == "$commit" ]]; then
    echo "fetched $name @ ${commit:0:12}"
  else
    echo "ERROR   $name: expected $commit, got $got" >&2; fail=1
  fi
done < "$PINS"
exit $fail
