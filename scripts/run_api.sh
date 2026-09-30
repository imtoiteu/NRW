#!/usr/bin/env bash
# Optional REST API (secondary to the browser UI). Default: localhost only.
#   export NRW_APP_PASSWORD='...'      # or a separate export NRW_API_TOKEN='...'
#   NRW_API_ADDRESS=0.0.0.0 NRW_API_PORT=8502 ./scripts/run_api.sh    # to expose remotely
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if [[ -z "${NRW_API_TOKEN:-${NRW_APP_PASSWORD:-}}" ]]; then
  echo "ERROR: set NRW_API_TOKEN or NRW_APP_PASSWORD" >&2; exit 2
fi
exec "${NRW_PYTHON:-$ROOT/.venv/bin/python}" -m uvicorn nrw.web.rest:app --app-dir "$ROOT/src" \
  --host "${NRW_API_ADDRESS:-127.0.0.1}" --port "${NRW_API_PORT:-8502}"
