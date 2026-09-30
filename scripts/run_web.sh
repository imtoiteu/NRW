#!/usr/bin/env bash
# Launch the NRW web application.
#
#   export NRW_APP_PASSWORD='a strong password'   # required, never stored in Git
#   export NRW_APP_PORT=8501                       # optional (default 8501)
#   export NRW_APP_ADDRESS=0.0.0.0                 # optional (default 0.0.0.0 = all interfaces)
#   ./scripts/run_web.sh
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if [[ -z "${NRW_APP_PASSWORD:-}" ]]; then
  echo "ERROR: NRW_APP_PASSWORD is not set. Example: export NRW_APP_PASSWORD='choose-a-strong-password'" >&2
  exit 2
fi
if (( ${#NRW_APP_PASSWORD} < 8 )); then
  echo "ERROR: NRW_APP_PASSWORD must have at least 8 characters." >&2
  exit 2
fi
PORT="${NRW_APP_PORT:-8501}"
ADDR="${NRW_APP_ADDRESS:-0.0.0.0}"
PY="${NRW_PYTHON:-$ROOT/.venv/bin/python}"
if [[ ! -x "$PY" ]]; then
  echo "ERROR: Python interpreter $PY not found. Create the venv first (see docs/web_app.md)." >&2
  exit 2
fi
echo "Starting NRW web app on http://${ADDR}:${PORT}  (Ctrl+C to stop)"
exec "$PY" -m streamlit run "$ROOT/app.py" --server.address "$ADDR" --server.port "$PORT"
