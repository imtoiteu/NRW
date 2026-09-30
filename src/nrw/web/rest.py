"""Optional REST API (FastAPI) for automated processing.

    POST /api/v1/extract   multipart: file (.s2p/.csv) + form fields  -> JSON results + diagnostics
                           (optional file2 + thickness2_mm: two-thickness joint solution)
    GET  /api/v1/health    -> {"status": "ok", "software": ...}

Authentication: ``Authorization: Bearer <token>`` where the token is ``NRW_API_TOKEN`` if set,
otherwise ``NRW_APP_PASSWORD``. Failed attempts are rate-limited per client IP exactly as in
the web UI. Uses the same analysis layer as the browser application (no separate numerics).

Start with scripts/run_api.sh (default 127.0.0.1:8502).
"""
from __future__ import annotations

import json
import os

from fastapi import FastAPI, File, Form, Header, HTTPException, Request, UploadFile
from fastapi.responses import JSONResponse

from ..io import _jsonable, describe, load_measurement_bytes
from . import analysis as A
from . import export as X
from . import security
from .version import software_id

MAX_BYTES = 50 * 1024 * 1024
app = FastAPI(title="NRW material parameter extraction API", version="1",
              docs_url="/api/v1/docs", openapi_url="/api/v1/openapi.json", redoc_url=None)


def _token() -> str | None:
    return os.environ.get("NRW_API_TOKEN") or security.configured_password()


def _auth(request: Request, authorization: str | None):
    tok = _token()
    if not tok:
        raise HTTPException(503, "server not configured: set NRW_API_TOKEN or NRW_APP_PASSWORD")
    client = request.client.host if request.client else "unknown"
    cand = authorization[7:] if authorization and authorization.lower().startswith("bearer ") else ""
    ok, msg = security.attempt_login(cand, "api:" + client, expected=tok)
    if not ok:
        raise HTTPException(401 if "Incorrect" in msg else 429, msg)


@app.get("/api/v1/health")
def health():
    return {"status": "ok", "software": software_id()}


@app.post("/api/v1/extract")
async def extract(request: Request,
                  file: UploadFile = File(...),
                  thickness_mm: float = Form(...),
                  geometry: str = Form("coax"),
                  waveguide_a_mm: float | None = Form(None),
                  method: str = Form("auto"),
                  assume_nonmagnetic: bool = Form(False),
                  branch_method: str = Form("anchored"),
                  f_min_ghz: float | None = Form(None),
                  f_max_ghz: float | None = Form(None),
                  offset1_mm: float = Form(0.0),
                  offset2_mm: float = Form(0.0),
                  noise_sigma: float = Form(1e-3),
                  freq_unit: str | None = Form(None),
                  kk: bool = Form(True),
                  file2: UploadFile | None = File(None),
                  thickness2_mm: float | None = Form(None),
                  offset1_2_mm: float = Form(0.0),
                  offset2_2_mm: float = Form(0.0),
                  authorization: str | None = Header(None)):
    _auth(request, authorization)
    content = await file.read()
    content2 = await file2.read() if file2 is not None else None
    if len(content) > MAX_BYTES or (content2 is not None and len(content2) > MAX_BYTES):
        raise HTTPException(413, "file too large")
    if (content2 is None) != (thickness2_mm is None):
        raise HTTPException(422, "two-thickness mode needs both file2 and thickness2_mm")
    try:
        data = load_measurement_bytes(file.filename or "upload.s2p", content, freq_unit=freq_unit)
        cfg = A.AnalysisConfig(thickness_mm=thickness_mm, geometry=geometry, waveguide_a_mm=waveguide_a_mm,
                               f_min_ghz=f_min_ghz, f_max_ghz=f_max_ghz, method=method,
                               assume_nonmagnetic=assume_nonmagnetic, branch_method=branch_method,
                               offset1_mm=offset1_mm, offset2_mm=offset2_mm, noise_sigma=noise_sigma, kk=kk,
                               two_thickness=content2 is not None, thickness2_mm=thickness2_mm,
                               offset1_2_mm=offset1_2_mm, offset2_2_mm=offset2_2_mm)
        data2 = (load_measurement_bytes(file2.filename or "upload2.s2p", content2, freq_unit=freq_unit)
                 if content2 is not None else None)
        out = A.run_analysis(data, cfg, data_summary=describe(data), data2=data2)
    except ValueError as e:
        raise HTTPException(422, str(e)) from e
    body = {"metadata": X.metadata(out),
            "results": {c: out.table[c].tolist() for c in out.table.columns},
            "diagnostics": {c: out.diagnostics[c].tolist() for c in out.diagnostics.columns}}
    return JSONResponse(json.loads(json.dumps(_jsonable(body))))
