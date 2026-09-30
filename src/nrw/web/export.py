"""Export of results: CSV, JSON, diagnostics, PNG figures, HTML report and a ZIP bundle.

Everything is built in memory (bytes) and handed to the browser for download; nothing is
written to a web-served location.
"""
from __future__ import annotations

import base64
import html
import io
import json
import zipfile
from dataclasses import asdict


from ..io import _jsonable
from .analysis import BRANCH_METHODS, METHODS
from .version import software_id

CONVENTION = ("exp(+jωt) time dependence; εr = ε′ − jε″ and μr = μ′ − jμ″, so a passive medium has "
              "ε″ ≥ 0 and μ″ ≥ 0 (Im εr ≤ 0, Im μr ≤ 0). Column *_imag = Im(x) ≤ 0 for passive media; "
              "*_loss = −Im(x) ≥ 0. No value is clipped or smoothed.")
KK_DISCLAIMER = ("The Kramers–Kronig check is a finite-band consistency diagnostic. It is not a mathematical "
                 "proof of global Kramers–Kronig compliance, which cannot be established from band-limited data.")


def results_csv(outcome) -> bytes:
    return outcome.table.to_csv(index=False, float_format="%.10g").encode()


def diagnostics_csv(outcome) -> bytes:
    return outcome.diagnostics.to_csv(index=False, float_format="%.10g").encode()


def metadata(outcome) -> dict:
    res = outcome.result
    cfg = outcome.config
    f = res.frequency
    return {
        "software": software_id(),
        "source_file": outcome.source,
        "analysis_timestamp_utc": outcome.timestamp_utc,
        "runtime_s": round(outcome.runtime_s, 3),
        "convention": CONVENTION,
        "configuration": asdict(cfg),
        "method_used": res.method,
        "method_label": METHODS.get(cfg.method, (cfg.method,))[0],
        "branch_method": res.branch_method,
        "branch_method_label": BRANCH_METHODS.get(res.branch_method, (res.branch_method,))[0],
        "frequency_interval_hz": [float(f[0]), float(f[-1])],
        "n_frequencies": int(f.size),
        "geometry": res.metadata.get("geometry"),
        "thickness_m": res.metadata.get("thickness_m"),
        "assessment": outcome.assessment.to_dict(),
        "warnings": list(res.warnings),
        "physics_checks": res.passivity.summary() if res.passivity is not None else None,
        "causality": {q: {"passed": r.passed, "statement": r.statement, "assumptions": r.assumptions,
                          "max_abs_normalized_residual": r.max_abs_normalized_residual,
                          "rms_normalized_residual": r.rms_normalized_residual, "threshold": r.threshold}
                      for q, r in (res.causality or {}).items()},
        "kk_disclaimer": KK_DISCLAIMER,
        "uncertainty": ({"spec": res.uncertainty["spec"], "n_failed_draws": res.uncertainty["n_failed_draws"]}
                        if res.uncertainty else None),
        "input_summary": outcome.data_summary,
    }


def results_json(outcome) -> bytes:
    doc = {"metadata": metadata(outcome), "data": {c: outcome.table[c].tolist() for c in outcome.table.columns}}
    return json.dumps(_jsonable(doc), indent=1).encode()


def _img(png: bytes, alt: str) -> str:
    return (f'<figure><img alt="{html.escape(alt)}" src="data:image/png;base64,'
            f'{base64.b64encode(png).decode()}"/><figcaption>{html.escape(alt)}</figcaption></figure>')


def report_html(outcome, figures: dict) -> bytes:
    """Self-contained HTML report (inline CSS and base64 PNGs)."""
    md = metadata(outcome)
    a = outcome.assessment
    res = outcome.result
    color = {"GREEN": "#1a7f37", "YELLOW": "#9a6700", "RED": "#cf222e"}[a.level]
    rows = "".join(
        f"<tr><td>{html.escape(r['severity'])}</td><td>{html.escape(r['rule'])}</td><td>{html.escape(r['message'])}"
        f"</td><td>{'' if r['n_points'] is None else r['n_points']}</td></tr>" for r in a.reasons) or \
        "<tr><td colspan=4>no diagnostic fired</td></tr>"
    checks = ""
    if res.passivity is not None:
        checks = "".join(f"<tr><td>{html.escape(c.name)}</td><td>{c.status}</td><td>{html.escape(c.description)}</td>"
                         f"<td>{c.n_violations}</td></tr>" for c in res.passivity.checks)
    kk = "".join(f"<li><b>{q}</b>: {html.escape(r.statement)}</li>" for q, r in (res.causality or {}).items()) \
        or "<li>not run / not available</li>"
    warns = "".join(f"<li>{html.escape(w)}</li>" for w in res.warnings) or "<li>none</li>"
    cfg = md["configuration"]
    sw = md["software"]
    info = [
        ("Source file", md["source_file"]), ("Analysis time (UTC)", md["analysis_timestamp_utc"]),
        ("Software", f"nrw {sw['package_version']}, commit {sw['commit']}" + (" (uncommitted changes)" if sw["dirty"] else "")),
        ("Sample thickness", f"{cfg['thickness_mm']:g} mm"),
        ("Geometry", "coaxial line (TEM)" if cfg["geometry"] == "coax" else f"rectangular waveguide TE10, a = {cfg['waveguide_a_mm']:g} mm"),
        ("Method", f"{md['method_label']} → {md['method_used']}"),
        ("Branch method", md["branch_method_label"]),
        ("Frequency interval", f"{md['frequency_interval_hz'][0]/1e9:.6g}–{md['frequency_interval_hz'][1]/1e9:.6g} GHz "
                               f"({md['n_frequencies']} points)"),
        ("Reference-plane offsets", f"port 1: {cfg['offset1_mm']:g} mm, port 2: {cfg['offset2_mm']:g} mm"),
        ("Assumed S-parameter noise σ", f"{cfg['noise_sigma']:g}"),
        ("Uncertainty", "Monte Carlo: " + json.dumps(md["uncertainty"]["spec"]) if md["uncertainty"] else "not estimated"),
    ]
    info_rows = "".join(f"<tr><th>{html.escape(k)}</th><td>{html.escape(str(v))}</td></tr>" for k, v in info)
    figs = "".join(_img(v, k) for k, v in figures.items())
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>NRW analysis report</title>
<style>
body{{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;max-width:980px;margin:24px auto;padding:0 16px;color:#1f2328;background:#fff}}
h1{{font-size:22px;margin-bottom:4px}} h2{{font-size:17px;margin-top:28px;border-bottom:1px solid #d0d7de;padding-bottom:4px}}
table{{border-collapse:collapse;width:100%;font-size:13px}} td,th{{border:1px solid #d0d7de;padding:4px 8px;text-align:left;vertical-align:top}}
th{{background:#f6f8fa;width:28%}} .level{{display:inline-block;padding:4px 12px;border-radius:6px;color:#fff;background:{color};font-weight:600}}
figure{{margin:12px 0}} img{{max-width:100%;border:1px solid #eaeef2}} figcaption{{font-size:12px;color:#57606a}}
.note{{font-size:13px;color:#57606a}}
</style></head><body>
<h1>NRW Material Parameter Extraction — analysis report</h1>
<p class="note">Complex permittivity and permeability retrieved from S-parameters with the tested <code>nrw</code> package.</p>
<h2>Summary</h2>
<p>Diagnostic status: <span class="level">{a.level}</span></p>
<p class="note">The status summarises documented diagnostic rules (src/nrw/assessment.py, rules v{a.rules_version}). It is not an
accuracy estimate: GREEN means that no diagnostic fired, not that the values are proven correct. Systematic errors
such as a wrong thickness, air gaps or calibration errors are not detectable by these diagnostics.</p>
<table>{info_rows}</table>
<h2>Convention</h2><p>{html.escape(CONVENTION)}</p>
<h2>Diagnostic findings</h2>
<table><tr><th>severity</th><th>rule</th><th>finding</th><th>points</th></tr>{rows}</table>
<h2>Warnings from the extraction</h2><ul>{warns}</ul>
<h2>Physical checks</h2>
<table><tr><th>check</th><th>status</th><th>criterion</th><th>violations</th></tr>{checks}</table>
<h2>Finite-band causality consistency</h2><p class="note">{html.escape(KK_DISCLAIMER)}</p><ul>{kk}</ul>
<h2>Figures</h2>{figs}
<p class="note">Generated by nrw {html.escape(sw['package_version'])} (commit {html.escape(sw['commit'])}).</p>
</body></html>"""
    return doc.encode()


def zip_bundle(outcome, figures: dict, input_name: str | None = None, input_bytes: bytes | None = None,
               comparison: dict | None = None) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("extracted_parameters.csv", results_csv(outcome))
        z.writestr("extracted_parameters.json", results_json(outcome))
        z.writestr("diagnostics.csv", diagnostics_csv(outcome))
        z.writestr("analysis_report.html", report_html(outcome, figures))
        z.writestr("configuration.json", json.dumps(_jsonable(metadata(outcome)), indent=1))
        for name, png in figures.items():
            z.writestr(f"figures/{name}", png)
        if comparison:
            for name, o in comparison.items():
                if hasattr(o, "table"):
                    z.writestr(f"comparison/{name}.csv", o.table.to_csv(index=False, float_format="%.10g"))
        if input_name and input_bytes is not None:
            z.writestr(f"input/{input_name}", input_bytes)
        z.writestr("README.txt", "NRW analysis export\n\n" + CONVENTION + "\n\n" + KK_DISCLAIMER + "\n\n"
                   "extracted_parameters.csv  per-frequency results and status\n"
                   "diagnostics.csv           Gamma, T, gamma, z, branch, Cramer-Rao bounds, solver data, flags\n"
                   "analysis_report.html      self-contained report\n"
                   "configuration.json        settings, software commit, assessment, warnings\n")
    return buf.getvalue()

