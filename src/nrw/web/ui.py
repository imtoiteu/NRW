"""Streamlit user interface. Numerical work is delegated to nrw.web.analysis -> nrw core.

Layout: login gate -> header -> data source -> tabs
(Data · Configure & extract · Results · Diagnostics · Passivity & causality · Compare methods ·
Uncertainty · Export).
"""
from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

from ..io import describe, load_measurement_bytes
from . import analysis as A
from . import export as X
from . import plots as P
from . import security, state
from .version import software_id, software_label

APP_TITLE = "NRW Material Parameter Extraction Laboratory"
SUBTITLE = "Complex permittivity and permeability retrieval from measured S-parameters"
DEMO_DIR = Path(__file__).resolve().parents[3] / "data" / "synthetic"
WAVEGUIDES = {  # name: (a in mm, recommended TE10 band in GHz)
    "WR-284": (72.136, "2.60–3.95"), "WR-187": (47.549, "3.95–5.85"), "WR-137": (34.849, "5.85–8.20"),
    "WR-90": (22.860, "8.20–12.40"), "WR-62": (15.799, "12.40–18.00"),
}
LEVEL_UI = {"GREEN": ("success", "GREEN — extraction appears numerically and physically consistent"),
            "YELLOW": ("warning", "YELLOW — extraction completed, but one or more diagnostics require attention"),
            "RED": ("error", "RED — the dataset/configuration does not support a trustworthy extraction")}


# ------------------------------------------------------------------------------------------
def main():
    st.set_page_config(page_title=APP_TITLE, page_icon=None, layout="wide", initial_sidebar_state="expanded")
    ss = st.session_state
    state.init(ss)
    problem = security.password_problem()
    if problem:
        _header()
        st.error("The application is locked because it is not configured securely.")
        st.info(problem)
        st.stop()
    if not ss.get("authenticated"):
        _login()
        st.stop()
    _sidebar()
    _header()
    _data_source()
    if "data" not in ss:
        st.info("Upload an .s2p or .csv file, or load synthetic demonstration data, to begin.")
        return
    tabs = st.tabs(["Data", "Configure & extract", "Results", "Diagnostics", "Passivity & causality",
                    "Compare methods", "Uncertainty", "Export"])
    with tabs[0]:
        _tab_data()
    with tabs[1]:
        _tab_configure()
    outcome = _current_outcome()
    with tabs[2]:
        _tab_results(outcome)
    with tabs[3]:
        _tab_diagnostics(outcome)
    with tabs[4]:
        _tab_passivity(outcome)
    with tabs[5]:
        _tab_compare()
    with tabs[6]:
        _tab_uncertainty()
    with tabs[7]:
        _tab_export(outcome)


# ------------------------------------------------------------------------------------------
def _header():
    st.markdown(f"## {APP_TITLE}")
    st.markdown(f"{SUBTITLE}")
    st.caption("Upload → Validate → Configure → Extract → Diagnose → Export  ·  " + software_label())


def _client_id() -> str:
    try:
        ip = st.context.ip_address
    except Exception:  # noqa: BLE001
        ip = None
    return ip or "unknown-client"


def _login():
    _header()
    st.markdown("#### Sign in")
    with st.form("login", clear_on_submit=True):
        pw = st.text_input("Access password", type="password")
        ok = st.form_submit_button("Sign in")
    if ok:
        good, msg = security.attempt_login(pw, _client_id())
        if good:
            st.session_state["authenticated"] = True
            st.rerun()
        else:
            st.error(msg)
    st.caption("Uploads and processing are available only after sign-in. Data stay in your browser session "
               "on this server and are not shared with other users.")


def _sidebar():
    ss = st.session_state
    with st.sidebar:
        st.markdown("**Session**")
        st.caption(software_label())
        if "data" in ss:
            st.caption(f"Loaded: {ss.get('input_name', '')}")
        if st.button("Delete session data", help="Removes the uploaded file, results and prepared exports from "
                                                  "this session (server memory)."):
            state.clear_session_data(ss)
            ss.pop("fig_memo", None)
            st.rerun()
        if st.button("Sign out"):
            state.logout(ss)
            st.rerun()
        st.divider()
        st.markdown("**Convention**")
        st.caption("exp(+jωt): εr = ε′ − jε″, μr = μ′ − jμ″. A passive medium has ε″ ≥ 0 and μ″ ≥ 0. "
                   "Values are never clipped, smoothed or sign-flipped.")
        st.markdown("**Documentation**")
        st.caption("docs/web_app.md · docs/final_method.md · docs/benchmark_report.md")


# ------------------------------------------------------------------------------------------
def _load(name: str, content: bytes, freq_unit=None, demo_meta=None):
    ss = st.session_state
    try:
        data = load_measurement_bytes(name, content, freq_unit=freq_unit)
    except ValueError as e:
        st.error(f"Could not read {name}: {e}")
        return False
    state.clear_session_data(ss)
    ss["data"] = data
    ss["input_name"] = name
    ss["input_bytes"] = content
    ss["input_summary"] = describe(data)
    ss["demo_meta"] = demo_meta
    return True


def _data_source():
    ss = st.session_state
    with st.expander("Data source", expanded="data" not in ss):
        src = st.radio("Source", ["Upload measured data", "Synthetic demonstration data"], horizontal=True,
                       label_visibility="collapsed")
        if src == "Upload measured data":
            c1, c2 = st.columns([3, 1])
            with c2:
                unit = st.selectbox("CSV frequency unit (only if not in the header)",
                                    ["(from header)", "Hz", "kHz", "MHz", "GHz"])
            with c1:
                up = st.file_uploader("Touchstone .s2p or CSV", type=["s2p", "csv", "txt", "tsv"],
                                      help="Touchstone v1 (RI/MA/DB; Hz…GHz) or a CSV with a header, e.g. "
                                           "'freq_GHz,S11_real,S11_imag,S21_real,S21_imag' or "
                                           "'Frequency (GHz),S11_dB,S11_phase_deg,…'.")
            if up is not None and st.button("Load file", type="primary"):
                fu = None if unit.startswith("(") else unit.lower()
                if _load(up.name, up.getvalue(), fu):
                    st.rerun()
        else:
            try:
                manifest = json.loads((DEMO_DIR / "manifest.json").read_text())
            except OSError:
                manifest = {}
            if not manifest:
                st.warning("No demonstration data found (run examples/make_synthetic_s2p.py).")
                return
            name = st.selectbox("Synthetic demonstration data set", list(manifest),
                                format_func=lambda n: f"{n} — {manifest[n]['description']}")
            st.caption("**Synthetic demonstration data** generated by the independent forward model "
                       "(nrw.forward) with 1e-3 complex noise. Not an experimental measurement.")
            if st.button("Load demo data", type="primary"):
                if _load(name, (DEMO_DIR / name).read_bytes(), demo_meta=manifest[name]):
                    st.rerun()


# ------------------------------------------------------------------------------------------
def _tab_data():
    ss = st.session_state
    data, summ = ss["data"], ss["input_summary"]
    if ss.get("demo_meta"):
        st.info("Synthetic demonstration data — not an experimental measurement. "
                f"Generating model: ε = {ss['demo_meta']['eps_model']}, μ = {ss['demo_meta']['mu_model']}, "
                f"thickness {ss['demo_meta']['thickness_mm']} mm.")
    c = st.columns(4)
    c[0].metric("Points", summ["n_points"])
    c[1].metric("Frequency range", f"{summ['f_min_hz']/1e9:.4g}–{summ['f_max_hz']/1e9:.4g} GHz")
    c[2].metric("Median spacing", f"{summ['df_median_hz']/1e6:.4g} MHz")
    c[3].metric("Uniform grid", "yes" if summ["uniform_grid"] else "no")
    rows = [("File", summ["source"]),
            ("Parameters present", ", ".join(k.upper() for k in summ["measured"])),
            ("Reference impedance (file)", f"{summ['z0_ohm']:g} Ω"),
            ("Frequency spacing min / max", f"{summ['df_min_hz']/1e6:.6g} / {summ['df_max_hz']/1e6:.6g} MHz"),
            ("|S11| min / max", f"{summ['s11_abs_min']:.4g} / {summ['s11_abs_max']:.4g}"),
            ("|S21| min / max", f"{summ['s21_abs_min']:.4g} / {summ['s21_abs_max']:.4g}")]
    st.table(pd.DataFrame(rows, columns=["property", "value"]).set_index("property"))
    problems = _validation_messages(summ)
    for p in problems:
        st.error(p)
    for n in summ["notes"]:
        st.caption(f"Parser note: {n}")
    unwrap = st.toggle("Show unwrapped phase (display only)", value=False)
    c1, c2 = st.columns(2)
    c1.plotly_chart(P.sparam_magnitude(data), width="stretch")
    c2.plotly_chart(P.sparam_phase(data, unwrap=unwrap), width="stretch")
    if not problems:
        st.success("Input validation passed. Continue with 'Configure & extract'.")


def _validation_messages(summ) -> list[str]:
    out = []
    if not summ["finite"]:
        out.append("The file contains non-finite values (NaN/Inf).")
    if not summ["strictly_increasing"]:
        out.append("Frequencies are not strictly increasing (duplicates or unsorted rows). Fix the file; the "
                   "application does not re-sort data silently.")
    if summ["n_points"] < 5:
        out.append("At least 5 frequency points are required.")
    if "s11" not in summ["measured"] or "s21" not in summ["measured"]:
        out.append("S11 and S21 are required.")
    if summ["s11_abs_max"] > 1.05 or summ["s21_abs_max"] > 1.05:
        out.append("|S| exceeds 1.05 somewhere: check calibration/normalisation (a passive sample cannot have gain). "
                   "Extraction is still possible but will be flagged.")
    return out


# ------------------------------------------------------------------------------------------
def _tab_configure():
    ss = st.session_state
    data, summ = ss["data"], ss["input_summary"]
    if _validation_messages(summ) and (not summ["finite"] or not summ["strictly_increasing"]
                                       or summ["n_points"] < 5):
        st.error("Input validation failed (see the Data tab). Extraction is disabled.")
        return
    demo = ss.get("demo_meta") or {}
    fmin, fmax = summ["f_min_hz"] / 1e9, summ["f_max_hz"] / 1e9
    with st.form("config"):
        st.markdown("##### Sample and measurement cell")
        c = st.columns(3)
        thickness = c[0].number_input("Sample thickness (mm)", min_value=0.0, value=float(demo.get("thickness_mm", 2.0)),
                                      step=0.01, format="%.4f", help="Length along the propagation direction.")
        geom = c[1].selectbox("Measurement cell", ["Coaxial air line (TEM)", "Rectangular waveguide (TE10)"])
        wg = c[2].selectbox("Waveguide size (TE10)", list(WAVEGUIDES) + ["custom"], index=3,
                            format_func=lambda k: f"{k} (a = {WAVEGUIDES[k][0]} mm, {WAVEGUIDES[k][1]} GHz)"
                            if k in WAVEGUIDES else "custom a")
        a_custom = c[2].number_input("Custom broad-wall width a (mm)", min_value=0.0, value=22.86, step=0.01)
        st.markdown("##### Frequency window")
        win = st.slider("Analysed band (GHz)", min_value=float(fmin), max_value=float(fmax),
                        value=(float(fmin), float(fmax)), step=max((fmax - fmin) / 1000, 1e-6), format="%.4f")
        st.markdown("##### Extraction method")
        method = st.radio("Method", list(A.METHODS), format_func=lambda k: A.METHODS[k][0], horizontal=False)
        nonmag = st.checkbox("The sample is known to be non-magnetic (μr = 1)", value=False,
                             help="Used by AUTO. Only tick this if μr = 1 is justified independently; a μ = 1 method "
                                  "on a magnetic sample returns a biased permittivity.")
        with st.expander("What do the methods do?"):
            for k, (lab, desc) in A.METHODS.items():
                st.markdown(f"**{lab}** — {desc}")
        st.markdown("##### Phase branch")
        branch = st.selectbox("Branch selection", list(A.BRANCH_METHODS), format_func=lambda k: A.BRANCH_METHODS[k][0])
        prior = st.number_input("Prior estimate of εr·μr (only for 'From a prior estimate')", min_value=0.0, value=0.0)
        with st.expander("What do the branch methods do?"):
            for k, (lab, desc) in A.BRANCH_METHODS.items():
                st.markdown(f"**{lab}** — {desc}")
            st.caption("numpy.unwrap alone never guarantees the physical branch (docs/benchmark_report.md §4).")
        st.markdown("##### Reference planes")
        c = st.columns(3)
        off1 = c[0].number_input("Offset port 1 (mm)", value=0.0, step=0.01, format="%.4f",
                                 help="Air-filled distance from the port-1 calibration plane to the sample face. "
                                      "De-embedded exactly: S11/R1², S21/(R1R2), R = exp(−γ0 d).")
        off2 = c[1].number_input("Offset port 2 (mm)", value=0.0, step=0.01, format="%.4f")
        holder = c[2].number_input("Holder length (mm, optional, NIST μ = 1 method only)", min_value=0.0, value=0.0,
                                   help="If > 0 and the Baker-Jarvis/NIST method is used, the reference-plane "
                                        "invariant form needs only L_air − L; offsets are then ignored.")
        with st.expander("Advanced"):
            noise = st.number_input("Assumed S-parameter noise σ (rms, linear)", min_value=1e-6, value=1e-3,
                                    format="%.2e", help="Used for the Cramér–Rao identifiability flags, the noise "
                                                        "floor and noise-aware physical checks.")
            avg = st.checkbox("Average ports for the NRW start value (reciprocal, symmetric sample)", value=True)
            reg = st.number_input("Regularisation weight (frequency-regularised method only)", min_value=0.0,
                                  value=1e-2, format="%.2e")
            kk = st.checkbox("Run finite-band causality (KK) diagnostic", value=True)
        go = st.form_submit_button("Run extraction", type="primary")
    if not go:
        if ss.get("current_key"):
            st.caption("Results for the last run are shown in the other tabs. Change settings and press "
                       "'Run extraction' to recompute.")
        return
    cfg = A.AnalysisConfig(
        thickness_mm=float(thickness), geometry="coax" if geom.startswith("Coax") else "waveguide",
        waveguide_a_mm=(a_custom if wg == "custom" else WAVEGUIDES[wg][0]) if geom.startswith("Rect") else None,
        f_min_ghz=None if win[0] <= fmin else win[0], f_max_ghz=None if win[1] >= fmax else win[1],
        method=method, assume_nonmagnetic=nonmag, branch_method=branch,
        prior_eps_mu=prior if branch == "prior" else None, offset1_mm=off1, offset2_mm=off2,
        holder_length_mm=holder if holder > 0 else None, noise_sigma=noise, average_ports=avg,
        regularization=reg, kk=kk)
    errs = cfg.validate()
    if errs:
        for e in errs:
            st.error(e)
        return
    if cfg.geometry == "waveguide":
        fc = 299_792_458.0 / (2 * cfg.waveguide_a_mm * 1e-3) / 1e9
        lo = win[0]
        if lo <= fc:
            st.error(f"TE10 cutoff of the empty guide is {fc:.4g} GHz; the band must start above it. Narrow the "
                     "frequency window.")
            return
        if win[1] > 2 * fc:
            st.warning(f"The band extends above 2·fc = {2*fc:.4g} GHz where the empty guide supports TE20; the "
                       "single-mode TE10 model is then invalid.")
    _run(cfg)


def _run(cfg, progress=None):
    ss = st.session_state
    key = cfg.key(A.data_hash(ss["data"]))
    out = state.recall(ss, key)
    if out is None:
        with st.spinner("Running extraction (core nrw package)…"):
            try:
                out = A.run_analysis(ss["data"], cfg, progress=progress, data_summary=ss["input_summary"])
            except ValueError as e:
                st.error(f"Extraction not possible: {e}")
                return None
        state.remember(ss, key, out)
    ss.pop("fig_memo", None)
    if not cfg.uncertainty:
        ss["current_key"] = key
        ss.pop("exports", None)
    st.success(f"Extraction finished in {out.runtime_s:.1f} s — see Results / Diagnostics.")
    return key


def _current_outcome():
    ss = st.session_state
    base = state.recall(ss, ss["current_key"]) if ss.get("current_key") else None
    mc = state.recall(ss, ss["mc_key"]) if ss.get("mc_key") else None
    if mc is not None and base is not None and _same_except_unc(mc.config, base.config):
        return mc
    return base


def _same_except_unc(a, b):
    da, db = asdict(a), asdict(b)
    da.pop("uncertainty")
    db.pop("uncertainty")
    return da == db


def _fig(outcome, name, builder, *args):
    """Session-local memo of figures per outcome (Streamlit re-executes the script on every
    interaction; figures of an unchanged outcome need not be rebuilt)."""
    ss = st.session_state
    memo = ss.setdefault("fig_memo", {})
    oid = id(outcome)
    if memo.get("_owner") != oid:
        memo.clear()
        memo["_owner"] = oid
    k = (name,) + tuple(str(a) for a in args)
    if k not in memo:
        memo[k] = builder(outcome, *args)
    return memo[k]


def _need(outcome) -> bool:
    if outcome is None:
        st.info("Run an extraction first ('Configure & extract').")
        return False
    return True


def _assessment_banner(outcome):
    a = outcome.assessment
    kind, text = LEVEL_UI[a.level]
    getattr(st, kind)(text)
    st.caption("The status is derived from documented diagnostic rules (src/nrw/assessment.py, docs/web_app.md). "
               "It is not an accuracy estimate: undetected systematic errors (thickness, air gaps, calibration) "
               "are possible even when GREEN.")
    if a.reasons:
        with st.expander(f"Diagnostic findings ({len(a.reasons)})", expanded=a.level != "GREEN"):
            for r in a.reasons:
                pts = f" — {r['n_points']} points ({100 * r['fraction']:.1f} %)" if r["n_points"] is not None else ""
                st.markdown(f"- **{r['severity']} {r['rule']}**: {r['message']}{pts}")


# ------------------------------------------------------------------------------------------
def _tab_results(outcome):
    if not _need(outcome):
        return
    res = outcome.result
    _assessment_banner(outcome)
    st.caption(f"Method: **{A.METHODS.get(outcome.config.method, ('',))[0]}** → `{res.method}` · branch: "
               f"`{res.branch_method}` · {len(res.frequency)} frequencies · shaded: red = unreliable, "
               "amber = attention. Curves are not smoothed and flagged points are not removed.")
    if res.method in A.EPS_ONLY:
        st.warning("μr = 1 was ASSUMED by this method; μ is not extracted.")
    keys = ["eps_real", "eps_loss"] + ([] if res.method in A.EPS_ONLY else ["mu_real", "mu_loss"])
    cols = st.columns(2)
    for i, k in enumerate(keys):
        cols[i % 2].plotly_chart(_fig(outcome, "param", P.parameter_figure, k), width="stretch")
    st.plotly_chart(_fig(outcome, "tan", P.loss_tangent_figure), width="stretch")
    st.markdown("##### Results table")
    st.dataframe(outcome.table, width="stretch", height=360, hide_index=True)


def _tab_diagnostics(outcome):
    if not _need(outcome):
        return
    res = outcome.result
    _assessment_banner(outcome)
    if res.warnings:
        st.markdown("##### Warnings")
        for w in res.warnings:
            st.warning(w)
    counts = {k: int(np.count_nonzero(v)) for k, v in res.flags.items() if np.shape(v) == res.frequency.shape}
    st.markdown("##### Flag counts (frequencies)")
    st.dataframe(pd.DataFrame([counts]), hide_index=True, width="stretch")
    c1, c2 = st.columns(2)
    c1.plotly_chart(_fig(outcome, "branch", P.branch_figure), width="stretch")
    c2.plotly_chart(_fig(outcome, "margin", P.margin_figure), width="stretch")
    c1, c2 = st.columns(2)
    c1.plotly_chart(_fig(outcome, "diag", P.diagnostic_figure, "gamma_T"), width="stretch")
    c2.plotly_chart(_fig(outcome, "diag", P.diagnostic_figure, "electrical_length"), width="stretch")
    c1, c2 = st.columns(2)
    c1.plotly_chart(_fig(outcome, "diag", P.diagnostic_figure, "identifiability"), width="stretch")
    c2.plotly_chart(_fig(outcome, "diag", P.diagnostic_figure, "solver"), width="stretch")
    c1, c2 = st.columns(2)
    c1.plotly_chart(_fig(outcome, "diag", P.diagnostic_figure, "impedance"), width="stretch")
    c2.plotly_chart(_fig(outcome, "diag", P.diagnostic_figure, "propagation"), width="stretch")
    st.markdown("##### Diagnostics table")
    st.dataframe(outcome.diagnostics, width="stretch", height=320, hide_index=True)


def _tab_passivity(outcome):
    if not _need(outcome):
        return
    res = outcome.result
    st.markdown("##### Passivity")
    st.caption("Measured result: the extracted values exactly as computed (never clipped). Physically constrained "
               "interpretation: 'pass' (loss ≥ 0), 'ambiguous' (negative but within 2σ of zero, i.e. consistent with "
               "a passive material given the uncertainty), 'fail' (below −2σ: inconsistent with a passive material; "
               "the data or the extraction are suspect there).")
    if res.passivity is not None:
        rows = [{"check": c.name, "status": c.status, "criterion": c.description, "violations": c.n_violations}
                for c in res.passivity.checks]
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        st.caption("status: exact = follows from passivity/linearity; assumption = holds only for reciprocal, "
                   "homogeneous, symmetric samples; heuristic = plausibility only.")
    t = outcome.table
    summary = {q: t[f"passivity_{q}"].value_counts().to_dict() for q in ("eps", "mu")}
    st.dataframe(pd.DataFrame(summary).fillna(0).astype(int, errors="ignore"), width="stretch")
    c1, c2 = st.columns(2)
    c1.plotly_chart(_fig(outcome, "pas", P.passivity_figure, "eps_loss"), width="stretch")
    if res.method not in A.EPS_ONLY:
        c2.plotly_chart(_fig(outcome, "pas", P.passivity_figure, "mu_loss"), width="stretch")
    st.markdown("##### Finite-band causality consistency")
    st.info(X.KK_DISCLAIMER)
    if not res.causality:
        st.warning("The KK diagnostic was not run or had too few well-conditioned points.")
    for q, r in res.causality.items():
        (st.success if r.passed else st.warning)(r.statement)
        with st.expander(f"Assumptions of the {q} diagnostic"):
            for s in r.assumptions:
                st.markdown(f"- {s}")
            st.caption(f"max |residual|/σ = {r.max_abs_normalized_residual:.3g}, rms = {r.rms_normalized_residual:.3g}; "
                       "ill-conditioned (flagged) frequencies are excluded from the test.")
        a, b = _fig(outcome, "kk", P.kk_figures, q)
        c1, c2 = st.columns(2)
        c1.plotly_chart(a, width="stretch")
        c2.plotly_chart(b, width="stretch")


def _tab_compare():
    ss = st.session_state
    base = state.recall(ss, ss["current_key"]) if ss.get("current_key") else None
    if not _need(base):
        return
    st.caption("Runs the core methods on the same data, frequency window, geometry, reference planes and branch "
               "method as the current extraction. Unstable regions are shown, not hidden. μ is not plotted for "
               "μ = 1 methods (it is assumed, not extracted).")
    choices = st.multiselect("Methods", ["nrw", "iterative_eps_mu", "baker_jarvis", "nrw_eps", "regularized_eps_mu"],
                             default=["nrw", "iterative_eps_mu", "baker_jarvis"],
                             format_func=lambda k: A.METHODS[k][0])
    key = ("cmp", ss["current_key"], tuple(choices))
    if st.button("Compare extraction methods", type="primary"):
        with st.spinner("Running methods…"):
            ss["comparison"] = (key, A.compare_methods(ss["data"], base.config, methods=choices))
    if not ss.get("comparison") or ss["comparison"][0] != key:
        return
    outs = ss["comparison"][1]
    for m, o in outs.items():
        if isinstance(o, Exception):
            st.error(f"{A.METHODS[m][0]}: {o}")
    labels = {k: A.METHODS[k][0] for k in A.METHODS}
    lv = {m: o.assessment.level for m, o in outs.items() if hasattr(o, "assessment")}
    st.dataframe(pd.DataFrame([{"method": labels[m], "status": lv.get(m, "error"),
                                "unreliable points": int(np.sum(o.assessment.point_status == "unreliable"))
                                if hasattr(o, "assessment") else None} for m, o in outs.items()]),
                 hide_index=True, width="stretch")
    cols = st.columns(2)
    for i, k in enumerate(["eps_real", "eps_loss", "mu_real", "mu_loss"]):
        cols[i % 2].plotly_chart(P.comparison_figure(outs, k, labels), width="stretch")
    dis = A.method_disagreement(outs)
    if not dis.empty:
        import plotly.graph_objects as go
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=dis["frequency_GHz"], y=dis["eps_rel_spread"], name="ε: max relative deviation"))
        fig.add_trace(go.Scatter(x=dis["frequency_GHz"], y=dis["mu_rel_spread"], name="μ: max relative deviation",
                                 line=dict(dash="dash")))
        fig.update_yaxes(type="log")
        st.plotly_chart(P._layout(fig, "Disagreement between simultaneous ε/μ methods", "relative deviation"),
                        width="stretch")


def _tab_uncertainty():
    ss = st.session_state
    base = state.recall(ss, ss["current_key"]) if ss.get("current_key") else None
    if not _need(base):
        return
    st.caption("Monte Carlo propagation through the complete pipeline (including branch selection) with the "
               "same configuration as the current extraction. Inputs are 1σ standard uncertainties.")
    with st.form("mc"):
        c = st.columns(3)
        s11m = c[0].number_input("|S11| relative (%)", min_value=0.0, value=0.5, step=0.1)
        s11p = c[1].number_input("S11 phase (deg)", min_value=0.0, value=0.5, step=0.1)
        noise = c[2].number_input("additive complex noise σ (linear)", min_value=0.0, value=0.0, format="%.2e")
        c = st.columns(3)
        s21m = c[0].number_input("|S21| relative (%)", min_value=0.0, value=0.5, step=0.1)
        s21p = c[1].number_input("S21 phase (deg)", min_value=0.0, value=0.5, step=0.1)
        seed = c[2].number_input("random seed", min_value=0, value=12345, step=1)
        c = st.columns(3)
        dL = c[0].number_input("thickness (µm)", min_value=0.0, value=5.0, step=1.0)
        dRP = c[1].number_input("reference plane, each port (µm)", min_value=0.0, value=5.0, step=1.0)
        nsel = c[2].selectbox("Monte Carlo draws", [100, 500, 1000, "custom"])
        ncustom = st.number_input("custom number of draws (20–5000)", min_value=20, max_value=5000, value=200)
        go = st.form_submit_button("Run uncertainty analysis", type="primary")
    n = int(ncustom if nsel == "custom" else nsel)
    per = max(base.runtime_s * 0.6, 0.05)
    st.caption(f"Estimated duration: about {n * per / 60:.1f} min (≈ {per:.2f} s per draw on this server).")
    if go:
        spec = {"s11_mag_rel": s11m / 100, "s11_phase_deg": s11p, "s21_mag_rel": s21m / 100, "s21_phase_deg": s21p,
                "s_noise": noise, "thickness_m": dL * 1e-6, "ref_plane_m": dRP * 1e-6, "n_draws": n, "seed": int(seed)}
        cfg = A.AnalysisConfig(**{**asdict(base.config), "uncertainty": spec})
        bar = st.progress(0.0, text="Monte Carlo…")
        t0 = time.time()

        def prog(done, total):
            el = time.time() - t0
            eta = el / done * (total - done)
            bar.progress(done / total, text=f"Monte Carlo draw {done}/{total} — about {eta:.0f} s remaining")

        key = _run(cfg, progress=prog)
        if key:
            ss["mc_key"] = key
            ss.pop("exports", None)
            bar.progress(1.0, text="Monte Carlo finished")
    out = _current_outcome()
    if out is not None and out.result.uncertainty:
        u = out.result.uncertainty
        st.success(f"Uncertainty available: {u['n_draws']} draws ({u['n_failed_draws']} failed). Bands are shown in "
                   "Results and in the plots below; the exported table contains percentiles.")
        st.caption("In the benchmark the 95 % Monte Carlo intervals covered the truth in 90–93 % of cases "
                   "(docs/benchmark_report.md §7): treat them as approximate.")
        keys = ["eps_real", "eps_loss"] + ([] if out.result.method in A.EPS_ONLY else ["mu_real", "mu_loss"])
        cols = st.columns(2)
        for i, k in enumerate(keys):
            cols[i % 2].plotly_chart(_fig(out, "param", P.parameter_figure, k), width="stretch", key=f"unc_{k}")


def _tab_export(outcome):
    ss = st.session_state
    if not _need(outcome):
        return
    c = st.columns(4)
    dpi = c[0].selectbox("Figure resolution (dpi)", [300, 600, 1200], index=0)
    show_unc = c[1].checkbox("Uncertainty bands in figures", value=True)
    show_flags = c[2].checkbox("Shade unreliable regions", value=True)
    include_input = c[3].checkbox("Include input file in ZIP", value=True)
    key = (id(outcome), dpi, show_unc, show_flags, include_input)
    if st.button("Prepare export files", type="primary") or (ss.get("exports") and ss["exports"][0] == key):
        if not ss.get("exports") or ss["exports"][0] != key:
            with st.spinner("Rendering figures and report…"):
                figs = P.publication_figures(outcome, dpi=dpi, show_uncertainty=show_unc, show_flags=show_flags)
                cmp = ss["comparison"][1] if ss.get("comparison") else None
                files = {
                    "extracted_parameters.csv": X.results_csv(outcome),
                    "extracted_parameters.json": X.results_json(outcome),
                    "diagnostics.csv": X.diagnostics_csv(outcome),
                    "analysis_report.html": X.report_html(outcome, figs),
                    "nrw_results.zip": X.zip_bundle(outcome, figs, ss.get("input_name") if include_input else None,
                                                    ss.get("input_bytes") if include_input else None, cmp),
                }
                ss["exports"] = (key, files, figs)
        _, files, figs = ss["exports"]
        mime = {"csv": "text/csv", "json": "application/json", "html": "text/html", "zip": "application/zip"}
        stem = Path(ss.get("input_name", "result")).stem
        st.markdown("##### Downloads")
        cols = st.columns(len(files))
        for col, (name, blob) in zip(cols, files.items()):
            fname = f"{stem}_{name}" if name != "nrw_results.zip" else f"{stem}_nrw_results.zip"
            col.download_button(name, blob, file_name=fname, mime=mime[name.rsplit('.', 1)[1]],
                                width="stretch")
        st.markdown("##### Publication figures (PNG)")
        cols = st.columns(3)
        for i, (name, png) in enumerate(figs.items()):
            with cols[i % 3]:
                st.image(png, caption=name, width="stretch")
                st.download_button(f"Download {name}", png, file_name=f"{stem}_{name}", mime="image/png",
                                   key=f"dl_{name}")
