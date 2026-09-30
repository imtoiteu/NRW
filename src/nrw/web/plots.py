"""Figures for the web application.

* Interactive figures: Plotly (hover, zoom). Problematic frequency regions are shaded
  (red = unreliable, amber = attention) and never removed from the curves.
* Publication figures: Matplotlib via nrw.plotting's style, PNG at a chosen dpi (>= 300).

Colours: validated categorical palette of nrw.plotting (identity also carried by line style /
marker symbol and legends, never colour alone).
"""
from __future__ import annotations

import io

import numpy as np
import plotly.graph_objects as go

from ..plotting import SERIES, STYLES, apply_style, plt

TRUTHLESS_NOTE = "measured-data extraction (no ground truth)"
SHADE = {"unreliable": "rgba(214, 39, 40, 0.13)", "attention": "rgba(237, 161, 0, 0.13)"}
DASH = ["solid", "dash", "dashdot", "dot", "longdash"]
SYMBOL = {"pass": "circle", "ambiguous": "diamond-open", "fail": "x", "undetermined": "square-open"}
PASS_COLOR = {"pass": SERIES[2], "ambiguous": SERIES[3], "fail": "#d62728", "undetermined": "#6b6b66"}

LABELS = {
    "eps_real": "ε′ (relative permittivity, real part)",
    "eps_loss": "ε″ = −Im(εr)  (loss, ≥ 0 if passive)",
    "mu_real": "μ′ (relative permeability, real part)",
    "mu_loss": "μ″ = −Im(μr)  (loss, ≥ 0 if passive)",
}
SHORT = {"eps_real": "ε′", "eps_loss": "ε″", "mu_real": "μ′", "mu_loss": "μ″"}


def _layout(fig, title, ytitle, xtitle="Frequency (GHz)", height=330):
    fig.update_layout(title=dict(text=title, x=0, font=dict(size=14)), height=height,
                      margin=dict(l=60, r=20, t=45, b=45), template="plotly_white",
                      legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, font=dict(size=11)),
                      hovermode="x unified")
    fig.update_xaxes(title_text=xtitle, showline=True, linecolor="#8a8a84", gridcolor="#eceae4")
    fig.update_yaxes(title_text=ytitle, showline=True, linecolor="#8a8a84", gridcolor="#eceae4")
    return fig


def _runs(mask):
    mask = np.asarray(mask, bool)
    if not mask.any():
        return []
    idx = np.flatnonzero(np.diff(np.concatenate([[0], mask.astype(int), [0]])))
    return list(zip(idx[::2], idx[1::2]))


def _merge_runs(runs, gap=2):
    """Merge runs separated by <= ``gap`` points (display only; exact flags stay in the table)."""
    out = []
    for a, b in runs:
        if out and a - out[-1][1] <= gap:
            out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def shade_status(fig, f_ghz, status):
    """Shade contiguous runs of 'unreliable' / 'attention' frequencies.

    All rectangles are added in ONE layout update (Plotly's add_vrect re-validates the whole
    layout per call, which is quadratic in the number of regions). 'unreliable' is drawn on
    top of 'attention'. Runs separated by <= 2 points are merged for display only.
    """
    if status is None or len(f_ghz) < 2:
        return fig
    status = np.asarray(status)
    df = float(np.median(np.diff(f_ghz)))
    shapes = list(fig.layout.shapes or [])
    for key in ("attention", "unreliable"):
        for a, b in _merge_runs(_runs(status == key)):
            shapes.append(dict(type="rect", xref="x", yref="paper", x0=float(f_ghz[a] - df / 2),
                               x1=float(f_ghz[b - 1] + df / 2), y0=0, y1=1, fillcolor=SHADE[key],
                               line=dict(width=0), layer="below"))
    fig.update_layout(shapes=shapes)
    return fig


def _legend_for_shading(fig, status):
    if status is None:
        return
    for key, name in (("unreliable", "unreliable region"), ("attention", "attention region")):
        if np.any(np.asarray(status) == key):
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=name,
                                     marker=dict(size=12, symbol="square", color=SHADE[key].replace("0.13", "0.5"))))


# ------------------------------------------------------------------------------------------
# Raw data
# ------------------------------------------------------------------------------------------
def sparam_magnitude(data, f_window=None):
    f = data.frequency / 1e9
    fig = go.Figure()
    for i, k in enumerate(("s11", "s21", "s12", "s22")):
        if k not in data.measured:
            continue
        with np.errstate(divide="ignore"):
            y = 20 * np.log10(np.abs(getattr(data, k)))
        fig.add_trace(go.Scatter(x=f, y=y, name=f"|{k.upper()}|", line=dict(color=SERIES[i], dash=DASH[i], width=1.6)))
    if f_window:
        _window_lines(fig, f_window)
    return _layout(fig, "Magnitude", "|S| (dB)")


def sparam_phase(data, unwrap=False, f_window=None):
    f = data.frequency / 1e9
    fig = go.Figure()
    for i, k in enumerate(("s11", "s21", "s12", "s22")):
        if k not in data.measured:
            continue
        ph = np.angle(getattr(data, k))
        if unwrap:
            ph = np.unwrap(ph)
        fig.add_trace(go.Scatter(x=f, y=np.degrees(ph), name=f"∠{k.upper()}",
                                 line=dict(color=SERIES[i], dash=DASH[i], width=1.4)))
    if f_window:
        _window_lines(fig, f_window)
    title = "Phase (unwrapped with numpy.unwrap - display only, not a branch decision)" if unwrap else "Phase (wrapped)"
    return _layout(fig, title, "phase (deg)")


def _window_lines(fig, f_window):
    lo, hi = f_window
    for x in (lo, hi):
        if x is not None:
            fig.add_vline(x=x, line=dict(color="#6b6b66", dash="dot", width=1))


# ------------------------------------------------------------------------------------------
# Results
# ------------------------------------------------------------------------------------------
def _values(res, key):
    return {"eps_real": res.epsilon_real, "eps_loss": res.epsilon_loss, "mu_real": res.mu_real,
            "mu_loss": res.mu_loss}[key]


def parameter_figure(outcome, key, show_status=True):
    res = outcome.result
    f = res.frequency / 1e9
    fig = go.Figure()
    status = outcome.assessment.point_status if show_status else None
    shade_status(fig, f, status)
    if res.uncertainty:
        u = res.uncertainty[key]
        fig.add_trace(go.Scatter(x=np.concatenate([f, f[::-1]]), y=np.concatenate([u["p97.5"], u["p2.5"][::-1]]),
                                 fill="toself", fillcolor="rgba(42,120,214,0.18)", line=dict(width=0),
                                 name="95 % Monte Carlo interval", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=f, y=_values(res, key), name=f"{SHORT[key]} ({res.method})",
                             line=dict(color=SERIES[0], width=1.8)))
    if key.endswith("loss"):
        fig.add_hline(y=0, line=dict(color="#111111", width=1))
    _legend_for_shading(fig, status)
    _robust_range(fig, _values(res, key), res.uncertainty[key] if res.uncertainty else None)
    return _layout(fig, LABELS[key], SHORT[key] + "  (–)")


def _robust_range(fig, y, u=None):
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    if y.size < 3:
        return
    lo, hi = np.percentile(y, [1, 99])
    if u is not None:
        lo = min(lo, np.nanpercentile(u["p2.5"], 2))
        hi = max(hi, np.nanpercentile(u["p97.5"], 98))
    span = max(hi - lo, 1e-3 * max(abs(hi), 1))
    ymin, ymax = np.min(y), np.max(y)
    # keep the full data range unless a few extreme spikes would flatten everything else;
    # clipped spikes are still in the data table and are marked by the shading
    if (ymax - ymin) > 6 * span:
        fig.update_yaxes(range=[lo - 0.3 * span, hi + 0.3 * span])
        fig.add_annotation(text="axis range limited to the 1-99 % quantiles; spikes extend beyond (see table)",
                           xref="paper", yref="paper", x=0, y=1.0, showarrow=False, font=dict(size=10, color="#6b6b66"),
                           xanchor="left", yanchor="bottom")


def loss_tangent_figure(outcome):
    res = outcome.result
    f = res.frequency / 1e9
    fig = go.Figure()
    shade_status(fig, f, outcome.assessment.point_status)
    fig.add_trace(go.Scatter(x=f, y=res.tan_delta_eps, name="tan δε = ε″/ε′", line=dict(color=SERIES[0])))
    if res.method not in ("baker_jarvis", "nrw_eps"):
        fig.add_trace(go.Scatter(x=f, y=res.tan_delta_mu, name="tan δμ = μ″/μ′",
                                 line=dict(color=SERIES[1], dash="dash")))
    _robust_range(fig, np.concatenate([res.tan_delta_eps, res.tan_delta_mu]))
    return _layout(fig, "Loss tangents", "tan δ  (–)")


def branch_figure(outcome):
    res = outcome.result
    f = res.frequency / 1e9
    m = res.selected_branch
    fig = go.Figure()
    shade_status(fig, f, outcome.assessment.point_status)
    fig.add_trace(go.Scatter(x=f, y=m, mode="lines", line=dict(shape="hvh", color=SERIES[0], width=2),
                             name=f"selected branch m ({res.branch_method})"))
    ch = np.flatnonzero(np.diff(m) != 0) + 1
    if ch.size:
        fig.add_trace(go.Scatter(x=f[ch], y=m[ch], mode="markers", name="branch change",
                                 marker=dict(symbol="diamond", size=9, color=SERIES[1])))
    dis = res.flags.get("branch_disagreement")
    if dis is not None and np.any(dis):
        fig.add_trace(go.Scatter(x=f[dis], y=m[dis], mode="markers", name="trackers disagree",
                                 marker=dict(symbol="x", size=9, color="#d62728")))
    fig.update_yaxes(dtick=1)
    return _layout(fig, "Phase-branch index m (βL = −arg T + 2πm)", "m")


def margin_figure(outcome):
    res = outcome.result
    f = res.frequency / 1e9
    fig = go.Figure()
    if res.branch_margin is not None:
        fig.add_trace(go.Scatter(x=f, y=res.branch_margin, name="ambiguity margin", line=dict(color=SERIES[0])))
        fig.add_hline(y=0.5, line=dict(color="#6b6b66", dash="dot"), annotation_text="weak below 0.5")
    if res.uncertainty and "branch_change_fraction" in res.uncertainty:
        fig.add_trace(go.Scatter(x=f, y=res.uncertainty["branch_change_fraction"], name="MC branch-change fraction",
                                 line=dict(color=SERIES[1], dash="dash")))
    fig.update_yaxes(range=[-0.05, 1.05])
    return _layout(fig, "Branch confidence", "fraction / margin", height=260)


def diagnostic_figure(outcome, which):
    """which: 'gamma_T', 'electrical_length', 'identifiability', 'solver', 'impedance', 'propagation'."""
    res = outcome.result
    base = res.nrw_baseline
    f = res.frequency / 1e9
    fig = go.Figure()
    shade_status(fig, f, outcome.assessment.point_status)
    if which == "gamma_T":
        fig.add_trace(go.Scatter(x=f, y=np.abs(res.reflection_coefficient), name="|Γ| (interface reflection)",
                                 line=dict(color=SERIES[0])))
        fig.add_trace(go.Scatter(x=f, y=np.abs(res.transmission), name="|T| (sample transmission)",
                                 line=dict(color=SERIES[1], dash="dash")))
        fig.add_hline(y=1, line=dict(color="#111111", width=1))
        return _layout(fig, "Reflection coefficient Γ and transmission T (NRW)", "magnitude")
    if which == "electrical_length" and base is not None:
        el = base.electrical_length_over_pi
        fig.add_trace(go.Scatter(x=f, y=el, name="βL/π", line=dict(color=SERIES[0])))
        for k in range(1, int(np.nanmax(el)) + 1 if np.isfinite(np.nanmax(el)) else 1):
            fig.add_hline(y=k, line=dict(color="#6b6b66", dash="dot", width=1))
        return _layout(fig, "Electrical length βL/π (integers = half-wavelength points)", "βL/π")
    if which == "identifiability":
        idf = res.identifiability
        for i, (k, lab) in enumerate((("eps_real_std", "σ(ε′) CRLB"), ("eps_loss_std", "σ(ε″) CRLB"),
                                      ("mu_real_std", "σ(μ′) CRLB"), ("mu_loss_std", "σ(μ″) CRLB"))):
            if k in idf:
                fig.add_trace(go.Scatter(x=f, y=idf[k], name=lab, line=dict(color=SERIES[i], dash=DASH[i])))
        fig.update_yaxes(type="log")
        return _layout(fig, f"Cramér–Rao bound at noise σ = {idf.get('noise_sigma_assumed', float('nan')):g}",
                       "standard deviation")
    if which == "solver":
        conv = res.convergence
        if "residual_rms" in conv:
            fig.add_trace(go.Scatter(x=f, y=conv["residual_rms"], name="fit residual rms |ΔS|",
                                     line=dict(color=SERIES[0])))
        if "jacobian_cond" in conv:
            fig.add_trace(go.Scatter(x=f, y=conv["jacobian_cond"], name="Jacobian condition number",
                                     line=dict(color=SERIES[1], dash="dash")))
        idf = res.identifiability
        if "jacobian_cond" in idf:
            fig.add_trace(go.Scatter(x=f, y=idf["jacobian_cond"], name="model Jacobian condition (CRLB)",
                                     line=dict(color=SERIES[2], dash="dot")))
        fig.update_yaxes(type="log")
        return _layout(fig, "Solver residual and conditioning", "value")
    if which == "impedance":
        fig.add_trace(go.Scatter(x=f, y=res.impedance.real, name="Re z", line=dict(color=SERIES[0])))
        fig.add_trace(go.Scatter(x=f, y=res.impedance.imag, name="Im z", line=dict(color=SERIES[1], dash="dash")))
        _robust_range(fig, np.concatenate([res.impedance.real, res.impedance.imag]))
        return _layout(fig, "Normalised wave impedance z = μγ0/γ", "z (–)")
    if which == "propagation":
        fig.add_trace(go.Scatter(x=f, y=res.propagation_constant.real, name="α (Np/m)", line=dict(color=SERIES[0])))
        fig.add_trace(go.Scatter(x=f, y=res.propagation_constant.imag, name="β (rad/m)",
                                 line=dict(color=SERIES[1], dash="dash")))
        return _layout(fig, "Propagation constant γ = α + jβ", "1/m")
    raise ValueError(which)


def passivity_figure(outcome, key):
    """Loss term with its uncertainty band; markers classify pass / ambiguous / fail."""
    from ..assessment import passivity_classes

    res = outcome.result
    f = res.frequency / 1e9
    which = key.split("_")[0]
    y = _values(res, key)
    sig = (res.uncertainty[key]["sigma_percentile"] if res.uncertainty
           else res.identifiability.get(f"{which}_loss_std"))
    cls = passivity_classes(y, sig)
    fig = go.Figure()
    if sig is not None:
        fig.add_trace(go.Scatter(x=np.concatenate([f, f[::-1]]), y=np.concatenate([y + 2 * sig, (y - 2 * sig)[::-1]]),
                                 fill="toself", fillcolor="rgba(42,120,214,0.15)", line=dict(width=0),
                                 name="±2σ " + ("(Monte Carlo)" if res.uncertainty else "(Cramér–Rao)"),
                                 hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=f, y=y, name=f"measured {SHORT[key]} (unaltered)", line=dict(color=SERIES[0], width=1.4)))
    for c in ("ambiguous", "fail", "undetermined"):
        m = cls == c
        if m.any():
            fig.add_trace(go.Scatter(x=f[m], y=np.where(np.isfinite(y[m]), y[m], 0), mode="markers", name=c,
                                     marker=dict(symbol=SYMBOL[c], size=8, color=PASS_COLOR[c])))
    fig.add_hline(y=0, line=dict(color="#111111", width=1.2))
    _robust_range(fig, y)
    return _layout(fig, f"Passivity of {SHORT[key]} (exp(+jωt): passive ⇒ {SHORT[key]} ≥ 0)", SHORT[key])


def kk_figures(outcome, q):
    r = outcome.result.causality.get(q)
    if r is None:
        return None, None
    f = r.frequency / 1e9
    a = go.Figure()
    a.add_trace(go.Scatter(x=f, y=r.real_measured, name=f"Re({q}) extracted", line=dict(color=SERIES[0])))
    a.add_trace(go.Scatter(x=f, y=r.real_predicted, name="finite-band KK + best passive out-of-band tail",
                           line=dict(color=SERIES[1], dash="dash")))
    _layout(a, f"{q}: finite-band causality consistency", f"Re({q})")
    b = go.Figure()
    b.add_trace(go.Scatter(x=f, y=r.normalized_residual, name="normalised residual", line=dict(color=SERIES[0])))
    for s in (-r.threshold, r.threshold):
        b.add_hline(y=s, line=dict(color="#d62728", dash="dot"))
    _layout(b, f"{q}: residual / σ (threshold ±{r.threshold:.3g})", "residual / σ", height=260)
    return a, b


def comparison_figure(outcomes: dict, key, method_labels: dict):
    fig = go.Figure()
    i = 0
    for name, o in outcomes.items():
        if not hasattr(o, "result"):
            continue
        res = o.result
        if key.startswith("mu") and res.method in ("baker_jarvis", "nrw_eps"):
            continue  # mu is assumed, not extracted
        fig.add_trace(go.Scatter(x=res.frequency / 1e9, y=_values(res, key), name=method_labels.get(name, name),
                                 line=dict(color=SERIES[i % len(SERIES)], dash=DASH[i % len(DASH)], width=1.5)))
        i += 1
    ys = [(_values(o.result, key)) for o in outcomes.values() if hasattr(o, "result")
          and not (key.startswith("mu") and o.result.method in ("baker_jarvis", "nrw_eps"))]
    if ys:
        _robust_range(fig, np.concatenate(ys))
    if key.endswith("loss"):
        fig.add_hline(y=0, line=dict(color="#111111", width=1))
    return _layout(fig, LABELS[key], SHORT[key])


# ------------------------------------------------------------------------------------------
# Publication figures (Matplotlib)
# ------------------------------------------------------------------------------------------
PUB_LABELS = {"eps_real": r"$\varepsilon'$", "eps_loss": r"$\varepsilon''$", "mu_real": r"$\mu'$",
              "mu_loss": r"$\mu''$"}


def _png(fig, dpi):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def _pub_axis(ax, res, key, status, show_unc, show_flags):
    f = res.frequency / 1e9
    y = _values(res, key)
    if show_flags and status is not None:
        for a, b in _runs(np.asarray(status) == "unreliable"):
            ax.axvspan(f[a], f[b - 1], color="#f3d6d6", lw=0, zorder=0)
    if show_unc and res.uncertainty:
        u = res.uncertainty[key]
        ax.fill_between(f, u["p2.5"], u["p97.5"], color=SERIES[0], alpha=0.2, lw=0, label="95 % interval (MC)")
    ax.plot(f, y, color=SERIES[0], lw=1.3, label="extracted")
    if key.endswith("loss"):
        ax.axhline(0, color="#111111", lw=0.7)
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel(PUB_LABELS[key])
    yy = y[np.isfinite(y)]
    if yy.size > 3:
        lo, hi = np.percentile(yy, [1, 99])
        span = max(hi - lo, 1e-3 * max(abs(hi), 1))
        if np.ptp(yy) > 6 * span:
            ax.set_ylim(lo - 0.3 * span, hi + 0.3 * span)


def publication_figures(outcome, dpi=300, show_uncertainty=True, show_flags=True) -> dict:
    """Return {filename: PNG bytes}: one figure per quantity, a 2x2 panel, loss tangents, branch."""
    apply_style()
    res = outcome.result
    status = outcome.assessment.point_status
    eps_only = res.method in ("baker_jarvis", "nrw_eps")
    keys = ["eps_real", "eps_loss"] + ([] if eps_only else ["mu_real", "mu_loss"])
    out = {}
    for k in keys:
        fig, ax = plt.subplots(figsize=(3.4, 2.6), constrained_layout=True)
        _pub_axis(ax, res, k, status, show_uncertainty, show_flags)
        if show_uncertainty and res.uncertainty:
            ax.legend(fontsize=7)
        out[f"{k}.png"] = _png(fig, dpi)
    n = len(keys)
    fig, axs = plt.subplots(n // 2, 2, figsize=(7.0, 2.6 * (n // 2)), constrained_layout=True, squeeze=False)
    for ax, k, tag in zip(axs.ravel(), keys, "abcd"):
        _pub_axis(ax, res, k, status, show_uncertainty, show_flags)
        ax.text(0.02, 0.95, f"({tag})", transform=ax.transAxes, va="top", fontsize=9)
    out["parameters_panel.png"] = _png(fig, dpi)
    fig, ax = plt.subplots(figsize=(3.4, 2.6), constrained_layout=True)
    f = res.frequency / 1e9
    ax.plot(f, res.tan_delta_eps, color=SERIES[0], label=r"tan$\,\delta_\varepsilon$")
    if not eps_only:
        ax.plot(f, res.tan_delta_mu, color=SERIES[1], ls=STYLES[1], label=r"tan$\,\delta_\mu$")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("loss tangent")
    ax.legend(fontsize=7)
    out["loss_tangents.png"] = _png(fig, dpi)
    fig, ax = plt.subplots(figsize=(3.4, 2.2), constrained_layout=True)
    ax.step(f, res.selected_branch, where="mid", color=SERIES[0])
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("branch index m")
    ax.yaxis.get_major_locator().set_params(integer=True)
    out["branch_index.png"] = _png(fig, dpi)
    return out
