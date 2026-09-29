"""Publication-quality figures (matplotlib, Agg backend).

Palette: validated categorical order (blue, orange, aqua, yellow, magenta). Ground truth
is always drawn in neutral black dashes; identity is never carried by colour alone
(legend + line style), and every figure has a CSV table next to it.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
import numpy as np  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
STYLES = ["-", "--", "-.", ":", (0, (5, 1, 1, 1))]
TRUTH = dict(color="#111111", ls=(0, (4, 2)), lw=1.2, label="ground truth")
INK_MUTED = "#6b6b66"
FLAG_FILL = "#e9e7e1"

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
    "legend.fontsize": 7.5, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#e6e5df", "grid.linewidth": 0.6, "lines.linewidth": 1.4, "axes.edgecolor": "#8a8a84",
    "xtick.color": "#4a4a45", "ytick.color": "#4a4a45", "legend.frameon": False,
})


def _shade(ax, f_ghz, mask, label=None):
    """Shade frequency spans where ``mask`` is True (e.g. ill-conditioned regions)."""
    mask = np.asarray(mask, bool)
    if not mask.any():
        return
    idx = np.flatnonzero(np.diff(np.concatenate([[0], mask.astype(int), [0]])))
    first = True
    for a, b in zip(idx[::2], idx[1::2]):
        lo = f_ghz[a]
        hi = f_ghz[min(b, len(f_ghz) - 1)] if b > a + 1 else f_ghz[a] + 1e-3
        ax.axvspan(lo, hi, color=FLAG_FILL, lw=0, zorder=0, label=label if first else None)
        first = False


def plot_result(res, outdir, truth=None, prefix=""):
    """Standard figure set for an ExtractionResult; returns list of written paths."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    f = res.frequency / 1e9
    paths = []
    unc = res.uncertainty
    sing = res.flags.get("singular_region", np.zeros(len(f), bool))
    quantities = [("eps_real", r"$\varepsilon'$", res.epsilon_real), ("eps_loss", r"$\varepsilon''$", res.epsilon_loss),
                  ("mu_real", r"$\mu'$", res.mu_real), ("mu_loss", r"$\mu''$", res.mu_loss)]
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.6), sharex=True, constrained_layout=True)
    for ax, (key, lab, y) in zip(axs.ravel(), quantities):
        _shade(ax, f, sing, "ill-conditioned (CRLB)")
        if unc:
            ax.fill_between(f, unc[key]["p2.5"], unc[key]["p97.5"], color=SERIES[0], alpha=0.18, lw=0,
                            label="95 % MC interval")
        ax.plot(f, y, color=SERIES[0], label=f"extracted ({res.method})")
        if truth is not None:
            tv = {"eps_real": truth[0].real, "eps_loss": -truth[0].imag, "mu_real": truth[1].real,
                  "mu_loss": -truth[1].imag}[key]
            ax.plot(f, tv, **TRUTH)
        ax.set_ylabel(lab)
        _robust_ylim(ax, y, unc[key] if unc else None, truth_vals=None)
    for ax in axs[1]:
        ax.set_xlabel("Frequency (GHz)")
    axs[0, 0].legend(loc="best")
    p = outdir / f"{prefix}parameters.png"
    fig.savefig(p)
    plt.close(fig)
    paths.append(p)

    # loss tangents
    fig, ax = plt.subplots(figsize=(3.5, 2.6), constrained_layout=True)
    _shade(ax, f, sing)
    ax.plot(f, res.tan_delta_eps, color=SERIES[0], label=r"tan$\delta_\varepsilon$")
    if res.method not in ("nrw_eps", "baker_jarvis"):
        ax.plot(f, res.tan_delta_mu, color=SERIES[1], ls="--", label=r"tan$\delta_\mu$")
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("loss tangent")
    _robust_ylim(ax, np.concatenate([res.tan_delta_eps, res.tan_delta_mu]))
    ax.legend()
    p = outdir / f"{prefix}loss_tangents.png"
    fig.savefig(p)
    plt.close(fig)
    paths.append(p)

    paths.append(plot_diagnostics(res, outdir / f"{prefix}diagnostics.png"))
    if res.causality:
        paths.append(plot_kk(res.causality, outdir / f"{prefix}kk_diagnostic.png"))
    return paths


def _robust_ylim(ax, y, u=None, truth_vals=None, pad=0.15):
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    if y.size == 0:
        return
    lo, hi = np.percentile(y, [2, 98])
    if u is not None:
        lo = min(lo, np.nanpercentile(u["p2.5"], 5))
        hi = max(hi, np.nanpercentile(u["p97.5"], 95))
    span = hi - lo if hi > lo else max(abs(hi), 1e-3)
    ax.set_ylim(lo - pad * span, hi + pad * span)


def plot_diagnostics(res, path):
    """Branch index, NRW singularity indicators, identifiability and passivity."""
    f = res.frequency / 1e9
    base = res.nrw_baseline
    fig, axs = plt.subplots(4, 1, figsize=(6.0, 6.6), sharex=True, constrained_layout=True)
    ax = axs[0]
    h = ax.step(f, res.selected_branch, where="mid", color=SERIES[0], label="selected branch m")
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_ylim(res.selected_branch.min() - 0.6, res.selected_branch.max() + 0.6)
    if res.branch_margin is not None:
        ax2 = ax.twinx()
        h += ax2.plot(f, res.branch_margin, color=SERIES[1], ls="--", lw=1, label="ambiguity margin (right)")
        ax2.set_ylim(-0.05, 1.3)
        ax2.set_ylabel("margin", color=INK_MUTED)
        ax2.spines["right"].set_visible(True)
        ax2.grid(False)
    ax.set_ylabel("branch m")
    ax.legend(h, [x.get_label() for x in h], loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0))
    ax.set_title(f"Phase branch ({res.branch_method})", loc="left")

    ax = axs[1]
    ax.semilogy(f, np.abs(base.X) if base is not None else np.nan, color=SERIES[0], label="|X|")
    ax.semilogy(f, 1 - np.abs(res.reflection_coefficient) + 1e-16, color=SERIES[1], ls="--", label="1 - |Γ|")
    ax.semilogy(f, np.abs(res.transmission), color=SERIES[2], ls="-.", label="|T|")
    ax.set_ylabel("NRW quantities")
    ax.legend(ncol=3, loc="best")
    el = base.electrical_length_over_pi if base is not None else None
    if el is not None:
        for k in range(1, int(np.nanmax(el)) + 1):
            idx = np.flatnonzero(np.diff(np.sign(el - k)))
            for i in idx:
                for a in axs:
                    a.axvline(f[i], color=INK_MUTED, lw=0.6, ls=":", zorder=0)
    ax.set_title("NRW singularity indicators (dotted: βL = kπ)", loc="left")

    ax = axs[2]
    idf = res.identifiability
    ax.semilogy(f, idf.get("eps_real_std", np.full(len(f), np.nan)), color=SERIES[0], label="CRLB σ(ε')")
    if "mu_real_std" in idf:
        ax.semilogy(f, idf["mu_real_std"], color=SERIES[1], ls="--", label="CRLB σ(μ')")
    if res.uncertainty:
        ax.semilogy(f, res.uncertainty["eps_real"]["sigma_percentile"], color=SERIES[2], ls="-.",
                    label="Monte Carlo σ(ε')")
    ax.set_ylabel("std. uncertainty")
    ax.legend(loc="best")
    ax.set_title(f"Identifiability at noise σ = {idf.get('noise_sigma_assumed', np.nan):g}", loc="left")

    ax = axs[3]
    ax.axhline(0, color="#111111", lw=0.8)
    ax.plot(f, res.epsilon_loss, color=SERIES[0], label="ε''")
    if res.method not in ("nrw_eps", "baker_jarvis"):
        ax.plot(f, res.mu_loss, color=SERIES[1], ls="--", label="μ''")
    absorbed = None
    if res.passivity is not None:
        try:
            absorbed = res.passivity["absorbed_power_port1"].metric
        except KeyError:
            absorbed = None
    if absorbed is not None:
        ax.plot(f, absorbed, color=SERIES[2], ls="-.", label="1-|S11|²-|S21|²")
    ax.set_ylabel("passivity (≥ 0)")
    ax.set_xlabel("Frequency (GHz)")
    ax.legend(loc="best", ncol=3)
    _robust_ylim(ax, np.concatenate([res.epsilon_loss, res.mu_loss, absorbed if absorbed is not None else []]))
    ax.set_title("Passivity: ε'', μ'' and absorbed power must be ≥ 0 (exp(+jωt))", loc="left")
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_kk(caus: dict, path):
    fig, axs = plt.subplots(len(caus), 2, figsize=(7.0, 2.4 * len(caus)), constrained_layout=True, squeeze=False)
    for row, (q, r) in enumerate(caus.items()):
        f = r.frequency / 1e9
        ax = axs[row, 0]
        ax.plot(f, r.real_measured, color=SERIES[0], label=f"Re({q}) extracted")
        ax.plot(f, r.real_predicted, color=SERIES[1], ls="--", label="finite-band KK + passive tail")
        ax.set_ylabel(f"Re({q})")
        ax.legend(loc="best")
        ax.set_title(f"{q}: {'consistent' if r.passed else 'NOT consistent'} (finite-band diagnostic)", loc="left")
        ax = axs[row, 1]
        ax.plot(f, r.normalized_residual, color=SERIES[0])
        thr = float(r.assumptions[-1].split("<=")[-1]) if r.assumptions else 3.0
        for s in (-thr, thr):
            ax.axhline(s, color=INK_MUTED, ls=":", lw=0.8)
        ax.set_ylabel("residual / σ")
        for a in axs[row]:
            a.set_xlabel("Frequency (GHz)")
    fig.savefig(path)
    plt.close(fig)
    return path
