"""Phase-branch selection study: when does "continuous unwrap from low frequency" work?

Strategies: principal (m=0), unwrap (from f_min with m=0), group_delay (Weir, pointwise),
anchored (segment-wise unwrap + band-level group-delay anchoring), dp (dynamic
programming on eps*mu continuity + group delay).

Conditions varied one at a time from a baseline (eps = 10 - 0.1j, mu = 1, L = 10 mm,
2-18 GHz, 401 points, noise 1e-3): start frequency, grid density, noise, thickness,
in-band dispersion, transmission near the noise floor, missing band.

Outputs results/studies/branch_*.csv, results/figures/branch_*.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchlib import GEOM, ROOT  # noqa: E402

from nrw.branch import select_branch  # noqa: E402
from nrw.classical import nrw  # noqa: E402
from nrw.forward import SampleSpec  # noqa: E402
from nrw.materials import Debye, Lorentz, as_model  # noqa: E402
from nrw.perturb import complex_gaussian  # noqa: E402
from nrw.plotting import SERIES, STYLES, plt  # noqa: E402

OUT = ROOT / "results" / "studies"
FIG = ROOT / "results" / "figures"
STRATS = ["principal", "unwrap", "group_delay", "anchored", "dp"]

BASE = dict(eps=10 - 0.1j, mu=1.0, L=10e-3, fmin=2e9, fmax=18e9, N=401, sigma=1e-3, gap=None)
CONDITIONS = [
    ("baseline", {}),
    ("start_8GHz", {"fmin": 8e9}),
    ("start_8GHz_thick", {"fmin": 8e9, "L": 30e-3}),
    ("grid_101", {"N": 101}),
    ("grid_33", {"N": 33}),
    ("grid_17", {"N": 17}),
    ("grid_17_thick", {"N": 17, "L": 30e-3}),
    ("grid_17_60mm_aliased", {"N": 17, "L": 60e-3}),
    ("noise_1e-2", {"sigma": 1e-2}),
    ("noise_3e-2", {"sigma": 3e-2}),
    ("thick_30mm", {"L": 30e-3}),
    ("thick_60mm", {"L": 60e-3}),
    ("dispersive_lorentz_inband", {"eps": Lorentz(4.0, 6.0, 10e9, 1.0e9), "L": 10e-3}),
    ("dispersive_debye_below_band", {"eps": Debye(4.0, 40.0, 0.5e9), "L": 10e-3}),
    ("lossy_S21_near_floor", {"eps": Debye(10.0, 30.0, 3e9), "L": 10e-3, "sigma": 3e-3}),
    ("missing_band_8_12GHz", {"gap": (8e9, 12e9)}),
    ("missing_band_thick", {"gap": (8e9, 12e9), "L": 30e-3}),
]


def dataset(cfg, seed):
    f = np.linspace(cfg["fmin"], cfg["fmax"], cfg["N"])
    if cfg["gap"] is not None:
        f = f[(f < cfg["gap"][0]) | (f > cfg["gap"][1])]
    spec = SampleSpec(eps=cfg["eps"], mu=cfg["mu"], length_m=cfg["L"], geometry=GEOM)
    s = spec.sparams(f)
    rng = np.random.default_rng(seed)
    if cfg["sigma"]:
        s = [complex_gaussian(x, cfg["sigma"], rng) for x in s]
    eps_t = as_model(cfg["eps"])(f)
    mu_t = as_model(cfg["mu"])(f)
    beta = GEOM.gamma(f, eps_t, mu_t).imag
    mtrue = np.floor((beta * cfg["L"] + np.pi) / (2 * np.pi)).astype(int)
    return f, s, mtrue, eps_t


def main(n_seeds=10):
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, over in CONDITIONS:
        cfg = dict(BASE, **over)
        for seed in range(n_seeds):
            f, s, mtrue, eps_t = dataset(cfg, seed)
            a11, a21 = (s[0] + s[3]) / 2, (s[1] + s[2]) / 2
            T = nrw(f, a11, a21, cfg["L"], GEOM).transmission
            psi_u = np.unwrap(-np.angle(T))
            max_step = float(np.max(np.abs(np.diff(-GEOM.gamma(f, eps_t, cfg['mu']).imag * cfg['L']))))
            for st in STRATS:
                with np.errstate(all="ignore"):
                    b = select_branch(st, f, T, cfg["L"], GEOM)
                wrong = b.m != mtrue
                r = nrw(f, a11, a21, cfg["L"], GEOM, branch=b.m)
                rel = np.abs(r.eps_mu_only - eps_t) / np.abs(eps_t)
                rows.append({"condition": name, "strategy": st, "seed": seed, "n_freq": len(f),
                             "m_true_min": int(mtrue.min()), "m_true_max": int(mtrue.max()),
                             "true_max_phase_step_rad": max_step,
                             "branch_error_rate": float(np.mean(wrong)), "any_branch_error": bool(wrong.any()),
                             "median_rel_err_eps_mu1": float(np.nanmedian(rel)),
                             "discontinuities": b.n_discontinuities})
        print(name, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "branch_raw.csv", index=False)
    summ = df.pivot_table(index="condition", columns="strategy", values="any_branch_error", aggfunc="mean")
    summ = summ.reindex([c for c, _ in CONDITIONS])[STRATS]
    summ.to_csv(OUT / "branch_summary_any_error_fraction.csv")
    summ2 = df.pivot_table(index="condition", columns="strategy", values="branch_error_rate", aggfunc="mean")
    summ2.reindex([c for c, _ in CONDITIONS])[STRATS].to_csv(OUT / "branch_summary_error_rate.csv")
    info = df.groupby("condition")[["n_freq", "m_true_min", "m_true_max", "true_max_phase_step_rad"]].first()
    info.reindex([c for c, _ in CONDITIONS]).to_csv(OUT / "branch_conditions.csv")
    print(summ.to_string())
    branch_figures()


def branch_figures():
    """Branch-diagnostic plots for two instructive failures of plain unwrapping."""
    for name in ("start_8GHz_thick", "grid_17_thick", "dispersive_lorentz_inband"):
        cfg = dict(BASE, **dict(CONDITIONS)[name])
        f, s, mtrue, eps_t = dataset(cfg, 0)
        a11, a21 = (s[0] + s[3]) / 2, (s[1] + s[2]) / 2
        T = nrw(f, a11, a21, cfg["L"], GEOM).transmission
        fig, axs = plt.subplots(2, 1, figsize=(6.0, 4.4), sharex=True, constrained_layout=True)
        axs[0].step(f / 1e9, mtrue, where="mid", color="#111111", ls=(0, (4, 2)), lw=1.2, label="true m")
        for i, st in enumerate(["unwrap", "group_delay", "anchored", "dp"]):
            b = select_branch(st, f, T, cfg["L"], GEOM)
            off = (i - 1.5) * 0.06
            axs[0].step(f / 1e9, b.m + off, where="mid", color=SERIES[i], ls=STYLES[i], label=st)
            r = nrw(f, a11, a21, cfg["L"], GEOM, branch=b.m)
            axs[1].plot(f / 1e9, r.eps_mu_only.real, color=SERIES[i], ls=STYLES[i], label=st,
                        marker="o" if len(f) < 40 else None, ms=3)
        axs[1].plot(f / 1e9, eps_t.real, color="#111111", ls=(0, (4, 2)), lw=1.2, label="truth")
        axs[0].set_ylabel("branch m")
        axs[0].legend(ncol=5, fontsize=6.5, loc="upper left")
        axs[0].set_title(f"Branch selection: {name}", loc="left")
        axs[1].set_ylabel("ε' (μ = 1)")
        axs[1].set_xlabel("Frequency (GHz)")
        lo, hi = np.percentile(eps_t.real, [0, 100])
        axs[1].set_ylim(lo - 0.6 * (hi - lo + 1), hi + 0.6 * (hi - lo + 1))
        fig.savefig(FIG / f"branch_{name}.png")
        plt.close(fig)


if __name__ == "__main__":
    main()
