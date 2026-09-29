"""Half-wavelength (NRW-singularity) study: algorithmic vs intrinsic instability.

For low-loss samples crossing beta L = k pi, compares
  * classical NRW (simultaneous eps, mu),
  * iterative simultaneous (eps, mu) least squares (all four S-parameters),
  * NRW-derived eps with mu = 1 (non-iterative),
  * Baker-Jarvis-type iterative eps-only (mu = 1),
  * two-thickness simultaneous (eps, mu),
against the Cramer-Rao lower bound (CRLB) for the same noise.

Outputs results/studies/singularity_*.csv and results/figures/singularity_*.png
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
from nrw.conditioning import crlb  # noqa: E402
from nrw.forward import SampleSpec  # noqa: E402
from nrw.iterative import solve_eps_mu_multithickness, solve_eps_mu_pointwise, solve_eps_only  # noqa: E402
from nrw.perturb import complex_gaussian  # noqa: E402
from nrw.plotting import SERIES, STYLES, plt  # noqa: E402

OUT = ROOT / "results" / "studies"
FIG = ROOT / "results" / "figures"


def run(eps_t, mu_t, L, L2, sigma=1e-3, n_seeds=20, N=801, label="dielectric"):
    f = np.linspace(2e9, 18e9, N)
    clean = SampleSpec(eps=eps_t, mu=mu_t, length_m=L, geometry=GEOM).sparams(f)
    clean2 = SampleSpec(eps=eps_t, mu=mu_t, length_m=L2, geometry=GEOM).sparams(f)
    E = {k: [] for k in ("nrw", "iter_eps_mu", "nrw_eps_mu1", "bj_eps", "two_thickness")}
    M = {k: [] for k in E}
    for seed in range(n_seeds):
        rng = np.random.default_rng(1000 + seed)
        s = [complex_gaussian(x, sigma, rng) for x in clean]
        s2 = [complex_gaussian(x, sigma, rng) for x in clean2]
        a11, a21 = (s[0] + s[3]) / 2, (s[1] + s[2]) / 2
        base = nrw(f, a11, a21, L, GEOM)
        b = select_branch("anchored", f, base.transmission, L, GEOM)
        r = nrw(f, a11, a21, L, GEOM, branch=b.m)
        it = solve_eps_mu_pointwise(f, s[0], s[1], L, GEOM, r.eps, r.mu, s12=s[2], s22=s[3], multistart=False)
        bj = solve_eps_only(f, s[0], s[1], L, GEOM, r.eps_mu_only, s12=s[2], s22=s[3], formulation="eq2.24+2.25")
        # two-thickness initialised from the (independent) second-sample NRW solution where
        # the first is ill-conditioned: use the eps-only estimate as a neutral starting point
        e0 = np.where(np.abs(a11) > 5 * sigma, r.eps, bj.eps)
        u0 = np.where(np.abs(a11) > 5 * sigma, r.mu, 1.0)
        tt = solve_eps_mu_multithickness(
            f, [{"s11": s[0], "s21": s[1], "s12": s[2], "s22": s[3], "length_m": L},
                {"s11": s2[0], "s21": s2[1], "s12": s2[2], "s22": s2[3], "length_m": L2}], GEOM, e0, u0)
        for k, (e, m) in {"nrw": (r.eps, r.mu), "iter_eps_mu": (it.eps, it.mu),
                          "nrw_eps_mu1": (r.eps_mu_only, np.ones(N)), "bj_eps": (bj.eps, bj.mu),
                          "two_thickness": (tt.eps, tt.mu)}.items():
            E[k].append(e)
            M[k].append(m)
    et = np.full(N, eps_t, complex)
    mt = np.full(N, mu_t, complex)
    std2, _ = crlb(f, GEOM, et, mt, L, sigma)
    std1, _ = crlb(f, GEOM, et, mt, L, sigma, params="eps")
    el = GEOM.gamma(f, et, mt).imag * L / np.pi
    rows = []
    per_f = {"f_hz": f, "beta_L_over_pi": el, "crlb_eps_real_std_eps_mu": std2[:, 0],
             "crlb_eps_real_std_eps_only": std1[:, 0], "crlb_mu_real_std": std2[:, 2]}
    near = np.abs(el - np.round(el)) < 0.05
    near &= np.round(el) >= 1
    for k in E:
        e = np.array(E[k])
        m = np.array(M[k])
        err = np.abs(e - et)
        rmse_f = np.sqrt(np.nanmean(np.abs(e.real - eps_t.real) ** 2, axis=0))
        per_f[f"rmse_eps_real_{k}"] = rmse_f
        rows.append({"study": label, "algorithm": k, "sigma": sigma, "n_seeds": n_seeds,
                     "rmse_eps_real_all": float(np.sqrt(np.nanmean((e.real - eps_t.real) ** 2))),
                     "rmse_eps_real_near_halfwave": float(np.sqrt(np.nanmean((e[:, near].real - eps_t.real) ** 2))),
                     "rmse_eps_real_away": float(np.sqrt(np.nanmean((e[:, ~near].real - eps_t.real) ** 2))),
                     "max_abs_err_eps": float(np.nanmax(err)),
                     "rmse_mu_real_all": float(np.sqrt(np.nanmean((m.real - np.real(mu_t)) ** 2))),
                     "max_abs_err_mu": float(np.nanmax(np.abs(m - mt))),
                     "mean_abs_bias_eps_real_away": float(np.abs(np.nanmean(e[:, ~near].real - eps_t.real))),
                     "nonfinite": int(np.count_nonzero(~np.isfinite(e)))})
    rows.append({"study": label, "algorithm": "CRLB_eps_mu", "rmse_eps_real_all": float(np.sqrt(np.mean(std2[:, 0] ** 2))),
                 "rmse_eps_real_near_halfwave": float(np.sqrt(np.mean(std2[near, 0] ** 2))),
                 "rmse_eps_real_away": float(np.sqrt(np.mean(std2[~near, 0] ** 2)))})
    rows.append({"study": label, "algorithm": "CRLB_eps_only", "rmse_eps_real_all": float(np.sqrt(np.mean(std1[:, 0] ** 2))),
                 "rmse_eps_real_near_halfwave": float(np.sqrt(np.mean(std1[near, 0] ** 2))),
                 "rmse_eps_real_away": float(np.sqrt(np.mean(std1[~near, 0] ** 2)))})
    return pd.DataFrame(rows), pd.DataFrame(per_f)


def figure(per_f, label, title):
    f = per_f["f_hz"] / 1e9
    fig, axs = plt.subplots(2, 1, figsize=(6.5, 5.0), sharex=True, constrained_layout=True)
    ax = axs[0]
    names = [("nrw", "classical NRW (ε, μ)"), ("iter_eps_mu", "iterative LSQ (ε, μ)"),
             ("two_thickness", "two-thickness LSQ (ε, μ)"), ("nrw_eps_mu1", "NRW ε with μ = 1"),
             ("bj_eps", "Baker-Jarvis-type iterative ε (μ = 1)")]
    for i, (k, lab) in enumerate(names):
        ax.semilogy(f, per_f[f"rmse_eps_real_{k}"], color=SERIES[i], ls=STYLES[i], label=lab)
    ax.semilogy(f, per_f["crlb_eps_real_std_eps_mu"], color="#111111", lw=0.9, ls=(0, (1, 1)),
                label="CRLB, simultaneous (ε, μ)")
    ax.semilogy(f, per_f["crlb_eps_real_std_eps_only"], color="#6b6b66", lw=0.9, ls=(0, (3, 1, 1, 1)),
                label="CRLB, ε only (μ = 1)")
    ax.set_ylabel("RMSE of ε' (20 noise draws)")
    ax.legend(fontsize=6.5, ncol=2, loc="upper left")
    ax.set_title(title, loc="left")
    ax = axs[1]
    ax.plot(f, per_f["beta_L_over_pi"], color=SERIES[0])
    for k in range(1, int(per_f["beta_L_over_pi"].max()) + 1):
        ax.axhline(k, color="#6b6b66", lw=0.6, ls=":")
    ax.set_ylabel("βL / π")
    ax.set_xlabel("Frequency (GHz)")
    FIG.mkdir(parents=True, exist_ok=True)
    p = FIG / f"singularity_{label}.png"
    fig.savefig(p)
    plt.close(fig)
    return p


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frames = []
    for label, eps_t, mu_t, L, L2, title in [
        ("dielectric", 4 - 0.004j, 1.0 + 0j, 20e-3, 13e-3, "Low-loss dielectric ε = 4 − 0.004j, μ = 1, L = 20 mm, σ = 1e-3"),
        ("magnetic", 4 - 0.004j, 1.5 - 0.003j, 20e-3, 13e-3,
         "Low-loss magnetic ε = 4 − 0.004j, μ = 1.5 − 0.003j, L = 20 mm, σ = 1e-3"),
    ]:
        df, per_f = run(eps_t, mu_t, L, L2, label=label)
        per_f.to_csv(OUT / f"singularity_{label}_per_frequency.csv", index=False)
        frames.append(df)
        figure(per_f, label, title)
        print(df.to_string(index=False))
    pd.concat(frames).to_csv(OUT / "singularity_summary.csv", index=False)


if __name__ == "__main__":
    main()
