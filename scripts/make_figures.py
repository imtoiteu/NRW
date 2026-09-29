"""Comparison figures: ground truth vs classical NRW vs best external code vs final workflow.

"Best external" = the external implementation that extracts both eps and mu with the lowest
mean failure rate in results/benchmarks/summary_by_algorithm.csv (run run_benchmarks.py
first); if that file is missing, the external curve is omitted.

Output: results/figures/compare_<case>_<scenario>.png and the full diagnostic figure set of
the final workflow in results/figures/final_<case>/.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchlib import CASE_BY_NAME, INTERNAL, ROOT, external_algorithms, make_dataset  # noqa: E402

from nrw import extract_material_parameters  # noqa: E402
from nrw.plotting import SERIES, STYLES, TRUTH, _shade, apply_style, plot_result, plt  # noqa: E402

FIG = ROOT / "results" / "figures"
PANELS = [("ε'", lambda e, m: e.real), ("ε''", lambda e, m: -e.imag), ("μ'", lambda e, m: m.real),
          ("μ''", lambda e, m: -m.imag)]


def best_external():
    p = ROOT / "results" / "benchmarks" / "summary_by_algorithm.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    raw = pd.read_csv(ROOT / "results" / "benchmarks" / "raw.csv")
    eps_mu = raw.groupby("algorithm")["mu_assumed_1"].max()
    ext = df[df.algorithm.str.startswith("ext:")]
    ext = ext[ext.algorithm.map(lambda a: not bool(eps_mu.get(a, True)))]
    if ext.empty:
        return None
    return ext.sort_values(["mean_failure_rate", "median_rel_err_eps"]).iloc[0].algorithm


def final_workflow(d, case):
    method = "iterative_eps_mu" if case.magnetic else "baker_jarvis"
    return extract_material_parameters(d.f, d.s11, d.s21, d.L_assumed, s12=d.s12, s22=d.s22, method=method,
                                       noise_sigma=d.info.get("sigma", 1e-3), kk=True)


def compare(case_name, scenario, ext_name, seed=0):
    case = CASE_BY_NAME[case_name]
    d = make_dataset(case, scenario, seed)
    curves = [("classical NRW, m = 0 (textbook)", INTERNAL["nrw_principal"](d)),
              ("classical NRW + anchored branch", INTERNAL["nrw_anchored"](d))]
    if ext_name:
        try:
            curves.append((f"best external: {ext_name[4:]}", external_algorithms()[ext_name](d)))
        except Exception as e:  # noqa: BLE001
            print("external failed:", e)
    with np.errstate(all="ignore"):
        res = final_workflow(d, case)
    curves.append((f"final workflow ({res.method})", {"eps": res.epsilon_complex, "mu": res.mu_complex}))
    f = d.f / 1e9
    apply_style()
    npan = 4 if case.magnetic else 2
    fig, axs = plt.subplots(1 if npan == 2 else 2, 2, figsize=(7.0, 2.7 if npan == 2 else 4.8),
                            constrained_layout=True, squeeze=False)
    for ax, (lab, fn) in zip(axs.ravel(), PANELS[:npan]):
        _shade(ax, f, res.singular_region, "flagged ill-conditioned")
        tv = fn(d.eps_true, d.mu_true)
        span = np.ptp(tv) + 0.05 * np.max(np.abs(tv)) + 0.02
        lo, hi = np.min(tv) - 1.2 * span, np.max(tv) + 1.2 * span
        for i, (name, out) in enumerate(curves):
            y = fn(np.asarray(out["eps"]), np.asarray(out["mu"]))
            off = np.mean(~((y >= lo) & (y <= hi)))
            label = name + (f" [{off:.0%} of points off scale]" if off > 0.02 else "")
            ax.plot(f, y, color=SERIES[i], ls=STYLES[i], lw=1.2, label=label,
                    marker="o" if len(f) < 60 else None, ms=2.5)
        ax.plot(f, tv, **TRUTH)
        ax.set_ylabel(lab)
        ax.set_ylim(lo, hi)
        ax.legend(fontsize=5.5, loc="best")
    for ax in axs[-1]:
        ax.set_xlabel("Frequency (GHz)")
    fig.suptitle(f"{case_name} — {scenario} (SYNTHETIC)", fontsize=8.5, x=0.02, ha="left")
    FIG.mkdir(parents=True, exist_ok=True)
    p = FIG / f"compare_{case_name}_{scenario}.png"
    fig.savefig(p)
    plt.close(fig)
    out = FIG / f"final_{case_name}_{scenario}"
    truth = (d.eps_true, d.mu_true)
    plot_result(res, out, truth=truth)
    print("wrote", p)


def main():
    ext = best_external()
    print("best external:", ext)
    for case, sc in [("magnetodielectric_composite", "noise_1e-3"), ("lowloss_dielectric", "noise_1e-3"),
                     ("thick_dielectric", "high_start"), ("lowloss_magnetodielectric_thick", "noise_1e-3"),
                     ("magnetic_lorentz", "coarse_grid"), ("very_lossy_magnetic", "noise_1e-3")]:
        compare(case, sc, ext)


if __name__ == "__main__":
    main()
