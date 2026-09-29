"""Power and limits of the finite-band KK consistency diagnostic.

For causal ground truths and for extractions corrupted in known ways, report how often the
diagnostic passes. A good diagnostic passes correct extractions (low false-alarm rate)
and fails corrupted ones (power). Uncertainty = Cramer-Rao std of the real and loss parts at the noise level
(the linearised uncertainty the software uses when no Monte Carlo is run).

Output: results/studies/kk_power.csv, results/figures/kk_examples.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchlib import CASE_BY_NAME, GEOM, ROOT, make_dataset  # noqa: E402

from nrw.branch import select_branch  # noqa: E402
from nrw.classical import nrw  # noqa: E402
from nrw.conditioning import crlb  # noqa: E402
from nrw.iterative import solve_eps_mu_pointwise, solve_eps_only  # noqa: E402
from nrw.kk import kk_consistency  # noqa: E402
from nrw.materials import Debye, Lorentz  # noqa: E402
from nrw.plotting import SERIES, plt  # noqa: E402


def extract(d, magnetic, branch_offset=0):
    a11, a21 = (d.s11 + d.s22) / 2, (d.s21 + d.s12) / 2
    base = nrw(d.f, a11, a21, d.L_assumed, GEOM)
    b = select_branch("anchored", d.f, base.transmission, d.L_assumed, GEOM)
    r = nrw(d.f, a11, a21, d.L_assumed, GEOM, branch=b.m + branch_offset)
    if magnetic:
        it = solve_eps_mu_pointwise(d.f, d.s11, d.s21, d.L_assumed, GEOM, r.eps, r.mu, s12=d.s12, s22=d.s22,
                                    multistart=False)
        return it.eps, it.mu, r
    it = solve_eps_only(d.f, d.s11, d.s21, d.L_assumed, GEOM, r.eps_mu_only, s12=d.s12, s22=d.s22,
                        formulation="eq2.24+2.25")
    return it.eps, it.mu, r


def sigma_for(d, eps, mu, magnetic, noise):
    std, _ = crlb(d.f, GEOM, eps, mu, d.L_assumed, noise, params="eps_mu" if magnetic else "eps")
    se = (std[:, 0], std[:, 1])
    sm = (std[:, 2], std[:, 3]) if magnetic else None
    return se, sm


def test(f, val, sig, q):
    sr, si = sig
    ok = np.isfinite(val) & np.isfinite(sr) & np.isfinite(si)
    r = kk_consistency(f[ok], val[ok], sigma=sr[ok], sigma_imag=si[ok], quantity=q)
    return r.passed, r.max_abs_normalized_residual


def main(n_seeds=5):
    rows = []
    cases = ["debye_dielectric", "lorentz_dielectric", "magnetodielectric_composite", "magnetic_lorentz",
             "lowloss_magnetodielectric_thick"]
    variants = {
        "truth": ("noise_free", 0), "correct_extraction_noise1e-3": ("noise_1e-3", 0),
        "wrong_branch_plus1": ("noise_1e-3", 1), "thickness_error_1pct": ("thickness_error", 0),
        "ref_plane_error_50um": ("ref_plane_error", 0), "air_gap_10um": ("air_gap", 0),
        "holder_mismatch_51ohm": ("holder_mismatch", 0),
    }
    for cname in cases:
        case = CASE_BY_NAME[cname]
        mag = case.magnetic
        for vname, (sc, boff) in variants.items():
            seeds = range(n_seeds) if sc == "noise_1e-3" else [0]
            for seed in seeds:
                d = make_dataset(case, sc, seed)
                if vname == "truth":
                    eps, mu = d.eps_true, d.mu_true
                else:
                    with np.errstate(all="ignore"):
                        eps, mu, _ = extract(d, mag, boff)
                se, sm = sigma_for(d, d.eps_true, d.mu_true, mag, 1e-3)
                for q, val, sg in (("eps", eps, se), ("mu", mu, sm)):
                    if q == "mu" and not mag:
                        continue
                    try:
                        passed, mx = test(d.f, val, sg, q)
                    except ValueError:
                        passed, mx = False, np.nan
                    rel = np.nanmedian(np.abs(val - (d.eps_true if q == "eps" else d.mu_true))
                                       / np.abs(d.eps_true if q == "eps" else d.mu_true))
                    rows.append({"case": cname, "variant": vname, "quantity": q, "seed": seed, "passed": passed,
                                 "max_norm_residual": mx, "median_rel_error_of_input": rel})
        print(cname, flush=True)
    # synthetic non-causal / non-physical inputs
    f = np.linspace(2e9, 18e9, 401)
    synth = {
        "noncausal_mixed_models": Debye(2, 8, 1e9)(f).real - 1j * (-Lorentz(3, 2, 10e9, 2e9)(f).imag),
        "active_conjugate": np.conj(Lorentz(3, 2, 10e9, 2e9)(f)),
        "constant_complex_4-0.4j": np.full(len(f), 4 - 0.4j),
        "real_part_plus_0.3(f/18GHz)^2": Lorentz(3, 2, 10e9, 2e9)(f) + 0.3 * (f / 18e9) ** 2,
        "real_part_minus_0.3(f/18GHz)^2": Lorentz(3, 2, 10e9, 2e9)(f) - 0.3 * (f / 18e9) ** 2,
        "real_part_offset_+0.5": Lorentz(3, 2, 10e9, 2e9)(f) + 0.5,
    }
    for name, val in synth.items():
        for rel in (1e-2, 1e-3):
            r = kk_consistency(f, val, sigma=rel * np.abs(val.real), sigma_imag=rel * np.abs(val.real))
            rows.append({"case": "synthetic", "variant": f"{name} (sigma={rel:g}|eps'|)", "quantity": "eps",
                         "seed": 0, "passed": r.passed, "max_norm_residual": r.max_abs_normalized_residual})
    df = pd.DataFrame(rows)
    out = ROOT / "results" / "studies"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "kk_power.csv", index=False)
    summ = df.groupby(["case", "variant", "quantity"]).agg(pass_rate=("passed", "mean"),
                                                           median_max_norm_residual=("max_norm_residual", "median"),
                                                           median_rel_error_of_input=("median_rel_error_of_input", "median"),
                                                           n=("passed", "size")).reset_index()
    summ.to_csv(out / "kk_power_summary.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(summ.to_string(index=False))
    figure()


def figure():
    case = CASE_BY_NAME["magnetodielectric_composite"]
    d = make_dataset(case, "noise_1e-3", 0)
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8), constrained_layout=True)
    for ax, (boff, title) in zip(axs, [(0, "correct branch"), (1, "branch offset +1")]):
        with np.errstate(all="ignore"):
            eps, mu, _ = extract(d, True, boff)
        se, _ = sigma_for(d, d.eps_true, d.mu_true, True, 1e-3)
        r = kk_consistency(d.f, eps, sigma=se[0], sigma_imag=se[1], quantity="eps")
        ax.plot(d.f / 1e9, r.real_measured, color=SERIES[0], label="Re(ε) extracted")
        ax.plot(d.f / 1e9, r.real_predicted, color=SERIES[1], ls="--", label="KK from Im(ε) + passive tail")
        ax.plot(d.f / 1e9, d.eps_true.real, color="#111111", ls=(0, (4, 2)), lw=1.1, label="ground truth")
        ax.set_title(f"{title}: {'pass' if r.passed else 'FAIL'} (max |res|/σ = {r.max_abs_normalized_residual:.3g})",
                     loc="left")
        ax.set_xlabel("Frequency (GHz)")
        ax.set_ylabel("ε'")
    axs[0].legend(fontsize=6.5)
    p = ROOT / "results" / "figures" / "kk_examples.png"
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
