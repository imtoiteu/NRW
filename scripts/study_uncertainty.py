"""Monte Carlo uncertainty: coverage check and comparison with the Cramer-Rao bound.

For a well-conditioned magnetic case and a low-loss dielectric crossing half-wavelength points:
repeat R times: draw a synthetic "measurement" with errors from the SAME spec used by the
Monte Carlo (S magnitude/phase, thickness, reference planes), run the full pipeline with MC
uncertainty, and record whether the truth lies inside the reported 95 % interval.
Coverage near 95 % supports the MC intervals; large deviations are reported as such.

Output: results/studies/uncertainty_coverage.csv, results/figures/uncertainty_bands.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchlib import GEOM, ROOT  # noqa: E402

from nrw import extract_material_parameters  # noqa: E402
from nrw.forward import SampleSpec  # noqa: E402
from nrw.materials import Debye, Lorentz, as_model  # noqa: E402
from nrw.perturb import shift_reference_planes  # noqa: E402
from nrw.plotting import SERIES, TRUTH, plt  # noqa: E402
from nrw.uncertainty import UncertaintySpec, _perturb_pair  # noqa: E402

SPEC = UncertaintySpec(s11_mag_rel=0.005, s11_phase_deg=0.5, s21_mag_rel=0.005, s21_phase_deg=0.5,
                       thickness_m=10e-6, ref_plane_m=10e-6, n_draws=60, seed=7)
CASES = {
    "magnetodielectric_composite_2mm": (Debye(6.0, 4.0, 5e9), Lorentz(1.0, 2.5, 4e9, 6e9), 2e-3, "iterative_eps_mu"),
    "lowloss_dielectric_10mm_bj": (Debye(4.0, 0.2, 30e9), 1.0, 10e-3, "baker_jarvis"),
    "lowloss_dielectric_10mm_eps_mu": (Debye(4.0, 0.2, 30e9), 1.0, 10e-3, "iterative_eps_mu"),
}


def one_measurement(eps, mu, L, f, rng):
    s11, s21, s12, s22 = SampleSpec(eps=eps, mu=mu, length_m=L, geometry=GEOM).sparams(f)
    a11, a21 = _perturb_pair(s11, s21, SPEC, rng)
    a22, a12 = _perturb_pair(s22, s12, SPEC, rng)
    d1, d2 = SPEC.ref_plane_m * rng.standard_normal(2)
    a11, a21, a12, a22 = shift_reference_planes(f, a11, a21, a12, a22, GEOM, d1, d2)
    L_meas = L + SPEC.thickness_m * rng.standard_normal()  # the value the experimenter believes
    return a11, a21, a12, a22, L_meas


def main(R=20, N=81):
    f = np.linspace(2e9, 18e9, N)
    rows = []
    example = {}
    for name, (eps_m, mu_m, L, method) in CASES.items():
        et, mt = as_model(eps_m)(f), as_model(mu_m)(f)
        for rep in range(R):
            rng = np.random.default_rng(500 + rep)
            a11, a21, a12, a22, Lm = one_measurement(eps_m, mu_m, L, f, rng)
            spec = UncertaintySpec(**{**SPEC.to_dict(), "seed": 1000 + rep})
            with np.errstate(all="ignore"):
                res = extract_material_parameters(f, a11, a21, Lm, s12=a12, s22=a22, method=method,
                                                  uncertainty=spec, kk=False, noise_sigma=2e-3)
            u = res.uncertainty
            for q, tv in (("eps_real", et.real), ("eps_loss", -et.imag), ("mu_real", np.real(mt) * np.ones(N)),
                          ("mu_loss", -np.imag(mt) * np.ones(N))):
                if method == "baker_jarvis" and q.startswith("mu"):
                    continue
                inside = (tv >= u[q]["p2.5"]) & (tv <= u[q]["p97.5"])
                sing = res.singular_region
                rows.append({"case": name, "method": method, "rep": rep, "quantity": q,
                             "coverage_95_all": float(np.mean(inside)),
                             "coverage_95_well_conditioned": float(np.mean(inside[~sing])) if (~sing).any() else np.nan,
                             "coverage_95_ill_conditioned": float(np.mean(inside[sing])) if sing.any() else np.nan,
                             "median_mc_sigma": float(np.nanmedian(u[q]["sigma_percentile"])),
                             "max_mc_sigma": float(np.nanmax(u[q]["sigma_percentile"])),
                             "median_abs_error": float(np.nanmedian(np.abs(
                                 {"eps_real": res.epsilon_real, "eps_loss": res.epsilon_loss,
                                  "mu_real": res.mu_real, "mu_loss": res.mu_loss}[q] - tv))),
                             "branch_unstable_points": int(np.count_nonzero(res.flags["mc_branch_unstable"])),
                             "n_singular_points": int(np.count_nonzero(sing))})
            if rep == 0:
                example[name] = (res, et, mt)
        print(name, flush=True)
    df = pd.DataFrame(rows)
    out = ROOT / "results" / "studies"
    df.to_csv(out / "uncertainty_coverage.csv", index=False)
    summ = df.groupby(["case", "quantity"]).agg(
        coverage_all=("coverage_95_all", "mean"), coverage_well=("coverage_95_well_conditioned", "mean"),
        coverage_ill=("coverage_95_ill_conditioned", "mean"), median_mc_sigma=("median_mc_sigma", "median"),
        median_abs_error=("median_abs_error", "median"), reps=("rep", "size")).reset_index()
    summ.to_csv(out / "uncertainty_coverage_summary.csv", index=False)
    print(summ.to_string(index=False))
    figure(example)


def figure(example):
    fig, axs = plt.subplots(len(example), 2, figsize=(7.0, 2.3 * len(example)), constrained_layout=True,
                            squeeze=False)
    for row, (name, (res, et, mt)) in enumerate(example.items()):
        f = res.frequency / 1e9
        for ax, (q, lab, y, tv) in zip(axs[row], [("eps_real", "ε'", res.epsilon_real, et.real),
                                                   ("eps_loss", "ε''", res.epsilon_loss, -et.imag)]):
            u = res.uncertainty[q]
            sing = res.singular_region
            for k in np.flatnonzero(sing):
                ax.axvspan(f[max(k - 1, 0)], f[min(k + 1, len(f) - 1)], color="#e9e7e1", lw=0, zorder=0)
            ax.fill_between(f, u["p2.5"], u["p97.5"], color=SERIES[0], alpha=0.2, lw=0, label="95 % MC interval")
            ax.plot(f, y, color=SERIES[0], label=f"extracted ({res.method})")
            ax.plot(f, tv, **TRUTH)
            ax.set_ylabel(lab)
            lo, hi = np.nanpercentile(np.concatenate([tv, u["p2.5"], u["p97.5"]]), [3, 97])
            pad = 0.3 * (hi - lo + 1e-3)
            ax.set_ylim(lo - pad, hi + pad)
            ax.set_title(name, loc="left", fontsize=7.5)
        axs[row, 0].legend(fontsize=6)
    for ax in axs[-1]:
        ax.set_xlabel("Frequency (GHz)")
    fig.savefig(ROOT / "results" / "figures" / "uncertainty_bands.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
