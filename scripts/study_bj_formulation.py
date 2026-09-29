"""Choice of NIST eps-only formulation and weight beta (TN 1341 eqs. 2.38, 2.25, 2.24+2.25).

Non-magnetic benchmark cases x noise levels x seeds; metric = RMSE of eps' and eps''.
Also a reference-plane scenario: planes displaced asymmetrically by unknown amounts with the
holder length known (only the invariant forms are then applicable without de-embedding).

Output: results/studies/bj_formulation.csv (+ summary printed)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchlib import CASES, GEOM, ROOT, make_dataset  # noqa: E402

from nrw.branch import select_branch  # noqa: E402
from nrw.classical import nrw  # noqa: E402
from nrw.iterative import solve_eps_only  # noqa: E402
from nrw.perturb import shift_reference_planes  # noqa: E402

VARIANTS = [("eq2.38", 0.0), ("eq2.38", 1.0), ("eq2.38", 10.0), ("eq2.25", 0.0), ("eq2.24+2.25", 0.0),
            ("direct", 0.0)]


def main(n_seeds=5):
    rows = []
    cases = [c for c in CASES if not c.magnetic]
    for case in cases:
        for sc in ("noise_1e-3", "noise_1e-2", "phase_noise", "mag_noise"):
            for seed in range(n_seeds):
                d = make_dataset(case, sc, seed)
                for planes in ("at_faces", "unknown_asymmetric"):
                    s11, s21, s12, s22 = d.s11, d.s21, d.s12, d.s22
                    if planes == "unknown_asymmetric":
                        # true planes 3 mm and 1 mm from the faces; only L_air - L = 4 mm is known
                        s11, s21, s12, s22 = shift_reference_planes(d.f, s11, s21, s12, s22, GEOM, 3e-3, 1e-3)
                    # initial estimate: NRW/Boughriet with branch (after nominal symmetric de-embedding)
                    g0 = GEOM.gamma0(d.f)
                    rr = np.exp(-g0 * (2e-3 if planes != "at_faces" else 0.0))
                    a11 = (s11 / rr**2 + s22 / rr**2) / 2
                    a21 = (s21 + s12) / 2 / rr**2
                    base = nrw(d.f, a11, a21, d.L_assumed, GEOM)
                    b = select_branch("anchored", d.f, base.transmission, d.L_assumed, GEOM)
                    e0 = nrw(d.f, a11, a21, d.L_assumed, GEOM, branch=b.m).eps_mu_only
                    for form, beta in VARIANTS:
                        kw = dict(formulation=form, beta=beta)
                        if planes == "unknown_asymmetric":
                            if form in ("eq2.38", "direct"):
                                kw.update(d1=2e-3, d2=2e-3)  # best available guess: symmetric split
                            else:
                                kw.update(l_air_minus_l=4e-3)
                        with np.errstate(all="ignore"):
                            it = solve_eps_only(d.f, s11, s21, d.L_assumed, GEOM, e0, s12=s12, s22=s22, **kw)
                        err = it.eps - d.eps_true
                        rows.append({"case": case.name, "scenario": sc, "seed": seed, "planes": planes,
                                     "formulation": form, "beta": beta,
                                     "rmse_eps_real": float(np.sqrt(np.nanmean(err.real ** 2))),
                                     "rmse_eps_loss": float(np.sqrt(np.nanmean(err.imag ** 2))),
                                     "max_abs_err": float(np.nanmax(np.abs(err))),
                                     "not_converged": int(np.count_nonzero(~it.converged))})
        print(case.name, flush=True)
    df = pd.DataFrame(rows)
    out = ROOT / "results" / "studies"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "bj_formulation.csv", index=False)
    df["variant"] = df.formulation + " beta=" + df.beta.astype(str)
    summ = df.groupby(["planes", "case", "variant"])[["rmse_eps_real", "rmse_eps_loss"]].median().unstack("variant")
    summ.to_csv(out / "bj_formulation_summary.csv")
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.float_format", "{:.2e}".format):
        print(df.groupby(["planes", "variant"])[["rmse_eps_real", "rmse_eps_loss", "max_abs_err", "not_converged"]]
              .median().to_string())
        print(summ["rmse_eps_real"].to_string())


if __name__ == "__main__":
    main()
