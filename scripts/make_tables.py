"""Generate the markdown tables of docs/benchmark_report.md from the result CSVs.

Every table in the report sits between markers
    <!-- TABLE:name -->  ...  <!-- /TABLE:name -->
and is regenerated here, so no number in the report is typed by hand.
Tables are also written to results/tables/<name>.md.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_benchmarks import EXT_SCENARIOS  # noqa: E402

B = ROOT / "results" / "benchmarks"
S = ROOT / "results" / "studies"
T = ROOT / "results" / "tables"
REPORT = ROOT / "docs" / "benchmark_report.md"


def fmt(x, digits=3):
    if isinstance(x, (bool, np.bool_)):
        return "yes" if x else "no"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(xf):
        return "–"
    if xf == 0:
        return "0"
    if abs(xf) >= 1e4 or abs(xf) < 1e-3:
        return f"{xf:.{digits - 1}e}"
    return f"{xf:.{digits}g}"


def md(df: pd.DataFrame, index=False) -> str:
    if index:
        df = df.reset_index()
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(fmt(v) for v in r.values) + " |")
    return "\n".join(out)


def algorithm_table(raw: pd.DataFrame, magnetic: bool) -> pd.DataFrame:
    """Comparison on the scenario set shared by all algorithms (fair internal vs external)."""
    from benchlib import CASE_BY_NAME
    sub = raw[raw.scenario.isin(EXT_SCENARIOS) & (raw.seed < 3)]
    sub = sub[sub.case.map(lambda c: CASE_BY_NAME[c].magnetic) == magnetic]
    g = sub.groupby("algorithm")
    t = pd.DataFrame({
        "mu = 1 assumed": g.mu_assumed_1.max(),
        "runs": g.size(),
        "median rel. err ε": g.median_rel_err_eps.median(),
        "median rel. err μ": g.median_rel_err_mu.median(),
        "mean failure rate": g.failure_rate.mean(),
        "runs with failures": g.failure_rate.apply(lambda x: np.mean(x > 0)),
        "mean branch error rate": g.branch_error_rate.mean(),
        "passivity viol. (>3σ CRLB)": g.passivity_violations.sum(),
        "non-finite values": g.nonfinite_count.sum(),
        "crashed runs": g.error.apply(lambda x: int(np.sum(x.fillna("").astype(str).str.len() > 0))),
        "median runtime (s)": g.runtime_s.median(),
    })
    if magnetic:
        t = t.drop(columns=[])
    return t.sort_values(["mu = 1 assumed", "mean failure rate", "median rel. err ε"])


def main():
    T.mkdir(parents=True, exist_ok=True)
    tables = {}
    raw = pd.read_csv(B / "raw.csv")
    tables["alg_nonmagnetic"] = md(algorithm_table(raw, False), index=True)
    tables["alg_magnetic"] = md(algorithm_table(raw, True), index=True)
    internal = raw[~raw.algorithm.str.startswith("ext:")]
    sc = internal.pivot_table(index="algorithm", columns="scenario", values="failure_rate", aggfunc="mean")
    tables["internal_by_scenario"] = md(sc.round(4), index=True)
    from benchlib import CASE_BY_NAME
    mag = internal[internal.case.map(lambda c: CASE_BY_NAME[c].magnetic)]
    tables["internal_by_case_magnetic"] = md(
        mag.pivot_table(index="algorithm", columns="case", values="failure_rate", aggfunc="mean").round(4), index=True)
    non = internal[~internal.case.map(lambda c: CASE_BY_NAME[c].magnetic)]
    tables["internal_by_case_nonmagnetic"] = md(
        non.pivot_table(index="algorithm", columns="case", values="failure_rate", aggfunc="mean").round(4), index=True)
    rm = raw[raw.algorithm.isin(["nrw_anchored", "iterative_eps_mu", "regularized_eps_mu", "bj_eps_iterative"])]
    rmse = rm.groupby(["algorithm", "case"])[["rmse_eps_real", "rmse_eps_loss", "rmse_mu_real", "rmse_mu_loss"]].median()
    tables["rmse_by_case"] = md(rmse, index=True)

    if (S / "singularity_summary.csv").exists():
        s = pd.read_csv(S / "singularity_summary.csv")
        tables["singularity"] = md(s[["study", "algorithm", "rmse_eps_real_all", "rmse_eps_real_near_halfwave",
                                      "rmse_eps_real_away", "max_abs_err_eps", "rmse_mu_real_all",
                                      "max_abs_err_mu"]])
    if (S / "branch_summary_any_error_fraction.csv").exists():
        a = pd.read_csv(S / "branch_summary_any_error_fraction.csv", index_col=0)
        info = pd.read_csv(S / "branch_conditions.csv", index_col=0)
        a.insert(0, "N", info["n_freq"])
        a.insert(1, "m range", info["m_true_min"].astype(str) + "–" + info["m_true_max"].astype(str))
        a.insert(2, "max true step (rad)", info["true_max_phase_step_rad"])
        tables["branch_any"] = md(a, index=True)
        e = pd.read_csv(S / "branch_summary_error_rate.csv", index_col=0)
        tables["branch_rate"] = md(e, index=True)
    if (S / "bj_formulation.csv").exists():
        d = pd.read_csv(S / "bj_formulation.csv")
        d["variant"] = d.formulation + " (β=" + d.beta.astype(str) + ")"
        t = d.groupby(["planes", "variant"])[["rmse_eps_real", "rmse_eps_loss", "max_abs_err", "not_converged"]].median()
        tables["bj_formulation"] = md(t, index=True)
    if (S / "kk_power_summary.csv").exists():
        k = pd.read_csv(S / "kk_power_summary.csv")
        piv = k[k.case != "synthetic"].pivot_table(index=["variant"], columns=["case", "quantity"], values="pass_rate")
        piv.columns = [f"{c[0]} ({c[1]})" for c in piv.columns]
        tables["kk_power"] = md(piv, index=True)
        syn = k[k.case == "synthetic"][["variant", "pass_rate", "median_max_norm_residual"]]
        tables["kk_synthetic"] = md(syn)
    if (S / "uncertainty_coverage_summary.csv").exists():
        u = pd.read_csv(S / "uncertainty_coverage_summary.csv")
        tables["uncertainty"] = md(u)
    for name, txt in tables.items():
        (T / f"{name}.md").write_text(txt + "\n")
    if REPORT.exists():
        rep = REPORT.read_text()
        for name, txt in tables.items():
            pat = re.compile(rf"(<!-- TABLE:{name} -->)(.*?)(<!-- /TABLE:{name} -->)", re.S)
            rep = pat.sub(lambda m: f"{m.group(1)}\n{txt}\n{m.group(3)}", rep)
        missing = [m for m in re.findall(r"<!-- TABLE:(\w+) -->", rep) if m not in tables]
        if missing:
            print("WARNING: no data for tables", missing)
        REPORT.write_text(rep)
    print("tables:", ", ".join(tables))


if __name__ == "__main__":
    main()
