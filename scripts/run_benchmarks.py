"""Main benchmark: every algorithm x case x scenario (x seeds for stochastic scenarios).

Usage:
    python scripts/run_benchmarks.py [--quick] [--no-external] [--seeds 5]

Outputs (results/benchmarks/):
    raw.csv                      one row per (algorithm, case, scenario, seed)
    summary_by_algorithm.csv     aggregated metrics per algorithm
    summary_by_scenario.csv      median failure rate per algorithm x scenario
    summary_by_case.csv          median relative error per algorithm x case
    environment.json             package versions, seeds, timestamp
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchlib import CASES, INTERNAL, ROOT, SCENARIOS, STOCHASTIC, external_algorithms, run_one  # noqa: E402

OUT = ROOT / "results" / "benchmarks"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="subset for smoke testing")
    ap.add_argument("--no-external", action="store_true")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--only", default=None, help="comma-separated algorithm names")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    algs = dict(INTERNAL)
    if not args.no_external:
        algs.update(external_algorithms())
    if args.only:
        keep = set(args.only.split(","))
        algs = {k: v for k, v in algs.items() if k in keep}
    cases = CASES[:3] if args.quick else CASES
    scenarios = SCENARIOS[:3] if args.quick else SCENARIOS
    rows = []
    for alg, fn in algs.items():
        for case in cases:
            for sc in scenarios:
                seeds = range(args.seeds) if sc in STOCHASTIC else [0]
                for seed in seeds:
                    row, _, _ = run_one(alg, fn, case, sc, seed)
                    rows.append(row)
        print(f"{alg:28s} done ({sum(r['algorithm'] == alg for r in rows)} runs)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "raw.csv", index=False)

    agg = df.groupby("algorithm").agg(
        median_rel_err_eps=("median_rel_err_eps", "median"),
        median_rel_err_mu=("median_rel_err_mu", "median"),
        mean_failure_rate=("failure_rate", "mean"),
        runs_with_any_failure=("failure_rate", lambda x: float(np.mean(x > 0))),
        passivity_violations=("passivity_violations", "sum"),
        passivity_violations_raw=("passivity_violations_raw", "sum"),
        unphysical_count=("unphysical_count", "sum"),
        nonfinite_count=("nonfinite_count", "sum"),
        mean_branch_error_rate=("branch_error_rate", "mean"),
        median_runtime_s=("runtime_s", "median"),
        crashed_runs=("error", lambda x: int(np.sum(x.astype(str).str.len() > 0))),
        runs=("error", "size"),
    ).reset_index()
    agg.to_csv(OUT / "summary_by_algorithm.csv", index=False)
    df.pivot_table(index="algorithm", columns="scenario", values="failure_rate", aggfunc="mean").to_csv(
        OUT / "summary_by_scenario.csv")
    df.pivot_table(index="algorithm", columns="case", values="failure_rate", aggfunc="mean").to_csv(
        OUT / "summary_by_case.csv")
    import scipy
    env = {"timestamp_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(),
           "numpy": np.__version__, "scipy": scipy.__version__, "pandas": pd.__version__,
           "platform": platform.platform(), "seeds_per_stochastic_scenario": args.seeds,
           "algorithms": list(algs), "cases": [c.name for c in cases], "scenarios": scenarios,
           "failure_definition": "fraction of frequencies with |eps-eps_true|/|eps_true| > 10% or "
                                 "|mu-mu_true|/|mu_true| > 10% or non-finite"}
    (OUT / "environment.json").write_text(json.dumps(env, indent=2))
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(agg.to_string(index=False))


if __name__ == "__main__":
    main()
