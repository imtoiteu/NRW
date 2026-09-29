"""Smoke test for every adapter in research/adapters.

Synthetic data from the closed-form single-slab model of docs/conventions.md
(exp(+jwt), gamma = sqrt(kc^2 - k0^2 eps mu) principal root, z = mu gamma0/gamma,
Gamma = (z-1)/(z+1), T = exp(-gamma L), reference planes at the sample faces,
reciprocal & symmetric slab so S12 = S21, S22 = S11).

Cases (required): coax, 2-18 GHz, 401 pts, L = 5 mm
  A: eps = 4 - 0.1j,  mu = 1
  B: eps = 10 - 0.5j, mu = 2 - 0.3j
Extra (informative):
  wg_A:       WR-90 waveguide a = 22.86 mm, 8.2-12.4 GHz, 401 pts, L = 5 mm,
              case-A material (only for adapters that claim waveguide support);
  coax_B_L12: case-B material with L = 12 mm (m = 0 at 2 GHz, beta*L ~ 20 rad at
              18 GHz, i.e. several phase wraps) - probes multi-wrap tracking;
  coax_B_L30: case-B material with L = 30 mm (beta*L ~ 5.6 rad at 2 GHz, i.e. the
              correct branch index is m = 1 at the lowest frequency) - probes
              whether branch selection works beyond "unwrap from principal".

Metric: max over frequency of |x_ext - x_true| / |x_true| for eps, and for mu
when the adapter returns mu. Also the median, and the fraction of points with
relative error < 1e-6 (i.e. "exact" points).

Run:  .venv/bin/python -m research.adapters.smoke_test   (from repo root)
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    __package__ = "research.adapters"

from research.adapters import ADAPTERS, available  # noqa: E402

C0 = 299_792_458.0


def slab_sparams(freq_hz, eps, mu, L, geometry):
    f = np.asarray(freq_hz, float)
    k0 = 2 * np.pi * f / C0
    kc = 0.0 if geometry["kind"] == "coax" else np.pi / geometry["a_m"]
    eps = np.broadcast_to(np.asarray(eps, complex), f.shape)
    mu = np.broadcast_to(np.asarray(mu, complex), f.shape)
    gam = np.sqrt(kc**2 - k0**2 * eps * mu + 0j)
    gam = np.where(gam.real < 0, -gam, gam)  # principal root with Re >= 0
    gam0 = np.sqrt(kc**2 - k0**2 + 0j)
    gam0 = np.where(gam0.real < 0, -gam0, gam0)
    gam0 = np.where((gam0.real == 0) & (gam0.imag < 0), -gam0, gam0)  # = +j beta0
    z = mu * gam0 / gam
    G = (z - 1) / (z + 1)
    T = np.exp(-gam * L)
    D = 1 - G**2 * T**2
    s11 = G * (1 - T**2) / D
    s21 = T * (1 - G**2) / D
    return s11, s21, s21.copy(), s11.copy()


CASES = {
    "coax_A": dict(geometry={"kind": "coax"}, f=(2e9, 18e9), eps=4 - 0.1j, mu=1 + 0j),
    "coax_B": dict(geometry={"kind": "coax"}, f=(2e9, 18e9), eps=10 - 0.5j, mu=2 - 0.3j),
    "wg_A": dict(geometry={"kind": "waveguide", "a_m": 0.02286}, f=(8.2e9, 12.4e9),
                 eps=4 - 0.1j, mu=1 + 0j),
    # informative: m = 0 at 2 GHz but beta*L grows to ~20 rad (several 2*pi wraps)
    "coax_B_L12": dict(geometry={"kind": "coax"}, f=(2e9, 18e9), eps=10 - 0.5j, mu=2 - 0.3j,
                       L=12e-3),
    # informative: thick sample, phase branch m >= 1 already at the lowest frequency
    "coax_B_L30": dict(geometry={"kind": "coax"}, f=(2e9, 18e9), eps=10 - 0.5j, mu=2 - 0.3j,
                       L=30e-3),
}
L = 5e-3  # default sample length
N = 401


def relerr(x, ref):
    x = np.asarray(x, complex)
    return np.abs(x - ref) / np.abs(ref)


def run_case(name, fn, case):
    f = np.linspace(*case["f"], N)
    Lc = case.get("L", L)
    s11, s21, s12, s22 = slab_sparams(f, case["eps"], case["mu"], Lc, case["geometry"])
    t0 = time.time()
    try:
        out = fn(f, s11, s21, s12, s22, Lc, dict(case["geometry"]))
    except NotImplementedError as e:
        return {"status": "unsupported", "msg": str(e)}
    except Exception as e:  # report, do not crash the whole smoke test
        return {"status": "error", "msg": f"{type(e).__name__}: {e}",
                "tb": traceback.format_exc(limit=3)}
    dt = time.time() - t0
    assert set(out) >= {"eps", "mu", "notes"}, "adapter must return eps, mu, notes"
    eps = np.asarray(out["eps"], complex)
    assert eps.shape == f.shape, f"eps shape {eps.shape} != {f.shape}"
    r = {"status": "ok", "seconds": round(dt, 2)}
    e = relerr(eps, case["eps"])
    finite = np.isfinite(e)
    r["eps_max"] = float(np.max(e)) if finite.all() else float("inf")
    r["eps_median"] = float(np.nanmedian(np.where(finite, e, np.nan)))
    r["eps_frac_exact"] = float(np.mean(e < 1e-6))
    r["eps_nonfinite"] = int((~finite).sum())
    r["eps_passive_violations"] = int(np.sum(eps.imag > 1e-12))
    if out["mu"] is not None:
        mu = np.asarray(out["mu"], complex)
        assert mu.shape == f.shape
        m = relerr(mu, case["mu"])
        fin = np.isfinite(m)
        r["mu_max"] = float(np.max(m)) if fin.all() else float("inf")
        r["mu_median"] = float(np.nanmedian(np.where(fin, m, np.nan)))
        r["mu_frac_exact"] = float(np.mean(m < 1e-6))
    else:
        r["mu_max"] = None
    return r


def main(argv=None):
    names = (argv or sys.argv[1:]) or list(ADAPTERS)
    avail = available()
    results = {}
    for name in names:
        results[name] = {}
        if name not in avail:
            results[name]["import"] = "unavailable"
            print(f"{name:28s} UNAVAILABLE")
            continue
        for cname, case in CASES.items():
            r = run_case(name, ADAPTERS[name], case)
            results[name][cname] = r
            if r["status"] == "ok":
                mu_s = "" if r["mu_max"] is None else f" mu_max={r['mu_max']:.3g} mu_med={r['mu_median']:.3g}"
                print(f"{name:28s} {cname:7s} eps_max={r['eps_max']:.3g} eps_med={r['eps_median']:.3g} "
                      f"exact={r['eps_frac_exact']:.2f} passive_viol={r['eps_passive_violations']}{mu_s} ({r['seconds']}s)")
            else:
                print(f"{name:28s} {cname:7s} {r['status'].upper()}: {r['msg'][:150]}")
    out = Path(__file__).resolve().parent / "smoke_results.json"
    out.write_text(json.dumps(results, indent=1))
    print(f"wrote {out}")
    return results


if __name__ == "__main__":
    main()
