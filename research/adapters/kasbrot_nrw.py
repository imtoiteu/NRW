"""Adapter: kasbrot123/nrw_method_eps_mu (Python, NO LICENSE -> all rights reserved;
used here only for local evaluation, code is not redistributed).

Calls ``nrw.nrw(S11, S21, f, d, lambda_cut)`` from the external repo unmodified
(n_best=None -> its per-frequency group-delay branch search over n in [-10, 10]).
Uses only S11/S21 (forward direction); S12/S22 are ignored, as in the original.

Geometry: lambda_cut = inf (coax/TEM; 1/inf = 0 in its formulas) or 2a (TE10).
Sign convention: the code is written for exp(-iwt) data. Evidence in the code:
tau_meas = -(1/2pi) d(unwrap arg T)/df and the branch is chosen by minimising
|Re tau_calc + Re tau_meas| with tau_calc > 0 for a normal material, which
requires d(arg T)/df > 0, i.e. T = exp(+j beta L) (exp(-iwt)); its demo also
uses eps = 8.4(1+0.05j) (Im > 0 for loss). With exp(+jwt) input the group-delay
search picks the wrong branch everywhere (verified: n_best=0 on the same data is
exact). The adapter therefore passes conj(S11), conj(S21) and returns
conj(eps), conj(mu).
"""
from __future__ import annotations

import contextlib
import io
import warnings

import numpy as np

from ._common import add_path, as_arrays, cutoff_wavelength, ext

_NAME = "nrw_method_eps_mu"


def _import():
    import matplotlib
    matplotlib.use("Agg", force=False)
    # Load by file path under a unique name: the external module is called "nrw.py" and
    # would otherwise collide with (or be shadowed by) this project's own `nrw` package.
    import importlib.util
    root = ext(_NAME)
    add_path(root)  # for any sibling imports inside the external module
    spec = importlib.util.spec_from_file_location("kasbrot123_nrw_external", root / "nrw.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    mod = _import()
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    lam_c = cutoff_wavelength(geometry)
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings(), \
            np.errstate(all="ignore"):
        warnings.simplefilter("ignore")
        eps, mu, _n = mod.nrw(np.conj(s11), np.conj(s21), f, thickness_m, lam_c)
    return {"eps": np.conj(np.asarray(eps, complex)), "mu": np.conj(np.asarray(mu, complex)),
            "notes": "kasbrot123 nrw(): NRW + per-point group-delay branch pick (n in -10..10), "
                     "S11/S21 only; exp(-iwt) code: inputs/outputs conjugated."}
