"""Adapter: kasbrot123/nrw_method_eps_mu (Python, NO LICENSE -> all rights reserved;
used here only for local evaluation, code is not redistributed).

Calls ``nrw.nrw(S11, S21, f, d, lambda_cut)`` from the external repo unmodified
(n_best=None -> its per-frequency group-delay branch search over n in [-10, 10]).
Uses only S11/S21 (forward direction); S12/S22 are ignored, as in the original.

Geometry: lambda_cut = inf (coax/TEM; 1/inf = 0 in its formulas) or 2a (TE10).
Sign convention: the code derives eps, mu from T and Gamma only through
ln(1/T) and a principal square root; no conjugation is applied anywhere, so
exp(+jwt) input produces exp(+jwt) output. Returned unchanged (the smoke test
checks the sign).
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
    add_path(ext(_NAME))
    import nrw  # noqa: the external module
    return nrw


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    mod = _import()
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    lam_c = cutoff_wavelength(geometry)
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings(), \
            np.errstate(all="ignore"):
        warnings.simplefilter("ignore")
        eps, mu, _n = mod.nrw(s11, s21, f, thickness_m, lam_c)
    return {"eps": np.asarray(eps, complex), "mu": np.asarray(mu, complex),
            "notes": "kasbrot123 nrw(): NRW + per-point group-delay branch pick (n in -10..10), "
                     "S11/S21 only; no sign flip."}
