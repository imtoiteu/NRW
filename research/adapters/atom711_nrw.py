"""Adapter: atom711/nrw-material-extraction (Python, NO LICENSE -> all rights
reserved; local evaluation only, code not redistributed).

Calls ``NRWExtractionEngine.execute_inversion_pipeline(s11, s21)`` unmodified.
The engine's constructor builds its own ``linspace`` frequency grid; the only
interface edit is that after construction we overwrite the instance attributes
``freqs``, ``omega`` and ``k_0`` with the caller's grid (the pipeline uses only
``self.k_0`` and ``self.d``).

Method: eps = -(ln P / (k0 d))^2 with P the NRW transmission term and phase
np.unwrap'ed from the principal value at the first frequency; assumes mu = 1
and TEM (no cutoff term) -> waveguide raises NotImplementedError.
Output: (eps', tan d) with eps'' = -Im(eps); we rebuild eps = eps'(1 - j tan d),
which is our convention (the engine itself uses eps = eps' - j eps'').
"""
from __future__ import annotations

import numpy as np

from ._common import add_path, as_arrays, ext, require_tem

_NAME = "nrw-material-extraction"


def _import():
    import matplotlib
    matplotlib.use("Agg", force=False)
    add_path(ext(_NAME) / "src")
    from NRWExtractionEngine import NRWExtractionEngine
    return NRWExtractionEngine


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    require_tem(geometry, "atom711_nrw")
    Engine = _import()
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    eng = Engine(f_start=f[0], f_stop=f[-1], num_points=len(f), thickness=thickness_m)
    eng.freqs = f
    eng.omega = 2 * np.pi * f
    eng.k_0 = eng.omega / eng.c
    with np.errstate(all="ignore"):
        eps_p, tand = eng.execute_inversion_pipeline(s11, s21)
    eps = eps_p * (1 - 1j * tand)
    return {"eps": np.asarray(eps, complex), "mu": None,
            "notes": "atom711 NRWExtractionEngine: eps=-(lnP/k0d)^2, mu=1 assumed, TEM only, "
                     "np.unwrap phase; returned (eps', tan d) recombined as eps'(1-j tan d)."}
