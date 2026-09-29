"""Adapter: marcdmc/meep-metamaterials ``retrieval/retrieval.py`` (Python, MIT).

Calls ``eff_parameters(freqs, d, S11, S21)`` unmodified with its defaults
(branch=0, continuity=False). This is the Smith/Chen-style metamaterial
retrieval: z = sqrt(((1+S11)^2-S21^2)/((1-S11)^2-S21^2)) (principal root),
e^{ink d} = S21/(1 - S11 (z-1)/(z+1)), n = [Im ln + 2 pi m - j Re ln]/(k d).

Interface translation:
  * The module imports ``meep`` (FDTD simulator, not pip-installable) and, via
    ``meep_metamaterials.metamaterials``, meep.mpb/meep.geom/IPython/pandas.
    ``eff_parameters`` uses none of them; permissive stub modules are inserted
    for meep, meep.mpb, meep.geom and IPython(.display) only while the module
    is being imported, and removed from sys.modules afterwards.
  * MEEP units: k = 2 pi * freqs with freqs in 1/length; we pass freqs = f/c0
    [1/m] and d in metres.
  * Time convention: MEEP/physics exp(-i w t) (S21 ~ e^{+ink d}, Im n > 0).
    Ours is exp(+jwt). We therefore pass conj(S11), conj(S21) and return
    conj(eps), conj(mu). No other change.
  * TEM/normal incidence only (no cutoff term) -> waveguide NotImplementedError.
"""
from __future__ import annotations

import importlib.util
import sys

import numpy as np

from ._common import C0, add_path, as_arrays, ext, require_tem, temporary_stubs

_NAME = "meep-metamaterials"
_mod = None


def _import():
    global _mod
    if _mod is not None:
        return _mod
    import matplotlib
    matplotlib.use("Agg", force=False)
    root = ext(_NAME)
    add_path(root)
    p = root / "meep_metamaterials" / "retrieval" / "retrieval.py"
    with temporary_stubs("meep", "meep.mpb", "meep.geom", "IPython", "IPython.display"):
        spec = importlib.util.spec_from_file_location("meep_mm_retrieval", p)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        # the package modules imported while stubs were active are dropped too
        for k in [k for k in sys.modules if k.startswith("meep_metamaterials")]:
            sys.modules.pop(k, None)
    _mod = m
    return m


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    require_tem(geometry, "meep_retrieval")
    m = _import()
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    with np.errstate(all="ignore"):
        r = m.eff_parameters(f / C0, thickness_m, np.conj(s11), np.conj(s21))
    return {"eps": np.conj(np.asarray(r["eps"], complex)),
            "mu": np.conj(np.asarray(r["mu"], complex)),
            "notes": "meep-metamaterials eff_parameters(branch=0, continuity=False); "
                     "exp(-iwt) code: inputs and outputs conjugated."}
