"""Adapter: dmtang-jpg/onlyEE-GUI ``core.py`` (Python, NO LICENSE -> all rights
reserved; local evaluation only, code not redistributed).

``core.py`` depends only on numpy and is loaded by file path (the package
directory name contains a hyphen, and gui.py needs PyQt, which we never import).

* ``extract``        -> core.method_II(f, s11, s21, d): simultaneous eps, mu via
                         mu = sqrt(z^2 * n^2), eps = sqrt(n^2 / z^2) with the
                         square-root branch taken from an np.unwrap'ed phase.
* ``extract_method1``-> core.method_I(f, s11, s21, d): eps only (mu = 1), T from
                         the "T-quadratic" X_T = (S21^2 - S11^2 + 1)/(2 S21).
Both are TEM only (c = 3e8 hard-coded, no cutoff term) -> waveguide raises
NotImplementedError. No de-embedding is applied (planes already at faces).
Sign: T = exp(-gamma d) with gamma = ln(1/T)/d, eps = (gamma c/(j w))^2, i.e. the
exp(+jwt) convention; returned unchanged (smoke test checks the sign).
"""
from __future__ import annotations

import importlib.util

import numpy as np

from ._common import as_arrays, ext, require_tem

_NAME = "onlyEE-GUI"
_mod = None


def _import():
    global _mod
    if _mod is None:
        p = ext(_NAME) / "core.py"
        spec = importlib.util.spec_from_file_location("onlyee_core", p)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        _mod = m
    return _mod


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    require_tem(geometry, "onlyee")
    core = _import()
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    eps, mu = core.method_II(f, s11, s21, thickness_m)
    return {"eps": np.asarray(eps, complex), "mu": np.asarray(mu, complex),
            "notes": "onlyEE core.method_II (NRW, unwrap-continuous sqrt branches), TEM, c=3e8."}


def extract_method1(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    require_tem(geometry, "onlyee_method1")
    core = _import()
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    eps = core.method_I(f, s11, s21, thickness_m)
    return {"eps": np.asarray(eps, complex), "mu": None,
            "notes": "onlyEE core.method_I (mu=1, T-quadratic, unwrap), TEM, c=3e8."}
