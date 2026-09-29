"""Registry of thin adapters around surveyed open-source NRW-type codes.

Every adapter has the signature

    extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict

freq_hz: 1-D float ndarray (Hz); s-params: complex ndarrays referenced to the
sample faces; thickness_m: float; geometry: {"kind": "coax"} or
{"kind": "waveguide", "a_m": 0.02286}. Returns {"eps": complex ndarray,
"mu": complex ndarray or None, "notes": str} in OUR convention
(exp(+jwt), eps = eps' - j eps'', passive => Im eps <= 0).

Adapters that do not support a geometry raise NotImplementedError.
The external code must first be fetched with scripts/fetch_external.sh.
See docs/github_survey.md for what each external code does (and its defects).
"""
from __future__ import annotations

import importlib

# name -> (module, function, requires, geometries)
_SPEC = {
    "permittivitycalc_nrw": ("permittivitycalc", "extract", "python", ("coax", "waveguide")),
    "permittivitycalc_nni": ("permittivitycalc", "extract_nni", "python", ("coax", "waveguide")),
    "kasbrot_nrw": ("kasbrot_nrw", "extract", "python", ("coax", "waveguide")),
    "atom711_nrw": ("atom711_nrw", "extract", "python", ("coax",)),
    "onlyee_method2": ("onlyee", "extract", "python", ("coax",)),
    "onlyee_method1": ("onlyee", "extract_method1", "python", ("coax",)),
    "sparam_inversion_nrw": ("sparam_inversion", "extract", "python+torch", ("coax", "waveguide")),
    "meep_retrieval": ("meep_retrieval", "extract", "python", ("coax",)),
    "musho_cpp": ("musho_cpp", "extract", "g++", ("coax",)),
    "scatter_opt_nrw": ("scatter_opt", "extract", "octave", ("coax", "waveguide")),
    "scatter_opt_gd": ("scatter_opt", "extract_groupdelay", "octave", ("coax", "waveguide")),
    "shikouzhong_nrw": ("shikouzhong_nrw", "extract", "octave", ("coax", "waveguide")),
    "njchorda_nrw": ("njchorda_nrw", "extract", "octave", ("coax", "waveguide")),
    "bj_tem_5g": ("bj_tem_5g", "extract", "octave", ("coax",)),
    "ermakov_wg": ("ermakov_wg", "extract", "octave", ("coax", "waveguide")),
}

GEOMETRIES = {k: v[3] for k, v in _SPEC.items()}


def _lazy(modname, fn):
    def call(*a, **k):
        m = importlib.import_module(f"{__name__}.{modname}")
        return getattr(m, fn)(*a, **k)
    call.__name__ = f"{modname}.{fn}"
    return call


ADAPTERS = {name: _lazy(mod, fn) for name, (mod, fn, _req, _g) in _SPEC.items()}


def _check(name) -> bool:
    import shutil
    mod, fn, req, _ = _SPEC[name]
    try:
        m = importlib.import_module(f"{__name__}.{mod}")
        getattr(m, fn)
        if req == "octave" and shutil.which("octave-cli") is None:
            return False
        if req == "g++" and shutil.which("g++") is None:
            return False
        if req == "python+torch":
            importlib.import_module("torch")
        # external repo present? (each module defines _NAME)
        from ._common import ext
        ext(m._NAME)
        if hasattr(m, "_import") and req.startswith("python"):
            m._import()
        return True
    except Exception:
        return False


def available() -> list[str]:
    """Names of adapters whose module, external checkout and runtime all load."""
    return [n for n in _SPEC if _check(n)]
