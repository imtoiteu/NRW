"""Adapter: bumbeneciconstantinbogdan/S-Param-Inversion ``nrw.py`` (Python/PyTorch, MIT).

Calls ``nrw.nrw_inverse(d, a, frequencies, S11, S21)`` unmodified (the classical
NRW baseline the repo benchmarks its MLP against; the MLP itself is not used).
Inputs are converted to complex128 torch tensors and back.

Geometry: TE10 with kt = pi/a. For coax we pass a = inf so kt = 0 (TEM); the
formulas are otherwise identical. Principal log only (no unwrap / no branch
index), so results are expected to fail once the electrical length exceeds pi.
Sign: P = exp(-j beta d), eps = eps' - j eps'' (exp(+jwt)); returned unchanged.
"""
from __future__ import annotations

import numpy as np

from ._common import as_arrays, ext

_NAME = "S-Param-Inversion"


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    # NB: the module name "nrw" clashes with kasbrot123's nrw.py; load by path.
    import importlib.util
    import torch
    p = ext(_NAME) / "nrw.py"
    spec = importlib.util.spec_from_file_location("sparam_inversion_nrw", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    a = float("inf") if geometry["kind"] == "coax" else float(geometry["a_m"])
    if geometry["kind"] not in ("coax", "waveguide"):
        raise ValueError(geometry)
    t = lambda x: torch.as_tensor(x, dtype=torch.complex128)
    eps, mu = mod.nrw_inverse(thickness_m, a, torch.as_tensor(f, dtype=torch.float64),
                              t(s11), t(s21))
    return {"eps": eps.numpy().astype(complex), "mu": mu.numpy().astype(complex),
            "notes": "S-Param-Inversion nrw_inverse (NRW, principal log, m=0 only); coax via a=inf."}
