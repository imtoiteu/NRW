"""Adapter: shikouzhong/Nicholson-Ross-Weir-Algorithm (MATLAB, MIT) under Octave.

Calls the function ``NRW(newS11, newS21, f, sample, k, kc)`` (NRW.m)
unmodified, forward (S11,S21) and reverse (S22,S12), and averages the two
exactly as the repo's driver script ``Untitled.m`` does. The driver's
reference-plane search (``rotate.m`` + Global-Optimization-Toolbox ``ga``) is
skipped because our S-parameters are already referenced to the sample faces.
k = 2 pi f / c with c = 299792458 (the driver uses c = 3.0e8); kc = 0 for coax
(as in the driver: fc = 0) or pi/a for TE10.
Sign: NRW.m returns er = Re eps, ei = -Im eps, ur, ui = -Im mu, i.e. eps'' and
mu'' positive for passive exp(+jwt) data; we rebuild eps = er - j ei,
mu = ur - j ui (our convention).
"""
from __future__ import annotations

import numpy as np

from ._common import C0, as_arrays, ext, run_octave

_NAME = "Nicholson-Ross-Weir-Algorithm"


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    root = ext(_NAME)
    f, s11, s21, s12, s22 = as_arrays(freq_hz, s11, s21, s12, s22)
    kc = 0.0 if geometry["kind"] == "coax" else np.pi / float(geometry["a_m"])
    code = """
k = 2*pi*f/c0;
[er1, ei1, ur1, ui1] = NRW(s11, s21, f, L, k, kc);
[er2, ei2, ur2, ui2] = NRW(s22, s12, f, L, k, kc);
er = 0.5*(er1+er2); ei = 0.5*(ei1+ei2); ur = 0.5*(ur1+ur2); ui = 0.5*(ui1+ui2);
save('-v7', outfile, 'er', 'ei', 'ur', 'ui');
"""
    r = run_octave(code, {"f": f, "s11": s11, "s21": s21, "s12": s12, "s22": s22,
                          "L": thickness_m, "kc": kc, "c0": C0}, [root])
    eps = r["er"] - 1j * r["ei"]
    mu = r["ur"] - 1j * r["ui"]
    return {"eps": np.asarray(eps, complex), "mu": np.asarray(mu, complex),
            "notes": "shikouzhong NRW.m (Octave), fwd/rev average; first-jump-only +2pi branch fix."}
