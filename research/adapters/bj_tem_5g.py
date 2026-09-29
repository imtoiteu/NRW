"""Adapter: tomaasaybr02/5G-New-Material-EM-Characterization-Techniques,
``MATLAB/NRW_BJ_TEM_2.m`` (MATLAB, MIT) under Octave.

Calls ``[n, mu, eps, Z] = NRW_BJ_TEM_2(Spar_sample, Spar_air, freq, c, L, Lair, NN)``
unmodified. This is the Baker-Jarvis (NIST TN 1355-R) reference-plane-invariant
closed form: A = S11 S22/(S21 S12), B = e^{2 gamma0 (Lair-L)} (S21 S12 - S11 S22),
Gamma^2 from a quadratic, then P from S21/S21_air, n via unwrap(angle(1/P))+2 pi NN.

Interface translation:
  * Spar_sample columns are [S11 S12 S21 S22] (the function reads S11=(:,1),
    S12=(:,2), S21=(:,3), S22=(:,4)).
  * The method needs the empty-holder S21 (Spar_air(:,3)) and the holder length
    Lair. With reference planes at the sample faces the "holder" is the sample
    length itself: Lair = L and S21_air = exp(-j k0 L) (ideal empty line),
    built with the same c that is passed to the function (c0 = 299792458).
  * NN = 0. TEM only (gamma_1 = j w / c) -> waveguide NotImplementedError.
Sign: n = ln(1/P)/(j k0 L) -> Im n < 0 for passive exp(+jwt) data; returned
unchanged.
"""
from __future__ import annotations

import numpy as np

from ._common import C0, as_arrays, ext, require_tem, run_octave

_NAME = "5G-New-Material-EM-Characterization-Techniques"


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    require_tem(geometry, "bj_tem_5g")
    root = ext(_NAME) / "MATLAB"
    f, s11, s21, s12, s22 = as_arrays(freq_hz, s11, s21, s12, s22)
    s21_air = np.exp(-1j * 2 * np.pi * f / C0 * thickness_m)
    code = """
Spar_sample = [s11(:) s12(:) s21(:) s22(:)];
Spar_air = [zeros(size(s21a(:))) s21a(:) s21a(:) zeros(size(s21a(:)))];
[n, mu, eps, Z] = NRW_BJ_TEM_2(Spar_sample, Spar_air, f(:), c0, L, L, 0);
save('-v7', outfile, 'eps', 'mu');
"""
    r = run_octave(code, {"f": f, "s11": s11, "s21": s21, "s12": s12, "s22": s22,
                          "s21a": s21_air, "L": thickness_m, "c0": C0}, [root])
    return {"eps": np.atleast_1d(r["eps"]).astype(complex),
            "mu": np.atleast_1d(r["mu"]).astype(complex),
            "notes": "5G repo NRW_BJ_TEM_2 (Baker-Jarvis ref-plane-invariant closed form, Octave); "
                     "Lair=L, ideal empty-line S21_air=exp(-j k0 L), NN=0."}
