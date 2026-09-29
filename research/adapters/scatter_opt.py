"""Adapter: jdhuang-csm/scatter-opt (MATLAB, BSD-3-Clause) run under GNU Octave.

Calls ``scatter_opt.nrw_extract(Sdata, L, L_air, lambda_c, ...)`` unmodified.

Interface translation (no edits to the external .m files):
  * ``Sdata`` is passed as a struct with column fields freq, s11, s21 (the
    MATLAB code expects a table; nrw_extract only uses dot-indexing on it).
  * Octave has no ``table`` class, so a 15-line shim ``array2table.m``
    (research/adapters/octave/shims) returns a struct with the same field
    names; the shim is placed on the Octave path, not in the external tree.
  * Default branch selection in the repo is 'rationalfit' (MATLAB RF Toolbox
    ``rationalfit`` - unavailable in Octave). ``extract`` uses
    BranchSelectMethod='manual', BaseBranch=0 (+ its own phase unwrapping of
    angle(1/T) inside nrw_lambda); ``extract_groupdelay`` uses its
    'GroupDelay' method (pure MATLAB) over branches -2..2.
  * lambda_c = Inf for coax (its formulas contain 1/lambda_c^2 -> 0) or 2a.
Sign: T = exp(-gamma L) algebra; eps = eps' - j eps'' for passive input.
Returned unchanged (smoke test checks sign).
"""
from __future__ import annotations

import numpy as np

from ._common import as_arrays, cutoff_wavelength, ext, run_octave

_NAME = "scatter-opt"


def _run(freq_hz, s11, s21, thickness_m, geometry, method_args: str):
    root = ext(_NAME)
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    lam_c = cutoff_wavelength(geometry)
    code = f"""
Sdata = struct('freq', f(:), 's11', s11(:), 's21', s21(:));
out = scatter_opt.nrw_extract(Sdata, L, 0, lam_c, {method_args});
eps = out.eps; mu = out.mu;
save('-v7', outfile, 'eps', 'mu');
"""
    r = run_octave(code, {"f": f, "s11": s11, "s21": s21, "L": thickness_m,
                          "lam_c": lam_c if np.isfinite(lam_c) else np.inf}, [root])
    return np.atleast_1d(r["eps"]).astype(complex), np.atleast_1d(r["mu"]).astype(complex)


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    eps, mu = _run(freq_hz, s11, s21, thickness_m, geometry,
                   "'BranchSelectMethod','manual','BaseBranch',0")
    return {"eps": eps, "mu": mu,
            "notes": "scatter-opt nrw_extract (Octave), manual BaseBranch=0 + unwrap(angle(1/T))."}


def extract_groupdelay(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    eps, mu = _run(freq_hz, s11, s21, thickness_m, geometry,
                   "'BranchSelectMethod','GroupDelay','BranchLimits',[-2 2]")
    return {"eps": eps, "mu": mu,
            "notes": "scatter-opt nrw_extract (Octave), GroupDelay base-branch selection "
                     "(lowest quarter of band) + unwrap."}
