"""Adapter: boivinalex/permittivitycalc (Python, MIT).

Calls ``permittivitycalc.sparam_data.AirlineData`` unmodified.

* ``extract``      -> AirlineData(..., nrw=True): NRW, simultaneous eps & mu.
* ``extract_nni``  -> AirlineData(..., nrw=False): Boughriet NNI, mu := 1.

Interface translation only:
  * AirlineData expects a METAS-VNA-Tools-like array
    [f, |S11|, deg S11, |S21|, deg S21, |S12|, deg S12, |S22|, deg S22] and the
    sample length in **cm**; that is what we build.
  * The package reports ``avg_*`` = mean of the forward (S11,S21) and reverse
    (S22,S12) results; we return those (they are its headline output).
  * tkinter (used only by file-dialog helpers) is absent from the venv's
    Python, so a stub module is inserted before import.
  * Geometry: the package hard-codes the module global LAM_C = inf (coax). For
    waveguide we temporarily set ``sparam_data.LAM_C = 2a`` in cm (the formulas
    in _permittivity_calc are written for a finite LAM_C) and restore it.
  * Sign convention: the package forms eps_r = dielec + 1j*lossfac (see its
    density-normalisation code), and returns ``lossfac = ep_r.imag``; the
    extraction formulas are written with T = exp(-gamma L) (exp(+jwt)), so for
    passive data lossfac < 0 already.  We therefore return
    eps = dielec + 1j*lossfac and mu = mu_real + 1j*mu_imag *without* sign flip.
    (The smoke test confirms Im eps <= 0 for passive input.)
"""
from __future__ import annotations

import contextlib
import io
import sys
import warnings

import numpy as np

from ._common import add_path, as_arrays, cutoff_wavelength, ext, stub_module

_NAME = "permittivitycalc"


def _import():
    if "tkinter" not in sys.modules:
        try:
            import tkinter  # noqa: F401
        except Exception:
            tk = stub_module("tkinter")
            fd = stub_module("tkinter.filedialog",
                             askopenfilename=lambda *a, **k: None,
                             askdirectory=lambda *a, **k: None)
            tk.filedialog = fd
    import matplotlib
    matplotlib.use("Agg", force=False)
    add_path(ext(_NAME))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from permittivitycalc import sparam_data
    return sparam_data


def _run(freq_hz, s11, s21, s12, s22, thickness_m, geometry, nrw):
    sd = _import()
    f, s11, s21, s12, s22 = as_arrays(freq_hz, s11, s21, s12, s22)
    cols = [f]
    for s in (s11, s21, s12, s22):
        cols += [np.abs(s), np.degrees(np.angle(s))]
    data = np.column_stack(cols)
    lam_c_cm = cutoff_wavelength(geometry) * 100.0
    old = sd.LAM_C
    sd.LAM_C = lam_c_cm
    try:
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            a = sd.AirlineData(L=thickness_m * 100.0, airline="custom", dataArray=data,
                               file=None, nrw=nrw)
    finally:
        sd.LAM_C = old
    eps = np.asarray(a.avg_dielec, float) + 1j * np.asarray(a.avg_lossfac, float)
    mu = None
    if nrw:
        mu = np.asarray(a.avg_mu_real, float) + 1j * np.asarray(a.avg_mu_imag, float)
    return eps, mu


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    eps, mu = _run(freq_hz, s11, s21, s12, s22, thickness_m, geometry, nrw=True)
    return {"eps": eps, "mu": mu,
            "notes": "permittivitycalc AirlineData(nrw=True): NRW, np.unwrap(angle T), "
                     "Gamma root test uses complex '<=' (lexicographic, see survey); "
                     "average of forward/reverse. No sign flip applied."}


def extract_nni(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    eps, _ = _run(freq_hz, s11, s21, s12, s22, thickness_m, geometry, nrw=False)
    return {"eps": eps, "mu": None,
            "notes": "permittivitycalc AirlineData(nrw=False): Boughriet NNI (mu=1), "
                     "average of forward/reverse. No sign flip applied."}
