"""Adapter: njchorda/Nicolson-Ross-Weir-Materal-Parameter-Extraction (MATLAB
script, MIT) under Octave.

``materialParameterExtraction.m`` is a *script* (clear/clc, reads a Touchstone
file with a third-party SPARAMS reader, plots). At run time the adapter reads
the pinned file and extracts, verbatim, the lines between the markers
``%% NRW Calculations`` and ``%% Plot results``; these are wrapped in a
temporary function file. Exact edits (nothing inside the calculation block is
changed):
  1. added header ``function [era, ura] = njchorda_nrw(S, a, L, c0)``, i.e. the
     ``%% File Input`` section (SPARAMS(filepath)) and the ``%% Dimensions``
     section (a = 22.86e-3; L = 6.3e-3) are replaced by function arguments;
  2. ``c0 = physconst('lightspeed')`` (a MATLAB toolbox function) is replaced by
     passing c0 = 299792458 (the value physconst returns);
  3. the plotting / printing sections after ``%% Plot results`` are dropped.
The calculation block hard-codes ``lamc = 2.*a``; for coax we pass a = Inf so
lamc = Inf (1/lamc^2 = 0).
Sign: era, ura are complex with Im < 0 for passive exp(+jwt) data (the script
itself plots |imag|); returned unchanged.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from ._common import C0, as_arrays, ext, run_octave

_NAME = "Nicolson-Ross-Weir-Materal-Parameter-Extraction"
_START, _END = "%% NRW Calculations", "%% Plot results"


def wrapped_source() -> str:
    src = (ext(_NAME) / "materialParameterExtraction.m").read_text()
    i, j = src.index(_START), src.index(_END)
    block = src[i:j]
    return ("function [era, ura] = njchorda_nrw(S, a, L, c0)\n"
            "% AUTO-WRAPPED at run time from njchorda materialParameterExtraction.m (MIT).\n"
            + block + "\nend\n")


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    a = np.inf if geometry["kind"] == "coax" else float(geometry["a_m"])
    with tempfile.TemporaryDirectory(prefix="njchorda_") as td:
        Path(td, "njchorda_nrw.m").write_text(wrapped_source())
        code = """
S = struct('f', f(:), 'S11', s11(:), 'S21', s21(:));
[era, ura] = njchorda_nrw(S, a, L, c0);
save('-v7', outfile, 'era', 'ura');
"""
        r = run_octave(code, {"f": f, "s11": s11, "s21": s21, "a": a, "L": thickness_m,
                              "c0": C0}, [Path(td)])
    return {"eps": np.atleast_1d(r["era"]).astype(complex),
            "mu": np.atleast_1d(r["ura"]).astype(complex),
            "notes": "njchorda NRW script block (Octave), principal branch m=0, no unwrap."}
