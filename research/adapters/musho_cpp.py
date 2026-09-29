"""Adapter: Dr-Musho-Research-Group/NRW_Permittivity_Extraction (C++17, NO
LICENSE -> all rights reserved; local evaluation only, code not redistributed).

The unmodified ``src/nicolson_ross_weir_permittivity.cpp`` is compiled with
``g++ -O2 -std=c++17`` into research/external/.build/ (gitignored) and run as a
subprocess on a CSV (Frequency,S11_real,S11_imag,S21_real,S21_imag). We read
the ``epsilon_r_real, epsilon_r_imag`` columns (the raw complex result, not the
sign-flipped "*_plot" columns).

Method: NRW eps-only with a user-supplied mu (we pass mu = 1 + 0j), principal
log (no unwrap, no branch index), and -- by design of that code -- a finite
cutoff wavelength lambda_c = pi (a + b)/2 derived from the coax radii even
for TEM. The geometry dict may carry "inner_radius_m"/"outer_radius_m"; if not,
APC-7 radii (1.52 mm / 3.50 mm) are used, as noted in `notes`. Coax only.
Sign: T = exp(-gamma L) algebra, eps = eps' - j eps''; returned unchanged.
"""
from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from ._common import EXTERNAL, as_arrays, ext, require_tem

_NAME = "NRW_Permittivity_Extraction"


def _binary() -> Path:
    src = ext(_NAME) / "src" / "nicolson_ross_weir_permittivity.cpp"
    out = EXTERNAL / ".build" / "nicolson_ross_weir_permittivity"
    if not out.exists() or out.stat().st_mtime < src.stat().st_mtime:
        out.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["g++", "-O2", "-std=c++17", "-o", str(out), str(src)],
                       check=True, capture_output=True, text=True)
    return out


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    require_tem(geometry, "musho_cpp")
    exe = _binary()
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    ri = float(geometry.get("inner_radius_m", 1.52e-3))
    ro = float(geometry.get("outer_radius_m", 3.50e-3))
    with tempfile.TemporaryDirectory(prefix="nrw_cpp_") as td:
        td = Path(td)
        with open(td / "in.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Frequency", "S11_real", "S11_imag", "S21_real", "S21_imag"])
            for row in zip(f, s11.real, s11.imag, s21.real, s21.imag):
                w.writerow([repr(float(x)) for x in row])
        proc = subprocess.run(
            [str(exe), "--input", "in.csv", "--output", "out.csv",
             "--sample-length-m", repr(thickness_m), "--mu-r-real", "1.0", "--mu-r-imag", "0.0",
             "--inner-radius-m", repr(ri), "--outer-radius-m", repr(ro)],
            cwd=td, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr[-1500:])
        with open(td / "out.csv") as fh:
            rows = list(csv.DictReader(fh))
    eps = np.array([float(r["epsilon_r_real"]) + 1j * float(r["epsilon_r_imag"]) for r in rows])
    return {"eps": eps, "mu": None,
            "notes": f"Musho C++ NRW eps-only, mu=1 given, principal log, lambda_c=pi(a+b)/2 "
                     f"with a={ri} m, b={ro} m."}
