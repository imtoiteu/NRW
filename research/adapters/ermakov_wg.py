"""Adapter: sergeyermakov239/WaveGuide_Measurements (MATLAB live script, NO
LICENSE -> all rights reserved; local evaluation only, nothing redistributed)
under Octave.

The extraction function ``read_data_phase_corr(filename, d, num_2pi)`` lives
inside ``right_measurements.mlx`` (a zip). At run time the adapter unzips the
pinned file, pulls the function text out of matlab/document.xml and applies
these exact edits (checked to match exactly once each):
  1. signature -> ``read_data_phase_corr_w(freq, S11, S21, d, num_2pi, a_in)``;
  2. the five lines reading the .s2p via RF Toolbox (``read(rfdata.data,..)``,
     ``freq=data.Freq``, ``S_mat=...``, ``S11=reshape(..)``, ``S21=reshape(..)``)
     are deleted (the data arrive as arguments);
  3. ``a = 23.1*10^(-3);`` -> ``a = a_in;`` (a = Inf for coax -> f0 = 0);
  4. ``L=270*10^(-3);`` -> ``L=d;`` and ``L2=130*10^(-3);`` -> ``L2=0;`` so that
     L1 = L - d - L2 = 0: the author's fixture de-embedding becomes identity
     because our planes are at the sample faces.
Everything else (c = 2.998e8, cumulative phase tracking, eps formula) is kept.
num_2pi = 0. eps only (mu = 1 assumed). Output eps has Im < 0 for passive data
(eps = f0^2/f^2 - (c gamma/(2 pi f))^2, gamma = ln(1/T)/d); returned unchanged.
"""
from __future__ import annotations

import html
import re
import tempfile
import zipfile
from pathlib import Path

import numpy as np

from ._common import as_arrays, ext, run_octave

_NAME = "WaveGuide_Measurements"


def _mlx_text() -> str:
    with zipfile.ZipFile(ext(_NAME) / "right_measurements.mlx") as z:
        s = z.read("matlab/document.xml").decode("utf8")
    paras = re.findall(r"<w:p>(.*?)</w:p>", s, flags=re.S)
    out = []
    for p in paras:
        t = "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", p, flags=re.S))
        t = html.unescape(t)
        t = t.replace("<![CDATA[", "").replace("]]>", "")
        out.append(t)
    return "\n".join(out)


def wrapped_source() -> str:
    txt = _mlx_text()
    i = txt.index("function [eps,freq,S11,S21] = read_data_phase_corr(")
    j = txt.index("function [freq,S11,S12,S21,S22] = get_s(")
    fn = txt[i:j]
    subs = [
        ("function [eps,freq,S11,S21] = read_data_phase_corr(filename,d,num_2pi)",
         "function [eps,freq,S11,S21] = read_data_phase_corr_w(freq,S11,S21,d,num_2pi,a_in)"),
        ("data=read(rfdata.data,filename);", ""),
        ("freq=data.Freq;", ""),
        ("S_mat=data.S_Parameters;", ""),
        ("S11=reshape(S_mat(1,1,:),[],1);", ""),
        ("S21=reshape(S_mat(2,1,:),[],1);", ""),
        ("a = 23.1*10^(-3);", "a = a_in;"),
        ("L=270*10^(-3);", "L=d;"),
        ("L2=130*10^(-3);", "L2=0;"),
    ]
    for old, new in subs:
        n = fn.count(old)
        if n != 1:
            raise RuntimeError(f"ermakov wrapper: expected exactly one {old!r}, found {n}")
        fn = fn.replace(old, new)
    return "% AUTO-WRAPPED at run time from right_measurements.mlx (see ermakov_wg.py)\n" + fn


def extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> dict:
    f, s11, s21 = as_arrays(freq_hz, s11, s21)
    a = np.inf if geometry["kind"] == "coax" else float(geometry["a_m"])
    with tempfile.TemporaryDirectory(prefix="ermakov_") as td:
        Path(td, "read_data_phase_corr_w.m").write_text(wrapped_source())
        code = """
[eps, fr, S11o, S21o] = read_data_phase_corr_w(f(:), s11(:), s21(:), L, 0, a);
save('-v7', outfile, 'eps');
"""
        r = run_octave(code, {"f": f, "s11": s11, "s21": s21, "a": a, "L": thickness_m},
                       [Path(td)])
    return {"eps": np.atleast_1d(r["eps"]).astype(complex), "mu": None,
            "notes": "ermakov read_data_phase_corr (Octave, wrapped from .mlx), mu=1, "
                     "cumulative phase tracking, c=2.998e8, num_2pi=0."}
