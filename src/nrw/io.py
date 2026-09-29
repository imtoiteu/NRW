"""Input/output: Touchstone .s2p (v1.x) and CSV.

Touchstone 2-port data lines are ordered  f  S11  S21  S12  S22  (note S21 before S12).
Supported option-line formats: RI, MA (magnitude, angle in degrees), DB (dB, degrees);
frequency units Hz/kHz/MHz/GHz; parameter S only. The reference resistance must equal
the holder's characteristic impedance (normally 50 ohm for a 7 mm line); a different
value is reported but not renormalised (the extraction assumes S normalised to the
empty holder).

A minimal own parser is used so that the core has no hard dependency; if scikit-rf is
installed, ``read_touchstone(..., backend="skrf")`` uses it instead (cross-checked in
tests).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_UNITS = {"hz": 1.0, "khz": 1e3, "mhz": 1e6, "ghz": 1e9}


@dataclass
class TwoPortData:
    frequency: np.ndarray
    s11: np.ndarray
    s21: np.ndarray
    s12: np.ndarray
    s22: np.ndarray
    z0: float = 50.0
    source: str = ""
    comments: list = field(default_factory=list)


def _to_complex(a, b, fmt):
    if fmt == "RI":
        return a + 1j * b
    if fmt == "MA":
        return a * np.exp(1j * np.deg2rad(b))
    if fmt == "DB":
        return 10 ** (a / 20) * np.exp(1j * np.deg2rad(b))
    raise ValueError(f"unsupported Touchstone format {fmt}")


def read_touchstone(path, backend="builtin") -> TwoPortData:
    path = Path(path)
    if backend == "skrf":
        import skrf  # optional dependency

        ntw = skrf.Network(str(path))
        if ntw.nports != 2:
            raise ValueError(f"{path}: expected a 2-port network, got {ntw.nports} ports")
        s = ntw.s
        return TwoPortData(ntw.f.copy(), s[:, 0, 0], s[:, 1, 0], s[:, 0, 1], s[:, 1, 1],
                           float(np.real(ntw.z0[0, 0])), str(path))
    if path.suffix.lower() not in (".s2p", ".ts"):
        raise ValueError(f"{path}: expected a .s2p file")
    unit, param, fmt, z0 = "ghz", "S", "MA", 50.0  # Touchstone defaults
    comments, values = [], []
    seen_option = False
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("!"):
            comments.append(line[1:].strip())
            continue
        if "!" in line:
            line = line.split("!", 1)[0].strip()
        if line.startswith("["):
            raise ValueError(f"{path}: Touchstone 2.0 keyword files are not supported by the builtin reader; "
                             "use backend='skrf'")
        if line.startswith("#"):
            if seen_option:
                continue  # Touchstone: only the first option line is significant
            seen_option = True
            toks = line[1:].upper().split()
            i = 0
            while i < len(toks):
                t = toks[i]
                if t.lower() in _UNITS:
                    unit = t.lower()
                elif t in ("S", "Y", "Z", "H", "G"):
                    param = t
                elif t in ("RI", "MA", "DB"):
                    fmt = t
                elif t == "R":
                    z0 = float(toks[i + 1])
                    i += 1
                i += 1
            continue
        values.extend(float(x) for x in line.split())
    if param != "S":
        raise ValueError(f"{path}: only S-parameters are supported (found {param})")
    if len(values) % 9:
        raise ValueError(f"{path}: number of values ({len(values)}) is not a multiple of 9 for a 2-port file")
    arr = np.asarray(values, float).reshape(-1, 9)
    f = arr[:, 0] * _UNITS[unit]
    s = [_to_complex(arr[:, 1 + 2 * k], arr[:, 2 + 2 * k], fmt) for k in range(4)]
    # Touchstone v1 2-port order: S11 S21 S12 S22
    return TwoPortData(f, s[0], s[1], s[2], s[3], z0, str(path), comments)


def write_touchstone(path, f, s11, s21, s12, s22, z0=50.0, comments=()):
    """Write a Touchstone v1 .s2p file in RI format, frequency in Hz."""
    path = Path(path)
    lines = [f"! {c}" for c in comments]
    lines.append(f"# HZ S RI R {z0:g}")
    for k in range(len(f)):
        vals = [f[k]]
        for s in (s11, s21, s12, s22):
            vals += [s[k].real, s[k].imag]
        lines.append(" ".join(f"{v:.12e}" for v in vals))
    path.write_text("\n".join(lines) + "\n")


def read_csv_sparams(path) -> TwoPortData:
    """CSV with header: f_hz,s11_re,s11_im,s21_re,s21_im[,s12_re,s12_im,s22_re,s22_im]."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    names = data.dtype.names
    need = ["f_hz", "s11_re", "s11_im", "s21_re", "s21_im"]
    missing = [n for n in need if n not in names]
    if missing:
        raise ValueError(f"{path}: missing CSV columns {missing}")
    c = lambda k: data[f"{k}_re"] + 1j * data[f"{k}_im"]  # noqa: E731
    s11, s21 = c("s11"), c("s21")
    s12 = c("s12") if "s12_re" in names else s21.copy()
    s22 = c("s22") if "s22_re" in names else s11.copy()
    return TwoPortData(np.asarray(data["f_hz"], float), s11, s21, s12, s22, source=str(path))


def read_measurement(path, backend="builtin") -> TwoPortData:
    p = Path(path)
    if p.suffix.lower() == ".csv":
        return read_csv_sparams(p)
    return read_touchstone(p, backend=backend)


def _jsonable(x):
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, float):
        return None if not np.isfinite(x) else x
    if isinstance(x, complex):
        return {"re": _jsonable(x.real), "im": _jsonable(x.imag)}
    if isinstance(x, np.ndarray):
        return [_jsonable(v) for v in x.tolist()] if x.dtype != object else [_jsonable(v) for v in x]
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    return x


def dump_json(obj, path):
    Path(path).write_text(json.dumps(_jsonable(obj), indent=2))
