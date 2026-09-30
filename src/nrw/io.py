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

import csv
import io as _io
import json
import re
import tempfile
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
    # which parameters were actually present in the file; absent S12/S22 are filled with
    # S21/S11 copies for backwards compatibility but must not be treated as measured
    measured: tuple = ("s11", "s21", "s12", "s22")
    notes: list = field(default_factory=list)  # parser decisions the user should know about

    def s12_measured(self):
        return self.s12 if "s12" in self.measured else None

    def s22_measured(self):
        return self.s22 if "s22" in self.measured else None


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


_FREQ_UNIT_RE = re.compile(r"(ghz|mhz|khz|hz)$")


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _classify_s_column(key: str, p: str):
    """Classify a normalised column name for S-parameter ``p`` ('11', '21', ...).

    Returns one of 're', 'im', 'db', 'mag', 'phase_deg', 'phase_rad' or None.
    """
    tag = f"s{p}"
    if tag not in key:
        return None
    # guard against e.g. 's111' or 's211' style collisions
    i = key.index(tag)
    if i + len(tag) < len(key) and key[i + len(tag)].isdigit():
        return None
    rest = key[:i] + key[i + len(tag):]
    if rest in ("re", "real", "r") or rest.startswith("real") or rest.endswith("real"):
        return "re"
    if rest in ("im", "imag", "i", "imaginary") or rest.startswith("imag") or rest.endswith("imag"):
        return "im"
    if any(t in rest for t in ("phase", "ang", "arg", "phi")) or rest in ("deg", "rad"):
        return "phase_rad" if "rad" in rest else "phase_deg"
    if "db" in rest:
        return "db"
    if rest in ("mag", "abs", "lin", "magnitude", "maglin", "linmag", "linear", "magnitudelinear", "maglinear"):
        return "mag"
    return None


def parse_csv_text(text: str, source: str = "", freq_unit: str | None = None) -> TwoPortData:
    """Parse a CSV/TSV S-parameter table with a header row.

    Recognised columns (case, spaces, brackets and underscores ignored):

    * frequency: ``f``, ``freq``, ``frequency`` with optional unit suffix (Hz, kHz, MHz, GHz),
      e.g. ``f_hz``, ``Frequency (GHz)``. Without a unit in the header, ``freq_unit`` must be
      given (no guessing).
    * per S-parameter (S11, S21, S12, S22), either
      real/imag: ``S11_real``/``S11_imag``, ``s11_re``/``s11_im``, ``Re(S11)``/``Im(S11)``; or
      magnitude + phase: ``S11_mag`` (linear) or ``S11_dB``, with ``S11_phase_deg`` /
      ``S11_phase_rad`` / ``S11_phase`` (degrees assumed, recorded in ``notes``).

    S11 and S21 are required. Missing S12/S22 are recorded as not measured.
    Comment lines starting with '#' or '!' are skipped; the delimiter (comma, semicolon,
    tab, whitespace) is detected.
    """
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith(("#", "!"))]
    if len(lines) < 2:
        raise ValueError(f"{source}: CSV needs a header row and at least one data row")
    head = lines[0]
    if "," in head:
        delim = ","
    elif ";" in head:
        delim = ";"
    elif "\t" in head:
        delim = "\t"
    else:
        delim = None  # whitespace
    rows = list(csv.reader(lines, delimiter=delim)) if delim else [ln.split() for ln in lines]
    header = [h.strip() for h in rows[0]]
    keys = [_norm(h) for h in header]
    notes = []
    # frequency column
    fcol = None
    for j, k in enumerate(keys):
        if re.fullmatch(r"(f|freq|frequency)(hz|khz|mhz|ghz)?", k):
            fcol = j
            break
    if fcol is None:
        raise ValueError(f"{source}: no frequency column found. Expected a header such as 'f_hz', "
                         f"'freq_GHz' or 'Frequency (GHz)'; got {header}")
    m = _FREQ_UNIT_RE.search(keys[fcol])
    unit = m.group(1) if m and keys[fcol] not in ("f",) else None
    if unit is None:
        if freq_unit is None:
            raise ValueError(f"{source}: the frequency column '{header[fcol]}' has no unit. Rename it (e.g. "
                             "'freq_GHz') or specify the frequency unit explicitly.")
        unit = freq_unit.lower()
        notes.append(f"frequency unit '{unit}' taken from user setting (not in header)")
    if unit not in _UNITS:
        raise ValueError(f"{source}: unknown frequency unit {unit!r}")
    try:
        data = np.array([[float(x) for x in r] for r in rows[1:]], dtype=float)
    except ValueError as e:
        raise ValueError(f"{source}: non-numeric value in data rows ({e})") from e
    if data.shape[1] != len(header):
        raise ValueError(f"{source}: rows have {data.shape[1]} columns but the header has {len(header)}")
    f = data[:, fcol] * _UNITS[unit]
    out = {}
    for p in ("11", "21", "12", "22"):
        cols = {}
        for j, k in enumerate(keys):
            c = _classify_s_column(k, p)
            if c and c not in cols:
                cols[c] = j
        if "re" in cols and "im" in cols:
            out[p] = data[:, cols["re"]] + 1j * data[:, cols["im"]]
        elif ("db" in cols or "mag" in cols) and ("phase_deg" in cols or "phase_rad" in cols):
            mag = 10 ** (data[:, cols["db"]] / 20) if "db" in cols else data[:, cols["mag"]]
            if "phase_rad" in cols:
                ph = data[:, cols["phase_rad"]]
            else:
                ph = np.deg2rad(data[:, cols["phase_deg"]])
                if not re.search(r"deg", keys[cols["phase_deg"]]):
                    notes.append(f"S{p} phase column '{header[cols['phase_deg']]}' has no unit: degrees assumed")
            out[p] = mag * np.exp(1j * ph)
        elif cols:
            raise ValueError(f"{source}: S{p} columns found ({[header[j] for j in cols.values()]}) but not a "
                             "complete pair (real+imag, or magnitude/dB + phase)")
    missing = [f"S{p}" for p in ("11", "21") if p not in out]
    if missing:
        raise ValueError(f"{source}: required columns for {missing} not found. Expected e.g. 'S11_real, S11_imag' "
                         f"or 'S11_dB, S11_phase_deg'; header was {header}")
    measured = tuple(f"s{p}" for p in ("11", "21", "12", "22") if p in out)
    s12 = out.get("12", out["21"].copy())
    s22 = out.get("22", out["11"].copy())
    if "12" not in out or "22" not in out:
        notes.append("S12/S22 not present: a reciprocal, symmetric sample is assumed where they are needed")
    return TwoPortData(f, out["11"], out["21"], s12, s22, source=source, measured=measured, notes=notes)


def read_csv_sparams(path, freq_unit=None) -> TwoPortData:
    """Read a CSV S-parameter table (see ``parse_csv_text`` for accepted headers)."""
    return parse_csv_text(Path(path).read_text(), source=str(path), freq_unit=freq_unit)


def load_measurement_bytes(name: str, content: bytes, freq_unit: str | None = None) -> TwoPortData:
    """Parse an uploaded file (.s2p or .csv/.txt/.tsv) from bytes.

    Touchstone files are parsed with scikit-rf when installed and cross-checked against the
    builtin v1 parser: both must agree, and structural errors detected by the builtin parser
    (e.g. a value count that is not a multiple of 9) are fatal, because scikit-rf silently
    pads some malformed files. scikit-rf alone is used only for Touchstone 2.0 keyword files. The bytes are written to a private temporary
    directory only for the duration of parsing.
    """
    suffix = Path(name).suffix.lower()
    if suffix in (".csv", ".txt", ".tsv", ".dat"):
        try:
            text = content.decode("utf-8-sig")
        except UnicodeDecodeError:
            text = content.decode("latin-1")
        return parse_csv_text(text, source=name, freq_unit=freq_unit)
    if suffix != ".s2p":
        raise ValueError(f"{name}: unsupported file type {suffix!r}; use .s2p (Touchstone) or .csv")
    with tempfile.TemporaryDirectory(prefix="nrw_upload_") as td:
        p = Path(td) / "upload.s2p"
        p.write_bytes(content)
        builtin, err_b = None, None
        try:
            builtin = read_touchstone(p)
        except Exception as e:  # noqa: BLE001
            err_b = e
        try:
            import skrf  # noqa: F401
            have_skrf = True
        except ImportError:
            have_skrf = False
        if have_skrf:
            try:
                d = read_touchstone(p, backend="skrf")
            except Exception as e:  # noqa: BLE001
                if builtin is None:
                    raise ValueError(f"{name}: could not parse Touchstone file (scikit-rf: {e}; builtin: {err_b})") from e
                builtin.notes.append(f"scikit-rf could not parse the file ({e}); builtin parser used")
                builtin.source = name
                return builtin
            d.source = name
            if builtin is not None:
                same = (builtin.frequency.shape == d.frequency.shape and np.allclose(builtin.frequency, d.frequency)
                        and all(np.allclose(getattr(builtin, k), getattr(d, k), rtol=1e-9, atol=1e-12)
                                for k in ("s11", "s21", "s12", "s22")))
                if not same:
                    raise ValueError(f"{name}: scikit-rf and the builtin Touchstone reader disagree on this file; "
                                     "check the option line (# unit S format R z0) and the data layout")
                d.comments = builtin.comments
                d.notes.append("parsed with scikit-rf; builtin parser agrees")
                return d
            # scikit-rf silently pads/broadcasts some malformed files (e.g. a truncated data line),
            # so it is trusted alone only for Touchstone 2.0 keyword files, which the builtin
            # reader does not support. Structural errors found by the builtin reader are fatal.
            if "Touchstone 2.0" not in str(err_b):
                raise ValueError(f"{name}: malformed Touchstone file: {err_b}")
            d.notes.append("Touchstone 2.0 file parsed with scikit-rf only (builtin reader supports v1)")
            return d
        if builtin is None:
            raise ValueError(f"{name}: could not parse Touchstone file: {err_b}")
        builtin.source = name
        builtin.notes.append("parsed with builtin Touchstone reader (scikit-rf not installed)")
        return builtin


def describe(data: TwoPortData) -> dict:
    """Summary used by the data preview (all values plain Python types)."""
    f = np.asarray(data.frequency, float)
    df = np.diff(f)
    out = {
        "source": data.source,
        "n_points": int(f.size),
        "f_min_hz": float(f.min()) if f.size else float("nan"),
        "f_max_hz": float(f.max()) if f.size else float("nan"),
        "strictly_increasing": bool(np.all(df > 0)) if f.size > 1 else False,
        "df_min_hz": float(df.min()) if df.size else float("nan"),
        "df_max_hz": float(df.max()) if df.size else float("nan"),
        "df_median_hz": float(np.median(df)) if df.size else float("nan"),
        "uniform_grid": bool(df.size and np.max(np.abs(df - np.median(df))) <= 1e-6 * abs(np.median(df))),
        "measured": list(data.measured),
        "z0_ohm": float(data.z0),
        "notes": list(data.notes),
        "finite": bool(all(np.all(np.isfinite(getattr(data, k))) for k in ("s11", "s21", "s12", "s22"))
                       and np.all(np.isfinite(f))),
    }
    for k in ("s11", "s21", "s12", "s22"):
        a = np.abs(getattr(data, k))
        out[f"{k}_abs_min"] = float(np.nanmin(a)) if a.size else float("nan")
        out[f"{k}_abs_max"] = float(np.nanmax(a)) if a.size else float("nan")
    return out


def read_measurement(path, backend="builtin", freq_unit=None) -> TwoPortData:
    p = Path(path)
    if p.suffix.lower() in (".csv", ".txt", ".tsv", ".dat"):
        return read_csv_sparams(p, freq_unit=freq_unit)
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
