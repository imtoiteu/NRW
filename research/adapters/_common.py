"""Shared helpers for the thin adapters around third-party NRW implementations.

Nothing in here implements an extraction algorithm; it only moves data in and
out of the external code (paths, Octave subprocess, geometry bookkeeping).

Conventions (docs/conventions.md): exp(+jwt), eps = eps' - j eps'' (passive =>
Im eps <= 0), S-parameters referenced to the sample faces, SI units.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import numpy as np

C0 = 299_792_458.0

ROOT = Path(__file__).resolve().parents[2]
EXTERNAL = ROOT / "research" / "external"
OCTAVE_DIR = Path(__file__).resolve().parent / "octave"


class ExternalMissing(ImportError):
    """Raised when the pinned external repository has not been fetched."""


def ext(name: str) -> Path:
    p = EXTERNAL / name
    if not p.is_dir():
        raise ExternalMissing(
            f"{p} not found - run scripts/fetch_external.sh to clone the pinned commit"
        )
    return p


def add_path(p: Path) -> None:
    s = str(p)
    if s not in sys.path:
        sys.path.insert(0, s)


def cutoff_wavelength(geometry: dict) -> float:
    """lambda_c in metres: inf for coax/TEM, 2a for rectangular TE10."""
    kind = geometry.get("kind")
    if kind == "coax":
        return float("inf")
    if kind == "waveguide":
        return 2.0 * float(geometry["a_m"])
    raise ValueError(f"unknown geometry kind {kind!r}")


def require_tem(geometry: dict, who: str) -> None:
    if geometry.get("kind") != "coax":
        raise NotImplementedError(f"{who}: external code supports only TEM (coax/free space)")


def as_arrays(freq_hz, *s):
    f = np.asarray(freq_hz, dtype=float).ravel()
    out = [np.asarray(x, dtype=complex).ravel() for x in s]
    for x in out:
        if x.shape != f.shape:
            raise ValueError("S-parameter arrays must have the same length as freq_hz")
    return (f, *out)


def stub_module(name: str, **attrs) -> types.ModuleType:
    """Insert a permissive stub module (used only for GUI/simulator imports that
    the extraction routine itself never calls)."""

    class _Any:
        def __init__(self, *a, **k):
            pass

        def __call__(self, *a, **k):
            return _Any()

        def __getattr__(self, item):
            return _Any()

    class _StubModule(types.ModuleType):
        def __getattr__(self, item):  # module-level fallback
            return _Any

    m = _StubModule(name)
    for k, v in attrs.items():
        setattr(m, k, v)
    sys.modules[name] = m
    parent, _, child = name.rpartition(".")
    if parent and parent in sys.modules:
        setattr(sys.modules[parent], child, m)
    return m


class temporary_stubs:
    """Context manager: stub the given module names only for the duration of an
    import, then remove them again from sys.modules so that other libraries
    (e.g. torch, which probes for IPython) are not affected."""

    def __init__(self, *names):
        self.names = names
        self.added = []

    def __enter__(self):
        for n in self.names:
            if n in sys.modules:
                continue
            try:
                __import__(n)
                continue
            except Exception:
                pass
            stub_module(n)
            self.added.append(n)
        for n in self.added:  # wire children onto parents
            parent, _, child = n.rpartition(".")
            if parent in sys.modules:
                setattr(sys.modules[parent], child, sys.modules[n])
        return self

    def __exit__(self, *exc):
        for n in self.added:
            sys.modules.pop(n, None)
        return False


# --------------------------------------------------------------------------
# Octave bridge
# --------------------------------------------------------------------------

def octave_available() -> bool:
    return shutil.which("octave") is not None


def run_octave(code: str, inputs: dict, addpaths: list[Path], timeout: int = 300) -> dict:
    """Run `code` in octave-cli. `inputs` are saved to in.mat and loaded as
    variables; the code must `save('-v7', outfile, ...)` its results into the
    variable `outfile`. Returns the loaded dict (squeezed arrays)."""
    from scipy.io import loadmat, savemat

    if not octave_available():
        raise RuntimeError("octave not installed")
    with tempfile.TemporaryDirectory(prefix="nrw_oct_") as td:
        td = Path(td)
        infile, outfile = td / "in.mat", td / "out.mat"
        savemat(infile, {k: (np.asarray(v).reshape(-1, 1) if np.ndim(v) == 1 else v)
                         for k, v in inputs.items()})
        paths = "".join(f"addpath('{p}');" for p in [OCTAVE_DIR / "shims", *addpaths])
        script = (
            "warning('off','all');" + paths +
            f"load('{infile}'); outfile='{outfile}';\n" + code + "\n"
        )
        (td / "nrw_adapter_job.m").write_text(script)
        proc = subprocess.run(
            ["octave-cli", "--no-gui", "--quiet", "--no-window-system", str(td / "nrw_adapter_job.m")],
            cwd=td, capture_output=True, text=True, timeout=timeout,
        )
        if proc.returncode != 0 or not outfile.exists():
            raise RuntimeError(
                f"octave failed (rc={proc.returncode}):\nSTDOUT:\n{proc.stdout[-2000:]}\n"
                f"STDERR:\n{proc.stderr[-2000:]}"
            )
        d = loadmat(outfile)
        return {k: np.squeeze(v) for k, v in d.items() if not k.startswith("__")}
