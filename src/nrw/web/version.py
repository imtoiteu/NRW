"""Software identification shown in the UI and written into every report."""
from __future__ import annotations

import functools
import subprocess
from pathlib import Path

from .. import __version__

REPO = Path(__file__).resolve().parents[3]


@functools.lru_cache(maxsize=1)
def software_id() -> dict:
    """Package version and Git commit (with a 'dirty' marker if the working tree has changes)."""
    out = {"package_version": __version__, "commit": "unknown", "commit_short": "unknown", "dirty": None}
    try:
        c = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True,
                           timeout=5)
        if c.returncode == 0:
            out["commit"] = c.stdout.strip()
            out["commit_short"] = out["commit"][:10]
            d = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, timeout=5)
            out["dirty"] = bool(d.stdout.strip()) if d.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        pass
    return out


def software_label() -> str:
    s = software_id()
    dirty = " (uncommitted changes)" if s["dirty"] else ""
    return f"nrw {s['package_version']} · commit {s['commit_short']}{dirty}"
