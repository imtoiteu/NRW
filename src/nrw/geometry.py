"""Measurement-cell geometries and guided-wave quantities.

Conventions (docs/conventions.md): exp(+jwt) time dependence, fields ~ exp(-gamma z),
gamma = sqrt(kc^2 - k0^2 eps mu) with Re(gamma) >= 0 for passive media.

Supported cells
---------------
* ``Coax``: coaxial air line, TEM mode, kc = 0. A 7 mm (APC-7) line covers 2-18 GHz
  without higher-order modes, which is why it is the default cell in this project.
* ``RectangularWaveguide``: TE10 mode, kc = pi / a (sample fills the cross-section).

Only the cutoff wavenumber enters the extraction equations; the coaxial diameters are
kept for the optional air-gap model (``nrw.forward.coax_air_gap_effective``) and for
higher-order-mode checks.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

C0 = 299_792_458.0  # speed of light in vacuum, m/s (exact, SI)


def principal_sqrt_decaying(x: np.ndarray) -> np.ndarray:
    """Square root choosing the branch with Re >= 0 (decaying wave for exp(-gamma z)).

    On the branch cut (Re = 0 exactly, lossless propagation) the root with Im >= 0 is
    returned, i.e. gamma = +j beta with beta > 0 (forward-propagating wave).
    """
    r = np.sqrt(np.asarray(x, dtype=complex))
    flip = (r.real < 0) | ((r.real == 0) & (r.imag < 0))
    return np.where(flip, -r, r)


@dataclass(frozen=True)
class Geometry:
    """Base class. ``kc`` is the cutoff wavenumber of the operating mode (rad/m)."""

    kind: str = field(init=False, default="generic")

    @property
    def kc(self) -> float:  # pragma: no cover - abstract
        raise NotImplementedError

    def k0(self, f: np.ndarray) -> np.ndarray:
        return 2 * np.pi * np.asarray(f, dtype=float) / C0

    def gamma0(self, f: np.ndarray) -> np.ndarray:
        """Propagation constant of the empty (air/vacuum filled) cell."""
        k0 = self.k0(f)
        return principal_sqrt_decaying(self.kc**2 - k0**2 + 0j)

    def gamma(self, f: np.ndarray, eps: np.ndarray, mu: np.ndarray) -> np.ndarray:
        """Propagation constant of the filled cell (forward-model root, Re >= 0)."""
        k0 = self.k0(f)
        return principal_sqrt_decaying(self.kc**2 - k0**2 * np.asarray(eps) * np.asarray(mu))

    def z_norm(self, f, eps, mu, gamma=None) -> np.ndarray:
        """Wave impedance of the filled cell normalised to the empty cell: z = mu*gamma0/gamma."""
        g = self.gamma(f, eps, mu) if gamma is None else gamma
        return np.asarray(mu) * self.gamma0(f) / g

    def eps_from_gamma_mu(self, f, gamma, mu) -> np.ndarray:
        """Invert gamma^2 = kc^2 - k0^2 eps mu for eps."""
        k0 = self.k0(f)
        return (self.kc**2 - np.asarray(gamma) ** 2) / (k0**2 * np.asarray(mu))

    def cutoff_frequency(self) -> float:
        return self.kc * C0 / (2 * np.pi)

    def to_dict(self) -> dict:
        raise NotImplementedError  # pragma: no cover


@dataclass(frozen=True)
class Coax(Geometry):
    """Coaxial air line (TEM). Default diameters are the 7 mm (APC-7) line."""

    inner_diameter_m: float = 3.04e-3
    outer_diameter_m: float = 7.00e-3
    kind: str = field(init=False, default="coax")

    @property
    def kc(self) -> float:
        return 0.0

    def te11_cutoff_hz(self, eps_mu: float = 1.0) -> float:
        """Approximate TE11 cutoff, fc ~ c / (pi (a+b)/2 sqrt(eps mu)) (a, b radii).

        Above this frequency a sample with the given eps*mu can support the first
        higher-order mode; used only as a warning.
        """
        a = self.inner_diameter_m / 2
        b = self.outer_diameter_m / 2
        return C0 / (np.pi * (a + b) * np.sqrt(np.real(eps_mu)))

    def to_dict(self) -> dict:
        return {
            "kind": "coax",
            "inner_diameter_m": self.inner_diameter_m,
            "outer_diameter_m": self.outer_diameter_m,
        }


@dataclass(frozen=True)
class RectangularWaveguide(Geometry):
    """Rectangular waveguide, TE10 mode. ``a_m`` is the broad-wall width."""

    a_m: float = 22.86e-3  # WR-90
    b_m: float = 10.16e-3
    kind: str = field(init=False, default="waveguide")

    @property
    def kc(self) -> float:
        return np.pi / self.a_m

    def to_dict(self) -> dict:
        return {"kind": "waveguide", "a_m": self.a_m, "b_m": self.b_m}


def geometry_from_dict(d: dict | Geometry | None) -> Geometry:
    """Build a geometry from a dict such as {"kind": "coax"} or {"kind": "waveguide", "a_m": ...}."""
    if d is None:
        return Coax()
    if isinstance(d, Geometry):
        return d
    kind = str(d.get("kind", "coax")).lower()
    if kind in ("coax", "coaxial", "tem", "freespace", "free_space"):
        return Coax(
            inner_diameter_m=float(d.get("inner_diameter_m", 3.04e-3)),
            outer_diameter_m=float(d.get("outer_diameter_m", 7.00e-3)),
        )
    if kind in ("waveguide", "rectangular", "wr", "te10"):
        if "a_m" not in d:
            raise ValueError("rectangular waveguide requires the broad-wall width 'a_m'")
        return RectangularWaveguide(a_m=float(d["a_m"]), b_m=float(d.get("b_m", d["a_m"] / 2)))
    raise ValueError(f"unknown geometry kind {kind!r}")
