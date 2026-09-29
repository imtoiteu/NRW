"""Independent forward model: S-parameters of a layered fill in a uniform guide.

The model cascades ABCD (chain) matrices of uniform line sections, each described by
its propagation constant gamma_i and its wave impedance normalised to the empty guide,
z_i = mu_i gamma0 / gamma_i (TEM and TE modes; docs/conventions.md). The two-port is then
converted to S-parameters with unit (empty-guide) reference impedance at both ports.

This derivation deliberately does NOT use the closed-form NRW slab expressions, so that
inverse routines are not validated against their own algebra. The closed-form slab
formula is provided separately (``slab_sparams_closed_form``) and the two are checked
against each other in the test-suite.

Non-idealities that can be modelled here (all optional):
* reference-plane offsets (empty-guide sections on either side of the sample),
* fixture/line attenuation in those sections (conductor-loss-like alpha ~ sqrt(f)),
* axial air gaps (extra empty sections between sample faces and reference planes),
* a radial air gap in a coaxial holder via the quasi-static layered-capacitor model
  (``coax_air_gap_effective``; approximation, see its docstring),
* sample-holder impedance mismatch (the ports see a reference impedance z_ref != 1).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .geometry import Coax, Geometry
from .materials import as_model


@dataclass
class Layer:
    """A uniform section of length ``length_m`` filled with (eps, mu).

    ``eps`` and ``mu`` may be complex scalars, arrays matching the frequency grid, or
    callables / ``nrw.materials.Model`` objects evaluated on the frequency grid.
    ``extra_alpha`` (Np/m, array or scalar) adds attenuation to gamma (fixture loss).
    """

    length_m: float
    eps: object = 1.0
    mu: object = 1.0
    extra_alpha: object = 0.0


def _evaluate(x, f):
    if callable(x):
        return np.asarray(as_model(x)(f), dtype=complex)
    arr = np.asarray(x, dtype=complex)
    return np.broadcast_to(arr, np.shape(f)).astype(complex)


def abcd_section(f, geom: Geometry, layer: Layer):
    """ABCD matrix (shape (N,2,2)) of one section, impedances normalised to empty guide."""
    eps = _evaluate(layer.eps, f)
    mu = _evaluate(layer.mu, f)
    g = geom.gamma(f, eps, mu) + _evaluate(layer.extra_alpha, f)
    z = mu * geom.gamma0(f) / geom.gamma(f, eps, mu)
    gl = g * layer.length_m
    ch, sh = np.cosh(gl), np.sinh(gl)
    m = np.empty(np.shape(f) + (2, 2), dtype=complex)
    m[..., 0, 0] = ch
    m[..., 0, 1] = z * sh
    m[..., 1, 0] = sh / z
    m[..., 1, 1] = ch
    return m


def abcd_to_s(m, z_ref=1.0):
    """Convert ABCD to S with (normalised) real reference impedance ``z_ref`` at both ports."""
    a, b, c, d = m[..., 0, 0], m[..., 0, 1] / z_ref, m[..., 1, 0] * z_ref, m[..., 1, 1]
    den = a + b + c + d
    s11 = (a + b - c - d) / den
    s21 = 2.0 / den
    s12 = 2.0 * (a * d - b * c) / den
    s22 = (-a + b - c + d) / den
    return s11, s21, s12, s22


def cascade_sparams(f, geom: Geometry, layers: Sequence[Layer], z_ref=1.0):
    """S-parameters (S11, S21, S12, S22) of the cascade of ``layers`` (port 1 -> port 2)."""
    f = np.asarray(f, dtype=float)
    m = np.broadcast_to(np.eye(2, dtype=complex), f.shape + (2, 2)).copy()
    for layer in layers:
        m = m @ abcd_section(f, geom, layer)
    return abcd_to_s(m, z_ref)


def slab_sparams_closed_form(f, geom: Geometry, eps, mu, length_m):
    """Textbook closed-form single-slab S-parameters (reference planes at the faces).

    S11 = Gamma (1 - T^2) / (1 - Gamma^2 T^2),  S21 = T (1 - Gamma^2) / (1 - Gamma^2 T^2).
    Used only to cross-check ``cascade_sparams``; never used to generate benchmark data.
    """
    f = np.asarray(f, dtype=float)
    eps = _evaluate(eps, f)
    mu = _evaluate(mu, f)
    g = geom.gamma(f, eps, mu)
    z = mu * geom.gamma0(f) / g
    gam = (z - 1) / (z + 1)
    t = np.exp(-g * length_m)
    den = 1 - gam**2 * t**2
    return gam * (1 - t**2) / den, t * (1 - gam**2) / den


def coax_air_gap_effective(eps, mu, geom: Coax, sample_inner_d_m, sample_outer_d_m):
    """Quasi-static layered-capacitor/inductor model of radial air gaps in a coax holder.

    The sample occupies radii [r2, r3] of a line with conductor radii r1 < r2, r3 < r4
    (gaps [r1, r2] and [r3, r4] are air). Treating the radial layers as series
    capacitances and the magnetic energy as parallel inductances gives (e.g. the
    capacitor model used in NIST TN 1355-R; see docs/theory_review.md):

        ln(r4/r1) / eps_eff = ln(r2/r1) + ln(r3/r2)/eps + ln(r4/r3)
        mu_eff ln(r4/r1)    = ln(r2/r1) + mu ln(r3/r2) + ln(r4/r3)

    Returns (eps_eff, mu_eff): what a homogeneous-fill extraction would "see".
    This is an approximation (quasi-static, TEM-like), adequate for small gaps.
    """
    r1, r4 = geom.inner_diameter_m / 2, geom.outer_diameter_m / 2
    r2, r3 = sample_inner_d_m / 2, sample_outer_d_m / 2
    if not (r1 <= r2 < r3 <= r4):
        raise ValueError("require inner conductor <= sample inner < sample outer <= outer conductor")
    l_tot = np.log(r4 / r1)
    l_in, l_s, l_out = np.log(r2 / r1), np.log(r3 / r2), np.log(r4 / r3)
    eps = np.asarray(eps, dtype=complex)
    mu = np.asarray(mu, dtype=complex)
    eps_eff = l_tot / (l_in + l_s / eps + l_out)
    mu_eff = (l_in + mu * l_s + l_out) / l_tot
    return eps_eff, mu_eff


@dataclass
class SampleSpec:
    """Convenience description of a single-sample T/R measurement."""

    eps: object
    mu: object
    length_m: float
    geometry: Geometry
    offset1_m: float = 0.0  # reference plane 1 lies this far from sample face 1 (air)
    offset2_m: float = 0.0
    fixture_alpha_np_per_m_at_1ghz: float = 0.0  # conductor-like loss ~ sqrt(f)
    z_ref: float = 1.0

    def sparams(self, f):
        f = np.asarray(f, dtype=float)
        alpha = self.fixture_alpha_np_per_m_at_1ghz * np.sqrt(f / 1e9)
        layers = []
        if self.offset1_m:
            layers.append(Layer(self.offset1_m, 1.0, 1.0, alpha))
        layers.append(Layer(self.length_m, self.eps, self.mu))
        if self.offset2_m:
            layers.append(Layer(self.offset2_m, 1.0, 1.0, alpha))
        return cascade_sparams(f, self.geometry, layers, self.z_ref)

    def truth(self, f):
        return _evaluate(self.eps, np.asarray(f, float)), _evaluate(self.mu, np.asarray(f, float))
