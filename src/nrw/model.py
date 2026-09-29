"""Closed-form single-slab S-parameter model used *inside* the iterative inverse solvers.

(The benchmark data are generated with the independent ABCD cascade in nrw.forward.)

Reference planes may be offset from the sample faces by distances d1, d2 in the empty
guide; the measured quantities are then

    S11 = R1^2 S11_0,  S22 = R2^2 S11_0,  S21 = S12 = R1 R2 S21_0,   R_i = exp(-gamma0 d_i).

Two combinations are invariant to the *individual* plane positions and depend only on
d1 + d2 = L_air - L (holder length minus sample length):

    S21 S12 - S11 S22 = R^2 (T^2 - Gamma^2) / (1 - Gamma^2 T^2)
    (S21 + S12) / 2   = R   T (1 - Gamma^2) / (1 - Gamma^2 T^2),      R = exp(-gamma0 (d1+d2)).
"""
from __future__ import annotations

import numpy as np

from .geometry import Geometry


def slab(f, geom: Geometry, eps, mu, length_m):
    """Return (Gamma, T, S11_0, S21_0) at the sample faces."""
    g = geom.gamma(f, eps, mu)
    z = np.asarray(mu) * geom.gamma0(f) / g
    gam = (z - 1) / (z + 1)
    t = np.exp(-g * length_m)
    den = 1 - gam**2 * t**2
    return gam, t, gam * (1 - t**2) / den, t * (1 - gam**2) / den


def measured(f, geom: Geometry, eps, mu, length_m, d1=0.0, d2=0.0):
    """(S11, S21, S12, S22) with reference-plane offsets d1, d2."""
    _, _, s11, s21 = slab(f, geom, eps, mu, length_m)
    g0 = geom.gamma0(f)
    r1, r2 = np.exp(-g0 * d1), np.exp(-g0 * d2)
    return s11 * r1**2, s21 * r1 * r2, s21 * r1 * r2, s11 * r2**2


def invariants(f, geom: Geometry, eps, mu, length_m, l_air_minus_l=0.0):
    """(S21 S12 - S11 S22, (S21+S12)/2) predicted by the model."""
    gam, t, _, _ = slab(f, geom, eps, mu, length_m)
    r = np.exp(-geom.gamma0(f) * l_air_minus_l)
    den = 1 - gam**2 * t**2
    return r**2 * (t**2 - gam**2) / den, r * t * (1 - gam**2) / den
