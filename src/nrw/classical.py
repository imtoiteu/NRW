"""Classical Nicolson-Ross-Weir (NRW) closed-form extraction with explicit diagnostics.

Equations (Nicolson & Ross 1970; Weir 1974; notation of docs/conventions.md):

    X     = (S11^2 - S21^2 + 1) / (2 S11)
    Gamma = X -/+ sqrt(X^2 - 1)                (root with |Gamma| <= 1, i.e. Re z >= 0)
    T     = (S11 + S21 - Gamma) / (1 - (S11 + S21) Gamma)
    gamma L = -ln T = -ln|T| + j(-arg T + 2 pi m)          (branch index m)
    z     = (1 + Gamma) / (1 - Gamma)          (sample impedance / empty-guide impedance)
    mu    = z gamma / gamma0
    eps   = (kc^2 - gamma^2) / (k0^2 mu)

For TEM (kc = 0) these reduce to n = gamma/(j k0), mu = n z, eps = n / z.

Nothing is clipped or replaced: non-finite values propagate to the caller and every
suspicious condition is reported as a boolean flag array in ``NRWResult.flags``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .geometry import Geometry

# Default diagnostic thresholds. These are *flag* thresholds, not corrections; they are
# expressed relative to a user-specified S-parameter noise floor where that makes sense.
DEFAULT_NOISE_FLOOR = 1e-3  # |S| below which a measured S-parameter is considered noise-dominated
GAMMA_UNIT_CIRCLE_TOL = 1e-3  # 1-|Gamma| below which the two roots are nearly indistinguishable
BLOWUP_LIMIT = 1e4  # |eps| or |mu| above this is flagged as numerical blow-up


@dataclass
class NRWResult:
    frequency: np.ndarray
    X: np.ndarray
    reflection_coefficient: np.ndarray  # Gamma (selected root)
    reflection_coefficient_alt: np.ndarray  # the rejected root (1/Gamma)
    transmission: np.ndarray  # T
    ln_T_principal: np.ndarray  # principal complex log of T
    branch: np.ndarray  # integer m
    propagation_constant: np.ndarray  # gamma (1/m)
    impedance: np.ndarray  # z normalised to empty guide
    eps: np.ndarray
    mu: np.ndarray
    eps_mu_only: np.ndarray  # eps computed with mu forced to 1 (permittivity-only reading)
    electrical_length_over_pi: np.ndarray  # beta L / pi (half-wave condition at integers)
    flags: dict = field(default_factory=dict)

    @property
    def refractive_index(self):
        return np.sqrt(self.eps * self.mu)


def gamma_from_s(s11, s21, root_tol=GAMMA_UNIT_CIRCLE_TOL):
    """Return (X, Gamma, Gamma_alt, flags) from S11, S21 using the NRW X-parameter.

    The two roots satisfy Gamma * Gamma_alt = 1, so exactly one has |Gamma| <= 1 unless
    both lie on the unit circle. Choosing |Gamma| <= 1 is equivalent to Re z >= 0,
    an exact passivity condition for the semi-infinite sample (not a heuristic).
    """
    s11 = np.asarray(s11, dtype=complex)
    s21 = np.asarray(s21, dtype=complex)
    with np.errstate(divide="ignore", invalid="ignore"):
        X = (s11**2 - s21**2 + 1) / (2 * s11)
        root = np.sqrt(X**2 - 1)
        g1, g2 = X + root, X - root
    pick1 = np.abs(g1) <= np.abs(g2)
    gam = np.where(pick1, g1, g2)
    alt = np.where(pick1, g2, g1)
    flags = {
        "gamma_root_ambiguous": np.abs(np.abs(gam) - np.abs(alt)) < 2 * root_tol,
        "gamma_outside_unit_circle": np.abs(gam) > 1 + 1e-12,
    }
    return X, gam, alt, flags


def transmission_from_s(s11, s21, gam):
    with np.errstate(divide="ignore", invalid="ignore"):
        return (s11 + s21 - gam) / (1 - (s11 + s21) * gam)


def material_from_gamma_T(f, geom: Geometry, gam, T, m, length_m):
    """Compute (gamma, z, eps, mu, eps_mu1) for a given integer branch array m."""
    f = np.asarray(f, dtype=float)
    lnT = np.log(T)  # principal: -pi < Im <= pi
    gl = -lnT + 2j * np.pi * np.asarray(m)
    g = gl / length_m
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (1 + gam) / (1 - gam)
        g0 = geom.gamma0(f)
        mu = z * g / g0
        eps = geom.eps_from_gamma_mu(f, g, mu)
        eps_mu1 = geom.eps_from_gamma_mu(f, g, 1.0)
    return lnT, g, z, eps, mu, eps_mu1


def nrw(
    f,
    s11,
    s21,
    length_m,
    geom: Geometry,
    branch=None,
    noise_floor=DEFAULT_NOISE_FLOOR,
):
    """Classical NRW extraction.

    Parameters
    ----------
    branch : None, int or int array
        Phase-branch index m per frequency. ``None`` means m = 0 everywhere (the naive
        principal-value choice, which is only correct when beta L < pi at every point);
        use :mod:`nrw.branch` to determine it properly.
    noise_floor : float
        |S| level used for flagging noise-dominated S11 / S21.
    """
    f = np.asarray(f, dtype=float)
    X, gam, alt, flags = gamma_from_s(s11, s21)
    T = transmission_from_s(np.asarray(s11, complex), np.asarray(s21, complex), gam)
    m = np.zeros(f.shape, dtype=int) if branch is None else np.broadcast_to(np.asarray(branch, dtype=int), f.shape)
    lnT, g, z, eps, mu, eps_mu1 = material_from_gamma_T(f, geom, gam, T, m, length_m)

    # beta L / pi: integer values are the half-guided-wavelength (NRW-singular) points.
    el = np.imag(g) * length_m / np.pi
    s11a = np.abs(np.asarray(s11))
    flags.update(
        {
            # S11 in the noise: X = (...)/(2 S11) and hence Gamma and z are ill-determined.
            "s11_below_noise": s11a < noise_floor,
            # Near a half-wavelength resonance (low-loss samples): |beta L/pi - k| < 0.05.
            "near_half_wavelength": (np.abs(el - np.round(el)) < 0.05) & (np.round(el) >= 1),
            # S21 in the noise: arg T and hence the branch are ill-determined.
            "s21_below_noise": np.abs(np.asarray(s21)) < noise_floor,
            "nonfinite": ~(np.isfinite(eps) & np.isfinite(mu)),
            "blowup": (np.abs(eps) > BLOWUP_LIMIT) | (np.abs(mu) > BLOWUP_LIMIT),
            "T_gain": np.abs(T) > 1 + 1e-9,  # |T|>1: negative attenuation (non-passive)
        }
    )
    return NRWResult(
        frequency=f,
        X=X,
        reflection_coefficient=gam,
        reflection_coefficient_alt=alt,
        transmission=T,
        ln_T_principal=lnT,
        branch=np.asarray(m, dtype=int),
        propagation_constant=g,
        impedance=z,
        eps=eps,
        mu=mu,
        eps_mu_only=eps_mu1,
        electrical_length_over_pi=el,
        flags=flags,
    )
