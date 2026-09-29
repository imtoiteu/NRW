"""Local identifiability / conditioning of the T/R inverse problem.

For additive circular complex Gaussian noise of standard deviation sigma on each measured
S-parameter, the Fisher information for the real parameter vector p is
F = (2/sigma^2) J^T J where J is the real Jacobian of the stacked real/imaginary parts of
the S-parameters with respect to p. Its inverse is the Cramer-Rao lower bound (CRLB) for
any unbiased estimator: no algorithm - classical, iterative or "robust" - can do better
from the same single-sample data. This lets us separate *algorithmic* instability from
*intrinsic* ill-posedness (project analysis; see docs/theory_review.md).

Parameter sets:
* "eps_mu": p = (eps', eps'', mu', mu'')  (simultaneous extraction)
* "eps":    p = (eps', eps'') with mu = 1  (permittivity-only, Baker-Jarvis 1990 setting)
"""
from __future__ import annotations

import numpy as np

from . import model
from .geometry import Geometry


def _complex_derivs(f, geom, eps, mu, length_m, h=1e-7):
    """dS11/d eps, dS21/d eps, dS11/d mu, dS21/d mu (S analytic in eps and mu)."""
    eps = np.asarray(eps, complex)
    mu = np.asarray(mu, complex)
    _, _, s11, s21 = model.slab(f, geom, eps, mu, length_m)
    he = h * np.maximum(1, np.abs(eps))
    hm = h * np.maximum(1, np.abs(mu))
    _, _, s11e, s21e = model.slab(f, geom, eps + he, mu, length_m)
    _, _, s11m, s21m = model.slab(f, geom, eps, mu + hm, length_m)
    return (s11e - s11) / he, (s21e - s21) / he, (s11m - s11) / hm, (s21m - s21) / hm


def real_jacobian(f, geom: Geometry, eps, mu, length_m, params="eps_mu", use_s22_s12=True):
    """Real Jacobian, shape (N, n_obs, n_par).

    Observations: Re/Im of S11, S21 (and S22, S12 if ``use_s22_s12``; for a symmetric
    reciprocal sample these repeat S11, S21 and simply double the information).
    Parameters are eps' and eps'' with eps = eps' - j eps'' (and likewise mu), so
    d/d eps'' = -j d/d eps.
    """
    d11e, d21e, d11m, d21m = _complex_derivs(f, geom, eps, mu, length_m)
    cols = [d11e, -1j * d11e, d11m, -1j * d11m] if params == "eps_mu" else [d11e, -1j * d11e]
    cols21 = [d21e, -1j * d21e, d21m, -1j * d21m] if params == "eps_mu" else [d21e, -1j * d21e]
    rows = []
    for c11, c21 in zip(cols, cols21):
        blk = [c11.real, c11.imag, c21.real, c21.imag]
        if use_s22_s12:
            blk += [c11.real, c11.imag, c21.real, c21.imag]
        rows.append(np.stack(blk, axis=-1))
    return np.stack(rows, axis=-1)


def crlb(f, geom: Geometry, eps, mu, length_m, sigma, params="eps_mu", use_s22_s12=True):
    """Cramer-Rao standard deviations (N, n_par) and Jacobian condition numbers (N,).

    ``sigma``: rms of the complex noise on each S-parameter (E|n|^2 = sigma^2, so each
    real component has variance sigma^2/2).
    """
    J = real_jacobian(f, geom, eps, mu, length_m, params, use_s22_s12)
    F = np.einsum("nop,noq->npq", J, J) / (sigma**2 / 2)
    sv = np.linalg.svd(J, compute_uv=False)
    cond = sv[:, 0] / np.where(sv[:, -1] > 0, sv[:, -1], np.nan)
    with np.errstate(invalid="ignore"):
        try:
            cov = np.linalg.inv(F)
        except np.linalg.LinAlgError:  # pragma: no cover
            cov = np.full_like(F, np.inf)
    std = np.sqrt(np.clip(np.einsum("npp->np", cov), 0, None))
    return std, cond
