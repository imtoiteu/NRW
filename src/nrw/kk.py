r"""Finite-band Kramers-Kronig (KK) consistency diagnostics.

Convention: exp(+jwt), eps(w) = eps'(w) - j eps''(w), eps'' >= 0 (passive). For a causal
response chi(w) = eps(w) - eps_inf that is square-integrable (e.g. no DC conductivity
pole in chi'; see docs/theory_review.md), the KK relations read

    eps'(w) - eps_inf = (2/pi) P int_0^inf  w' eps''(w') / (w'^2 - w^2) dw'
    eps''(w)          = -(2 w/pi) P int_0^inf (eps'(w') - eps_inf) / (w'^2 - w^2) dw'

(identical in form to the exp(-iwt), eps = eps_1 + i eps_2 textbook statement because
eps'' here equals eps_2 there). The same holds for mu with mu_inf = 1.

Why a strict KK claim is impossible from 2-18 GHz data
------------------------------------------------------
The integrals run over (0, inf); measured data exist only on [w1, w2]. The
out-of-band loss contributes to eps' in-band a term that is *unknown* but, for a
passive medium (eps'' >= 0 everywhere), has a constrained form: loss at a frequency
W > w2 adds a non-negative multiple of W/(W^2 - w^2), loss at W < w1 adds a non-negative
multiple of W/(W^2 - w^2) (negative-valued in band). Hence finite-band data can only be
tested for *consistency with the existence of some passive out-of-band continuation*;
they can never prove that the material satisfies the global KK relations, and even a
failure is conditional on the passive-tail assumption and on the data uncertainty.

Diagnostic implemented here (singly-subtractive, anchored at w0 in band):

    eps'(w) - eps'(w0) = (2/pi) [I(w) - I(w0)] + tail(w) - tail(w0),
    I(w) = P int_{w1}^{w2} w' eps''(w')/(w'^2 - w^2) dw'   (exact for piecewise-linear eps'')
    tail(w) - tail(w0) = anchored KK transform of an out-of-band loss that is
    non-negative, continuous with the measured loss at the band edges and piecewise
    linear on log-spaced nodes (free coefficients c_i >= 0 at those nodes).

c_i are found by non-negative least squares. The normalised residual after the best
passive tail is the diagnostic. "Pass" means: consistent within the stated uncertainty,
bandwidth and assumptions - nothing more.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import nnls


def _pv_integral_piecewise_linear(x, y, w):
    """P int_{x0}^{xN} x' y(x') / (x'^2 - w^2) dx' for piecewise-linear y on nodes x.

    Uses the closed-form antiderivative on each segment: with y = p + q x',
    x'(p + q x')/(x'^2 - w^2) = q + p x'/(x'^2 - w^2) + q w^2/(x'^2 - w^2), whose
    antiderivative is q x' + (p/2) ln|x'^2 - w^2| + (q w/2) ln|(x' - w)/(x' + w)|.
    If w coincides with a node, the two adjacent segments carry equal and opposite
    ln|x' - w| -> -inf terms (their coefficients both equal y(w)/2 by continuity), so
    the principal value is obtained exactly by setting ln 0 := 0 at that node.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    w = np.atleast_1d(np.asarray(w, float))
    q = np.diff(y) / np.diff(x)  # (K,)
    p = y[:-1] - q * x[:-1]
    out = np.zeros(w.shape)
    for j, wj in enumerate(w):
        total = np.sum(q * np.diff(x))
        # terms in ln|x+w| (regular for w > 0) and ln|x-w| (singular at nodes)
        a_minus = (p + q * wj) / 2  # coefficient of ln|x - w| on each segment
        a_plus = (p - q * wj) / 2  # coefficient of ln|x + w|
        lxm = np.log(np.abs(x - wj), where=np.abs(x - wj) > 0, out=np.zeros_like(x))
        lxp = np.log(x + wj)
        total += np.sum(a_plus * (lxp[1:] - lxp[:-1]))
        # ln|x - w| terms: sum over segments of a_k (L_{k+1} - L_k), grouped by node
        contrib = a_minus * lxm[1:] - a_minus * lxm[:-1]
        total += np.sum(contrib)
        out[j] = total
    return out


def kk_real_from_imag_finite(f, loss, f_eval=None):
    """(2/pi) P int_{f1}^{f2} f' loss(f') / (f'^2 - f^2) df' (dimensionless; Hz units cancel)."""
    f = np.asarray(f, float)
    fe = f if f_eval is None else np.asarray(f_eval, float)
    return (2 / np.pi) * _pv_integral_piecewise_linear(f, loss, fe)


@dataclass
class KKResult:
    frequency: np.ndarray
    real_measured: np.ndarray
    real_predicted: np.ndarray  # anchored finite-band prediction including best passive tail
    residual: np.ndarray
    normalized_residual: np.ndarray  # residual / sigma
    tail_coefficients: np.ndarray
    anchor_index: int
    sigma: np.ndarray
    passed: bool
    statement: str
    assumptions: list = field(default_factory=list)
    max_abs_normalized_residual: float = np.nan
    rms_normalized_residual: float = np.nan


def _passive_tail_basis(f, loss, fa, n_tail, tail_span):
    """Basis for the in-band real-part contribution of a passive out-of-band loss.

    The out-of-band loss is modelled as a continuous, non-negative, piecewise-linear
    function on log-spaced nodes below f1 and above f2. It joins the measured loss
    continuously at the band edges (the edge values are fixed to the measured ones:
    a jump in loss would make the principal-value integral diverge at the edge) and its
    values at the out-of-band nodes are free non-negative coefficients.

    Returns (base, Phi): base = anchored finite contribution of the fixed edge ramps and
    in-band data; Phi[:, i] = anchored contribution of the i-th out-of-band hat.
    """
    f1, f2 = f[0], f[-1]
    nh = max(n_tail // 2, 2)
    low = np.geomspace(f1 / tail_span, f1, nh + 1)[:-1]
    high = np.geomspace(f2, f2 * tail_span, nh + 1)[1:]
    X = np.concatenate([[low[0] / 2], low, f, high, [high[-1] * 2]])  # zero end nodes
    i0 = 1 + len(low)
    y_base = np.zeros(len(X))
    y_base[i0 : i0 + len(f)] = loss
    ev = np.concatenate([f, [fa]])

    def anchored(y):
        v = kk_real_from_imag_finite(X, y, ev)
        return v[:-1] - v[-1]

    base = anchored(y_base)
    cols = []
    for j in list(range(1, 1 + len(low))) + list(range(i0 + len(f), i0 + len(f) + len(high))):
        y = np.zeros(len(X))
        y[j] = 1.0
        cols.append(anchored(y))
    return base, np.stack(cols, axis=1)


def kk_consistency(f, value, sigma=None, rel_tol=0.01, anchor="center", n_tail=24, tail_span=100.0,
                   threshold=3.0, quantity="eps"):
    """Finite-band, singly-subtractive KK consistency test with a passive out-of-band tail.

    value : complex eps (or mu) on the band, exp(+jwt) convention.
    sigma : standard uncertainty of Re(value) per frequency (e.g. from Monte Carlo). If
            None, ``rel_tol * |value|`` is used and the statement says so.
    Returns a KKResult whose ``statement`` is the only wording the software endorses.
    """
    f = np.asarray(f, float)
    v = np.asarray(value, complex)
    ok = np.isfinite(v)
    if ok.sum() < 5:
        raise ValueError("need at least 5 finite points for a KK diagnostic")
    if np.any(np.diff(f[ok]) <= 0):
        raise ValueError("frequencies must be strictly increasing")
    f, v = f[ok], v[ok]
    loss = -v.imag
    real = v.real
    sig = rel_tol * np.abs(v) if sigma is None else np.asarray(sigma, float)[ok]
    sig = np.maximum(sig, 1e-12)
    k0 = len(f) // 2 if anchor == "center" else int(anchor)
    fa = f[k0]
    base, Phi = _passive_tail_basis(f, loss, fa, n_tail, tail_span)
    pred0 = real[k0] + base
    colscale = np.max(np.abs(Phi), axis=0)
    colscale[colscale == 0] = 1
    A = (Phi / colscale) / sig[:, None]
    b = (real - pred0) / sig
    c, _ = nnls(A, b, maxiter=50 * A.shape[1])
    pred = pred0 + (Phi / colscale) @ c
    f1, f2 = f[0], f[-1]
    res = real - pred
    nres = res / sig
    mx = float(np.max(np.abs(nres)))
    rms = float(np.sqrt(np.mean(nres**2)))
    passed = bool(mx <= threshold)
    assumptions = [
        f"band limited to {f1/1e9:.4g}-{f2/1e9:.4g} GHz; out-of-band behaviour unknown",
        "out-of-band loss assumed passive (non-negative), continuous at the band edges and "
        f"piecewise-linear on {2 * max(n_tail // 2, 2)} log-spaced nodes within x{tail_span:g} of the band edges",
        f"{quantity}'' interpolated piecewise-linearly between measured points",
        f"singly-subtractive form anchored at {fa/1e9:.4g} GHz (removes the unknown {quantity}_inf)",
        "uncertainty: " + ("supplied sigma" if sigma is not None else f"assumed {rel_tol:g} x |{quantity}|"),
        f"pass threshold: max |residual|/sigma <= {threshold:g}",
    ]
    if passed:
        statement = ("KK consistency diagnostic passed within the available bandwidth and assumptions: "
                     f"Re({quantity}) is reproduced from Im({quantity}) by a finite-band subtractive KK "
                     "relation plus a passive out-of-band continuation, to within "
                     f"{threshold:g} standard uncertainties. This does NOT establish that the material "
                     "strictly satisfies the global Kramers-Kronig relations.")
    else:
        statement = ("KK consistency diagnostic FAILED within the available bandwidth and assumptions: no "
                     "passive out-of-band continuation reproduces Re from Im to within "
                     f"{threshold:g} standard uncertainties (max {mx:.3g}). Possible causes: extraction "
                     "error (e.g. wrong phase branch, half-wavelength artefacts), underestimated "
                     "uncertainty, systematic errors, or non-passive/ill-posed data.")
    return KKResult(f, real, pred, res, nres, c, k0, sig, passed, statement, assumptions, mx, rms)
