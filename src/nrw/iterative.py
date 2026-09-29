"""Iterative (nonlinear least-squares) solvers for the T/R problem.

All solvers fit the closed-form slab model (nrw.model) to measured S-parameters.
Unknowns are real vectors (eps', eps'' [, mu', mu'']) per frequency, with
eps = eps' - j eps'' (docs/conventions.md).

Important property (demonstrated in the benchmarks): because T = exp(-gamma L) is
periodic in beta with period 2 pi / L, *every* phase branch is an exact zero-residual
solution. An iterative solver therefore does not resolve the branch ambiguity; it
inherits the branch from its initial estimate. Branch selection (nrw.branch) must
precede or accompany any iterative refinement.

Convergence diagnostics are returned for every frequency; nothing is silently replaced.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix

from . import model
from .geometry import Geometry


def complex_lm(residual, x0, max_iter=100, xtol=1e-12, ftol=1e-15, lam0=1e-3, h=1e-7):
    """Vectorised Levenberg-Marquardt for residuals holomorphic in complex unknowns.

    residual(x) maps x of shape (N, P) (complex) to F of shape (N, M) (complex); each
    of the N rows (frequencies) is an independent problem. Because F is holomorphic,
    min sum |F|^2 over complex x is equivalent to the real least-squares problem in
    (Re x, Im x), and the Gauss-Newton step is  (J^H J + lam diag(J^H J)) dx = -J^H F
    with the complex Jacobian J = dF/dx (Wirtinger calculus). The real Jacobian has the
    same singular values as J (each doubled), so cond(J) is reported.

    Returns x, cost (= sum |F|^2 per row), n_iter, converged, cond.
    """
    x = np.array(x0, dtype=complex, copy=True)
    N, P = x.shape
    F = residual(x)
    cost = np.sum(np.abs(F) ** 2, axis=1)
    lam = np.full(N, lam0)
    active = np.isfinite(cost)
    conv = np.zeros(N, bool)
    nit = np.zeros(N, int)
    J = None
    for _ in range(max_iter):
        if not active.any():
            break
        hs = h * np.maximum(1.0, np.abs(x))
        J = np.empty(F.shape + (P,), complex)
        for p in range(P):
            xp = x.copy()
            xp[:, p] += hs[:, p]
            J[..., p] = (residual(xp) - F) / hs[:, p : p + 1]
        A = np.einsum("nmp,nmq->npq", J.conj(), J)
        g = np.einsum("nmp,nm->np", J.conj(), F)
        D = np.einsum("npp->np", A).real
        Ad = A + (lam[:, None] * np.maximum(D, 1e-300))[:, :, None] * np.eye(P)[None]
        with np.errstate(all="ignore"):
            try:
                dx = -np.linalg.solve(Ad, g[..., None])[..., 0]
            except np.linalg.LinAlgError:
                dx = np.full_like(x, np.nan)
        dx[~active] = 0
        xn = x + dx
        Fn = residual(xn)
        cn = np.sum(np.abs(Fn) ** 2, axis=1)
        better = active & np.isfinite(cn) & (cn <= cost)
        rel_step = np.max(np.abs(dx) / np.maximum(1.0, np.abs(x)), axis=1)
        rel_cost = (cost - cn) / np.maximum(cost, 1e-300)
        x[better] = xn[better]
        F[better] = Fn[better]
        done = better & ((rel_step < xtol) | (rel_cost < ftol) | (cn < 1e-30))
        cost = np.where(better, cn, cost)
        lam = np.where(better, lam / 3, lam * 4)
        nit += active
        conv |= done
        stuck = active & (lam > 1e12)
        active &= ~(done | stuck)
    if J is None:
        cond = np.full(N, np.nan)
    else:
        sv = np.linalg.svd(J, compute_uv=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            cond = sv[:, 0] / sv[:, -1]
    return x, cost, nit, conv, cond


@dataclass
class IterativeResult:
    frequency: np.ndarray
    eps: np.ndarray
    mu: np.ndarray
    converged: np.ndarray
    n_iterations: np.ndarray
    residual_rms: np.ndarray  # rms |S_model - S_meas| over the fitted S-parameters
    jacobian_cond: np.ndarray
    alternative_minimum: np.ndarray  # a different solution with comparable residual found
    method: str
    notes: list = field(default_factory=list)


def _stack(*arrs):
    return np.concatenate([np.concatenate([a.real.ravel(), a.imag.ravel()]) for a in arrs])


def _meas_list(s11, s21, s12, s22):
    out = [("s11", s11), ("s21", s21)]
    if s12 is not None:
        out.append(("s12", s12))
    if s22 is not None:
        out.append(("s22", s22))
    return out


def _model_list(f, geom, eps, mu, L, keys, d1, d2):
    m11, m21, m12, m22 = model.measured(f, geom, eps, mu, L, d1, d2)
    d = {"s11": m11, "s21": m21, "s12": m12, "s22": m22}
    return [d[k] for k in keys]


# ---------------------------------------------------------------------------------------
# Pointwise simultaneous eps, mu
# ---------------------------------------------------------------------------------------
def _pack_result(f, x, cost, nit, conv, cond, nobs, method, mu_fixed=False):
    eps = x[:, 0]
    mu = np.ones(len(f), complex) if mu_fixed else x[:, 1]
    res = np.sqrt(cost / nobs)
    return IterativeResult(f, eps, mu, conv, nit, res, cond, np.zeros(len(f), bool), method)


def _best_of(cands):
    """Element-wise choose the lowest-cost of several complex_lm outputs."""
    xs, costs, nits, convs, conds = (np.stack(z) for z in zip(*cands))
    c = np.where(np.isfinite(costs), costs, np.inf)
    i = np.argmin(c, axis=0)
    ar = np.arange(xs.shape[1])
    return xs[i, ar], costs[i, ar], nits[i, ar], convs[i, ar], conds[i, ar], i


def solve_eps_mu_pointwise(f, s11, s21, length_m, geom: Geometry, eps0, mu0, s12=None, s22=None,
                           d1=0.0, d2=0.0, multistart=True, max_iter=100):
    """Least-squares fit of (eps, mu) at each frequency to all supplied S-parameters.

    eps0, mu0: initial estimates per frequency (normally NRW with a vetted branch).
    multistart: additionally start from points displaced along the eps/mu trade-off
    valley (eps*mu fixed, mu scaled by 1.5 and 1/1.5). If a *different* solution with a
    residual within 10x of the best is found, ``alternative_minimum`` is flagged, which
    signals that eps and mu are not separately identifiable at that frequency.
    """
    f = np.asarray(f, float)
    meas = _meas_list(*(np.asarray(x, complex) if x is not None else None for x in (s11, s21, s12, s22)))
    keys = [k for k, _ in meas]
    y = np.stack([v for _, v in meas], axis=1)  # (N, M)

    def residual(x):
        mod = _model_list(f, geom, x[:, 0], x[:, 1], length_m, keys, d1, d2)
        return np.stack(mod, axis=1) - y

    e0 = np.broadcast_to(np.asarray(eps0, complex), f.shape)
    u0 = np.broadcast_to(np.asarray(mu0, complex), f.shape)
    x0 = np.stack([e0, u0], axis=1)
    x0 = np.where(np.isfinite(x0), x0, 1.0 + 0j)
    first = complex_lm(residual, x0, max_iter=max_iter)
    alt = np.zeros(len(f), bool)
    x, cost, nit, conv, cond = first
    if multistart:
        n2 = x[:, 0] * x[:, 1]
        others = []
        for s in (1.5, 1 / 1.5):
            xs = np.stack([n2 / (x[:, 1] * s), x[:, 1] * s], axis=1)
            others.append(complex_lm(residual, xs, max_iter=max_iter))
        for xo, co, _, cvo, _ in others:
            dx = np.max(np.abs(xo - x) / np.maximum(1, np.abs(x)), axis=1)
            alt |= cvo & (co <= np.maximum(10 * cost, 1e-24)) & (dx > 1e-3)
        x, cost, nit, conv, cond, _ = _best_of([first] + others)
    out = _pack_result(f, x, cost, nit, conv, cond, y.shape[1], "iterative_eps_mu")
    out.alternative_minimum = alt
    return out


# ---------------------------------------------------------------------------------------
# Permittivity-only (mu = 1) iterative solution: NIST TN 1341 / Baker-Jarvis et al. (1990)
# ---------------------------------------------------------------------------------------
def solve_eps_only(f, s11, s21, length_m, geom: Geometry, eps0, s12=None, s22=None, l_air_minus_l=None,
                   d1=0.0, d2=0.0, beta=0.0, formulation="eq2.38", max_iter=100, weight=None):
    """NIST (Baker-Jarvis et al. 1990; NIST TN 1341) permittivity-only iterative solution, mu_r = 1.

    Formulations (equation numbers of NIST TN 1341 [R4]; TN 1536 [R7] eqs. 97-98):

    "eq2.38" (default) - weighted transmission/reflection equation at the sample faces
        (known plane offsets d1, d2 are de-embedded first):
            1/2 {(S12 + S21) + beta (S11 + S22)} = [T(1 - G^2) + beta G (1 - T^2)] / (1 - T^2 G^2)
        beta = uncertainty(S21) / uncertainty(S11); TN 1341 recommends beta = 0 for low-loss
        samples (strong S21) and a large beta for high-loss samples.
    "eq2.25" - reference-plane-invariant determinant equation
            S11 S22 - S21 S12 = exp(-2 gamma0 (L_air - L)) (G^2 - T^2) / (1 - T^2 G^2),
        which needs only L_air - L (``l_air_minus_l``), not the individual plane positions.
    "eq2.24+2.25" - eq. (2.25) stacked with the invariant transmission equation (2.24),
            (S21 + S12)/2 = exp(-gamma0 (L_air - L)) T (1 - G^2) / (1 - T^2 G^2),
        solved in the least-squares sense (4 real equations, 2 unknowns).
    "direct" - least squares on S11, S21, S12, S22 with known d1, d2 (not a NIST form; control).

    Solver: vectorised Levenberg-Marquardt (complex_lm), i.e. damped Newton on the real and
    imaginary parts, started from the supplied estimate (normally the NRW / Boughriet value with
    a vetted branch). TN 1341 uses undamped Newton with frequency continuation; damping only
    enlarges the basin of convergence and does not change the solution.

    With mu = 1 the problem keeps full rank when S11 -> 0 at half-wavelength points (see
    nrw.conditioning), which is why the method does not diverge there - for non-magnetic
    samples only. Note: as with every T/R method, eps is defined only up to the phase branch;
    the solver stays on the branch of its initial estimate.
    """
    f = np.asarray(f, float)
    s11 = np.asarray(s11, complex)
    s21 = np.asarray(s21, complex)
    s12c = s21 if s12 is None else np.asarray(s12, complex)
    s22c = s11 if s22 is None else np.asarray(s22, complex)
    if weight is not None:  # backwards-compatible alias
        beta = weight
    if l_air_minus_l is None:
        l_air_minus_l = d1 + d2
    g0 = geom.gamma0(f)

    def gt(x):
        gam, t, _, _ = model.slab(f, geom, x[:, 0], 1.0, length_m)
        return gam, t

    if formulation == "eq2.38":
        r1, r2 = np.exp(-g0 * d1), np.exp(-g0 * d2)
        f11, f22, f21, f12 = s11 / r1**2, s22c / r2**2, s21 / (r1 * r2), s12c / (r1 * r2)
        y = (0.5 * ((f12 + f21) + beta * (f11 + f22)))[:, None]

        def residual(x):
            gam, t = gt(x)
            return ((t * (1 - gam**2) + beta * gam * (1 - t**2)) / (1 - t**2 * gam**2))[:, None] - y
    elif formulation in ("eq2.25", "eq2.24+2.25"):
        R = np.exp(-g0 * l_air_minus_l)
        det = (s11 * s22c - s21 * s12c)
        trn = (s21 + s12c) / 2

        def residual(x):
            gam, t = gt(x)
            den = 1 - t**2 * gam**2
            r_det = R**2 * (gam**2 - t**2) / den - det
            if formulation == "eq2.25":
                return r_det[:, None]
            return np.stack([r_det, R * t * (1 - gam**2) / den - trn], axis=1)
    elif formulation == "direct":
        y = np.stack([s11, s21, s12c, s22c], axis=1)

        def residual(x):
            return np.stack(model.measured(f, geom, x[:, 0], 1.0, length_m, d1, d2), axis=1) - y
    else:
        raise ValueError(f"unknown formulation {formulation!r}")
    e0 = np.broadcast_to(np.asarray(eps0, complex), f.shape)
    x0 = np.where(np.isfinite(e0), e0, 1.0 + 0j)[:, None]
    x, cost, nit, conv, cond = complex_lm(residual, x0, max_iter=max_iter)
    nobs = residual(x).shape[1]
    return _pack_result(f, x, cost, nit, conv, cond, nobs, f"nist_eps_only[{formulation}, beta={beta:g}]",
                        mu_fixed=True)


# ---------------------------------------------------------------------------------------
# Frequency-regularised simultaneous eps, mu (global fit)
# ---------------------------------------------------------------------------------------
def solve_eps_mu_regularized(f, s11, s21, length_m, geom: Geometry, eps0, mu0, s12=None, s22=None,
                             d1=0.0, d2=0.0, lam=1e-2, sigma=1e-3, max_nfev=50):
    """Global Tikhonov-regularised fit over all frequencies (4N real unknowns).

    minimise  sum_k |S_model(eps_k, mu_k) - S_meas,k|^2 / sigma^2
            + lam * sum_k |second difference of (eps', eps'', mu', mu'') along k|^2

    This is a *model assumption* (smoothness in frequency) and introduces bias at sharp
    features; it is offered to stabilise simultaneous extraction in half-wavelength
    regions where single-frequency data carry almost no information on eps/mu
    separately. The regularisation weight must be reported with any result.
    """
    f = np.asarray(f, float)
    meas = _meas_list(*(np.asarray(x, complex) if x is not None else None for x in (s11, s21, s12, s22)))
    keys = [k for k, _ in meas]
    N = len(f)
    nS = len(keys)
    y = np.stack([v for _, v in meas], axis=0)  # (nS, N)
    sq = np.sqrt(lam)

    def unpack(p):
        p = p.reshape(N, 4)
        return p[:, 0] - 1j * p[:, 1], p[:, 2] - 1j * p[:, 3], p

    def fun(p):
        e, u, P = unpack(p)
        mod = np.stack(_model_list(f, geom, e, u, length_m, keys, d1, d2), axis=0)
        r = (mod - y) / sigma  # (nS, N)
        data = np.concatenate([r.real.T, r.imag.T], axis=1).ravel()  # per-frequency blocks
        d2p = (P[2:] - 2 * P[1:-1] + P[:-2]) * sq
        return np.concatenate([data, d2p.ravel()])

    n_data = N * 2 * nS
    n_reg = max(N - 2, 0) * 4
    S = lil_matrix((n_data + n_reg, 4 * N), dtype=int)
    for k in range(N):
        S[k * 2 * nS : (k + 1) * 2 * nS, 4 * k : 4 * k + 4] = 1
    for k in range(N - 2):
        for j in range(4):
            row = n_data + 4 * k + j
            for kk in (k, k + 1, k + 2):
                S[row, 4 * kk + j] = 1
    e0 = np.broadcast_to(np.asarray(eps0, complex), f.shape)
    u0 = np.broadcast_to(np.asarray(mu0, complex), f.shape)
    p0 = np.stack([e0.real, -e0.imag, u0.real, -u0.imag], axis=1).ravel()
    r = least_squares(fun, p0, jac_sparsity=S.tocsr(), method="trf", x_scale="jac", max_nfev=max_nfev,
                      xtol=1e-10, ftol=1e-10)
    e, u, _ = unpack(r.x)
    mod = np.stack(_model_list(f, geom, e, u, length_m, keys, d1, d2), axis=0)
    res = np.sqrt(np.mean(np.abs(mod - y) ** 2, axis=0))
    out = IterativeResult(f, e, u, np.full(N, bool(r.success)), np.full(N, r.nfev), res,
                          np.full(N, np.nan), np.zeros(N, bool), f"regularized_eps_mu[lam={lam:g}]")
    out.notes.append(f"optimizer status {r.status}: {r.message}")
    return out


# ---------------------------------------------------------------------------------------
# Two (or more) sample lengths
# ---------------------------------------------------------------------------------------
def solve_eps_mu_multithickness(f, samples, geom: Geometry, eps0, mu0, max_iter=100):
    """Pointwise fit of one (eps, mu) to several samples of different length.

    samples: list of dicts {"s11", "s21", "length_m", optional "s12", "s22", "d1", "d2"}.
    With lengths whose half-wavelength points do not coincide, at least one sample has
    |S11| well above the noise at every frequency, restoring identifiability of eps/mu
    (two-sample idea discussed in NIST TN 1355-R; see docs/theory_review.md).
    """
    f = np.asarray(f, float)
    prepared = []
    ys = []
    for s in samples:
        meas = _meas_list(*(np.asarray(s[k], complex) if s.get(k) is not None else None
                            for k in ("s11", "s21", "s12", "s22")))
        prepared.append(([k for k, _ in meas], s["length_m"], s.get("d1", 0.0), s.get("d2", 0.0)))
        ys += [v for _, v in meas]
    y = np.stack(ys, axis=1)

    def residual(x):
        out = []
        for keys, L, d1, d2 in prepared:
            out += _model_list(f, geom, x[:, 0], x[:, 1], L, keys, d1, d2)
        return np.stack(out, axis=1) - y

    e0 = np.broadcast_to(np.asarray(eps0, complex), f.shape)
    u0 = np.broadcast_to(np.asarray(mu0, complex), f.shape)
    x0 = np.stack([e0, u0], axis=1)
    x0 = np.where(np.isfinite(x0), x0, 1.0 + 0j)
    x, cost, nit, conv, cond = complex_lm(residual, x0, max_iter=max_iter)
    return _pack_result(f, x, cost, nit, conv, cond, y.shape[1], "multithickness_eps_mu")
