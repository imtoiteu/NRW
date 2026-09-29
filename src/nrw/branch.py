r"""Phase-branch (integer ambiguity) determination for gamma L = -ln T.

With T = |T| e^{j theta}, theta in (-pi, pi], every integer m gives an admissible
solution beta L = psi + 2 pi m, psi = -theta (docs/conventions.md). The inverse problem
itself carries no information about m at a single frequency; m must come from
additional information:

* continuity in frequency (unwrapping)          -> fixes *changes* of m, not m itself;
* a known starting branch (e.g. m = 0 when the sample is electrically short, beta L < pi,
  at the lowest frequency)                       -> fixes the offset if that holds;
* group delay (Weir 1974): for weakly dispersive media the phase delay computed from a
  candidate m must equal the measured group delay -> fixes m pointwise (noisy) or
  the global offset (robust, when combined with unwrapping);
* continuity of the *material* quantity eps*mu between neighbouring frequencies
  (candidate-enumeration with dynamic programming) -> tolerates coarse frequency grids
  where adjacent phase steps exceed pi and unwrapping is ill-posed;
* a prior estimate of the refractive index       -> fixes m directly;
* a second sample of different length           -> see nrw.multithickness.

Observation used below (project analysis, see docs/theory_review.md): Weir's
group-delay criterion evaluated with a non-dispersive eps*mu is algebraically the
statement "phase delay = group delay", i.e. a local linear extrapolation of the phase
to f = 0 (TEM). It therefore misassigns m when the in-band dispersion makes
|f (tau_group - tau_phase)| approach 1/2 cycle.

All strategies return a :class:`BranchResult` exposing m(f) and diagnostics.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .geometry import C0, Geometry


@dataclass
class BranchResult:
    m: np.ndarray
    method: str
    # per-frequency ambiguity margin in [0, 1]: 1 = unambiguous, 0 = two candidates tie.
    margin: np.ndarray | None = None
    # indices k where the step k-1 -> k is unreliable for unwrapping (|dpsi| near/over pi,
    # or noise-dominated transmission).
    unreliable_steps: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    segments: list = field(default_factory=list)  # [(start, stop)] index ranges anchored independently
    notes: list = field(default_factory=list)

    @property
    def n_discontinuities(self) -> int:
        return int(np.count_nonzero(np.diff(self.m) != 0)) if self.m.size > 1 else 0


def principal_phase(T):
    """psi = -arg T in [-pi, pi)."""
    return -np.angle(T)


def candidate_eps_mu(f, geom: Geometry, T, length_m, m, candidates=False):
    """eps*mu and gamma for branch index m.

    If ``candidates`` is True, ``m`` is a 1-D array of K candidate integers and the
    outputs have shape (N, K); otherwise ``m`` broadcasts against T (shape (N,)).
    """
    lnT = np.log(np.asarray(T, complex))
    k0 = geom.k0(f)
    if candidates:
        lnT, k0 = lnT[:, None], k0[:, None]
    g = (-lnT + 2j * np.pi * np.asarray(m)) / length_m
    return (geom.kc**2 - g**2) / k0**2, g


def phase_delay_nondispersive(f, geom: Geometry, gamma, eps_mu, length_m):
    """Group delay predicted for a *non-dispersive* eps*mu (Weir 1974's calculated delay).

    From gamma^2 = kc^2 - k0^2 eps mu with d(eps mu)/d omega neglected,
    d gamma/d omega = -k0 eps mu / (c gamma), and tau = L d beta/d omega
    = -L Im(k0 eps mu / (c gamma)). For TEM this is the phase delay L n'/c.
    """
    k0 = geom.k0(f)
    if np.ndim(gamma) > np.ndim(k0):
        k0 = k0[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        return -length_m * np.imag(k0 * eps_mu / (C0 * gamma))


def measured_group_delay(f, T):
    """Group delay of T from finite differences of the *locally* unwrapped phase.

    Only adjacent-point phase steps are unwrapped (requires |dpsi| < pi between
    neighbours, i.e. delta_f * tau < 1/2). Returns (tau, max_abs_step).
    """
    psi_u = np.unwrap(principal_phase(T))
    tau = np.gradient(psi_u, f) / (2 * np.pi)
    steps = np.abs(np.diff(psi_u))
    return tau, steps


def _unreliable(f, T, noise_floor, step_limit):
    psi = principal_phase(T)
    d = np.abs(np.angle(np.exp(1j * np.diff(psi))))  # wrapped step magnitude in [0, pi]
    bad = np.flatnonzero(d > step_limit) + 1
    low = np.flatnonzero(np.abs(T) < noise_floor)
    return np.unique(np.concatenate([bad, low, low + 1]).clip(0, len(T) - 1)).astype(int)


def principal(f, T, length_m, geom, **_):
    return BranchResult(m=np.zeros(np.shape(f), dtype=int), method="principal",
                        notes=["m = 0 everywhere; valid only if beta L < pi across the band"])


def unwrap_from_start(f, T, length_m, geom, m0: int = 0, noise_floor=1e-3, step_limit=0.75 * np.pi, **_):
    """Conventional approach: np.unwrap of psi from the lowest frequency, with m(f0)=m0."""
    psi = principal_phase(T)
    psi_u = np.unwrap(psi) + 2 * np.pi * m0
    m = np.rint((psi_u - psi) / (2 * np.pi)).astype(int)
    return BranchResult(m=m, method="unwrap", unreliable_steps=_unreliable(f, T, noise_floor, step_limit),
                        notes=[f"assumes m = {m0} at the lowest frequency"])


def group_delay_pointwise(f, T, length_m, geom, m_max=None, **_):
    """Weir (1974): at each frequency choose m minimising |tau_calc(m) - tau_measured|."""
    f = np.asarray(f, float)
    tau_meas, _ = measured_group_delay(f, T)
    if m_max is None:
        m_max = _default_m_max(f, T, length_m, geom)
    cands = np.arange(0, m_max + 1)
    em, g = candidate_eps_mu(f, geom, T, length_m, cands, candidates=True)
    tau_c = phase_delay_nondispersive(f, geom, g, em, length_m)
    cost = np.abs(tau_c - tau_meas[:, None])
    cost = np.where(np.isfinite(cost), cost, np.inf)
    order = np.argsort(cost, axis=1)
    best = order[:, 0]
    c1 = np.take_along_axis(cost, order[:, :1], 1)[:, 0]
    c2 = np.take_along_axis(cost, order[:, 1:2], 1)[:, 0] if cands.size > 1 else np.full_like(c1, np.inf)
    margin = np.clip(1 - c1 / np.where(c2 > 0, c2, np.inf), 0, 1)
    return BranchResult(m=cands[best].astype(int), method="group_delay", margin=margin)


def _default_m_max(f, T, length_m, geom, n_max=50.0):
    """Upper bound on m from an assumed maximum refractive index (|n| <= n_max)."""
    return int(np.ceil(np.max(geom.k0(f)) * n_max * length_m / (2 * np.pi))) + 1


def _segments(n, breaks):
    edges = [0] + sorted(set(int(b) for b in breaks if 0 < b < n)) + [n]
    return [(a, b) for a, b in zip(edges[:-1], edges[1:]) if b > a]


def anchored_unwrap(f, T, length_m, geom, noise_floor=1e-3, step_limit=0.75 * np.pi, m_max=None,
                    weights=None, **_):
    """Unwrap within reliable segments; fix each segment's integer offset by group delay.

    Within a segment the unwrapped phase determines changes of m exactly; the single
    unknown integer offset is chosen to minimise the weighted median of
    |tau_calc - tau_measured| over the *whole segment* (robust to noise, unlike the
    pointwise rule). Segments are split at unreliable steps (|dpsi| > step_limit or
    |T| below the noise floor) so that one bad step cannot corrupt the rest of the band.
    """
    f = np.asarray(f, float)
    psi = principal_phase(T)
    unrel = _unreliable(f, T, noise_floor, step_limit)
    segs = _segments(len(f), unrel)
    m = np.zeros(len(f), dtype=int)
    margin = np.zeros(len(f))
    if m_max is None:
        m_max = _default_m_max(f, T, length_m, geom)
    w_all = np.abs(T) if weights is None else np.asarray(weights, float)
    notes = []
    for a, b in segs:
        idx = slice(a, b)
        psi_u = np.unwrap(psi[idx])
        rel = np.rint((psi_u - psi[idx]) / (2 * np.pi)).astype(int)
        if b - a < 3:
            m[idx] = rel  # too short to estimate a group delay; offset 0 by default
            notes.append(f"segment {a}:{b} too short to anchor; offset 0 assumed")
            continue
        tau_meas = np.gradient(psi_u, f[idx]) / (2 * np.pi)
        best_cost, best_m0, second = np.inf, 0, np.inf
        for m0 in range(-int(rel.min()), m_max + 1):
            mm = rel + m0
            em, g = candidate_eps_mu(f[idx], geom, T[idx], length_m, mm)
            tau_c = phase_delay_nondispersive(f[idx], geom, g, em, length_m)
            r = np.abs(tau_c - tau_meas)
            ok = np.isfinite(r)
            if not ok.any():
                continue
            cost = _weighted_median(r[ok], w_all[idx][ok])
            if cost < best_cost:
                second, best_cost, best_m0 = best_cost, cost, m0
            elif cost < second:
                second = cost
        m[idx] = rel + best_m0
        margin[idx] = 1 - best_cost / second if np.isfinite(second) and second > 0 else 1.0
    return BranchResult(m=m, method="anchored_unwrap", margin=margin, unreliable_steps=unrel,
                        segments=segs, notes=notes)


def _weighted_median(x, w):
    o = np.argsort(x)
    cw = np.cumsum(w[o])
    return x[o][np.searchsorted(cw, 0.5 * cw[-1])]


def dynamic_programming(f, T, length_m, geom, m_max=None, noise_floor=1e-3, rel_smooth=0.05,
                        gd_weight=1.0, prior_eps_mu=None, prior_rel_sigma=0.5, **_):
    """Global branch path by dynamic programming (Viterbi) over candidate integers.

    Cost of a path m_0..m_{N-1}:
        sum_k  gd_weight * ((tau_calc(m_k) - tau_meas_k) / s_tau)^2        (group delay)
             + sum_k ((eps_mu_k(m_k) - eps_mu_{k-1}(m_{k-1})) / s_k)^2      (continuity)
             + optional prior term on eps*mu.
    s_k = rel_smooth * |eps_mu| * max(1, df_k / median(df)) allows larger changes across
    gaps. Continuity is imposed on the *material* quantity eps*mu, not on the phase,
    so the method does not require |dpsi| < pi between neighbours.
    The group-delay term uses the locally unwrapped phase and is down-weighted
    automatically where |T| is noise-dominated.
    """
    f = np.asarray(f, float)
    n = len(f)
    if m_max is None:
        m_max = _default_m_max(f, T, length_m, geom)
    cands = np.arange(0, m_max + 1)
    em, g = candidate_eps_mu(f, geom, T, length_m, cands, candidates=True)  # (n, K)
    tau_meas, _ = measured_group_delay(f, T)
    tau_c = phase_delay_nondispersive(f, geom, g, em, length_m)
    scale_tau = np.maximum(np.abs(tau_meas), length_m / C0)  # relative mismatch
    node = gd_weight * ((tau_c - tau_meas[:, None]) / scale_tau[:, None]) ** 2
    rel = np.clip(np.abs(T) / noise_floor, 0, 1) ** 2  # confidence in measured group delay
    node = node * rel[:, None]
    if prior_eps_mu is not None:
        pe = np.broadcast_to(np.asarray(prior_eps_mu, complex), f.shape)
        node = node + (np.abs(em - pe[:, None]) / (prior_rel_sigma * np.abs(pe[:, None]))) ** 2
    node = np.where(np.isfinite(node), node, 1e12)
    # physical admissibility (exact for passive ordinary media): beta > 0 at m >= 0 is
    # automatic; reject candidates with Re(eps mu) <= 0 only via cost (not hard) since
    # plasma-like/negative-index media exist.
    df = np.diff(f)
    dfr = np.maximum(1.0, df / np.median(df)) if n > 1 else np.ones(0)
    cost = node[0].copy()
    back = np.zeros((n, len(cands)), dtype=int)
    for k in range(1, n):
        s = rel_smooth * dfr[k - 1] * np.maximum(np.abs(em[k])[None, :], np.abs(em[k - 1])[:, None])
        trans = (np.abs(em[k][None, :] - em[k - 1][:, None]) / s) ** 2  # (K_prev, K)
        trans = np.where(np.isfinite(trans), trans, 1e12)
        tot = cost[:, None] + trans
        back[k] = np.argmin(tot, axis=0)
        cost = tot[back[k], np.arange(len(cands))] + node[k]
    path = np.zeros(n, dtype=int)
    path[-1] = int(np.argmin(cost))
    srt = np.sort(cost)
    final_margin = 1 - srt[0] / srt[1] if len(srt) > 1 and srt[1] > 0 else 1.0
    for k in range(n - 1, 0, -1):
        path[k - 1] = back[k, path[k]]
    margin = np.full(n, final_margin)
    return BranchResult(m=cands[path].astype(int), method="dp", margin=margin,
                        unreliable_steps=_unreliable(f, T, noise_floor, 0.75 * np.pi),
                        notes=[f"path-cost margin (best vs runner-up end state) = {final_margin:.3g}"])


def from_prior(f, T, length_m, geom, prior_eps_mu=None, **_):
    """m nearest to the phase implied by a prior eps*mu estimate (e.g. datasheet value)."""
    if prior_eps_mu is None:
        raise ValueError("branch method 'prior' requires prior_eps_mu")
    pe = np.broadcast_to(np.asarray(prior_eps_mu, complex), np.shape(f))
    beta_p = np.real(np.sqrt(geom.k0(f) ** 2 * pe - geom.kc**2 + 0j))
    psi = principal_phase(T)
    m = np.rint((beta_p * length_m - psi) / (2 * np.pi)).astype(int)
    frac = np.abs((beta_p * length_m - psi) / (2 * np.pi) - m)
    return BranchResult(m=np.maximum(m, 0), method="prior", margin=1 - 2 * frac)


METHODS = {
    "principal": principal,
    "unwrap": unwrap_from_start,
    "continuous": unwrap_from_start,  # alias: "continuous unwrap from low frequency"
    "group_delay": group_delay_pointwise,
    "anchored": anchored_unwrap,
    "anchored_unwrap": anchored_unwrap,
    "dp": dynamic_programming,
    "prior": from_prior,
}


def select_branch(method: str, f, T, length_m, geom: Geometry, **kw) -> BranchResult:
    try:
        fn = METHODS[method]
    except KeyError as e:
        raise ValueError(f"unknown branch method {method!r}; choose from {sorted(METHODS)}") from e
    return fn(np.asarray(f, float), np.asarray(T, complex), float(length_m), geom, **kw)
