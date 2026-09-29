"""Physical-consistency diagnostics for measured S-parameters and extracted eps, mu.

Each check is labelled with its epistemic status (docs/final_method.md, docs/theory_review.md):

* ``exact``      - follows from passivity/linearity of the sample itself; a violation
                   larger than the uncertainty means the data or the extraction is wrong.
* ``assumption`` - holds only under a modelling assumption (reciprocal, homogeneous,
                   symmetric sample and fixture); a violation means the assumption fails.
* ``heuristic``  - plausibility rule of thumb; a violation is a prompt to inspect,
                   not proof of error (e.g. eps' < 1 is physical above a resonance).

Sign conventions: exp(+jwt), eps = eps' - j eps''. For a passive, isotropic, linear medium
the time-averaged absorbed power density is (w/2)(eps0 eps''|E|^2 + mu0 mu''|H|^2) >= 0
for *all* field amplitudes, hence eps'' >= 0 and mu'' >= 0 individually.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Check:
    name: str
    status: str  # exact | assumption | heuristic
    description: str
    violated: np.ndarray  # bool per frequency
    metric: np.ndarray | None = None

    @property
    def n_violations(self) -> int:
        return int(np.count_nonzero(self.violated))

    def to_dict(self):
        return {"name": self.name, "status": self.status, "description": self.description,
                "n_violations": self.n_violations}


@dataclass
class PhysicsReport:
    checks: list = field(default_factory=list)

    def __getitem__(self, name) -> Check:
        for c in self.checks:
            if c.name == name:
                return c
        raise KeyError(name)

    def summary(self) -> dict:
        return {c.name: c.to_dict() for c in self.checks}

    @property
    def n_exact_violations(self) -> int:
        return int(sum(c.n_violations for c in self.checks if c.status == "exact"))


def _tol(x, sigma, k):
    return k * np.asarray(sigma) if sigma is not None else np.full(np.shape(x), 0.0)


def check_network(s11, s21, s12=None, s22=None, tol=1e-3, noise_sigma=None, k_sigma=3.0):
    """Checks on the measured two-port itself (before any extraction).

    Measurement noise of rms sigma perturbs |S|^2 by about 2|S| sigma, so for noisy data a
    violation is declared only beyond max(tol, 2 k_sigma sigma).
    """
    if noise_sigma is not None:
        tol = max(tol, 2 * k_sigma * noise_sigma)
    s11 = np.asarray(s11, complex)
    s21 = np.asarray(s21, complex)
    s12 = s21 if s12 is None else np.asarray(s12, complex)
    s22 = s11 if s22 is None else np.asarray(s22, complex)
    S = np.stack([np.stack([s11, s12], -1), np.stack([s21, s22], -1)], -2)
    # passive network <=> I - S^H S positive semidefinite <=> largest singular value <= 1
    smax = np.linalg.svd(S, compute_uv=False)[:, 0]
    absorbed = 1 - np.abs(s11) ** 2 - np.abs(s21) ** 2
    recip = np.abs(s21 - s12) / np.maximum(np.abs(s21), 1e-12)
    sym = np.abs(s11 - s22) / np.maximum(np.abs(s11), np.abs(s22)).clip(1e-12)
    return [
        Check("network_passivity", "exact", f"largest singular value of S <= 1 + {tol:.2g}", smax > 1 + tol, smax),
        Check("absorbed_power_port1", "exact", "1 - |S11|^2 - |S21|^2 >= -tol", absorbed < -tol, absorbed),
        Check("reciprocity", "assumption", "|S21 - S12|/|S21| <= 5 % (linear, non-gyrotropic sample; stable VNA)",
              recip > 0.05, recip),
        Check("port_symmetry", "assumption",
              "|S11 - S22|/max|S11|,|S22| <= 20 % (homogeneous sample, symmetric planes)", sym > 0.2, sym),
    ]


def check_material(eps, mu, gamma=None, z=None, eps_sigma=None, mu_sigma=None, k_sigma=2.0,
                   eps_mu_fixed=False, n_max=50.0, alpha_tol=None):
    """Checks on extracted parameters.

    eps_sigma / mu_sigma: per-frequency standard uncertainty of eps'' / mu'' (if known);
    a passivity violation is only declared when eps'' < -k_sigma * sigma.
    """
    eps = np.asarray(eps, complex)
    mu = np.asarray(mu, complex)
    epp = -eps.imag  # eps''
    mpp = -mu.imag
    out = [
        Check("eps_passivity", "exact", "eps'' >= 0 (exp(+jwt) convention), within k*sigma",
              epp < -_tol(epp, eps_sigma, k_sigma), epp),
    ]
    if not eps_mu_fixed:
        out.append(Check("mu_passivity", "exact", "mu'' >= 0, within k*sigma", mpp < -_tol(mpp, mu_sigma, k_sigma), mpp))
    if gamma is not None:
        g = np.asarray(gamma)
        tol_a = 0.0 if alpha_tol is None else np.asarray(alpha_tol)
        out.append(Check("attenuation_nonnegative", "exact", "Re(gamma) >= 0 (|T| <= 1), beyond noise",
                         g.real < -tol_a, g.real))
    if z is not None:
        zz = np.asarray(z)
        out.append(Check("impedance_real_part", "exact", "Re(z) >= 0 (|Gamma| <= 1)", zz.real < 0, zz.real))
    finite = np.isfinite(eps) & np.isfinite(mu)
    out.append(Check("finite", "exact", "eps and mu finite", ~finite))
    n = np.sqrt(eps * mu)
    out.append(Check("refractive_index_plausible", "heuristic", f"|n| <= {n_max}", np.abs(n) > n_max, np.abs(n)))
    out.append(Check("eps_real_ge_1", "heuristic", "eps' >= 1 (violated legitimately above resonances/plasmas)",
                     eps.real < 1, eps.real))
    if not eps_mu_fixed:
        out.append(Check("mu_real_positive", "heuristic", "mu' > 0 (violated legitimately above FMR)",
                         mu.real <= 0, mu.real))
    return out


def continuity_check(f, x, name, k=8.0):
    """Heuristic: flag jumps |x_k - x_{k-1}| larger than k * (robust local scale)."""
    x = np.asarray(x, complex)
    d = np.abs(np.diff(x))
    finite = np.isfinite(d)
    if finite.sum() < 3:
        return Check(f"{name}_continuity", "heuristic", "no isolated jumps", np.zeros(len(x), bool))
    med = np.median(d[finite])
    mad = np.median(np.abs(d[finite] - med)) * 1.4826
    scale = max(med + mad, 1e-12 * max(1.0, np.nanmax(np.abs(x))))
    jump = np.concatenate([[False], d > k * scale])
    return Check(f"{name}_continuity", "heuristic", f"|delta {name}| <= {k:g} x robust step scale", jump,
                 np.concatenate([[0.0], d]))


def run_all(f, s11, s21, s12, s22, eps, mu, gamma=None, z=None, eps_sigma=None, mu_sigma=None,
            eps_mu_fixed=False, noise_sigma=None, length_m=None):
    """All checks. With ``noise_sigma`` (and ``length_m``), exact checks on measured quantities
    tolerate deviations up to 3 sigma of the measurement noise: |S|^2 by 6 sigma and the
    attenuation by 3 sigma / (|S21| L) (since d ln|T| ~ d|S21| / |S21|)."""
    rep = PhysicsReport()
    rep.checks += check_network(s11, s21, s12, s22, noise_sigma=noise_sigma)
    alpha_tol = None
    if noise_sigma is not None and length_m:
        alpha_tol = 3 * noise_sigma / (np.maximum(np.abs(np.asarray(s21)), 1e-12) * length_m)
    rep.checks += check_material(eps, mu, gamma, z, eps_sigma, mu_sigma, eps_mu_fixed=eps_mu_fixed,
                                 alpha_tol=alpha_tol)
    rep.checks.append(continuity_check(f, eps, "eps"))
    if not eps_mu_fixed:
        rep.checks.append(continuity_check(f, mu, "mu"))
    return rep
