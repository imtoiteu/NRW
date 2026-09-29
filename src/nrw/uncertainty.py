"""Monte Carlo propagation of measurement uncertainty through the *entire* pipeline.

Every draw perturbs the inputs and re-runs the complete extraction, including phase-branch
selection, so that branch instabilities appear as multimodal output distributions
(reported as ``branch_change_fraction``) instead of being hidden by a linearisation.
Robust summaries (median, percentile intervals, MAD-based sigma) are reported because
the output distributions near half-wavelength points are heavy-tailed.

A linearised (Cramer-Rao) bound for white noise is available separately in
nrw.conditioning and is reported next to the Monte Carlo result for comparison.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .perturb import shift_reference_planes


@dataclass
class UncertaintySpec:
    """Standard uncertainties (1 sigma) of the inputs. Zero disables a component.

    Magnitude uncertainties are relative (e.g. 0.005 = 0.5 %); phases in degrees.
    ``s_noise`` is an additive circular complex Gaussian noise rms on every S-parameter.
    Reference-plane uncertainties are independent for the two ports (m).
    """

    s11_mag_rel: float = 0.0
    s11_phase_deg: float = 0.0
    s21_mag_rel: float = 0.0
    s21_phase_deg: float = 0.0
    s_noise: float = 0.0
    thickness_m: float = 0.0
    ref_plane_m: float = 0.0
    n_draws: int = 200
    seed: int = 12345

    def to_dict(self):
        return asdict(self)


def _perturb_pair(s_refl, s_trans, spec, rng):
    n = s_refl.shape
    out_r = s_refl * (1 + spec.s11_mag_rel * rng.standard_normal(n)) * np.exp(
        1j * np.deg2rad(spec.s11_phase_deg) * rng.standard_normal(n))
    out_t = s_trans * (1 + spec.s21_mag_rel * rng.standard_normal(n)) * np.exp(
        1j * np.deg2rad(spec.s21_phase_deg) * rng.standard_normal(n))
    if spec.s_noise:
        out_r = out_r + spec.s_noise / np.sqrt(2) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
        out_t = out_t + spec.s_noise / np.sqrt(2) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return out_r, out_t


def monte_carlo(extract_fn, f, s11, s21, s12, s22, thickness_m, geom, spec: UncertaintySpec,
                nominal_branch=None):
    """Run ``extract_fn(f, s11, s21, s12, s22, thickness_m) -> (eps, mu, branch)`` on perturbed data.

    Returns a dict of per-frequency summaries for eps', eps'', mu', mu''.
    """
    rng = np.random.default_rng(spec.seed)
    f = np.asarray(f, float)
    N = len(f)
    E = np.full((spec.n_draws, N), np.nan + 0j)
    M = np.full((spec.n_draws, N), np.nan + 0j)
    B = np.zeros((spec.n_draws, N), dtype=int)
    failures = 0
    for i in range(spec.n_draws):
        a11, a21 = _perturb_pair(s11, s21, spec, rng)
        a22, a12 = _perturb_pair(s22, s12, spec, rng)
        if spec.ref_plane_m:
            d1, d2 = spec.ref_plane_m * rng.standard_normal(2)
            a11, a21, a12, a22 = shift_reference_planes(f, a11, a21, a12, a22, geom, d1, d2)
        L = thickness_m + spec.thickness_m * rng.standard_normal()
        try:
            e, m, b = extract_fn(f, a11, a21, a12, a22, L)
        except Exception:  # noqa: BLE001 - a failed draw is counted, not hidden
            failures += 1
            continue
        E[i], M[i], B[i] = e, m, b
    out = {"n_draws": spec.n_draws, "n_failed_draws": failures, "spec": spec.to_dict()}
    for name, arr in (("eps_real", E.real), ("eps_loss", -E.imag), ("mu_real", M.real), ("mu_loss", -M.imag)):
        with np.errstate(invalid="ignore"):
            med = np.nanmedian(arr, axis=0)
            lo, hi = np.nanpercentile(arr, [2.5, 97.5], axis=0)
            q16, q84 = np.nanpercentile(arr, [15.87, 84.13], axis=0)
            mad = 1.4826 * np.nanmedian(np.abs(arr - med), axis=0)
        out[name] = {"median": med, "p2.5": lo, "p97.5": hi, "p15.87": q16, "p84.13": q84,
                     "sigma_robust": mad, "sigma_percentile": (q84 - q16) / 2}
    out["nonfinite_fraction"] = np.mean(~np.isfinite(E) | ~np.isfinite(M), axis=0)
    if nominal_branch is not None:
        out["branch_change_fraction"] = np.mean(B != np.asarray(nominal_branch)[None, :], axis=0)
    return out
