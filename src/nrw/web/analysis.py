"""Analysis layer of the web application.

Everything numerical is delegated to the tested core (``nrw.extract_material_parameters``,
``nrw.assessment.assess``). This module only
* turns UI settings into a validated, hashable configuration,
* restricts the frequency window,
* builds the result / diagnostics tables,
* runs several core methods for the comparison view.
"""
from __future__ import annotations

import hashlib
import json
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from ..api import ExtractionResult, extract_material_parameters
from ..assessment import Assessment, assess, passivity_classes
from ..io import TwoPortData
from ..uncertainty import UncertaintySpec

# ------------------------------------------------------------------------------------------
# Method and branch catalogue (only what is implemented and tested in the core)
# ------------------------------------------------------------------------------------------
METHODS = {
    "auto": ("AUTO (recommended)",
             "Evidence-based default (docs/final_method.md): if the sample is declared non-magnetic, the NIST "
             "mu = 1 iterative method; otherwise simultaneous eps/mu least squares started from NRW, with "
             "ill-conditioned (half-wavelength) points flagged via the Cramer-Rao bound."),
    "nrw": ("Classical NRW (closed form)",
            "Nicolson-Ross (1970) / Weir (1974) explicit equations for eps and mu. Diverges near beta L = k pi "
            "for low-loss samples because S11 -> 0 makes the eps/mu split unidentifiable (not hidden: flagged)."),
    "iterative_eps_mu": ("Iterative eps & mu least squares",
                         "Fits the slab model to S11, S21, S12, S22 at every frequency (damped Newton), started "
                         "from NRW on the selected branch. Same information limit as NRW near half-wave points; "
                         "reports convergence, residual, condition number and alternative minima."),
    "baker_jarvis": ("Baker-Jarvis / NIST iterative (mu = 1)",
                     "NIST TN 1341 / Baker-Jarvis et al. (1990) permittivity-only method, reference-plane "
                     "invariant eqs. (2.24)+(2.25). Assumes mu_r = 1: stable through half-wave points for "
                     "NON-MAGNETIC samples only; on magnetic samples it returns a biased effective permittivity."),
    "nrw_eps": ("Non-iterative eps with mu = 1 (Boughriet)",
                "Closed-form permittivity from the propagation constant with mu_r = 1 (Boughriet et al. 1997, as "
                "restated in NIST TN 1536). Non-magnetic samples only."),
    "regularized_eps_mu": ("Frequency-regularised eps & mu",
                           "Global Tikhonov fit with a smoothness penalty across frequency. A MODEL ASSUMPTION "
                           "that can bias sharp features; use only for comparison and report the weight."),
}
EPS_ONLY = ("baker_jarvis", "nrw_eps")

BRANCH_METHODS = {
    "anchored": ("Automatic robust tracking (recommended)",
                 "Unwraps within reliable segments and anchors each segment's integer offset by the band-median "
                 "group delay; cross-checked by an independent dynamic-programming tracker (disagreements flagged)."),
    "dp": ("Branch enumeration / optimisation (dynamic programming)",
           "Chooses the branch path minimising eps*mu discontinuity plus group-delay mismatch over all candidate "
           "integers (Viterbi). Needed when the frequency grid is too coarse for unwrapping."),
    "unwrap": ("Continuous unwrap from the lowest frequency (m = 0 there)",
               "Classical approach. Correct ONLY if the sample is electrically shorter than half a guided "
               "wavelength at the lowest frequency and the grid resolves the phase; it does not guarantee the "
               "physical branch."),
    "group_delay": ("Pointwise group delay (Weir 1974)",
                    "Chooses m independently at each frequency by matching Weir's non-dispersive group delay. "
                    "Noise-sensitive (numerical derivative); for comparison."),
    "prior": ("From a prior estimate of eps*mu",
              "Uses a user-supplied approximate eps_r*mu_r (e.g. datasheet) to fix the branch."),
    "principal": ("Principal value only (m = 0)",
                  "Textbook choice; valid only if beta L < pi across the whole band. Diagnostic use only."),
}


@dataclass
class AnalysisConfig:
    thickness_mm: float
    geometry: str = "coax"  # coax | waveguide
    waveguide_a_mm: float | None = None
    f_min_ghz: float | None = None
    f_max_ghz: float | None = None
    method: str = "auto"
    assume_nonmagnetic: bool = False
    branch_method: str = "anchored"
    prior_eps_mu: float | None = None
    offset1_mm: float = 0.0
    offset2_mm: float = 0.0
    holder_length_mm: float | None = None
    noise_sigma: float = 1e-3
    average_ports: bool = True
    regularization: float = 1e-2
    kk: bool = True
    uncertainty: dict | None = None  # UncertaintySpec fields

    def validate(self) -> list[str]:
        errs = []
        if not (np.isfinite(self.thickness_mm) and self.thickness_mm > 0):
            errs.append("Sample thickness must be a positive number (mm).")
        if self.geometry not in ("coax", "waveguide"):
            errs.append("Geometry must be 'coax' or 'waveguide'.")
        if self.geometry == "waveguide" and not (self.waveguide_a_mm and self.waveguide_a_mm > 0):
            errs.append("Rectangular waveguide needs the broad-wall width a (mm).")
        if self.method not in METHODS:
            errs.append(f"Unknown method {self.method!r}.")
        if self.branch_method not in BRANCH_METHODS:
            errs.append(f"Unknown branch method {self.branch_method!r}.")
        if self.branch_method == "prior" and not (self.prior_eps_mu and self.prior_eps_mu > 0):
            errs.append("Branch method 'prior' needs a positive estimate of eps_r*mu_r.")
        if self.f_min_ghz is not None and self.f_max_ghz is not None and self.f_min_ghz >= self.f_max_ghz:
            errs.append("Frequency window: minimum must be below maximum.")
        if self.holder_length_mm is not None and self.holder_length_mm < self.thickness_mm:
            errs.append("Holder length must be at least the sample thickness.")
        if not (self.noise_sigma > 0):
            errs.append("Noise level must be positive.")
        return errs

    def key(self, data_hash: str) -> str:
        blob = json.dumps({"cfg": asdict(self), "data": data_hash}, sort_keys=True, default=str)
        return hashlib.sha256(blob.encode()).hexdigest()[:16]

    def geometry_arg(self):
        if self.geometry == "waveguide":
            return {"kind": "waveguide", "a_m": self.waveguide_a_mm * 1e-3}
        return {"kind": "coax"}


def data_hash(data: TwoPortData) -> str:
    h = hashlib.sha256()
    for a in (data.frequency, data.s11, data.s21, data.s12, data.s22):
        h.update(np.ascontiguousarray(a).tobytes())
    h.update(",".join(data.measured).encode())
    return h.hexdigest()[:16]


def window(data: TwoPortData, cfg: AnalysisConfig):
    """Return (f, s11, s21, s12|None, s22|None) inside the configured frequency window."""
    f = np.asarray(data.frequency, float)
    sel = np.ones(f.size, bool)
    if cfg.f_min_ghz is not None:
        sel &= f >= cfg.f_min_ghz * 1e9 * (1 - 1e-12)
    if cfg.f_max_ghz is not None:
        sel &= f <= cfg.f_max_ghz * 1e9 * (1 + 1e-12)
    s12 = data.s12_measured()
    s22 = data.s22_measured()
    return (f[sel], data.s11[sel], data.s21[sel], None if s12 is None else s12[sel],
            None if s22 is None else s22[sel])


@dataclass
class Outcome:
    config: AnalysisConfig
    result: ExtractionResult
    assessment: Assessment
    table: pd.DataFrame
    diagnostics: pd.DataFrame
    runtime_s: float
    timestamp_utc: str
    source: str
    data_summary: dict = field(default_factory=dict)


def run_analysis(data: TwoPortData, cfg: AnalysisConfig, progress=None, data_summary=None) -> Outcome:
    """Run the core extraction for ``cfg`` on ``data`` (raises ValueError on invalid settings)."""
    errs = cfg.validate()
    if errs:
        raise ValueError(" ".join(errs))
    f, s11, s21, s12, s22 = window(data, cfg)
    if f.size < 5:
        raise ValueError("Fewer than 5 frequencies in the selected window.")
    unc = None
    if cfg.uncertainty:
        unc = UncertaintySpec(**cfg.uncertainty)
    t0 = time.perf_counter()
    res = extract_material_parameters(
        frequency=f, s11=s11, s21=s21, s12=s12, s22=s22, thickness=cfg.thickness_mm * 1e-3,
        geometry=cfg.geometry_arg(), method=cfg.method, branch_method=cfg.branch_method,
        reference_plane_offsets=(cfg.offset1_mm * 1e-3, cfg.offset2_mm * 1e-3),
        holder_length=None if cfg.holder_length_mm is None else cfg.holder_length_mm * 1e-3,
        assume_nonmagnetic=cfg.assume_nonmagnetic, noise_sigma=cfg.noise_sigma, average_ports=cfg.average_ports,
        regularization=cfg.regularization, prior_eps_mu=cfg.prior_eps_mu, uncertainty=unc, kk=cfg.kk,
        uncertainty_progress=progress)
    dt = time.perf_counter() - t0
    a = assess(res)
    return Outcome(cfg, res, a, results_table(res, a), diagnostics_table(res), dt,
                   datetime.now(timezone.utc).isoformat(timespec="seconds"), data.source, data_summary or {})


def _kk_metric(res: ExtractionResult, q: str):
    out = np.full(len(res.frequency), np.nan)
    r = res.causality.get(q) if res.causality else None
    if r is None:
        return out
    idx = np.searchsorted(res.frequency, r.frequency)
    ok = (idx < len(res.frequency))
    idx = idx[ok]
    match = np.isclose(res.frequency[idx], r.frequency[ok], rtol=1e-12, atol=0)
    out[idx[match]] = np.asarray(r.normalized_residual)[ok][match]
    return out


def results_table(res: ExtractionResult, a: Assessment) -> pd.DataFrame:
    """Per-frequency results. Loss terms are reported both as Im(x) (<= 0 for passive media,
    exp(+jwt)) and as x'' = -Im(x) (>= 0); nothing is clipped."""
    f = res.frequency
    eps_only = res.method in EPS_ONLY
    eps, mu = res.epsilon_complex, res.mu_complex
    nan = np.full(len(f), np.nan)
    sig_e = res.uncertainty["eps_loss"]["sigma_percentile"] if res.uncertainty else res.identifiability.get("eps_loss_std")
    pe = passivity_classes(-eps.imag, sig_e)
    if eps_only:
        pm = np.full(len(f), "not extracted (mu = 1 assumed)", dtype=object)
    else:
        sig_m = (res.uncertainty["mu_loss"]["sigma_percentile"] if res.uncertainty
                 else res.identifiability.get("mu_loss_std"))
        pm = passivity_classes(-mu.imag, sig_m)
    rank = {"pass": 0, "ambiguous": 1, "fail": 2, "undetermined": 3}
    worst = np.array([max((x for x in (p1, p2) if x in rank), key=lambda v: rank[v], default="pass")
                      for p1, p2 in zip(pe, pm)], dtype=object)
    with np.errstate(divide="ignore", invalid="ignore"):
        tab = {
            "frequency_Hz": f,
            "frequency_GHz": f / 1e9,
            "epsilon_real": eps.real,
            "epsilon_imag": eps.imag,
            "epsilon_loss": -eps.imag,
            "mu_real": nan if eps_only else mu.real,
            "mu_imag": nan if eps_only else mu.imag,
            "mu_loss": nan if eps_only else -mu.imag,
            "tan_delta_e": res.tan_delta_eps,
            "tan_delta_m": nan if eps_only else res.tan_delta_mu,
            "branch_index": res.selected_branch,
            "branch_margin": res.branch_margin if res.branch_margin is not None else nan,
            "passivity_eps": pe,
            "passivity_mu": pm,
            "passivity_flag": worst,
            "causality_metric_eps": _kk_metric(res, "eps"),
            "causality_metric_mu": nan if eps_only else _kk_metric(res, "mu"),
            "crlb_std_eps_real": res.identifiability.get("eps_real_std", nan),
            "crlb_std_mu_real": nan if eps_only else res.identifiability.get("mu_real_std", nan),
        }
    if res.uncertainty:
        for q in ("eps_real", "eps_loss", "mu_real", "mu_loss"):
            if eps_only and q.startswith("mu"):
                continue
            u = res.uncertainty[q]
            tab[f"{q}_p2.5"] = u["p2.5"]
            tab[f"{q}_median"] = u["median"]
            tab[f"{q}_p97.5"] = u["p97.5"]
            tab[f"{q}_sigma"] = u["sigma_percentile"]
        tab["mc_branch_change_fraction"] = res.uncertainty.get("branch_change_fraction", nan)
    tab["method"] = res.method
    tab["status"] = a.point_status
    tab["warnings"] = a.point_reasons
    return pd.DataFrame(tab)


def diagnostics_table(res: ExtractionResult) -> pd.DataFrame:
    base = res.nrw_baseline
    f = res.frequency
    d = {
        "frequency_Hz": f,
        "frequency_GHz": f / 1e9,
        "gamma_refl_re": res.reflection_coefficient.real,
        "gamma_refl_im": res.reflection_coefficient.imag,
        "gamma_refl_abs": np.abs(res.reflection_coefficient),
        "T_re": res.transmission.real,
        "T_im": res.transmission.imag,
        "T_abs": np.abs(res.transmission),
        "alpha_Np_per_m": res.propagation_constant.real,
        "beta_rad_per_m": res.propagation_constant.imag,
        "z_re": res.impedance.real,
        "z_im": res.impedance.imag,
        "branch_index": res.selected_branch,
        "branch_margin": res.branch_margin if res.branch_margin is not None else np.full(len(f), np.nan),
    }
    if base is not None:
        d["nrw_X_abs"] = np.abs(base.X)
        d["beta_L_over_pi"] = base.electrical_length_over_pi
    for k, v in res.identifiability.items():
        if isinstance(v, np.ndarray) and v.shape == f.shape:
            d[f"crlb_{k}"] = v
    for k, v in res.convergence.items():
        v = np.asarray(v)
        if v.shape == f.shape:
            d[f"solver_{k}"] = v.astype(int) if v.dtype == bool else v
    for k, v in res.flags.items():
        if np.shape(v) == f.shape:
            d[f"flag_{k}"] = np.asarray(v, int)
    return pd.DataFrame(d)


def compare_methods(data: TwoPortData, cfg: AnalysisConfig, methods=None) -> dict:
    """Run several core methods on the same data/configuration (no KK, no Monte Carlo)."""
    methods = methods or ["nrw", "iterative_eps_mu", "baker_jarvis", "nrw_eps"]
    out = {}
    for m in methods:
        c = AnalysisConfig(**{**asdict(cfg), "method": m, "kk": False, "uncertainty": None})
        try:
            out[m] = run_analysis(data, c)
        except ValueError as e:
            out[m] = e
    return out


def method_disagreement(outcomes: dict) -> pd.DataFrame:
    """Relative spread of eps and mu across the simultaneous (eps, mu) methods."""
    rows = {}
    sims = {k: v for k, v in outcomes.items() if isinstance(v, Outcome) and k not in EPS_ONLY}
    if len(sims) < 2:
        return pd.DataFrame()
    f = next(iter(sims.values())).result.frequency
    E = np.array([v.result.epsilon_complex for v in sims.values()])
    M = np.array([v.result.mu_complex for v in sims.values()])
    with np.errstate(invalid="ignore", divide="ignore"):
        rows["frequency_GHz"] = f / 1e9
        rows["eps_rel_spread"] = (np.nanmax(np.abs(E - np.nanmedian(E.real, 0) - 1j * np.nanmedian(E.imag, 0)), 0)
                                  / np.abs(np.nanmedian(E.real, 0)))
        rows["mu_rel_spread"] = (np.nanmax(np.abs(M - np.nanmedian(M.real, 0) - 1j * np.nanmedian(M.imag, 0)), 0)
                                 / np.abs(np.nanmedian(M.real, 0)))
    return pd.DataFrame(rows)
