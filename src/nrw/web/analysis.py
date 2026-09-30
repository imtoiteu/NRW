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

from ..api import ExtractionResult, extract_material_parameters, extract_multi_thickness
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
JOINT_METHOD = "multi_thickness_eps_mu"
JOINT_LABEL = ("Two-thickness joint ε & μ least squares",
               "Fits one (ε, μ) at every frequency to the S-parameters of BOTH samples of the same material "
               "(nrw.extract_multi_thickness). Each sample is first solved alone (NRW + branch + ε/μ least squares); "
               "the better-conditioned single-sample solution starts the joint fit. Because the half-wavelength "
               "points of different lengths do not coincide, the ε/μ separation stays identifiable except where "
               "both samples resonate at once (NIST TN 1355-R two-length idea; docs/benchmark_report.md §3).")

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
    # two-thickness mode (second sample of the same material, same frequency grid)
    two_thickness: bool = False
    thickness2_mm: float | None = None
    offset1_2_mm: float = 0.0
    offset2_2_mm: float = 0.0

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
        if self.two_thickness:
            if not (self.thickness2_mm is not None and np.isfinite(self.thickness2_mm) and self.thickness2_mm > 0):
                errs.append("Two-thickness mode needs a positive thickness for sample 2 (mm).")
            elif abs(self.thickness2_mm - self.thickness_mm) <= 1e-6 * max(self.thickness_mm, 1e-9):
                errs.append("The two samples must have different thicknesses.")
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


def combined_hash(data: TwoPortData, data2: TwoPortData | None, cfg: "AnalysisConfig") -> str:
    """Cache key component covering every dataset the configuration actually uses."""
    h = data_hash(data)
    if cfg.two_thickness and data2 is not None:
        h = hashlib.sha256((h + data_hash(data2)).encode()).hexdigest()[:16]
    return h


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


def grids_compatible(f1, f2, rtol=1e-9) -> str | None:
    """None if two frequency grids are identical, otherwise a human-readable reason.

    Measured S-parameters are never interpolated onto another grid: both samples must be
    measured with the same VNA sweep (or restricted to a common window of identical points).
    """
    f1, f2 = np.asarray(f1, float), np.asarray(f2, float)
    if f1.shape != f2.shape:
        return (f"the samples have different numbers of frequency points in the analysed window ({f1.size} vs "
                f"{f2.size})")
    if not np.allclose(f1, f2, rtol=rtol, atol=0):
        k = int(np.argmax(np.abs(f1 - f2)))
        return f"the frequency grids differ (e.g. {f1[k]/1e9:.9g} vs {f2[k]/1e9:.9g} GHz at point {k})"
    return None


def run_analysis(data: TwoPortData, cfg: AnalysisConfig, progress=None, data_summary=None,
                 data2: TwoPortData | None = None) -> Outcome:
    """Run the core extraction for ``cfg`` (raises ValueError on invalid settings).

    With ``cfg.two_thickness`` the second sample ``data2`` is required and the joint solver
    ``nrw.extract_multi_thickness`` is used; ``cfg.method`` is then not used.
    """
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
    if cfg.two_thickness:
        if data2 is None:
            raise ValueError("Two-thickness mode is selected but no second sample is loaded.")
        f2, t11, t21, t12, t22 = window(data2, cfg)
        why = grids_compatible(f, f2)
        if why:
            raise ValueError(f"Two-thickness mode needs both samples on the same frequency grid: {why}. Measured "
                             "data are not interpolated; re-measure with the same sweep or choose a window where "
                             "the points coincide.")
        meas = [dict(s11=s11, s21=s21, s12=s12, s22=s22, thickness=cfg.thickness_mm * 1e-3,
                     reference_plane_offsets=(cfg.offset1_mm * 1e-3, cfg.offset2_mm * 1e-3)),
                dict(s11=t11, s21=t21, s12=t12, s22=t22, thickness=cfg.thickness2_mm * 1e-3,
                     reference_plane_offsets=(cfg.offset1_2_mm * 1e-3, cfg.offset2_2_mm * 1e-3))]
        try:
            res = extract_multi_thickness(f, meas, geometry=cfg.geometry_arg(), branch_method=cfg.branch_method,
                                          noise_sigma=cfg.noise_sigma, kk=cfg.kk, uncertainty=unc,
                                          uncertainty_progress=progress, prior_eps_mu=cfg.prior_eps_mu)
        except ValueError as e:  # InputError is a ValueError
            raise ValueError(str(e)) from e
        source = f"{data.source} + {data2.source}"
    else:
        res = extract_material_parameters(
            frequency=f, s11=s11, s21=s21, s12=s12, s22=s22, thickness=cfg.thickness_mm * 1e-3,
            geometry=cfg.geometry_arg(), method=cfg.method, branch_method=cfg.branch_method,
            reference_plane_offsets=(cfg.offset1_mm * 1e-3, cfg.offset2_mm * 1e-3),
            holder_length=None if cfg.holder_length_mm is None else cfg.holder_length_mm * 1e-3,
            assume_nonmagnetic=cfg.assume_nonmagnetic, noise_sigma=cfg.noise_sigma,
            average_ports=cfg.average_ports, regularization=cfg.regularization, prior_eps_mu=cfg.prior_eps_mu,
            uncertainty=unc, kk=cfg.kk, uncertainty_progress=progress)
        source = data.source
    dt = time.perf_counter() - t0
    a = assess(res)
    return Outcome(cfg, res, a, results_table(res, a), diagnostics_table(res), dt,
                   datetime.now(timezone.utc).isoformat(timespec="seconds"), source, data_summary or {})


def single_sample_outcomes(outcome: "Outcome") -> dict:
    """Wrap the single-sample solutions of a two-thickness outcome (no recomputation)."""
    res = outcome.result
    out = {}
    if not res.sample_results:
        return out
    Ls = res.metadata.get("thicknesses_m", [])
    for i, r in enumerate(res.sample_results):
        a = assess(r)
        label = f"sample {i + 1} alone (L = {Ls[i] * 1e3:g} mm)" if i < len(Ls) else f"sample {i + 1} alone"
        out[label] = Outcome(outcome.config, r, a, results_table(r, a), diagnostics_table(r), 0.0,
                             outcome.timestamp_utc, outcome.source, outcome.data_summary)
    return out


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
    if res.sample_results:
        for i, r in enumerate(res.sample_results):
            tab[f"branch_index_sample{i + 1}"] = r.selected_branch
            tab[f"ill_conditioned_sample{i + 1}_alone"] = np.asarray(r.flags.get("singular_region",
                                                                               np.zeros(len(f), bool)), bool)
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
    """Run several single-sample core methods on ``data`` with the same settings (no KK, no
    Monte Carlo). In two-thickness mode this compares methods on sample 1 only."""
    methods = methods or ["nrw", "iterative_eps_mu", "baker_jarvis", "nrw_eps"]
    out = {}
    for m in methods:
        c = AnalysisConfig(**{**asdict(cfg), "method": m, "kk": False, "uncertainty": None, "two_thickness": False})
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
