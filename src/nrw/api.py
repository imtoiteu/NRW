"""User-facing extraction API.

    result = extract_material_parameters(frequency=f, s11=..., s21=..., thickness=2e-3,
                                         method="iterative_eps_mu", branch_method="anchored")

Pipeline (docs/final_method.md gives the full rationale):
 1. validate input (finite, strictly increasing positive frequency, positive thickness,
    above waveguide cutoff);
 2. reference planes: de-embed known offsets (or use holder length for the
    reference-plane-invariant Baker-Jarvis-type formulation);
 3. classical NRW algebra -> Gamma, T (|Gamma| <= 1 root: exact passivity condition);
 4. phase branch m(f) from the selected strategy (default "anchored": segment-wise
    unwrapping with group-delay anchoring), with ambiguity margins;
 5. material solution by the chosen method, initialised from the vetted NRW branch;
 6. diagnostics: NRW singularity flags, local identifiability (Cramer-Rao bound at the
    stated noise level), physical checks (exact / assumption / heuristic), finite-band
    KK consistency;
 7. optional Monte Carlo uncertainty through the whole pipeline.
Nothing is clipped, smoothed or NaN-replaced; problems are reported in ``warnings`` and
per-frequency ``flags``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import __version__
from .branch import select_branch
from .classical import nrw
from .conditioning import crlb
from .geometry import Geometry, geometry_from_dict
from .iterative import solve_eps_mu_pointwise, solve_eps_mu_regularized, solve_eps_only
from .kk import kk_consistency
from .physics import run_all
from .uncertainty import UncertaintySpec, monte_carlo

METHODS = ("nrw", "nrw_eps", "baker_jarvis", "iterative_eps_mu", "regularized_eps_mu", "auto")
EPS_ONLY_METHODS = ("nrw_eps", "baker_jarvis")


class InputError(ValueError):
    pass


@dataclass
class ExtractionResult:
    frequency: np.ndarray
    epsilon_complex: np.ndarray
    mu_complex: np.ndarray
    selected_branch: np.ndarray
    branch_margin: np.ndarray | None
    reflection_coefficient: np.ndarray
    transmission: np.ndarray
    propagation_constant: np.ndarray
    impedance: np.ndarray
    method: str
    branch_method: str
    convergence: dict = field(default_factory=dict)
    flags: dict = field(default_factory=dict)
    identifiability: dict = field(default_factory=dict)
    passivity: object = None  # nrw.physics.PhysicsReport
    causality: dict = field(default_factory=dict)  # {"eps": KKResult, "mu": KKResult}
    uncertainty: dict | None = None
    warnings: list = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    nrw_baseline: object = None  # nrw.classical.NRWResult

    # --- convenience views --------------------------------------------------------
    @property
    def epsilon_real(self):
        return self.epsilon_complex.real

    @property
    def epsilon_loss(self):
        return -self.epsilon_complex.imag

    @property
    def mu_real(self):
        return self.mu_complex.real

    @property
    def mu_loss(self):
        return -self.mu_complex.imag

    @property
    def tan_delta_eps(self):
        return self.epsilon_loss / self.epsilon_real

    @property
    def tan_delta_mu(self):
        return self.mu_loss / self.mu_real

    @property
    def singular_region(self):
        return self.flags.get("singular_region", np.zeros(len(self.frequency), bool))

    def table(self) -> dict:
        cols = {
            "f_hz": self.frequency,
            "eps_real": self.epsilon_real,
            "eps_loss": self.epsilon_loss,
            "mu_real": self.mu_real,
            "mu_loss": self.mu_loss,
            "tan_delta_eps": self.tan_delta_eps,
            "tan_delta_mu": self.tan_delta_mu,
            "branch_m": self.selected_branch,
            "branch_margin": self.branch_margin if self.branch_margin is not None else np.full(len(self.frequency), np.nan),
            "gamma_refl_re": self.reflection_coefficient.real,
            "gamma_refl_im": self.reflection_coefficient.imag,
            "T_re": self.transmission.real,
            "T_im": self.transmission.imag,
            "alpha_np_per_m": self.propagation_constant.real,
            "beta_rad_per_m": self.propagation_constant.imag,
            "z_re": self.impedance.real,
            "z_im": self.impedance.imag,
        }
        for k, v in self.identifiability.items():
            if isinstance(v, np.ndarray) and v.shape == self.frequency.shape:
                cols[f"crlb_{k}"] = v
        if self.uncertainty:
            for q in ("eps_real", "eps_loss", "mu_real", "mu_loss"):
                u = self.uncertainty[q]
                cols[f"{q}_p2.5"] = u["p2.5"]
                cols[f"{q}_p97.5"] = u["p97.5"]
                cols[f"{q}_sigma_robust"] = u["sigma_robust"]
            if "branch_change_fraction" in self.uncertainty:
                cols["mc_branch_change_fraction"] = self.uncertainty["branch_change_fraction"]
        for k, v in self.flags.items():
            cols[f"flag_{k}"] = np.asarray(v, int)
        return cols

    def summary(self) -> dict:
        out = {
            "method": self.method,
            "branch_method": self.branch_method,
            "n_frequencies": int(len(self.frequency)),
            "f_min_hz": float(self.frequency[0]),
            "f_max_hz": float(self.frequency[-1]),
            "branch_range": [int(self.selected_branch.min()), int(self.selected_branch.max())],
            "flag_counts": {k: int(np.count_nonzero(v)) for k, v in self.flags.items()},
            "convergence": {k: (int(np.count_nonzero(v)) if np.asarray(v).dtype == bool else None)
                            for k, v in self.convergence.items()},
            "physics": self.passivity.summary() if self.passivity is not None else None,
            "causality": {k: {"passed": r.passed, "statement": r.statement, "assumptions": r.assumptions,
                              "max_abs_normalized_residual": r.max_abs_normalized_residual}
                          for k, r in self.causality.items()},
            "warnings": self.warnings,
            "metadata": self.metadata,
        }
        if self.uncertainty:
            out["uncertainty"] = {"spec": self.uncertainty["spec"], "n_failed_draws": self.uncertainty["n_failed_draws"]}
        return out

    def save(self, outdir, plots=True):
        """Write results.csv, summary.json and (optionally) figures into ``outdir``."""
        from .io import dump_json

        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        cols = self.table()
        names = list(cols)
        data = np.column_stack([np.asarray(cols[n], float) for n in names])
        np.savetxt(outdir / "results.csv", data, delimiter=",", header=",".join(names), comments="", fmt="%.10g")
        dump_json(self.summary(), outdir / "summary.json")
        if plots:
            from .plotting import plot_result

            plot_result(self, outdir)
        return outdir


def _validate(f, s11, s21, s12, s22, thickness, geom: Geometry):
    f = np.asarray(f, dtype=float)
    if f.ndim != 1 or f.size < 3:
        raise InputError("frequency must be a 1-D array with at least 3 points")
    arrs = []
    for name, s in (("s11", s11), ("s21", s21), ("s12", s12), ("s22", s22)):
        if s is None:
            arrs.append(None)
            continue
        a = np.asarray(s, dtype=complex)
        if a.shape != f.shape:
            raise InputError(f"{name} has shape {a.shape}, expected {f.shape}")
        if not np.all(np.isfinite(a)):
            raise InputError(f"{name} contains non-finite values at indices {np.flatnonzero(~np.isfinite(a))[:10]}")
        arrs.append(a)
    if not np.all(np.isfinite(f)) or np.any(f <= 0):
        raise InputError("frequencies must be finite and positive")
    if np.any(np.diff(f) <= 0):
        raise InputError("frequency must be strictly increasing (sort the data and remove duplicates explicitly)")
    if not (np.isfinite(thickness) and thickness > 0):
        raise InputError(f"thickness must be positive and finite, got {thickness!r}")
    fc = geom.cutoff_frequency()
    if fc > 0 and np.any(f <= fc):
        raise InputError(f"frequencies at or below the empty-guide cutoff {fc/1e9:.4g} GHz are not supported")
    return f, arrs


def _core(f, s11, s21, s12, s22, L, geom, method, branch_method, opts):
    """One pass of the pipeline without uncertainty/KK. Returns a dict."""
    both = s12 is not None and s22 is not None
    a11 = (s11 + s22) / 2 if (both and opts["average_ports"]) else s11
    a21 = (s21 + s12) / 2 if (both and opts["average_ports"]) else s21
    base = nrw(f, a11, a21, L, geom, noise_floor=opts["noise_floor"])
    br = select_branch(branch_method, f, base.transmission, L, geom, noise_floor=opts["noise_floor"],
                       **({"prior_eps_mu": opts["prior_eps_mu"]} if opts["prior_eps_mu"] is not None else {}))
    r = nrw(f, a11, a21, L, geom, branch=br.m, noise_floor=opts["noise_floor"])
    # independent cross-check of the branch (different information: eps*mu continuity)
    other = "dp" if branch_method != "dp" else "anchored"
    br_check = select_branch(other, f, base.transmission, L, geom, noise_floor=opts["noise_floor"])
    conv = {}
    if method == "nrw":
        eps, mu = r.eps, r.mu
    elif method == "nrw_eps":
        eps, mu = r.eps_mu_only, np.ones_like(r.eps)
    elif method == "baker_jarvis":
        it = solve_eps_only(f, s11, s21, L, geom, r.eps_mu_only, s12=s12, s22=s22,
                            l_air_minus_l=opts["l_air_minus_l"], formulation=opts["bj_formulation"],
                            beta=opts["bj_beta"])
        eps, mu = it.eps, it.mu
        conv = {"converged": it.converged, "n_iterations": it.n_iterations, "residual_rms": it.residual_rms,
                "jacobian_cond": it.jacobian_cond}
    elif method == "iterative_eps_mu":
        it = solve_eps_mu_pointwise(f, s11, s21, L, geom, r.eps, r.mu, s12=s12, s22=s22,
                                    multistart=opts["multistart"])
        eps, mu = it.eps, it.mu
        conv = {"converged": it.converged, "n_iterations": it.n_iterations, "residual_rms": it.residual_rms,
                "jacobian_cond": it.jacobian_cond, "alternative_minimum": it.alternative_minimum}
    elif method == "regularized_eps_mu":
        it = solve_eps_mu_regularized(f, s11, s21, L, geom, r.eps, r.mu, s12=s12, s22=s22,
                                      lam=opts["regularization"], sigma=opts["noise_sigma"])
        eps, mu = it.eps, it.mu
        conv = {"converged": it.converged, "residual_rms": it.residual_rms}
    else:  # pragma: no cover - validated earlier
        raise ValueError(method)
    return {"eps": eps, "mu": mu, "nrw": r, "branch": br, "branch_check": br_check, "convergence": conv}


def extract_material_parameters(frequency, s11, s21, thickness, s12=None, s22=None, geometry="coax",
                                method="auto", branch_method="anchored", reference_plane_offsets=(0.0, 0.0),
                                holder_length=None, assume_nonmagnetic=False, noise_sigma=1e-3,
                                noise_floor=None, average_ports=True, multistart=True, bj_formulation="eq2.24+2.25", bj_beta=0.0,
                                regularization=1e-2, prior_eps_mu=None, uncertainty=None, kk=True,
                                kk_threshold=None, singular_rel_threshold=0.01,
                                uncertainty_progress=None) -> ExtractionResult:
    """Extract complex eps_r and mu_r from two-port T/R S-parameters.

    Parameters
    ----------
    frequency : array (Hz), strictly increasing.
    s11, s21 [, s12, s22] : complex arrays, normalised to the empty holder impedance.
    thickness : sample length along the propagation direction (m).
    geometry : "coax", "waveguide" (then pass a dict with "a_m") or a nrw.geometry object.
    method : "auto" | "nrw" | "nrw_eps" | "baker_jarvis" | "iterative_eps_mu" | "regularized_eps_mu".
        "auto" = "baker_jarvis" if ``assume_nonmagnetic`` else "iterative_eps_mu".
    branch_method : see nrw.branch.METHODS (default "anchored").
    reference_plane_offsets : (d1, d2) in m, distance of the calibration planes from the
        sample faces (air-filled). They are de-embedded for all methods except
        "baker_jarvis" when ``holder_length`` is given (then only L_air - L is needed).
    noise_sigma : rms complex noise assumed on each S-parameter; used for the Cramer-Rao
        identifiability bound and as the data scale of the regularised solver.
    noise_floor : |S| below which a parameter is considered noise-dominated
        (default 3 * noise_sigma).
    bj_formulation, bj_beta : NIST TN 1341 formulation for method "baker_jarvis" ("eq2.24+2.25"
        default: reference-plane invariant; "eq2.25", "eq2.38" with weight beta, "direct").
        The default was chosen from results/studies/bj_formulation.csv (docs/benchmark_report.md).
    uncertainty : None, True (default spec) or nrw.uncertainty.UncertaintySpec.
    uncertainty_progress : optional callable(done, total) called after each Monte Carlo draw.
    singular_rel_threshold : a frequency is flagged ``singular_region`` (ill-conditioned) when the
        Cramer-Rao std of eps' at ``noise_sigma`` exceeds this fraction of eps' (default 1 %,
        about 10x the level of well-conditioned points in the benchmark cases).
    """
    geom = geometry_from_dict({"kind": geometry} if isinstance(geometry, str) else geometry)
    if method not in METHODS:
        raise InputError(f"unknown method {method!r}; choose from {METHODS}")
    if method == "auto":
        method = "baker_jarvis" if assume_nonmagnetic else "iterative_eps_mu"
    f, (s11, s21, s12, s22) = _validate(frequency, s11, s21, s12, s22, thickness, geom)
    L = float(thickness)
    warn = []
    if noise_floor is None:
        noise_floor = 3 * noise_sigma

    # --- reference planes ---------------------------------------------------------
    d1, d2 = (float(x) for x in reference_plane_offsets)
    l_air_minus_l = None
    if method == "baker_jarvis" and holder_length is not None:
        l_air_minus_l = float(holder_length) - L
        if l_air_minus_l < 0:
            raise InputError("holder_length must be >= thickness")
        if d1 or d2:
            warn.append("holder_length given: reference_plane_offsets ignored (invariant formulation)")
        d1 = d2 = 0.0
    elif d1 or d2:
        g0 = geom.gamma0(f)
        r1, r2 = np.exp(-g0 * d1), np.exp(-g0 * d2)
        s11 = s11 / r1**2
        s21 = s21 / (r1 * r2)
        s12 = None if s12 is None else s12 / (r1 * r2)
        s22 = None if s22 is None else s22 / r2**2
    if holder_length is not None and method != "baker_jarvis":
        warn.append("holder_length is only used by method='baker_jarvis'; planes de-embedded with offsets instead")
    if method in EPS_ONLY_METHODS:
        warn.append(f"method '{method}' assumes mu_r = 1 (non-magnetic sample); mu is not extracted")

    opts = dict(noise_floor=noise_floor, noise_sigma=noise_sigma, average_ports=average_ports,
                multistart=multistart, bj_formulation=bj_formulation, bj_beta=bj_beta, regularization=regularization,
                prior_eps_mu=prior_eps_mu, l_air_minus_l=l_air_minus_l)
    core = _core(f, s11, s21, s12, s22, L, geom, method, branch_method, opts)
    eps, mu, base, br = core["eps"], core["mu"], core["nrw"], core["branch"]

    # --- diagnostics ----------------------------------------------------------------
    flags = {k: np.asarray(v, bool) for k, v in base.flags.items()}
    eps_only = method in EPS_ONLY_METHODS
    params = "eps" if eps_only else "eps_mu"
    ok = np.isfinite(eps) & np.isfinite(mu)
    ident = {}
    std = np.full((len(f), 2 if eps_only else 4), np.nan)
    cond = np.full(len(f), np.nan)
    if ok.any():
        std[ok], cond[ok] = crlb(f[ok], geom, eps[ok], mu[ok], L, noise_sigma, params=params,
                                 use_s22_s12=s12 is not None)
    ident["eps_real_std"], ident["eps_loss_std"] = std[:, 0], std[:, 1]
    if not eps_only:
        ident["mu_real_std"], ident["mu_loss_std"] = std[:, 2], std[:, 3]
    ident["jacobian_cond"] = cond
    ident["noise_sigma_assumed"] = noise_sigma
    rel = std[:, 0] / np.maximum(np.abs(eps.real), 1e-12)
    flags["singular_region"] = ~(rel <= singular_rel_threshold)  # includes NaN
    flags["branch_ambiguous"] = (br.margin < 0.5) if br.margin is not None else np.zeros(len(f), bool)
    flags["branch_unreliable_step"] = np.isin(np.arange(len(f)), br.unreliable_steps)
    flags["branch_disagreement"] = core["branch_check"].m != br.m
    if "converged" in core["convergence"]:
        flags["not_converged"] = ~np.asarray(core["convergence"]["converged"], bool)
    if "alternative_minimum" in core["convergence"]:
        flags["alternative_minimum"] = np.asarray(core["convergence"]["alternative_minimum"], bool)
    if eps_only:
        # eps-only methods: NRW's own S11-driven flags do not apply to the eps estimate
        flags.pop("near_half_wavelength", None)

    # --- uncertainty ------------------------------------------------------------------
    unc = None
    if uncertainty:
        spec = uncertainty if isinstance(uncertainty, UncertaintySpec) else UncertaintySpec(
            s11_mag_rel=0.005, s11_phase_deg=0.5, s21_mag_rel=0.005, s21_phase_deg=0.5,
            thickness_m=5e-6, ref_plane_m=5e-6, n_draws=100)
        mc_opts = dict(opts, multistart=False)

        def run(ff, a11, a21, a12, a22, LL):
            c = _core(ff, a11, a21, a12, a22, LL, geom, method, branch_method, mc_opts)
            return c["eps"], c["mu"], c["branch"].m

        unc = monte_carlo(run, f, s11, s21, s12 if s12 is not None else s21, s22 if s22 is not None else s11,
                          L, geom, spec, nominal_branch=br.m, progress=uncertainty_progress)
        flags["mc_branch_unstable"] = unc["branch_change_fraction"] > 0.05

    eps_sig = unc["eps_loss"]["sigma_percentile"] if unc else None
    mu_sig = unc["mu_loss"]["sigma_percentile"] if unc else None
    phys = run_all(f, s11, s21, s12, s22, eps, mu, base.propagation_constant, base.impedance,
                   eps_sig, mu_sig, eps_mu_fixed=eps_only, noise_sigma=noise_sigma, length_m=L)

    caus = {}
    if kk:
        for q, val, key in (("eps", eps, "eps_real"), ("mu", mu, "mu_real")):
            if q == "mu" and eps_only:
                continue
            good = np.isfinite(val) & ~flags["singular_region"]
            if good.sum() < 5:
                warn.append(f"KK diagnostic for {q} skipped: fewer than 5 well-conditioned points")
                continue
            sig = sig_i = None
            lkey = "eps_loss" if q == "eps" else "mu_loss"
            col = 0 if q == "eps" else 2
            if unc:
                sig = np.maximum(unc[key]["sigma_percentile"], 1e-6 * np.abs(val.real))[good]
                sig_i = np.maximum(unc[lkey]["sigma_percentile"], 1e-6 * np.abs(val.real))[good]
            elif np.isfinite(std[good, col]).all():
                sig = np.maximum(std[good, col], 1e-4 * np.abs(val.real[good]))
                sig_i = np.maximum(std[good, col + 1], 1e-4 * np.abs(val.real[good]))
            try:
                caus[q] = kk_consistency(f[good], val[good], sigma=sig, sigma_imag=sig_i, threshold=kk_threshold,
                                         quantity=q)
            except ValueError as e:
                warn.append(f"KK diagnostic for {q} failed to run: {e}")

    # --- warnings ---------------------------------------------------------------------
    nsing = int(flags["singular_region"].sum())
    if nsing:
        warn.append(f"{nsing} frequencies in ill-conditioned (half-wavelength-type) regions: Cramer-Rao "
                    f"std(eps') > {singular_rel_threshold:.0%} of eps' at noise {noise_sigma:g}; treat these values "
                    "as unreliable")
    if flags["branch_ambiguous"].any():
        warn.append(f"phase-branch ambiguity margin < 0.5 at {int(flags['branch_ambiguous'].sum())} frequencies")
    if flags["branch_disagreement"].any():
        warn.append(f"branch strategies '{branch_method}' and '{core['branch_check'].method}' disagree at "
                    f"{int(flags['branch_disagreement'].sum())} frequencies: the phase branch there is not "
                    "determined by the data; supply a prior (branch_method='prior'), a second thickness, or "
                    "denser/lower-frequency data")
    if br.unreliable_steps.size:
        warn.append(f"{br.unreliable_steps.size} unreliable phase steps (|dpsi| > 0.75 pi or |T| at noise floor); "
                    f"branch anchored independently in {len(br.segments) or 1} segment(s)")
    if phys.n_exact_violations:
        bad = [c.name for c in phys.checks if c.status == "exact" and c.n_violations]
        warn.append(f"exact physical constraints violated: {bad}")
    for c in phys.checks:
        if c.status == "assumption" and c.n_violations:
            warn.append(f"model assumption '{c.name}' violated at {c.n_violations} frequencies ({c.description})")
    if eps_only:
        good = ~flags.get("s11_below_noise", np.zeros(len(f), bool)) & np.isfinite(base.mu)
        if good.any():
            dev = np.nanmedian(np.abs(base.mu[good] - 1))
            if dev > 0.05:
                warn.append(f"NRW indicates |mu - 1| ~ {dev:.2g} (median): the mu = 1 assumption of '{method}' "
                            "may be invalid")
    if geom.kind == "coax":
        n2 = np.nanmax(np.abs((eps * mu).real)) if np.isfinite(eps * mu).any() else 1.0
        fte = geom.te11_cutoff_hz(max(n2, 1.0))
        if f[-1] > fte:
            warn.append(f"the filled line can support TE11 above ~{fte/1e9:.3g} GHz (eps*mu ~ {n2:.3g}); it is "
                        "excited only by sample asymmetry/eccentricity or gaps, but resonant artefacts there "
                        "invalidate the TEM model")

    return ExtractionResult(
        frequency=f,
        epsilon_complex=eps,
        mu_complex=mu,
        selected_branch=br.m,
        branch_margin=br.margin,
        reflection_coefficient=base.reflection_coefficient,
        transmission=base.transmission,
        propagation_constant=base.propagation_constant,
        impedance=base.impedance,
        method=method,
        branch_method=branch_method,
        convergence=core["convergence"],
        flags=flags,
        identifiability=ident,
        passivity=phys,
        causality=caus,
        uncertainty=unc,
        warnings=warn,
        metadata={"nrw_version": __version__, "thickness_m": L, "geometry": geom.to_dict(),
                  "reference_plane_offsets_m": [d1, d2], "holder_length_m": holder_length,
                  "noise_sigma": noise_sigma, "noise_floor": noise_floor,
                  "singular_rel_threshold": singular_rel_threshold, "average_ports": average_ports,
                  "regularization": regularization if method == "regularized_eps_mu" else None,
                  "bj_formulation": bj_formulation if method == "baker_jarvis" else None,
                  "bj_beta": bj_beta if method == "baker_jarvis" else None,
                  "conventions": "exp(+jwt); eps = eps' - j eps''; passive => eps'', mu'' >= 0"},
        nrw_baseline=base,
    )


def extract_multi_thickness(frequency, measurements, geometry="coax", branch_method="anchored", noise_sigma=1e-3,
                            singular_rel_threshold=0.01, kk=True):
    """Simultaneous (eps, mu) from two or more samples of different length (same material).

    measurements: list of dicts {"s11", "s21", optional "s12", "s22", "thickness" (m),
    optional "reference_plane_offsets"}. Each sample is first processed by the single-sample
    pipeline (NRW + branch); at every frequency the single-sample solution with the smaller
    Cramer-Rao std of eps' initialises a joint least-squares fit to all samples
    (nrw.iterative.solve_eps_mu_multithickness). Identifiability is lost only where all
    samples are near a half-wavelength condition simultaneously; the joint Cramer-Rao bound
    (sum of the Fisher informations) is reported and flagged as in the single-sample case.
    """
    from .conditioning import real_jacobian
    from .iterative import solve_eps_mu_multithickness

    if len(measurements) < 2:
        raise InputError("extract_multi_thickness needs at least two measurements")
    Ls = [float(m["thickness"]) for m in measurements]
    if len(set(Ls)) < len(Ls):
        raise InputError("sample lengths must differ")
    singles = [extract_material_parameters(frequency, m["s11"], m["s21"], m["thickness"], s12=m.get("s12"),
                                           s22=m.get("s22"), geometry=geometry, method="iterative_eps_mu",
                                           branch_method=branch_method, noise_sigma=noise_sigma,
                                           reference_plane_offsets=m.get("reference_plane_offsets", (0.0, 0.0)),
                                           multistart=False, kk=False)
               for m in measurements]
    f = singles[0].frequency
    geom = geometry_from_dict({"kind": geometry} if isinstance(geometry, str) else geometry)
    stds = np.stack([np.where(np.isfinite(r.identifiability["eps_real_std"]), r.identifiability["eps_real_std"],
                              np.inf) for r in singles])
    pick = np.argmin(stds, axis=0)
    e0 = np.choose(pick, [r.epsilon_complex for r in singles])
    u0 = np.choose(pick, [r.mu_complex for r in singles])
    samples = []
    for m, r in zip(measurements, singles):
        d1, d2 = m.get("reference_plane_offsets", (0.0, 0.0))
        samples.append({"s11": np.asarray(m["s11"], complex), "s21": np.asarray(m["s21"], complex),
                        "s12": None if m.get("s12") is None else np.asarray(m["s12"], complex),
                        "s22": None if m.get("s22") is None else np.asarray(m["s22"], complex),
                        "length_m": float(m["thickness"]), "d1": d1, "d2": d2})
    it = solve_eps_mu_multithickness(f, samples, geom, e0, u0)
    F = np.zeros((len(f), 4, 4))
    ok = np.isfinite(it.eps) & np.isfinite(it.mu)
    for smp in samples:  # Fisher information adds over independent samples
        J = real_jacobian(f[ok], geom, it.eps[ok], it.mu[ok], smp["length_m"], use_s22_s12=smp["s12"] is not None)
        F[ok] += np.einsum("nop,noq->npq", J, J) / (noise_sigma**2 / 2)
    std = np.full((len(f), 4), np.nan)
    with np.errstate(all="ignore"):
        for k in np.flatnonzero(ok):
            try:
                std[k] = np.sqrt(np.clip(np.diag(np.linalg.inv(F[k])), 0, None))
            except np.linalg.LinAlgError:
                pass
    flags = {"singular_region": ~(std[:, 0] / np.maximum(np.abs(it.eps.real), 1e-12) <= singular_rel_threshold),
             "not_converged": ~it.converged}
    for i, r in enumerate(singles):
        flags[f"sample{i}_branch_disagreement"] = r.flags["branch_disagreement"]
    phys = run_all(f, samples[0]["s11"], samples[0]["s21"], samples[0]["s12"], samples[0]["s22"], it.eps, it.mu,
                   noise_sigma=noise_sigma)
    caus = {}
    if kk:
        good = ok & ~flags["singular_region"]
        if good.sum() >= 5:
            caus["eps"] = kk_consistency(f[good], it.eps[good], sigma=std[good, 0], sigma_imag=std[good, 1],
                                         quantity="eps")
            caus["mu"] = kk_consistency(f[good], it.mu[good], sigma=std[good, 2], sigma_imag=std[good, 3],
                                        quantity="mu")
    warn = [w for r in singles for w in r.warnings if "branch" in w]
    if flags["singular_region"].any():
        warn.append(f"{int(flags['singular_region'].sum())} frequencies remain ill-conditioned for all samples")
    return ExtractionResult(
        frequency=f, epsilon_complex=it.eps, mu_complex=it.mu, selected_branch=singles[0].selected_branch,
        branch_margin=singles[0].branch_margin, reflection_coefficient=singles[0].reflection_coefficient,
        transmission=singles[0].transmission, propagation_constant=singles[0].propagation_constant,
        impedance=singles[0].impedance, method="multi_thickness_eps_mu", branch_method=branch_method,
        convergence={"converged": it.converged, "residual_rms": it.residual_rms, "jacobian_cond": it.jacobian_cond},
        flags=flags, identifiability={"eps_real_std": std[:, 0], "eps_loss_std": std[:, 1], "mu_real_std": std[:, 2],
                                      "mu_loss_std": std[:, 3], "noise_sigma_assumed": noise_sigma},
        passivity=phys, causality=caus, warnings=warn,
        metadata={"nrw_version": __version__, "thicknesses_m": Ls, "geometry": geom.to_dict(),
                  "noise_sigma": noise_sigma,
                  "conventions": "exp(+jwt); eps = eps' - j eps''; passive => eps'', mu'' >= 0"},
        nrw_baseline=singles[0].nrw_baseline)
