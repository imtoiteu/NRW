"""Traffic-light assessment of an extraction (GREEN / YELLOW / RED) from documented rules.

The assessment is a *summary of diagnostics*, not an accuracy estimate. GREEN means that
no documented diagnostic fired; it does not prove the extracted values are correct (e.g.
uncorrected thickness, air-gap or calibration errors are invisible to every diagnostic here;
docs/benchmark_report.md section 6).

Per-frequency status
--------------------
``unreliable`` if any of: non-finite/blow-up value; ill-conditioned (Cramer-Rao) point;
branch trackers disagree; Monte Carlo branch unstable; |S21| at the noise floor;
solver not converged; significant passivity violation (loss < -k sigma) or negative
attenuation beyond noise.
``attention`` if any of: near a half-wavelength point; |S11| at the noise floor; weak branch
margin or unreliable phase step; alternative least-squares minimum; ambiguous Gamma root;
negative loss within uncertainty.
``ok`` otherwise.

Overall level (fractions refer to the analysed frequencies)
----------------------------------------------------------
RED (the data/configuration do not support a trustworthy extraction) if any:
  R1 non-finite or blown-up values at > 10 %;
  R2 measured network violates passivity (sigma_max(S) > 1 beyond noise) at > 10 %;
  R3 phase branch not determined (tracker disagreement or MC branch instability) at > 10 %;
  R4 |S21| below the noise floor at > 25 %;
  R5 ill-conditioned (Cramer-Rao) at > 50 %;
  R6 a mu = 1 method was used although NRW indicates |mu - 1| > 0.05;
  R7 significant material passivity violation at > 10 %.
YELLOW if not RED and any: an ``unreliable`` or ``attention`` point exists; the KK
diagnostic failed or could not be run; a model assumption (reciprocity, port symmetry) is
violated.
GREEN otherwise.

INFO notes (shown, but not affecting the level): I1 the filled coaxial line can support TE11 in
part of the band (depends on eps*mu and line size, not on the data; excited only by sample
asymmetry/gaps, which no diagnostic here can detect).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

RULES_VERSION = "1.0"
LIMITS = {"R1_nonfinite": 0.10, "R2_network_passivity": 0.10, "R3_branch": 0.10, "R4_s21_floor": 0.25,
          "R5_ill_conditioned": 0.50, "R7_material_passivity": 0.10}

UNRELIABLE_FLAGS = {
    "nonfinite": "non-finite value",
    "blowup": "numerical blow-up (|eps| or |mu| > 1e4)",
    "singular_region": "Potential half-wavelength instability (ill-conditioned, Cramer-Rao)",
    "branch_disagreement": "Branch selection uncertain (independent trackers disagree)",
    "mc_branch_unstable": "Branch selection uncertain (Monte Carlo draws change branch)",
    "s21_below_noise": "Transmission magnitude too small for reliable inversion",
    "not_converged": "iterative solver did not converge",
}
ATTENTION_FLAGS = {
    "near_half_wavelength": "near a half-wavelength point (beta L ~ k pi)",
    "s11_below_noise": "|S11| at the noise floor (reflection phase unreliable)",
    "branch_ambiguous": "Branch selection uncertain (weak ambiguity margin)",
    "branch_unreliable_step": "unreliable phase step (branch re-anchored)",
    "alternative_minimum": "alternative least-squares minimum (eps/mu trade-off)",
    "gamma_root_ambiguous": "|Gamma| ~ 1: reflection root ambiguous",
    "T_gain": "|T| > 1 (apparent gain)",
}


@dataclass
class Assessment:
    level: str  # GREEN | YELLOW | RED
    reasons: list = field(default_factory=list)  # dicts: severity, rule, message, n_points, fraction
    point_status: np.ndarray | None = None  # "ok" | "attention" | "unreliable" per frequency
    point_reasons: list = field(default_factory=list)  # per-frequency list of short reasons
    rules_version: str = RULES_VERSION

    def to_dict(self):
        return {"level": self.level, "rules_version": self.rules_version, "reasons": self.reasons,
                "n_ok": int(np.sum(self.point_status == "ok")),
                "n_attention": int(np.sum(self.point_status == "attention")),
                "n_unreliable": int(np.sum(self.point_status == "unreliable"))}


def passivity_classes(loss, sigma=None, k=2.0):
    """Classify a loss term (eps'' or mu'', exp(+jwt): passive => >= 0) per frequency.

    'pass' if loss >= 0; 'ambiguous' if -k sigma <= loss < 0 (consistent with passivity
    within uncertainty); 'fail' if loss < -k sigma; 'undetermined' if loss is not finite.
    Without sigma every negative value is 'fail'. The values themselves are never altered.
    """
    loss = np.asarray(loss, float)
    sig = np.zeros_like(loss) if sigma is None else np.nan_to_num(np.asarray(sigma, float), nan=0.0)
    out = np.where(loss >= 0, "pass", np.where(loss >= -k * sig, "ambiguous", "fail")).astype(object)
    out[~np.isfinite(loss)] = "undetermined"
    return out


def _material_sigma(res, which):
    if res.uncertainty:
        return res.uncertainty[f"{which}_loss"]["sigma_percentile"]
    return res.identifiability.get(f"{which}_loss_std")


def assess(res) -> Assessment:
    """Assess an ``ExtractionResult`` with the rules documented in this module."""
    n = len(res.frequency)
    flags = {k: np.asarray(v, bool) for k, v in res.flags.items() if np.shape(v) == (n,)}
    eps_only = res.method in ("nrw_eps", "baker_jarvis")
    unrel = np.zeros(n, bool)
    att = np.zeros(n, bool)
    point_reasons = [[] for _ in range(n)]
    for name, msg in UNRELIABLE_FLAGS.items():
        if name in flags:
            unrel |= flags[name]
            for i in np.flatnonzero(flags[name]):
                point_reasons[i].append(name)
    for name, msg in ATTENTION_FLAGS.items():
        if name in flags and not (eps_only and name in ("near_half_wavelength", "s11_below_noise")):
            att |= flags[name]
            for i in np.flatnonzero(flags[name]):
                point_reasons[i].append(name)
    # passivity of the material (never altered, only classified)
    pv_fail = np.zeros(n, bool)
    pv_amb = np.zeros(n, bool)
    for which, val in (("eps", res.epsilon_loss), ("mu", res.mu_loss)):
        if which == "mu" and eps_only:
            continue
        cls = passivity_classes(val, _material_sigma(res, which))
        pv_fail |= cls == "fail"
        pv_amb |= cls == "ambiguous"
        for i in np.flatnonzero(cls == "fail"):
            point_reasons[i].append(f"{which}_passivity_fail")
        for i in np.flatnonzero(cls == "ambiguous"):
            point_reasons[i].append(f"{which}_passivity_ambiguous")
    unrel |= pv_fail
    att |= pv_amb
    status = np.where(unrel, "unreliable", np.where(att, "attention", "ok")).astype(object)

    reasons = []

    def add(sev, rule, msg, mask=None):
        cnt = int(np.count_nonzero(mask)) if mask is not None else None
        reasons.append({"severity": sev, "rule": rule, "message": msg, "n_points": cnt,
                        "fraction": (cnt / n) if (cnt is not None and n) else None})

    frac = lambda m: float(np.mean(m)) if n else 0.0  # noqa: E731
    nonfin = flags.get("nonfinite", np.zeros(n, bool)) | flags.get("blowup", np.zeros(n, bool))
    if frac(nonfin) > LIMITS["R1_nonfinite"]:
        add("RED", "R1", "non-finite or blown-up values at more than 10 % of frequencies", nonfin)
    net = np.zeros(n, bool)
    if res.passivity is not None:
        try:
            net = np.asarray(res.passivity["network_passivity"].violated, bool)
        except KeyError:
            pass
    if frac(net) > LIMITS["R2_network_passivity"]:
        add("RED", "R2", "measured S-matrix is not passive (sigma_max(S) > 1 beyond noise) at more than 10 % of "
                         "frequencies: check calibration, normalisation and port assignment", net)
    br = flags.get("branch_disagreement", np.zeros(n, bool)) | flags.get("mc_branch_unstable", np.zeros(n, bool))
    if frac(br) > LIMITS["R3_branch"]:
        add("RED", "R3", "Branch selection uncertain at more than 10 % of frequencies: the phase branch is not "
                         "determined by the data (use a prior, a second thickness, or denser/lower-frequency data)", br)
    s21 = flags.get("s21_below_noise", np.zeros(n, bool))
    if frac(s21) > LIMITS["R4_s21_floor"]:
        add("RED", "R4", "Transmission magnitude too small for reliable inversion at more than 25 % of frequencies "
                         "(sample too thick/lossy for this dynamic range)", s21)
    sing = flags.get("singular_region", np.zeros(n, bool))
    if frac(sing) > LIMITS["R5_ill_conditioned"]:
        add("RED", "R5", "ill-conditioned (Cramer-Rao) at more than half of the frequencies", sing)
    if any("mu = 1 assumption" in w for w in res.warnings):
        add("RED", "R6", "a mu = 1 method was used but NRW indicates a magnetic sample: the permittivity is biased "
                         "(use a simultaneous eps/mu method)")
    if frac(pv_fail) > LIMITS["R7_material_passivity"]:
        add("RED", "R7", "significant passivity violations (loss < -2 sigma) at more than 10 % of frequencies", pv_fail)

    if unrel.any():
        add("YELLOW", "Y1", "frequencies flagged unreliable (see per-frequency status)", unrel)
    if att.any():
        add("YELLOW", "Y2", "frequencies requiring attention (see per-frequency status)", att)
    for q in ("eps", "mu"):
        if q == "mu" and eps_only:
            continue
        r = res.causality.get(q) if res.causality else None
        if r is None:
            add("YELLOW", "Y3", f"finite-band KK diagnostic for {q} not available (not run or too few "
                                "well-conditioned points)")
        elif not r.passed:
            add("YELLOW", "Y3", f"finite-band KK consistency diagnostic for {q} FAILED "
                                f"(max normalised residual {r.max_abs_normalized_residual:.3g} > {r.threshold:.3g})")
    if res.passivity is not None:
        for c in res.passivity.checks:
            if c.status == "assumption" and c.n_violations:
                add("YELLOW", "Y4", f"model assumption violated: {c.name} ({c.description})",
                    np.asarray(c.violated, bool))
    if any("TE11" in w for w in res.warnings):
        add("INFO", "I1", "the filled coaxial line can support TE11 in part of the band; it is excited only by "
                          "sample asymmetry, eccentricity or gaps, which the data diagnostics cannot detect - look "
                          "for narrow resonant artefacts")
    sev = {r["severity"] for r in reasons}
    level = "RED" if "RED" in sev else ("YELLOW" if "YELLOW" in sev else "GREEN")
    return Assessment(level, reasons, status, [",".join(p) for p in point_reasons])
