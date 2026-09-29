# Final recommended extraction workflow

This document specifies the workflow implemented by `nrw.extract_material_parameters`
(single sample) and `nrw.extract_multi_thickness` (two or more lengths) so that it can be
reproduced independently. Equation numbers [Rn] refer to `docs/references.md`, and theory
sections (§) to `docs/theory_review.md`. Conventions are fixed in `docs/conventions.md`:
exp(+jωt), ε = ε′ − jε″, fields ∝ exp(−γz), T = exp(−γL), z = μγ0/γ, Γ = (z − 1)/(z + 1).

The evidence behind each choice is in `docs/benchmark_report.md` (BR) and in the study
outputs under `results/studies/`.

---

## 0. Decision summary

| sample | recommended method | why (evidence) |
|---|---|---|
| known **non-magnetic** (μ<sub>r</sub> = 1 can be justified independently) | NIST ε-only iterative method [R3, R4], invariant form TN 1341 eqs. (2.24)+(2.25) (`method="baker_jarvis"`) | Reaches its Cramér–Rao bound across half-wave points (BR §3). The invariant form is immune to an unknown split of the reference planes (BR §5). |
| magnetic or **unknown** | simultaneous (ε, μ) least squares on all four S-parameters, started from NRW with a vetted branch (`method="iterative_eps_mu"`), **with CRLB-based flagging** | No single-sample method beats the (ε, μ) CRLB. NRW/LSQ attain it (BR §3), so the only honest treatment of half-wave regions is to flag them and inflate the uncertainty. |
| magnetic, low loss, half-wave regions matter | **two or more sample lengths** (`extract_multi_thickness`) | Restores identifiability except where all samples resonate at once (BR §3; TN 1355-R [R6]). |
| any | branch: group-delay-anchored segment unwrapping, cross-checked by dynamic-programming tracking | Survives thick samples, a high start frequency and gaps where plain unwrapping fails over the whole band (BR §4). |

Do **not** use `baker_jarvis`/`nrw_eps` on magnetic samples: they return the effective
quantity εμ/1, biased by a factor μ (BR §3, magnetic case).

---

## 1. Measurement import

* `.s2p` Touchstone v1 (RI/MA/DB; Hz…GHz). The 2-port order is `S11 S21 S12 S22`. Also CSV.
  Reader: `nrw.io.read_measurement` (builtin; scikit-rf optional and cross-checked in tests).
* The reference resistance must equal the empty-holder impedance. It is reported, not
  renormalised.

## 2. Validation (`nrw.api._validate`; errors are raised, never repaired)

* 1-D frequency array with ≥ 3 points, finite, **strictly increasing**, > 0. Non-monotonic or
  duplicated frequencies are rejected, not sorted silently.
* S-parameters finite and of matching shape; thickness finite and > 0.
* Waveguide: every frequency above the empty-guide cutoff.
* Coax: a warning when the filled line could support TE11 (f > c/(π(a+b)√(εμ))). That mode is
  excited only by asymmetry, but it invalidates the TEM model if excited.
* Network checks before extraction (`nrw.physics.check_network`):
  * exact: largest singular value of S ≤ 1, and 1 − |S11|² − |S21|² ≥ 0;
  * model assumptions: reciprocity |S21 − S12|/|S21| ≤ 5 %, and port symmetry.

## 3. Reference-plane handling

* Known offsets d1, d2 (planes in air, at distance d_i from the sample faces) are de-embedded
  exactly: S11/R1², S21/(R1R2), S12/(R1R2), S22/R2², with R_i = exp(−γ0 d_i).
* For **non-magnetic** samples with uncertain individual plane positions, give the holder
  length L_air instead. The NIST invariant equations (2.24)+(2.25) [R4] then depend only on
  L_air − L. In the formulation study (BR §5), displacing the planes asymmetrically by
  3 mm/1 mm left this form unchanged, while four-S-parameter direct least squares degraded
  about 10×.
* An uncorrected plane error is **not** detectable from the extracted values by any exact
  test. The KK diagnostic detected a 50 µm port-1 error in some cases only (BR §6), so
  plane positions must come from the calibration.

## 4. Preprocessing

* No smoothing, filtering, clipping or NaN replacement.
* If S12 and S22 are present, the NRW *initialisation* uses the port averages
  (S11 + S22)/2 and (S21 + S12)/2 (a reciprocal, symmetric sample is assumed and checked in
  step 2). The iterative solvers then fit all four S-parameters individually.
* Frequencies may be non-uniform or have gaps. Branch tracking accounts for them (step 7).

## 5. Initial extraction (classical NRW, `nrw.classical.nrw`)

X = (S11² − S21² + 1)/(2S11); Γ = X ∓ √(X² − 1), taking the root with |Γ| ≤ 1. Because the
two roots satisfy Γ·Γ′ = 1, this choice is the exact passivity condition Re z ≥ 0.
Then T = (S11 + S21 − Γ)/(1 − (S11 + S21)Γ). Intermediate quantities (X, both Γ roots, T,
ln T, βL/π) are returned.

Flags:
* `s11_below_noise`: |S11| < noise floor (default 3σ);
* `near_half_wavelength`: |βL/π − k| < 0.05 with k ≥ 1;
* `gamma_root_ambiguous`: both roots have |Γ| ≈ 1;
* `s21_below_noise`;
* `T_gain`: |T| > 1;
* `nonfinite`;
* `blowup`: |ε| or |μ| > 10⁴.

## 6. Phase handling

γL = −ln T = −ln|T| + j(−arg T + 2πm). The principal log is taken, and the integer branch m
is determined separately (step 7) and **reported per frequency** (`selected_branch`). An
iterative solver cannot resolve m: every branch reproduces the data exactly
(`tests/test_iterative.py::test_iterative_eps_mu_branch_inherited`).

## 7. Branch determination (`nrw.branch`)

1. **Segmenting.** The band is split into segments at "unreliable steps", which are any of:
   * a wrapped phase step |Δψ| > 0.75π;
   * |T| below the noise floor;
   * a gap where the step predicted from the band-median group delay, 2πτΔf, exceeds 0.75π.
2. **Unwrapping within segments.** Within each segment the phase is unwrapped, which
   determines the *changes* of m.
3. **Anchoring each segment.** The single unknown integer offset of each segment is the one
   minimising the |T|-weighted median of |τ_calc(m) − τ_meas| over the segment.
   * τ_calc is Weir's non-dispersive group delay [R2], written exactly for lossy media as
     τ = −L·Im(k0·εμ/(c·γ)).
   * τ_meas is the finite-difference group delay of the locally unwrapped phase of T.
   * Using the band median, rather than choosing m pointwise as in [R2], makes the choice
     robust to noise (BR §4: pointwise selection had 12–48 % wrong points at σ = 1e-2 to 3e-2).
4. **Cross-check.** A dynamic-programming (Viterbi [R44]) tracker is run independently. It
   uses continuity of εμ between neighbouring frequencies plus the group-delay term, and does
   not need |Δψ| < π. Frequencies where the two trackers disagree are flagged
   `branch_disagreement` and produce a warning.
   * On aliased grids (step > π) only the DP tracker is correct (BR §4). In that case,
     re-run with `branch_method="dp"`, or supply a prior (`branch_method="prior"`), a
     second thickness, or denser data.
5. **Margin.** Each segment reports an ambiguity margin (1 − best/second-best cost). A
   margin below 0.5 is flagged `branch_ambiguous`.

Known limitation (BR §4): with strong in-band dispersion, where |S21| reaches the noise
floor over a resonance (ε′ even negative), **no phase-only strategy** is reliable after
the resonance. Only prior knowledge or a second thickness helps.

## 8. Robust material solution (`nrw.iterative`)

* **Non-magnetic** (`baker_jarvis`): complex damped-Newton / Levenberg–Marquardt [R42] on
  TN 1341 eqs. (2.24)+(2.25) [R4]. Newton on real and imaginary parts is the published
  solver; damping only enlarges the convergence basin.
  * Unknown: ε (2 real); μ = 1 fixed.
  * Start value: the Boughriet/NRW μ = 1 value ε = (k_c² − γ²)/k0² on the vetted branch.
  * Convergence: relative step < 1e-12, relative cost decrease < 1e-15, or cost < 1e-30;
    at most 100 iterations. Convergence, iteration count, residual rms and Jacobian
    condition number are reported per frequency.
  * Alternatives: eq. (2.38) with weight β, eq. (2.25) alone, or direct least squares
    (`bj_formulation`).
* **Magnetic/unknown** (`iterative_eps_mu`): the same solver on S11, S21, S12, S22 with 4 real
  unknowns, started from NRW (ε, μ).
  * A multistart probes the ε/μ trade-off valley: εμ fixed, μ scaled by 1.5 and 1/1.5. A
    different solution with residual within 10× of the best is flagged
    `alternative_minimum`.
  * For a symmetric sample with averaged ports NRW is already the exact solution. The LSQ
    step matters when S11 ≠ S22 or S21 ≠ S12, i.e. with unequal errors at the two ports.
* **Multi-thickness** (`extract_multi_thickness`): the joint LSQ over all samples, started at
  each frequency from the single-sample solution with the smaller CRLB.
* **Regularised** (`regularized_eps_mu`): Tikhonov second-difference smoothing over frequency.
  This is a *model assumption* that introduces bias at sharp features. It is not in the
  default workflow, and the weight must be reported if used.

## 9. Singular-region treatment

The Cramér–Rao standard deviations of (ε′, ε″, μ′, μ″) are computed at the solution for the
stated noise level (`noise_sigma`) [R40]. Frequencies where CRLB σ(ε′)/ε′ > 1 % (default)
are flagged `singular_region`. The values are **reported, not removed or interpolated**,
and are excluded from the KK diagnostic.

The threshold is about 10× the level of well-conditioned points in the benchmark. Near
half-wave points the (ε, μ) CRLB rises by more than 20× (BR §3), which is why flagging, not
reformulation, is the remedy.

## 10. Physical validation (`nrw.physics`)

| check | status | criterion |
|---|---|---|
| ε″ ≥ 0, μ″ ≥ 0 | exact (passive, isotropic, linear medium; §6 of theory) | violation only if below −2σ (MC σ when available) |
| Re γ ≥ 0 (\|T\| ≤ 1) | exact | |
| Re z ≥ 0 (\|Γ\| ≤ 1) | exact | |
| σ_max(S) ≤ 1, absorbed power ≥ 0 | exact (passive network) | tolerance 1e-3 |
| S21 = S12 | assumption (reciprocal sample, stable VNA) | 5 % |
| S11 = S22 | assumption (homogeneous sample, symmetric planes) | 20 % |
| \|n\| ≤ 50, ε′ ≥ 1, μ′ > 0, continuity | heuristic | flag only |

## 11. Causal-consistency diagnostic (`nrw.kk.kk_consistency`)

The real part is predicted from the loss by the singly-subtractive finite-band KK integral.
The integral is computed exactly for piecewise-linear data using an analytic kernel matrix,
and is anchored at the band centre so that ε∞ cancels. To it is added the best
**passive out-of-band continuation**: a non-negative, piecewise-linear loss on 24
log-spaced nodes below 2 GHz and 24 above 18 GHz, joined continuously to the measured loss
at the band edges, with coefficients fitted by NNLS [R41].

Residuals are normalised by σ_tot² = σ′² + σ″² + σ_num², where σ_num is the self-estimated
quadrature error from a half-resolution recomputation. "Pass" means
max |residual|/σ_tot ≤ Φ⁻¹(1 − α/2N), a Bonferroni threshold [R43] with α = 0.01.

* **Power:** a wrong branch fails in 100 % of the tested runs (BR §6).
* **Calibration:** ground truth passes in all cases, and correct noisy extractions pass in
  38/40 runs.
* **Blind spots:** in the magnetic (ε, μ) cases, 1 % thickness errors, 10 µm air gaps and a
  51 Ω holder all **pass**; so do air gaps in the dielectric cases, a constant offset in ε′,
  and any real-part trend that a passive absorption above 18 GHz could produce.

The software's pass statement reads "…passed within the available bandwidth and
assumptions… does NOT establish that the material strictly satisfies the global
Kramers–Kronig relations".

## 12. Uncertainty estimation (`nrw.uncertainty`)

Monte Carlo (JCGM 101 [R33]) **through the whole pipeline, including branch selection**.

* Inputs perturbed:
  * |S11|, |S21| (relative), and the S11/S21 phases (degrees), independently for the two
    port pairs;
  * optional additive complex noise;
  * the sample length;
  * each reference plane independently.
* Defaults for the CLI: 0.5 %, 0.5°, 5 µm, 5 µm, 200 draws, seed 12345.
* Output per frequency: median, 2.5/97.5 and 15.87/84.13 percentiles, MAD-σ, the fraction of
  draws whose branch differs from the nominal one (`mc_branch_unstable` if > 5 %), and the
  fraction of non-finite draws.
* The CRLB is reported alongside as a lower bound for white noise.
* Coverage of the MC 95 % intervals is quantified in BR §7.

## 13. Result export

`ExtractionResult.save(dir)` writes:
* `results.csv`: every per-frequency quantity, including ε′, ε″, μ′, μ″, tanδ, m, margin, Γ,
  T, α, β, z, CRLB σ, MC percentiles and all flags;
* `summary.json`: settings, versions, flag counts, physics checks with their status, KK
  statements with assumptions, and warnings;
* figures: parameters with MC bands, loss tangents, diagnostics (branch, NRW indicators,
  identifiability, passivity), and KK.

The CLI prints warnings and the KK statements.

---

## Reproducing this workflow without the package

Steps 5–8 need only the equations above. The published NIST equations are TN 1341 (2.24),
(2.25) and (2.38) [R4] (also TN 1536 eqs. 97–98 [R7]). Step 9 needs the 8×4 real Jacobian
of (Re, Im) S11, S21, S12, S22 with respect to (ε′, ε″, μ′, μ″), with Fisher information
F = 2JᵀJ/σ². Step 11 needs the kernel of `nrw.kk.kk_kernel_matrix` (closed-form
antiderivative in its docstring) and an NNLS solver.
