# Paper-ready methodology

Text for the Methods section of a materials/electromagnetics paper using this software.
Bracketed items [ … ] must be filled with the actual measurement details. Reference tags
[Rn] map to `docs/references.md`. Every quantitative statement here comes from
`docs/benchmark_report.md` (BR), which is generated from `results/`.

---

## Methods — Retrieval of complex permittivity and permeability

**Measurement.** S-parameters (S11, S21, S12, S22) of [material] samples of length
[L ± u(L)] were measured from 2 to 18 GHz in a [7 mm coaxial air line of length L_air] with
a [VNA model] after [TRL/SOLT] calibration. The reference planes were [at the sample faces /
at distances d1, d2 from the faces, de-embedded analytically].

**Constitutive-parameter retrieval.** We use the exp(+jωt) convention, with
ε_r = ε′ − jε″ and μ_r = μ′ − jμ″; passive media then have ε″, μ″ ≥ 0.

1. *Initial solution.* The Nicolson–Ross–Weir (NRW) equations [R1, R2] give the interface
   reflection coefficient Γ (root |Γ| ≤ 1, equivalent to Re z ≥ 0) and the sample transmission
   factor T = exp(−γL).
2. *Phase branch.* The branch index m in γL = −ln T + j2πm was determined as follows.
   * The phase of T was unwrapped within segments free of unreliable steps.
   * Each segment's integer offset was fixed by requiring the non-dispersive delay computed
     from the candidate branch to match the measured group delay in the band median [R2].
   * The result was cross-checked by an independent dynamic-programming tracker based on the
     continuity of εμ.
   * The branch index m(f) is reported. [No / N] frequencies had disagreeing trackers.
3. *Material solution.*
   * *(Non-magnetic reference samples)* ε_r was obtained, with μ_r = 1 assumed, by solving the
     reference-plane-invariant NIST transmission/reflection equations [R4, eqs. (2.24)–(2.25);
     R3] with damped Newton iteration, starting from the NRW value.
   * *(Magnetic samples)* ε_r and μ_r were obtained simultaneously by least-squares fitting of
     the single-slab model to all four S-parameters at each frequency, starting from the NRW
     solution on the selected branch.
   * [If two lengths were measured: ε_r and μ_r were fitted jointly to both samples
     (L1 = …, L2 = …), which restores the ε/μ separation near the half-wavelength frequencies
     of either sample [R6].]
4. *Conditioning.* For a single low-loss sample, the ε/μ separation is ill-conditioned where
   the sample length approaches an integer multiple of half the guided wavelength
   (βL → kπ) [R6, R7].
   * We therefore computed the Cramér–Rao bound [R40] of (ε′, ε″, μ′, μ″) at the measured
     noise level ([σ]).
   * Frequencies where the bound on ε′ exceeded 1 % are marked in the figures [and excluded
     from the averages in Table X].
5. *Physical checks.*
   * Passivity: ε″ ≥ 0 and μ″ ≥ 0 within the uncertainty; |Γ| ≤ 1; |T| ≤ 1; the largest
     singular value of the measured S-matrix ≤ 1.
   * Reciprocity and port symmetry were checked as model assumptions.
6. *Finite-band causal consistency.* The real part was predicted from the measured loss by a
   singly-subtractive Kramers–Kronig relation [R26, R27], evaluated over 2–18 GHz and
   complemented by the best non-negative (passive) out-of-band loss continuation.
   * Residuals were compared with the combined uncertainty using a Bonferroni-corrected
     threshold (α = 0.01) [R43].
   * [ε_r and μ_r passed / did not pass] this diagnostic.
   * The test establishes consistency with causality and passivity *within the measured band
     and under the stated assumptions*. It is not a proof that the global Kramers–Kronig
     relations hold [R24].
7. *Uncertainty.* Standard uncertainties were propagated by Monte Carlo (JCGM 101 [R33];
   [N] draws) through the complete retrieval, including branch selection. The inputs
   perturbed were:
   * |S11|, |S21| [± …%];
   * the phases of S11 and S21 [± …°];
   * the sample length [± … µm];
   * the reference-plane positions [± … µm].

   We report medians with 95 % intervals.

**Software.** All steps are implemented in the open-source package `nrw` [repository URL,
commit hash]. Its test suite includes forward/inverse recovery with an independent
transfer-matrix model, branch-crossing, noise and singular-region tests.

---

## Claims that should NOT be made

The ratings use the scale of `docs/terminology_verdicts.md`. The primary-source reasoning is
in `docs/theory_review.md` §5. The experimental evidence below is from this project (BR).

### 1. "The NIST-modified NRW method eliminates half-wave divergence."

**Verdict: incorrect for simultaneous ε and μ. Conditionally correct only for the NIST
permittivity-only method with μ_r = 1.**

* "NIST-modified NRW" is not a defined NIST term. The NIST iterative method of Baker-Jarvis
  et al. [R3, R4] avoids the divergence **because it sets μ_r = 1** ("minimizes instabilities
  … by setting the permeability to unity", [R4]).
* NIST states that simultaneous (ε, μ) solutions are unstable at nλ/2 [R6, R7].
* *Our experiment* (BR §3; low-loss sample, L = 20 mm, σ = 10⁻³, 20 noise draws):
  * Near βL = kπ, the RMSE of ε′ is 0.131 for NRW and 0.131 for iterative (ε, μ) least
    squares. The Cramér–Rao bound for *any* unbiased single-sample (ε, μ) estimator is 0.136.
    Away from these points the RMSE is 0.006. The divergence is therefore a loss of
    information in the data, not a defect of the NRW algebra, and no reformulation of the
    same data removes it.
  * The μ = 1 NIST method gives an RMSE of 9e-4 at its own bound.
  * On a magnetic sample (μ = 1.5) the μ = 1 method is biased by ≈ 2.0 in ε′.
  * A second sample length reduces the near-half-wave RMSE of ε′ to 0.006.

**Replacement:** "For non-magnetic samples, ε_r was retrieved with the NIST
transmission/reflection permittivity algorithm (μ_r = 1 assumed) [Baker-Jarvis et al. 1990],
which is not affected by the half-wavelength instability of the NRW equations. When ε_r and
μ_r are retrieved simultaneously from a single sample, the separation of ε and μ is
intrinsically ill-conditioned near L = nλ_g/2. These frequencies are identified from the
Cramér–Rao bound and flagged [or: resolved with a second sample length]."

### 2. "Continuous phase unwrapping guarantees the correct physical phase branch."

**Verdict: incorrect as a guarantee. The procedure is conditionally correct.**

* Unwrapping only propagates the branch assumed at the first frequency.
* *Our experiment* (BR §4; ε = 10 − 0.1j unless stated, 10 noise draws per condition).
  Plain unwrap from the lowest frequency with m = 0 assigned a wrong branch **in every
  draw, over essentially the whole band** when:
  * the sample was thick (30 or 60 mm), so βL > π at 2 GHz;
  * the band started at 8 GHz;
  * a 4 GHz gap in the data was combined with a thick sample;
  * the grid was coarse enough to alias (17 points, 60 mm).
* Group-delay-anchored segment unwrapping fixed all of these except the aliased grid, which
  only the dynamic-programming tracker handled.
* With a strongly dispersive resonance that drives |S21| to the noise floor, **all**
  phase-based strategies failed in part of the band.

**Replacement:** "The branch index m(f) of the transmission phase was determined by
segment-wise unwrapping anchored to the measured group delay and cross-checked by an
independent continuity-based tracker. The lowest frequency satisfies βL [< π / = …], the
frequency step satisfies Δf·τ_g < 0.5, and the selected branch reproduces the measured group
delay. The branch index is reported."

### 3. "The retrieved parameters strictly satisfy the Kramers–Kronig relations."

**Verdict: incorrect / unverifiable from 2–18 GHz data.**

* The KK integrals run over (0, ∞). Out-of-band absorption adds unknown in-band
  contributions, so finite-band data can only be *consistent* with causality [R24].
* *Our experiment* (BR §6). The finite-band diagnostic has real power: a phase-branch error
  was rejected in 100 % of runs, and an uncorrected 50 µm reference-plane error was rejected
  in μ for two of three magnetic cases and in ε for the dielectric cases. But the following **passed**:
  * a 1 % thickness error (magnetic cases);
  * 10 µm air gaps;
  * a 51 Ω holder mismatch (magnetic cases);
  * a constant offset of +0.5 in ε′;
  * an added rising trend +0.3(f/18 GHz)² in ε′.

  These are indistinguishable from a passive continuation outside the band.
* The software therefore never prints "satisfies the Kramers–Kronig relations".

**Replacement:** "The retrieved constitutive parameters were checked for passivity and for
finite-band causal consistency using a singly-subtractive Kramers–Kronig diagnostic with a
passive out-of-band continuation. Within the 2–18 GHz band and the stated uncertainty, no
inconsistency was detected." (Use only if the diagnostic passed. Otherwise report the
failure and the frequencies involved.)

### 4. "…using the Nicolson–Ross–Weir method … a NIST/Baker–Jarvis-type robust formulation."

**Verdict: ambiguous. "Robust" is too strong for single-sample (ε, μ).** Name the equations
actually solved (NRW closed form; NIST μ = 1 iterative; (ε, μ) least squares; two-length),
how the branch was chosen, and how the half-wave regions were treated. See
`docs/terminology_verdicts.md` row 4 for replacement wording.

### 5. Further claims to avoid

* "The method is exact / has no uncertainty at resonance." Near βL = kπ the uncertainty of
  single-sample (ε, μ) is large *by the Cramér–Rao bound*.
* "Passivity (ε″ ≥ 0) proves the result is correct." It is necessary, not sufficient: wrong
  branches can remain passive. In our benchmark the textbook m = 0 NRW produced large errors
  whose losses stayed mostly positive.
* "A 1 % thickness error causes a 1 % error in ε." In a transmission-dominated retrieval
  with μ = 1, the measured phase fixes βL, so ε ∝ (βL)²/(k0 L)² ∝ L⁻². A 1 % length error then
  gives 1.8–2.0 % in ε′ for the dielectric cases, and 1.2–1.4 % for the thin and near-unity cases
  (BR §2, `thickness_error` scenario, NIST μ = 1 method).
* "The Kramers–Kronig test validates the thickness/air gap/calibration." It does not
  (item 3).
