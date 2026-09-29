# Theory review: transmission/reflection retrieval of ε_r and μ_r

Scope: two-port transmission/reflection (T/R) measurements of a homogeneous, isotropic, linear
sample filling the cross-section of a coaxial airline (TEM) or rectangular waveguide (TE10),
roughly 2–18 GHz. Citations `[Rn]` refer to `docs/references.md`, which records how each source
was verified. Statements marked *(project analysis)* are our own derivations and are not
attributed to any paper.

---

## 1. Conventions, propagation, S-parameters, geometry

### 1.1 Time dependence and sign conventions

We use the conventions in `docs/conventions.md`: phasors ∝ exp(+jωt), ε_r = ε' − jε'',
μ_r = μ' − jμ''. This is the convention of the NIST reports [R4], [R6], [R7] (they write ε*_R,
μ*_R for the complex quantities). The physics literature (Landau & Lifshitz [R29], most
metamaterial papers such as [R12], [R18]–[R20]) uses exp(−iωt), where ε = ε' + iε''. For real ω the
two conventions are related by complex conjugation, so the **numbers ε', ε'' ≥ 0 are identical in
both**. Only the signs of imaginary parts in formulas change. Equations copied from exp(−iωt)
papers therefore need j → −i, i.e. complex conjugation.

### 1.2 Propagation

In a uniform guide filled with (ε_r, μ_r), fields vary as exp(−γz) with

$$\gamma=\sqrt{k_c^2-k_0^2\,\varepsilon_r\mu_r},\qquad k_0=\omega/c_0,\qquad
\gamma_0=\sqrt{k_c^2-k_0^2},$$

where k_c = 0 (coax TEM) or k_c = π/a (TE10). The branch of the square root is Re γ ≥ 0. With
λ_c = 2π/k_c and λ_0 = 2π/k_0, this is the NIST form
γ = j(2π/λ_0)·sqrt(ε_r μ_r − (λ_0/λ_c)²) ([R6] eq. (2.36)–(2.37) context; [R11] eq. (5)).

The sample transmission factor is T = exp(−γL). NIST writes it as z [R4] or Z [R6], [R7].
**Warning:** NIST's lower-case z is the transmission factor, not the impedance. In this project z
always means the normalized wave impedance.

### 1.3 Impedance and interface reflection

The normalized wave impedance for TE and TEM modes is z = μ_r γ_0/γ, and the air→sample
interface reflection is Γ = (z − 1)/(z + 1). This is identical to NIST's definition

$$\Gamma=\frac{\mu/\gamma-\mu_0/\gamma_0}{\mu/\gamma+\mu_0/\gamma_0}$$

([R6] eq. (2.37); [R7] eq. (113)). The NIST form holds for TE and TEM modes; for TM modes,
replace μ/γ by γ/ε ([R7] text after (113)). For TEM, z = sqrt(μ_r/ε_r) ([R6] eq. (2.38)).

### 1.4 S-parameters of the sample section

The sample of length L sits in an airline, with L_1 and L_2 the air sections between the
calibration planes and the sample faces (L_air = L_1 + L + L_2):

$$S_{11}=R_1^2\frac{\Gamma(1-T^2)}{1-\Gamma^2T^2},\quad
S_{22}=R_2^2\frac{\Gamma(1-T^2)}{1-\Gamma^2T^2},\quad
S_{21}=S_{12}=R_1R_2\frac{T(1-\Gamma^2)}{1-\Gamma^2T^2},\quad R_i=e^{-\gamma_0L_i}.$$

Source: [R6] eqs. (2.31)–(2.35); [R4] eqs. (2.12)–(2.16); [R7] §7.5. This matches
`conventions.md`. The empty holder gives S21⁰ = R1 R2 exp(−γ0 L) ([R4] eq. (2.20)). This
formulation assumes a single mode, a perfect sample (no gaps, faces normal to the axis, homogeneous
and isotropic) and a reciprocal sample (S12 = S21) ([R6] §2.1.3). NIST notes that 2π/λ_c must be
known accurately in waveguide ([R4] App. E), and recommends a mode-filter section of about two
free-space wavelengths between the coax-to-waveguide adapter and the holder ([R6] §2.3).

### 1.5 Reference-plane-invariant combinations

The following combinations are independent of L_1 and L_2 separately; they depend only on
L_air − L ([R4] eqs. (2.21)–(2.25); [R6] eqs. (2.51)–(2.52); [R7] eqs. (97)–(98)):

$$S_{11}S_{22}-S_{21}S_{12}=e^{-2\gamma_0(L_{air}-L)}\,\frac{\Gamma^2-T^2}{1-\Gamma^2T^2},\qquad
\frac{S_{21}+S_{12}}{2}=e^{-\gamma_0(L_{air}-L)}\,\frac{T(1-\Gamma^2)}{1-\Gamma^2T^2}.$$

Other invariants from [R4] eqs. (2.21)–(2.24) are |S11|, |S21|, S11S22/(S12S21) (coax only), and
S21/S21⁰ = e^{γ0 L} T(1−Γ²)/(1−T²Γ²), which needs an empty-holder measurement.

---

## 2. Classical NRW (Nicolson–Ross 1970; Weir 1974)

Nicolson and Ross [R1] worked in the time domain with a TEM fixture. Weir [R2] gave the
frequency-domain network-analyzer version for waveguide and TEM lines. The explicit equations
below are as reproduced by NIST ([R4] eqs. (2.26)–(2.35); [R6] eqs. (2.40)–(2.50); [R7]
eqs. (112)–(124)). The reference planes must first be rotated to the sample faces ([R6] §2.2.1).

1. $X=\dfrac{1-V_1V_2}{V_1-V_2}=\dfrac{S_{11}^2-S_{21}^2+1}{2S_{11}}$ with V1 = S21 + S11 and
   V2 = S21 − S11 ([R6] (2.41)–(2.43); the second form is [R11] eq. (7)).
2. $\Gamma=X\pm\sqrt{X^2-1}$. Choose the root with |Γ| ≤ 1 ([R6] (2.40)). Both roots satisfy
   Γ² − 2XΓ + 1 = 0, so their product is 1: exactly one has |Γ| < 1 unless |Γ| = 1
   *(project analysis)*. In our convention |Γ| ≤ 1 ⇔ Re z ≥ 0 (passive half-space).
3. $T=\dfrac{S_{11}+S_{21}-\Gamma}{1-(S_{11}+S_{21})\Gamma}$ ([R6] (2.44)).
4. Propagation constant, with branch index m (`conventions.md`):
   $$\gamma L=-\ln T=-\ln|T|+j(-\arg T+2\pi m).$$
   NIST writes $1/\Lambda^2=-[\ln(1/T)/(2\pi L)]^2$, i.e. γ = j2π/Λ ([R6] (2.45)).
5. Permeability and permittivity (general TE10/TEM form, [R6] (2.46)–(2.47)):
   $$\mu_r=\frac{1+\Gamma}{(1-\Gamma)\Lambda\sqrt{1/\lambda_0^2-1/\lambda_c^2}}
   =z\,\frac{\gamma}{\gamma_0},\qquad
   \varepsilon_r=\frac{\lambda_0^2}{\mu_r}\left[\frac{1}{\lambda_c^2}+\frac{1}{\Lambda^2}\right]
   =\frac{k_c^2-\gamma^2}{k_0^2\,\mu_r}.$$
   The equivalence of the two forms is *(project analysis)*: γ0 = j2π/λ_0g with
   1/λ_0g² = 1/λ0² − 1/λc². For **TEM** (k_c = 0) they reduce to
   $$\mu_r=z\,\frac{\gamma}{\gamma_0},\qquad \varepsilon_r=\frac{\gamma}{\gamma_0}\,\frac{1}{z},\qquad
   n=\sqrt{\varepsilon_r\mu_r}=\frac{\gamma}{\gamma_0},\quad z=\sqrt{\mu_r/\varepsilon_r}.$$
   An equivalent TEM-only closed form is
   z = ±sqrt{[(1+S11)² − S21²]/[(1−S11)² − S21²]} ([R11] eq. (16)). For **TE10**,
   μ_r = zγ/γ0 still holds, but ε_r ≠ γ/(γ0 z) because of the k_c² term.
6. Weir's group-delay branch selection ([R6] (2.48)–(2.50); [R7] (123)–(124)). Compute
   $$\tau_{calc}=L\frac{d}{df}\sqrt{\frac{\varepsilon_r\mu_r f^2}{c^2}-\frac{1}{\lambda_c^2}},$$
   for candidate n = 0, ±1, …, compare it with τ_meas = −(1/2π) dφ/df (φ = phase of T), and
   keep the n that agrees. Remarks:
   - [R6] eq. (2.48) prints a leading minus sign, while [R4] (2.34) and [R7] (123) do not. It is a
     sign-convention inconsistency between the NIST reports.
   - [R7] adds: "Usually we fit the measured group delay around a specific point by a
     least-squares polynomial" (§7.5.4).
   - NIST notes that because of this ambiguity the NRW solution is not explicit "in the strict
     sense" ([R6] §2.2.1; [R7] §7.5.1).

---

## 3. Singularity near L = nλ_g/2: an identifiability problem

### 3.1 What NIST says (verified wording)

* [R6] §2.2.1: "Where there is no loss in the sample under test, the NRW solution is divergent at
  integral multiples of one-half wavelength in the sample. This occurs because the phase of S11
  cannot be accurately measured for small |S11|. Also in this limit both of the scattering
  equations reduce to the relation Z² → 1, which is only a relation for the phase velocity and
  therefore solutions for ε*_R and μ*_R are not separable." [R7] §7.5.1 repeats this.
* [R6] Ch. 3 intro: "for low-loss materials both S21 and S11 become equations for the phase
  velocity, and the permittivity and permeability therefore enter as a product."
* [R11] §2.1 (NIST authors, 2014) shows on measured nylon data that the refractive index n stays
  smooth through the λ/2 resonance, while the wave impedance ς = z has a resonant artifact. Since
  ε_r = n/ς and μ_r = nς, ε and μ show mirror-image artifacts.

### 3.2 Mechanism *(project analysis)*

At a lossless half-wave condition, βL = mπ and α = 0, so T² = 1. Then for **any** Γ:

$$S_{11}=\frac{\Gamma(1-T^2)}{1-\Gamma^2T^2}=0,\qquad S_{21}=\frac{T(1-\Gamma^2)}{1-\Gamma^2}=T=\pm1 .$$

The S-matrix is therefore [[0, T], [T, 0]], independent of z. The consequences are:

* **Strict non-identifiability.** The data determine T, hence γ and the product ε_r μ_r
  (= (k_c² − γ²)/k_0²), but carry no information on z. The split ε/μ (TEM) or μγ0/γ (TE10) is
  unidentifiable at that frequency. No algorithm working on single-length, single-frequency
  two-port data can recover it, explicit or iterative. This is a property of the measurement, not
  of the NRW formula.
* **Near-resonance ill-conditioning.** With T fixed, ∂S11/∂Γ = (1−T²)(1+Γ²T²)/(1−Γ²T²)². Hence
  δΓ ≈ δS11/|1−T²| to leading order. Near βL = mπ we have |1 − T²| ≈ |1 − e^{−2αL}e^{−2jβL}|,
  whose minimum is 1 − e^{−2αL} ≈ 2αL. The noise amplification on Γ, and hence on z, ε and μ
  separately, is bounded by about 1/(2αL). It is finite for lossy samples and diverges as the loss
  goes to zero. This matches the NIST observation that the divergence is a low-loss phenomenon.
* **Where the NRW formula breaks.** Near the resonance X = (S11² − S21² + 1)/(2S11) →
  (1 − S21²)/(2S11) is a ratio of two small quantities, and Γ ≈ 1/(2X) ≈ S11/(1 − S21²). The phase
  of S11 sets the phase of Γ, and that phase is poorly measured when |S11| is small. [R11] saw the
  impedance start to oscillate already at |S11| ≈ 0.139 (|S21| ≈ 0.966) on nylon, and attributes it
  both to the small values and to the VNA noise floor.
* **What stays well determined.** T ≈ S21 (up to O(S11·Γ)) and hence n or εμ. This is why
  S21-only formulas with μ_r assumed (the NIST μ_r = 1 iterative method with β = 0, and the
  Boughriet formula) show no λ/2 divergence (§4, §5). They supply the missing information, μ_r,
  by assumption.

The half-wave condition refers to the guided wavelength in the sample,
λ_g = 2π/β = 2π/Im γ. For TE10 this is not λ_0/√(ε'μ').

---

## 4. Baker-Jarvis / NIST formulations

"Baker-Jarvis method" covers several distinct algorithms in the NIST reports. They differ in
unknowns, assumptions and λ/2 behaviour.

### 4.1 NIST permittivity-only iterative method (TN 1341 / Baker-Jarvis et al. 1990)

* **Source.** [R4] §2.2 and App. E (July 1990). The journal version [R3] has the same abstract; its
  full text was not accessible to us. NIST's 2005 restatement is [R7] §7.3.1.
* **Assumption.** μ_r = 1: "This procedure minimizes instabilities in determination of
  permittivity by setting the permeability to unity" ([R4] Ch. 1; also [R4] p. 10: "In this report
  μ*_R = 1").
* **Unknowns.** ε'_r and ε''_r (2 real). Optionally also L and L_air, when using the determinant
  equation ([R4] §2.2: "four real equations for the four unknowns: ε'_R, ε''_R, L, L_air").
* **Equation 1** — weighted combination; S-parameters must first be rotated to the sample faces
  using R1, R2 ([R4] eq. (2.38); [R7] eq. (97)):
  $$\tfrac12\{[S_{12}+S_{21}]+\beta[S_{11}+S_{22}]\}=\frac{T(1-\Gamma^2)+\beta\,\Gamma(1-T^2)}{1-T^2\Gamma^2}.$$
  β is a real constant weighting reflection against transmission. For low-loss materials "the
  S21 signal is strong and so we can set β = 0", while for high-loss materials "a large value of β
  is appropriate". "A general relation for β is given by the ratio of the uncertainty in S21
  divided by the uncertainty in S11" ([R4] §2.2; [R7] §7.3.1, where it is printed β1). [R4]
  Figs. 2.3–2.4 compare β = 0, 1, → ∞ for PTFE in X-band; "the best results are produced by using
  only S21 and S12 data" for low-loss PTFE.
* **Equation 2** — reference-plane invariant, used "when the reference plane positions are
  uncertain" ([R4] eq. (2.25); [R7] eq. (98)):
  $$S_{11}S_{22}-S_{21}S_{12}=e^{-2\gamma_0(L_{air}-L)}\frac{\Gamma^2-T^2}{1-T^2\Gamma^2}.$$
  "Equation (2.25) works well for both low-loss and high-loss materials" ([R4] §2.2). Also usable:
  (2.24), S21/S21⁰ = e^{γ0L}T(1−Γ²)/(1−T²Γ²), which "does not require the reference plane positions
  … but only the length of the waveguide and the sample length" ([R4] App. E).
* **Solver.** Newton–Raphson on the real and imaginary parts treated as separate real equations,
  x_k = x_{k−1} − J⁻¹(x_{k−1})F(x_{k−1}), with the Jacobian computed numerically
  ([R4] App. E, eqs. (6.177)–(6.180)).
* **Initial guess.** "we use the solutions to the Nicolson-Ross-Weir equations as a starting value
  and then use the previously obtained permittivity at one frequency as the initial guess for the
  next frequency" ([R4] §2.2).
* **λ/2 behaviour.** Stable for non-magnetic samples. [R4]'s abstract states that "the minimum
  uncertainty for transmission/reflection measurements of nonmagnetic materials occurs at integral
  multiples of one-half wavelength." It is **not** stable once μ is also unknown (§4.3).
* **Vendor naming.** Keysight's "NIST Precision" ("Reflection/Transmission Epsilon Precision")
  returns ε_r only, uses S11, S21, S12, S22, and cites [R3] as its reference ([R34], [R35]).
  Chalapat et al. describe it the same way: "The method of Baker-Jarvis is tailored for dielectric
  materials, thus one takes μr = 1" ([R12] footnote 1).

**Implementer's summary (TN 1341 / TN 1536 form).** At each frequency, solve F(ε') = 0 with
F = LHS(meas) − RHS(Γ(ε), T(ε)) for ε = ε' − jε''. Use μ_r = 1, Γ = (γ0 − γ)/(γ0 + γ) (from μ = 1),
T = exp(−γL) and γ = sqrt(k_c² − k0²ε). Use either Eq. 1 with a chosen β (β = 0 means
transmission-only, which has no λ/2 problem) or Eq. 2 (the determinant, with L_air, and optionally
L, as extra unknowns). Solve by Newton on [Re, Im]. Start from NRW at the first frequency, then
continue from the previous frequency. With β = 0 and a lossless sample, S21 alone still has
multiple solutions (γ is defined modulo 2πj/L). The continuation from the previous frequency is
what selects the branch.

### 4.2 Explicit μ_r = 1 formula attributed to Boughriet et al. (1997)

[R7] §7.3.2 (verified) states: "there are two independent transmission-line expressions for the
permittivity … Boughriet [91] has studied the second equation for permittivity. This equation is
relatively stable for low-loss materials":

$$\varepsilon_r=\frac{\lambda_0^2}{\mu_r}\left[\frac{1}{\lambda_c^2}-\left[\frac{1}{2\pi L}\ln\frac1T\right]^2\right]\quad(\mu_r=1),\qquad\text{[R7] eq. (99)}$$

where T comes from the NRW steps 1–3 and the branch of ln is chosen as in NRW. The alternative
expression, "relatively stable for high-loss materials", is [R7] eqs. (100)–(101):
ε_r = (c²/ω²)[((z+1)/(z−1))² + (2π/λ_c)²] with z = (Γ+1)/(Γ−1). The printed (100) appears to use
z as a propagation-type quantity; we have not validated it and do not recommend implementing it
from [R7] without an independent derivation.

The same content is often written in "effective parameter" notation (secondary source [R36]; the
original paper [R10] was not read):
μ_eff = 1, ε_eff = (λ_0g/Λ)², and
$$\varepsilon_r=\Bigl(1-\frac{\lambda_0^2}{\lambda_c^2}\Bigr)\varepsilon_{eff}+\frac{\lambda_0^2}{\lambda_c^2}\frac{1}{\mu_{eff}}.$$
With 1/λ_0g² = 1/λ0² − 1/λc², this is algebraically identical to [R7] eq. (99) when μ_r = 1
*(project analysis)*. The general identity with μ_eff = μ_r = zγ/γ0 and ε_eff = (λ_0g/Λ)²/μ_eff
reproduces NRW exactly *(project analysis)*.

Setting μ_eff = 1 makes ε depend only on T. T stays well conditioned at the resonance (§3.2), so
the formula is stable for low-loss **non-magnetic** samples. That is what [R10]'s abstract
claims ("no divergence … for low-loss materials"; "applicable to permittivity measurements"). It
is **closed-form (non-iterative)**, but it still needs the ln branch m.

### 4.3 NIST formulations for simultaneous ε and μ ([R6], [R7])

| Method (NIST name) | Unknowns / data | Solver | λ/2 behaviour (as stated by NIST) |
|---|---|---|---|
| NRW [R6 §2.2.1, R7 §7.5.1] | ε, μ; S11, S21 at sample faces | explicit + branch search | "divergent" for no loss at nλ/2 |
| "2-port solution where position is determined solely by L_air and L" [R6 §2.2.2] = "Modified Nicolson-Ross: Reference-Plane Invariant Algorithm" [R7 §7.5.2–7.5.4] | ε, μ; det S and (S21+S12)/2, eqs. (2.51)–(2.52) / (97)–(98) | Newton iteration ("works well if good initial guesses are available") **or** explicit solution (2.53)–(2.57) / (125)–(128) | unstable: "When both the permittivity and permeability are solved for in materials with both low-loss magnetic and dielectric properties, the solutions … will be unstable at integral multiples of one-half wavelength in the material" ([R7] §7.3.1; same in [R4] §2.2). [R6] Table 5.1: "Full 2-port, one sample, Magnetic … Unstable at nλm/2" |
| Two samples of different length [R6 §2.2.3] | ε, μ; S21 of lengths L and α1L, eqs. (2.58)–(2.59) | iterative | "unstable … if … both materials resonate at a certain frequency simultaneously" (§2.2.3); Table 5.1 calls it "Stable over all frequencies" (an internal inconsistency — the text is the more precise statement) |
| Short-circuit line, two positions / two samples [R6 Ch. 4] | ε, μ; one-port S11 | iterative | two-position: "Unstable at nλm/2"; two-sample: stable except simultaneous resonance (Table 5.2) |
| Optimized (multi-frequency) solution [R6 Ch. 3; R8; R9] | coefficients of a causal pole/Laurent model for ε(ω), μ(ω) + small corrections to L_1, L_2, λ_c; all frequencies at once | nonlinear least squares (ODRPACK, Levenberg–Marquardt type), initial guess from Stuchly or Weir | "stable over all frequencies" ([R6] abstract); "Very stable", weak point "non-global minima" (Table 5.1); "the existence of alternative local optima … required some testing" (§3.4.2) |

**Explicit reference-plane-invariant solution** ([R6] (2.53)–(2.57); [R7] (125)–(128)).
Let x = (S21S12 − S11S22)e^{2γ0(L_air−L)} and y = [(S21 + S12)/2]e^{γ0(L_air−L)}. Then:

$$T=\frac{x+1}{2y}\pm\sqrt{\Bigl(\frac{x+1}{2y}\Bigr)^2-1},\qquad
\Gamma_2=\pm\sqrt{\frac{x-T^2}{xT^2-1}},$$

$$\mu_r=-\frac{1+\Gamma_2}{1-\Gamma_2}\frac{1}{\gamma_0L}(\ln T+2\pi jn),\qquad
\varepsilon_r=\frac{c^2}{\omega^2}\Bigl[\Bigl(\frac{2\pi}{\lambda_c}\Bigr)^2-\frac{1}{L^2}(\ln T+2\pi jn)^2\Bigr]\Big/\mu_r .$$

Notes on this solution:

* We checked x = (T² − Γ²)/(1 − Γ²T²) and y = T(1 − Γ²)/(1 − Γ²T²). Hence (x+1)/y = T + 1/T, so
  the two T roots are reciprocal and |T| ≤ 1 selects one *(project analysis)*.
* The ± sign of Γ2 is not reference-plane invariant. NIST resolves it by comparison with
  Γ3 = [α(T² − 1) ± sqrt(α²T⁴ + 2T²(2S11 − α²) + α²)]/(2S11T²), where α = e^{−2γ0L1}, which "needs
  an estimate of L1" ([R6] (2.55)).
* In TEM, −Γ ↔ z⁻¹, which swaps ε and μ ([R12] §II makes the same observation).
* At T² = 1, x = 1 and Γ2² = (x − T²)/(xT² − 1) = 0/0. This is the identifiability singularity
  of §3.2 in explicit form *(project analysis)*.

[R11] (NIST, 2014) confirms for the iterative version: "the material parameters calculated from
these equations generate similar artifacts when simultaneously measuring both ε_r and μ_r of a
low-loss magnetic material. However, … the BJ method does not have the resonant issue if the
method uses only Eq. (13) [the (S21+S12)/2 equation] … assuming that μ_r = 1". [R11] also notes
that the earlier fixes ([R10], [R13], [R12]) "remain limited to fixing artifacts for only
dielectric material measurements". Its own proposal regresses the input impedance around the
resonance under a non-dispersive, locally flat characteristic-impedance assumption.

---

## 5. Answers to the terminology questions

**(1) Are "NIST-modified NRW method" and "Baker–Jarvis algorithm" equivalent terms?** No. Neither
is a uniquely defined algorithm.

* "Baker-Jarvis (BJ) method" in NIST's own later usage ([R11] §2.1) means the *iterative
  reference-plane-invariant* solution of det S and (S21 + S12)/2 ([R6] (2.51)–(2.52)). It can be
  run for (ε, μ) or, with μ_r = 1 and transmission only, for ε alone.
* The instrument name "NIST Precision" refers to the **permittivity-only (μ_r = 1)** algorithm of
  [R3]/[R4] ([R34], [R35]).
* The phrase closest to "NIST-modified NRW" in NIST documents is the section title "Modified
  Nicolson-Ross: Reference-Plane Invariant Algorithm" ([R7] §7.5.2). It denotes the (ε, μ)
  reference-plane-invariant equations (97)–(98), solved iteratively or explicitly, and NIST states
  they are unstable at nλ/2 for low-loss magnetic materials.

A text using either term must say which variant (ε-only vs ε, μ; weighted Eq. (2.38) vs
determinant; explicit vs Newton).

What "NIST Precision" / the 1990 method solves (verified via [R4]; journal text [R3] not read):

* Unknowns ε', ε'' (optionally also L, L_air), with **μ_r = 1**.
* Equations: ½{(S12 + S21) + β(S11 + S22)} = [T(1−Γ²) + βΓ(1−T²)]/(1 − T²Γ²) at sample-face
  planes, and/or the invariant S11S22 − S21S12 = e^{−2γ0(L_air−L)}(Γ² − T²)/(1 − T²Γ²), and/or
  S21/S21⁰.
* β is the ratio of S21 uncertainty to S11 uncertainty; β = 0 for low loss.
* Newton–Raphson on four or two real equations with a numerical Jacobian.
* Initial guess from NRW, then frequency continuation.

**(2) Which Baker-Jarvis/NIST method is appropriate for simultaneous ε and μ?** Single-sample
two-port data can be reduced by NRW or by the reference-plane-invariant (ε, μ) solution
([R6] §2.2.2; [R7] §7.5.2–7.5.4). Both are unstable at nλ_g/2 for low-loss samples.

NIST's remedies for magnetic materials are:

* two samples of different lengths ([R6] §2.2.3, eqs. (2.58)–(2.59)), chosen so that they do not
  resonate simultaneously;
* short-circuit-line two-sample measurements ([R6] Ch. 4);
* the multi-frequency causal-model regression ("optimized solution", [R6] Ch. 3, [R8]), which
  NIST's discussion ranks with the reference-plane-invariant solution as "the most accurate and
  stable methods. However, local minima have to be avoided" ([R6] Ch. 5).

For single-frequency work NIST notes that minimum-uncertainty magnetic measurements come from one
measurement at an electric-field maximum and one at a magnetic-field maximum ([R6] Ch. 5, SCL
context). The permittivity-only NIST iterative method is **not** appropriate for magnetic samples.

**(3) Iterative vs closed-form.**

* Closed-form apart from the ln branch: NRW [R1], [R2]; the explicit reference-plane-invariant
  (ε, μ) solution [R6] (2.53)–(2.57); Boughriet's μ = 1 formula [R10] / [R7] (99); Chalapat et al.
  [R12] (explicit, reference-plane invariant, group-delay based).
* Iterative: NIST ε-only Newton method [R3], [R4]; iterative (ε, μ) reference-plane-invariant
  solution [R6] §2.2.2; two-sample and SCL solutions [R6]; multi-frequency ODRPACK regression [R6]
  Ch. 3, [R8].

All closed-form methods still need a discrete branch decision ("Many researchers think of the NRW
solution as an explicit solution; however, due to the phase ambiguity, it is not in the strict
sense", [R6] §2.2.1).

**(4) Exact cause of the NRW singularity.**

* For a lossless sample at βL = mπ, T² = 1. S11 vanishes for every impedance and S21 = T, so the
  data fix εμ (via γ) but not z. That is loss of identifiability of the ε/μ split (§3.2).
* Near that point Γ is obtained as a ratio of two vanishing quantities. Noise in the phase of the
  small S11 is amplified by roughly 1/|1 − T²| ≥ 1/(2αL) *(project analysis)*.
* This is NIST's explanation: small |S11| gives an ill-measured phase, and both equations reduce to
  T² → 1, "a relation for the phase velocity", making ε and μ "not separable" ([R6] §2.2.1,
  Ch. 3; [R7] §7.5.1; [R11] §2.1).
* It is a property of single-length two-port data, not a defect of the particular NRW algebra.
  Reformulations using the same data inherit it, as shown for the iterative BJ method by [R11]
  and stated by [R7] §7.3.1.

**(5) How NIST methods mitigate it, and whether for simultaneous ε, μ.**

* Adding information (μ_r = 1) removes it, via NIST ε-only iterative [R3], [R4] or the Boughriet
  formula [R10], [R7] (99). This works for non-magnetic samples only.
* For (ε, μ), NIST mitigates it only by adding data or prior structure:
  - two lengths ([R6] §2.2.3), still unstable at simultaneous resonances;
  - two SCL samples ([R6] Ch. 4);
  - multi-frequency regression with a causal model that bridges the resonance ([R6] Ch. 3, [R8]);
  - local regression near λ/2 assuming a non-dispersive material ([R11]).
* The reference-plane-invariant formulation alone does **not** remove it ([R7] §7.3.1, Table 5.1
  of [R6], [R11]).

**(6) How to resolve the phase ambiguity.**

* **Group-delay comparison** (Weir [R2] as reproduced in [R6] (2.48)–(2.50), [R7] (123)–(124)).
  Strictly it compares the measured *group* delay with a delay computed from each candidate's ε_m
  μ_m. If d(ε_m μ_m)/df is computed consistently from the same data, the comparison becomes an
  identity for every m. The method discriminates only because the computed delay is effectively
  the *phase* delay (dispersion neglected), so it assumes weak dispersion *(project analysis)*.
  [R11]: "If the material sample possesses a dispersive characteristic, the correct branch cut
  will remain ambiguous". [R7]: "For dispersionless materials the group delay is constant with
  frequency."
* **Two sample lengths** ([R6] §2.3: "Another approach to root selection is the measurement of two
  samples of differing lengths"). [R17] concludes this is the only reviewed approach "not
  requiring an a priori knowledge about wavelength and dispersive properties of the material, nor
  … measurements over very wide frequency ranges".
* **Low-frequency start + continuity** ([R6] §2.3: at low frequencies "the roots are usually more
  widely spaced and therefore root selection is simplified"). [R17] §5 notes that one "can select
  the proper value of p by obtaining the low frequency values … and then change the value to
  ensure their continuity. However, it is doubtful whether this process can be automated for
  general structures … since branch crossings and other special features may obscure the choice
  of correct branches."
* **Group-velocity data in a reference-plane-invariant explicit method** ([R12]).
* **Kramers–Kronig reconstruction of β from α** ([R15] for dispersive materials; [R16] for
  metamaterials). [R17] cautions that truncating the KK integral has unquantified error, and that
  spatial-dispersion artifacts can invalidate it.
* **Amplitude-based / single-valued reformulations for low-loss dielectrics** ([R13]; [R14]
  exists, content unverified).
* **Stepwise NRW without branch seeking** for samples thicker than λ/2 ([R22], abstract only).
* **Metamaterial-retrieval branch/sign rules** ([R18]–[R20]). [R17] shows that the signs of z and
  n are dependent and "both sign combinations lead to the same permittivity and permeability". The
  only genuine ambiguity is the log branch.

**(7)–(8) Is a low-frequency-start continuous unwrap sufficient? Failure conditions.**
*(project analysis; consistent with [R6] §2.3, [R17] §5, [R11])*

Unwrapping arg T along frequency, starting from a correctly known branch m₀ at f_min, gives the
physical branch if and only if all of the following hold:

1. **Initial branch known.** Starting with m = 0 is correct only if βL ≤ π at f_min, i.e. the
   sample is shorter than λ_g/2 at the lowest frequency. For TEM this means
   f_min < c/(2L·Re n) (e.g. L = 10 mm, n = 5 gives f_min < 3 GHz). Otherwise m₀ must come from a
   prior (index estimate accurate to better than ±c/(2 f_min L)), from group delay (requires weak
   dispersion: roughly f·|τ_phase − τ_group| < ½), or from a second length.
   - In rectangular waveguide the lowest usable frequency is above the empty-guide cutoff. For a
     multi-band 2–18 GHz campaign each band (WR-284 … WR-62) starts at a high f_min. Each band
     must be anchored independently; continuity cannot carry the branch across bands.
2. **Sampling fine enough.** The phase change between consecutive points must be below π:
   Δf·τ_g < ½, where τ_g is the group delay of the quantity being unwrapped. Unwrapping raw S21
   before de-embedding means τ_g includes the full airline (L_air), not just the sample.
3. **Phase noise small.** σ_θ ≈ σ_S/|T| must be well below π/2 for every step. This fails for very
   lossy or long samples (|S21| near the VNA floor), near absorption/magnetic resonances (e.g.
   ferrite FMR in-band) and where higher-order-mode resonances notch S21 ([R6] Ch. 3 and 5
   discuss higher-mode resonances).
4. **No gaps.** Masked frequencies, band edges, or points discarded near λ/2 must not hide a 2π
   step. Continuity cannot bridge an unobserved interval longer than 1/(2τ_g).
5. **No genuine branch changes that continuity would suppress.** In the metamaterial example of [R17] §4, the retrieved
   parameters stay continuous only if the branch index is switched mid-band (p = 0 → 1 at about
   11.95 GHz), and [R17] warns that "branch crossings" can obscure the choice. A naive
   continuity rule can lock onto the wrong branch.

Continuous unwrapping from the low end is therefore a *heuristic that is correct under stated
conditions*, not a guarantee. It must be accompanied by a check such as a second length, group
delay vs phase delay, forward-model residuals, or passivity.

**(9) Legitimate checks.**

* *Passivity of the retrieved medium*: ε'' ≥ 0, μ'' ≥ 0 (exp(+jωt), derived in §6), Re γ ≥ 0 and
  |Γ| ≤ 1 (Re z ≥ 0), all within uncertainty. [R7] §7.5.5 states the root conditions
  Re(jωμ/γ) ≥ 0 and Re(γ) ≥ 0 for NIM extraction.
  - *Caveat*: for *effective* parameters of inhomogeneous metamaterials, one of ε'', μ'' can be
    negative although the structure is passive ([R21]). For homogeneous materials both must be
    ≥ 0.
* *Passivity of the measured network*: I − SᴴS positive semidefinite, i.e. the singular values of
  S are ≤ 1 (standard network theory; cf. [R28]).
* *Reciprocity*: S12 = S21 within VNA uncertainty (the model assumes it, [R6] §2.1.3).
* *Symmetry*: S11 R1⁻² = S22 R2⁻² for a homogeneous symmetric sample. The ratio
  S11/S22 = e^{−2γ0(L1−L2)} gives a position check.
* *Forward-model consistency*: re-simulate S from the retrieved (ε, μ) and compare with all four
  measured S-parameters.
* *Forward/reverse and two-length agreement*.
* *Group-delay vs phase-delay consistency* for the chosen branch.
* *Airline length* from the empty line ([R4] eq. (2.20); [R12] eq. (10)).
* *Causality*: consistency with KK in finite-band form (see (10)).

**(10) Can 2–18 GHz data "strictly satisfy the Kramers–Kronig relations"?** No, not as a verifiable
statement.

* The KK integrals run over 0 < ω < ∞ ([R4] App. F; [R6] App. D: "Direct solution requires
  complete data over full spectrum for one component").
* Band-limited data can at most be shown *consistent* with causality up to the unknown
  out-of-band contribution. That contribution appears in-band as a constant (ε_∞) plus a smooth,
  slowly varying function.
* Tools that make this quantitative:
  - singly or multiply subtractive KK with anchor points ([R4] eqs. (6.183)–(6.184) give the
    once-subtracted form; [R26]; commonly attributed to [R27]);
  - finite-range dispersion bounds ([R24]: "Unlike the Kramers-Kronig dispersion relations, they
    only require measurements over a finite frequency range");
  - fitting a model that is causal by construction (sum of Debye/Lorentz poles, NIST's causal
    Laurent/pole model [R6] Ch. 3, [R8]). The *model* then satisfies KK exactly, and the claim
    becomes "the data are reproduced by a causal model to within the stated residuals".
* KK-based branch selection ([R15], [R16]) inherits the truncation problem ([R17] §5.2.3).
* Defensible statements: "the retrieved spectra are consistent with causality within uncertainty
  over 2–18 GHz, as tested by …", or "the retrieved spectra are well represented by a causal
  (KK-consistent) dispersion model".

---

## 6. Passivity (derivation in the exp(+jωt) convention)

From ∇×E = −jωμ0μ_r H and ∇×H = jωε0ε_r E (no sources):

$$-\tfrac12\nabla\cdot(\mathbf E\times\mathbf H^*)=\tfrac12\bigl(j\omega\mu_0\mu_r|\mathbf H|^2-j\omega\varepsilon_0\varepsilon_r^*|\mathbf E|^2\bigr).$$

Taking the real part gives the time-averaged power absorbed per unit volume,

$$p=\tfrac{\omega}{2}\bigl(\varepsilon_0\varepsilon''|\mathbf E|^2+\mu_0\mu''|\mathbf H|^2\bigr).$$

A passive medium has p ≥ 0 for all fields, hence ε'' ≥ 0 and μ'' ≥ 0 (Im ε_r ≤ 0 and Im μ_r ≤ 0
in this convention) *(project analysis, standard)*. In the exp(−iωt) convention of Landau &
Lifshitz [R29] the same statement reads Im ε > 0. Passivity does **not** constrain the signs of ε'
or μ'. For a passive half-space, Re z ≥ 0 ⇔ |Γ| ≤ 1, and α = Re γ ≥ 0.

---

## 7. Causality and Kramers–Kronig in the exp(+jωt), ε = ε' − jε'' convention

Linear response is D(t) = ε0[ε_∞E(t) + ∫₀^∞ G(τ)E(t−τ)dτ] with real G, so
ε_r(ω) = ε_∞ + ∫₀^∞ G(τ)e^{−jωτ}dτ (compare [R6] eqs. (3.2)–(3.3)). This expression is analytic
in the **lower** half ω-plane (Im ω < 0) and satisfies ε_r(−ω) = ε_r(ω)*. Hence ε' is even and ε''
is odd. Integrating (ε_r − ε_∞)/(ω' − ω) along the real axis closed in the lower half-plane gives,
for a non-conducting dielectric *(project analysis)*:

$$\varepsilon'(\omega)-\varepsilon_\infty=\frac{2}{\pi}\,\mathrm{P}\!\!\int_0^\infty\frac{\omega'\,\varepsilon''(\omega')}{\omega'^2-\omega^2}\,d\omega',\qquad
\varepsilon''(\omega)=-\frac{2\omega}{\pi}\,\mathrm{P}\!\!\int_0^\infty\frac{\varepsilon'(\omega')-\varepsilon_\infty}{\omega'^2-\omega^2}\,d\omega'.$$

These are numerically **the same formulas as in the exp(−iωt) physics convention**, because ε' and
ε'' are the same numbers in both. The same relations hold for μ_r with μ_∞. For a conductor, add
σ_dc/(ε0ω) to ε'' and exclude it from the first integral.

**Debye check** (numerical, scratch script): take ε = ε_∞ + Δ/(1 + jωτ), so
ε' = ε_∞ + Δ/(1+ω²τ²) and ε'' = Δωτ/(1+ω²τ²) ≥ 0.

* At ω = 0 the first relation gives (2/π)∫₀^∞ Δτ/(1+ω'²τ²)dω' = Δ, which is correct.
* Numerical quadrature at ωτ = 0.2, 1.0, 2.4 (ε_∞ = 2, Δ = 3) reproduces ε' and ε'' to 8 digits.
  The + sign in front of the first integral is required.

**Discrepancy in the NIST reports.** [R4] eq. (6.181) and [R6] eq. (D.8) print the first relation
as ε'(ω) − ε_∞ = **−**(2/π)∫₀^∞[θε''(θ) − ωε''(ω)]/(θ² − ω²)dθ. The ωε''(ω) subtraction term
integrates to zero in principal value and is harmless. The leading minus sign, however, fails the
Debye test with ε'' ≥ 0. It is either a typographical sign error or an unstated opposite sign
convention for ε''. **Do not implement the printed sign.** The second relation ([R4] (6.182),
[R6] (D.9)) agrees with ours: P∫₀^∞dθ/(θ² − ω²) = 0 makes subtracting ε'(ω) equivalent to
subtracting ε_∞.

Once-subtracted forms (reference frequency ω₀, [R4] eqs. (6.183)–(6.184)) reduce sensitivity to
out-of-band data. The integrals over (−∞, ∞) use the even/odd symmetry. We did not re-derive their
printed signs; derive them from the relations above before use.

**Magnetic caveat.** [R29] (Ch. IX; § numbers not verified by us) discusses that μ(ω) loses
physical meaning at sufficiently high (optical) frequencies. At 2–18 GHz μ_r is meaningful, but
"μ_∞" must be read as the high-frequency limit of the *magnetic* dispersion in the model, not a
literal ω → ∞ value.

**Finite-band limitations.** Out-of-band poles contribute in-band terms that are indistinguishable
from a constant plus a slowly varying background. A finite-band test can reject gross
non-causality, such as a wrong branch giving a ±2π/L jump in β, but cannot certify causality
([R24]; [R17] §5.2.3; see §5 (10)).

---

## 8. Uncertainty

* **NIST analytic propagation.** [R6] §2.5, eqs. (2.67)–(2.68), for one sample at one position:
  $$\frac{\Delta\varepsilon'_R}{\varepsilon'_R}=\frac{1}{\varepsilon'_R}\sqrt{\sum_{\alpha}\Bigl[\Bigl(\frac{\partial\varepsilon'_R}{\partial|S_\alpha|}\Delta|S_\alpha|\Bigr)^2+\Bigl(\frac{\partial\varepsilon'_R}{\partial\theta_\alpha}\Delta\theta_\alpha\Bigr)^2\Bigr]+\Bigl(\frac{\partial\varepsilon'_R}{\partial L}\Delta L\Bigr)^2+\Bigl(\frac{\partial\varepsilon'_R}{\partial d}\Delta d\Bigr)^2}$$
  with α ∈ {11, 21}, d the air gap, and analogous forms for ε'' and μ. The partial derivatives come
  from implicit differentiation of the S-parameter equations ([R6] eqs. (2.69)–(2.81)). This is a
  first-order (Taylor) analysis that assumes small, independent input uncertainties ([R6] §2.5).
  [R6] lists as error sources the S-parameter magnitude and phase, gaps, holder dimensions, sample
  length, line losses and connector mismatch, and reference-plane positions.
* **Qualitative NIST conclusions.**
  - For low-loss materials, longer samples reduce uncertainty but risk over-moding.
  - For lossy materials, the optimum is about one attenuation length.
  - SCL: keep |S11| > −20 dB.
  - The ε-only method has minimum uncertainty at nλ/2 for non-magnetic samples.
  ([R6] Ch. 5; [R4] abstract; [R7] Fig. 25.)
* **Near λ/2 in (ε, μ) retrieval** the linearized uncertainty of the ε/μ split grows like
  1/|1 − T²| (§3.2). Linear propagation becomes unreliable where the Jacobian is nearly singular.
  A Monte Carlo propagation per JCGM 101 [R33] (GUM framework [R32]) is more defensible there, and
  results inside a declared neighbourhood of nλ_g/2 should be flagged or masked, not reported as
  measured.
* **Boughriet** [R10] (abstract) reports ε accuracy "similar to that obtained with a more recently
  proposed iterative technique". Its uncertainty equations were not read.
* **Interlaboratory evidence**: [R31] compares T/R permittivity results from eleven organizations
  against cavity references (abstract).

---

## 9. Reference-plane and air-gap effects

### 9.1 Reference planes

* NRW needs the S-parameters rotated to the sample faces, so errors δL_i produce phase errors
  2β0 δL_i in S_ii that grow with frequency ([R12] eq. (11); [R6] §2.2.1).
* The remedies are:
  - reference-plane-invariant combinations (§1.5; [R4], [R6], [R7], [R12]);
  - estimating L_1, L_2 from the data at a frequency where ε is known ([R4] eqs. (2.36)–(2.37));
  - letting a regression fit small length corrections (β_L1, β_L2) within their physical
    tolerance ([R6] Ch. 3 eqs. (3.14)–(3.16)).
* The explicit invariant (ε, μ) solution still needs an L_1 estimate to fix the sign of Γ
  ([R6] (2.55)).

### 9.2 Air gaps

* Gaps lower the apparent ε', especially in coax near the centre conductor ([R4] §2.3; [R6] §2.3.1
  "The effect of the air gaps is to measure values of the material parameters that are lower than
  the actual values"; [R7] §7.4.1).
* NIST capacitor model for coax dielectrics ([R6] App. C.2.1, eqs. (C.16)–(C.22)). Radii: R1 inner
  conductor, R2 sample inner, R3 sample outer, R4 outer conductor. Define L1 = ln(R2/R1) + ln(R4/R3),
  L2 = ln(R3/R2), L3 = ln(R4/R1). The approximate correction is
  $$\varepsilon'_{cR}=\varepsilon'_{mR}\frac{L_2}{L_3-\varepsilon'_{mR}L_1},\qquad
  \tan\delta_c=\tan\delta_m\Bigl[1+\varepsilon'_{mR}\frac{L_1}{L_2}\Bigr],$$
  which "breaks down when ε'_mR ≥ L3/L1". The exact forms are (C.16)–(C.17).
* Rectangular-waveguide E-plane gap: (C.23)–(C.24).
* Magnetic coax (series-inductance model, (C.33)–(C.34)):
  $$\mu'_{cR}=\frac{\mu'_{mR}\ln(R_4/R_1)-[\ln(R_2/R_1)+\ln(R_4/R_3)]}{\ln(R_3/R_2)},\qquad
  \mu''_{cR}=\mu''_{mR}\frac{\ln(R_4/R_1)}{\ln(R_3/R_2)} .$$
  The permeability correction is much smaller than the permittivity correction ([R6] C.3).
* NIST warns that literature gap formulas "generally under-correct for the real part of the
  permittivity and over-correct for the imaginary part" ([R6] §2.5).
* Mitigation: conducting paint, indium–gallium alloy or conducting grease in the gap ([R6] §C.5).
  Combining T/R with resonant measurements is another option ([R39], metadata only).
* Gap-correction heritage is cited by NIST to Westphal [R38] (not independently verified).

---

## 10. Finite-band and other limitations

* **Single mode only.** High-permittivity or inhomogeneous samples excite higher-order modes, whose
  resonances contaminate the data ([R6] Ch. 3, 5; [R7] §7.4.3).
* **Low-frequency end.** Short samples have poor sensitivity at low frequency ([R12] abstract:
  "the uncertainties are larger at low frequencies and at the Fabry–Pérot resonances").
* **Waveguide bands.** 2–18 GHz in rectangular guide needs several bands. Branch and continuity
  arguments are per band (§5 (7)).
* **KK / causality over 2–18 GHz is not closed** (§7).
* **Homogeneity.** For structured or composite samples the retrieved parameters are *effective*
  and may violate homogeneous-medium sign rules ([R19], [R21]).

---

## 11. Summary for implementation

1. **NRW (TEM/TE10)**: §2 equations, with the branch m reported and chosen by the documented
   strategy of §5 (6)–(8).
2. **NIST ε-only iterative** (μ_r = 1): Newton on Eq. (2.38)/(97) with weight β (default 0), or on
   the determinant (2.25)/(98). Initial guess from NRW, then frequency continuation.
3. **Boughriet NNI** (μ_r = 1): ε_r = λ0²[1/λc² − (ln(1/T)/(2πL))²], i.e. ε_r = −(λ0/2πL)²(ln T)² +
   (λ0/λc)² with the branch m. This is closed form.
4. **(ε, μ) reference-plane invariant** ([R6] (2.51)–(2.57)): expect and flag λ/2 ill-conditioning.
   Report the condition indicator 1/|1 − T²| or the Jacobian condition number.
5. **Masking and uncertainty**: mask or flag nλ_g/2 neighbourhoods for (ε, μ). Use Monte Carlo
   uncertainty near them.
6. **Checks**: passivity, reciprocity, and finite-band causality consistency, not "strict KK".
