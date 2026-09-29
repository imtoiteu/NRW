# Survey of open-source NRW / transmission–reflection extraction codes

*Survey date: 2026-09-29. All statements about a repository come from reading
its code at the pinned commit (`research/external_pins.csv`); the smoke-test
numbers come from running it through the adapters in `research/adapters/`.
Machine-readable version: `research/github_survey.csv`.*

Conventions used for all comparisons: `docs/conventions.md` (exp(+jωt),
ε = ε′ − jε″, passive ⇒ Im ε ≤ 0, planes at the sample faces).

---

## 1. Search method

### 1.1 GitHub repository search (REST, unauthenticated)

Endpoint `https://api.github.com/search/repositories?q=<query>&per_page=30`
(the code-search endpoint needs authentication, so it was not used), with
≥7 s between requests (10 req/min limit). One request hit the limit and was
repeated. Exact queries, in the order issued on 2026-09-29, with hit counts:

| # | query | hits | | # | query | hits |
|---|---|---|---|---|---|---|
| 1 | nicolson ross weir python | 0 | | 27 | material characterization waveguide | 2 |
| 2 | Nicolson Ross Weir | 2 | | 28 | transmission reflection permittivity | rate-limited → #41 |
| 3 | nicolson-ross-weir | 2 | | 29 | dielectric permittivity waveguide | 3 |
| 4 | NRW permittivity permeability | 2 | | 30 | free space permittivity measurement | 0 |
| 5 | NRW permittivity | 4 | | 31 | effective parameter retrieval | 3 |
| 6 | S parameter material parameter extraction | 4 | | 32 | S-parameter retrieval metamaterial | 0 |
| 7 | Baker Jarvis permittivity | 0 | | 33 | complex permittivity measurement vna | 0 |
| 8 | NIST transmission reflection dielectric | 0 | | 34 | material characterization vna | 0 |
| 9 | NIST NRW | 0 | | 35 | permittivity coaxial | 7 |
| 10 | phase unwrap NRW | 0 | | 36 | effective medium retrieval S parameters | 0 |
| 11 | complex permittivity S11 S21 | 0 | | 37 | permittivity from s-parameters | 5 |
| 12 | material parameter retrieval touchstone | 0 | | 41 | transmission reflection permittivity | 2 |
| 13 | waveguide permittivity permeability extraction | 0 | | 42 | permittivitycalc | 1 |
| 14 | coaxial transmission reflection method | 0 | | 43 | airline permittivity | 0 |
| 15 | metamaterial parameter retrieval | 1 | | 44 | dielectric measurement s-parameters python | 0 |
| 16 | permittivity extraction matlab | 4 | | 45 | material parameters s-parameters | 12 |
| 17 | epsilon mu retrieval S-parameters julia | 0 | | 46 | NRW algorithm | 6 |
| 18 | permittivity extraction | 14 | | 47 | Ross Weir | 4 |
| 19 | permittivity measurement | 16 | | 48 | retrieval permittivity permeability | 0 |
| 20 | dielectric characterization S-parameters | 0 | | 49 | eps mu extraction | 0 |
| 21 | NRW method | 8 | | 50 | dielectric properties coaxial airline | 0 |
| 22 | NRW extraction | 4 | | 51 | ASTM D5568 | 1 |
| 23 | nicholson ross weir | 2 | | 61 | S-parameter retrieval | 5 |
| 24 | Nicholson-Ross-Weir | 2 | | 62 | effective parameters metamaterial | 1 |
| 25 | permittivity permeability extraction | 1 | | 63 | metamaterial retrieval | 1 |
| 26 | permittivity permeability S-parameters | 1 | | 64 | electromagnetic material characterization | 2 |
| | | | | 65 | permeability extraction | 5 |
| | | | | 66–68 | transmission line permittivity / dielectric constant waveguide / free-space measurement permittivity | 0 |
| | | | | 69 | complex permittivity | 16 |
| | | | | 70 | language:julia permittivity | 2 |

(Query numbering has gaps where batches were restarted; no other queries were
issued.) Note that GitHub's repository search matches name/description/topics
only, so repositories that implement NRW without saying so in their metadata
are invisible to it — hence the web search below.

### 1.2 Web search

WebSearch queries (2026-09-29): "Nicolson Ross Weir python github permittivity
permeability S-parameters"; "github NRW material parameter extraction waveguide
S-parameters matlab"; "github Baker-Jarvis iterative permittivity transmission
reflection code"; "MATLAB File Exchange Nicolson Ross Weir permittivity
permeability extraction"; "github metamaterial effective parameter retrieval
S-parameters python Smith Chen retrieval"; "julia OR python package
"permittivity" "permeability" "transmission/reflection" open source code github
coaxial airline"; "github S-parameter retrieval epsilon mu phase unwrapping
branch python notebook".

New candidates from the web: `jdhuang-csm/scatter-opt` (not found by the
repository search), `marcdmc/meep-metamaterials`, and confirmation of
`permittivitycalc` (also on PyPI). The MATLAB File Exchange / ResearchGate /
figshare hits ("Matlab code of NRW algorithm…", PLOS figshare item) are not on
GitHub, have no version control, and were not pursued. MathWorks Answers
threads contain snippets only.

### 1.3 Known-package checks

* **scikit-rf 2.1.0** (installed in the venv): no NRW / material-extraction
  routine (it has Touchstone I/O, calibration and media models only).
* **"PyMaterials", "mwpermittivity"**: no such GitHub repositories or PyPI
  packages were found by the searches above; they do not appear to exist.
* **Julia**: only `Permittivity.jl` (dispersion models, no S-parameter
  inversion) and `LorentzDrudeMetals.jl`; nothing relevant.
* Open-ended coaxial-probe codes (`tyoon124/PyOECP`, `vna-dielectric-permittivity`,
  `DielectricPermittivityMeasurementTool`, the Stuchly–Stuchly MATLAB repo) are
  reflection-only methods and therefore out of scope; two were cloned for the
  record, the others were only screened.

### 1.4 Inclusion criteria

A repository was cloned if its metadata or the web suggested two-port
transmission/reflection extraction of ε (and possibly μ) from S-parameters, or
closely related retrieval (metamaterial retrieval, waveguide T/R). 22 were
cloned; 12 contain a relevant, runnable extraction routine and received an
adapter; 10 are recorded as screened-out (no code, different method, or
out of scope).

---

## 2. Candidate table

"Max rel. err." is max over the 401 frequencies of |x − x_true|/|x_true|,
taking the larger of ε and μ (ε only if the code returns no μ), on the two
required smoke-test cases (coax, 2–18 GHz, L = 5 mm): **A** ε = 4 − 0.1j, μ = 1;
**B** ε = 10 − 0.5j, μ = 2 − 0.3j. "by design" = the code assumes μ = 1.

| repository (adapter) | lang. | licence | method | ε/μ | geometry | branch handling | A | B |
|---|---|---|---|---|---|---|---|---|
| boivinalex/permittivitycalc (`permittivitycalc_nrw`) | Py | MIT | NRW | ε, μ | coax | np.unwrap | 0.025 | 0.099 |
| ″ (`permittivitycalc_nni`) | Py | MIT | NNI (Boughriet) | ε | coax | np.unwrap | 4.6e-14 | 1.0 (by design) |
| jdhuang-csm/scatter-opt (`scatter_opt_nrw`) | MATLAB | BSD-3 | NRW + branch integer | ε, μ | coax/WG | unwrap + base branch | 5.9e-13 | 5.9e-13 |
| ″ (`scatter_opt_gd`) | MATLAB | BSD-3 | NRW + group-delay base branch | ε, μ | coax/WG | GD + unwrap | 5.9e-13 | 5.9e-13 |
| kasbrot123/nrw_method_eps_mu (`kasbrot_nrw`) | Py | none | NRW + per-point GD search | ε, μ | coax/WG | GD argmin | 4.4e-15 | 1.4e-15 |
| shikouzhong/Nicholson-Ross-Weir-Algorithm (`shikouzhong_nrw`) | MATLAB | MIT | NRW fwd/rev avg | ε, μ | coax/WG | first jump only | 4.7e-15 | 1.6e-15 |
| tomaasaybr02/5G-New-Material… (`bj_tem_5g`) | MATLAB | MIT | Baker-Jarvis ref.-plane-invariant closed form | ε, μ | TEM | unwrap + NN | 8.2e-15 | 3.5e-15 |
| njchorda/Nicolson-Ross-Weir-… (`njchorda_nrw`) | MATLAB | MIT | NRW (TN 1536) | ε, μ | WG (coax via a=∞) | none (m = 0) | 0.34 | 1.0 |
| sergeyermakov239/WaveGuide_Measurements (`ermakov_wg`) | MATLAB | none | T-only, μ = 1 | ε | WG (coax via a=∞) | cumulative phase | 5.0e-5 | 1.0 (by design) |
| atom711/nrw-material-extraction (`atom711_nrw`) | Py | none | T-only, μ = 1 | ε | TEM | np.unwrap | 2.1e-15 | 1.0 (by design) |
| dmtang-jpg/onlyEE-GUI (`onlyee_method2`) | Py | none | NRW (z², n² with continuous √) | ε, μ | TEM | np.unwrap | 6.9e-4 | 6.9e-4 |
| ″ (`onlyee_method1`) | Py | none | T-quadratic, μ = 1 | ε | TEM | np.unwrap | 1.4e-3 | 1.0 (by design) |
| bumbeneciconstantinbogdan/S-Param-Inversion (`sparam_inversion_nrw`) | Py/torch | MIT | NRW (baseline for an MLP) | ε, μ | WG (coax via a=∞) | none | 2.0 | 2.0 |
| marcdmc/meep-metamaterials (`meep_retrieval`) | Py | MIT | Smith/Chen retrieval | ε, μ | TEM | fixed m = 0 | 2.0 | 2.0 |
| Dr-Musho-Research-Group/NRW_Permittivity_Extraction (`musho_cpp`) | C++ | none | NRW, μ given | ε | coax (with TE11 λc!) | none | 90 | 37 |

Informative extra cases (all in `research/adapters/smoke_results.json`):

| adapter | wg_A (WR-90, 8.2–12.4 GHz, ε=4−0.1j) | coax_B_L12 (several wraps) | coax_B_L30 (m = 1 at f_min) |
|---|---|---|---|
| permittivitycalc_nrw | 0.030 | 0.099 | 1.0 |
| permittivitycalc_nni | 1.3 | 1.0 | 1.0 |
| scatter_opt_nrw | 1.6e-12 | 5.9e-13 | 1.0 |
| scatter_opt_gd | **8.0** | 5.9e-13 | 5.9e-13 |
| kasbrot_nrw | 1.8e-15 | 1.1e-15 | 1.3e-15 |
| shikouzhong_nrw | 1.6e-15 | **0.80** | 1.0 |
| bj_tem_5g | n/a (TEM) | 3.9e-15 | 1.1 |
| njchorda_nrw | 1.8e-15 | 1.0 | 1.0 |
| ermakov_wg | 5.0e-5 | 1.0 (μ≠1) | 1.0 |
| onlyee_method2 | n/a | 6.9e-4 | 1.9 |
| sparam_inversion_nrw | 3.7e-10 | 2.0 | 1.3 |
| meep_retrieval | n/a | 2.0 | 1.3 |
| musho_cpp | n/a | 37 | 35 |

Screened out (cloned, no adapter): `senglongyu/S2RFP` (CPW multiline
propagation constant, different method), `Viceduroz/waveguide-dielectric-optimization`
(|S21|-only fit of one constant ε for an 8-slab stack), `mikkokataja/2D_metamaterial_homogenization`
(field-averaging homogenization, no S-parameters), `lsito/ParameterRetrieval`,
`anahithambardzumyan/YPI-…`, `ZiadHatab/dielectric-permittivity-extraction`,
`DrTissot/Exact_HF` (README/docs only, no code at the pinned commit),
`Vladmg1101/Solid-dielectric-materials` (PDF only), and the two open-ended
probe codes.

---

## 3. Per-repository review notes

### 3.1 boivinalex/permittivitycalc — Python, MIT, last commit 2021-01-26
The only packaged, documented, tested (CI config, ~100 unittest cases) Python
code in the survey; on PyPI. `sparam_data.AirlineData._permittivity_calc`
implements NRW (`nrw=True`) and Boughriet's NNI (`nrw=False`) for coax
airlines (VAL/PAL/GAL), forward (S11,S21) and reverse (S22,S12) separately,
then averages. `iter_data.AirlineIter` fits multi-pole Cole-Cole ε (and
optionally μ) directly to the S-parameters with lmfit/emcee (a
"modified Baker-Jarvis" approach after Houtz et al. 2016); its forward model
uses T = exp(−γL), ε = k′ − jk″ (our convention). Input is the METAS VNA
Tools II text export with uncertainties; bundled data are real measurements
(rexolite, serpentine). Extras: washer de-embedding with a synthetic PTFE
model, boundary (air-gap) correction, bulk-density normalisation,
Boughriet/Baker-Jarvis uncertainty propagation.

Defects found (code-level, confirmed numerically):
1. **Conjugated propagation constant.** `ln_1_T = 1j*t_phase_unwrap + np.log(1/abs(t))`
   uses +arg T; for exp(+jωt) data ln(1/T) = ln(1/|T|) − j arg T. In NNI mode
   this just conjugates ε (so the reported `lossfac` is ε″ > 0 and the result
   is exact after `eps = dielec − j·lossfac`). In NRW mode the conjugated γ is
   combined with the *unconjugated* Γ, giving μ = z·n* and ε = n*/z: for a
   μ = 1 sample it returns ε = |ε| (real) and |μ| = 1 with a spurious phase —
   wrong for every lossy material (case A error 0.025, i.e. ε″ lost; case B
   0.099). The package's own `_de_embed` and `iter_data` models use the correct
   sign, so the package is internally inconsistent.
2. **Γ root test** `sign[np.absolute(gam1<=1)]` compares complex numbers with
   `<=` (NumPy lexicographic order: real part first) instead of |Γ| ≤ 1. For
   ordinary dielectrics (Re Γ < 0) it happens to select correctly; over
   uniformly drawn |Γ| < 1 it selects 1/Γ in 25 % of cases (Re Γ > 0 region,
   e.g. μ > ε samples).
3. **NNI cutoff typo** `(lam_0**2/LAM_C*2)` instead of `/LAM_C**2`: harmless
   for coax (LAM_C = ∞), wrong for waveguide (wg_A error 1.3).
4. Uncertainty code uses `new_t = 1j*phase + |t|`, not a representation of T.
5. Tests only assert `IsNotNone`; there is no numerical check of any result.
Branch: `np.unwrap` from the principal value at f_min (fails when m ≥ 1 at f_min, L30).

### 3.2 jdhuang-csm/scatter-opt — MATLAB, BSD-3, last commit 2024-02-15
The most complete algorithmic toolbox: NRW and NNI with an explicit branch
integer (`nrw_lambda`: ln|1/T| + j(unwrap(angle(1/T)) + 2πp)), base-branch
selection by GroupDelay, Kramers–Kronig (Szabó 2010) or rational fitting of γ
(RF Toolbox `rationalfit`, the default), sample-position optimisation
(`OptimizeL` → R1, R2 de-embedding), geometric-mean S11/S22 option, and
Domich/Baker-Jarvis-style nonlinear least-squares fitting of pole models
(`lm_*`, `de_*`, `lsqcurvefit`). Author disclaimer: "hastily developed",
untouched since 2020. With the correct base branch it is exact (6e-13, from
`c_vac = 299792458.000176`).

Defects: (1) **the GroupDelay test is degenerate**: `tcalc = -L*gradient(sqrt(eps.*mu./lambda_0.^2 - 1/lambda_c.^2))`
equals −L·d(1/Λ_p)/df and the branch offset p/L is frequency-independent, so
every branch without a sign flip gets the same error. Confirmed: for the
WR-90 case branches 0, 1, 2 all have mean error 0.010644 and branch 2 was
selected (error 8.0); it worked on coax only by tie-breaking luck.
(2) The KK method uses `k0 = 2*Sdata.freq./c_vac` (missing π).
(3) `nrw_ref_coef` treats Γ = 0 exactly as "no root". (4) Default branch
method and `table` usage require MATLAB RF Toolbox / MATLAB-only classes.
`pre_gamma_0`/`rev_gamma_0` use c_lab = 299704644.54 (speed in air) for the
empty line, while NRW uses c_vac — a deliberate but undocumented mix.

### 3.3 kasbrot123/nrw_method_eps_mu — Python, no licence, 2025-08-17
Compact `nrw()` with a per-frequency search over n ∈ [−10, 10] minimising
|Re τ_calc + Re τ_meas|; a gmpy2 high-precision variant for stability studies.
**Written for exp(−iωt)**: the criterion needs d(arg T)/df > 0 (its demo uses
ε = 8.4(1 + 0.05j)). Fed VNA-convention data directly it chooses wrong
branches everywhere (case A error 13); the adapter conjugates in and out,
after which all cases are exact, including L30. Caveat: as in scatter-opt,
τ_calc = L·d(1/Λ_n)/df is analytically independent of n (we measured
τ_calc = 0.448 ns for n = −1…3 at mid-band); the choice is made by the sign
flips of the principal square root and by round-off, and the branch is
chosen independently at every frequency. It was right in all our cases but
this is not a sound discriminator. `nrw_reverse` picks the |Γ| > 1 root and
returns NaN for `lambda_cut = inf`.

### 3.4 shikouzhong/Nicholson-Ross-Weir-Algorithm — MATLAB, MIT, 2017-04-24
`NRW.m` is a clean NRW function (Re(1/Λ) ≥ 0 root, returns ε″, μ″ positive);
the driver averages forward and reverse and finds the sample position with a
genetic algorithm (`ga`) that equalises the S11/S22 phases. Real measured coax
data (holder 6.95 mm, 4.74 mm and 2.00 mm samples) are included.
Defect: `for m = index + 1 : length(tempLog)` with a vector `index` uses only
the first element, so exactly 2π is added after the first wrap and later
wraps are not handled (coax_B_L12: 40 % of points correct). Cannot start on
m ≥ 1. Driver uses c = 3.0e8.

### 3.5 tomaasaybr02/5G-New-Material-EM-Characterization-Techniques — MATLAB, MIT, 2025-03-19
Master's-project code; `NRW_BJ_TEM_{2,3,medio}.m` implement the
Baker-Jarvis (TN 1355-R) **reference-plane-invariant** closed form: A =
S11S22/(S21S12), B = e^{2γ0(L_air−L)}(S21S12 − S11S22), Γ² from a quadratic,
T from S21/S21_air, n from unwrap + user integer NN. The only surveyed code
that uses all four S-parameters plus an empty-holder measurement. Exact on
TEM cases with m = 0 at f_min. Defects: Γ's sign is fixed by the principal
square root and Z = (1 − Γ)/(1 + Γ) — implicitly assumes Re Γ_true < 0 (true
for ε > μ; would swap ε and μ for μ-dominated samples); TEM only; a hard-coded
length uncertainty (0.1433 mm) in the sensitivity formula; many near-duplicate
scripts.

### 3.6 njchorda/Nicolson-Ross-Weir-Materal-Parameter-Extraction — MATLAB script, MIT, 2026-03-02
Readable WR-90 script following TN 1536 p. 112, with an HFSS-simulated PTFE
example. **No branch handling** (`angle(Z1) + -2*pi*0`); τ_calc/τ_meas are
computed but never used. Exact while βL < π (wg_A), wrong and non-passive
beyond (coax A: 76 non-passive points). `magTanD = imag/real` vs
`tanD = −imag/real` — inconsistent signs.

### 3.7 sergeyermakov239/WaveGuide_Measurements — MATLAB Live Script, no licence, 2025-04-13
ε-only (μ = 1) WR-guide extraction with fixture de-embedding and cumulative
phase tracking; real Agilent E8362C measurements (empty guide + 4 samples,
10.6–70.2 mm). README claims Luukkonen's stepwise NRW, but the code is plain
phase accumulation plus a manual `num_2pi`. c = 2.998e8 (5e-5 bias).
All geometry hard-coded in the function.

### 3.8 atom711/nrw-material-extraction — Python, no licence, 2026-09-22
Educational ε-only extraction ε = −(ln P/(k0 d))² (μ = 1, TEM) with unwrap;
exact for case A. Synthetic data only.

### 3.9 dmtang-jpg/onlyEE-GUI — Python/PyQt, no licence, 2026-05-19
GUI with VNA SCPI control. `core.method_II` computes μ = √(z²n²), ε = √(n²/z²)
with square-root branches from an unwrapped phase; `method_I` (μ = 1) uses
the T-quadratic. Automatic sample-position estimation from the S11/S22 phase
difference. c = 3e8 hard-coded (6.9e-4 bias). The continuous-√ approach
inherits a wrong initial branch at every point (L30: all 401 points
non-passive). `load_s2p` ignores the Touchstone option line (MA assumed).

### 3.10 bumbeneciconstantinbogdan/S-Param-Inversion — Python/PyTorch, MIT, 2026-07-19
NRW baseline `nrw_inverse` (principal log only) used to motivate an MLP
inverse model; the baseline fails as soon as βL > π, which overstates NRW's
weakness relative to the MLP. 1e-10 residual from μ0 = 4π·10⁻⁷ vs c.

### 3.11 marcdmc/meep-metamaterials — Python, MIT, 2023-01-10
Smith/Chen-type retrieval for MEEP simulations, exp(−iωt) convention (adapter
conjugates). Default fixed branch m = 0; optional "continuity" mode checks
the derivative of the base branch only with a unit-dependent threshold.
Importing it requires the MEEP stack (stubbed in the adapter during import).

### 3.12 Dr-Musho-Research-Group/NRW_Permittivity_Extraction — C++17, no licence, 2026-06-03
Well-commented CLI, but it deliberately inserts λc = π(a + b)/2 (TE11
cut-off) into the TEM NRW formula, giving huge low-frequency bias (case A
max error 90); no unwrapping; c = 3e8; the "regression test" baseline
ε′ ≈ 12097 is accepted as passing.

---

## 4. State of practice (synthesis)

* **Algorithms.** 11 of 12 runnable codes are the classic NRW Γ/T inversion
  (or its μ = 1 T-only reduction); one implements the Baker-Jarvis
  reference-plane-invariant closed form; one (scatter-opt) and one
  (permittivitycalc) add model-based least-squares fitting. **No surveyed
  code implements the NIST iterative per-frequency Baker-Jarvis solver**
  (Newton on the TN 1355-R S-parameter equations); "Baker-Jarvis" in these
  repositories means either the invariant closed form or a global pole-model fit.
* **Branch (phase) ambiguity** is the dominant failure mode. Practice ranges
  from none (m = 0: njchorda, S-Param-Inversion, Musho, meep default) through
  `np.unwrap` from the principal value at f_min (permittivitycalc, atom711,
  onlyEE, scatter-opt manual, 5G) to group-delay selection (kasbrot,
  scatter-opt). Both group-delay implementations differentiate 1/Λ_m, which
  makes the criterion branch-independent; only scatter-opt's KK and
  rational-fit options are genuinely branch-discriminating, and neither runs
  outside MATLAB (RF Toolbox) — and KK has a missing-π bug. Nobody reports the
  branch index m or a confidence in it.
* **Γ root selection** is |Γ| ≤ 1 everywhere except permittivitycalc
  (lexicographic complex comparison) and the 5G code (principal √ of Γ²).
* **Half-wavelength (|S11| → 0) instability** is acknowledged only by
  permittivitycalc (NNI alternative, resonance-aware averaging) and
  scatter-opt (NNI, S smoothing); no code flags or down-weights points near
  n·λg/2 in NRW mode.
* **Reference planes.** Handled by: position optimisation (scatter-opt
  OptimizeL, shikouzhong `ga`, onlyEE phase-difference), fixed fixture lengths
  (ermakov), washer de-embedding (permittivitycalc), or invariance (5G).
* **S-parameter use.** Most use S11/S21 only; forward/reverse averaging in
  permittivitycalc and shikouzhong; geometric-mean S11/S22 option in
  scatter-opt; all four in 5G.
* **Sign conventions** are a real hazard: two codes are written for exp(−iωt)
  (kasbrot, meep), one mixes conventions internally (permittivitycalc NRW),
  and several report ε″ with an ad-hoc sign flip or `abs()`.
* **Constants.** c = 3e8 or 2.998e8 in four codes (biases of 1e-4–1e-3 in ε).
* **Validation.** Only permittivitycalc (unittest, CI configs) and the Musho
  C++ code (a Windows .bat regression script against its own previous output)
  have automated tests; neither checks against a known truth. No code ships a synthetic round-trip test;
  measured example data exist in permittivitycalc, shikouzhong, ermakov and
  (provenance unclear) scatter-opt and Musho.
* **Touchstone.** Native readers only in onlyEE (MA only); others rely on
  MATLAB RF Toolbox, a third-party reader, or bespoke text formats.

---

## 5. Licensing constraints

| licence | repositories | what we may do |
|---|---|---|
| MIT | permittivitycalc, njchorda, shikouzhong, 5G-New-Material, S-Param-Inversion, meep-metamaterials, vna-dielectric-permittivity | reuse/adapt with copyright + licence notice |
| BSD-3-Clause | scatter-opt, ZiadHatab (no code) | reuse with notice; no endorsement use of names |
| GPL-3.0 | Exact_HF (no code) | copyleft — do not mix into MIT `src/` |
| custom permissive text | DielectricPermittivityMeasurementTool | read the text before any reuse |
| **none** | nrw_method_eps_mu, nrw-material-extraction, onlyEE-GUI, WaveGuide_Measurements, NRW_Permittivity_Extraction, S2RFP, waveguide-dielectric-optimization, 2D_metamaterial_homogenization, ParameterRetrieval, YPI, Solid-dielectric-materials | **all rights reserved**: may be read and run locally for evaluation; must not be copied, modified-and-distributed, or vendored |

Consequences for this project: **no external code is copied into `src/`**.
The external repositories are not committed (`research/external/` is
gitignored) and are re-created at the pinned commits by
`scripts/fetch_external.sh`. Adapters call the external code in place; the
two script-type MATLAB sources (njchorda, ermakov) are wrapped *at run time*
from the pinned checkout into a temporary directory, so no third-party text is
stored in this repository. Our own implementation in `src/nrw` should be
written from the primary literature (NRW 1970/1974, Baker-Jarvis TN 1355-R /
TN 1536, Boughriet 1997), not from any of these codes.

---

## 6. Adapters and smoke test (reproduction)

```bash
scripts/fetch_external.sh                       # clone all pins
uv pip install -p .venv/bin/python uncertainties lmfit emcee corner seaborn tqdm gmpy2
uv pip install -p .venv/bin/python torch --index-url https://download.pytorch.org/whl/cpu
apt-get install -y octave                       # 8.4.0 installed in ~8.5 min here
.venv/bin/python -m research.adapters.smoke_test   # writes research/adapters/smoke_results.json
```

`research/adapters/__init__.py` exposes `ADAPTERS` (15 entries, 12
repositories) and `available()`; every adapter has the signature
`extract(freq_hz, s11, s21, s12, s22, thickness_m, geometry) -> {"eps", "mu", "notes"}`
and returns our convention. Adapters raise `NotImplementedError` for
geometries the external code cannot represent (`GEOMETRIES` in the registry).
Interface-only changes made to run the external code (all documented in each
adapter's docstring):

| adapter | what was needed besides calling the function |
|---|---|
| permittivitycalc_* | METAS-style mag/deg array, L in cm; tkinter stubbed during import; LAM_C module global set to 2a (cm) for waveguide; ε = dielec − j·lossfac |
| kasbrot_nrw, meep_retrieval | conjugate inputs and outputs (exp(−iωt) codes); meep: MEEP/IPython stubbed during import, freqs = f/c0 |
| atom711_nrw | overwrite `freqs/omega/k_0` of the engine with the caller's grid |
| scatter_opt_* | Octave; struct instead of `table`; shims `array2table.m`, `height.m`, no-op `import.m` + `+scatter_opt` on the path; branch method 'manual' (m = 0) or 'GroupDelay' |
| shikouzhong_nrw | Octave; fwd/rev average as in the driver; no `ga` plane search; ε = er − j·ei |
| njchorda_nrw | Octave; calc block extracted verbatim at run time and wrapped in a function (inputs replace file/dimension sections; physconst → 299792458) |
| bj_tem_5g | Octave; L_air = L, ideal empty-line S21_air = exp(−jk0L), NN = 0 |
| ermakov_wg | Octave; function extracted from the .mlx at run time; 4 exact text edits (signature, RF-Toolbox file read removed, a from geometry, fixture lengths → 0) |
| onlyee_* | load `core.py` by path |
| sparam_inversion_nrw | torch tensors; coax via a = ∞ |
| musho_cpp | compile with g++, CSV in/out; μ = 1 passed; APC-7 radii default |

The smoke results above show that, once the interface is right, **four codes
are exact to round-off on both required cases** (scatter-opt, kasbrot after
conjugation, shikouzhong, 5G Baker-Jarvis) and that every failure traces to a
specific branch, sign-convention, constant or model-assumption issue listed in
Section 3.
