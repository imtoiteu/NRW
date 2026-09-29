# Conventions used throughout this repository

These conventions are fixed. Every module, document and figure uses them.

## Time dependence and constitutive parameters

* Phasor time dependence: **exp(+jωt)** (engineering / NIST convention, as in
  Baker-Jarvis et al., NIST TN 1355-R and TN 1536).
* ε_r = ε' − jε'' and μ_r = μ' − jμ''.
* With exp(+jωt), a **passive** isotropic medium has ε'' ≥ 0 and μ'' ≥ 0
  (i.e. Im ε_r ≤ 0, Im μ_r ≤ 0). Derivation: docs/theory_review.md §Passivity.
* Loss tangents: tanδ_ε = ε''/ε', tanδ_μ = μ''/μ'.

## Propagation

* Fields in a uniform guide vary as exp(−γz), γ = α + jβ, with α ≥ 0 for a
  passive medium (principal square root with Re γ ≥ 0 in the forward model).
* γ = sqrt(k_c² − k0² ε_r μ_r), k0 = ω/c0, with k_c = 0 for coaxial TEM and
  free space at normal incidence, k_c = π/a for rectangular TE10.
* Empty-guide propagation constant: γ0 = sqrt(k_c² − k0²) (= jβ0 above cutoff).
* Sample transmission factor: T = exp(−γL), L = sample length (m).

## Impedance and reflection

* Normalized sample wave impedance (relative to the empty guide):
  z = μ_r γ0 / γ (TE modes and TEM; for TEM this reduces to sqrt(μ_r/ε_r)).
* Interface reflection coefficient (air→sample): Γ = (z − 1)/(z + 1).
* Passivity of the half-space ⇔ Re z ≥ 0 ⇔ |Γ| ≤ 1.

## S-parameters

* Reference planes at the sample faces unless offsets are declared:
  S11 = Γ(1 − T²)/(1 − Γ²T²), S21 = T(1 − Γ²)/(1 − Γ²T²).
* Offsets: a port-i reference plane located a distance L_i (m) in the empty
  guide *away* from the sample face multiplies S_ii by R_i² and S_21 by R_1R_2,
  with R_i = exp(−γ0 L_i).
* Units: frequency in Hz internally (CLI accepts GHz), lengths in metres.

## Branch index

* With T = |T| exp(jθ), θ = arg T ∈ (−π, π]:
  γL = −ln T = −ln|T| + j(−θ + 2πm), i.e. αL = −ln|T| and βL = −θ + 2πm.
* The integer **m** is the phase-branch index reported by the software. For an
  ordinary (positive-index) passive sample βL > −π, so m ≥ 0; m = 0 means the
  electrical length βL lies in (−π, π].
