# nrw — robust extraction of complex ε<sub>r</sub> and μ<sub>r</sub> from two-port T/R S-parameters

Research software and technical report on retrieving complex permittivity
ε<sub>r</sub> = ε′ − jε″ and permeability μ<sub>r</sub> = μ′ − jμ″ from measured
S11, S21, S12, S22 of a sample in a coaxial air line (TEM) or rectangular waveguide
(TE10), aimed at the 2–18 GHz range.

The package implements the classical Nicolson–Ross–Weir (NRW) equations, the NIST
permittivity-only iterative method (Baker-Jarvis et al. 1990 / NIST TN 1341), the
Boughriet non-iterative (μ = 1) formula, pointwise and regularised iterative
simultaneous (ε, μ) least squares, and two-thickness extraction. It adds explicit
phase-branch tracking, local identifiability (Cramér–Rao) analysis, passivity checks
and a finite-band Kramers–Kronig *consistency* diagnostic. It returns diagnostics
instead of silently returning questionable numbers.

**All results in this repository come from synthetic ground-truth data.** No measured
data are included (see `data/raw/README.md`).

## Main findings (details: `docs/benchmark_report.md`)

* The NRW divergence near βL = kπ is an **identifiability limit of single-sample
  simultaneous (ε, μ) data**, not only a formula artefact. Classical NRW and iterative
  least squares both follow the Cramér–Rao bound (CRLB) up to the same spike (RMSE of ε′
  near half-wave 0.131 vs a CRLB of 0.136; `results/studies/singularity_summary.csv`).
  No re-formulation removes it.
* The NIST iterative method avoids the spike **because it assumes μ<sub>r</sub> = 1**.
  On a non-magnetic sample it reaches its own CRLB (RMSE ≈ 9e-4). On a magnetic sample it
  returns an effective ε with a bias of ≈ 2 in ε′ (μ = 1.5 case).
* A second sample thickness restores identifiability of (ε, μ): RMSE of ε′ near
  half-wave falls from 0.13 to 0.006.
* "Continuous unwrapping from the lowest frequency" fails over the **whole band** whenever
  βL > π at the lowest frequency (thick samples, high start frequency). Plain unwrapping
  also fails across gaps and on aliased grids. Group-delay-anchored segment unwrapping
  and dynamic-programming branch tracking fix these cases
  (`results/studies/branch_summary_*.csv`). No phase-only method survives a band where
  |S21| reaches the noise floor.
* A finite-band KK diagnostic can reveal gross errors (a wrong branch, uncorrected
  reference planes). It cannot detect scale errors such as sample thickness, air gaps or
  holder mismatch, and it can never prove global KK compliance.

## Install

```bash
uv venv .venv -p 3.12            # or: python3.12 -m venv .venv
uv pip install -p .venv/bin/python -e ".[dev]"
# exact versions used for the shipped results:
uv pip install -p .venv/bin/python -r requirements-lock.txt --extra-index-url https://download.pytorch.org/whl/cpu
```

Optional, only to re-run the third-party comparison: `scripts/fetch_external.sh`
(clones the pinned repositories into the gitignored `research/external/`), GNU Octave ≥ 8,
`g++`, and `uv pip install -p .venv/bin/python -e ".[external]"`.

## Run the tests

```bash
.venv/bin/python -m pytest            # add --cov=nrw for coverage
```

## Python API

```python
from nrw import extract_material_parameters

res = extract_material_parameters(
    frequency=f_hz, s11=s11, s21=s21, s12=s12, s22=s22,   # complex arrays
    thickness=2.0e-3,                                      # m
    geometry="coax",                                       # or {"kind": "waveguide", "a_m": 22.86e-3}
    method="iterative_eps_mu",   # "auto" | "nrw" | "nrw_eps" | "baker_jarvis" | "iterative_eps_mu" | "regularized_eps_mu"
    branch_method="anchored",    # "principal" | "unwrap" | "group_delay" | "anchored" | "dp" | "prior"
    noise_sigma=1e-3,            # rms complex noise per S-parameter (for identifiability flags)
    uncertainty=True,            # Monte Carlo through the whole pipeline
)
res.epsilon_complex, res.mu_complex, res.selected_branch, res.flags["singular_region"]
res.causality["eps"].statement   # the only KK wording the software endorses
res.save("results/my_sample")    # results.csv, summary.json, figures
```

`method="baker_jarvis"` is the NIST ε-only method (μ<sub>r</sub> = 1 **assumed**). Use it
only for non-magnetic samples; the software warns when NRW indicates μ ≠ 1.

## Processing a measurement (.s2p)

```bash
# magnetic or unknown sample, planes at the sample faces
.venv/bin/python -m nrw.cli --input data/raw/sample01.s2p --thickness-mm 2.0 \
    --method iterative-eps-mu --uncertainty --output results/sample01

# non-magnetic sample, reference-plane-invariant NIST formulation (holder length known)
.venv/bin/python -m nrw.cli --input data/raw/sample02.s2p --thickness-mm 10.0 \
    --method baker-jarvis --holder-length-mm 50.0 --uncertainty --output results/sample02

# waveguide (WR-90), planes offset from the sample faces
.venv/bin/python -m nrw.cli --input data/raw/x_band.s2p --thickness-mm 5.0 --geometry waveguide \
    --waveguide-a-mm 22.86 --offset1-mm 3.0 --offset2-mm 3.0 --output results/x_band
```

`python -m nrw.cli --help` lists all options: frequency window, noise level,
Monte Carlo spec, seed, branch method, and Touchstone backend.

Example on shipped **synthetic** data (the output is committed in `results/examples/`):

```bash
.venv/bin/python examples/make_synthetic_s2p.py
.venv/bin/python -m nrw.cli --input data/synthetic/synthetic_magnetic_composite_2mm.s2p \
    --thickness-mm 2.0 --method iterative-eps-mu --uncertainty --mc-draws 100 \
    --output results/examples/magnetic_composite_2mm
.venv/bin/python -m nrw.cli --input data/synthetic/synthetic_lowloss_dielectric_10mm.s2p \
    --thickness-mm 10.0 --method baker-jarvis --uncertainty --mc-draws 100 \
    --output results/examples/lowloss_dielectric_10mm
```

## Reproduce benchmarks, tables and figures

```bash
.venv/bin/python scripts/run_benchmarks.py --seeds 5 --jobs 4   # results/benchmarks/*.csv (add --no-external to skip third-party codes)
.venv/bin/python scripts/study_singularity.py                   # half-wavelength / CRLB study
.venv/bin/python scripts/study_branch.py                        # phase-branch failure map
.venv/bin/python scripts/study_bj_formulation.py                # NIST eps-only formulation / beta choice
.venv/bin/python scripts/study_kk.py                            # KK diagnostic power
.venv/bin/python scripts/study_uncertainty.py                   # Monte Carlo coverage
.venv/bin/python scripts/make_figures.py                        # comparison figures
.venv/bin/python scripts/make_tables.py                         # markdown tables used in docs/benchmark_report.md
```

All random draws are seeded. `results/benchmarks/environment.json` records the versions
and settings used. Runtime on a shared 4-core machine: the internal benchmark takes
minutes; the full external comparison takes about an hour because of the Octave-based codes.

## Repository layout

| path | content |
|---|---|
| `src/nrw/` | package: `geometry`, `materials`, `forward` (independent ABCD model), `model`, `classical` (NRW), `branch`, `iterative` (NIST ε-only, (ε, μ) LSQ, regularised, multi-thickness), `conditioning` (CRLB), `physics`, `kk`, `uncertainty`, `io`, `api`, `cli`, `plotting` |
| `tests/` | pytest suite |
| `scripts/` | benchmark harness and studies |
| `examples/` | synthetic data generator |
| `docs/` | theory review, references, survey, benchmark report, final method, paper methodology, conventions, terminology verdicts |
| `research/` | open-source survey (CSV), pinned external repositories, adapters |
| `results/` | benchmark tables, study outputs, figures, example outputs |
| `data/raw/` | location for future measured data (format spec inside) |

## Conventions

exp(+jωt); ε<sub>r</sub> = ε′ − jε″; passive ⇒ ε″, μ″ ≥ 0; fields ∝ exp(−γz);
T = exp(−γL); z = μγ0/γ; Γ = (z − 1)/(z + 1); branch index m with βL = −arg T + 2πm.
See `docs/conventions.md`.

## Licence

MIT for this repository's own code. Third-party code is **not** redistributed; see
`research/external/README.md` and `docs/github_survey.md` for the licences of the surveyed
projects.
