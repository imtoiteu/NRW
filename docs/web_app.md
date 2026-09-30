# NRW Material Parameter Extraction Laboratory — web application

A browser front end for the tested `nrw` package. Every number shown comes from the core
(`nrw.extract_material_parameters`, `nrw.assessment`, `nrw.kk`, `nrw.uncertainty`). The web
layer (`src/nrw/web/`) contains no extraction algorithm of its own.

```
app.py                     Streamlit entry point (3 lines)
src/nrw/web/ui.py          page layout and widgets
src/nrw/web/analysis.py    configuration → core call → result/diagnostic tables, method comparison
src/nrw/web/plots.py       interactive (Plotly) and publication (Matplotlib) figures
src/nrw/web/export.py      CSV / JSON / PNG / HTML report / ZIP (built in memory)
src/nrw/web/security.py    password gate, rate limiting
src/nrw/web/state.py       session-scoped storage and "delete session data"
src/nrw/web/version.py     package version + Git commit shown in the UI and reports
src/nrw/web/rest.py        optional FastAPI endpoint
src/nrw/assessment.py      GREEN/YELLOW/RED rules (core, also usable from Python)
.streamlit/config.toml     server settings (no secrets)
scripts/run_web.sh         launcher (0.0.0.0:8501 by default)
scripts/run_api.sh         optional REST API launcher (127.0.0.1:8502 by default)
```

---

## 1. Install

```bash
cd /root/imtoiteu/NRW
uv venv .venv -p 3.12                                  # once (or: python3.12 -m venv .venv)
uv pip install -p .venv/bin/python -e ".[web]"         # core + streamlit, plotly, pandas, scikit-rf, fastapi
```

## 2. Set the password (never stored in Git)

The application reads the password **only** from the environment variable `NRW_APP_PASSWORD`
(at least 8 characters). If it is unset, the app shows a "locked" page and accepts no uploads.

A convenient way to keep it out of shell history and out of the repository:

```bash
mkdir -p ~/.config/nrw && chmod 700 ~/.config/nrw
python3 -c "import secrets; print(secrets.token_urlsafe(18))" > ~/.config/nrw/app_password
chmod 600 ~/.config/nrw/app_password
cat ~/.config/nrw/app_password          # the password to type into the browser
export NRW_APP_PASSWORD="$(cat ~/.config/nrw/app_password)"
```

To change it, write a new value to that file and restart the server.

## 3. Start / stop

```bash
export NRW_APP_PASSWORD="$(cat ~/.config/nrw/app_password)"
export NRW_APP_PORT=8501            # optional, default 8501
./scripts/run_web.sh                # binds 0.0.0.0:$NRW_APP_PORT; Ctrl+C stops it
```

Equivalent direct command:

```bash
.venv/bin/python -m streamlit run app.py --server.address 0.0.0.0 --server.port 8501
```

### Persistent run inside tmux

```bash
tmux new-session -d -s nrw-web \
  'export NRW_APP_PASSWORD="$(cat ~/.config/nrw/app_password)"; cd /root/imtoiteu/NRW && ./scripts/run_web.sh'
tmux ls                              # shows the session
tmux attach -t nrw-web               # view the log; detach again with Ctrl+b then d
tmux kill-session -t nrw-web         # stop the server
```

After changing code, restart the server: the configuration disables automatic reloading
(`fileWatcherType = "none"`).

### Find the VPS IP and open the browser

```bash
curl -s https://api.ipify.org; echo       # public IPv4 of the VPS
hostname -I                               # all local addresses
curl -s -o /dev/null -w '%{http_code}\n' http://127.0.0.1:8501/_stcore/health   # 200 = server up
```

Then open `http://<VPS_IP>:8501` from any browser.

### Firewall

Check before changing anything:

```bash
ufw status              # "inactive" = UFW does not filter
iptables -S INPUT | head
```

If UFW is **active**, allow the port with `ufw allow 8501/tcp`. A hosting-provider firewall
(cloud security group) may also need a rule for TCP 8501.

## 4. Security model and limitations

* Password gate before any upload or processing. Constant-time comparison. Per-IP rate
  limiting: 5 failures within 5 min lock that IP for 5 min, and every failure costs 1 s.
* The password is never written to disk by the application, never logged, and never kept in
  session state (only an "authenticated" flag).
* Uploaded files and results stay in the user's **own Streamlit session** in server memory.
  * No global cache is used for user data.
  * Static file serving is disabled, and downloads are generated on demand for the session
    that owns the data.
  * Touchstone parsing writes the upload to a private temporary directory that is deleted
    immediately afterwards.
* **Delete session data** (sidebar) removes the upload, results and prepared exports.
  **Sign out** also ends the authentication.
* Tracebacks are not shown to remote users (`client.showErrorDetails = false`).
* **Limitation: plain HTTP.** The password and data travel unencrypted. This is acceptable
  only for short tests. For regular use, either
  * tunnel over SSH: `ssh -L 8501:localhost:8501 root@<VPS_IP>`, then open
    `http://localhost:8501` locally (you can then bind the app to `127.0.0.1` with
    `NRW_APP_ADDRESS=127.0.0.1`); or
  * put a TLS reverse proxy (e.g. Caddy or nginx with Let's Encrypt) in front of it.
* A single shared password is not per-user access control.

## 5. Workflow in the browser

**Upload → Validate → Configure → Extract → Diagnose → Export.** The header shows the package
version and Git commit.

1. **Data source.** Upload `.s2p`/`.csv`, or *Load demo data* (clearly labelled **synthetic
   demonstration data**, produced by `examples/make_synthetic_s2p.py` with the independent
   forward model; never measured data).
2. **Data tab.** Shows:
   * file name, frequency range, number of points, spacing (min/median/max), whether the grid
     is uniform, which of S11/S21/S12/S22 are present, and |S11|, |S21| min/max;
   * |S| in dB and phase, wrapped or unwrapped (the unwrapped view is for display only);
   * parser notes.

   Extraction is disabled until validation passes: finite values, strictly increasing
   frequencies, ≥ 5 points, S11 and S21 present.
3. **Configure & extract.** Settings:
   * thickness (mm);
   * cell: coaxial TEM, or rectangular waveguide TE10 with a WR preset or custom a. Frequencies
     at or below the TE10 cutoff are refused, and a warning is shown above the TE20 cutoff;
   * frequency window;
   * method and "known non-magnetic" flag;
   * branch method, and a prior εμ when needed;
   * reference-plane offsets (mm) and optional holder length;
   * advanced: noise level, port averaging, regularisation weight, KK on/off.
4. **Results.**
   * GREEN/YELLOW/RED status with its findings;
   * ε′, ε″, μ′, μ″ plotted separately, with unreliable regions shaded red and attention
     regions amber;
   * loss tangents;
   * an interactive table.
5. **Diagnostics.**
   * warnings and flag counts;
   * branch index m(f) with branch changes marked and tracker disagreements crossed;
   * branch margin and Monte Carlo branch-change fraction;
   * |Γ| and |T|, electrical length βL/π with the half-wavelength lines;
   * Cramér–Rao bounds, solver residual and condition numbers;
   * impedance and propagation constant;
   * the full diagnostics table.
6. **Passivity & causality.**
   * physical checks with their status (exact / assumption / heuristic);
   * ε″ and μ″ **exactly as extracted**, with a ±2σ band. Each frequency is classified *pass*
     (≥ 0), *ambiguous* (negative, within 2σ), *fail* (< −2σ) or *undetermined*;
   * the finite-band KK diagnostic statement, its assumptions, and the residual plot.
7. **Compare methods.** Overlays the selected methods on the same data and settings, with a
   status per method and the relative disagreement between the simultaneous ε/μ methods.
8. **Uncertainty (optional).** Monte Carlo through the whole pipeline.
   * Inputs: |S11|, |S21| (%), phases (°), additive noise, thickness (µm) and plane position
     (µm).
   * 100 / 500 / 1000 / custom draws. The progress bar shows the estimated remaining time.
9. **Export.**
   * `extracted_parameters.csv`, `extracted_parameters.json`, `diagnostics.csv`;
   * `analysis_report.html`, a self-contained report with source file, timestamp, commit,
     thickness, geometry, method, branch method, frequency interval, warnings, diagnostic
     summary and figures;
   * a ZIP with everything, including the input file (optional) and the comparison tables;
   * publication PNGs at 300/600/1200 dpi.

## 6. Supported inputs

| format | details |
|---|---|
| Touchstone `.s2p` (v1) | `# <Hz/kHz/MHz/GHz> S <RI/MA/DB> R <z0>`. Parsed with scikit-rf **and** the builtin reader; the two must agree, and structural errors are rejected (scikit-rf alone silently pads malformed files). Touchstone 2.0 keyword files: scikit-rf only. |
| CSV / TSV / semicolon / whitespace tables | Header row required. Frequency column `f`, `freq`, `frequency` with a unit in the header (`freq_GHz`, `Frequency (Hz)`), or a unit chosen in the UI. Per S-parameter either real/imag (`S11_real`/`S11_imag`, `s11_re`/`s11_im`, `Re(S11)`/`Im(S11)`) or magnitude + phase (`S11_mag` linear or `S11_dB`, with `S11_phase_deg`/`S11_phase_rad`/`S11_phase`; degrees assumed without a unit, and noted). S11 and S21 required; S12/S22 optional. |

The S-parameters must be normalised to the characteristic impedance of the empty sample
holder.

## 7. Geometries, methods, branch strategies

* **Geometries:**
  * coaxial line, TEM (k_c = 0; default 7 mm line);
  * rectangular waveguide, TE10 (k_c = π/a).

  Both use the same slab equations with the correct γ = √(k_c² − k0²εμ) and z = μγ0/γ
  (docs/conventions.md). The waveguide path is covered by unit tests but not by the benchmark
  grid.
* **Methods** (docs/final_method.md, docs/benchmark_report.md):
  * AUTO: the NIST μ = 1 method if the sample is declared non-magnetic, otherwise iterative
    ε & μ;
  * Classical NRW;
  * Iterative ε & μ least squares;
  * Baker-Jarvis / NIST iterative (μ = 1; permittivity only);
  * Non-iterative ε with μ = 1 (Boughriet);
  * Frequency-regularised ε & μ (model assumption; comparison only).

  The two μ = 1 methods do not report μ. Classical NRW is never labelled Baker-Jarvis.
* **Branch strategies:**
  * automatic robust tracking (anchored segment unwrapping, recommended, with an independent
    DP cross-check);
  * dynamic programming;
  * continuous unwrap from the lowest frequency;
  * pointwise Weir group delay;
  * from a prior εμ;
  * principal value.

  No strategy is claimed to guarantee the physical branch.

## 8. Status rules (GREEN / YELLOW / RED)

The rules are implemented in `src/nrw/assessment.py` (rules version 1.0). The status summarises
diagnostics. **It is not an accuracy percentage and does not prove the values correct.**
Wrong thickness, air gaps and calibration errors are invisible to these diagnostics
(docs/benchmark_report.md §6).

Per frequency:
* **unreliable**: non-finite or blow-up values; ill-conditioned by the Cramér–Rao bound
  ("Potential half-wavelength instability"); branch trackers disagree or Monte Carlo draws
  change the branch ("Branch selection uncertain"); |S21| at the noise floor ("Transmission
  magnitude too small for reliable inversion"); no convergence; loss < −2σ.
* **attention**: near βL = kπ; |S11| at the noise floor; weak branch margin or unreliable
  phase step; alternative least-squares minimum; |Γ| ≈ 1; |T| > 1; loss negative but within
  2σ.

Overall:
* **RED** if any of:
  * R1: non-finite values at > 10 % of frequencies;
  * R2: measured S-matrix not passive at > 10 %;
  * R3: branch undetermined at > 10 %;
  * R4: |S21| at the noise floor at > 25 %;
  * R5: ill-conditioned at > 50 %;
  * R6: a μ = 1 method used on a sample that NRW indicates is magnetic;
  * R7: significant passivity violations at > 10 %.
* **YELLOW** if not RED and any of: unreliable or attention points (Y1, Y2); the KK
  diagnostic failed or was unavailable (Y3); a model assumption (reciprocity, port symmetry)
  is violated beyond noise (Y4).
* **GREEN** otherwise.
* **INFO** notes are shown but do not change the level. I1: the filled coaxial line can
  support TE11 in part of the band. This depends on εμ and the line size, not on the data;
  TE11 is excited only by asymmetry or gaps.

## 9. Optional REST API

```bash
export NRW_APP_PASSWORD="$(cat ~/.config/nrw/app_password)"   # or NRW_API_TOKEN
./scripts/run_api.sh                                          # 127.0.0.1:8502 (localhost only)
curl -s http://127.0.0.1:8502/api/v1/health
curl -s -H "Authorization: Bearer $NRW_APP_PASSWORD" \
     -F file=@data/synthetic/synthetic_magnetic_composite_2mm.s2p -F thickness_mm=2.0 \
     -F method=auto -F branch_method=anchored http://127.0.0.1:8502/api/v1/extract | head -c 600
```

Form fields: `thickness_mm` (required), `geometry` (`coax`/`waveguide`), `waveguide_a_mm`,
`method`, `assume_nonmagnetic`, `branch_method`, `f_min_ghz`, `f_max_ghz`, `offset1_mm`,
`offset2_mm`, `noise_sigma`, `freq_unit`, `kk`.

The response is JSON: `metadata` (settings, commit, assessment, warnings, KK statements),
`results` (the same columns as `extracted_parameters.csv`) and `diagnostics`. The interactive
schema is at `/api/v1/docs`. To expose the API remotely, set `NRW_API_ADDRESS=0.0.0.0` (same
HTTP caveat as above).

## 10. Performance

* Extraction runs only when "Run extraction" is pressed. Results are memoised per session,
  keyed by a hash of the data and settings, so revisiting tabs or re-running identical
  settings does not recompute.
* Typical times on this VPS: 801 points in 5–12 s including the KK diagnostic.
* Monte Carlo is separate and optional. It costs roughly 0.6 × the extraction time per draw.

## 11. Known limitations

* Single-sample extraction only in the UI. The two-thickness solver
  (`nrw.extract_multi_thickness`), the recommended remedy for half-wave regions of magnetic
  samples, is available from Python but not yet from the browser.
* Homogeneous, isotropic slab model filling the cross-section. Air gaps are not corrected:
  the result is the effective parameters of the gapped line.
* The UI's reference-plane correction assumes air-filled, lossless offset sections.
* Plain HTTP (see §4). Sessions and results live in server memory and are lost on restart.
