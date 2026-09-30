"""Tests for the web-facing helpers (nrw.web.*), the flexible I/O and the assessment rules."""
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pytest

from conftest import make
from nrw.assessment import assess, passivity_classes
from nrw.io import describe, load_measurement_bytes, parse_csv_text, write_touchstone
from nrw.materials import Debye, Lorentz
from nrw.perturb import complex_gaussian

pytest.importorskip("pandas")
DEMO = Path(__file__).resolve().parents[1] / "data" / "synthetic"


# ------------------------------------------------------------------------------------------ security
def test_password_only_from_env(monkeypatch):
    from nrw.web import security
    monkeypatch.delenv("NRW_APP_PASSWORD", raising=False)
    assert security.password_problem() is not None
    assert not security.check_password("anything")
    monkeypatch.setenv("NRW_APP_PASSWORD", "short")
    assert "shorter" in security.password_problem()
    monkeypatch.setenv("NRW_APP_PASSWORD", "long-enough-password")
    assert security.password_problem() is None
    assert security.check_password("long-enough-password")
    assert not security.check_password("long-enough-passwor")


def test_login_lockout(monkeypatch):
    from nrw.web import security
    security._reset_for_tests()
    monkeypatch.setenv("NRW_APP_PASSWORD", "long-enough-password")
    for _ in range(security.MAX_FAILURES):
        ok, _ = security.attempt_login("bad", "1.2.3.4", sleep=False)
        assert not ok
    ok, msg = security.attempt_login("long-enough-password", "1.2.3.4", sleep=False)
    assert not ok and "Too many" in msg  # locked even with the right password
    ok, _ = security.attempt_login("long-enough-password", "5.6.7.8", sleep=False)
    assert ok  # other clients unaffected
    security._reset_for_tests()


# ------------------------------------------------------------------------------------------ input
@pytest.fixture
def sdata(coax):
    f = np.linspace(2e9, 18e9, 81)
    s = make(f, Debye(6, 4, 5e9)(f), Lorentz(1, 2.5, 4e9, 6e9)(f), 2e-3, coax)
    return f, s


def _csv(header, rows):
    return header + "\n" + "\n".join(",".join(f"{v:.12g}" for v in r) for r in rows)


def test_csv_real_imag_and_db_phase_agree(sdata):
    f, s = sdata
    ri = _csv("freq_GHz,S11_real,S11_imag,S21_real,S21_imag",
              [[f[k] / 1e9, s[0][k].real, s[0][k].imag, s[1][k].real, s[1][k].imag] for k in range(len(f))])
    dbp = _csv("Frequency (Hz);|S11| (dB);S11 phase (deg);S21_dB;S21_phase_deg".replace(",", ";"),
               [[f[k], 20 * np.log10(abs(s[0][k])), np.degrees(np.angle(s[0][k])), 20 * np.log10(abs(s[1][k])),
                 np.degrees(np.angle(s[1][k]))] for k in range(len(f))]).replace(",", ";")
    mag_rad = _csv("f_mhz,S11_mag,S11_phase_rad,S21_mag,S21_phase_rad,S12_mag,S12_phase_rad,S22_mag,S22_phase_rad",
                   [[f[k] / 1e6] + sum(([abs(x[k]), np.angle(x[k])] for x in s), []) for k in range(len(f))])
    a = parse_csv_text(ri, "a.csv")
    b = parse_csv_text(dbp, "b.csv")
    c = parse_csv_text(mag_rad, "c.csv")
    for d in (a, b, c):
        assert np.allclose(d.frequency, f, rtol=1e-10)
        assert np.allclose(d.s11, s[0], atol=1e-9) and np.allclose(d.s21, s[1], atol=1e-9)
    assert a.measured == ("s11", "s21") and a.s12_measured() is None
    assert c.measured == ("s11", "s21", "s12", "s22") and np.allclose(c.s22_measured(), s[3], atol=1e-9)


def test_csv_errors_are_explicit():
    with pytest.raises(ValueError, match="no unit"):
        parse_csv_text("freq,S11_re,S11_im,S21_re,S21_im\n1,0,0,1,0\n2,0,0,1,0", "x.csv")
    d = parse_csv_text("freq,S11_re,S11_im,S21_re,S21_im\n1,0,0,1,0\n2,0,0,1,0", "x.csv", freq_unit="GHz")
    assert d.frequency[1] == 2e9 and any("user setting" in n for n in d.notes)
    with pytest.raises(ValueError, match="no frequency column"):
        parse_csv_text("x,S11_re,S11_im\n1,0,0\n", "x.csv")
    with pytest.raises(ValueError, match="S21"):
        parse_csv_text("f_hz,S11_re,S11_im\n1,0,0\n2,0,0\n", "x.csv")
    with pytest.raises(ValueError, match="complete pair"):
        parse_csv_text("f_hz,S11_re,S11_im,S21_re\n1,0,0,1\n", "x.csv")
    with pytest.raises(ValueError, match="non-numeric"):
        parse_csv_text("f_hz,S11_re,S11_im,S21_re,S21_im\n1,a,0,1,0\n", "x.csv")


def test_load_bytes_s2p_uses_skrf_and_crosschecks(tmp_path, sdata):
    f, s = sdata
    p = tmp_path / "m.s2p"
    write_touchstone(p, f, *s)
    d = load_measurement_bytes("m.s2p", p.read_bytes())
    assert np.allclose(d.s21, s[1]) and np.allclose(d.s12, s[2])
    pytest.importorskip("skrf")
    assert any("agrees" in n for n in d.notes)
    with pytest.raises(ValueError, match="unsupported"):
        load_measurement_bytes("m.xlsx", b"123")
    with pytest.raises(ValueError):
        load_measurement_bytes("bad.s2p", b"# GHz S RI R 50\n1 2 3\n")
    summ = describe(d)
    assert summ["n_points"] == 81 and summ["uniform_grid"] and summ["strictly_increasing"]


def test_demo_manifest_complete():
    man = json.loads((DEMO / "manifest.json").read_text())
    assert len(man) >= 4
    for name, meta in man.items():
        assert (DEMO / name).exists() and meta["synthetic"] and meta["thickness_mm"] > 0
        assert "SYNTHETIC" in (DEMO / name).read_text()[:300]


# ------------------------------------------------------------------------------------------ analysis
@pytest.fixture(scope="module")
def composite_outcome():
    from nrw.web.analysis import AnalysisConfig, run_analysis
    d = load_measurement_bytes("c.s2p", (DEMO / "synthetic_magnetic_composite_2mm.s2p").read_bytes())
    return run_analysis(d, AnalysisConfig(thickness_mm=2.0, f_min_ghz=2.0, f_max_ghz=10.0))


def test_config_validation():
    from nrw.web.analysis import AnalysisConfig
    assert AnalysisConfig(thickness_mm=0).validate()
    assert AnalysisConfig(thickness_mm=1, geometry="waveguide").validate()
    assert AnalysisConfig(thickness_mm=1, branch_method="prior").validate()
    assert AnalysisConfig(thickness_mm=1, f_min_ghz=5, f_max_ghz=4).validate()
    assert AnalysisConfig(thickness_mm=1, method="magic").validate()
    assert not AnalysisConfig(thickness_mm=1).validate()
    a, b = AnalysisConfig(thickness_mm=1), AnalysisConfig(thickness_mm=1.1)
    assert a.key("h") != b.key("h") and a.key("h") == AnalysisConfig(thickness_mm=1).key("h")


def test_outcome_tables_and_signs(composite_outcome):
    o = composite_outcome
    t = o.table
    assert t["frequency_GHz"].min() >= 2.0 and t["frequency_GHz"].max() <= 10.0 + 1e-9
    for c in ("frequency_Hz", "frequency_GHz", "epsilon_real", "epsilon_imag", "mu_real", "mu_imag", "tan_delta_e",
              "tan_delta_m", "branch_index", "passivity_flag", "causality_metric_eps", "method", "status", "warnings"):
        assert c in t.columns
    assert np.allclose(t["epsilon_imag"], -t["epsilon_loss"]) and (t["epsilon_loss"] > 0).all()
    assert np.allclose(t["tan_delta_e"], t["epsilon_loss"] / t["epsilon_real"])
    assert o.assessment.level == "GREEN"
    assert {"gamma_refl_abs", "T_abs", "beta_L_over_pi", "flag_singular_region"} <= set(o.diagnostics.columns)


def test_eps_only_method_does_not_report_mu(composite_outcome):
    from nrw.web.analysis import AnalysisConfig, run_analysis
    d = load_measurement_bytes("c.s2p", (DEMO / "synthetic_magnetic_composite_2mm.s2p").read_bytes())
    o = run_analysis(d, AnalysisConfig(thickness_mm=2.0, method="baker_jarvis", f_max_ghz=6, kk=False))
    assert o.table["mu_real"].isna().all()
    assert o.assessment.level == "RED" and any(r["rule"] == "R6" for r in o.assessment.reasons)


def test_assessment_red_when_branch_undetermined(coax):
    from nrw import extract_material_parameters
    f = np.linspace(2e9, 18e9, 17)  # aliased grid: trackers disagree
    s = make(f, 10 - 0.1j, 1, 60e-3, coax)
    r = extract_material_parameters(f, s[0], s[1], 60e-3, s12=s[2], s22=s[3], method="nrw_eps", kk=False)
    a = assess(r)
    assert a.level == "RED" and any(x["rule"] == "R3" for x in a.reasons)


def test_assessment_half_wave_is_flagged_not_hidden():
    from nrw.web.analysis import AnalysisConfig, run_analysis
    d = load_measurement_bytes("h.s2p", (DEMO / "synthetic_lowloss_magnetic_20mm_halfwave.s2p").read_bytes())
    o = run_analysis(d, AnalysisConfig(thickness_mm=20.0, kk=False))
    assert o.assessment.level in ("YELLOW", "RED")
    assert (o.table["status"] == "unreliable").sum() > 0
    assert len(o.table) == len(d.frequency)  # nothing removed


def test_passivity_classes():
    c = passivity_classes(np.array([0.1, -0.01, -1.0, np.nan]), np.array([0.1, 0.1, 0.1, 0.1]))
    assert list(c) == ["pass", "ambiguous", "fail", "undetermined"]


# ------------------------------------------------------------------------------------------ plots/export
def test_plots_and_exports(composite_outcome):
    from nrw.web import export as X
    from nrw.web import plots as P
    o = composite_outcome
    for k in ("eps_real", "eps_loss", "mu_real", "mu_loss"):
        P.parameter_figure(o, k)
        P.passivity_figure(o, k) if k.endswith("loss") else None
    for w in ("gamma_T", "electrical_length", "identifiability", "solver", "impedance", "propagation"):
        P.diagnostic_figure(o, w)
    P.branch_figure(o)
    P.kk_figures(o, "eps")
    figs = P.publication_figures(o, dpi=300)
    assert {"eps_real.png", "mu_loss.png", "parameters_panel.png"} <= set(figs)
    assert all(v[:8] == b"\x89PNG\r\n\x1a\n" for v in figs.values())
    js = json.loads(X.results_json(o))
    assert js["metadata"]["software"]["commit"] and "exp(+jωt)" in js["metadata"]["convention"]
    html = X.report_html(o, figs).decode()
    assert "not a mathematical proof" in html
    # the phrase may only appear inside the software's own negation
    assert "strictly satisf" not in html.replace("does NOT establish that the material strictly satisfies", "")
    z = zipfile.ZipFile(io.BytesIO(X.zip_bundle(o, figs, "in.s2p", b"abc")))
    names = set(z.namelist())
    assert {"extracted_parameters.csv", "extracted_parameters.json", "diagnostics.csv", "analysis_report.html",
            "input/in.s2p", "figures/eps_real.png"} <= names


# ------------------------------------------------------------------------------------------ REST API
def test_rest_api_auth_and_extract(monkeypatch):
    pytest.importorskip("fastapi")
    pytest.importorskip("httpx")
    from fastapi.testclient import TestClient

    from nrw.web import rest, security
    security._reset_for_tests()
    monkeypatch.setattr(security, "FAIL_DELAY_S", 0.0)
    monkeypatch.setenv("NRW_APP_PASSWORD", "long-enough-password")
    monkeypatch.delenv("NRW_API_TOKEN", raising=False)
    c = TestClient(rest.app)
    assert c.get("/api/v1/health").json()["status"] == "ok"
    content = (DEMO / "synthetic_magnetic_composite_2mm.s2p").read_bytes()
    form = {"thickness_mm": "2.0", "f_max_ghz": "6", "kk": "false"}
    r = c.post("/api/v1/extract", files={"file": ("a.s2p", content)}, data=form)
    assert r.status_code == 401
    r = c.post("/api/v1/extract", files={"file": ("a.s2p", content)}, data=form,
               headers={"Authorization": "Bearer long-enough-password"})
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["metadata"]["assessment"]["level"] in ("GREEN", "YELLOW", "RED")
    assert len(body["results"]["epsilon_real"]) == len(body["results"]["frequency_Hz"])
    r = c.post("/api/v1/extract", files={"file": ("a.s2p", content)}, data={"thickness_mm": "-1"},
               headers={"Authorization": "Bearer long-enough-password"})
    assert r.status_code == 422
    security._reset_for_tests()


# ------------------------------------------------------------------------------------------ UI smoke
def test_streamlit_app_gate_and_demo_flow(monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    from nrw.web import security
    security._reset_for_tests()
    app = str(Path(__file__).resolve().parents[1] / "app.py")
    monkeypatch.delenv("NRW_APP_PASSWORD", raising=False)
    at = AppTest.from_file(app, default_timeout=300).run()
    assert at.error and "locked" in at.error[0].value and not at.file_uploader  # nothing usable
    monkeypatch.setenv("NRW_APP_PASSWORD", "long-enough-password")
    monkeypatch.setattr(security, "FAIL_DELAY_S", 0.0)
    at = AppTest.from_file(app, default_timeout=300).run()
    assert not at.tabs  # no data pages before login
    at.text_input[0].input("long-enough-password")
    at.button[0].click()
    at.run()
    assert at.session_state["authenticated"]
    at.radio[0].set_value("Synthetic demonstration data").run()
    [b for b in at.button if b.label == "Load demo data"][0].click()
    at.run()
    assert "data" in at.session_state and not at.exception
    [b for b in at.button if b.label == "Run extraction"][0].click()
    at.run()
    assert not at.exception and at.session_state["current_key"]
    [b for b in at.button if b.label == "Delete session data"][0].click()
    at.run()
    assert "data" not in at.session_state
    security._reset_for_tests()


def test_streamlit_compare_uncertainty_export_flow(monkeypatch):
    """Regression: Results and Uncertainty tabs show the same memoised figures (needs unique keys);
    the whole compare -> Monte Carlo -> export path must run without exceptions."""
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    from nrw.web import security
    security._reset_for_tests()
    monkeypatch.setenv("NRW_APP_PASSWORD", "long-enough-password")
    monkeypatch.setattr(security, "FAIL_DELAY_S", 0.0)
    at = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=600).run()
    at.text_input[0].input("long-enough-password")
    at.button[0].click()
    at.run()
    at.radio[0].set_value("Synthetic demonstration data").run()
    [b for b in at.button if b.label == "Load demo data"][0].click()
    at.run()
    [b for b in at.button if b.label == "Run extraction"][0].click()
    at.run()
    [b for b in at.button if b.label == "Compare extraction methods"][0].click()
    at.run()
    assert not at.exception and "comparison" in at.session_state
    [s for s in at.selectbox if s.label == "Monte Carlo draws"][0].set_value("custom")
    [n for n in at.number_input if n.label.startswith("custom number")][0].set_value(20)
    [b for b in at.button if b.label == "Run uncertainty analysis"][0].click()
    at.run()
    assert not at.exception and "mc_key" in at.session_state
    [b for b in at.button if b.label == "Prepare export files"][0].click()
    at.run()
    assert not at.exception
    files = at.session_state["exports"][1]
    z = zipfile.ZipFile(io.BytesIO(files["nrw_results.zip"]))
    assert "figures/parameters_panel.png" in z.namelist() and "comparison/nrw.csv" in z.namelist()
    t = files["extracted_parameters.csv"].decode().splitlines()[0]
    assert "eps_real_p2.5" in t  # Monte Carlo percentiles exported
    security._reset_for_tests()


# ------------------------------------------------------------------------------------------ two-thickness
PAIR = ("synthetic_lowloss_magnetic_20mm_halfwave.s2p", "synthetic_lowloss_magnetic_13mm_halfwave.s2p")


def _pair():
    return [load_measurement_bytes(n, (DEMO / n).read_bytes()) for n in PAIR]


def test_two_thickness_analysis_beats_single_samples():
    from nrw.web.analysis import AnalysisConfig, run_analysis, single_sample_outcomes
    d1, d2 = _pair()
    cfg = AnalysisConfig(thickness_mm=20.0, two_thickness=True, thickness2_mm=13.0, kk=False)
    o = run_analysis(d1, cfg, data2=d2)
    assert o.result.method == "multi_thickness_eps_mu" and o.result.sample_results
    singles = single_sample_outcomes(o)
    assert len(singles) == 2
    eps_t, mu_t = 4 - 0.004j, 1.5 - 0.003j
    err_joint = np.nanmax(np.abs(o.result.epsilon_complex - eps_t))
    for s in singles.values():
        assert np.nanmax(np.abs(s.result.epsilon_complex - eps_t)) > 5 * err_joint
        assert s.result.flags["singular_region"].sum() > o.result.flags["singular_region"].sum()
    assert np.nanmax(np.abs(o.result.mu_complex - mu_t)) < 0.05
    assert {"branch_index_sample1", "branch_index_sample2", "ill_conditioned_sample1_alone"} <= set(o.table.columns)
    assert (o.assessment.point_status == "unreliable").sum() == 0


def test_two_thickness_refuses_mismatched_grid_and_equal_thickness():
    from dataclasses import replace

    from nrw.web.analysis import AnalysisConfig, run_analysis
    d1, d2 = _pair()
    d2b = replace(d2, frequency=d2.frequency[::2], s11=d2.s11[::2], s21=d2.s21[::2], s12=d2.s12[::2],
                  s22=d2.s22[::2])
    with pytest.raises(ValueError, match="same frequency grid"):
        run_analysis(d1, AnalysisConfig(thickness_mm=20.0, two_thickness=True, thickness2_mm=13.0, kk=False),
                     data2=d2b)
    with pytest.raises(ValueError, match="different thicknesses"):
        run_analysis(d1, AnalysisConfig(thickness_mm=20.0, two_thickness=True, thickness2_mm=20.0), data2=d2)
    with pytest.raises(ValueError, match="no second sample"):
        run_analysis(d1, AnalysisConfig(thickness_mm=20.0, two_thickness=True, thickness2_mm=13.0))


def test_two_thickness_monte_carlo_and_export():
    from nrw.web import export as X
    from nrw.web.analysis import AnalysisConfig, run_analysis
    d1, d2 = _pair()
    cfg = AnalysisConfig(thickness_mm=20.0, two_thickness=True, thickness2_mm=13.0, f_max_ghz=6.0, kk=False,
                         uncertainty={"s21_phase_deg": 0.5, "thickness_m": 5e-6, "n_draws": 6, "seed": 3})
    calls = []
    o = run_analysis(d1, cfg, data2=d2, progress=lambda i, n: calls.append(i))
    u = o.result.uncertainty
    assert calls == list(range(1, 7)) and u["n_draws"] == 6 and u["n_failed_draws"] == 0
    assert np.all(u["eps_real"]["p2.5"] <= u["eps_real"]["p97.5"])
    md = X.metadata(o)
    assert md["two_thickness"] and np.allclose(md["thicknesses_m"], [0.02, 0.013])
    z = zipfile.ZipFile(io.BytesIO(X.zip_bundle(o, {}, "a.s2p", b"1", extra_inputs={"b.s2p": b"2"})))
    assert {"input/a.s2p", "input/b.s2p"} <= set(z.namelist())
    assert "sample 2: 13 mm" in X.report_html(o, {}).decode()


def test_rest_api_two_thickness(monkeypatch):
    pytest.importorskip("fastapi")
    pytest.importorskip("httpx")
    from fastapi.testclient import TestClient

    from nrw.web import rest, security
    security._reset_for_tests()
    monkeypatch.setattr(security, "FAIL_DELAY_S", 0.0)
    monkeypatch.setenv("NRW_APP_PASSWORD", "long-enough-password")
    monkeypatch.delenv("NRW_API_TOKEN", raising=False)
    c = TestClient(rest.app)
    h = {"Authorization": "Bearer long-enough-password"}
    files = {"file": (PAIR[0], (DEMO / PAIR[0]).read_bytes()), "file2": (PAIR[1], (DEMO / PAIR[1]).read_bytes())}
    r = c.post("/api/v1/extract", files=files, data={"thickness_mm": "20", "thickness2_mm": "13", "f_max_ghz": "6",
                                                     "kk": "false"}, headers=h)
    assert r.status_code == 200, r.text
    assert r.json()["metadata"]["method_used"] == "multi_thickness_eps_mu"
    r = c.post("/api/v1/extract", files=files, data={"thickness_mm": "20"}, headers=h)
    assert r.status_code == 422  # file2 without thickness2_mm
    security._reset_for_tests()


def test_streamlit_two_thickness_flow(monkeypatch):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    from nrw.web import security
    security._reset_for_tests()
    monkeypatch.setenv("NRW_APP_PASSWORD", "long-enough-password")
    monkeypatch.setattr(security, "FAIL_DELAY_S", 0.0)
    at = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=900).run()
    at.text_input[0].input("long-enough-password")
    at.button[0].click()
    at.run()
    at.radio[0].set_value("Synthetic demonstration data").run()
    at.selectbox[0].set_value(PAIR[0]).run()
    [b for b in at.button if b.label == "Load demo data"][0].click()
    at.run()
    [b for b in at.button if b.label.startswith("Load paired synthetic demo sample")][0].click()
    at.run()
    assert "data2" in at.session_state and not at.exception
    assert [c for c in at.checkbox if c.label.startswith("Use the joint two-thickness")][0].value
    [b for b in at.button if b.label == "Run extraction"][0].click()
    at.run()
    assert not at.exception
    key = at.session_state["current_key"]
    assert at.session_state["results"][key].result.method == "multi_thickness_eps_mu"
    [b for b in at.button if b.label == "Prepare export files"][0].click()
    at.run()
    assert not at.exception
    z = zipfile.ZipFile(io.BytesIO(at.session_state["exports"][1]["nrw_results.zip"]))
    assert {f"input/{PAIR[0]}", f"input/{PAIR[1]}"} <= set(z.namelist())
    [b for b in at.button if b.label == "Remove second sample"][0].click()
    at.run()
    assert "data2" not in at.session_state and "data" in at.session_state
    security._reset_for_tests()


def test_two_thickness_detects_inconsistent_samples(coax):
    from nrw import extract_multi_thickness
    f = np.linspace(2e9, 18e9, 161)
    rng = np.random.default_rng(4)

    def meas(eps, L):
        s = [complex_gaussian(x, 1e-3, rng) for x in make(f, eps, 1.5 - 0.003j, L, coax)]
        return dict(s11=s[0], s21=s[1], s12=s[2], s22=s[3], thickness=L)

    ok = extract_multi_thickness(f, [meas(4 - 0.004j, 20e-3), meas(4 - 0.004j, 13e-3)], kk=False)
    bad = extract_multi_thickness(f, [meas(4 - 0.004j, 20e-3), meas(4.2 - 0.004j, 13e-3)], kk=False)
    assert ok.flags["samples_inconsistent"].mean() < 0.05
    assert not any(r["rule"] == "R8" for r in assess(ok).reasons)
    a = assess(bad)
    assert bad.flags["samples_inconsistent"].mean() > 0.5
    assert a.level == "RED" and any(r["rule"] == "R8" for r in a.reasons)
    assert any("not consistent with one" in w for w in bad.warnings)
