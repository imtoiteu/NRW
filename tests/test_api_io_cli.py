import json

import numpy as np
import pytest

from conftest import make
from nrw import extract_material_parameters
from nrw.api import InputError
from nrw.cli import main
from nrw.io import read_csv_sparams, read_touchstone, write_touchstone
from nrw.materials import Debye, Lorentz
from nrw.perturb import complex_gaussian


@pytest.fixture
def magnetic(coax):
    f = np.linspace(2e9, 18e9, 201)
    eps, mu = Debye(6, 4, 5e9)(f), Lorentz(1, 2.5, 4e9, 6e9)(f)
    rng = np.random.default_rng(11)
    s = [complex_gaussian(x, 5e-4, rng) for x in make(f, eps, mu, 2e-3, coax)]
    return f, s, eps, mu


def test_api_recovers_magnetic(magnetic):
    f, s, eps, mu = magnetic
    r = extract_material_parameters(f, s[0], s[1], 2e-3, s12=s[2], s22=s[3], method="iterative_eps_mu",
                                    noise_sigma=5e-4)
    good = ~r.singular_region
    assert good.mean() > 0.9
    assert np.median(np.abs(r.epsilon_complex[good] - eps[good]) / np.abs(eps[good])) < 0.01
    assert np.median(np.abs(r.mu_complex[good] - mu[good]) / np.abs(mu[good])) < 0.01
    assert r.causality["eps"].passed
    assert set(["frequency", "epsilon_complex", "mu_complex"]).issubset(r.__dict__)
    assert np.allclose(r.epsilon_loss, -r.epsilon_complex.imag)


@pytest.mark.parametrize("method", ["nrw", "nrw_eps", "baker_jarvis", "iterative_eps_mu", "regularized_eps_mu"])
def test_api_methods_run(coax, method):
    f = np.linspace(2e9, 18e9, 81)
    s = make(f, 4 - 0.05j, 1, 6e-3, coax)
    r = extract_material_parameters(f, s[0], s[1], 6e-3, s12=s[2], s22=s[3], method=method, kk=False)
    assert np.nanmedian(np.abs(r.epsilon_complex - (4 - 0.05j))) < 1e-3


def test_api_auto_nonmagnetic_uses_bj(coax):
    f = np.linspace(2e9, 18e9, 81)
    s = make(f, 4 - 0.05j, 1, 6e-3, coax)
    r = extract_material_parameters(f, s[0], s[1], 6e-3, assume_nonmagnetic=True)
    assert r.method == "baker_jarvis"
    assert any("mu_r = 1" in w for w in r.warnings)


def test_api_warns_bj_on_magnetic(magnetic):
    f, s, *_ = magnetic
    r = extract_material_parameters(f, s[0], s[1], 2e-3, method="baker_jarvis", kk=False)
    assert any("mu = 1 assumption" in w for w in r.warnings)


def test_api_reference_plane_deembedding(coax):
    from nrw.perturb import shift_reference_planes
    f = np.linspace(2e9, 18e9, 81)
    s = shift_reference_planes(f, *make(f, 4 - 0.05j, 1.5 - 0.1j, 6e-3, coax), coax, 2e-3, 5e-3)
    r = extract_material_parameters(f, s[0], s[1], 6e-3, s12=s[2], s22=s[3], method="nrw",
                                    reference_plane_offsets=(2e-3, 5e-3), kk=False)
    assert np.allclose(r.epsilon_complex, 4 - 0.05j, rtol=1e-6)
    rb = extract_material_parameters(f, s[0], s[1], 6e-3, s12=s[2], s22=s[3], method="baker_jarvis",
                                     holder_length=13e-3, kk=False)
    assert rb.metadata["holder_length_m"] == 13e-3


@pytest.mark.parametrize("bad", ["nonmonotonic", "zero_thickness", "negative_thickness", "nan", "shape",
                                 "below_cutoff", "method"])
def test_api_rejects_malformed(coax, bad):
    f = np.linspace(2e9, 18e9, 21)
    s11, s21, _, _ = make(f, 4, 1, 5e-3, coax)
    kw = dict(frequency=f, s11=s11, s21=s21, thickness=5e-3)
    if bad == "nonmonotonic":
        kw["frequency"] = f[::-1].copy()
    elif bad == "zero_thickness":
        kw["thickness"] = 0.0
    elif bad == "negative_thickness":
        kw["thickness"] = -1e-3
    elif bad == "nan":
        s11 = s11.copy()
        s11[3] = np.nan
        kw["s11"] = s11
    elif bad == "shape":
        kw["s21"] = s21[:-1]
    elif bad == "below_cutoff":
        kw["geometry"] = {"kind": "waveguide", "a_m": 22.86e-3}
    elif bad == "method":
        kw["method"] = "nist_magic"
    with pytest.raises(InputError):
        extract_material_parameters(**kw)


def test_touchstone_roundtrip_and_formats(tmp_path, coax):
    f = np.linspace(2e9, 18e9, 11)
    s = make(f, 4 - 0.1j, 1, 5e-3, coax)
    p = tmp_path / "a.s2p"
    write_touchstone(p, f, s[0], s[1], s[2] * 0.99, s[3] * 0.98, comments=["test"])
    d = read_touchstone(p)
    assert np.allclose(d.frequency, f) and np.allclose(d.s11, s[0]) and np.allclose(d.s21, s[1])
    assert np.allclose(d.s12, s[2] * 0.99) and np.allclose(d.s22, s[3] * 0.98)  # S21 before S12 ordering
    # MA and DB in GHz, with inline comments
    lines = ["! header", "# GHz S MA R 50"]
    lines2 = ["# GHz S DB R 50"]
    for k in range(len(f)):
        vals = [f[k] / 1e9]
        vals2 = [f[k] / 1e9]
        for x in (s[0][k], s[1][k], s[2][k], s[3][k]):
            vals += [abs(x), np.degrees(np.angle(x))]
            vals2 += [20 * np.log10(abs(x)), np.degrees(np.angle(x))]
        lines.append(" ".join(map(str, vals)) + " ! c")
        lines2.append(" ".join(map(str, vals2)))
    (tmp_path / "ma.s2p").write_text("\n".join(lines))
    (tmp_path / "db.s2p").write_text("\n".join(lines2))
    for name in ("ma.s2p", "db.s2p"):
        d = read_touchstone(tmp_path / name)
        assert np.allclose(d.s21, s[1], rtol=1e-6) and np.allclose(d.frequency, f)
    skrf = pytest.importorskip("skrf")  # noqa: F841
    d2 = read_touchstone(p, backend="skrf")
    assert np.allclose(d2.s12, d.s12 if False else s[2] * 0.99)


def test_touchstone_rejects_bad(tmp_path):
    (tmp_path / "bad.s2p").write_text("# GHz S RI R 50\n1 2 3\n")
    with pytest.raises(ValueError):
        read_touchstone(tmp_path / "bad.s2p")
    (tmp_path / "y.s2p").write_text("# GHz Y RI R 50\n" + " ".join(["1"] * 9) + "\n")
    with pytest.raises(ValueError):
        read_touchstone(tmp_path / "y.s2p")
    (tmp_path / "x.txt").write_text("")
    with pytest.raises(ValueError):
        read_touchstone(tmp_path / "x.txt")


def test_csv_reader(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("f_hz,s11_re,s11_im,s21_re,s21_im\n1e9,0.1,0.2,0.5,-0.1\n2e9,0.2,0.1,0.4,-0.2\n")
    d = read_csv_sparams(p)
    assert np.allclose(d.s12, d.s21) and d.frequency[1] == 2e9


def test_cli_end_to_end(tmp_path, coax):
    f = np.linspace(2e9, 18e9, 61)
    s = make(f, Debye(4, 1, 8e9)(f), 1, 5e-3, coax)
    p = tmp_path / "m.s2p"
    write_touchstone(p, f, *s)
    out = tmp_path / "out"
    rc = main(["--input", str(p), "--thickness-mm", "5", "--method", "baker-jarvis", "--output", str(out),
               "--no-plots"])
    assert rc == 0
    summ = json.loads((out / "summary.json").read_text())
    assert summ["method"] == "baker_jarvis"
    data = np.genfromtxt(out / "results.csv", delimiter=",", names=True)
    assert np.allclose(data["eps_real"], Debye(4, 1, 8e9)(f).real, rtol=1e-6)


def test_result_save_with_plots(tmp_path, magnetic):
    f, s, *_ = magnetic
    r = extract_material_parameters(f, s[0], s[1], 2e-3, s12=s[2], s22=s[3], noise_sigma=5e-4)
    r.save(tmp_path / "o")
    for name in ("results.csv", "summary.json", "parameters.png", "diagnostics.png", "kk_diagnostic.png"):
        assert (tmp_path / "o" / name).exists()
