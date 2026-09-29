import numpy as np
import pytest

from conftest import make
from nrw.kk import kk_consistency, kk_real_from_imag_finite
from nrw.materials import Debye, Lorentz
from nrw.physics import check_material, check_network, continuity_check
from nrw.uncertainty import UncertaintySpec, monte_carlo


def test_kk_integral_exact_over_wide_band():
    m = Lorentz(3, 2, 10e9, 2e9)
    f = np.linspace(1e6, 2e12, 200001)
    fe = np.array([2e9, 5e9, 15e9])
    pred = 3 + kk_real_from_imag_finite(f, -m(f).imag, fe)
    assert np.allclose(pred, m(fe).real, atol=2e-3)


def test_kk_sign_convention_debye():
    # Debye with exp(+jwt): eps'' > 0 and eps' - eps_inf = (2/pi) P int w' eps''/(w'^2 - w^2)
    m = Debye(2, 5, 3e9)
    f = np.geomspace(1e3, 1e15, 400001)
    pred = 2 + kk_real_from_imag_finite(f, -m(f).imag, np.array([3e9]))
    assert np.isclose(pred[0], m(3e9).real, rtol=1e-3)


@pytest.mark.parametrize("model", [Lorentz(3, 2, 10e9, 2e9), Debye(2, 8, 1e9), Debye(2, 8, 40e9),
                                   Debye(2, 8, 6e9, sigma_dc=0.05), Lorentz(1, 4, 3e9, 2e9)])
def test_kk_passes_causal_models(band, model):
    r = kk_consistency(band, model(band))
    assert r.passed, r.max_abs_normalized_residual
    assert "does NOT establish" in r.statement


def test_kk_fails_noncausal(band):
    mixed = Debye(2, 8, 1e9)(band).real - 1j * (-Lorentz(3, 2, 10e9, 2e9)(band).imag)
    assert not kk_consistency(band, mixed).passed
    active = np.conj(Lorentz(3, 2, 10e9, 2e9)(band))
    assert not kk_consistency(band, active).passed


def test_kk_input_validation(band):
    with pytest.raises(ValueError):
        kk_consistency(band[:4], np.ones(4, complex))
    with pytest.raises(ValueError):
        kk_consistency(band[::-1], Debye(2, 3, 1e9)(band))


def test_network_checks(coax, band):
    s = make(band, 5 - 0.5j, 1, 5e-3, coax)
    checks = {c.name: c for c in check_network(*s)}
    assert checks["network_passivity"].n_violations == 0
    assert checks["reciprocity"].n_violations == 0
    bad = check_network(s[0] * 1.5, s[1] * 1.5, s[2], s[3])
    assert {c.name: c for c in bad}["network_passivity"].n_violations > 0
    nonrec = {c.name: c for c in check_network(s[0], s[1], s[2] * 1.3, s[3])}
    assert nonrec["reciprocity"].status == "assumption" and nonrec["reciprocity"].n_violations > 0


def test_material_checks_statuses():
    eps = np.array([4 - 0.1j, 4 + 0.1j, 0.5 - 0.1j])
    mu = np.array([1 - 0j, 1 - 0.01j, 1 + 0.2j])
    c = {x.name: x for x in check_material(eps, mu)}
    assert c["eps_passivity"].violated.tolist() == [False, True, False]
    assert c["mu_passivity"].violated.tolist() == [False, False, True]
    assert c["eps_real_ge_1"].status == "heuristic"
    # with uncertainty the small violation is tolerated
    c2 = {x.name: x for x in check_material(eps, mu, eps_sigma=np.array([0.1, 0.1, 0.1]))}
    assert c2["eps_passivity"].n_violations == 0


def test_continuity_detects_spike(band):
    x = np.full(len(band), 4 - 0.1j) + 1e-4 * np.arange(len(band))
    x[200] += 3
    assert continuity_check(band, x, "eps").violated[200]


def test_monte_carlo_seeded_and_summaries(coax):
    f = np.linspace(2e9, 18e9, 41)
    s = make(f, 4 - 0.1j, 1, 5e-3, coax)

    def ex(ff, a11, a21, a12, a22, L):
        from nrw.classical import nrw
        r = nrw(ff, a11, a21, L, coax)
        return r.eps_mu_only, np.ones(len(ff)), r.branch

    spec = UncertaintySpec(s21_phase_deg=0.5, thickness_m=5e-6, n_draws=30, seed=1)
    a = monte_carlo(ex, f, *s, 5e-3, coax, spec, nominal_branch=np.zeros(len(f), int))
    b = monte_carlo(ex, f, *s, 5e-3, coax, spec, nominal_branch=np.zeros(len(f), int))
    assert np.array_equal(a["eps_real"]["median"], b["eps_real"]["median"])
    assert np.all(a["eps_real"]["p2.5"] <= a["eps_real"]["p97.5"])
    assert np.all(a["eps_real"]["sigma_robust"] > 0)
