import numpy as np
import pytest

from conftest import make
from nrw.branch import METHODS, measured_group_delay, select_branch
from nrw.classical import gamma_from_s, nrw
from nrw.materials import Debye, Lorentz
from nrw.perturb import complex_gaussian


def true_branch(f, geom, eps, mu, L):
    beta = geom.gamma(f, eps, mu).imag
    return np.floor((beta * L + np.pi) / (2 * np.pi)).astype(int)


@pytest.mark.parametrize("eps,mu,L", [(4 - 0.1j, 1, 5e-3), (10 - 0.5j, 2 - 0.3j, 10e-3), (25 - 2j, 1, 20e-3)])
def test_nrw_exact_recovery_with_true_branch(coax, band, eps, mu, L):
    s11, s21, _, _ = make(band, eps, mu, L, coax)
    m = true_branch(band, coax, eps, mu, L)
    r = nrw(band, s11, s21, L, coax, branch=m)
    good = np.abs(s11) > 1e-6  # exactly at a lossless resonance S11 = 0
    assert np.allclose(r.eps[good], eps, rtol=1e-8)
    assert np.allclose(r.mu[good], mu, rtol=1e-8)
    assert np.all(np.abs(r.reflection_coefficient) <= 1)


def test_nrw_waveguide(wr90):
    f = np.linspace(8.2e9, 12.4e9, 201)
    eps, mu = Debye(6, 2, 5e9)(f), Lorentz(1, 1, 7e9, 4e9)(f)
    s11, s21, _, _ = make(f, eps, mu, 6e-3, wr90)
    b = select_branch("anchored", f, nrw(f, s11, s21, 6e-3, wr90).transmission, 6e-3, wr90)
    r = nrw(f, s11, s21, 6e-3, wr90, branch=b.m)
    assert np.allclose(r.eps, eps, rtol=1e-8) and np.allclose(r.mu, mu, rtol=1e-8)


def test_principal_branch_fails_for_thick_sample(coax, band):
    s11, s21, _, _ = make(band, 10 - 0.5j, 2 - 0.3j, 10e-3, coax)
    r = nrw(band, s11, s21, 10e-3, coax)  # m = 0 everywhere
    assert np.max(np.abs(r.eps - (10 - 0.5j))) > 5


def test_gamma_root_is_passive():
    X, g, alt, flags = gamma_from_s(np.array([0.3 + 0.2j]), np.array([0.5 - 0.4j]))
    assert np.abs(g[0]) <= 1 <= np.abs(alt[0])
    assert np.isclose(g[0] * alt[0], 1)


def test_singularity_flags_near_half_wavelength(coax):
    f = np.linspace(2e9, 18e9, 1601)
    L, eps = 20e-3, 4 - 1e-4j
    s11, s21, _, _ = make(f, eps, 1, L, coax)
    r = nrw(f, s11, s21, L, coax, branch=true_branch(f, coax, eps, 1, L), noise_floor=1e-3)
    assert r.flags["near_half_wavelength"].any()
    assert r.flags["s11_below_noise"].any()
    # flagged S11 minima coincide with integer beta L / pi
    k = np.argmin(np.abs(s11))
    assert abs(r.electrical_length_over_pi[k] - round(r.electrical_length_over_pi[k])) < 0.01


def test_nrw_does_not_hide_nan(coax):
    f = np.array([2e9, 3e9, 4e9])
    r = nrw(f, np.zeros(3, complex), np.ones(3, complex) * 0.9, 1e-3, coax)
    assert r.flags["nonfinite"].all() or r.flags["blowup"].any() or not np.all(np.isfinite(r.eps))


@pytest.mark.parametrize("method", ["unwrap", "anchored", "dp", "group_delay"])
def test_branch_methods_noise_free(coax, band, method):
    eps, mu, L = 10 - 0.5j, 2 - 0.3j, 10e-3
    s11, s21, _, _ = make(band, eps, mu, L, coax)
    T = nrw(band, s11, s21, L, coax).transmission
    b = select_branch(method, band, T, L, coax)
    assert np.array_equal(b.m, true_branch(band, coax, eps, mu, L))


def test_unwrap_fails_when_first_point_not_short(coax):
    """'continuous unwrap from low frequency' assumes beta L < pi at f_min."""
    f = np.linspace(8e9, 18e9, 401)
    eps, L = 30 - 0.3j, 15e-3  # beta L / pi ~ 4.4 at 8 GHz
    s11, s21, _, _ = make(f, eps, 1, L, coax)
    T = nrw(f, s11, s21, L, coax).transmission
    mt = true_branch(f, coax, eps, 1, L)
    assert not np.array_equal(select_branch("unwrap", f, T, L, coax).m, mt)
    assert np.array_equal(select_branch("anchored", f, T, L, coax).m, mt)


def test_unwrap_fails_on_coarse_grid_dp_survives(coax):
    """Phase step per point > pi breaks unwrapping; continuity of eps*mu does not need it."""
    f = np.linspace(2e9, 18e9, 17)  # 1 GHz steps
    eps, L = 16 - 0.2j, 40e-3  # tau ~ 0.53 ns -> 0.53 cycles per 1 GHz step
    s11, s21, _, _ = make(f, eps, 1, L, coax)
    T = nrw(f, s11, s21, L, coax).transmission
    mt = true_branch(f, coax, eps, 1, L)
    assert not np.array_equal(select_branch("unwrap", f, T, L, coax).m, mt)
    assert np.array_equal(select_branch("dp", f, T, L, coax).m, mt)


def test_prior_branch(coax, band):
    eps, L = 10 - 0.5j, 10e-3
    s11, s21, _, _ = make(band, eps, 1, L, coax)
    T = nrw(band, s11, s21, L, coax).transmission
    b = select_branch("prior", band, T, L, coax, prior_eps_mu=11.0)
    assert np.array_equal(b.m, true_branch(band, coax, eps, 1, L))
    with pytest.raises(ValueError):
        select_branch("prior", band, T, L, coax)


def test_unknown_branch_method(coax, band):
    with pytest.raises(ValueError):
        select_branch("magic", band, np.ones(len(band), complex), 1e-3, coax)
    assert "anchored" in METHODS


def test_measured_group_delay_matches_nL_over_c(coax, band):
    L, n = 10e-3, 3.0
    s11, s21, _, _ = make(band, n**2, 1, L, coax)
    T = nrw(band, s11, s21, L, coax).transmission
    tau, _ = measured_group_delay(band, T)
    assert np.allclose(tau[1:-1], n * L / 299792458.0, rtol=1e-6)


def test_anchored_robust_to_noise(coax, band):
    rng = np.random.default_rng(3)
    eps, mu, L = 10 - 0.5j, 2 - 0.3j, 10e-3
    s = [complex_gaussian(x, 3e-3, rng) for x in make(band, eps, mu, L, coax)]
    T = nrw(band, s[0], s[1], L, coax).transmission
    b = select_branch("anchored", band, T, L, coax)
    assert np.mean(b.m == true_branch(band, coax, eps, mu, L)) > 0.99
