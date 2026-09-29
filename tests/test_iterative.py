import numpy as np
import pytest
from scipy.optimize import least_squares

from conftest import make
from nrw import model
from nrw.classical import nrw
from nrw.conditioning import crlb
from nrw.iterative import (complex_lm, solve_eps_mu_multithickness, solve_eps_mu_pointwise,
                           solve_eps_mu_regularized, solve_eps_only)
from nrw.materials import Debye, Lorentz
from nrw.perturb import complex_gaussian, shift_reference_planes


def test_complex_lm_matches_scipy(coax):
    f = np.array([3e9, 7e9, 11e9])
    eps_t = 5 - 0.3j
    _, _, s11, s21 = model.slab(f, coax, eps_t, 1.0, 8e-3)

    def res(x):
        _, _, m11, m21 = model.slab(f, coax, x[:, 0], 1.0, 8e-3)
        return np.stack([m11 - s11, m21 - s21], axis=1)

    x, cost, nit, conv, cond = complex_lm(res, np.full((3, 1), 4.5 - 0.1j))
    assert conv.all() and np.allclose(x[:, 0], eps_t, rtol=1e-9)
    for k in range(3):
        def fr(p, k=k):
            _, _, m11, m21 = model.slab(f[k:k + 1], coax, p[0] - 1j * p[1], 1.0, 8e-3)
            r = np.concatenate([m11 - s11[k], m21 - s21[k]])
            return np.concatenate([r.real, r.imag])
        sp = least_squares(fr, [4.5, 0.1], method="lm", xtol=1e-14)
        assert np.allclose(sp.x[0] - 1j * sp.x[1], x[k, 0], rtol=1e-8)


def test_eps_only_invariant_reference_plane(coax, band):
    L, eps = 10e-3, 4 - 0.02j
    s = make(band, eps, 1, L, coax)
    # planes displaced by unknown individual amounts but known total air length
    s = shift_reference_planes(band, *s, coax, 3e-3, 1e-3)
    r = solve_eps_only(band, s[0], s[1], L, coax, 3.8, s12=s[2], s22=s[3], l_air_minus_l=4e-3,
                       formulation="eq2.24+2.25")
    assert r.converged.all()
    assert np.allclose(r.eps, eps, rtol=1e-7)
    for form, kw in (("eq2.25", {"l_air_minus_l": 4e-3}), ("eq2.38", {"d1": 2e-3, "d2": 2e-3})):
        r = solve_eps_only(band, s[0], s[1], L, coax, 3.8, s12=s[2], s22=s[3], formulation=form, **kw)
        assert np.allclose(r.eps, eps, rtol=1e-7), form  # eq2.38 with beta=0 needs only d1 + d2


def test_eps_only_stable_at_half_wavelength(coax):
    f = np.linspace(2e9, 18e9, 641)
    L, eps = 20e-3, 4 - 0.004j
    rng = np.random.default_rng(1)
    s = [complex_gaussian(x, 2e-3, rng) for x in make(f, eps, 1, L, coax)]
    base = nrw(f, s[0], s[1], L, coax)
    from nrw.branch import select_branch
    b = select_branch("anchored", f, base.transmission, L, coax)
    r0 = nrw(f, s[0], s[1], L, coax, branch=b.m)
    bj = solve_eps_only(f, s[0], s[1], L, coax, r0.eps_mu_only, s12=s[2], s22=s[3])
    assert np.max(np.abs(bj.eps - eps)) < 0.03  # stays bounded
    assert np.max(np.abs(r0.eps - eps)) > 0.2  # classical simultaneous NRW spikes


def test_iterative_eps_mu_branch_inherited(coax, band):
    """Every branch is an exact solution: LM keeps the branch of its initial guess."""
    L, eps, mu = 10e-3, 10 - 0.5j, 2 - 0.3j
    s = make(band, eps, mu, L, coax)
    base = nrw(band, s[0], s[1], L, coax, branch=0)
    it = solve_eps_mu_pointwise(band, s[0], s[1], L, coax, base.eps, base.mu, multistart=False)
    assert it.converged.all() and np.nanmax(it.residual_rms) < 1e-9
    assert np.max(np.abs(it.eps - eps)) > 1  # wrong branch reproduces data exactly


def test_crlb_blows_up_only_for_simultaneous(coax):
    f = np.linspace(2e9, 18e9, 801)
    L, eps = 20e-3, 4 - 0.004j
    std2, _ = crlb(f, coax, np.full(len(f), eps), np.ones(len(f)), L, 1e-3)
    std1, _ = crlb(f, coax, np.full(len(f), eps), np.ones(len(f)), L, 1e-3, params="eps")
    assert std2[:, 0].max() / np.median(std2[:, 0]) > 50
    assert std1[:, 0].max() / np.median(std1[:, 0]) < 20


def test_multithickness_resolves_half_wavelength(coax):
    f = np.linspace(2e9, 18e9, 401)
    eps, mu = 4 - 0.004j, 1.0
    rng = np.random.default_rng(5)
    samples = []
    for L in (20e-3, 13e-3):
        s = [complex_gaussian(x, 1e-3, rng) for x in make(f, eps, mu, L, coax)]
        samples.append({"s11": s[0], "s21": s[1], "s12": s[2], "s22": s[3], "length_m": L})
    r = solve_eps_mu_multithickness(f, samples, coax, eps * 1.01, 1.0)
    single = solve_eps_mu_pointwise(f, samples[0]["s11"], samples[0]["s21"], 20e-3, coax, eps * 1.01, 1.0,
                                    s12=samples[0]["s12"], s22=samples[0]["s22"], multistart=False)
    assert np.max(np.abs(r.eps - eps)) < 0.25 * np.max(np.abs(single.eps - eps))


def test_regularized_runs_and_improves_singular_points(coax):
    f = np.linspace(2e9, 18e9, 161)
    L = 20e-3
    eps, mu = Debye(4, 0.3, 20e9)(f), Lorentz(1.0, 0.3, 30e9, 20e9)(f)
    rng = np.random.default_rng(7)
    s = [complex_gaussian(x, 1e-3, rng) for x in make(f, eps, mu, L, coax)]
    it = solve_eps_mu_pointwise(f, s[0], s[1], L, coax, eps, mu, s12=s[2], s22=s[3], multistart=False)
    rg = solve_eps_mu_regularized(f, s[0], s[1], L, coax, it.eps, it.mu, s12=s[2], s22=s[3], lam=1.0, sigma=1e-3)
    assert np.max(np.abs(rg.eps - eps)) < np.max(np.abs(it.eps - eps))


def test_eps_only_bad_formulation(coax, band):
    with pytest.raises(ValueError):
        solve_eps_only(band, np.zeros(len(band)), np.ones(len(band)), 1e-3, coax, 2.0, formulation="x")
