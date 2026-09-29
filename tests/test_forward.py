import numpy as np
import pytest

from nrw.forward import Layer, cascade_sparams, coax_air_gap_effective, slab_sparams_closed_form
from nrw.geometry import Coax, RectangularWaveguide, geometry_from_dict, principal_sqrt_decaying
from nrw.materials import Constant, Debye, Lorentz, Sum
from nrw.perturb import shift_reference_planes


@pytest.mark.parametrize("geom,f", [(Coax(), np.linspace(2e9, 18e9, 50)),
                                    (RectangularWaveguide(22.86e-3), np.linspace(8.2e9, 12.4e9, 50))])
def test_cascade_matches_closed_form(geom, f):
    eps, mu = Debye(5, 3, 4e9)(f), Lorentz(1, 1.5, 6e9, 3e9)(f)
    s11, s21, s12, s22 = cascade_sparams(f, geom, [Layer(7e-3, eps, mu)])
    c11, c21 = slab_sparams_closed_form(f, geom, eps, mu, 7e-3)
    assert np.allclose(s11, c11, atol=1e-13) and np.allclose(s21, c21, atol=1e-13)
    assert np.allclose(s12, s21, atol=1e-13)  # reciprocity
    assert np.allclose(s22, s11, atol=1e-13)  # symmetry


def test_air_sample_is_transparent(coax, band):
    s11, s21, _, _ = cascade_sparams(band, coax, [Layer(0.01, 1.0, 1.0)])
    assert np.allclose(s11, 0, atol=1e-14)
    assert np.allclose(s21, np.exp(-coax.gamma0(band) * 0.01), atol=1e-14)


def test_passive_forward_energy(coax, band):
    s11, s21, _, _ = cascade_sparams(band, coax, [Layer(5e-3, 10 - 2j, 3 - 1j)])
    assert np.all(np.abs(s11) ** 2 + np.abs(s21) ** 2 < 1)


def test_lossless_conserves_energy(coax, band):
    s11, s21, _, _ = cascade_sparams(band, coax, [Layer(5e-3, 6.0, 2.0)])
    assert np.allclose(np.abs(s11) ** 2 + np.abs(s21) ** 2, 1, atol=1e-12)


def test_offsets_equal_plane_shift(coax, band):
    eps = 4 - 0.1j
    base = cascade_sparams(band, coax, [Layer(5e-3, eps, 1)])
    off = cascade_sparams(band, coax, [Layer(2e-3, 1, 1), Layer(5e-3, eps, 1), Layer(3e-3, 1, 1)])
    shifted = shift_reference_planes(band, *base, coax, 2e-3, 3e-3)
    for a, b in zip(off, shifted):
        assert np.allclose(a, b, atol=1e-13)


def test_sqrt_branch():
    r = principal_sqrt_decaying(np.array([-4 + 0j, 4 + 0j, -4 - 1e-3j, -4 + 1e-3j]))
    assert np.all(r.real >= 0)
    assert r[0] == 2j


def test_materials_passive_causal_sign(band):
    for m in (Debye(2, 5, 3e9, sigma_dc=0.01), Lorentz(3, 2, 10e9, 1e9), Sum(2, [Debye(2, 5, 3e9), Lorentz(3, 2, 10e9, 1e9)])):
        assert np.all(-m(band).imag >= 0)  # eps'' >= 0 for exp(+jwt)
    assert Constant(3 - 1j)(band).shape == band.shape


def test_air_gap_model_limits():
    g = Coax()
    e, m = coax_air_gap_effective(10 - 1j, 2 - 0.5j, g, g.inner_diameter_m, g.outer_diameter_m)
    assert np.isclose(e, 10 - 1j) and np.isclose(m, 2 - 0.5j)
    e, m = coax_air_gap_effective(10 - 1j, 2 - 0.5j, g, g.inner_diameter_m + 20e-6, g.outer_diameter_m)
    assert e.real < 10  # air gap lowers the apparent permittivity (series capacitance)
    assert m.real < 2
    with pytest.raises(ValueError):
        coax_air_gap_effective(4, 1, g, 1e-3, 8e-3)


def test_geometry_from_dict():
    assert geometry_from_dict({"kind": "coax"}).kc == 0
    assert np.isclose(geometry_from_dict({"kind": "waveguide", "a_m": 0.02286}).cutoff_frequency(), 6.557e9, rtol=1e-3)
    with pytest.raises(ValueError):
        geometry_from_dict({"kind": "waveguide"})
    with pytest.raises(ValueError):
        geometry_from_dict({"kind": "stripline"})
