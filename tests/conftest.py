import numpy as np
import pytest

from nrw.forward import SampleSpec
from nrw.geometry import Coax, RectangularWaveguide


@pytest.fixture
def coax():
    return Coax()


@pytest.fixture
def wr90():
    return RectangularWaveguide(a_m=22.86e-3)


@pytest.fixture
def band():
    return np.linspace(2e9, 18e9, 401)


def make(f, eps, mu, L, geom, **kw):
    return SampleSpec(eps=eps, mu=mu, length_m=L, geometry=geom, **kw).sparams(f)
