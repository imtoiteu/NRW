"""nrw: robust extraction of complex permittivity and permeability from two-port
transmission/reflection S-parameters.

Conventions: exp(+jwt) time dependence, eps_r = eps' - j eps'', mu_r = mu' - j mu''
(passive media: eps'', mu'' >= 0). See docs/conventions.md.
"""
__version__ = "0.1.0"

from .api import ExtractionResult, extract_material_parameters  # noqa: E402
from .geometry import Coax, RectangularWaveguide  # noqa: E402

__all__ = ["extract_material_parameters", "ExtractionResult", "Coax", "RectangularWaveguide", "__version__"]
