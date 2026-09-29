"""Causal, passive dispersion models used as synthetic ground truth.

All models use exp(+jwt): eps = eps' - j eps'' with eps'' >= 0 for passive media.
Every model below is a sum of terms that are analytic in the lower half of the complex
omega plane (poles in the upper half plane for exp(+jwt)), hence causal, and has a
non-negative loss part for real omega > 0, hence passive.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np


class Model:
    """A frequency-dependent complex relative parameter (eps_r or mu_r)."""

    def __call__(self, f: np.ndarray) -> np.ndarray:  # pragma: no cover - abstract
        raise NotImplementedError

    @property
    def high_frequency_limit(self) -> float:  # value at f -> infinity
        raise NotImplementedError  # pragma: no cover


@dataclass(frozen=True)
class Constant(Model):
    """Frequency-independent value. Note: a constant *loss* is not exactly KK-consistent
    over an infinite band (it is the classical counter-example); use it only for
    extraction tests, not for causality tests."""

    value: complex

    def __call__(self, f):
        return np.full(np.shape(f), complex(self.value))

    @property
    def high_frequency_limit(self):
        return complex(self.value).real


@dataclass(frozen=True)
class Debye(Model):
    """eps(f) = eps_inf + delta / (1 + j f/f_relax)  (+ optional DC conductivity term)."""

    eps_inf: float
    delta: float
    f_relax: float
    sigma_dc: float = 0.0  # S/m, adds -j sigma/(w eps0); passive, causal, singular at DC

    def __call__(self, f):
        f = np.asarray(f, dtype=float)
        val = self.eps_inf + self.delta / (1 + 1j * f / self.f_relax)
        if self.sigma_dc:
            eps0 = 8.8541878128e-12
            val = val - 1j * self.sigma_dc / (2 * np.pi * f * eps0)
        return val

    @property
    def high_frequency_limit(self):
        return self.eps_inf


@dataclass(frozen=True)
class Lorentz(Model):
    """eps(f) = eps_inf + delta f0^2 / (f0^2 - f^2 + j f gamma_f).

    ``gamma_f`` is the damping (full width) in Hz. With exp(+jwt) the loss part is
    delta f0^2 f gamma_f / ((f0^2-f^2)^2 + f^2 gamma_f^2) >= 0.
    """

    eps_inf: float
    delta: float
    f0: float
    gamma_f: float

    def __call__(self, f):
        f = np.asarray(f, dtype=float)
        return self.eps_inf + self.delta * self.f0**2 / (self.f0**2 - f**2 + 1j * f * self.gamma_f)

    @property
    def high_frequency_limit(self):
        return self.eps_inf


@dataclass(frozen=True)
class Sum(Model):
    """eps_inf + sum of the dispersive parts of several models."""

    eps_inf: float
    terms: Sequence[Model] = field(default_factory=tuple)

    def __call__(self, f):
        out = np.full(np.shape(f), complex(self.eps_inf))
        for t in self.terms:
            out = out + (t(f) - t.high_frequency_limit)
        return out

    @property
    def high_frequency_limit(self):
        return self.eps_inf


@dataclass(frozen=True)
class Callable_(Model):
    fn: Callable[[np.ndarray], np.ndarray]
    limit: float = 1.0

    def __call__(self, f):
        return np.asarray(self.fn(np.asarray(f, dtype=float)), dtype=complex)

    @property
    def high_frequency_limit(self):
        return self.limit


def as_model(x) -> Model:
    if isinstance(x, Model):
        return x
    if callable(x):
        return Callable_(x)
    return Constant(complex(x))
