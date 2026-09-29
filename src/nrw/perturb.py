"""Measurement-error models used to stress-test extraction algorithms.

All functions take an explicit ``numpy.random.Generator`` so experiments are seeded.
None of these are used by the extraction code itself.
"""
from __future__ import annotations

import numpy as np


def complex_gaussian(s, sigma, rng):
    """Additive circular complex Gaussian noise, E|n|^2 = sigma^2 (VNA trace-noise-like)."""
    s = np.asarray(s, dtype=complex)
    return s + sigma / np.sqrt(2) * (rng.standard_normal(s.shape) + 1j * rng.standard_normal(s.shape))


def magnitude_noise(s, rel_sigma, rng):
    """Multiplicative magnitude noise: |S| -> |S| (1 + rel_sigma N(0,1))."""
    s = np.asarray(s, dtype=complex)
    return s * (1 + rel_sigma * rng.standard_normal(s.shape))


def phase_noise(s, sigma_rad, rng):
    s = np.asarray(s, dtype=complex)
    return s * np.exp(1j * sigma_rad * rng.standard_normal(s.shape))


def shift_reference_planes(f, s11, s21, s12, s22, geom, d1_m, d2_m):
    """Move reference planes by d1, d2 (m) *away* from the sample in the empty guide.

    S11 *= R1^2, S22 *= R2^2, S21, S12 *= R1 R2 with R_i = exp(-gamma0 d_i).
    Negative d moves the plane into the sample side (pure phase de-embedding).
    """
    g0 = geom.gamma0(f)
    r1, r2 = np.exp(-g0 * d1_m), np.exp(-g0 * d2_m)
    return s11 * r1**2, s21 * r1 * r2, s12 * r1 * r2, s22 * r2**2
