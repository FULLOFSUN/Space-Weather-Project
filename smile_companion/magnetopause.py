"""Magnetopause and bow-shock models driven by upstream solar wind.

Shue et al. (1998), JGR 103(A8), 17691:
    r(theta) = r0 * (2 / (1 + cos theta))^alpha
    r0    = (10.22 + 1.29 tanh(0.184 (Bz + 8.14))) * Dp^(-1/6.6)     [R_E]
    alpha = (0.58 - 0.007 Bz) * (1 + 0.024 ln Dp)
with Bz [nT] (GSM) and dynamic pressure Dp [nPa].
"""
from __future__ import annotations

import numpy as np


def dynamic_pressure_npa(n_cc, v_kms):
    """Solar-wind dynamic pressure [nPa], OMNI convention (4% alpha particles):
    P = 2e-6 * n * v^2 with n in cm^-3 and v in km/s."""
    return 2.0e-6 * np.asarray(n_cc, float) * np.asarray(v_kms, float) ** 2


def shue_params(bz_nt, dp_npa):
    bz = np.asarray(bz_nt, float)
    dp = np.clip(np.asarray(dp_npa, float), 0.05, None)
    r0 = (10.22 + 1.29 * np.tanh(0.184 * (bz + 8.14))) * dp ** (-1.0 / 6.6)
    alpha = (0.58 - 0.007 * bz) * (1.0 + 0.024 * np.log(dp))
    return r0, alpha


def standoff_re(bz_nt, dp_npa):
    """Subsolar magnetopause distance r0 [R_E]."""
    return shue_params(bz_nt, dp_npa)[0]


def surface(bz_nt: float, dp_npa: float, theta_max_deg: float = 150.0, n: int = 300):
    """Magnetopause curve in a plane containing the Sun-Earth line (X, rho) [R_E]."""
    r0, alpha = shue_params(bz_nt, dp_npa)
    th = np.radians(np.linspace(-theta_max_deg, theta_max_deg, n))
    r = r0 * (2.0 / (1.0 + np.cos(th))) ** alpha
    return r * np.cos(th), r * np.sin(th)


def radius_at(theta_rad, bz_nt, dp_npa):
    r0, alpha = shue_params(bz_nt, dp_npa)
    denom = np.maximum(1.0 + np.cos(theta_rad), 1e-12)   # theta = 180 deg -> very large, not inf
    return r0 * (2.0 / denom) ** alpha


def position_relative_to_magnetopause(x_re, y_re, z_re, bz_nt, dp_npa,
                                      max_theta_deg: float = 120.0):
    """Signed distance [R_E] of a point from the magnetopause along its radial
    direction: positive = outside (magnetosheath / solar wind).

    Returns NaN beyond `max_theta_deg` from the Sun-Earth line, where the model
    surface flares into the tail and a radial distance is not meaningful."""
    x, y, z = (np.asarray(v, float) for v in (x_re, y_re, z_re))
    r = np.sqrt(x * x + y * y + z * z)
    theta = np.arccos(np.clip(x / np.where(r > 0, r, 1), -1, 1))
    d = r - radius_at(theta, bz_nt, dp_npa)
    return np.where(np.degrees(theta) <= max_theta_deg, d, np.nan)


def bow_shock_standoff_re(r0_mp, mach_ms: float = 8.0, gamma: float = 5 / 3):
    """Farris & Russell (1994) subsolar bow-shock distance from the magnetopause
    standoff, for a magnetosonic Mach number M."""
    m2 = mach_ms ** 2
    return r0_mp * (1 + 1.1 * ((gamma - 1) * m2 + 2) / ((gamma + 1) * (m2 - 1)))
