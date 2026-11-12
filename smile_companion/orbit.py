"""SMILE orbit: TLE handling, propagation, and frame conversions.

Propagation uses SGP4/SDP4 (``sgp4`` package) when a TLE is available, and a
two-body Kepler propagator otherwise (offline demo, or no ``sgp4`` installed).
SGP4 returns TEME coordinates; for a dashboard-level picture TEME is treated
as inertial (the difference from J2000 is well under a degree).
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd

from .config import MU_KM3_S2, NOMINAL_ORBIT, RE_KM

try:
    from sgp4.api import Satrec, jday
    HAVE_SGP4 = True
except ImportError:  # pragma: no cover - exercised only without the package
    HAVE_SGP4 = False

DEG = np.pi / 180.0


# ---------------------------------------------------------------------------
# Time helpers
# ---------------------------------------------------------------------------
def julian_date(t: pd.DatetimeIndex) -> np.ndarray:
    t = pd.DatetimeIndex(t).tz_localize(None) if pd.DatetimeIndex(t).tz else pd.DatetimeIndex(t)
    return t.to_julian_date().to_numpy()


def gmst_rad(jd: np.ndarray) -> np.ndarray:
    """Greenwich mean sidereal time (IAU 1982), radians."""
    T = (jd - 2451545.0) / 36525.0
    g = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T ** 2 - T ** 3 / 38710000.0
    return np.mod(g, 360.0) * DEG


# ---------------------------------------------------------------------------
# Orbital elements and TLEs
# ---------------------------------------------------------------------------
@dataclass
class Elements:
    a_km: float
    e: float
    inc_deg: float
    raan_deg: float
    argp_deg: float
    M0_deg: float
    epoch: datetime

    @property
    def period_h(self) -> float:
        return 2 * np.pi * np.sqrt(self.a_km ** 3 / MU_KM3_S2) / 3600.0

    @property
    def perigee_alt_km(self) -> float:
        return self.a_km * (1 - self.e) - RE_KM

    @property
    def apogee_alt_km(self) -> float:
        return self.a_km * (1 + self.e) - RE_KM


def nominal_elements(epoch: datetime | None = None) -> Elements:
    o = NOMINAL_ORBIT
    rp, ra = RE_KM + o["perigee_alt_km"], RE_KM + o["apogee_alt_km"]
    return Elements(a_km=(rp + ra) / 2, e=(ra - rp) / (ra + rp), inc_deg=o["inclination_deg"],
                    raan_deg=o["raan_deg"], argp_deg=o["argp_deg"], M0_deg=o["mean_anomaly_deg"],
                    epoch=epoch or datetime.now(timezone.utc).replace(microsecond=0))


def parse_tle(text: str) -> tuple[str, str, str]:
    """Return (name, line1, line2) from a 2- or 3-line TLE block."""
    lines = [ln.rstrip() for ln in text.strip().splitlines() if ln.strip()]
    l1 = next(ln for ln in lines if ln.startswith("1 "))
    l2 = next(ln for ln in lines if ln.startswith("2 "))
    name = lines[0].strip() if not lines[0].startswith(("1 ", "2 ")) else "UNKNOWN"
    return name, l1, l2


def tle_elements(line1: str, line2: str) -> Elements:
    """Mean elements from a TLE (used by the Kepler fallback)."""
    yy, doy = int(line1[18:20]), float(line1[20:32])
    year = 2000 + yy if yy < 57 else 1900 + yy
    epoch = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1)
    n_rev_day = float(line2[52:63])
    n_rad_s = n_rev_day * 2 * np.pi / 86400.0
    return Elements(a_km=(MU_KM3_S2 / n_rad_s ** 2) ** (1 / 3), e=float("0." + line2[26:33].strip()),
                    inc_deg=float(line2[8:16]), raan_deg=float(line2[17:25]),
                    argp_deg=float(line2[34:42]), M0_deg=float(line2[43:51]), epoch=epoch)


def fetch_tle(norad_id: int, timeout: float = 15.0) -> str:
    import requests
    from .config import URLS, USER_AGENT
    r = requests.get(URLS["tle"].format(norad=norad_id), timeout=timeout,
                     headers={"User-Agent": USER_AGENT})
    r.raise_for_status()
    if "1 " not in r.text:
        raise ValueError(f"No TLE returned for NORAD {norad_id}: {r.text[:80]!r}")
    return r.text


# ---------------------------------------------------------------------------
# Propagation
# ---------------------------------------------------------------------------
def _solve_kepler(M: np.ndarray, e: float) -> np.ndarray:
    E = np.where(e > 0.8, np.pi * np.ones_like(M), M)
    for _ in range(50):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < 1e-12:
            break
    return E


def propagate_kepler(el: Elements, times: pd.DatetimeIndex) -> np.ndarray:
    """Two-body inertial positions [km], shape (N, 3)."""
    epoch = pd.Timestamp(el.epoch).tz_convert("UTC") if pd.Timestamp(el.epoch).tz else pd.Timestamp(el.epoch, tz="UTC")
    tt = pd.DatetimeIndex(times)
    tt = tt.tz_convert("UTC") if tt.tz else tt.tz_localize("UTC")
    dt = (tt - epoch).total_seconds().to_numpy()
    n = np.sqrt(MU_KM3_S2 / el.a_km ** 3)
    E = _solve_kepler(np.mod(el.M0_deg * DEG + n * dt, 2 * np.pi), el.e)
    nu = 2 * np.arctan2(np.sqrt(1 + el.e) * np.sin(E / 2), np.sqrt(1 - el.e) * np.cos(E / 2))
    r = el.a_km * (1 - el.e * np.cos(E))
    # perifocal -> inertial (3-1-3 rotation by RAAN, i, argp)
    O, i, w = el.raan_deg * DEG, el.inc_deg * DEG, el.argp_deg * DEG
    u = w + nu
    x = r * (np.cos(O) * np.cos(u) - np.sin(O) * np.sin(u) * np.cos(i))
    y = r * (np.sin(O) * np.cos(u) + np.cos(O) * np.sin(u) * np.cos(i))
    z = r * (np.sin(u) * np.sin(i))
    return np.column_stack([x, y, z])


def propagate_tle(line1: str, line2: str, times: pd.DatetimeIndex) -> np.ndarray:
    """SGP4/SDP4 positions [km] (TEME), falling back to Kepler without ``sgp4``."""
    if not HAVE_SGP4:
        return propagate_kepler(tle_elements(line1, line2), times)
    sat = Satrec.twoline2rv(line1, line2)
    tt = pd.DatetimeIndex(times)
    tt = tt.tz_convert("UTC").tz_localize(None) if tt.tz else tt
    jd, fr = np.array([jday(t.year, t.month, t.day, t.hour, t.minute,
                             t.second + t.microsecond * 1e-6) for t in tt]).T
    err, r, _ = sat.sgp4_array(jd, fr)
    r = np.asarray(r, float)
    r[np.asarray(err) != 0] = np.nan
    return r


# ---------------------------------------------------------------------------
# Frames
# ---------------------------------------------------------------------------
def sun_unit_vector(jd: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Low-precision solar direction in the equatorial inertial frame
    (Astronomical Almanac, ~0.01 deg) and the obliquity [rad]."""
    n = jd - 2451545.0
    L = (280.460 + 0.9856474 * n) * DEG
    g = (357.528 + 0.9856003 * n) * DEG
    lam = L + (1.915 * np.sin(g) + 0.020 * np.sin(2 * g)) * DEG
    eps = (23.439 - 0.0000004 * n) * DEG
    s = np.column_stack([np.cos(lam), np.cos(eps) * np.sin(lam), np.sin(eps) * np.sin(lam)])
    return s, eps


def eci_to_gse(r_eci: np.ndarray, jd: np.ndarray) -> np.ndarray:
    """Inertial -> GSE: X toward the Sun, Z ecliptic north, Y = Z x X (dusk)."""
    X, eps = sun_unit_vector(jd)
    Z = np.column_stack([np.zeros_like(eps), -np.sin(eps), np.cos(eps)])
    Y = np.cross(Z, X)
    return np.column_stack([(r_eci * X).sum(1), (r_eci * Y).sum(1), (r_eci * Z).sum(1)])


def eci_to_geodetic(r_eci: np.ndarray, jd: np.ndarray) -> pd.DataFrame:
    """Sub-satellite latitude/longitude (spherical Earth) and altitude."""
    th = gmst_rad(jd)
    x = np.cos(th) * r_eci[:, 0] + np.sin(th) * r_eci[:, 1]
    y = -np.sin(th) * r_eci[:, 0] + np.cos(th) * r_eci[:, 1]
    z = r_eci[:, 2]
    rn = np.linalg.norm(r_eci, axis=1)
    return pd.DataFrame({"lat_deg": np.degrees(np.arcsin(z / rn)),
                         "lon_deg": np.degrees(np.arctan2(y, x)),
                         "alt_km": rn - RE_KM})


def track(times: pd.DatetimeIndex, tle: tuple[str, str] | None = None,
          elements: Elements | None = None) -> pd.DataFrame:
    """Full SMILE track table: inertial, GSE (in R_E) and geodetic columns."""
    times = pd.DatetimeIndex(times)
    if tle is not None:
        r = propagate_tle(tle[0], tle[1], times)
    else:
        r = propagate_kepler(elements or nominal_elements(), times)
    jd = julian_date(times)
    gse = eci_to_gse(r, jd) / RE_KM
    geo = eci_to_geodetic(r, jd)
    df = pd.DataFrame({"time": times, "x_eci_km": r[:, 0], "y_eci_km": r[:, 1], "z_eci_km": r[:, 2],
                       "x_gse_re": gse[:, 0], "y_gse_re": gse[:, 1], "z_gse_re": gse[:, 2]})
    df["r_re"] = np.linalg.norm(r, axis=1) / RE_KM
    return pd.concat([df, geo], axis=1)
