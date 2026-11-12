"""Offline tests: physics, frames and response parsers. Run with  pytest -q"""
import json

import numpy as np
import pandas as pd
import pytest

from smile_companion import magnetopause as mp
from smile_companion import orbit, spaceweather
from smile_companion.config import RE_KM


def _checksum(line: str) -> str:
    s = sum(int(c) if c.isdigit() else (1 if c == "-" else 0) for c in line[:68])
    return line[:68] + str(s % 10)


def synthetic_tle(inc=70.0, raan=10.0, ecc=0.8363, argp=270.0, M=180.0, n=0.4743):
    """A TLE built from SMILE's NOMINAL elements, for tests only (not real data)."""
    l1 = _checksum("1 69123U 26109A   26307.50000000  .00000000  00000-0  00000-0 0  999 ")
    l2 = _checksum(f"2 69123 {inc:8.4f} {raan:8.4f} {f'{ecc:.7f}'[2:]} {argp:8.4f} {M:8.4f} {n:11.8f}    1 ")
    return f"SMILE\n{l1}\n{l2}\n"


# ---------------------------------------------------------------- magnetopause
def test_shue_quiet_standoff_about_10_re():
    assert float(mp.standoff_re(0.0, 2.0)) == pytest.approx(10.25, abs=0.05)


def test_magnetopause_compresses_with_pressure_and_southward_bz():
    assert mp.standoff_re(0, 10) < mp.standoff_re(0, 2)
    assert mp.standoff_re(-10, 2) < mp.standoff_re(5, 2)


def test_dynamic_pressure_omni_convention():
    assert float(mp.dynamic_pressure_npa(5, 400)) == pytest.approx(1.6)


def test_point_outside_dayside_magnetopause():
    assert mp.position_relative_to_magnetopause(15, 0, 0, 0, 2) > 0
    assert mp.position_relative_to_magnetopause(5, 0, 0, 0, 2) < 0
    assert np.isnan(mp.position_relative_to_magnetopause(-15, 0, 0, 0, 2))


# ---------------------------------------------------------------- orbit
def test_nominal_orbit_geometry():
    el = orbit.nominal_elements()
    assert el.perigee_alt_km == pytest.approx(5000, abs=1)
    assert el.apogee_alt_km == pytest.approx(121182, abs=1)
    assert 49 < el.period_h < 52


def test_kepler_track_reaches_perigee_and_apogee():
    el = orbit.nominal_elements(epoch=pd.Timestamp("2026-11-01", tz="UTC").to_pydatetime())
    t = pd.date_range("2026-11-01", periods=2000, freq=pd.Timedelta(hours=el.period_h) / 1999, tz="UTC")
    trk = orbit.track(t, elements=el)
    assert trk.alt_km.min() == pytest.approx(5000, rel=0.01)
    assert trk.alt_km.max() == pytest.approx(121182, rel=0.001)
    assert trk.lat_deg.abs().max() <= 70.0 + 1e-6
    # perigee over the southern hemisphere for argp = 270 deg
    assert trk.loc[trk.alt_km.idxmin(), "lat_deg"] < -60


def test_tle_parsing_and_elements():
    name, l1, l2 = orbit.parse_tle(synthetic_tle())
    el = orbit.tle_elements(l1, l2)
    assert name == "SMILE"
    assert el.inc_deg == pytest.approx(70.0)
    assert el.e == pytest.approx(0.8363)
    assert el.perigee_alt_km == pytest.approx(5000, rel=0.01)
    assert el.apogee_alt_km == pytest.approx(121182, rel=0.01)
    assert el.epoch.year == 2026 and el.epoch.month == 11


def test_tle_track_runs_with_or_without_sgp4():
    _, l1, l2 = orbit.parse_tle(synthetic_tle())
    t = pd.date_range("2026-11-03", periods=50, freq="1h", tz="UTC")
    trk = orbit.track(t, tle=(l1, l2))
    assert trk[["x_gse_re", "z_gse_re", "alt_km"]].notna().all().all()


def test_gse_x_points_to_sun():
    jd = orbit.julian_date(pd.DatetimeIndex(["2026-06-21 12:00"]))
    s, _ = orbit.sun_unit_vector(jd)
    gse = orbit.eci_to_gse(s * 1e8, jd)
    assert gse[0, 0] == pytest.approx(1e8, rel=1e-9)
    # summer solstice: Sun near declination +23.4 deg
    assert np.degrees(np.arcsin(s[0, 2])) == pytest.approx(23.44, abs=0.1)


# ---------------------------------------------------------------- parsers
SWPC_PLASMA = json.dumps([["time_tag", "density", "speed", "temperature"],
                          ["2026-11-03 10:00:00.000", "5.2", "412.1", "90000"],
                          ["2026-11-03 10:01:00.000", None, "415.0", "91000"]])
SWPC_MAG = json.dumps([["time_tag", "bx_gsm", "by_gsm", "bz_gsm", "lon_gsm", "lat_gsm", "bt"],
                       ["2026-11-03 10:00:00.000", "1.0", "-2.0", "-4.5", "300", "-30", "5.0"]])
OMNI_CSV = ("2023-04-23T18:00:00.000Z,550.2,12.3,-15.1,7.2,-120\n"
            "2023-04-23T18:01:00.000Z,99999.9,999.99,9999.99,99.99,99999\n")
GFZ_KP = json.dumps({"meta": {}, "datetime": ["2023-04-23T18:00:00Z", "2023-04-23T21:00:00Z"],
                     "Kp": [8.333, 7.667], "status": ["def", "def"]})


def test_parse_swpc():
    p = spaceweather.parse_swpc_table(SWPC_PLASMA)
    m = spaceweather.parse_swpc_table(SWPC_MAG)
    assert p["speed"].iloc[0] == pytest.approx(412.1) and np.isnan(p["density"].iloc[1])
    assert m["bz_gsm"].iloc[0] == pytest.approx(-4.5)
    assert str(p.index.tz) == "UTC"


def test_parse_omni_fill_values():
    df = spaceweather.parse_omni_hapi_csv(OMNI_CSV)
    assert df["flow_speed"].iloc[0] == pytest.approx(550.2)
    assert df.iloc[1].isna().all()


def test_parse_gfz_kp_and_storm_level():
    kp = spaceweather.parse_gfz_kp(GFZ_KP)
    assert kp["kp"].iloc[0] == pytest.approx(8.333)
    assert spaceweather.storm_level(kp["kp"].iloc[0]) == "G4 severe"
    assert spaceweather.storm_level(3.0) == "quiet / unsettled"
