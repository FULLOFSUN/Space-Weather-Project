#!/usr/bin/env python3
"""Quick-look figures without the dashboard.

    python scripts/make_figures.py --demo                       # offline, synthetic
    python scripts/make_figures.py --start 2023-04-23 --end 2023-04-25   # real OMNI + Kp
    python scripts/make_figures.py --start 2026-11-01 --end 2026-11-03 --live-tle
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from smile_companion import figures, orbit, spaceweather  # noqa: E402
from smile_companion.config import SMILE_NORAD_ID  # noqa: E402
from smile_companion.demo import synthetic_storm  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true", help="synthetic data, no network")
    ap.add_argument("--start")
    ap.add_argument("--end")
    ap.add_argument("--live-tle", action="store_true", help="use SMILE's latest CelesTrak TLE")
    ap.add_argument("--out", default="docs/figures")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    if a.demo:
        sw = synthetic_storm()
        kp = sw.attrs["kp"]
        note = " — SYNTHETIC DEMO"
    else:
        if not (a.start and a.end):
            ap.error("give --start and --end, or use --demo")
        sw = spaceweather.fetch_omni(a.start, a.end)
        kp = spaceweather.fetch_kp(a.start, a.end)
        note = f" — {sw.attrs['source']}"

    tle = None
    if a.live_tle:
        _, l1, l2 = orbit.parse_tle(orbit.fetch_tle(SMILE_NORAD_ID))
        tle = (l1, l2)
        el = orbit.tle_elements(l1, l2)
    else:
        el = orbit.nominal_elements(epoch=sw.index[0].to_pydatetime())
    times = pd.date_range(sw.index[0], periods=400,
                          freq=pd.Timedelta(hours=el.period_h) / 399)
    trk = orbit.track(times, tle=tle, elements=el)

    quiet_bz, quiet_dp = 2.0, 1.5
    i_min = sw["r0_mp_re"].idxmin()
    storm_bz, storm_dp = float(sw.loc[i_min, "bz_nt"]), float(sw.loc[i_min, "pdyn_npa"])
    orbit_note = "" if tle else " — nominal orbit, arbitrary phase"
    f1 = figures.orbit_vs_magnetopause(trk, i_min, (quiet_bz, quiet_dp), (storm_bz, storm_dp), orbit_note)
    f1.savefig(out / "orbit_magnetopause.png", dpi=150)
    f2 = figures.solar_wind_panel(sw, kp, note)
    f2.savefig(out / "solar_wind.png", dpi=150)
    print(f"period {el.period_h:.1f} h, perigee {el.perigee_alt_km:,.0f} km, apogee {el.apogee_alt_km:,.0f} km")
    print(f"minimum standoff {sw['r0_mp_re'].min():.2f} R_E at {i_min:%Y-%m-%d %H:%M} UTC")
    print("wrote", out / "orbit_magnetopause.png", "and", out / "solar_wind.png")


if __name__ == "__main__":
    main()
