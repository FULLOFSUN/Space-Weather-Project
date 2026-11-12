"""Upstream solar wind (NOAA SWPC real-time, NASA OMNI) and Kp (GFZ Potsdam).

Every fetcher has a matching ``parse_*`` function that works on the raw
response text, so parsing is unit-tested without network access.
"""
from __future__ import annotations

import io
import json
from datetime import datetime

import numpy as np
import pandas as pd

from .config import OMNI_DATASET, OMNI_FILL, OMNI_PARAMS, URLS, USER_AGENT
from .magnetopause import dynamic_pressure_npa, standoff_re


def _get(url: str, timeout: float = 30.0) -> str:
    import requests
    r = requests.get(url, timeout=timeout, headers={"User-Agent": USER_AGENT})
    r.raise_for_status()
    return r.text


def _iso(t) -> str:
    return pd.Timestamp(t).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# NOAA SWPC real-time solar wind (DSCOVR/ACE at L1), last 2 h ... 7 days
# Files are JSON arrays whose first row is the header.
# ---------------------------------------------------------------------------
def parse_swpc_table(text: str) -> pd.DataFrame:
    rows = json.loads(text)
    df = pd.DataFrame(rows[1:], columns=rows[0])
    df["time_tag"] = pd.to_datetime(df["time_tag"], utc=True)
    for c in df.columns.drop("time_tag"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.set_index("time_tag").sort_index()


def fetch_swpc(window: str = "1-day") -> pd.DataFrame:
    """window in {'2-hour', '6-hour', '1-day', '3-day', '7-day'}."""
    plasma = parse_swpc_table(_get(URLS["swpc_plasma"].format(window=window)))
    mag = parse_swpc_table(_get(URLS["swpc_mag"].format(window=window)))
    df = plasma[["density", "speed"]].join(mag[["bz_gsm", "bt"]], how="outer")
    df = df.rename(columns={"density": "n_cc", "speed": "v_kms", "bz_gsm": "bz_nt", "bt": "bt_nt"})
    return derive(df.resample("1min").mean(), source="NOAA SWPC real-time (L1, not time-shifted)")


# ---------------------------------------------------------------------------
# NASA OMNI 1-min (bow-shock time-shifted) through the CDAWeb HAPI server
# ---------------------------------------------------------------------------
def parse_omni_hapi_csv(text: str, params=OMNI_PARAMS) -> pd.DataFrame:
    df = pd.read_csv(io.StringIO(text), header=None, names=["time"] + list(params), comment="#")
    df["time"] = pd.to_datetime(df["time"], utc=True)
    df = df.set_index("time")
    for p in params:
        df[p] = pd.to_numeric(df[p], errors="coerce")
        if p in OMNI_FILL:
            df.loc[df[p] >= OMNI_FILL[p], p] = np.nan
    return df


def fetch_omni(start, end, dataset: str = OMNI_DATASET) -> pd.DataFrame:
    url = URLS["omni_hapi"].format(dataset=dataset, params=",".join(OMNI_PARAMS),
                                   start=_iso(start), end=_iso(end))
    raw = parse_omni_hapi_csv(_get(url, timeout=120))
    df = raw.rename(columns={"flow_speed": "v_kms", "proton_density": "n_cc", "BZ_GSM": "bz_nt",
                             "Pressure": "pdyn_npa", "SYM_H": "symh_nt"})
    return derive(df, source=f"NASA OMNI ({dataset}, bow-shock nose)")


# ---------------------------------------------------------------------------
# Kp index (GFZ Potsdam, definitive + nowcast)
# ---------------------------------------------------------------------------
def parse_gfz_kp(text: str) -> pd.DataFrame:
    d = json.loads(text)
    df = pd.DataFrame({"time": pd.to_datetime(d["datetime"], utc=True), "kp": d["Kp"]})
    if "status" in d:
        df["status"] = d["status"]
    return df.set_index("time")


def fetch_kp(start, end) -> pd.DataFrame:
    return parse_gfz_kp(_get(URLS["gfz_kp"].format(start=_iso(start), end=_iso(end))))


def storm_level(kp: float) -> str:
    """NOAA G-scale from Kp."""
    if not np.isfinite(kp) or kp < 5:
        return "quiet / unsettled" if np.isfinite(kp) else "n/a"
    return {5: "G1 minor", 6: "G2 moderate", 7: "G3 strong", 8: "G4 severe"}.get(int(kp), "G5 extreme")


# ---------------------------------------------------------------------------
def derive(df: pd.DataFrame, source: str = "") -> pd.DataFrame:
    """Add dynamic pressure (if missing) and Shue magnetopause standoff."""
    df = df.copy()
    if "pdyn_npa" not in df or df["pdyn_npa"].isna().all():
        df["pdyn_npa"] = dynamic_pressure_npa(df["n_cc"], df["v_kms"])
    df["r0_mp_re"] = standoff_re(df["bz_nt"], df["pdyn_npa"])
    df.attrs["source"] = source
    return df


def latest_valid(df: pd.DataFrame, col: str):
    s = df[col].dropna()
    return (s.index[-1], float(s.iloc[-1])) if len(s) else (None, float("nan"))
