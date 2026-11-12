"""Synthetic inputs for offline development.

Everything here is MADE-UP data shaped like a moderate CME-driven storm.
It exists so the dashboard and figures run without internet access; it is
never mixed with real data and the app labels it as synthetic.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .spaceweather import derive


def synthetic_storm(start: str = "2026-11-01", days: float = 3.0, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    t = pd.date_range(start, periods=int(days * 1440), freq="1min", tz="UTC")
    h = (t - t[0]).total_seconds().to_numpy() / 3600.0
    shock, ejecta_end = 20.0, 44.0
    after = h >= shock
    in_cme = after & (h < ejecta_end)
    v = np.where(after, 650 - 120 * np.clip((h - shock) / 40, 0, 1), 380.0)
    n = np.where(after, 4 + 18 * np.exp(-(h - shock) / 3.0), 5.0)
    bz = np.where(in_cme, -12 * np.sin(np.pi * (h - shock) / (ejecta_end - shock)), 1.0)
    df = pd.DataFrame({
        "v_kms": v + rng.normal(0, 8, h.size),
        "n_cc": np.clip(n + rng.normal(0, 0.6, h.size), 0.5, None),
        "bz_nt": bz + rng.normal(0, 1.2, h.size),
    }, index=t)
    df.index.name = "time"
    out = derive(df, source="SYNTHETIC DEMO DATA (not observations)")
    # a crude Kp proxy so the panels have something to show
    kp3 = out["bz_nt"].clip(upper=0).abs().resample("3h").mean() * 0.45 \
        + out["v_kms"].resample("3h").mean() / 250
    out.attrs["kp"] = pd.DataFrame({"kp": np.clip(np.round(kp3 * 3) / 3, 0, 9)})
    return out
