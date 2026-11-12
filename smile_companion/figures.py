"""Static matplotlib quick-looks (used by scripts/make_figures.py and the README)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import magnetopause as mp

C_ORBIT, C_QUIET, C_STORM, C_NOW = "#1f6fb4", "#7f7f7f", "#d1495b", "#111111"


def orbit_vs_magnetopause(trk: pd.DataFrame, now: pd.Timestamp, quiet: tuple[float, float],
                          storm: tuple[float, float], title_note: str = ""):
    """SMILE track in the GSE X-Z plane with quiet and storm magnetopause curves.
    quiet/storm = (Bz [nT], Pdyn [nPa])."""
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.add_patch(plt.Circle((0, 0), 1, color="#3a7bd5", zorder=3))
    for (bz, dp), col, lab in ((quiet, C_QUIET, "quiet"), (storm, C_STORM, "storm")):
        x, rho = mp.surface(bz, dp)
        r0 = float(mp.standoff_re(bz, dp))
        ax.plot(x, rho, color=col, lw=1.8, label=f"magnetopause, {lab}: Bz={bz:+.0f} nT, "
                                                    f"Pdyn={dp:.1f} nPa, r0={r0:.1f} RE")
        ax.plot(x, -rho, color=col, lw=1.8)
    ax.plot(trk["x_gse_re"], trk["z_gse_re"], color=C_ORBIT, lw=1.2, label="SMILE orbit (one period)")
    i = int(np.argmin(np.abs((trk["time"] - now).dt.total_seconds())))
    ax.plot(trk["x_gse_re"].iloc[i], trk["z_gse_re"].iloc[i], "o", color=C_NOW, ms=7, zorder=5,
            label=f"SMILE at {now:%Y-%m-%d %H:%M} UTC")
    ax.annotate("Sun →", xy=(0.97, 0.5), xycoords="axes fraction", ha="right", color="#b07d00")
    ax.set(xlim=(-25, 25), ylim=(-25, 25), aspect="equal", xlabel="X GSE [R_E]", ylabel="Z GSE [R_E]",
           title=f"SMILE and the magnetopause (GSE X–Z){title_note}")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left", fontsize=7.5, frameon=True)
    fig.tight_layout()
    return fig


def solar_wind_panel(sw: pd.DataFrame, kp: pd.DataFrame | None = None, note: str = ""):
    rows = [("v_kms", "Speed\n[km/s]"), ("n_cc", "Density\n[cm⁻³]"), ("bz_nt", "Bz GSM\n[nT]"),
            ("pdyn_npa", "Pdyn\n[nPa]"), ("r0_mp_re", "Standoff r0\n[R_E]")]
    extra = kp is not None and len(kp)
    fig, axes = plt.subplots(len(rows) + bool(extra), 1, figsize=(9, 9), sharex=True)
    for ax, (col, lab) in zip(axes, rows):
        ax.plot(sw.index, sw[col], lw=0.8, color=C_ORBIT)
        ax.set_ylabel(lab, fontsize=9)
        ax.grid(alpha=0.3)
    axes[2].axhline(0, color="k", lw=0.6)
    axes[4].axhline(6.6, color=C_STORM, lw=0.8, ls="--")
    axes[4].text(sw.index[0], 6.7, " geosynchronous orbit (6.6 R_E)", color=C_STORM, fontsize=8, va="bottom")
    if extra:
        ax = axes[-1]
        colors = np.where(kp["kp"] >= 5, C_STORM, C_QUIET)
        ax.bar(kp.index, kp["kp"], width=0.11, align="edge", color=colors)
        ax.set(ylabel="Kp", ylim=(0, 9))
        ax.axhline(5, color=C_STORM, lw=0.6, ls="--")
        ax.grid(alpha=0.3)
    axes[0].set_title(f"Upstream solar wind and magnetopause response{note}")
    fig.autofmt_xdate()
    fig.tight_layout()
    return fig
