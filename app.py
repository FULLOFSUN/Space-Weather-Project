"""SMILE Companion dashboard.

    streamlit run app.py

Modes
-----
Live        : SMILE's latest TLE (CelesTrak) + real-time solar wind (NOAA SWPC) + Kp (GFZ)
Case study  : any past interval with NASA OMNI 1-min data + Kp (e.g. the April 2023 G4 storm)
Offline demo: nominal orbit + SYNTHETIC solar wind, for development without internet
"""
from __future__ import annotations

from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from smile_companion import magnetopause as mp
from smile_companion import orbit, spaceweather
from smile_companion.config import SMILE_NORAD_ID
from smile_companion.demo import synthetic_storm

st.set_page_config(page_title="SMILE Companion", page_icon="🛰️", layout="wide")

BLUE, RED, GREY = "#1f6fb4", "#d1495b", "#8a8a8a"


# --------------------------------------------------------------------------- data
@st.cache_data(ttl=3600, show_spinner="Fetching SMILE TLE from CelesTrak…")
def get_tle(norad: int):
    return orbit.parse_tle(orbit.fetch_tle(norad))


@st.cache_data(ttl=300, show_spinner="Fetching real-time solar wind (NOAA SWPC)…")
def get_swpc(window: str):
    return spaceweather.fetch_swpc(window)


@st.cache_data(ttl=24 * 3600, show_spinner="Fetching OMNI 1-min data (CDAWeb)…")
def get_omni(start, end):
    return spaceweather.fetch_omni(start, end)


@st.cache_data(ttl=1800, show_spinner="Fetching Kp (GFZ Potsdam)…")
def get_kp(start, end):
    return spaceweather.fetch_kp(start, end)


# --------------------------------------------------------------------------- sidebar
st.sidebar.title("SMILE Companion")
mode = st.sidebar.radio("Data mode", ["Live", "Case study", "Offline demo"])
norad = st.sidebar.number_input("NORAD ID", value=SMILE_NORAD_ID, step=1)

notes = []
tle = None
try:
    if mode == "Offline demo":
        sw = synthetic_storm()
        kp = sw.attrs["kp"]
        now = sw["r0_mp_re"].idxmin()
        notes.append("⚠️ **Synthetic demo data** — solar wind and Kp are made up; the orbit is "
                     "SMILE's nominal science orbit with an arbitrary phase.")
    elif mode == "Live":
        window = st.sidebar.selectbox("Solar-wind window", ["6-hour", "1-day", "3-day", "7-day"], index=1)
        sw = get_swpc(window)
        now = pd.Timestamp(datetime.now(timezone.utc))
        kp = get_kp(now - timedelta(days=3), now)
    else:
        c1, c2 = st.sidebar.columns(2)
        d0 = c1.date_input("Start", value=datetime(2023, 4, 23))
        d1 = c2.date_input("End", value=datetime(2023, 4, 25))
        sw = get_omni(pd.Timestamp(d0, tz="UTC"), pd.Timestamp(d1, tz="UTC"))
        kp = get_kp(pd.Timestamp(d0, tz="UTC"), pd.Timestamp(d1, tz="UTC"))
        naive = lambda t: t.tz_convert(None).to_pydatetime()  # noqa: E731  slider wants naive UTC
        sel = st.sidebar.slider("Time in interval", min_value=naive(sw.index[0]),
                                max_value=naive(sw.index[-1]),
                                value=naive(sw["r0_mp_re"].idxmin()), format="MM-DD HH:mm")
        now = pd.Timestamp(sel, tz="UTC")
        if sw.index[0] < pd.Timestamp("2026-05-19", tz="UTC"):
            notes.append("ℹ️ This interval is **before SMILE's launch** (19 May 2026): the orbit "
                         "shown is today's orbit replayed, to show what SMILE *would* have seen.")

    if mode != "Offline demo":
        _, l1, l2 = get_tle(int(norad))
        tle = (l1, l2)
except Exception as exc:  # network problems, bad dates, API changes
    st.error(f"Could not load data for **{mode}** mode: {exc}\n\nTry *Offline demo* mode.")
    st.stop()

el = orbit.tle_elements(*tle) if tle else orbit.nominal_elements(epoch=sw.index[0].to_pydatetime())
period = pd.Timedelta(hours=el.period_h)
when = now
times = pd.date_range(when - period / 2, when + period / 2, periods=500)
trk = orbit.track(times, tle=tle, elements=None if tle else el)
here = orbit.track(pd.DatetimeIndex([when]), tle=tle, elements=None if tle else el).iloc[0]

# upstream conditions at `now` (nearest valid sample within 30 min)
near = sw.loc[(sw.index > now - pd.Timedelta("30min")) & (sw.index <= now + pd.Timedelta("1min"))]
cur = near.dropna(subset=["bz_nt", "pdyn_npa"]).iloc[-1] if len(near.dropna(subset=["bz_nt", "pdyn_npa"])) \
    else sw.dropna(subset=["bz_nt", "pdyn_npa"]).iloc[-1]
bz, dp = float(cur["bz_nt"]), float(cur["pdyn_npa"])
r0 = float(mp.standoff_re(bz, dp))
dist_mp = float(mp.position_relative_to_magnetopause(here.x_gse_re, here.y_gse_re, here.z_gse_re, bz, dp))
kp_now = float(kp["kp"].loc[:now].iloc[-1]) if len(kp.loc[:now]) else float("nan")

# --------------------------------------------------------------------------- header
st.title("🛰️ SMILE Companion")
st.caption("Where ESA–CAS SMILE is, and what the solar wind is doing to the magnetosphere it images. "
           f"Solar wind: {sw.attrs.get('source', '')}. Orbit: "
           f"{'CelesTrak TLE, SGP4' if tle and orbit.HAVE_SGP4 else 'CelesTrak TLE, two-body' if tle else 'nominal elements'}.")
for n_ in notes:
    st.info(n_)

k = st.columns(6)
k[0].metric("SMILE altitude", f"{here.alt_km:,.0f} km", f"{here.r_re:.1f} R_E from centre", delta_color="off")
if np.isfinite(dist_mp):
    k[1].metric("vs magnetopause", "outside" if dist_mp > 0 else "inside", f"{dist_mp:+.1f} R_E", delta_color="off")
else:
    k[1].metric("vs magnetopause", "nightside", "tail region, no distance", delta_color="off")
k[2].metric("Solar-wind speed", f"{cur['v_kms']:.0f} km/s")
k[3].metric("IMF Bz (GSM)", f"{bz:+.1f} nT", "southward" if bz < 0 else "northward", delta_color="off")
k[4].metric("Standoff r0 (Shue 98)", f"{r0:.1f} R_E", f"Pdyn {dp:.1f} nPa", delta_color="off")
k[5].metric("Kp", f"{kp_now:.1f}" if np.isfinite(kp_now) else "–", spaceweather.storm_level(kp_now), delta_color="off")
st.caption(f"Conditions at {now:%Y-%m-%d %H:%M} UTC")

# --------------------------------------------------------------------------- orbit views
left, right = st.columns([1.05, 1])
with left:
    st.subheader("Orbit vs magnetopause (GSE X–Z)")
    fig = go.Figure()
    for (b, p, col, name) in ((2.0, 1.5, GREY, "quiet (Bz +2 nT, 1.5 nPa)"),
                              (bz, dp, RED, f"now (Bz {bz:+.1f} nT, {dp:.1f} nPa)")):
        x, rho = mp.surface(b, p)
        fig.add_trace(go.Scatter(x=np.r_[x, np.nan, x], y=np.r_[rho, np.nan, -rho], mode="lines",
                                 line=dict(color=col, width=2), name=f"magnetopause, {name}"))
    fig.add_trace(go.Scatter(x=trk.x_gse_re, y=trk.z_gse_re, mode="lines", line=dict(color=BLUE),
                             name="SMILE ±½ orbit", customdata=trk.time.dt.strftime("%m-%d %H:%M"),
                             hovertemplate="%{customdata}<br>X %{x:.1f} R_E, Z %{y:.1f} R_E"))
    fig.add_trace(go.Scatter(x=[here.x_gse_re], y=[here.z_gse_re], mode="markers",
                             marker=dict(size=11, color="black"), name="SMILE now"))
    fig.add_shape(type="circle", x0=-1, y0=-1, x1=1, y1=1, fillcolor="#3a7bd5", line_width=0)
    fig.update_layout(height=520, xaxis=dict(range=[-25, 25], title="X GSE [R_E] (Sun →)"),
                      yaxis=dict(range=[-25, 25], title="Z GSE [R_E]", scaleanchor="x"),
                      legend=dict(orientation="h", y=-0.2), margin=dict(t=10))
    st.plotly_chart(fig, use_container_width=True)

with right:
    st.subheader("Ground track")
    g = go.Figure(go.Scattergeo(lat=trk.lat_deg, lon=trk.lon_deg, mode="lines", line=dict(color=BLUE),
                                name="±½ orbit"))
    g.add_trace(go.Scattergeo(lat=[here.lat_deg], lon=[here.lon_deg], mode="markers",
                              marker=dict(size=11, color="black"), name="now"))
    g.update_layout(height=300, margin=dict(t=10, b=0), geo=dict(projection_type="natural earth",
                                                                 showland=True, landcolor="#e8e8e8"))
    st.plotly_chart(g, use_container_width=True)
    st.subheader("Altitude")
    a = go.Figure(go.Scatter(x=trk.time, y=trk.alt_km / 1e3, line=dict(color=BLUE)))
    a.add_vline(x=when, line_dash="dot")
    a.update_layout(height=200, margin=dict(t=10), yaxis_title="Altitude [10³ km]")
    st.plotly_chart(a, use_container_width=True)

# --------------------------------------------------------------------------- solar wind
st.subheader("Upstream driver and magnetospheric response")
rows = [("v_kms", "Speed [km/s]"), ("n_cc", "Density [cm⁻³]"), ("bz_nt", "Bz GSM [nT]"),
        ("pdyn_npa", "Pdyn [nPa]"), ("r0_mp_re", "Standoff r0 [R_E]")]
has_symh = "symh_nt" in sw and sw["symh_nt"].notna().any()
n_rows = len(rows) + 1 + int(has_symh)
fig = make_subplots(rows=n_rows, cols=1, shared_xaxes=True, vertical_spacing=0.015)
for i, (col, lab) in enumerate(rows, start=1):
    fig.add_trace(go.Scattergl(x=sw.index, y=sw[col], line=dict(color=BLUE, width=1), name=lab,
                               showlegend=False), row=i, col=1)
    fig.update_yaxes(title_text=lab, row=i, col=1, title_font_size=11)
fig.add_hline(y=0, line_color="black", line_width=0.6, row=3, col=1)
fig.add_hline(y=6.6, line_color=RED, line_dash="dash", row=5, col=1,
              annotation_text="geosynchronous orbit", annotation_position="bottom right")
r = len(rows) + 1
if has_symh:
    fig.add_trace(go.Scattergl(x=sw.index, y=sw["symh_nt"], line=dict(color=RED, width=1),
                               showlegend=False), row=r, col=1)
    fig.update_yaxes(title_text="SYM-H [nT]", row=r, col=1, title_font_size=11)
    r += 1
fig.add_trace(go.Bar(x=kp.index + pd.Timedelta("90min"), y=kp["kp"], width=3 * 3600e3 * 0.9,
                     marker_color=np.where(kp["kp"] >= 5, RED, GREY), showlegend=False), row=r, col=1)
fig.update_yaxes(title_text="Kp", range=[0, 9], row=r, col=1, title_font_size=11)
fig.add_vline(x=now, line_dash="dot")
fig.update_layout(height=170 * n_rows, margin=dict(t=10))
st.plotly_chart(fig, use_container_width=True)

# --------------------------------------------------------------------------- LEO impact
st.subheader("LEO impact (bring your own satellite)")
st.write("Upload a CSV with a `time` column and one numeric column (e.g. semi-major-axis residuals "
         "or drag acceleration from `leo-storm-decay`) to overlay it on Kp.")
up = st.file_uploader("CSV", type="csv")
if up is not None:
    df = pd.read_csv(up, parse_dates=["time"]).set_index("time")
    col = st.selectbox("Column", [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])])
    f2 = make_subplots(specs=[[{"secondary_y": True}]])
    f2.add_trace(go.Bar(x=kp.index, y=kp["kp"], name="Kp", marker_color=GREY, opacity=0.5), secondary_y=True)
    f2.add_trace(go.Scatter(x=df.index, y=df[col], name=col, line=dict(color=BLUE)), secondary_y=False)
    f2.update_yaxes(title_text="Kp", range=[0, 9], secondary_y=True)
    st.plotly_chart(f2, use_container_width=True)

st.caption("Models: Shue et al. (1998) magnetopause; SGP4/SDP4 orbit (Vallado et al. 2006). "
           "Data: CelesTrak, NOAA SWPC, NASA OMNI via CDAWeb HAPI, GFZ Potsdam Kp.")
