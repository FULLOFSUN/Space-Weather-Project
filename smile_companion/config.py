"""Mission constants and data-source URLs."""

SMILE_NORAD_ID = 69123            # COSPAR 2026-109A, launched 19 May 2026 on Vega-C

# Nominal science orbit (ESA SMILE operations page): ~5,000 km perigee,
# ~121,182 km (19 RE) apogee, 70 deg inclination for a Vega-C launch, perigee
# over the southern hemisphere. Used ONLY in offline/demo mode; live mode
# uses the latest TLE from CelesTrak.
NOMINAL_ORBIT = {
    "perigee_alt_km": 5000.0,
    "apogee_alt_km": 121182.0,
    "inclination_deg": 70.0,
    "raan_deg": 0.0,              # arbitrary in demo mode
    "argp_deg": 270.0,            # perigee at the southernmost point
    "mean_anomaly_deg": 180.0,    # start at apogee
}

RE_KM = 6371.0                    # mean Earth radius used for R_E units
MU_KM3_S2 = 398600.4418

URLS = {
    "tle": "https://celestrak.org/NORAD/elements/gp.php?CATNR={norad}&FORMAT=TLE",
    "swpc_plasma": "https://services.swpc.noaa.gov/products/solar-wind/plasma-{window}.json",
    "swpc_mag": "https://services.swpc.noaa.gov/products/solar-wind/mag-{window}.json",
    "omni_hapi": ("https://cdaweb.gsfc.nasa.gov/hapi/data?id={dataset}"
                  "&parameters={params}&time.min={start}&time.max={end}&format=csv"),
    "gfz_kp": "https://kp.gfz-potsdam.de/app/json/?start={start}&end={end}&index=Kp",
}

OMNI_DATASET = "OMNI_HRO_1MIN"
OMNI_PARAMS = ["flow_speed", "proton_density", "BZ_GSM", "Pressure", "SYM_H"]
# OMNI fill values (anything at or above these is missing)
OMNI_FILL = {"flow_speed": 99999.0, "proton_density": 999.0, "BZ_GSM": 9999.0,
             "Pressure": 99.0, "SYM_H": 99999.0}

USER_AGENT = "smile-companion/0.1 (research dashboard)"
