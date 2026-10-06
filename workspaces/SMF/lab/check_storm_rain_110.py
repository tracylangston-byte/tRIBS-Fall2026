"""
check_storm_rain_110.py   (version 1, read-only: it writes nothing)

PURPOSE
  For 1 Aug to 29 Sep 2014, put three things side by side for each day that matters:
    - how much rain each of the two rain gauges (SMF, SMPHQ) recorded that day,
    - the highest real streamflow at the SMF gauge that day,
    - which day is one of the four storms we plan to use (12 Aug, 19 Aug, 8 Sep, 27 Sep).
  It answers: when the model makes no flow on a storm date, is that the model (for example channel loss
  absorbing the storm) or is there simply little rain in the forcing files for that day?

RUN (from the lab/ folder)
  python check_storm_rain_110.py

NOTES
  - Rain files are read as rates in mm/hr, one row per 15 minutes (4 rows share each hour), so a day's total is
    sum(R) x 0.25 h. If that assumption is wrong the totals are off by a constant factor; the peak-intensity column
    and the inch totals let you check it against the storm list (8 Sep: about 3.56 in at SMF).
  - Only days with 5 mm or more of rain at either gauge, or 1 cfs or more of real flow, or one of the four storms
    are listed.
"""
import sys
from pathlib import Path

import pandas as pd

MET = Path("../init_data/met")
RAIN = {"SMF": MET / "precip_SMF_2014-2014.mdf", "SMPHQ": MET / "precip_SMPHQ_2014-2014.mdf"}
OBS = MET / "SMF_Observations_1993-2025.xlsx"
START, END = pd.Timestamp("2014-08-01"), pd.Timestamp("2014-09-29")
STORMS = {"2014-08-12": "storm 1 (12 Aug)", "2014-08-19": "storm 2 (19 Aug)",
          "2014-09-08": "storm 3 (8 Sep)", "2014-09-27": "storm 4 (27 Sep)"}
ROW_HR = 0.25
MM_PER_IN = 25.4


def day_rain(path):
    if not path.exists():
        sys.exit("Cannot find %s. Run this from the lab/ folder." % path)
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    rcol = [c for c in df.columns if c.upper().startswith("R")][0]
    d = pd.to_datetime(dict(year=df["Year"], month=df["Month"], day=df["Day"]))
    out = pd.DataFrame({"date": d, "R": df[rcol]})
    g = out.groupby("date")["R"]
    return pd.DataFrame({"mm": g.sum() * ROW_HR, "peak": g.max()})


rain = {k: day_rain(p) for k, p in RAIN.items()}

if not OBS.exists():
    sys.exit("Cannot find %s. Run this from the lab/ folder." % OBS)
obs = pd.read_excel(OBS, sheet_name="Discharge", skiprows=6)
obs["date"] = pd.to_datetime(obs["Date"].astype(str)).dt.normalize()
q = obs.groupby("date")["cfs"].agg(["max", "count"])

days = pd.date_range(START, END, freq="D")
t = pd.DataFrame(index=days)
t["SMF mm"] = rain["SMF"]["mm"].reindex(days).fillna(0.0)
t["SMPHQ mm"] = rain["SMPHQ"]["mm"].reindex(days).fillna(0.0)
t["SMF peak mm/hr"] = rain["SMF"]["peak"].reindex(days).fillna(0.0)
t["real max cfs"] = q["max"].reindex(days)
t["gauge readings"] = q["count"].reindex(days).fillna(0).astype(int)
t["note"] = [STORMS.get(d.strftime("%Y-%m-%d"), "") for d in days]

keep = (t["SMF mm"] >= 5) | (t["SMPHQ mm"] >= 5) | (t["real max cfs"].fillna(0) >= 1) | (t["note"] != "")
show = t[keep].copy()
show.index = show.index.strftime("%d %b")
show["real max cfs"] = show["real max cfs"].map(lambda v: "no data" if pd.isna(v) else "%.1f" % v)
print("Forcing and real gauge, 1 Aug to 29 Sep 2014 (rain in mm per calendar day, flow in cfs)")
print(show.to_string(float_format=lambda v: "%.1f" % v))
print()
print("Rain on the four storm dates in inches (SMF / SMPHQ), for matching against the storm list:")
for d, lab in STORMS.items():
    r = t.loc[pd.Timestamp(d)]
    print("  %-18s SMF %.2f in   SMPHQ %.2f in   real max %s cfs"
          % (lab, r["SMF mm"] / MM_PER_IN, r["SMPHQ mm"] / MM_PER_IN,
             "no data" if pd.isna(r["real max cfs"]) else "%.1f" % r["real max cfs"]))
print()
print("Whole window: SMF %.0f mm, SMPHQ %.0f mm; days with real flow of 1 cfs or more: %d"
      % (t["SMF mm"].sum(), t["SMPHQ mm"].sum(), int((t["real max cfs"].fillna(0) >= 1).sum())))
