#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wy2017_slide_data.py   (v1.0)
=============================
READ-ONLY. Pulls the small tables needed for the "Original vs Re-tuned, whole water year 2017" slide out of result files
that already exist. It runs no model, changes no input file and overwrites nothing.

WHAT IT READS (from workspaces/Sycamore/lab/, or change the folders with --results-dir / --obs-file)
  results/W5_d012_k10_base_nocl_Outlet.qout                     Original, outlet flow
  results/W5_d012_k10_f0p15_as300_ksx1p03_nocl_Outlet.qout      Re-tuned, outlet flow
  results/<run>8759_00.rft and <run>8759_00.mrf                 runoff made on the land, by type, and the basin rain (both runs)
  ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv        USGS gauge 09510200

WHAT IT WRITES (a new folder, wy2017_slide_data/ by default)
  wy2017_daily_flow.csv      one row per day, Oct 1 2016 - Sep 30 2017: daily mean and daily highest hourly flow (m3/s) for the gauge,
                             Original and Re-tuned, how many gauge hours each day has, and the basin rain (mm per day)
  wy2017_runoff_by_type.csv  one row per run and month, plus a "year" row: Hsrf, Sbsrf, Psrf, Satsrf and their total (million m3)
  wy2017_checks.txt          the checks printed below, saved so you can attach them

RUN IT from workspaces/Sycamore/lab/:
    python wy2017_slide_data.py

WORDS
  Hsrf   infiltration excess     Sbsrf  saturation excess     Psrf  perched return flow     Satsrf  groundwater exfiltration
  runoff made   runoff GENERATED on the land, all cells added, BEFORE routing (not flow at the outlet)
  hourly mean   the model and gauge files are averaged into hours first, as in the project's other scripts
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

START = pd.Timestamp("2016-10-01 00:00")
END = pd.Timestamp("2017-10-01 00:00")
RUNS = [("Original", "W5_d012_k10_base_nocl"),
        ("Re-tuned", "W5_d012_k10_f0p15_as300_ksx1p03_nocl")]
TYPES = ["Hsrf", "Sbsrf", "Psrf", "Satsrf"]
# Totals already on record (W5_Run_Report_v1, W5_Kchan_Report_v1, Model_State_Output_v1), for the checks at the end.
KNOWN = {"gauge year": 18.07, "Original delivered": 38.62, "Re-tuned delivered": 38.38,
         "Original made": 38.76, "Re-tuned made": 38.52}

lines = []


def log(msg=""):
    print(msg, flush=True)
    lines.append(msg)


def find(folder, name):
    """Look for a file in the given folder first, then anywhere below the current folder."""
    p = Path(folder) / name
    if p.is_file():
        return p
    hits = sorted(Path(".").rglob(name))
    if hits:
        log(f"  (did not find {p}; using {hits[0]})")
        return hits[0]
    sys.exit(f"STOPPED: cannot find {name}\nRun this from workspaces/Sycamore/lab/, or use --results-dir to point at the folder that holds it.")


def read_qout(path):
    raw = pd.read_csv(path).iloc[:, :2]
    raw.columns = ["t", "Q"]
    s = pd.Series(raw.Q.to_numpy(float), index=START + pd.to_timedelta(raw.t.to_numpy(), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_hourly_table(path, need):
    d = pd.read_csv(path)
    missing = [c for c in need if c not in d.columns]
    if missing:
        sys.exit(f"STOPPED: {path} has no column(s) {missing}. Columns found: {list(d.columns)}")
    d = d[d["Time_hr"] <= 8759.5].copy()                       # the last row (hour 8760) is a filler row
    d.index = START + pd.to_timedelta(d["Time_hr"].to_numpy(float), unit="h")   # label = END of the hour
    return d


def hourly(q, a, b):
    s = q[(q.index >= pd.Timestamp(a)) & (q.index < pd.Timestamp(b))]
    return s.resample("60min").mean().dropna()


def vol(h):
    return float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--obs-file", default="../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv")
    ap.add_argument("--out-dir", default="wy2017_slide_data")
    a = ap.parse_args()
    out = Path(a.out_dir)
    if out.exists() and any(out.iterdir()):
        sys.exit(f"STOPPED: {out}/ already has files in it; this script does not overwrite. Use --out-dir with a new name.")
    out.mkdir(parents=True, exist_ok=True)

    log("wy2017_slide_data.py v1.0 (read-only)")
    gpath = find(Path(a.obs_file).parent, Path(a.obs_file).name)
    g = pd.read_csv(gpath, parse_dates=["datetime"]).set_index("datetime")["q_cms"]
    q, mrf, rft = {}, {}, {}
    for lab, rid in RUNS:
        q[lab] = read_qout(find(a.results_dir, f"{rid}_Outlet.qout"))
        mrf[lab] = read_hourly_table(find(a.results_dir, f"{rid}8759_00.mrf"), ["Time_hr", "MAP_mm_hr", "Srf_m3_s"])
        rft[lab] = read_hourly_table(find(a.results_dir, f"{rid}8759_00.rft"), ["Time_hr"] + [f"{t}_m3_s" for t in TYPES])
    log(f"  gauge file: {gpath}")

    # ---- 1. daily flow ---------------------------------------------------------------------------------
    days = pd.date_range(START, END - pd.Timedelta(days=1), freq="D")
    rows = []
    hg, h1, h2 = hourly(g, START, END), hourly(q["Original"], START, END), hourly(q["Re-tuned"], START, END)
    rain_h = mrf["Original"]["MAP_mm_hr"]
    for d0 in days:
        d1 = d0 + pd.Timedelta(days=1)
        sel = lambda h: h[(h.index >= d0) & (h.index < d1)]
        gd, od, rd = sel(hg), sel(h1), sel(h2)
        rn = rain_h[(rain_h.index > d0) & (rain_h.index <= d1)]          # hour-END labels: the hours of this day
        rows.append(dict(date=d0.strftime("%Y-%m-%d"),
                         gauge_hours=len(gd), gauge_mean_m3s=gd.mean() if len(gd) else np.nan, gauge_max_m3s=gd.max() if len(gd) else np.nan,
                         original_mean_m3s=od.mean(), original_max_m3s=od.max(),
                         retuned_mean_m3s=rd.mean(), retuned_max_m3s=rd.max(), rain_mm=rn.sum()))
    D = pd.DataFrame(rows)
    D.to_csv(out / "wy2017_daily_flow.csv", index=False, float_format="%.5f")

    # ---- 2. runoff made by type, by month --------------------------------------------------------------
    months = pd.period_range("2016-10", "2017-09", freq="M")
    trows = []
    for lab, _ in RUNS:
        d = rft[lab]
        hs = d.copy()
        hs.index = hs.index - pd.Timedelta(hours=1)                     # label by the START of each hour (calendar months)
        per = hs.index.to_period("M")
        sums = {t: hs[f"{t}_m3_s"].groupby(per).sum() * 3600 / 1e6 for t in TYPES}
        for p in months:
            v = {t: float(sums[t].get(p, 0.0)) for t in TYPES}
            trows.append(dict(run=lab, month=str(p), **{t: v[t] for t in TYPES}, total=sum(v.values())))
        yv = {t: float(d[f"{t}_m3_s"].sum() * 3600 / 1e6) for t in TYPES}
        trows.append(dict(run=lab, month="year", **yv, total=sum(yv.values())))
    T = pd.DataFrame(trows)
    T.to_csv(out / "wy2017_runoff_by_type.csv", index=False, float_format="%.5f")

    # ---- 3. checks -------------------------------------------------------------------------------------
    log("\nCHECKS (the 'on record' numbers come from the project's earlier reports; small differences in the last digit are fine)")
    got = {"gauge year": vol(hg), "Original delivered": vol(h1), "Re-tuned delivered": vol(h2)}
    for lab, _ in RUNS:
        got[f"{lab} made"] = float(mrf[lab]["Srf_m3_s"].sum() * 3600 / 1e6)
    for k, want in KNOWN.items():
        have = got[k]
        log(f"  {k:20s} here {have:7.2f} M m3   on record {want:6.2f}   {'ok' if abs(have - want) <= 0.02 else 'DIFFERS'}")
    for lab, _ in RUNS:
        tot = float(T[(T.run == lab) & (T.month == "year")]["total"].iloc[0])
        srf = got[f"{lab} made"]
        log(f"  {lab:20s} the four types add to {tot:7.3f} M m3; the .mrf total runoff (Srf) is {srf:7.3f}   {'ok' if abs(tot - srf) <= 0.02 else 'DIFFERS'}")
    log("\nMONTHLY OUTLET VOLUME, M m3 (sums of hourly mean flows)")
    log(f"  {'month':9s} {'gauge':>8s} {'Original':>9s} {'Re-tuned':>9s}")
    for p in months:
        a0, b0 = p.start_time, p.end_time.ceil("D")
        log(f"  {str(p):9s} {vol(hourly(g, a0, b0)):8.3f} {vol(hourly(q['Original'], a0, b0)):9.3f} {vol(hourly(q['Re-tuned'], a0, b0)):9.3f}")
    log("\nRUNOFF MADE ON THE LAND, WHOLE YEAR, M m3 (before routing)")
    log(f"  {'run':9s} {'Hsrf':>8s} {'Sbsrf':>8s} {'Psrf':>8s} {'Satsrf':>8s} {'total':>8s}")
    for lab, _ in RUNS:
        r = T[(T.run == lab) & (T.month == "year")].iloc[0]
        log(f"  {lab:9s} {r.Hsrf:8.3f} {r.Sbsrf:8.3f} {r.Psrf:8.3f} {r.Satsrf:8.3f} {r.total:8.3f}")
    (out / "wy2017_checks.txt").write_text("\n".join(lines) + "\n")
    log(f"\nWrote {out}/wy2017_daily_flow.csv, {out}/wy2017_runoff_by_type.csv and {out}/wy2017_checks.txt")
    log("Next step: attach those three files to the chat.")


if __name__ == "__main__":
    main()
