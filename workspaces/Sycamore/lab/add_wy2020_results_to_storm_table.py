#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
add_wy2020_results_to_storm_table.py   (v1.0, 2026-11-06)
=========================================================

Fills the WY2020 gauge and model result columns of storm_table.csv and writes the result as a NEW file.
Your storm_table.csv is never changed.

WHAT IT FILLS (WY2020 winter storms only: the 7 on the 20 mm list plus Dec 3-4, Dec 7-9 and Jan 21)
  Gauge hours, Window (d)                      hours of gauge data / hours in the window, and the window length
  Gauge volume, Gauge runoff ratio             outlet volume at the USGS gauge (M m3) and volume / (basin rain x 0.4246 M m3 per mm)
  Model volume, Model minus gauge, Model runoff ratio    the same for the Re-tuned run (all hours, not only gauge hours)
  Original, loss off / Kchan 0.5 / Kchan 0.2   Original-setting volumes (M m3) in the same window
  Gauge / Model volume, storm window           the scored storm window (see below): storms of the features run only
  Water table depth                            model basin-mean depth to the water table at 00:00 on the first rain day (Re-tuned)
  Root-zone moisture, Outlet flow at start     the same moment; storms of the features run only

  These are the same definitions used for WY2017 (model_state_analysis_wy2017.py):
    Event window   00:00 on the first rain day to 00:00 four days after the last rain day, or the next event's first day if that is sooner.
    Storm window   00:00 on the first rain day to 00:00 four days after the last rain day (does NOT stop at the next event).
    Volume         hourly mean flow x 3,600 s, added over the hours that have data (so a gauge volume with missing hours is a minimum).
    Features-run storms: Nov 19-22, Nov 28-30, Dec 7-9 (the near miss), Dec 23-28, Feb 10-11, Feb 22-23, Mar 10-13, Mar 18-19.

WHAT IT LEAVES ALONE
  Every other cell, including all rain columns, the pre-registered columns and the Note column (except as listed below).
  The two September storms (no gauge record before Oct 1). The temperature and model-snow columns (parked analysis, WY2017 only).
  Where the gauge has NO data in a window (Dec 3-4 and Dec 7-9 fall in the Dec 2-17 dropout), the gauge cells stay empty
  and a short sentence is added to that row's Note. Rows whose gauge has some hours missing get the same kind of note.

CHECKS (printed)
  1. The seven scored storm windows (gauge, Re-tuned, Original, Kchan 0.5, Kchan 0.2) are compared with the WY2020 write-up: MATCH or MISMATCH.
  2. The model's own basin rain in the .mrf file is compared with the table's Basin rain for each storm (shows the .mrf is the right file and time-aligned).
  3. After writing, every cell that was not meant to change is compared with the original: it must be identical.
  If a check says MISMATCH, do not use the new file until it is explained.

WHAT IT READS (read-only)
  storm_table.csv                                       short headers (year, first_day, ...) OR the readable ones (Water year, First rain day, ...)
  results/Sycamore_WY2020_Retuned_Outlet.qout, ..._Original_Outlet.qout, ..._Original_kchan0p5_Outlet.qout, ..._Original_kchan0p2_Outlet.qout
                                                        (Time_hr 0 = Sep 1 2019 00:00)
  ../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv                columns datetime, q_cms
  results/Sycamore_WY2020_Retuned*.mrf                  the Re-tuned run's basin-average file (for water table depth and moisture). If it is missing,
                                                        those three columns stay empty and the script says so.
WHAT IT WRITES
  storm_table_with_wy2020.csv   (next to storm_table.csv unless you set --out)

COMMANDS (from workspaces/Sycamore/lab/)
  python add_wy2020_results_to_storm_table.py --dry-run     # reads everything, lists the windows, writes nothing
  python add_wy2020_results_to_storm_table.py               # writes storm_table_with_wy2020.csv
  python add_wy2020_results_to_storm_table.py --mrf results/<file>.mrf     # if the .mrf file has another name

THEN (your step, by hand): look at the printed checks, open the new file, and only if you are happy copy it over storm_table.csv.
NOTE: build_storm_table_sycamore.py (the script that first makes storm_table.csv) still leaves these columns empty. If you rebuild the
table with it, run this script again afterward.

WORDS
  runoff ratio   outlet volume / rain volume on the basin (1 mm over the basin = 0.4246 M m3).   M m3   million cubic meters.
  .mrf           the basin-average output tRIBS writes each hour (rain, water table, soil moisture, ...).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUN_START = pd.Timestamp("2019-09-01 00:00")
MM_TO_MM3 = 424.6 / 1000.0            # 1 mm of rain over the basin = 0.4246 M m3
AFTER_DAYS = 3
CHECK_TOL = 0.002                     # M m3 (the write-up gives 3 decimals)
RAIN_TOL_MM = 0.6                     # .mrf rain vs table rain, mm over a whole storm (rounding and the 1-hour label shift)

RUN_STEMS = {"retuned": "Sycamore_WY2020_Retuned", "original": "Sycamore_WY2020_Original",
             "kchan05": "Sycamore_WY2020_Original_kchan0p5", "kchan02": "Sycamore_WY2020_Original_kchan0p2"}
RUN_LABEL = {"retuned": "Re-tuned", "original": "Original", "kchan05": "Kchan 0.5", "kchan02": "Kchan 0.2"}

# first days of the storms the features run covered (their gauge/model storm-window and start-state columns are filled)
FEATURES_FIRST_DAYS = {"2019-11-19", "2019-11-28", "2019-12-07", "2019-12-23", "2020-02-10", "2020-02-22", "2020-03-10", "2020-03-18"}

# WY2020 write-up, scored storm windows (M m3): gauge, Re-tuned, Original, Kchan 0.5, Kchan 0.2
EXPECTED = {
    "Nov 19-22": (0.015, 2.800, 3.165, 3.019, 3.091), "Nov 28-30": (1.473, 3.235, 3.538, 3.341, 3.433),
    "Dec 23-28": (4.069, 0.218, 0.154, 0.044, 0.083), "Feb 10-11": (0.386, 1.767, 1.765, 1.641, 1.700),
    "Feb 22-23": (4.131, 4.839, 4.981, 4.825, 4.899), "Mar 10-13": (5.626, 0.875, 0.493, 0.340, 0.413),
    "Mar 18-19": (3.317, 1.031, 0.841, 0.683, 0.775),
}
EXPECTED_ORDER = ["gauge", "retuned", "original", "kchan05", "kchan02"]

CONTEXT = {  # key -> accepted headers
    "year": ["year", "Water year"], "event": ["event", "Event"], "first_day": ["first_day", "First rain day"],
    "last_day": ["last_day", "Last rain day"], "status": ["status", "List status"], "basin_mm": ["basin_mm", "Basin rain (mm)"],
}
NOTE_HEADERS = ["note", "Note"]
TARGETS = {
    "gauge_hours": ["gauge_hours", "Gauge hours"],
    "window_days": ["window_days", "Window (d)"],
    "gauge_Mm3": ["gauge_Mm3", "Gauge volume (M m3)"],
    "gauge_ratio_pct": ["gauge_ratio_pct", "Gauge runoff ratio (%)"],
    "model_Mm3": ["model_Mm3", "Model volume (M m3)"],
    "model_minus_gauge_Mm3": ["model_minus_gauge_Mm3", "Model minus gauge (M m3)"],
    "model_ratio_pct": ["model_ratio_pct", "Model runoff ratio (%)"],
    "gauge_storm_Mm3": ["gauge_storm_Mm3", "Gauge volume, storm window (M m3)"],
    "model_storm_Mm3": ["model_storm_Mm3", "Model volume, storm window (M m3)"],
    "asdrawn_nocl_Mm3": ["asdrawn_nocl_Mm3", "Original, loss off (M m3)"],
    "asdrawn_k05_Mm3": ["asdrawn_k05_Mm3", "Original, Kchan 0.5 (M m3)"],
    "asdrawn_k02_Mm3": ["asdrawn_k02_Mm3", "Original, Kchan 0.2 (M m3)"],
    "water_table_depth_mm": ["water_table_depth_mm", "Water table depth (mm)"],
    "soil_moisture_root": ["soil_moisture_root", "Root-zone moisture (fraction)"],
    "outlet_flow_start_m3s": ["outlet_flow_start_m3s", "Outlet flow at start (m3/s)"],
}


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    sys.exit("\nSTOPPED: " + msg + "\n")


def need(path: Path, what: str) -> Path:
    if not Path(path).is_file():
        die("cannot find the {} file:\n    {}\n  Run this from workspaces/Sycamore/lab/, or point to it with the matching option (see the top of this script).".format(what, path))
    return Path(path)


# ------------------------------------------------------------------------------------------------ readers and numbers
def read_qout(path: Path) -> pd.Series:
    raw = pd.read_csv(path, float_precision="round_trip").iloc[:, :2]
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    s = pd.Series(raw["Q"].to_numpy(float), index=RUN_START + pd.to_timedelta(raw["Time_hr"].to_numpy(float), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_obs(path: Path) -> pd.Series:
    o = pd.read_csv(path, float_precision="round_trip")
    for c in ("datetime", "q_cms"):
        if c not in o.columns:
            die("{} has no column '{}'. Columns found: {}".format(path, c, list(o.columns)[:10]))
    o["datetime"] = pd.to_datetime(o["datetime"])
    s = pd.Series(pd.to_numeric(o["q_cms"], errors="coerce").to_numpy(float), index=o["datetime"])
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_mrf(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, float_precision="round_trip")
    if "Time_hr" not in d.columns:
        die("{} has no 'Time_hr' column. Columns found: {}".format(path, list(d.columns)[:12]))
    if "RainMin" in d.columns:
        d = d[pd.to_numeric(d["RainMin"], errors="coerce") < 9000].copy()       # tRIBS adds a filler row at the end (RainMin 9999.99)
    d.index = RUN_START + pd.to_timedelta(d["Time_hr"].to_numpy(float), unit="h")      # label = END of the hour
    return d


def hourly(q: pd.Series, a, b) -> pd.Series:
    s = q[(q.index >= a) & (q.index < b)]
    return s.resample("60min").mean().dropna()


def vol(h: pd.Series) -> float:
    return float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")


def fnum(v, nd: int = 4) -> str:
    """Number -> text in the style of the existing table (4 decimals at most, no trailing zeros); empty for missing."""
    if v is None or not np.isfinite(v):
        return ""
    return repr(round(float(v), nd))


def pick_header(raw: pd.DataFrame, options: list):
    hit = [c for c in options if c in raw.columns]
    return hit[0] if hit else None


# ------------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fill the WY2020 result columns of storm_table.csv into a new file (see the top of this script).")
    ap.add_argument("--root", default=str(HERE), help="the lab folder (default: the folder this script is in)")
    ap.add_argument("--table", default=None, help="storm_table.csv (default <root>/storm_table.csv)")
    ap.add_argument("--results", default=None, help="folder with the WY2020 outlet files (default <root>/results)")
    ap.add_argument("--obs", default=None, help="USGS WY2020 gauge file (default <root>/../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv)")
    ap.add_argument("--mrf", default=None, help="the Re-tuned run's .mrf file (default: results/Sycamore_WY2020_Retuned*.mrf)")
    ap.add_argument("--out", default=None, help="the new file (default <root>/storm_table_with_wy2020.csv)")
    ap.add_argument("--force", action="store_true", help="overwrite WY2020 result cells that already have values (default: stop)")
    ap.add_argument("--dry-run", action="store_true", help="read everything and list the windows; write nothing")
    a = ap.parse_args(argv)

    root = Path(a.root)
    table = Path(a.table) if a.table else root / "storm_table.csv"
    results = Path(a.results) if a.results else root / "results"
    obs_path = Path(a.obs) if a.obs else root.parent / "init_data" / "obs" / "sycamore_09510200_obs_flow_WY2020.csv"
    out_path = Path(a.out) if a.out else root / "storm_table_with_wy2020.csv"
    if out_path.resolve() == table.resolve():
        die("--out is the same file as the table. This script only writes a NEW file.")

    log("\nadd_wy2020_results_to_storm_table: reading files")
    raw = pd.read_csv(need(table, "storm table"), dtype=str, keep_default_na=False, encoding="utf-8-sig")
    original = raw.copy()
    ctx, tgt, missing = {}, {}, []
    for k, opts in CONTEXT.items():
        h = pick_header(raw, opts)
        (ctx.__setitem__(k, h) if h else missing.append((k, opts)))
    for k, opts in TARGETS.items():
        h = pick_header(raw, opts)
        (tgt.__setitem__(k, h) if h else missing.append((k, opts)))
    if missing:
        die("the table is missing column(s): {}\n  All columns found: {}".format(["/".join(o) for _, o in missing], list(raw.columns)))
    note_h = pick_header(raw, NOTE_HEADERS)
    log("  table: {} rows, {} columns; header style: {}".format(len(raw), raw.shape[1], "short" if ctx["year"] == "year" else "readable"))

    y20 = raw[ctx["year"]].str.contains("2020")
    spin = raw[ctx["status"]].str.strip().str.startswith("Spin")
    winter_idx = list(raw.index[y20 & ~spin])
    if not winter_idx:
        die("no WY2020 winter rows found in the table.")
    ev = []
    for i in winter_idx:
        first, last = pd.to_datetime(raw.at[i, ctx["first_day"]], errors="coerce"), pd.to_datetime(raw.at[i, ctx["last_day"]], errors="coerce")
        if pd.isna(first) or pd.isna(last):
            die("could not read the dates of '{}' ({} / {}).".format(raw.at[i, ctx["event"]], raw.at[i, ctx["first_day"]], raw.at[i, ctx["last_day"]]))
        bm = pd.to_numeric(raw.at[i, ctx["basin_mm"]], errors="coerce")
        ev.append({"idx": i, "name": raw.at[i, ctx["event"]].replace("–", "-").strip(), "first": first, "last": last,
                   "basin_mm": float(bm) if pd.notna(bm) else float("nan"),
                   "listed": raw.at[i, ctx["status"]].strip().startswith("On 20")})
    ev.sort(key=lambda d: d["first"])
    for k, e in enumerate(ev):
        e["lo"] = e["first"]
        e["hi_storm"] = e["last"] + pd.Timedelta(days=1 + AFTER_DAYS)
        e["hi"] = min(e["hi_storm"], ev[k + 1]["first"]) if k + 1 < len(ev) else e["hi_storm"]
        e["features"] = e["first"].strftime("%Y-%m-%d") in FEATURES_FIRST_DAYS

    # refuse to overwrite results that are already there
    already = [(raw.at[e["idx"], ctx["event"]], k) for e in ev for k, h in tgt.items() if raw.at[e["idx"], h].strip() != ""]
    if already and not a.force:
        die("{} WY2020 result cell(s) already have values (for example {}). This script fills empty cells only.\n"
            "  Run it on the table that has them blank, or add --force to overwrite.".format(len(already), already[0]))

    log("  WY2020 winter storms: {} ({} on the 20 mm list)".format(len(ev), sum(1 for e in ev if e["listed"])))

    # data
    gauge = read_obs(need(obs_path, "WY2020 gauge"))
    q = {"retuned": read_qout(need(results / (RUN_STEMS["retuned"] + "_Outlet.qout"), "Re-tuned outlet"))}
    log("  gauge: {} readings, {:%b %d %Y} to {:%b %d %Y}".format(gauge.notna().sum(), gauge.index.min(), gauge.index.max()))
    log("  Re-tuned outlet: {} rows".format(len(q["retuned"])))
    for key in ("original", "kchan05", "kchan02"):
        p = results / (RUN_STEMS[key] + "_Outlet.qout")
        if p.is_file():
            q[key] = read_qout(p)
            log("  {} outlet: {} rows".format(RUN_LABEL[key], len(q[key])))
        else:
            log("  NOT FOUND (those columns stay empty): {}".format(p))
    mrf = None
    mrf_path = Path(a.mrf) if a.mrf else None
    if mrf_path is None:
        found = sorted(results.glob("Sycamore_WY2020_Retuned*.mrf"))
        mrf_path = found[0] if found else None
    if mrf_path is not None and mrf_path.is_file():
        mrf = read_mrf(mrf_path)
        need_cols = [c for c in ("MDGW_mm", "MSMRt_[]", "MAP_mm_hr") if c not in mrf.columns]
        if need_cols:
            die("{} is missing column(s) {}. Columns found: {}".format(mrf_path, need_cols, list(mrf.columns)))
        log("  Re-tuned .mrf: {} ({} rows)".format(mrf_path.name, len(mrf)))
    else:
        log("  NOT FOUND (water table depth, root-zone moisture and outlet flow at start stay empty): results/Sycamore_WY2020_Retuned*.mrf")

    log("\nWindows:")
    log("  {:11s} {:>16s}  {:>16s}  {:>16s}  {}".format("storm", "event window", "to", "storm window to", "fills"))
    for e in ev:
        log("  {:11s} {:%Y-%m-%d %H:%M}  {:%Y-%m-%d %H:%M}  {:%Y-%m-%d %H:%M}  {}{}".format(
            e["name"], e["lo"], e["hi"], e["hi_storm"], "event" + ("+storm window+start state" if e["features"] else ""),
            "" if e["listed"] or e["features"] else " (not on the list)"))
    if a.dry_run:
        log("\nDRY RUN: all inputs were read and the windows listed. Nothing was written.")
        return 0

    # numbers
    notes_added, rows_out = [], []
    bad_writeup, bad_rain = [], []
    for e in ev:
        i, lo, hi, hib = e["idx"], e["lo"], e["hi"], e["hi_storm"]
        hours_exp = int((hi - lo) / pd.Timedelta(hours=1))
        gh = hourly(gauge, lo, hi)
        g = vol(gh)
        m = vol(hourly(q["retuned"], lo, hi))
        bm = e["basin_mm"]
        vals = {"gauge_hours": "{}/{}".format(len(gh), hours_exp), "window_days": str(int((hi - lo) / pd.Timedelta(days=1)))}
        vals["gauge_Mm3"] = fnum(g)
        vals["gauge_ratio_pct"] = fnum(100.0 * g / (MM_TO_MM3 * bm)) if np.isfinite(g) and np.isfinite(bm) and bm > 0 else ""
        vals["model_Mm3"] = fnum(m)
        vals["model_minus_gauge_Mm3"] = fnum(m - g) if np.isfinite(g) else ""
        vals["model_ratio_pct"] = fnum(100.0 * m / (MM_TO_MM3 * bm)) if np.isfinite(bm) and bm > 0 else ""
        for key, col in (("original", "asdrawn_nocl_Mm3"), ("kchan05", "asdrawn_k05_Mm3"), ("kchan02", "asdrawn_k02_Mm3")):
            vals[col] = fnum(vol(hourly(q[key], lo, hi))) if key in q else ""
        gs = ms = float("nan")
        if e["features"]:
            gs = vol(hourly(gauge, lo, hib))
            ms = vol(hourly(q["retuned"], lo, hib))
            vals["gauge_storm_Mm3"], vals["model_storm_Mm3"] = fnum(gs), fnum(ms)
        else:
            vals["gauge_storm_Mm3"] = vals["model_storm_Mm3"] = ""
        wt = sm = of = float("nan")
        if mrf is not None:
            row_ = mrf[mrf.index <= lo]
            if len(row_):
                wt = float(row_["MDGW_mm"].iloc[-1])
                if e["features"]:
                    sm = float(row_["MSMRt_[]"].iloc[-1])
            qq = q["retuned"][q["retuned"].index <= lo]
            if e["features"] and len(qq):
                of = float(qq.iloc[-1])
            # check 2: the model's own rain vs the table's basin rain over the storm days (hours labeled by their END)
            r = mrf["MAP_mm_hr"]
            mm_rain = float(r[(r.index > e["first"]) & (r.index <= e["last"] + pd.Timedelta(days=1))].sum())
            if np.isfinite(bm) and abs(mm_rain - bm) > RAIN_TOL_MM:
                bad_rain.append("{}: model rain {:.2f} mm vs table {:.2f} mm".format(e["name"], mm_rain, bm))
        vals["water_table_depth_mm"] = fnum(wt)
        vals["soil_moisture_root"] = fnum(sm)
        vals["outlet_flow_start_m3s"] = fnum(of)
        for k, v in vals.items():
            raw.at[i, tgt[k]] = v

        # notes where the gauge record is incomplete
        if note_h is not None and len(gh) < 0.97 * hours_exp:
            old = raw.at[i, note_h]
            if not any(w in old for w in ("ropout", "Gauge missing", "gap")):
                add = ("No gauge data in this window; the gauge columns are empty." if len(gh) == 0 else
                       "Gauge has data for {} of {} h in this window, so the gauge volume is a minimum and model / gauge is not like for like.".format(len(gh), hours_exp))
                raw.at[i, note_h] = (old + " " + add).strip()
                notes_added.append(e["name"])

        # check 1 (write-up): the scored storm window of the listed storms, all five series
        if e["listed"] and e["name"] in EXPECTED:
            got = {"gauge": vol(hourly(gauge, lo, hib))}
            for key in ("retuned", "original", "kchan05", "kchan02"):
                got[key] = vol(hourly(q[key], lo, hib)) if key in q else float("nan")
            for key, exp in zip(EXPECTED_ORDER, EXPECTED[e["name"]]):
                if key in ("original", "kchan05", "kchan02") and key not in q:
                    continue
                if not (np.isfinite(got[key]) and abs(got[key] - exp) <= CHECK_TOL):
                    bad_writeup.append("{} {}: here {} vs write-up {:.3f}".format(e["name"], RUN_LABEL.get(key, key), fnum(got[key], 3), exp))
            e["checked"] = True
        rows_out.append((e["name"], hours_exp, len(gh), g, m, vals["asdrawn_nocl_Mm3"], vals["asdrawn_k05_Mm3"], vals["asdrawn_k02_Mm3"], gs, ms, wt))

    log("\nResults written to the new table (volumes in M m3; event window unless marked)")
    log("  {:11s} {:>9s} {:>8s} {:>9s} {:>8s} {:>8s} {:>8s}   {:>10s} {:>10s}   {:>9s}".format(
        "storm", "gauge h", "gauge", "Re-tuned", "Orig", "K 0.5", "K 0.2", "gauge stm", "model stm", "water tbl"))
    for nm, he, gh_, g, m, o, k5, k2, gs, ms, wt in rows_out:
        log("  {:11s} {:>9s} {:>8s} {:>9s} {:>8s} {:>8s} {:>8s}   {:>10s} {:>10s}   {:>9s}".format(
            nm, "{}/{}".format(gh_, he), fnum(g, 3) or "-", fnum(m, 3) or "-", fnum(float(o), 3) if o else "-", fnum(float(k5), 3) if k5 else "-",
            fnum(float(k2), 3) if k2 else "-", fnum(gs, 3) or "-", fnum(ms, 3) or "-", fnum(wt, 1) or "-"))

    n_checked = sum(1 for e in ev if e.get("checked"))
    log("")
    if n_checked:
        if bad_writeup:
            log("CHECK 1 (write-up): MISMATCH ({} value(s) differ by more than {:g} M m3):".format(len(bad_writeup), CHECK_TOL))
            for x in bad_writeup:
                log("  " + x)
        else:
            log("CHECK 1 (write-up): MATCH. Scored-window volumes of {} listed storms (gauge and the runs found) equal the write-up values (within {:g} M m3).".format(n_checked, CHECK_TOL))
    if mrf is not None:
        if bad_rain:
            log("CHECK 2 (model rain vs table rain): DIFFERENCE for {} storm(s):".format(len(bad_rain)))
            for x in bad_rain:
                log("  " + x)
        else:
            log("CHECK 2 (model rain vs table rain): MATCH. The .mrf basin rain equals the table's Basin rain for every storm (within {:g} mm).".format(RAIN_TOL_MM))
    else:
        log("CHECK 2 skipped: no .mrf file.")
    if notes_added:
        log("Notes added to: " + ", ".join(notes_added))

    # check 3: nothing else changed
    changed_cols = set(tgt.values()) | ({note_h} if note_h else set())
    changed_rows = set(e["idx"] for e in ev)
    n_same = n_diff = 0
    for col in raw.columns:
        for i in raw.index:
            same = raw.at[i, col] == original.at[i, col]
            if col in changed_cols and i in changed_rows:
                continue
            n_same += int(same)
            n_diff += int(not same)
    if n_diff:
        die("CHECK 3 failed: {} cell(s) outside the target cells changed. Nothing was written.".format(n_diff))
    n_filled = sum(1 for e in ev for k, h in tgt.items() if raw.at[e["idx"], h].strip() != "")
    log("CHECK 3 (nothing else changed): MATCH. {} other cells are identical to the original table.".format(n_same))
    log("Filled {} cells in {} WY2020 rows.".format(n_filled, len(ev)))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        raw.to_csv(f, index=False)
    log("\nWrote: {}\n  (storm_table.csv was not changed.)".format(out_path))
    return 1 if (bad_writeup or bad_rain) else 0


if __name__ == "__main__":
    sys.exit(main())
