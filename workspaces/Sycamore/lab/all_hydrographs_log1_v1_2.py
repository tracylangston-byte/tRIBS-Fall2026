#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
all_hydrographs_log1_v1_2.py   (v1.2, 2026-11-07)
=================================================

A NEW file. all_hydrographs_log1.py (v1.1) is not changed. Run with no new options, this script does what v1.1 does
(WY2017 and WY2020, gauge and the Re-tuned model) and writes the same figure and CSV names.

WHAT IS NEW IN v1.2
  --wy2020-run     which model run(s) to draw for the WY2020 storms. One or more, separated by commas, up to 4. Default: retuned.
                   Short names:  retuned  original  kchan05  kchan02
                   Any other name is a run name from results/ (the part before _Outlet.qout). The prefix Sycamore_WY2020_ may be left off:
                       S1_Ks1p0   means   results/Sycamore_WY2020_S1_Ks1p0_Outlet.qout
                   Examples:   --wy2020-run S1_Ks1p0                  (S1 instead of Re-tuned)
                               --wy2020-run retuned,S1_Ks1p0          (Re-tuned and S1 together on every WY2020 panel)
  --model-label    the names shown on the figure, in the same order, separated by commas (so no commas inside a name). Keep them short
                   (about 18 characters or fewer): they are used in the legend and in the number panel under each storm.
                   Example:   --model-label "Re-tuned,S1 Ks x1.0"
                   If you leave it off, retuned is called Re-tuned and S1 to S4 get names such as "S1 Ks x1.0".
  The WY2017 panels always show the Re-tuned model (the screens were not run on WY2017). Each run keeps its own color and line style:
  Re-tuned is always red and dashed; the other runs are blue, orange, green and purple.

WHAT CHANGES IN THE CHECK
  The check against the storm table / WY2020 write-up applies to the GAUGE volumes (always) and to the RE-TUNED volumes (when Re-tuned
  is drawn). Any other run has no reference numbers, so it is listed as "not checked". It is compared with the gauge and with Re-tuned,
  not with the table. A MISMATCH on a gauge or Re-tuned volume still means a wrong file, window or run; do not use the figure.

WHERE THE FILES GO
  With the default run the names are the v1.1 names (fig_all_hydrographs_log1.png, all_hydrographs_log1.csv in all_hydro_figs/).
  With any other WY2020 run the run names are added to the file names, for example fig_all_hydrographs_log1_retuned_S1_Ks1p0.png, so the
  Re-tuned figure and CSV are never overwritten.

ONE FIGURE with the 12 storms side by side in time order (6 on top, 6 below): WY2017 first, then WY2020.
Lines drawn: the USGS gauge (black) and the model run(s).

EACH COLUMN HAS THREE PARTS
  top      hourly radar rain over the basin (mm/h, bars hang from the top)   <- basin mean, the same rain the model received
  middle   outlet flow in CFS on a LOG axis: gauge (black) and the model run(s)
  bottom   a panel of numbers for the storm, in this order (rows 1, 6 to 9 come from storm_table.csv; the rest are worked out here).
           With ONE run (the v1.1 layout) there are 10 rows. With more runs, rows 4 and 5 are repeated for each run, named with its label.
             1  Storm rain (mm)            basin rain over the storm days
             2  Runoff ratio (% of rain)   GAUGE volume / rain volume x 100. Rain volume = storm rain (mm) x 0.4246 M m3 per mm (the basin area).
                                           100% would mean every drop of rain left the basin as streamflow. Same definition as the table's gauge runoff ratio.
             3  Volume at gauge (M m3)     USGS gauge volume over the storm window (M m3 = million cubic meters)
             4  tRIBS volume (M m3)        model volume over the same window (one row per run; the label is the run's name when there are several)
             5  tRIBS / gauge volume       model volume / gauge volume (1.00 = perfect, above 1 = model made too much, below 1 = too little).
                                           Shown as n/a when the gauge volume is below --min-gauge-vol (default 0.05 M m3),
                                           because dividing by a near-zero gauge volume gives a meaningless number.
                                           COLOR FLAG on this value (the symbol is a second cue, so it also works in black and white):
                                             green  check mark     within 25% of 1.00 (0.75 to 1.25); --flag-tol 0.10 would mean 10%
                                             orange up triangle    model volume more than 25% HIGH (ratio above 1.25)
                                             blue   down triangle  model volume more than 25% LOW  (ratio below 0.75)
                                           no flag when the ratio is n/a.
             6  Rain above 4 mm/h (mm)     the part of the rain falling faster than 4 mm/h (the table's intensity number)
             7  Radar:gauge, 3 gauges      radar rain / gauge rain, using all 3 rain gauges
             8  Radar:gauge, 2 gauges      the same, without the Mt. Ord gauge
             9  Dry days before            days without rain before the storm
            10  Rain, prior 30 days (mm)   basin radar rain in the 720 hours before 00:00 on the first rain day (the storm itself is not counted);
                                           n/a if fewer than 97% of those hours are in the rain file

THE WINDOW OF EACH PANEL (same rule as the scoring)
  00:00 on the first rain day  to  00:00 four days after the last rain day.
  Before it, EXTRA DAYS are drawn on a gray band for context (default 1 day; --pad-days 0 to 7 changes it).
  Nothing in the gray band is part of the storm window or the volumes.

PANEL NAMES
  "Jan 19-24" is the table's storm "Jan 19-21" (its scored window runs through Jan 24). "Jan 23-24" is a separate storm inside that window,
  so you will see it in both panels. The numbers in the bottom panel of "Jan 19-24" are those of the Jan 19-21 row of the table.

SCALES
  All panels share ONE flow axis (CFS, log) and ONE rain axis, so storms can be compared by eye. Use --free-y to give each panel its own.
  Flow of zero (or below --log-min-cfs, default 0.1 CFS) cannot be drawn on a log axis, so the line stops there. That is not missing data.

UNITS   flow in CFS (cubic feet per second) = m3/s x 35.3147. Volumes are in M m3 = million cubic meters.
        Volumes (rows 2 to 5) are added up over the storm window, so for a storm whose window is longer than the table's event window
        (for example Nov 28-30) the runoff ratio can differ a little from the table's "Gauge runoff ratio (%)" column.

CHECK (printed, nothing is changed)
  Volume = hourly mean flow x 3,600 s over the window (the scoring rule). It is compared with
    WY2017: the table's "Gauge volume, storm window" and "Model volume, storm window" columns (if the table has them), and
    WY2020: the numbers in the WY2020 write-up (gauge and Re-tuned) when the table does not have them.
  It prints MATCH or MISMATCH. Only the gauge and the Re-tuned model are checked (see above).

WHAT IT READS (read-only; no model is run)
  storm_table.csv                                 (short headers such as "first_day" OR the readable ones such as "First rain day")
  storm_features_out/hourly_basin_wy2017.csv and hourly_basin_wy2020.csv      columns interval_start_mst, basin_mean_mm
  results/W5_d012_k10_f0p15_as300_ksx1p03_nocl_Outlet.qout     WY2017 Re-tuned (Time_hr 0 = Oct 1 2016 00:00)
  results/Sycamore_WY2020_<run>_Outlet.qout                    each WY2020 run you ask for (Time_hr 0 = Sep 1 2019 00:00)
  ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv and ..._WY2020.csv    columns datetime, q_cms
WHAT IT WRITES (in all_hydro_figs/ unless you set --out-dir)
  fig_all_hydrographs_log1[_runs].png   the figure
  all_hydrographs_log1[_runs].csv       per storm: window, peaks (CFS), peak rain, and ALL the numbers from the bottom panel
                                        (the raw gauge and model volumes are always kept; a ratio column is blank where the figure says n/a).
                                        Re-tuned columns are named retuned_...; every other run gets columns named <run>_peak_cfs,
                                        <run>_Mm3, <run>_over_gauge_volume and <run>_volume_flag. A column is blank for a panel that does not draw that run.

COMMANDS (from workspaces/Sycamore/lab/)
  python all_hydrographs_log1_v1_2.py --dry-run                 # checks the files and lists the 12 windows; writes nothing
  python all_hydrographs_log1_v1_2.py                           # same as v1.1
  python all_hydrographs_log1_v1_2.py --wy2020-run retuned,S1_Ks1p0 --model-label "Re-tuned,S1 Ks x1.0" --cols 4 \\
         --storms "WY2020:Nov 19-22;WY2020:Nov 28-30;WY2020:Dec 23-28;WY2020:Feb 10-11;WY2020:Feb 22-23;WY2020:Mar 10-13;WY2020:Mar 18-19"
  python all_hydrographs_log1_v1_2.py --pad-days 4              # 4 extra days before each storm (0 to 7)
  python all_hydrographs_log1_v1_2.py --free-y                  # each panel gets its own flow and rain scale
  python all_hydrographs_log1_v1_2.py --log-min-cfs 0.01        # draw smaller flows
  python all_hydrographs_log1_v1_2.py --min-gauge-vol 0.1       # show the tRIBS / gauge ratio only when the gauge volume is at least 0.1 M m3
  python all_hydrographs_log1_v1_2.py --flag-tol 0.10           # green flag only when the volume is within 10% (default 0.25 = 25%)
  python all_hydrographs_log1_v1_2.py --storms "WY2017:Dec 21-25;WY2017:Dec 30-Jan 2;WY2020:Nov 19-22" --cols 3     # other storms (table names)
  Storm spec: YEAR:Event  or  YEAR:Event=Label   (several separated by ;). The year is WY2017 or WY2020. Event is the table's name.

WORDS
  hyetograph  rain bars against time.     log axis  each gridline is 10 times the one below.
  CFS         cubic feet per second.      window    the stretch of time a storm's volume is added up over.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CFS_PER_CMS = 35.3146667
AFTER_DAYS = 3                 # window ends at 00:00, (1 + 3) days after the last rain day's 00:00
FIG_DIR = "all_hydro_figs"
FIG_NAME = "fig_all_hydrographs_log1.png"
CSV_NAME = "all_hydrographs_log1.csv"
MATCH_TOL = 0.005              # M m3
BASIN_MM_TO_MM3 = 0.4246       # M m3 of water per mm of rain over the whole basin (basin area 424.6 km2)
MIN_GAUGE_VOL = 0.05           # M m3; below this the tRIBS / gauge ratio is shown as n/a (see --min-gauge-vol)
FLAG_TOL = 0.25                # tRIBS / gauge volume within +-25% of 1.00 gets the green flag (see --flag-tol)
# flag name -> (symbol, color). Symbols are the second cue so the flag never depends on color alone.
FLAG_STYLE = {"within": ("✓", "#3aa66b"), "high": ("▲", "#e8871e"), "low": ("▼", "#3b82d4")}
PRIOR_DAYS = 30                # days of rain counted before the first rain day
PRIOR_HOURS = PRIOR_DAYS * 24
PRIOR_MIN_HOURS = int(math.ceil(0.97 * PRIOR_HOURS))   # at least 97% of those hours must be in the rain file

INK, INK_2, GRID = "#0b0b0b", "#52514e", "#d9d8d3"
RAIN_BLUE, PAD_GRAY = "#8db4e2", "#ebeae6"
GAUGE_COLOR = "#0b0b0b"

# ---- model runs ------------------------------------------------------------------------------------------------
RETUNED = "retuned"            # the key of the Re-tuned run (the only model run that is checked against the table / write-up)
MAX_RUNS = 4
WY2017_STEM = "W5_d012_k10_f0p15_as300_ksx1p03_nocl"
WY2020_PREFIX = "Sycamore_WY2020_"
# short name -> (results stem, default label)
WY2020_SHORT = {
    "retuned": (WY2020_PREFIX + "Retuned", "Re-tuned"),
    "original": (WY2020_PREFIX + "Original", "Original"),
    "kchan05": (WY2020_PREFIX + "Original_kchan0p5", "Orig. Kchan 0.5"),
    "kchan02": (WY2020_PREFIX + "Original_kchan0p2", "Orig. Kchan 0.2"),
}
NICE_LABELS = {"S1_Ks1p0": "S1 Ks x1.0", "S2_Ks0p5": "S2 Ks x0.5", "S3_Bed0p5": "S3 bedrock +0.5 m", "S4_GW0p5": "S4 start GW 0.5"}
RETUNED_STYLE = ("#d62728", "--", 1.8)                                      # red, dashed (as in v1.1)
OTHER_STYLES = [("#2a78d6", "-", 1.7), ("#e8871e", "-.", 1.7), ("#3aa66b", ":", 2.0), ("#7b4fc4", "-", 1.5)]

YEARS = {
    "WY2017": {"start": "2016-10-01 00:00", "hourly": "hourly_basin_wy2017.csv", "obs": "sycamore_09510200_obs_flow_WY2017.csv"},
    "WY2020": {"start": "2019-09-01 00:00", "hourly": "hourly_basin_wy2020.csv", "obs": "sycamore_09510200_obs_flow_WY2020.csv"},
}
# (year, event name in the storm table, label on the panel)
DEFAULT_STORMS = [
    ("WY2017", "Dec 21-25", "Dec 21-25"), ("WY2017", "Jan 19-21", "Jan 19-24"), ("WY2017", "Jan 23-24", "Jan 23-24"),
    ("WY2017", "Feb 18-19", "Feb 18-19"), ("WY2017", "Feb 27-28", "Feb 27-28"),
    ("WY2020", "Nov 19-22", "Nov 19-22"), ("WY2020", "Nov 28-30", "Nov 28-30"), ("WY2020", "Dec 23-28", "Dec 23-28"),
    ("WY2020", "Feb 10-11", "Feb 10-11"), ("WY2020", "Feb 22-23", "Feb 22-23"),
    ("WY2020", "Mar 10-13", "Mar 10-13"), ("WY2020", "Mar 18-19", "Mar 18-19"),
]
# WY2020 write-up volumes (M m3): gauge and Re-tuned, scored storm windows
EXPECTED_WY2020 = {
    "Nov 19-22": (0.015, 2.800), "Nov 28-30": (1.473, 3.235), "Dec 23-28": (4.069, 0.218), "Feb 10-11": (0.386, 1.767),
    "Feb 22-23": (4.131, 4.839), "Mar 10-13": (5.626, 0.875), "Mar 18-19": (3.317, 1.031),
}
# storm_table.csv comes with short headers or readable headers; both are accepted
COLS = {
    "year": ["year", "Water year"],
    "event": ["event", "Event"],
    "first_day": ["first_day", "First rain day"],
    "last_day": ["last_day", "Last rain day"],
    "dry_days": ["dry_days_before", "Dry days before (d)"],
    "rain_mm": ["basin_mm", "Basin rain (mm)"],
    "x4_mm": ["excess_4_mm", "Rain above 4 mm/h (mm)"],
    "r3": ["ratio_3g", "Radar / gauge (3) (ratio)"],
    "r2": ["ratio_2g", "Radar / gauge (2) (ratio)"],
    "gauge_storm": ["gauge_storm_Mm3", "Gauge volume, storm window (M m3)"],
    "model_storm": ["model_storm_Mm3", "Model volume, storm window (M m3)"],
}
REQUIRED = ["year", "event", "first_day", "last_day"]


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    sys.exit("\nSTOPPED: " + msg + "\n")


def need(path: Path, what: str) -> Path:
    if not Path(path).is_file():
        die("cannot find the {} file:\n    {}\n  Run this from workspaces/Sycamore/lab/, or point to the file with the matching option "
            "(see the top of this script).".format(what, path))
    return Path(path)


def tag_of(key: str) -> str:
    """A run key made safe for a column or file name."""
    return re.sub(r"[^A-Za-z0-9_]+", "_", key).strip("_")


# ------------------------------------------------------------------------------------------------ the model runs
def parse_wy2020_runs(spec: str, labels):
    names = [s.strip() for s in (spec or RETUNED).split(",") if s.strip()]
    if not names:
        die("--wy2020-run needs at least one run, for example  --wy2020-run retuned,S1_Ks1p0")
    if len(names) > MAX_RUNS:
        die("--wy2020-run takes at most {} runs (you gave {}).".format(MAX_RUNS, len(names)))
    runs, seen = [], set()
    for nm in names:
        low = nm.lower()
        if low in WY2020_SHORT:
            key = low
            stem, label = WY2020_SHORT[low]
        else:
            stem = nm if nm.startswith("Sycamore_") else WY2020_PREFIX + nm
            if stem.endswith("_Outlet.qout"):
                stem = stem[:-len("_Outlet.qout")]
            key = stem[len(WY2020_PREFIX):] if stem.startswith(WY2020_PREFIX) else stem
            label = NICE_LABELS.get(key, key)
            for sk, (sstem, slabel) in WY2020_SHORT.items():      # a full stem that is one of the short names counts as that name
                if stem == sstem:
                    key, label = sk, slabel
        if key in seen:
            die("--wy2020-run lists '{}' twice.".format(nm))
        seen.add(key)
        runs.append({"key": key, "stem": stem, "label": label})
    if labels:
        labs = [s.strip() for s in labels.split(",")]
        if len(labs) != len(runs) or any(not x for x in labs):
            die("--model-label needs one name for each run, in the same order, separated by commas.\n"
                "  You gave {} name(s) for {} run(s): runs = {}".format(len([x for x in labs if x]), len(runs), ", ".join(r["key"] for r in runs)))
        for r, lab in zip(runs, labs):
            r["label"] = lab
    others = iter(OTHER_STYLES)
    for r in runs:
        r["color"], r["ls"], r["lw"] = RETUNED_STYLE if r["key"] == RETUNED else next(others)
    return runs


def need_run(results: Path, stem: str, what: str) -> Path:
    p = results / (stem + "_Outlet.qout")
    if p.is_file():
        return p
    have = sorted(f.name[:-len("_Outlet.qout")] for f in results.glob("*_Outlet.qout")) if results.is_dir() else []
    msg = "cannot find the outlet file for {}:\n    {}".format(what, p)
    if have:
        msg += "\n  Outlet files found in {}:\n    ".format(results) + "\n    ".join(have[:25])
        msg += "\n  Give the run name without _Outlet.qout (the prefix {} may be left off).".format(WY2020_PREFIX)
    else:
        msg += "\n  Run this from workspaces/Sycamore/lab/, or point to the folder with --results."
    die(msg)


# ------------------------------------------------------------------------------------------------ reading
def read_table(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    t = pd.DataFrame(index=raw.index)
    missing = []
    for want, options in COLS.items():
        hit = [c for c in options if c in raw.columns]
        if hit:
            t[want] = raw[hit[0]]
        elif want in REQUIRED:
            missing.append(want)
        else:
            t[want] = np.nan
    if missing:
        die("{} is missing column(s) for {} (looked for {}).\n  All columns found: {}".format(
            path, missing, [COLS[m] for m in missing], list(raw.columns)))
    t["year"] = t["year"].astype(str).str.strip()
    t["event"] = t["event"].astype(str).str.replace("–", "-").str.replace("—", "-").str.strip()
    t["first_day"] = pd.to_datetime(t["first_day"], errors="coerce")
    t["last_day"] = pd.to_datetime(t["last_day"], errors="coerce")
    for c in ("dry_days", "rain_mm", "x4_mm", "r3", "r2", "gauge_storm", "model_storm"):
        t[c] = pd.to_numeric(t[c], errors="coerce")
    return t


def read_qout(path: Path, start: str) -> pd.Series:
    raw = pd.read_csv(path).iloc[:, :2]
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    idx = pd.Timestamp(start) + pd.to_timedelta(raw["Time_hr"].to_numpy(float), unit="h")
    s = pd.Series(raw["Q"].to_numpy(float), index=idx)
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_obs(path: Path) -> pd.Series:
    o = pd.read_csv(path)
    for c in ("datetime", "q_cms"):
        if c not in o.columns:
            die("{} has no column '{}'. Columns found: {}".format(path, c, list(o.columns)[:10]))
    o["datetime"] = pd.to_datetime(o["datetime"])
    s = pd.Series(pd.to_numeric(o["q_cms"], errors="coerce").to_numpy(float), index=o["datetime"])
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_rain(path: Path) -> pd.Series:
    h = pd.read_csv(path)
    for c in ("interval_start_mst", "basin_mean_mm"):
        if c not in h.columns:
            die("{} has no column '{}'. Columns found: {}".format(path, c, list(h.columns)[:12]))
    h["interval_start_mst"] = pd.to_datetime(h["interval_start_mst"])
    s = pd.Series(pd.to_numeric(h["basin_mean_mm"], errors="coerce").to_numpy(float), index=h["interval_start_mst"])
    return s[~s.index.duplicated(keep="last")].sort_index()


# ------------------------------------------------------------------------------------------------ numbers
def volume(q: pd.Series, lo, hi):
    """Million m3 from the hourly mean flow x 3,600 s over [lo, hi), counting only hours with data. Returns (volume, hours)."""
    s = q[(q.index >= lo) & (q.index < hi)]
    h = s.resample("60min").mean().dropna()
    return (float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")), int(len(h))


def peak(q: pd.Series, lo, hi) -> float:
    s = q[(q.index >= lo) & (q.index < hi)].dropna()
    return float(s.max()) if len(s) else float("nan")


def logsafe(y, floor: float) -> np.ndarray:
    y = np.asarray(y, float).copy()
    y[~(y > floor)] = np.nan
    return y


def fmt(v, nd: int = 2) -> str:
    return "n/a" if v is None or not np.isfinite(v) else "{:.{n}f}".format(v, n=nd)


def fmt_pct(v) -> str:
    """Percent with a % sign; two decimals below 1% so small values do not all print as 0.0%."""
    if v is None or not np.isfinite(v):
        return "n/a"
    return "{:.2f}%".format(v) if abs(v) < 1.0 else "{:.1f}%".format(v)


def fmt_ratio(v) -> str:
    """A ratio with about three digits: 0.43, 4.6, 46, 460."""
    if v is None or not np.isfinite(v):
        return "n/a"
    if abs(v) < 10:
        return "{:.2f}".format(v)
    if abs(v) < 100:
        return "{:.1f}".format(v)
    return "{:.0f}".format(v)


def volume_flag(ratio, tol):
    """'within' if tRIBS/gauge volume is inside 1 +- tol, 'high' if the model made more, 'low' if less, None if no ratio."""
    if ratio is None or not np.isfinite(ratio):
        return None
    if ratio > 1.0 + tol + 1e-9:
        return "high"
    if ratio < 1.0 - tol - 1e-9:
        return "low"
    return "within"


def parse_storms(spec: str):
    out = []
    for part in [p.strip() for p in spec.split(";") if p.strip()]:
        if ":" not in part:
            die("--storms item '{}' needs a year, like WY2017:Dec 21-25".format(part))
        yr, rest = part.split(":", 1)
        yr = yr.strip().upper()
        yr = yr if yr.startswith("WY") else "WY" + yr
        if yr not in YEARS:
            die("--storms year '{}' must be WY2017 or WY2020.".format(yr))
        ev, _, lab = rest.partition("=")
        ev = ev.strip().replace("–", "-")
        out.append((yr, ev, lab.strip() or ev))
    return out


# ------------------------------------------------------------------------------------------------ drawing
def tidy(ax) -> None:
    ax.grid(alpha=0.35, color=GRID, which="major", lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=7.5)


def flow_tick(v, _pos) -> str:
    return "{:,.0f}".format(v) if v >= 1 else "{:g}".format(v)


def number_rows(it):
    """The rows of the number panel for one storm: a list of (label, value text), and {row index: run} for the rows that get a flag."""
    runs = it["runs"]
    single_retuned = len(runs) == 1 and runs[0]["key"] == RETUNED     # the v1.1 wording
    rows = [("Storm rain (mm)", fmt(it["rain_mm"], 1)),
            ("Runoff ratio (% of rain)", fmt_pct(it["runoff_pct"])),
            ("Volume at gauge (M m³)", fmt(it["gauge_vol"], 3))]
    flag_rows = {}
    for r in runs:
        lab_v, lab_r = ("tRIBS volume (M m³)", "tRIBS / gauge volume") if single_retuned else \
                       ("{} vol (M m³)".format(r["label"]), "{} / gauge".format(r["label"]))
        rows.append((lab_v, fmt(r["vol"], 3)))
        flag_rows[len(rows)] = r
        rows.append((lab_r, fmt_ratio(r["ratio"])))
    rows += [("Rain above 4 mm/h (mm)", fmt(it["x4_mm"], 2)),
             ("Radar:gauge, 3 gauges", fmt(it["r3"], 2)),
             ("Radar:gauge, 2 gauges", fmt(it["r2"], 2)),
             ("Dry days before", fmt(it["dry_days"], 0)),
             ("Rain, prior 30 days (mm)", fmt(it["prior_mm"], 1))]
    return rows, flag_rows


def make_figure(out: Path, items: list, a, flow_top: float, rain_top: float) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    n = len(items)
    cols = max(1, min(a.cols, n))
    blocks = int(math.ceil(n / float(cols)))
    n_max = max(len(number_rows(it)[0]) for it in items)           # 10 rows with one run; two more for each extra run
    info_ratio = 2.3 * max(1.0, n_max / 10.0)                       # taller number panel when there are more rows
    block_h = 6.1 * (0.8 + 2.8 + 0.55 + info_ratio) / 6.45
    step = 0.1 if n_max <= 10 else 0.9 / (n_max - 1)

    # legend entries: the gauge, then each distinct run once, in the order first drawn
    run_handles, seen = [], set()
    for it in items:
        for r in it["runs"]:
            sig = (r["key"], r["label"])
            if sig not in seen:
                seen.add(sig)
                run_handles.append(Line2D([], [], color=r["color"], lw=r["lw"], ls=r["ls"], label="{} model".format(r["label"])))
    handles = [Line2D([], [], color=GAUGE_COLOR, lw=2.2, label="USGS gauge 09510200")] + run_handles + [
               Patch(color=RAIN_BLUE, label="hourly radar rain, basin mean"),
               Patch(color=PAD_GRAY, label="extra day(s) before the storm window ({:d})".format(a.pad_days)),
               Patch(color=FLAG_STYLE["within"][1], label="{} tRIBS / gauge volume within {:g}% of 1.00".format(FLAG_STYLE["within"][0], 100 * a.flag_tol)),
               Patch(color=FLAG_STYLE["high"][1], label="{} model volume more than {:g}% high".format(FLAG_STYLE["high"][0], 100 * a.flag_tol)),
               Patch(color=FLAG_STYLE["low"][1], label="{} model volume more than {:g}% low".format(FLAG_STYLE["low"][0], 100 * a.flag_tol))]
    leg_rows = int(math.ceil(len(handles) / 4.0))
    fig_w, fig_h = 3.55 * cols + 1.1, block_h * blocks + 1.9

    # footnote: two paragraphs, wrapped to the figure width (about 15 characters per inch at this font size). In the 6-column figure
    # nothing wraps, so that figure is the same as v1.1's. In a narrower figure the wrapped lines need more room, so the legend and
    # the panels are moved up by the extra lines.
    foot_paras = [
        "Flow of zero or below {:g} CFS cannot be drawn on a log axis, so a line stops there (not missing data). 'Jan 19-24' is the table's Jan 19-21 storm "
        "(window runs through Jan 24). 2 gauges = without the Mt. Ord gauge.".format(a.log_min_cfs),
        "Runoff ratio = gauge volume / (storm rain x {:.4f} M m³ per mm). Volumes are over the storm window. tRIBS / gauge volume: 1 = perfect, above 1 = model made too much, "
        "below 1 = too little; n/a = gauge volume under {:g} M m³ (or not computed). Prior 30 days = basin radar rain before the first rain day.".format(
            BASIN_MM_TO_MM3, a.min_gauge_vol)]
    wrap_w = max(60, int((fig_w - 0.6) * 15))
    foot_lines = []
    for para in foot_paras:
        foot_lines += textwrap.wrap(para, width=wrap_w) or [""]
    foot = "\n".join(foot_lines)
    foot_extra = max(0, len(foot_lines) - 2) * 0.135 / fig_h          # 0.135 inch per extra line
    bottom = 0.105 + 0.022 * max(0, leg_rows - 2) + foot_extra

    fig = plt.figure(figsize=(fig_w, fig_h))
    # a short figure (one band of panels) keeps the same space in INCHES for the title lines as the 12-storm figure has, so the
    # subtitle does not run into the panel titles; a tall figure uses the v1.1 fractions unchanged
    short = fig_h < 12.0
    top_frac = 0.912 if not short else 1.0 - 1.241 / fig_h
    y_title = 0.985 if not short else 1.0 - 0.2115 / fig_h
    y_sub = 0.947 if not short else 1.0 - 0.747 / fig_h
    outer = fig.add_gridspec(blocks, 1, left=0.05, right=0.995, top=top_frac, bottom=bottom, hspace=0.17)

    for b in range(blocks):
        inner = outer[b].subgridspec(4, cols, height_ratios=[0.8, 2.8, 0.55, info_ratio], hspace=0.06, wspace=0.10)
        for c in range(cols):
            k = b * cols + c
            if k >= n:
                continue
            it = items[k]
            ax_r = fig.add_subplot(inner[0, c])
            ax_q = fig.add_subplot(inner[1, c], sharex=ax_r)
            ax_i = fig.add_subplot(inner[3, c])
            lo, hi, lo_plot = it["lo"], it["hi"], it["lo_plot"]

            # rain
            r = it["rain"]
            ax_r.axvspan(lo_plot, lo, color=PAD_GRAY, lw=0, zorder=0)
            ax_r.bar(r.index, r.values, width=1.0 / 24.0, align="edge", color=RAIN_BLUE, zorder=2)
            top_r = rain_top if not a.free_y else max(1.0, math.ceil(float(np.nanmax(r.values))) if len(r) and np.isfinite(r.values).any() else 1.0)
            ax_r.set_ylim(top_r, 0)
            ax_r.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
            ax_r.set_xlim(lo_plot, hi)
            ax_r.set_title("{}, {:%Y}".format(it["label"], it["first"]), loc="left", fontsize=10.5, color=INK, pad=4)
            ax_r.tick_params(labelbottom=False)

            # flow
            ax_q.axvspan(lo_plot, lo, color=PAD_GRAY, lw=0, zorder=0)
            g = it["gauge_cfs"]
            ax_q.plot(g.index, logsafe(g.values, a.log_min_cfs), color=GAUGE_COLOR, lw=2.2, label="USGS gauge 09510200", zorder=3)
            for z, run in enumerate(it["runs"]):
                m = run["cfs"]
                ax_q.plot(m.index, logsafe(m.values, a.log_min_cfs), color=run["color"], lw=run["lw"], ls=run["ls"],
                          label="{} model".format(run["label"]), zorder=4 + z)
            ax_q.set_yscale("log")
            if a.free_y:
                peaks = [it["gauge_peak_cfs"]] + [run["peak_cfs"] for run in it["runs"]]
                top = max([x for x in peaks if np.isfinite(x)] + [a.log_min_cfs * 10]) * 2.5
            else:
                top = flow_top
            ax_q.set_ylim(a.log_min_cfs, top)
            ax_q.yaxis.set_major_formatter(FuncFormatter(flow_tick))
            loc = mdates.AutoDateLocator(minticks=3, maxticks=5)
            ax_q.xaxis.set_major_locator(loc)
            ax_q.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
            ax_q.set_xlim(lo_plot, hi)
            for ax in (ax_r, ax_q):
                tidy(ax)
            if c == 0 or a.free_y:
                ax_r.set_ylabel("Rain\n(mm/h)", fontsize=8, color=INK_2)
                ax_q.set_ylabel("Outlet flow (CFS), log", fontsize=8, color=INK_2)
            else:
                ax_r.tick_params(labelleft=False)
                ax_q.tick_params(labelleft=False)

            # numbers
            ax_i.set_xticks([])
            ax_i.set_yticks([])
            for s in ax_i.spines.values():
                s.set_color(GRID)
            lines, flag_rows = number_rows(it)
            for j, (lab, val) in enumerate(lines):
                y = 0.95 - j * step
                ax_i.text(0.04, y, lab, transform=ax_i.transAxes, ha="left", va="center", fontsize=8, color=INK_2)
                kw = {}
                run_here = flag_rows.get(j)                # flag_rows is keyed by the row number of each run's ratio row
                if run_here is not None and run_here.get("flag"):
                    sym, col = FLAG_STYLE[run_here["flag"]]
                    val = sym + " " + val
                    kw = {"bbox": dict(boxstyle="round,pad=0.25", fc=col, ec="none")}
                ax_i.text(0.96, y, val, transform=ax_i.transAxes, ha="right", va="center", fontsize=8.5, color=INK,
                          fontweight="bold", **kw)
            # thin lines between the three groups: water volumes | rain intensity and radar | antecedent conditions
            g1_end = 2 + 2 * len(it["runs"])
            for after_row in (g1_end, g1_end + 3):
                yy = 0.95 - (after_row + 0.5) * step
                ax_i.plot([0.02, 0.98], [yy, yy], transform=ax_i.transAxes, color=GRID, lw=0.8, clip_on=False)

    labels = []
    for it in items:
        for r in it["runs"]:
            if r["label"] not in labels:
                labels.append(r["label"])
    who = "Re-tuned model" if labels == ["Re-tuned"] else "models ({})".format(", ".join(labels))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.040 + foot_extra), ncol=4, frameon=False, fontsize=9.5)
    title = "Lower Sycamore Wash: gauge and {}, storm by storm (outlet flow in CFS, log axis)".format(who)
    if len(title) > int((fig_w - 0.5) * 9.5):          # too long for a narrow figure with several runs: the legend names the runs anyway
        title = "Lower Sycamore Wash: gauge and models, storm by storm (outlet flow in CFS, log axis)"
    fig.suptitle(title, x=0.012, y=y_title, ha="left", fontsize=14, color=INK)
    years_txt = " then ".join(sorted(set(it["year"] for it in items)))       # 'WY2017 then WY2020', or just the one year drawn
    sub = ("{}, in time order. {} flow and rain scales in every panel. "
           "Storm window: 00:00 on the first rain day to 00:00 four days after the last rain day.").format(
               years_txt, "Separate" if a.free_y else "Shared")
    sub = "\n".join(textwrap.wrap(sub, width=max(40, int((fig_w - 0.4) * 15))))       # no wrapping in the 6-column figure
    fig.text(0.012, y_sub, sub, fontsize=9.5, color=INK_2, ha="left", va="top")
    fig.text(0.012, 0.006, foot, fontsize=7.8, color=INK_2, ha="left", va="bottom")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=a.dpi)
    plt.close(fig)
    return str(out)


# ------------------------------------------------------------------------------------------------ main
def pad_type(s: str) -> int:
    try:
        v = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError("--pad-days must be a whole number from 0 to 7")
    if v < 0 or v > 7:
        raise argparse.ArgumentTypeError("--pad-days must be from 0 to 7 (you gave {})".format(v))
    return v


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="12 storm hydrographs in one figure: gauge and model run(s), flow in CFS on a log axis (see the top of this file).")
    ap.add_argument("--root", default=str(HERE), help="the lab folder (default: the folder this script is in)")
    ap.add_argument("--table", default=None, help="storm_table.csv (default <root>/storm_table.csv)")
    ap.add_argument("--hourly-dir", default=None, help="folder with hourly_basin_wy2017.csv and ..._wy2020.csv (default <root>/storm_features_out)")
    ap.add_argument("--results", default=None, help="folder with the *_Outlet.qout files (default <root>/results)")
    ap.add_argument("--obs-dir", default=None, help="folder with the two USGS files (default <root>/../init_data/obs)")
    ap.add_argument("--out-dir", default=None, help="where the figure and CSV go (default <root>/" + FIG_DIR + ")")
    ap.add_argument("--wy2020-run", default=RETUNED,
                    help="model run(s) to draw for the WY2020 storms, comma separated, up to {}: retuned, original, kchan05, kchan02, or a run name "
                         "such as S1_Ks1p0 (default: retuned)".format(MAX_RUNS))
    ap.add_argument("--model-label", default=None,
                    help="names to show for those runs, same order, comma separated, short (for example \"Re-tuned,S1 Ks x1.0\")")
    ap.add_argument("--storms", default=None, help="other storms, 'YEAR:Event=Label;YEAR:Event' (default: the 12 asked for)")
    ap.add_argument("--cols", type=int, default=6, help="panels per row (default 6)")
    ap.add_argument("--pad-days", type=pad_type, default=1, help="extra days drawn before each storm window, 0 to 7 (default 1)")
    ap.add_argument("--log-min-cfs", type=float, default=0.1, help="lowest flow on the log axis, CFS (default 0.1)")
    ap.add_argument("--free-y", action="store_true", help="each panel gets its own flow and rain scale (default: shared)")
    ap.add_argument("--flag-tol", type=float, default=FLAG_TOL,
                    help="the tRIBS / gauge volume ratio gets the green flag if it is within 1 +- this fraction (default {:g} = 25%%)".format(FLAG_TOL))
    ap.add_argument("--min-gauge-vol", type=float, default=MIN_GAUGE_VOL,
                    help="the tRIBS / gauge volume ratio is shown only if the gauge volume is at least this many M m3 (default {:g})".format(MIN_GAUGE_VOL))
    ap.add_argument("--dpi", type=int, default=130)
    ap.add_argument("--dry-run", action="store_true", help="check the files and list the windows; write nothing")
    a = ap.parse_args(argv)
    if a.log_min_cfs <= 0:
        die("--log-min-cfs must be above zero.")
    if not (0 < a.flag_tol < 1):
        die("--flag-tol must be between 0 and 1 (0.25 means within 25%).")
    if a.min_gauge_vol <= 0:
        die("--min-gauge-vol must be above zero (a ratio cannot be made with a gauge volume of zero).")

    root = Path(a.root)
    table = Path(a.table) if a.table else root / "storm_table.csv"
    hourly_dir = Path(a.hourly_dir) if a.hourly_dir else root / "storm_features_out"
    results = Path(a.results) if a.results else root / "results"
    obs_dir = Path(a.obs_dir) if a.obs_dir else root.parent / "init_data" / "obs"
    out_dir = Path(a.out_dir) if a.out_dir else root / FIG_DIR
    storms = parse_storms(a.storms) if a.storms else DEFAULT_STORMS

    wy20_runs = parse_wy2020_runs(a.wy2020_run, a.model_label)
    retuned20 = [r for r in wy20_runs if r["key"] == RETUNED]
    wy17_label = retuned20[0]["label"] if retuned20 else "Re-tuned"
    runs_by_year = {
        "WY2017": [{"key": RETUNED, "stem": WY2017_STEM, "label": wy17_label,
                    "color": RETUNED_STYLE[0], "ls": RETUNED_STYLE[1], "lw": RETUNED_STYLE[2]}],
        "WY2020": wy20_runs,
    }

    log("\nall_hydrographs_log1_v1_2: reading files")
    tab = read_table(need(table, "storm table"))
    log("  storm table: {} rows".format(len(tab)))
    need_years = sorted(set(s[0] for s in storms))
    if "WY2020" not in need_years and (a.wy2020_run != RETUNED or a.model_label):
        log("  note: no WY2020 storms are selected, so --wy2020-run and --model-label are not used.")
    rain, gauge, model = {}, {}, {}
    for yr in need_years:
        info = YEARS[yr]
        rain[yr] = read_rain(need(hourly_dir / info["hourly"], "hourly rain " + yr))
        gauge[yr] = read_obs(need(obs_dir / info["obs"], "USGS gauge " + yr))
        log("  {}: rain {} h ({:%b %d %Y} to {:%b %d %Y}); gauge {} readings ({:%b %d %Y} to {:%b %d %Y})".format(
            yr, len(rain[yr]), rain[yr].index.min(), rain[yr].index.max(), gauge[yr].notna().sum(), gauge[yr].index.min(), gauge[yr].index.max()))
        for r in runs_by_year[yr]:
            p = need_run(results, r["stem"], "{} ({})".format(r["label"], yr))
            model[(yr, r["key"])] = read_qout(p, info["start"])
            log("      model run {} [{}]: {} rows".format(r["label"], r["stem"], len(model[(yr, r["key"])])))

    # one table row per storm
    items = []
    for yr, ev, label in storms:
        rows = tab[(tab["year"].str.contains(yr[2:])) & (tab["event"] == ev)]
        if rows.empty:
            names = list(tab.loc[tab["year"].str.contains(yr[2:]), "event"])
            die("storm '{}' is not in the {} rows of {}.\n  Names found: {}".format(ev, yr, table, ", ".join(names)))
        r = rows.iloc[0]
        if pd.isna(r["first_day"]) or pd.isna(r["last_day"]):
            die("the first/last rain day of '{}' could not be read as a date in {}.".format(ev, table))
        first, last = pd.Timestamp(r["first_day"]), pd.Timestamp(r["last_day"])
        lo = first
        hi = last + pd.Timedelta(days=1 + AFTER_DAYS)
        # basin radar rain in the 30 days (720 hours) before 00:00 on the first rain day; the storm itself is not counted
        pr = rain[yr][(rain[yr].index >= first - pd.Timedelta(days=PRIOR_DAYS)) & (rain[yr].index < first)].dropna()
        prior_hours = int(len(pr))
        prior_mm = float(pr.sum()) if prior_hours >= PRIOR_MIN_HOURS else float("nan")
        it = {"year": yr, "event": ev, "label": label, "first": first, "last": last, "lo": lo, "hi": hi,
              "lo_plot": lo - pd.Timedelta(days=a.pad_days),
              "dry_days": float(r["dry_days"]), "rain_mm": float(r["rain_mm"]), "x4_mm": float(r["x4_mm"]),
              "r3": float(r["r3"]), "r2": float(r["r2"]),
              "prior_mm": prior_mm, "prior_hours": prior_hours,
              "tab_gauge": float(r["gauge_storm"]), "tab_model": float(r["model_storm"]),
              "runs": [dict(x) for x in runs_by_year[yr]]}
        items.append(it)
    items.sort(key=lambda d: d["first"])

    log("\nWindows ({} storms, {} extra day(s) before each):".format(len(items), a.pad_days))
    for it in items:
        log("  {:10s} {:9s} window {:%Y-%m-%d %H:%M} to {:%Y-%m-%d %H:%M}   drawn from {:%Y-%m-%d}   model: {}".format(
            it["label"], it["year"], it["lo"], it["hi"], it["lo_plot"], ", ".join(r["label"] for r in it["runs"])))
        if it["lo_plot"] < gauge[it["year"]].index.min() or it["lo_plot"] < rain[it["year"]].index.min():
            log("    NOTE: the extra days start before the gauge or rain record; the first part will be empty.")
        log("    rain in the {} days before the first rain day: {} mm ({} of {} hours in the rain file)".format(
            PRIOR_DAYS, fmt(it["prior_mm"], 1), it["prior_hours"], PRIOR_HOURS))
        if not np.isfinite(it["prior_mm"]):
            log("    NOTE: under 97% of those {} hours are in the rain file, so 'Rain, prior 30 days' will show n/a for this storm.".format(PRIOR_HOURS))
    if a.dry_run:
        log("\nDRY RUN: all inputs were read and the windows listed. Nothing was written.")
        return 0

    # slices, peaks, check
    bad = []
    n_runs_max = max(len(it["runs"]) for it in items)
    run_labels = []
    for it in items:
        for r in it["runs"]:
            if r["label"] not in run_labels:
                run_labels.append(r["label"])
    log("\nVolumes over the window (M m3), scoring rule; peaks in CFS are description only")
    head = "  {:10s} {:>9s} ".format("storm", "gauge") + " ".join("{:>10.10s}".format(x) for x in run_labels) + "   " + \
           "{:>10s} ".format("gauge pk") + " ".join("{:>10.10s}".format(x) for x in run_labels) + "   {:>9s} ".format("runoff %") + \
           " ".join("{:>12.12s}".format(x[:9] + " /g") for x in run_labels)
    log(head)
    n_gauge_checked, n_retuned_checked, unchecked_runs = 0, 0, []
    for it in items:
        yr = it["year"]
        g_cms = gauge[yr]
        lo, hi, lo_plot = it["lo"], it["hi"], it["lo_plot"]
        it["gauge_cfs"] = g_cms[(g_cms.index >= lo_plot) & (g_cms.index <= hi)] * CFS_PER_CMS
        it["rain"] = rain[yr][(rain[yr].index >= lo_plot) & (rain[yr].index < hi)]
        it["gauge_peak_cfs"] = peak(g_cms, lo, hi) * CFS_PER_CMS
        it["rain_peak"] = float(it["rain"].max()) if len(it["rain"]) else float("nan")
        gv, gh = volume(g_cms, lo, hi)
        it["gauge_vol"] = gv
        # runoff ratio (gauge): volume / (rain depth x basin area), in percent
        rain_vol = it["rain_mm"] * BASIN_MM_TO_MM3
        it["runoff_pct"] = 100.0 * gv / rain_vol if (np.isfinite(gv) and np.isfinite(rain_vol) and rain_vol > 0) else float("nan")
        for run in it["runs"]:
            m_cms = model[(yr, run["key"])]
            run["cfs"] = m_cms[(m_cms.index >= lo_plot) & (m_cms.index <= hi)] * CFS_PER_CMS
            run["peak_cfs"] = peak(m_cms, lo, hi) * CFS_PER_CMS
            mv, mh = volume(m_cms, lo, hi)
            run["vol"] = mv
            # tRIBS / gauge volume: only when the gauge volume is big enough for a ratio to mean something
            run["ratio"] = mv / gv if (np.isfinite(gv) and np.isfinite(mv) and gv >= a.min_gauge_vol) else float("nan")
            run["flag"] = volume_flag(run["ratio"], a.flag_tol)
        cells = {r["label"]: r for r in it["runs"]}
        line = "  {:10s} {:>9s} ".format(it["label"], fmt(gv, 3))
        line += " ".join("{:>10s}".format(fmt(cells[x]["vol"], 3) if x in cells else "") for x in run_labels) + "   "
        line += "{:>10s} ".format(fmt(it["gauge_peak_cfs"], 1))
        line += " ".join("{:>10s}".format(fmt(cells[x]["peak_cfs"], 1) if x in cells else "") for x in run_labels) + "   "
        line += "{:>9s} ".format(fmt_pct(it["runoff_pct"]))
        line += " ".join("{:>12s}".format(fmt_ratio(cells[x]["ratio"]) if x in cells else "") for x in run_labels)
        log(line)
        # check: the gauge always, the Re-tuned run when it is drawn
        exp = None
        if np.isfinite(it["tab_gauge"]) and np.isfinite(it["tab_model"]):
            exp = (it["tab_gauge"], it["tab_model"], "the table")
        elif yr == "WY2020" and it["event"] in EXPECTED_WY2020:
            exp = EXPECTED_WY2020[it["event"]] + ("the WY2020 write-up",)
        if exp is not None:
            n_gauge_checked += 1
            if not (np.isfinite(gv) and abs(gv - exp[0]) <= MATCH_TOL):
                bad.append("{} gauge: here {} vs {} {:.3f}".format(it["label"], fmt(gv, 3), exp[2], exp[0]))
            for run in it["runs"]:
                if run["key"] == RETUNED:
                    n_retuned_checked += 1
                    if not (np.isfinite(run["vol"]) and abs(run["vol"] - exp[1]) <= MATCH_TOL):
                        bad.append("{} {}: here {} vs {} {:.3f}".format(it["label"], run["label"], fmt(run["vol"], 3), exp[2], exp[1]))
        for run in it["runs"]:
            if run["key"] != RETUNED and run["label"] not in unchecked_runs:
                unchecked_runs.append(run["label"])
        if gh < 0.97 * int((hi - lo) / pd.Timedelta(hours=1)):
            log("    note: the gauge has data for {} of {} hours in this window; its volume counts only those hours.".format(
                gh, int((hi - lo) / pd.Timedelta(hours=1))))
    if n_gauge_checked:
        if bad:
            log("\nCHECK: MISMATCH ({} value(s) differ by more than {:g} M m3):".format(len(bad), MATCH_TOL))
            for x in bad:
                log("  " + x)
            log("  The figure is still drawn, but do not use it until this is explained.")
        else:
            log("\nCHECK: MATCH. Gauge volumes for {} of {} storms and Re-tuned volumes for {} panels equal the table / write-up values (within {:g} M m3).".format(
                n_gauge_checked, len(items), n_retuned_checked, MATCH_TOL))
    if n_gauge_checked < len(items):
        log("  ({} storm(s) have nothing to check against and are not checked.)".format(len(items) - n_gauge_checked))
    if unchecked_runs:
        log("  Not checked (no reference values exist for them): {}. Compare these with the gauge and with Re-tuned.".format(", ".join(unchecked_runs)))

    # shared scales
    tops = [x for it in items for x in [it["gauge_peak_cfs"]] + [r["peak_cfs"] for r in it["runs"]] if np.isfinite(x)]
    flow_top = (max(tops) if tops else 100.0) * 2.5
    rmax = [it["rain_peak"] for it in items if np.isfinite(it["rain_peak"])]
    rain_top = max(1.0, float(math.ceil(max(rmax)))) if rmax else 1.0

    # file names: the v1.1 names when the only WY2020 run is Re-tuned, otherwise the run names are added
    used20 = "WY2020" in need_years
    suffix = "" if (not used20 or [r["key"] for r in wy20_runs] == [RETUNED]) else "_" + "_".join(tag_of(r["key"]) for r in wy20_runs)
    fig_path = make_figure(out_dir / FIG_NAME.replace(".png", suffix + ".png"), items, a, flow_top, rain_top)

    # CSV: the v1.1 columns for Re-tuned, then four columns for each other run
    extra_keys = []
    for it in items:
        for r in it["runs"]:
            if r["key"] != RETUNED and r["key"] not in extra_keys:
                extra_keys.append(r["key"])
    flag_text = {"within": "within {:g}%".format(100 * a.flag_tol), "high": "model too high", "low": "model too low"}
    any_retuned = any(r["key"] == RETUNED for it in items for r in it["runs"])
    out_rows = []
    for it in items:
        by_key = {r["key"]: r for r in it["runs"]}
        rt = by_key.get(RETUNED)
        row = {"storm": it["label"], "table_event": it["event"], "year": it["year"], "window_start": it["lo"], "window_end": it["hi"],
               "extra_days_before": a.pad_days, "gauge_peak_cfs": it["gauge_peak_cfs"]}
        if any_retuned:
            row["retuned_peak_cfs"] = rt["peak_cfs"] if rt else np.nan
        row["peak_rain_mm_h"] = it["rain_peak"]
        row["gauge_Mm3"] = it["gauge_vol"]
        if any_retuned:
            row["retuned_Mm3"] = rt["vol"] if rt else np.nan
        row["storm_rain_mm"] = it["rain_mm"]
        row["gauge_runoff_ratio_pct"] = it["runoff_pct"]
        if any_retuned:
            row["retuned_over_gauge_volume"] = rt["ratio"] if rt else np.nan
            row["volume_flag"] = flag_text.get(rt["flag"], "") if rt else ""
        row.update({"rain_above_4mm_h_mm": it["x4_mm"], "radar_gauge_3": it["r3"], "radar_gauge_2": it["r2"],
                    "dry_days_before": it["dry_days"], "rain_prior_30d_mm": it["prior_mm"], "prior_30d_hours_in_file": it["prior_hours"]})
        for k in extra_keys:
            r = by_key.get(k)
            t = tag_of(k)
            row[t + "_peak_cfs"] = r["peak_cfs"] if r else np.nan
            row[t + "_Mm3"] = r["vol"] if r else np.nan
            row[t + "_over_gauge_volume"] = r["ratio"] if r else np.nan
            row[t + "_volume_flag"] = flag_text.get(r["flag"], "") if r else ""
        out_rows.append(row)
    csv = out_dir / CSV_NAME.replace(".csv", suffix + ".csv")
    pd.DataFrame(out_rows).to_csv(csv, index=False)
    log("\nWrote:\n  {}\n  {}".format(fig_path, csv))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
