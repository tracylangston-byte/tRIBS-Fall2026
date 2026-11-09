#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
all_hydrographs_log1.py   (v1.1, 2026-11-06)
===========================================

ONE figure with the 12 storms side by side in time order (6 on top, 6 below): WY2017 first, then WY2020.
Only two lines are drawn: the USGS gauge (black) and the Re-tuned model (red, dashed).

EACH COLUMN HAS THREE PARTS
  top      hourly radar rain over the basin (mm/h, bars hang from the top)   <- basin mean, the same rain the model received
  middle   outlet flow in CFS on a LOG axis: gauge (black) and Re-tuned model (red dashed)
  bottom   a panel of 10 numbers for the storm, in this order (rows 1, 6 to 9 come from storm_table.csv; the rest are worked out here):
             1  Storm rain (mm)            basin rain over the storm days
             2  Runoff ratio (% of rain)   GAUGE volume / rain volume x 100. Rain volume = storm rain (mm) x 0.4246 M m3 per mm (the basin area).
                                           100% would mean every drop of rain left the basin as streamflow. Same definition as the table's gauge runoff ratio.
             3  Volume at gauge (M m3)     USGS gauge volume over the storm window (M m3 = million cubic meters)
             4  tRIBS volume (M m3)        Re-tuned model volume over the same window
             5  tRIBS / gauge volume       model volume / gauge volume (the usual direction: 1.00 = perfect, above 1 = model made too much,
                                           below 1 = too little). Shown as n/a when the gauge volume is below --min-gauge-vol (default 0.05 M m3),
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
    WY2020: the numbers in the WY2020 write-up (gauge and Re-tuned).
  It prints MATCH or MISMATCH for each. A MISMATCH means a wrong file, wrong window or wrong run was used: do not use the figure until that is explained.

WHAT IT READS (read-only; no model is run)
  storm_table.csv                                 (short headers such as "first_day" OR the readable ones such as "First rain day")
  storm_features_out/hourly_basin_wy2017.csv and hourly_basin_wy2020.csv      columns interval_start_mst, basin_mean_mm
  results/W5_d012_k10_f0p15_as300_ksx1p03_nocl_Outlet.qout     WY2017 Re-tuned (Time_hr 0 = Oct 1 2016 00:00)
  results/Sycamore_WY2020_Retuned_Outlet.qout                   WY2020 Re-tuned (Time_hr 0 = Sep 1 2019 00:00)
  ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv and ..._WY2020.csv    columns datetime, q_cms
WHAT IT WRITES (in all_hydro_figs/ unless you set --out-dir)
  fig_all_hydrographs_log1.png     the figure
  all_hydrographs_log1.csv         per storm: window, peaks (CFS), peak rain, and ALL the numbers from the bottom panel
                                   (the raw gauge and Re-tuned volumes are always kept; the ratio column is blank where the figure says n/a)

COMMANDS (from workspaces/Sycamore/lab/)
  python all_hydrographs_log1.py --dry-run                 # checks the files and lists the 12 windows; writes nothing
  python all_hydrographs_log1.py                           # the figure, 1 extra day before each storm
  python all_hydrographs_log1.py --pad-days 4              # 4 extra days before each storm (0 to 7)
  python all_hydrographs_log1.py --free-y                  # each panel gets its own flow and rain scale
  python all_hydrographs_log1.py --log-min-cfs 0.01        # draw smaller flows
  python all_hydrographs_log1.py --min-gauge-vol 0.1       # show the tRIBS / gauge ratio only when the gauge volume is at least 0.1 M m3
  python all_hydrographs_log1.py --flag-tol 0.10           # green flag only when the volume is within 10% (default 0.25 = 25%)
  python all_hydrographs_log1.py --storms "WY2017:Dec 21-25;WY2017:Dec 30-Jan 2;WY2020:Nov 19-22" --cols 3     # other storms (table names)
  Storm spec: YEAR:Event  or  YEAR:Event=Label   (several separated by ;). The year is WY2017 or WY2020. Event is the table's name.

WORDS
  hyetograph  rain bars against time.     log axis  each gridline is 10 times the one below.
  CFS         cubic feet per second.      window    the stretch of time a storm's volume is added up over.
"""
from __future__ import annotations

import argparse
import math
import sys
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
MODEL_COLOR = "#d62728"        # red, as asked

YEARS = {
    "WY2017": {"start": "2016-10-01 00:00", "stem": "W5_d012_k10_f0p15_as300_ksx1p03_nocl",
               "hourly": "hourly_basin_wy2017.csv", "obs": "sycamore_09510200_obs_flow_WY2017.csv"},
    "WY2020": {"start": "2019-09-01 00:00", "stem": "Sycamore_WY2020_Retuned",
               "hourly": "hourly_basin_wy2020.csv", "obs": "sycamore_09510200_obs_flow_WY2020.csv"},
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
    fig = plt.figure(figsize=(3.55 * cols + 1.1, 6.1 * blocks + 1.9))
    outer = fig.add_gridspec(blocks, 1, left=0.05, right=0.995, top=0.912, bottom=0.105, hspace=0.17)

    for b in range(blocks):
        inner = outer[b].subgridspec(4, cols, height_ratios=[0.8, 2.8, 0.55, 2.3], hspace=0.06, wspace=0.10)
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
            g, m = it["gauge_cfs"], it["model_cfs"]
            ax_q.plot(g.index, logsafe(g.values, a.log_min_cfs), color=GAUGE_COLOR, lw=2.2, label="USGS gauge 09510200", zorder=3)
            ax_q.plot(m.index, logsafe(m.values, a.log_min_cfs), color=MODEL_COLOR, lw=1.8, ls="--", label="Re-tuned model", zorder=4)
            ax_q.set_yscale("log")
            if a.free_y:
                top = max([x for x in (it["gauge_peak_cfs"], it["model_peak_cfs"]) if np.isfinite(x)] + [a.log_min_cfs * 10]) * 2.5
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
            lines = [("Storm rain (mm)", fmt(it["rain_mm"], 1)),
                     ("Runoff ratio (% of rain)", fmt_pct(it["runoff_pct"])),
                     ("Volume at gauge (M m³)", fmt(it["gauge_vol"], 3)),
                     ("tRIBS volume (M m³)", fmt(it["model_vol"], 3)),
                     ("tRIBS / gauge volume", fmt_ratio(it["vol_ratio"])),
                     ("Rain above 4 mm/h (mm)", fmt(it["x4_mm"], 2)),
                     ("Radar:gauge, 3 gauges", fmt(it["r3"], 2)),
                     ("Radar:gauge, 2 gauges", fmt(it["r2"], 2)),
                     ("Dry days before", fmt(it["dry_days"], 0)),
                     ("Rain, prior 30 days (mm)", fmt(it["prior_mm"], 1))]
            step = 0.1
            for j, (lab, val) in enumerate(lines):
                y = 0.95 - j * step
                ax_i.text(0.04, y, lab, transform=ax_i.transAxes, ha="left", va="center", fontsize=8, color=INK_2)
                kw = {}
                if j == 4 and it.get("vol_flag"):          # the tRIBS / gauge volume row gets a colored flag
                    sym, col = FLAG_STYLE[it["vol_flag"]]
                    val = sym + " " + val
                    kw = {"bbox": dict(boxstyle="round,pad=0.25", fc=col, ec="none")}
                ax_i.text(0.96, y, val, transform=ax_i.transAxes, ha="right", va="center", fontsize=8.5, color=INK,
                          fontweight="bold", **kw)
            # thin lines between the three groups: water volumes | rain intensity and radar | antecedent conditions
            for after_row in (4, 7):
                yy = 0.95 - (after_row + 0.5) * step
                ax_i.plot([0.02, 0.98], [yy, yy], transform=ax_i.transAxes, color=GRID, lw=0.8, clip_on=False)

    handles = [Line2D([], [], color=GAUGE_COLOR, lw=2.2, label="USGS gauge 09510200"),
               Line2D([], [], color=MODEL_COLOR, lw=1.8, ls="--", label="Re-tuned model"),
               Patch(color=RAIN_BLUE, label="hourly radar rain, basin mean"),
               Patch(color=PAD_GRAY, label="extra day(s) before the storm window ({:d})".format(a.pad_days)),
               Patch(color=FLAG_STYLE["within"][1], label="{} tRIBS / gauge volume within {:g}% of 1.00".format(FLAG_STYLE["within"][0], 100 * a.flag_tol)),
               Patch(color=FLAG_STYLE["high"][1], label="{} model volume more than {:g}% high".format(FLAG_STYLE["high"][0], 100 * a.flag_tol)),
               Patch(color=FLAG_STYLE["low"][1], label="{} model volume more than {:g}% low".format(FLAG_STYLE["low"][0], 100 * a.flag_tol))]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.040), ncol=4, frameon=False, fontsize=9.5)
    fig.suptitle("Lower Sycamore Wash: gauge and Re-tuned model, storm by storm (outlet flow in CFS, log axis)",
                 x=0.012, y=0.985, ha="left", fontsize=14, color=INK)
    sub = ("WY2017 then WY2020, in time order. {} flow and rain scales in every panel. "
           "Storm window: 00:00 on the first rain day to 00:00 four days after the last rain day.").format(
               "Separate" if a.free_y else "Shared")
    fig.text(0.012, 0.947, sub, fontsize=9.5, color=INK_2, ha="left", va="top")
    foot = ("Flow of zero or below {:g} CFS cannot be drawn on a log axis, so a line stops there (not missing data). 'Jan 19-24' is the table's Jan 19-21 storm "
            "(window runs through Jan 24). 2 gauges = without the Mt. Ord gauge.\n"
            "Runoff ratio = gauge volume / (storm rain x {:.4f} M m³ per mm). Volumes are over the storm window. tRIBS / gauge volume: 1 = perfect, above 1 = model made too much, "
            "below 1 = too little; n/a = gauge volume under {:g} M m³ (or not computed). Prior 30 days = basin radar rain before the first rain day.").format(
                a.log_min_cfs, BASIN_MM_TO_MM3, a.min_gauge_vol)
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
    ap = argparse.ArgumentParser(description="12 storm hydrographs in one figure: gauge and Re-tuned model, flow in CFS on a log axis (see the top of this file).")
    ap.add_argument("--root", default=str(HERE), help="the lab folder (default: the folder this script is in)")
    ap.add_argument("--table", default=None, help="storm_table.csv (default <root>/storm_table.csv)")
    ap.add_argument("--hourly-dir", default=None, help="folder with hourly_basin_wy2017.csv and ..._wy2020.csv (default <root>/storm_features_out)")
    ap.add_argument("--results", default=None, help="folder with the *_Outlet.qout files (default <root>/results)")
    ap.add_argument("--obs-dir", default=None, help="folder with the two USGS files (default <root>/../init_data/obs)")
    ap.add_argument("--out-dir", default=None, help="where the figure and CSV go (default <root>/" + FIG_DIR + ")")
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

    log("\nall_hydrographs_log1: reading files")
    tab = read_table(need(table, "storm table"))
    log("  storm table: {} rows".format(len(tab)))
    need_years = sorted(set(s[0] for s in storms))
    rain, gauge, model = {}, {}, {}
    for yr in need_years:
        info = YEARS[yr]
        rain[yr] = read_rain(need(hourly_dir / info["hourly"], "hourly rain " + yr))
        gauge[yr] = read_obs(need(obs_dir / info["obs"], "USGS gauge " + yr))
        model[yr] = read_qout(need(results / (info["stem"] + "_Outlet.qout"), "Re-tuned outlet " + yr), info["start"])
        log("  {}: rain {} h ({:%b %d %Y} to {:%b %d %Y}); gauge {} readings ({:%b %d %Y} to {:%b %d %Y}); Re-tuned {} rows".format(
            yr, len(rain[yr]), rain[yr].index.min(), rain[yr].index.max(), gauge[yr].notna().sum(), gauge[yr].index.min(), gauge[yr].index.max(),
            len(model[yr])))

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
              "tab_gauge": float(r["gauge_storm"]), "tab_model": float(r["model_storm"])}
        items.append(it)
    items.sort(key=lambda d: d["first"])

    log("\nWindows ({} storms, {} extra day(s) before each):".format(len(items), a.pad_days))
    for it in items:
        log("  {:10s} {:9s} window {:%Y-%m-%d %H:%M} to {:%Y-%m-%d %H:%M}   drawn from {:%Y-%m-%d}".format(
            it["label"], it["year"], it["lo"], it["hi"], it["lo_plot"]))
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
    log("\nVolumes over the window (M m3), scoring rule; peaks in CFS are description only")
    log("  {:10s} {:>9s} {:>10s}   {:>10s} {:>10s}   {:>9s} {:>11s}".format(
        "storm", "gauge", "Re-tuned", "gauge pk", "model pk", "runoff %", "tRIBS/gauge"))
    for it in items:
        yr = it["year"]
        g_cms, m_cms = gauge[yr], model[yr]
        lo, hi, lo_plot = it["lo"], it["hi"], it["lo_plot"]
        it["gauge_cfs"] = g_cms[(g_cms.index >= lo_plot) & (g_cms.index <= hi)] * CFS_PER_CMS
        it["model_cfs"] = m_cms[(m_cms.index >= lo_plot) & (m_cms.index <= hi)] * CFS_PER_CMS
        it["rain"] = rain[yr][(rain[yr].index >= lo_plot) & (rain[yr].index < hi)]
        it["gauge_peak_cfs"] = peak(g_cms, lo, hi) * CFS_PER_CMS
        it["model_peak_cfs"] = peak(m_cms, lo, hi) * CFS_PER_CMS
        it["rain_peak"] = float(it["rain"].max()) if len(it["rain"]) else float("nan")
        gv, gh = volume(g_cms, lo, hi)
        mv, mh = volume(m_cms, lo, hi)
        it["gauge_vol"], it["model_vol"] = gv, mv
        # runoff ratio (gauge): volume / (rain depth x basin area), in percent
        rain_vol = it["rain_mm"] * BASIN_MM_TO_MM3
        it["runoff_pct"] = 100.0 * gv / rain_vol if (np.isfinite(gv) and np.isfinite(rain_vol) and rain_vol > 0) else float("nan")
        # tRIBS / gauge volume: only when the gauge volume is big enough for a ratio to mean something
        it["vol_ratio"] = mv / gv if (np.isfinite(gv) and np.isfinite(mv) and gv >= a.min_gauge_vol) else float("nan")
        it["vol_flag"] = volume_flag(it["vol_ratio"], a.flag_tol)
        log("  {:10s} {:>9s} {:>10s}   {:>10s} {:>10s}   {:>9s} {:>11s}".format(
            it["label"], fmt(gv, 3), fmt(mv, 3), fmt(it["gauge_peak_cfs"], 1), fmt(it["model_peak_cfs"], 1),
            fmt_pct(it["runoff_pct"]), fmt_ratio(it["vol_ratio"])))
        exp = None
        if np.isfinite(it["tab_gauge"]) and np.isfinite(it["tab_model"]):
            exp = (it["tab_gauge"], it["tab_model"], "the table")
        elif yr == "WY2020" and it["event"] in EXPECTED_WY2020:
            exp = EXPECTED_WY2020[it["event"]] + ("the WY2020 write-up",)
        if exp is not None:
            it["checked"] = True
            if not (np.isfinite(gv) and abs(gv - exp[0]) <= MATCH_TOL):
                bad.append("{} gauge: here {} vs {} {:.3f}".format(it["label"], fmt(gv, 3), exp[2], exp[0]))
            if not (np.isfinite(mv) and abs(mv - exp[1]) <= MATCH_TOL):
                bad.append("{} Re-tuned: here {} vs {} {:.3f}".format(it["label"], fmt(mv, 3), exp[2], exp[1]))
        else:
            it["checked"] = False
        if gh < 0.97 * int((hi - lo) / pd.Timedelta(hours=1)):
            log("    note: the gauge has data for {} of {} hours in this window; its volume counts only those hours.".format(
                gh, int((hi - lo) / pd.Timedelta(hours=1))))
    n_checked = sum(1 for it in items if it["checked"])
    if n_checked:
        if bad:
            log("\nCHECK: MISMATCH ({} value(s) differ by more than {:g} M m3):".format(len(bad), MATCH_TOL))
            for x in bad:
                log("  " + x)
            log("  The figure is still drawn, but do not use it until this is explained.")
        else:
            log("\nCHECK: MATCH. Gauge and Re-tuned volumes for {} of {} storms equal the table / write-up values (within {:g} M m3).".format(
                n_checked, len(items), MATCH_TOL))
    if n_checked < len(items):
        log("  ({} storm(s) have nothing to check against and are not checked.)".format(len(items) - n_checked))

    # shared scales
    tops = [x for it in items for x in (it["gauge_peak_cfs"], it["model_peak_cfs"]) if np.isfinite(x)]
    flow_top = (max(tops) if tops else 100.0) * 2.5
    rmax = [it["rain_peak"] for it in items if np.isfinite(it["rain_peak"])]
    rain_top = max(1.0, float(math.ceil(max(rmax)))) if rmax else 1.0

    fig_path = make_figure(out_dir / FIG_NAME, items, a, flow_top, rain_top)
    rows = [{"storm": it["label"], "table_event": it["event"], "year": it["year"], "window_start": it["lo"], "window_end": it["hi"],
             "extra_days_before": a.pad_days, "gauge_peak_cfs": it["gauge_peak_cfs"], "retuned_peak_cfs": it["model_peak_cfs"],
             "peak_rain_mm_h": it["rain_peak"], "gauge_Mm3": it["gauge_vol"], "retuned_Mm3": it["model_vol"],
             "storm_rain_mm": it["rain_mm"], "gauge_runoff_ratio_pct": it["runoff_pct"],
             "retuned_over_gauge_volume": it["vol_ratio"],
             "volume_flag": {"within": "within {:g}%".format(100 * a.flag_tol), "high": "model too high", "low": "model too low"}.get(it["vol_flag"], ""),
             "rain_above_4mm_h_mm": it["x4_mm"], "radar_gauge_3": it["r3"], "radar_gauge_2": it["r2"],
             "dry_days_before": it["dry_days"], "rain_prior_30d_mm": it["prior_mm"], "prior_30d_hours_in_file": it["prior_hours"]}
            for it in items]
    csv = out_dir / CSV_NAME
    pd.DataFrame(rows).to_csv(csv, index=False)
    log("\nWrote:\n  {}\n  {}".format(fig_path, csv))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
