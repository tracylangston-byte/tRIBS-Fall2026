#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_wy2020_season_sycamore.py   (v1.0, 2026-11-09)
===================================================

The WHOLE WY2020 run, Sep 1 2019 to May 31 2020, in one figure: the USGS gauge against up to FOUR model runs.
A NEW file. all_hydrographs_log1_v1_3.py and plot_wy2020_hydrographs_sycamore.py are not changed.

WHAT THE FIGURE HAS (one PNG)
  Rows 1 to 3   Sep to Nov, Dec to Feb, Mar to May. Each row has
                  top     hourly radar rain over the basin (mm/h, bars hang from the top; the same rain the model received)
                  bottom  outlet flow in CFS on a LOG axis: the gauge (thick black) and each model run
                All three rows share ONE rain scale and ONE flow scale, so the rows can be compared by eye (--free-y gives each row its own flow scale).
                Gray bands are the seven scored storm windows, named at the top, plus the Sep 23 storm (not scored).
  Row 4 left    CUMULATIVE VOLUME: all the water that has left the outlet since Sep 1, in M m3, for the gauge and each run.
                The legend gives each run's total for the season and its ratio to the gauge's total.
  Row 4 right   RUN MINUS GAUGE: the same cumulative volumes with the gauge's taken away. A line going UP means the run is sending out more
                water than the gauge saw in that stretch; DOWN means less. Zero means the run has matched the gauge so far.

THE RUNS (--runs, up to four, in the order you give them; the gauge is always drawn)
  Names are the same as in all_hydrographs_log1_v1_3.py:
    retuned   original   kchan05   kchan02            short names
    S5_Ks1p0_Bed2p0                                    any run in results/ (the part before _Outlet.qout; the prefix Sycamore_WY2020_ may be left off,
                                                       and a run saved with no such prefix, such as a W5 run, is found too)
  --labels gives the names shown on the figure, in the same order, separated by commas (so no commas inside a name).
  If you leave --labels off, Re-tuned is called Re-tuned and the screens get names such as "S5 bedrock +2 m", "S6 Ks x0.5", "S7 As 500", "S8 f x0.02".
  Re-tuned is always red and dashed; the other runs are blue, orange, green and purple, with different line styles as well (so no run depends on color alone).

HOW THE NUMBERS ARE WORKED OUT (the scoring rule)
  Flow = hourly mean of the readings; volume = hourly mean x 3,600 s. Hours with no reading add nothing (the gauge has 2 such hours).
  A storm window is 00:00 on the first rain day to 00:00 four days after the last rain day (the table's rule). Model rows are placed at Sep 1 00:00 + Time_hr and averaged by clock hour,
  exactly as the scoring and edge-check scripts do, so the volumes here can be compared with theirs.
  The figure plots hourly means, so a short peak can look a little lower than in the 15-minute data.

WHAT IS PRINTED (and saved in a CSV)
  1  the volume in each of the seven scored storms (gauge, each run, run / gauge)
  2  the volume in each month, Sep to May, and the season total
  3  the volume OUTSIDE the seven scored windows (the same number the edge check calls "outside the windows")
  A CHECK: the gauge and Re-tuned volumes (when Re-tuned is drawn) are compared with the WY2020 write-up for the seven storms, and with the
  monthly numbers of the Nov 7 edge check for every month. It prints MATCH or MISMATCH. A MISMATCH means a wrong file, a shifted clock or a
  wrong run; do not use the figure until it is explained. Other runs have no reference numbers, so they are not checked.

WHAT IT READS (read-only; no model is run)
  storm_table.csv                                    the seven storms' first and last rain days (short or readable headers both work)
  storm_features_out/hourly_basin_wy2020.csv         columns interval_start_mst, basin_mean_mm
  results/Sycamore_WY2020_<run>_Outlet.qout          first two columns Time_hr, flow (m3/s); Time_hr 0 = Sep 1 2019 00:00
  ../init_data/obs/sycamore_09510200_obs_flow_Sep2019_to_May2020.csv    columns datetime, q_cms
        This is the combined gauge file made by add_sep2019_gauge_sycamore.py. The old ..._WY2020.csv starts on Oct 1 and has no September,
        so the script stops if it is given a gauge file with no September readings.
WHAT IT WRITES (in season_figs/ unless you set --out-dir; an existing file with the same name is replaced and the script says so)
  fig_season_<runs>.png        the figure
  season_volumes_<runs>.csv    every number printed in 1 to 3

COMMANDS (from workspaces/Sycamore/lab/)
  python plot_wy2020_season_sycamore.py --dry-run                         # checks the files, lists the runs and windows; writes nothing
  python plot_wy2020_season_sycamore.py                                   # Re-tuned only
  python plot_wy2020_season_sycamore.py --runs retuned,S5_Ks1p0_Bed2p0,S7_Ks1p0_Bed2p0_As500,S8_Ks1p0_Bed2p0_f0p02 --labels "Re-tuned,S5 bedrock +2 m,S7 As 500,S8 f x0.02"
  python plot_wy2020_season_sycamore.py --runs S5_Ks1p0_Bed2p0,S6_Ks0p5_Bed2p0,S7_Ks1p0_Bed2p0_As500,S8_Ks1p0_Bed2p0_f0p02
  python plot_wy2020_season_sycamore.py --free-y                          # each row gets its own flow scale
  python plot_wy2020_season_sycamore.py --extra-windows ""                # no extra gray band for the Sep 23 storm

WORDS
  hyetograph  a bar chart of rain against time.   hydrograph  flow at one place against time.
  log axis    each gridline is 10 times the one below it; it shows small flows next to big ones. Flow of zero cannot be drawn on a log axis,
              so a line stops there (--log-min-cfs, default 0.1 CFS). That is not missing data.
  cumulative  a running total: each point is everything up to that time added up.
UNITS   flow in CFS (cubic feet per second) = m3/s x 35.3147. Volumes in M m3 = million cubic meters.

VERSION HISTORY
  v1.0 (2026-11-09)  first version. Tested on a stand-in lab (the same file names and formats as the real ones), not on your own files.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUN_START = pd.Timestamp("2019-09-01 00:00")          # Time_hr 0 in the WY2020 outlet files
SEASON_END = pd.Timestamp("2020-06-01 00:00")         # the run ends May 31 2020 23:00 (6,575 h)
AFTER_DAYS = 3                                        # window ends 3 days after the day AFTER the last rain day
CFS_PER_CMS = 35.3147
MAX_RUNS = 4
FIG_DIR = "season_figs"
GAUGE_FILE = "sycamore_09510200_obs_flow_Sep2019_to_May2020.csv"
MATCH_TOL = 0.002                                     # the write-up and edge-check numbers are given to 3 or 4 decimals (M m3)
MIN_RATIO_GAUGE = 0.05                                # run / gauge is shown only if the gauge volume is at least this many M m3

INK, INK_2, GRID = "#0b0b0b", "#52514e", "#d9d8d3"
RAIN_BLUE, BAND_GRAY, BAND_SEP = "#8db4e2", "#ebeae6", "#f1e6d6"
GAUGE_COLOR = "#0b0b0b"
RETUNED = "retuned"
RETUNED_STYLE = ("#d62728", "--", 1.8)                # red, dashed (as in the storm figure)
OTHER_STYLES = [("#2a78d6", "-", 1.7), ("#e8871e", "-.", 1.7), ("#3aa66b", ":", 2.0), ("#7b4fc4", "-", 1.5)]

WY2020_PREFIX = "Sycamore_WY2020_"
WY2020_SHORT = {
    "retuned": (WY2020_PREFIX + "Retuned", "Re-tuned"),
    "original": (WY2020_PREFIX + "Original", "Original"),
    "kchan05": (WY2020_PREFIX + "Original_kchan0p5", "Orig. Kchan 0.5"),
    "kchan02": (WY2020_PREFIX + "Original_kchan0p2", "Orig. Kchan 0.2"),
}
SCREEN_LABELS = {
    "S1": "S1 Ks x1.0", "S2": "S2 Ks x0.5", "S3": "S3 bedrock +0.5 m", "S4": "S4 start GW 0.5",
    "S5": "S5 bedrock +2 m", "S6": "S6 Ks x0.5", "S7": "S7 As 500", "S8": "S8 f x0.02",
}

# The seven scored storms: name -> (first rain day, last rain day), used only if storm_table.csv has no usable row for it.
SCORED = ["Nov 19-22", "Nov 28-30", "Dec 23-28", "Feb 10-11", "Feb 22-23", "Mar 10-13", "Mar 18-19"]
STORM_DAYS_FALLBACK = {
    "Nov 19-22": ("2019-11-19", "2019-11-22"), "Nov 28-30": ("2019-11-28", "2019-11-30"), "Dec 23-28": ("2019-12-23", "2019-12-28"),
    "Feb 10-11": ("2020-02-10", "2020-02-11"), "Feb 22-23": ("2020-02-22", "2020-02-23"), "Mar 10-13": ("2020-03-10", "2020-03-13"),
    "Mar 18-19": ("2020-03-18", "2020-03-19"),
}
# The Sep 23 storm is not in storm_table.csv (the table starts Oct 1). Same window as the storm figure's September panel.
DEFAULT_EXTRA = "Sep 23 (not scored)=2019-09-23:2019-09-23"

# Reference numbers for the CHECK (M m3), scoring rule. Gauge and Re-tuned only.
EXPECTED_STORMS = {   # WY2020 write-up: (gauge, Re-tuned)
    "Nov 19-22": (0.015, 2.800), "Nov 28-30": (1.473, 3.235), "Dec 23-28": (4.069, 0.218), "Feb 10-11": (0.386, 1.767),
    "Feb 22-23": (4.131, 4.839), "Mar 10-13": (5.626, 0.875), "Mar 18-19": (3.317, 1.031),
}
EXPECTED_MONTHS = {   # Nov 7 edge check, monthly outlet table, all hours: (gauge, Re-tuned); months by the start of each hour
    "2019-09": (3.4952, 29.4115), "2019-10": (0.0497, 0.1044), "2019-11": (1.0745, 5.7818), "2019-12": (5.9273, 0.8034),
    "2020-01": (2.4398, 0.8061), "2020-02": (5.6316, 6.8090), "2020-03": (12.0199, 2.2584), "2020-04": (1.7729, 0.1661),
    "2020-05": (0.3216, 0.0249),
}
MONTHS = list(EXPECTED_MONTHS)

COLS = {
    "year": ["year", "Water year"],
    "event": ["event", "Event"],
    "first_day": ["first_day", "First rain day"],
    "last_day": ["last_day", "Last rain day"],
}


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
    """A run key made safe for a file name; the screens (S5_Ks1p0_Bed2p0) are shortened to S5."""
    m = re.match(r"^(S\d+)_", key)
    if m:
        return m.group(1)
    return re.sub(r"[^A-Za-z0-9_]+", "_", key).strip("_")


# ------------------------------------------------------------------------------------------------ the runs
def run_label(key: str) -> str:
    m = re.match(r"^(S\d+)(_|$)", key)
    if m and m.group(1) in SCREEN_LABELS:
        return SCREEN_LABELS[m.group(1)]
    return key


def parse_runs(spec: str, labels):
    names = [s.strip() for s in (spec or RETUNED).split(",") if s.strip()]
    if not names:
        die("--runs needs at least one run, for example  --runs retuned,S5_Ks1p0_Bed2p0")
    if len(names) > MAX_RUNS:
        die("--runs takes at most {} runs (you gave {}).".format(MAX_RUNS, len(names)))
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
            label = run_label(key)
            for sk, (sstem, slabel) in WY2020_SHORT.items():       # a full stem that is one of the short names counts as that name
                if stem == sstem:
                    key, label = sk, slabel
        if key in seen:
            die("--runs lists '{}' twice.".format(nm))
        seen.add(key)
        runs.append({"key": key, "stem": stem, "label": label, "tag": tag_of(key)})
    if len({r["tag"] for r in runs}) != len(runs):
        die("two of the runs would get the same short name in the file names: {}".format(", ".join(r["tag"] for r in runs)))
    if labels:
        labs = [s.strip() for s in labels.split(",")]
        if len(labs) != len(runs) or any(not x for x in labs):
            die("--labels needs one name for each run, in the same order, separated by commas.\n"
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
    if stem.startswith(WY2020_PREFIX):                      # a run saved without the prefix (for example a W5 run) is found too
        bare = results / (stem[len(WY2020_PREFIX):] + "_Outlet.qout")
        if bare.is_file():
            return bare
    have = sorted(f.name[:-len("_Outlet.qout")] for f in results.glob("*_Outlet.qout")) if results.is_dir() else []
    msg = "cannot find the outlet file for {}:\n    {}".format(what, p)
    if have:
        msg += "\n  Outlet files found in {}:\n    ".format(results) + "\n    ".join(have[:25])
        msg += "\n  Give the run name without _Outlet.qout (the prefix {} may be left off).".format(WY2020_PREFIX)
    else:
        msg += "\n  Run this from workspaces/Sycamore/lab/, or point to the folder with --results."
    die(msg)
    return p


# ------------------------------------------------------------------------------------------------ reading
def read_qout(path: Path) -> pd.Series:
    raw = pd.read_csv(path).iloc[:, :2]
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    idx = RUN_START + pd.to_timedelta(raw["Time_hr"].to_numpy(float), unit="h")
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


def read_scored_storms(path: Path):
    """The seven scored storms as a list of dicts (name, first, last). Dates come from storm_table.csv; a storm with no usable row uses the built-in days."""
    raw = pd.read_csv(path)
    pick = {}
    for want, options in COLS.items():
        hit = [c for c in options if c in raw.columns]
        if not hit:
            die("{} has no column for '{}' (looked for {}).\n  Columns found: {}".format(path, want, options, list(raw.columns)))
        pick[want] = hit[0]
    t = pd.DataFrame({k: raw[v] for k, v in pick.items()})
    t = t[t["year"].astype(str).str.contains("2020")].copy()
    t["event"] = t["event"].astype(str).str.replace("–", "-").str.replace("—", "-").str.strip()
    t["first_day"] = pd.to_datetime(t["first_day"], errors="coerce")
    t["last_day"] = pd.to_datetime(t["last_day"], errors="coerce")
    out, notes = [], []
    for name in SCORED:
        row = t[t["event"] == name]
        ok = (not row.empty and pd.notna(row.iloc[0]["first_day"]) and pd.notna(row.iloc[0]["last_day"])
              and pd.Timestamp("2019-09-01") <= row.iloc[0]["first_day"] <= row.iloc[0]["last_day"] <= pd.Timestamp("2020-06-01"))
        if ok:
            a, b = pd.Timestamp(row.iloc[0]["first_day"]), pd.Timestamp(row.iloc[0]["last_day"])
        else:
            a, b = (pd.Timestamp(x) for x in STORM_DAYS_FALLBACK[name])
            notes.append(name)
        out.append({"name": name, "first": a, "last": b, "scored": True})
    return out, notes


def parse_extra(spec: str):
    """'NAME=YYYY-MM-DD:YYYY-MM-DD;NAME2=...' -> list of window dicts (first rain day : last rain day)."""
    out = []
    for part in [p for p in (spec or "").split(";") if p.strip()]:
        m = re.match(r"^\s*(.+?)\s*=\s*(\d{4}-\d{2}-\d{2})\s*:\s*(\d{4}-\d{2}-\d{2})\s*$", part)
        if not m:
            die("--extra-windows: could not read '{}'. Use  NAME=FIRST_DAY:LAST_DAY  with dates like 2019-09-23, separated by ';'.".format(part))
        a, b = pd.Timestamp(m.group(2)), pd.Timestamp(m.group(3))
        if b < a:
            die("--extra-windows: the last day is before the first day in '{}'.".format(part))
        out.append({"name": m.group(1), "first": a, "last": b, "scored": False})
    return out


def add_window(w: dict) -> dict:
    w["lo"] = w["first"]
    w["hi"] = w["last"] + pd.Timedelta(days=1 + AFTER_DAYS)
    return w


# ------------------------------------------------------------------------------------------------ numbers
def hourly_mean(q: pd.Series) -> pd.Series:
    """Hourly mean over the whole season; an hour with no reading is NaN."""
    h = q.resample("60min").mean()
    return h[(h.index >= RUN_START) & (h.index < SEASON_END)]


def vol_per_hour(h: pd.Series) -> pd.Series:
    """M m3 in each hour (hourly mean flow x 3,600 s); an hour with no reading adds nothing."""
    return (h * 3600.0 / 1e6).fillna(0.0)


def window_volume(v: pd.Series, lo, hi) -> float:
    return float(v[(v.index >= lo) & (v.index < hi)].sum())


def month_volumes(v: pd.Series) -> dict:
    g = v.groupby(v.index.strftime("%Y-%m")).sum()
    return {m: float(g.get(m, np.nan)) for m in MONTHS}


def outside_volume(v: pd.Series, windows) -> float:
    keep = pd.Series(True, index=v.index)
    for w in windows:
        keep &= ~((v.index >= w["lo"]) & (v.index < w["hi"]))
    return float(v[keep].sum())


def fmt_v(x: float) -> str:
    if not np.isfinite(x):
        return "n/a"
    return "{:.3f}".format(x) if abs(x) < 0.1 else "{:.2f}".format(x)


def ratio(x: float, g: float) -> float:
    return x / g if (np.isfinite(x) and np.isfinite(g) and g >= MIN_RATIO_GAUGE) else float("nan")


def fmt_r(r: float) -> str:
    return "n/a" if not np.isfinite(r) else "{:.2f}".format(r)


def run_checks(vols_h: dict, storms: list, run_keys: list, tol: float = MATCH_TOL):
    """Compare gauge (always) and Re-tuned (if drawn) with the write-up (seven storms) and the edge check (months).
    vols_h = {'gauge': hourly volume series, key: hourly volume series}. Returns (list of mismatch strings, number of values checked)."""
    bad, n = [], 0
    for key in ["gauge"] + ([RETUNED] if RETUNED in run_keys else []):
        col = 0 if key == "gauge" else 1
        for w in storms:
            exp = EXPECTED_STORMS[w["name"]][col]
            got = window_volume(vols_h[key], w["lo"], w["hi"])
            n += 1
            if not (np.isfinite(got) and abs(got - exp) <= tol):
                bad.append("storm {} {}: here {:.4f} vs write-up {:.3f}".format(w["name"], key, got, exp))
        mv = month_volumes(vols_h[key])
        for m in MONTHS:
            exp = EXPECTED_MONTHS[m][col]
            n += 1
            if not (np.isfinite(mv[m]) and abs(mv[m] - exp) <= tol):
                bad.append("month {} {}: here {:.4f} vs edge check {:.4f}".format(m, key, mv[m], exp))
    return bad, n


# ------------------------------------------------------------------------------------------------ drawing
def tidy(ax) -> None:
    ax.grid(alpha=0.35, color=GRID, which="major", lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=8.5)


def logsafe(y, floor: float) -> np.ndarray:
    y = np.asarray(y, float).copy()
    y[y <= floor] = np.nan
    return y


def make_figure(out: Path, runs, gauge_h, run_h, cum, rain, storms, extras, args, check_msg: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    blocks = [(pd.Timestamp("2019-09-01"), pd.Timestamp("2019-12-01"), "Sep to Nov 2019"),
              (pd.Timestamp("2019-12-01"), pd.Timestamp("2020-03-01"), "Dec 2019 to Feb 2020"),
              (pd.Timestamp("2020-03-01"), pd.Timestamp("2020-06-01"), "Mar to May 2020")]
    windows = storms + extras

    fig = plt.figure(figsize=(15.5, 17.0))
    outer = fig.add_gridspec(4, 1, height_ratios=[2.9, 2.9, 2.9, 3.0], hspace=0.30, left=0.065, right=0.965, top=0.925, bottom=0.125)

    flow_all = [gauge_h] + [run_h[r["key"]] for r in runs]
    top_cfs = max(np.nanmax(logsafe(h.values * CFS_PER_CMS, args.log_min_cfs)) for h in flow_all if np.isfinite(h.values).any())
    ymax = 10 ** np.ceil(np.log10(top_cfs * 1.8))
    rain_max = float(np.nanmax(rain.values)) if len(rain) else 1.0
    rain_lim = max(6.0, np.ceil(rain_max / 6.0) * 6.0)

    for i, (a, b, title) in enumerate(blocks):
        inner = outer[i].subgridspec(2, 1, height_ratios=[0.45, 2.45], hspace=0.05)
        ax_r = fig.add_subplot(inner[0])
        ax_q = fig.add_subplot(inner[1], sharex=ax_r)
        r = rain[(rain.index >= a) & (rain.index < b)]
        ax_r.bar(r.index, r.values, width=pd.Timedelta(hours=1), align="edge", color=RAIN_BLUE, zorder=2)
        ax_r.set_ylim(rain_lim, 0)
        ax_r.set_yticks([0, rain_lim / 2.0])                      # no tick at the bottom edge: its label would sit on the flow axis below
        ax_r.set_ylabel("Rain\n(mm/h)", color=INK_2, fontsize=8.5)
        plt.setp(ax_r.get_xticklabels(), visible=False)
        ax_r.set_title(title, loc="left", fontsize=11.5, color=INK, pad=4)
        if not args.no_shade:
            for w in windows:
                if w["hi"] > a and w["lo"] < b:
                    col = BAND_GRAY if w["scored"] else BAND_SEP
                    ax_r.axvspan(w["lo"], w["hi"], color=col, lw=0, zorder=0)
                    ax_q.axvspan(w["lo"], w["hi"], color=col, lw=0, zorder=0)
        for run in runs:
            h = run_h[run["key"]]
            h = h[(h.index >= a) & (h.index <= b)]
            ax_q.plot(h.index, logsafe(h.values * CFS_PER_CMS, args.log_min_cfs), color=run["color"], ls=run["ls"], lw=run["lw"], zorder=3)
        g = gauge_h[(gauge_h.index >= a) & (gauge_h.index <= b)]
        ax_q.plot(g.index, logsafe(g.values * CFS_PER_CMS, args.log_min_cfs), color=GAUGE_COLOR, lw=2.2, zorder=5)
        ax_q.set_yscale("log")
        if args.free_y:
            tops = [np.nanmax(logsafe(h[(h.index >= a) & (h.index <= b)].values * CFS_PER_CMS, args.log_min_cfs))
                    for h in flow_all if np.isfinite(h[(h.index >= a) & (h.index <= b)].values).any()]
            ax_q.set_ylim(args.log_min_cfs, 10 ** np.ceil(np.log10(max(tops) * 1.8)) if tops else ymax)
        else:
            ax_q.set_ylim(args.log_min_cfs, ymax)
        ax_q.set_xlim(a, b)
        ax_q.set_ylabel("Outlet flow (CFS), log", color=INK_2, fontsize=9)
        ax_q.xaxis.set_major_locator(mdates.MonthLocator())
        ax_q.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
        ax_q.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=[8, 15, 22]))
        ax_q.grid(which="minor", axis="x", alpha=0.25, color=GRID, lw=0.5)
        for ax in (ax_r, ax_q):
            tidy(ax)
        if not args.no_shade:
            for w in windows:
                if w["hi"] > a and w["lo"] < b:
                    x = max(w["lo"], a) + pd.Timedelta(hours=8)
                    ax_q.text(x, 0.985, w["name"], transform=ax_q.get_xaxis_transform(), ha="left", va="top", fontsize=8,
                              color=INK_2, zorder=6)

    # row 4: cumulative volume (left) and run minus gauge (right)
    inner4 = outer[3].subgridspec(1, 2, wspace=0.16)
    ax_c = fig.add_subplot(inner4[0])
    ax_d = fig.add_subplot(inner4[1])
    for ax in (ax_c, ax_d):
        if not args.no_shade:
            for w in windows:
                ax.axvspan(w["lo"], w["hi"], color=(BAND_GRAY if w["scored"] else BAND_SEP), lw=0, zorder=0)
        ax.set_xlim(RUN_START, SEASON_END)
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        tidy(ax)
    gtot = float(cum["gauge"].iloc[-1])
    ax_c.plot(cum["gauge"].index, cum["gauge"].values, color=GAUGE_COLOR, lw=2.2, zorder=5)
    for run in runs:
        c = cum[run["key"]]
        tot = float(c.iloc[-1])
        ax_c.plot(c.index, c.values, color=run["color"], ls=run["ls"], lw=run["lw"], zorder=3)
        run["total_text"] = "{} M m³ ({}× gauge)".format(fmt_v(tot), fmt_r(ratio(tot, gtot)))
        d = c - cum["gauge"]
        ax_d.plot(d.index, d.values, color=run["color"], ls=run["ls"], lw=run["lw"], zorder=3)
    ax_d.axhline(0, color=GAUGE_COLOR, lw=1.4, zorder=4)
    ax_c.set_ylim(bottom=0)
    ax_c.set_ylabel("Water out of the outlet since Sep 1 (M m³)", color=INK_2, fontsize=9)
    ax_d.set_ylabel("Run minus gauge, cumulative (M m³)", color=INK_2, fontsize=9)
    ax_c.set_title("Cumulative volume, whole season (season totals in the key below)", loc="left", fontsize=11.5, color=INK, pad=4)
    ax_d.set_title("Run minus gauge (up = run has sent out more water)", loc="left", fontsize=11.5, color=INK, pad=4)

    names = ", ".join(r["label"] for r in runs)
    fig.suptitle("Lower Sycamore Wash, Sep 1 2019 to May 31 2020: gauge and model ({}), outlet flow in CFS (log axis)".format(names),
                 x=0.012, y=0.985, ha="left", fontsize=14, color=INK)
    fig.text(0.012, 0.953, "Hourly mean flow. One flow scale and one rain scale (0 to {:g} mm/h, bars hang from the top) in rows 1 to 3{}. Gray bands: the seven scored storm windows (00:00 on the first rain day "
             "to 00:00 four days after the last); tan band: the Sep 23 storm, which is not scored.".format(rain_lim, " EXCEPT that each row has its own flow scale (--free-y)" if args.free_y else ""),
             fontsize=8.8, color=INK_2, ha="left", va="top")
    handles = [Line2D([], [], color=GAUGE_COLOR, lw=2.2, label="USGS gauge 09510200: season total {} M m³".format(fmt_v(gtot)))]
    handles += [Line2D([], [], color=r["color"], ls=r["ls"], lw=r["lw"], label="{} model: season total {}".format(r["label"], r["total_text"])) for r in runs]
    handles += [Line2D([], [], color=RAIN_BLUE, lw=6, label="hourly radar rain, basin mean")]
    if not args.no_shade:
        handles += [Patch(color=BAND_GRAY, label="scored storm window")]
        if extras:
            handles += [Patch(color=BAND_SEP, label="Sep 23 storm (not scored)")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.040), ncol=3, frameon=False, fontsize=9.5)
    foot = ("Volume = hourly mean flow x 3,600 s; hours with no reading add nothing. Flow of zero or below {:g} CFS cannot be drawn on a log axis, so a line stops there (not missing data). "
            "{}").format(args.log_min_cfs, check_msg)
    fig.text(0.012, 0.006, foot, fontsize=7.6, color=INK_2, ha="left", va="bottom", wrap=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="The whole WY2020 run (Sep 1 to May 31): gauge against up to four model runs (see the top of this file).")
    ap.add_argument("--root", default=str(HERE), help="the lab folder (default: the folder this script is in)")
    ap.add_argument("--table", default=None, help="storm_table.csv (default <root>/storm_table.csv)")
    ap.add_argument("--hourly", default=None, help="hourly_basin_wy2020.csv (default <root>/storm_features_out/hourly_basin_wy2020.csv)")
    ap.add_argument("--results", default=None, help="folder with the *_Outlet.qout files (default <root>/results)")
    ap.add_argument("--obs", default=None, help="the combined gauge file (default <root>/../init_data/obs/" + GAUGE_FILE + ")")
    ap.add_argument("--out-dir", default=None, help="where the figure and CSV go (default <root>/" + FIG_DIR + ")")
    ap.add_argument("--runs", default=RETUNED, help="up to four model runs, comma separated: retuned, original, kchan05, kchan02, or a run name such as S5_Ks1p0_Bed2p0 (default: retuned)")
    ap.add_argument("--labels", default=None, help="names to show for those runs, same order, comma separated, short (for example \"Re-tuned,S5 bedrock +2 m\")")
    ap.add_argument("--free-y", action="store_true", help="each of the three rows gets its own flow scale (default: one shared scale)")
    ap.add_argument("--log-min-cfs", type=float, default=0.1, help="lowest flow on the log axes, CFS (default 0.1)")
    ap.add_argument("--extra-windows", default=DEFAULT_EXTRA,
                    help="extra gray bands that are not scored, 'NAME=FIRST_DAY:LAST_DAY;NAME2=...' (default: the Sep 23 storm; \"\" for none)")
    ap.add_argument("--no-shade", action="store_true", help="draw no storm bands")
    ap.add_argument("--skip-checks", action="store_true", help="do not compare the gauge and Re-tuned volumes with the write-up and the edge check")
    ap.add_argument("--dry-run", action="store_true", help="check the files and list the runs and windows; write nothing")
    ap.add_argument("--dpi", type=int, default=150)
    a = ap.parse_args(argv)

    root = Path(a.root)
    table = Path(a.table) if a.table else root / "storm_table.csv"
    hourly = Path(a.hourly) if a.hourly else root / "storm_features_out" / "hourly_basin_wy2020.csv"
    results = Path(a.results) if a.results else root / "results"
    obs_path = Path(a.obs) if a.obs else root.parent / "init_data" / "obs" / GAUGE_FILE
    out_dir = Path(a.out_dir) if a.out_dir else root / FIG_DIR

    runs = parse_runs(a.runs, a.labels)
    log("\nWY2020 whole-season figure: reading files")
    storms_raw, notes = read_scored_storms(need(table, "storm table"))
    storms = [add_window(w) for w in storms_raw]
    extras = [add_window(w) for w in parse_extra(a.extra_windows)]
    rain = read_rain(need(hourly, "hourly basin rain"))
    if not obs_path.is_file() and obs_path.name == GAUGE_FILE:
        die("cannot find the combined gauge file:\n    {}\n  It is made by add_sep2019_gauge_sycamore.py (it holds September 2019 to May 2020). "
            "The old ..._WY2020.csv starts on Oct 1 and cannot be used here.".format(obs_path))
    gauge = read_obs(need(obs_path, "gauge"))
    if gauge.dropna().empty or gauge.dropna().index.min() > pd.Timestamp("2019-09-05"):
        die("the gauge file {} has no readings in early September 2019 (it starts {}). This figure needs the combined file that runs from Sep 1 "
            "(add_sep2019_gauge_sycamore.py makes it).".format(obs_path, gauge.dropna().index.min() if not gauge.dropna().empty else "never"))
    log("  rain: {} hours, {:%b %d %Y} to {:%b %d %Y}".format(len(rain), rain.index.min(), rain.index.max()))
    log("  gauge: {} readings, {:%b %d %Y} to {:%b %d %Y}   ({})".format(gauge.notna().sum(), gauge.index.min(), gauge.index.max(), obs_path.name))
    if notes:
        log("  NOTE: storm_table.csv had no usable WY2020 row for {}; the built-in storm days were used for them.".format(", ".join(notes)))
    series = {}
    for r in runs:
        p = need_run(results, r["stem"], "{} ({})".format(r["label"], r["key"]))
        q = read_qout(p)
        series[r["key"]] = q
        log("  {}: {} rows, {:%b %d %Y} to {:%b %d %Y %H:%M}   ({})".format(r["label"], len(q), q.index.min(), q.index.max(), p.name))
        if q.index.max() < pd.Timestamp("2020-05-31 22:00"):
            log("  WARNING: {} stops at {:%b %d %Y %H:%M}; the run should reach May 31 23:00. Its volumes are for a shorter run.".format(r["label"], q.index.max()))

    log("\nRuns drawn: " + ", ".join("{} [{}]".format(r["label"], r["key"]) for r in runs))
    log("Windows (scored storms = 00:00 on the first rain day to 00:00 four days after the last rain day):")
    for w in storms + extras:
        log("  {:22s} {:%Y-%m-%d %H:%M} to {:%Y-%m-%d %H:%M}   {}".format(w["name"], w["lo"], w["hi"], "scored" if w["scored"] else "not scored"))
    if a.dry_run:
        log("\nDRY RUN: all inputs were read and the runs and windows listed. Nothing was written.")
        return 0

    # ---- numbers
    gauge_h = hourly_mean(gauge * 1.0)
    run_h = {r["key"]: hourly_mean(series[r["key"]]) for r in runs}
    vols_h = {"gauge": vol_per_hour(gauge_h)}
    vols_h.update({k: vol_per_hour(h) for k, h in run_h.items()})
    cum = {k: v.cumsum() for k, v in vols_h.items()}
    keys = ["gauge"] + [r["key"] for r in runs]
    cols = ["gauge"] + [r["label"] for r in runs]

    def line(label, vals):
        gv = vals[0]
        s = "  {:24s}{:>9s}".format(label, fmt_v(gv))
        for v in vals[1:]:
            s += "{:>14s}".format("{} ({})".format(fmt_v(v), fmt_r(ratio(v, gv))))
        return s

    head = "  {:24s}{:>9s}".format("", "gauge") + "".join("{:>14s}".format(c[:13]) for c in cols[1:])
    rows = []
    log("\n1. Volume in each scored storm (M m3), run (run / gauge)")
    log(head)
    for w in storms:
        vals = [window_volume(vols_h[k], w["lo"], w["hi"]) for k in keys]
        log(line(w["name"], vals))
        rows.append(("storm: " + w["name"], w["lo"], w["hi"], vals))
    vals = [sum(window_volume(vols_h[k], w["lo"], w["hi"]) for w in storms) for k in keys]
    log(line("seven storms", vals))
    rows.append(("seven storms together", None, None, vals))
    log("\n2. Volume in each month (M m3; months by the start of each hour)")
    log(head)
    mv = {k: month_volumes(vols_h[k]) for k in keys}
    for m in MONTHS:
        vals = [mv[k][m] for k in keys]
        log(line(m, vals))
        rows.append(("month: " + m, pd.Timestamp(m + "-01"), pd.Timestamp(m + "-01") + pd.offsets.MonthBegin(1), vals))
    vals = [float(vols_h[k].sum()) for k in keys]
    log(line("season total", vals))
    rows.append(("season total (Sep 1 to May 31)", RUN_START, SEASON_END, vals))
    log("\n3. Volume OUTSIDE the seven scored windows, Sep to May (M m3; the Sep 23 storm is NOT left out)")
    log(head)
    vals = [outside_volume(vols_h[k], storms) for k in keys]
    log(line("outside the windows", vals))
    rows.append(("outside the seven scored windows", None, None, vals))

    # ---- checks
    check_msg = "CHECK not run (--skip-checks)."
    bad_checks = []
    if not a.skip_checks:
        bad_checks, n_checked = run_checks(vols_h, storms, [r["key"] for r in runs])
        who = "gauge" + (" and Re-tuned" if RETUNED in [r["key"] for r in runs] else "")
        if bad_checks:
            log("\nCHECK against the write-up (seven storms) and the Nov 7 edge check (months), {}: MISMATCH ({} of {} values differ by more than {:g} M m3):".format(who, len(bad_checks), n_checked, MATCH_TOL))
            for b in bad_checks:
                log("  " + b)
            log("  The figure was still drawn, but do not use it until this is explained.")
            check_msg = "CHECK: MISMATCH for the {} volumes (see the printed list); do not use this figure.".format(who)
        else:
            log("\nCHECK against the write-up (seven storms) and the Nov 7 edge check (months), {}: MATCH. All {} values agree within {:g} M m3.".format(who, n_checked, MATCH_TOL))
            check_msg = "CHECK: the {} volumes in the seven storms and in each month equal the write-up and the Nov 7 edge check (within {:g} M m3). Other runs are not checked.".format(who, MATCH_TOL)

    # ---- files
    tags = "_".join(r["tag"] for r in runs)
    fig_path = out_dir / "fig_season_{}.png".format(tags)
    csv_path = out_dir / "season_volumes_{}.csv".format(tags)
    replaced = [p.name for p in (fig_path, csv_path) if p.exists()]
    out_dir.mkdir(parents=True, exist_ok=True)
    recs = []
    for name, lo, hi, vals in rows:
        rec = {"period": name, "start": lo, "end": hi, "gauge_Mm3": vals[0]}
        for r, v in zip(runs, vals[1:]):
            rec["{}_Mm3".format(r["tag"])] = v
            rec["{}_over_gauge".format(r["tag"])] = ratio(v, vals[0])
        recs.append(rec)
    df = pd.DataFrame(recs)
    num = [c for c in df.columns if c.endswith("_Mm3") or c.endswith("_over_gauge")]
    df[num] = df[num].round(5)
    df.to_csv(csv_path, index=False)
    make_figure(fig_path, runs, gauge_h, run_h, cum, rain, storms, extras, a, check_msg)
    log("\nWrote:\n  {}\n  {}".format(fig_path, csv_path))
    if replaced:
        log("  (replaced existing file(s): {})".format(", ".join(replaced)))
    return 1 if bad_checks else 0


if __name__ == "__main__":
    sys.exit(main())
