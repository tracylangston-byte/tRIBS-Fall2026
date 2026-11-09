#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_wy2020_hydrographs_sycamore.py   (v1.0, 2026-11-04)
=========================================================

WY2020 storm hydrographs: the USGS gauge against up to FOUR model runs, one figure per storm.

THE FOUR RUNS (all drawn against the gauge, all on the same axes)
  Re-tuned: loss off      Sycamore_WY2020_Retuned_Outlet.qout            solid line, purple
  Original: loss off      Sycamore_WY2020_Original_Outlet.qout           dashed line, orange
  Original: Kchan 0.5     Sycamore_WY2020_Original_kchan0p5_Outlet.qout  dotted line, blue
  Original: Kchan 0.2     Sycamore_WY2020_Original_kchan0p2_Outlet.qout  dash-dot line, green
  The colors match the WY2017 figures. The line styles differ as well, so the runs can be told apart without color.
  The gauge is the thick black line. Use --runs to draw only some of the four (for example --runs retuned,original).

WHAT EACH FIGURE HAS
  One PNG per storm, two columns: LEFT = linear flow axis, RIGHT = log flow axis.
    top row     hyetograph: rain the model received (basin mean, mm/h; hangs from the top)
    bottom row  flow at the outlet in m3/s: the gauge and the model runs
  The linear column shows the peaks. The log column shows small flows (the gauge in Nov 19-22 is about 0.1 m3/s at its
  highest, while the model is above 50 m3/s, so on a linear axis the gauge looks flat). A pale band marks the day
  before the scored window, for context only.
  The title states the gauge and model volumes (million m3) for the scored window.

ALSO WRITTEN
  fig_wy2020_contact_sheet.png        every storm on one page (linear axes) so you can choose which ones to use
  wy2020_hydrograph_windows.csv       window, rain, volume, peak and peak time for the gauge and each run
  A volume table is printed. For the seven scored storms it is checked against the numbers in the WY2020 write-up
  (gauge, Re-tuned, Original, Kchan 0.5, Kchan 0.2) and says MATCH or MISMATCH for each.

WINDOWS (same rule the scoring used, so the volumes here equal the volumes in the write-up)
  scored window = 00:00 on the first rain day to 00:00 four days after the last rain day, as listed in storm_table.csv.
  Volume = hourly mean flow x 3,600 s, summed over the hours that have data (the same rule as the scoring).
  Default storms = the seven that were scored (WY2020 rows with status "On 20 mm list").

UNITS
  Flow is in m3/s (what the model and gauge files use). Volumes are in million m3 (M m3).
  Peaks and peak timing were NOT pre-registered and are NOT scored. They are printed and saved for description only.

WHAT IT READS (read-only; nothing is changed, no model is run)
  storm_table.csv                                   columns year, event, first_day, last_day, status (and basin_mm, used as a rain check)
  storm_features_out/hourly_basin_wy2020.csv        columns interval_start_mst, basin_mean_mm
  results/<run>_Outlet.qout                         first two columns: Time_hr, flow (m3/s); Time_hr 0 = Sep 1 2019 00:00
  ../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv   columns datetime, q_cms
WHAT IT WRITES  (in wy2020_hydro_figs/ unless you set --out-dir)
  fig_wy2020_<storm>.png, fig_wy2020_contact_sheet.png, wy2020_hydrograph_windows.csv

COMMANDS (from workspaces/Sycamore/lab/)
  python plot_wy2020_hydrographs_sycamore.py --dry-run                 # checks the files and lists the windows; writes nothing
  python plot_wy2020_hydrographs_sycamore.py                           # the seven scored storms, all four runs
  python plot_wy2020_hydrographs_sycamore.py --storms "Nov 19-22,Dec 23-28"      # only these (names as in storm_table.csv)
  python plot_wy2020_hydrographs_sycamore.py --storms all              # every WY2020 storm from October on (September storms are spin-up)
  python plot_wy2020_hydrographs_sycamore.py --runs retuned,original   # only two runs
  python plot_wy2020_hydrographs_sycamore.py --runs original,kchan05,kchan02 --storms "Nov 19-22"     # the channel-loss comparison
  python plot_wy2020_hydrographs_sycamore.py --custom "2019-11-19 00:00" "2019-12-02 00:00" --custom-name "Nov 19 to Dec 1"
  python plot_wy2020_hydrographs_sycamore.py --no-contact              # skip the contact sheet

WORDS
  hyetograph  a bar chart of rain against time.   hydrograph  flow at one place against time.
  log axis    each gridline is 10 times the one below it; it shows small flows next to big ones.
  hour labels the model labels an hour by its END (as the scoring does); the rain file labels by the START, and the rain bars
              are drawn over the hour they cover.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUN_START = pd.Timestamp("2019-09-01 00:00")        # Time_hr 0 in the WY2020 outlet files
AFTER_DAYS = 3                                      # scored window ends 3 days after the day AFTER the last rain day
FIG_DIR = "wy2020_hydro_figs"
CSV_NAME = "wy2020_hydrograph_windows.csv"
CONTACT_NAME = "fig_wy2020_contact_sheet.png"
MATCH_TOL = 0.002                                   # write-up volumes are given to 3 decimals (M m3)

INK, INK_2, GRID = "#0b0b0b", "#52514e", "#d9d8d3"
RAIN_BLUE, PAD_GRAY = "#8db4e2", "#f0efec"
GAUGE_COLOR = "#0b0b0b"

# key, label, color, line style, outlet file stem. Colors are the WY2017 figure colors (checked for color-blind separation);
# the line styles are the second cue, so no run depends on color alone.
RUNS = [
    ("retuned", "Re-tuned: loss off", "#6b3fa0", "-", "Sycamore_WY2020_Retuned"),
    ("original", "Original: loss off", "#eb6834", "--", "Sycamore_WY2020_Original"),
    ("kchan05", "Original: Kchan 0.5", "#2a78d6", ":", "Sycamore_WY2020_Original_kchan0p5"),
    ("kchan02", "Original: Kchan 0.2", "#4a9d6b", "-.", "Sycamore_WY2020_Original_kchan0p2"),
]

# Scored storm-window volumes (M m3) as written in the WY2020 validation write-up. Used only as a check.
EXPECTED = {
    "Nov 19-22": {"gauge": 0.015, "retuned": 2.800, "original": 3.165, "kchan05": 3.019, "kchan02": 3.091},
    "Nov 28-30": {"gauge": 1.473, "retuned": 3.235, "original": 3.538, "kchan05": 3.341, "kchan02": 3.433},
    "Dec 23-28": {"gauge": 4.069, "retuned": 0.218, "original": 0.154, "kchan05": 0.044, "kchan02": 0.083},
    "Feb 10-11": {"gauge": 0.386, "retuned": 1.767, "original": 1.765, "kchan05": 1.641, "kchan02": 1.700},
    "Feb 22-23": {"gauge": 4.131, "retuned": 4.839, "original": 4.981, "kchan05": 4.825, "kchan02": 4.899},
    "Mar 10-13": {"gauge": 5.626, "retuned": 0.875, "original": 0.493, "kchan05": 0.340, "kchan02": 0.413},
    "Mar 18-19": {"gauge": 3.317, "retuned": 1.031, "original": 0.841, "kchan05": 0.683, "kchan02": 0.775},
}


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    sys.exit("\nSTOPPED: " + msg + "\n")


# ------------------------------------------------------------------------------------------------ reading
def need(path: Path, what: str) -> Path:
    if not Path(path).is_file():
        die("cannot find the {} file:\n    {}\n  Run this from workspaces/Sycamore/lab/, or point to the file with the matching "
            "option (see the top of this script).".format(what, path))
    return Path(path)


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


# storm_table.csv comes in two header styles (short names, or the readable names in the formatted copy). Both are accepted.
STORM_COLUMN_NAMES = {
    "year": ["year", "Water year"],
    "event": ["event", "Event"],
    "first_day": ["first_day", "First rain day"],
    "last_day": ["last_day", "Last rain day"],
    "status": ["status", "List status"],
    "basin_mm": ["basin_mm", "Basin rain (mm)"],          # optional (only used for a rain cross-check)
}
# Used ONLY if the date columns cannot be read as full dates (for example "Nov 19" with no year). A note is printed when used.
STORM_DAYS_FALLBACK = {
    "Sep 10": ("2019-09-10", "2019-09-10"), "Sep 23": ("2019-09-23", "2019-09-23"),
    "Nov 19-22": ("2019-11-19", "2019-11-22"), "Nov 28-30": ("2019-11-28", "2019-11-30"),
    "Dec 3-4": ("2019-12-03", "2019-12-04"), "Dec 7-9": ("2019-12-07", "2019-12-09"),
    "Dec 23-28": ("2019-12-23", "2019-12-28"), "Jan 21": ("2020-01-21", "2020-01-21"),
    "Feb 10-11": ("2020-02-10", "2020-02-11"), "Feb 22-23": ("2020-02-22", "2020-02-23"),
    "Mar 10-13": ("2020-03-10", "2020-03-13"), "Mar 18-19": ("2020-03-18", "2020-03-19"),
}


def read_storms(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    t = pd.DataFrame()
    missing = []
    for want, options in STORM_COLUMN_NAMES.items():
        hit = [c for c in options if c in raw.columns]
        if hit:
            t[want] = raw[hit[0]]
        elif want != "basin_mm":
            missing.append(want)
    if missing:
        die("{} is missing column(s) for {} (looked for {}).\n  All columns found: {}".format(
            path, missing, [STORM_COLUMN_NAMES[m] for m in missing], list(raw.columns)))
    t = t[t["year"].astype(str).str.contains("2020")].copy()
    if t.empty:
        die("no WY2020 rows in {}. The water-year column holds: {}".format(path, sorted(set(raw.iloc[:, 0].astype(str)))[:6]))
    # en dashes in event names (Nov 19–22) would not match the plain hyphen used elsewhere
    t["event"] = t["event"].astype(str).str.replace("–", "-").str.replace("—", "-").str.strip()
    t["status"] = t["status"].astype(str).str.strip()
    t["first_day"] = pd.to_datetime(t["first_day"], errors="coerce")
    t["last_day"] = pd.to_datetime(t["last_day"], errors="coerce")
    lo_ok, hi_ok = pd.Timestamp("2019-09-01"), pd.Timestamp("2020-09-30")
    bad = t["first_day"].isna() | t["last_day"].isna() | (t["first_day"] < lo_ok) | (t["last_day"] > hi_ok)
    if bad.any():
        unknown = [e for e in t.loc[bad, "event"] if e not in STORM_DAYS_FALLBACK]
        if unknown:
            die("the dates in {} could not be read for: {}.\n  Raw first/last day values: {}".format(
                path, unknown, raw.loc[t.index[bad], [c for c in STORM_COLUMN_NAMES["first_day"] + STORM_COLUMN_NAMES["last_day"] if c in raw.columns]].head(3).values.tolist()))
        for i in t.index[bad]:
            a, b = STORM_DAYS_FALLBACK[t.at[i, "event"]]
            t.at[i, "first_day"], t.at[i, "last_day"] = pd.Timestamp(a), pd.Timestamp(b)
        log("  NOTE: the dates in storm_table.csv were not full dates for {} storm(s); the built-in WY2020 storm days were used for them.".format(int(bad.sum())))
    return t.reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ numbers
def volume(q: pd.Series, lo, hi):
    """Million m3 from hourly mean flow x 3,600 s over [lo, hi), counting only hours with data. Returns (volume, hours)."""
    s = q[(q.index >= lo) & (q.index < hi)]
    h = s.resample("60min").mean().dropna()
    return (float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")), int(len(h))


def peak(q: pd.Series, lo, hi):
    s = q[(q.index >= lo) & (q.index < hi)].dropna()
    if not len(s):
        return float("nan"), pd.NaT
    return float(s.max()), s.idxmax()


def fmt_v(v: float) -> str:
    if not np.isfinite(v):
        return "n/a"
    return "{:.3f}".format(v) if v < 0.1 else "{:.2f}".format(v)


def logsafe(y) -> np.ndarray:
    y = np.asarray(y, float).copy()
    y[y <= 0] = np.nan
    return y


# ------------------------------------------------------------------------------------------------ drawing
def tidy(ax) -> None:
    ax.grid(alpha=0.35, color=GRID, which="major", lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=8.5)


def draw_flow(ax, gauge: pd.Series, runs: list, lo_plot, lo, hi, scale: str, log_min: float, lw: float = 1.8, gauge_lw: float = 2.4):
    """Gauge (thick black, on top) and the runs, in [lo_plot, hi]. runs = list of (label, color, style, series)."""
    ax.axvspan(lo_plot, lo, color=PAD_GRAY, lw=0, zorder=0)
    tops = []
    for label, color, ls, q in runs:
        w = q[(q.index >= lo_plot) & (q.index <= hi)]
        y = logsafe(w.values) if scale == "log" else w.values
        ax.plot(w.index, y, color=color, ls=ls, lw=lw, label=label, zorder=3)
        if np.isfinite(logsafe(w.values)).any():
            tops.append(np.nanmax(logsafe(w.values)))
    g = gauge[(gauge.index >= lo_plot) & (gauge.index <= hi)]
    gy = logsafe(g.values) if scale == "log" else g.values
    ax.plot(g.index, gy, color=GAUGE_COLOR, lw=gauge_lw, label="USGS gauge 09510200", zorder=5)
    if np.isfinite(logsafe(g.values)).any():
        tops.append(np.nanmax(logsafe(g.values)))
    if scale == "log":
        ax.set_yscale("log")
        ax.set_ylim(log_min, max((max(tops) * 3) if tops else 1.0, log_min * 100))
    else:
        ax.set_ylim(bottom=0)
    ax.set_xlim(lo_plot, hi)


def make_storm_figure(out: Path, name: str, lo, hi, lo_plot, rain: pd.Series, gauge, runs, vols, rain_mm, args) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 7.6), sharex="col",
                             gridspec_kw={"height_ratios": [0.7, 2.5], "hspace": 0.07, "wspace": 0.17})
    r = rain[(rain.index >= lo_plot) & (rain.index < hi)]
    for j, scale in enumerate(("linear", "log")):
        ax_r, ax_q = axes[0, j], axes[1, j]
        ax_r.axvspan(lo_plot, lo, color=PAD_GRAY, lw=0, zorder=0)
        ax_r.bar(r.index, r.values, width=pd.Timedelta(hours=1), align="edge", color=RAIN_BLUE, zorder=2)
        ax_r.invert_yaxis()
        ax_r.set_xlim(lo_plot, hi)
        draw_flow(ax_q, gauge, runs, lo_plot, lo, hi, scale, args.log_min)
        ax_q.set_xlabel("Phoenix time (MST)", color=INK_2, fontsize=9)
        loc = mdates.AutoDateLocator()
        ax_q.xaxis.set_major_locator(loc)
        ax_q.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
        for ax in (ax_r, ax_q):
            tidy(ax)
        ax_r.set_title("Linear flow axis" if scale == "linear" else "Log flow axis", loc="left", fontsize=10, color=INK_2)
        if j == 0:
            ax_r.set_ylabel("Rain (mm/h)\nbasin mean", color=INK_2, fontsize=9)
            ax_q.set_ylabel("Flow at the outlet (m$^3$/s)", color=INK_2, fontsize=9)
        else:
            ax_q.set_ylabel("Flow at the outlet (m$^3$/s), log scale", color=INK_2, fontsize=9)
    if pd.Timestamp(lo_plot) < pd.Timestamp(lo):
        axes[1, 0].text(0.004, 0.97, "before the scored window", transform=axes[1, 0].transAxes, fontsize=7.5, color=INK_2, va="top")

    gv = vols["gauge"]
    model_vals = [vols[k] for k, _, _, _ in [(r_[0], 0, 0, 0) for r_ in RUNS] if k in vols and k != "gauge" and np.isfinite(vols[k])]
    if model_vals:
        title = "{}: gauge {} M m³, model {} to {} M m³".format(name, fmt_v(gv), fmt_v(min(model_vals)), fmt_v(max(model_vals)))
    else:
        title = "{}: gauge {} M m³".format(name, fmt_v(gv))
    fig.suptitle(title, x=0.01, y=0.985, ha="left", fontsize=14, color=INK)
    sub = "Rain {:.1f} mm over the scored window ({:%b %d %H:%M} to {:%b %d %H:%M}).  Volumes in M m³: gauge {}".format(
        rain_mm, lo, hi, fmt_v(gv))
    for k, label, _, _, _ in RUNS:
        if k in vols and k != "gauge":
            sub += ";  {} {}".format(label, fmt_v(vols[k]))
    fig.text(0.01, 0.945, sub, fontsize=8.8, color=INK_2, ha="left", va="top")

    handles = [Line2D([], [], color=GAUGE_COLOR, lw=2.4, label="USGS gauge 09510200")]
    handles += [Line2D([], [], color=c, ls=ls, lw=1.8, label=lab) for lab, c, ls, _ in runs]
    handles += [Line2D([], [], color=RAIN_BLUE, lw=6, label="rain the model received")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.045), ncol=len(handles), frameon=False, fontsize=9)
    foot = ("Scored window: 00:00 on the first rain day to 00:00 four days after the last rain day; volumes use hourly means (the scoring rule). "
            "Peaks and timing were not pre-registered or scored; read them as description only.")
    fig.text(0.01, 0.005, foot, fontsize=7.4, color=INK_2, ha="left", va="bottom")
    fig.subplots_adjust(left=0.07, right=0.985, top=0.86, bottom=0.17)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)
    return str(out)


def make_contact_sheet(out: Path, items: list, runs, args) -> str:
    """items = list of dicts with name, lo, hi, lo_plot, gauge, vols. One small linear panel per storm."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    n = len(items)
    ncol = 4 if n > 3 else n
    nrow = int(np.ceil((n + 1) / float(ncol)))          # one extra cell for the legend
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.9 * ncol, 3.1 * nrow + 0.6), squeeze=False)
    flat = list(axes.ravel())
    for ax, it in zip(flat, items):
        draw_flow(ax, it["gauge"], runs, it["lo_plot"], it["lo"], it["hi"], "linear", args.log_min, lw=1.4, gauge_lw=2.0)
        loc = mdates.AutoDateLocator(minticks=3, maxticks=4)
        ax.xaxis.set_major_locator(loc)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
        tidy(ax)
        ax.set_title("{}\ngauge {} M m³".format(it["name"], fmt_v(it["vols"]["gauge"])), fontsize=9.5, color=INK, loc="left")
        ax.set_ylabel("m$^3$/s", color=INK_2, fontsize=8.5)
    legend_ax = flat[n] if n < len(flat) else None
    for ax in flat[n:]:
        ax.axis("off")
    handles = [Line2D([], [], color=GAUGE_COLOR, lw=2.4, label="USGS gauge 09510200")]
    handles += [Line2D([], [], color=c, ls=ls, lw=1.8, label=lab) for lab, c, ls, _ in runs]
    if legend_ax is not None:
        legend_ax.legend(handles=handles, loc="center", frameon=False, fontsize=9.5)
    else:
        fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, fontsize=9)
    fig.suptitle("WY2020 storms: gauge against the model runs (linear flow axes, each panel on its own scale)", x=0.01, y=0.995, ha="left", fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)
    return str(out)


# ------------------------------------------------------------------------------------------------ main
def slug(name: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in name).strip("_")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="WY2020 hydrographs: gauge against up to four model runs (see the top of this file).")
    ap.add_argument("--root", default=str(HERE), help="the lab folder (default: the folder this script is in)")
    ap.add_argument("--table", default=None, help="storm_table.csv (default <root>/storm_table.csv)")
    ap.add_argument("--hourly", default=None, help="hourly_basin_wy2020.csv (default <root>/storm_features_out/hourly_basin_wy2020.csv)")
    ap.add_argument("--results", default=None, help="folder with the *_Outlet.qout files (default <root>/results)")
    ap.add_argument("--obs", default=None, help="USGS gauge file for WY2020 (default <root>/../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv)")
    ap.add_argument("--out-dir", default=None, help="where the PNGs and the CSV go (default <root>/" + FIG_DIR + ")")
    ap.add_argument("--runs", default="retuned,original,kchan05,kchan02", help="comma list of: retuned, original, kchan05, kchan02 (default all four)")
    ap.add_argument("--storms", default=None, help="'scored' (the seven scored storms; default), 'all' (every WY2020 storm from October on), or a comma list of names, e.g. \"Nov 19-22,Dec 23-28\"")
    ap.add_argument("--custom", nargs=2, metavar=("START", "END"), default=None, help="one extra window, for example --custom \"2019-11-19 00:00\" \"2019-12-02 00:00\"")
    ap.add_argument("--custom-name", default="Custom window")
    ap.add_argument("--pad-days", type=float, default=1.0, help="days drawn before the window, for context only (default 1)")
    ap.add_argument("--log-min", type=float, default=0.001, help="lowest value on the log axes, m3/s (default 0.001)")
    ap.add_argument("--no-contact", action="store_true", help="skip the contact sheet")
    ap.add_argument("--dry-run", action="store_true", help="check the files and list the windows; write nothing")
    ap.add_argument("--dpi", type=int, default=150)
    a = ap.parse_args(argv)

    root = Path(a.root)
    table = Path(a.table) if a.table else root / "storm_table.csv"
    hourly = Path(a.hourly) if a.hourly else root / "storm_features_out" / "hourly_basin_wy2020.csv"
    results = Path(a.results) if a.results else root / "results"
    obs_path = Path(a.obs) if a.obs else root.parent / "init_data" / "obs" / "sycamore_09510200_obs_flow_WY2020.csv"
    out_dir = Path(a.out_dir) if a.out_dir else root / FIG_DIR

    want = [k.strip() for k in a.runs.split(",") if k.strip()]
    known = [r[0] for r in RUNS]
    bad = [k for k in want if k not in known]
    if bad or not want:
        die("--runs must be a comma list of {} (got {}).".format(", ".join(known), bad or "nothing"))
    chosen = [r for r in RUNS if r[0] in want]

    log("\nWY2020 hydrographs: reading files")
    storms = read_storms(need(table, "storm table"))
    rain = read_rain(need(hourly, "hourly basin rain"))
    gauge = read_obs(need(obs_path, "WY2020 gauge"))
    log("  storm table: {} WY2020 rows".format(len(storms)))
    log("  rain: {} hours, {:%b %d %Y} to {:%b %d %Y}".format(len(rain), rain.index.min(), rain.index.max()))
    log("  gauge: {} readings, {:%b %d %Y} to {:%b %d %Y}".format(gauge.notna().sum(), gauge.index.min(), gauge.index.max()))

    runs_found, runs_missing = [], []
    for key, label, color, ls, stem in chosen:
        p = results / "{}_Outlet.qout".format(stem)
        (runs_found if p.is_file() else runs_missing).append((key, label, color, ls, stem, p))
    for key, label, _, _, _, p in runs_missing:
        log("  NOT FOUND (left out): {}  ->  {}".format(label, p))
    if not runs_found:
        die("none of the chosen outlet files were found in {}.".format(results))
    series = {}
    for key, label, color, ls, stem, p in runs_found:
        series[key] = read_qout(p)
        q = series[key]
        log("  {}: {} rows, {:%b %d %Y} to {:%b %d %Y}".format(label, len(q), q.index.min(), q.index.max()))
    runs = [(label, color, ls, series[key]) for key, label, color, ls, _, _ in runs_found]
    run_keys = [key for key, _, _, _, _, _ in runs_found]

    # which windows
    windows = []                       # dict(name, lo, hi, scored)
    spec = a.storms
    if spec is None:
        spec = "scored" if a.custom is None else ""
    if spec == "scored":
        sel = storms[storms["status"].astype(str).str.startswith("On 20")]
        if sel.empty:
            die("no WY2020 row has a status starting with 'On 20'. The statuses found are: {}\n  Use --storms with names instead.".format(sorted(set(storms["status"]))))
        if len(sel) != 7:
            log("  NOTE: found {} storms with status starting 'On 20' in storm_table.csv; the write-up scored 7.".format(len(sel)))
    elif spec == "all":
        sel = storms[storms["first_day"] >= pd.Timestamp("2019-10-01")]      # September storms are spin-up
    elif spec:
        names = [s.strip() for s in spec.split(",") if s.strip()]
        absent = [n for n in names if n not in set(storms["event"])]
        if absent:
            die("these storms are not in the WY2020 rows of {}: {}\n  Names found: {}".format(table, absent, ", ".join(storms["event"])))
        sel = storms[storms["event"].isin(names)]
        sel = sel.set_index("event").loc[names].reset_index()
    else:
        sel = storms.iloc[0:0]
    for _, r in sel.iterrows():
        lo = pd.Timestamp(r["first_day"])
        hi = pd.Timestamp(r["last_day"]) + pd.Timedelta(days=1 + AFTER_DAYS)
        windows.append({"name": r["event"], "lo": lo, "hi": hi, "first": lo, "last": pd.Timestamp(r["last_day"]),
                        "table_mm": float(r["basin_mm"]) if "basin_mm" in storms.columns and pd.notna(r.get("basin_mm")) else float("nan"),
                        "check": r["event"] in EXPECTED and str(r["status"]).startswith("On 20")})
    if a.custom:
        lo, hi = pd.Timestamp(a.custom[0]), pd.Timestamp(a.custom[1])
        if not hi > lo:
            die("--custom END must be later than START.")
        windows.append({"name": a.custom_name, "lo": lo, "hi": hi, "first": None, "last": None, "table_mm": float("nan"), "check": False})
    if not windows:
        die("no storms selected.")

    log("\nWindows ({}):".format(len(windows)))
    for w in windows:
        log("  {:16s} {:%Y-%m-%d %H:%M} to {:%Y-%m-%d %H:%M}".format(w["name"], w["lo"], w["hi"]))
    if a.dry_run:
        log("\nDRY RUN: all inputs were read and the windows listed. Nothing was written.")
        return 0

    # numbers
    rows, items, bad_checks, warn = [], [], [], []
    log("\nVolumes (M m3), scored-window rule; peaks in m3/s are description only (not scored)")
    head = "  {:14s}{:>9s}".format("storm", "gauge") + "".join("{:>11s}".format(k) for k in run_keys)
    log(head)
    for w in windows:
        lo, hi = w["lo"], w["hi"]
        gv, gh = volume(gauge, lo, hi)
        vols = {"gauge": gv}
        gp, gpt = peak(gauge, lo, hi)
        hrs = int((hi - lo) / pd.Timedelta(hours=1))
        rain_win = float(rain[(rain.index >= lo) & (rain.index < hi)].sum())
        rows.append({"storm": w["name"], "window_start": lo, "window_end": hi, "run_key": "gauge", "series": "gauge", "volume_Mm3": gv, "hours_with_data": gh,
                     "hours_in_window": hrs, "peak_m3s": gp, "peak_time": gpt, "peak_lag_h_vs_gauge": 0.0, "rain_mm_window": rain_win})
        for key, label, _, _, _, _ in runs_found:
            v, vh = volume(series[key], lo, hi)
            pk, pkt = peak(series[key], lo, hi)
            vols[key] = v
            lag = float((pkt - gpt) / pd.Timedelta(hours=1)) if pd.notna(pkt) and pd.notna(gpt) else float("nan")
            rows.append({"storm": w["name"], "window_start": lo, "window_end": hi, "run_key": key, "series": label, "volume_Mm3": v, "hours_with_data": vh,
                         "hours_in_window": hrs, "peak_m3s": pk, "peak_time": pkt, "peak_lag_h_vs_gauge": lag, "rain_mm_window": rain_win})
        log("  {:14s}{:>9s}".format(w["name"], fmt_v(gv)) + "".join("{:>11s}".format(fmt_v(vols[k])) for k in run_keys))
        if gh < 0.97 * hrs:
            warn.append("{}: the gauge has data for {} of {} hours; the volume is scored on the hours with data (as in the write-up). Dec 23-28 is expected to show this.".format(w["name"], gh, hrs))
        if w["first"] is not None and np.isfinite(w["table_mm"]):
            storm_rain = float(rain[(rain.index >= w["first"]) & (rain.index < w["last"] + pd.Timedelta(days=1))].sum())
            if abs(storm_rain - w["table_mm"]) > 0.1:
                warn.append("{}: rain over the storm days is {:.2f} mm here but {:.2f} mm in storm_table.csv".format(w["name"], storm_rain, w["table_mm"]))
        if w["check"]:
            exp = EXPECTED[w["name"]]
            for key in ["gauge"] + run_keys:
                if key in exp and not (np.isfinite(vols[key]) and abs(vols[key] - exp[key]) <= MATCH_TOL):
                    bad_checks.append("{} {}: here {} vs write-up {:.3f}".format(w["name"], key, fmt_v(vols[key]), exp[key]))
        # peaks
        w["vols"], w["rain_mm"] = vols, rain_win
        w["lo_plot"] = lo - pd.Timedelta(days=a.pad_days)

    short = {"gauge": "gauge", "retuned": "Re-tuned", "original": "Original", "kchan05": "Kchan 0.5", "kchan02": "Kchan 0.2"}
    log("\nPeaks (m3/s) and peak time relative to the gauge (hours; positive = model later). Description only, not scored:")
    for w in windows:
        parts = []
        for r in rows:
            if r["storm"] == w["name"]:
                lag = ""
                if r["run_key"] != "gauge" and np.isfinite(r["peak_lag_h_vs_gauge"]):
                    lag = " ({:+.1f} h)".format(r["peak_lag_h_vs_gauge"])
                parts.append("{} {:.3g}{}".format(short[r["run_key"]], r["peak_m3s"], lag))
        log("  {:16s}".format(w["name"]) + "; ".join(parts))

    n_check = sum(1 for w in windows if w["check"])
    if n_check:
        if bad_checks:
            log("\nCHECK against the write-up: MISMATCH ({} value(s) differ by more than {:g} M m3):".format(len(bad_checks), MATCH_TOL))
            for b in bad_checks:
                log("  " + b)
            log("  The figures were still drawn, but do not use them until this is explained.")
        else:
            log("\nCHECK against the write-up: MATCH. All volumes for {} scored storm(s) equal the write-up values (within {:g} M m3).".format(n_check, MATCH_TOL))
    for x in warn:
        log("WARNING: " + x)

    # figures
    written = []
    for w in windows:
        fn = out_dir / "fig_wy2020_{}.png".format(slug(w["name"]))
        written.append(make_storm_figure(fn, w["name"], w["lo"], w["hi"], w["lo_plot"], rain, gauge, runs, w["vols"], w["rain_mm"], a))
    if not a.no_contact and len(windows) > 1:
        written.append(make_contact_sheet(out_dir / CONTACT_NAME, [{"name": w["name"], "lo": w["lo"], "hi": w["hi"], "lo_plot": w["lo_plot"], "gauge": gauge, "vols": w["vols"]} for w in windows], runs, a))
    csv = out_dir / CSV_NAME
    pd.DataFrame(rows).to_csv(csv, index=False)
    written.append(str(csv))
    log("\nWrote:\n  " + "\n  ".join(written))
    return 1 if bad_checks else 0


if __name__ == "__main__":
    sys.exit(main())
