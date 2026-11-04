#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_w5_rain_sycamore.py   (v1.0, 2026-11-03)
=============================================

A figure that puts four things for the whole water year 2017 (Oct 2016 - Sep 2017) on one time axis:

  TOP     monthly outlet volume (M m3): USGS gauge, the FIRST W5 run (re-tuned, loss off), and the three
          "draw 12 / k10 as drawn" runs (loss off, Kchan 0.5, Kchan 0.2).
  MIDDLE  a hyetograph: daily radar rain (grid mean, mm/day) hanging down from the top, with the daily mean
          of the two good rain gauges (Sycamore Upper and Iron Dike) as black dots.
  BOTTOM  the radar / gauge ratio for each month, two ways: all 3 gauges, and 2 gauges (Mt. Ord left out).
          1.0 means the radar and the gauges agree. Above 1.0 the radar read higher than the gauges.

The July - September flow panel of the old fig_w5_kchan.png is NOT drawn here.

WORDS
  hyetograph     a bar chart of rain against time, drawn upside down (bars hang from the top) so it sits
                 above a flow plot the way weather and flow are usually shown.
  ratio          radar rain divided by gauge rain over the same days. The radar value is the mean of the 3 x 3
                 radar cells around each gauge (the "3x3 mean" from the radar check), so it compares like with like.
  grid mean      the average over every cell of the 32 x 49 radar grid. This is what the model sees on average.

WHAT IT READS (nothing is changed, nothing is run)
  results/<run id>_Outlet.qout          the finished W5 outlet files
  ../init_data/obs/...WY2017.csv        the USGS gauge file (through run_tail_screen_sycamore.py)
  radar_check_out/radar_gauge_daily.csv, radar_check_out_aprsep_all3/radar_gauge_daily.csv
                                        daily rain tables from the Stage 0 radar check (Oct-Mar and Apr-Sep)
WHAT IT WRITES
  w5_kchan_figs/fig_w5_kchan_rain.png    the figure (the old fig_w5_kchan.png is not touched)
  w5_rain_ratio_months.csv               the monthly numbers behind the bottom panel

COMMANDS (run from lab/)
  python plot_w5_rain_sycamore.py --table-only     # print the monthly rain table only; needs no model output
  python plot_w5_rain_sycamore.py                  # table plus figure
Options: --rain-dirs a,b   --out path.png   --exclude "Mt. Ord"   --min-gauge-mm 5
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

DEFAULT_RAIN_DIRS = ["radar_check_out", "radar_check_out_aprsep_all3"]   # Oct-Mar first, then Apr-Sep
FIG_DIR = "w5_kchan_figs"
OUT_NAME = "fig_w5_kchan_rain.png"
MONTHS_CSV = "w5_rain_ratio_months.csv"
MIN_GAUGE_MM = 5.0              # a monthly ratio needs at least this much gauge rain (summed over the gauges), as in the radar check
EXCLUDE = "Mt. Ord"             # the gauge left out of the "2 gauges" ratio (reads low in cold storms)
WY_FIRST, WY_LAST = "2016-10", "2017-09"
GAUGE_GRAY, INK, INK_2, GRID = "#8a8984", "#0b0b0b", "#52514e", "#d9d8d3"
RAIN_BLUE = "#8db4e2"
RATIO3_COLOR, RATIO2_COLOR = "#333333", "#c0392b"
RETUNED_COLOR = "#6b3fa0"
RUN_COLORS = ["#2a78d6", "#eb6834", "#4a9d6b"]      # as drawn: loss off, Kchan 0.5, Kchan 0.2 (same colors as the old figure)
NAN = float("nan")


def log(msg: str = "") -> None:
    print(msg, flush=True)


# --------------------------------------------------------------------------------------------
# Rain tables
# --------------------------------------------------------------------------------------------
def load_daily(lab: Path, dirs: list):
    """Stack the radar_gauge_daily.csv files (first folder wins when the same day is in two). Returns (frame, notes)."""
    frames, notes = [], []
    for name in dirs:
        p = lab / name / "radar_gauge_daily.csv"
        if not p.is_file():
            notes.append(f"missing: {p}")
            continue
        d = pd.read_csv(p, index_col=0)
        d.index = pd.to_datetime(d.index).normalize()
        d.index.name = "day"
        frames.append(d)
        notes.append(f"{name}: {len(d)} days, {d.index.min():%b %d %Y} to {d.index.max():%b %d %Y}")
    if not frames:
        return None, notes
    d = pd.concat(frames)
    d = d[~d.index.duplicated(keep="first")].sort_index()
    d = d[(d.index >= pd.Timestamp(WY_FIRST + "-01")) & (d.index < (pd.Period(WY_LAST, "M") + 1).start_time)]
    return d, notes


def gauge_names(d: pd.DataFrame) -> list:
    out = []
    for c in d.columns:
        if c.endswith(" | gauge_mm"):
            g = c[: -len(" | gauge_mm")]
            if f"{g} | radar_3x3_mean_mm" in d.columns and g not in out:
                out.append(g)
    return out


def monthly_table(d: pd.DataFrame, exclude: str = EXCLUDE, min_mm: float = MIN_GAUGE_MM) -> pd.DataFrame:
    """One row per month: radar grid-mean rain, and for 3 gauges and for 2 gauges the mean gauge total, the mean radar 3x3
    total at the same gauges over the same days, and their ratio (blank if the gauges caught under min_mm in all)."""
    names = gauge_names(d)
    two = [g for g in names if exclude not in g]
    rows = []
    for p in pd.period_range(WY_FIRST, WY_LAST, freq="M"):
        m = d[(d.index >= p.start_time) & (d.index < (p + 1).start_time)]
        row = {"month": p.strftime("%b %Y"), "days_with_data": len(m), "days_in_month": p.days_in_month,
               "radar_grid_mean_mm": float(m["grid_mean_radar_mm"].sum()) if len(m) else NAN}
        for tag, group in (("3g", names), ("2g", two)):
            gs, rs = [], []
            for g in group:
                pair = m[[f"{g} | gauge_mm", f"{g} | radar_3x3_mean_mm"]].dropna()
                if len(pair):
                    gs.append(float(pair.iloc[:, 0].sum()))
                    rs.append(float(pair.iloc[:, 1].sum()))
            row[f"gauge_mm_{tag}"] = float(np.mean(gs)) if gs else NAN
            row[f"radar3x3_mm_{tag}"] = float(np.mean(rs)) if rs else NAN
            row[f"ratio_{tag}"] = (sum(rs) / sum(gs)) if (gs and sum(gs) >= min_mm) else NAN
        rows.append(row)
    return pd.DataFrame(rows)


def daily_two_gauge_mean(d: pd.DataFrame, exclude: str = EXCLUDE) -> pd.Series:
    two = [g for g in gauge_names(d) if exclude not in g]
    return d[[f"{g} | gauge_mm" for g in two]].mean(axis=1)


# --------------------------------------------------------------------------------------------
# The figure
# --------------------------------------------------------------------------------------------
def make_figure(out: Path, vol_series: list, tab: pd.DataFrame, d: pd.DataFrame, exclude: str = EXCLUDE) -> str:
    """vol_series: list of (label, color, [12 monthly volumes in M m3]) in the order they are drawn (the gauge first)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    periods = list(pd.period_range(WY_FIRST, WY_LAST, freq="M"))
    starts = [p.start_time for p in periods]
    mids = [s + pd.Timedelta(days=(p.days_in_month - 1) / 2.0) for s, p in zip(starts, periods)]
    s_num = mdates.date2num(pd.DatetimeIndex(starts).to_pydatetime())
    m_num = mdates.date2num(pd.DatetimeIndex(mids).to_pydatetime())
    end_num = mdates.date2num((periods[-1] + 1).start_time.to_pydatetime())

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(11.5, 11), sharex=True, gridspec_kw={"height_ratios": [3, 2, 2.1], "hspace": 0.1})

    # ---- top: monthly volumes ------------------------------------------------------------------------
    n = len(vol_series)
    bw = 4.6
    offs = (np.arange(n) - (n - 1) / 2.0) * bw
    top = 0.0
    for i, (lab_, col, vals) in enumerate(vol_series):
        v = np.asarray(vals, dtype=float)
        ax1.bar(m_num + offs[i], v, width=bw, color=col, label=lab_)
        if np.isfinite(v).any():
            top = max(top, float(np.nanmax(v)))
    ax1.set_ylim(0, max(top, 1e-6) * 1.4)
    ax1.set_ylabel("Monthly outlet volume (M m$^3$)")
    ax1.set_title("Whole water year 2017: outlet volume, radar rain, and radar / gauge ratio", loc="left")
    ax1.legend(frameon=False, ncol=2, loc="upper left", fontsize=8.5)

    # ---- middle: hyetograph --------------------------------------------------------------------------
    x = mdates.date2num(d.index.to_pydatetime())
    rain = d["grid_mean_radar_mm"].fillna(0.0).to_numpy(dtype=float)
    ax2.bar(x, rain, width=1.0, align="edge", color=RAIN_BLUE, label="radar rain, grid mean (all cells)")
    g2 = daily_two_gauge_mean(d, exclude)
    wet = (g2 >= 0.1).to_numpy()          # skip trace amounts so a dry month is not a line of dots
    ax2.plot(x[wet] + 0.5, g2.to_numpy()[wet], "o", ms=2.8, color=INK, label=f"gauges, mean of 2 (without {exclude})")
    ymax = max(float(np.nanmax(rain)) if len(rain) else 1.0, float(np.nanmax(g2.to_numpy())) if g2.notna().any() else 1.0, 1.0) * 1.12
    ax2.set_ylim(ymax, 0)                      # upside down: bars hang from the top
    ax2.set_ylabel("Rain (mm per day)")
    ax2.legend(frameon=False, ncol=2, loc="lower right", fontsize=8.5)

    # ---- bottom: monthly radar / gauge ratio ---------------------------------------------------------
    r3, r2 = tab["ratio_3g"].to_numpy(dtype=float), tab["ratio_2g"].to_numpy(dtype=float)
    ax3.axhline(1.0, color=INK_2, lw=1, ls=(0, (4, 3)))
    ax3.plot(m_num, r3, "-o", color=RATIO3_COLOR, lw=1.6, ms=5, label="3 gauges (with Mt. Ord)")
    ax3.plot(m_num, r2, "--s", color=RATIO2_COLOR, lw=1.6, ms=5, label=f"2 gauges (without {exclude})")
    for xi, a, b in zip(m_num, r3, r2):
        a_above = not (np.isfinite(a) and np.isfinite(b) and a < b)      # the higher of the two gets its number above, the lower below
        if np.isfinite(a):
            ax3.annotate(f"{a:.2f}", (xi, a), xytext=(0, 7 if a_above else -13), textcoords="offset points", ha="center", fontsize=7.5, color=RATIO3_COLOR)
        if np.isfinite(b):
            ax3.annotate(f"{b:.2f}", (xi, b), xytext=(0, -13 if a_above and np.isfinite(a) else 7), textcoords="offset points", ha="center", fontsize=7.5, color=RATIO2_COLOR)
    finite = np.concatenate([r3[np.isfinite(r3)], r2[np.isfinite(r2)]])
    ax3.set_ylim(0, max(2.0, float(finite.max()) * 1.2) if finite.size else 2.0)
    ax3.set_ylabel("Radar / gauge rain, by month\n(1.0 = they agree)")
    ax3.legend(frameon=False, ncol=2, loc="upper left", fontsize=8.5)

    # ---- shared x axis -------------------------------------------------------------------------------
    for ax in (ax1, ax2, ax3):
        ax.set_xticks(s_num, minor=True)
        ax.grid(True, which="minor", axis="x", alpha=0.25)
        ax.grid(True, which="major", axis="y", alpha=0.25)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    ax3.set_xticks(m_num)
    ax3.set_xticklabels([p.strftime("%b") for p in periods])
    ax3.tick_params(axis="x", which="major", length=0)
    ax3.set_xlim(s_num[0] - 2, end_num + 2)
    foot = ("Ratio = radar 3x3-cell mean at each gauge / gauge rain, over the same days (summed over the gauges). A month with under "
            f"{MIN_GAUGE_MM:g} mm of gauge rain has no ratio.\nHyetograph bars use the whole radar grid, so they are not the same numbers as the ratio.")
    fig.text(0.01, 0.005, foot, fontsize=8, color=INK_2, va="bottom", ha="left")
    fig.subplots_adjust(left=0.09, right=0.985, top=0.95, bottom=0.085)
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return f"wrote {out}"


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def print_table(tab: pd.DataFrame, notes: list) -> None:
    log("RAIN BY MONTH (mm). Ratio = radar 3x3 mean at the gauges / gauge rain, same days.")
    for n in notes:
        log("  " + n)
    show = tab.rename(columns={"days_with_data": "days", "days_in_month": "of", "radar_grid_mean_mm": "radar grid mm",
                               "gauge_mm_3g": "gauge mm (3)", "radar3x3_mm_3g": "radar mm (3)", "ratio_3g": "ratio 3 gauges",
                               "gauge_mm_2g": "gauge mm (2)", "radar3x3_mm_2g": "radar mm (2)", "ratio_2g": "ratio 2 gauges"})
    log(show.to_string(index=False, float_format=lambda x: f"{x:.2f}", na_rep="-"))
    short = tab[tab["days_with_data"] < tab["days_in_month"]]
    for r in short.itertuples():
        log(f"  WARNING: {r.month}: only {r.days_with_data} of {r.days_in_month} days have rain data; its totals are low by the missing days.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Monthly volume, daily radar rain (hyetograph) and radar/gauge ratio for water year 2017 (see the top of this file).")
    ap.add_argument("--rain-dirs", default=",".join(DEFAULT_RAIN_DIRS), help="folders holding radar_gauge_daily.csv, Oct-Mar first then Apr-Sep (default %(default)s)")
    ap.add_argument("--out", default=None, help=f"figure file (default {FIG_DIR}/{OUT_NAME})")
    ap.add_argument("--exclude", default=EXCLUDE, help="text in the name of the gauge left out of the '2 gauges' ratio (default %(default)s)")
    ap.add_argument("--min-gauge-mm", type=float, default=MIN_GAUGE_MM)
    ap.add_argument("--table-only", action="store_true", help="print the monthly rain table and write the csv, but make no figure (needs no model output)")
    a = ap.parse_args(argv)
    os.chdir(HERE)

    dirs = [x.strip() for x in a.rain_dirs.split(",") if x.strip()]
    d, notes = load_daily(HERE, dirs)
    if d is None:
        log("No radar_gauge_daily.csv found in: " + ", ".join(dirs))
        for n in notes:
            log("  " + n)
        return 2
    if not gauge_names(d):
        log("The daily tables have no '<gauge> | gauge_mm' columns, so there is nothing to compare. Columns found: " + ", ".join(d.columns[:8]))
        return 2
    tab = monthly_table(d, a.exclude, a.min_gauge_mm)
    print_table(tab, notes)
    tab.to_csv(HERE / MONTHS_CSV, index=False)
    log(f"  Wrote {MONTHS_CSV}.")
    if a.table_only:
        return 0

    import run_w5_kchan_sycamore as K              # brings the tested W4/W5 scripts; importing builds and runs nothing
    T, W = K.T, K.W
    K.set_plan(K.DEFAULT_KCHAN)
    gauge_q = T.read_obs(HERE / T.INIT_DATA / T.OBS_FILE)
    runs = [("first W5 run: re-tuned, loss off (Ks x1.03, f x0.15, As 300)", RETUNED_COLOR, W.W5_ID)]
    for i, j in enumerate(K.planned()):
        runs.append((f"as drawn (f x0.5874, As 1): {K.label_of(j)}", RUN_COLORS[i % len(RUN_COLORS)], j["run_id"]))
    vol_series, gauge_vol = [], None
    for label, color, rid in runs:
        qf = HERE / "results" / f"{rid}_Outlet.qout"
        if not qf.is_file():
            log(f"  skipped (outlet file not found): {qf.name}")
            continue
        rows = W.month_rows(T.read_qout(qf, T.RUN_START_TS), gauge_q)
        if gauge_vol is None:
            gauge_vol = [r["gauge_Mm3"] for r in rows]
        vol_series.append((label, color, [r["model_Mm3"] for r in rows]))
    if not vol_series:
        log("No finished W5 outlet files were found in results/, so no figure.")
        return 2
    vol_series.insert(0, ("USGS gauge", GAUGE_GRAY, gauge_vol))
    out = Path(a.out) if a.out else HERE / FIG_DIR / OUT_NAME
    log("  Figure: " + make_figure(out, vol_series, tab, d, a.exclude))
    return 0


if __name__ == "__main__":
    sys.exit(main())
