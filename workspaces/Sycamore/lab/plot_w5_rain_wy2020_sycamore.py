#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_w5_rain_wy2020_sycamore.py   (v1.1, 2026-11-04)
====================================================

The water year 2020 version of plot_w5_rain_sycamore.py. Same three panels, same look, for the scored window
(October 2019 - May 2020). September 2019 is model spin-up and has no gauge flow, so it is not drawn.

  TOP     monthly outlet volume (M m3): USGS gauge, the Re-tuned run, and the Original run (both with channel loss off),
          plus (v1.1) the Original with channel loss on at Kchan 0.5 and 0.2 mm/hr when their outlet files are found.
  MIDDLE  a hyetograph: daily radar rain (grid mean, mm/day) hanging down from the top, with the daily mean of the two
          good rain gauges (Sycamore Upper and Iron Dike) as black dots.
  BOTTOM  the radar / gauge ratio for each month, two ways: all 3 gauges, and 2 gauges (Mt. Ord left out).
          1.0 means the radar and the gauges agree. Above 1.0 the radar read higher than the gauges.

Both model settings were frozen before WY2020 was run. Nothing here changes a setting; it only draws what the runs gave.

WORDS
  hyetograph     a bar chart of rain against time, drawn upside down (bars hang from the top) so it sits above a flow plot.
  ratio          radar rain divided by gauge rain over the same days (radar = mean of the 3 x 3 cells around each gauge).
  grid mean      the average over every cell of the 32 x 49 radar grid.

WHAT IT READS (nothing is changed, nothing is run)
  results/Sycamore_WY2020_Retuned_Outlet.qout, results/Sycamore_WY2020_Original_Outlet.qout   the finished WY2020 runs
  ../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv                                       the USGS gauge file
  radar_check_out_wy2020/radar_gauge_daily.csv                                                 daily rain table from the WY2020 radar check
WHAT IT WRITES
  w5_kchan_figs/fig_w5_rain_wy2020.png     the figure
  w5_rain_ratio_months_wy2020.csv          the monthly rain numbers behind the bottom panel
  wy2020_monthly_volumes.csv               the monthly volumes behind the top panel

COMMANDS (run from lab/)
  python plot_w5_rain_wy2020_sycamore.py --table-only    # print the monthly tables only (no figure)
  python plot_w5_rain_wy2020_sycamore.py                 # tables plus figure
Options: --rain-dirs a,b  --obs path  --retuned path  --original path  --out path.png  --exclude "Mt. Ord"  --min-gauge-mm 5
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

DEFAULT_RAIN_DIRS = ["radar_check_out_wy2020"]
DEFAULT_OBS = "../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv"
DEFAULT_RETUNED = "results/Sycamore_WY2020_Retuned_Outlet.qout"
DEFAULT_ORIGINAL = "results/Sycamore_WY2020_Original_Outlet.qout"
DEFAULT_KCHAN05 = "results/Sycamore_WY2020_Original_kchan0p5_Outlet.qout"      # v1.1: the two Original runs with channel loss on (optional)
DEFAULT_KCHAN02 = "results/Sycamore_WY2020_Original_kchan0p2_Outlet.qout"
RUN_START = pd.Timestamp("2019-09-01 00:00")      # Time_hr 0 in the outlet files
FIG_DIR = "w5_kchan_figs"
OUT_NAME = "fig_w5_rain_wy2020.png"
MONTHS_CSV = "w5_rain_ratio_months_wy2020.csv"
VOLUMES_CSV = "wy2020_monthly_volumes.csv"
MIN_GAUGE_MM = 5.0
EXCLUDE = "Mt. Ord"
WY_FIRST, WY_LAST = "2019-10", "2020-05"
GAUGE_GRAY, INK, INK_2 = "#8a8984", "#0b0b0b", "#52514e"
RAIN_BLUE = "#8db4e2"
RATIO3_COLOR, RATIO2_COLOR = "#333333", "#c0392b"
RETUNED_COLOR, ORIGINAL_COLOR = "#6b3fa0", "#eb6834"
KCHAN05_COLOR, KCHAN02_COLOR = "#2a78d6", "#4a9d6b"        # same colors as the Kchan 0.5 / 0.2 runs in the WY2017 figures
NAN = float("nan")


def log(msg: str = "") -> None:
    print(msg, flush=True)


def periods():
    return list(pd.period_range(WY_FIRST, WY_LAST, freq="M"))


# --------------------------------------------------------------------------------------------
# Rain tables (same method as the WY2017 script)
# --------------------------------------------------------------------------------------------
def load_daily(lab: Path, dirs: list):
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
    names = gauge_names(d)
    two = [g for g in names if exclude not in g]
    rows = []
    for p in periods():
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
# Flow files and monthly volumes
# --------------------------------------------------------------------------------------------
def read_qout(path: Path) -> pd.Series:
    raw = pd.read_csv(path, header=0).iloc[:, :2]
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    s = pd.Series(raw["Q"].to_numpy(float), index=RUN_START + pd.to_timedelta(raw["Time_hr"].to_numpy(), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_obs(path: Path) -> pd.Series:
    o = pd.read_csv(path, parse_dates=["datetime"]).set_index("datetime")
    s = pd.to_numeric(o["q_cms"], errors="coerce")
    return s[~s.index.duplicated(keep="last")].sort_index()


def month_volumes(q: pd.Series) -> pd.DataFrame:
    """Monthly volume = sum of hourly mean flow x 3600 (M m3), same rule as the storm table. Also the hours counted."""
    rows = []
    for p in periods():
        s = q[(q.index >= p.start_time) & (q.index < (p + 1).start_time)]
        h = s.resample("60min").mean().dropna()
        rows.append({"month": p.strftime("%b %Y"), "Mm3": float(h.sum() * 3600.0 / 1e6) if len(h) else NAN,
                     "hours": int(len(h)), "hours_expected": int(p.days_in_month * 24)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------
# The figure
# --------------------------------------------------------------------------------------------
def make_figure(out: Path, vol_series: list, tab: pd.DataFrame, d: pd.DataFrame, exclude: str = EXCLUDE) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    pers = periods()
    starts = [p.start_time for p in pers]
    mids = [s + pd.Timedelta(days=(p.days_in_month - 1) / 2.0) for s, p in zip(starts, pers)]
    s_num = mdates.date2num(pd.DatetimeIndex(starts).to_pydatetime())
    m_num = mdates.date2num(pd.DatetimeIndex(mids).to_pydatetime())
    end_num = mdates.date2num((pers[-1] + 1).start_time.to_pydatetime())

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(11.5, 11), sharex=True, gridspec_kw={"height_ratios": [3, 2, 2.1], "hspace": 0.1})

    # ---- top: monthly volumes ------------------------------------------------------------------------
    n = len(vol_series)
    bw = min(7.0, 27.0 / n)                      # narrower bars when there are more than three series, so months do not overlap
    offs = (np.arange(n) - (n - 1) / 2.0) * bw
    top = 0.0
    for i, (lab_, col, vals) in enumerate(vol_series):
        v = np.asarray(vals, dtype=float)
        ax1.bar(m_num + offs[i], v, width=bw, color=col, label=lab_)
        if np.isfinite(v).any():
            top = max(top, float(np.nanmax(v)))
    ax1.set_ylim(0, max(top, 1e-6) * 1.4)
    ax1.set_ylabel("Monthly outlet volume (M m$^3$)")
    ax1.set_title("Water year 2020 (Oct 2019 - May 2020): outlet volume, radar rain, and radar / gauge ratio", loc="left")
    ax1.legend(frameon=False, ncol=1, loc="upper right", fontsize=8.5)

    # ---- middle: hyetograph --------------------------------------------------------------------------
    x = mdates.date2num(d.index.to_pydatetime())
    rain = d["grid_mean_radar_mm"].fillna(0.0).to_numpy(dtype=float)
    ax2.bar(x, rain, width=1.0, align="edge", color=RAIN_BLUE, label="radar rain, grid mean (all cells)")
    g2 = daily_two_gauge_mean(d, exclude)
    wet = (g2 >= 0.1).to_numpy()
    ax2.plot(x[wet] + 0.5, g2.to_numpy()[wet], "o", ms=2.8, color=INK, label=f"gauges, mean of 2 (without {exclude})")
    ymax = max(float(np.nanmax(rain)) if len(rain) else 1.0, float(np.nanmax(g2.to_numpy())) if g2.notna().any() else 1.0, 1.0) * 1.12
    ax2.set_ylim(ymax, 0)
    ax2.set_ylabel("Rain (mm per day)")
    ax2.legend(frameon=False, ncol=2, loc="lower right", fontsize=8.5)

    # ---- bottom: monthly radar / gauge ratio ---------------------------------------------------------
    r3, r2 = tab["ratio_3g"].to_numpy(dtype=float), tab["ratio_2g"].to_numpy(dtype=float)
    ax3.axhline(1.0, color=INK_2, lw=1, ls=(0, (4, 3)))
    ax3.plot(m_num, r3, "-o", color=RATIO3_COLOR, lw=1.6, ms=5, label="3 gauges (with Mt. Ord)")
    ax3.plot(m_num, r2, "--s", color=RATIO2_COLOR, lw=1.6, ms=5, label=f"2 gauges (without {exclude})")
    for xi, a, b in zip(m_num, r3, r2):
        a_above = not (np.isfinite(a) and np.isfinite(b) and a < b)
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
    ax3.set_xticklabels([p.strftime("%b") for p in pers])
    ax3.tick_params(axis="x", which="major", length=0)
    ax3.set_xlim(s_num[0] - 2, end_num + 2)
    foot = ("Both model settings were frozen before WY2020 was run (nothing was tuned on this year). Model and gauge volumes = sum of hourly mean flow x 3600.\n"
            "Ratio = radar 3x3-cell mean at each gauge / gauge rain, over the same days (summed over the gauges).\n"
            f"A month with under {MIN_GAUGE_MM:g} mm of gauge rain has no ratio. Hyetograph bars use the whole radar grid, so they are not the same numbers as the ratio.")
    fig.text(0.01, 0.005, foot, fontsize=8, color=INK_2, va="bottom", ha="left")
    fig.subplots_adjust(left=0.09, right=0.985, top=0.95, bottom=0.105)
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return f"wrote {out}"


# --------------------------------------------------------------------------------------------
def print_rain_table(tab: pd.DataFrame, notes: list) -> None:
    log("RAIN BY MONTH (mm). Ratio = radar 3x3 mean at the gauges / gauge rain, same days.")
    for n in notes:
        log("  " + n)
    show = tab.rename(columns={"days_with_data": "days", "days_in_month": "of", "radar_grid_mean_mm": "radar grid mm",
                               "gauge_mm_3g": "gauge mm (3)", "radar3x3_mm_3g": "radar mm (3)", "ratio_3g": "ratio 3 gauges",
                               "gauge_mm_2g": "gauge mm (2)", "radar3x3_mm_2g": "radar mm (2)", "ratio_2g": "ratio 2 gauges"})
    log(show.to_string(index=False, float_format=lambda x: f"{x:.2f}", na_rep="-"))
    for r in tab[tab["days_with_data"] < tab["days_in_month"]].itertuples():
        log(f"  WARNING: {r.month}: only {r.days_with_data} of {r.days_in_month} days have rain data; its totals are low by the missing days.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="WY2020 monthly volume, daily radar rain and radar/gauge ratio (see the top of this file).")
    ap.add_argument("--rain-dirs", default=",".join(DEFAULT_RAIN_DIRS), help="folder(s) holding radar_gauge_daily.csv (default %(default)s)")
    ap.add_argument("--obs", default=DEFAULT_OBS, help="USGS gauge file (default %(default)s)")
    ap.add_argument("--retuned", default=DEFAULT_RETUNED, help="Re-tuned outlet file (default %(default)s)")
    ap.add_argument("--original", default=DEFAULT_ORIGINAL, help="Original outlet file (default %(default)s)")
    ap.add_argument("--kchan05", default=DEFAULT_KCHAN05, help="Original with Kchan 0.5 outlet file; left out if not found (default %(default)s)")
    ap.add_argument("--kchan02", default=DEFAULT_KCHAN02, help="Original with Kchan 0.2 outlet file; left out if not found (default %(default)s)")
    ap.add_argument("--out", default=None, help=f"figure file (default {FIG_DIR}/{OUT_NAME})")
    ap.add_argument("--exclude", default=EXCLUDE, help="text in the name of the gauge left out of the '2 gauges' ratio (default %(default)s)")
    ap.add_argument("--min-gauge-mm", type=float, default=MIN_GAUGE_MM)
    ap.add_argument("--table-only", action="store_true", help="print the monthly tables and write the csv files, but make no figure")
    a = ap.parse_args(argv)
    os.chdir(HERE)

    dirs = [x.strip() for x in a.rain_dirs.split(",") if x.strip()]
    d, notes = load_daily(HERE, dirs)
    if d is None:
        log("No radar_gauge_daily.csv found in: " + ", ".join(dirs))
        for n in notes:
            log("  " + n)
        log("It comes from the WY2020 radar check (output folder radar_check_out_wy2020). Use --rain-dirs to point at the right folder.")
        return 2
    if not gauge_names(d):
        log("The daily table has no '<gauge> | gauge_mm' columns, so there is nothing to compare. Columns found: " + ", ".join(d.columns[:8]))
        return 2
    tab = monthly_table(d, a.exclude, a.min_gauge_mm)
    print_rain_table(tab, notes)
    tab.to_csv(HERE / MONTHS_CSV, index=False)
    log(f"  Wrote {MONTHS_CSV}.\n")

    missing = [p for p in (a.obs, a.retuned, a.original) if not (HERE / p).is_file()]
    if missing:
        log("Cannot draw the volume panel; file(s) not found (run from lab/, or pass the path):")
        for p in missing:
            log("  " + p)
        return 2
    series = [("USGS gauge", GAUGE_GRAY, month_volumes(read_obs(HERE / a.obs))),
              ("Re-tuned (Ks x1.03, f x0.15, As 300), channel loss off", RETUNED_COLOR, month_volumes(read_qout(HERE / a.retuned))),
              ("Original (Ks x0.187, f x0.5874, As 1), channel loss off", ORIGINAL_COLOR, month_volumes(read_qout(HERE / a.original)))]
    extra = []
    for path_, lab_, col_, name_ in ((a.kchan05, "Original, channel loss on (Kchan 0.5 mm/hr)", KCHAN05_COLOR, "kchan05_Mm3"),
                                     (a.kchan02, "Original, channel loss on (Kchan 0.2 mm/hr)", KCHAN02_COLOR, "kchan02_Mm3")):
        if (HERE / path_).is_file():
            series.append((lab_, col_, month_volumes(read_qout(HERE / path_))))
            extra.append(name_)
        else:
            log(f"  NOTE: {path_} not found; that run is left out.")
    vol = pd.DataFrame({"month": series[0][2]["month"], "gauge_Mm3": series[0][2]["Mm3"], "retuned_Mm3": series[1][2]["Mm3"],
                        "original_Mm3": series[2][2]["Mm3"], **{n_: series[3 + i][2]["Mm3"] for i, n_ in enumerate(extra)}, "gauge_hours": series[0][2]["hours"], "hours_expected": series[0][2]["hours_expected"]})
    log("MONTHLY OUTLET VOLUME (M m3)")
    log(vol.to_string(index=False, float_format=lambda x: f"{x:.3f}", na_rep="-"))
    short = vol[vol["gauge_hours"] < vol["hours_expected"]]
    for r in short.itertuples():
        log(f"  NOTE: {r.month}: gauge has {r.gauge_hours} of {r.hours_expected} hours; its volume is low by the missing hours.")
    vol.to_csv(HERE / VOLUMES_CSV, index=False)
    log(f"  Wrote {VOLUMES_CSV}.")
    if a.table_only:
        return 0
    out = Path(a.out) if a.out else HERE / FIG_DIR / OUT_NAME
    log("  Figure: " + make_figure(out, [(l, c, v["Mm3"].tolist()) for l, c, v in series], tab, d, a.exclude))
    return 0


if __name__ == "__main__":
    sys.exit(main())
