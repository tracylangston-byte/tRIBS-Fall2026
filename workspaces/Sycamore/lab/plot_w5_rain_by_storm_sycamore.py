#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_w5_rain_by_storm_sycamore.py   (v1.0, 2026-11-03)
======================================================

The same three things as plot_w5_rain_sycamore.py (outlet volume, rain, radar / gauge ratio), but cut apart BY STORM:
one column per storm, oldest to newest. Two figures:

  WY2017  (5 winter storms)   rows: VOLUME  /  RAIN  /  RATIO
  WY2020  (7 winter storms)   rows:            RAIN  /  RATIO       (no volume row: WY2020 flow and model results stay SEALED)

THE ROWS
  VOLUME  outlet volume (M m3) over the storm window for the USGS gauge, the first W5 run (re-tuned, loss off) and the three
          "as drawn" runs (loss off, Kchan 0.5, Kchan 0.2). The window runs from midnight at the start of the first rain day to
          3 days after the end of the last rain day (--after-days), because the flow arrives after the rain. WY2017 only.
  RAIN    daily radar rain (mm/day, mean of the whole radar grid) hanging from the top, with the daily mean of the two good gauges
          (Sycamore Upper and Iron Dike) as black dots. The day before and the day after the storm are drawn paler. The rain scale is
          the same in every column of one figure, so a big storm looks big.
  RATIO   radar / gauge rain for the storm days: all 3 gauges (black circle) and the 2 good gauges, Mt. Ord left out (red square).
          The radar value is the mean of the 3 x 3 radar cells around each gauge, summed over the same days as the gauge.
          1.0 means radar and gauges agree; above 1.0 the radar read higher. A storm with under 5 mm of gauge rain has no ratio.

THE STORMS (the proposed list: runs of consecutive days with at least 0.5 mm of grid-mean radar rain, total 20 mm or more).
  Tiers by grid-mean total: very big 60 mm or more, medium 25 to 60, small 20 to 25. The list is PROPOSED, not yet accepted.

WHAT IT READS (nothing is changed, nothing is run)
  WY2017 rain:  radar_check_out/radar_gauge_daily.csv (and radar_check_out_octdec/ if the first is missing days)
  WY2020 rain:  radar_check_out_wy2020/radar_gauge_daily.csv
  WY2017 flow:  results/<run id>_Outlet.qout and the USGS gauge file (through run_tail_screen_sycamore.py)
  WY2020 flow and model results are never opened.
WHAT IT WRITES
  w5_kchan_figs/fig_w5_by_storm_wy2017.png, fig_w5_by_storm_wy2020.png, and w5_rain_by_storm.csv

COMMANDS (run from lab/)
  python plot_w5_rain_by_storm_sycamore.py --table-only   # the storm table only (rain only, no model output needed)
  python plot_w5_rain_by_storm_sycamore.py                # table plus both figures (the 2017 volumes need the W5 outlet files)
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
import plot_w5_rain_sycamore as P                    # noqa: E402  (its colors, gauge-name finder and constants; importing runs nothing)

NAN = float("nan")
FIG_DIR = P.FIG_DIR
CSV_NAME = "w5_rain_by_storm.csv"
AFTER_DAYS = 3
MIN_GAUGE_MM = P.MIN_GAUGE_MM
EXCLUDE = P.EXCLUDE

# (name, first rain day, last rain day), local Arizona days, inclusive
STORMS = {
    "WY2017": [("Dec 21-25", "2016-12-21", "2016-12-25"), ("Dec 30-Jan 2", "2016-12-30", "2017-01-02"),
               ("Jan 19-21", "2017-01-19", "2017-01-21"), ("Feb 18-19", "2017-02-18", "2017-02-19"),
               ("Feb 27-28", "2017-02-27", "2017-02-28")],
    "WY2020": [("Nov 19-22", "2019-11-19", "2019-11-22"), ("Nov 28-30", "2019-11-28", "2019-11-30"),
               ("Dec 23-28", "2019-12-23", "2019-12-28"), ("Feb 10-11", "2020-02-10", "2020-02-11"),
               ("Feb 22-23", "2020-02-22", "2020-02-23"), ("Mar 10-13", "2020-03-10", "2020-03-13"),
               ("Mar 18-19", "2020-03-18", "2020-03-19")],
}
RAIN_DIRS = {"WY2017": ["radar_check_out", "radar_check_out_octdec"], "WY2020": ["radar_check_out_wy2020"]}


def log(msg: str = "") -> None:
    print(msg, flush=True)


def tier(total_mm: float) -> str:
    if not np.isfinite(total_mm):
        return "no data"
    return "very big" if total_mm >= 60 else ("medium" if total_mm >= 25 else ("small" if total_mm >= 20 else "below 20 mm"))


# --------------------------------------------------------------------------------------------
# Rain
# --------------------------------------------------------------------------------------------
def load_daily(lab: Path, dirs: list):
    """Stack radar_gauge_daily.csv files (first folder wins when a day is in two). Returns (frame or None, notes)."""
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
    return d[~d.index.duplicated(keep="first")].sort_index(), notes


def storm_stats(d: pd.DataFrame, name: str, first: str, last: str, exclude: str = EXCLUDE, min_mm: float = MIN_GAUGE_MM) -> dict:
    a, b = pd.Timestamp(first), pd.Timestamp(last)
    days = pd.date_range(a, b, freq="D")
    m = d[(d.index >= a) & (d.index <= b)]
    names = P.gauge_names(d)
    two = [g for g in names if exclude not in g]
    row = {"storm": name, "first_day": a, "last_day": b, "days_expected": len(days), "days_with_data": len(m),
           "grid_total_mm": float(m["grid_mean_radar_mm"].sum()) if len(m) else NAN}
    row["tier"] = tier(row["grid_total_mm"])
    row["wettest_day"] = m["grid_mean_radar_mm"].idxmax() if len(m) and m["grid_mean_radar_mm"].notna().any() else pd.NaT
    row["wettest_day_mm"] = float(m["grid_mean_radar_mm"].max()) if len(m) else NAN
    for tag, group in (("3g", names), ("2g", two)):
        gs, rs = [], []
        for g in group:
            cols = [f"{g} | gauge_mm", f"{g} | radar_3x3_mean_mm"]
            if not all(c in m.columns for c in cols):
                continue
            pair = m[cols].dropna()
            if len(pair):
                gs.append(float(pair.iloc[:, 0].sum()))
                rs.append(float(pair.iloc[:, 1].sum()))
        row[f"gauge_mm_{tag}"] = float(np.mean(gs)) if gs else NAN
        row[f"radar3x3_mm_{tag}"] = float(np.mean(rs)) if rs else NAN
        row[f"ratio_{tag}"] = (sum(rs) / sum(gs)) if (gs and sum(gs) >= min_mm) else NAN
    return row


# --------------------------------------------------------------------------------------------
# Volumes (WY2017 only)
# --------------------------------------------------------------------------------------------
def storm_windows(storms: list, after_days: int = AFTER_DAYS) -> list:
    out = []
    for name, first, last in storms:
        a = pd.Timestamp(first)
        b = pd.Timestamp(last) + pd.Timedelta(days=1 + after_days)
        out.append((name, a, b))
    return out


def compute_volumes(storms: list, run_q: list, gauge_q: pd.Series, hourly, vol_mm3, after_days: int = AFTER_DAYS) -> list:
    """run_q: [(label, color, flow series)]. Returns one dict per storm: gauge volume, hours of gauge record, model volumes."""
    rows = []
    for name, a, b in storm_windows(storms, after_days):
        gh = hourly(gauge_q, a, b)
        row = {"storm": name, "window_start": a, "window_end": b, "hours_expected": int((b - a) / pd.Timedelta(hours=1)),
               "gauge_hours": int(len(gh)), "gauge_Mm3": float(vol_mm3(gh)), "runs": []}
        for label, color, q in run_q:
            row["runs"].append(float(vol_mm3(hourly(q, a, b))))
        rows.append(row)
    return rows


# --------------------------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------------------------
def make_year_figure(out: Path, year: str, stats: list, d: pd.DataFrame, vol_rows: list | None, run_meta: list,
                     exclude: str = EXCLUDE, after_days: int = AFTER_DAYS) -> str:
    """vol_rows None = no volume row (WY2020). run_meta: [(label, color)] for the model runs, in the order of each vol row's 'runs'."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    n = len(stats)
    with_vol = vol_rows is not None
    nrows = 3 if with_vol else 2
    ratios = [2.2, 2.4, 1.7] if with_vol else [2.4, 1.7]
    fig, axes = plt.subplots(nrows, n, figsize=(max(9.0, 2.35 * n + 1.1), 3.1 * nrows + 1.0), squeeze=False,
                             gridspec_kw={"height_ratios": ratios, "hspace": 0.28, "wspace": 0.12})
    r = 0
    names_gauges = P.gauge_names(d)
    two = [g for g in names_gauges if exclude not in g]

    # ---- volume row ----------------------------------------------------------------------------------
    if with_vol:
        vmax = max([v["gauge_Mm3"] if np.isfinite(v["gauge_Mm3"]) else 0.0 for v in vol_rows] + [x for v in vol_rows for x in v["runs"] if np.isfinite(x)] + [0.01])
        series_n = 1 + len(run_meta)
        bw = 0.8 / series_n
        for i, v in enumerate(vol_rows):
            ax = axes[r, i]
            vals = [v["gauge_Mm3"]] + list(v["runs"])
            cols = [P.GAUGE_GRAY] + [c for _, c in run_meta]
            for k, (val, col) in enumerate(zip(vals, cols)):
                ax.bar(-0.4 + bw * (k + 0.5), val, width=bw, color=col)
            ax.set_xlim(-0.55, 0.55)
            ax.set_ylim(0, vmax * 1.12)
            ax.set_xticks([])
            if v["gauge_hours"] < 0.97 * v["hours_expected"]:
                ax.text(0.5, 0.97, f"gauge record: {v['gauge_hours']} of {v['hours_expected']} h\n(volume is a minimum)", transform=ax.transAxes,
                        ha="center", va="top", fontsize=7, color=P.INK_2)
            if i == 0:
                ax.set_ylabel("Outlet volume (M m$^3$)\nover the storm window")
            else:
                ax.tick_params(labelleft=False)
        r += 1

    # ---- rain row ------------------------------------------------------------------------------------
    pads = []
    for s in stats:
        days = pd.date_range(s["first_day"] - pd.Timedelta(days=1), s["last_day"] + pd.Timedelta(days=1), freq="D")
        pads.append(days)
    rain_max = 1.0
    for s, days in zip(stats, pads):
        m = d.reindex(days)
        rain_max = max(rain_max, float(np.nanmax(m["grid_mean_radar_mm"].to_numpy(dtype=float))) if m["grid_mean_radar_mm"].notna().any() else 0.0)
        g2 = m[[f"{g} | gauge_mm" for g in two]].mean(axis=1) if two else pd.Series(dtype=float)
        if len(g2) and g2.notna().any():
            rain_max = max(rain_max, float(np.nanmax(g2.to_numpy(dtype=float))))
    for i, (s, days) in enumerate(zip(stats, pads)):
        ax = axes[r, i]
        m = d.reindex(days)
        x = np.arange(len(days))
        inside = (days >= s["first_day"]) & (days <= s["last_day"])
        bars = m["grid_mean_radar_mm"].fillna(0.0).to_numpy(dtype=float)
        ax.bar(x, bars, width=0.85, color=[P.RAIN_BLUE if k else "#cfe0f3" for k in inside])
        g2 = m[[f"{g} | gauge_mm" for g in two]].mean(axis=1) if two else pd.Series(np.nan, index=days)
        gv = g2.to_numpy(dtype=float)
        ok = np.isfinite(gv) & (gv >= 0.1)
        ax.plot(x[ok], gv[ok], "o", ms=4.2, color=P.INK)
        ax.set_ylim(rain_max * 1.12, 0)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{t.day}" for t in days], fontsize=7.5)
        ax.set_xlabel(f"{days[0]:%b %Y}" if days[0].month == days[-1].month else f"{days[0]:%b} / {days[-1]:%b %Y}", fontsize=8, color=P.INK_2)
        if i == 0:
            ax.set_ylabel("Rain (mm per day)")
        else:
            ax.tick_params(labelleft=False)
        ax.set_title(f"{s['storm']}\n{s['grid_total_mm']:.1f} mm, {s['tier']}", fontsize=9.5, loc="center", color=P.INK)
    # the titles of the first row belong above the first row of the figure: move them there
    top_row = 0
    for i, s in enumerate(stats):
        t = axes[r, i].get_title()
        axes[r, i].set_title("")
        axes[top_row, i].set_title(t, fontsize=9.5, color=P.INK)
    r += 1

    # ---- ratio row -----------------------------------------------------------------------------------
    allr = [s[k] for s in stats for k in ("ratio_3g", "ratio_2g") if np.isfinite(s[k])]
    ymax = max(2.0, 1.25 * max(allr)) if allr else 2.0
    for i, s in enumerate(stats):
        ax = axes[r, i]
        ax.axhline(1.0, color=P.INK_2, lw=1, ls=(0, (4, 3)))
        for xpos, key, mk, col, dy in ((-0.2, "ratio_3g", "o", P.RATIO3_COLOR, 0.0), (0.2, "ratio_2g", "s", P.RATIO2_COLOR, 0.0)):
            v = s[key]
            if np.isfinite(v):
                ax.plot([xpos], [v], mk, ms=8, color=col, mec="white", mew=0.8)
                below = v < 1.0
                ax.annotate(f"{v:.2f}", (xpos, v), xytext=(0, -15 if below else 9), textcoords="offset points", ha="center", fontsize=8.5, color=col)
            else:
                ax.annotate("no ratio", (xpos, 0.05 * ymax), ha="center", fontsize=7, color=P.INK_2, rotation=90, va="bottom")
        ax.set_xlim(-0.6, 0.6)
        ax.set_ylim(0, ymax)
        ax.set_xticks([])
        if i == 0:
            ax.set_ylabel("Radar / gauge rain\n(1.0 = they agree)")
        else:
            ax.tick_params(labelleft=False)
    # ---- shared look and legend ----------------------------------------------------------------------
    for ax in axes.ravel():
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
        ax.grid(True, axis="y", alpha=0.25)
        ax.set_axisbelow(True)
    handles = []
    if with_vol:
        handles += [Patch(color=P.GAUGE_GRAY, label="USGS gauge")] + [Patch(color=c, label=lab) for lab, c in run_meta]
    handles += [Patch(color=P.RAIN_BLUE, label="radar rain, grid mean (all cells)"),
                Line2D([], [], marker="o", ls="", color=P.INK, label=f"gauges, mean of 2 (without {exclude})"),
                Line2D([], [], marker="o", ls="", color=P.RATIO3_COLOR, label="ratio, 3 gauges (with Mt. Ord)"),
                Line2D([], [], marker="s", ls="", color=P.RATIO2_COLOR, label=f"ratio, 2 gauges (without {exclude})")]
    fig.suptitle(f"{year}: storm by storm" + ("" if with_vol else "  (rain only: WY2020 flow and model results stay sealed)"), x=0.01, ha="left", fontsize=12)
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.028), ncol=5 if with_vol else 4, frameon=False, fontsize=8.5)
    foot = ("Ratio = radar 3x3-cell mean at each gauge / gauge rain over the storm days, summed over the gauges (under "
            f"{MIN_GAUGE_MM:g} mm of gauge rain: no ratio).\nPaler bars are the day before and the day after."
            + (f" Volume window: 00:00 on the first rain day to {after_days} days after the last rain day." if with_vol else ""))
    fig.text(0.01, 0.004, foot, fontsize=7.5, color=P.INK_2, va="bottom", ha="left")
    fig.subplots_adjust(left=0.07, right=0.99, top=0.88 if not with_vol else 0.9, bottom=0.20 if with_vol else 0.22)
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return f"wrote {out}"


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def print_table(year: str, stats: list, notes: list) -> None:
    log(f"\n{year}: STORMS (rain, mm). Ratio = radar 3x3 mean at the gauges / gauge rain, over the storm days.")
    for n in notes:
        log("  " + n)
    df = pd.DataFrame(stats)[["storm", "days_with_data", "days_expected", "grid_total_mm", "tier", "wettest_day_mm",
                               "gauge_mm_3g", "radar3x3_mm_3g", "ratio_3g", "gauge_mm_2g", "radar3x3_mm_2g", "ratio_2g"]]
    df = df.rename(columns={"days_with_data": "days", "days_expected": "of", "grid_total_mm": "grid mm", "wettest_day_mm": "wettest day mm",
                            "gauge_mm_3g": "gauge mm (3)", "radar3x3_mm_3g": "radar mm (3)", "ratio_3g": "ratio 3 gauges",
                            "gauge_mm_2g": "gauge mm (2)", "radar3x3_mm_2g": "radar mm (2)", "ratio_2g": "ratio 2 gauges"})
    log(df.to_string(index=False, float_format=lambda x: f"{x:.2f}", na_rep="-"))
    for s in stats:
        if s["days_with_data"] < s["days_expected"]:
            log(f"  WARNING: {s['storm']}: only {s['days_with_data']} of {s['days_expected']} storm days have rain data; its totals are low.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Outlet volume, rain and radar/gauge ratio, one column per storm (see the top of this file).")
    ap.add_argument("--after-days", type=int, default=AFTER_DAYS, help="days after the last rain day included in the volume window (default %(default)s)")
    ap.add_argument("--exclude", default=EXCLUDE, help="text in the name of the gauge left out of the '2 gauges' ratio (default %(default)s)")
    ap.add_argument("--min-gauge-mm", type=float, default=MIN_GAUGE_MM)
    ap.add_argument("--table-only", action="store_true", help="print the storm table and write the csv; make no figure (needs no model output)")
    ap.add_argument("--dirs-wy2017", default=",".join(RAIN_DIRS["WY2017"]), help="WY2017 daily-rain folders, first wins (default %(default)s)")
    ap.add_argument("--dirs-wy2020", default=",".join(RAIN_DIRS["WY2020"]), help="WY2020 daily-rain folder(s) (default %(default)s)")
    a = ap.parse_args(argv)
    os.chdir(HERE)

    data, all_stats, csv_rows = {}, {}, []
    for year, dirs_text in (("WY2017", a.dirs_wy2017), ("WY2020", a.dirs_wy2020)):
        d, notes = load_daily(HERE, [x.strip() for x in dirs_text.split(",") if x.strip()])
        if d is None or not P.gauge_names(d):
            log(f"{year}: no usable radar_gauge_daily.csv in {dirs_text}")
            for n in notes:
                log("  " + n)
            return 2
        stats = [storm_stats(d, nm, f, l, a.exclude, a.min_gauge_mm) for nm, f, l in STORMS[year]]
        print_table(year, stats, notes)
        data[year], all_stats[year] = d, stats
        for s in stats:
            csv_rows.append({"year": year, **s})

    vols = None
    run_meta = []
    if not a.table_only:
        import run_w5_kchan_sycamore as K                 # brings the tested W4/W5 scripts; importing builds and runs nothing
        T, W = K.T, K.W
        K.set_plan(K.DEFAULT_KCHAN)
        gauge_q = T.read_obs(HERE / T.INIT_DATA / T.OBS_FILE)
        runs = [("first W5 run: re-tuned, loss off", P.RETUNED_COLOR, W.W5_ID)]
        for i, j in enumerate(K.planned()):
            runs.append((f"as drawn: {K.label_of(j)}", P.RUN_COLORS[i % len(P.RUN_COLORS)], j["run_id"]))
        run_q = []
        for label, color, rid in runs:
            qf = HERE / "results" / f"{rid}_Outlet.qout"
            if not qf.is_file():
                log(f"  skipped (outlet file not found): {qf.name}")
                continue
            run_q.append((label, color, T.read_qout(qf, T.RUN_START_TS)))
        if not run_q:
            log("No finished W5 outlet files in results/, so the WY2017 volume row cannot be drawn.")
            return 2
        vols = compute_volumes(STORMS["WY2017"], run_q, gauge_q, W.hourly, W.vol_mm3, a.after_days)
        run_meta = [(lab, col) for lab, col, _ in run_q]
        log("\nWY2017: OUTLET VOLUME per storm window, M m3 (window = first rain day 00:00 to %d days after the last rain day)" % a.after_days)
        vt = pd.DataFrame([{"storm": v["storm"], "gauge hours": f"{v['gauge_hours']} of {v['hours_expected']}", "USGS gauge": v["gauge_Mm3"],
                            **{lab: x for (lab, _), x in zip(run_meta, v["runs"])}} for v in vols])
        log(vt.to_string(index=False, float_format=lambda x: f"{x:.3f}", na_rep="-"))
        for row, v in zip([r for r in csv_rows if r["year"] == "WY2017"], vols):
            row["gauge_Mm3"] = v["gauge_Mm3"]
            for (lab, _), x in zip(run_meta, v["runs"]):
                row[f"model_Mm3 | {lab}"] = x

    pd.DataFrame(csv_rows).to_csv(HERE / CSV_NAME, index=False)
    log(f"\n  Wrote {CSV_NAME}.")
    if a.table_only:
        return 0
    for year in ("WY2017", "WY2020"):
        out = HERE / FIG_DIR / f"fig_w5_by_storm_{year.lower()}.png"
        log("  Figure: " + make_year_figure(out, year, all_stats[year], data[year], vols if year == "WY2017" else None, run_meta, a.exclude, a.after_days))
    return 0


if __name__ == "__main__":
    sys.exit(main())
