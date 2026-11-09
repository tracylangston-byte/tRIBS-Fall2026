#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_storm_by_storm_all_sycamore.py   (v1.3, 2026-11-04)
=========================================================

The "storm by storm" figure for EVERY storm in the storm table (storm_table.csv), not only the 20 mm list.
It replaces plot_w5_rain_by_storm_sycamore.py's two figures (5 and 7 storms). Same three rows:

  VOLUME                 outlet volume in M m3 for the USGS gauge and the model runs. WY2017: Re-tuned and the three Original runs
                         (loss off, Kchan 0.5, Kchan 0.2). WY2020: Re-tuned and Original (loss off) and, new in v1.3, the Original with Kchan 0.5 and Kchan 0.2 (same colors as WY2017). The window runs from 00:00 on the first rain day to
                         3 days after the last rain day, or the next storm's first day if that comes sooner (the
                         table's "event window"; no flow is counted in two panels). Above each panel: gauge volume / Re-tuned volume (M m3);
                         WY2020 also shows the Original volume (gauge / Re-tuned / Original).
                         A star (*) replaces the gauge bar and number when the gauge recorded NO flow in the window
                         (volume under 0.005 M m3). WY2020 September storms (before the USGS file starts on Oct 1) have no gauge record: the panel says
                         "no gauge record" and shows the model only. They are spin-up and do not set the scale of the other panels.
  RAIN                   daily radar rain (mm per day, mean of the whole radar grid), from the hourly basin rain files.
  RATIO                  radar / gauge rain over the storm days (3 gauges, and 2 gauges without Mt. Ord), taken from the storm
                         table. The table has them only for the storms on the 20 mm list; other storms show "not computed".
                         (No gauge dots in this version: they need the daily gauge files, which were not available.)

The WY2020 volume row was added after the WY2020 scoring plan was logged and the runs were scored (nothing here is tuned on WY2020).
A year with more than 8 storms is cut into two bands of panels (same scales in both bands).

NAMES (log entry 21:10): Re-tuned = Ks x1.03, f x0.15, As 300 (the frozen setting); Original = the "as drawn" setting.

WHAT IT READS  (nothing is changed, nothing is run)
  storm_table.csv                                    the storm list (year, event, first_day, last_day, grid_mm, status)
  hourly_basin_wy2017.csv, hourly_basin_wy2020.csv   hourly radar rain (grid_mean_mm), summed to local days
  results/<run>_Outlet.qout (4 WY2017 files, 2 WY2020 files) + the two USGS files
WHAT IT WRITES
  <out-dir>/fig_storm_by_storm_all_wy2017.png, fig_storm_by_storm_all_wy2020.png, storm_by_storm_all.csv

COMMAND (from lab/):   python plot_storm_by_storm_all_sycamore.py
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plot_w5_rain_sycamore as P                          # noqa: E402  colors and gauge-name finder

AFTER_DAYS = 3
NO_FLOW_MM3 = 0.005                                        # gauge volume below this = "no gauge flow"
RUN_START = {"WY2017": "2016-10-01 00:00", "WY2020": "2019-09-01 00:00"}      # Time_hr 0 of each run
RUNS_BY_YEAR = {  # label, color, outlet file stem
    "WY2017": [
        ("Re-tuned: loss off", P.RETUNED_COLOR, "W5_d012_k10_f0p15_as300_ksx1p03_nocl"),
        ("Original: loss off", P.RUN_COLORS[0], "W5_d012_k10_base_nocl"),
        ("Original: Kchan 0.5", P.RUN_COLORS[1], "W5_d012_k10_base_kchan0p5"),
        ("Original: Kchan 0.2", P.RUN_COLORS[2], "W5_d012_k10_base_kchan0p2"),
    ],
    "WY2020": [
        ("Re-tuned: loss off", P.RETUNED_COLOR, "Sycamore_WY2020_Retuned"),
        ("Original: loss off", P.RUN_COLORS[0], "Sycamore_WY2020_Original"),
        ("Original: Kchan 0.5", P.RUN_COLORS[1], "Sycamore_WY2020_Original_kchan0p5"),      # v1.3: the two channel-loss comparison runs
        ("Original: Kchan 0.2", P.RUN_COLORS[2], "Sycamore_WY2020_Original_kchan0p2"),
    ],
}


def log(msg=""):
    print(msg, flush=True)


# ------------------------------------------------------------------------------------------------ data
def load_daily(path: Path) -> pd.DataFrame:
    """Hourly basin rain file -> daily grid-mean rain (mm per day). Hours are labeled by the START of the hour (Arizona time)."""
    h = pd.read_csv(path, parse_dates=["interval_start_mst"])
    day = h["interval_start_mst"].dt.normalize()
    d = h.groupby(day)["grid_mean_mm"].sum().to_frame("grid_mean_radar_mm")
    d["hours"] = h.groupby(day)["grid_mean_mm"].size()
    d.index.name = "day"
    log(f"  {path.name}: {len(d)} days, {d.index.min():%b %d %Y} to {d.index.max():%b %d %Y}")
    return d


def read_qout(path: Path, start: str) -> pd.Series:
    raw = pd.read_csv(path, sep=",", header=0).iloc[:, :2]
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    s = pd.Series(raw["Q"].to_numpy(float), index=pd.Timestamp(start) + pd.to_timedelta(raw["Time_hr"].to_numpy(), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_obs(path: Path) -> pd.Series:
    o = pd.read_csv(path, parse_dates=["datetime"]).set_index("datetime")
    s = pd.to_numeric(o["q_cms"], errors="coerce")
    return s[~s.index.duplicated(keep="last")].sort_index()


def volume(q: pd.Series, a, b):
    s = q[(q.index >= a) & (q.index < b)]
    h = s.resample("60min").mean().dropna()
    return (float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")), int(len(h))


# ------------------------------------------------------------------------------------------------ figure
def make_year_figure(out: Path, year: str, rows: list, d: pd.DataFrame, vols, runs, exclude="Mt. Ord"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    n = len(rows)
    per = n if n <= 8 else math.ceil(n / 2)
    bands = [list(range(i, min(i + per, n))) for i in range(0, n, per)]
    with_vol = vols is not None
    nr = 3 if with_vol else 2
    ratios = [2.0, 2.4, 1.6] if with_vol else [2.4, 1.6]
    band_h = 6.2 if with_vol else 4.0
    fig = plt.figure(figsize=(2.15 * per + 1.2, band_h * len(bands) + 2.0))
    hr = []
    for b_ in range(len(bands)):
        hr += ratios + ([0.6] if b_ < len(bands) - 1 else [])           # a spacer row between bands
    gs = fig.add_gridspec(len(hr), per, height_ratios=hr, left=0.07, right=0.99, top=0.94, bottom=0.115 if with_vol else 0.17,
                          hspace=0.5, wspace=0.12)

    # shared scales (same in every band)
    pads = {i: pd.date_range(rows[i]["first"] - pd.Timedelta(days=1), rows[i]["last"] + pd.Timedelta(days=1), freq="D") for i in range(n)}
    rain_max = 1.0
    for i in range(n):
        m = d.reindex(pads[i])
        rain_max = max(rain_max, float(np.nanmax(m["grid_mean_radar_mm"].to_numpy(float))) if m["grid_mean_radar_mm"].notna().any() else 0)
    vmax = 0.01
    if with_vol:
        scored = [v for v, st in zip(vols, rows) if st["status_short"] != "spin-up"] or vols      # September spin-up storms do not set the scale
        vmax = max([v["gauge"] for v in scored if np.isfinite(v["gauge"])] + [x for v in scored for x in v["runs"] if np.isfinite(x)] + [0.01])
    ratio_vals = [r[k] for r in rows for k in ("ratio_3g", "ratio_2g") if np.isfinite(r[k])]
    ymax_r = max(2.0, 1.25 * max(ratio_vals)) if ratio_vals else 2.0

    for b, idx in enumerate(bands):
        off = b * (nr + 1)
        axes = np.array([[fig.add_subplot(gs[off + r_, c_]) for c_ in range(per)] for r_ in range(nr)])
        for c in range(per):
            if c >= len(idx):
                for r_ in range(nr):
                    axes[r_, c].axis("off")
                continue
            i = idx[c]
            s = rows[i]
            r = 0
            if with_vol:
                ax = axes[r, c]
                v = vols[i]
                vals = [v["gauge"]] + list(v["runs"])
                cols = [P.GAUGE_GRAY] + [col for _, col, _ in runs]
                bw = 0.8 / len(vals)
                for k, (val, col) in enumerate(zip(vals, cols)):
                    if k == 0 and v["no_record"]:
                        continue                                                           # no gauge file for these days: draw nothing
                    if k == 0 and v["no_flow"]:
                        ax.text(-0.4 + bw * 0.5, vmax * 0.015, "*", ha="center", va="bottom", fontsize=17, color=P.INK, fontweight="bold")
                        continue
                    ax.bar(-0.4 + bw * (k + 0.5), val, width=bw, color=col)
                ax.set_xlim(-0.55, 0.55)
                ax.set_ylim(0, vmax * 1.28)
                if any(np.isfinite(x) and x > vmax * 1.28 for x in vals):
                    ax.text(0.5, 0.74, "bars off scale", transform=ax.transAxes, ha="center", va="top", fontsize=6.3, color=P.INK_2, bbox=dict(facecolor="white", edgecolor="none", pad=1.2), zorder=5)
                ax.set_xticks([])
                gtxt = "n/a" if v["no_record"] else ("*" if v["no_flow"] else f"{v['gauge']:.2f}")
                if year == "WY2020" and len(runs) >= 4:                                   # line 1: gauge / Re-tuned / Original; line 2: Kchan 0.5 / Kchan 0.2
                    txt = " / ".join([gtxt] + [f"{x:.2f}" for x in v["runs"][:2]]) + "\n" + " / ".join(f"{x:.2f}" for x in v["runs"][2:4])
                else:
                    shown = v["runs"][:2] if len(runs) <= 2 else v["runs"][:1]           # WY2017: Re-tuned only
                    txt = " / ".join([gtxt] + [f"{x:.2f}" for x in shown])
                ax.text(0.5, 0.985, txt, transform=ax.transAxes, ha="center", va="top", fontsize=7.8, color=P.INK, bbox=dict(facecolor="white", edgecolor="none", pad=1.2), zorder=5)
                if v["no_record"]:
                    ax.text(0.5, 0.86, "no gauge record", transform=ax.transAxes, ha="center", va="top", fontsize=6.3, color=P.INK_2, bbox=dict(facecolor="white", edgecolor="none", pad=1.2), zorder=5)
                elif v["gauge_hours"] < 0.97 * v["hours_expected"]:
                    ax.text(0.5, 0.86, f"gauge {v['gauge_hours']} of {v['hours_expected']} h", transform=ax.transAxes, ha="center", va="top",
                            fontsize=6.3, color=P.INK_2, bbox=dict(facecolor="white", edgecolor="none", pad=1.2), zorder=5)
                if c == 0:
                    ax.set_ylabel("Outlet volume (M m$^3$)\nover the storm window")
                else:
                    ax.tick_params(labelleft=False)
                r += 1
            # rain
            ax = axes[r, c]
            days = pads[i]
            m = d.reindex(days)
            x = np.arange(len(days))
            inside = (days >= s["first"]) & (days <= s["last"])
            ax.bar(x, m["grid_mean_radar_mm"].fillna(0.0).to_numpy(float), width=0.85, color=[P.RAIN_BLUE if k else "#cfe0f3" for k in inside])
            ax.set_ylim(rain_max * 1.12, 0)
            ax.set_xticks(x)
            ax.set_xticklabels([f"{t.day}" for t in days], fontsize=7)
            ax.set_xlabel(f"{days[0]:%b %Y}" if days[0].month == days[-1].month else f"{days[0]:%b} / {days[-1]:%b %Y}", fontsize=7.5, color=P.INK_2)
            if c == 0:
                ax.set_ylabel("Rain (mm per day)")
            else:
                ax.tick_params(labelleft=False)
            axes[0, c].set_title(f"{s['event']}\n{s['grid_mm']:.1f} mm, {s['status_short']}", fontsize=9, color=P.INK)
            r += 1
            # ratio
            ax = axes[r, c]
            ax.axhline(1.0, color=P.INK_2, lw=1, ls=(0, (4, 3)))
            for xpos, key, mk, col in ((-0.2, "ratio_3g", "o", P.RATIO3_COLOR), (0.2, "ratio_2g", "s", P.RATIO2_COLOR)):
                val = s[key]
                if np.isfinite(val):
                    ax.plot([xpos], [val], mk, ms=7, color=col, mec="white", mew=0.8)
                    ax.annotate(f"{val:.2f}", (xpos, val), xytext=(0, -14 if val < 1 else 8), textcoords="offset points", ha="center", fontsize=8, color=col)
                else:
                    if xpos < 0:
                        ax.text(0.5, 0.5, "not computed", transform=ax.transAxes, ha="center", va="center", fontsize=6.8, color=P.INK_2, rotation=90)
            ax.set_xlim(-0.6, 0.6)
            ax.set_ylim(0, ymax_r)
            ax.set_xticks([])
            if c == 0:
                ax.set_ylabel("Radar / gauge rain\n(1.0 = they agree)")
            else:
                ax.tick_params(labelleft=False)
        for ax in axes.ravel():
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            ax.grid(True, axis="y", alpha=0.25)
            ax.set_axisbelow(True)

    handles = []
    if with_vol:
        handles += [Patch(color=P.GAUGE_GRAY, label="USGS gauge")] + [Patch(color=c, label=l) for l, c, _ in runs]
    handles += [Patch(color=P.RAIN_BLUE, label="radar rain, grid mean (all cells)"),
                Line2D([], [], marker="o", ls="", color=P.RATIO3_COLOR, label="ratio, 3 gauges (with Mt. Ord)"),
                Line2D([], [], marker="s", ls="", color=P.RATIO2_COLOR, label=f"ratio, 2 gauges (without {exclude})")]
    fig.suptitle(f"{year}: storm by storm, all {n} storms in the storm table" + ("" if with_vol else "  (rain only)"),
                 x=0.01, y=0.992, ha="left", fontsize=12)
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.035), ncol=5 if with_vol else 4, frameon=False, fontsize=8.3)
    foot = ("Ratio = radar 3x3-cell mean at each gauge / gauge rain over the storm days, summed over the gauges; the storm table has it only for the 20 mm list. "
            "Paler bars = day before and day after.\n")
    if with_vol:
        foot += (f"Volume window: 00:00 on the first rain day to {AFTER_DAYS} days after the last rain day, or the next storm's first day if sooner.\n"
                 f"Numbers above the bars: gauge / Re-tuned" + (" / Original (WY2020: second line = Original with Kchan 0.5 / Kchan 0.2)" if year == "WY2020" else "") + " volume (M m3). "
                 f"* = USGS gauge recorded no flow in the window (under {NO_FLOW_MM3:g} M m3).")
    fig.text(0.01, 0.004, foot, fontsize=7.2, color=P.INK_2, va="bottom", ha="left")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return str(out)


# ------------------------------------------------------------------------------------------------ the storm table
# The table can come in three shapes: my storm_table.csv (short names such as "event"), a CSV made from the spreadsheet (readable
# headers such as "Grid rain (mm)"), or the spreadsheet itself (Sycamore_Storm_Table.xlsx, sheet "Storm table", headers on row 5).
HEADER_TO_KEY = {"water year": "year", "event": "event", "first rain day": "first_day", "last rain day": "last_day",
                 "list status": "status", "grid rain": "grid_mm", "radar / gauge (3)": "ratio_3g", "radar / gauge (2)": "ratio_2g"}


def read_table(path: str) -> pd.DataFrame:
    if str(path).lower().endswith((".xlsx", ".xlsm")):
        t = pd.read_excel(path, sheet_name="Storm table", header=4)
    else:
        t = pd.read_csv(path)
        if "event" not in t.columns and "Event" not in t.columns:                 # the header may sit lower down; look for it
            raise SystemExit(f"{path}: no 'event' or 'Event' column. Columns found: {list(t.columns)[:8]} ...")
    new = {}
    for c in t.columns:
        label = str(c).split("\n")[0]
        label = re.sub(r"\s*\([^)]*\)\s*$", "", label).strip().lower() if label.lower() not in HEADER_TO_KEY else label.strip().lower()
        if label in HEADER_TO_KEY and HEADER_TO_KEY[label] not in new.values():
            new[c] = HEADER_TO_KEY[label]
    t = t.rename(columns=new)
    missing = [k for k in ("year", "event", "first_day", "last_day", "status", "ratio_3g", "ratio_2g") if k not in t.columns]
    if missing:
        raise SystemExit(f"{path}: these columns could not be found: {missing}. Columns found: {list(t.columns)[:10]} ...")
    return t.dropna(subset=["event"])


# ------------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Storm-by-storm figure for every storm in storm_table.csv (see the top of this file).")
    ap.add_argument("--table", default="storm_table.csv", help="storm_table.csv, a CSV made from the spreadsheet, or Sycamore_Storm_Table.xlsx")
    ap.add_argument("--root", default=str(HERE), help="folder that holds results/")
    ap.add_argument("--hourly-wy2017", default="hourly_basin_wy2017.csv")
    ap.add_argument("--hourly-wy2020", default="hourly_basin_wy2020.csv")
    ap.add_argument("--results", default=None, help="folder with the four *_Outlet.qout files (default <root>/results)")
    ap.add_argument("--obs", default=None, help="USGS gauge file, WY2017 (default <root>/../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv)")
    ap.add_argument("--obs-wy2020", default=None, help="USGS gauge file, WY2020 (default <root>/../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv)")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--no-volumes", action="store_true")
    a = ap.parse_args(argv)
    root = Path(a.root)
    res = Path(a.results) if a.results else root / "results"
    obs_by_year = {"WY2017": Path(a.obs) if a.obs else root.parent / "init_data" / "obs" / "sycamore_09510200_obs_flow_WY2017.csv",
                   "WY2020": Path(a.obs_wy2020) if a.obs_wy2020 else root.parent / "init_data" / "obs" / "sycamore_09510200_obs_flow_WY2020.csv"}
    out_dir = Path(a.out_dir) if a.out_dir else root / P.FIG_DIR

    tab = read_table(a.table)
    out_rows = []
    for year, hp in (("WY2017", a.hourly_wy2017), ("WY2020", a.hourly_wy2020)):
        log(f"\n{year}: rain")
        d = load_daily(Path(hp))
        t = tab[tab["year"] == year].reset_index(drop=True)
        rows = []
        for _, r in t.iterrows():
            first, last = pd.Timestamp(r["first_day"]), pd.Timestamp(r["last_day"])
            m = d.reindex(pd.date_range(first, last, freq="D"))
            st = {"event": r["event"], "first": first, "last": last,
                  "grid_mm": float(m["grid_mean_radar_mm"].sum()),
                  "ratio_3g": float(r["ratio_3g"]) if pd.notna(r["ratio_3g"]) else float("nan"),
                  "ratio_2g": float(r["ratio_2g"]) if pd.notna(r["ratio_2g"]) else float("nan")}
            s_ = str(r["status"])
            st["status_short"] = "20 mm list" if s_.startswith("On 20") else ("spin-up" if s_.startswith("Spin") else "below 20 mm")
            if m["grid_mean_radar_mm"].isna().any():
                log(f"  WARNING: {r['event']}: some storm days have no rain data")
            if abs(st["grid_mm"] - float(r["grid_mm"])) > 0.05:
                log(f"  WARNING: {r['event']}: grid rain {st['grid_mm']:.2f} mm here vs {float(r['grid_mm']):.2f} mm in the table")
            rows.append(st)
        vols = None
        runs = [r_ for r_ in RUNS_BY_YEAR[year] if year != "WY2020" or "kchan" not in r_[2] or (res / f"{r_[2]}_Outlet.qout").is_file()]
        if len(runs) < len(RUNS_BY_YEAR[year]):
            log(f"  {year}: Kchan outlet file(s) not found in {res}: those bars are left out")
        need = [obs_by_year[year]] + [res / f"{stem}_Outlet.qout" for _, _, stem in runs]
        absent = [str(x) for x in need if not Path(x).is_file()]
        if absent and not a.no_volumes:
            log(f"  {year}: no volume row (file not found): " + "; ".join(absent))
        if not a.no_volumes and not absent:
            gq = read_obs(obs_by_year[year])
            qs = [(lab, col, read_qout(res / f"{stem}_Outlet.qout", RUN_START[year])) for lab, col, stem in runs]
            vols = []
            for k_, st in enumerate(rows):
                lo = st["first"]
                hi = st["last"] + pd.Timedelta(days=1 + AFTER_DAYS)
                if k_ + 1 < len(rows):                                    # stop at the next storm's first day, so no flow is counted twice
                    hi = min(hi, rows[k_ + 1]["first"])
                g, gh = volume(gq, lo, hi)
                no_record = gh == 0                                       # the USGS file has no hours in this window (WY2020 September)
                g = 0.0 if not np.isfinite(g) else g
                v = {"gauge": g, "gauge_hours": gh, "hours_expected": int((hi - lo) / pd.Timedelta(hours=1)),
                     "no_record": no_record, "no_flow": (g < NO_FLOW_MM3) and not no_record, "runs": [volume(q, lo, hi)[0] for _, _, q in qs]}
                vols.append(v)
            log(f"\n{year} volumes (M m3), window = first rain day 00:00 to 3 days after the last rain day (or the next storm's first day)")
            for st, v in zip(rows, vols):
                log(f"  {st['event']:14s} gauge {v['gauge']:.4f} ({v['gauge_hours']}/{v['hours_expected']} h){'  *' if v['no_flow'] else ''}{'  no gauge record' if v['no_record'] else ''}  "
                    + "  ".join(f"{x:.4f}" for x in v["runs"]))
        for i, st in enumerate(rows):
            rec = {"year": year, "event": st["event"], "first_day": st["first"].date(), "last_day": st["last"].date(), "grid_mm": st["grid_mm"],
                   "ratio_3g": st["ratio_3g"], "ratio_2g": st["ratio_2g"]}
            if vols is not None:
                rec.update({"gauge_Mm3": float("nan") if vols[i]["no_record"] else vols[i]["gauge"], "gauge_no_flow": vols[i]["no_flow"],
                            "gauge_no_record": vols[i]["no_record"]})
                rec.update({lab: x for (lab, _, _), x in zip(runs, vols[i]["runs"])})
            out_rows.append(rec)
        log("  wrote " + make_year_figure(out_dir / f"fig_storm_by_storm_all_{year.lower()}.png", year, rows, d, vols, runs))
    pd.DataFrame(out_rows).to_csv(out_dir / "storm_by_storm_all.csv", index=False)
    log("  wrote " + str(out_dir / "storm_by_storm_all.csv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
