#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_candidate_anchors_sycamore_v2.py   (v2.1, 2026-11-01)
=======================================================

v2.1: HYETOGRAPHS. A hyetograph is a bar plot of rainfall against time, drawn upside down (bars hang from the top) so that
  the rain "falls" onto the hydrograph below it. Every hydrograph panel in candidate_panels.png and every storm panel in
  candidate_storms.png now has one on top. The rain is the basin-mean rain the MODEL is given (column MAP_mm_hr of the
  run's .mrf file, read the same way score_run_sycamore.py reads it; each hour's value is placed at the START of that hour).
  It is the same for every run, because only Ks, f and PsiB change between these runs. The script checks this when it can.
  The rain totals (mm) are printed on the screen and shown in the numbers box (panels) or the panel titles (storms).
  If no .mrf file is found, the pictures are drawn as before with a note where the rain would be.

v2.0 (the earlier version is plot_candidate_anchors_sycamore.py, v1.1, and is unchanged):
  * the default candidates are draws 12, 123 and 25;
  * in candidate_panels.png every draw has its own ROW with two panels of the same Jan 18 - Mar 12 picture: LOG flow
    on the left (small flows and recession tails) and LINEAR flow on the right (the real size of the big peaks);
  * draw 99 is no longer drawn anywhere (not in the panels, the storms picture or the cumulative picture); the two
    current anchors, 85 and 5, are still the gray references in the storms and cumulative pictures;
  * extra empty space is left above the data in every panel, and the legend and the numbers box sit in that space,
    so they never cover the flood-week peak;
  * the pictures go to plot_out_candidates_v2/ so v1's pictures are not overwritten.
  To get the old eight-draw set, use  --candidates 89,82,25,107,52,123,12,94

v1.1: the choice rule changed. Feb 18-19 is left out of the choice (it is only about 5% of the Oct-Mar water and
SRP cares about the big flows), the flood week may be up to 1.12 times the gauge volume, and the table now says,
for each run, whether it meets the NEW RULE (see below). Default candidates are now the eight draws listed below.

WHAT THIS SCRIPT DOES
---------------------
It lays the existing W2 model runs of several CANDIDATE soil anchors over the USGS 09510200 gauge for Jan 18 to
Mar 12 2017, so you can look at their hydrographs and choose the replacement anchor. It runs no model, changes
none of your files, and only READS the gauge file and the kept *_Outlet.qout files in results/.

It is built from plot_w4_gauge_slowlimb_sycamore.py (same file readers, same colors for the gauge). The slow-flow
cut is NOT used here: we decided to score on the whole flow inside each storm window and to state the missing slow
tail plainly. That script stays as a stand-alone diagnostic picture.

VOCABULARY
----------
  anchor       one (Ks, f, PsiB) soil setting that the routing sweep (cv, flowexp, n) will be run on.
  candidate    a draw we are considering as the replacement anchor.
  reference    a draw drawn in gray only for comparison: the two current anchors (85 and 5).
  draw         one row of the 150-draw Latin hypercube; "d089" in a run name means draw 89.
  hydrograph   a plot of flow (m3/s) against time.
  volume ratio model volume / gauge volume over a window. 1.00 = same water, below 1 = model too dry.
  peak error   (model peak / gauge peak - 1) as a percent. Negative = model peak too low.

THE NEW RULE (decided Oct 31 2026; the numbers are proposals until written into the design doc)
-----------------------------------------------------------------------------------------------
  Choose the anchor on the two BIG storms: the flood week (Jan 18-25) and Feb 27-28. Feb 18-19 is shown but not scored.
    flood week   volume ratio 0.90 to 1.12
    Feb 27-28    volume ratio 0.90 to 1.10
  The table columns flood_ok, feb27_ok and both_ok say whether a run meets each part. Peaks and arrival are SHOWN
  but not gated here: these runs all use the baseline routing (cv 6), and routing is what moves the peak (in Stage 2,
  cv 3 to 15 moved the flood-week peak from -50% to +50% without changing volume) and arrival (n). So the soil anchor is
  chosen on volume, and the peaks are fixed afterward by the routing sweep.

WHICH DRAWS (set in the lists at the top; change them with the options below)
-----------------------------------------------------------------------------
  Existing W2 runs (already in results/):
    89    best flood week and best Feb 27-28 of the 13 W2 runs, Feb 18-19 far too big
    82    the middle compromise
  New W2 runs, picked from the 150 W1 draws because their flood week is 1.06 to 1.12 and PsiB is low, which the
  13 W2 runs suggest brings Feb 27-28 up:
    25, 107, 52   flood week inside 0.90 to 1.10
    123, 12, 94   flood week 1.11 to 1.12 (inside the new 1.12 limit)
  The new draws are drawn only after their W2 runs exist; until then the script skips them with a message.
  References (gray): 85 and 5.
  No setting fixes the storm balance (Stage 2 result), so every candidate is a compromise. The pictures are there to
  show which compromise you prefer to look at.

PICTURES AND FILES IT WRITES (to --out-dir, default plot_out_candidates_v2/)
-------------------------------------------------------------------------
  candidate_panels.png   one ROW per candidate, Jan 18 - Mar 12: gauge (dark), the candidate (color).
                         LEFT panel = LOG flow scale, RIGHT panel = LINEAR flow scale. The numbers box (same
                         numbers as the CSV) is on the left panel.
  candidate_storms.png   three panels (flood week Jan 18-25, Feb 18-21, Feb 27 - Mar 3), LINEAR scale, all candidates
                         and the gray references together, with the gauge in black.
  candidate_cumulative.png  running total of water (million m3) from Jan 18: how fast each run "fills the bucket"
                         against the gauge. The gap at the end is the missing water; the steps are the storms.
  candidate_table.csv    the numbers behind the pictures (volume ratio and peak error per storm, whole window ratio).
  A table is also printed on the screen.

  Volumes here use a simple rule (model flow interpolated onto the gauge's 15-minute times, counted only where the
  gauge has data, rectangle sum), so they can differ from score_run_sycamore.py by about 1%. To see the difference
  add  --w2-summary w2_summary.csv  and the script prints the scorer's flood-week ratio next to its own.

HOW TO RUN IT (from workspaces/Sycamore/lab/)
---------------------------------------------
  python -u plot_candidate_anchors_sycamore.py
  python -u plot_candidate_anchors_sycamore.py --w2-summary w2_summary.csv
  python -u plot_candidate_anchors_sycamore.py --candidates 89,82,35 --refs 85,5
  python -u plot_candidate_anchors_sycamore.py --candidates 89,W2_ks0p2075_cv6_psib0p4317_f0p8215_d089_nocl
  python -u plot_candidate_anchors_sycamore.py --results-dir results --obs ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv

  A draw can be given as its number (89) or as a full run id. For a number the script looks for
  results/W2_*_dNNN_nocl_Outlet.qout (and uses a W1_ file if there is no W2 one: it then only covers Jan 18-25).
  A draw that cannot be found is skipped with a message; the others are still drawn.

It takes a few seconds, runs in the foreground, and needs pandas, numpy and matplotlib (already in your environment).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")                      # draw to files only: no window is needed in a Codespace
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import pandas as pd

VERSION = "2.1"
DEFAULT_OBS = "../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv"

CANDIDATES = [12, 123, 25]

# The new choice rule (Oct 31 2026): the two big storms only. Feb 18-19 is not scored.
FLOOD_VOL_BAND = (0.90, 1.12)
FEB27_VOL_BAND = (0.90, 1.10)
REFERENCES = {85: "current anchor", 5: "current anchor"}
REF_STYLES = {85: (0, (6, 2)), 5: (0, (4, 1, 1, 1)), 99: (0, (1, 1.5))}     # dashed, dash-dot, dotted
DEFAULT_REF_STYLE = (0, (3, 3))

RUN_STARTS = {"W1": "2017-01-18 00:00", "W2": "2017-01-18 00:00",
              "W3": "2016-10-01 00:00", "W4": "2016-10-01 00:00", "W5": "2016-10-01 00:00"}

SPAN = ("2017-01-18 00:00", "2017-03-12 00:00")                 # the W2 window
# Windows (local Arizona time, [start, end)). Same as the scorer's flood week and the earlier plot script's storm windows.
EVENTS = {
    "flood": dict(label="Flood week (Jan 18-25)", win=("2017-01-18 00:00", "2017-01-25 00:00")),
    "feb18": dict(label="Feb 18-19 storm (not scored under the new rule)", win=("2017-02-18 00:00", "2017-02-21 00:00")),
    "feb27": dict(label="Feb 27-28 storm", win=("2017-02-27 00:00", "2017-03-03 00:00")),
}

C_GAUGE = "#1b3a5c"
C_REF = "#7f7f7f"
CAND_COLORS = ["#d55e00", "#009e73", "#cc79a7", "#0072b2", "#e69f00", "#56b4e9", "#8c564b", "#000000"]
LOG_FLOOR = 0.05                    # m3/s: the bottom of the log-scale panels
LOG_TOP = 3000.0                    # m3/s: top of the log panels. The gauge peak is about 220, so the top part is empty space for the legend and numbers box
C_RAIN = "#5b8fc7"                  # hyetograph bars
RAIN_LIMIT_FACTOR = 1.6             # the rain axis reaches this many times the highest hourly rain (the bars hang from the top)
LIN_HEADROOM = 1.55                 # the linear panels reach this many times the highest flow drawn

_trap = getattr(np, "trapezoid", None) or np.trapz     # numpy 2 renamed trapz


# ----------------------------------------------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------------------------------------------
def die(msg: str) -> None:
    print(f"\nERROR: {msg}\n", file=sys.stderr)
    sys.exit(2)


def ts(x) -> pd.Timestamp:
    return pd.Timestamp(x)


def between(s: pd.Series, a, b) -> pd.Series:
    """The part of a series with a <= time < b (explicit comparisons: text slicing would include the whole last day)."""
    a, b = ts(a), ts(b)
    return s[(s.index >= a) & (s.index < b)]


# ----------------------------------------------------------------------------------------------------------------
# reading files
# ----------------------------------------------------------------------------------------------------------------
def read_gauge(path: str) -> pd.Series:
    p = Path(path)
    if not p.exists():
        die(f"gauge file not found: {p}   (use --obs to point to sycamore_09510200_obs_flow_WY2017.csv)")
    df = pd.read_csv(p, parse_dates=["datetime"]).set_index("datetime")
    if "q_cms" not in df.columns:
        die(f"{p} has no q_cms column; columns found: {list(df.columns)}")
    s = pd.to_numeric(df["q_cms"], errors="coerce")
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_qout(path: Path, start) -> pd.Series:
    """Routed outlet flow (m3/s) from a tRIBS *_Outlet.qout file: a header line, then Time_hr and Q."""
    with open(path) as fh:
        fh.readline()
        second = fh.readline()
    if "," in second:
        raw = pd.read_csv(path, sep=",", header=0)
    else:
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None)
    raw = raw.iloc[:, :2].copy()
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    s = pd.Series(raw["Q"].to_numpy(float), index=ts(start) + pd.to_timedelta(raw["Time_hr"].to_numpy(), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def find_run(item: str, results_dir: Path):
    """item = a draw number ('89') or a full run id. Returns (run_id, qout path, window name) or None."""
    item = str(item).strip()
    if not results_dir.is_dir():
        print(f"  skipped {item}: results folder not found: {results_dir}")
        return None
    if re.fullmatch(r"\d+", item):
        n = int(item)
        for pre in ("W2", "W1"):
            hits = sorted(results_dir.glob(f"**/{pre}_*_d{n:03d}_nocl_Outlet.qout"))
            if len(hits) == 1:
                rid = hits[0].name[: -len("_Outlet.qout")]
                return rid, hits[0], pre
            if len(hits) > 1:
                print(f"  skipped {item}: {len(hits)} files match {pre}_*_d{n:03d}_nocl_Outlet.qout; give the full run id")
                return None
        print(f"  skipped {item}: no W2_ or W1_ outlet file for draw {n} under {results_dir}")
        return None
    hits = sorted(results_dir.glob(f"**/{item}_Outlet.qout"))
    if len(hits) != 1:
        print(f"  skipped {item}: expected exactly one {item}_Outlet.qout under {results_dir}, found {len(hits)}")
        return None
    return item, hits[0], item[:2].upper()


def parse_params(run_id: str) -> dict:
    """Ks, f and PsiB multipliers from a run id like W2_ks0p2524_cv6_psib0p5743_f1p165_d099_nocl."""
    def grab(key):
        m = re.search(rf"_{key}(\d+p?\d*)", run_id)
        return float(m.group(1).replace("p", ".")) if m else np.nan
    m = re.search(r"_d(\d+)_", run_id)
    return dict(draw=int(m.group(1)) if m else -1, ks=grab("ks"), f=grab("f"), psib=grab("psib"))


def read_mrf_rain(path: Path, start) -> pd.Series:
    """Basin-mean rain (mm per hour) from a tRIBS .mrf file. The file labels each hour by its END, so the time is moved back one
    hour: the value then sits at the START of the hour it describes (as in score_run_sycamore.py). The unfilled last row is dropped."""
    df = pd.read_csv(path)
    if "Time_hr" not in df.columns or "MAP_mm_hr" not in df.columns:
        raise ValueError(f"no Time_hr / MAP_mm_hr columns; columns found: {list(df.columns)}")
    if "RainMin_mm_hr" in df.columns:
        df = df[pd.to_numeric(df["RainMin_mm_hr"], errors="coerce") < 9000]
    df = df.apply(pd.to_numeric, errors="coerce").dropna(subset=["Time_hr", "MAP_mm_hr"])
    idx = ts(start) + pd.to_timedelta(df["Time_hr"].to_numpy(float), unit="h") - pd.Timedelta(hours=1)
    s = pd.Series(df["MAP_mm_hr"].to_numpy(float), index=idx)
    return s[~s.index.duplicated(keep="last")].sort_index()


def load_rain(runs, results_dir: Path):
    """Rain for the pictures: from the run whose .mrf file covers the longest time. Every other run's .mrf is compared with it
    where they overlap (the rain should be the same for all of them). Returns (series or None, text saying where it came from)."""
    found = []
    for r in runs:
        hits = sorted(results_dir.glob(f"**/{r['run_id']}[0-9]*.mrf"))
        if not hits:
            continue
        start = RUN_STARTS.get(r["window"])
        try:
            found.append((r, hits[-1], read_mrf_rain(hits[-1], start)))
        except Exception as exc:                      # a damaged or unexpected file must not stop the pictures
            print(f"  rain: could not read {hits[-1].name}: {exc}")
    if not found:
        return None, "no .mrf file found next to the outlet files"
    found.sort(key=lambda t: len(t[2]), reverse=True)
    main_run, main_path, main = found[0]
    for r, path, other in found[1:]:
        common = other.index.intersection(main.index)
        if len(common) == 0:
            continue
        tm, to = float(main.loc[common].sum()), float(other.loc[common].sum())
        if tm > 0 and abs(to / tm - 1) > 0.005:
            print(f"  WARNING: the rain in {path.name} differs from the rain in {main_path.name} by {100 * (to / tm - 1):+.1f}% "
                  f"over their common hours. The pictures use {main_path.name}.")
    return main, f"{main_path.name} (basin-mean rain, {main.index.min():%b %d} to {main.index.max():%b %d})"


def rain_total_mm(rain, a, b) -> float:
    """Basin-mean rain depth (mm) in [a, b): the sum of the hourly values (each is mm per hour for one hour)."""
    if rain is None:
        return float("nan")
    r = between(rain, a, b)
    return float(r.sum()) if len(r) else float("nan")


def draw_hyeto(ax, rain, a, b, ymax):
    """An upside-down bar plot of hourly rain on ax between a and b."""
    ax.grid(alpha=0.25)
    ax.tick_params(labelbottom=False, labelsize=8)
    if rain is None or len(between(rain, a, b)) == 0:
        ax.text(0.5, 0.5, "no rain data (.mrf file not found)", transform=ax.transAxes, ha="center", va="center", fontsize=8, color="0.4")
        ax.set_yticks([])
        return
    r = between(rain, a, b)
    ax.fill_between(r.index, 0, r.to_numpy(float), step="post", color=C_RAIN, linewidth=0)
    ax.set_ylim(ymax * RAIN_LIMIT_FACTOR, 0)              # upside down: 0 at the top
    ax.set_ylabel("rain\n(mm/hr)", fontsize=8)


# ----------------------------------------------------------------------------------------------------------------
# numbers
# ----------------------------------------------------------------------------------------------------------------
def on_gauge_times(model: pd.Series, gauge: pd.Series, a, b):
    """Model and gauge flow on the gauge's own 15-minute times in [a, b), counted only where the gauge has data.
    Returns (model array, gauge array, index), both with NaN where the gauge is blank or the model has no coverage."""
    g = between(gauge, a, b)
    m = model.reindex(model.index.union(g.index)).interpolate(method="time", limit_area="inside").reindex(g.index)
    m = m.where(g.notna())
    return m, g, g.index


def rect_volume_m3(q: pd.Series, step_s: float = 900.0) -> float:
    return float(np.nansum(q.to_numpy(float)) * step_s)


def event_numbers(model: pd.Series, gauge: pd.Series, win) -> dict:
    m, g, idx = on_gauge_times(model, gauge, *win)
    covered = m.notna().sum()
    if covered < 0.9 * g.notna().sum() or g.notna().sum() == 0:
        return dict(vol_ratio=np.nan, peak_err_pct=np.nan, vol_model_Mm3=np.nan, vol_gauge_Mm3=np.nan)
    vm = rect_volume_m3(m) / 1e6
    vg = rect_volume_m3(g.where(m.notna())) / 1e6
    pm = float(between(model, *win).max())
    pg = float(g.max())
    return dict(vol_ratio=vm / vg if vg > 0 else np.nan,
                peak_err_pct=(pm / pg - 1.0) * 100.0 if pg > 0 else np.nan,
                vol_model_Mm3=vm, vol_gauge_Mm3=vg)


def cumulative(model: pd.Series, gauge: pd.Series, a, b):
    """Running totals of water (million m3) from time a, on the gauge's 15-minute times.
    Returns (model total, gauge total). The gauge total depends on the gauge alone (blanks count as zero water), so it
    never depends on which run is drawn. The model total is counted only where the gauge has data, and it is left blank
    (NaN) after the end of a run that stops early (a one-week W1 run), so it is not drawn as a flat line."""
    g = between(gauge, a, b)
    m, _, _ = on_gauge_times(model, gauge, a, b)
    cg = (g.fillna(0.0) * 900.0).cumsum() / 1e6
    cm = (m.fillna(0.0) * 900.0).cumsum() / 1e6
    last = model.index.max()
    cm = cm.where(cm.index <= last)
    return cm, cg


# ----------------------------------------------------------------------------------------------------------------
# pictures
# ----------------------------------------------------------------------------------------------------------------
def fmt_date_axis(ax):
    ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MO, interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)


def pct(x):
    return "n/a" if not np.isfinite(x) else f"{x:+.0f}%"


def rat(x):
    return "n/a" if not np.isfinite(x) else f"{x:.2f}"


def draw_panels(gauge, cands, refs, rain, out_png, dpi):
    """One row per candidate: LOG flow on the left, LINEAR flow on the right. Above every hydrograph: its hyetograph."""
    n = len(cands)
    fig_h = 4.8 * n
    fig = plt.figure(figsize=(15, fig_h))
    outer = GridSpec(n, 2, figure=fig, hspace=0.22, wspace=0.14, left=0.065, right=0.99,
                     top=1 - 1.25 / fig_h, bottom=0.95 / fig_h)
    a, b = ts(SPAN[0]), ts(SPAN[1])
    g = between(gauge, a, b)
    rain_max = float(between(rain, a, b).max()) if (rain is not None and len(between(rain, a, b))) else 1.0
    first = None
    for i, c in enumerate(cands):
        tot = {k: rain_total_mm(rain, *ev["win"]) for k, ev in EVENTS.items()}
        for j, scale in enumerate(("log", "linear")):
            inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i, j], height_ratios=[1, 3.6], hspace=0.07)
            axh = fig.add_subplot(inner[0], sharex=first)
            ax = fig.add_subplot(inner[1], sharex=axh)
            first = first or axh
            draw_hyeto(axh, rain, a, b, rain_max)
            floor = LOG_FLOOR if scale == "log" else 0.0
            ax.plot(g.index, np.maximum(g.fillna(0.0).to_numpy(float), floor), color=C_GAUGE, lw=1.4, label="gauge")
            s = between(c["series"], a, b)
            ax.plot(s.index, np.maximum(s.to_numpy(float), floor), color=c["color"], lw=1.6, label=f"draw {c['draw']}")
            if scale == "log":
                ax.set_yscale("log")
                ax.set_ylim(LOG_FLOOR, LOG_TOP)
                ax.set_ylabel("flow (m3/s), LOG scale", fontsize=9)
            else:
                top = max(float(np.nanmax(g.to_numpy(float))), float(np.nanmax(s.to_numpy(float))))
                ax.set_ylim(0, top * LIN_HEADROOM)
                ax.set_ylabel("flow (m3/s), LINEAR scale", fontsize=9)
            ax.set_xlim(a, b)
            ax.grid(alpha=0.25)
            ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
            if i == n - 1:
                fmt_date_axis(ax)
            else:
                ax.tick_params(labelbottom=False)
            if j == 0:
                axh.set_title(f"draw {c['draw']}   Ks {c['ks']:.3g}   f {c['f']:.3g}   PsiB {c['psib']:.3g}"
                              + ("" if c["window"] == "W2" else f"   [{c['window']} run: covers a week only]"),
                              fontsize=9.5, loc="left")
                t = c["nums"]
                rain_line = ("" if rain is None else
                             f"\nrain (mm): flood {tot['flood']:.0f}, Feb 18-19 {tot['feb18']:.0f}, Feb 27-28 {tot['feb27']:.0f}")
                ax.text(0.17, 0.97,
                        f"flood week: vol {rat(t['flood']['vol_ratio'])}, peak {pct(t['flood']['peak_err_pct'])}\n"
                        f"Feb 18-19*: vol {rat(t['feb18']['vol_ratio'])}, peak {pct(t['feb18']['peak_err_pct'])}   (*not scored)\n"
                        f"Feb 27-28:  vol {rat(t['feb27']['vol_ratio'])}, peak {pct(t['feb27']['peak_err_pct'])}" + rain_line,
                        transform=ax.transAxes, va="top", ha="left", fontsize=8, family="monospace",
                        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.9))
            else:
                axh.set_title("same run, linear flow scale", fontsize=9.5, loc="left")
    fig.suptitle("Candidate anchors against the gauge, Jan 18 - Mar 12 2017 (W2 runs, routing at baseline cv 6)\n"
                 "top of each panel: rain the model is given (basin mean, mm/hr)   left: log flow   right: linear flow   "
                 "vol = model / gauge volume (1.00 = same water); peak = model peak vs gauge peak", fontsize=11, y=0.995)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)


def draw_storms(gauge, cands, refs, rain, out_png, dpi):
    fig = plt.figure(figsize=(17, 6.8))
    outer = GridSpec(1, 3, figure=fig, wspace=0.22, left=0.05, right=0.99, top=0.84, bottom=0.13)
    for k, (key, ev) in enumerate(EVENTS.items()):
        a, b = ts(ev["win"][0]), ts(ev["win"][1])
        inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[k], height_ratios=[1, 3.8], hspace=0.07)
        axh = fig.add_subplot(inner[0])
        ax = fig.add_subplot(inner[1], sharex=axh)
        rw = between(rain, a, b) if rain is not None else None
        draw_hyeto(axh, rain, a, b, float(rw.max()) if (rw is not None and len(rw)) else 1.0)
        g = between(gauge, a, b)
        for r in refs:
            s = between(r["series"], a, b)
            if len(s):
                ax.plot(s.index, s.to_numpy(float), color=C_REF, lw=1.0, ls=REF_STYLES.get(r["draw"], DEFAULT_REF_STYLE),
                        label=f"ref {r['draw']} ({REFERENCES.get(r['draw'], 'reference')})")
        for c in cands:
            s = between(c["series"], a, b)
            if len(s):
                ax.plot(s.index, s.to_numpy(float), color=c["color"], lw=1.5, label=f"draw {c['draw']}")
        ax.plot(g.index, g.to_numpy(float), color="black", lw=2.2, label="gauge", zorder=5)
        ax.set_xlim(a, b)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.4)                 # empty space above the highest line, for the legend
        tot = rain_total_mm(rain, a, b)
        axh.set_title(ev["label"] + ("" if not np.isfinite(tot) else f"     rain {tot:.0f} mm"), fontsize=10.5, loc="left")
        ax.set_ylabel("flow (m3/s), linear scale", fontsize=9)
        ax.grid(alpha=0.25)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
        if k == 0:
            ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
    fig.suptitle("The three storms: rain the model is given (top), candidates (colors), references (gray) and gauge (black)", fontsize=11.5, y=0.97)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)


def draw_cumulative(gauge, cands, refs, out_png, dpi):
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    a, b = ts(SPAN[0]), ts(SPAN[1])
    any_series = cumulative(cands[0]["series"], gauge, a, b)[1]
    for r in refs:
        cm, cg = cumulative(r["series"], gauge, a, b)
        ax.plot(cm.index, cm.to_numpy(float), color=C_REF, lw=1.0, ls=REF_STYLES.get(r["draw"], DEFAULT_REF_STYLE),
                label=f"ref {r['draw']} ({REFERENCES.get(r['draw'], 'reference')})")
    for c in cands:
        cm, cg = cumulative(c["series"], gauge, a, b)
        ax.plot(cm.index, cm.to_numpy(float), color=c["color"], lw=1.5, label=f"draw {c['draw']}")
    ax.plot(any_series.index, any_series.to_numpy(float), color="black", lw=2.4, label="gauge", zorder=5)
    ax.set_xlim(a, b)
    ax.set_ylabel("running total of water since Jan 18 (million m3)")
    ax.set_title("Running total of water: each step up is a storm. The gap at the end is the water the model misses.",
                 fontsize=10.5, loc="left")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.9)
    fmt_date_axis(ax)
    fig.tight_layout()
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="Lay candidate anchor runs over the gauge, Jan-Mar 2017.")
    ap.add_argument("--obs", default=DEFAULT_OBS, help=f"gauge CSV (default {DEFAULT_OBS})")
    ap.add_argument("--results-dir", default="results", help="folder holding the *_Outlet.qout files (default results)")
    ap.add_argument("--candidates", default=",".join(str(c) for c in CANDIDATES),
                    help="comma list of draw numbers or full run ids to treat as candidates")
    ap.add_argument("--refs", default=",".join(str(r) for r in REFERENCES),
                    help="comma list of draw numbers drawn in gray for comparison (default 85,5)")
    ap.add_argument("--out-dir", default="plot_out_candidates_v2", help="folder for the pictures and CSV")
    ap.add_argument("--w2-summary", default=None, help="optional w2_summary.csv: prints the scorer's flood-week ratio next to ours")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    print(f"plot_candidate_anchors_sycamore_v2.py v{VERSION}")
    gauge = read_gauge(args.obs)
    rd = Path(args.results_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    def load(items):
        loaded = []
        for it in [x.strip() for x in items.split(",") if x.strip()]:
            hit = find_run(it, rd)
            if hit is None:
                continue
            rid, path, pre = hit
            start = RUN_STARTS.get(pre)
            if start is None:
                print(f"  skipped {it}: cannot tell when {rid} starts (run ids must begin W1_ ... W5_)")
                continue
            series = read_qout(path, start)
            p = parse_params(rid)
            nums = {k: event_numbers(series, gauge, ev["win"]) for k, ev in EVENTS.items()}
            span = event_numbers(series, gauge, SPAN)
            loaded.append(dict(run_id=rid, path=str(path), window=pre, series=series, nums=nums, span=span, **p))
            print(f"  read: draw {p['draw']:>3}  {rid}  ({series.index.min():%b %d} to {series.index.max():%b %d})")
        return loaded

    print("Candidates:")
    cands = load(args.candidates)
    print("References:")
    refs = load(args.refs)
    if not cands:
        die("no candidate runs could be read. Check --results-dir (run this from workspaces/Sycamore/lab/).")
    for i, c in enumerate(cands):
        c["color"] = CAND_COLORS[i % len(CAND_COLORS)]

    print("Rain:")
    rain, rain_note = load_rain(cands + refs, rd)
    print(f"  {rain_note}")
    if rain is not None:
        print("  basin-mean rain the model is given (mm): " + ", ".join(
            f"{ev['label'].split(' (')[0]} {rain_total_mm(rain, *ev['win']):.1f}" for ev in EVENTS.values())
            + f", whole window {rain_total_mm(rain, *SPAN):.1f}")

    # ---- the table ----
    rows = []
    for role, group in (("candidate", cands), ("reference", refs)):
        for r in group:
            row = dict(role=role, draw=r["draw"], run_id=r["run_id"], window=r["window"], ks=r["ks"], f=r["f"], psib=r["psib"])
            for k in EVENTS:
                row[f"{k}_vol_ratio"] = r["nums"][k]["vol_ratio"]
                row[f"{k}_peak_err_pct"] = r["nums"][k]["peak_err_pct"]
            row["whole_window_vol_ratio"] = r["span"]["vol_ratio"]
            rows.append(row)
    table = pd.DataFrame(rows)

    # ---- does each run meet the new rule? (blank = not enough of the run to tell) ----
    def within(x, lo, hi):
        return (x >= lo) & (x <= hi)

    flood_ok = within(table["flood_vol_ratio"], *FLOOD_VOL_BAND)
    feb27_ok = within(table["feb27_vol_ratio"], *FEB27_VOL_BAND)
    have = table["feb27_vol_ratio"].notna() & table["flood_vol_ratio"].notna()
    def yesno(ok, known):
        return np.where(known, np.where(ok, "yes", "no"), "")
    table["flood_ok"] = yesno(flood_ok, table["flood_vol_ratio"].notna())
    table["feb27_ok"] = yesno(feb27_ok, table["feb27_vol_ratio"].notna())
    table["both_ok"] = yesno(flood_ok & feb27_ok, have)

    if args.w2_summary:
        p = Path(args.w2_summary)
        if p.exists():
            sc = pd.read_csv(p, usecols=["run_id", "flood_week__vol_ratio", "flood_week__peak_err_pct", "span__vol_ratio"])
            table = table.merge(sc.rename(columns={"flood_week__vol_ratio": "scorer_flood_vol_ratio",
                                                   "flood_week__peak_err_pct": "scorer_flood_peak_err_pct",
                                                   "span__vol_ratio": "scorer_whole_vol_ratio"}), on="run_id", how="left")
            bad = table[(table["scorer_flood_vol_ratio"].notna())
                        & ((table["scorer_flood_vol_ratio"] - table["flood_vol_ratio"]).abs() > 0.03)]
            if len(bad):
                print("\nWARNING: flood-week volume ratio differs from the scorer's by more than 0.03 for draws "
                      f"{bad['draw'].tolist()}; check that the files and windows are what you expect.")
        else:
            print(f"  --w2-summary file not found: {p} (skipped)")

    table.to_csv(out / "candidate_table.csv", index=False, float_format="%.4f")

    show = ["role", "draw", "ks", "f", "psib", "flood_vol_ratio", "flood_peak_err_pct", "feb18_vol_ratio", "feb18_peak_err_pct",
            "feb27_vol_ratio", "feb27_peak_err_pct", "whole_window_vol_ratio", "flood_ok", "feb27_ok", "both_ok"]
    if "scorer_flood_vol_ratio" in table.columns:
        show += ["scorer_flood_vol_ratio"]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print("\nNumbers (vol = model / gauge; peak = percent vs gauge peak).")
        print(f"New rule (volume only): flood week vol {FLOOD_VOL_BAND[0]}-{FLOOD_VOL_BAND[1]}, Feb 27-28 vol "
              f"{FEB27_VOL_BAND[0]}-{FEB27_VOL_BAND[1]}. Peaks are shown but not gated (routing fixes them later). "
              f"Feb 18-19 is shown but not scored.")
        print(table[show].round(3).to_string(index=False))

    draw_panels(gauge, cands, refs, rain, out / "candidate_panels.png", args.dpi)
    draw_storms(gauge, cands, refs, rain, out / "candidate_storms.png", args.dpi)
    draw_cumulative(gauge, cands, refs, out / "candidate_cumulative.png", args.dpi)
    print(f"\nWROTE {out / 'candidate_panels.png'}\n      {out / 'candidate_storms.png'}\n      "
          f"{out / 'candidate_cumulative.png'}\n      {out / 'candidate_table.csv'}")


if __name__ == "__main__":
    main()
