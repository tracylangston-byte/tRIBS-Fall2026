#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rain_runoff_storms_sycamore.py   (v1.0)
=======================================

QUESTION IT ANSWERS
-------------------
The model makes 3 to 8 m3/s of flow in the fall (Oct 8, Nov 4, Dec 22) where the gauge sees almost nothing.
Is that because the RAIN in those storms was big (then a streambed that soaks up water would be the fix),
or because the model turns a modest amount of rain into runoff too easily (then the soil settings or the
rain input are the suspects)?

HOW IT ANSWERS
--------------
It finds every storm in the run window (October to March) from the basin-average rain the model was given, then,
for each storm, puts side by side:

  * the rain: total (mm), heaviest hour (mm/hr), and how wet the days before it were (7-day and 30-day rain)
  * the model's water table depth at the start of the storm (from the .mrf; a bigger number is deeper)
  * the gauge outflow and the model outflow over the storm plus two days of recession (M m3 and m3/s peak)
  * the RUNOFF RATIO = outflow volume / rain volume (rain volume = basin-average rain depth x basin area)

It does NOT run tRIBS and does not change any model file. It reads files that are already in results/ and the
gauge file. Output: one table on screen, one CSV, and one figure (rain per storm against outflow per storm).

WHAT IT NEEDS
-------------
  results/<run>_Outlet.qout   the routed flow of the run (the scorer's reader is used)
  results/<run>NNNN.mrf       the basin-average rain (MAP_mm_hr) and water table (MDGW_mm)
  ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv   the gauge file (the scorer's default)
  Nothing else: the file readers and the volume method below are copied from score_run_sycamore.py (v1.0), so the numbers
  agree with the scorer's, but this script does not import it and works in any folder.

USAGE (from workspaces/Sycamore/lab/)
-------------------------------------
  python rain_runoff_storms_sycamore.py --dry-run
  python rain_runoff_storms_sycamore.py
  python rain_runoff_storms_sycamore.py --also W4_ks0p2524_cv6_psib0p5743_f1p165_d099_cc2   # add a Kchan run's outflow
  python rain_runoff_storms_sycamore.py --run <another run id> --window W4

WORDS USED
----------
  basin-average rain   the average over the whole watershed of the radar rain the model was given (the .mrf column MAP_mm_hr)
  runoff ratio         outflow volume divided by rain volume; 0.1 means 10% of the rain left the basin as streamflow
  antecedent rain      rain that fell BEFORE the storm; it sets how wet the ground is when the storm arrives
  water table (MDGW)   mean depth to the water table in the model (mm); a BIGGER number means a DEEPER water table
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402

VERSION = "1.0"
HOUR = pd.Timedelta(hours=1)

DEFAULT_RUN = "W4_ks0p2524_cv6_psib0p5743_f1p165_d099_nocl"   # draw 99, W4, channel loss off
DEFAULT_WINDOW = "W4"

RAIN_DAY_MM = 2.0        # a day counts as a rainy day if the basin-average rain is at least this (mm)
MERGE_GAP_DAYS = 2       # rainy days separated by up to this many dry days belong to the same storm
TAIL_DAYS = 2            # the storm window runs this many days past the last rainy day (to catch the recession)
MIN_STORM_MM = 5.0       # list a storm if its rain is at least this (mm) ...
LIST_IF_FLOW_ABOVE = 0.13  # ... or if the gauge or any run goes above this (m3/s), the design's false-flow threshold
FALL_END = "2017-01-01"  # storms that start before this date are tagged "fall", the rest "winter"
EXCLUDE = [("2017-01-01 00:00", "2017-01-03 00:00")]   # the gauge outage, left out of both sides as the scorer does


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    print(f"\nERROR: {msg}\n", file=sys.stderr, flush=True)
    sys.exit(2)


# ------------------------------------------------------------------------------------------------------------------
# readers and the volume method: copied from score_run_sycamore.py v1.0 so the numbers agree with the scorer's
# ------------------------------------------------------------------------------------------------------------------
BASIN_KM2 = 424.6        # watershed area from the build (424.6151 km2), as in the scorer
MAX_STEP_H = 3.0         # a gap between two flow samples longer than this counts as a blank (no volume)
DEFAULT_OBS = "../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv"
RUN_STARTS = {"W1": "2017-01-18 00:00", "W2": "2017-01-18 00:00", "W3": "2016-10-01 00:00",
              "W4": "2016-10-01 00:00", "W5": "2016-10-01 00:00", "V": "2019-09-01 00:00"}


def ts(x) -> pd.Timestamp:
    return pd.Timestamp(x)


def run_window_from_id(run_id: str):
    import re
    m = re.match(r"^(W[1-5]|V)(?=[_\-\s.]|$)", run_id, flags=re.IGNORECASE)
    return m.group(1).upper() if m else None


def read_qout(path: Path, start) -> pd.Series:
    """Routed outlet flow (m3/s). Comma-delimited with a header line; an older whitespace layout is also accepted."""
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


def read_obs(path: str) -> pd.Series:
    if not Path(path).exists():
        die(f"gauge file not found: {path} (use --obs)")
    obs = pd.read_csv(path, parse_dates=["datetime"]).set_index("datetime")
    s = pd.to_numeric(obs["q_cms"], errors="coerce")
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_mrf(path: Path, start) -> pd.DataFrame:
    sim = pd.read_csv(path)
    sim = sim[sim["RainMin_mm_hr"] < 9000]                       # drop the unfilled last row
    sim.index = ts(start) + pd.to_timedelta(sim["Time_hr"], unit="h")
    return sim


def find_run_files(run_id: str, results_dir: Path):
    qout = sorted(results_dir.glob(f"**/{run_id}_Outlet.qout"))
    if len(qout) != 1:
        return None, None, f"expected exactly one {run_id}_Outlet.qout under {results_dir}, found {len(qout)}"
    mrf = sorted(results_dir.glob(f"**/{run_id}[0-9]*.mrf"))
    return qout[0], (mrf[-1] if mrf else None), ""


class Flow:
    """A flow series with trapezoid steps. A step is blank (no volume) if either end is NaN or the gap is over MAX_STEP_H."""

    def __init__(self, q: pd.Series):
        q = q[~q.index.duplicated(keep="last")].sort_index()
        self.q = q
        t = q.index.values.astype("datetime64[s]").astype("int64")
        v = q.to_numpy(float)
        self.t0, self.t1 = t[:-1], t[1:]
        dt = np.diff(t)
        self.valid = ~np.isnan(v[1:]) & ~np.isnan(v[:-1]) & (dt <= MAX_STEP_H * 3600)
        self.vol = np.where(self.valid, 0.5 * (v[1:] + v[:-1]) * dt, 0.0)

    def _keep(self, a, b, exclude):
        a_s, b_s = ts(a).value // 10**9, ts(b).value // 10**9
        keep = (self.t0 >= a_s) & (self.t0 < b_s)
        for ea, eb in exclude:
            keep &= ~((self.t0 >= ts(ea).value // 10**9) & (self.t0 < ts(eb).value // 10**9))
        return keep

    def volume(self, a, b, exclude=(), match_to: "Flow | None" = None) -> float:
        """Volume (m3) of the steps that START in [a, b), minus excluded times. With match_to, only steps whose midpoint
        falls in a valid (non-blank) step of that other series count."""
        keep = self._keep(a, b, exclude) & self.valid
        if match_to is not None and len(match_to.t0):
            mid = 0.5 * (self.t0 + self.t1)
            i = np.searchsorted(match_to.t0, mid, side="right") - 1
            ok = (i >= 0)
            ii = np.clip(i, 0, len(match_to.t0) - 1)
            ok &= (mid < match_to.t1[ii]) & match_to.valid[ii]
            keep &= ok
        return float(self.vol[keep].sum())

    def has_data(self, a, b, exclude=()) -> bool:
        return bool((self._keep(a, b, exclude) & self.valid).any())

    def window(self, a, b, exclude=()) -> pd.Series:
        s = self.q[(self.q.index >= ts(a)) & (self.q.index < ts(b))]
        for ea, eb in exclude:
            s = s[~((s.index >= ts(ea)) & (s.index < ts(eb)))]
        return s


def fmt(x, nd=2) -> str:
    try:
        return "-" if (x is None or not math.isfinite(float(x))) else f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return "-"


# ------------------------------------------------------------------------------------------------------------------
# finding storms from the rain
# ------------------------------------------------------------------------------------------------------------------
def find_storms(rain_hr: pd.Series):
    """Group rainy days into storms. Returns a list of (first_day, last_day) Timestamps (midnight)."""
    daily = rain_hr.resample("D").sum()
    rainy = daily[daily >= RAIN_DAY_MM].index
    storms = []
    for d in rainy:
        if storms and (d - storms[-1][1]).days - 1 <= MERGE_GAP_DAYS:
            storms[-1][1] = d
        else:
            storms.append([d, d])
    return [(a, b) for a, b in storms]


def storm_label(a: pd.Timestamp, b: pd.Timestamp) -> str:
    if a == b:
        return f"{a:%b %d}"
    if a.month == b.month:
        return f"{a:%b %d}-{b:%d}"
    return f"{a:%b %d}-{b:%b %d}"


def build_table(run_id, window, results_dir, obs_path, also, basin_km2):
    start = ts(RUN_STARTS[window])
    qpath, mpath, err = find_run_files(run_id, results_dir)
    if qpath is None:
        die(err)
    if mpath is None:
        die(f"no .mrf file found for {run_id} under {results_dir}. The rain comes from that file.")
    q = read_qout(qpath, start)
    mrf = read_mrf(mpath, start)
    gauge = read_obs(obs_path)
    M, G = Flow(q), Flow(gauge)

    rain_hr = mrf["MAP_mm_hr"].copy()
    rain_hr.index = rain_hr.index - HOUR            # the .mrf labels each hour by its END; use the hour's START
    water = mrf["MDGW_mm"]

    extra = {}
    for rid in also:
        qp, _m, e2 = find_run_files(rid, results_dir)
        if qp is None:
            die(f"--also run: {e2}")
        extra[rid] = Flow(read_qout(qp, ts(RUN_STARTS[run_window_from_id(rid) or window])))

    excl = [tuple(e) for e in EXCLUDE]
    rows = []
    for a, b in find_storms(rain_hr):
        w0, w1 = a, b + pd.Timedelta(days=1 + TAIL_DAYS)
        if w0 < start or w1 > q.index[-1] + HOUR:       # a storm the run does not fully cover is skipped
            continue
        r = rain_hr[(rain_hr.index >= w0) & (rain_hr.index < w1)]
        rain_mm = float(r.sum())
        rain_vol = rain_mm * basin_km2 / 1000.0          # mm over km2 -> M m3
        prior7 = float(rain_hr[(rain_hr.index >= w0 - pd.Timedelta(days=7)) & (rain_hr.index < w0)].sum())
        prior30 = float(rain_hr[(rain_hr.index >= w0 - pd.Timedelta(days=30)) & (rain_hr.index < w0)].sum())
        wt0 = float(water.asof(w0 + HOUR)) if len(water) else float("nan")
        gv = G.volume(w0, w1, excl) / 1e6 if G.has_data(w0, w1, excl) else float("nan")
        mv = M.volume(w0, w1, excl, match_to=G) / 1e6
        gq, mq = G.window(w0, w1, excl), M.window(w0, w1, excl)
        row = {
            "storm": storm_label(a, b), "period": "fall" if w0 < pd.Timestamp(FALL_END) else "winter",
            "window_start": w0, "window_end": w1,
            "rain_mm": rain_mm, "rain_volume_Mm3": rain_vol,
            "heaviest_hour_mm": float(r.max()) if len(r) else float("nan"),
            "rain_prior_7d_mm": prior7, "rain_prior_30d_mm": prior30, "water_table_at_start_mm": wt0,
            "gauge_vol_Mm3": gv, "gauge_peak_m3s": float(gq.max()) if gq.notna().any() else float("nan"),
            "model_vol_Mm3": mv, "model_peak_m3s": float(mq.max()) if mq.notna().any() else float("nan"),
        }
        row["gauge_runoff_ratio"] = gv / rain_vol if (rain_vol > 0 and math.isfinite(gv)) else float("nan")
        row["model_runoff_ratio"] = mv / rain_vol if rain_vol > 0 else float("nan")
        for rid, F in extra.items():
            tag = rid.split("_")[-1]
            ev = F.volume(w0, w1, excl, match_to=G) / 1e6
            row[f"{tag}_vol_Mm3"] = ev
            row[f"{tag}_runoff_ratio"] = ev / rain_vol if rain_vol > 0 else float("nan")
            fq = F.window(w0, w1, excl)
            row[f"{tag}_peak_m3s"] = float(fq.max()) if fq.notna().any() else float("nan")
        rows.append(row)
    df = pd.DataFrame(rows)
    keep = ((df["rain_mm"] >= MIN_STORM_MM) | (df["gauge_peak_m3s"] >= LIST_IF_FLOW_ABOVE) | (df["model_peak_m3s"] >= LIST_IF_FLOW_ABOVE)) if len(df) else []
    return (df[keep].reset_index(drop=True) if len(df) else df), {"qpath": qpath, "mpath": mpath, "run_first": q.index[0], "run_last": q.index[-1]}


# ------------------------------------------------------------------------------------------------------------------
# printing and plotting
# ------------------------------------------------------------------------------------------------------------------
def print_tables(df, run_id, extra_tags):
    log(f"\nSTORMS FROM THE RAIN  (run {run_id})")
    log(f"  A storm = rainy days (basin rain >= {RAIN_DAY_MM:g} mm/day), merged if separated by {MERGE_GAP_DAYS} dry days or fewer; "
        f"the window runs {TAIL_DAYS} days past the last rainy day to catch the recession.")
    log(f"  Listed if rain >= {MIN_STORM_MM:g} mm, or the gauge or model passes {LIST_IF_FLOW_ABOVE} m3/s. Jan 1-2 (gauge outage) is left out of both sides.")
    log("")
    hdr = (f"  {'storm':<14}{'rain':>7}{'peak':>7}{'prior7':>8}{'prior30':>9}{'wtable':>8} |{'gauge':>7}{'model':>7}{'g.peak':>8}{'m.peak':>8} |{'g/rain':>8}{'m/rain':>8}")
    for t in extra_tags:
        hdr += f"{t + ' vol':>12}{t + '/rain':>11}"
    log(hdr)
    log(f"  {'':<14}{'mm':>7}{'mm/hr':>7}{'mm':>8}{'mm':>9}{'mm':>8} |{'M m3':>7}{'M m3':>7}{'m3/s':>8}{'m3/s':>8} |{'ratio':>8}{'ratio':>8}")
    for period in ("fall", "winter"):
        part = df[df["period"] == period]
        if part.empty:
            continue
        log(f"  --- {period} ({'before' if period == 'fall' else 'from'} Jan 1) ---")
        for _, r in part.iterrows():
            line = (f"  {r['storm']:<14}{fmt(r['rain_mm'], 1):>7}{fmt(r['heaviest_hour_mm'], 1):>7}{fmt(r['rain_prior_7d_mm'], 1):>8}"
                    f"{fmt(r['rain_prior_30d_mm'], 1):>9}{fmt(r['water_table_at_start_mm'], 0):>8} |"
                    f"{fmt(r['gauge_vol_Mm3'], 3):>7}{fmt(r['model_vol_Mm3'], 3):>7}{fmt(r['gauge_peak_m3s'], 2):>8}{fmt(r['model_peak_m3s'], 2):>8} |"
                    f"{fmt(r['gauge_runoff_ratio'], 3):>8}{fmt(r['model_runoff_ratio'], 3):>8}")
            for t in extra_tags:
                line += f"{fmt(r[t + '_vol_Mm3'], 3):>12}{fmt(r[t + '_runoff_ratio'], 3):>11}"
            log(line)
    log("")
    log("  rain = basin-average rain over the storm window; peak = heaviest single hour; prior7 / prior30 = rain in the 7 / 30 days BEFORE the window;")
    log("  wtable = model mean depth to water table at the start (bigger = deeper); g/rain and m/rain = gauge or model outflow volume / rain volume.")
    log("  The gauge volume skips any blank stretch, so it can only be too small.")


def print_summary(df):
    log("\nSUMMARY (medians; storms with rain >= 10 mm only, so tiny events do not dominate)")
    big = df[df["rain_mm"] >= 10.0]
    for period in ("fall", "winter"):
        part = big[big["period"] == period]
        if part.empty:
            log(f"  {period}: no storm with 10 mm or more")
            continue
        log(f"  {period:<7} {len(part):>2} storms   median rain {fmt(part['rain_mm'].median(), 1)} mm   "
            f"median prior-7-day rain {fmt(part['rain_prior_7d_mm'].median(), 1)} mm   "
            f"model runoff ratio {fmt(part['model_runoff_ratio'].median(), 3)}   gauge runoff ratio {fmt(part['gauge_runoff_ratio'].median(), 3)}")
    fall = df[df["period"] == "fall"]
    if len(fall):
        worst = fall.sort_values("model_vol_Mm3", ascending=False).head(4)
        log(f"  The {len(worst)} fall storms with the most model outflow:")
        for _, r in worst.iterrows():
            log(f"    {r['storm']:<14} rain {fmt(r['rain_mm'], 1):>6} mm (heaviest hour {fmt(r['heaviest_hour_mm'], 1)} mm/hr)   "
                f"model {fmt(r['model_vol_Mm3'], 3)} M m3 (peak {fmt(r['model_peak_m3s'], 2)} m3/s)   gauge {fmt(r['gauge_vol_Mm3'], 3)} M m3")
    log("\nHOW TO READ IT")
    log("  * Fall storms with RAIN and heaviest-hour values like the winter storms, but a model runoff ratio near the winter one: the model treats fall and winter alike.")
    log("    The gauge does not, so something is different in the real basin in the fall (dry streambed, dry soil, rain that was smaller than the radar says).")
    log("  * Fall storms with MUCH BIGGER rain than any winter storm: the model may be right to make runoff, and the gauge is the surprise.")
    log("  * Fall storms with SMALL rain but a big model runoff ratio: the model is making runoff from too little rain (soil settings or the rain input).")
    log("  * Compare the water-table and prior-rain columns: if fall is much drier before the storm, the model is not using that dryness to hold water back.")
    log("  This is a pattern check, not a verdict: a basin average smooths a patchy storm, and each row is one storm.")


def make_figure(df, out_png, run_id, dpi, basin_km2=BASIN_KM2):
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    floor = 0.001                                   # flows below this (M m3), including zero, are drawn ON this line
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    colors = {"fall": "#c0392b", "winter": "#1f77b4"}
    d = df[df["rain_mm"] > 0]
    rmin, rmax = float(d["rain_mm"].min()), float(d["rain_mm"].max())
    xs = np.array([rmin * 0.8, rmax * 1.25])
    for frac, lab in ((0.01, "1% of the rain"), (0.1, "10%"), (0.5, "50%")):
        ys = xs * basin_km2 / 1000.0 * frac
        ax.plot(xs, ys, color="gray", lw=0.9, ls="--", zorder=1)
        ax.annotate(lab, (xs[1], ys[1]), xytext=(-4, -11), textcoords="offset points", ha="right", fontsize=8, color="gray")
    for period in ("fall", "winter"):
        p = d[d["period"] == period]
        if p.empty:
            continue
        mv = p["model_vol_Mm3"].clip(lower=floor)
        gv = p["gauge_vol_Mm3"].clip(lower=floor)
        ax.scatter(p["rain_mm"], mv, s=70, color=colors[period], edgecolor="black", linewidth=0.6, zorder=3, label=f"model, {period}")
        ax.scatter(p["rain_mm"], gv, s=70, marker="^", facecolor="white", edgecolor=colors[period], linewidth=1.8, zorder=3, label=f"gauge, {period}")
        for (_, r), y in zip(p.iterrows(), mv):
            if y >= 0.1 and (period == "fall" or y >= 1.0):
                ax.annotate(r["storm"], (r["rain_mm"], y), xytext=(6, 4), textcoords="offset points", fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(rmin * 0.8, rmax * 1.25)
    ax.set_ylim(floor / 2.5, max(float(d[["model_vol_Mm3", "gauge_vol_Mm3"]].max().max()) * 3.0, 1.0))
    ticks = [t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500) if rmin * 0.8 <= t <= rmax * 1.25]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:g}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:g}"))
    ax.set_xlabel("Basin-average rain in the storm window (mm)")
    ax.set_ylabel("Outflow volume in the storm window (million m3)")
    ax.set_title("Rain against outflow, storm by storm\nHigher on the page at the same rain = a larger share of the rain leaves as flow", fontsize=10.5)
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, frameon=False)
    fig.text(0.01, 0.005, f"run {run_id}\nDashed lines: outflow equal to 1%, 10% and 50% of the rain volume. "
             f"Flows under {floor:g} M m3 (including zero) are drawn on the bottom line.", fontsize=7, color="gray")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Rain and runoff ratio, storm by storm (see the text at the top of this file).",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", default=DEFAULT_RUN, help=f"run id whose rain and flow to use (default {DEFAULT_RUN})")
    p.add_argument("--window", default=None, help="W1..W5 (default: read from the run id)")
    p.add_argument("--also", default="", help="extra run ids (comma separated) whose outflow is added to the table, e.g. the Kchan runs")
    p.add_argument("--results-dir", default="results")
    p.add_argument("--obs", default=None, help=f"gauge file (default {DEFAULT_OBS})")
    p.add_argument("--out-dir", default="rain_runoff_out")
    p.add_argument("--basin-km2", type=float, default=None, help="basin area in km2 (default 424.6, as in the scorer)")
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--dry-run", action="store_true", help="find the files and the storms, print the table, write nothing")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)
    window = (a.window or run_window_from_id(a.run) or DEFAULT_WINDOW).upper()
    if window not in RUN_STARTS:
        die(f"unknown window {window!r}")
    obs = a.obs or DEFAULT_OBS
    basin = a.basin_km2 or BASIN_KM2
    also = [s.strip() for s in a.also.split(",") if s.strip()]
    results_dir = Path(a.results_dir)
    log(f"rain_runoff_storms_sycamore.py v{VERSION}   run {a.run}   window {window}   basin {basin:g} km2   gauge {obs}")
    df, info = build_table(a.run, window, results_dir, obs, also, basin)
    log(f"  outlet file {info['qpath']}\n  .mrf file   {info['mpath']}\n  run covers {info['run_first']:%b %d %Y %H:%M} to {info['run_last']:%b %d %Y %H:%M}")
    if df.empty:
        die("no storms found. Is the .mrf rain all zero, or does the run not cover October to March?")
    tags = [r.split("_")[-1] for r in also]
    print_tables(df, a.run, tags)
    print_summary(df)
    if a.dry_run:
        log("\nDRY RUN: nothing was written.")
        return 0
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv = out / "rain_runoff_storms.csv"
    png = out / "rain_vs_outflow_by_storm.png"
    df.to_csv(csv, index=False)
    make_figure(df, png, a.run, a.dpi, basin)
    log(f"\nWROTE {csv}\n      {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
