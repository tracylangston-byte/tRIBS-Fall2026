#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_w4_gauge_slowlimb_sycamore.py   (v1.0, 2026-10-31)
=======================================================

WHAT THIS SCRIPT DOES
---------------------
It draws the USGS 09510200 gauge record for the whole W4 window (Oct 1 2016 to Apr 1 2017) and shows where the
TENTATIVE CUT would fall: the part of each big storm's water that we are calling "slow flow" and might take out
before comparing the model to the gauge. It runs no model, changes none of your files, and only reads the gauge
file (and, if you ask, kept model outlet files to lay over the gauge).

VOCABULARY (the words used in the pictures)
-------------------------------------------
  hydrograph   a plot of flow (m3/s) against time
  recession    the falling part of a hydrograph after a peak
  tail         the long, slow end of a recession (flow of about 1 to 5 m3/s that lasts for days)
  slow flow    the part of a storm's flow that sits on the tail's line. ORANGE (an area on the linear picture, a
               line on the log pictures, because on a log scale an area down to the axis exaggerates it).
  quick flow   everything above the slow flow, the sharp flood part. Shaded BLUE.
  fit window   the days of the late tail that the line is fitted to. Shaded GRAY.
  half-life    how many hours it takes the slow flow to fall to half. Longer = a slower tail.

THE TENTATIVE RULE (the "line-fit rule", set in advance from the gauge alone, no model results used)
-----------------------------------------------------------------------------------------------------
  1. For each storm, take the days of its late tail (the fit window) and fit a straight line to the LOG of the flow
     against time. A straight line on a log scale is an exponential decay.
  2. Extend that line back to the storm's peak.
  3. Slow flow = the line, but never more than the flow the gauge actually recorded.
     Before the peak, slow flow = the flow just before the storm started (the starting base flow).
  4. Quick flow = recorded flow minus slow flow.
  Note the orange line STEPS UP at each peak (from the starting base flow to the line). That step is how the rule is
  written, not something the gauge measured. A smoother rule could ramp up over the rise; it would change the volume
  very little because the rise lasts hours, but it is a choice we should make on purpose.
  The picture's percentage is the slow flow's share of the storm's gauge volume over the storm window.

  This is a tentative rule. The size of the cut depends on which days you fit. For the January flood the script also
  draws two other fit windows so you can see how much the answer moves (about 14% to 28%).

PICTURES AND FILES IT WRITES (to --out-dir, default plot_out_w4_gauge/)
-----------------------------------------------------------------------
  w4_gauge_overview.png     top: the whole W4 gauge record (linear scale). Bottom: the same on a LOG scale, which
                            is the only way to see a 1 m3/s tail beside a 221 m3/s flood.
  w4_storm_zooms.png        one panel per storm (log scale): recorded flow, slow flow, quick flow, fit window, line.
  w4_storm_cut_summary.csv  the numbers behind the pictures.

LAYING MODEL RUNS OVER THE GAUGE (optional, later)
--------------------------------------------------
  python plot_w4_gauge_slowlimb_sycamore.py --overlay W4_ks0p2524_cv6_psib0p5743_f1p165_d099_nocl
  python plot_w4_gauge_slowlimb_sycamore.py --overlay RUN_A,RUN_B --results-dir results

  It looks for results/<run>_Outlet.qout (searched in sub-folders too) and draws it as a dashed colored line on every
  panel. A run id that starts with W1_ or W2_ is taken to start Jan 18 2017; W3_, W4_ or W5_ start Oct 1 2016.
  A run that is missing is skipped with a message; the gauge picture is still made.

HOW TO RUN IT (from workspaces/Sycamore/lab/)
---------------------------------------------
  python -u plot_w4_gauge_slowlimb_sycamore.py
  python -u plot_w4_gauge_slowlimb_sycamore.py --obs ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv
  python -u plot_w4_gauge_slowlimb_sycamore.py --out-dir my_plots --dpi 200

It takes a few seconds, runs in the foreground, and needs pandas, numpy and matplotlib (already in your environment).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")                      # draw to files only: no window is needed in a Codespace
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

VERSION = "1.0"
DEFAULT_OBS = "../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv"

W4_START = pd.Timestamp("2016-10-01 00:00")
W4_END = pd.Timestamp("2017-04-01 00:00")
OUTAGE = (pd.Timestamp("2017-01-01 00:00"), pd.Timestamp("2017-01-03 00:00"))   # gauge outage, Jan 1-2
FALSE_FLOW_THRESHOLD = 0.13        # m3/s: the gauge's largest Oct-Dec flow, the scorer's false-flow threshold
MIN_FIT_Q = 0.05                   # m3/s: flows below this are left out of the line fit
LOG_FLOOR = 0.05                   # m3/s: the bottom of the log-scale panels

# One entry per storm. All times are local Arizona time, [start, end).
#   vol  = the storm window the volume is counted over (no recession tail), same as the scorer's event windows
#   fit  = the late-tail days the line is fitted to
#   zoom = what the zoom panel shows
STORMS = {
    "jan": dict(label="January flood (main storm)",
                vol=("2017-01-18 00:00", "2017-01-23 12:00"),
                fit=("2017-01-26 00:00", "2017-02-03 00:00"),
                zoom=("2017-01-17 00:00", "2017-02-04 00:00")),
    "feb18": dict(label="Feb 18-19 storm",
                  vol=("2017-02-18 00:00", "2017-02-21 00:00"),
                  fit=("2017-02-22 12:00", "2017-02-26 12:00"),
                  zoom=("2017-02-17 00:00", "2017-02-27 00:00")),
    "feb27": dict(label="Feb 27-28 storm",
                  vol=("2017-02-27 00:00", "2017-03-03 00:00"),
                  fit=("2017-03-04 00:00", "2017-03-12 00:00"),
                  zoom=("2017-02-26 00:00", "2017-03-12 00:00")),
}
# Other fit windows for the January flood only, to show how much the cut depends on the days chosen.
ALT_JAN_FITS = {
    "fit Jan 27 - Feb 3": ("2017-01-27 00:00", "2017-02-03 00:00"),
    "fit Jan 26 - Feb 1": ("2017-01-26 00:00", "2017-02-01 00:00"),
    "fit Jan 28 - Feb 4": ("2017-01-28 00:00", "2017-02-04 00:00"),
}
ALT_STYLES = ["-.", ":", (0, (5, 2))]

RUN_STARTS = {"W1": "2017-01-18 00:00", "W2": "2017-01-18 00:00",
              "W3": "2016-10-01 00:00", "W4": "2016-10-01 00:00", "W5": "2016-10-01 00:00"}

C_GAUGE = "#1b3a5c"       # dark blue: the gauge
C_QUICK = "#8fb8de"       # light blue: quick flow
C_SLOW = "#e69f00"        # orange: slow flow
C_FIT = "#7f7f7f"         # gray: fit window
C_LINE = "#000000"        # black dashed: the fitted line
OVERLAY_COLORS = ["#d55e00", "#009e73", "#cc79a7", "#56b4e9", "#8c564b", "#7f7f7f"]

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
    """The part of a series with a <= time < b. (Plain text slicing like s['2017-02-18':'2017-02-21'] would include
    the whole of Feb 21, so this uses explicit comparisons instead.)"""
    a, b = ts(a), ts(b)
    return s[(s.index >= a) & (s.index < b)]


def hours_since_start(index) -> np.ndarray:
    return np.asarray((index - W4_START) / pd.Timedelta(hours=1), dtype=float)


def volume_m3(s: pd.Series) -> float:
    """Volume (m3) of a flow series (m3/s) by the trapezoid rule; blanks are dropped first."""
    s = s.dropna()
    if len(s) < 2:
        return 0.0
    t = np.asarray((s.index - s.index[0]).total_seconds(), dtype=float)
    return float(_trap(s.to_numpy(float), t))


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
    s = s[~s.index.duplicated(keep="last")].sort_index()
    return s


def read_qout(path: Path, start) -> pd.Series:
    """Routed outlet flow (m3/s) from a tRIBS *_Outlet.qout file: comma-delimited, a header line, Time_hr then Q."""
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


def run_start_from_id(run_id: str):
    prefix = run_id[:2].upper()
    return RUN_STARTS.get(prefix)


def load_overlays(arg: str | None, results_dir: str) -> list:
    """Returns [(run_id, series), ...] for the runs that could be read."""
    out = []
    if not arg:
        return out
    rd = Path(results_dir)
    for run_id in [r.strip() for r in arg.split(",") if r.strip()]:
        start = run_start_from_id(run_id)
        if start is None:
            print(f"  overlay skipped: {run_id}: cannot tell when it starts (run ids must begin W1_ ... W5_)")
            continue
        hits = sorted(rd.glob(f"**/{run_id}_Outlet.qout")) if rd.is_dir() else []
        if len(hits) != 1:
            print(f"  overlay skipped: {run_id}: expected exactly one {run_id}_Outlet.qout under {rd}, found {len(hits)}")
            continue
        out.append((run_id, read_qout(hits[0], start)))
        print(f"  overlay read: {run_id}  ({hits[0]})")
    return out


# ----------------------------------------------------------------------------------------------------------------
# the line-fit rule
# ----------------------------------------------------------------------------------------------------------------
def fit_line(gauge: pd.Series, fit_window) -> tuple:
    """Fit ln(flow) = c + k * t over the fit window (t in hours since Oct 1 2016). Returns (k, c, n points)."""
    f = between(gauge, *fit_window).dropna()
    f = f[f > MIN_FIT_Q]
    if len(f) < 10:
        raise ValueError(f"only {len(f)} usable gauge points in the fit window {fit_window}")
    k, c = np.polyfit(hours_since_start(f.index), np.log(f.to_numpy(float)), 1)
    return float(k), float(c), int(len(f))


def slow_flow(gauge: pd.Series, storm: dict, k: float, c: float, start, end):
    """Slow flow (m3/s) over [start, end) on the gauge's own times. Returns (slow series, gauge series, peak time,
    starting base flow, the extended line as a series)."""
    w = between(gauge, start, end).dropna()
    vol_part = between(gauge, *storm["vol"]).dropna()
    base = float(vol_part.iloc[:4].mean())                  # flow just as the storm window opens
    tpk = vol_part.idxmax()                                 # time of the storm's peak
    line = np.exp(c + k * hours_since_start(w.index))
    before = np.minimum(w.to_numpy(float), base)
    after = np.minimum(w.to_numpy(float), np.maximum(line, base))
    slow = pd.Series(np.where(w.index < tpk, before, after), index=w.index)
    return slow, w, tpk, base, pd.Series(line, index=w.index)


def storm_numbers(gauge: pd.Series, storm: dict, fit_window=None) -> dict:
    fw = fit_window or storm["fit"]
    k, c, n = fit_line(gauge, fw)
    slow, w, tpk, base, line = slow_flow(gauge, storm, k, c, *storm["zoom"])
    vol_gauge = volume_m3(between(w, *storm["vol"])) / 1e6
    vol_slow = volume_m3(between(slow, *storm["vol"])) / 1e6
    peak = float(between(w, *storm["vol"]).max())
    slow_pk = float(slow.asof(tpk))
    return dict(k=k, c=c, n_fit=n, slow=slow, gauge=w, tpk=tpk, base=base, line=line, fit_window=fw,
                vol_gauge=vol_gauge, vol_slow=vol_slow, share=vol_slow / vol_gauge if vol_gauge > 0 else np.nan,
                half_life_h=float(np.log(2) / -k) if k < 0 else np.nan, peak=peak, slow_at_peak=slow_pk,
                quick_peak=peak - slow_pk)


# ----------------------------------------------------------------------------------------------------------------
# pictures
# ----------------------------------------------------------------------------------------------------------------
def style_time_axis(ax):
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.grid(True, alpha=0.25, linewidth=0.6)


def draw_overview(gauge, nums, overlays, out_path, dpi):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8.5), sharex=True, gridspec_kw={"height_ratios": [1, 1]})
    g = gauge[(gauge.index >= W4_START) & (gauge.index < W4_END)]
    for ax, log in ((ax1, False), (ax2, True)):
        gp = g.where(g > 0) if log else g                 # a log axis cannot show zero flow
        ax.plot(gp.index, gp.values, color=C_GAUGE, linewidth=1.1, label="gauge (USGS 09510200)")
        ax.axvspan(*OUTAGE, color="#d62728", alpha=0.12, linewidth=0,
                   label="gauge outage (Jan 1-2, blank in the record)")
        for i, (name, nm) in enumerate(nums.items()):
            ax.axvspan(*[ts(x) for x in nm["fit_window"]], color=C_FIT, alpha=0.20, linewidth=0,
                       label="fit window (days the line is fitted to)" if i == 0 else None)
            s, w = nm["slow"], nm["gauge"]
            ax.fill_between(w.index, s.values, w.values,
                            color=C_QUICK, alpha=0.9, linewidth=0, label="quick flow (kept)" if i == 0 else None)
            if log:     # on a log scale an area down to the axis exaggerates the slow flow, so draw a line
                ax.plot(s.index, s.values, color=C_SLOW, linewidth=2.2,
                        label="slow flow (the tentative cut; a line on this log scale)" if i == 0 else None)
            else:       # on the linear scale the orange area is proportional to the volume that would be cut
                ax.fill_between(s.index, 0, s.values, color=C_SLOW, alpha=0.9, linewidth=0,
                                label="slow flow (the tentative cut)" if i == 0 else None)
        for j, (rid, q) in enumerate(overlays):
            qq = q[(q.index >= W4_START) & (q.index < W4_END)]
            qq = qq.where(qq > 0) if log else qq
            ax.plot(qq.index, qq.values, color=OVERLAY_COLORS[j % len(OVERLAY_COLORS)], linestyle="--",
                    linewidth=0.9, label=f"model: {rid}")
        ax.set_xlim(W4_START, W4_END)
        style_time_axis(ax)
    # top panel: label each storm at its peak
    label_pos = {"jan": (4, 225, "left"), "feb18": (-16, 125, "left"), "feb27": (30, 92, "right")}   # (days from the peak, flow, alignment)
    for name, nm in nums.items():
        dx, ly, ha = label_pos.get(name, (4, 100, "left"))
        ax1.annotate(f"{STORMS[name]['label']}\npeak {nm['peak']:.0f} m3/s; slow flow is {nm['share']:.0%} of\n"
                     f"the storm's volume", xy=(nm["tpk"], nm["peak"]),
                     xytext=(nm["tpk"] + pd.Timedelta(days=dx), ly), textcoords="data",
                     arrowprops=dict(arrowstyle="-", color="#555555", linewidth=0.7),
                     fontsize=8, va="center", ha=ha)
    ax1.set_ylim(0, 260)
    ax1.set_ylabel("flow (m3/s), linear scale")
    ax1.set_title("W4 gauge record (Oct 1 2016 to Apr 1 2017) with the tentative slow-flow cut", loc="left", fontsize=12)
    # bottom panel: log scale, the false-flow threshold, the fitted lines
    ax2.set_yscale("log")
    ax2.set_ylim(LOG_FLOOR, 400)
    ax2.axhline(FALSE_FLOW_THRESHOLD, color="#d62728", linestyle=":", linewidth=1.2)
    ax2.text(W4_START + pd.Timedelta(days=2), FALSE_FLOW_THRESHOLD * 1.15,
             f"{FALSE_FLOW_THRESHOLD} m3/s = largest gauge flow in Oct-Dec (the model makes false flow above this)",
             color="#d62728", fontsize=8, va="bottom")
    for i, (name, nm) in enumerate(nums.items()):
        after = nm["line"][nm["line"].index >= nm["tpk"]]
        ax2.plot(after.index, after.values, color=C_LINE, linestyle="--", linewidth=1.0,
                 label="fitted line, extended back to the peak" if i == 0 else None)
    ax2.set_ylabel("flow (m3/s), LOG scale (zero flow cannot be drawn)")
    ax2.set_xlabel("date (Arizona local time)")
    h, l = ax2.get_legend_handles_labels()
    ax2.legend(h, l, fontsize=8, loc="center left", ncol=1, framealpha=0.9)   # the empty Oct-Dec part of the plot
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def draw_zooms(gauge, nums, alt_jan, overlays, out_path, dpi):
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.2))
    for ax, (name, nm) in zip(axes, nums.items()):
        st = STORMS[name]
        z0, z1 = [ts(x) for x in st["zoom"]]
        ax.axvspan(*[ts(x) for x in st["vol"]], color=C_QUICK, alpha=0.10, linewidth=0,
                   label="storm window the volume is counted over")
        ax.axvspan(*[ts(x) for x in nm["fit_window"]], color=C_FIT, alpha=0.22, linewidth=0, label="fit window")
        g = between(gauge, z0, z1)
        ax.plot(g.index, g.where(g > 0).values, color=C_GAUGE, linewidth=1.3, label="gauge")
        s, w = nm["slow"], nm["gauge"]
        ax.fill_between(w.index, s.values, w.values, color=C_QUICK, alpha=0.9, linewidth=0, label="quick flow (kept)")
        ax.plot(s.index, s.values, color=C_SLOW, linewidth=2.4, label="slow flow (cut; a line on a log scale)")
        after = nm["line"][nm["line"].index >= nm["tpk"]]
        ax.plot(after.index, after.values, color=C_LINE, linestyle="--", linewidth=1.1, label="fitted line")
        if name == "jan":
            for (lab, fw), ls in zip(ALT_JAN_FITS.items(), ALT_STYLES):
                a = alt_jan[lab]
                ln = a["line"][a["line"].index >= a["tpk"]]
                ax.plot(ln.index, ln.values, color="#444444", linestyle=ls, linewidth=1.0,
                        label=f"{lab}: cut would be {a['share']:.0%}")
        for j, (rid, q) in enumerate(overlays):
            qq = between(q, z0, z1)
            ax.plot(qq.index, qq.where(qq > 0).values, color=OVERLAY_COLORS[j % len(OVERLAY_COLORS)],
                    linestyle="--", linewidth=1.1, label=f"model: {rid}")
        ax.set_yscale("log")
        ax.set_ylim(LOG_FLOOR, max(nm["peak"] * 2.5, 5))
        ax.set_xlim(z0, z1)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=2 if (z1 - z0).days > 12 else 1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.tick_params(axis="x", labelrotation=45, labelsize=8)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.set_title(f"{st['label']}\nslow flow = {nm['share']:.0%} of the storm's gauge volume "
                     f"({nm['vol_slow']:.2f} of {nm['vol_gauge']:.2f} M m3)\n"
                     f"tail half-life {nm['half_life_h']:.0f} h;  peak {nm['peak']:.0f} m3/s, "
                     f"slow flow at the peak {nm['slow_at_peak']:.1f}", fontsize=9, loc="left")
        ax.set_ylabel("flow (m3/s), LOG scale")
        ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    fig.suptitle("Each storm: what the line-fit rule would call slow flow (orange) and quick flow (blue)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Draw the W4 gauge record and the tentative slow-flow cut "
                                             "(see the text at the top of this file).")
    ap.add_argument("--obs", default=DEFAULT_OBS, help=f"gauge file (default {DEFAULT_OBS})")
    ap.add_argument("--out-dir", default="plot_out_w4_gauge", help="folder for the pictures and the CSV")
    ap.add_argument("--overlay", default=None, help="run id, or several separated by commas, to lay over the gauge")
    ap.add_argument("--results-dir", default="results", help="folder holding the kept *_Outlet.qout files (default results)")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args(argv)

    gauge = read_gauge(args.obs)
    n_blank = int(gauge.isna().sum())
    print(f"plot_w4_gauge_slowlimb_sycamore.py v{VERSION}")
    print(f"  gauge: {gauge.index[0]:%b %d %Y} to {gauge.index[-1]:%b %d %Y}, {int(gauge.notna().sum())} values, "
          f"{n_blank} blank")
    overlays = load_overlays(args.overlay, args.results_dir)

    nums = {}
    for name, st in STORMS.items():
        try:
            nums[name] = storm_numbers(gauge, st)
        except ValueError as exc:
            die(f"{st['label']}: {exc}")
    alt_jan = {}
    for lab, fw in ALT_JAN_FITS.items():
        alt_jan[lab] = storm_numbers(gauge, STORMS["jan"], fit_window=fw)

    print("\nTENTATIVE CUT (line-fit rule; volume = trapezoid rule over the storm window, M m3 = million m3)")
    print(f"  {'storm':<28}{'gauge vol':>10}{'slow vol':>10}{'share':>8}{'half-life h':>13}{'peak':>8}{'slow at peak':>14}{'quick peak':>12}")
    for name, nm in nums.items():
        print(f"  {STORMS[name]['label']:<28}{nm['vol_gauge']:>10.3f}{nm['vol_slow']:>10.3f}{nm['share']:>8.1%}"
              f"{nm['half_life_h']:>13.0f}{nm['peak']:>8.1f}{nm['slow_at_peak']:>14.1f}{nm['quick_peak']:>12.1f}")
    print("\nHOW MUCH THE JANUARY CUT DEPENDS ON THE DAYS FITTED")
    print(f"  {'fit window':<24}{'share':>8}{'slow vol':>10}{'slow at peak':>14}")
    print(f"  {'Jan 26 - Feb 3 (rule)':<24}{nums['jan']['share']:>8.1%}{nums['jan']['vol_slow']:>10.3f}{nums['jan']['slow_at_peak']:>14.1f}")
    for lab, a in alt_jan.items():
        print(f"  {lab:<24}{a['share']:>8.1%}{a['vol_slow']:>10.3f}{a['slow_at_peak']:>14.1f}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    draw_overview(gauge, nums, overlays, out_dir / "w4_gauge_overview.png", args.dpi)
    draw_zooms(gauge, nums, alt_jan, overlays, out_dir / "w4_storm_zooms.png", args.dpi)

    rows = []
    for name, nm in nums.items():
        rows.append(dict(storm=STORMS[name]["label"], fit_window=f"{nm['fit_window'][0]} to {nm['fit_window'][1]}",
                         storm_window=f"{STORMS[name]['vol'][0]} to {STORMS[name]['vol'][1]}",
                         gauge_vol_Mm3=round(nm["vol_gauge"], 4), slow_vol_Mm3=round(nm["vol_slow"], 4),
                         slow_share=round(nm["share"], 4), k_per_h=round(nm["k"], 6),
                         half_life_h=round(nm["half_life_h"], 1), gauge_peak_m3s=round(nm["peak"], 2),
                         slow_at_peak_m3s=round(nm["slow_at_peak"], 2), quick_peak_m3s=round(nm["quick_peak"], 2),
                         base_flow_m3s=round(nm["base"], 3), n_fit_points=nm["n_fit"], rule="main"))
    for lab, a in alt_jan.items():
        rows.append(dict(storm=STORMS["jan"]["label"], fit_window=f"{a['fit_window'][0]} to {a['fit_window'][1]}",
                         storm_window=f"{STORMS['jan']['vol'][0]} to {STORMS['jan']['vol'][1]}",
                         gauge_vol_Mm3=round(a["vol_gauge"], 4), slow_vol_Mm3=round(a["vol_slow"], 4),
                         slow_share=round(a["share"], 4), k_per_h=round(a["k"], 6),
                         half_life_h=round(a["half_life_h"], 1), gauge_peak_m3s=round(a["peak"], 2),
                         slow_at_peak_m3s=round(a["slow_at_peak"], 2), quick_peak_m3s=round(a["quick_peak"], 2),
                         base_flow_m3s=round(a["base"], 3), n_fit_points=a["n_fit"], rule=f"sensitivity: {lab}"))
    pd.DataFrame(rows).to_csv(out_dir / "w4_storm_cut_summary.csv", index=False)
    print(f"\nWROTE {out_dir / 'w4_gauge_overview.png'}\n      {out_dir / 'w4_storm_zooms.png'}\n"
          f"      {out_dir / 'w4_storm_cut_summary.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
