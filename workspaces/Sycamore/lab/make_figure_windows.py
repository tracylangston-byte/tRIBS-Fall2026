#!/usr/bin/env python3
"""Figure [WINDOWS]: calibration and scoring windows on the Sycamore Creek gauge record.

Two stacked panels:
  (a) Water year 2017: the six calibration windows (bars) under the gauge hydrograph,
      with the five scored storm windows shaded on the hydrograph.
  (b) Water year 2020: lead-in month, scored window, left-out summer, and the seven
      scored storm windows.

Usage:  python3 make_figure_windows.py WY2017.csv WY2020.csv OUTPUT_DIR
Inputs: gauge CSVs with columns datetime,q_cms,q_cfs,tier,note (15-minute record).
Outputs: Figure_WINDOWS.png (300 dpi) and Figure_WINDOWS.svg in OUTPUT_DIR.

Conventions: blue only for water; gauge symbol/line blue; timed storms vermillion
(Okabe-Ito #D55E00); volume-only storms and calibration bars gray; plain sans-serif
text at 8 pt or larger (Arial if installed, Liberation Sans otherwise).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

# ----------------------------------------------------------------------------
# Style
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.6,
    "svg.fonttype": "none",       # keep text as text so it can be restyled later
    "axes.unicode_minus": False,
    "hatch.linewidth": 0.45,
})

BLUE = "#0072B2"      # measured gauge flow (water)
LBLUE = "#56B4E9"     # estimated gauge flow (water, lighter)
VERM = "#D55E00"      # timed storms
GRAY = "#7F7F7F"      # volume-only storms
BAR = "#6B6B6B"       # calibration/run bars
BAR_LIGHT = "#C8C8C8"
INK = "#222222"
GRID = "#E3E3E3"

FLOOR = 1e-4          # zero flow is plotted here (m3/s) so it can appear on a log axis
YLIM = (6e-5, 3e3)

# ----------------------------------------------------------------------------
# Windows (from the Section 4.4.2 outline, D5, and the timing handoff, section 3)
# ----------------------------------------------------------------------------
T = pd.Timestamp

# (name, start, end, hours, start_state, extra)
CAL_WINDOWS = [
    ("Flood week",           T("2017-01-18"), T("2017-01-25"),       "168 h",   "0.82", ""),
    ("Late February",        T("2017-02-21"), T("2017-03-03"),       "240 h",   "0.82", ""),
    ("January to mid-March", T("2017-01-18"), T("2017-03-12"),       "1,272 h", "0.82", ""),
    ("October to January",   T("2016-10-01"), T("2017-01-25"),       "2,784 h", "0.95", "start-state test only"),
    ("Six-month",            T("2016-10-01"), T("2017-04-01"),       "4,368 h", "0.95", ""),
    ("Whole water year",     T("2016-10-01"), T("2017-09-30 23:00"), "8,759 h", "0.95", ""),
]

# (label, start, end, timed?) - scoring windows as used for the timing/shape scores
STORMS_2017 = [
    ("Dec 21–25",     T("2016-12-21"), T("2016-12-29"), False),
    ("Dec 30–Jan 2",  T("2016-12-30"), T("2017-01-06"), False),
    ("Jan 19–21",     T("2017-01-19"), T("2017-01-23"), True),
    ("Feb 18–19",     T("2017-02-18"), T("2017-02-23"), False),
    ("Feb 27–28",     T("2017-02-27"), T("2017-03-04"), True),
]
STORMS_2020 = [
    ("Nov 19–22",     T("2019-11-19"), T("2019-11-26"), False),
    ("Nov 28–30",     T("2019-11-28"), T("2019-12-04"), False),
    ("Dec 23–28",     T("2019-12-23"), T("2020-01-01"), True),
    ("Feb 10–11",     T("2020-02-10"), T("2020-02-15"), False),
    ("Feb 22–23",     T("2020-02-22"), T("2020-02-27"), True),
    ("Mar 10–13",     T("2020-03-10"), T("2020-03-17"), True),
    ("Mar 18–19",     T("2020-03-18"), T("2020-03-23"), True),
]

WY17_AXIS = (T("2016-10-01"), T("2017-10-01"))
WY20_AXIS = (T("2019-09-01"), T("2020-10-01"))
LEAD_IN = (T("2019-09-01"), T("2019-10-01"))
SCORED_20 = (T("2019-10-01"), T("2020-06-01"))       # 1 Oct 2019 to 31 May 2020 (end of the 31 May 23:00 step)
LEFT_OUT = (T("2020-06-01"), T("2020-10-01"))


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def load(path):
    d = pd.read_csv(path, parse_dates=["datetime"])
    q = d["q_cms"].astype(float)
    d["q_plot"] = np.where(q > 0, q, np.where(q.isna(), np.nan, FLOOR))
    return d


def draw_hydrograph(ax, d):
    """Blue line for the gauge record; estimated stretches in lighter blue; gaps break the line."""
    x = d["datetime"].to_numpy()
    y = d["q_plot"].to_numpy()
    ax.plot(x, y, color=BLUE, lw=0.6, solid_joinstyle="round", zorder=3)
    est = (d["tier"] == "estimated").to_numpy()
    # include one neighbor on each side so the lighter segment joins the darker one
    near = est | np.roll(est, 1) | np.roll(est, -1)
    ye = np.where(near, y, np.nan)
    ax.plot(x, ye, color=LBLUE, lw=0.6, zorder=4)
    # the lighter segments overdraw at most one 15-minute step at each end


def style_hydro_axis(ax, xlim, month_grid=True):
    ax.set_xlim(*xlim)
    ax.set_yscale("log")
    ax.set_ylim(*YLIM)
    ticks = [FLOOR, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]
    ax.set_yticks(ticks)
    ax.set_yticklabels(["0", "0.001", "0.01", "0.1", "1", "10", "100", "1000"])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_ylabel("Gauge discharge (m³/s)")
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", which="both", labelbottom=False, length=0)
    ax.tick_params(axis="y", length=2.5, width=0.6)
    if month_grid:
        add_month_grid(ax, xlim)


def add_month_grid(ax, xlim):
    months = pd.date_range(xlim[0], xlim[1], freq="MS")
    for m in months:
        ax.axvline(m, color=GRID, lw=0.5, zorder=0)


def month_labels(ax, xlim, year_at):
    """Month names centered between month starts; year printed under the months in year_at."""
    starts = pd.date_range(xlim[0], xlim[1], freq="MS")
    mids = [s + (s + pd.offsets.MonthBegin(1) - s) / 2 for s in starts[:-1]]
    labels = []
    for s in starts[:-1]:
        lab = s.strftime("%b")
        if (s.year, s.month) in year_at:
            lab += "\n" + str(s.year)
        labels.append(lab)
    ax.set_xticks(mdates.date2num(mids), minor=True)
    ax.set_xticklabels(labels, minor=True)
    ax.set_xticks(mdates.date2num(starts))
    ax.set_xticklabels([])
    ax.tick_params(axis="x", which="major", length=3, width=0.6)
    ax.tick_params(axis="x", which="minor", length=0, pad=2)


def shade_storms(ax, storms):
    for lab, a, b, timed in storms:
        ax.axvspan(a, b, color=VERM if timed else GRAY, alpha=0.20 if timed else 0.18,
                   lw=0, zorder=1)


def storm_labels(ax, storms, width_in, xlim, min_sep_in=0.14):
    """Rotated labels above the axes, nudged apart so they never touch, with short leaders."""
    span_days = (xlim[1] - xlim[0]).total_seconds() / 86400.0
    in_per_day = width_in / span_days
    centers = np.array([((a - xlim[0]) + (b - xlim[0])).total_seconds() / 2 / 86400.0 * in_per_day
                        for _, a, b, _ in storms])
    pos = centers.copy()
    for _ in range(400):
        moved = False
        for i in range(1, len(pos)):
            gap = pos[i] - pos[i - 1]
            if gap < min_sep_in - 1e-9:
                push = (min_sep_in - gap) / 2
                pos[i - 1] -= push
                pos[i] += push
                moved = True
        if not moved:
            break
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for (lab, a, b, timed), c_in, p_in in zip(storms, centers, pos):
        c = xlim[0] + pd.Timedelta(days=c_in / in_per_day)
        p = xlim[0] + pd.Timedelta(days=p_in / in_per_day)
        ax.plot([c, c, p], [1.0, 1.015, 1.035], color="#9A9A9A", lw=0.5, transform=tr,
                clip_on=False, zorder=5)
        ax.text(p, 1.04, lab, rotation=90, ha="center", va="bottom", fontsize=8, color=INK,
                transform=tr, clip_on=False)


def lane_axis(fig, rect, xlim, n_lanes):
    ax = fig.add_axes(rect)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.5, n_lanes - 0.5)
    ax.invert_yaxis()
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    add_month_grid(ax, xlim)
    ax.set_axisbelow(True)
    return ax


# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
def build(path17, path20, outdir):
    d17, d20 = load(path17), load(path20)

    W = 6.5
    L, R = 0.78, 0.12
    PW = W - L - R                     # plot width (in)

    # vertical layout, in inches, built from the bottom up
    x_lab = 0.40                       # month labels under the lane axis
    lane_h = 0.17                      # one lane
    hyd_h = 1.55
    lab_h = 0.80                       # room for rotated storm labels
    title_h = 0.24
    gap = 0.12

    y = 0.05
    b_lane_bot = y + x_lab
    b_lane_h = 0.30
    b_hyd_bot = b_lane_bot + b_lane_h + 0.05
    b_title_y = b_hyd_bot + hyd_h + lab_h + 0.06

    a_lane_bot = b_title_y + title_h + gap + x_lab
    a_lane_h = lane_h * len(CAL_WINDOWS)
    a_hyd_bot = a_lane_bot + a_lane_h + 0.05
    a_title_y = a_hyd_bot + hyd_h + lab_h + 0.06
    H = a_title_y + title_h + 0.02

    fig = plt.figure(figsize=(W, H), dpi=300)

    def rect(bottom, height):
        return [L / W, bottom / H, PW / W, height / H]

    # ---------------- panel (a)
    axa = fig.add_axes(rect(a_hyd_bot, hyd_h))
    shade_storms(axa, STORMS_2017)
    draw_hydrograph(axa, d17)
    style_hydro_axis(axa, WY17_AXIS)
    storm_labels(axa, STORMS_2017, PW, WY17_AXIS)

    # gauge outage 1 Jan 2017 (24.5 h missing), derived from the 'missing' tier
    miss = d17[d17["tier"] == "missing"]["datetime"]
    if len(miss):
        t0, t1 = miss.min(), miss.max() + pd.Timedelta(minutes=15)
        axa.axvspan(t0, t1, ymin=0, ymax=1, facecolor="none", edgecolor=INK, hatch="////",
                    lw=0.6, zorder=2)
        hrs = (t1 - t0).total_seconds() / 3600
        axa.annotate("gauge outage\n1 Jan\n(%.1f h)" % hrs, xy=(t0, 4e-3),
                     xytext=(-6, 0), textcoords="offset points", ha="right", va="center",
                     fontsize=8, color=INK, zorder=6)

    lanes_a = lane_axis(fig, rect(a_lane_bot, a_lane_h), WY17_AXIS, len(CAL_WINDOWS))
    for i, (name, s, e, hrs, start, extra) in enumerate(CAL_WINDOWS):
        lanes_a.broken_barh([(mdates.date2num(s), mdates.date2num(e) - mdates.date2num(s))],
                            (i - 0.31, 0.62), facecolors=BAR, lw=0, zorder=3)
        txt = "%s, %s, start %s" % (name, hrs, start)
        if extra:
            txt += " (%s)" % extra
        if name == "Whole water year":
            lanes_a.text(s + pd.Timedelta(days=3), i, txt, ha="left", va="center", fontsize=8,
                         color="white", zorder=4)
        else:
            lanes_a.text(e + pd.Timedelta(days=4), i, txt, ha="left", va="center", fontsize=8,
                         color=INK, zorder=4)
    month_labels(lanes_a, WY17_AXIS, year_at={(2016, 10), (2017, 1)})
    lanes_a.tick_params(axis="x", which="major", length=3, width=0.6)

    # ---------------- panel (b)
    axb = fig.add_axes(rect(b_hyd_bot, hyd_h))
    for span in (LEAD_IN, LEFT_OUT):      # unscored time: hatched, same as the lane below
        axb.axvspan(*span, facecolor="none", edgecolor="#B8B8B8", hatch="////", lw=0, zorder=0)
    shade_storms(axb, STORMS_2020)
    draw_hydrograph(axb, d20)
    style_hydro_axis(axb, WY20_AXIS)
    storm_labels(axb, STORMS_2020, PW, WY20_AXIS)
    # the 23 Sep 2019 storm falls in the lead-in month and is not scored
    axb.text(LEAD_IN[1] + pd.Timedelta(days=2), 60, "23 Sep storm\n(lead-in,\nnot scored)",
             ha="left", va="center", fontsize=8, color=INK, zorder=6, linespacing=1.15)

    lanes_b = lane_axis(fig, rect(b_lane_bot, b_lane_h), WY20_AXIS, 1)

    def seg(ax, a, b, fc, hatch=None, ec="none"):
        ax.broken_barh([(mdates.date2num(a), mdates.date2num(b) - mdates.date2num(a))],
                       (-0.38, 0.76), facecolors=fc, edgecolors=ec, hatch=hatch, lw=0.6, zorder=3)

    seg(lanes_b, *LEAD_IN, "white", hatch="////", ec="#9A9A9A")
    seg(lanes_b, *SCORED_20, BAR)
    seg(lanes_b, *LEFT_OUT, "white", hatch="////", ec="#9A9A9A")
    lanes_b.text(LEAD_IN[0] + (LEAD_IN[1] - LEAD_IN[0]) / 2, 0, "lead-in", ha="center",
                 va="center", fontsize=8, color=INK, zorder=4,
                 bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
    lanes_b.text(SCORED_20[0] + (SCORED_20[1] - SCORED_20[0]) / 2, 0,
                 "scored window, 1 Oct 2019 – 31 May 2020", ha="center", va="center",
                 fontsize=8, color="white", zorder=4)
    lanes_b.text(LEFT_OUT[0] + (LEFT_OUT[1] - LEFT_OUT[0]) / 2, 0, "left out: June 2020 Bush Fire",
                 ha="center", va="center", fontsize=8, color=INK, zorder=4,
                 bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
    month_labels(lanes_b, WY20_AXIS, year_at={(2019, 9), (2020, 1)})

    # ---------------- legend (a), upper right where the hydrograph is empty
    handles = [
        Line2D([0], [0], color=BLUE, lw=1.0, label="gauge, measured"),
        Line2D([0], [0], color=LBLUE, lw=1.0, label="gauge, estimated"),
        Patch(facecolor=VERM, alpha=0.30, lw=0, label="storm window scored for timing and volume"),
        Patch(facecolor=GRAY, alpha=0.30, lw=0, label="storm window scored for volume only"),
        Line2D([0], [0], color="none", label="0 on the axis = zero flow, plotted at 0.0001 m\u00b3/s"),
    ]
    axa.legend(handles=handles, loc="upper right", frameon=False, handlelength=1.6,
               borderaxespad=0.2, labelspacing=0.2)

    # ---------------- panel titles
    fig.text(0.012, (a_title_y + 0.02) / H,
             "(a)  Water year 2017: six calibration windows", fontsize=9, fontweight="bold",
             color=INK, va="bottom")
    fig.text(0.012, (b_title_y + 0.02) / H,
             "(b)  Water year 2020: one validation run (1 Sep 2019 start at 0.95, 6,575 h)",
             fontsize=9, fontweight="bold", color=INK, va="bottom")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "Figure_WINDOWS.png", dpi=300, facecolor="white")
    fig.savefig(outdir / "Figure_WINDOWS.svg", facecolor="white")
    print("figure size (in):", round(W, 2), "x", round(H, 2))
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    build(sys.argv[1], sys.argv[2], sys.argv[3])
