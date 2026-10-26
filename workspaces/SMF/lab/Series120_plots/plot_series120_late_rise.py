"""
plot_series120_late_rise.py
===========================
Shows the small extra rise of flow, hours after the main flood, that some volume-matched Series 120 runs
make (about 07:00 to 08:00 on 13 Aug 2014) next to the RAIN the model was fed, so you can see for yourself
what comes before it.

Background. In plot_series120_hydrograph_shape.py the top row zooms to the main flood, so this small late rise
is left out of the picture (it is only mentioned in a footnote and in the runs table). This script is the picture
of the part that was left out.

It uses the same three storms and the same volume-matched runs as the hydrograph-shape script
(tRIBS 6.0.0, South Mountain Fan, synthetic truth Ks 7.0 / f 0.012):

    0.8x storm   (rain forcing multiplied by 0.8)
    1.0x storm   (rain as recorded)
    1.25x storm  (rain forcing multiplied by 1.25)

THE TWO FIGURES (both PNG)
    1  fig_series120_late_rise_1_night.png
       One column per storm, three rows, all on the same clock (21:00 on 12 Aug to 09:00 on 13 Aug):
         top     the RAIN the model was fed: both gauges (SMF and SMPHQ), in 15-minute blocks, already
                 multiplied by the storm's scale. The vertical axis is the same in every column.
         middle  the flow of every volume-matched run (colour = f group) and the synthetic truth (black).
                 The vertical axis is the same in every column, so you can see the late rise grow with the storm.
         bottom  the same flows with the vertical axis cut short so the quiet hours show: the overnight tail
                 that high-f runs hold a little above the truth, and the smaller humps before the late rise.
                 The late rise itself runs off the top of this row; that is on purpose.
       The grey band marks the strongest rain block after midnight (both gauges added together).

    2  fig_series120_late_rise_2_size_vs_f.png
       The same runs as numbers, one marker per run, by f:
         top     how far the run is above the truth at its worst after 05:00 on 13 Aug (m3/s)
         bottom  the extra water the run puts out after midnight compared with the truth, as a percent of the
                 truth's total volume (positive = more water after midnight than the truth)
       This figure draws the SIZE of the late rise as a continuous number, so it does not depend on any
       cut-off for deciding which runs "have" one.

It also writes two small tables you can open in Excel:
    series120_late_rise_runs.csv       one row per volume-matched run per storm
    series120_late_rise_summary.csv    one row per storm and f group (counts and medians)

WORDS USED
    late rise            how far a run's flow is above the truth's, at its largest, after 05:00 on 13 Aug.
                         A run with a late rise of at least 0.10 m3/s counts as "having" one (--bump_min).
                         That line is a choice. It decides the COUNTS ("14 of 16 runs") but not the figures.
    start of the rise    the lowest point of (run minus truth) in the two hours before the late rise peaks.
    extra volume         sum of (run minus truth) over every 5-minute step from midnight on, times 300 s,
                         shown as a percent of the truth's total volume.
    strongest rain block the 15-minute block after midnight with the most rain, both gauges added together.

WHAT THIS SCRIPT DOES NOT DO
    It does not run tRIBS and does not change any file you already have. It only READS the results CSVs, the
    hydrograph files and the two rain files, and WRITES PNGs and two small CSVs into its own folder. It is safe
    to run while a sweep is going.

    It also does not say WHY the late rise happens. It shows what is in the rain forcing and what the flow does.
    The reason (what the water is doing inside the model) needs a look at the model's own soil-moisture and
    runoff output.

HOW TO RUN (from the lab/ folder, like your other scripts)
    python plot_series120_late_rise.py
(or  python Series110_plots/plot_series120_late_rise.py  if you keep the scripts in that subfolder)
This file borrows helper functions from plot_series120_hydrograph_shape.py, so that file has to sit in the same
folder as this one.

WHERE THE INPUTS COME FROM
    calibration_work/03_comparisons/summary_tables/lhs_results_synth_Ks_f_120*.csv    results CSVs (one per storm)
    calibration_work/03_comparisons/csv_exports/<run_id>_compare_obs_sim.csv           the hydrographs
    init_data/met/precip_SMF_2014-2014.mdf and precip_SMPHQ_2014-2014.mdf              the rain forcing
and the output goes to
    calibration_work/03_comparisons/sensitivity_plots/series120_late_rise/

OPTIONS (all optional)
    --legs 080,100,125     which storms to use (default all three)
    --csv_080 PATH --csv_100 PATH --csv_125 PATH   use a different results file
    --summary_dir PATH     look for the results CSVs in a different folder
    --csv_dir PATH         look for the hydrograph files in a different folder
    --met_dir PATH         look for the two rain files in a different folder
    --out_dir PATH         save somewhere else
    --window_start "2014-08-12 21:00"   --window_end "2014-08-13 09:00"   the clock window drawn in figure 1
    --late_from "2014-08-13 05:00"      the late rise is looked for from this time on
    --bump_min 0.10        a run "has" a late rise if it gets at least this many m3/s above the truth
    --feasible_pct 2       "volume-matched" means within this many percent of the truth volume
    --f_low 0.007 --f_high 0.021       the f group cut-offs (same as the hydrograph-shape script)
    --watermark TEXT       a faint word across both figures (e.g. PRELIMINARY)
    --dpi 150              picture sharpness (72 to 600)

THINGS TO KEEP IN MIND WHEN READING IT
    * Synthetic, noise-free data: the truth is a model run, not a gauge. Real gauge error and rain-gauge error
      are not in these pictures.
    * The rain files hold a rate (mm/hr) four times an hour. A block is labelled by its START time.
    * The 0.8x and 1.25x storms use the rain multiplied by 0.8 and 1.25. A zero stays zero, so WHEN it rains is
      the same in all three storms; only the amounts change.
    * Two rain gauges are two points. SMPHQ sits outside the basin. They say what rain was in the forcing, not what
      fell on each part of the basin.
    * Along the zero-bias ridge, Ks rises as f rises, so f is also a stand-in for "where on the ridge".
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # draw to a file; no window needed (Codespaces has no screen)
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

try:
    import plot_series120_hydrograph_shape as S      # shared helpers: loaders, checks, colours, layout
except ImportError:
    sys.exit("This script borrows helper functions from plot_series120_hydrograph_shape.py.\n"
             "Put that file in the same folder as this one (the folder this file is in), then run it again.")

# ======================================================================
# SETTINGS you might want to change
# ======================================================================
WINDOW_START = "2014-08-12 21:00"
WINDOW_END = "2014-08-13 09:00"
LATE_FROM = "2014-08-13 05:00"      # look for the late rise from this time on
MIDNIGHT = "2014-08-13 00:00"       # "extra volume" and "strongest rain block" count from here
BUMP_MIN_M3S = 0.10                 # a run "has" a late rise if it gets at least this far above the truth
FOOT_LOOKBACK_MIN = 120             # the start of the rise = lowest (run minus truth) in this many minutes before its peak
BLOCK_MIN = 15                      # the rain files hold a rate four times an hour
NO_DATA = 9999                      # the forcing's NO_DATA flag is 9999.99
OUT_SUBFOLDER = "series120_late_rise"
FIG1 = "fig_series120_late_rise_1_night.png"
FIG2 = "fig_series120_late_rise_2_size_vs_f.png"
RUNS_CSV = "series120_late_rise_runs.csv"
SUMMARY_CSV = "series120_late_rise_summary.csv"
GAUGES = [("SMF", "precip_SMF_2014-2014.mdf"), ("SMPHQ", "precip_SMPHQ_2014-2014.mdf")]

# Rain bars use two neutral inks (a solid dark one and a light hatched one) so they cannot be mistaken for the
# blue / grey / red f groups; the hatch is a second cue besides lightness.
RAIN_STYLE = {
    "SMF":   dict(facecolor="#52514e", edgecolor="#52514e", hatch=None),
    "SMPHQ": dict(facecolor="#e4e2d9", edgecolor="#52514e", hatch="////"),
}
BLOCK_SHADE = "#c3c2b7"


class RainProblem(Exception):
    """The rain files cannot be used. Printed as a plain message."""


# ======================================================================
# Rain forcing
# ======================================================================
def read_forcing(path):
    """Rain rate (mm/hr) indexed by the start time of each block."""
    path = Path(path)
    if not path.exists():
        raise RainProblem(f"Could not find {path}")
    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]
    if "Year" not in df.columns:                       # an older, whitespace-delimited file
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = [str(c).strip() for c in df.columns]
    rcols = [c for c in df.columns if c.upper().startswith("R")]
    if not rcols or not {"Year", "Month", "Day", "Hour"} <= set(df.columns):
        raise RainProblem(f"{path.name}: expected the columns Year, Month, Day, Hour and R_mm; found {list(df.columns)}.")
    rate = pd.to_numeric(df[rcols[0]], errors="coerce")
    rate = rate.where(rate < NO_DATA, np.nan)
    grp = df.groupby(["Year", "Month", "Day", "Hour"], sort=False)
    minute = grp.cumcount() * 60.0 / grp["Year"].transform("size")
    t = pd.to_datetime(dict(year=df["Year"], month=df["Month"], day=df["Day"], hour=df["Hour"])) \
        + pd.to_timedelta(minute, unit="m")
    s = pd.Series(rate.to_numpy(float), index=pd.DatetimeIndex(t)).sort_index()
    return s[~s.index.duplicated(keep="first")]


def load_rain(met_dir, t0, t1):
    """{gauge: unscaled rain rate (mm/hr) on 15-minute blocks from t0 up to (not including) t1}, plus the block times."""
    blocks = pd.date_range(t0, t1 - pd.Timedelta(minutes=BLOCK_MIN), freq=f"{BLOCK_MIN}min")
    out = {}
    for g, fname in GAUGES:
        out[g] = read_forcing(Path(met_dir) / fname).reindex(blocks)
    return blocks, out


def strongest_block(rain, scale, after):
    """Start time of the 15-minute block after `after` with the most rain (both gauges added), or None."""
    tot = None
    for g in rain:
        v = (rain[g] * scale).fillna(0.0)
        tot = v if tot is None else tot + v
    sub = tot[tot.index >= after]
    if sub.empty or float(sub.max()) <= 0:
        return None
    return sub.idxmax()


# ======================================================================
# The late rise, per run
# ======================================================================
def analyse_run(idx, truth, sim, late_from, midnight, bump_min):
    """Numbers for one run. d = run minus truth, on the 5-minute grid."""
    d = sim - truth
    after = np.flatnonzero(np.asarray(idx >= late_from))
    if len(after) == 0:
        return None
    ip = int(after[int(np.argmax(d[after]))])                      # where the run is furthest above the truth
    lb = int(round(FOOT_LOOKBACK_MIN / S.STEP_MIN))
    lo = max(0, ip - lb)
    seg = d[lo:ip + 1]
    i_foot = lo + int(np.flatnonzero(seg == seg.min())[-1])        # the lowest point just before it
    day = np.asarray(idx >= midnight)
    step_s = S.STEP_MIN * 60.0
    extra = float(d[day].sum()) * step_s
    vol = float(truth.sum()) * step_s
    P = float(truth.max())
    return dict(
        late_excess_m3s=float(d[ip]),
        late_excess_pct_of_truth_peak=100.0 * float(d[ip]) / P,
        late_peak_flow_m3s=float(sim[ip]),
        late_peak_time=idx[ip],
        rise_start_time=idx[i_foot],
        rise_from_start_m3s=float(d[ip] - d[i_foot]),
        extra_volume_after_midnight_m3=extra,
        extra_volume_after_midnight_pct=100.0 * extra / vol,
        has_late_rise=bool(d[ip] >= bump_min),
        above_1pct_of_peak=bool(d[ip] >= 0.01 * P),
    )


# ======================================================================
# Figure 1: the night, rain over flow
# ======================================================================
def figure_night(keys, matched, hyd, data, rain, blocks, strongest, an, csv_dir, f_low, f_high, args, win, out_png):
    n = len(keys)
    W = max(5.4 * n + 1.3, S.MIN_FIG_W)
    left_in = 1.3
    left = left_in / W
    lo, hi = win
    late_from = pd.Timestamp(args.late_from)
    n_all = sum(len(matched[k]) for k in keys)
    n_ok = sum(len(hyd[k]["runs"]) for k in keys)
    warn = S.partial_text(keys, data)
    sub_text = S.wrap_text(
        f"Top: the rain the model was fed (both gauges, 15-minute blocks, already multiplied by the storm's scale). "
        f"Middle: flow of each volume-matched run (colour = f group) and the synthetic truth (black), {lo:%H:%M} on {lo:%d %b} to "
        f"{hi:%H:%M} on {hi:%d %b}. Bottom: the same flows with the vertical axis cut short so the quiet hours show; the late "
        f"rise runs off the top there. Grey band = the strongest rain block after midnight (both gauges added).", W, left_in, 11)
    notes = []
    rt = S.routing_text(data, keys)
    if rt:
        notes.append(rt)
    spreads = [hyd[k]["truth_spread"] for k in keys if hyd[k]["truth"] is not None]
    notes.append(f"{n_ok} of {n_all} hydrograph files read from {S._short(csv_dir)} and checked against the results CSVs "
                 f"(largest difference between the truth curves of one storm {(max(spreads) if spreads else 0.0):.1g} of the peak).")
    notes.append(f"Late rise = how far a run is above the truth at its largest after {late_from:%H:%M} on {late_from:%d %b}; "
                 f"a run with at least {args.bump_min:g} m³/s counts as having one. That line is a choice and sets the counts only; "
                 f"figure 2 shows the sizes.")
    notes.append("Rain bars are the forcing the run was given. Two gauge points (SMPHQ is outside the basin) say what rain was in the "
                 "forcing, not what fell on each part of the basin. Noise-free synthetic data: real gauge and rain-gauge error are not included.")
    notes_text = S.wrap_text("\n".join(notes), W, left_in, 9.6)
    top_in = S.top_block_in(sub_text, warn, 3)
    bottom_in, notes_h = S.bottom_block_in(notes_text, 9.6)
    H = top_in + 8.6 + bottom_in
    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor(S.SURFACE)
    gs = fig.add_gridspec(3, n, height_ratios=[0.8, 1.3, 1.0], left=left, right=0.985, top=1 - top_in / H,
                          bottom=bottom_in / H, hspace=0.42, wspace=0.20)

    # shared vertical ranges, so the three storms can be compared by eye
    rain_max = max([float(np.nanmax(np.nan_to_num(rain[g].to_numpy(float)) * S.LEGS[k]["scale"]))
                    for k in keys for g in rain] + [1.0])
    rain_lim = S.nice_ceil(rain_max * 1.12, 1.0)
    quiet_lo, quiet_hi = lo + pd.Timedelta(minutes=60), late_from + pd.Timedelta(minutes=60)
    mid_top, quiet_top = 0.2, 0.1
    for k in keys:
        h = hyd[k]
        if h["truth"] is None:
            continue
        idx = h["idx"]
        after = np.asarray((idx >= pd.Timestamp(MIDNIGHT)) & (idx <= hi))
        qmask = np.asarray((idx >= quiet_lo) & (idx <= quiet_hi))
        for s in list(h["runs"].values()) + [h["truth"]]:
            if after.any():
                mid_top = max(mid_top, float(s[after].max()))
            if qmask.any():
                quiet_top = max(quiet_top, float(s[qmask].max()))
    mid_lim = S.nice_ceil(mid_top * 1.12, 0.2)
    quiet_lim = S.nice_ceil(quiet_top * 1.15, 0.05)

    bw = 7.0 / 1440.0                                    # bar width in days (7 of the 15 minutes)
    xb = mdates.date2num(blocks.to_pydatetime())
    for col, k in enumerate(keys):
        spec, h, m = S.LEGS[k], hyd[k], matched[k]
        ax1 = fig.add_subplot(gs[0, col])
        ax2 = fig.add_subplot(gs[1, col], sharex=ax1)
        ax3 = fig.add_subplot(gs[2, col], sharex=ax1)
        for ax in (ax1, ax2, ax3):
            S.style_axes(ax)
        used = [i for i in h["runs"]]
        n_rise = {b: sum(1 for i in used if an[k][i]["has_late_rise"] and m.loc[i, "f_band"] == b) for b in S.BAND_ORDER}
        n_band = {b: sum(1 for i in used if m.loc[i, "f_band"] == b) for b in S.BAND_ORDER}
        tot_rise = sum(n_rise.values())
        line1 = f"{len(m)} runs within {args.feasible_pct:g}% of the truth volume" + \
                ("" if len(used) == len(m) else f"  ({len(used)} plotted)")
        line2 = f"late rise of {args.bump_min:g} m³/s or more: {tot_rise} of {len(used)} runs"
        line3 = f"low f {n_rise['low']} of {n_band['low']}, near {n_rise['near']} of {n_band['near']}, high f {n_rise['high']} of {n_band['high']}"
        S.panel_heading(ax1, f"{spec['name']}  ({spec['note']})", f"{line1}\n{line2}\n{line3}", sub_lines=3)

        # row 1: rain
        for g, off in (("SMF", 0.5 / 1440.0), ("SMPHQ", 7.5 / 1440.0)):
            v = np.nan_to_num(rain[g].to_numpy(float)) * spec["scale"]
            st = RAIN_STYLE[g]
            ax1.bar(xb + off, v, width=bw, align="edge", facecolor=st["facecolor"], edgecolor=st["edgecolor"],
                    hatch=st["hatch"], linewidth=0.5, zorder=3)
        ax1.set_ylim(0, rain_lim)
        S.set_time_axis(ax1, lo, hi)
        ax1.set_ylabel("Rain rate (mm/hr)", color=S.INK_SOFT, fontsize=11.5)

        if h["truth"] is None:
            for ax in (ax2, ax3):
                ax.set_xticks([])
                ax.set_yticks([])
                ax.text(0.5, 0.5, "no usable hydrograph files\nfor these runs", ha="center", va="center",
                        transform=ax.transAxes, color=S.INK_SOFT, fontsize=11)
            continue
        idx, truth = h["idx"], h["truth"]
        x = idx.to_numpy()
        order = sorted(h["runs"], key=lambda i: S.BAND_STYLE[m.loc[i, "f_band"]]["z"])
        for i in order:
            st = S.BAND_STYLE[m.loc[i, "f_band"]]
            for ax in (ax2, ax3):
                ax.plot(x, h["runs"][i], color=st["color"], linewidth=st["lw"], alpha=0.85, zorder=st["z"], solid_capstyle="round")
        for ax in (ax2, ax3):
            ax.plot(x, truth, color=S.INK, linewidth=1.5, zorder=6)
        for ax, lim in ((ax2, mid_lim), (ax3, quiet_lim)):
            S.set_time_axis(ax, lo, hi)
            ax.set_ylim(0, lim)
        ax2.set_ylabel("Flow (m³/s)", color=S.INK_SOFT, fontsize=11.5)
        ax3.set_ylabel("Flow, axis cut short (m³/s)", color=S.INK_SOFT, fontsize=11.5)
        ax3.set_xlabel(S.time_label(lo, hi), color=S.INK_SOFT, fontsize=10.5)

        # the strongest rain block after midnight, in all three rows
        t_blk = strongest[k]
        if t_blk is not None and lo <= t_blk < hi:
            for ax in (ax1, ax2, ax3):
                ax.axvspan(t_blk, t_blk + pd.Timedelta(minutes=BLOCK_MIN), color=BLOCK_SHADE, alpha=0.55, linewidth=0, zorder=1)
            ax2.annotate(f"strongest rain\nafter midnight\n({t_blk:%H:%M})", xy=(t_blk, 0.97), xycoords=("data", "axes fraction"),
                         xytext=(-6, 0), textcoords="offset points", ha="right", va="top", fontsize=9, color=S.INK_SOFT)

    S.fig_title(fig, "The small late rise, against the rain: tRIBS 6.0.0, Series 120", sub_text, left, warn=warn)
    handles = [Line2D([], [], color=S.INK, linewidth=2.0, label=f"synthetic truth (Ks {S.TRUTH_KS:g}, f {S.TRUTH_F:g})")]
    handles += S.band_handles(f_low, f_high, "line")
    for g in ("SMF", "SMPHQ"):
        st = RAIN_STYLE[g]
        handles.append(Patch(facecolor=st["facecolor"], edgecolor=st["edgecolor"], hatch=st["hatch"], linewidth=0.5,
                             label=f"rain, {g} gauge"))
    handles.append(Patch(facecolor=BLOCK_SHADE, alpha=0.55, linewidth=0, label="strongest rain block after midnight"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=10.5, labelcolor=S.INK_SOFT,
               bbox_to_anchor=(0.5, (0.15 + notes_h + 0.20) / H))
    fig.text(left, (0.15 + notes_h) / H, notes_text, fontsize=9.6, color=S.INK_SOFT, ha="left", va="top", linespacing=1.45)
    S.add_watermark(fig, args.watermark)
    S.save(fig, out_png, args.dpi)


# ======================================================================
# Figure 2: size of the late rise against f
# ======================================================================
def figure_size(keys, matched, hyd, data, an, f_low, f_high, args, out_png):
    n = len(keys)
    W = max(5.4 * n + 1.3, S.MIN_FIG_W)
    left_in = 2.0
    left = left_in / W
    warn = S.partial_text(keys, data)
    late_from = pd.Timestamp(args.late_from)
    sub_text = S.wrap_text(
        f"Each marker is one run within {args.feasible_pct:g}% of its storm's truth volume (the same runs as figure 1). Top: how far the run is above "
        f"the truth at its worst after {late_from:%H:%M} on {late_from:%d %b}. Bottom: the extra water the run puts out after midnight, as a percent of "
        f"the truth's total volume. Gold dashed line = the truth f ({S.TRUTH_F:g}).", W, left_in, 11)
    notes = [
        f"The dotted line in the top row is {args.bump_min:g} m³/s, the line used to count a run as having a late rise. A run that matches the truth exactly sits at zero "
        f"in both rows, because the truth's own small response after {late_from:%H:%M} (if it has one) is subtracted out.",
        "Volume-matched means the WHOLE event's volume is within the tolerance of the truth's. Extra water after midnight must therefore be balanced by less water "
        "before midnight: the run has about the same total, spread differently.",
        "Noise-free synthetic data. The dotted grey vertical lines are the f cut-offs between the three groups; they are choices, not measurements.",
    ]
    notes_text = S.wrap_text("\n".join(notes), W, left_in, 9.4)
    top_in = S.top_block_in(sub_text, warn, 2)
    bottom_in, notes_h = S.bottom_block_in(notes_text, 9.4)
    H = top_in + 8.2 + bottom_in
    rows = [("late_excess_m3s", "Largest amount above the truth\nafter 05:00 (m³/s)"),
            ("extra_volume_after_midnight_pct", "Extra volume after midnight\n(% of the truth's volume)")]
    fig, axes = plt.subplots(2, n, figsize=(W, H), sharex=True, sharey="row", squeeze=False)
    fig.patch.set_facecolor(S.SURFACE)
    fig.subplots_adjust(left=left, right=0.985, top=1 - top_in / H, bottom=bottom_in / H, hspace=0.12, wspace=0.07)
    for col, k in enumerate(keys):
        spec, m, h = S.LEGS[k], matched[k], hyd[k]
        used = list(h["runs"])
        table = pd.DataFrame({i: an[k][i] for i in used}).T if used else pd.DataFrame()
        counts = {b: int((m.loc[used, "f_band"] == b).sum()) for b in S.BAND_ORDER} if used else {b: 0 for b in S.BAND_ORDER}
        for r, (colname, ylabel) in enumerate(rows):
            ax = axes[r, col]
            S.style_axes(ax)
            ax.set_xscale("log")
            ax.set_xlim(0.003 / 1.07, 0.05 * 1.10)
            ax.set_xticks(S.F_TICKS)
            ax.set_xticklabels([f"{t:g}" for t in S.F_TICKS])
            ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            for edge in (f_low, f_high):
                ax.axvline(edge, color=S.AXIS, linewidth=1.0, linestyle=(0, (2, 3)), zorder=1)
            ax.axvline(S.TRUTH_F, color=S.GOLD, linewidth=1.6, linestyle=(0, (5, 3)), zorder=2)
            if r == 0:
                ax.axhline(args.bump_min, color=S.INK_SOFT, linewidth=1.0, linestyle=(0, (1, 2)), zorder=2)
            else:
                ax.axhline(0, color=S.INK, linewidth=1.1, zorder=2)
            for b in S.BAND_ORDER:
                ids = [i for i in used if m.loc[i, "f_band"] == b]
                if not ids:
                    continue
                st = S.BAND_STYLE[b]
                ax.scatter(m.loc[ids, "f_RS_abs"], table.loc[ids, colname].astype(float), s=58, marker=st["marker"],
                           color=st["color"], edgecolors=S.SURFACE, linewidths=0.7, zorder=4)
            if col == 0:
                ax.set_ylabel(ylabel, color=S.INK_SOFT, fontsize=10.5)
            else:
                ax.tick_params(labelleft=False)
            if r < 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("f (1/mm, log scale)", color=S.INK_SOFT, fontsize=11.5)
        S.panel_heading(axes[0, col], f"{spec['name']}  ({spec['note']})",
                        f"{len(used)} runs within {args.feasible_pct:g}% of the truth volume\n"
                        f"{counts['low']} low f, {counts['near']} near, {counts['high']} high f", sub_lines=2)
    S.fig_title(fig, "How big the late rise is, by f and storm size: tRIBS 6.0.0, Series 120", sub_text, left, warn=warn)
    handles = S.band_handles(f_low, f_high, "marker")
    handles.append(Line2D([], [], color=S.GOLD, linewidth=1.8, linestyle=(0, (5, 3)), label=f"truth f ({S.TRUTH_F:g})"))
    handles.append(Line2D([], [], color=S.INK_SOFT, linewidth=1.2, linestyle=(0, (1, 2)), label=f"{args.bump_min:g} m³/s line (counts as a late rise)"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=10.5, labelcolor=S.INK_SOFT,
               bbox_to_anchor=(0.5, (0.15 + notes_h + 0.20) / H))
    fig.text(left, (0.15 + notes_h) / H, notes_text, fontsize=9.4, color=S.INK_SOFT, ha="left", va="top", linespacing=1.45)
    S.add_watermark(fig, args.watermark)
    S.save(fig, out_png, args.dpi)


# ======================================================================
# Tables and console summary
# ======================================================================
def build_tables(keys, matched, hyd, an, strongest):
    rows = []
    for k in keys:
        m = matched[k]
        for i, r in m.iterrows():
            base = {"storm": S.LEGS[k]["name"], "run_id": r["run_id"], "Ks_mult": r["Ks_mult"], "f_RS_abs": r["f_RS_abs"],
                    "f_group": r["f_band"], "pbias_pct": r["pbias_pct"], "peak_error_pct": r["peak_error_pct"]}
            a = an[k].get(i)
            if a is None:
                base["hydrograph_file"] = hyd[k]["status"][i] or "not read"
            else:
                blk = strongest[k]
                lag = (a["late_peak_time"] - blk).total_seconds() / 60.0 if blk is not None else np.nan
                base.update(hydrograph_file="ok", **{c: a[c] for c in (
                    "has_late_rise", "above_1pct_of_peak", "late_excess_m3s", "late_excess_pct_of_truth_peak", "late_peak_flow_m3s")},
                    late_peak_time=str(a["late_peak_time"]), rise_start_time=str(a["rise_start_time"]),
                    rise_from_start_m3s=a["rise_from_start_m3s"],
                    extra_volume_after_midnight_m3=a["extra_volume_after_midnight_m3"],
                    extra_volume_after_midnight_pct=a["extra_volume_after_midnight_pct"],
                    strongest_rain_block_start=(str(blk) if blk is not None else ""),
                    minutes_from_strongest_rain_to_late_peak=lag)
            rows.append(base)
    runs = pd.DataFrame(rows)
    srows = []
    ok = runs[runs["hydrograph_file"] == "ok"]
    for (storm, grp), g in ok.groupby(["storm", "f_group"], sort=False):
        w = g[g["has_late_rise"].astype(bool)]
        srows.append({"storm": storm, "f_group": grp, "n_runs": len(g), "n_with_late_rise": len(w),
                      "n_above_1pct_of_peak": int(g["above_1pct_of_peak"].astype(bool).sum()),
                      "late_excess_m3s_median_all": g["late_excess_m3s"].median(),
                      "late_excess_m3s_max": g["late_excess_m3s"].max(),
                      "extra_volume_pct_median": g["extra_volume_after_midnight_pct"].median(),
                      "extra_volume_pct_max": g["extra_volume_after_midnight_pct"].max(),
                      "minutes_rain_to_late_peak_median_with_rise": w["minutes_from_strongest_rain_to_late_peak"].median() if len(w) else np.nan})
    summ = pd.DataFrame(srows)
    order = {b: j for j, b in enumerate(S.BAND_ORDER)}
    if len(summ):
        summ = summ.sort_values(["storm", "f_group"], key=lambda s: s.map(order) if s.name == "f_group" else s,
                                kind="stable").reset_index(drop=True)
    return runs, summ


# ======================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--legs", default="080,100,125", help="comma-separated storms to use (default 080,100,125)")
    for k in S.LEGS:
        ap.add_argument(f"--csv_{k}", type=Path, default=None)
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--csv_dir", type=Path, default=None, help="folder holding the <run_id>_compare_obs_sim.csv files")
    ap.add_argument("--met_dir", type=Path, default=None, help="folder holding the two precip_*.mdf rain files")
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--window_start", default=WINDOW_START)
    ap.add_argument("--window_end", default=WINDOW_END)
    ap.add_argument("--late_from", default=LATE_FROM)
    ap.add_argument("--bump_min", type=float, default=BUMP_MIN_M3S)
    ap.add_argument("--feasible_pct", type=float, default=S.FEASIBLE_PCT)
    ap.add_argument("--f_low", type=float, default=S.F_LOW)
    ap.add_argument("--f_high", type=float, default=S.F_HIGH)
    ap.add_argument("--watermark", default="", help="a faint word across both figures, e.g. PRELIMINARY")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")
    if not (0.1 <= args.feasible_pct <= 20.0):
        ap.error("--feasible_pct must be between 0.1 and 20.")
    if not (0.003 < args.f_low < args.f_high < 0.05):
        ap.error("--f_low and --f_high must satisfy 0.003 < f_low < f_high < 0.05.")
    if not (0.0 < args.bump_min <= 10.0):
        ap.error("--bump_min must be more than 0 and at most 10 (m3/s).")
    try:
        lo, hi = pd.Timestamp(args.window_start), pd.Timestamp(args.window_end)
        late_from = pd.Timestamp(args.late_from)
    except Exception:
        ap.error('The times must look like "2014-08-12 21:00".')
    if not (hi - lo >= pd.Timedelta(hours=3)) or not (lo < late_from < hi):
        ap.error("Need --window_end at least 3 hours after --window_start, and --late_from between them.")

    keys = [k.strip() for k in args.legs.split(",") if k.strip()]
    bad = [k for k in keys if k not in S.LEGS]
    if bad or len(set(keys)) != len(keys) or not keys:
        ap.error(f"--legs must be a comma-separated list drawn from {list(S.LEGS)} with no repeats (got {args.legs!r}).")
    keys = [k for k in S.LEGS if k in keys]

    root = Path.cwd().parent
    calib_dir = root / "calibration_work"
    summary_dir = args.summary_dir or (calib_dir / "03_comparisons" / "summary_tables")
    csv_dir = args.csv_dir or (calib_dir / "03_comparisons" / "csv_exports")
    met_dir = args.met_dir or (root / "init_data" / "met")
    out_dir = args.out_dir or (calib_dir / "03_comparisons" / "sensitivity_plots" / OUT_SUBFOLDER)
    out_res = out_dir.resolve()
    for guarded in (csv_dir, summary_dir, met_dir):
        if out_res == guarded.resolve():
            print(f"The output folder {out_dir} is the same as an input folder ({guarded}); refusing to write there.")
            raise SystemExit(1)

    paths = {k: (getattr(args, f"csv_{k}") or summary_dir / S.LEGS[k]["csv"]) for k in keys}
    missing = [p for p in paths.values() if not Path(p).exists()]
    if missing:
        print("Could not find these results file(s):")
        for p in missing:
            print(f"    {p}")
        print("Run this from the lab/ folder, or point to the file(s) with --csv_080 / --csv_100 / --csv_125 or --summary_dir.")
        raise SystemExit(1)

    print(f"Loading {len(keys)} storm(s) from {summary_dir}:")
    try:
        data = {k: S.load_leg(k, paths[k]) for k in keys}
    except S.DataProblem as exc:
        print(f"\nSTOPPED: {exc}")
        raise SystemExit(1)
    matched = {k: S.pick_matched(data[k], args.feasible_pct, args.f_low, args.f_high) for k in keys}
    print(f"\nVolume-matched runs (|PBIAS| < {args.feasible_pct:g}%):")
    for k in keys:
        c = matched[k]["f_band"].value_counts()
        print(f"  {S.LEGS[k]['name']:<12} {len(matched[k]):3d} runs   low f {int(c.get('low', 0))}, near {int(c.get('near', 0))}, high f {int(c.get('high', 0))}")
    empty = [S.LEGS[k]["name"] for k in keys if len(matched[k]) == 0]
    if empty:
        print(f"\nSTOPPED: no volume-matched runs for {', '.join(empty)}. Try a larger --feasible_pct.")
        raise SystemExit(1)

    # ---- the rain ----
    print(f"\nReading the rain forcing from {met_dir}:")
    try:
        blocks, rain = load_rain(met_dir, lo, hi)
    except RainProblem as exc:
        print(f"\nSTOPPED: {exc}\nRun this from the lab/ folder, or point to the folder with --met_dir.")
        raise SystemExit(1)
    for g in rain:
        s = rain[g]
        n_nan = int(s.isna().sum())
        if n_nan == len(s):
            print(f"\nSTOPPED: the {g} rain file has no values at all between {lo:%d %b %H:%M} and {hi:%d %b %H:%M}. "
                  f"Check --window_start and --window_end against the dates in the file.")
            raise SystemExit(1)
        wet = s[s > 0]
        print(f"  {g:<6} {float(np.nansum(s.to_numpy())) * BLOCK_MIN / 60.0:5.1f} mm in the window (unscaled record); "
              + (f"first rain block {wet.index[0]:%H:%M}, last {wet.index[-1]:%H:%M}" if len(wet) else "no rain in the window")
              + (f"; {n_nan} block(s) missing" if n_nan else ""))
    midnight = pd.Timestamp(MIDNIGHT)
    strongest = {k: strongest_block(rain, S.LEGS[k]["scale"], midnight) for k in keys}

    # ---- the hydrograph files ----
    print(f"\nReading the hydrograph files from {csv_dir}:")
    if not csv_dir.exists():
        print(f"\nSTOPPED: the folder {csv_dir} does not exist. Point to the right folder with --csv_dir.")
        raise SystemExit(1)
    hyd = {}
    for k in keys:
        hyd[k] = S.load_hydrographs(matched[k], csv_dir)
        st = pd.Series(hyd[k]["status"])
        n_ok = int((st == "ok").sum())
        n_nf = int((st == "file not found").sum())
        print(f"  {S.LEGS[k]['name']:<12} {n_ok:3d} of {len(st)} usable   ({n_nf} not found, {len(st) - n_ok - n_nf} left out)")
        shown = 0
        for i, s in enumerate(hyd[k]["status"]):
            if s not in ("ok", "file not found") and shown < 8:
                print(f"      {matched[k].loc[i, 'run_id']}: {s}")
                shown += 1
    if sum(len(hyd[k]["runs"]) for k in keys) == 0:
        print("\nSTOPPED: no usable hydrograph files were found, so there is nothing to plot.")
        print("  Check where your <run_id>_compare_obs_sim.csv files are and pass that folder with --csv_dir.")
        raise SystemExit(1)
    for k in keys:
        if hyd[k]["truth"] is not None and not (hyd[k]["idx"][0] <= lo and hyd[k]["idx"][-1] >= late_from):
            print(f"\nSTOPPED: the hydrographs of the {S.LEGS[k]['name']} run from {hyd[k]['idx'][0]:%d %b %H:%M} to "
                  f"{hyd[k]['idx'][-1]:%d %b %H:%M}, which does not cover the window you asked for. Check --window_start and --late_from.")
            raise SystemExit(1)

    # ---- the numbers ----
    an = {k: {} for k in keys}
    for k in keys:
        h = hyd[k]
        if h["truth"] is None:
            continue
        for i, sim in h["runs"].items():
            a = analyse_run(h["idx"], h["truth"], sim, late_from, midnight, args.bump_min)
            if a is None:
                print(f"\nSTOPPED: {S.LEGS[k]['name']}: no data after {late_from:%d %b %H:%M}. Check --late_from.")
                raise SystemExit(1)
            an[k][i] = a

    runs, summ = build_tables(keys, matched, hyd, an, strongest)

    out_dir.mkdir(parents=True, exist_ok=True)
    figure_night(keys, matched, hyd, data, rain, blocks, strongest, an, csv_dir, args.f_low, args.f_high, args, (lo, hi), out_dir / FIG1)
    figure_size(keys, matched, hyd, data, an, args.f_low, args.f_high, args, out_dir / FIG2)
    runs.to_csv(out_dir / RUNS_CSV, index=False)
    summ.to_csv(out_dir / SUMMARY_CSV, index=False)

    # ---- console summary ----
    print("\nStrongest rain block after midnight (both gauges added; the rain as the model was fed, i.e. scaled for the storm):")
    for k in keys:
        b = strongest[k]
        if b is None:
            print(f"  {S.LEGS[k]['name']:<12} no rain after midnight in the window")
            continue
        parts = ", ".join(f"{g} {float(rain[g][b]) * S.LEGS[k]['scale']:.1f} mm/hr" for g in rain if pd.notna(rain[g][b]))
        print(f"  {S.LEGS[k]['name']:<12} block starting {b:%H:%M}: {parts}")
    print(f"\nLate rise of the volume-matched runs (a run 'has' one if it gets at least {args.bump_min:g} m3/s above the truth after {late_from:%H:%M}):")
    show = summ[["storm", "f_group", "n_runs", "n_with_late_rise", "n_above_1pct_of_peak", "late_excess_m3s_max",
                 "extra_volume_pct_median", "minutes_rain_to_late_peak_median_with_rise"]].copy()
    show.columns = ["storm", "f group", "runs", f"rise >= {args.bump_min:g}", "rise >= 1% of peak", "biggest rise m3/s",
                    "median extra vol %", "min: rain to late peak"]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print("\n  'rise >= 1% of peak' is the line the hydrograph-shape script's footnote used; the counts depend on the line, which is why figure 2 shows sizes.")
    print("  'extra vol %' = run minus truth, summed from midnight on, as a percent of the truth's volume; 'min: rain to late peak' = minutes from the start")
    print("  of the strongest rain block after midnight to the run's late peak (median over the runs that have a late rise).")
    print(f"\nSaved to {out_dir}:")
    for w in (FIG1, FIG2, RUNS_CSV, SUMMARY_CSV):
        print(f"    {w}")
    print("\nREMINDER: synthetic, noise-free data, and two rain-gauge points. This shows what is in the rain forcing and what the flow does; "
          "it does not show why the water leaves the model late.")


if __name__ == "__main__":
    main()
