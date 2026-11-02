#!/usr/bin/env python3
"""
plot_routing_hydrographs_sycamore.py -- hydrographs from the Stage 5 routing sweep (rt1): the 60 runs at draw 12 and the 30 at draw 25.

Nothing is run or changed. It reads the flow files that already exist (results/<run>_Outlet.qout), the sweep summary (rt1_summary.csv)
and the gauge file, and draws pictures.

For EACH soil setting (draw 12, then draw 25) it makes two figures:

  hydro_draw<NN>_overview.png     four panels
      1. flood week Jan 18-25, linear flow axis. Gray = sweep runs that miss a closing criterion (arrival within 1 h, peak within 15%);
         blue = runs that meet BOTH; thick green = baseline routing (cv 6, flowexp 0.3, n 0.04); orange = k10 (the W1B winner);
         purple = k19 (the best-shaped passing run); black = gauge. The dotted lines mark the gauge peak +/-15%.
      2. the same week on a log flow axis, so the tail and the small flows show.
      3. zoom on the rising limb, Jan 20 18:00 - Jan 21 12:00. The dotted horizontal line is 20% of the gauge peak: the "arrival" is the
         time a curve first rises through that line, so a curve that crosses it earlier than the gauge's is early.
      4. cumulative volume (million m3), model and gauge.
  hydro_draw<NN>_by_parameter.png  three panels, all the rising-limb zoom, colored by cv, then flowexp, then n
      (one parameter per panel; gauge black, baseline green). It shows which knob moves the arrival time, the peak and the shape.

and ONE figure that compares the two soil settings:

  hydro_replication_12_vs_25.png    flood week and rising-limb zoom for the same routing points (baseline, k10, k19) at both soil settings:
      solid = draw 12, dashed = draw 25. The text under the title gives the largest difference in flow between the two.

ADDED (second storm and whole season), in the same folder routing_hydrographs/:

  hydro_w1b_draw<NN>.png   the W1B second-storm check (Feb 21 - Mar 3 runs; the Feb 27-28 storm is scored), one figure per soil setting:
      the 8 chosen sweep points + the baseline routing, from rt1b_summary.csv. Four panels: the storm (linear, with the gauge peak +/-15%),
      the whole run on a log axis, the rising-limb zoom (its time window is found from the gauge: 6 h before to 18 h after the gauge's
      own arrival; the dashed line is 20% of the gauge peak), and cumulative volume. Blue = passed BOTH storms (W1 and W1B); gray = did not;
      green = baseline routing; orange = k10; purple = k19; black = gauge.
  hydro_w4_draw12_season.png   the whole W4 season (Oct 1 2016 - Apr 1 2017) for draw 12: baseline routing (green), the k10 routing (orange)
      and the old W4 baseline (Ks 0.2, f 1, PsiB 0.5; thin gray dashed), against the gauge (black). Six panels: Oct-Dec (the gauge is dry, so any
      model flow is false flow; the legend gives each run's Oct-Dec volume), the whole season on a log axis, cumulative volume, the flood
      week Jan 18-25, Feb 18-21, and Feb 27 - Mar 3.
  hydro_w4_draw12_rising_limbs.png   three rising-limb zooms (flood week, Feb 18-19, Feb 27-28), one per storm. The window starts 3 h before the
      earliest arrival among the curves and runs to 15 h after the gauge's. The panel title gives each run's arrival error in hours
      (model minus gauge; negative = early), found with the scorer's own arrival function, so the numbers should match the scorer's table.

Run it from workspaces/Sycamore/lab/:
    python plot_routing_hydrographs_sycamore.py                       # everything below
    python plot_routing_hydrographs_sycamore.py --only sweep          # the W1 sweep figures only (also: w1b, w4; comma separated)
    python plot_routing_hydrographs_sycamore.py --anchors 12          # only draw 12 (sweep and W1B figures)
    python plot_routing_hydrographs_sycamore.py --label-k 10,19       # which sweep points are drawn as named, colored lines (default)

Needs, in the same folder: score_run_sycamore.py (its file readers are used), rt1_summary.csv, results/ with the sweep's Outlet.qout files,
and ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv. The W1B figures also need rt1b_summary.csv and the W1B runs in results/; the W4
figures need the three W4 flow files named below (options --w4-new, --w4-routing-baseline, --w4-baseline change the run ids).
Writes: routing_hydrographs/*.png (the folder is created if it is not there).
Prints what it drew and which files, if any, were missing. A run whose flow file is missing is left out of the picture, never invented; a group
whose summary or main files are missing is skipped with a message and the other groups still run.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import score_run_sycamore as sc

START = pd.Timestamp(sc.RUN_WINDOWS["W1"][0])
END = pd.Timestamp(sc.RUN_WINDOWS["W1"][1])
ZOOM = (pd.Timestamp("2017-01-20 18:00"), pd.Timestamp("2017-01-21 12:00"))
ARRIVAL_H, PEAK_PCT = 1.0, 15.0           # the closing criteria
LOG_FLOOR = 1e-3
C_FAIL, C_PASS, C_BASE = "#b8b8b8", "#3b7dc4", "#1b9e3e"
C_LABEL = {10: "#e8710a", 19: "#7b3fa0"}  # k10 orange, k19 purple; other labeled points get the gray-brown below
C_OTHER_LABEL = "#8c564b"
ANCHOR_TEXT = {12: "draw 12 (Ks x0.187, f x0.5874, PsiB x0.2939)", 25: "draw 25 (Ks x0.1852, f x0.5055, PsiB x0.3358)"}


def cum_volume(q: pd.Series, a, b) -> pd.Series:
    """Cumulative volume (million m3), trapezoid rule on the series' own time steps."""
    s = q[(q.index >= a) & (q.index <= b)].dropna()
    t = (s.index - s.index[0]).total_seconds().to_numpy()
    return pd.Series(sc.cumulative(t, s.to_numpy(float)) / 1e6, index=s.index)


def load_flow(results: Path, run_id: str, missing: list, start=None):
    p = results / f"{run_id}_Outlet.qout"
    if not p.is_file():
        missing.append(p.name)
        return None
    return sc.read_qout(p, START if start is None else start)


def clip(q: pd.Series, a, b) -> pd.Series:
    return q[(q.index >= a) & (q.index <= b)]


def setup_axis(ax, title, ylabel="flow (m3/s)", logy=False):
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)
    ax.tick_params(axis="x", rotation=30)
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(LOG_FLOOR, None)


def plot_series(ax, q, a, b, logy=False, **kw):
    s = clip(q, a, b)
    y = np.maximum(s.to_numpy(float), LOG_FLOOR) if logy else s.to_numpy(float)
    ax.plot(s.index, y, **kw)


def meets_both(df: pd.DataFrame) -> pd.Series:
    return (df["flood_week__arrival_err_h"].abs() <= ARRIVAL_H) & (df["flood_week__peak_err_pct"].abs() <= PEAK_PCT)


def line_style(row, label_ks):
    """(color, linewidth, alpha, zorder, name or None) for one summary row."""
    k = int(row["k"])
    if row["kind"] == "baseline":
        return C_BASE, 2.8, 1.0, 6, "baseline routing (cv 6, flowexp 0.3, n 0.04)"
    if k in label_ks:
        col = C_LABEL.get(k, C_OTHER_LABEL)
        return col, 2.2, 1.0, 7, f"k{k} (cv {row['cv']:g}, flowexp {row['flowexp']:g}, n {row['n']:g})"
    if row["meets"]:
        return C_PASS, 1.1, 0.75, 4, None
    return C_FAIL, 0.8, 0.6, 2, None


def draw_anchor(anchor, df, results, gauge, out_dir, label_ks, missing):
    df = df.copy()
    df["meets"] = meets_both(df)
    sweep = df[df["kind"] == "sweep"]
    n_pass, n_all = int(sweep["meets"].sum()), len(sweep)
    g = clip(gauge, START, END)
    gpk = float(g.max())
    level = 0.2 * gpk
    flows, n_drawn = {}, 0
    for _, r in df.iterrows():
        q = load_flow(results, r["run_id"], missing)
        if q is not None:
            flows[(r["kind"], int(r["k"]))] = q
            n_drawn += 1

    # ---------------- overview ----------------
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    a1, a2, a3, a4 = axes.ravel()
    for _, r in df.iterrows():
        q = flows.get((r["kind"], int(r["k"])))
        if q is None:
            continue
        col, lw, al, z, name = line_style(r, label_ks)
        for ax, (a, b, logy) in zip((a1, a2, a3), ((START, END, False), (START, END, True), (ZOOM[0], ZOOM[1], False))):
            plot_series(ax, q, a, b, logy=logy, color=col, lw=lw, alpha=al, zorder=z)
        c = cum_volume(q, START, END)
        a4.plot(c.index, c.to_numpy(), color=col, lw=lw, alpha=al, zorder=z)
    for ax, (a, b, logy) in zip((a1, a2, a3), ((START, END, False), (START, END, True), (ZOOM[0], ZOOM[1], False))):
        plot_series(ax, g, a, b, logy=logy, color="black", lw=2.0, zorder=9)
    gc = cum_volume(g, START, END)
    a4.plot(gc.index, gc.to_numpy(), color="black", lw=2.0, zorder=9)
    for ax in (a1, a3):
        ax.axhline(gpk * (1 + PEAK_PCT / 100), color="black", ls=":", lw=0.9, zorder=1)
        ax.axhline(gpk * (1 - PEAK_PCT / 100), color="black", ls=":", lw=0.9, zorder=1)
    a3.axhline(level, color="black", ls="--", lw=0.9, zorder=1)
    a3.text(ZOOM[0] + pd.Timedelta(minutes=20), level * 1.04, f"arrival level = 20% of the gauge peak ({level:.1f} m3/s)", fontsize=8)
    a1.text(0.01, 0.97, f"dotted: gauge peak {gpk:.0f} m3/s +/-{PEAK_PCT:g}%", transform=a1.transAxes, va="top", fontsize=8)
    setup_axis(a1, "Flood week Jan 18-25 (linear flow)")
    setup_axis(a2, "Flood week Jan 18-25 (log flow)", logy=True)
    setup_axis(a3, "Rising limb zoom, Jan 20 18:00 - Jan 21 12:00")
    setup_axis(a4, f"Cumulative volume, {START:%b %d} - {END:%b %d}", ylabel="million m3")
    handles = [Line2D([], [], color="black", lw=2.0, label=f"gauge (peak {gpk:.0f} m3/s)"),
               Line2D([], [], color=C_BASE, lw=2.8, label="baseline routing (cv 6, flowexp 0.3, n 0.04)")]
    for k in sorted(label_ks):
        rr = df[(df["kind"] == "sweep") & (df["k"] == k)]
        if len(rr) and ("sweep", k) in flows:
            r = rr.iloc[0]
            handles.append(Line2D([], [], color=C_LABEL.get(k, C_OTHER_LABEL), lw=2.2,
                                  label=f"k{k} (cv {r['cv']:g}, flowexp {r['flowexp']:g}, n {r['n']:g})"))
    handles += [Line2D([], [], color=C_PASS, lw=1.4, label=f"meets both criteria ({n_pass} of {n_all})"),
                Line2D([], [], color=C_FAIL, lw=1.4, label=f"misses a criterion ({n_all - n_pass} of {n_all})")]
    a1.legend(handles=handles, loc="upper right", fontsize=8)
    fig.suptitle(f"Stage 5 routing sweep on W1, {ANCHOR_TEXT.get(anchor, 'draw ' + str(anchor))}: {n_all} runs + baseline vs the gauge\n"
                 f"criteria: first arrival within {ARRIVAL_H:g} h and peak within {PEAK_PCT:g}%", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p1 = out_dir / f"hydro_draw{anchor:02d}_overview.png"
    fig.savefig(p1, dpi=120)
    plt.close(fig)

    # ---------------- by parameter ----------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for ax, (col, name, cmap) in zip(axes, (("cv", "cv (hillslope velocity coefficient)", "viridis"),
                                           ("flowexp", "flowexp (hillslope velocity exponent)", "plasma"),
                                           ("n", "n (channel roughness)", "cividis"))):
        vals = sweep[col].to_numpy(float)
        norm = matplotlib.colors.Normalize(vmin=vals.min(), vmax=vals.max())
        cm = plt.get_cmap(cmap)
        for _, r in sweep.iterrows():
            q = flows.get(("sweep", int(r["k"])))
            if q is not None:
                plot_series(ax, q, ZOOM[0], ZOOM[1], color=cm(norm(r[col])), lw=1.1, alpha=0.85, zorder=3)
        base = df[df["kind"] == "baseline"]
        if len(base) and ("baseline", 0) in flows:
            plot_series(ax, flows[("baseline", 0)], ZOOM[0], ZOOM[1], color=C_BASE, lw=2.8, zorder=6)
        plot_series(ax, g, ZOOM[0], ZOOM[1], color="black", lw=2.0, zorder=9)
        ax.axhline(level, color="black", ls="--", lw=0.9, zorder=1)
        setup_axis(ax, f"colored by {name}")
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cm), ax=ax, pad=0.02)
        cb.set_label(f"{col} (baseline: {float(base[col].iloc[0]):g})" if len(base) else col)
    axes[0].legend(handles=[Line2D([], [], color="black", lw=2.0, label="gauge"),
                            Line2D([], [], color=C_BASE, lw=2.8, label="baseline routing"),
                            Line2D([], [], color="black", ls="--", lw=0.9, label="arrival level (20% of gauge peak)")],
                   loc="upper left", fontsize=8)
    fig.suptitle(f"Rising limb, Jan 20 18:00 - Jan 21 12:00, {ANCHOR_TEXT.get(anchor, 'draw ' + str(anchor))}: one parameter per panel "
                 f"({len(sweep)} sweep runs, all three parameters change together)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p2 = out_dir / f"hydro_draw{anchor:02d}_by_parameter.png"
    fig.savefig(p2, dpi=120)
    plt.close(fig)

    print(f"  draw {anchor}: {n_drawn} of {len(df)} runs drawn; {n_pass} of {n_all} sweep runs meet both criteria "
          f"(k = {', '.join(str(int(k)) for k in sweep[sweep['meets']]['k'])})")
    print(f"    wrote {p1}")
    print(f"    wrote {p2}")
    return flows


def draw_replication(flows_by_anchor, summ, out_dir, label_ks, gauge):
    if not all(a in flows_by_anchor for a in (12, 25)):
        print("  replication figure skipped: it needs both draw 12 and draw 25")
        return
    g = clip(gauge, START, END)
    pts = [("baseline", 0, C_BASE, "baseline routing")] + [("sweep", k, C_LABEL.get(k, C_OTHER_LABEL), f"k{k}") for k in sorted(label_ks)]
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
    diffs = []
    for kind, k, col, name in pts:
        q12, q25 = flows_by_anchor[12].get((kind, k)), flows_by_anchor[25].get((kind, k))
        if q12 is None or q25 is None:
            print(f"  replication: {name} skipped (a flow file is missing)")
            continue
        for ax, (a, b) in zip(axes, ((START, END), ZOOM)):
            plot_series(ax, q12, a, b, color=col, lw=2.0, ls="-", zorder=5)
            plot_series(ax, q25, a, b, color=col, lw=2.0, ls="--", zorder=5)
        t = q12.index.union(q25.index)
        a_ = q12.reindex(t).interpolate(method="time", limit_area="inside")
        b_ = q25.reindex(t).interpolate(method="time", limit_area="inside")
        d = (a_ - b_).abs()
        diffs.append(f"{name}: largest difference {float(d.max()):.1f} m3/s ({100 * float(d.max()) / float(q12.max()):.1f}% of the draw-12 peak)")
    for ax, (a, b), title in zip(axes, ((START, END), ZOOM), ("Flood week Jan 18-25", "Rising limb zoom, Jan 20 18:00 - Jan 21 12:00")):
        plot_series(ax, g, a, b, color="black", lw=2.2, zorder=9)
        setup_axis(ax, title)
    gpk = float(g.max())
    axes[1].axhline(0.2 * gpk, color="black", ls="--", lw=0.9)
    handles = [Line2D([], [], color="black", lw=2.2, label="gauge"),
               Line2D([], [], color="gray", lw=2.0, ls="-", label="draw 12 (solid)"),
               Line2D([], [], color="gray", lw=2.0, ls="--", label="draw 25 (dashed)")]
    for kind, k, col, name in pts:
        handles.append(Line2D([], [], color=col, lw=2.4, label=name))
    axes[0].legend(handles=handles, loc="upper right", fontsize=9)
    fig.suptitle("Replication: the same routing points at the two soil settings (draw 12 solid, draw 25 dashed)\n" + "; ".join(diffs), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    p = out_dir / "hydro_replication_12_vs_25.png"
    fig.savefig(p, dpi=120)
    plt.close(fig)
    print(f"  wrote {p}")
    for d in diffs:
        print("    " + d)


# =============================================================================================================
# Second storm (W1B) and the whole season (W4)
# =============================================================================================================
B_START = pd.Timestamp("2017-02-21 00:00")                 # W1B runs start here (240 h; Feb 27-28 is scored)
W4_START = pd.Timestamp(sc.RUN_WINDOWS["W4"][0])
W4_END = pd.Timestamp(sc.RUN_WINDOWS["W4"][1])
EVENTS = sc.PROFILES["wy2017"]["events"]                   # the scorer's own storm windows
W4_NEW = "W4_ks0p187_cv6p39_psib0p2939_f0p5874_d012_fe0p2973_n0p067_nocl"
W4_ROUTING_BASELINE = "W4_ks0p187_cv6_psib0p2939_f0p5874_d012_nocl"
W4_BASELINE = "Sycamore_octmar_ks0p2_cv6_psib0p5_nocl"
C_NEW, C_OLD = "#e8710a", "#7f7f7f"


def arrival_info(gauge: pd.Series, name: str):
    """(level, gauge arrival time) for a scorer event window: level = 20% of the gauge peak in that window (as the scorer does)."""
    a, b = EVENTS[name]
    ev = gauge[(gauge.index >= pd.Timestamp(a)) & (gauge.index < pd.Timestamp(b))].dropna()
    level = 0.2 * float(ev.max())
    tg, flag = sc.arrival(ev, level)
    return level, (tg if flag == "ok" else None), float(ev.max())


def model_arrival(q: pd.Series, name: str, level: float):
    a, b = EVENTS[name]
    ev = q[(q.index >= pd.Timestamp(a)) & (q.index < pd.Timestamp(b))].dropna()
    tm, flag = sc.arrival(ev, level)
    return tm if flag == "ok" else None


def draw_w1b(anchor, df, results, gauge, out_dir, label_ks, missing):
    """One figure for the W1B check at one soil setting. df = the rt1b_summary rows of that anchor."""
    df = df.copy()
    df["meets"] = df["passes"].astype(str).str.lower().eq("true")
    sweep = df[df["kind"] == "sweep"]
    n_pass, n_all = int(sweep["meets"].sum()), len(sweep)
    ev_a, ev_b = (pd.Timestamp(x) for x in EVENTS["feb27_28"])
    end = ev_b
    level, tg, gpk = arrival_info(gauge, "feb27_28")
    zoom = (tg - pd.Timedelta(hours=6), tg + pd.Timedelta(hours=18)) if tg is not None else (ev_a, ev_a + pd.Timedelta(hours=24))
    g = clip(gauge, B_START, end)
    flows, n_drawn = {}, 0
    for _, r in df.iterrows():
        q = load_flow(results, r["run_id"], missing, start=B_START)
        if q is not None:
            flows[(r["kind"], int(r["k"]))] = q
            n_drawn += 1
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    a1, a2, a3, a4 = axes.ravel()
    spans = ((ev_a, end, False), (B_START, end, True), (zoom[0], zoom[1], False))
    for _, r in df.iterrows():
        q = flows.get((r["kind"], int(r["k"])))
        if q is None:
            continue
        col, lw, al, z, name = line_style(r, label_ks)
        for ax, (a, b, logy) in zip((a1, a2, a3), spans):
            plot_series(ax, q, a, b, logy=logy, color=col, lw=lw, alpha=al, zorder=z)
        c = cum_volume(q, B_START, end)
        a4.plot(c.index, c.to_numpy(), color=col, lw=lw, alpha=al, zorder=z)
    for ax, (a, b, logy) in zip((a1, a2, a3), spans):
        plot_series(ax, g, a, b, logy=logy, color="black", lw=2.0, zorder=9)
    gc = cum_volume(g, B_START, end)
    a4.plot(gc.index, gc.to_numpy(), color="black", lw=2.0, zorder=9)
    for ax in (a1, a3):
        ax.axhline(gpk * (1 + PEAK_PCT / 100), color="black", ls=":", lw=0.9, zorder=1)
        ax.axhline(gpk * (1 - PEAK_PCT / 100), color="black", ls=":", lw=0.9, zorder=1)
    a3.axhline(level, color="black", ls="--", lw=0.9, zorder=1)
    a3.text(zoom[0] + pd.Timedelta(minutes=20), level * 1.04, f"arrival level = 20% of the gauge peak ({level:.1f} m3/s)", fontsize=8)
    a1.text(0.01, 0.97, f"dotted: gauge peak {gpk:.0f} m3/s +/-{PEAK_PCT:g}%", transform=a1.transAxes, va="top", fontsize=8)
    setup_axis(a1, f"Storm window {ev_a:%b %d} - {ev_b:%b %d} (linear flow)")
    setup_axis(a2, f"Whole W1B run {B_START:%b %d} - {end:%b %d} (log flow)", logy=True)
    setup_axis(a3, f"Rising limb zoom, {zoom[0]:%b %d %H:%M} - {zoom[1]:%b %d %H:%M} (6 h before to 18 h after the gauge's arrival)")
    setup_axis(a4, f"Cumulative volume, {B_START:%b %d} - {end:%b %d}", ylabel="million m3")
    handles = [Line2D([], [], color="black", lw=2.0, label=f"gauge (peak {gpk:.0f} m3/s)"),
               Line2D([], [], color=C_BASE, lw=2.8, label="baseline routing (cv 6, flowexp 0.3, n 0.04)")]
    for k in sorted(label_ks):
        rr = df[(df["kind"] == "sweep") & (df["k"] == k)]
        if len(rr) and ("sweep", k) in flows:
            r = rr.iloc[0]
            handles.append(Line2D([], [], color=C_LABEL.get(k, C_OTHER_LABEL), lw=2.2,
                                  label=f"k{k} (cv {r['cv']:g}, flowexp {r['flowexp']:g}, n {r['n']:g})"))
    handles += [Line2D([], [], color=C_PASS, lw=1.4, label=f"passed both storms ({n_pass} of {n_all})"),
                Line2D([], [], color=C_FAIL, lw=1.4, label=f"did not ({n_all - n_pass} of {n_all})")]
    a1.legend(handles=handles, loc="upper right", fontsize=8)
    fig.suptitle(f"W1B second-storm check (Feb 27-28), {ANCHOR_TEXT.get(anchor, 'draw ' + str(anchor))}: the {n_all} points chosen from the sweep + baseline vs the gauge\n"
                 f"criteria: first arrival within {ARRIVAL_H:g} h and peak within {PEAK_PCT:g}%; passed = met both on W1 and on W1B", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / f"hydro_w1b_draw{anchor:02d}.png"
    fig.savefig(p, dpi=120)
    plt.close(fig)
    print(f"  W1B draw {anchor}: {n_drawn} of {len(df)} runs drawn; {n_pass} of {n_all} passed both storms (k = "
          f"{', '.join(str(int(k)) for k in sweep[sweep['meets']]['k'])})")
    print(f"    wrote {p}")


def volume_between(q: pd.Series, a, b) -> float:
    c = cum_volume(q, a, b)
    return float(c.iloc[-1]) if len(c) else float("nan")


def draw_w4(results, gauge, out_dir, ids: dict, missing):
    """The whole-season figure and the three rising-limb zooms for W4 at draw 12."""
    runs = [("base_old", ids["baseline"], C_OLD, 1.2, "--", "old W4 baseline (Ks 0.2, f 1, PsiB 0.5; cv 6)"),
            ("route0", ids["routing_baseline"], C_BASE, 2.2, "-", "draw 12, baseline routing (cv 6, 0.3, 0.04)"),
            ("new", ids["new"], C_NEW, 2.2, "-", "draw 12, k10 routing (cv 6.39, 0.2973, 0.067)")]
    flows = {}
    for key, rid, *_ in runs:
        q = load_flow(results, rid, missing, start=W4_START)
        if q is not None:
            flows[key] = q
    if not flows:
        print("  W4: none of the three W4 flow files were found, skipped")
        return
    g = clip(gauge, W4_START, W4_END)
    od_a, od_b = pd.Timestamp("2016-10-01"), pd.Timestamp("2017-01-01")

    def labeled(key, label):
        q = flows.get(key)
        return label if q is None else f"{label}; Oct-Dec {volume_between(q, od_a, od_b):.2f} M m3"

    fig, axes = plt.subplots(2, 3, figsize=(19, 10))
    a1, a2, a3, a4, a5, a6 = axes.ravel()
    wins = [(a1, od_a, od_b, False), (a2, W4_START, W4_END, True), (a4, *(pd.Timestamp(x) for x in EVENTS["flood_week"]), False),
            (a5, pd.Timestamp("2017-02-18"), pd.Timestamp("2017-02-21"), False), (a6, *(pd.Timestamp(x) for x in EVENTS["feb27_28"]), False)]
    for key, rid, col, lw, ls, label in runs:
        q = flows.get(key)
        if q is None:
            continue
        for ax, a, b, logy in wins:
            plot_series(ax, q, a, b, logy=logy, color=col, lw=lw, ls=ls, zorder=5)
        c = cum_volume(q, W4_START, W4_END)
        a3.plot(c.index, c.to_numpy(), color=col, lw=lw, ls=ls, zorder=5)
    for ax, a, b, logy in wins:
        plot_series(ax, g, a, b, logy=logy, color="black", lw=2.0, zorder=9)
    gc = cum_volume(g, W4_START, W4_END)
    a3.plot(gc.index, gc.to_numpy(), color="black", lw=2.0, zorder=9)
    for ax, name in ((a4, "flood_week"), (a6, "feb27_28")):
        _, _, gpk = arrival_info(gauge, name)
        ax.axhline(gpk * (1 + PEAK_PCT / 100), color="black", ls=":", lw=0.9, zorder=1)
        ax.axhline(gpk * (1 - PEAK_PCT / 100), color="black", ls=":", lw=0.9, zorder=1)
    setup_axis(a1, "Oct 1 - Dec 31: the gauge is dry, so any model flow is false flow")
    setup_axis(a2, "Whole season Oct 1 - Apr 1 (log flow)", logy=True)
    setup_axis(a3, "Cumulative volume Oct 1 - Apr 1", ylabel="million m3")
    setup_axis(a4, "Flood week Jan 18-25 (dotted: gauge peak +/-15%)")
    setup_axis(a5, "Feb 18-21")
    setup_axis(a6, "Feb 27 - Mar 3 (dotted: gauge peak +/-15%)")
    handles = [Line2D([], [], color="black", lw=2.0, label="gauge")]
    for key, rid, col, lw, ls, label in runs:
        if key in flows:
            handles.append(Line2D([], [], color=col, lw=lw, ls=ls, label=labeled(key, label)))
    a1.legend(handles=handles, loc="upper right", fontsize=8)
    fig.suptitle("W4 whole season (Oct 1 2016 - Apr 1 2017), draw 12 (Ks x0.187, f x0.5874, PsiB x0.2939; Kchan off): baseline routing vs the k10 routing", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "hydro_w4_draw12_season.png"
    fig.savefig(p, dpi=120)
    plt.close(fig)
    print(f"  wrote {p}")
    for key, rid, col, lw, ls, label in runs:
        if key in flows:
            print(f"    {labeled(key, label)}")

    # ---- rising limbs ----
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.5))
    for ax, (name, title) in zip(axes, (("flood_week", "Flood week"), ("feb18_19", "Feb 18-19"), ("feb27_28", "Feb 27-28"))):
        level, tg, gpk = arrival_info(gauge, name)
        if tg is None:
            ax.set_title(f"{title}: no gauge arrival found")
            continue
        arrivals = {}
        for key, rid, col, lw, ls, label in runs:
            q = flows.get(key)
            if q is not None:
                arrivals[key] = model_arrival(q, name, level)
        starts = [tg] + [t for t in arrivals.values() if t is not None]
        a, b = min(starts) - pd.Timedelta(hours=3), tg + pd.Timedelta(hours=15)
        for key, rid, col, lw, ls, label in runs:
            q = flows.get(key)
            if q is not None:
                plot_series(ax, q, a, b, color=col, lw=lw, ls=ls, zorder=5)
        plot_series(ax, g, a, b, color="black", lw=2.0, zorder=9)
        ax.axhline(level, color="black", ls="--", lw=0.9, zorder=1)
        ax.axvline(tg, color="black", ls=":", lw=0.9, zorder=1)
        parts = []
        for key, rid, col, lw, ls, label in runs:
            t = arrivals.get(key, None)
            if key in flows:
                parts.append(f"{label.split(' (')[0].replace('draw 12, ', '')}: " + (f"{(t - tg) / pd.Timedelta(hours=1):+.2f} h" if t is not None else "no arrival"))
        setup_axis(ax, f"{title}: arrival error (model minus gauge)\n" + "; ".join(parts))
    axes[0].legend(handles=[Line2D([], [], color="black", lw=2.0, label="gauge")] +
                   [Line2D([], [], color=col, lw=lw, ls=ls, label=label) for key, rid, col, lw, ls, label in runs if key in flows] +
                   [Line2D([], [], color="black", ls="--", lw=0.9, label="arrival level (20% of the gauge peak)"),
                    Line2D([], [], color="black", ls=":", lw=0.9, label="gauge arrival time")], loc="upper left", fontsize=7)
    fig.suptitle("W4 rising limbs, draw 12: where each curve rises through 20% of the gauge peak (the scorer's arrival)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = out_dir / "hydro_w4_draw12_rising_limbs.png"
    fig.savefig(p, dpi=120)
    plt.close(fig)
    print(f"  wrote {p}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Hydrographs from the Stage 5 routing sweep (see the top of this file).")
    ap.add_argument("--summary", default="rt1_summary.csv")
    ap.add_argument("--results", default="results")
    ap.add_argument("--obs", default="../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv")
    ap.add_argument("--out", default="routing_hydrographs")
    ap.add_argument("--anchors", default="12,25", help="draw numbers, comma separated (default 12,25)")
    ap.add_argument("--label-k", default="10,19", help="sweep points drawn as named, colored lines (default 10,19)")
    ap.add_argument("--only", default="sweep,w1b,w4", help="which groups to draw: sweep, w1b, w4 (comma separated; default all three)")
    ap.add_argument("--w1b-summary", default="rt1b_summary.csv", help="the W1B check's summary (default rt1b_summary.csv)")
    ap.add_argument("--w4-new", default=W4_NEW, help="run id of the W4 run with the k10 routing")
    ap.add_argument("--w4-routing-baseline", default=W4_ROUTING_BASELINE, help="run id of the W4 run at draw 12 with baseline routing")
    ap.add_argument("--w4-baseline", default=W4_BASELINE, help="run id of the old W4 baseline (Ks 0.2, f 1, PsiB 0.5)")
    a = ap.parse_args(argv)

    groups = {x.strip() for x in a.only.split(",") if x.strip()}
    bad = groups - {"sweep", "w1b", "w4"}
    if bad:
        print(f"--only: unknown group(s) {sorted(bad)}; choose from sweep, w1b, w4.")
        return 1
    try:
        anchors = [int(x) for x in a.anchors.replace(" ", "").split(",") if x]
        label_ks = {int(x) for x in a.label_k.replace(" ", "").split(",") if x}
    except ValueError:
        print("--anchors and --label-k must be numbers separated by commas.")
        return 1
    gauge = sc.read_obs(a.obs)
    out_dir = Path(a.out)
    out_dir.mkdir(exist_ok=True)
    results = Path(a.results)
    missing: list = []
    print(f"Hydrographs; flow files from {results}/; pictures go to {out_dir}/")

    if "sweep" in groups:
        print(f"\nThe W1 routing sweep ({a.summary})")
        if not Path(a.summary).is_file():
            print(f"  {a.summary} not found, sweep figures skipped. Run this from workspaces/Sycamore/lab/.")
        else:
            summ = pd.read_csv(a.summary).copy()
            need = ["anchor", "k", "kind", "cv", "flowexp", "n", "run_id", "flood_week__arrival_err_h", "flood_week__peak_err_pct"]
            miss = [c for c in need if c not in summ.columns]
            if miss:
                print(f"  {a.summary} lacks columns {miss}, sweep figures skipped.")
            else:
                summ = summ[summ["status"].eq("done")] if "status" in summ.columns else summ
                flows_by_anchor = {}
                for anchor in anchors:
                    df = summ[summ["anchor"] == anchor]
                    if df.empty:
                        print(f"  draw {anchor}: no rows in {a.summary}, skipped")
                        continue
                    flows_by_anchor[anchor] = draw_anchor(anchor, df, results, gauge, out_dir, label_ks, missing)
                if len(flows_by_anchor) == 2:
                    draw_replication(flows_by_anchor, summ, out_dir, label_ks, gauge)

    if "w1b" in groups:
        print(f"\nThe W1B second-storm check ({a.w1b_summary})")
        if not Path(a.w1b_summary).is_file():
            print(f"  {a.w1b_summary} not found, W1B figures skipped.")
        else:
            bs = pd.read_csv(a.w1b_summary).copy()
            need = ["anchor", "k", "kind", "cv", "flowexp", "n", "run_id", "passes"]
            miss = [c for c in need if c not in bs.columns]
            if miss:
                print(f"  {a.w1b_summary} lacks columns {miss}, W1B figures skipped.")
            else:
                for anchor in anchors:
                    df = bs[bs["anchor"] == anchor]
                    if df.empty:
                        print(f"  W1B draw {anchor}: no rows in {a.w1b_summary}, skipped")
                        continue
                    draw_w1b(anchor, df, results, gauge, out_dir, label_ks, missing)

    if "w4" in groups:
        print("\nW4, the whole season at draw 12")
        draw_w4(results, gauge, out_dir, {"new": a.w4_new, "routing_baseline": a.w4_routing_baseline, "baseline": a.w4_baseline}, missing)

    if missing:
        print(f"\n{len(missing)} flow file(s) were missing and those runs are left out: " + ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else ""))
    print("\nDone. Nothing was run or changed except the pictures in " + str(out_dir) + "/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
