#!/usr/bin/env python
"""Figure: hourly basin-mean radar rain above the observed hydrograph for the five WY2017 winter events.

Inputs: init_data/obs/sycamore_09510200_obs_flow_WY2017.csv (repo, 15-min, MST)
        out/basin_mean_hourly_rain_wy2017.csv (from load_radar.py: MRMS grids area-weighted over the basin outline)
Rain bars: each bar is one hour, drawn from its start time to one hour later (Arizona time, MST).
"""
import argparse, textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patheffects as pe
import srp_common as c

ORANGE = "#eb6834"      # categorical slot 2: rain (flow is slot 1, blue)
EVENTS = c.EVENTS
HALO = [pe.withStroke(linewidth=3, foreground=c.SURFACE)]
SUBTITLE = "For the Jan 21, Jan 24 and Feb 28 storms, the wettest radar hour comes 4 to 7 hours before the flow peak"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(c.OUT / "fig1_rain_and_flow.png"))
    args = ap.parse_args()
    c.style()

    d = pd.read_csv(c.INIT / "obs" / "sycamore_09510200_obs_flow_WY2017.csv", parse_dates=["datetime"]).set_index("datetime")
    q = d.q_cms
    rain = pd.read_csv(c.OUT / "basin_mean_hourly_rain_wy2017.csv", parse_dates=["interval_start_mst"]).set_index("interval_start_mst").basin_mean_mm
    rmax = max(rain.loc[t0:t1].max() for _, _, t0, t1 in EVENTS)
    rtop = 8.0 if rmax <= 8 else float(np.ceil(rmax))

    peaks = {}
    for lab, name, t0, t1 in EVENTS:
        w = d.loc[t0:t1]; i = w.q_cms.idxmax()
        peaks[lab] = (i, w.q_cms[i], w.q_cfs[i])

    fig = plt.figure(figsize=(13, 10.2))
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.7], left=0.065, right=0.985, top=0.845, bottom=0.165, hspace=0.50)
    ax0 = fig.add_subplot(outer[0])
    inner = outer[1].subgridspec(2, 5, height_ratios=[0.55, 1], hspace=0.07, wspace=0.28)

    # ---- overview ----
    ax0.plot(q.index, q.values, color=c.BLUE, lw=1.3, solid_joinstyle="round")
    for lab, name, t0, t1 in EVENTS:
        ax0.axvspan(pd.Timestamp(t0), pd.Timestamp(t1), color=c.GRID, alpha=0.55, lw=0)
        t, v, _ = peaks[lab]
        ax0.annotate(lab, (t, v), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=10,
                     fontweight="bold", color=c.INK, path_effects=HALO)
    ax0.set_xlim(pd.Timestamp("2016-10-01"), pd.Timestamp("2017-10-01")); ax0.set_ylim(0, 240)
    ax0.xaxis.set_major_locator(mdates.MonthLocator()); ax0.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax0.set_ylabel("Discharge (m³/s)")
    ax0.set_xlabel("Water year 2017, October 2016 to September 2017 (Arizona time, MST)")
    ax0.grid(axis="y"); ax0.set_axisbelow(True)
    ax0.set_title("Whole year (shaded: the five events below)", loc="left", pad=6, fontsize=10, color=c.INK2, fontweight="normal")

    # ---- five events: rain row above flow row ----
    for k, (lab, name, t0, t1) in enumerate(EVENTS):
        axr = fig.add_subplot(inner[0, k])
        ax = fig.add_subplot(inner[1, k], sharex=axr)
        T0, T1 = pd.Timestamp(t0), pd.Timestamp(t1)

        rw = rain.loc[t0:t1]
        axr.bar(rw.index, rw.values, width=pd.Timedelta(hours=1), align="edge", color=ORANGE, lw=0)
        axr.set_ylim(rtop, 0)
        axr.set_yticks([0, 2, 4, 6][: int(rtop // 2) + 1])
        axr.grid(axis="y"); axr.set_axisbelow(True)
        axr.tick_params(labelbottom=False, labelsize=8.5)
        axr.set_title(f"{lab}  {name}", loc="left", pad=22)
        axr.text(0.99, 0.06, f"{rw.sum():.0f} mm total", transform=axr.transAxes, ha="right", va="bottom",
                 fontsize=8.5, color=c.INK2)
        if k == 0:
            axr.set_ylabel("Rain (mm/hr)", fontsize=9)

        w = d.loc[t0:t1]
        ax.plot(w.index, w.q_cms, color=c.BLUE, lw=1.6, solid_joinstyle="round")
        t, v, cfs = peaks[lab]
        ymax = w.q_cms.max() * 1.28
        ax.set_ylim(0, ymax); ax.set_xlim(T0, T1)
        ax.plot([t], [v], marker="o", ms=6, mfc=c.BLUE, mec=c.SURFACE, mew=1.5, ls="none", zorder=5)
        ax.xaxis.set_major_locator(mdates.DayLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[6, 12, 18]))
        ax.grid(axis="y"); ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelsize=8.5)
        plt.setp(ax.get_xticklabels(), ha="center")
        if lab == "1":
            axr.text(0, 1.045, "recession only; USGS-estimated", transform=axr.transAxes, fontsize=8.5, color=c.INK2, va="bottom")
            gap = d.index[d.tier == "missing"]; gap = gap[gap < "2017-01-03"]
            ax.axvspan(gap.min(), gap.max(), color=c.GRID, alpha=0.9, lw=0)
            ax.text(gap.min() + (gap.max() - gap.min()) / 2, ymax * 0.55, "USGS outage,\nno record", fontsize=8.5,
                    color=c.INK2, ha="center", va="center")
            txt = f"{v:.1f} m³/s\n({cfs:,.0f} cfs)\n{t:%b %d, %H:%M}"
            note = "first reading"
        else:
            axr.text(0, 1.045, f"flow peak {t:%b %d, %H:%M}", transform=axr.transAxes, fontsize=8.5, color=c.INK2, va="bottom")
            txt = (f"{v:.1f}" if v < 100 else f"{v:.0f}") + f" m³/s\n({cfs:,.0f} cfs)"
            note = None
        right = (t - T0) < pd.Timedelta(hours=40)
        ax.annotate(txt, (t, v), xytext=(9 if right else -9, 0), textcoords="offset points",
                    ha="left" if right else "right", va="center", fontsize=8.5, color=c.INK, path_effects=HALO)
        if k == 0:
            ax.set_ylabel("Discharge (m³/s)")

    fig.text(0.065, 0.958, "Rain and observed flow, Sycamore Creek near Fort McDowell (USGS 09510200), WY2017",
             fontsize=15, fontweight="bold", va="center")
    fig.text(0.065, 0.926, SUBTITLE, fontsize=10.5, color=c.INK2, va="center")
    paras = ["Rain: MRMS gauge-corrected hourly radar rain, area-weighted mean over the basin outline. Each bar is one hour of rain "
             "ending at the radar archive time stamp (assumed to mark the end of the hour); the rain scale is shared across the five panels, and \"mm total\" is the basin-mean rain summed over the window shown. "
             "Flow: USGS 15-minute record, Arizona time (MST); a break in the line is a gap with no record. Each flow panel has its own scale.",
             "Radar has a low bias at the gauges: in a hand check of four widespread storms, 72-hour radar totals were 69-87% of the mean of three FCDMC gauges."]
    fig.text(0.065, 0.012, "\n".join(textwrap.fill(p, 150) for p in paras), fontsize=8.5, color=c.INK2, va="bottom", linespacing=1.5)
    fig.savefig(args.out)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
