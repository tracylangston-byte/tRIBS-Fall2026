#!/usr/bin/env python
"""Slide: event table for the five WY2017 winter storms (rain, peak flow, lag, runoff, runoff ratio) + whole year.

Needs out/basin_mean_hourly_rain_wy2017.csv (run load_radar.py first). Numbers come from event_table.py.
Writes out/fig5_event_table.png (16:9 slide) and out/event_table_wy2017.csv.
"""
import argparse, textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import srp_common as c
import event_table as E

ORANGE = "#eb6834"      # rain (categorical slot 2); flow is blue (slot 1)


def fmt_q(v):
    return f"{v:.1f}" if v < 100 else f"{v:,.0f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(c.OUT / "fig5_event_table.png"))
    args = ap.parse_args()
    c.style()
    ev, yr = E.compute()
    ev.to_csv(c.OUT / "event_table_wy2017.csv", index=False)

    # monthly rain and flow (daily means) for the subtitle
    q, r, cfs = E.load()
    m3 = cfs * E.CFS_TO_CMS * 86400
    flow_m = m3.resample("MS").sum()
    jfm_flow = flow_m.iloc[3:6].sum() / flow_m.sum() * 100
    jfm_rain = r.resample("MS").sum().iloc[3:6].sum() / r.sum() * 100
    share2 = ev.runoff_mm[1] / yr["flow_mm"] * 100
    # how much the ratios move if the end of direct runoff is cut off 1.5 to 5 days after the peak
    sens = {k: [] for k in ("2", "5")}
    for nd in (1.5, E.N_DAYS, 3.5, 5.0):
        e2, _ = E.compute(nd)
        for k in sens:
            sens[k].append(float(e2.runoff_ratio[e2.storm == k].iloc[0]) * 100)
    rng = {k: (min(v), max(v)) for k, v in sens.items()}
    shares = [E.compute(nd)[0].runoff_mm[1] / yr["flow_mm"] * 100 for nd in (1.5, E.N_DAYS, 3.5, 5.0)]

    W, H = 13.333, 7.5
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # ---- headline ----
    head = f"The Jan 21 storm peaked at {yr['peak_cfs']:,.0f} cfs and carried about a quarter of the year's flow"
    ax.text(0.05, 0.93, "\n".join(textwrap.wrap(head, 66)), fontsize=23, fontweight="bold", va="center", linespacing=1.25)
    sub = (f"Over the water year, about {yr['ratio']*100:.0f}% of the rain became streamflow. January to March brought "
           f"{jfm_rain:.0f}% of the rain and {jfm_flow:.0f}% of the flow; July brought 24% of the rain and almost none.")
    ax.text(0.05, 0.815, "\n".join(textwrap.wrap(sub, 125)), fontsize=13, color=c.INK2, va="center", linespacing=1.4)

    # ---- columns ----
    X = dict(storm=0.05, rain_n=0.225, rain_b=0.235, wet=0.375, q_n=0.485, q_b=0.495, lag=0.645, ro_n=0.795, ro_b=0.805, ratio=0.925)
    y_head, y0, dy = 0.725, 0.645, 0.077
    heads = [("storm", "Storm", "left"), ("rain_n", "Rain (mm)\nbasin mean", "right"), ("wet", "Wettest hour\n(mm/hr)", "center"),
             ("q_n", "Peak flow\n(m³/s)", "right"), ("lag", "Wettest hour\nto peak (h)", "center"),
             ("ro_n", "Direct runoff\n(mm)", "right"), ("ratio", "Runoff\nratio", "center")]
    for k, label, ha in heads:
        ax.text(X[k], y_head, label, fontsize=12, fontweight="bold", color=c.INK2, ha=ha, va="center", linespacing=1.25)
    ax.plot([0.05, 0.97], [y_head - 0.04, y_head - 0.04], color=c.INK2, lw=1.2)

    maxr = ev.rain_mm.max() * 1.05; maxq = ev.peak_cms.max() * 1.05; maxro = ev.runoff_mm.max() * 1.05
    BW = {"rain": 0.08, "q": 0.095, "ro": 0.07}; BH = 0.026

    def bar(x, y, val, vmax, wmax, color):
        ax.add_patch(plt.Rectangle((x, y - BH / 2), wmax * val / vmax, BH, color=color, lw=0))

    for i, row in ev.iterrows():
        y = y0 - i * dy
        ax.text(X["storm"], y, row.storm, fontsize=15, fontweight="bold", color=c.INK2, ha="left", va="center")
        ax.text(X["storm"] + 0.022, y, row["name"].replace("-", "–"), fontsize=15, fontweight="bold", ha="left", va="center")
        ax.text(X["rain_n"], y, f"{row.rain_mm:.0f}", fontsize=15, ha="right", va="center")
        bar(X["rain_b"], y, row.rain_mm, maxr, BW["rain"], ORANGE)
        ax.text(X["wet"], y + 0.012, f"{row.wettest_hr_mm:.1f}", fontsize=15, ha="center", va="center")
        ax.text(X["wet"], y - 0.022, f"{row.wettest_hr_start:%b %-d, %H:%M}", fontsize=9.5, color=c.MUTED, ha="center", va="center")
        if row.peak_recorded:
            ax.text(X["q_n"], y + 0.012, fmt_q(row.peak_cms), fontsize=15, ha="right", va="center")
            ax.text(X["q_n"], y - 0.022, f"{row.peak_cfs:,.0f} cfs", fontsize=9.5, color=c.MUTED, ha="right", va="center")
            bar(X["q_b"], y, row.peak_cms, maxq, BW["q"], c.BLUE)
            lag = row.lag_h
            ax.text(X["lag"], y + (0.012 if row.storm == "4" else 0), f"{lag:.1f}" if lag < 10 else f"{lag:.0f}", fontsize=15, ha="center", va="center")
            if row.storm == "4":
                ax.text(X["lag"], y - 0.022, "slow response", fontsize=9.5, color=c.MUTED, ha="center", va="center")
            ax.text(X["ro_n"], y, f"{row.runoff_mm:.1f}", fontsize=15, ha="right", va="center")
            bar(X["ro_b"], y, row.runoff_mm, maxro, BW["ro"], c.BLUE)
            ax.text(X["ratio"], y, f"{row.runoff_ratio*100:.0f}%", fontsize=16, fontweight="bold", ha="center", va="center")
        else:
            ax.text(X["q_n"], y, "–", fontsize=15, color=c.MUTED, ha="right", va="center")
            ax.text(X["q_b"], y, "USGS outage:\npeak not recorded", fontsize=10.5, color=c.INK2, style="italic", ha="left", va="center", linespacing=1.25)
            for k in ("lag", "ro_n", "ratio"):
                ax.text(X[k], y, "–", fontsize=15, color=c.MUTED, ha="center" if k != "ro_n" else "right", va="center")
        ax.plot([0.05, 0.97], [y - dy / 2, y - dy / 2], color=c.GRID, lw=1.0)

    # ---- whole-year row ----
    y = y0 - 5 * dy
    ax.add_patch(plt.Rectangle((0.05, y - dy / 2 + 0.002), 0.92, dy - 0.004, color="#f1f0ec", lw=0, zorder=0))
    ax.plot([0.05, 0.97], [y + dy / 2, y + dy / 2], color=c.INK2, lw=1.2)
    ax.text(X["storm"], y, "Whole year", fontsize=15, fontweight="bold", ha="left", va="center")
    ax.text(X["rain_n"], y, f"{yr['rain_mm']:.0f}", fontsize=15, fontweight="bold", ha="right", va="center")
    ax.text(X["wet"], y + 0.012, f"{yr['wettest_hr_mm']:.1f}", fontsize=15, fontweight="bold", ha="center", va="center")
    ax.text(X["wet"], y - 0.022, f"{yr['wettest_hr_start']:%b %-d, %H:%M}", fontsize=9.5, color=c.MUTED, ha="center", va="center")
    ax.text(X["q_n"], y + 0.012, fmt_q(yr["peak_cms"]), fontsize=15, fontweight="bold", ha="right", va="center")
    ax.text(X["q_n"], y - 0.022, f"{yr['peak_cfs']:,.0f} cfs", fontsize=9.5, color=c.MUTED, ha="right", va="center")
    ax.text(X["lag"], y + 0.012, "–", fontsize=15, color=c.MUTED, ha="center", va="center")
    ax.text(X["ro_n"], y + 0.012, f"{yr['flow_mm']:.1f}", fontsize=15, fontweight="bold", ha="right", va="center")
    ax.text(X["ro_n"], y - 0.022, f"all flow, {yr['flow_acft']:,.0f} acre-ft", fontsize=9.5, color=c.MUTED, ha="right", va="center")
    ax.text(X["ratio"], y, f"{yr['ratio']*100:.0f}%", fontsize=16, fontweight="bold", ha="center", va="center")
    ax.plot([0.05, 0.97], [y - dy / 2, y - dy / 2], color=c.INK2, lw=1.2)

    # ---- notes ----
    notes = [
        "Rain: MRMS gauge-corrected hourly radar, area-weighted mean over the basin outline, summed over each storm's 72-hour window in Figure 1. The radar reads low at the gauges "
        "(69-87% of the mean of three FCDMC gauges in four storms), so rain is probably underestimated and runoff ratios overestimated. Radar hour labels are assumed to mark the end "
        "of the hour (a one-hour error shifts the lags by an hour).",
        f"Direct runoff: flow above a straight baseline drawn from the start of the rise to {E.N_DAYS:.1f} days after the peak (Linsley's rule, 0.827 A^0.2 days, A = 164 mi²), "
        f"as depth over the {E.AREA_KM2:.1f} km² basin. Cutting it off 1.5 to 5 days after the peak gives Jan 21 {rng['2'][0]:.0f}-{rng['2'][1]:.0f}%, Feb 28 {rng['5'][0]:.0f}-{rng['5'][1]:.0f}%, "
        f"and Jan 21's share of the year's flow {min(shares):.0f}-{max(shares):.0f}% ({share2:.0f}% as shown).",
        "Jan 1-2: a USGS gauge outage hid the rise and the peak (USGS-estimated daily means: 168 cfs on Jan 1, 232 cfs on Jan 2), so none is shown. "
        "Whole year: USGS daily means. Times are MST.",
    ]
    ax.text(0.05, 0.2, "\n".join(textwrap.fill(n, 195) for n in notes), fontsize=9, color=c.INK2, va="top", linespacing=1.45)
    fig.savefig(args.out, dpi=200)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
