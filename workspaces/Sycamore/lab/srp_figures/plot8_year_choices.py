#!/usr/bin/env python
"""Slide: the water-year screen (WY2016 to WY2026) and why WY2017 was chosen.

Every number is typed in from the project doc "Choosing the Water Year: Data Pulled and How We Screened It"
(Step 3 table, Step 5, Step 6; doc dated 2026-10-09). Nothing is read from model files, so this script has no inputs.
Peaks are counted on USGS daily-mean flow at gauge 09510200, so a peak is not the same as a storm.
Writes out/fig8_year_choices.png (16:9 slide) and out/year_screen_wy2016_2026.csv.
"""
import argparse, textwrap
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import srp_common as c

PRE = ("#cde2fb", "#184f95")          # light fill, dark text: before the June 2020 Bush Fire
POST = ("#fbe0d3", "#8f3411")         # after it

# (water year, fire label, fire colours, approved share, cool-season peaks >=50 cfs, monsoon peaks >=50 cfs,
#  largest daily mean cfs, its date, peaks flagged e (of all peaks), what counts against it)
ROWS = [
    ("WY2016", "Pre-fire", PRE, "100%", 1, 0, 347, "Jan 7", 0, "Only one peak above 50 cfs"),
    ("WY2017", "Pre-fire", PRE, "100%", 5, 0, 1920, "Jan 21", 1, "Winter-only; the Jan 21 peak dominates; Jan 1-2 flow is estimated or missing"),
    ("WY2018", "Pre-fire", PRE, "100%", 0, 2, 161, "Aug 8", 2, "Dry year; both peaks are estimated"),
    ("WY2019", "Pre-fire", PRE, "100%", 4, 1, 1160, "Sep 23", 3, "Three of its five peaks are estimated by USGS"),
    ("WY2020", "Pre-fire to June", PRE, "100%", 7, 1, 754, "Mar 13", 1, "Pre-fire only Oct to May; wet basin (March median 87 cfs)"),
    ("WY2021", "Post-fire", POST, "100%", 0, 5, 1310, "Aug 18", 1, "No cool-season flow; first year after the fire"),
    ("WY2022", "Post-fire", POST, "100%", 4, 0, 667, "Jan 1", 0, "Peaks packed into one week (Dec 25 to Jan 1)"),
    ("WY2023", "Post-fire", POST, "100%", 6, 0, 643, "Mar 22", 1, "Post-fire; wet basin (March median 112 cfs)"),
    ("WY2024", "Post-fire", POST, "100%", 3, 0, 382, "Apr 1", 0, "Small peaks only (largest 382 cfs)"),
    ("WY2025", "Post-fire", POST, "85%", 0, 1, 155, "Sep 26", 0, "Provisional from Aug 6, 2025; one peak (155 cfs)"),
    ("WY2026", "Post-fire", POST, "0%", 2, 2, 155, "Oct 12, 2025", 1, "All provisional; 19 days missing (Jul 17 to Aug 4)"),
]
# the e-count denominators in the doc's Step 3 table (cool + monsoon peaks)
E_OF = {"WY2016": 1, "WY2017": 5, "WY2018": 2, "WY2019": 5, "WY2020": 8, "WY2021": 5, "WY2022": 4,
        "WY2023": 6, "WY2024": 3, "WY2025": 1, "WY2026": 4}


def wrap(s, n):
    return "\n".join(textwrap.wrap(s, n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(c.OUT / "fig8_year_choices.png"))
    args = ap.parse_args()
    c.style()

    for r in ROWS:                                           # the table must agree with itself
        assert r[4] + r[5] == E_OF[r[0]], r[0]
    pd.DataFrame([dict(water_year=r[0], fire_status=r[1], approved_days=r[3], cool_season_peaks_ge50=r[4], monsoon_peaks_ge50=r[5],
                       largest_daily_mean_cfs=r[6], largest_date=r[7], peaks_flagged_e=f"{r[8]} of {E_OF[r[0]]}", against=r[9]) for r in ROWS]
                 ).to_csv(c.OUT / "year_screen_wy2016_2026.csv", index=False)

    W, H = 13.333, 7.5
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    head = "WY2017 is the only full pre-fire year with five winter peaks, and only one of them is estimated"
    ax.text(0.05, 0.935, wrap(head, 58), fontsize=21, fontweight="bold", va="center", linespacing=1.25)
    ax.text(0.05, 0.848, "Eleven water years screened on USGS daily flow at gauge 09510200 and on fire history. Counts are daily-mean peaks, not storms.",
            fontsize=12.5, color=c.INK2, va="center")

    X = dict(yr=0.05, fire=0.122, appr=0.215, win=0.285, mon=0.375, big=0.445, bigl=0.514, est=0.612, why=0.69)
    y_head = 0.792
    for k, label in (("yr", "Water year"), ("fire", "Fire status"), ("appr", "Approved\ndays"), ("win", "Winter\npeaks ≥50 cfs"),
                     ("mon", "Monsoon\npeaks"), ("big", "Largest daily\nmean (cfs)"), ("est", "Estimated\npeaks"), ("why", "What counts against it")):
        ax.text(X[k], y_head, label, fontsize=10.5, fontweight="bold", color=c.INK2, ha="left", va="center", linespacing=1.15)
    ax.plot([0.05, 0.97], [y_head - 0.031, y_head - 0.031], color=c.INK2, lw=1.2)

    dy = 0.0452
    y0 = y_head - 0.031 - dy / 2
    vmax, barw = 1920.0, 0.062
    DX = 0.0108
    for i, (yr, flab, fcol, appr, nwin, nmon, big, bigd, ne, why) in enumerate(ROWS):
        y = y0 - i * dy
        chosen = yr == "WY2017"
        if chosen:
            ax.add_patch(Rectangle((0.05, y - dy / 2 + 0.001), 0.92, dy - 0.002, color="#e8f1fc", lw=0, zorder=0))
            ax.add_patch(Rectangle((0.05, y - dy / 2 + 0.001), 0.004, dy - 0.002, color=c.BLUE, lw=0, zorder=1))
            ax.text(X["yr"] + 0.008, y + 0.0085, yr, fontsize=11.5, fontweight="bold", ha="left", va="center")
            ax.text(X["yr"] + 0.008, y - 0.0105, "chosen", fontsize=8.5, fontweight="bold", color=c.BLUE, ha="left", va="center")
        else:
            ax.text(X["yr"] + 0.008, y, yr, fontsize=11.5, fontweight="bold", ha="left", va="center")
        ax.text(X["fire"], y, flab, fontsize=9.5, fontweight="bold", color=fcol[1], ha="left", va="center",
                bbox=dict(boxstyle="round,pad=0.3,rounding_size=0.5", fc=fcol[0], ec="none"))
        ax.text(X["appr"], y, appr, fontsize=10, ha="left", va="center")
        kw = dict(marker="o", ms=7.5, mec=c.BLUE, mew=1.5, ls="none", clip_on=False)
        for k in range(nwin):
            ax.plot([X["win"] + 0.004 + k * DX], [y], mfc=c.BLUE, **kw)
        for k in range(nmon):
            ax.plot([X["mon"] + 0.004 + k * DX], [y], mfc=c.SURFACE, **kw)
        if nwin == 0:
            ax.text(X["win"], y, "none", fontsize=9, color=c.MUTED, ha="left", va="center")
        if nmon == 0:
            ax.text(X["mon"], y, "none", fontsize=9, color=c.MUTED, ha="left", va="center")
        ax.add_patch(Rectangle((X["big"], y - 0.007), barw * big / vmax, 0.014, color=c.BLUE, lw=0))
        ax.text(X["bigl"], y, f"{big:,} ({bigd})", fontsize=9.3, ha="left", va="center")
        ax.text(X["est"], y, f"{ne} of {E_OF[yr]}", fontsize=10, ha="left", va="center")
        ax.text(X["why"], y, wrap(why, 58), fontsize=9.3, fontweight="bold" if chosen else "normal", ha="left", va="center", linespacing=1.2)
        ax.plot([0.05, 0.97], [y - dy / 2, y - dy / 2], color=c.GRID, lw=1.0)

    # ---------------- what the choice costs ----------------
    ax.text(0.05, 0.232, "What choosing WY2017 costs", fontsize=10.5, fontweight="bold", color=c.INK2, va="center")
    cost = ("Winter-only: nearly all runoff comes in January to March, so the model is tested on winter storms, not monsoon flash floods (a question for Josh). "
            "One dominant event: Jan 21 (1,920 cfs daily mean, 7,800 cfs at 01:35) can dominate a fit score, so storms should also be scored one at a time. "
            "Jan 1 to 00:30 on Jan 2 has no 15-minute readings; those 98 rows are left blank and not scored.")
    ax.text(0.05, 0.213, wrap(cost, 185), fontsize=9.3, color=c.INK2, va="top", linespacing=1.3)

    notes = [
        "Peak: a daily mean of at least 50 cfs that rises at least half its height above nearby days; one storm can count as several peaks or merge into one. Estimated: flagged e by USGS.",
        "Filled dot: cool-season peak. Hollow dot: monsoon peak (Jun 15 to Sep 30). Peaks of 10 to 50 cfs are not shown. WY2016 is the first year the Iowa State radar archive covers safely.",
        "Source: project doc \"Choosing the Water Year\" (Oct 9, 2026), Steps 2 to 6; USGS gauge 09510200, daily means.",
    ]
    ax.text(0.05, 0.105, "\n".join(notes), fontsize=9, color=c.INK2, va="top", linespacing=1.5)
    fig.savefig(args.out, dpi=200)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
