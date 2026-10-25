#!/usr/bin/env python
"""Slide: wildfire history and why the model uses a pre-fire year (WY2017).

Every number is typed in from the project doc "Choosing the Water Year: Data Pulled and How We Screened It"
(Step 5, Step 6 and the Sources list; doc dated 2026-10-09). Nothing is read from model files, so this script has no
inputs. The burned areas come from a Forest Service report and have NOT been clipped to the basin boundary.
Writes out/fig7_fire_history.png (16:9 slide).
"""
import argparse, textwrap
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import srp_common as c

PRE = ("#cde2fb", "#184f95")          # light fill, dark text: before the June 2020 Bush Fire
POST = ("#fbe0d3", "#8f3411")         # after it
SEV = [("Unburned", 9.4, "#e6e5e1"), ("Low", 41.3, "#f6cbb0"), ("Moderate", 44.7, "#eb6834"), ("High", 4.6, "#8f3411")]
HUC = [("Upper Sycamore", 15007, 39506), ("Middle Sycamore", 11388, 32885), ("Lower Sycamore", 4997, 28327)]

# (water year, status text, status colours, what the screen found)
YEARS = [
    ("WY2016", "Pre-fire", PRE, "Only one peak above 50 cfs"),
    ("WY2017", "Pre-fire", PRE, "Chosen. Every day approved; five winter peaks, 190 to 1,920 cfs (daily mean)"),
    ("WY2018", "Pre-fire", PRE, "Dry year"),
    ("WY2019", "Pre-fire", PRE, "Three of its five peaks are estimated by USGS"),
    ("WY2020", "Pre-fire to June", PRE, "Seven winter peaks, but only Oct to May is pre-fire"),
    ("WY2021", "Post-fire", POST, "First year after the fire; no cool-season flow"),
    ("WY2022", "Post-fire", POST, "Four peaks above 50 cfs, packed into Dec 25 to Jan 1"),
    ("WY2023", "Post-fire", POST, "Six winter peaks; the Bullet Fire came after the winter"),
    ("WY2024", "Post-fire", POST, "Small peaks only (largest 382 cfs)"),
    ("WY2025", "Post-fire", POST, "Provisional or incomplete; fires not checked"),
    ("WY2026", "Post-fire", POST, "Provisional or incomplete record"),
]


def wrap(s, n):
    return "\n".join(textwrap.wrap(s, n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(c.OUT / "fig7_fire_history.png"))
    args = ap.parse_args()
    c.style()

    burned = sum(b for _, b, _ in HUC)
    basin_acres = 424.0e6 / 4046.856                      # the doc compares against the 424 km2 basin
    share = burned / basin_acres * 100
    assert abs(sum(p for _, p, _ in SEV) - 100) < 0.05
    assert 29 < share < 31, share

    W, H = 13.333, 7.5
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    head = "The 2020 Bush Fire burned part of this basin, so the model uses a pre-fire year: WY2017"
    ax.text(0.05, 0.935, wrap(head, 58), fontsize=21, fontweight="bold", va="center", linespacing=1.25)
    ax.text(0.05, 0.848, wrap("A pre-fire year means one set of land-cover and soil conditions describes the whole run. "
                              "Years after June 2020 describe a partly burned basin.", 150), fontsize=12.5, color=c.INK2, va="center", linespacing=1.4)

    # ---------------- left: the year-by-year screen ----------------
    XL = dict(yr=0.05, st=0.125, note=0.245)
    y_head = 0.792
    for k, label in (("yr", "Water year"), ("st", "Fire status"), ("note", "What the flow screen found")):
        ax.text(XL[k], y_head, label, fontsize=11.5, fontweight="bold", color=c.INK2, ha="left", va="center")
    ax.plot([0.05, 0.60], [y_head - 0.025, y_head - 0.025], color=c.INK2, lw=1.2)
    dy = 0.046
    y0 = y_head - 0.025 - dy / 2
    for i, (yr, st, col, note) in enumerate(YEARS):
        y = y0 - i * dy
        chosen = yr == "WY2017"
        if chosen:
            ax.add_patch(Rectangle((0.05, y - dy / 2 + 0.001), 0.55, dy - 0.002, color="#e8f1fc", lw=0, zorder=0))
            ax.add_patch(Rectangle((0.05, y - dy / 2 + 0.001), 0.004, dy - 0.002, color=c.BLUE, lw=0, zorder=1))
        ax.text(XL["yr"] + 0.008, y, yr, fontsize=11.5, fontweight="bold", ha="left", va="center")
        ax.text(XL["st"], y, st, fontsize=9.5, fontweight="bold", color=col[1], ha="left", va="center",
                bbox=dict(boxstyle="round,pad=0.3,rounding_size=0.5", fc=col[0], ec="none"))
        ax.text(XL["note"], y, wrap(note, 58), fontsize=9.8, fontweight="bold" if chosen else "normal", ha="left", va="center", linespacing=1.2)
        ax.plot([0.05, 0.60], [y - dy / 2, y - dy / 2], color=c.GRID, lw=1.0)

    # ---------------- right: what the Bush Fire burned ----------------
    X0, XW = 0.645, 0.325
    y = 0.792
    ax.text(X0, y, "The Bush Fire, June 13 to about July 6, 2020", fontsize=12.5, fontweight="bold", va="center")
    y -= 0.034
    ax.text(X0, y, wrap("195,107 acres in a July 21, 2020 burned-area memo; 193,455 acres on Wikipedia. Tonto National Forest.", 62),
            fontsize=9.5, color=c.INK2, va="top", linespacing=1.3)
    # burn severity bar
    y -= 0.075
    ax.text(X0, y, "Burn severity, whole fire", fontsize=10.5, fontweight="bold", color=c.INK2, va="center")
    y -= 0.032
    x = X0
    for name, pct, colr in SEV:
        ax.add_patch(Rectangle((x, y - 0.012), XW * pct / 100, 0.024, color=colr, lw=0))
        x += XW * pct / 100
    y -= 0.032
    for (name, pct, colr), dx in zip(SEV, [0.0, 0.093, 0.1635, 0.259]):
        ax.add_patch(Rectangle((X0 + dx, y - 0.007), 0.009, 0.014, color=colr, lw=0))
        ax.text(X0 + dx + 0.012, y, f"{name} {pct:.1f}%", fontsize=9, color=c.INK2, va="center", ha="left")
    # burned area by HUC-12
    y -= 0.058
    ax.text(X0, y, "Burned acres, Sycamore Creek sub-watersheds", fontsize=10.5, fontweight="bold", color=c.INK2, va="center")
    vmax = max(t for _, _, t in HUC)
    for name, b, t in HUC:
        y -= 0.056
        ax.text(X0, y, f"{name}: {b:,} of {t:,} ({b / t * 100:.0f}%)", fontsize=9.5, va="center")
        ax.add_patch(Rectangle((X0, y - 0.027), XW * t / vmax, 0.014, color="#e6e5e1", lw=0))
        ax.add_patch(Rectangle((X0, y - 0.027), XW * b / vmax, 0.014, color="#eb6834", lw=0))
    y -= 0.065
    ax.text(X0, y, wrap(f"Together about {burned:,} acres, roughly {share:.0f}% of the 424 km² basin at most. Treat that as an upper bound: part of Lower "
                        f"Sycamore lies below the gauge, and the report's own percent column does not match its acres (Upper: 38% from the acres, 41.5% in its table).", 66),
            fontsize=9.3, color=c.INK2, va="top", linespacing=1.3)
    # ---------------- left, under the table: fires after WY2017 ----------------
    ax.text(0.05, 0.232, "Fires after WY2017, sizes as reported (sources disagree)", fontsize=10.5, fontweight="bold", color=c.INK2, va="center")
    ax.text(0.05, 0.213, wrap("Bullet, Jun 2023, about 3,355 acres, Lower Sycamore area. Diamond, Jul 2023, about 5,060, and Sand Stone, Jul 2024, 13,000 to 27,000, "
                              "both near Sunflower; the articles do not show either inside the basin, which is not the same as confirmed outside. "
                              "Sycamore, Jun 2026, about 12,000, in a different watershed.", 112),
            fontsize=9.3, color=c.INK2, va="top", linespacing=1.3)

    # ---------------- notes ----------------
    notes = [
        "Not measured: how much of the burn falls inside the model boundary. Next step: intersect the National Interagency Fire Center perimeter history with the basin polygon in ArcGIS Pro.",
        "Not checked: the 2012 Sunflower Fire (about 17,500 acres in the upper Sycamore area, per an older note) and any fires in WY2025.",
        "Sources: USFS Bush Fire Watershed Specialist Report (Jul 2020) and BAER funding request; Wikipedia; Fox10, KJZZ and ABC15. From the project doc \"Choosing the Water Year\" (Oct 9, 2026).",
    ]
    ax.text(0.05, 0.105, "\n".join(notes), fontsize=9, color=c.INK2, va="top", linespacing=1.5)
    fig.savefig(args.out, dpi=200)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
