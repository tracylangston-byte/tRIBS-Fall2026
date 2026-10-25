#!/usr/bin/env python
"""Figure: total rainfall over the basin, WY2017 (MRMS gauge-corrected radar, 1 km cells) plus monthly basin-mean depth.

Inputs: out/radar_wy2017_cache.npz (from load_radar.py), out/basin_mean_hourly_rain_wy2017.csv, the basin outline shapefile.
"""
import argparse, textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
import srp_common as c

EXT = (440000.0, 472000.0, 3720000.0, 3769000.0)
CMAP = LinearSegmentedColormap.from_list("blues", c.BLUE_RAMP)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boundary", default=str(c.BOUNDARY))
    ap.add_argument("--out", default=str(c.OUT / "fig4_rainfall.png"))
    args = ap.parse_args()
    c.style()

    z = np.load(c.OUT / "radar_wy2017_cache.npz")
    total = z["rain"].sum(axis=0)                  # mm over the year, per 1 km cell
    cover = z["cover"]
    w = cover / cover.sum()
    basin_total = float((total * w).sum())
    inside = cover >= 0.5
    lo, hi = float(total[inside].min()), float(total[inside].max())
    rings = c.load_boundary(args.boundary)

    h = pd.read_csv(c.OUT / "basin_mean_hourly_rain_wy2017.csv", parse_dates=["interval_start_mst"]).set_index("interval_start_mst")
    monthly = h.basin_mean_mm.resample("MS").sum()
    assert abs(monthly.sum() - basin_total) < 0.5, (monthly.sum(), basin_total)
    share = 100 * (monthly["2017-07-01"] + monthly["2017-01-01"]) / monthly.sum()

    fig = plt.figure(figsize=(12, 8.8))
    ax = fig.add_axes([0.06, 0.27, 0.38, 0.58])
    im = ax.imshow(total, cmap=CMAP, vmin=lo, vmax=hi, extent=EXT, interpolation="nearest", origin="upper")
    compound = MPath.make_compound_path(*[MPath(r) for r in rings])
    clip = PathPatch(compound, transform=ax.transData, facecolor="none", edgecolor="none")
    ax.add_patch(clip); im.set_clip_path(clip)
    for r in rings:
        ax.plot(r[:, 0], r[:, 1], color=c.INK, lw=1.2)
    ax.plot(*c.GAUGE_XY, marker="^", ms=8, mfc=c.INK, mec=c.SURFACE, mew=1.5, ls="none")
    allp = np.vstack(rings)
    ax.set_xlim(allp[:, 0].min() - 1000, allp[:, 0].max() + 1000); ax.set_ylim(allp[:, 1].min() - 1000, allp[:, 1].max() + 1000)
    ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(lambda v, p: f"{v/1000:.0f}"); ax.yaxis.set_major_formatter(lambda v, p: f"{v/1000:.0f}")
    ax.set_xticks([450000, 455000, 460000, 465000])
    ax.set_xlabel("Easting (km, UTM zone 12N)"); ax.set_ylabel("Northing (km)")
    for s in ax.spines.values():
        s.set_visible(True); s.set_color(c.GRID)
    ax.set_title("Total rain, Oct 2016 to Sep 2017", loc="left", pad=8)
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0, 0.185, pos.width, 0.014])
    fig.canvas.draw(); pos = ax.get_position(); cax.set_position([pos.x0, 0.185, pos.width, 0.014])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False); cb.ax.tick_params(length=0); cb.set_label("mm per 1 km radar cell", color=c.INK2, labelpad=4)

    ax2 = fig.add_axes([0.58, 0.27, 0.39, 0.58])
    labs = [m.strftime("%b") for m in monthly.index]
    ax2.bar(range(12), monthly.values, width=0.62, color=c.BLUE, lw=0)
    ax2.set_xticks(range(12)); ax2.set_xticklabels(labs, fontsize=9)
    ax2.set_ylim(0, monthly.max() * 1.15)
    ax2.set_ylabel("Basin-mean rain (mm)"); ax2.grid(axis="y"); ax2.set_axisbelow(True)
    ax2.spines["left"].set_visible(False); ax2.tick_params(axis="y", length=0); ax2.tick_params(axis="x", length=0)
    for i in monthly.reset_index(drop=True).nlargest(2).index:
        ax2.annotate(f"{monthly.iloc[i]:.0f}", (i, monthly.iloc[i]), xytext=(0, 4), textcoords="offset points",
                     ha="center", fontsize=9.5, color=c.INK)
    ax2.set_title("Basin-mean rain by month", loc="left", pad=8)

    fig.text(0.06, 0.955, "Rainfall over Lower Sycamore Wash, water year 2017", fontsize=15, fontweight="bold", va="center")
    fig.text(0.06, 0.92, f"Basin-mean total {basin_total:.0f} mm ({basin_total/25.4:.1f} in); July and January together brought {share:.0f}% of it",
             fontsize=10.5, color=c.INK2, va="center")
    paras = ["MRMS gauge-corrected hourly radar rain on 1 km cells, summed over 8,760 hours (Arizona time). Basin mean is area-weighted over the basin "
             f"outline; cell totals inside the basin run from {lo:.0f} to {hi:.0f} mm. The 17 hours missing from the radar archive (0.19%) count as zero rain.",
             "Radar has a low bias at the gauges: in a hand check of four widespread storms, 72-hour radar totals were 69-87% of the mean of three FCDMC gauges. "
             "Triangle = USGS 09510200 (basin outlet)."]
    fig.text(0.06, 0.012, "\n".join(textwrap.fill(p, 150) for p in paras), fontsize=8.5, color=c.INK2, va="bottom", linespacing=1.5)
    fig.savefig(args.out)
    print("wrote", args.out)
    print(f"basin-mean total {basin_total:.1f} mm; cell range inside basin {lo:.0f}-{hi:.0f} mm; Jul+Jan share {share:.1f}%")
    print(monthly.round(1).to_string())


if __name__ == "__main__":
    main()
