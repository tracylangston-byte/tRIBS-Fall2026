#!/usr/bin/env python
"""Figure: land cover map of Lower Sycamore Wash (8 tRIBS classes) with basin areas.

Inputs (read-only, from the repo): init_data/LandUse.asc, init_data/LandUse_crosswalk.csv
Optional: --boundary path/to/Sycamore_boundary.shp  (NAD83 / UTM 12N) -> clips the map to the basin
Areas in the table come from the crosswalk's BASIN_CELLS column (10 m cells inside the basin outline
used when the crosswalk was written), so they do not depend on the map being clipped.
"""
import argparse
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap
import srp_common as c

CLASSES = {1: "Woodland", 2: "Desert scrub", 3: "Chaparral and oak shrub", 4: "Bare rock",
           5: "Recently burned or logged", 6: "Developed", 7: "Riparian and water", 8: "Grassland"}
COLORS = {1: "#2f6b3f", 2: "#d9b56c", 3: "#7d8f2a", 4: "#bdbab2",
          5: "#5b3a2e", 6: "#c2415d", 7: "#2a78d6", 8: "#dfe6a3"}
BLOCK = 3   # 10 m cells -> 30 m map cells (majority class of each 3 x 3 block)


def majority_downsample(cls, k):
    ny, nx = (cls.shape[0] // k) * k, (cls.shape[1] // k) * k
    cls = cls[:ny, :nx]
    best = np.zeros((ny // k, nx // k), dtype=np.int32)
    best_n = np.full(best.shape, -1, dtype=np.int32)
    for cid in CLASSES:
        n = (cls == cid).reshape(ny // k, k, nx // k, k).sum(axis=(1, 3))
        take = n > best_n
        best[take], best_n[take] = cid, n[take]
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boundary", default=str(c.BOUNDARY) if c.BOUNDARY.exists() else None)
    ap.add_argument("--out", default=str(c.OUT / "fig2_landcover.png"))
    args = ap.parse_args()
    c.style()

    hdr, a, (L, R, B, T) = c.read_asc(c.INIT / "LandUse.asc")
    cls = np.nan_to_num(a, nan=0).astype(np.int16)
    assert set(np.unique(cls)) <= set(CLASSES) | {0}, np.unique(cls)
    small = majority_downsample(cls, BLOCK)
    cs = hdr["cellsize"] * BLOCK
    ext = (L, L + small.shape[1] * cs, T - small.shape[0] * cs, T)

    rings = c.load_boundary(args.boundary) if args.boundary else None
    img = small.astype(float)
    if rings:
        img[~c.grid_mask(rings, ext, small.shape, ss=4)] = np.nan
    img[img == 0] = np.nan

    cw = pd.read_csv(c.INIT / "LandUse_crosswalk.csv")
    g = cw.groupby("TRIBS_ID").BASIN_CELLS.sum()
    km2 = g * (hdr["cellsize"] ** 2) / 1e6
    pct = 100 * g / g.sum()

    fig = plt.figure(figsize=(11, 7.6))
    ax = fig.add_axes([0.06, 0.09, 0.44, 0.78])
    cmap = ListedColormap([COLORS[i] for i in CLASSES])
    ax.imshow(img, cmap=cmap, vmin=0.5, vmax=8.5, extent=ext, interpolation="nearest", origin="upper")
    if rings:
        for r in rings:
            ax.plot(r[:, 0], r[:, 1], color=c.INK, lw=1.2)
    ax.plot(*c.GAUGE_XY, marker="^", ms=9, mfc=c.INK, mec=c.SURFACE, mew=1.5, ls="none")
    ax.annotate("USGS 09510200\n(basin outlet)", c.GAUGE_XY, xytext=(10, -6), textcoords="offset points",
                fontsize=8.5, color=c.INK, va="top",
                path_effects=[pe.withStroke(linewidth=3, foreground=c.SURFACE)])
    if rings:      # zoom to the basin outline (plus 1 km)
        allp = np.vstack(rings)
        ax.set_xlim(allp[:, 0].min() - 1000, allp[:, 0].max() + 1000)
        ax.set_ylim(allp[:, 1].min() - 1000, allp[:, 1].max() + 1000)
    else:
        ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(lambda v, p: f"{v/1000:.0f}")
    ax.yaxis.set_major_formatter(lambda v, p: f"{v/1000:.0f}")
    ax.set_xlabel("Easting (km, UTM zone 12N)"); ax.set_ylabel("Northing (km)")
    for s in ax.spines.values():
        s.set_visible(True); s.set_color(c.GRID)

    # legend table on the right
    tx = fig.add_axes([0.54, 0.09, 0.42, 0.78]); tx.axis("off")
    tx.set_xlim(0, 1); tx.set_ylim(0, 1)
    tx.text(0, 0.97, "Class", fontsize=9.5, color=c.INK2, fontweight="bold", va="top")
    tx.text(0.80, 0.97, "km²", fontsize=9.5, color=c.INK2, fontweight="bold", va="top", ha="right")
    tx.text(0.97, 0.97, "% of basin", fontsize=9.5, color=c.INK2, fontweight="bold", va="top", ha="right")
    y = 0.91
    for cid, name in CLASSES.items():
        tx.add_patch(plt.Rectangle((0, y - 0.034), 0.05, 0.034, color=COLORS[cid], lw=0, transform=tx.transAxes))
        tx.text(0.075, y - 0.017, name, fontsize=10.5, color=c.INK, va="center")
        tx.text(0.80, y - 0.017, f"{km2[cid]:.1f}", fontsize=10.5, color=c.INK, va="center", ha="right")
        tx.text(0.97, y - 0.017, f"{pct[cid]:.1f}", fontsize=10.5, color=c.INK, va="center", ha="right")
        y -= 0.058
    tx.plot([0, 0.97], [y + 0.012] * 2, color=c.GRID, lw=1)
    tx.text(0.075, y - 0.025, "Basin total", fontsize=10.5, color=c.INK, fontweight="bold", va="center")
    tx.text(0.80, y - 0.025, f"{km2.sum():.1f}", fontsize=10.5, color=c.INK, fontweight="bold", va="center", ha="right")
    paras = ["Source: LANDFIRE 2016 Remap Existing Vegetation Type, grouped into 8 tRIBS classes "
             "on the 10 m model grid.",
             "The map is drawn at 30 m (majority class per 3 x 3 cells), so narrow riparian and road "
             "cells are smoothed out of the picture but are counted in the areas."]
    if not rings:
        paras.append("DRAFT: basin outline not yet applied. The map shows the full model grid "
                     "(basin plus about 5 km buffer); the areas are for the basin only.")
    note = "\n\n".join(textwrap.fill(p, 66) for p in paras)
    tx.text(0, y - 0.11, note, fontsize=8.5, color=c.INK2, va="top", linespacing=1.4)

    fig.text(0.06, 0.955, "Land cover, Lower Sycamore Wash", fontsize=15, fontweight="bold", va="center")
    fig.text(0.06, 0.915, "Woodland and desert scrub cover about 82% of the basin", fontsize=10.5,
             color=c.INK2, va="center")
    fig.savefig(args.out)
    print("wrote", args.out)
    print((pd.DataFrame({"class": pd.Series(CLASSES), "km2": km2.round(2), "pct": pct.round(1)})).to_string())


if __name__ == "__main__":
    main()
