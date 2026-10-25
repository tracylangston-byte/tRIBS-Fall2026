#!/usr/bin/env python
"""Figure: SOLUS-derived soil inputs for Lower Sycamore Wash - surface Ks, decay parameter f, depth to bedrock.

Inputs (read-only, from the repo): init_data/solus/Ks_0-5cm.asc, init_data/solus/f.asc, init_data/bedrock_depth_m.asc
Optional: --boundary path/to/Sycamore_boundary.shp  (NAD83 / UTM 12N) -> blanks cells outside the basin.
Each grid is drawn with its OWN header (the bedrock grid is 444 x 569 cells of ~100 m; Ks and f are
443 x 568 cells of 101 m), so nothing here assumes the grids line up cell for cell.
The map frame is the land cover grid extent (basin plus about 5 km buffer).
"""
import argparse
import textwrap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
import srp_common as c

FRAME = (441494.35, 470774.35, 3720135.19, 3768105.19)     # LandUse.asc extent (L, R, B, T)
CMAP = LinearSegmentedColormap.from_list("blues", c.BLUE_RAMP)


def window(a, ext, frame):
    """Cells of grid `a` (extent ext) whose centres fall inside `frame`."""
    ny, nx = a.shape
    xs = ext[0] + (np.arange(nx) + 0.5) * (ext[1] - ext[0]) / nx
    ys = ext[3] - (np.arange(ny) + 0.5) * (ext[3] - ext[2]) / ny
    X, Y = np.meshgrid(xs, ys)
    return (X >= frame[0]) & (X <= frame[1]) & (Y >= frame[2]) & (Y <= frame[3]), X, Y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boundary", default=str(c.BOUNDARY) if c.BOUNDARY.exists() else None)
    ap.add_argument("--out", default=str(c.OUT / "fig3_soils.png"))
    args = ap.parse_args()
    c.style()
    rings = c.load_boundary(args.boundary) if args.boundary else None

    panels = [
        ("Surface Ks (0-5 cm)", "solus/Ks_0-5cm.asc", "mm/hr", "linear"),
        ("Decay parameter f", "solus/f.asc", "1/mm (log scale)", "log"),
        ("Depth to bedrock", "bedrock_depth_m.asc", "m", "linear"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 8.6))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.85, bottom=0.33, wspace=0.18)
    stats = {}
    ims = []
    for ax, (title, rel, unit, scale) in zip(axes, panels):
        hdr, a, ext = c.read_asc(c.INIT / rel)
        inwin, X, Y = window(a, ext, FRAME)
        a = a.copy()
        if rings:
            a[~c.grid_mask(rings, ext, a.shape, ss=8)] = np.nan
        vals = a[inwin & np.isfinite(a)]
        lo, hi = np.percentile(vals, [2, 98])
        if scale == "log":
            norm = LogNorm(vmin=max(lo, 1e-4), vmax=hi)
        else:
            norm = Normalize(vmin=0 if rel.startswith("bedrock") else lo, vmax=hi)
        im = ax.imshow(a, cmap=CMAP, norm=norm, extent=ext, interpolation="nearest", origin="upper")
        if rings:
            for r in rings:
                ax.plot(r[:, 0], r[:, 1], color=c.INK, lw=1.2)
        ax.plot(*c.GAUGE_XY, marker="^", ms=8, mfc=c.INK, mec=c.SURFACE, mew=1.5, ls="none")
        if rings:
            allp = np.vstack(rings)
            ax.set_xlim(allp[:, 0].min() - 1000, allp[:, 0].max() + 1000)
            ax.set_ylim(allp[:, 1].min() - 1000, allp[:, 1].max() + 1000)
        else:
            ax.set_xlim(FRAME[0], FRAME[1]); ax.set_ylim(FRAME[2], FRAME[3])
        ax.set_aspect("equal")
        ax.set_title(title, loc="left", pad=8)
        ax.xaxis.set_major_formatter(lambda v, p: f"{v/1000:.0f}")
        ax.yaxis.set_major_formatter(lambda v, p: f"{v/1000:.0f}")
        ax.set_xticks([450000, 460000, 470000])
        ax.set_xlabel("Easting (km)")
        if ax is axes[0]:
            ax.set_ylabel("Northing (km, UTM zone 12N)")
        for s in ax.spines.values():
            s.set_visible(True); s.set_color(c.GRID)
        ims.append((ax, im, unit))
        at_min = float(np.mean(vals == vals.min())) * 100
        stats[title] = (float(vals.min()), float(np.median(vals)), float(vals.max()), lo, hi,
                        int(np.isnan(a[inwin]).sum()), at_min)

    fig.canvas.draw()                     # settles the equal-aspect map boxes
    for ax, im, unit in ims:
        pos = ax.get_position()
        cax = fig.add_axes([pos.x0, 0.255, pos.width, 0.014])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal")
        cb.outline.set_visible(False); cb.ax.tick_params(length=0)
        cb.set_label(unit, color=c.INK2, labelpad=4)
        if im.norm.__class__.__name__ == "LogNorm":
            cb.set_ticks([1e-4, 1e-3]); cb.ax.minorticks_off()

    fig.text(0.07, 0.945, "Soil inputs from SOLUS100, Lower Sycamore Wash", fontsize=15, fontweight="bold", va="center")
    fig.text(0.07, 0.905, "Model inputs for tRIBS: infiltration capacity, how fast it falls with depth, and how deep "
             "soil sits above bedrock", fontsize=10.5, color=c.INK2, va="center")
    f_floor = stats["Decay parameter f"][6]
    paras = ["Ks and f are not measurements. Ks comes from the Rosetta v3 pedotransfer function applied to SOLUS100 "
             "sand, silt, clay and bulk density; f is fitted to Ks over 0-30 cm.",
             f"Color scales span the 2nd to 98th percentile of the cells shown. {f_floor:.0f}% of the f cells shown sit "
             "at 0.0001, the smallest value in the grid (palest color). The bedrock scale starts at 0 m (rock at the "
             "surface). " + ("Cells outside the outline are blank. " if rings else "Blank bedrock cells have no SOLUS value. ") + "Triangle = USGS 09510200 (basin outlet)."]
    if not rings:
        paras.append("DRAFT: basin outline not yet applied. Maps show the full model frame (basin plus about 5 km buffer), "
                     "so the color ranges include cells outside the basin.")
    foot = "\n".join(textwrap.fill(p, 145) for p in paras)
    fig.text(0.07, 0.02, foot, fontsize=8.5, color=c.INK2, va="bottom", linespacing=1.5)
    fig.savefig(args.out)
    print("wrote", args.out)
    for k, v in stats.items():
        print(f"{k}: min {v[0]:.4g}  median {v[1]:.4g}  max {v[2]:.4g}  (2-98% {v[3]:.3g}-{v[4]:.3g}; NaN in frame {v[5]}; at min {v[6]:.1f}%)")


if __name__ == "__main__":
    main()
