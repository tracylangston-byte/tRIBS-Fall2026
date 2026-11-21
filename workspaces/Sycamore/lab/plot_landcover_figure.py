"""
Land cover map of the Lower Sycamore basin (thesis figure for Section 4.3.3).

One map of the eight land cover classes in LandUse.asc, clipped to the basin outline,
with a legend that lists each class and its share of the basin, plus the boundary and
gauge symbols, a scale bar and a north arrow. The shares are counted from the grid, not
typed in: compare the printed table with Appendix E (about 46.7, 35.4, 10.3, 2.9, 2.6,
0.8, 0.7 and 0.7 percent).

Usage in the notebook (working folder lab/):

    import glob, geopandas as gpd
    from plot_landcover_figure import plot_landcover
    ws = gpd.read_file(glob.glob('data/preprocessing/*boundary*.shp')[0])
    fig = plot_landcover(ws, out='Sycamore_landcover_figure.png')

The script looks for LandUse.asc in ../init_data and init_data. If yours is elsewhere,
pass the path: plot_landcover(ws, landuse_file='path/to/LandUse.asc').

Colors: every figure color in this thesis other than the slide-deck schematics has its own
palette. These eight were chosen to stay distinguishable for each pair of classes in
normal vision (smallest difference 15.2 on the OKLab x 100 scale) and under simulated
red-blind and green-blind vision (11.4 and 11.1), checked with a palette validator.
Bare rock is gray and woodland is dark green by design, so the legend is needed to read
the map; the classes cannot be told apart by color against the white page alone.
The outlet gauge is a red circle with a black edge and a white ring (red is not used by
any class), the same red as in the mesh figure. GAUGE_MS sets its size.

The grid is drawn cell by cell without smoothing. At 6.5 in wide, one printed pixel
covers about 35 m, so narrow washes (riparian class) can look broken even though every
10 m cell is counted in the shares.
"""
import glob
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.plot import plotting_extent
from shapely.geometry import Point

GRID_DIRS = ['../init_data', 'init_data', '../init_data/landcover', 'init_data/landcover']
LANDUSE_NAMES = ['LandUse.asc', 'landuse.asc', 'LandUse*.asc']
GAUGE = (449721, 3728390)          # USGS 09510200, UTM 12N (EPSG:26912)
# GAUGE_COLOR = '#d62728'            # same red as the mesh figure
GAUGE_COLOR = 'white'            
GAUGE_MS = 12                       # marker size in points; raise it to make the gauge larger

# (class ID in LandUse.asc, legend name, color)
CLASSES = [
    (1, 'Woodland',                  '#1e5631'),
    (2, 'Desert scrub',              '#b89560'),
    (3, 'Chaparral and oak shrub',   '#6e7b2a'),
    (4, 'Bare rock',                 '#bfc5cf'),
    (5, 'Recently burned or logged', '#4a2a6a'),
    (6, 'Developed',                 '#cc3d8e'),
    (7, 'Riparian and water',        '#00a6d6'),
    (8, 'Grassland',                 '#c8e06a'),
]

FONT_SETTINGS = {
    'font.family': 'Liberation Sans',
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Liberation Sans',
    'mathtext.it': 'Liberation Sans:italic',
    'mathtext.bf': 'Liberation Sans:bold',
    'mathtext.sf': 'Liberation Sans',
    'mathtext.tt': 'Liberation Sans',
    'mathtext.cal': 'Liberation Sans',
}


def _find_landuse():
    for d in GRID_DIRS:
        for name in LANDUSE_NAMES:
            hits = sorted(glob.glob(os.path.join(d, name)))
            if hits:
                return hits[0]
    hits = sorted(glob.glob('../**/LandUse.asc', recursive=True))
    if hits:
        print(f'LandUse.asc found by wider search: {hits[0]}')
        return hits[0]
    raise FileNotFoundError(f'LandUse.asc not found in {GRID_DIRS}; pass landuse_file=...')


def _read_landuse(path, ws):
    """Class grid clipped to the basin: NaN outside the outline and at nodata cells."""
    with rasterio.open(path) as src:
        crs = src.crs if src.crs else ws.crs
        outline = ws.to_crs(crs).geometry
        data, tr = mask(src, outline, crop=True, nodata=0, filled=True)
        res = src.res
    arr = data[0].astype('float32')
    arr[(arr < 1) | (arr > len(CLASSES))] = np.nan
    return arr, plotting_extent(arr, tr), crs, res


def plot_landcover(ws, landuse_file=None, gauge=GAUGE, out='Sycamore_landcover_figure.png',
                   dpi=300, width_in=6.5, max_height_in=8.8):
    """ws: basin outline (GeoDataFrame)."""
    path = landuse_file or _find_landuse()
    arr, ext, crs, res = _read_landuse(path, ws)

    # class table: counted from the grid
    total = int(np.isfinite(arr).sum())
    cell_km2 = res[0] * res[1] / 1e6
    shares = {}
    print(f'{os.path.basename(path)}: cell size {res[0]:g} x {res[1]:g} m, '
          f'{total:,} cells in basin = {total * cell_km2:,.1f} km2')
    print(f'{"ID":<3} {"Class":<27} {"cells":>10} {"km2":>8} {"share":>7}')
    for cid, name, _ in CLASSES:
        n = int((arr == cid).sum())
        shares[cid] = 100.0 * n / total if total else 0.0
        print(f'{cid:<3} {name:<27} {n:>10,} {n * cell_km2:>8.1f} {shares[cid]:>6.1f}%')
        if n == 0:
            print(f'    note: no cells of class {cid} inside the outline')

    ws_g = ws.to_crs(crs)
    xmin, ymin, xmax, ymax = ws_g.total_bounds
    pad = 0.01 * max(xmax - xmin, ymax - ymin)
    xmin, xmax, ymin, ymax = xmin - pad, xmax + pad, ymin - pad, ymax + pad
    gx, gy = gpd.GeoSeries([Point(gauge)], crs='EPSG:26912').to_crs(crs).iloc[0].coords[0]

    # layout in inches: map on the left, legend and keys on the right
    mx, gap = 0.08, 0.18
    aspect = (ymax - ymin) / (xmax - xmin)
    map_w = 3.95
    map_w = min(map_w, (max_height_in - 2 * mx) / aspect)
    map_h = map_w * aspect
    leg_w = width_in - 2 * mx - gap - map_w
    body_h = max(map_h, 3.9)                       # the legend needs about this much room
    W, H = width_in, 2 * mx + body_h

    cmap = ListedColormap([c for _, _, c in CLASSES]).with_extremes(bad=(1, 1, 1, 0))
    norm = BoundaryNorm(np.arange(0.5, len(CLASSES) + 1.5, 1.0), cmap.N)

    with plt.rc_context(FONT_SETTINGS):
        fig = plt.figure(figsize=(W, H), facecolor='white')

        ax = fig.add_axes([mx / W, (H - mx - map_h) / H, map_w / W, map_h / H])
        ax.imshow(arr, extent=ext, origin='upper', cmap=cmap, norm=norm, interpolation='nearest')
        ws_g.boundary.plot(ax=ax, color='black', linewidth=0.6)
        # white ring under a red circle with a black edge, so it shows on dark woodland too
        ax.plot(gx, gy, marker='o', ms=GAUGE_MS + 4, mfc='white', mec='white', mew=0,
                ls='none', zorder=4)
        ax.plot(gx, gy, marker='o', ms=GAUGE_MS, mfc=GAUGE_COLOR, mec='black', mew=0.9,
                ls='none', zorder=5)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        # legend and keys: one axes in inches, origin at its lower left
        lx = mx + map_w + gap
        lax = fig.add_axes([lx / W, mx / H, leg_w / W, body_h / H])
        lax.set_xlim(0, leg_w)
        lax.set_ylim(0, body_h)
        lax.axis('off')
        fs = 7.5
        top = body_h - 0.04
        lax.text(0.0, top, 'Land cover class', fontsize=fs, fontweight='bold', va='center',
                 color='black')
        lax.text(leg_w - 0.02, top, 'Share of basin', fontsize=fs, fontweight='bold', va='center',
                 ha='right', color='black')
        row_h = 0.25
        y = top - 0.30
        for cid, name, color in CLASSES:
            lax.add_patch(Rectangle((0.0, y - 0.07), 0.22, 0.14, facecolor=color,
                                    edgecolor='#4d4d4d', linewidth=0.4))
            lax.text(0.32, y, name, fontsize=fs, va='center', color='black')
            lax.text(leg_w - 0.02, y, f'{shares[cid]:.1f}%', fontsize=fs, va='center',
                     ha='right', color='black')
            y -= row_h

        y1 = y - 0.12
        y2 = y1 - 0.26
        lax.plot([0.02, 0.20], [y1, y1], color='black', lw=0.8)
        lax.text(0.32, y1, 'Basin boundary', fontsize=fs, va='center', color='black')
        lax.plot([0.11], [y2], marker='o', ms=GAUGE_MS * 0.8, mfc=GAUGE_COLOR, mec='black',
                 mew=0.9, ls='none')
        lax.text(0.32, y2, 'Outlet gauge', fontsize=fs, va='center', color='black')

        m_per_in = (xmax - xmin) / map_w                  # map meters per inch
        km = next(k for k in (20, 10, 5, 2, 1) if k * 1000 / m_per_in <= 0.62 * leg_w)
        L = km * 1000 / m_per_in
        yb = 0.16
        lax.plot([0.02, 0.02 + L], [yb, yb], color='black', lw=1.0, solid_capstyle='butt')
        for xe in (0.02, 0.02 + L):
            lax.plot([xe, xe], [yb - 0.04, yb + 0.04], color='black', lw=1.0)
        lax.text(0.02 + L / 2, yb + 0.07, f'{km} km', fontsize=fs, ha='center', va='bottom',
                 color='black')
        xn = leg_w - 0.18
        lax.annotate('', xy=(xn, 0.50), xytext=(xn, 0.10),
                     arrowprops=dict(arrowstyle='-|>', color='black', lw=0.8, mutation_scale=7))
        lax.text(xn, 0.53, 'N', fontsize=fs, ha='center', va='bottom', color='black')

        fig.savefig(out, dpi=dpi, facecolor='white')
    print(f'Saved {out}  ({W:.2f} x {H:.2f} in at {dpi} dpi)')
    return fig
