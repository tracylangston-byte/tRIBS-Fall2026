"""Check whether the soil grids cover the whole Lower Sycamore basin.

The bedrock grid reaches all 42,449 basin cells, but the soil grids hold only 41,617.
This lists the basin cells (100 m cell centers inside the outline) that each soil grid
does not reach, and says whether they lie off the edge of the grid or inside it but
empty. It reads files only and changes nothing.

Usage in the notebook (working folder lab/):

    import glob, geopandas as gpd
    from soil_grid_gap_check import check_soil_grid_gap
    ws = gpd.read_file(glob.glob('data/preprocessing/*boundary*.shp')[0])
    check_soil_grid_gap(ws)
"""
import glob
import os

import numpy as np
import matplotlib.pyplot as plt
import rasterio
from rasterio.features import geometry_mask
from rasterio.transform import rowcol, xy
from rasterio.warp import transform as warp_xy

plt.rcParams['font.family'] = 'Liberation Sans'

REF = 'bedrock_depth_m.asc'                      # the grid that reaches the whole basin
SOIL = {'Ks': 'Ks_0-5cm.asc', 'theta_s': 'theta_s_0-5cm.asc', 'f': 'f.asc'}
DIRS = ['../init_data/solus', '../init_data', 'init_data/solus', 'init_data']


def _find(name):
    for d in DIRS:
        hits = glob.glob(os.path.join(d, name))
        if hits:
            return hits[0]
    raise FileNotFoundError(f'{name} not found in any of {DIRS}')


def check_soil_grid_gap(ws, ref=REF, soil=SOIL, out='soil_grid_gap_check.png'):
    """ws: basin outline (GeoDataFrame). Prints a report and saves a map of the missing cells."""
    with rasterio.open(_find(ref)) as r:
        ref_crs = r.crs if r.crs else ws.crs
        inside = geometry_mask(ws.to_crs(ref_crs).geometry, out_shape=r.shape,
                               transform=r.transform, invert=True)
        rows, cols = np.nonzero(inside)
        xs, ys = xy(r.transform, rows, cols)
        xs, ys = np.array(xs), np.array(ys)
        print(f'{ref}: {len(xs):,} basin cells (centers inside the outline); '
              f'grid {r.width} x {r.height} cells, cell size {r.res[0]:g}, '
              f'x {r.bounds.left:.0f}-{r.bounds.right:.0f}, y {r.bounds.bottom:.0f}-{r.bounds.top:.0f}')

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ws.to_crs(ref_crs).boundary.plot(ax=ax, color='0.4', linewidth=0.6)
    miss_any = np.zeros(len(xs), dtype=bool)

    for label, fname in soil.items():
        with rasterio.open(_find(fname)) as s:
            px, py = warp_xy(ref_crs, s.crs, xs, ys) if s.crs else (xs, ys)
            rr, cc = rowcol(s.transform, px, py)
            rr, cc = np.array(rr), np.array(cc)
            off = (rr < 0) | (rr >= s.height) | (cc < 0) | (cc >= s.width)
            band = s.read(1).astype('float64')
            if s.nodata is not None:
                band[band == s.nodata] = np.nan
            val = np.full(len(xs), np.nan)
            val[~off] = band[rr[~off], cc[~off]]
            empty = ~off & ~np.isfinite(val)
            miss = off | empty
            print(f'{fname}: grid {s.width} x {s.height}, cell size {s.res[0]:g}, '
                  f'x {s.bounds.left:.0f}-{s.bounds.right:.0f}, y {s.bounds.bottom:.0f}-{s.bounds.top:.0f}; '
                  f'basin cells missing {miss.sum():,} (off the grid {off.sum():,}, '
                  f'inside but empty {empty.sum():,})')
            if miss.any():
                print(f'    missing cells span x {xs[miss].min():.0f}-{xs[miss].max():.0f}, '
                      f'y {ys[miss].min():.0f}-{ys[miss].max():.0f}')
            miss_any |= miss

    ax.scatter(xs[~miss_any], ys[~miss_any], s=0.3, color='0.85', marker='s', linewidths=0)
    ax.scatter(xs[miss_any], ys[miss_any], s=1.5, color='black', marker='s', linewidths=0)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f'Basin cells the soil grids do not reach (black): {miss_any.sum():,}',
                 fontsize=9, loc='left')
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print('Saved', out)
    return xs[miss_any], ys[miss_any]
