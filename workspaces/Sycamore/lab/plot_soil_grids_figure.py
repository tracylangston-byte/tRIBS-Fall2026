"""
Soil parameter and depth-to-bedrock maps for the Lower Sycamore basin (thesis figure).

Seven maps in table order (Ks, theta_s, theta_r, m, psi_b, f, depth to bedrock), each
clipped to the basin outline with its own horizontal color bar, plus a key panel
(outline, gauge, north arrow, scale bar).

Usage in the notebook (working folder lab/):

    import glob, geopandas as gpd
    ws = gpd.read_file(glob.glob('data/preprocessing/*boundary*.shp')[0])
    fig = plot_soil_grids(ws, out='Sycamore_soil_grids_figure.png')
"""
import glob
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import FuncFormatter, MaxNLocator
import geopandas as gpd
import rasterio
from rasterio.features import geometry_mask
from rasterio.mask import mask
from rasterio.plot import plotting_extent
from shapely.geometry import Point

# Folders searched for the grids (relative to lab/)
GRID_DIRS = ['../init_data/solus', '../init_data', 'init_data/solus', 'init_data']
GAUGE = (449721, 3728390)          # outlet gauge, UTM 12N (EPSG:26912)

# (file, panel title, color ramp), in the same order as the soils table
PANELS = [
    ('Ks_0-5cm.asc',        r'$K_s$ (mm/hr)',        'soil'),
    ('theta_s_0-5cm.asc',   r'$\theta_s$',           'soil'),
    ('theta_r_0-5cm.asc',   r'$\theta_r$',           'soil'),
    ('m_0-5cm.asc',         r'$m$',                  'soil'),
    ('psib_0-5cm.asc',      r'$\psi_b$ (mm)',        'soil'),
    ('f.asc',               r'$f$ (1/mm)',           'soil'),
    ('bedrock_depth_m.asc', 'Depth to bedrock (m)',  'rock'),
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


def _find(name, dirs=GRID_DIRS):
    for d in dirs:
        hits = glob.glob(os.path.join(d, name))
        if hits:
            return hits[0]
    raise FileNotFoundError(f'{name} not found in any of {dirs}')


def _cmap(kind):
    """Soil parameters: viridis (perceptually uniform, colorblind-safe, prints in grayscale).
    Bedrock depth: a separate yellow-orange-brown ramp so it does not read as a soil parameter."""
    if kind == 'soil':
        cm = plt.cm.viridis
    else:
        cm = LinearSegmentedColormap.from_list('rock', plt.cm.YlOrBr(np.linspace(0.15, 1.0, 256)))
    return cm.with_extremes(bad=(1, 1, 1, 0))      # cells outside the basin stay transparent


def _read_clipped(path, ws):
    """Grid clipped to the basin: NaN outside the outline and at nodata cells."""
    with rasterio.open(path) as src:
        crs = src.crs if src.crs else ws.crs
        outline = ws.to_crs(crs).geometry
        data, tr = mask(src, outline, crop=True, nodata=np.nan, filled=True)
        arr = data[0].astype('float64')
        if src.nodata is not None:
            arr[arr == src.nodata] = np.nan
    # cells whose centers fall inside the outline, whether or not they hold data
    inside = int(geometry_mask(outline, out_shape=arr.shape, transform=tr, invert=True).sum())
    return arr, plotting_extent(arr, tr), crs, inside


def _fmt(x, _):
    return f'{x:.3g}'.replace('-', '−')


def plot_soil_grids(ws, gauge=GAUGE, out='Sycamore_soil_grids_figure.png', dpi=300,
                    width_in=6.5, max_height_in=8.8, clip=(2, 98),
                    dtb_file='bedrock_depth_m.asc'):
    """ws: basin outline (GeoDataFrame). clip: percentiles that set each color scale.
    dtb_file: use 'bedrock_depth_plus1m.asc' to map the depth used in the model."""
    panels = [(dtb_file if f.startswith('bedrock') else f, t, k) for f, t, k in PANELS]

    layers = []
    for fname, title, kind in panels:
        arr, ext, crs, inside = _read_clipped(_find(fname), ws)
        print(f'{fname:<24} cells inside outline {inside:>7,}   valid {int(np.isfinite(arr).sum()):>7,}')
        layers.append((arr, ext, title, kind))

    ws_g = ws.to_crs(crs)
    xmin, ymin, xmax, ymax = ws_g.total_bounds
    pad = 0.01 * max(xmax - xmin, ymax - ymin)
    xmin, xmax, ymin, ymax = xmin - pad, xmax + pad, ymin - pad, ymax + pad
    gx, gy = gpd.GeoSeries([Point(gauge)], crs='EPSG:26912').to_crs(crs).iloc[0].coords[0]

    # layout in inches: 4 columns x 2 rows; the 8th cell is the key
    ncol, nrow = 4, 2
    mx, gap, row_gap = 0.08, 0.14, 0.12
    title_h, cb_gap, cb_h, tick_h = 0.20, 0.10, 0.08, 0.22
    aspect = (ymax - ymin) / (xmax - xmin)
    pw = (width_in - 2 * mx - (ncol - 1) * gap) / ncol
    fixed_h = 2 * mx + nrow * (title_h + cb_gap + cb_h + tick_h) + (nrow - 1) * row_gap
    pw = min(pw, (max_height_in - fixed_h) / (nrow * aspect))     # keep it on one page
    ph = pw * aspect
    W = 2 * mx + ncol * pw + (ncol - 1) * gap
    row_h = title_h + ph + cb_gap + cb_h + tick_h
    H = 2 * mx + nrow * row_h + (nrow - 1) * row_gap

    with plt.rc_context(FONT_SETTINGS):
        fig = plt.figure(figsize=(W, H), facecolor='white')

        def cell_origin(i):
            r, c = divmod(i, ncol)
            x0 = mx + c * (pw + gap)
            y_panel = H - mx - r * (row_h + row_gap) - title_h - ph
            return x0, y_panel

        for i, (arr, ext, title, kind) in enumerate(layers):
            x0, yp = cell_origin(i)
            ax = fig.add_axes([x0 / W, yp / H, pw / W, ph / H])
            lo, hi = np.nanpercentile(arr, clip)
            if not hi > lo:
                lo, hi = np.nanmin(arr), np.nanmax(arr)
            if not hi > lo:
                hi = lo + 1e-9
            im = ax.imshow(arr, extent=ext, origin='upper', cmap=_cmap(kind),
                           norm=Normalize(lo, hi), interpolation='nearest')
            ws_g.boundary.plot(ax=ax, color='#4d4d4d', linewidth=0.5)
            ax.plot(gx, gy, marker='o', ms=3.2, mfc='white', mec='black', mew=0.7,
                    ls='none', zorder=5)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_aspect('equal', adjustable='box')
            ax.axis('off')
            ax.set_title(f'({chr(97 + i)}) {title}', fontsize=8, loc='left', pad=2, color='black')

            cax = fig.add_axes([x0 / W, (yp - cb_gap - cb_h) / H, pw / W, cb_h / H])
            cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
            cb.outline.set_linewidth(0.4)
            cb.ax.xaxis.set_major_locator(MaxNLocator(3))
            cb.ax.xaxis.set_major_formatter(FuncFormatter(_fmt))
            cb.ax.tick_params(labelsize=7, length=2, width=0.4, pad=1.5, labelcolor='black')

        # key panel (8th cell)
        # the key also uses the space where the other cells have a color bar
        x0, yp = cell_origin(len(layers))
        extra = cb_gap + cb_h + tick_h - 0.02
        kh = ph + extra
        kax = fig.add_axes([x0 / W, (yp - extra) / H, pw / W, kh / H])
        kax.set_xlim(0, pw)
        kax.set_ylim(0, kh)
        kax.axis('off')
        y1, y2 = kh - 0.10, kh - 0.30
        kax.plot([0.02, 0.28], [y1, y1], color='#4d4d4d', lw=0.8)
        kax.text(0.36, y1, 'Basin boundary', fontsize=7, va='center', color='black')
        kax.plot([0.15], [y2], marker='o', ms=3.2, mfc='white', mec='black', mew=0.7, ls='none')
        kax.text(0.36, y2, 'Outlet gauge', fontsize=7, va='center', color='black')

        m_per_in = (xmax - xmin) / pw                       # map meters per inch
        km = next(k for k in (20, 10, 5, 2, 1) if k * 1000 / m_per_in <= 0.65 * pw)
        L = km * 1000 / m_per_in
        yb = 0.12
        kax.plot([0.02, 0.02 + L], [yb, yb], color='black', lw=1.0, solid_capstyle='butt')
        for xe in (0.02, 0.02 + L):
            kax.plot([xe, xe], [yb - 0.04, yb + 0.04], color='black', lw=1.0)
        kax.text(0.02 + L / 2, yb + 0.07, f'{km} km', fontsize=7, ha='center', va='bottom', color='black')
        xn = pw - 0.15
        kax.annotate('', xy=(xn, 0.42), xytext=(xn, 0.08),
                     arrowprops=dict(arrowstyle='-|>', color='black', lw=0.8, mutation_scale=7))
        kax.text(xn, 0.45, 'N', fontsize=7, ha='center', va='bottom', color='black')

        fig.savefig(out, dpi=dpi, facecolor='white')
    print(f'Saved {out}  ({W:.2f} x {H:.2f} in at {dpi} dpi)')
    return fig
