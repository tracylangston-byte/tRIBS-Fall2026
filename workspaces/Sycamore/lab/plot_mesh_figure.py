"""Figure of the Lower Sycamore Wash mesh for Section 4.3.

Two panels: (a) the whole basin, (b) a zoom at the gauge so the 100 m stream
spacing and the hillslope triangles can be seen. Run it as one notebook cell.

Inputs (GeoDataFrames in the DEM's CRS, EPSG:26912):
  tri : mesh triangles (the 'Sycamore_mesh_triangles' QA shapefile)
  ws  : basin polygon  (the 'Sycamore_boundary' shapefile)
  st  : stream lines   (the 'Sycamore_stream' shapefile)

Example (notebook folder is lab/):
  import glob, geopandas as gpd
  tri = gpd.read_file(glob.glob('data/preprocessing/*mesh*triangles*.shp')[0])
  ws  = gpd.read_file(glob.glob('data/preprocessing/*boundary*.shp')[0])
  st  = gpd.read_file(glob.glob('data/preprocessing/*stream*.shp')[0])
  fig = plot_mesh(tri, ws, st, out='Sycamore_mesh_figure.png')

The footnote under the legend is computed from the files, not typed in: basin
area from the boundary, triangle count from the triangle file, and node count
as the number of distinct triangle corners. Check them against the mesh of
record (15,792 nodes, 30,317 triangles).
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

GAUGE = (449721, 3728390)   # USGS 09510200, UTM 12N (EPSG 26912)
ZOOM_HALF_WIDTH = 1500      # m; the zoom panel is 3 km x 3 km

plt.rcParams['font.family'] = 'Liberation Sans'

RED = '#d62728'
BLUE = '#1f6fb4'


def _scalebar(ax, length_m, label, x0, y0):
    ax.plot([x0, x0 + length_m], [y0, y0], color='k', lw=2, solid_capstyle='butt')
    ax.annotate(label, (x0 + length_m / 2, y0), xytext=(0, -6), textcoords='offset points',
                ha='center', va='top', fontsize=9)


def _footnote(tri, ws):
    area = ws.area.sum() / 1e6
    corners = np.round(np.array([c for g in tri.geometry for c in g.exterior.coords[:3]]), 2)
    nodes = len(np.unique(corners, axis=0))
    text = f'Basin {area:,.1f} km²  ·  {nodes:,} nodes  ·  {len(tri):,} triangles'
    epsg = ws.crs.to_epsg() if ws.crs is not None else None
    if epsg == 26912:
        text += '  ·  UTM Zone 12N (EPSG:26912)'
    elif epsg:
        text += f'  ·  EPSG:{epsg}'
    return text


def plot_mesh(tri, ws, st, gauge=GAUGE, half=ZOOM_HALF_WIDTH,
              out='Sycamore_mesh_figure.png', dpi=300, footnote=True):
    gx, gy = gauge
    xmin, ymin, xmax, ymax = ws.total_bounds
    xlo, xhi = xmin - 500, xmax + 500
    ylo, yhi = ymin - 4000, ymax + 500           # extra room below for the scale bar

    # Axes are laid out in inches so panel (a) fits its own shape and the gap
    # between the panels stays small: the height is fixed, and panel (a)'s
    # width follows the basin's aspect ratio.
    H = 5.2                                       # panel height, in
    wa = H * (xhi - xlo) / (yhi - ylo)            # panel (a) width, in
    wb = H                                        # panel (b) is square
    left, gap, right, bottom, top = 0.1, 0.3, 0.1, 0.85, 0.35
    fig_w = left + wa + gap + wb + right
    fig_h = bottom + H + top
    fig = plt.figure(figsize=(fig_w, fig_h))
    a = fig.add_axes([left / fig_w, bottom / fig_h, wa / fig_w, H / fig_h])
    b = fig.add_axes([(left + wa + gap) / fig_w, bottom / fig_h, wb / fig_w, H / fig_h])

    # (a) whole basin
    tri.boundary.plot(ax=a, color='0.45', linewidth=0.12)
    ws.boundary.plot(ax=a, color='k', linewidth=0.8)
    st.plot(ax=a, color=BLUE, linewidth=0.7)
    a.plot(gx, gy, marker='o', markerfacecolor=RED, markersize=5, markeredgecolor='k',
           markeredgewidth=0.7, linestyle='none', zorder=5)
    a.add_patch(Rectangle((gx - half, gy - half), 2 * half, 2 * half, fill=False,
                          edgecolor=RED, linewidth=1))
    a.set_xlim(xlo, xhi)
    a.set_ylim(ylo, yhi)
    _scalebar(a, 10000, '10 km', xlo + 800, ymin - 2200)
    a.annotate('N', xy=(0.92, 0.97), xytext=(0.92, 0.86), xycoords='axes fraction',
               textcoords='axes fraction', ha='center', va='bottom', fontsize=10,
               arrowprops=dict(arrowstyle='-|>', color='k', lw=1.2))
    a.set_title('(a) Whole basin', fontsize=10, loc='left')

    # (b) zoom at the gauge
    tri.boundary.plot(ax=b, color='0.35', linewidth=0.4)
    st.plot(ax=b, color=BLUE, linewidth=1.4)
    b.plot(gx, gy, marker='o', markerfacecolor=RED, markersize=8, markeredgecolor='k',
           markeredgewidth=0.7, linestyle='none', zorder=5)
    b.set_xlim(gx - half, gx + half)
    b.set_ylim(gy - half, gy + half)
    _scalebar(b, 500, '500 m', gx - half + 150, gy - half + 280)
    b.set_title('(b) Zoom at the gauge', fontsize=10, loc='left')

    for ax in (a, b):
        ax.set_aspect('equal', adjustable='box')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(''); ax.set_ylabel('')
        for s in ax.spines.values():
            s.set_linewidth(0.6)

    handles = [Line2D([0], [0], color='0.45', lw=1, label='Triangle edges'),
               Line2D([0], [0], color=BLUE, lw=1.5, label='Stream network'),
               Line2D([0], [0], color='k', lw=1, label='Basin boundary'),
               Line2D([0], [0], marker='o', color='w', markerfacecolor=RED,
                      markeredgecolor='k', markeredgewidth=0.7, markersize=6,
                      label='USGS 09510200'),
               Rectangle((0, 0), 1, 1, fill=False, edgecolor=RED, linewidth=1,
                         label='Area shown in (b)')]
    fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False, fontsize=8.5,
               handlelength=1.6, columnspacing=1.3,
               bbox_to_anchor=(0.5, 0.36 / fig_h))
    if footnote:
        fig.text(0.5, 0.10 / fig_h, _footnote(tri, ws), ha='center', va='bottom',
                 fontsize=8.5, color='0.25')

    fig.savefig(out, dpi=dpi)
    print('wrote', out, '| triangles:', len(tri))
    return fig
