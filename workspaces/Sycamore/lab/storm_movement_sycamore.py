#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
storm_movement_sycamore.py
==========================
Where was the rain, hour by hour, and which way did it move?
(Lower Sycamore Wash, hourly MRMS radar grids, one stretch of time at a time.)

VERSION 2: the time labels under the overview map no longer run together on short windows, and the "first"/"last"
labels are no longer cut off at the east edge. Nothing else changed (same options, same CSV columns, same numbers).
For widespread frontal rain the centroid track is NOT a good tool (it is dragged by the edge of the radar grid);
use arrival_time_sycamore.py for that, and keep this script for the hourly animation, hyetograph and CSV.

WHAT IT MAKES  (three files per run, in --out-dir, all starting with the --tag)
------------------------------------------------------------------------------
  <tag>_hourly_summary.csv   one row per hour: basin-mean rain, how much of the basin was wet, where the
                             rain "centre of mass" (centroid) was, whether the storm was touching the edge of
                             the radar grid, how fast and which way it moved since the previous hour, and how
                             far it was from the USGS gauge.
  <tag>_overview.png         the whole stretch on one page: a map of the total rain with the centroid track
                             on top (colour = time), plus the hyetograph and centroid-to-gauge distance.
  <tag>_animation.gif        the hourly frames played as a movie, with a fixed colour scale, the basin
                             outline, the gauge, a short trail behind the centroid and a time cursor on
                             the hyetograph.
  <tag>_run_info.txt         the settings used and a short numeric summary (also printed on screen).

QUICK START  (run from the lab/ folder of workspaces/Sycamore)
--------------------------------------------------------------
    # 0. check that the script itself works (uses made-up storms, touches none of your data)
    python storm_movement_sycamore.py self-test

    # 1. look before you leap: lists the files it will read and writes NOTHING
    python storm_movement_sycamore.py run --dry-run

    # 2. the default run = the January 2017 flood week (Jan 18 00:00 to Jan 25 00:00 MST)
    python storm_movement_sycamore.py run

    # 3. any other stretch, e.g. only the Feb 27-28 storm (use whole hours, Arizona time)
    python storm_movement_sycamore.py run --start "2017-02-27 09:00" --end "2017-03-02 09:00" --tag feb28

WORDS USED
----------
  hourly frame   one radar grid = the rain that fell during ONE hour, in mm per hour (so a cell value of 2.0
                 means 2.0 mm fell in that hour). 32 x 49 cells, each 1 km x 1 km, EPSG:26912.
  centroid       the rain-weighted average position of the rain: add up (rain x position) over every wet cell
                 and divide by the total rain. Think of the rain as weights on a flat board; the centroid is
                 where the board would balance.
  spread         how far the rain is scattered around its centroid (rain-weighted root-mean-square distance, km).
                 A rising spread with a still centroid means the storm is growing, not moving.
  bearing        the direction the centroid moved, in degrees clockwise from north, pointing TO where it went:
                 0 = north, 90 = east, 180 = south, 270 = west. A bearing of 225 means "moving toward the
                 south-west".
  hyetograph     a bar chart of rain against time.

HOW THE CENTROID IS MADE  (so you can explain it and change it)
---------------------------------------------------------------
  1. Only cells with at least --cell-min mm/hr (default 0.1) count. This keeps faint radar noise from dragging
     the centroid around.
  2. An hour only gets a centroid if at least --min-wet-cells cells (default 20 = 20 km2) pass that test.
     Hours with less rain are left blank (NaN) on purpose: with almost no rain "where the storm is" has no
     meaning.
  3. The centroid uses the WHOLE radar grid (32 x 49 km), not just the basin, so the storm can be seen
     arriving from outside. A basin-only centroid is dragged around by the basin edge.
  4. Hours that MRMS had no data for (listed as "filled_zero_missing_hour" in manifest.csv, which sits one
     folder above rain/) are marked and never get a centroid, because a zero frame looks like a storm that
     vanished. None fall in January-March 2017.
  5. Speed and bearing come from the centroid at hour t minus the centroid at hour t-1, and only when BOTH
     hours have a centroid. Each step is a noisy 1-hour difference: look at the pattern over several hours,
     not at single rows.
  6. EDGE WARNING. If a storm is partly outside the radar grid (arriving from the west, or leaving), the grid
     only holds part of it and the centroid is pulled toward the inside of the grid. In a test the centroid
     was accurate to under 100 m while the storm stayed more than about 12 km from every edge, and was
     1-13 km off as the storm crossed the edge. Hours where more than --edge-frac (default 2%) of the wet
     rain sits in the outermost ring of cells are marked edge_flag = True in the CSV, drawn as hollow markers
     in the animation and as x marks in the overview, and left out of the speed statistics in the summary.
     This is a warning, not a correction.

THE BASIN MEAN uses the same area weighting as srp_figures/load_radar.py (share of each 1 km cell inside the
basin outline), so basin-mean values here should match basin_mean_hourly_rain_wy2017.csv.

LIMITS  (read before interpreting)
----------------------------------
  * The frames are one hour apart. A fast thunderstorm can cross the whole grid inside one frame, so cell-by-
    cell tracking is not possible with this data. Winter fronts last many hours and are the target here.
  * The model (tRIBS) is given these same hourly frames, so this shows the rain pattern the model saw.
  * The grid is only 32 x 49 km. Rain that falls outside it is invisible, so a storm's first and last hours
    (see the edge warning above) are the least reliable part of its track.
  * Over a window that holds several separate storm pulses (like Jan 18-25) the first-to-last displacement is
    NOT a storm velocity. Use the hourly table and the track, or run one storm at a time.
  * Radar read about 69-87% of the gauges in widespread storms, so trust the pattern more than the depths.

CONVENTIONS
-----------
  * Arizona time (MST, UTC-7, no daylight saving). --start and --end are the rain interval [start, end).
  * File names are radar_MMDDYYYYHH.asc, where HH is the hour that ENDS (the hour ending at midnight is written
    as 00 of the NEXT day). So the Jan 18-25 window reads the 168 files with labels Jan 18 01 .. Jan 25 00.
  * The gauge is USGS 09510200 at x 449721, y 3728390 (EPSG:26912), the value used in srp_figures/srp_common.py.
    Distances are straight-line, not along the channel.

NEEDS: numpy, pandas, matplotlib, Pillow, and ONE of pyshp (import shapefile) or geopandas to read the basin
outline (your Codespace already has geopandas for the notebooks).
"""
from __future__ import annotations

import argparse
import contextlib
import io
import math
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HOUR = pd.Timedelta(hours=1)
DEFAULT_RAIN_DIR = "~/mrms_wk_wy2017/rain"
DEFAULT_BOUNDARY = "data/preprocessing/Sycamore_boundary.shp"
GAUGE_XY = (449721.0, 3728390.0)               # USGS 09510200, EPSG:26912 (srp_common.py / data_sources.md)
FILLED_STATUS = ("filled_zero_missing_hour",)  # manifest.csv status for hours MRMS had no data
DEFAULT_START, DEFAULT_END = "2017-01-18 00:00", "2017-01-25 00:00"

CSV_COLUMNS = [
    "label_end_mst", "interval_start_mst", "grid_mean_mm_hr", "basin_mean_mm_hr", "basin_max_mm_hr",
    "basin_wet_frac", "grid_wet_cells", "filled_by_zero", "centroid_valid", "edge_rain_frac", "edge_flag",
    "centroid_x_m", "centroid_y_m", "spread_km", "dist_to_gauge_km", "step_dx_km", "step_dy_km", "step_speed_kmh", "step_bearing_deg",
]

# colours (same family as srp_figures/srp_common.py)
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
BLUE, NAVY, RED = "#2a78d6", "#0C234B", "#AB0520"
BLUE_LIGHT = "#bcd6f6"


def die(msg: str):
    raise SystemExit("ERROR: " + msg)


# ================================================================================================================
# files, grids, basin outline
# ================================================================================================================
def file_for_label(rain_dir: Path, prefix: str, label: pd.Timestamp) -> Path:
    """radar_MMDDYYYYHH.asc for the hour that ENDS at `label` (00 of the next day for the hour ending midnight)."""
    return Path(rain_dir) / f"{prefix}{label:%m%d%Y%H}.asc"


def parse_window(start: str, end: str):
    s, e = pd.Timestamp(start), pd.Timestamp(end)
    for name, t in (("--start", s), ("--end", e)):
        if t.minute or t.second or t.microsecond:
            die(f"{name} must be a whole hour (like '2017-01-18 00:00'), got {t}")
    if e <= s:
        die("--end must be later than --start")
    n = int(round((e - s) / HOUR))
    labels = pd.DatetimeIndex([s + HOUR * (i + 1) for i in range(n)])
    return s, e, labels


def read_ascii(path: Path):
    """ESRI ASCII grid -> (header dict with lower-case keys, 2-D float array, first row = NORTH row). NODATA -> NaN."""
    with open(path) as fh:
        hdr = {}
        while True:
            pos = fh.tell()
            parts = fh.readline().split()
            if len(parts) == 2 and parts[0].replace("_", "").isalpha():
                hdr[parts[0].lower()] = float(parts[1])
            else:
                fh.seek(pos)
                break
        body = fh.read().split()
    need = ("ncols", "nrows", "cellsize")
    if any(k not in hdr for k in need) or not (("xllcorner" in hdr or "xllcenter" in hdr) and
                                                ("yllcorner" in hdr or "yllcenter" in hdr)):
        raise ValueError(f"{Path(path).name}: header is missing ncols/nrows/cellsize/xll/yll")
    nx, ny = int(hdr["ncols"]), int(hdr["nrows"])
    if len(body) != nx * ny:
        raise ValueError(f"{Path(path).name}: expected {nx * ny} values, found {len(body)}")
    a = np.array(body, dtype="float64").reshape(ny, nx)
    nd = hdr.get("nodata_value")
    if nd is not None:
        a[a == nd] = np.nan
    return hdr, a


def grid_geometry(hdr: dict) -> dict:
    nx, ny, cs = int(hdr["ncols"]), int(hdr["nrows"]), float(hdr["cellsize"])
    if "xllcorner" in hdr:
        L, B = hdr["xllcorner"], hdr["yllcorner"]
    else:
        L, B = hdr["xllcenter"] - cs / 2, hdr["yllcenter"] - cs / 2
    R, T = L + nx * cs, B + ny * cs
    xc = L + (np.arange(nx) + 0.5) * cs
    yc = T - (np.arange(ny) + 0.5) * cs            # row 0 is the NORTH row
    X, Y = np.meshgrid(xc, yc)
    return dict(nx=nx, ny=ny, cs=cs, L=L, R=R, B=B, T=T, X=X, Y=Y,
                sig=(nx, ny, round(L, 3), round(B, 3), cs))


def load_stack(files):
    """Read all grids -> array (n, ny, nx) with NaN where nodata; every file must share one header."""
    hdr0, a0 = read_ascii(files[0])
    geo = grid_geometry(hdr0)
    stack = np.empty((len(files), geo["ny"], geo["nx"]), dtype="float32")
    stack[0] = a0
    for k in range(1, len(files)):
        h, a = read_ascii(files[k])
        if grid_geometry(h)["sig"] != geo["sig"]:
            die(f"{files[k].name} has a different grid header from {files[0].name}: "
                f"{grid_geometry(h)['sig']} vs {geo['sig']}")
        stack[k] = a
    return geo, stack


def _rings_pyshp(path: Path):
    import shapefile                                   # pyshp
    sf = shapefile.Reader(str(path))
    rings = []
    for shp in sf.shapes():
        pts = np.array(shp.points)
        if len(pts) == 0:
            continue
        parts = list(shp.parts) + [len(pts)]
        for i in range(len(parts) - 1):
            rings.append(pts[parts[i]:parts[i + 1]])
    return rings


def _rings_geopandas(path: Path):
    import geopandas as gpd
    rings = []
    for geom in gpd.read_file(path).geometry:
        if geom is None or geom.is_empty or geom.geom_type not in ("Polygon", "MultiPolygon"):
            continue
        for p in (list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]):
            rings.append(np.asarray(p.exterior.coords)[:, :2])
            for h in p.interiors:
                rings.append(np.asarray(h.coords)[:, :2])
    return rings


def load_boundary_rings(path, backend: str = "auto"):
    """Basin outline as a list of (N, 2) arrays in EPSG:26912. Refuses an outline whose .prj is not UTM 12N."""
    path = Path(path).expanduser()
    if not path.exists():
        die(f"basin outline not found: {path}  (run from lab/, or pass --boundary path/to/Sycamore_boundary.shp)")
    prj = path.with_suffix(".prj")
    if prj.exists():
        w = prj.read_text()
        if "UTM zone 12N" not in w and "UTM_Zone_12N" not in w and "26912" not in w:
            die(f"{prj.name} is not NAD83 / UTM zone 12N; reproject the outline first.")
    else:
        print(f"  warning: no .prj next to {path.name}; assuming EPSG:26912")
    order = ["pyshp", "geopandas"] if backend == "auto" else [backend]
    for b in order:
        try:
            rings = _rings_pyshp(path) if b == "pyshp" else _rings_geopandas(path)
        except ImportError:
            continue
        if not rings:
            die(f"no polygon found in {path.name}")
        return rings
    die("cannot read the basin outline: install pyshp (pip install pyshp) or geopandas")


def grid_cover(rings, geo: dict, ss: int = 50):
    """Fraction (0-1) of each cell that lies inside the basin outline. Same method as srp_common.grid_cover:
    the polygon is drawn on a grid `ss` times finer (even-odd rule) and averaged back down."""
    from PIL import Image, ImageDraw
    L, R, B, T, nx, ny = geo["L"], geo["R"], geo["B"], geo["T"], geo["nx"], geo["ny"]
    fine = np.zeros((ny * ss, nx * ss), dtype=bool)
    for r in rings:
        if len(r) < 3:
            continue
        px = (r[:, 0] - L) / (R - L) * nx * ss
        py = (T - r[:, 1]) / (T - B) * ny * ss
        img = Image.new("1", (nx * ss, ny * ss), 0)
        ImageDraw.Draw(img).polygon(list(zip(px.tolist(), py.tolist())), fill=1)
        fine ^= np.asarray(img, dtype=bool)
    return fine.reshape(ny, ss, nx, ss).mean(axis=(1, 3))


def read_filled_keys(manifest: Path):
    """Set of 'YYYY-MM-DD HH:MM' labels MRMS had no data for, or None if the manifest is missing/unreadable."""
    manifest = Path(manifest).expanduser()
    if not manifest.exists():
        return None
    try:
        m = pd.read_csv(manifest)
        lab = pd.to_datetime(m["label_local"])
        sel = m["status"].isin(FILLED_STATUS)
    except Exception:
        return None
    return {f"{t:%Y-%m-%d %H:%M}" for t in lab[sel]}


# ================================================================================================================
# the numbers
# ================================================================================================================
def compass(deg: float) -> str:
    if deg is None or not np.isfinite(deg):
        return "n/a"
    return ["N", "NE", "E", "SE", "S", "SW", "W", "NW"][int(((deg % 360) + 22.5) // 45) % 8]


def compute_table(labels, rain, geo, cover, filled_keys, cell_min, min_wet, edge_thresh=0.02):
    """One row per hour. `rain` is (n, ny, nx) mm/hr with NaN already replaced by 0."""
    X, Y = geo["X"], geo["Y"]
    wb = cover / cover.sum()
    inside = cover >= 0.5
    gx, gy = GAUGE_XY
    ring = np.zeros(X.shape, dtype=bool)                  # the outermost ring of cells of the radar grid
    ring[0, :] = ring[-1, :] = ring[:, 0] = ring[:, -1] = True
    rows, prev = [], None
    for k, lab in enumerate(labels):
        r = rain[k].astype("float64")
        wet = r >= cell_min
        n_wet = int(wet.sum())
        filled = f"{lab:%Y-%m-%d %H:%M}" in filled_keys
        valid = bool(n_wet >= min_wet and not filled)
        cx = cy = spread = dist = np.nan
        w_all = np.where(wet, r, 0.0)
        tot = float(w_all.sum())
        edge_frac = float(w_all[ring].sum() / tot) if tot > 0 else np.nan
        edge_flag = bool(np.isfinite(edge_frac) and edge_frac > edge_thresh)
        if valid:
            w = w_all
            S = w.sum()
            cx, cy = float((w * X).sum() / S), float((w * Y).sum() / S)
            spread = math.sqrt(float((w * ((X - cx) ** 2 + (Y - cy) ** 2)).sum() / S)) / 1000.0
            dist = math.hypot(cx - gx, cy - gy) / 1000.0
        dx = dy = sp = br = np.nan
        if valid and prev is not None:                    # previous HOUR also had a centroid
            dx, dy = (cx - prev[0]) / 1000.0, (cy - prev[1]) / 1000.0
            sp = math.hypot(dx, dy)                       # frames are exactly 1 h apart -> km per hour
            if sp > 1e-9:
                br = (math.degrees(math.atan2(dx, dy)) + 360.0) % 360.0
        prev = (cx, cy) if valid else None
        rows.append({
            "label_end_mst": lab, "interval_start_mst": lab - HOUR,
            "grid_mean_mm_hr": float(r.mean()), "basin_mean_mm_hr": float((r * wb).sum()),
            "basin_max_mm_hr": float(r[inside].max()) if inside.any() else np.nan,
            "basin_wet_frac": float((wet * cover).sum() / cover.sum()), "grid_wet_cells": n_wet,
            "filled_by_zero": filled, "centroid_valid": valid, "edge_rain_frac": edge_frac, "edge_flag": edge_flag,
            "centroid_x_m": cx, "centroid_y_m": cy,
            "spread_km": spread, "dist_to_gauge_km": dist, "step_dx_km": dx, "step_dy_km": dy,
            "step_speed_kmh": sp, "step_bearing_deg": br,
        })
    return pd.DataFrame(rows, columns=CSV_COLUMNS)


def summarize(df: pd.DataFrame, area_km2: float, a) -> list[str]:
    n = len(df)
    L = []
    bm = df["basin_mean_mm_hr"]
    L.append(f"window (rain intervals, MST): {df['interval_start_mst'].iloc[0]:%Y-%m-%d %H:%M} to "
             f"{df['label_end_mst'].iloc[-1]:%Y-%m-%d %H:%M}  ({n} hourly frames)")
    L.append(f"basin area from 1 km coverage: {area_km2:.1f} km2")
    L.append(f"basin-mean rain over the window: {bm.sum():.1f} mm   (grid-mean {df['grid_mean_mm_hr'].sum():.1f} mm)")
    top = df.loc[bm.idxmax()]
    L.append(f"wettest basin-mean hour: {top['interval_start_mst']:%b %d %H:%M}-{top['label_end_mst']:%H:%M} "
             f"= {bm.max():.2f} mm/hr;  hours with basin mean >= {a.cell_min:g} mm/hr: {int((bm >= a.cell_min).sum())}")
    nf = int(df["filled_by_zero"].sum())
    if nf:
        L.append(f"hours MRMS had no data (zero-filled): {nf}  -> {', '.join(f'{t:%b %d %H:%M}' for t in df.loc[df['filled_by_zero'], 'label_end_mst'])}")
    v = df[df["centroid_valid"]]
    L.append(f"hours with a centroid: {len(v)} of {n}  (cell-min {a.cell_min:g} mm/hr, at least {a.min_wet_cells} wet cells)")
    n_edge = int(v["edge_flag"].sum())
    if len(v):
        L.append(f"  of those, {n_edge} are edge-flagged (more than {a.edge_frac:.0%} of the wet rain in the outermost ring of "
                 f"grid cells: storm partly off the grid, centroid pulled inward)")
    if len(v) >= 2:
        f, l = v.iloc[0], v.iloc[-1]
        dx, dy = (l["centroid_x_m"] - f["centroid_x_m"]) / 1000, (l["centroid_y_m"] - f["centroid_y_m"]) / 1000
        dist = math.hypot(dx, dy)
        br = (math.degrees(math.atan2(dx, dy)) + 360) % 360 if dist > 1e-9 else float("nan")
        L.append(f"first -> last centroid ({f['label_end_mst']:%b %d %H:%M} -> {l['label_end_mst']:%b %d %H:%M}): "
                 f"{dist:.1f} km toward {compass(br)} ({br:.0f} deg)   [a straight line between two points; see LIMITS]")
        i = v["dist_to_gauge_km"].idxmin()
        L.append(f"centroid distance to gauge: {f['dist_to_gauge_km']:.1f} km at the start, {l['dist_to_gauge_km']:.1f} km at the end, "
                 f"closest {v.loc[i, 'dist_to_gauge_km']:.1f} km at {v.loc[i, 'label_end_mst']:%b %d %H:%M}")
        clear = ~(df["edge_flag"] | df["edge_flag"].shift(1, fill_value=False))      # neither end of the step touches the edge
        sp = df.loc[clear, "step_speed_kmh"].dropna()
        if len(sp):
            L.append(f"hour-to-hour centroid steps with both hours clear of the grid edge: {len(sp)} usable; median {sp.median():.1f} km/h, "
                     f"90th percentile {sp.quantile(0.9):.1f} km/h, largest {sp.max():.1f} km/h  (noisy: single hours jump)")
        else:
            L.append("no hour-to-hour step has both hours clear of the grid edge, so no speed statistics are given (see the CSV)")
        if bool(f["edge_flag"]) or bool(l["edge_flag"]):
            L.append("  WARNING: the first and/or last centroid is edge-flagged, so the first -> last line above is biased")
    else:
        L.append("fewer than two hours have a centroid, so no track or speed was computed")
    return L


# ================================================================================================================
# figures
# ================================================================================================================
def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "axes.edgecolor": GRID, "axes.spines.top": False, "axes.spines.right": False,
        "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False,
    })
    return plt


def _tnum(series):
    import matplotlib.dates as mdates
    return mdates.date2num(series.to_numpy("datetime64[ns]"))


def _time_ticks(lo, hi, max_ticks=5, two_line=True):
    """Tick positions (matplotlib date numbers) and labels for a time colour bar between `lo` and `hi`, chosen so the
    labels never run together: the smallest 'round' step (1, 2, 3, 6, 12 hours, 1, 2, 3, 7 days ...) that gives at most
    `max_ticks` ticks, with ticks on whole multiples of the step counted from midnight of the first day."""
    import matplotlib.dates as mdates
    day0 = math.floor(lo)
    pos, step = [], 1
    for step in (1, 2, 3, 6, 12, 24, 48, 72, 168, 336, 720):
        d = step / 24.0
        k0 = math.ceil((lo - day0) / d - 1e-9)
        k1 = math.floor((hi - day0) / d + 1e-9)
        pos = [day0 + k * d for k in range(k0, k1 + 1)]
        if len(pos) <= max_ticks:
            break
    if len(pos) < 2:                                   # window shorter than one step: label the two ends instead
        pos, step = [lo, hi], 1
    fmt = ("%b %d\n%H:%M" if two_line else "%b %d %H:%M") if step < 24 else "%b %d"
    return pos, [mdates.num2date(p).strftime(fmt) for p in pos]


def _km_axes(ax):
    from matplotlib.ticker import FuncFormatter
    f = FuncFormatter(lambda v, _: f"{v / 1000:.0f}")
    ax.xaxis.set_major_formatter(f)
    ax.yaxis.set_major_formatter(f)
    ax.set_xlabel("UTM 12N easting (km)")
    ax.set_ylabel("UTM 12N northing (km)")


def _draw_outline(ax, rings):
    for r in rings:
        ax.plot(r[:, 0], r[:, 1], color=NAVY, lw=1.3, zorder=4)


def make_overview(df, total_mm, geo, rings, out_png, title):
    plt = _mpl()
    import matplotlib.dates as mdates
    fig = plt.figure(figsize=(11.5, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.55], hspace=0.32, wspace=0.28)
    axm = fig.add_subplot(gs[:, 0])
    axh = fig.add_subplot(gs[0, 1])
    axd = fig.add_subplot(gs[1, 1], sharex=axh)
    ext = [geo["L"], geo["R"], geo["B"], geo["T"]]
    im = axm.imshow(total_mm, extent=ext, origin="upper", cmap="Blues", vmin=0, vmax=max(float(total_mm.max()), 1.0))
    cb = fig.colorbar(im, ax=axm, fraction=0.046, pad=0.02)
    cb.set_label("rain over the window (mm)")
    _draw_outline(axm, rings)
    axm.plot(*GAUGE_XY, marker="*", ms=13, color=RED, mec="white", mew=0.8, ls="none", zorder=6)
    axm.annotate("gauge", xy=GAUGE_XY, xytext=(8, -4), textcoords="offset points", fontsize=8, color=RED, zorder=9)
    v = df[df["centroid_valid"]]
    if len(v):
        tn = _tnum(v["label_end_mst"])
        full_t = _tnum(df["label_end_mst"])
        # lines only between hours that are consecutive AND both valid
        for i in range(1, len(df)):
            if df["centroid_valid"].iat[i] and df["centroid_valid"].iat[i - 1]:
                axm.plot(df["centroid_x_m"].iloc[i - 1:i + 1], df["centroid_y_m"].iloc[i - 1:i + 1],
                         color=INK2, lw=0.8, alpha=0.7, zorder=5)
        sc = axm.scatter(v["centroid_x_m"], v["centroid_y_m"], c=tn, cmap="viridis", s=34, zorder=7,
                         vmin=full_t.min(), vmax=full_t.max(), edgecolor="white", linewidth=0.6)
        cbt = fig.colorbar(sc, ax=axm, orientation="horizontal", fraction=0.04, pad=0.10)
        tk_pos, tk_lab = _time_ticks(float(full_t.min()), float(full_t.max()))     # few, well-spaced labels (never crowd)
        cbt.set_ticks(tk_pos)
        cbt.set_ticklabels(tk_lab)
        cbt.ax.tick_params(labelsize=7.5)
        cbt.set_label("time of rain centroid (hour ending, MST)")
        ve = v[v["edge_flag"]]
        if len(ve):
            axm.scatter(ve["centroid_x_m"], ve["centroid_y_m"], marker="x", s=26, color=INK, zorder=8, linewidths=1.0,
                        label="x = storm touching the grid edge\n(centroid biased inward)")
            axm.legend(loc="lower right", fontsize=7.5)
        for idx, word, va in ((0, "first", "bottom"), (-1, "last", "top")):
            px, py = float(v["centroid_x_m"].iloc[idx]), float(v["centroid_y_m"].iloc[idx])
            east = px > geo["L"] + 0.7 * (geo["R"] - geo["L"])               # near the east edge: write the word to the left
            axm.text(px, py, f"{word}  " if east else f"  {word}", fontsize=8, color=INK, zorder=9, va=va,
                     ha="right" if east else "left",
                     bbox=dict(boxstyle="round,pad=0.12", fc=SURFACE, ec="none", alpha=0.85))
    else:
        axm.text(0.5, 0.5, "no hour had enough rain\nfor a centroid", transform=axm.transAxes, ha="center")
    axm.set_xlim(geo["L"], geo["R"])
    axm.set_ylim(geo["B"], geo["T"])
    _km_axes(axm)
    axm.set_title("Total rain and where the rain centroid went")
    # hyetograph
    x0 = _tnum(df["interval_start_mst"])
    axh.bar(x0, df["basin_mean_mm_hr"], width=1 / 24, align="edge", color=BLUE, lw=0)
    for t in x0[df["filled_by_zero"].to_numpy()]:
        axh.axvspan(t, t + 1 / 24, color=MUTED, alpha=0.5, lw=0)
    axh.set_ylabel("basin-mean rain (mm/hr)")
    axh.set_title("Hyetograph (grey band = MRMS had no data that hour)")
    axh.set_axisbelow(True)
    axh.grid(color=GRID)
    axh.tick_params(labelbottom=False)
    axh.xaxis.offsetText.set_visible(False)
    # distance to gauge
    xe = _tnum(df["label_end_mst"])
    axd.plot(xe, df["dist_to_gauge_km"], marker="o", ms=3.5, lw=1.0, color=NAVY)
    de = df[df["edge_flag"] & df["centroid_valid"]]
    if len(de):
        axd.plot(_tnum(de["label_end_mst"]), de["dist_to_gauge_km"], ls="none", marker="x", ms=7, color=INK)
    axd.set_ylabel("centroid to gauge (km)")
    axd.set_title("Rain centroid distance to the gauge (down = toward the outlet)")
    axd.set_axisbelow(True)
    axd.grid(color=GRID)
    loc = mdates.AutoDateLocator()
    axd.xaxis.set_major_locator(loc)
    axd.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
    fig.suptitle(title, fontsize=12, fontweight="bold", x=0.02, ha="left")
    fig.savefig(out_png, dpi=170, bbox_inches="tight")
    plt.close(fig)


def make_gif(df, rain, geo, rings, out_gif, a, title):
    plt = _mpl()
    from matplotlib import colors
    from matplotlib.animation import FuncAnimation, PillowWriter
    n = len(df)
    fig = plt.figure(figsize=(10, 4.9), dpi=80)
    gs = fig.add_gridspec(1, 2, width_ratios=[0.9, 1.5], wspace=0.45, left=0.08, right=0.97, bottom=0.2, top=0.86)
    axm, axh = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    cmap = plt.get_cmap("Blues").copy()
    cmap.set_bad(alpha=0.0)
    norm = colors.PowerNorm(gamma=0.5, vmin=0.0, vmax=a.vmax)
    ext = [geo["L"], geo["R"], geo["B"], geo["T"]]
    axm.set_facecolor("#f1f0ec")
    im = axm.imshow(np.ma.masked_less(rain[0], 0.05), extent=ext, origin="upper", cmap=cmap, norm=norm, interpolation="nearest")
    cb = fig.colorbar(im, ax=axm, fraction=0.05, pad=0.03)
    cb.set_ticks([t for t in (0.1, 0.5, 1, 2, 4, 6, 8, 12) if t <= a.vmax])
    cb.set_label("mm/hr")
    _draw_outline(axm, rings)
    axm.plot(*GAUGE_XY, marker="*", ms=12, color=RED, mec="white", mew=0.8, ls="none", zorder=6)
    axm.annotate("gauge", xy=GAUGE_XY, xytext=(8, -4), textcoords="offset points", fontsize=8, color=RED, zorder=9)
    trail, = axm.plot([], [], color=NAVY, lw=1.2, marker="o", ms=3, zorder=7)
    now, = axm.plot([], [], marker="o", ms=9, color="#f2a900", mec=INK, mew=1.0, ls="none", zorder=8)
    axm.set_xlim(geo["L"], geo["R"])
    axm.set_ylim(geo["B"], geo["T"])
    _km_axes(axm)
    ttl = axm.set_title("")
    x0 = _tnum(df["interval_start_mst"])
    bars = axh.bar(x0, df["basin_mean_mm_hr"], width=1 / 24, align="edge", color=BLUE_LIGHT, lw=0)
    for t in x0[df["filled_by_zero"].to_numpy()]:
        axh.axvspan(t, t + 1 / 24, color=MUTED, alpha=0.5, lw=0)
    cursor = axh.axvline(x0[0], color=RED, lw=1.2)
    axh.set_ylabel("basin-mean rain (mm/hr)")
    axh.set_title("Basin-mean hyetograph")
    axh.set_axisbelow(True)
    axh.grid(color=GRID)
    import matplotlib.dates as mdates
    loc = mdates.AutoDateLocator()
    axh.xaxis.set_major_locator(loc)
    axh.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
    fig.suptitle(title, fontsize=11, fontweight="bold", x=0.02, ha="left")
    cx, cy, valid = df["centroid_x_m"].to_numpy(), df["centroid_y_m"].to_numpy(), df["centroid_valid"].to_numpy()
    edge = df["edge_flag"].to_numpy()
    state = {"prev": 0, "printed": -1}

    def update(k):
        im.set_data(np.ma.masked_less(rain[k], 0.05))
        lo = max(0, k - a.trail_hours)
        tx = np.where(valid[lo:k + 1], cx[lo:k + 1], np.nan)
        ty = np.where(valid[lo:k + 1], cy[lo:k + 1], np.nan)
        trail.set_data(tx, ty)
        if valid[k]:
            now.set_data([cx[k]], [cy[k]])
            now.set_markerfacecolor("none" if edge[k] else "#f2a900")
            extra = f"centroid {df['dist_to_gauge_km'].iat[k]:.1f} km from gauge" + ("  (touches grid edge)" if edge[k] else "")
        else:
            now.set_data([], [])
            extra = "MRMS had no data this hour (zero-filled)" if df["filled_by_zero"].iat[k] else "too little rain for a centroid"
        t0, t1 = df["interval_start_mst"].iat[k], df["label_end_mst"].iat[k]
        ttl.set_text(f"{t0:%a %b %d %H:%M} to {t1:%H:%M} MST")
        axm.set_xlabel(f"basin mean {df['basin_mean_mm_hr'].iat[k]:.2f} mm/hr\n{extra}", fontsize=8)
        bars[state["prev"]].set_color(BLUE_LIGHT)
        bars[k].set_color(NAVY)
        state["prev"] = k
        cursor.set_xdata([x0[k] + 0.5 / 24, x0[k] + 0.5 / 24])
        if k % 24 == 0 and state["printed"] != k:
            state["printed"] = k
            print(f"    animation frame {k + 1}/{n}")
        return []

    anim = FuncAnimation(fig, update, frames=n, interval=1000 / a.fps, blit=False)
    anim.save(str(out_gif), writer=PillowWriter(fps=a.fps))
    plt.close(fig)


# ================================================================================================================
# the run command
# ================================================================================================================
def cmd_run(a) -> int:
    start, end, labels = parse_window(a.start, a.end)
    rain_dir = Path(a.rain_dir).expanduser()
    files = [file_for_label(rain_dir, a.prefix, lab) for lab in labels]
    tag = a.tag or f"{start:%Y%m%d}_{end:%Y%m%d}"
    out_dir = Path(a.out_dir).expanduser()
    manifest = Path(a.manifest).expanduser() if a.manifest else rain_dir.parent / "manifest.csv"
    n = len(labels)
    print("storm_movement_sycamore: " + ("DRY RUN (nothing will be written)" if a.dry_run else "run"))
    print(f"  window: rain from {start:%Y-%m-%d %H:%M} to {end:%Y-%m-%d %H:%M} MST  ->  {n} hourly files, labels "
          f"{labels[0]:%Y-%m-%d %H:%M} .. {labels[-1]:%Y-%m-%d %H:%M} (end of hour)")
    print(f"  rain folder: {rain_dir}")
    missing = [f for f in files if not f.exists()]
    print(f"  files present: {n - len(missing)} of {n}" + ("" if not missing else f"   <-- {len(missing)} MISSING, first: {missing[0].name}"))
    if missing:
        if a.dry_run:
            print("  (a real run would stop here)")
        else:
            die(f"{len(missing)} radar files are missing (first {missing[0].name}, last {missing[-1].name}). "
                f"Check --rain-dir, --prefix and the window.")
    geo = None
    if len(missing) < n:
        first = next(f for f in files if f.exists())
        hdr, _ = read_ascii(first)
        geo = grid_geometry(hdr)
        print(f"  grid: {geo['nx']} x {geo['ny']} cells of {geo['cs']:g} m, x {geo['L']:.0f}-{geo['R']:.0f}, y {geo['B']:.0f}-{geo['T']:.0f}")
    rings = load_boundary_rings(a.boundary)
    cover = None
    if geo is not None:
        cover = grid_cover(rings, geo)
        area = float(cover.sum()) * geo["cs"] ** 2 / 1e6
        print(f"  basin outline: {a.boundary}  ->  {area:.1f} km2 on the radar grid")
        if cover.sum() <= 0:
            die("the basin outline does not overlap the radar grid (wrong coordinate system?)")
        if a.expect_area_km2 and abs(area - a.expect_area_km2) > 0.03 * a.expect_area_km2:
            print(f"  WARNING: expected about {a.expect_area_km2:g} km2 (the DEM-derived basin is 424.6 km2); is this the right outline?")
    filled = read_filled_keys(manifest)
    if filled is None:
        print(f"  manifest: not found at {manifest} -> zero-filled hours cannot be flagged "
              f"(none are expected in Jan-Mar 2017; use --manifest to point at it)")
        filled = set()
    else:
        here = sum(1 for lab in labels if f"{lab:%Y-%m-%d %H:%M}" in filled)
        print(f"  manifest: {manifest}  ->  {len(filled)} zero-filled hours in the whole file, {here} inside this window")
    print(f"  outputs ({out_dir}/):  {tag}_hourly_summary.csv  {tag}_overview.png  "
          + ("" if a.no_gif else f"{tag}_animation.gif  ") + f"{tag}_run_info.txt")
    if a.dry_run:
        print("\nDry run finished: nothing was computed or written.")
        return 0

    print(f"\nreading {n} grids ...")
    geo, stack = load_stack(files)
    n_nan = int(np.isnan(stack).sum())
    if n_nan:
        print(f"  note: {n_nan} NODATA cells were set to 0")
    rain = np.clip(np.nan_to_num(stack, nan=0.0), 0.0, None)
    area = float(cover.sum()) * geo["cs"] ** 2 / 1e6
    df = compute_table(labels, rain, geo, cover, filled, a.cell_min, a.min_wet_cells, a.edge_frac)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path, png_path, gif_path, txt_path = (out_dir / f"{tag}_{s}" for s in
                                              ("hourly_summary.csv", "overview.png", "animation.gif", "run_info.txt"))
    out = df.copy()
    out["label_end_mst"] = out["label_end_mst"].dt.strftime("%Y-%m-%d %H:%M")
    out["interval_start_mst"] = out["interval_start_mst"].dt.strftime("%Y-%m-%d %H:%M")
    out.round({"grid_mean_mm_hr": 4, "basin_mean_mm_hr": 4, "basin_max_mm_hr": 3, "basin_wet_frac": 4, "edge_rain_frac": 4,
               "centroid_x_m": 1, "centroid_y_m": 1, "spread_km": 3, "dist_to_gauge_km": 3,
               "step_dx_km": 3, "step_dy_km": 3, "step_speed_kmh": 3, "step_bearing_deg": 1}).to_csv(csv_path, index=False)
    title = f"Lower Sycamore Wash: hourly MRMS rain, {start:%b %d} to {end:%b %d %Y} (MST)"
    print("drawing the overview figure ...")
    make_overview(df, rain.sum(axis=0), geo, rings, png_path, title)
    if not a.no_gif:
        print("drawing the animation (about a second per 10 frames) ...")
        make_gif(df, rain, geo, rings, gif_path, a, title)
    lines = summarize(df, area, a)
    settings = (f"settings: rain-dir {rain_dir} | boundary {a.boundary} | cell-min {a.cell_min:g} | min-wet-cells {a.min_wet_cells} | "
                f"edge-frac {a.edge_frac:g} | trail-hours {a.trail_hours} | vmax {a.vmax:g} | fps {a.fps} | manifest {'used' if filled else 'not used / nothing flagged'}")
    txt_path.write_text("\n".join(lines + ["", settings]) + "\n")
    print("\n" + "\n".join("  " + s for s in lines))
    print(f"\nwrote to {out_dir}/ :")
    for p in (csv_path, png_path, txt_path) + (() if a.no_gif else (gif_path,)):
        print(f"  {p.name}  ({p.stat().st_size / 1024:.0f} KB)")
    return 0


# ================================================================================================================
# self-test (made-up storms with known answers; touches none of your data)
# ================================================================================================================
def _write_asc(path: Path, arr, nodata=-9999):
    ny, nx = arr.shape
    with open(path, "w") as f:
        f.write(f"ncols {nx}\nnrows {ny}\nxllcorner 440000\nyllcorner 3720000\ncellsize 1000\nNODATA_value {nodata}\n")
        np.savetxt(f, arr, fmt="%.3f")


def _write_test_boundary(path: Path, x0, y0, x1, y1):
    ring = [(x0, y0), (x0, y1), (x1, y1), (x1, y0), (x0, y0)]          # clockwise = outer ring
    try:
        import shapefile
        w = shapefile.Writer(str(path))
        w.field("id", "N")
        w.poly([ring])
        w.record(1)
        w.close()
    except ImportError:
        import geopandas as gpd
        from shapely.geometry import Polygon
        gpd.GeoDataFrame({"id": [1]}, geometry=[Polygon(ring)], crs="EPSG:26912").to_file(path)
    path.with_suffix(".prj").write_text(
        'PROJCS["NAD_1983_UTM_Zone_12N",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",'
        'SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
        'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],'
        'PARAMETER["Central_Meridian",-111.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],'
        'UNIT["Meter",1.0]]')


def cmd_self_test(_a) -> int:
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("  PASS  " if cond else "  FAIL  ") + msg)
        ok = ok and bool(cond)

    tmp = Path(tempfile.mkdtemp(prefix="stormmove_selftest_"))
    try:
        print("storm_movement_sycamore SELF-TEST (made-up storms; your data is not touched)")
        rain_dir = tmp / "wk" / "rain"
        rain_dir.mkdir(parents=True)
        # ---- 1. names and window --------------------------------------------------------------------------
        s, e, labels = parse_window("2017-01-19 18:00", "2017-01-20 07:00")
        check(len(labels) == 13 and labels[0] == pd.Timestamp("2017-01-19 19:00") and labels[-1] == pd.Timestamp("2017-01-20 07:00"),
              "window 18:00 -> 07:00 gives 13 hourly labels, first 19:00 and last 07:00 (end of hour)")
        check(file_for_label(Path("x"), "radar_", pd.Timestamp("2017-01-20 00:00")).name == "radar_0120201700.asc",
              "the hour ending at midnight is named with hour 00 of the NEXT day (radar_0120201700.asc)")
        check(file_for_label(Path("x"), "radar_", pd.Timestamp("2017-01-18 01:00")).name == "radar_0118201701.asc",
              "label 01:00 on Jan 18 is radar_0118201701.asc")
        full = len(parse_window("2017-01-18 00:00", "2017-01-25 00:00")[2])
        check(full == 168, "the default Jan 18-25 window reads 168 files")
        # ---- 2. synthetic grids ---------------------------------------------------------------------------
        geo0 = grid_geometry({"ncols": 32, "nrows": 49, "xllcorner": 440000, "yllcorner": 3720000, "cellsize": 1000})
        X, Y = geo0["X"], geo0["Y"]
        sigma, v_kmh, brg = 2500.0, 1.6, 225.0
        comp = v_kmh * 1000 * math.sqrt(0.5)                             # metres per hour in each axis (toward SW)
        x_start, y_start = 460000.0, 3758000.0

        def blob(k):
            xk, yk = x_start - comp * k, y_start - comp * k
            return xk, yk, 5.0 * np.exp(-((X - xk) ** 2 + (Y - yk) ** 2) / (2 * sigma ** 2))

        block = np.zeros_like(X)
        inside = (X > 452000) & (X < 462000) & (Y > 3744000) & (Y < 3754000)
        block[inside] = 4.0
        frames = {}
        for k in range(13):
            if k == 8:
                frames[k] = np.full_like(X, 0.05)                         # faint, uniform: below the 0.1 floor
            elif k == 12:
                frames[k] = block                                         # uniform block exactly over the test basin
            else:
                frames[k] = blob(k)[2]
        for k, lab in enumerate(labels):
            _write_asc(file_for_label(rain_dir, "radar_", lab), frames[k])
        (tmp / "wk" / "manifest.csv").write_text(
            "label_local,status\n" + "\n".join(
                f"{lab:%Y-%m-%d %H:%M},{'filled_zero_missing_hour' if k == 9 else 'ok'}" for k, lab in enumerate(labels)) + "\n")
        shp = tmp / "test_basin.shp"
        _write_test_boundary(shp, 452000, 3744000, 462000, 3754000)
        # ---- 3. reading -----------------------------------------------------------------------------------
        hdr, arr = read_ascii(file_for_label(rain_dir, "radar_", labels[0]))
        g = grid_geometry(hdr)
        check(arr.shape == (49, 32) and g["T"] == 3769000 and abs(g["Y"][0, 0] - 3768500) < 1e-6 and abs(g["X"][0, -1] - 471500) < 1e-6,
              "grid is 32 x 49; row 0 is the NORTH row (y 3768500) and the last column is the EAST column (x 471500)")
        k_peak = int(np.nanargmax(arr))
        check(abs(g["X"].ravel()[k_peak] - blob(0)[0]) < 1000 and abs(g["Y"].ravel()[k_peak] - blob(0)[1]) < 1000,
              "the brightest cell of the first frame sits where the blob was put (so the grid is not flipped)")
        # ---- 4. boundary backends -------------------------------------------------------------------------
        covers = {}
        for b in ("pyshp", "geopandas"):
            try:
                rings = load_boundary_rings(shp, backend=b)
                covers[b] = grid_cover(rings, g)
            except SystemExit:
                print(f"  SKIP  boundary reader '{b}' is not installed here")
        check(len(covers) >= 1, "at least one boundary reader (pyshp or geopandas) works")
        if len(covers) == 2:
            check(np.allclose(covers["pyshp"], covers["geopandas"]), "pyshp and geopandas give the same coverage grid")
        cov = next(iter(covers.values()))
        check(abs(cov.sum() - 100) < 2.0, f"a 10 km x 10 km test basin covers {cov.sum():.2f} cells (about 100)")
        # ---- 5. dry run writes nothing --------------------------------------------------------------------
        out1 = tmp / "out1"
        base = ["--rain-dir", str(rain_dir), "--boundary", str(shp), "--start", "2017-01-19 18:00", "--end", "2017-01-20 07:00",
                "--out-dir", str(out1), "--tag", "t", "--expect-area-km2", "100"]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = main(["run", "--dry-run"] + base)
        check(rc == 0 and not out1.exists() and "files present: 13 of 13" in buf.getvalue() and "inside this window" in buf.getvalue(),
              "--dry-run finds 13 of 13 files and the manifest, and writes nothing")
        # ---- 6. missing files stop a real run --------------------------------------------------------------
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                main(["run", "--no-gif"] + base[:6] + ["--end", "2017-01-20 08:00"] + base[8:])
            stopped = False
        except SystemExit as ex:
            stopped = "missing" in str(ex)
        check(stopped, "a window that needs a file that is not there stops with a 'missing' error")
        # ---- 7. full run ----------------------------------------------------------------------------------
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = main(["run", "--fps", "4"] + base)
        csvp, pngp, gifp, txtp = (out1 / f"t_{s}" for s in ("hourly_summary.csv", "overview.png", "animation.gif", "run_info.txt"))
        check(rc == 0 and all(p.exists() and p.stat().st_size > 0 for p in (csvp, pngp, gifp, txtp)),
              "a full run writes the CSV, the overview PNG, the GIF and the run-info text")
        try:
            from PIL import Image
            check(Image.open(gifp).n_frames == 13, "the GIF has one frame per hour (13)")
        except Exception as ex:                                            # pragma: no cover
            check(False, f"could not open the GIF: {ex}")
        df = pd.read_csv(csvp, parse_dates=["label_end_mst", "interval_start_mst"])
        check(list(df.columns) == CSV_COLUMNS and len(df) == 13, "the CSV has the documented columns and 13 rows")
        # ---- 8. the numbers --------------------------------------------------------------------------------
        errs = []
        for k in (0, 1, 4, 7, 10, 11):
            xk, yk, _ = blob(k)
            errs.append(max(abs(df["centroid_x_m"].iat[k] - xk), abs(df["centroid_y_m"].iat[k] - yk)))
        check(max(errs) < 250, f"centroids land on the true blob centres (worst error {max(errs):.0f} m, limit 250 m)")
        sp = df.loc[[1, 2, 5, 7, 11], "step_speed_kmh"]
        bb = df.loc[[1, 2, 5, 7, 11], "step_bearing_deg"]
        check(float((sp - v_kmh).abs().max()) < 0.1, f"hour-to-hour speed is {v_kmh} km/h (got {sp.min():.3f}-{sp.max():.3f})")
        check(float((bb - brg).abs().max()) < 3, f"bearing is {brg:.0f} deg = toward the SW (got {bb.min():.1f}-{bb.max():.1f})")
        check(not df["centroid_valid"].iat[8] and np.isnan(df["centroid_x_m"].iat[8]),
              "the faint uniform hour (0.05 mm/hr, under the 0.1 floor) gets no centroid")
        check(abs(df["basin_mean_mm_hr"].iat[8] - 0.05) < 1e-4 and abs(df["grid_mean_mm_hr"].iat[8] - 0.05) < 1e-4,
              "a uniform field gives the same basin mean and grid mean (area weights add to 1)")
        check(bool(df["filled_by_zero"].iat[9]) and not df["centroid_valid"].iat[9],
              "the manifest-flagged hour is marked filled and gets no centroid even though it contains rain")
        check(np.isnan(df["step_speed_kmh"].iat[10]) and np.isfinite(df["step_speed_kmh"].iat[11]),
              "no speed is computed across a gap (hour 10 blank), and speeds resume at hour 11")
        check(df["step_speed_kmh"].isna().iat[0], "the first hour has no step (no previous hour)")
        want_d = math.hypot(457000 - GAUGE_XY[0], 3749000 - GAUGE_XY[1]) / 1000
        check(abs(df["centroid_x_m"].iat[12] - 457000) < 1 and abs(df["centroid_y_m"].iat[12] - 3749000) < 1,
              "a uniform 10 x 10 km block has its centroid at the block centre (457000, 3749000)")
        check(abs(df["dist_to_gauge_km"].iat[12] - want_d) < 0.01, f"distance to the gauge is {want_d:.3f} km")
        check(abs(df["basin_mean_mm_hr"].iat[12] - 4.0) < 0.12 and df["basin_wet_frac"].iat[12] > 0.95,
              f"rain only over the test basin gives a basin mean near 4.0 (got {df['basin_mean_mm_hr'].iat[12]:.3f}) and wet fraction near 1")
        check(abs(df["spread_km"].iat[12] - math.sqrt(2 * (10 ** 2 - 1) / 12)) < 0.02,
              "spread of a uniform 10 x 10 km block matches the exact value")
        check(not bool(df.loc[[0, 1, 4, 7, 10, 11, 12], "edge_flag"].any()), "storms well inside the grid are not edge-flagged")
        # edge behaviour, checked directly on in-memory grids: a storm 2 km from the west edge vs one in the middle
        one = lambda xk: (5.0 * np.exp(-((X - xk) ** 2 + (Y - 3745000.0) ** 2) / (2 * sigma ** 2)))[None, :, :]
        ones = np.ones((49, 32))
        t_edge = compute_table(pd.DatetimeIndex([labels[0]]), one(442000.0), geo0, ones, set(), 0.1, 20, 0.02)
        t_mid = compute_table(pd.DatetimeIndex([labels[0]]), one(456000.0), geo0, ones, set(), 0.1, 20, 0.02)
        check(bool(t_edge["edge_flag"].iat[0]) and not bool(t_mid["edge_flag"].iat[0]),
              f"a storm 2 km from the west edge is edge-flagged ({t_edge['edge_rain_frac'].iat[0]:.1%} of its rain in the outer ring); one in the middle is not")
        check(abs(t_edge["centroid_x_m"].iat[0] - 442000.0) > 300 and abs(t_mid["centroid_x_m"].iat[0] - 456000.0) < 50,
              f"...and the flagged storm's centroid really is biased ({t_edge['centroid_x_m'].iat[0] - 442000.0:+.0f} m) while the central one is not")
        # time-axis labels on the overview colour bar: few, inside the window, never repeated, for short and long windows
        import matplotlib.dates as mdates
        lo = float(mdates.date2num(pd.Timestamp("2017-01-20 07:00").to_pydatetime()))
        for hrs in (3, 26, 168, 24 * 40):
            pos, lab = _time_ticks(lo, lo + hrs / 24)
            inside_ok = all(lo - 1e-9 <= p <= lo + hrs / 24 + 1e-9 for p in pos)
            check(2 <= len(pos) <= 5 and inside_ok and len(set(lab)) == len(lab),
                  f"time-axis labels for a {hrs}-hour window: {len(pos)} ticks, all inside the window, none repeated")
        info = txtp.read_text()
        check("hours with a centroid: 11 of 13" in info, "run_info.txt reports 11 of 13 hours with a centroid")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\nSELF-TEST " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


# ================================================================================================================
# command line
# ================================================================================================================
def build_parser():
    p = argparse.ArgumentParser(description="Rain movement over the Lower Sycamore basin from hourly MRMS grids "
                                            "(see the text at the top of this file).",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="make the CSV, overview figure and animation for one window")
    r.add_argument("--rain-dir", default=DEFAULT_RAIN_DIR, help=f"folder with radar_MMDDYYYYHH.asc (default {DEFAULT_RAIN_DIR})")
    r.add_argument("--prefix", default="radar_", help="file-name start before MMDDYYYYHH (default radar_)")
    r.add_argument("--boundary", default=DEFAULT_BOUNDARY, help=f"basin outline shapefile in EPSG:26912 (default {DEFAULT_BOUNDARY})")
    r.add_argument("--manifest", default=None, help="manifest.csv from the radar build (default: one folder above --rain-dir)")
    r.add_argument("--start", default=DEFAULT_START, help=f"start of the first rain hour, MST, whole hour (default {DEFAULT_START})")
    r.add_argument("--end", default=DEFAULT_END, help=f"end of the last rain hour, MST, whole hour (default {DEFAULT_END})")
    r.add_argument("--tag", default=None, help="start of the output file names (default: the two dates)")
    r.add_argument("--out-dir", default="storm_movement_out", help="folder for the outputs (default storm_movement_out)")
    r.add_argument("--cell-min", type=float, default=0.1, help="a cell counts as wet from this mm/hr (default 0.1)")
    r.add_argument("--min-wet-cells", type=int, default=20, help="wet cells needed before an hour gets a centroid (default 20)")
    r.add_argument("--edge-frac", type=float, default=0.02, help="mark an hour edge_flag when more than this share of the wet rain is in the outermost ring of grid cells (default 0.02)")
    r.add_argument("--trail-hours", type=int, default=6, help="hours of centroid trail drawn in the animation (default 6)")
    r.add_argument("--vmax", type=float, default=8.0, help="top of the animation colour scale in mm/hr (default 8)")
    r.add_argument("--fps", type=int, default=6, help="animation frames per second (default 6)")
    r.add_argument("--no-gif", action="store_true", help="skip the animation (the CSV and PNG are quick)")
    r.add_argument("--expect-area-km2", type=float, default=424.6, help="warn if the outline's area on the grid differs by >3%% (0 = off)")
    r.add_argument("--dry-run", action="store_true", help="check the files and settings, write nothing")
    sub.add_parser("self-test", help="check the script on made-up storms")
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    return cmd_run(a) if a.cmd == "run" else cmd_self_test(a)


if __name__ == "__main__":
    sys.exit(main())
