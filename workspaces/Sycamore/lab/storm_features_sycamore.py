#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
storm_features_sycamore.py   (v1.0, 2026-11-03)
===============================================
A feature table for every storm on the list, from the hourly radar grids ALONE (no model output, no stream flow, so it is safe to
run on WY2020 while WY2020 stays sealed).

WHY
  The by-storm figures showed that storm TOTAL depth does not tell us how much the stream (or the model) will give: Dec 21-25
  (48 mm over five days) gave no flow, Jan 19-21 (56 mm, 41 mm of it in one day) gave the flood. This script measures the other
  things about each storm that the totals hide: how hard it rained, where, how the rain moved, and how wet the basin already was.

WHAT IT MAKES  (all in --out-dir, default storm_features_out/)
  storm_features.csv            one row per storm (WY2017 and WY2020 together), every number listed below
  storm_features_report.txt     the same numbers as readable tables, plus a rule check (does the 0.5 mm / 20 mm rule find the storms on the list?)
  hourly_basin_wy2017.csv       hour-by-hour basin-mean rain and intensity numbers for the whole grid record (also useful on its own)
  hourly_basin_wy2020.csv       same for WY2020 (rain only)
  fig_storm_intensity_wy2017.png, fig_storm_intensity_wy2020.png     one column per storm: hourly basin rain (bars) and the share of the
                                basin raining harder than the soil can take (line)

THE NUMBERS  (one row per storm; the storm = whole local days from the first to the last rain day)
  Intensity (item 1)
    basin_mm               rain over the basin, the same area-weighted mean tRIBS is given (mm)
    grid_mm                mean over the whole 32 x 49 km radar grid (mm; matches the earlier tables)
    peak_1h, peak_3h, peak_6h, peak_24h   largest basin-mean rain in any 1, 3, 6, 24 hours (mm in that span)
    wettest_day_mm         the wettest local day (mm)
    hours_wet              hours with basin-mean rain of at least 0.1 mm
    hours_ge_<T>, depth_ge_<T>   hours when the BASIN-MEAN rain was at least T mm/h (T = 2.2 and 4 by default), and the rain that fell in them
    excess_<T>_mm          "rain the soil could not take": for every cell and hour, the part of the rain above T mm/h, averaged over the basin
                           and added up (mm). This is the better intensity number, because the soil sees each cell's rain, not the basin mean.
    excess_<T>_share       excess_<T>_mm divided by basin_mm
    peak_cell_mm_hr        the largest rain rate in any basin cell in any hour
    top6h_share            the share of the storm's rain that fell in its six wettest hours
    t50_frac               how far through the storm (0 = start, 1 = end) half of the rain had fallen (small = front-loaded)
  Where the rain fell and how it moved (item 2)
    basin_over_grid        basin_mm / grid_mm
    spatial_cv             how uneven the rain was across the basin (standard deviation / mean of the cell totals, area weighted)
    frac_basin_ge_20mm, frac_basin_ge_40mm   share of the basin area whose storm total reached 20 / 40 mm
    rain_dist_km, area_dist_km, dist_skew_km   rain-weighted mean distance of the basin's rain from the stream gauge, the plain
                           area-weighted mean distance of the basin, and the difference (negative = more rain than usual near the outlet)
    centroid_bearing, centroid_dist_km   where the storm-total rain "balances" relative to the gauge (bearing in degrees clockwise from north)
    arrival_dir, arrival_speed_kmh, arrival_r2, arrival_spread_h   which way the rain MOVED across the basin, from the time each cell had
                           half its rain (a straight-front fit; the speed is only given when the fit is decent, R2 of 0.3 or more);
                           arrival_spread_h = time between the 10th and 90th percentile cell (how long the storm took to cross the basin)
  Starting conditions that come from rain alone (item 3, rain side)
    prev7_mm, prev30_mm, season_mm    basin rain in the 7 and 30 days before the first rain day, and since Oct 1
    dry_days_before        number of days in a row, just before the storm, with less than 2 mm of basin rain (the dry spell)
    prev_event, prev_event_mm, prev_gap_days    the previous rain event of 10 mm or more (grid mean): its dates, size and the gap
    lookback_incomplete    True when the radar record does not reach back far enough for prev30 (spin-up storms)
  Data quality (item 7, radar side)
    radar_filled_hours, radar_missing_hours   hours in the storm (and the 24 hours before) when MRMS had no data (zero-filled) or the file is missing
  Rain uncertainty (item 8; only when --ratios is given)
    ratio_3g, ratio_2g     radar / gauge rain for the storm, from w5_rain_by_storm.csv
    basin_mm_if_gauges_2g  basin_mm divided by ratio_2g: what the basin would have received if the two good gauges were right and the
                           radar were off by the same factor everywhere (a bracket, not a correction)

THE STORMS
  The built-in list is the PROPOSED list (not yet accepted): the 12 winter storms (5 in WY2017, 7 in WY2020), plus the sub-20 mm
  near-misses and the Sept 2019 spin-up storms as context rows (column kind = listed / near-miss / spin-up). Use --storms-file to
  give your own list (columns: year,name,first_day,last_day,kind).

COMMANDS  (run from lab/; nothing is changed, nothing is run except this script)
  python storm_features_sycamore.py self-test                  # checks the script on made-up storms with known answers
  python storm_features_sycamore.py run --dry-run              # lists the folders and files it would read, writes nothing
  python storm_features_sycamore.py run --ratios w5_rain_by_storm.csv
  python storm_features_sycamore.py run --years wy2020         # one year only

WORDS
  hourly grid     one radar file = the rain that fell in ONE hour, mm per cell; file names are radar_MMDDYYYYHH.asc, HH = the hour that ENDS
  area weighting  the basin mean counts each 1 km cell by the share of the cell inside the basin outline
  centroid        the point where the rain would balance if it were weights on a board
  front           the leading edge of a band of rain as it crosses the land
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

VERSION = 1
HOUR = pd.Timedelta(hours=1)
DAY = pd.Timedelta(days=1)
GAUGE_XY = (449721.0, 3728390.0)               # USGS 09510200, EPSG:26912
DEFAULT_BOUNDARY = "data/preprocessing/Sycamore_boundary.shp"
FILLED_STATUS = ("filled_zero_missing_hour",)
INK, INK_2, GRID_C, RAIN_LIGHT, RAIN_DARK, RED = "#0b0b0b", "#52514e", "#d9d8d3", "#bcd6f6", "#2a78d6", "#c0392b"

YEARS = {
    "WY2017": dict(start="2016-10-01 00:00", end="2017-10-01 00:00", rain_dir="~/mrms_wk_wy2017/rain", wy_start="2016-10-01"),
    "WY2020": dict(start="2019-09-01 00:00", end="2020-06-01 00:00", rain_dir="~/mrms_wk_wy2020/rain", wy_start="2019-10-01"),
}
# (name, first rain day, last rain day, kind)
STORMS = {
    "WY2017": [("Dec 15-17", "2016-12-15", "2016-12-17", "near-miss"),
               ("Dec 21-25", "2016-12-21", "2016-12-25", "listed"),
               ("Dec 30-Jan 2", "2016-12-30", "2017-01-02", "listed"),
               ("Jan 19-21", "2017-01-19", "2017-01-21", "listed"),
               ("Jan 23-24", "2017-01-23", "2017-01-24", "near-miss"),
               ("Feb 18-19", "2017-02-18", "2017-02-19", "listed"),
               ("Feb 27-28", "2017-02-27", "2017-02-28", "listed")],
    "WY2020": [("Sep 10", "2019-09-10", "2019-09-10", "spin-up"),
               ("Sep 23", "2019-09-23", "2019-09-23", "spin-up"),
               ("Nov 19-22", "2019-11-19", "2019-11-22", "listed"),
               ("Nov 28-30", "2019-11-28", "2019-11-30", "listed"),
               ("Dec 7-8", "2019-12-07", "2019-12-08", "near-miss"),
               ("Dec 23-28", "2019-12-23", "2019-12-28", "listed"),
               ("Feb 10-11", "2020-02-10", "2020-02-11", "listed"),
               ("Feb 22-23", "2020-02-22", "2020-02-23", "listed"),
               ("Mar 10-13", "2020-03-10", "2020-03-13", "listed"),
               ("Mar 18-19", "2020-03-18", "2020-03-19", "listed")],
}


def die(msg: str):
    raise SystemExit("ERROR: " + msg)


def log(msg: str = "") -> None:
    print(msg, flush=True)


# ================================================================================================================
# reading the grids and the basin outline (same conventions as storm_movement_sycamore.py)
# ================================================================================================================
def file_for_label(rain_dir: Path, label: pd.Timestamp) -> Path:
    return Path(rain_dir) / f"radar_{label:%m%d%Y%H}.asc"


def read_ascii(path: Path):
    """ESRI ASCII grid -> (header dict, 2-D float array with the NORTH row first). NODATA -> NaN."""
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
    return dict(nx=nx, ny=ny, cs=cs, L=L, R=R, B=B, T=T, X=X, Y=Y, sig=(nx, ny, round(L, 3), round(B, 3), cs))


def load_boundary_rings(path):
    """Basin outline as a list of (N, 2) arrays in EPSG:26912 (geopandas or pyshp)."""
    path = Path(path).expanduser()
    if not path.exists():
        die(f"basin outline not found: {path}  (run from lab/, or pass --boundary path/to/Sycamore_boundary.shp)")
    prj = path.with_suffix(".prj")
    if prj.exists():
        w = prj.read_text()
        if "UTM zone 12N" not in w and "UTM_Zone_12N" not in w and "26912" not in w:
            die(f"{prj.name} is not NAD83 / UTM zone 12N; reproject the outline first.")
    rings = []
    try:
        import geopandas as gpd
        for geom in gpd.read_file(path).geometry:
            if geom is None or geom.is_empty or geom.geom_type not in ("Polygon", "MultiPolygon"):
                continue
            for p in (list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]):
                rings.append(np.asarray(p.exterior.coords)[:, :2])
                for h in p.interiors:
                    rings.append(np.asarray(h.coords)[:, :2])
    except ImportError:
        try:
            import shapefile
        except ImportError:
            die("cannot read the basin outline: install geopandas or pyshp")
        sf = shapefile.Reader(str(path))
        for shp in sf.shapes():
            pts = np.array(shp.points)
            parts = list(shp.parts) + [len(pts)]
            for i in range(len(parts) - 1):
                rings.append(pts[parts[i]:parts[i + 1]])
    if not rings:
        die(f"no polygon found in {path.name}")
    return rings


def grid_cover(rings, geo: dict, ss: int = 50):
    """Share (0-1) of each cell inside the basin outline (outline drawn on a grid ss times finer, even-odd rule)."""
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


def load_year(rain_dir: Path, start: str, end: str):
    """Read every hourly grid between start and end (whole hours). Returns (labels, stack mm, geo, missing flags).
    A missing file becomes an all-zero hour and is flagged; a grid with a different header stops the run."""
    s, e = pd.Timestamp(start), pd.Timestamp(end)
    n = int(round((e - s) / HOUR))
    labels = pd.DatetimeIndex([s + HOUR * (i + 1) for i in range(n)])
    files = [file_for_label(rain_dir, lab) for lab in labels]
    present = [f.exists() for f in files]
    if not any(present):
        return labels, None, None, np.ones(n, dtype=bool)
    first = files[present.index(True)]
    hdr0, _ = read_ascii(first)
    geo = grid_geometry(hdr0)
    stack = np.zeros((n, geo["ny"], geo["nx"]), dtype="float32")
    for k, f in enumerate(files):
        if not present[k]:
            continue
        h, a = read_ascii(f)
        if grid_geometry(h)["sig"] != geo["sig"]:
            die(f"{f.name} has a different grid header from {first.name}")
        stack[k] = np.clip(np.nan_to_num(a, nan=0.0), 0.0, None)
        if k % 2000 == 0 and k:
            log(f"    read {k} of {n} hourly grids ...")
    return labels, stack, geo, ~np.array(present)


# ================================================================================================================
# hourly numbers
# ================================================================================================================
def tag_of(T: float) -> str:
    return f"{T:g}".replace(".", "p")


def hourly_table(labels, rain, geo, cover, filled_keys, missing, thresholds, cell_min=0.1):
    """One row per hour: basin-mean and grid-mean rain, wettest basin cell, wet share, and for each threshold T the basin share
    raining at least T and the basin-mean 'excess' (rain above T). `rain` is (n, ny, nx) mm per hour."""
    wb = cover / cover.sum()
    inside = cover >= 0.5
    r64 = rain.astype("float64")
    df = pd.DataFrame({
        "label_end_mst": labels, "interval_start_mst": labels - HOUR,
        "grid_mean_mm": r64.mean(axis=(1, 2)),
        "basin_mean_mm": np.tensordot(r64, wb, axes=([1, 2], [0, 1])),
        "basin_max_mm": r64[:, inside].max(axis=1) if inside.any() else np.nan,
        "basin_wet_frac": np.tensordot((r64 >= cell_min).astype("float64"), wb, axes=([1, 2], [0, 1])),
    })
    for T in thresholds:
        df[f"frac_basin_ge_{tag_of(T)}"] = np.tensordot((r64 >= T).astype("float64"), wb, axes=([1, 2], [0, 1]))
        df[f"basin_excess_{tag_of(T)}_mm"] = np.tensordot(np.maximum(r64 - T, 0.0), wb, axes=([1, 2], [0, 1]))
    df["filled_zero"] = [f"{lab:%Y-%m-%d %H:%M}" in filled_keys for lab in labels]
    df["missing_file"] = np.asarray(missing, dtype=bool)
    return df


def daily_sum(df: pd.DataFrame, col: str) -> pd.Series:
    s = df.set_index("interval_start_mst")[col]
    return s.groupby(s.index.floor("D")).sum()


def find_events(daily: pd.Series, min_day: float = 0.5, min_total: float = 20.0):
    """Runs of consecutive days with at least min_day mm; keep runs whose total is at least min_total. -> [(first, last, total)]"""
    out, run = [], []
    days = pd.date_range(daily.index.min(), daily.index.max(), freq="D")
    d = daily.reindex(days).fillna(0.0)
    for day, v in d.items():
        if v >= min_day:
            run.append((day, v))
        elif run:
            if sum(x for _, x in run) >= min_total:
                out.append((run[0][0], run[-1][0], float(sum(x for _, x in run))))
            run = []
    if run and sum(x for _, x in run) >= min_total:
        out.append((run[0][0], run[-1][0], float(sum(x for _, x in run))))
    return out


# ================================================================================================================
# per-storm numbers
# ================================================================================================================
def arrival_fit(rain_win, geo, cover, min_total_mm=2.0):
    """Which way did the rain move across the basin? Half-way time of every basin cell, then a straight-front fit
    t = c + a*east + b*north (weighted by the cell's rain). Direction = where the rain moved TO (degrees clockwise from north)."""
    nan = dict(arrival_dir=np.nan, arrival_speed_kmh=np.nan, arrival_r2=np.nan, arrival_spread_h=np.nan, arrival_cells=0)
    tot = rain_win.sum(axis=0)
    sel = (cover >= 0.5) & (tot >= min_total_mm)
    m = int(sel.sum())
    if m < 30:
        return {**nan, "arrival_cells": m}
    r = rain_win[:, sel].astype("float64")
    cum = np.cumsum(r, axis=0)
    half = 0.5 * tot[sel].astype("float64")
    ge = cum >= half[None, :] - 1e-9
    idx = ge.argmax(axis=0)
    cols = np.arange(m)
    prev = np.where(idx > 0, cum[np.maximum(idx - 1, 0), cols], 0.0)
    rk = r[idx, cols]
    frac = np.where(rk > 0, (half - prev) / np.where(rk > 0, rk, 1.0), 0.5)
    t = idx + np.clip(frac, 0.0, 1.0)                       # hours since the window start (hour k covers k to k+1)
    x = geo["X"][sel] / 1000.0
    y = geo["Y"][sel] / 1000.0
    w = tot[sel].astype("float64")
    A = np.column_stack([np.ones(m), x - x.mean(), y - y.mean()])
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(A * sw[:, None], t * sw, rcond=None)
    pred = A @ coef
    tbar = np.average(t, weights=w)
    ss_tot = float(np.sum(w * (t - tbar) ** 2))
    r2 = 1.0 - float(np.sum(w * (t - pred) ** 2)) / ss_tot if ss_tot > 1e-9 else np.nan
    gx, gy = float(coef[1]), float(coef[2])                 # hours per km toward the east and the north
    g = math.hypot(gx, gy)
    bearing = (math.degrees(math.atan2(gx, gy)) + 360.0) % 360.0 if g > 1e-9 else np.nan
    speed = 1.0 / g if g > 1e-9 else np.nan
    p10, p90 = np.percentile(t, [10, 90])
    if not (np.isfinite(r2) and r2 >= 0.3):
        speed = np.nan                                      # a weak fit: do not quote a speed
    return dict(arrival_dir=bearing, arrival_speed_kmh=speed, arrival_r2=r2, arrival_spread_h=float(p90 - p10), arrival_cells=m)


def storm_features(hours: pd.DataFrame, rain, geo, cover, name, first, last, kind, thresholds, ctx: dict) -> dict:
    """All numbers for one storm. `hours` is the hourly table of the whole record, `rain` the matching (n, ny, nx) stack (mm/h)."""
    a, b = pd.Timestamp(first), pd.Timestamp(last) + DAY
    sel = ((hours["interval_start_mst"] >= a) & (hours["interval_start_mst"] < b)).to_numpy()
    h = hours[sel]
    row = dict(storm=name, kind=kind, first_day=a, last_day=pd.Timestamp(last), window_hours=int(sel.sum()))
    if not sel.any():
        row["note"] = "storm is outside the radar record"
        return row
    bm = h["basin_mean_mm"].to_numpy()
    total = float(bm.sum())
    row["basin_mm"] = total
    row["grid_mm"] = float(h["grid_mean_mm"].sum())
    roll = pd.Series(bm).rolling
    for k in (1, 3, 6, 24):
        row[f"peak_{k}h"] = float(roll(min(k, len(bm))).sum().max())
    daily = pd.Series(bm, index=h["interval_start_mst"].dt.floor("D").to_numpy()).groupby(level=0).sum()
    row["wettest_day_mm"] = float(daily.max())
    row["hours_wet"] = int((bm >= 0.1).sum())
    for T in thresholds:
        tg = tag_of(T)
        hit = bm >= T
        row[f"hours_ge_{tg}"] = int(hit.sum())
        row[f"depth_ge_{tg}"] = float(bm[hit].sum())
        ex = float(h[f"basin_excess_{tg}_mm"].sum())
        row[f"excess_{tg}_mm"] = ex
        row[f"excess_{tg}_share"] = ex / total if total > 0 else np.nan
        row[f"peak_frac_basin_ge_{tg}"] = float(h[f"frac_basin_ge_{tg}"].max())
    row["peak_cell_mm_hr"] = float(h["basin_max_mm"].max())
    row["top6h_share"] = float(np.sort(bm)[::-1][:6].sum() / total) if total > 0 else np.nan
    if total > 0:
        cum = np.cumsum(bm) / total
        row["t50_frac"] = float((int(np.argmax(cum >= 0.5)) + 0.5) / len(bm))
    else:
        row["t50_frac"] = np.nan
    # ---- where the rain fell ----
    win = rain[sel]
    tot = win.sum(axis=0).astype("float64")
    wb = cover / cover.sum()
    mean_cell = float((wb * tot).sum())
    row["basin_over_grid"] = total / row["grid_mm"] if row["grid_mm"] > 0 else np.nan
    if mean_cell > 0:
        var = float((wb * (tot - mean_cell) ** 2).sum())
        row["spatial_cv"] = math.sqrt(var) / mean_cell
    else:
        row["spatial_cv"] = np.nan
    row["frac_basin_ge_20mm"] = float((wb * (tot >= 20.0)).sum())
    row["frac_basin_ge_40mm"] = float((wb * (tot >= 40.0)).sum())
    gx, gy = GAUGE_XY
    dist = np.hypot(geo["X"] - gx, geo["Y"] - gy) / 1000.0
    area_d = float((wb * dist).sum())
    rain_d = float((wb * tot * dist).sum() / (wb * tot).sum()) if (wb * tot).sum() > 0 else np.nan
    row["area_dist_km"], row["rain_dist_km"], row["dist_skew_km"] = area_d, rain_d, rain_d - area_d
    if tot.sum() > 0:
        cx = float((tot * geo["X"]).sum() / tot.sum())
        cy = float((tot * geo["Y"]).sum() / tot.sum())
        row["centroid_dist_km"] = math.hypot(cx - gx, cy - gy) / 1000.0
        row["centroid_bearing"] = (math.degrees(math.atan2(cx - gx, cy - gy)) + 360.0) % 360.0
    else:
        row["centroid_dist_km"] = row["centroid_bearing"] = np.nan
    row.update(arrival_fit(win, geo, cover))
    # ---- starting conditions from rain alone ----
    hb = hours["basin_mean_mm"].to_numpy()
    st = hours["interval_start_mst"]
    rec0 = st.min()
    for k in (7, 30):
        lo = a - k * DAY
        m = ((st >= lo) & (st < a)).to_numpy()
        row[f"prev{k}_mm"] = float(hb[m].sum()) if lo >= rec0 else np.nan
    row["lookback_incomplete"] = bool(a - 30 * DAY < rec0)
    wy0 = pd.Timestamp(ctx["wy_start"])
    row["season_mm"] = float(hb[((st >= wy0) & (st < a)).to_numpy()].sum()) if (a >= wy0 and rec0 <= wy0) else np.nan
    dly = ctx["daily_basin"]
    wet_days = dly[(dly.index < a) & (dly >= ctx["wet_day_mm"])]
    row["dry_days_before"] = float((a - (wet_days.index.max() + DAY)) / DAY) if len(wet_days) else np.nan
    prev = [e for e in ctx["events"] if e[1] < a]
    if prev:
        e = prev[-1]
        row["prev_event"] = (f"{e[0]:%b} {e[0].day}" if e[0] == e[1] else
                             f"{e[0]:%b} {e[0].day}-{e[1].day}" if e[0].month == e[1].month else
                             f"{e[0]:%b} {e[0].day}-{e[1]:%b} {e[1].day}")
        row["prev_event_mm"] = e[2]
        row["prev_gap_days"] = float((a - (e[1] + DAY)) / DAY)
    else:
        row["prev_event"], row["prev_event_mm"], row["prev_gap_days"] = "", np.nan, np.nan
    # ---- quality ----
    pre = ((st >= a - DAY) & (st < b)).to_numpy()
    row["radar_filled_hours"] = int(hours["filled_zero"].to_numpy()[pre].sum())
    row["radar_missing_hours"] = int(hours["missing_file"].to_numpy()[pre].sum())
    return row


# ================================================================================================================
# figure and report
# ================================================================================================================
def make_figure(out: Path, year: str, rows: list, hours: pd.DataFrame, thresholds: list) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    rows = [r for r in rows if "basin_mm" in r]
    if not rows:
        return "no storms in the record, no figure"
    n = len(rows)
    T0 = thresholds[0]
    tg = tag_of(T0)
    fig, axes = plt.subplots(2, n, figsize=(max(9.0, 2.3 * n + 1.0), 6.4), squeeze=False,
                             gridspec_kw={"height_ratios": [2.3, 1.3], "hspace": 0.32, "wspace": 0.12})
    ymax = 0.5
    wins = []
    for r in rows:
        a, b = r["first_day"] - 12 * HOUR, r["last_day"] + DAY + 12 * HOUR
        m = (hours["interval_start_mst"] >= a) & (hours["interval_start_mst"] < b)
        wins.append((a, b, hours[m]))
        ymax = max(ymax, float(hours[m]["basin_mean_mm"].max()))
    for i, (r, (a, b, h)) in enumerate(zip(rows, wins)):
        ax, ax2 = axes[0, i], axes[1, i]
        x = mdates.date2num(h["interval_start_mst"].to_numpy("datetime64[ns]"))
        bm = h["basin_mean_mm"].to_numpy()
        ax.bar(x, bm, width=1 / 24, align="edge", color=[RAIN_DARK if v >= T0 else RAIN_LIGHT for v in bm], lw=0)
        for T in thresholds:
            ax.axhline(T, color=RED, lw=0.8, ls=(0, (4, 3)))
        ax.set_ylim(0, ymax * 1.12)
        ax.set_title(f"{r['storm']}\n{r['basin_mm']:.1f} mm, peak {r['peak_1h']:.1f} mm/h", fontsize=9.5, color=INK)
        ax2.plot(x + 0.5 / 24, h[f"frac_basin_ge_{tg}"].to_numpy(), color=INK, lw=1.0, drawstyle="steps-mid")
        ax2.set_ylim(0, 1.0)
        for a_ in (ax, ax2):
            a_.set_xlim(mdates.date2num(a.to_datetime64()), mdates.date2num(b.to_datetime64()))
            a_.xaxis.set_major_locator(mdates.DayLocator())
            a_.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
            a_.tick_params(labelsize=7.5)
            for s_ in ("top", "right"):
                a_.spines[s_].set_visible(False)
            a_.grid(True, axis="y", alpha=0.25)
            a_.set_axisbelow(True)
        ax.tick_params(labelbottom=False)
        ax2.set_xlabel(f"{r['first_day']:%b %Y}", fontsize=8, color=INK_2)
        if i == 0:
            ax.set_ylabel("Basin-mean rain (mm per hour)")
            ax2.set_ylabel(f"Share of basin raining\nat least {T0:g} mm/h")
        else:
            ax.tick_params(labelleft=False)
            ax2.tick_params(labelleft=False)
    leg = [Patch(color=RAIN_DARK, label=f"hour with basin-mean rain at least {T0:g} mm/h"), Patch(color=RAIN_LIGHT, label="hours below that rate"),
           plt.Line2D([], [], color=RED, ls=(0, (4, 3)), label="thresholds: " + ", ".join(f"{T:g}" for T in thresholds) + " mm/h")]
    fig.suptitle(f"{year}: how hard it rained, storm by storm (radar, area-weighted over the basin)", x=0.01, ha="left", fontsize=12)
    fig.legend(handles=leg, loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.07, right=0.99, top=0.88, bottom=0.14)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return f"wrote {out}"


def report_tables(df: pd.DataFrame, thresholds: list) -> list:
    L = []
    tg = [tag_of(T) for T in thresholds]

    def show(title, cols, note=""):
        L.append("")
        L.append(title)
        if note:
            L.append("  " + note)
        use = ["year", "storm", "kind"] + [c for c in cols if c in df.columns]
        L.append(df[use].to_string(index=False, float_format=lambda x: f"{x:.2f}", na_rep="-"))

    show("1. INTENSITY (mm; rain per hour in mm/h)",
         ["basin_mm", "grid_mm", "peak_1h", "peak_3h", "peak_6h", "peak_24h", "wettest_day_mm", "hours_wet", "peak_cell_mm_hr", "top6h_share", "t50_frac"],
         "peak_Nh = the most rain in any N hours (mm).  top6h_share = share of the storm's rain in its six wettest hours.  t50_frac = how far through the storm half the rain had fallen.")
    show("1b. RAIN ABOVE WHAT THE SOIL CAN TAKE",
         sum([[f"hours_ge_{t}", f"depth_ge_{t}", f"excess_{t}_mm", f"excess_{t}_share", f"peak_frac_basin_ge_{t}"] for t in tg], []),
         "hours_ge/depth_ge use the basin mean; excess_ adds up, cell by cell, the part of the rain above the threshold (the better number).")
    show("2. WHERE THE RAIN FELL AND HOW IT MOVED",
         ["basin_over_grid", "spatial_cv", "frac_basin_ge_20mm", "frac_basin_ge_40mm", "rain_dist_km", "area_dist_km", "dist_skew_km",
          "centroid_dist_km", "centroid_bearing", "arrival_dir", "arrival_speed_kmh", "arrival_r2", "arrival_spread_h"],
         "dist_skew_km below zero = more rain than usual near the outlet.  Bearings are degrees clockwise from north (0 N, 90 E, 180 S, 270 W).  "
         "arrival_dir = where the rain moved TO; a speed is given only when R2 is 0.3 or more.")
    show("3. STARTING CONDITIONS FROM RAIN ALONE (basin rain, mm)",
         ["prev7_mm", "prev30_mm", "season_mm", "dry_days_before", "prev_event", "prev_event_mm", "prev_gap_days", "lookback_incomplete"],
         "season_mm = rain since Oct 1.  prev_event = the last rain event of 10 mm or more (grid mean) before the storm.  dry_days_before: days under 2 mm of basin rain in a row just before the storm.")
    show("7. RADAR DATA QUALITY", ["window_hours", "radar_filled_hours", "radar_missing_hours"],
         "hours inside the storm and the 24 hours before it when MRMS had no data (zero-filled) or the file is missing.")
    if "ratio_2g" in df.columns:
        show("8. RAIN UNCERTAINTY (radar / gauge, from the ratios file)",
             ["basin_mm", "ratio_3g", "ratio_2g", "basin_mm_if_gauges_2g", "grid_mm_check"],
             "basin_mm_if_gauges_2g = basin_mm / ratio_2g, a bracket only.  grid_mm_check = grid_mm minus the total in the ratios file (should be about 0).")
    return L


# ================================================================================================================
# the run command
# ================================================================================================================
def parse_threshold_list(text: str) -> list:
    try:
        v = sorted({float(x) for x in text.split(",") if x.strip()})
    except ValueError:
        die(f"--thresholds must be numbers separated by commas, got {text!r}")
    if not v or any(x <= 0 for x in v):
        die("--thresholds needs at least one positive number")
    return v


def load_storms(path, years):
    if not path:
        return {y: list(STORMS[y]) for y in years}
    d = pd.read_csv(path)
    need = {"year", "name", "first_day", "last_day"}
    if not need <= set(d.columns):
        die(f"--storms-file needs the columns {sorted(need)} (and optionally kind)")
    out = {}
    for y in years:
        sub = d[d["year"].astype(str).str.upper() == y]
        out[y] = [(r["name"], str(r["first_day"]), str(r["last_day"]), str(r["kind"]) if "kind" in d.columns else "listed") for _, r in sub.iterrows()]
    return out


def run_year(year, spec, storms, a, thresholds, boundary_rings):
    rain_dir = Path(spec["rain_dir"]).expanduser()
    log(f"\n{year}: rain folder {rain_dir}")
    labels, stack, geo, missing = load_year(rain_dir, spec["start"], spec["end"])
    n = len(labels)
    log(f"  {n - int(missing.sum())} of {n} hourly grids present ({spec['start']} to {spec['end']})")
    if stack is None:
        log("  none found; skipping this year (check --rain-dir-" + year.lower() + ")")
        return None
    cover = grid_cover(boundary_rings, geo)
    area = float(cover.sum()) * geo["cs"] ** 2 / 1e6
    log(f"  grid {geo['nx']} x {geo['ny']} cells; basin outline covers {area:.1f} km2 of it")
    if a.expect_area_km2 and abs(area - a.expect_area_km2) > 0.03 * a.expect_area_km2:
        log(f"  WARNING: expected about {a.expect_area_km2:g} km2; is --boundary the right outline?")
    filled = read_filled_keys(Path(a.manifest).expanduser() if a.manifest else rain_dir.parent / "manifest.csv")
    log("  manifest: " + ("not found, zero-filled hours cannot be flagged" if filled is None else f"{len(filled)} zero-filled hours in the record"))
    hours = hourly_table(labels, stack, geo, cover, filled or set(), missing, thresholds, a.cell_min)
    ctx = dict(wy_start=spec["wy_start"], wet_day_mm=a.wet_day_mm, daily_basin=daily_sum(hours, "basin_mean_mm"),
               events=find_events(daily_sum(hours, "grid_mean_mm"), 0.5, a.prev_event_mm))
    rows = [{"year": year, **storm_features(hours, stack, geo, cover, nm, f, l, k, thresholds, ctx)} for nm, f, l, k in storms]
    ev20 = find_events(daily_sum(hours, "grid_mean_mm"), 0.5, 20.0)
    return dict(hours=hours, rows=rows, events20=ev20, area=area)


def cmd_run(a) -> int:
    years = [y.strip().upper() for y in a.years.split(",") if y.strip()]
    bad = [y for y in years if y not in YEARS]
    if bad:
        die(f"unknown year {bad}; use wy2017 and/or wy2020")
    thresholds = parse_threshold_list(a.thresholds)
    storms = load_storms(a.storms_file, years)
    specs = {}
    for y in years:
        specs[y] = dict(YEARS[y])
        specs[y]["rain_dir"] = getattr(a, f"rain_dir_{y.lower()}") or YEARS[y]["rain_dir"]
        w = getattr(a, f"window_{y.lower()}")
        if w:
            parts = [p.strip() for p in w.split(",")]
            if len(parts) != 2:
                die(f"--window-{y.lower()} must be 'start,end' like '2016-10-01 00:00,2017-10-01 00:00'")
            specs[y]["start"], specs[y]["end"] = parts
    out_dir = Path(a.out_dir).expanduser()
    log("storm_features_sycamore: " + ("DRY RUN (nothing is read in full, nothing is written)" if a.dry_run else "run"))
    log(f"  thresholds {thresholds} mm/h, cell-min {a.cell_min:g}, wet day {a.wet_day_mm:g} mm, previous event {a.prev_event_mm:g} mm")
    if a.dry_run:
        for y in years:
            sp = specs[y]
            rd = Path(sp["rain_dir"]).expanduser()
            s, e = pd.Timestamp(sp["start"]), pd.Timestamp(sp["end"])
            labs = [s + HOUR * (i + 1) for i in range(int(round((e - s) / HOUR)))]
            have = sum(file_for_label(rd, lab).exists() for lab in labs)
            log(f"  {y}: {rd}  {have} of {len(labs)} files found;  storms: {len(storms[y])}")
        bp = Path(a.boundary).expanduser()
        log(f"  basin outline: {bp}  ({'found' if bp.exists() else 'NOT FOUND'})")
        log(f"  would write to {out_dir}/")
        return 0
    rings = load_boundary_rings(a.boundary)
    results = {}
    for y in years:
        r = run_year(y, specs[y], storms[y], a, thresholds, rings)
        if r is not None:
            results[y] = r
    if not results:
        die("no radar grids were found for any year")
    out_dir.mkdir(parents=True, exist_ok=True)
    allrows = [row for y in results for row in results[y]["rows"]]
    df = pd.DataFrame(allrows)
    if a.ratios:
        rp = Path(a.ratios).expanduser()
        if rp.exists():
            rt = pd.read_csv(rp)
            keep = [c for c in ("year", "storm", "ratio_3g", "ratio_2g", "grid_total_mm") if c in rt.columns]
            df = df.merge(rt[keep], on=["year", "storm"], how="left")
            if "ratio_2g" in df.columns:
                df["basin_mm_if_gauges_2g"] = df["basin_mm"] / df["ratio_2g"]
            if "grid_total_mm" in df.columns:
                df["grid_mm_check"] = df["grid_mm"] - df["grid_total_mm"]
                df = df.drop(columns=["grid_total_mm"])
        else:
            log(f"  note: --ratios file not found: {rp}")
    csv_path = out_dir / "storm_features.csv"
    df.to_csv(csv_path, index=False)
    for y in results:
        results[y]["hours"].to_csv(out_dir / f"hourly_basin_{y.lower()}.csv", index=False, float_format="%.5f")
    lines = ["storm_features_sycamore v%d: rain-only features of each storm (radar, area-weighted over the basin)" % VERSION,
             f"thresholds {thresholds} mm/h; the storm = whole local days from the first to the last rain day"]
    lines += report_tables(df, thresholds)
    lines += ["", "RULE CHECK: events found by the rule (consecutive days of at least 0.5 mm grid-mean rain, total 20 mm or more)"]
    for y in results:
        listed = {(pd.Timestamp(f), pd.Timestamp(l)) for nm, f, l, k in storms[y] if k == "listed"}
        lines.append(f"  {y} (basin outline {results[y]['area']:.1f} km2):")
        for e in results[y]["events20"]:
            tag = "on the list" if (e[0], e[1]) in listed else ("NOT on the list" + ("  (summer: left out by decision)" if e[0].month in (6, 7, 8, 9) else ""))
            lines.append(f"    {e[0]:%b %d %Y} to {e[1]:%b %d %Y}   {e[2]:6.1f} mm   {tag}")
        have = {(e[0], e[1]) for e in results[y]["events20"]}
        for f, l in sorted(listed - have):
            lines.append(f"    LISTED but not found by the rule: {f:%b %d} to {l:%b %d}")
    report = "\n".join(lines)
    (out_dir / "storm_features_report.txt").write_text(report + "\n")
    log("\n" + report)
    if not a.no_figures:
        for y in results:
            rws = [r for r in results[y]["rows"] if r.get("kind") == "listed"]
            log("  Figure: " + make_figure(out_dir / f"fig_storm_intensity_{y.lower()}.png", y, rws, results[y]["hours"], thresholds))
    log(f"\nWrote to {out_dir}/ : storm_features.csv, storm_features_report.txt, " + ", ".join(f"hourly_basin_{y.lower()}.csv" for y in results))
    return 0


# ================================================================================================================
# self-test: made-up storms with known answers
# ================================================================================================================
def _write_asc(path: Path, arr, nodata=-9999):
    ny, nx = arr.shape
    with open(path, "w") as f:
        f.write(f"ncols {nx}\nnrows {ny}\nxllcorner 440000\nyllcorner 3720000\ncellsize 1000\nNODATA_value {nodata}\n")
        np.savetxt(f, arr, fmt="%.3f")


def _write_test_boundary(path: Path, x0, y0, x1, y1):
    import geopandas as gpd
    from shapely.geometry import Polygon
    gpd.GeoDataFrame({"id": [1]}, geometry=[Polygon([(x0, y0), (x0, y1), (x1, y1), (x1, y0)])], crs="EPSG:26912").to_file(path)


def cmd_self_test(_a) -> int:
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("  PASS  " if cond else "  FAIL  ") + msg)
        ok = ok and bool(cond)

    tmp = Path(tempfile.mkdtemp(prefix="stormfeat_selftest_"))
    try:
        print("storm_features_sycamore SELF-TEST (made-up storms; your data is not touched)")
        geo = grid_geometry({"ncols": 32, "nrows": 49, "xllcorner": 440000, "yllcorner": 3720000, "cellsize": 1000})
        X, Y = geo["X"], geo["Y"]
        ones = np.ones((49, 32))
        thresholds = [2.2, 4.0]
        # ---- 1. a uniform storm with a known hourly pattern over a known day ----------------------------
        start = pd.Timestamp("2017-01-10 00:00")
        n = 24 * 20
        labels = pd.DatetimeIndex([start + HOUR * (i + 1) for i in range(n)])
        pat = np.zeros(n)
        day = 24 * 12                                         # Jan 22 00:00 is hour 288; storm on Jan 22-23 (days 12-13)
        pat[day + 6:day + 10] = [1.0, 6.0, 3.0, 1.0]          # a burst on the first day
        pat[day + 24 + 3:day + 24 + 5] = [2.0, 2.0]           # and a light spell the next day
        rain = (pat[:, None, None] * ones[None, :, :]).astype("float32")
        hours = hourly_table(labels, rain, geo, ones, set(), np.zeros(n, dtype=bool), thresholds)
        ctx = dict(wy_start="2016-10-01", wet_day_mm=2.0, daily_basin=daily_sum(hours, "basin_mean_mm"),
                   events=find_events(daily_sum(hours, "grid_mean_mm"), 0.5, 10.0))
        r = storm_features(hours, rain, geo, ones, "t", "2017-01-22", "2017-01-23", "listed", thresholds, ctx)
        check(abs(r["basin_mm"] - 15.0) < 1e-6 and abs(r["grid_mm"] - 15.0) < 1e-6, "uniform storm: basin and grid totals are both 15.0 mm")
        check(abs(r["peak_1h"] - 6.0) < 1e-6 and abs(r["peak_3h"] - 10.0) < 1e-6 and abs(r["peak_6h"] - 11.0) < 1e-6,
              f"peak 1 h = 6.0, peak 3 h = 10.0, peak 6 h = 11.0 (got {r['peak_1h']:.1f}, {r['peak_3h']:.1f}, {r['peak_6h']:.1f})")
        check(abs(r["wettest_day_mm"] - 11.0) < 1e-6 and r["hours_wet"] == 6, "wettest day 11.0 mm and 6 wet hours")
        check(r["hours_ge_2p2"] == 2 and abs(r["depth_ge_2p2"] - 9.0) < 1e-6 and r["hours_ge_4"] == 1,
              "2 hours at or above 2.2 mm/h holding 9.0 mm of rain; 1 hour at or above 4")
        check(abs(r["excess_2p2_mm"] - (3.8 + 0.8)) < 1e-6 and abs(r["excess_4_mm"] - 2.0) < 1e-6,
              f"excess over 2.2 mm/h = 3.8 + 0.8 = 4.6 mm; over 4 mm/h = 2.0 mm (got {r['excess_2p2_mm']:.2f}, {r['excess_4_mm']:.2f})")
        check(abs(r["top6h_share"] - 1.0) < 1e-9 and r["window_hours"] == 48, "all the rain falls in the six wettest hours; the window is 48 hours")
        check(abs(r["spatial_cv"]) < 1e-9 and abs(r["frac_basin_ge_20mm"]) < 1e-12, "uniform rain has zero spatial variation and no cell reaches 20 mm")
        check(abs(r["prev7_mm"]) < 1e-9 and r["lookback_incomplete"] is True and np.isnan(r["prev30_mm"]),
              "a storm with a dry week before it: prev7 = 0; a 30-day lookback that starts before the record is flagged, not guessed")
        # ---- 2. antecedent numbers ------------------------------------------------------------------------
        pat2 = pat.copy()
        pat2[24 * 8:24 * 8 + 4] = 2.5                         # 10 mm on Jan 18 (day 8), four days before the storm
        pat2[24 * 1:24 * 1 + 2] = 3.0                         # 6 mm on Jan 11, eleven days before
        rain2 = (pat2[:, None, None] * ones[None, :, :]).astype("float32")
        hours2 = hourly_table(labels, rain2, geo, ones, set(), np.zeros(n, dtype=bool), thresholds)
        ctx2 = dict(wy_start="2016-10-01", wet_day_mm=2.0, daily_basin=daily_sum(hours2, "basin_mean_mm"),
                    events=find_events(daily_sum(hours2, "grid_mean_mm"), 0.5, 5.0))
        r2 = storm_features(hours2, rain2, geo, ones, "t", "2017-01-22", "2017-01-23", "listed", thresholds, ctx2)
        check(abs(r2["prev7_mm"] - 10.0) < 1e-6, f"rain in the 7 days before the storm = 10.0 mm (got {r2['prev7_mm']})")
        check(np.isnan(r2["prev30_mm"]) and r2["lookback_incomplete"], "a 30-day lookback that reaches before the record is flagged incomplete, not guessed")
        check(abs(r2["dry_days_before"] - 3.0) < 1e-9, f"dry spell before the storm = 3 days (Jan 19, 20, 21; got {r2['dry_days_before']})")
        check(r2["prev_event"] == "Jan 18" and abs(r2["prev_event_mm"] - 10.0) < 1e-6 and abs(r2["prev_gap_days"] - 3.0) < 1e-9,
              f"previous event Jan 18, 10.0 mm, gap 3 dry days (got {r2['prev_event']}, {r2['prev_event_mm']}, {r2['prev_gap_days']})")
        # ---- 3. a moving band: direction and speed -------------------------------------------------------
        v = 4.0                                               # km per hour toward the east
        nh = 30
        band = np.zeros((nh, 49, 32), dtype="float32")
        for k in range(nh):
            xk = 428000.0 + v * 1000.0 * k
            band[k] = 6.0 * np.exp(-((X - xk) ** 2) / (2 * 3000.0 ** 2))
        fit = arrival_fit(band, geo, ones)
        check(abs(fit["arrival_dir"] - 90.0) < 3.0, f"a band moving east is reported as moving toward 90 deg (got {fit['arrival_dir']:.1f})")
        check(abs(fit["arrival_speed_kmh"] - v) < 0.4, f"...at {v} km/h (got {fit['arrival_speed_kmh']:.2f})")
        check(fit["arrival_r2"] > 0.95, f"...with an R2 near 1 (got {fit['arrival_r2']:.3f})")
        band_sw = np.zeros((nh + 20, 49, 32), dtype="float32")
        comp = v * 1000.0 * math.sqrt(0.5)
        for k in range(nh + 20):
            xk, yk = 462000.0 - comp * k * 1.0 + 4000.0, 3760000.0 - comp * k
            band_sw[k] = 6.0 * np.exp(-(((X - xk) ** 2 + (Y - yk) ** 2) * 0.5) / (2 * 3000.0 ** 2) * 2)
        fit2 = arrival_fit(band_sw, geo, ones)
        check(np.isfinite(fit2["arrival_dir"]) and abs(((fit2["arrival_dir"] - 225.0) + 180) % 360 - 180) < 15.0,
              f"a blob moving south-west is reported as moving toward about 225 deg (got {fit2['arrival_dir']:.0f})")
        flat = np.ones((12, 49, 32), dtype="float32")         # everything arrives at once
        fit3 = arrival_fit(flat, geo, ones)
        check(np.isnan(fit3["arrival_speed_kmh"]), "rain that arrives everywhere at once gets no speed (weak or no fit)")
        # ---- 4. distance skew: rain only near the gauge -------------------------------------------------
        near = np.zeros((6, 49, 32), dtype="float32")
        d_g = np.hypot(X - GAUGE_XY[0], Y - GAUGE_XY[1])
        near[:, d_g < 6000.0] = 5.0
        lab6 = pd.DatetimeIndex([pd.Timestamp("2017-01-22 00:00") + HOUR * (i + 1) for i in range(6)])
        h6 = hourly_table(lab6, near, geo, ones, set(), np.zeros(6, dtype=bool), thresholds)
        ctx6 = dict(wy_start="2016-10-01", wet_day_mm=2.0, daily_basin=daily_sum(h6, "basin_mean_mm"), events=[])
        r6 = storm_features(h6, near, geo, ones, "near", "2017-01-22", "2017-01-22", "listed", thresholds, ctx6)
        check(r6["rain_dist_km"] < 4.5 and r6["dist_skew_km"] < -5.0, f"rain only within 6 km of the gauge: rain-weighted distance {r6['rain_dist_km']:.1f} km, skew {r6['dist_skew_km']:.1f}")
        check(r6["spatial_cv"] > 2.0, f"concentrated rain gives a large spatial variation (CV {r6['spatial_cv']:.1f})")
        # ---- 5. events rule ------------------------------------------------------------------------------
        dly = pd.Series([0, 3, 30, 0.6, 0, 0, 12, 12, 0, 25], index=pd.date_range("2017-02-01", periods=10))
        ev = find_events(dly, 0.5, 20.0)
        check(len(ev) == 3 and abs(ev[0][2] - 33.6) < 1e-9 and ev[0][0] == pd.Timestamp("2017-02-02") and ev[0][1] == pd.Timestamp("2017-02-04"),
              "the day-run rule groups Feb 2-4 (3 + 30 + 0.6 = 33.6 mm), then 24 mm and 25 mm")
        # ---- 6. end to end on small made-up files ------------------------------------------------------------
        rd = tmp / "wk" / "rain"
        rd.mkdir(parents=True)
        s0 = pd.Timestamp("2016-12-14 00:00")
        nn = 24 * 14
        labs = [s0 + HOUR * (i + 1) for i in range(nn)]
        rng = np.random.default_rng(3)
        for i, lab in enumerate(labs):
            fld = np.zeros((49, 32))
            if 24 * 7 <= i < 24 * 9:                          # Dec 21-22: 2 mm/h everywhere for six hours a day
                if (i % 24) in range(8, 14):
                    fld[:] = 2.0 + rng.uniform(0, 0.2, fld.shape)
            if i == 5:
                continue                                      # one file is missing on purpose
            _write_asc(file_for_label(rd, lab), fld)
        (tmp / "wk" / "manifest.csv").write_text("label_local,status\n" + f"{labs[30]:%Y-%m-%d %H:%M},filled_zero_missing_hour\n")
        shp = tmp / "basin.shp"
        _write_test_boundary(shp, 450000, 3735000, 460000, 3745000)
        sfile = tmp / "storms.csv"
        sfile.write_text("year,name,first_day,last_day,kind\nWY2017,Dec 21-22,2016-12-21,2016-12-22,listed\n")
        out = tmp / "out"
        argv = ["run", "--years", "wy2017", "--rain-dir-wy2017", str(rd), "--window-wy2017", f"{s0:%Y-%m-%d %H:%M},{s0 + HOUR * nn:%Y-%m-%d %H:%M}",
                "--boundary", str(shp), "--storms-file", str(sfile), "--out-dir", str(out), "--expect-area-km2", "100"]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = main(argv)
        txt = buf.getvalue()
        csvf = out / "storm_features.csv"
        check(rc == 0 and csvf.exists() and (out / "hourly_basin_wy2017.csv").exists() and (out / "storm_features_report.txt").exists()
              and (out / "fig_storm_intensity_wy2017.png").exists(), "an end-to-end run writes the feature CSV, the hourly CSV, the report and the figure")
        d = pd.read_csv(csvf)
        check(len(d) == 1 and abs(d["basin_mm"].iat[0] - (2.1 * 12)) < 1.0, f"end to end: 12 wet hours near 2.1 mm/h give about 25 mm (got {d['basin_mm'].iat[0]:.1f})")
        check(d["radar_missing_hours"].iat[0] == 0, "the missing file (hour 6 of the record) is outside the storm window, so it is not counted for the storm")
        hb = pd.read_csv(out / "hourly_basin_wy2017.csv")
        check(int(hb["missing_file"].sum()) == 1 and int(hb["filled_zero"].sum()) == 1, "the hourly CSV flags 1 missing file and 1 zero-filled hour")
        check("LISTED but not found by the rule" not in txt and "on the list" in txt, "the rule check finds the listed storm")
        # dry run writes nothing
        out2 = tmp / "out2"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = main(["run", "--dry-run", "--years", "wy2017", "--rain-dir-wy2017", str(rd), "--window-wy2017",
                       f"{s0:%Y-%m-%d %H:%M},{s0 + HOUR * nn:%Y-%m-%d %H:%M}", "--boundary", str(shp), "--out-dir", str(out2)])
        check(rc == 0 and not out2.exists() and "335 of 336 files found" in buf.getvalue(), "--dry-run counts 335 of 336 files and writes nothing")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\nSELF-TEST " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


# ================================================================================================================
# command line
# ================================================================================================================
def build_parser():
    p = argparse.ArgumentParser(description="Rain-only storm features from the hourly radar grids (see the text at the top of this file).",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="build the feature table, report, hourly CSVs and figures")
    r.add_argument("--years", default="wy2017,wy2020", help="which years (default wy2017,wy2020)")
    r.add_argument("--rain-dir-wy2017", default=None, help=f"WY2017 grids (default {YEARS['WY2017']['rain_dir']})")
    r.add_argument("--rain-dir-wy2020", default=None, help=f"WY2020 grids (default {YEARS['WY2020']['rain_dir']})")
    r.add_argument("--window-wy2017", default=None, help="'start,end' of the WY2017 record to read (default whole water year)")
    r.add_argument("--window-wy2020", default=None, help="'start,end' of the WY2020 record to read (default Sep 1 2019 to Jun 1 2020)")
    r.add_argument("--boundary", default=DEFAULT_BOUNDARY, help=f"basin outline in EPSG:26912 (default {DEFAULT_BOUNDARY})")
    r.add_argument("--manifest", default=None, help="manifest.csv from the radar build (default: one folder above each rain folder)")
    r.add_argument("--storms-file", default=None, help="CSV with year,name,first_day,last_day[,kind] to use instead of the built-in list")
    r.add_argument("--ratios", default=None, help="w5_rain_by_storm.csv, to add the radar/gauge ratios and check the grid totals")
    r.add_argument("--thresholds", default="2.2,4", help="rain rates in mm/h that the soil can take (default 2.2,4)")
    r.add_argument("--cell-min", type=float, default=0.1, help="a cell counts as wet from this mm/h (default 0.1)")
    r.add_argument("--wet-day-mm", type=float, default=2.0, help="a 'wet day' has at least this much basin rain (default 2)")
    r.add_argument("--prev-event-mm", type=float, default=10.0, help="the 'previous event' is a rain event of at least this many mm (default 10)")
    r.add_argument("--out-dir", default="storm_features_out", help="folder for the outputs (default storm_features_out)")
    r.add_argument("--no-figures", action="store_true", help="skip the figures")
    r.add_argument("--expect-area-km2", type=float, default=424.6, help="warn if the outline's area on the grid differs by more than 3 percent (0 = off)")
    r.add_argument("--dry-run", action="store_true", help="list what would be read and written; read and write nothing else")
    sub.add_parser("self-test", help="check the script on made-up storms")
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    return cmd_run(a) if a.cmd == "run" else cmd_self_test(a)


if __name__ == "__main__":
    sys.exit(main())
