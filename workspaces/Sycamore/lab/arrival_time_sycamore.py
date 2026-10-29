#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
arrival_time_sycamore.py
========================
When did the rain arrive at each part of the Lower Sycamore basin, and which way did it move?
(Hourly MRMS radar grids, ONE widespread storm at a time.)

WHY THIS EXISTS
---------------
storm_movement_sycamore.py follows the rain "centre of mass" (centroid). That works for a compact storm that fits
inside the 32 x 49 km radar grid. A winter front fills the whole grid and runs off its edges, so the centroid mostly
tracks where the grid stops. This script asks a different question that does not have that problem: for EVERY 1 km
cell, at what clock time had the cell received half of its own storm total? Mapping that single time per cell shows
which side of the basin got its rain first and how fast the rain spread across.

WHAT IT MAKES  (five files per run, in --out-dir, all starting with the --tag)
------------------------------------------------------------------------------
  <tag>_arrival_overview.png    one page: (left) map of the half-way time of every cell with an arrow for the fitted
                                direction of movement; (top right) the same times plotted against distance along that
                                arrow, to show how straight the pattern really is; (middle and bottom right) basin
                                rain over time for the east half vs the west half, and the north half vs the south half
  <tag>_arrival_cells.csv       one row per radar cell: position, share inside the basin, window total, and the times
                                (hours after the window start) when it had 25%, 50% and 75% of its rain; load it in a
                                notebook to make your own maps (row 0 = the NORTH row)
  <tag>_arrival_fits.csv        the straight-front fits (direction, speed, R2) for 2 areas x 3 fractions
  <tag>_arrival_halves.csv      hourly mean rain over the west / east / north / south halves of the basin
  <tag>_arrival_run_info.txt    the settings used and a plain-language summary with warnings (also printed on screen)

QUICK START  (run from the lab/ folder of workspaces/Sycamore; storm_movement_sycamore.py version 2 must be in the
-----------   same folder because this script borrows its file reading and basin-outline code)
    # 0. check the script itself (made-up storms with known answers; touches none of your data)
    python arrival_time_sycamore.py self-test

    # 1. look before you leap: lists the files it will read and writes NOTHING
    python arrival_time_sycamore.py run --dry-run

    # 2. the default run = the main January 2017 storm (Jan 20 06:00 to Jan 21 08:00 MST)
    python arrival_time_sycamore.py run

    # 3. any other single storm (whole hours, Arizona time), e.g. the smaller Jan 23-24 one
    python arrival_time_sycamore.py run --start "2017-01-23 12:00" --end "2017-01-25 00:00" --tag jan23

WORDS USED
----------
  half-way time   for ONE cell: add up its hourly rain over the window to get its window total, then find the moment
                  the running total first reaches half of that total. (Also called the half-mass time or the arrival
                  time of the storm at that cell.) Inside an hour the rain is assumed steady, so times fall between
                  the hour marks (for example 18:20).
  front           the leading edge of a band of rain as it moves across the land.
  straight-front  the best-fitting flat, tilted surface through the map of half-way times:
  fit             time = c + a x east + b x north (least squares, each cell counted in proportion to its rain).
                  The direction in which time rises fastest is the direction the rain MOVES TOWARD; the speed is
                  1 / (how many hours the time changes per km).
  bearing         direction in degrees clockwise from north: 0 = north, 90 = east, 180 = south, 270 = west. The
                  script reports where the rain moves TO and, in words, where it comes FROM (the opposite way): moving
                  toward the SW (225 degrees) is the same as arriving from the NE.
  R2              how much of the cell-to-cell variation in half-way time the straight front explains. 1 = perfect,
                  0 = none. Rule of thumb used here (not a law): 0.6 or more = good, 0.3 to 0.6 = moderate, under 0.3
                  = weak (do not quote the speed).
  misfit          the typical gap (hours) between a cell's real half-way time and the fitted front.
  half of basin   the basin outline is split by a straight north-south (or east-west) line into two halves of EQUAL
                  AREA (cells cut by the line are shared out in proportion). The hyetograph of each half is its mean
                  rain per hour.
  lag             how much later one half of the basin reached its half-way time than the other half. A second,
                  rougher check slides one hyetograph over the other, one whole hour at a time, and keeps the shift
                  where they match best (cross-correlation).

HOW THE NUMBERS ARE MADE
------------------------
  1. Read the hourly grids for the window and add them up cell by cell (running total).
  2. Cells with less than --min-total mm (default 5) in the window are left out: the timing of almost no rain is noise.
  3. Half-way time per cell = where the cell's running total crosses 25%, 50% and 75% of its window total, found by
     straight-line interpolation inside the hour. 50% is the one that is mapped; 25% and 75% are used only to test
     whether the answer depends on which fraction you pick.
  4. Straight-front fit to the 50% times of the basin cells (those at least half inside the outline). The same fit is
     repeated for the whole radar grid and for 25% and 75%, giving 6 fits. If their directions all agree (within 30
     degrees) the direction is called stable.
  5. West/east and north/south halves: half-way time of each half's hyetograph and their difference, plus the
     cross-correlation lag (hourly resolution) as a check.

LIMITS  (read before interpreting)
----------------------------------
  * Use ONE storm per window. If the window holds two separate bursts, "half of the rain" is a meaningless moment.
    The script warns when the first or last hour of the window still holds more than 3% of the basin rain (the storm
    probably continues past the window: widen --start/--end).
  * Frames are one hour apart, so differences of well under an hour between places are at the limit of the data. The
    script says so when the fitted time difference across the basin is under 1 hour.
  * A straight front is a simplification. Rain bands that curve, speed up, or break into cells give a low R2.
  * Radar timing has not been checked against gauges (the Stage 0 radar-gauge check is on hold). Trust the pattern more
    than the depths.
  * tRIBS is given these same hourly grids, so this shows the timing the model saw.
  * An hour that MRMS had no data for (zero-filled, listed in manifest.csv) shifts every cell's timing after it. The
    script lists any such hour inside the window. None fall in January-March 2017.

CONVENTIONS
-----------
  * Arizona time (MST, UTC-7, no daylight saving). --start and --end are the rain interval [start, end).
  * File names are radar_MMDDYYYYHH.asc, HH = the hour that ENDS (the hour ending midnight is 00 of the NEXT day).
  * EPSG:26912 (UTM 12N). The gauge is USGS 09510200 at x 449721, y 3728390.

NEEDS: numpy, pandas, matplotlib, Pillow, ONE of pyshp / geopandas, and storm_movement_sycamore.py (version 2).
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

try:
    from storm_movement_sycamore import (
        BLUE, DEFAULT_BOUNDARY, DEFAULT_RAIN_DIR, GAUGE_XY, GRID, HOUR, INK, INK2, MUTED, RED, SURFACE,
        _km_axes, _mpl, _time_ticks, _write_asc, _write_test_boundary, compass, die, file_for_label, grid_cover,
        grid_geometry, load_boundary_rings, load_stack, parse_window, read_ascii, read_filled_keys)
except ImportError as exc:                                      # pragma: no cover
    raise SystemExit(
        "ERROR: arrival_time_sycamore.py needs storm_movement_sycamore.py (VERSION 2, the one dated after the "
        "colour-bar fix) in the same folder.\n       Python said: " + str(exc))

ORANGE = "#eb6834"                      # categorical slot 2; validated together with BLUE (adjacent pair, light mode)
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
             "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]       # one-hue ramp, light -> dark
FRACTIONS = (0.25, 0.5, 0.75)
DEFAULT_START, DEFAULT_END = "2017-01-20 06:00", "2017-01-21 08:00"                    # the main January 2017 storm
EDGE_SHARE_WARN = 0.03                  # warn if the first or last hour holds more than this share of the basin rain
STABLE_DEG = 30.0                       # directions of the six fits within this many degrees = "stable"
MAX_SPEED_KMH = 150.0                   # a fitted speed above this is called "no usable direction"


# ================================================================================================================
# the numbers
# ================================================================================================================
def cumulative(rain):
    """(n, ny, nx) hourly depths -> (n+1, ny, nx) running totals in mm; index 0 = window start, index k = END of hour k."""
    cum = np.zeros((rain.shape[0] + 1,) + tuple(rain.shape[1:]), dtype="float64")
    np.cumsum(rain, axis=0, dtype="float64", out=cum[1:])
    return cum


def arrival_hours(cum, frac, min_total):
    """Hours after the window start when each cell had received `frac` of its own window total (NaN if the total is
    below `min_total`). Rain is assumed steady inside each hour, so the running total is a straight line between hour
    marks and the crossing is found by linear interpolation."""
    total = cum[-1]
    ok = total >= max(float(min_total), 1e-9)
    target = frac * total
    first = (cum >= target[None, ...]).argmax(axis=0)                 # first hour mark whose running total reaches it
    first = np.maximum(first, 1)
    lo = np.take_along_axis(cum, (first - 1)[None, ...], axis=0)[0]
    hi = np.take_along_axis(cum, first[None, ...], axis=0)[0]
    span = hi - lo
    f = np.divide(target - lo, span, out=np.ones_like(span), where=span > 0)
    t = (first - 1) + np.clip(f, 0.0, 1.0)
    return np.where(ok, t, np.nan)


def arrival_1d(series, frac=0.5):
    """Half-way (or other fraction) time, in hours after the start, of one hourly series. NaN if there is no rain."""
    cum = cumulative(np.asarray(series, dtype="float64")[:, None, None])
    return float(arrival_hours(cum, frac, 0.0)[0, 0])


def fit_plane(x_km, y_km, t_h, w):
    """Weighted least-squares straight front  t = c + a*(x-xm) + b*(y-ym)  (x east, y north, km; t hours).
    The rain moves toward the direction (a, b): bearing = clockwise from north, speed = 1/|(a, b)| km/h."""
    x_km, y_km, t_h, w = (np.asarray(v, dtype="float64") for v in (x_km, y_km, t_h, w))
    n = len(t_h)
    out = dict(ok=False, n=n, bearing_deg=float("nan"), speed_kmh=float("nan"), r2=float("nan"), rms_h=float("nan"),
               across_h=float("nan"), grad_h_per_km=float("nan"), a=float("nan"), b=float("nan"), c=float("nan"),
               xm=float("nan"), ym=float("nan"))
    if n < 10 or w.sum() <= 0:
        return out
    xm, ym = float(np.average(x_km, weights=w)), float(np.average(y_km, weights=w))
    A = np.column_stack([np.ones(n), x_km - xm, y_km - ym])
    sw = np.sqrt(w)
    coef, _, rank, _ = np.linalg.lstsq(A * sw[:, None], t_h * sw, rcond=None)
    if rank < 3:
        return out
    c, a, b = (float(v) for v in coef)
    resid = t_h - A @ coef
    tbar = float(np.average(t_h, weights=w))
    sst = float((w * (t_h - tbar) ** 2).sum())
    ssr = float((w * resid ** 2).sum())
    g = math.hypot(a, b)
    out.update(ok=True, c=c, a=a, b=b, xm=xm, ym=ym, grad_h_per_km=g,
               r2=(1.0 - ssr / sst) if sst > 0 else float("nan"), rms_h=math.sqrt(ssr / float(w.sum())))
    if g > 1e-9:
        out["bearing_deg"] = (math.degrees(math.atan2(a, b)) + 360.0) % 360.0
        out["speed_kmh"] = 1.0 / g
        s = (x_km - xm) * (a / g) + (y_km - ym) * (b / g)             # km along the direction of movement
        out["across_h"] = g * float(np.percentile(s, 95) - np.percentile(s, 5))
    return out


def has_direction(fit) -> bool:
    return bool(fit["ok"] and np.isfinite(fit["bearing_deg"]) and fit["speed_kmh"] <= MAX_SPEED_KMH)


def fit_verdict(fit) -> str:
    if not fit["ok"]:
        return "NO FIT: too few cells with enough rain"
    r2 = fit["r2"]
    if not np.isfinite(r2) or r2 < 0.3:
        return "WEAK FIT (R2 under 0.3): the timing is not described by one straight front; do not quote a speed or direction"
    if r2 < 0.6:
        return "MODERATE FIT (R2 0.3 to 0.6): a straight front explains part of the timing; quote the direction only as a tendency"
    return "GOOD FIT (R2 0.6 or more): a straight front describes the timing well"


def heading_words(fit) -> str:
    if not fit["ok"] or not np.isfinite(fit["bearing_deg"]):
        return "no direction could be fitted"
    if fit["speed_kmh"] > MAX_SPEED_KMH:
        return f"no usable direction: the half-way times barely change across the area (fitted speed over {MAX_SPEED_KMH:.0f} km/h)"
    b = fit["bearing_deg"]
    return (f"rain moving toward the {compass(b)} ({b:.0f} deg), i.e. arriving from the {compass((b + 180) % 360)}, "
            f"at about {fit['speed_kmh']:.1f} km/h")


def bearing_spread(fits) -> float:
    """Largest angle (degrees) between any fit's direction and the average direction of all usable fits."""
    b = [f["bearing_deg"] for f in fits.values() if has_direction(f)]
    if len(b) < 2:
        return float("nan")
    r = np.radians(b)
    m = math.degrees(math.atan2(float(np.sin(r).mean()), float(np.cos(r).mean())))
    return max(abs((x - m + 180.0) % 360.0 - 180.0) for x in b)


def half_weights(cover, axis):
    """Split the basin into two halves of equal area with a straight line. `cover` = share (0-1) of each cell inside the
    outline. Returns two weight grids that add up to `cover`: axis 'x' -> (west, east); axis 'y' -> (north, south).
    Cells the line cuts are shared out in proportion."""
    line = cover.sum(axis=0) if axis == "x" else cover.sum(axis=1)      # per column (west->east) or per row (north->south)
    half = float(line.sum()) / 2.0
    cum = np.cumsum(line)
    j = int(np.argmax(cum >= half - 1e-12))
    prev = float(cum[j - 1]) if j > 0 else 0.0
    frac = (half - prev) / float(line[j]) if line[j] > 0 else 0.0
    share = np.zeros(len(line))
    share[:j] = 1.0
    share[j] = min(max(frac, 0.0), 1.0)
    first = cover * (share[None, :] if axis == "x" else share[:, None])
    return first, cover - first


def half_series(rain, w):
    """Hourly mean rain (mm/hr) over the area described by the weight grid `w`."""
    return (rain * w[None, ...]).sum(axis=(1, 2)) / float(w.sum())


def xcorr_lag(a, b, max_lag):
    """Whole-hour shift that makes series b best match series a. Positive = b happens LATER than a.
    Returns (lag, correlation); (nan, nan) if a series is flat or the overlap is too short."""
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    n = len(a)
    best = (float("nan"), float("nan"))
    for lag in range(-int(max_lag), int(max_lag) + 1):
        x, y = (a[:n - lag], b[lag:]) if lag >= 0 else (a[-lag:], b[:n + lag])
        if len(x) < 6 or x.std() == 0 or y.std() == 0:
            continue
        c = float(np.corrcoef(x, y)[0, 1])
        if not np.isfinite(best[1]) or c > best[1] + 1e-12:
            best = (lag, c)
    return best


def analyze(rain, geo, cover, min_total=5.0, max_lag=6):
    """Every number the script reports. `rain` is (n, ny, nx) mm per hour with NaN already replaced by 0."""
    cum = cumulative(rain)
    total = cum[-1]
    Xk, Yk = geo["X"] / 1000.0, geo["Y"] / 1000.0
    inside = cover >= 0.5
    t = {p: arrival_hours(cum, p, min_total) for p in FRACTIONS}
    fits = {}
    for area, mask in (("basin", inside), ("grid", np.ones(inside.shape, dtype=bool))):
        for p in FRACTIONS:
            m = mask & np.isfinite(t[p])
            fits[(area, p)] = fit_plane(Xk[m], Yk[m], t[p][m], total[m])
    basin_mean = (rain * (cover / cover.sum())[None, ...]).sum(axis=(1, 2))
    tot_b = float(basin_mean.sum())
    halves = {}
    for axis, names in (("x", ("west", "east")), ("y", ("north", "south"))):
        w1, w2 = half_weights(cover, axis)
        s1, s2 = half_series(rain, w1), half_series(rain, w2)
        t1, t2 = arrival_1d(s1), arrival_1d(s2)
        lag, corr = xcorr_lag(s1, s2, max_lag)
        halves[axis] = dict(names=names, s1=s1, s2=s2, t1=t1, t2=t2, diff_h=t2 - t1, lag_h=lag, lag_corr=corr)
    return dict(n=rain.shape[0], total=total, t=t, fits=fits, inside=inside, basin_mean=basin_mean, halves=halves,
                first_share=float(basin_mean[0] / tot_b) if tot_b > 0 else float("nan"),
                last_share=float(basin_mean[-1] / tot_b) if tot_b > 0 else float("nan"),
                n_cells=int(np.isfinite(t[0.5]).sum()), n_basin=int(inside.sum()),
                n_cells_basin=int((inside & np.isfinite(t[0.5])).sum()))


def clock(start, hours) -> str:
    return "n/a" if not np.isfinite(hours) else f"{start + pd.Timedelta(hours=float(hours)):%b %d %H:%M}"


def later_sentence(name1, t1, name2, t2) -> str:
    d = t2 - t1
    if not np.isfinite(d):
        return "no half-way time could be found (not enough rain)"
    if abs(d) < 0.25:
        return f"the two halves were within 15 minutes (difference {d:+.2f} h)"
    later, earlier = (name2, name1) if d > 0 else (name1, name2)
    return f"the {later} half was {abs(d):.1f} h LATER than the {earlier} half"


def summarize(res, labels, a, area_km2, filled_in_window) -> list[str]:
    start = labels[0] - HOUR
    L = []
    n = res["n"]
    L.append(f"window (rain intervals, MST): {start:%Y-%m-%d %H:%M} to {labels[-1]:%Y-%m-%d %H:%M}  ({n} hourly frames)")
    L.append(f"basin area from 1 km coverage: {area_km2:.1f} km2")
    L.append(f"basin-mean rain over the window: {res['basin_mean'].sum():.1f} mm")
    L.append(f"cells with at least {a.min_total:g} mm in the window: {res['n_cells']} of {res['total'].size} on the grid "
             f"({res['n_cells_basin']} of {res['n_basin']} basin cells)")
    fs, ls = res["first_share"], res["last_share"]
    L.append(f"WINDOW CHECK: the first hour holds {fs:.1%} and the last hour {ls:.1%} of the basin rain")
    if (np.isfinite(fs) and fs > EDGE_SHARE_WARN) or (np.isfinite(ls) and ls > EDGE_SHARE_WARN):
        L.append(f"  WARNING: more than {EDGE_SHARE_WARN:.0%} falls in the first or last hour, so the storm probably continues past "
                 f"the window; the half-way times are biased. Widen --start / --end and run again.")
    if filled_in_window:
        L.append(f"  WARNING: MRMS had no data (zero-filled) for {len(filled_in_window)} hour(s) inside the window: "
                 f"{', '.join(filled_in_window)}. Timing after those hours is shifted.")
    t50 = res["t"][0.5]
    sel = res["inside"] & np.isfinite(t50)
    L.append("")
    L.append("WHEN THE RAIN ARRIVED  (half-way time = when a cell had received half of its window total; basin cells)")
    if sel.sum() >= 2:
        p10, p90 = (float(v) for v in np.percentile(t50[sel], [10, 90]))
        L.append(f"  the earliest 10% of basin cells reached it by {clock(start, p10)}, the latest 10% by {clock(start, p90)}  "
                 f"(a spread of {p90 - p10:.1f} h)")
    main = res["fits"][("basin", 0.5)]
    L.append(f"  straight-front fit to the basin cells: {heading_words(main)}")
    if has_direction(main):
        L.append(f"    fit quality: R2 {main['r2']:.2f}, typical misfit {main['rms_h']:.2f} h")
        L.append(f"    {fit_verdict(main)}")
        L.append(f"    the fitted half-way time changes by about {main['across_h']:.1f} h from one side of the basin to the other")
        if main["across_h"] < 1.0:
            L.append("    NOTE: under 1 hour across the basin. With hourly frames the storm reached the whole basin "
                     "essentially at once; the direction above is at the limit of what the data can show.")
    elif main["ok"]:
        L.append("    the half-way times are almost the same everywhere in the basin: with hourly frames the storm reached "
                 "the whole basin essentially at once, so no direction of movement is reported")
    L.append("")
    L.append("DOES THE ANSWER DEPEND ON HOW IT IS MEASURED?  (same straight-front fit, different cells and fractions)")
    L.append(f"    {'cells':<6} {'at':>4}   {'toward':<9} {'km/h':>6}   {'R2':>5}   {'cells used':>10}")
    for (area, p), f in res["fits"].items():
        if f["ok"] and np.isfinite(f["bearing_deg"]):
            if f["speed_kmh"] <= MAX_SPEED_KMH:
                dir_txt, sp = f"{compass(f['bearing_deg'])} ({f['bearing_deg']:.0f})", f"{f['speed_kmh']:6.1f}"
            else:
                dir_txt, sp = "none", f">{MAX_SPEED_KMH:.0f}".rjust(6)
            L.append(f"    {area:<6} {int(round(p * 100)):>3d}%   {dir_txt:<9} {sp}   {f['r2']:5.2f}   {f['n']:10d}")
        else:
            L.append(f"    {area:<6} {int(round(p * 100)):>3d}%   {'no fit':<9} {'':>6}   {'':>5}   {f['n']:10d}")
    sp_deg = bearing_spread(res["fits"])
    if np.isfinite(sp_deg):
        L.append(f"  direction across the usable fits: largest departure from the average {sp_deg:.0f} deg -> "
                 + ("STABLE (within 30 deg)" if sp_deg <= STABLE_DEG else
                    "NOT STABLE (over 30 deg): the direction depends on how it is measured; treat it as unresolved"))
    for axis, title in (("x", "EAST vs WEST halves of the basin (equal areas)"), ("y", "NORTH vs SOUTH halves of the basin (equal areas)")):
        h = res["halves"][axis]
        n1, n2 = h["names"]
        L.append("")
        L.append(title)
        L.append(f"  {n1} half reached half its rain at {clock(start, h['t1'])}; {n2} half at {clock(start, h['t2'])}: "
                 + later_sentence(n1, h["t1"], n2, h["t2"]))
        lag = h["lag_h"]
        if np.isfinite(lag):
            word = "no shift" if lag == 0 else (f"{n2} later by {int(lag)} h" if lag > 0 else f"{n1} later by {int(-lag)} h")
            L.append(f"  cross-correlation check (whole hours): best match at {word}  (correlation {h['lag_corr']:.2f})")
            if np.isfinite(h["diff_h"]) and abs(lag - h["diff_h"]) > 1.0:
                L.append("  NOTE: the two lag measures differ by more than an hour, which usually means more than one burst of rain; "
                         "trust neither number fully.")
        else:
            L.append("  cross-correlation check: not available (a flat series or too few hours)")
    return L


# ================================================================================================================
# figure
# ================================================================================================================
def _outline_halo(ax, rings):
    for r in rings:
        ax.plot(r[:, 0], r[:, 1], color="white", lw=3.6, zorder=4, solid_capstyle="round")
        ax.plot(r[:, 0], r[:, 1], color=INK, lw=1.3, zorder=5)


def _half_panel(ax, h, s0, nh, title, ylabel):
    x = s0 + np.arange(nh + 1) / 24.0
    for ser, col, nm, tt in ((h["s1"], BLUE, h["names"][0], h["t1"]), (h["s2"], ORANGE, h["names"][1], h["t2"])):
        ax.step(x, np.append(ser, ser[-1]), where="post", color=col, lw=1.8, label=f"{nm} half")
        if np.isfinite(tt):
            ax.axvline(s0 + tt / 24.0, ymax=0.7, color=col, lw=1.2, ls="--")
    ymax = float(max(np.max(h["s1"]), np.max(h["s2"]), 1e-6))
    ax.set_ylim(0, ymax * 1.45)
    ax.set_xlim(s0, s0 + nh / 24.0)
    ax.set_title(title)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_axisbelow(True)
    ax.grid(color=GRID)
    ax.legend(loc="upper right", fontsize=8)
    n1, n2 = h["names"]
    lag = h["lag_h"]
    lag_txt = "n/a" if not np.isfinite(lag) else f"{int(lag):+d} h"
    ax.text(0.01, 0.97, f"dashed = half-way time.  {later_sentence(n1, h['t1'], n2, h['t2'])}\n"
                        f"cross-correlation shift ({n2} minus {n1}): {lag_txt}",
            transform=ax.transAxes, va="top", ha="left", fontsize=7.5, color=INK2,
            bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none", alpha=0.9))


def make_figure(res, geo, rings, labels, a, out_png, title):
    plt = _mpl()
    import matplotlib.dates as mdates
    import matplotlib.patheffects as pe
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrowPatch
    start = labels[0] - HOUR
    s0 = float(mdates.date2num(start.to_pydatetime()))
    nh = res["n"]
    t50 = res["t"][0.5]
    valid = np.isfinite(t50)
    fit = res["fits"][("basin", 0.5)]
    Xk, Yk = geo["X"] / 1000.0, geo["Y"] / 1000.0
    fig = plt.figure(figsize=(11.8, 8.6))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.6], hspace=0.5, wspace=0.26)
    axm = fig.add_subplot(gs[:, 0])
    axt = fig.add_subplot(gs[0, 1])
    axe = fig.add_subplot(gs[1, 1])
    axn = fig.add_subplot(gs[2, 1], sharex=axe, sharey=axe)
    # ---- map of half-way times --------------------------------------------------------------------------------
    cmap = LinearSegmentedColormap.from_list("blue_seq", BLUE_RAMP).copy()
    cmap.set_bad(alpha=0.0)
    axm.set_facecolor("#f1f0ec")
    tnum = s0 + t50 / 24.0
    vmin, vmax = (float(v) for v in np.nanpercentile(tnum[valid], [1, 99]))
    if vmax - vmin < 1.0 / 24.0:
        vmin, vmax = vmin - 0.5 / 24.0, vmax + 0.5 / 24.0
    ext = [geo["L"], geo["R"], geo["B"], geo["T"]]
    im = axm.imshow(np.ma.masked_invalid(tnum), extent=ext, origin="upper", cmap=cmap, vmin=vmin, vmax=vmax,
                    interpolation="nearest")
    cb = fig.colorbar(im, ax=axm, orientation="horizontal", fraction=0.045, pad=0.11)
    pos, lab = _time_ticks(vmin, vmax, max_ticks=4, two_line=True)
    cb.set_ticks(pos)
    cb.set_ticklabels(lab)
    cb.ax.tick_params(labelsize=7.5)
    cb.set_label("time the cell had received half its rain (MST)")
    _outline_halo(axm, rings)
    axm.plot(*GAUGE_XY, marker="*", ms=13, color=RED, mec="white", mew=0.8, ls="none", zorder=6)
    axm.annotate("gauge", xy=GAUGE_XY, xytext=(9, -5), textcoords="offset points", fontsize=8, color=INK, zorder=9,
                 bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none", alpha=0.9))
    if has_direction(fit) and fit["r2"] >= 0.3:
        g = fit["grad_h_per_km"]
        ux, uy = fit["a"] / g, fit["b"] / g
        cx, cy = fit["xm"] * 1000.0, fit["ym"] * 1000.0
        half_len = 6000.0
        arr = FancyArrowPatch((cx - ux * half_len, cy - uy * half_len), (cx + ux * half_len, cy + uy * half_len),
                              arrowstyle="-|>", mutation_scale=18, lw=2.4, color=INK, zorder=8)
        arr.set_path_effects([pe.withStroke(linewidth=5.0, foreground="white")])
        axm.add_patch(arr)
        axm.text(cx, cy - 7500.0, f"toward {compass(fit['bearing_deg'])}", fontsize=8, ha="center", va="top", color=INK,
                 zorder=9, bbox=dict(boxstyle="round,pad=0.15", fc=SURFACE, ec="none", alpha=0.85))
    axm.text(0.02, 0.985, f"grey: under {a.min_total:g} mm (no timing)\ncolours: 1st-99th percentile", transform=axm.transAxes,
             va="top", fontsize=7.5, color=INK2, zorder=9, bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none", alpha=0.85))
    axm.set_xlim(geo["L"], geo["R"])
    axm.set_ylim(geo["B"], geo["T"])
    _km_axes(axm)
    axm.set_title("When each cell had received half its rain")
    # ---- half-way time against distance along the arrow ---------------------------------------------------------
    axt.set_title("Half-way time vs distance along the arrow")
    if has_direction(fit):
        g = fit["grad_h_per_km"]
        ux, uy = fit["a"] / g, fit["b"] / g
        s_km = (Xk - fit["xm"]) * ux + (Yk - fit["ym"]) * uy
        ins = res["inside"] & valid
        out = valid & ~res["inside"]
        axt.scatter(s_km[out], tnum[out], s=9, color=MUTED, alpha=0.45, linewidths=0, zorder=2)
        axt.scatter(s_km[ins], tnum[ins], s=11, color=BLUE, alpha=0.55, linewidths=0, zorder=3)
        sl = np.array([float(s_km[ins].min()), float(s_km[ins].max())]) if ins.any() else np.array([-10.0, 10.0])
        axt.plot(sl, s0 + (fit["c"] + g * sl) / 24.0, color=INK, lw=1.8, ls="--", zorder=4)
        axt.set_xlabel(f"distance along the arrow, toward {compass(fit['bearing_deg'])} (km)")
        lo_t, hi_t = float(np.nanmin(tnum[valid])), float(np.nanmax(tnum[valid]))
        pad = max(0.02, 0.05 * (hi_t - lo_t))
        axt.set_ylim(lo_t - pad, hi_t + pad)
        yp, yl = _time_ticks(lo_t - pad, hi_t + pad, max_ticks=6, two_line=False)
        axt.set_yticks(yp)
        axt.set_yticklabels(yl)
        axt.tick_params(labelsize=8)
        axt.text(0.02, 0.97, f"dashed line = fitted straight front\n{heading_words(fit).replace('rain moving ', '')}\n"
                             f"R2 {fit['r2']:.2f}, typical misfit {fit['rms_h']:.2f} h",
                 transform=axt.transAxes, va="top", ha="left", fontsize=7.5, color=INK2,
                 bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none", alpha=0.9))
        axt.legend(handles=[Line2D([], [], marker="o", ls="none", color=BLUE, alpha=0.7, ms=5, label="basin cells"),
                            Line2D([], [], marker="o", ls="none", color=MUTED, alpha=0.7, ms=5, label="other radar cells")],
                   loc="lower right", fontsize=8)
    else:
        axt.text(0.5, 0.5, "no usable direction\n(too few cells, or the half-way times\nbarely change across the area)",
                 transform=axt.transAxes, ha="center", va="center", color=INK2)
        axt.set_xticks([])
        axt.set_yticks([])
    axt.set_axisbelow(True)
    axt.grid(color=GRID)
    # ---- the halves ------------------------------------------------------------------------------------------------
    _half_panel(axe, res["halves"]["x"], s0, nh, "West vs east half of the basin", "mean rain over the half (mm/hr)")
    _half_panel(axn, res["halves"]["y"], s0, nh, "North vs south half of the basin", "mean rain over the half (mm/hr)")
    xp, xl = _time_ticks(s0, s0 + nh / 24.0, max_ticks=6, two_line=True)
    axn.set_xticks(xp)
    axn.set_xticklabels(xl)
    axe.tick_params(labelbottom=False)
    axn.set_xlabel("time (hour starting, MST)")
    fig.suptitle(title, fontsize=12, fontweight="bold", x=0.02, ha="left")
    fig.savefig(out_png, dpi=170, bbox_inches="tight")
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
    names = ("arrival_cells.csv", "arrival_fits.csv", "arrival_halves.csv", "arrival_run_info.txt") + (() if a.no_png else ("arrival_overview.png",))
    print("arrival_time_sycamore: " + ("DRY RUN (nothing will be written)" if a.dry_run else "run"))
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
    geo = cover = None
    rings = load_boundary_rings(a.boundary)
    if len(missing) < n:
        hdr, _ = read_ascii(next(f for f in files if f.exists()))
        geo = grid_geometry(hdr)
        print(f"  grid: {geo['nx']} x {geo['ny']} cells of {geo['cs']:g} m, x {geo['L']:.0f}-{geo['R']:.0f}, y {geo['B']:.0f}-{geo['T']:.0f}")
        cover = grid_cover(rings, geo)
        area = float(cover.sum()) * geo["cs"] ** 2 / 1e6
        print(f"  basin outline: {a.boundary}  ->  {area:.1f} km2 on the radar grid")
        if cover.sum() <= 0:
            die("the basin outline does not overlap the radar grid (wrong coordinate system?)")
        if a.expect_area_km2 and abs(area - a.expect_area_km2) > 0.03 * a.expect_area_km2:
            print(f"  WARNING: expected about {a.expect_area_km2:g} km2 (the DEM-derived basin is 424.6 km2); is this the right outline?")
    filled = read_filled_keys(manifest)
    if filled is None:
        print(f"  manifest: not found at {manifest} -> zero-filled hours cannot be flagged (none are expected in Jan-Mar 2017)")
        filled = set()
    in_window = [f"{lab:%b %d %H:%M}" for lab in labels if f"{lab:%Y-%m-%d %H:%M}" in filled]
    if filled:
        print(f"  manifest: {manifest}  ->  {len(filled)} zero-filled hours in the whole file, {len(in_window)} inside this window")
    print(f"  outputs ({out_dir}/):  " + "  ".join(f"{tag}_{s}" for s in names))
    if a.dry_run:
        print("\nDry run finished: nothing was computed or written.")
        return 0

    print(f"\nreading {n} grids ...")
    geo, stack = load_stack(files)
    n_nan = int(np.isnan(stack).sum())
    if n_nan:
        print(f"  note: {n_nan} NODATA cells were set to 0")
    rain = np.clip(np.nan_to_num(stack, nan=0.0), 0.0, None).astype("float64")
    area = float(cover.sum()) * geo["cs"] ** 2 / 1e6
    res = analyze(rain, geo, cover, a.min_total, a.max_lag)
    if res["n_cells"] == 0:
        die(f"no cell received at least {a.min_total:g} mm in this window; check --start/--end or lower --min-total")
    out_dir.mkdir(parents=True, exist_ok=True)
    p_cells, p_fits, p_halves, p_txt, p_png = (out_dir / f"{tag}_{s}" for s in
                                               ("arrival_cells.csv", "arrival_fits.csv", "arrival_halves.csv",
                                                "arrival_run_info.txt", "arrival_overview.png"))
    start_ts = labels[0] - HOUR
    rr, cc = np.indices(res["total"].shape)
    cells = pd.DataFrame({"row": rr.ravel(), "col": cc.ravel(), "x_m": geo["X"].ravel(), "y_m": geo["Y"].ravel(),
                          "basin_cover": cover.ravel().round(4), "total_mm": res["total"].ravel().round(3),
                          "t25_h": res["t"][0.25].ravel().round(3), "t50_h": res["t"][0.5].ravel().round(3),
                          "t75_h": res["t"][0.75].ravel().round(3)})
    cells["t50_mst"] = (start_ts + pd.to_timedelta(cells["t50_h"], unit="h")).dt.strftime("%Y-%m-%d %H:%M")
    cells.to_csv(p_cells, index=False)
    rows = []
    for (area_name, p), f in res["fits"].items():
        rows.append({"cells": area_name, "fraction": p, "fit_ok": f["ok"], "n_cells": f["n"],
                     "toward_bearing_deg": f["bearing_deg"], "speed_kmh": f["speed_kmh"], "r2": f["r2"], "misfit_h": f["rms_h"],
                     "across_area_h": f["across_h"]})
    pd.DataFrame(rows).round(4).to_csv(p_fits, index=False)
    hh = pd.DataFrame({"label_end_mst": labels.strftime("%Y-%m-%d %H:%M"), "interval_start_mst": (labels - HOUR).strftime("%Y-%m-%d %H:%M"),
                       "basin_mean_mm_hr": res["basin_mean"],
                       "west_mm_hr": res["halves"]["x"]["s1"], "east_mm_hr": res["halves"]["x"]["s2"],
                       "north_mm_hr": res["halves"]["y"]["s1"], "south_mm_hr": res["halves"]["y"]["s2"]})
    hh.round(4).to_csv(p_halves, index=False)
    lines = summarize(res, labels, a, area, in_window)
    settings = (f"settings: rain-dir {rain_dir} | boundary {a.boundary} | min-total {a.min_total:g} mm | max-lag {a.max_lag} h | "
                f"manifest {'used' if filled else 'not used / nothing flagged'}")
    p_txt.write_text("\n".join(lines + ["", settings]) + "\n")
    if not a.no_png:
        print("drawing the figure ...")
        make_figure(res, geo, rings, labels, a, p_png,
                    f"Lower Sycamore Wash: how the rain arrived, {start:%b %d %H:%M} to {end:%b %d %H:%M %Y} (MST)")
    print("\n" + "\n".join("  " + s for s in lines))
    print(f"\nwrote to {out_dir}/ :")
    for p in (p_cells, p_fits, p_halves, p_txt) + (() if a.no_png else (p_png,)):
        print(f"  {p.name}  ({p.stat().st_size / 1024:.0f} KB)")
    return 0


# ================================================================================================================
# self-test (made-up storms with known answers; touches none of your data)
# ================================================================================================================
def _plane_wave(X, Y, n, heading, speed, tau0, sigma, amp, x0=456000.0, y0=3744000.0):
    """Hourly depths (n, ny, nx) of a rain pulse (Gaussian in time, spread `sigma` h) that sweeps across the grid toward
    `heading` (degrees) at `speed` km/h. Cell's pulse is centred at tau = tau0 + distance_along_heading / speed, so its
    true half-way time is exactly tau. `amp` = each cell's window total in mm (number or grid). Returns (depths, tau)."""
    b = math.radians(heading)
    tau = tau0 + ((X - x0) * math.sin(b) + (Y - y0) * math.cos(b)) / 1000.0 / speed
    erf = np.vectorize(math.erf)
    k = np.arange(n + 1, dtype="float64")[:, None, None]
    cdf = 0.5 * (1.0 + erf((k - tau[None]) / (sigma * math.sqrt(2.0))))
    return amp * (cdf[1:] - cdf[:-1]), tau


def _angdiff(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def cmd_self_test(_a) -> int:
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("  PASS  " if cond else "  FAIL  ") + msg)
        ok = ok and bool(cond)

    tmp = Path(tempfile.mkdtemp(prefix="arrivaltest_"))
    try:
        print("arrival_time_sycamore SELF-TEST (made-up storms; your data is not touched)")
        # ---- 1. half-way time of one cell: exact cases ------------------------------------------------------------
        def one(series, frac, min_total=0.0):
            cum = cumulative(np.asarray(series, dtype=float)[:, None, None])
            return float(arrival_hours(cum, frac, min_total)[0, 0])
        check(abs(one([0, 0, 2, 2, 0, 0], 0.5) - 3.0) < 1e-9, "rain [0,0,2,2,0,0]: half the rain has fallen exactly 3.0 h after the start")
        check(abs(one([0, 4, 0, 0], 0.5) - 1.5) < 1e-9 and abs(one([0, 4, 0, 0], 0.25) - 1.25) < 1e-9 and abs(one([0, 4, 0, 0], 0.75) - 1.75) < 1e-9,
              "rain [0,4,0,0]: 25%, 50%, 75% are reached at 1.25, 1.5, 1.75 h (steady rain inside the hour)")
        check(abs(one([3, 1], 0.5) - 2.0 / 3.0) < 1e-9, "rain [3,1]: half (2 mm) is reached 0.667 h in")
        check(np.isnan(one([0, 0, 0], 0.5)) and np.isnan(one([1, 1], 0.5, min_total=5.0)),
              "a cell with no rain, or less than --min-total, gets no half-way time")
        # ---- 2. halves of the basin ---------------------------------------------------------------------------------
        cov = np.zeros((10, 8))
        cov[:, 3:6] = 1.0                                               # 3 columns x 10 rows
        w_w, w_e = half_weights(cov, "x")
        check(abs(w_w.sum() - 15.0) < 1e-9 and abs(w_e.sum() - 15.0) < 1e-9 and abs(w_w[:, 4].sum() - 5.0) < 1e-9 and w_w[:, 3].sum() == 10.0
              and np.allclose(w_w + w_e, cov) and w_w.min() >= 0 and w_e.min() >= 0,
              "the east-west split gives two halves of equal area; the middle column is shared 50/50")
        cov2 = np.zeros((10, 8))
        cov2[2:8, :] = 1.0
        w_n, w_s = half_weights(cov2, "y")
        check(abs(w_n.sum() - w_s.sum()) < 1e-9 and w_n[2:5].sum() == 24.0 and w_s[5:8].sum() == 24.0 and w_n[7].sum() == 0.0,
              "the north-south split puts row 0 (the NORTH row) in the north half")
        # ---- 3. straight-front fit: exact recovery ----------------------------------------------------------------------
        gx, gy = np.meshgrid(np.arange(0.0, 20.0), np.arange(0.0, 30.0))
        gx, gy = gx.ravel(), gy.ravel()
        wts = np.ones_like(gx)
        for brg, spd in ((225.0, 12.0), (90.0, 20.0), (0.0, 10.0), (180.0, 8.0), (315.0, 15.0)):
            tt = (gx * math.sin(math.radians(brg)) + gy * math.cos(math.radians(brg))) / spd + 4.0
            f = fit_plane(gx, gy, tt, wts)
            check(f["ok"] and _angdiff(f["bearing_deg"], brg) < 1e-4 and abs(f["speed_kmh"] - spd) < 1e-6 and f["r2"] > 0.999999,
                  f"a front moving toward {brg:.0f} deg at {spd:g} km/h is recovered exactly (got {f['bearing_deg']:.3f} deg, {f['speed_kmh']:.4f} km/h)")
        rng = np.random.default_rng(3)
        fr = fit_plane(gx, gy, rng.uniform(0, 10, gx.size), wts)
        check(fr["ok"] and fr["r2"] < 0.3 and fit_verdict(fr).startswith("WEAK"), f"random timing gives a weak fit (R2 {fr['r2']:.3f}) and a WEAK warning")
        check(not fit_plane(gx[:5], gy[:5], np.arange(5.0), np.ones(5))["ok"], "fewer than 10 cells gives no fit")
        # ---- 4. cross-correlation lag ----------------------------------------------------------------------------------
        bump = np.zeros(24)
        bump[4:7] = (1, 3, 1)
        late = np.roll(bump, 3)
        lag_a, _ = xcorr_lag(bump, late, 6)
        lag_b, _ = xcorr_lag(late, bump, 6)
        check(lag_a == 3 and lag_b == -3, "the second series 3 h later gives lag +3 (and -3 when swapped)")
        check(np.isnan(xcorr_lag(np.zeros(24), bump, 6)[0]), "a flat series gives no lag")
        # ---- 5. synthetic storm sweeping the grid ------------------------------------------------------------------------
        geo0 = grid_geometry({"ncols": 32, "nrows": 49, "xllcorner": 440000, "yllcorner": 3720000, "cellsize": 1000})
        shp = tmp / "test_basin.shp"
        _write_test_boundary(shp, 448000, 3726000, 468000, 3758000)           # 20 km x 32 km = 640 km2
        rings = load_boundary_rings(shp)
        cover = grid_cover(rings, geo0)
        check(abs(cover.sum() - 640.0) < 2.0, f"the 20 x 32 km test basin covers {cover.sum():.1f} cells (about 640)")
        n = 24
        rain, tau = _plane_wave(geo0["X"], geo0["Y"], n, 225.0, 12.0, 12.0, 1.5, 30.0)
        res = analyze(rain, geo0, cover, 5.0, 6)
        err = float(np.nanmax(np.abs(res["t"][0.5] - tau)))
        check(err < 0.1, f"half-way times match the true centre times of the pulses (worst error {err:.3f} h, limit 0.1 h)")
        f = res["fits"][("basin", 0.5)]
        check(has_direction(f) and _angdiff(f["bearing_deg"], 225.0) < 1.0 and abs(f["speed_kmh"] - 12.0) < 0.4 and f["r2"] > 0.99,
              f"storm sweeping toward the SW at 12 km/h: recovered {f['bearing_deg']:.1f} deg, {f['speed_kmh']:.2f} km/h, R2 {f['r2']:.4f}")
        check(all(has_direction(ff) and _angdiff(ff["bearing_deg"], 225.0) < 2.0 for ff in res["fits"].values()) and bearing_spread(res["fits"]) < 2.0,
              f"all 6 fits (basin/grid x 25/50/75%) agree on the direction (largest departure {bearing_spread(res['fits']):.2f} deg)")
        bx, by = geo0["X"][res["inside"]] / 1000.0, geo0["Y"][res["inside"]] / 1000.0
        s_true = (bx - bx.mean()) * math.sin(math.radians(225.0)) + (by - by.mean()) * math.cos(math.radians(225.0))
        want_across = (float(np.percentile(s_true, 95)) - float(np.percentile(s_true, 5))) / 12.0
        check(abs(f["across_h"] - want_across) < 0.15 and want_across > 1.0,
              f"fitted time difference across the basin is {f['across_h']:.2f} h (expected {want_across:.2f} h from the true speed)")
        slope = math.sin(math.radians(225.0)) / 12.0                           # hours per km along x (negative: east is earlier)
        hx, hy = res["halves"]["x"], res["halves"]["y"]
        want_x = slope * 10.0                                                  # half centres are 10 km apart east-west
        want_y = (math.cos(math.radians(225.0)) / 12.0) * -16.0                # south half's centre is 16 km south of the north half's
        check(abs(hx["diff_h"] - want_x) < 0.15 and hx["diff_h"] < 0,
              f"east half arrives earlier than west by {-hx['diff_h']:.2f} h (expected {-want_x:.2f} h; sign: east minus west is negative)")
        check(abs(hy["diff_h"] - want_y) < 0.15 and hy["diff_h"] > 0,
              f"south half arrives later than north by {hy['diff_h']:.2f} h (expected {want_y:.2f} h)")
        check(abs(hx["lag_h"] - want_x) <= 1.0 and abs(hy["lag_h"] - want_y) <= 1.0,
              f"cross-correlation lags agree within one hour (east-west {hx['lag_h']:+.0f} h, north-south {hy['lag_h']:+.0f} h)")
        check(res["first_share"] < 0.005 and res["last_share"] < 0.005, "a storm that fits inside the window has almost no rain in its first/last hour")
        # a different storm: toward the east at 20 km/h, rain amount varying across the grid
        amp = 25.0 + 15.0 * (geo0["X"] - 440000.0) / 32000.0
        rain2, tau2 = _plane_wave(geo0["X"], geo0["Y"], n, 90.0, 20.0, 12.0, 1.0, amp)
        res2 = analyze(rain2, geo0, cover, 5.0, 6)
        f2 = res2["fits"][("basin", 0.5)]
        check(has_direction(f2) and _angdiff(f2["bearing_deg"], 90.0) < 1.0 and abs(f2["speed_kmh"] - 20.0) < 0.8,
              f"storm sweeping east at 20 km/h with uneven rain depth: recovered {f2['bearing_deg']:.1f} deg, {f2['speed_kmh']:.2f} km/h")
        check(res2["halves"]["x"]["diff_h"] > 0.3 and res2["halves"]["x"]["lag_h"] >= 0, "for an eastward storm the east half is the LATER one (positive east minus west)")
        # window that cuts the storm
        cut = analyze(rain[3:15], geo0, cover, 5.0, 6)
        check(cut["first_share"] > EDGE_SHARE_WARN or cut["last_share"] > EDGE_SHARE_WARN, "a window that cuts the storm in two trips the first/last-hour warning")
        # rain that reaches the whole basin at once must not report a direction or a GOOD FIT
        rain_fast, _ = _plane_wave(geo0["X"], geo0["Y"], n, 225.0, 400.0, 12.0, 1.5, 30.0)
        res_fast = analyze(rain_fast, geo0, cover, 5.0, 6)
        txt = "\n".join(summarize(res_fast, parse_window("2017-01-20 00:00", "2017-01-21 00:00")[2],
                                 build_parser().parse_args(["run"]), 640.0, []))
        check(not has_direction(res_fast["fits"][("basin", 0.5)]) and "at once" in txt and "GOOD FIT" not in txt
              and "toward the" not in txt and "none" in txt,
              "rain reaching the whole basin at once (400 km/h) reports no direction and no GOOD FIT, and says 'at once'")
        # ---- 6. a full run through the files ---------------------------------------------------------------------------------
        rain_dir = tmp / "wk" / "rain"
        rain_dir.mkdir(parents=True)
        s, e, labels = parse_window("2017-01-20 00:00", "2017-01-21 00:00")
        for k, lab in enumerate(labels):
            _write_asc(file_for_label(rain_dir, "radar_", lab), rain[k])
        (tmp / "wk" / "manifest.csv").write_text("label_local,status\n" + "\n".join(
            f"{lab:%Y-%m-%d %H:%M},{'filled_zero_missing_hour' if k == 5 else 'ok'}" for k, lab in enumerate(labels)) + "\n")
        out1 = tmp / "out1"
        base = ["--rain-dir", str(rain_dir), "--boundary", str(shp), "--start", "2017-01-20 00:00", "--end", "2017-01-21 00:00",
                "--out-dir", str(out1), "--tag", "t", "--expect-area-km2", "640"]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = main(["run", "--dry-run"] + base)
        check(rc == 0 and not out1.exists() and "files present: 24 of 24" in buf.getvalue(), "--dry-run finds 24 of 24 files and writes nothing")
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                main(["run", "--no-png", "--rain-dir", str(rain_dir), "--boundary", str(shp), "--start", "2017-01-20 00:00",
                      "--end", "2017-01-21 01:00", "--out-dir", str(out1)])
            stopped = False
        except SystemExit as ex:
            stopped = "missing" in str(ex)
        check(stopped, "a window that needs a file that is not there stops with a 'missing' error")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = main(["run"] + base)
        names = [f"t_{x}" for x in ("arrival_overview.png", "arrival_cells.csv", "arrival_fits.csv", "arrival_halves.csv", "arrival_run_info.txt")]
        check(rc == 0 and all((out1 / nme).exists() and (out1 / nme).stat().st_size > 0 for nme in names), "a full run writes the PNG, three CSVs and the run-info text")
        info = (out1 / "t_arrival_run_info.txt").read_text()
        check("toward the SW" in info and "arriving from the NE" in info and "GOOD FIT" in info and "STABLE" in info,
              "run_info.txt says the rain moved toward the SW (arriving from the NE) with a GOOD FIT and a STABLE direction")
        check("zero-filled" in info and "Jan 20 06:00" in info, "the zero-filled hour in the manifest is listed as a warning")
        cells = pd.read_csv(out1 / "t_arrival_cells.csv")
        check(len(cells) == 32 * 49 and {"row", "col", "x_m", "y_m", "basin_cover", "total_mm", "t25_h", "t50_h", "t75_h", "t50_mst"} <= set(cells.columns),
              "the cells CSV has one row per radar cell (1568) and the documented columns")
        c0 = cells[(cells["row"] == 0) & (cells["col"] == 0)].iloc[0]
        check(abs(c0["y_m"] - 3768500) < 1e-6 and abs(c0["x_m"] - 440500) < 1e-6, "row 0 / col 0 is the NORTH-WEST cell (x 440500, y 3768500)")
        fits = pd.read_csv(out1 / "t_arrival_fits.csv")
        halves = pd.read_csv(out1 / "t_arrival_halves.csv")
        check(len(fits) == 6 and len(halves) == 24 and list(halves.columns[-4:]) == ["west_mm_hr", "east_mm_hr", "north_mm_hr", "south_mm_hr"],
              "the fits CSV has 6 rows and the halves CSV has 24 hourly rows")
        # a run whose window cuts the storm prints the warning
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = main(["run", "--no-png", "--rain-dir", str(rain_dir), "--boundary", str(shp), "--start", "2017-01-20 03:00",
                       "--end", "2017-01-20 15:00", "--out-dir", str(out1), "--tag", "cut", "--expect-area-km2", "640"])
        check("WARNING: more than 3% falls in the first or last hour" in buf.getvalue(), "a window that cuts the storm prints the widen-the-window warning")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\nSELF-TEST " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


# ================================================================================================================
# command line
# ================================================================================================================
def build_parser():
    p = argparse.ArgumentParser(description="When the rain arrived at each part of the Lower Sycamore basin, from hourly MRMS "
                                            "grids (see the text at the top of this file).",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="make the map, figure and tables for one storm window")
    r.add_argument("--rain-dir", default=DEFAULT_RAIN_DIR, help=f"folder with radar_MMDDYYYYHH.asc (default {DEFAULT_RAIN_DIR})")
    r.add_argument("--prefix", default="radar_", help="file-name start before MMDDYYYYHH (default radar_)")
    r.add_argument("--boundary", default=DEFAULT_BOUNDARY, help=f"basin outline shapefile in EPSG:26912 (default {DEFAULT_BOUNDARY})")
    r.add_argument("--manifest", default=None, help="manifest.csv from the radar build (default: one folder above --rain-dir)")
    r.add_argument("--start", default=DEFAULT_START, help=f"start of the first rain hour, MST, whole hour (default {DEFAULT_START})")
    r.add_argument("--end", default=DEFAULT_END, help=f"end of the last rain hour, MST, whole hour (default {DEFAULT_END})")
    r.add_argument("--tag", default=None, help="start of the output file names (default: the two dates)")
    r.add_argument("--out-dir", default="storm_movement_out", help="folder for the outputs (default storm_movement_out)")
    r.add_argument("--min-total", type=float, default=5.0, help="leave out cells with less than this many mm in the window (default 5)")
    r.add_argument("--max-lag", type=int, default=6, help="largest shift, in hours, tried by the cross-correlation check (default 6)")
    r.add_argument("--expect-area-km2", type=float, default=424.6, help="warn if the outline's area on the grid differs by >3%% (0 = off)")
    r.add_argument("--no-png", action="store_true", help="skip the figure (the CSVs and text are quick)")
    r.add_argument("--dry-run", action="store_true", help="check the files and settings, write nothing")
    sub.add_parser("self-test", help="check the script on made-up storms")
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    return cmd_run(a) if a.cmd == "run" else cmd_self_test(a)


if __name__ == "__main__":
    sys.exit(main())
