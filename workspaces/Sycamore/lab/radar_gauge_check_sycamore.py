#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
radar_gauge_check_sycamore.py   (v1.0, 2026-10-27)
==================================================

STAGE 0 of the Lower Sycamore calibration design: does the hourly MRMS radar rain that the model is given
agree with the rain gauges inside the watershed?

The design doc says (Stage 0, "Radar check"):
    For Mt. Ord (2000), Sycamore Cr Upper (82200) and Iron Dike (83000), get hourly or 15-minute totals.
    Take the radar value at each gauge's 1 km cell and its 3x3 neighbours, using end-of-hour labels.
    Compare event totals and hourly timing, and compute radar / gauge per storm.
    A consistent ratio across storms makes rain scale a calibrated knob with that range.
    A ratio that swings storm to storm keeps it a probe.

WHAT THIS SCRIPT DOES
---------------------
For every gauge you list and every storm window you list it:
  1. reads the radar grids the model uses (radar_MMDDYYYYHH.asc, 1 km, mm/hr, EPSG:26912);
  2. turns the gauge record (any time step finer than or equal to one hour) into hourly depths;
  3. takes the radar value in the gauge's own 1 km cell, and in the 3x3 block of cells around it;
  4. adds up gauge and radar over the SAME hours (an hour where either one is missing is left out of both);
  5. prints and saves: storm totals, radar / gauge ratios, peak-hour and centre-of-rain timing differences,
     a lag test (does the radar pattern line up with the gauge at lag 0, or one hour early or late?),
     a day-by-day table (to help you adjust the storm windows), and pictures;
  6. ends with a PROPOSED reading of the rain-scale question (consistent ratio = calibrated knob,
     swinging ratio = probe). The thresholds for "consistent" and "swings" are proposals, not rules.

WHAT YOU HAVE TO SUPPLY (this script does not have them)
---------------------------------------------------------
  * the gauge files and where each gauge is: x, y in EPSG:26912 metres, OR lon, lat in decimal degrees
    (a gauges CSV, see `template`);
  * optionally edited storm windows (a storms CSV). The built-in defaults come from the windows in
    storm_rain_cell.py and the compare cell. Edit them after reading the day-by-day table.

COMMANDS
--------
  python radar_gauge_check_sycamore.py template                 # writes gauges_template.csv and storms_template.csv
  python radar_gauge_check_sycamore.py run --dry-run  --radar-dir <dir> --gauges gauges.csv
  python radar_gauge_check_sycamore.py run            --radar-dir <dir> --gauges gauges.csv [--storms storms.csv]
  python radar_gauge_check_sycamore.py self-test                # builds synthetic data with known answers and checks them

TIME CONVENTIONS (read this, it is the easiest place to be off by an hour)
--------------------------------------------------------------------------
  * All clocks are Arizona local time (UTC-7, no daylight saving).
  * The Sycamore radar files are labeled with the END of the rain hour (default --radar-label end): the file
    radar_0121201713.asc holds the rain that fell 12:00-13:00. Use --radar-label start if a set is labeled
    by the start of the hour.
  * Inside this script every hour is labeled by when it STARTS (the same as storm_rain_cell.py, which moves
    the .mrf labels back one hour). A storm window [start, end) contains the hours that START inside it.
    Printed "hour ending" times add one hour back.
  * A gauge record time-stamped HH:00 exactly counts in the hour ENDING at HH:00.
  * Gauge file clocks: --gauge-utc-offset (or a utc_offset column in the gauges CSV). -7 means the gauge file
    is already on Arizona time. 0 means the gauge file is in UTC.
  * LAG SIGN: lag +1 means the radar value one hour LATER matches the gauge now, i.e. the radar is one hour LATE.

REQUIRES: numpy, pandas, matplotlib (pictures only). No pyproj needed (the lon/lat to UTM 12N conversion is built in).

WHAT IS NOT CHECKED
-------------------
  * Gauge catch error (wind, tipping-bucket under-count in heavy rain) is not corrected.
  * A gauge sees a point; a radar cell is a 1 km average. For a patchy convective storm they can differ a lot
    even if both are right. Three gauges cannot separate that from a radar bias.
  * MRMS gauge-corrected products use many gauges. Whether these FCDMC gauges are among them is not known here,
    so the check may not be fully independent.
  * Not yet run on the real radar files or the real gauge files (written without access to them).
    `self-test` checks the logic on synthetic data with known answers.
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
import tempfile
import textwrap
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = "1.0"
HOUR = pd.Timedelta(hours=1)

# ----------------------------------------------------------------------------------------------------------------
# DEFAULT STORM WINDOWS   (name, first hour included, end not included; local clock; an hour belongs to the window
# that contains the time it STARTS). Taken from the windows in storm_rain_cell.py and the compare cell
# (flood week Jan 18-25, second Jan storm from Jan 23 12:00, Feb 18-20 and Feb 27-Mar 2 "first days").
# These are rain windows, so they leave out the long recession tails that the runoff windows include.
# Edit with --storms after reading the day-by-day table.
# ----------------------------------------------------------------------------------------------------------------
DEFAULT_STORMS = [
    ("Jan 18-23 flood storm",       "2017-01-18 00:00", "2017-01-23 12:00"),
    ("Jan 23-25 second Jan storm",  "2017-01-23 12:00", "2017-01-26 00:00"),
    ("Feb 18-20 storm",             "2017-02-18 00:00", "2017-02-21 00:00"),
    ("Feb 27-Mar 2 storm",          "2017-02-27 00:00", "2017-03-03 00:00"),
]

# Palette: first slots of the validated default categorical palette (blue, orange, aqua, yellow, magenta, green,
# violet, red); markers differ too, so a gauge is never identified by color alone.
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SERIES_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
INK = "#0b0b0b"
INK_2 = "#52514e"
GAUGE_BAR = "#8a8984"          # neutral gray for the observation
RADAR_BAND = "#9ec5f4"         # light end of the blue ramp, for the 3x3 min-max band
GRID_COLOR = "#d9d8d3"

LAGS = list(range(-3, 4))


def die(msg: str, code: int = 2) -> None:
    print(f"\nERROR: {msg}\n", file=sys.stderr)
    sys.exit(code)


# ================================================================================================================
# small helpers
# ================================================================================================================
def lonlat_to_utm(lon: float, lat: float, zone: int = 12):
    """WGS84 lon/lat (decimal degrees) to UTM easting/northing in metres, northern hemisphere.
    EPSG:26912 is NAD83 / UTM zone 12N; NAD83 and WGS84 differ by about a metre here, which is far below
    the 1 km size of a radar cell. Standard transverse-Mercator series (Snyder 1987)."""
    a, f, k0 = 6378137.0, 1 / 298.257223563, 0.9996
    e2 = f * (2 - f)
    ep2 = e2 / (1 - e2)
    lon0 = math.radians(-183 + 6 * zone)
    phi, lam = math.radians(lat), math.radians(lon)
    n = a / math.sqrt(1 - e2 * math.sin(phi) ** 2)
    t, c, aa = math.tan(phi) ** 2, ep2 * math.cos(phi) ** 2, math.cos(phi) * (lam - lon0)
    e4, e6 = e2 * e2, e2 * e2 * e2
    m = a * ((1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * phi
             - (3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * math.sin(2 * phi)
             + (15 * e4 / 256 + 45 * e6 / 1024) * math.sin(4 * phi)
             - (35 * e6 / 3072) * math.sin(6 * phi))
    x = k0 * n * (aa + (1 - t + c) * aa ** 3 / 6 + (5 - 18 * t + t * t + 72 * c - 58 * ep2) * aa ** 5 / 120) + 500000.0
    y = k0 * (m + n * math.tan(phi) * (aa * aa / 2 + (5 - t + 9 * c + 4 * c * c) * aa ** 4 / 24
                                       + (61 - 58 * t + t * t + 600 * c - 330 * ep2) * aa ** 6 / 720))
    return x, y


def read_ascii(path: Path):
    """ESRI ASCII grid -> (header dict with lower-case keys, 2-D array, rows top to bottom). NODATA becomes NaN."""
    with open(path) as fh:
        hdr = {}
        for _ in range(6):
            parts = fh.readline().split()
            if len(parts) != 2:
                raise ValueError("malformed header line")
            hdr[parts[0].lower()] = float(parts[1])
        body = fh.read().split()
    n = int(hdr["ncols"]) * int(hdr["nrows"])
    if len(body) != n:
        raise ValueError(f"expected {n} values, found {len(body)}")
    data = np.array(body, dtype="float64").reshape(int(hdr["nrows"]), int(hdr["ncols"]))
    data[data == hdr["nodata_value"]] = np.nan
    return hdr, data


def hour_starts(first: pd.Timestamp, end_excl: pd.Timestamp) -> pd.DatetimeIndex:
    return pd.date_range(first, end_excl - HOUR, freq="60min")


def fmt_end(ts) -> str:
    """Printable 'hour ending' time for an hour that STARTS at ts."""
    return "n.a." if pd.isna(ts) else f"{ts + HOUR:%b %d %H:%M}"


def safe_div(a, b):
    return a / b if (b is not None and np.isfinite(b) and b > 0 and np.isfinite(a)) else np.nan


def nanagg(arr, fn, axis):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return fn(arr, axis=axis)


# ================================================================================================================
# reading the configuration files
# ================================================================================================================
def read_storms(path: str | None):
    if not path:
        return [(n, pd.Timestamp(a), pd.Timestamp(b)) for n, a, b in DEFAULT_STORMS]
    df = pd.read_csv(path, comment="#", skipinitialspace=True, dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]
    for need in ("name", "start", "end"):
        if need not in df.columns:
            die(f"storms file {path} needs the columns name,start,end (found {list(df.columns)})")
    out = []
    for r in df.itertuples():
        a, b = pd.Timestamp(r.start), pd.Timestamp(r.end)
        if b <= a:
            die(f"storm {r.name!r}: end {b} is not after start {a}")
        out.append((str(r.name).strip(), a, b))
    if not out:
        die(f"storms file {path} has no rows")
    return out


GAUGE_DEFAULTS = {"units": "mm", "value_type": "increment", "time_col": "", "value_col": "", "scale": "1",
                  "no_record": "zero", "valid_from": "", "valid_to": "", "utc_offset": "",
                  "x": "", "y": "", "lon": "", "lat": ""}


def read_gauge_config(path: str, default_utc_offset: float, base_dir: Path | None = None):
    """Read the gauges CSV. Columns: name, file, then x,y (EPSG:26912 m) OR lon,lat, then optional:
    units (mm|in|cm), value_type (increment|cumulative), time_col, value_col, scale (multiplier applied first, e.g.
    the depth of one tip), no_record (zero|missing: what an hour with NO record means), valid_from, valid_to
    (local clock: the period the file really covers), utc_offset (clock of this gauge file)."""
    df = pd.read_csv(path, comment="#", skipinitialspace=True, dtype=str, keep_default_na=False)
    df.columns = [c.strip().lower() for c in df.columns]
    for need in ("name", "file"):
        if need not in df.columns:
            die(f"gauges file {path} needs the columns name and file (found {list(df.columns)})")
    base = base_dir or Path(path).resolve().parent
    gauges = []
    for i, r in df.iterrows():
        g = {k: str(r.get(k, v)).strip() if k in df.columns else v for k, v in GAUGE_DEFAULTS.items()}
        g["name"], g["file"] = str(r["name"]).strip(), str(r["file"]).strip()
        for k in GAUGE_DEFAULTS:                      # blank cell in the CSV -> default
            if g[k] == "":
                g[k] = GAUGE_DEFAULTS[k]
        if g["file"] == "":
            die(f"gauge {g['name']!r} has no file in {path}")
        p = Path(g["file"])
        g["path"] = p if p.is_absolute() else (base / p if (base / p).exists() else Path.cwd() / p)
        if g["x"] and g["y"]:
            g["x"], g["y"] = float(g["x"]), float(g["y"])
            g["location_from"] = "x,y"
        elif g["lon"] and g["lat"]:
            g["x"], g["y"] = lonlat_to_utm(float(g["lon"]), float(g["lat"]))
            g["location_from"] = "lon,lat"
        else:
            die(f"gauge {g['name']!r} has no location: fill x and y (EPSG:26912 m) or lon and lat in {path}")
        if g["units"] not in ("mm", "in", "cm"):
            die(f"gauge {g['name']!r}: units must be mm, in or cm, not {g['units']!r}")
        if g["value_type"] not in ("increment", "cumulative"):
            die(f"gauge {g['name']!r}: value_type must be increment or cumulative, not {g['value_type']!r}")
        if g["no_record"] not in ("zero", "missing"):
            die(f"gauge {g['name']!r}: no_record must be zero or missing, not {g['no_record']!r}")
        g["scale"] = float(g["scale"])
        g["utc_offset"] = float(g["utc_offset"]) if g["utc_offset"] != "" else float(default_utc_offset)
        gauges.append(g)
    if not gauges:
        die(f"gauges file {path} has no rows")
    return gauges


# ================================================================================================================
# gauge data -> hourly depth (mm), labeled by hour START, on the model clock
# ================================================================================================================
def load_gauge_hourly(g: dict, hours: pd.DatetimeIndex, model_utc_offset: float):
    """Return (hourly Series on `hours`, info dict). Hours the file does not cover are NaN."""
    if not g["path"].exists():
        die(f"gauge file not found for {g['name']!r}: {g['path']}")
    df = pd.read_csv(g["path"], comment="#", skipinitialspace=True)
    tcol = g["time_col"] or df.columns[0]
    vcol = g["value_col"] or df.columns[1]
    for c in (tcol, vcol):
        if c not in df.columns:
            die(f"gauge {g['name']!r}: column {c!r} is not in {g['path']} (columns: {list(df.columns)})")
    ts = pd.to_datetime(df[tcol], errors="coerce")
    if ts.dtype == object:                                  # mixed offsets in the strings
        ts = pd.to_datetime(df[tcol], errors="coerce", utc=True)
    if getattr(ts.dt, "tz", None) is not None:             # time zone written in the file: trust it
        ts = ts.dt.tz_convert("UTC").dt.tz_localize(None) + pd.Timedelta(hours=model_utc_offset)
    else:                                                   # naive: shift from the gauge's clock to the model clock
        ts = ts + pd.Timedelta(hours=model_utc_offset - g["utc_offset"])
    vals = pd.to_numeric(df[vcol], errors="coerce") * g["scale"] * {"mm": 1.0, "in": 25.4, "cm": 10.0}[g["units"]]
    ok = ts.notna() & vals.notna()
    s = pd.Series(vals[ok].to_numpy(float), index=pd.DatetimeIndex(ts[ok])).sort_index()
    n_records = len(s)
    if g["value_type"] == "cumulative":
        inc = s.diff()                                      # the FIRST record is only a baseline reading: it has no increment
        reset = inc < 0                                     # counter reset: the new value is the rain since the reset
        inc[reset] = s[reset]
        s, n_resets = inc.dropna(), int(reset.sum())
    else:
        n_resets = 0
    info = {"n_records": n_records, "n_resets": n_resets,
            "first_record": s.index.min() if len(s) else pd.NaT, "last_record": s.index.max() if len(s) else pd.NaT,
            "n_negative": int((s < 0).sum()) if len(s) else 0}
    if not len(s):
        return pd.Series(np.nan, index=hours), info
    s = s.clip(lower=0)
    # hour ending at HH:00 holds records in (HH-1:00, HH:00]; label = START of that hour
    h = s.resample("60min", closed="right", label="left").sum(min_count=1)
    vf = pd.Timestamp(g["valid_from"]) if g["valid_from"] else h.index[0]
    vt = pd.Timestamp(g["valid_to"]) if g["valid_to"] else h.index[-1] + HOUR
    covered = hour_starts(vf.floor("60min"), vt.ceil("60min"))
    h = h.reindex(covered)
    if g["no_record"] == "zero":
        h = h.fillna(0.0)
    out = h.reindex(hours)
    info["valid_from"], info["valid_to"] = covered[0], covered[-1] + HOUR
    return out, info


# ================================================================================================================
# radar
# ================================================================================================================
def load_radar(radar_dir: Path, prefix: str, ext: str, hours: pd.DatetimeIndex, label_kind: str):
    """Stack of hourly grids, shape (len(hours), nrows, ncols), mm/hr (= mm per hour). Missing files give NaN."""
    stack, hdr0, missing, bad = None, None, [], []
    for i, t0 in enumerate(hours):
        t_label = t0 + HOUR if label_kind == "end" else t0
        path = radar_dir / f"{prefix}{t_label:%m%d%Y%H}.{ext}"
        if not path.exists():
            missing.append(t0)
            continue
        try:
            hdr, data = read_ascii(path)
        except (ValueError, OSError) as exc:
            bad.append((path.name, str(exc)))
            continue
        if hdr0 is None:
            hdr0 = hdr
            stack = np.full((len(hours), data.shape[0], data.shape[1]), np.nan, dtype="float32")
        elif hdr != hdr0:
            bad.append((path.name, "header differs from the first file"))
            continue
        stack[i] = data
    if hdr0 is None:
        die(f"no readable radar files like {prefix}MMDDYYYYHH.{ext} in {radar_dir} for "
            f"{hours[0]:%Y-%m-%d %H:%M} to {hours[-1] + HOUR:%Y-%m-%d %H:%M} "
            f"(label kind: {label_kind}; first name tried: {prefix}{(hours[0] + (HOUR if label_kind == 'end' else 0 * HOUR)):%m%d%Y%H}.{ext})")
    return stack, hdr0, missing, bad


def locate_cell(x: float, y: float, hdr: dict):
    cs = hdr["cellsize"]
    col = int(math.floor((x - hdr["xllcorner"]) / cs))
    row = int(hdr["nrows"]) - 1 - int(math.floor((y - hdr["yllcorner"]) / cs))
    inside = 0 <= col < int(hdr["ncols"]) and 0 <= row < int(hdr["nrows"])
    return row, col, inside


# ================================================================================================================
# the comparison
# ================================================================================================================
class GaugeSeries:
    """Everything about one gauge on the hourly index: gauge depth, radar cell, radar 3x3 block."""

    def __init__(self, g, info, hours, gauge_hr, stack, hdr):
        self.name, self.cfg, self.info = g["name"], g, info
        self.hours = hours
        row, col, inside = locate_cell(g["x"], g["y"], hdr)
        self.inside, self.row, self.col = inside, row, col
        nrows, ncols = int(hdr["nrows"]), int(hdr["ncols"])
        if inside:
            r0, r1, c0, c1 = max(0, row - 1), min(nrows, row + 2), max(0, col - 1), min(ncols, col + 2)
            block = stack[:, r0:r1, c0:c1].reshape(len(hours), -1).astype("float64")
            self.n_nb = block.shape[1]
            self.edge = self.n_nb < 9
            self.cell = pd.Series(stack[:, row, col].astype("float64"), index=hours)
            self.nb = block                                                    # (T, k) every neighbour hour by hour
            self.nb_mean = pd.Series(nanagg(block, np.nanmean, 1), index=hours)
            self.nb_min = pd.Series(nanagg(block, np.nanmin, 1), index=hours)
            self.nb_max = pd.Series(nanagg(block, np.nanmax, 1), index=hours)
        else:
            self.n_nb, self.edge = 0, False
            nan = pd.Series(np.nan, index=hours)
            self.cell = self.nb_mean = self.nb_min = self.nb_max = nan
            self.nb = np.full((len(hours), 1), np.nan)
        self.gauge = gauge_hr


def lag_table(gser: pd.Series, rser: pd.Series, mask: np.ndarray, min_pairs: int = 12):
    """Correlation of the gauge at hour t with the radar at hour t+k, for k in LAGS, over the hours in `mask`.
    k > 0 best = the radar is k hours LATE."""
    rows = []
    for k in LAGS:
        shifted = rser.copy()
        shifted.index = shifted.index - pd.Timedelta(hours=k)         # the radar value at t+k now sits at t
        shifted = shifted.reindex(gser.index)
        ok = mask & gser.notna().to_numpy() & shifted.notna().to_numpy()
        n = int(ok.sum())
        r = np.nan
        if n >= min_pairs:
            a, b = gser.to_numpy()[ok], shifted.to_numpy()[ok]
            if a.std() > 0 and b.std() > 0:
                r = float(np.corrcoef(a, b)[0, 1])
        rows.append((k, r, n))
    return rows


def best_lag(rows):
    valid = [(k, r) for k, r, n in rows if np.isfinite(r)]
    if not valid:
        return np.nan, np.nan, np.nan
    best = max(valid, key=lambda kr: (round(kr[1], 6), -abs(kr[0])))
    r0 = next((r for k, r in valid if k == 0), np.nan)
    return best[0], best[1], r0


def one_storm_row(G: GaugeSeries, name, s0, s1, hours, grid_mean, model_rain, min_gauge_mm):
    m = ((hours >= s0) & (hours < s1))
    n_hours = int(m.sum())
    g, c, nbm = G.gauge[m], G.cell[m], G.nb_mean[m]
    ok = (g.notna() & c.notna()).to_numpy()
    n_ok = int(ok.sum())
    row = {"gauge": G.name, "storm": name, "window_start": s0, "window_end": s1,
           "hours_in_window": n_hours, "hours_used": n_ok,
           "coverage_pct": 100.0 * n_ok / n_hours if n_hours else np.nan}
    nan_keys = ["gauge_mm", "radar_cell_mm", "radar_3x3_mean_mm", "radar_3x3_min_cell_mm", "radar_3x3_max_cell_mm",
                "ratio_cell", "ratio_3x3_mean", "ratio_3x3_min_cell", "ratio_3x3_max_cell",
                "gauge_peak_hour_mm", "radar_cell_peak_hour_mm", "gauge_peak_hour_ending", "radar_cell_peak_hour_ending",
                "peak_hour_diff_h", "centre_of_rain_diff_h", "best_lag_h", "r_at_lag0", "r_at_best_lag",
                "grid_mean_radar_mm", "model_basin_mm"]
    for k in nan_keys:
        row[k] = np.nan
    if n_ok == 0:
        row["note"] = "no hour with both gauge and radar"
        return row
    gv, cv = g.to_numpy()[ok], c.to_numpy()[ok]
    idx = g.index[ok]
    gt, ct = float(gv.sum()), float(cv.sum())
    nbt = nbm.to_numpy()[ok]
    nbt = float(np.nansum(nbt))
    per_nb = np.nansum(G.nb[m.to_numpy() if hasattr(m, "to_numpy") else m][ok], axis=0)
    row.update(gauge_mm=gt, radar_cell_mm=ct, radar_3x3_mean_mm=nbt,
               radar_3x3_min_cell_mm=float(per_nb.min()), radar_3x3_max_cell_mm=float(per_nb.max()))
    enough = gt >= min_gauge_mm
    if enough:
        row.update(ratio_cell=ct / gt, ratio_3x3_mean=nbt / gt, ratio_3x3_min_cell=float(per_nb.min()) / gt,
                   ratio_3x3_max_cell=float(per_nb.max()) / gt)
    row["note"] = "" if enough else f"gauge total below {min_gauge_mm:g} mm: ratio not computed"
    if G.edge:
        row["note"] = (row["note"] + "; " if row["note"] else "") + "gauge cell is on the grid edge: fewer than 9 cells"
    if gt > 0:
        ig, ic = int(np.argmax(gv)), int(np.argmax(cv))
        row.update(gauge_peak_hour_mm=float(gv[ig]), radar_cell_peak_hour_mm=float(cv[ic]),
                   gauge_peak_hour_ending=idx[ig] + HOUR, radar_cell_peak_hour_ending=idx[ic] + HOUR,
                   peak_hour_diff_h=(idx[ic] - idx[ig]) / HOUR)
        if ct > 0:
            t = np.asarray((idx - s0) / HOUR, dtype=float) + 0.5       # middle of each hour
            row["centre_of_rain_diff_h"] = float((t * cv).sum() / cv.sum() - (t * gv).sum() / gv.sum())
    k, rb, r0 = best_lag(lag_table(G.gauge, G.cell, np.asarray(m)))
    row.update(best_lag_h=k, r_at_best_lag=rb, r_at_lag0=r0)
    gm = grid_mean[m].to_numpy()[ok]
    row["grid_mean_radar_mm"] = float(np.nansum(gm))
    if model_rain is not None:
        mr = model_rain.reindex(hours)[m].to_numpy()[ok]
        row["model_basin_mm"] = float(np.nansum(mr))
    return row


def rain_scale_reading(storm_df: pd.DataFrame, cv_max: float, spread_max: float, swing_min: float):
    """PROPOSED reading of the design's rule. Per storm, pool the gauges (sum of radar 3x3-mean totals over sum of
    gauge totals), then look at how much that ratio moves from storm to storm."""
    lines = []
    d = storm_df[storm_df["ratio_3x3_mean"].notna()]
    if d.empty:
        return ["Not enough rain at the gauges to judge (no storm reached the minimum gauge total).", ""], "undecided", None
    sums = d.groupby("storm")[["radar_3x3_mean_mm", "gauge_mm"]].sum()          # works on old and new pandas
    pooled = sums["radar_3x3_mean_mm"] / sums["gauge_mm"]
    pooled = pooled.reindex([s for s in storm_df["storm"].unique() if s in pooled.index])
    n_storms, n_gauges = len(pooled), d["gauge"].nunique()
    lines.append("Pooled radar / gauge ratio by storm (radar = mean of the 3x3 block; gauges added together):")
    for s, v in pooled.items():
        lines.append(f"    {s:<32} {v:6.2f}")
    if n_storms < 3 or n_gauges < 2:
        lines.append(f"Only {n_storms} storm(s) and {n_gauges} gauge(s) with usable totals: too few to call the ratio "
                     "consistent or swinging. Treat rain scale as a probe for now.")
        return lines + [""], "undecided", None
    mean, sd = float(pooled.mean()), float(pooled.std(ddof=1))
    cv, spread = sd / mean, float(pooled.max() / pooled.min())
    lines.append(f"Across {n_storms} storms: mean {mean:.2f}, spread (largest / smallest) {spread:.2f}, "
                 f"coefficient of variation {cv:.2f}.")
    if spread <= spread_max and cv <= cv_max:
        verdict = "consistent"
        lo, hi = 1 / pooled.max(), 1 / pooled.min()
        lines.append(f"PROPOSED READING: CONSISTENT (spread <= {spread_max:g} and CV <= {cv_max:g}). Rain scale can be a calibrated "
                     f"knob. The factor that would bring the radar to the gauges runs from {lo:.2f} to {hi:.2f} "
                     "(1 / ratio). That is the range the knob needs to cover, if the gauges are right and representative.")
    elif spread >= swing_min or cv > 1.5 * cv_max:
        verdict = "swings"
        lines.append(f"PROPOSED READING: SWINGS (spread >= {swing_min:g} or CV > {1.5 * cv_max:g}). Keep rain scale a probe: "
                     "one number cannot fix storms that disagree this much.")
        lines.append("A swinging ratio can come from the radar, from the gauges (point vs 1 km cell, catch error), or both. "
                     "Look at the hyetograph pictures before deciding which.")
    else:
        verdict = "borderline"
        lines.append("PROPOSED READING: BORDERLINE (between the consistent and swinging thresholds). Decide with the "
                     "pictures and the SRP team; storms with the lowest gauge totals carry the most noise.")
    lines.append(f"(Thresholds are proposals to be set before the first model run: consistent = spread <= {spread_max:g} and "
                 f"CV <= {cv_max:g}; swings = spread >= {swing_min:g} or CV > {1.5 * cv_max:g}.)")
    return lines + [""], verdict, pooled


# ================================================================================================================
# pictures
# ================================================================================================================
def style_axes(ax):
    ax.set_facecolor("white")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.grid(axis="y", color=GRID_COLOR, lw=0.7)
    ax.set_axisbelow(True)


def make_plots(out_dir: Path, gauges, storms, storm_df, lag_df, hours):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.dates as mdates
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib is not installed: skipping pictures)")
        return []
    made = []
    # ---- 1. one hyetograph picture per storm ---------------------------------------------------------------------
    for si, (sname, s0, s1) in enumerate(storms, start=1):
        m = (hours >= s0) & (hours < s1)
        if not m.any():
            continue
        fig, axes = plt.subplots(len(gauges), 1, figsize=(11, 2.5 * len(gauges) + 0.9), sharex=True, squeeze=False)
        wet = np.zeros(len(hours), dtype=bool)               # zoom on the wet part of the window (any gauge or radar cell)
        for G in gauges:
            wet |= np.nan_to_num(G.gauge.to_numpy(), nan=0.0) > 0.05
            wet |= np.nan_to_num(G.cell.to_numpy(), nan=0.0) > 0.05
        wet &= np.asarray(m)
        if wet.any():
            lo = max(s0, hours[wet][0] - 3 * HOUR)
            hi = min(s1, hours[wet][-1] + 4 * HOUR)
        else:
            lo, hi = s0, s1
        for ax, G in zip(axes[:, 0], gauges):
            style_axes(ax)
            t = hours[m]
            ax.bar(t, G.gauge[m].to_numpy(), width=1 / 24, align="edge", color=GAUGE_BAR, label="gauge, hourly")
            ax.fill_between(t, G.nb_min[m].to_numpy(), G.nb_max[m].to_numpy(), step="post", color=RADAR_BAND,
                            alpha=0.65, lw=0, label="radar, lowest to highest of the 3x3 cells")
            ax.step(t, G.cell[m].to_numpy(), where="post", color=SERIES_COLORS[0], lw=1.6, label="radar, the gauge's own cell")
            r = storm_df[(storm_df["gauge"] == G.name) & (storm_df["storm"] == sname)]
            if len(r) and r["hours_used"].iloc[0] > 0:
                r = r.iloc[0]
                txt = (f"{G.name}:  gauge {r.gauge_mm:.1f} mm   radar cell {r.radar_cell_mm:.1f} mm "
                       f"({'ratio %.2f' % r.ratio_cell if np.isfinite(r.ratio_cell) else 'no ratio'})   "
                       f"3x3 mean {r.radar_3x3_mean_mm:.1f} mm ({'%.2f' % r.ratio_3x3_mean if np.isfinite(r.ratio_3x3_mean) else 'no ratio'})")
            else:
                txt = f"{G.name}: no hours with both gauge and radar"
            ax.set_title(txt, loc="left", fontsize=10, color=INK)
            ax.set_ylabel("mm per hour", color=INK_2, fontsize=9)
        axes[0, 0].legend(loc="upper right", fontsize=8, frameon=False)
        axes[-1, 0].set_xlim(lo, hi)
        axes[-1, 0].xaxis.set_major_locator(mdates.AutoDateLocator(minticks=5, maxticks=9))
        axes[-1, 0].xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:%M"))
        fig.suptitle(f"{sname}: gauge vs radar. Each bar or step covers the hour that STARTS at its left edge.\n"
                     f"Shown: {lo:%b %d %H:%M} to {hi:%b %d %H:%M} (storm window {s0:%b %d %H:%M} to {s1:%b %d %H:%M}).",
                     x=0.01, ha="left", fontsize=10, color=INK)
        fig.tight_layout()
        p = out_dir / f"fig_storm{si}_hyetographs.png"
        fig.savefig(p, dpi=140, facecolor="white")
        plt.close(fig)
        made.append(p)
    # ---- 2. ratio by storm, and totals against totals ------------------------------------------------------------
    d = storm_df[storm_df["ratio_3x3_mean"].notna()]
    if not d.empty:
        fig, (axa, axb) = plt.subplots(1, 2, figsize=(11.5, 4.4), gridspec_kw={"width_ratios": [1.5, 1]})
        for ax in (axa, axb):
            style_axes(ax)
        names = [s[0] for s in storms]
        n_g = len(gauges)
        for gi, G in enumerate(gauges):
            dd = d[d["gauge"] == G.name]
            xs = np.array([names.index(s) for s in dd["storm"]], dtype=float) + (gi - (n_g - 1) / 2) * 0.12
            col, mk = SERIES_COLORS[gi % 8], SERIES_MARKERS[gi % 8]
            axa.vlines(xs, dd["ratio_3x3_min_cell"], dd["ratio_3x3_max_cell"], color=col, lw=1.2, alpha=0.6)
            axa.plot(xs, dd["ratio_3x3_mean"], mk, color=col, ms=7, mec="white", mew=1.0, label=G.name)
            axb.plot(dd["gauge_mm"], dd["radar_3x3_mean_mm"], mk, color=col, ms=7, mec="white", mew=1.0, label=G.name)
        axa.axhline(1.0, color=INK_2, lw=1, ls=(0, (4, 3)))
        axa.set_xticks(range(len(names)))
        axa.set_xticklabels([textwrap.fill(n, 14) for n in names], fontsize=8)
        axa.set_ylabel("radar / gauge storm total", fontsize=9, color=INK_2)
        axa.set_title("Ratio by storm: dot = 3x3 mean, bar = lowest to highest single cell", loc="left", fontsize=10, color=INK)
        handles, labels = axa.get_legend_handles_labels()
        hi = float(np.nanmax([d["gauge_mm"].max(), d["radar_3x3_mean_mm"].max()])) * 1.08
        axb.plot([0, hi], [0, hi], color=INK_2, lw=1, ls=(0, (4, 3)))
        axb.set_xlim(0, hi)
        axb.set_ylim(0, hi)
        axb.set_xlabel("gauge storm total (mm)", fontsize=9, color=INK_2)
        axb.set_ylabel("radar 3x3 mean storm total (mm)", fontsize=9, color=INK_2)
        axb.set_title("Totals (dashed = 1:1)", loc="left", fontsize=10, color=INK)
        axb.set_aspect("equal", adjustable="box")
        fig.tight_layout(rect=(0, 0.07, 1, 1))
        fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False, fontsize=9)
        p = out_dir / "fig_ratio_by_storm.png"
        fig.savefig(p, dpi=140, facecolor="white")
        plt.close(fig)
        made.append(p)
    # ---- 3. lag test ---------------------------------------------------------------------------------------------
    pooled = lag_df[lag_df["storm"] == "all storms pooled"]
    if not pooled.empty:
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
        for ax, series, title in zip(axes, ("cell", "3x3 mean"), ("Radar: the gauge's own cell", "Radar: mean of 3x3 cells")):
            style_axes(ax)
            for gi, G in enumerate(gauges):
                q = pooled[(pooled["gauge"] == G.name) & (pooled["radar_series"] == series)]
                ax.plot(q["lag_h"], q["r"], "-" + SERIES_MARKERS[gi % 8], color=SERIES_COLORS[gi % 8], lw=1.5, ms=6,
                        mec="white", mew=0.8, label=G.name)
            ax.axvline(0, color=INK_2, lw=1, ls=(0, (4, 3)))
            ax.set_xticks(LAGS)
            ax.set_xlabel("lag in hours (+ means the radar is LATE)", fontsize=9, color=INK_2)
            ax.set_title(title, loc="left", fontsize=10, color=INK)
        axes[0].set_ylabel("hourly correlation, gauge vs radar", fontsize=9, color=INK_2)
        axes[0].legend(frameon=False, fontsize=8)
        fig.tight_layout()
        p = out_dir / "fig_lag_test.png"
        fig.savefig(p, dpi=140, facecolor="white")
        plt.close(fig)
        made.append(p)
    return made


# ================================================================================================================
# the whole check
# ================================================================================================================
def read_mrf_rain(path: str, start: str):
    sim = pd.read_csv(path)
    sim = sim[sim["RainMin_mm_hr"] < 9000]
    idx = pd.Timestamp(start) + pd.to_timedelta(sim["Time_hr"], unit="h")        # .mrf labels each hour by its END
    return pd.Series(sim["MAP_mm_hr"].to_numpy(float), index=idx - HOUR)         # -> START labels


def run_check(args) -> dict:
    storms = read_storms(args.storms)
    first = min(s[1] for s in storms) - 3 * HOUR
    last = max(s[2] for s in storms) + 3 * HOUR
    if args.span_start:
        first = pd.Timestamp(args.span_start)
    if args.span_end:
        last = pd.Timestamp(args.span_end)
    hours = hour_starts(first.floor("60min"), last.ceil("60min"))
    gcfg = read_gauge_config(args.gauges, args.gauge_utc_offset)
    radar_dir = Path(args.radar_dir)
    if not radar_dir.is_dir():
        die(f"radar folder not found: {radar_dir}")

    print(f"radar_gauge_check_sycamore.py v{VERSION}")
    print(f"  radar: {radar_dir}   files {args.radar_prefix}MMDDYYYYHH.{args.radar_ext}   labels = {args.radar_label.upper()} of the hour")
    print(f"  hours loaded: {hours[0]:%Y-%m-%d %H:%M} to {hours[-1] + HOUR:%Y-%m-%d %H:%M} ({len(hours)} hours; each labeled by its START)")
    print(f"  gauges: {len(gcfg)}   storms: {len(storms)}   minimum gauge total for a ratio: {args.min_gauge_mm:g} mm")

    if args.dry_run:
        n_have = sum((radar_dir / f"{args.radar_prefix}{(t + HOUR if args.radar_label == 'end' else t):%m%d%Y%H}.{args.radar_ext}").exists() for t in hours)
        print(f"\nDRY RUN: radar files present {n_have} of {len(hours)} needed" + ("" if n_have == len(hours) else "  <-- some are missing"))
        if n_have:
            first_have = next(t for t in hours if (radar_dir / f"{args.radar_prefix}{(t + HOUR if args.radar_label == 'end' else t):%m%d%Y%H}.{args.radar_ext}").exists())
            t_label = first_have + HOUR if args.radar_label == "end" else first_have
            hdr, _ = read_ascii(radar_dir / f"{args.radar_prefix}{t_label:%m%d%Y%H}.{args.radar_ext}")
            print(f"  grid: {int(hdr['ncols'])} x {int(hdr['nrows'])} cells of {hdr['cellsize']:g} m, xll {hdr['xllcorner']:.0f}, yll {hdr['yllcorner']:.0f}")
            for g in gcfg:
                row, col, inside = locate_cell(g["x"], g["y"], hdr)
                print(f"  gauge {g['name']!r}: x {g['x']:.0f}, y {g['y']:.0f} (from {g['location_from']}) -> row {row}, col {col} "
                      f"{'(inside the grid)' if inside else '<-- OUTSIDE THE GRID'}")
        for g in gcfg:
            hr, info = load_gauge_hourly(g, hours, args.model_utc_offset)
            print(f"  gauge {g['name']!r}: {info['n_records']} records ({info['first_record']} to {info['last_record']}), "
                  f"{int(hr.notna().sum())} of {len(hours)} hours covered, hours with rain {int((hr > 0).sum())}, "
                  f"total {float(hr.sum()):.1f} mm")
        print("  storm windows:")
        for n, a, b in storms:
            print(f"    {n:<32} {a:%Y-%m-%d %H:%M} to {b:%Y-%m-%d %H:%M}  ({int((b - a) / HOUR)} h)")
        print("\nNothing was computed or written.")
        return {}

    stack, hdr, missing, bad = load_radar(radar_dir, args.radar_prefix, args.radar_ext, hours, args.radar_label)
    print(f"  grid: {int(hdr['ncols'])} x {int(hdr['nrows'])} cells of {hdr['cellsize']:g} m, xll {hdr['xllcorner']:.0f}, yll {hdr['yllcorner']:.0f}")
    if missing:
        print(f"  WARNING: {len(missing)} radar hours have no file (first {missing[0]:%b %d %H:%M}); they are left out of both gauge and radar.")
    if bad:
        print(f"  WARNING: {len(bad)} radar files could not be used, e.g. {bad[:2]}")
    grid_mean = pd.Series(nanagg(stack.reshape(len(hours), -1).astype("float64"), np.nanmean, 1), index=hours)

    model_rain = read_mrf_rain(args.mrf, args.mrf_start) if args.mrf else None
    gauges = []
    for g in gcfg:
        hr, info = load_gauge_hourly(g, hours, args.model_utc_offset)
        G = GaugeSeries(g, info, hours, hr, stack, hdr)
        if not G.inside:
            print(f"  WARNING: gauge {G.name!r} (x {g['x']:.0f}, y {g['y']:.0f}) is OUTSIDE the radar grid: skipped.")
            continue
        gauges.append(G)
        print(f"  gauge {G.name!r}: cell row {G.row}, col {G.col}; {info['n_records']} records"
              + (f", {info['n_resets']} counter resets" if info["n_resets"] else "")
              + f"; {int(hr.notna().sum())} of {len(hours)} hours covered; rain in the loaded span {float(hr.sum()):.1f} mm"
              + ("; EDGE CELL" if G.edge else ""))
        if info["n_negative"]:
            print(f"    note: {info['n_negative']} negative increments were set to 0")
    if not gauges:
        die("no gauge is inside the radar grid")

    # ---- storm table ----------------------------------------------------------------------------------------------------
    rows = []
    for G in gauges:
        for name, s0, s1 in storms:
            rows.append(one_storm_row(G, name, s0, s1, hours, grid_mean, model_rain, args.min_gauge_mm))
    storm_df = pd.DataFrame(rows)

    # ---- lag table ------------------------------------------------------------------------------------------------------
    lag_rows = []
    all_mask = np.zeros(len(hours), dtype=bool)
    for _, s0, s1 in storms:
        all_mask |= np.asarray((hours >= s0) & (hours < s1))
    for G in gauges:
        for series_name, ser in (("cell", G.cell), ("3x3 mean", G.nb_mean)):
            for label, mask in [("all storms pooled", all_mask)] + [(n, np.asarray((hours >= a) & (hours < b))) for n, a, b in storms]:
                for k, r, n in lag_table(G.gauge, ser, mask):
                    lag_rows.append({"gauge": G.name, "storm": label, "radar_series": series_name, "lag_h": k, "r": r, "n_hours": n})
    lag_df = pd.DataFrame(lag_rows)

    # ---- hourly pairs and daily table ----------------------------------------------------------------------------------
    pair_frames, daily_frames = [], []
    for G in gauges:
        df = pd.DataFrame({"hour_start": hours, "gauge": G.name, "gauge_mm": G.gauge.to_numpy(), "radar_cell_mm": G.cell.to_numpy(),
                           "radar_3x3_mean_mm": G.nb_mean.to_numpy(), "radar_3x3_min_mm": G.nb_min.to_numpy(),
                           "radar_3x3_max_mm": G.nb_max.to_numpy(), "in_storm_window": all_mask})
        df.insert(1, "hour_ending", df["hour_start"] + HOUR)
        pair_frames.append(df)
        both = df["gauge_mm"].notna() & df["radar_cell_mm"].notna()
        d = df[both].copy()
        d["day"] = d["hour_start"].dt.floor("D")
        dd = d.groupby("day")[["gauge_mm", "radar_cell_mm", "radar_3x3_mean_mm"]].sum()
        dd["hours_used"] = d.groupby("day").size()
        dd.insert(0, "gauge", G.name)
        daily_frames.append(dd.reset_index())
    pairs_df = pd.concat(pair_frames, ignore_index=True)
    daily_long = pd.concat(daily_frames, ignore_index=True)
    daily_wide = daily_long.pivot(index="day", columns="gauge", values=["gauge_mm", "radar_cell_mm", "radar_3x3_mean_mm"])
    daily_wide.columns = [f"{b} | {a}" for a, b in daily_wide.columns]
    gm_day = grid_mean.groupby(grid_mean.index.floor("D")).sum()
    daily_wide.insert(0, "grid_mean_radar_mm", gm_day.reindex(daily_wide.index))
    if model_rain is not None:
        mr = model_rain.reindex(hours)
        daily_wide.insert(1, "model_basin_mm", mr.groupby(mr.index.floor("D")).sum(min_count=1).reindex(daily_wide.index))

    # ---- save -----------------------------------------------------------------------------------------------------------
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    storm_df.to_csv(out_dir / "radar_gauge_storm_table.csv", index=False)
    lag_df.to_csv(out_dir / "radar_gauge_lag_test.csv", index=False)
    pairs_df.to_csv(out_dir / "radar_gauge_hourly_pairs.csv", index=False)
    daily_wide.to_csv(out_dir / "radar_gauge_daily.csv")

    # ---- report ---------------------------------------------------------------------------------------------------------
    L = []
    L.append(f"RADAR vs GAUGE CHECK (Stage 0)   radar_gauge_check_sycamore.py v{VERSION}")
    L.append(f"Radar: {radar_dir}   labels = {args.radar_label.upper()} of the hour; hours below are labeled by when they START.")
    L.append(f"Hours {hours[0]:%Y-%m-%d %H:%M} to {hours[-1] + HOUR:%Y-%m-%d %H:%M}; radar hours missing: {len(missing)}.")
    L.append("")
    L.append("1. STORM TOTALS  (hours where both gauge and radar have a value; radar / gauge ratio only if gauge total >= "
             f"{args.min_gauge_mm:g} mm)")
    show = ["gauge", "storm", "coverage_pct", "gauge_mm", "radar_cell_mm", "radar_3x3_mean_mm", "ratio_cell", "ratio_3x3_mean",
            "ratio_3x3_min_cell", "ratio_3x3_max_cell"]
    ren = {"coverage_pct": "hours used %", "gauge_mm": "gauge mm", "radar_cell_mm": "radar cell mm", "radar_3x3_mean_mm": "radar 3x3 mean mm",
           "ratio_cell": "ratio cell", "ratio_3x3_mean": "ratio 3x3 mean", "ratio_3x3_min_cell": "ratio lowest cell",
           "ratio_3x3_max_cell": "ratio highest cell"}
    L.append(storm_df[show].rename(columns=ren).to_string(index=False, float_format=lambda x: f"{x:.2f}", na_rep="-"))
    low_cov = storm_df[storm_df["coverage_pct"] < 90]
    if len(low_cov):
        L.append("  WARNING: these gauge/storm pairs use under 90% of the window's hours (gauge file or radar gap): " +
                 "; ".join(f"{r.gauge} / {r.storm} ({r.coverage_pct:.0f}%)" for r in low_cov.itertuples()))
    notes = storm_df[storm_df["note"].fillna("") != ""]
    for r in notes.itertuples():
        L.append(f"  note: {r.gauge} / {r.storm}: {r.note}")
    L.append("")
    L.append("2. TIMING  (radar minus gauge; + means the radar is later)")
    tshow = ["gauge", "storm", "gauge_peak_hour_mm", "radar_cell_peak_hour_mm", "peak_hour_diff_h", "centre_of_rain_diff_h", "best_lag_h",
             "r_at_lag0", "r_at_best_lag"]
    tren = {"gauge_peak_hour_mm": "gauge peak mm/h", "radar_cell_peak_hour_mm": "radar peak mm/h", "peak_hour_diff_h": "peak hour diff h",
            "centre_of_rain_diff_h": "centre of rain diff h", "best_lag_h": "best lag h", "r_at_lag0": "r at lag 0", "r_at_best_lag": "r at best lag"}
    L.append(storm_df[tshow].rename(columns=tren).to_string(index=False, float_format=lambda x: f"{x:.2f}", na_rep="-"))
    L.append("  Peak hour is the single wettest hour: for a storm with several similar hours it can jump around. The centre of rain")
    L.append("  (rain-weighted mean time) and the lag test below are steadier.")
    L.append("")
    L.append("3. LAG TEST, all storm hours pooled  (hourly correlation of gauge at t with radar at t+lag; best lag marked)")
    for G in gauges:
        for series_name in ("cell", "3x3 mean"):
            q = lag_df[(lag_df["gauge"] == G.name) & (lag_df["storm"] == "all storms pooled") & (lag_df["radar_series"] == series_name)]
            k, rb, r0 = best_lag([(a, b, c) for a, b, c in zip(q["lag_h"], q["r"], q["n_hours"])])
            cells = "  ".join(f"{int(a):+d}:{b:.2f}" if np.isfinite(b) else f"{int(a):+d}:n.a." for a, b in zip(q["lag_h"], q["r"]))
            flag = "" if k == 0 else ("   <-- not lag 0" if np.isfinite(k) else "")
            L.append(f"    {G.name:<28} radar {series_name:<9} {cells}   best lag {k:+.0f} h{flag}" if np.isfinite(k) else
                     f"    {G.name:<28} radar {series_name:<9} not enough hours for a correlation")
    lags_cell = []
    for G in gauges:
        q = lag_df[(lag_df["gauge"] == G.name) & (lag_df["storm"] == "all storms pooled") & (lag_df["radar_series"] == "3x3 mean")]
        k, _, _ = best_lag([(a, b, c) for a, b, c in zip(q["lag_h"], q["r"], q["n_hours"])])
        if np.isfinite(k):
            lags_cell.append(k)
    if lags_cell and all(k == 0 for k in lags_cell):
        L.append("  READING: every gauge lines up with the radar at lag 0 (3x3 mean), so the file labels and the gauge clock agree "
                 f"under the {args.radar_label.upper()}-of-hour convention.")
    elif lags_cell and len(set(lags_cell)) == 1:
        L.append(f"  READING: every gauge's best lag is {lags_cell[0]:+.0f} h. That points to a label or clock offset of about that size "
                 "(radar file labels, or the gauge file clock), not to a radar rain error. Check --radar-label and --gauge-utc-offset, then rerun.")
    elif lags_cell:
        L.append("  READING: the gauges do not agree on a best lag. Look at the hyetograph pictures: patchy storms moving across the "
                 "basin can make single gauges differ by an hour.")
    L.append("  This tests the radar file labels against the gauges only. It does not test how tRIBS reads those labels.")
    L.append("")
    L.append("4. RAIN-SCALE READING  (design rule: consistent ratio = calibrated knob; swinging ratio = probe)")
    reading, verdict, pooled = rain_scale_reading(storm_df, args.consistent_cv, args.consistent_spread, args.swing_spread)
    L.extend(reading)
    if args.mrf:
        L.append("5. MODEL RAIN vs RADAR GRID MEAN, per storm (context only: mm over the basin, from the .mrf, and mm averaged over every cell of the radar grid)")
        ms = storm_df.drop_duplicates("storm")[["storm"]].copy()
        for name, s0, s1 in storms:
            mm = (hours >= s0) & (hours < s1)
            ms.loc[ms["storm"] == name, "model basin mean mm"] = float(model_rain.reindex(hours)[mm].sum())
            ms.loc[ms["storm"] == name, "radar grid mean mm"] = float(grid_mean[mm].sum())
        L.append(ms.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
        L.append("  The radar grid also covers land outside the basin, so the two will not match exactly; the ratio between storms is what to look at.")
        L.append("")
    L.append("CAUTIONS")
    L.append("  * A gauge is a point and a radar cell is 1 km across. For patchy convective storms the two can differ even when both are right.")
    L.append("  * Tipping-bucket gauges under-count in very heavy rain and in wind; no catch correction is applied here.")
    L.append("  * MRMS gauge-corrected products use many gauges. If these FCDMC gauges feed the product, agreement is partly built in.")
    L.append("  * Three gauges and four storms is a small sample. Treat the proposed reading as a prompt for the SRP discussion, not a result.")
    L.append("")
    L.append("DAY BY DAY (days with 1 mm or more at any gauge or radar cell; use it to adjust the storm windows; full table in radar_gauge_daily.csv)")
    L.append("  Each gauge column reads: gauge mm / radar 3x3-mean mm. Days are labeled by when the hours START (midnight to midnight).")
    compact = pd.DataFrame(index=daily_wide.index)
    compact["radar grid mean mm"] = daily_wide["grid_mean_radar_mm"]
    if "model_basin_mm" in daily_wide.columns:
        compact["model basin mm"] = daily_wide["model_basin_mm"]
    keep = pd.Series(False, index=daily_wide.index)
    for G in gauges:
        gm, rm = daily_wide[f"{G.name} | gauge_mm"], daily_wide[f"{G.name} | radar_3x3_mean_mm"]
        keep |= (gm.fillna(0) >= 1.0) | (rm.fillna(0) >= 1.0)
        compact[G.name] = [("-" if pd.isna(a) and pd.isna(b) else f"{a:5.1f} / {b:5.1f}") for a, b in zip(gm, rm)]
    compact = compact[keep]
    compact.index = compact.index.strftime("%b %d")
    L.append(compact.to_string(float_format=lambda x: f"{x:.1f}", na_rep="-") if len(compact) else "  (no such days)")
    report = "\n".join(L)
    (out_dir / "radar_check_report.txt").write_text(report + "\n")
    print("\n" + report)
    made = [] if args.no_plots else make_plots(out_dir, gauges, storms, storm_df, lag_df, hours)
    print("\nFILES WRITTEN to " + str(out_dir.resolve()))
    for name in ("radar_check_report.txt", "radar_gauge_storm_table.csv", "radar_gauge_daily.csv", "radar_gauge_hourly_pairs.csv",
                 "radar_gauge_lag_test.csv"):
        print("  " + name)
    for p in made:
        print("  " + p.name)
    return {"storm_df": storm_df, "lag_df": lag_df, "daily": daily_wide, "pairs": pairs_df, "verdict": verdict, "pooled": pooled,
            "missing_hours": missing, "gauges": gauges}


# ================================================================================================================
# templates
# ================================================================================================================
GAUGE_TEMPLATE = """# gauges_template.csv  -- fill in the blanks, then save as gauges.csv.   Lines starting with # are ignored.
# name        what to call the gauge in tables and pictures
# file        the gauge record (CSV with a header). Relative paths are looked up next to this file, then where you run the script
# x,y         EPSG:26912 (UTM 12N) metres        OR        lon,lat   decimal degrees (negative longitude)  -- give one pair
# units       mm | in | cm    (depth units of the value column; default mm)
# value_type  increment (depth per record; default) | cumulative (a counter that grows; a drop is treated as a reset to zero)
# time_col, value_col   column names; default = first column and second column of the file
# scale       a multiplier applied first (for example the depth of one tip, if the file counts tips; default 1)
# no_record   zero | missing : what an hour with NO record means. Event-triggered ALERT-style logs only write when it rains,
#             so use zero (default). Use missing if a gap means the logger was down.
# valid_from, valid_to   the period the file really covers (local clock). Default: first to last record in the file.
# utc_offset  the clock of this gauge file: -7 = Arizona time (default -7), 0 = UTC
name,file,x,y,lon,lat,units,value_type,time_col,value_col,scale,no_record,valid_from,valid_to,utc_offset
Mt. Ord (2000),,,,,,in,increment,,,1,zero,,,-7
Sycamore Cr Upper (82200),,,,,,in,increment,,,1,zero,,,-7
Iron Dike (83000),,,,,,in,increment,,,1,zero,,,-7
"""

STORM_TEMPLATE = """# storms_template.csv  -- the built-in defaults, one row per storm. Edit and save as storms.csv, then use --storms storms.csv
# Local Arizona time. A window contains the hours that START inside [start, end). Use whole hours.
name,start,end
""" + "\n".join(f"{n},{a},{b}" for n, a, b in DEFAULT_STORMS) + "\n"


def do_template(out_dir: str) -> int:
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    (d / "gauges_template.csv").write_text(GAUGE_TEMPLATE)
    (d / "storms_template.csv").write_text(STORM_TEMPLATE)
    print(f"wrote {d / 'gauges_template.csv'} and {d / 'storms_template.csv'}")
    print("Fill in the gauge files and locations (the script does not have them), then run with --dry-run first.")
    return 0


# ================================================================================================================
# self-test: synthetic data with known answers
# ================================================================================================================
def _write_ascii(path: Path, arr: np.ndarray, xll: float, yll: float, cs: float):
    header = f"ncols {arr.shape[1]}\nnrows {arr.shape[0]}\nxllcorner {xll:g}\nyllcorner {yll:g}\ncellsize {cs:g}\nNODATA_value -9999\n"
    with open(path, "w") as fh:
        fh.write(header)
        np.savetxt(fh, arr, fmt="%.4f")


def _selftest_build(tmp: Path, radar_lag_h: int, drop_hours=()):
    """Synthetic world: 3 gauges, radar = 0.8 x gauge in each gauge's cell (so the true ratio is 0.8), optional lag."""
    rng = np.random.default_rng(11)
    storms = [(n, pd.Timestamp(a), pd.Timestamp(b)) for n, a, b in DEFAULT_STORMS]
    first, last = min(s[1] for s in storms) - 3 * HOUR, max(s[2] for s in storms) + 3 * HOUR
    hours = hour_starts(first, last)
    T = len(hours)
    nrows, ncols, xll, yll, cs = 32, 49, 440000.0, 3720000.0, 1000.0
    cells = [(10, 10), (20, 30), (8, 40)]
    true = np.zeros((3, T))
    for name, s0, s1 in storms:
        start = int((s0 - first) / HOUR) + 5
        dur = int(rng.integers(8, 13))
        shape = rng.gamma(2.0, 1.0, dur)
        shape /= shape.sum()
        for gi in range(3):
            total = 30.0 * (0.8 + 0.4 * rng.random())
            true[gi, start:start + dur] += shape * total
    radar_t = np.zeros_like(true)
    for gi in range(3):
        sh = np.roll(true[gi], radar_lag_h)
        if radar_lag_h > 0:
            sh[:radar_lag_h] = 0
        radar_t[gi] = 0.8 * sh
    rdir = tmp / "radar"
    rdir.mkdir()
    dropped = []
    for i, t0 in enumerate(hours):
        grid = np.full((nrows, ncols), 0.8 * true[:, i].mean() * 0.5)
        for gi, (r, c) in enumerate(cells):
            grid[r - 1:r + 2, c - 1:c + 2] = radar_t[gi, i] * (1 + 0.05 * rng.standard_normal((3, 3)))
            grid[r, c] = radar_t[gi, i]
        grid = np.clip(grid, 0, None)
        if t0 in drop_hours:
            dropped.append(t0)
            continue
        _write_ascii(rdir / f"radar_{t0 + HOUR:%m%d%Y%H}.asc", grid, xll, yll, cs)
    # gauge files: 4 records per wet hour at :15 :30 :45 :00(end); dry hours write nothing (event-triggered log)
    gdir = tmp / "gauges"
    gdir.mkdir()
    cfg = ["name,file,x,y,lon,lat,units,value_type,time_col,value_col,scale,no_record,valid_from,valid_to,utc_offset"]
    xs = [xll + (c + 0.5) * cs for r, c in cells]
    ys = [yll + (nrows - 1 - r + 0.5) * cs for r, c in cells]
    for gi in range(3):
        recs = []
        for i in np.nonzero(true[gi] > 0)[0]:
            for q in (15, 30, 45, 60):
                recs.append((hours[i] + pd.Timedelta(minutes=q), true[gi, i] / 4.0))
        times = pd.DatetimeIndex([a for a, _ in recs])
        vals = np.array([b for _, b in recs])
        if gi == 0:                                    # inches, increments, Arizona clock
            df = pd.DataFrame({"time": times.strftime("%Y-%m-%d %H:%M:%S"), "inches": vals / 25.4})
            cfg.append(f"Gauge A,gauges/a.csv,{xs[0]},{ys[0]},,,in,increment,time,inches,1,zero,{first:%Y-%m-%d %H:%M},{last:%Y-%m-%d %H:%M},-7")
            df.to_csv(gdir / "a.csv", index=False)
        elif gi == 1:                                  # mm, cumulative with a counter reset after the second storm, location given as lon/lat
            cum = np.cumsum(vals)
            reset_at = int(np.searchsorted(times, pd.Timestamp("2017-02-05 00:00")))
            cum2 = cum.copy()
            cum2[reset_at:] = cum[reset_at:] - cum[reset_at - 1]
            # a real counter has a reading before the first tip: the first record of a cumulative file is only a baseline
            times_b = times.insert(0, first)
            cum_b = np.insert(cum2, 0, 0.0)
            df = pd.DataFrame({"t": times_b.strftime("%Y-%m-%d %H:%M:%S"), "total_mm": cum_b})
            # give the location as lon/lat computed from the cell centre by inverting with a small search
            lon, lat = _utm_to_lonlat_search(xs[1], ys[1])
            cfg.append(f"Gauge B,gauges/b.csv,,,{lon:.7f},{lat:.7f},mm,cumulative,t,total_mm,1,zero,{first:%Y-%m-%d %H:%M},{last:%Y-%m-%d %H:%M},-7")
            df.to_csv(gdir / "b.csv", index=False)
        else:                                          # tips of 0.01 in, file clock UTC
            df = pd.DataFrame({"when_utc": (times + pd.Timedelta(hours=7)).strftime("%Y-%m-%dT%H:%M:%S"),
                               "tips": vals / (0.01 * 25.4)})
            cfg.append(f"Gauge C,gauges/c.csv,{xs[2]},{ys[2]},,,in,increment,when_utc,tips,0.01,zero,{first:%Y-%m-%d %H:%M},{last:%Y-%m-%d %H:%M},0")
            df.to_csv(gdir / "c.csv", index=False)
    (tmp / "gauges.csv").write_text("\n".join(cfg) + "\n")
    return {"hours": hours, "true": true, "storms": storms, "dropped": dropped}


def _utm_to_lonlat_search(x, y, lon_guess=-111.7, lat_guess=33.9):
    """Invert lonlat_to_utm by Newton iteration (used only by the self-test to exercise the lon/lat path)."""
    lon, lat = lon_guess, lat_guess
    for _ in range(30):
        x0, y0 = lonlat_to_utm(lon, lat)
        dx, dy = x - x0, y - y0
        if abs(dx) < 1e-4 and abs(dy) < 1e-4:
            break
        x1, y1 = lonlat_to_utm(lon + 1e-5, lat)
        x2, y2 = lonlat_to_utm(lon, lat + 1e-5)
        j = np.array([[(x1 - x0) / 1e-5, (x2 - x0) / 1e-5], [(y1 - y0) / 1e-5, (y2 - y0) / 1e-5]])
        step = np.linalg.solve(j, [dx, dy])
        lon, lat = lon + step[0], lat + step[1]
    return lon, lat


def _ns(**kw):
    base = dict(radar_dir=None, gauges=None, storms=None, out_dir=None, radar_prefix="radar_", radar_ext="asc", radar_label="end",
                gauge_utc_offset=-7.0, model_utc_offset=-7.0, min_gauge_mm=5.0, consistent_cv=0.20, consistent_spread=1.30,
                swing_spread=1.60, mrf=None, mrf_start=None, span_start=None, span_end=None, dry_run=False, no_plots=False)
    base.update(kw)
    return argparse.Namespace(**base)


def do_self_test() -> int:
    ok = True

    def check(cond, text):
        nonlocal ok
        print(("  PASS  " if cond else "  FAIL  ") + text)
        ok = ok and bool(cond)

    print("SELF-TEST: synthetic radar and gauges with known answers (nothing here is real data)")
    # (0) the built-in UTM conversion against values from pyproj (EPSG:4326 -> EPSG:26912) computed when this script was written
    refs = [((-111.0, 34.0), (500000.0, 3762155.977)), ((-111.7, 33.9), (435280.253, 3751288.777)),
            ((-112.1, 34.4), (398893.738, 3807057.001)), ((-110.5, 33.2), (546602.490, 3673570.663))]
    err = max(max(abs(lonlat_to_utm(*a)[0] - b[0]), abs(lonlat_to_utm(*a)[1] - b[1])) for a, b in refs)
    check(err < 0.05, f"lon/lat -> UTM 12N matches pyproj reference points to within {err:.4f} m")
    tmp = Path(tempfile.mkdtemp(prefix="radar_check_selftest_"))
    try:
        # (1) radar exactly 0.8 x gauge, no lag, three radar hours deleted
        drop = [pd.Timestamp("2017-01-20 14:00"), pd.Timestamp("2017-01-20 15:00"), pd.Timestamp("2017-02-19 03:00")]
        w1 = tmp / "world1"
        w1.mkdir()
        _selftest_build(w1, radar_lag_h=0, drop_hours=drop)
        res = run_check(_ns(radar_dir=str(w1 / "radar"), gauges=str(w1 / "gauges.csv"), out_dir=str(w1 / "out"), no_plots=False))
        sd = res["storm_df"]
        used = sd[sd["gauge_mm"] >= 5]
        check(len(used) == 12, f"12 gauge/storm pairs have a ratio (got {len(used)})")
        check(np.allclose(used["ratio_cell"], 0.8, atol=0.02), f"ratio in the gauge's own cell = 0.80 for every pair (range {used['ratio_cell'].min():.3f}-{used['ratio_cell'].max():.3f}), "
              "which also checks inch conversion, the cumulative counter reset, tip scaling, the UTC shift and the lon/lat path")
        check(len(res["missing_hours"]) == 3, "the 3 deleted radar hours were reported missing")
        gA = sd[(sd["gauge"] == "Gauge A") & (sd["storm"] == "Jan 18-23 flood storm")].iloc[0]
        check(gA["hours_used"] == gA["hours_in_window"] - 2, "the two deleted hours on Jan 20 were left out of both gauge and radar")
        check((sd["best_lag_h"].dropna() == 0).all(), "best lag is 0 for every storm and gauge")
        lg = res["lag_df"]
        pooled = lg[(lg["storm"] == "all storms pooled") & (lg["radar_series"] == "cell")]
        bl = [best_lag([(a, b, c) for a, b, c in zip(q["lag_h"], q["r"], q["n_hours"])])[0] for _, q in pooled.groupby("gauge")]
        check(all(k == 0 for k in bl), f"pooled best lag is 0 for all gauges (got {bl})")
        check(abs(sd["centre_of_rain_diff_h"].dropna()).max() < 0.1, "centre-of-rain difference is about 0 h")
        check(res["verdict"] == "consistent", f"proposed reading is 'consistent' (got {res['verdict']!r})")
        pf = [p.name for p in (w1 / "out").glob("*.png")]
        check(len(pf) >= 6, f"pictures were written ({len(pf)} png files)")
        # (2) the radar is one hour LATE
        w2 = tmp / "world2"
        w2.mkdir()
        _selftest_build(w2, radar_lag_h=1)
        res2 = run_check(_ns(radar_dir=str(w2 / "radar"), gauges=str(w2 / "gauges.csv"), out_dir=str(w2 / "out"), no_plots=True))
        lg2 = res2["lag_df"]
        pooled2 = lg2[(lg2["storm"] == "all storms pooled") & (lg2["radar_series"] == "cell")]
        bl2 = [best_lag([(a, b, c) for a, b, c in zip(q["lag_h"], q["r"], q["n_hours"])])[0] for _, q in pooled2.groupby("gauge")]
        check(all(k == 1 for k in bl2), f"a radar made one hour late gives best lag +1 for all gauges (got {bl2})")
        cd = res2["storm_df"]["centre_of_rain_diff_h"].dropna()
        check(abs(cd.mean() - 1.0) < 0.15, f"centre-of-rain difference is about +1 h for the late radar (mean {cd.mean():.2f})")
        # (3) --radar-label start is the same files read one hour earlier: the totals must change, i.e. the option is not ignored
        res3 = run_check(_ns(radar_dir=str(w1 / "radar"), gauges=str(w1 / "gauges.csv"), out_dir=str(w1 / "out3"), no_plots=True, radar_label="start"))
        lag3 = res3["lag_df"]
        p3 = lag3[(lag3["storm"] == "all storms pooled") & (lag3["radar_series"] == "cell")]
        bl3 = [best_lag([(a, b, c) for a, b, c in zip(q["lag_h"], q["r"], q["n_hours"])])[0] for _, q in p3.groupby("gauge")]
        check(all(k == 1 for k in bl3), f"reading end-labeled files as start-labeled shows up as best lag +1 (got {bl3})")
        # (4) dry run writes nothing and does not fail
        before = set(p.name for p in (w1 / "out").iterdir())
        run_check(_ns(radar_dir=str(w1 / "radar"), gauges=str(w1 / "gauges.csv"), out_dir=str(w1 / "out"), dry_run=True))
        after = set(p.name for p in (w1 / "out").iterdir())
        check(before == after, "--dry-run wrote nothing")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\nSELF-TEST " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


# ================================================================================================================
# command line
# ================================================================================================================
def build_parser():
    p = argparse.ArgumentParser(description="Stage 0: compare hourly MRMS radar rain with rain gauges inside the Lower Sycamore watershed "
                                            "(see the text at the top of this file).",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("template", help="write gauges_template.csv and storms_template.csv")
    t.add_argument("--out-dir", default=".")
    sub.add_parser("self-test", help="run the logic checks on synthetic data")
    r = sub.add_parser("run", help="run the check")
    r.add_argument("--radar-dir", required=True, help="folder holding radar_MMDDYYYYHH.asc")
    r.add_argument("--gauges", required=True, help="gauges CSV (see `template`)")
    r.add_argument("--storms", help="storms CSV (default: built-in windows)")
    r.add_argument("--out-dir", default="radar_check_out")
    r.add_argument("--radar-prefix", default="radar_")
    r.add_argument("--radar-ext", default="asc")
    r.add_argument("--radar-label", choices=["end", "start"], default="end", help="does the hour in the file name mark the END (default) or START of the rain hour?")
    r.add_argument("--gauge-utc-offset", type=float, default=-7.0, help="clock of the gauge files (-7 = Arizona time, 0 = UTC); a utc_offset column overrides")
    r.add_argument("--model-utc-offset", type=float, default=-7.0, help="clock of the model (Arizona = -7)")
    r.add_argument("--min-gauge-mm", type=float, default=5.0, help="smallest gauge storm total that gets a ratio (default 5 mm)")
    r.add_argument("--consistent-cv", type=float, default=0.20, help="PROPOSED: coefficient of variation of storm ratios for 'consistent'")
    r.add_argument("--consistent-spread", type=float, default=1.30, help="PROPOSED: largest / smallest storm ratio for 'consistent'")
    r.add_argument("--swing-spread", type=float, default=1.60, help="PROPOSED: largest / smallest storm ratio for 'swings'")
    r.add_argument("--mrf", help="optional: a model .mrf to compare the model's basin rain with the radar grid")
    r.add_argument("--mrf-start", help="start of that model run, e.g. '2016-10-01 00:00' (needed with --mrf)")
    r.add_argument("--span-start", help="override the first hour loaded (default: 3 h before the first storm)")
    r.add_argument("--span-end", help="override the end of the hours loaded (default: 3 h after the last storm)")
    r.add_argument("--no-plots", action="store_true")
    r.add_argument("--dry-run", action="store_true", help="check the inputs and print the plan; compute and write nothing")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd == "template":
        return do_template(args.out_dir)
    if args.cmd == "self-test":
        return do_self_test()
    if args.mrf and not args.mrf_start:
        die("--mrf needs --mrf-start (the model run's start, e.g. '2016-10-01 00:00')")
    run_check(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
