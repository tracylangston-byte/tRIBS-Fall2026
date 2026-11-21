#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
score_windows_sycamore.py   (v1.0, 2026-11-21)
==============================================

Scores the Lower Sycamore runs against the USGS 09510200 gauge by WINDOW, with storm windows found from the RAIN (not from the flow),
and writes tables of PBIAS, KGE', peak flow, timing and the runoff mechanisms. It also compares the new Sep-1-start WY2017 runs with the
old Oct-1-start WY2017 runs (the spin-up check). READ-ONLY: it opens result files and writes CSV and text files into an output folder.
It runs no model and never overwrites a file it did not write itself.

THREE COMMANDS
  windows     finds the 12 storm windows from the basin-mean rain (.mrf file), tests the rain threshold, the dry gap and R (days after the
              last rain) using RAIN AND GAUGE ONLY, checks the gauge files against the numbers on record, and writes storm_windows.csv
  score       scores every run in every window and writes the tables (needs storm_windows.csv from the windows command)
  self-test   builds a small fake world with known answers and checks every number this script produces

  python score_windows_sycamore.py windows --data bundle_folder --out-dir out
  python score_windows_sycamore.py score   --data bundle_folder --windows out/storm_windows.csv --out-dir out
  python score_windows_sycamore.py self-test

THE WINDOWS (all times are Arizona local time; a window is [start, end): the end hour is not included)
  OctMay2017, OctMay2020   Oct 1 00:00 to May 31 23:00 of the water year (September is spin-up and is never scored; the runs stop at hour
                           6551 / 6575, so the last hour of May 31 is never simulated and is left out for the gauge too)
  Both                     both of those together (the two runs of one setting, pooled)
  12 storms                start = the first hour whose basin-mean rain is at least THR mm/h after at least GAP dry hours;
                           end = the end of the LAST rain hour of the storm + R days, or the start of the next storm, whichever is first.
                           "Storm" for the cut = any rain cluster with at least MIN_NEXT mm of basin-mean rain. The rain is the same in every
                           run, so every run shares the same windows. THR, GAP, MIN_NEXT and R are chosen with the windows command and are
                           written in the header of storm_windows.csv.
  AllStorms                the 12 storm windows together (AllStorms_2017 and AllStorms_2020 are the same split by year)
  OFF_2017, OFF_2020, OFF_both   everything in Oct-May that is outside the 12 storm windows (this includes the smaller storms)
  Oct2016, OctDec2016      two extra windows for the spin-up check only (Oct 1 to Nov 1 and Oct 1 to Jan 1, 2016)

THE METRICS (model = the run, gauge = USGS; both turned into HOURLY MEANS first)
  Hourly mean         the average of the readings in a clock hour [hh:00, hh+1:00). The model writes every 3.75 minutes, the gauge every 15.
                      An hour with no reading is left out, never filled and never set to zero. Hours are named by their START.
  Matched hours       hours where BOTH the gauge and the model have a value. Every volume, PBIAS and KGE' uses matched hours, so a gauge
                      outage (for example Jan 1 2017) cannot make the model look too high. vol_model_allhours_Mm3 counts every model hour.
  Volume              sum of hourly means x 3,600 s, in M m3 (million cubic meters). 1 M m3 is about 811 acre-feet.
  PBIAS (%)           100 x (model volume - gauge volume) / gauge volume. MODEL MINUS GAUGE: POSITIVE MEANS THE MODEL IS TOO HIGH. This is
                      the opposite sign of the textbook definition (Gupta et al. 1999), which is gauge minus model.
  KGE'                Kling-Gupta efficiency, 2012 version: 1 - sqrt((r-1)^2 + (beta-1)^2 + (gamma-1)^2), where r = Pearson correlation,
                      beta = model mean / gauge mean (so PBIAS = 100 x (beta - 1)), gamma = (model sd / model mean) / (gauge sd / gauge mean).
                      1 is perfect. It is NOT gated: read it with beta next to it. The default basis is hourly means (--kge-basis 15min
                      interpolates the model onto the gauge's 15-minute times instead, as the old per-event scorer did).
  Gauge class         zero: gauge volume under 0.005 M m3 (no flow recorded); near-zero: under 0.05 M m3; ok: 0.05 M m3 or more.
                      For zero and near-zero windows the plain differences (volume in M m3, peak in cfs) are shown and every ratio-type
                      value (PBIAS, KGE', peak ratio, timing) is left blank with the reason in blank_reasons.
  Peak flow (storms)  the highest single reading in the window (gauge 15-minute, model 3.75-minute), in cfs (m3/s / 0.0283168).
                      peak_diff_cfs = model - gauge; peak_ratio = model / gauge; dt_peak_h = model peak time - gauge peak time (hours;
                      positive = the model peaks LATER).
  Timing (storms)     all in hours, MODEL MINUS GAUGE (positive = model later):
                      dt10_h, dt50_h, dt90_h   the time by which 10%, 50% and 90% of the window volume has passed (matched hours)
                      arr20_diff_h, arr50_diff_h   the time the flow first RISES THROUGH 20% (50%) of the GAUGE peak; blank if the model never
                                                  reaches that level (that blank is a result, not missing data)
                      lag_h   the whole-hour shift within +/-12 h at which the model hourly series correlates best with the gauge's
                              (positive = model later)
                      tail_vol_ratio   model volume / gauge volume from the gauge's t90 to the window end
                      Noise floor: the model (3.75-min) and gauge (15-min) sampling differ, so identical flow still shows dt10/50/90 of about 0.1 h.
  AllStorms row       PBIAS and KGE' are pooled over all matched hours of the 12 windows. Peak and timing values are the MEDIAN across the
                      storms whose gauge class is ok (n_ok_storms of them); peak_diff_cfs is the median over all storms that have both peaks.
  OFF rows            volume, PBIAS (only if the gauge class is ok) and the number of days whose highest flow is above 1 cfs and above
                      10 cfs, for the model and for the gauge (days with data only). No KGE', peak or timing.
  Mechanisms          from the .rft file: the runoff the model MADE on the land (all cells added, BEFORE routing, so not outlet flow):
                      Hsrf infiltration excess, Sbsrf subsurface stormflow, Psrf perched-water runoff, Satsrf saturation excess.
                      Volume = hourly value (m3/s) x 3,600 s, in M m3, over the window. Shares are percent of the four added. They are
                      described, not scored. rain_mm is the basin-mean rain from the .mrf file over the window. The .mrf and .rft files
                      label each hour by its END; this script moves the label back one hour, so an hour belongs to the window it falls in.

THE SPIN-UP CHECK
  The new runs (Sycamore_WY2017sep_*, Time_hr 0 = Sep 1 2016) are compared with the old Oct-1-start WY2017 runs for Re-tuned and Original.
  Both are put on one clock using each run's own start date (from its .in file, else from its name), then scored in exactly the same
  2017 windows. The check is the same table for both, plus the difference (new minus old).

WHAT IT READS (from --data, searched below it): *_Outlet.qout, the .mrf and .rft files, the .in files, and the gauge CSV files.
WHAT IT WRITES (to --out-dir): storm_windows.csv, windows_grid_test.csv, windows_R_test.csv (windows command);
  window_scores.csv, window_mechanisms.csv, spinup_check.csv, score_report.txt (score command). Every CSV starts with a line that marks
  it as written by this script and with comment lines (#) that define the columns; read it with comment='#'.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = "1.0"
SELF_MARK = "# written by score_windows_sycamore.py"
HOUR = pd.Timedelta(hours=1)
DAY = pd.Timedelta(days=1)
CFS = 0.0283168                 # m3/s in one cfs
BASIN_KM2 = 424.6
ZERO_GAUGE_VOL = 0.005          # M m3: below this the gauge recorded no flow
MIN_GAUGE_VOL = 0.05            # M m3: below this a ratio to the gauge is not shown
LOW_MODEL_RATIO = 0.25          # note only: model volume under this fraction of the gauge's
MAX_LAG_H, MIN_COMMON_H = 12, 12
MAX_STEP_H = 3.0                # a gap longer than this between two readings is not interpolated across
DAYS_ABOVE_CFS = (1.0, 10.0)
MECH_TYPES = ["Hsrf", "Sbsrf", "Psrf", "Satsrf"]
OLD_STEMS = {"Retuned": "W5_d012_k10_f0p15_as300_ksx1p03_nocl", "Original": "W5_d012_k10_base_nocl"}
DEFAULT_START = {"WY2017sep": "2016-09-01 00:00", "WY2020": "2019-09-01 00:00", "WY2017old": "2016-10-01 00:00"}
# Oct 1 00:00 to May 31 23:00 (end not included): the runs stop at hour 6551 (WY2017) and 6575 (WY2020), so the last hour of May 31 is never simulated
PERIODS = {2017: (pd.Timestamp("2016-10-01"), pd.Timestamp("2017-05-31 23:00")), 2020: (pd.Timestamp("2019-10-01"), pd.Timestamp("2020-05-31 23:00"))}
# thresholds tried by the windows command
GRID_THR = (0.02, 0.05, 0.1, 0.25, 0.5)
GRID_GAP = (6, 12, 24, 36)
R_TEST = (3, 4, 5, 7)

# every storm in the storm table (year, name, first rain day, last rain day, basin mm, one of the 12 scored storms)
STORM_TABLE = [
    (2017, "Oct 2", "2016-10-02", "2016-10-02", 10.537, False), (2017, "Nov 3", "2016-11-03", "2016-11-03", 15.910, False),
    (2017, "Nov 20-21", "2016-11-20", "2016-11-21", 16.349, False), (2017, "Nov 27", "2016-11-27", "2016-11-27", 8.710, False),
    (2017, "Dec 15-17", "2016-12-15", "2016-12-17", 17.694, False), (2017, "Dec 21-25", "2016-12-21", "2016-12-25", 47.803, True),
    (2017, "Dec 30-Jan 2", "2016-12-30", "2017-01-02", 48.178, True), (2017, "Jan 13-16", "2017-01-13", "2017-01-16", 16.854, False),
    (2017, "Jan 19-21", "2017-01-19", "2017-01-21", 63.400, True), (2017, "Jan 23-24", "2017-01-23", "2017-01-24", 18.058, False),
    (2017, "Feb 18-19", "2017-02-18", "2017-02-19", 46.410, True), (2017, "Feb 27-28", "2017-02-27", "2017-02-28", 41.500, True),
    (2017, "Mar 23", "2017-03-23", "2017-03-23", 9.806, False), (2017, "May 8-9", "2017-05-08", "2017-05-09", 6.590, False),
    (2020, "Nov 19-22", "2019-11-19", "2019-11-22", 71.545, True), (2020, "Nov 28-30", "2019-11-28", "2019-11-30", 44.039, True),
    (2020, "Dec 3-4", "2019-12-03", "2019-12-04", 5.067, False), (2020, "Dec 7-9", "2019-12-07", "2019-12-09", 20.898, False),
    (2020, "Dec 23-28", "2019-12-23", "2019-12-28", 75.196, True), (2020, "Jan 21", "2020-01-21", "2020-01-21", 8.763, False),
    (2020, "Feb 10-11", "2020-02-10", "2020-02-11", 30.319, True), (2020, "Feb 22-23", "2020-02-22", "2020-02-23", 52.637, True),
    (2020, "Mar 10-13", "2020-03-10", "2020-03-13", 64.518, True), (2020, "Mar 18-19", "2020-03-18", "2020-03-19", 24.610, True),
]
# gauge volumes (M m3) on record for the OLD window rule: 00:00 on the first rain day to 00:00 four days after the last rain day
# (WY2017: cut at the next storm in the storm table; WY2020: not cut). Used only to check that the gauge files are the right ones.
GAUGE_CHECK = {
    (2017, "Dec 21-25"): 0.0000, (2017, "Dec 30-Jan 2"): 0.8158, (2017, "Jan 19-21"): 5.5667, (2017, "Feb 18-19"): 1.2193,
    (2017, "Feb 27-28"): 3.8067, (2020, "Nov 19-22"): 0.015, (2020, "Nov 28-30"): 1.473, (2020, "Dec 23-28"): 4.069,
    (2020, "Feb 10-11"): 0.386, (2020, "Feb 22-23"): 4.131, (2020, "Mar 10-13"): 5.626, (2020, "Mar 18-19"): 3.317,
}
GAUGE_CHECK_TOL = 0.005


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    sys.exit("\nSTOPPED: " + msg + "\n")


def isnum(x) -> bool:
    return x is not None and isinstance(x, (int, float, np.integer, np.floating)) and bool(np.isfinite(x))


def fnum(x, nd=2, signed=False) -> str:
    if not isnum(x):
        return "-"
    return (("{:+." + str(nd) + "f}") if signed else ("{:." + str(nd) + "f}")).format(x)


def write_csv(path, header_lines, df) -> None:
    path = Path(path)
    if path.exists():
        with open(path, encoding="utf-8", errors="replace") as fh:
            first = fh.readline()
        if not first.startswith(SELF_MARK):
            die("{} already exists and was not written by this script. Move it, or choose another --out-dir.".format(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(SELF_MARK + " v{}\n".format(VERSION))
        for line in header_lines:
            fh.write("# " + line + "\n")
        df.to_csv(fh, index=False, float_format="%.6g")


def write_text(path, text) -> None:
    path = Path(path)
    if path.exists():
        with open(path, encoding="utf-8", errors="replace") as fh:
            first = fh.readline()
        if not first.startswith(SELF_MARK):
            die("{} already exists and was not written by this script. Move it, or choose another --out-dir.".format(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(SELF_MARK + " v{}\n".format(VERSION) + text, encoding="utf-8")


# ================================================================================================ reading the files
def read_qout(path: Path, start: pd.Timestamp) -> pd.Series:
    """Outlet flow (m3/s). tRIBS 6.0.0: comma-delimited with a header (Time_hr, Q); an older whitespace layout is also accepted."""
    with open(path) as fh:
        fh.readline()
        second = fh.readline()
    if "," in second:
        raw = pd.read_csv(path, sep=",", header=0)
    else:
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None)
    raw = raw.iloc[:, :2].copy()
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    secs = np.round(raw["Time_hr"].to_numpy(float) * 3600.0).astype("int64")
    s = pd.Series(raw["Q"].to_numpy(float), index=start + pd.to_timedelta(secs, unit="s"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_gauge(path: Path) -> pd.Series:
    o = pd.read_csv(path)
    for c in ("datetime", "q_cms"):
        if c not in o.columns:
            die("{} has no column '{}'. Columns found: {}".format(path, c, list(o.columns)[:10]))
    s = pd.Series(pd.to_numeric(o["q_cms"], errors="coerce").to_numpy(float), index=pd.to_datetime(o["datetime"]))
    return s[~s.index.duplicated(keep="last")].sort_index()


def _hour_start_index(time_hr, start) -> pd.DatetimeIndex:
    """The .mrf and .rft files label each hour by its END. Move the label back one hour so it names the hour's start."""
    secs = np.round(np.asarray(time_hr, float) * 3600.0).astype("int64")
    return start + pd.to_timedelta(secs, unit="s") - HOUR


def read_mrf(path: Path, start) -> pd.Series:
    d = pd.read_csv(path)
    if "MAP_mm_hr" not in d.columns or "Time_hr" not in d.columns:
        die("{} has no Time_hr / MAP_mm_hr columns. Columns found: {}".format(path, list(d.columns)))
    d = d.apply(pd.to_numeric, errors="coerce").dropna(subset=["Time_hr", "MAP_mm_hr"])
    d = d[(d.drop(columns=["Time_hr"]) < 9000).all(axis=1)]            # an unfilled last row is flagged 9999
    s = pd.Series(d["MAP_mm_hr"].to_numpy(float), index=_hour_start_index(d["Time_hr"], start))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_rft(path: Path, start) -> pd.DataFrame:
    d = pd.read_csv(path)
    want = ["Time_hr"] + ["{}_m3_s".format(t) for t in MECH_TYPES]
    miss = [c for c in want if c not in d.columns]
    if miss:
        die("{} has no column(s) {}. Columns found: {}".format(path, miss, list(d.columns)))
    d = d[want].apply(pd.to_numeric, errors="coerce").dropna()
    d = d[(d[want[1:]] < 9000).all(axis=1)]
    out = d[want[1:]].copy()
    out.index = _hour_start_index(d["Time_hr"], start)
    return out[~out.index.duplicated(keep="last")].sort_index()


def read_in(path: Path) -> dict:
    """KEY: line followed by the value line (the tRIBS input file layout)."""
    kv = {}
    lines = Path(path).read_text(errors="replace").splitlines()
    for i, ln in enumerate(lines[:-1]):
        m = re.match(r"^([A-Z0-9_]+):", ln.strip())
        if m:
            kv[m.group(1)] = lines[i + 1].strip()
    return kv


def parse_startdate(txt: str):
    m = re.match(r"^(\d\d)/(\d\d)/(\d{4})/(\d\d)/(\d\d)$", txt.strip())
    if not m:
        return None
    mo, d, y, hh, mi = (int(x) for x in m.groups())
    return pd.Timestamp(year=y, month=mo, day=d, hour=hh, minute=mi)


# ================================================================================================ series helpers
def hourly(s: pd.Series) -> pd.Series:
    s = s.dropna()
    return s.resample("60min").mean().dropna() if len(s) else s


def sel(x, lo, hi):
    return x[(x.index >= lo) & (x.index < hi)]


def vol_mm3(h) -> float:
    return float(np.sum(h) * 3600.0 / 1e6)


def gauge_class(v: float, n_matched: int) -> str:
    if n_matched == 0 or not isnum(v):
        return "no data"
    if v < ZERO_GAUGE_VOL:
        return "zero"
    if v < MIN_GAUGE_VOL:
        return "near-zero"
    return "ok"


# ================================================================================================ runs and gauges
def classify(stem: str, old_map: dict):
    """Which set, year and setting a run name belongs to (None for any other run)."""
    for setting, st in old_map.items():
        if stem == st:
            return {"set": "WY2017old", "year": 2017, "setting": setting + "_old", "base": setting}
    for prefix, year, sname in (("Sycamore_WY2017sep_", 2017, "WY2017sep"), ("Sycamore_WY2020_", 2020, "WY2020")):
        if stem.startswith(prefix):
            tail = stem[len(prefix):]
            if tail in ("Retuned", "Original"):
                return {"set": sname, "year": year, "setting": tail, "base": tail}
            m = re.match(r"^(S\d+)_", tail)
            if m:
                return {"set": sname, "year": year, "setting": m.group(1), "base": m.group(1)}
    return None


def setting_sort_key(s: str):
    base = s.replace("_old", "")
    order = {"Retuned": 0, "Original": 1}
    n = int(base[1:]) if re.match(r"^S\d+$", base) else 99
    return (1 if s.endswith("_old") else 0, order.get(base, 2), n)


class Gauge:
    def __init__(self, year, path, raw):
        self.year, self.path, self.raw = year, path, raw
        self.h = hourly(raw)


class RunYear:
    def __init__(self, stem, info, paths, notes):
        self.stem, self.info, self.paths = stem, info, paths
        default = pd.Timestamp(DEFAULT_START[info["set"]])
        start, src = default, "run name"
        self.runtime = None
        if paths.get("in"):
            kv = read_in(paths["in"])
            sd = parse_startdate(kv.get("STARTDATE", ""))
            if sd is not None:
                start, src = sd, ".in file"
                if sd != default:
                    notes.append("{}: STARTDATE in the .in file ({:%Y-%m-%d %H:%M}) differs from the name's default ({:%Y-%m-%d %H:%M}); the .in file is used.".format(stem, sd, default))
            try:
                self.runtime = float(kv.get("RUNTIME", ""))
            except ValueError:
                self.runtime = None
        self.start, self.start_source = start, src
        self.raw = read_qout(paths["qout"], start)
        self.h = hourly(self.raw)
        self.rain = read_mrf(paths["mrf"], start) if paths.get("mrf") else None
        self.rft = read_rft(paths["rft"], start) if paths.get("rft") else None
        if self.runtime is not None and len(self.raw):
            last_h = (self.raw.index[-1] - start) / HOUR
            if last_h < self.runtime - 1.0:
                notes.append("{}: the outlet file ends at hour {:.1f}, but RUNTIME in the .in file is {:g}: the run did not finish or the file is cut.".format(stem, last_h, self.runtime))


def find_companions(data: Path, outlet: Path, stem: str):
    pat = re.compile("^" + re.escape(stem) + r"(\d.*)?\.(mrf|rft)$")
    found = {"mrf": [], "rft": []}
    for folder in dict.fromkeys([outlet.parent, data]):
        it = folder.iterdir() if folder == outlet.parent else data.rglob("*")
        for p in it:
            m = pat.match(p.name)
            if m and p.is_file() and p not in found[m.group(2)]:
                found[m.group(2)].append(p)
        if found["mrf"] or found["rft"]:
            break
    ins = [p for p in data.rglob(stem + ".in") if p.is_file()]
    out = {"qout": outlet}
    for k in ("mrf", "rft"):
        if found[k]:
            out[k] = max(found[k], key=lambda p: p.stat().st_size)
    if ins:
        out["in"] = ins[0]
    return out


def discover(data: Path, old_map: dict, only=None):
    """{setting: {year: paths}} for every recognized run under data, plus the list of runs left out."""
    outs = sorted(p for p in data.rglob("*_Outlet.qout") if "_partial_" not in str(p))
    by_setting, ignored, seen = {}, [], set()
    for p in outs:
        stem = p.name[:-len("_Outlet.qout")]
        if stem in seen:
            continue
        seen.add(stem)
        info = classify(stem, old_map)
        if info is None or (only and info["setting"] not in only and info["base"] not in only):
            ignored.append(stem)
            continue
        by_setting.setdefault(info["setting"], {})[info["year"]] = (stem, info, find_companions(data, p, stem))
    return by_setting, ignored


def load_gauges(a, data: Path):
    cands = sorted(p for p in data.rglob("sycamore_09510200_obs_flow*.csv"))
    chosen = {}
    for year, opt, keys in ((2017, getattr(a, "gauge_2017", None), ("WY2017",)), (2020, getattr(a, "gauge_2020", None), ("WY2020", "Sep2019", "2020"))):
        if opt:
            fs = [Path(opt)]
        else:
            fs = [p for p in cands if any(k in p.name for k in keys)]
        if not fs:
            die("no gauge file found for WY{} under {} (looked for sycamore_09510200_obs_flow*{}*.csv). Use --gauge-{}.".format(year, data, keys[0], year))
        reads = [(p, read_gauge(p)) for p in fs]
        p0, s0 = reads[0]
        if len(reads) > 1:
            lo, hi = PERIODS[year]
            h0 = sel(hourly(s0), lo, hi)
            for p, s in reads[1:]:
                h1 = sel(hourly(s), lo, hi)
                both = h0.index.intersection(h1.index)
                dmax = float((h0[both] - h1[both]).abs().max()) if len(both) else float("nan")
                same_idx = len(h0.index.symmetric_difference(h1.index)) == 0
                log("  gauge files for WY{}: {} and {} {} over Oct-May (largest hourly difference {:.2g} m3/s; hours only in one: {}). Using {}.".format(
                    year, p0.name, p.name, "agree" if (isnum(dmax) and dmax < 1e-6 and same_idx) else "DIFFER", dmax,
                    len(h0.index.symmetric_difference(h1.index)), p0.name))
        chosen[year] = Gauge(year, p0, s0)
    return chosen


# ================================================================================================ the rain windows
def find_clusters(rain: pd.Series, thr: float, gap_h: int, lo, hi) -> pd.DataFrame:
    """Storms in the hourly basin-mean rain: hours at or above thr, joined while fewer than gap_h dry hours lie between them."""
    cols = ["start", "end", "total_mm", "peak_mm_h", "wet_hours"]
    r = sel(rain, lo, hi).dropna()
    wet = r[r >= thr]
    if wet.empty:
        return pd.DataFrame(columns=cols)
    t = wet.index
    secs = t.values.astype("datetime64[s]").astype("int64")
    dry_between = np.r_[np.inf, np.diff(secs) / 3600.0 - 1.0]
    gid = np.cumsum(dry_between >= gap_h)
    out = []
    for g in np.unique(gid):
        tt = t[gid == g]
        start, end = tt[0], tt[-1] + HOUR
        seg = r[(r.index >= start) & (r.index < end)]
        out.append((start, end, float(seg.sum()), float(seg.max()), int((seg >= thr).sum())))
    return pd.DataFrame(out, columns=cols)


def build_storm_windows(rain_by_year, thr, gap_h, R_days, min_next_mm):
    rows, warns = [], []
    for year in (2017, 2020):
        p0, p1 = PERIODS[year]
        cl = find_clusters(rain_by_year[year], thr, gap_h, p0, p1)
        for (_, name, f, l, mm_tab, named) in [s for s in STORM_TABLE if s[0] == year and s[5]]:
            r0, r1 = pd.Timestamp(f), pd.Timestamp(l) + DAY
            idx = [i for i, c in cl.iterrows() if c["start"] < r1 and c["end"] > r0]
            if not idx:
                warns.append("{} {}: no rain cluster found in the named days ({} to {}); no window".format(year, name, f, l))
                continue
            mine = cl.loc[idx]
            start, last_end = mine["start"].min(), mine["end"].max()
            others = cl[(~cl.index.isin(idx)) & (cl["total_mm"] >= min_next_mm)]
            inside = others[(others["start"] >= start) & (others["start"] < last_end)]
            if len(inside):
                warns.append("{} {}: a separate rain cluster of {:.1f} mm starts inside the storm ({:%b %d %H:%M}); it is not treated as a next storm".format(
                    year, name, float(inside["total_mm"].iloc[0]), inside["start"].iloc[0]))
            later = others[others["start"] >= last_end]
            nxt = later["start"].min() if len(later) else pd.NaT
            nxt_mm = float(later.loc[later["start"].idxmin(), "total_mm"]) if len(later) else np.nan
            end_r = last_end + R_days * DAY
            end = min(end_r, nxt) if pd.notna(nxt) else end_r
            end = min(end, p1)
            rows.append({"year": year, "name": name, "named_first_day": f, "named_last_day": l, "n_clusters_merged": len(idx),
                         "rain_start": start, "rain_last_end": last_end, "cluster_rain_mm": float(mine["total_mm"].sum()),
                         "table_basin_mm": mm_tab, "peak_rain_mm_h": float(mine["peak_mm_h"].max()),
                         "window_start": start, "window_end": end, "window_hours": int(round((end - start) / HOUR)),
                         "cut_by_next_storm": bool(pd.notna(nxt) and nxt < end_r), "next_storm_start": nxt, "next_storm_mm": nxt_mm,
                         "R_days": R_days})
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("window_start").reset_index(drop=True)
        df.insert(0, "storm_id", ["storm{:02d}".format(i + 1) for i in range(len(df))])
        for i in range(len(df) - 1):
            if df.loc[i, "window_end"] > df.loc[i + 1, "window_start"]:
                warns.append("windows of {} and {} overlap".format(df.loc[i, "name"], df.loc[i + 1, "name"]))
    return df, warns


def evaluate_params(rain_by_year, thr, gap_h, min_next_mm):
    named_clean = named_split = named_missing = small_found = small_total = unlisted = n_next = 0
    for year in (2017, 2020):
        p0, p1 = PERIODS[year]
        cl = find_clusters(rain_by_year[year], thr, gap_h, p0, p1)
        big = cl[cl["total_mm"] >= min_next_mm]
        n_next += len(big)
        listed = []
        for (_, name, f, l, mm, named) in [s for s in STORM_TABLE if s[0] == year]:
            r0, r1 = pd.Timestamp(f), pd.Timestamp(l) + DAY
            ov = cl[(cl["start"] < r1) & (cl["end"] > r0)]
            ovb = big[(big["start"] < r1) & (big["end"] > r0)]
            listed.append((r0, r1))
            if named:
                named_clean += int(len(ov) == 1)
                named_split += int(len(ov) > 1)
                named_missing += int(len(ov) == 0)
            else:
                small_total += 1
                small_found += int(len(ovb) >= 1)
        for _, c in big.iterrows():
            if not any(c["start"] < r1 and c["end"] > r0 for r0, r1 in listed):
                unlisted += 1
    return {"thr_mm_h": thr, "gap_h": gap_h, "named_one_cluster": named_clean, "named_split": named_split, "named_missing": named_missing,
            "table_small_storms_found": "{}/{}".format(small_found, small_total), "clusters_ge_min_next": n_next, "unlisted_clusters_ge_min_next": unlisted}


def gauge_hourly_stats(gh: pd.Series, lo, hi):
    h = sel(gh, lo, hi)
    return h, vol_mm3(h.to_numpy())


def r_test(storms: pd.DataFrame, gauges: dict, Rs=R_TEST):
    rows = []
    for _, s in storms.iterrows():
        year = int(s["year"])
        p1 = PERIODS[year][1]
        gh = gauges[year].h
        _, v7 = gauge_hourly_stats(gh, s["window_start"], min(s["rain_last_end"] + max(Rs) * DAY, s["next_storm_start"] if pd.notna(s["next_storm_start"]) else p1, p1))
        base = {"storm_id": s["storm_id"], "year": year, "name": s["name"]}
        gall = sel(gh, s["window_start"], min(s["rain_last_end"] + max(Rs) * DAY, p1))
        pk_h = float(gall.max()) if len(gall) else np.nan
        pk_t = gall.idxmax() if len(gall) else pd.NaT
        # hours after the last rain until the gauge hourly mean first falls to 5% of its peak (or 1 cfs, whichever is higher)
        recede = np.nan
        if len(gall) and pk_h / CFS >= 1.0:
            lim = max(0.05 * pk_h, 1.0 * CFS)
            nxt = s["next_storm_start"] if pd.notna(s["next_storm_start"]) else p1
            after = gall[(gall.index >= pk_t) & (gall.index < nxt)]
            low = after[after <= lim]
            recede = float((low.index[0] - s["rain_last_end"]) / HOUR) if len(low) else np.nan
        for R in Rs:
            end = min(s["rain_last_end"] + R * DAY, s["next_storm_start"] if pd.notna(s["next_storm_start"]) else p1, p1)
            hw, vol = gauge_hourly_stats(gh, s["window_start"], end)
            edge = float(hw.iloc[-1]) if len(hw) else np.nan
            rows.append({**base, "R_days": R, "window_end": end, "cut_by_next_storm": bool(pd.notna(s["next_storm_start"]) and s["next_storm_start"] < s["rain_last_end"] + R * DAY),
                         "gauge_vol_Mm3": vol, "gauge_hours": int(len(hw)), "window_hours": int(round((end - s["window_start"]) / HOUR)),
                         "pct_of_longest_window": 100.0 * vol / v7 if v7 > 0 else np.nan,
                         "gauge_peak_hourly_cfs": pk_h / CFS if isnum(pk_h) else np.nan,
                         "edge_flow_cfs": edge / CFS if isnum(edge) else np.nan,
                         "edge_pct_of_peak": 100.0 * edge / pk_h if (isnum(edge) and isnum(pk_h) and pk_h > 0) else np.nan,
                         "hours_after_last_rain_to_5pct_of_peak": recede})
    return pd.DataFrame(rows)


def old_rule_gauge_check(gauges: dict):
    """Gauge volume over the OLD window rule for the 12 storms, against the numbers on record. A check on the gauge FILES only."""
    out = []
    for year in (2017, 2020):
        tab = [s for s in STORM_TABLE if s[0] == year]
        for i, (_, name, f, l, mm, named) in enumerate(tab):
            if not named:
                continue
            lo = pd.Timestamp(f)
            hi = pd.Timestamp(l) + 4 * DAY
            if year == 2017 and i + 1 < len(tab):
                hi = min(hi, pd.Timestamp(tab[i + 1][2]))
            h = sel(gauges[year].h, lo, hi)
            v = vol_mm3(h.to_numpy())
            exp = GAUGE_CHECK[(year, name)]
            out.append((year, name, lo, hi, v, exp, abs(v - exp) <= GAUGE_CHECK_TOL, len(h), int(round((hi - lo) / HOUR))))
    return out


def load_rain(year, by_setting, notes):
    """Basin-mean hourly rain for a year from one run's .mrf; checks that the other runs of that year have the same rain."""
    pref = ["Retuned", "Original"] + sorted(by_setting, key=setting_sort_key)
    ref, ref_stem = None, None
    others = []
    for st in pref:
        ent = by_setting.get(st, {}).get(year)
        if not ent or "mrf" not in ent[2] or ent[1]["set"] == "WY2017old":
            continue
        stem, info, paths = ent
        run_start = DEFAULT_START[info["set"]]
        if paths.get("in"):
            sd = parse_startdate(read_in(paths["in"]).get("STARTDATE", ""))
            run_start = sd if sd is not None else run_start
        s = read_mrf(paths["mrf"], pd.Timestamp(run_start))
        if ref is None:
            ref, ref_stem = s, stem
        else:
            others.append((stem, s))
    if ref is None:
        # last resort: an old run
        for st in by_setting:
            ent = by_setting[st].get(year)
            if ent and "mrf" in ent[2]:
                stem, info, paths = ent
                run_start = pd.Timestamp(DEFAULT_START[info["set"]])
                ref, ref_stem = read_mrf(paths["mrf"], run_start), stem
                break
    if ref is None:
        die("no .mrf file found for WY{}: the storm windows need the basin-mean rain from one run's .mrf file.".format(year))
    worst = 0.0
    for stem, s in others:
        lo, hi = PERIODS[year]
        a, b = sel(ref, lo, hi), sel(s, lo, hi)
        idx = a.index.intersection(b.index)
        d = float((a[idx] - b[idx]).abs().max()) if len(idx) else float("nan")
        if len(a.index.symmetric_difference(b.index)):
            notes.append("WY{}: the rain of {} covers different hours than {}".format(year, stem, ref_stem))
        if isnum(d) and d > 1e-6:
            notes.append("WY{}: the rain in {} differs from {} by up to {:.3g} mm/h".format(year, stem, ref_stem, d))
        worst = max(worst, d if isnum(d) else 0.0)
    return ref, ref_stem, len(others), worst


# ================================================================================================ the windows command
def cmd_windows(a) -> int:
    data = Path(a.data)
    out = Path(a.out_dir)
    if not data.is_dir():
        die("--data folder not found: {}".format(data))
    notes = []
    log("\nscore_windows_sycamore.py v{}   windows command (rain and gauge only; no model result is used except the rain in the .mrf)".format(VERSION))
    by_setting, ignored = discover(data, OLD_STEMS)
    gauges = load_gauges(a, data)
    for y, g in gauges.items():
        log("  WY{} gauge: {} ({} readings, {:%b %d %Y} to {:%b %d %Y})".format(y, g.path.name, int(g.raw.notna().sum()), g.raw.index.min(), g.raw.index.max()))
    rain, info = {}, {}
    for year in (2017, 2020):
        rain[year], stem, n_other, worst = load_rain(year, by_setting, notes)
        info[year] = (stem, n_other, worst)
        log("  WY{} rain: from {}; {} other run(s) compared, largest difference {:.3g} mm/h".format(year, stem, n_other, worst))
    for n in notes:
        log("  NOTE: " + n)

    log("\n1. THE GAUGE FILES against the volumes on record (old window rule: 00:00 on the first rain day to 00:00 four days after the last rain day;")
    log("   WY2017 cut at the next storm in the storm table). This checks the FILES, not the new windows.")
    chk = old_rule_gauge_check(gauges)
    for (y, name, lo, hi, v, exp, ok, nh, nw) in chk:
        log("   {} {:<13} {:%b %d} to {:%b %d}  gauge {:7.3f}  on record {:7.3f}  {}{}".format(
            y, name, lo, hi, v, exp, "MATCH" if ok else "MISMATCH", "" if nh >= nw else "   (gauge has {} of {} hours)".format(nh, nw)))
    nbad = sum(1 for c in chk if not c[6])
    log("   -> {}".format("all 12 match" if nbad == 0 else "{} of 12 do NOT match: a gauge file or a window rule differs from the earlier scoring. Do not use the results until this is explained.".format(nbad)))

    log("\n2. RAIN THRESHOLD AND DRY GAP, tested on the 12 named storms (rain only). A storm is clean when exactly ONE rain cluster covers its named days.")
    log("   MIN_NEXT (smallest rain cluster that counts as a 'next storm') = {:g} mm".format(a.min_next))
    grid = pd.DataFrame([evaluate_params(rain, t, g, a.min_next) for t in GRID_THR for g in GRID_GAP])
    log(grid.to_string(index=False))

    log("\n3. THE 12 STORM WINDOWS with THR = {:g} mm/h, GAP = {} h, MIN_NEXT = {:g} mm, R = {} days".format(a.thr, a.gap, a.min_next, a.R))
    storms, warns = build_storm_windows(rain, a.thr, a.gap, a.R, a.min_next)
    for w in warns:
        log("   WARNING: " + w)
    show = storms[["storm_id", "year", "name", "n_clusters_merged", "rain_start", "rain_last_end", "cluster_rain_mm", "table_basin_mm", "window_end", "window_hours", "cut_by_next_storm"]].copy()
    log(show.to_string(index=False, float_format=lambda x: "{:.1f}".format(x)))
    if len(storms) != 12:
        log("   WARNING: {} of 12 storms have a window".format(len(storms)))

    log("\n4. HOW MANY DAYS AFTER THE LAST RAIN (R)? Gauge and rain only; windows are cut at the next storm in every case.")
    rt = r_test(storms, gauges)
    summ = []
    for R in R_TEST:
        sub = rt[(rt["R_days"] == R) & (rt["gauge_vol_Mm3"] >= MIN_GAUGE_VOL)]
        summ.append({"R_days": R, "storms_cut_by_next": int(rt[rt["R_days"] == R]["cut_by_next_storm"].sum()),
                     "median_pct_of_longest_window": float(sub["pct_of_longest_window"].median()) if len(sub) else np.nan,
                     "storms_still_above_10pct_of_peak_at_end": int((sub["edge_pct_of_peak"] > 10).sum()),
                     "storms_still_above_25pct_of_peak_at_end": int((sub["edge_pct_of_peak"] > 25).sum()), "storms_counted": len(sub)})
    log(pd.DataFrame(summ).to_string(index=False, float_format=lambda x: "{:.1f}".format(x)))
    log("   (only storms with a gauge volume of at least {:g} M m3 are counted; 'longest window' = R {} days, cut at the next storm)".format(MIN_GAUGE_VOL, max(R_TEST)))
    pv = rt[rt["R_days"].isin([3, 4, 5])].pivot(index="name", columns="R_days", values=["gauge_vol_Mm3", "edge_pct_of_peak"])
    pv.columns = ["{}_R{}".format(c[0].replace("gauge_vol_Mm3", "vol").replace("edge_pct_of_peak", "edge%peak"), c[1]) for c in pv.columns]
    first = rt[rt["R_days"] == 4].set_index("name")["hours_after_last_rain_to_5pct_of_peak"].rename("h_to_5%peak")
    log(pv.join(first).to_string(float_format=lambda x: "{:.2f}".format(x)))

    write_csv(out / "windows_grid_test.csv", ["THR (mm/h) and GAP (dry hours) tested on the 12 named storms; MIN_NEXT = {:g} mm".format(a.min_next),
                                              "named_one_cluster = storms covered by exactly one rain cluster; named_split = more than one; named_missing = none"], grid)
    write_csv(out / "windows_R_test.csv", ["gauge volume and right-edge flow of every storm window for R = 3, 4, 5 and 7 days after the last rain (cut at the next storm)",
                                           "edge_flow_cfs = last gauge hourly mean in the window; edge_pct_of_peak = that as a percent of the gauge hourly peak"], rt)
    hdr = ["storm windows found from the basin-mean hourly rain in the .mrf file (times are Arizona local time; window = [window_start, window_end))",
           "PARAMETERS thr_mm_h={:g} gap_h={} min_next_mm={:g} R_days={}".format(a.thr, a.gap, a.min_next, a.R),
           "window_start = first hour with rain >= thr after at least gap_h dry hours; window_end = end of the last rain hour + R_days, or the start of the next storm "
           "(a rain cluster of at least min_next_mm), whichever is first",
           "n_clusters_merged = rain clusters inside the storm's named days that were joined into one storm"]
    write_csv(out / "storm_windows.csv", hdr, storms)
    log("\nWROTE {}, {}, {}".format(out / "storm_windows.csv", out / "windows_grid_test.csv", out / "windows_R_test.csv"))
    return 1 if (len(storms) != 12 or nbad) else 0


# ================================================================================================ windows for scoring
class Window:
    def __init__(self, wid, kind, parts, name=""):
        self.wid, self.kind, self.parts, self.name = wid, kind, parts, name or wid

    @property
    def years(self):
        return sorted({y for y, _, _ in self.parts})


def complement(period, intervals):
    p0, p1 = period
    cur, out = p0, []
    for lo, hi in sorted(intervals):
        if lo > cur:
            out.append((cur, lo))
        cur = max(cur, hi)
    if cur < p1:
        out.append((cur, p1))
    return out


def build_windows(storms: pd.DataFrame):
    W = [Window("OctMay2017", "season", [(2017,) + PERIODS[2017]]), Window("OctMay2020", "season", [(2020,) + PERIODS[2020]]),
         Window("Both", "pooled", [(2017,) + PERIODS[2017], (2020,) + PERIODS[2020]])]
    st_parts = {2017: [], 2020: []}
    for _, s in storms.iterrows():
        y = int(s["year"])
        W.append(Window("{}_{}_{}".format(s["storm_id"], y, re.sub(r"[^A-Za-z0-9-]+", "", s["name"])), "storm", [(y, s["window_start"], s["window_end"])], s["name"]))
        st_parts[y].append((y, s["window_start"], s["window_end"]))
    W.append(Window("AllStorms", "all_storms", st_parts[2017] + st_parts[2020]))
    for y in (2017, 2020):
        W.append(Window("AllStorms_{}".format(y), "all_storms", st_parts[y]))
    off = {y: [(y, lo, hi) for lo, hi in complement(PERIODS[y], [(p[1], p[2]) for p in st_parts[y]])] for y in (2017, 2020)}
    W.append(Window("OFF_2017", "off", off[2017]))
    W.append(Window("OFF_2020", "off", off[2020]))
    W.append(Window("OFF_both", "off", off[2017] + off[2020]))
    W.append(Window("Oct2016", "season", [(2017, pd.Timestamp("2016-10-01"), pd.Timestamp("2016-11-01"))], "Oct 1 - Nov 1 2016 (spin-up check)"))
    W.append(Window("OctDec2016", "season", [(2017, pd.Timestamp("2016-10-01"), pd.Timestamp("2017-01-01"))], "Oct 1 - Jan 1 2016-17 (spin-up check)"))
    return W


# ================================================================================================ metric functions
def kge_prime(m: np.ndarray, g: np.ndarray):
    out = {"kge": np.nan, "r": np.nan, "beta": np.nan, "gamma": np.nan, "why": ""}
    if len(g) < 3:
        out["why"] = "fewer than 3 matched values"
        return out
    mg, gg = float(np.mean(m)), float(np.mean(g))
    if gg <= 0 or float(np.std(g)) == 0:
        out["why"] = "gauge flow is zero or constant"
        return out
    out["beta"] = mg / gg
    if float(np.std(m)) == 0 or mg == 0:
        out["why"] = "model flow is zero or constant (r and gamma undefined)"
        return out
    out["r"] = float(np.corrcoef(m, g)[0, 1])
    out["gamma"] = (float(np.std(m)) / mg) / (float(np.std(g)) / gg)
    out["kge"] = float(1 - np.sqrt((out["r"] - 1) ** 2 + (out["beta"] - 1) ** 2 + (out["gamma"] - 1) ** 2))
    return out


def cum_times(h: pd.Series, lo, fracs=(0.1, 0.5, 0.9)):
    """Hours after `lo` by which each fraction of the volume of the hourly means h has passed (linear within the hour)."""
    v = h.to_numpy(float) * 3600.0
    tot = v.sum()
    if not np.isfinite(tot) or tot <= 0:
        return None
    cum = np.cumsum(v)
    starts = ((h.index - lo) / HOUR).to_numpy(float)
    out = {}
    for f in fracs:
        i = min(int(np.searchsorted(cum, f * tot, side="left")), len(v) - 1)
        prev = cum[i - 1] if i > 0 else 0.0
        out[f] = float(starts[i] + ((f * tot - prev) / v[i] if v[i] > 0 else 0.0))
    return out


def arrival(raw: pd.Series, level: float):
    """First time the flow RISES THROUGH level (linear between readings). Returns (Timestamp or None, reason)."""
    s = raw.dropna()
    if s.empty:
        return None, "has no data"
    above = s.to_numpy() >= level
    rising = np.nonzero(above[1:] & ~above[:-1])[0] + 1
    if len(rising) == 0:
        return None, ("is already above the level at the window start" if above[0] else "never reaches the level")
    i = int(rising[0])
    t0, t1 = s.index[i - 1], s.index[i]
    q0, q1 = float(s.iloc[i - 1]), float(s.iloc[i])
    if (t1 - t0) > pd.Timedelta(hours=MAX_STEP_H) or q1 == q0:
        return t1, ""
    return t0 + (t1 - t0) * ((level - q0) / (q1 - q0)), ""


def best_lag(gh: pd.Series, mh: pd.Series):
    """Whole-hour shift within +/-MAX_LAG_H at which model hourly means correlate best with the gauge's. Positive = model later."""
    best, best_r = None, -np.inf
    for k in sorted(range(-MAX_LAG_H, MAX_LAG_H + 1), key=lambda x: (abs(x), x)):
        gs = gh.copy()
        gs.index = gs.index + k * HOUR
        df = pd.concat([mh.rename("m"), gs.rename("g")], axis=1, join="inner").dropna()
        if len(df) < MIN_COMMON_H or df["m"].std() == 0 or df["g"].std() == 0:
            continue
        r = float(df["m"].corr(df["g"]))
        if r > best_r + 1e-12:
            best, best_r = k, r
    return best, (best_r if best is not None else np.nan)


def interp_15min(model_raw: pd.Series, gauge_raw: pd.Series, lo, hi):
    g = sel(gauge_raw, lo, hi).dropna()
    mt = model_raw.index.values.astype("datetime64[s]").astype("int64").astype(float)
    gt = g.index.values.astype("datetime64[s]").astype("int64").astype(float)
    m = np.interp(gt, mt, model_raw.to_numpy(float), left=np.nan, right=np.nan)
    ok = np.isfinite(m)
    return m[ok], g.to_numpy(float)[ok]


def storm_extra(run, gauge, lo, hi, gh_m, mh_m, mh_win, gh_win, cls, vol_ratio):
    """Peak and timing for one storm window. Returns (values, reasons)."""
    out, why = {}, []
    gw = sel(gauge.raw, lo, hi).dropna()
    mw = sel(run.raw, lo, hi).dropna()
    if len(gw):
        out["peak_gauge_cfs"] = float(gw.max() / CFS)
    if len(mw):
        out["peak_model_cfs"] = float(mw.max() / CFS)
    if "peak_gauge_cfs" in out and "peak_model_cfs" in out:
        out["peak_diff_cfs"] = out["peak_model_cfs"] - out["peak_gauge_cfs"]
    if cls != "ok":
        why.append("peak_ratio, dt_peak_h, dt10/50/90_h, arr20/50_diff_h, lag_h, tail_vol_ratio: gauge class is {} (no gauge hydrograph to compare with)".format(cls))
        return out, why
    out["peak_ratio"] = out["peak_model_cfs"] / out["peak_gauge_cfs"] if out.get("peak_gauge_cfs", 0) > 0 and "peak_model_cfs" in out else np.nan
    if len(gw) and len(mw):
        out["dt_peak_h"] = float((mw.idxmax() - gw.idxmax()) / HOUR)
    tg, tm = cum_times(gh_m, lo), cum_times(mh_m, lo)
    if tg is None:
        why.append("dt10/50/90_h, tail_vol_ratio: gauge volume is zero")
    elif tm is None:
        why.append("dt10/50/90_h: the model made no flow in this window (nothing to time)")
    else:
        for f, c in ((0.1, "dt10_h"), (0.5, "dt50_h"), (0.9, "dt90_h")):
            out[c] = tm[f] - tg[f]
        t0 = lo + tg[0.9] * HOUR
        gtail = vol_mm3(gh_m[gh_m.index >= t0].to_numpy())
        mtail = vol_mm3(mh_m[mh_m.index >= t0].to_numpy())
        if gtail > 0:
            out["tail_vol_ratio"] = mtail / gtail
        else:
            why.append("tail_vol_ratio: the gauge has no volume after its t90")
    if len(gw):
        gpk = float(gw.max())
        for frac, key in ((0.2, "arr20_diff_h"), (0.5, "arr50_diff_h")):
            tg_, wg = arrival(gw, frac * gpk)
            tm_, wm = arrival(mw, frac * gpk)
            if tg_ is None:
                why.append("{}: the gauge {}".format(key, wg))
            elif tm_ is None:
                why.append("{}: the model {} ({:.0f}% of the gauge peak)".format(key, wm, 100 * frac))
            else:
                out[key] = float((tm_ - tg_) / HOUR)
    lag, lr = best_lag(gh_win, mh_win)
    if lag is None:
        why.append("lag_h: fewer than {} common hours or a flat series".format(MIN_COMMON_H))
    else:
        out["lag_h"] = float(lag)
    if isnum(vol_ratio) and vol_ratio < LOW_MODEL_RATIO:
        why.append("NOTE (not a blank): model volume is only {:.0%} of the gauge's, so the timing columns are not a like-for-like comparison".format(vol_ratio))
    return out, why


def score_window(win: Window, runs: dict, gauges: dict, basis: str):
    """One row of metrics for one setting in one window. runs = {year: RunYear}. None if the setting lacks a needed year."""
    yrs = win.years
    if any(y not in runs for y in yrs):
        return None
    row = {"window": win.wid, "kind": win.kind, "name": win.name, "years": "+".join(str(y) for y in yrs), "runs": " + ".join(runs[y].stem for y in yrs)}
    if len(win.parts) == 1:
        row["start"], row["end"] = win.parts[0][1], win.parts[0][2]
    M, G, Mk, Gk = [], [], [], []
    n_hours = gauge_hours = model_hours = 0
    vol_model_all = 0.0
    single = None
    for (y, lo, hi) in win.parts:
        run, g = runs[y], gauges[y]
        n_hours += int(round((hi - lo) / HOUR))
        gh, mh_all = sel(g.h, lo, hi), sel(run.h, lo, hi)
        gauge_hours += len(gh)
        model_hours += len(mh_all)
        vol_model_all += vol_mm3(mh_all.to_numpy())
        mh = run.h.reindex(gh.index)
        ok = mh.notna()
        M.append(mh[ok].to_numpy(float))
        G.append(gh[ok].to_numpy(float))
        if basis == "15min":
            a_, b_ = interp_15min(run.raw, g.raw, lo, hi)
            Mk.append(a_)
            Gk.append(b_)
        if len(win.parts) == 1:
            single = (run, g, lo, hi, gh[ok], mh[ok], mh_all, gh)
    m = np.concatenate(M) if M else np.array([])
    g_ = np.concatenate(G) if G else np.array([])
    vg, vm = vol_mm3(g_), vol_mm3(m)
    cls = gauge_class(vg, len(g_))
    row.update({"n_hours": n_hours, "gauge_hours": gauge_hours, "model_hours": model_hours, "matched_hours": int(len(g_)), "gauge_class": cls,
                "vol_gauge_Mm3": vg if len(g_) else np.nan, "vol_model_Mm3": vm if len(g_) else np.nan,
                "vol_diff_Mm3": (vm - vg) if len(g_) else np.nan, "vol_model_allhours_Mm3": vol_model_all})
    why, notes = [], []
    if gauge_hours < n_hours:
        notes.append("gauge has {} of {} hours".format(gauge_hours, n_hours))
    if model_hours < n_hours:
        notes.append("model has {} of {} hours".format(model_hours, n_hours))
    if cls == "ok":
        row["pbias_pct"] = 100.0 * (vm - vg) / vg
    else:
        why.append("pbias_pct: gauge class is {} (gauge volume {} M m3)".format(cls, fnum(vg, 3)))
    if win.kind != "off":
        if cls == "ok":
            if basis == "15min":
                k = kge_prime(np.concatenate(Mk), np.concatenate(Gk)) if Mk else kge_prime(np.array([]), np.array([]))
            else:
                k = kge_prime(m, g_)
            row.update({"kge_prime": k["kge"], "kge_r": k["r"], "kge_beta": k["beta"], "kge_gamma": k["gamma"]})
            if k["why"]:
                why.append("kge_prime: " + k["why"])
        else:
            why.append("kge_prime, r, beta, gamma: gauge class is {}".format(cls))
    if win.kind == "storm" and single is not None:
        run, g, lo, hi, gh_m, mh_m, mh_all, gh_win = single
        ex, w2 = storm_extra(run, g, lo, hi, gh_m, mh_m, mh_all, gh_win, cls, (vm / vg) if (cls == "ok" and vg > 0) else np.nan)
        row.update(ex)
        why += w2
    if win.kind in ("season", "pooled", "off"):
        for thr in DAYS_ABOVE_CFS:
            dm = dg = 0
            for (y, lo, hi) in win.parts:
                sm, sg = sel(runs[y].raw, lo, hi).dropna(), sel(gauges[y].raw, lo, hi).dropna()
                dm += int((sm.resample("D").max().dropna() > thr * CFS).sum()) if len(sm) else 0
                dg += int((sg.resample("D").max().dropna() > thr * CFS).sum()) if len(sg) else 0
            row["days_model_gt{:g}cfs".format(thr)], row["days_gauge_gt{:g}cfs".format(thr)] = dm, dg
        row["model_max_cfs"] = float(max((sel(runs[y].raw, lo, hi).max() for (y, lo, hi) in win.parts), default=np.nan) / CFS)
        row["gauge_max_cfs"] = float(max((sel(gauges[y].raw, lo, hi).max() for (y, lo, hi) in win.parts), default=np.nan) / CFS)
    row["notes"] = " | ".join(notes)
    row["blank_reasons"] = " | ".join(why)
    return row


MEDIAN_COLS = ["peak_ratio", "dt_peak_h", "dt10_h", "dt50_h", "dt90_h", "arr20_diff_h", "arr50_diff_h", "lag_h", "tail_vol_ratio"]


def fill_all_storm_medians(rows_by_window: dict, win: Window, storm_rows: list):
    row = rows_by_window.get(win.wid)
    if row is None:
        return
    okr = [r for r in storm_rows if r.get("gauge_class") == "ok"]
    counts = []
    for c in MEDIAN_COLS:
        vals = [r[c] for r in okr if isnum(r.get(c))]
        if vals:
            row[c] = float(np.median(vals))
        counts.append("{}:{}".format(c, len(vals)))
    pd_vals = [r["peak_diff_cfs"] for r in storm_rows if isnum(r.get("peak_diff_cfs"))]
    if pd_vals:
        row["peak_diff_cfs"] = float(np.median(pd_vals))
    row["n_ok_storms"] = len(okr)
    row["n_storms"] = len(storm_rows)
    row["notes"] = (row.get("notes", "") + " | " if row.get("notes") else "") + "peak and timing columns are MEDIANS across the storms with gauge class ok; storms used per column = " + ", ".join(counts)


def mech_row(win: Window, runs: dict):
    yrs = win.years
    if any(y not in runs for y in yrs):
        return None
    row = {"window": win.wid, "kind": win.kind, "name": win.name, "years": "+".join(str(y) for y in yrs), "runs": " + ".join(runs[y].stem for y in yrs)}
    tot = {t: 0.0 for t in MECH_TYPES}
    rain, have_rft, have_rain = 0.0, True, True
    for (y, lo, hi) in win.parts:
        run = runs[y]
        if run.rft is None:
            have_rft = False
        else:
            r = sel(run.rft, lo, hi)
            for t in MECH_TYPES:
                tot[t] += float(r["{}_m3_s".format(t)].sum() * 3600.0 / 1e6)
        if run.rain is None:
            have_rain = False
        else:
            rain += float(sel(run.rain, lo, hi).sum())
    row["rain_mm"] = rain if have_rain else np.nan
    row["rain_Mm3"] = rain * BASIN_KM2 / 1000.0 if have_rain else np.nan
    if have_rft:
        total = sum(tot.values())
        for t in MECH_TYPES:
            row["{}_Mm3".format(t)] = tot[t]
        row["made_total_Mm3"] = total
        for t in MECH_TYPES:
            row["{}_pct".format(t)] = 100.0 * tot[t] / total if total > 0 else np.nan
        row["made_over_rain_pct"] = 100.0 * total / row["rain_Mm3"] if have_rain and row["rain_Mm3"] > 0 else np.nan
    else:
        row["notes"] = "no .rft file for this run"
    return row


# ================================================================================================ the score command
SCORE_COLS = ["setting", "window", "kind", "name", "years", "runs", "start", "end", "n_hours", "gauge_hours", "model_hours", "matched_hours", "gauge_class",
              "vol_gauge_Mm3", "vol_model_Mm3", "vol_diff_Mm3", "vol_model_allhours_Mm3", "pbias_pct", "kge_prime", "kge_r", "kge_beta", "kge_gamma",
              "peak_gauge_cfs", "peak_model_cfs", "peak_diff_cfs", "peak_ratio", "dt_peak_h", "dt10_h", "dt50_h", "dt90_h", "arr20_diff_h", "arr50_diff_h",
              "lag_h", "tail_vol_ratio", "days_model_gt1cfs", "days_gauge_gt1cfs", "days_model_gt10cfs", "days_gauge_gt10cfs", "model_max_cfs", "gauge_max_cfs",
              "n_ok_storms", "n_storms", "notes", "blank_reasons"]
MECH_COLS = ["setting", "window", "kind", "name", "years", "runs", "rain_mm", "rain_Mm3", "Hsrf_Mm3", "Sbsrf_Mm3", "Psrf_Mm3", "Satsrf_Mm3", "made_total_Mm3",
             "Hsrf_pct", "Sbsrf_pct", "Psrf_pct", "Satsrf_pct", "made_over_rain_pct", "notes"]

DEF_LINES = [
    "Arizona local time. Hourly means (matched hours: both gauge and model have the hour). Volume = sum of hourly means x 3,600 s in M m3.",
    "pbias_pct = 100 x (model - gauge) / gauge: MODEL MINUS GAUGE, POSITIVE = MODEL TOO HIGH (opposite of the textbook Gupta sign). It is blank when gauge_class is not ok.",
    "kge_prime = Kling-Gupta efficiency 2012 = 1 - sqrt((r-1)^2 + (beta-1)^2 + (gamma-1)^2); beta = model mean / gauge mean; gamma = CV ratio. Not gated: read it with beta.",
    "gauge_class: zero = gauge volume < 0.005 M m3; near-zero = < 0.05 M m3; ok = 0.05 M m3 or more. Ratio-type values are blank (reason in blank_reasons) unless ok.",
    "peak_*_cfs = highest single reading in the storm window (cfs); dt_* and arr*_diff_h and lag_h = hours, model minus gauge (positive = model later).",
    "dt10/50/90_h = time by which 10/50/90% of window volume passed; arr20/50_diff_h = first rise through 20/50% of the GAUGE peak (blank = model never reaches it); tail_vol_ratio = model/gauge volume after the gauge t90.",
    "AllStorms rows: pbias and KGE' pooled over the storm windows; peak and timing columns are medians over storms with gauge_class ok. OFF rows: no KGE', peak or timing.",
    "days_*_gt1cfs / gt10cfs = number of days whose highest flow is above 1 / 10 cfs (gauge: days with data only).",
    "Noise floor: the model is written every 3.75 min and the gauge every 15 min, so even identical flow gives dt10/50/90 differences of about 0.1 h (6 min). Treat timing differences under about 0.25 h as zero.",
]


def fmt_df(df: pd.DataFrame, cols, widths=None, nd=2) -> str:
    d = df[[c for c in cols if c in df.columns]].copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].map(lambda x: "-" if not np.isfinite(x) else ("{:.%df}" % nd).format(x))
    return d.to_string(index=False)


def spinup_pairs(df: pd.DataFrame, mech: pd.DataFrame):
    rows = []
    wins = [w for w in df["window"].unique() if (w.startswith("OctMay2017") or w.startswith("storm") or w in ("AllStorms_2017", "OFF_2017", "Oct2016", "OctDec2016"))]
    wins = [w for w in wins if (w.startswith("storm") and "_2017_" in w) or not w.startswith("storm")]
    metrics = [("vol_model_allhours_Mm3", "model volume, every hour (M m3)"), ("pbias_pct", "PBIAS (%)"), ("kge_prime", "KGE'"), ("peak_model_cfs", "model peak (cfs)"),
               ("dt_peak_h", "peak time, model - gauge (h)"), ("dt50_h", "t50, model - gauge (h)"), ("model_max_cfs", "model highest flow (cfs)"),
               ("days_model_gt1cfs", "days with model flow above 1 cfs")]
    for base in ("Retuned", "Original"):
        new, old = base, base + "_old"
        for w in wins:
            a = df[(df["setting"] == new) & (df["window"] == w)]
            b = df[(df["setting"] == old) & (df["window"] == w)]
            if a.empty or b.empty:
                continue
            for col, label in metrics:
                if col not in df.columns:
                    continue
                x, y = a.iloc[0].get(col, np.nan), b.iloc[0].get(col, np.nan)
                if isnum(x) or isnum(y):
                    rows.append({"setting": base, "window": w, "metric": label, "new_Sep1_start": x, "old_Oct1_start": y, "new_minus_old": (x - y) if (isnum(x) and isnum(y)) else np.nan})
            ma = mech[(mech["setting"] == new) & (mech["window"] == w)]
            mb = mech[(mech["setting"] == old) & (mech["window"] == w)]
            if len(ma) and len(mb):
                for col, label in (("made_total_Mm3", "runoff made, all types (M m3)"), ("Hsrf_Mm3", "Hsrf made (M m3)"), ("Satsrf_Mm3", "Satsrf made (M m3)")):
                    x, y = ma.iloc[0].get(col, np.nan), mb.iloc[0].get(col, np.nan)
                    if isnum(x) or isnum(y):
                        rows.append({"setting": base, "window": w, "metric": label, "new_Sep1_start": x, "old_Oct1_start": y, "new_minus_old": (x - y) if (isnum(x) and isnum(y)) else np.nan})
    return pd.DataFrame(rows)


def cmd_score(a) -> int:
    data = Path(a.data)
    out = Path(a.out_dir)
    if not data.is_dir():
        die("--data folder not found: {}".format(data))
    if not a.windows or not Path(a.windows).is_file():
        die("give --windows storm_windows.csv (written by the windows command). Run: python score_windows_sycamore.py windows --data {} --out-dir {}".format(a.data, a.out_dir))
    storms = pd.read_csv(a.windows, comment="#", parse_dates=["window_start", "window_end", "rain_start", "rain_last_end"])
    if len(storms) != 12:
        log("  WARNING: {} has {} storms (12 expected)".format(a.windows, len(storms)))
    old_map = dict(OLD_STEMS)
    if a.old_retuned:
        old_map["Retuned"] = a.old_retuned
    if a.old_original:
        old_map["Original"] = a.old_original
    only = [s.strip() for s in a.settings.split(",")] if a.settings else None
    log("\nscore_windows_sycamore.py v{}   score command".format(VERSION))
    notes = []
    gauges = load_gauges(a, data)
    by_setting, ignored = discover(data, old_map, only)
    if not by_setting:
        die("no recognized runs under {}. Expected names starting Sycamore_WY2017sep_ or Sycamore_WY2020_ (and the two old W5 runs).".format(data))
    wins = build_windows(storms)
    settings = sorted(by_setting, key=setting_sort_key)
    log("  runs recognized: {}".format(", ".join("{} [{}]".format(s, "+".join(str(y) for y in sorted(by_setting[s]))) for s in settings)))
    if ignored:
        log("  outlet files not recognized and left out: {}".format(", ".join(ignored)))

    rows, mech_rows = [], []
    run_notes = []
    for st in settings:
        runs = {}
        for year, (stem, info, paths) in sorted(by_setting[st].items()):
            runs[year] = RunYear(stem, info, paths, run_notes)
            r = runs[year]
            log("  loaded {:<46} year {} start {:%Y-%m-%d %H:%M} (from {}), {} outlet rows, .mrf {}, .rft {}".format(
                stem, year, r.start, r.start_source, len(r.raw), "yes" if r.rain is not None else "NO", "yes" if r.rft is not None else "NO"))
        by_win, storm_rows = {}, {2017: [], 2020: [], "all": []}
        for win in wins:
            row = score_window(win, runs, gauges, a.kge_basis)
            if row is None:
                continue
            row["setting"] = st
            by_win[win.wid] = row
            if win.kind == "storm":
                storm_rows["all"].append(row)
                storm_rows[win.years[0]].append(row)
            mr = mech_row(win, runs)
            if mr is not None:
                mr["setting"] = st
                mech_rows.append(mr)
        for win in wins:
            if win.kind == "all_storms" and win.wid in by_win:
                fill_all_storm_medians(by_win, win, storm_rows["all"] if win.wid == "AllStorms" else storm_rows[win.years[0]])
        rows += [by_win[w.wid] for w in wins if w.wid in by_win]
    for n in run_notes:
        log("  NOTE: " + n)
    df = pd.DataFrame(rows)
    for c in SCORE_COLS:
        if c not in df.columns:
            df[c] = np.nan
    df = df[SCORE_COLS]
    mech = pd.DataFrame(mech_rows)
    for c in MECH_COLS:
        if c not in mech.columns:
            mech[c] = np.nan
    mech = mech[MECH_COLS]
    spin = spinup_pairs(df, mech)

    # ---- the report
    L = ["score_windows_sycamore.py v{}   window scores".format(VERSION), "windows file: {}   KGE' basis: {}".format(a.windows, a.kge_basis), ""]
    L += ["DEFINITIONS"] + ["  " + x for x in DEF_LINES] + [""]
    wf = pd.read_csv(a.windows, comment="#")
    with open(a.windows, encoding="utf-8") as fh:
        prm = [ln.strip() for ln in fh if ln.startswith("# PARAMETERS")]
    L += ["STORM WINDOWS ({})".format(prm[0][2:] if prm else "parameters not recorded"),
          fmt_df(wf, ["storm_id", "year", "name", "window_start", "window_end", "window_hours", "cut_by_next_storm", "cluster_rain_mm"], nd=1), ""]
    main_cols = ["window", "gauge_class", "vol_gauge_Mm3", "vol_model_Mm3", "vol_diff_Mm3", "pbias_pct", "kge_prime", "kge_r", "kge_beta", "kge_gamma"]
    pt_cols = ["window", "peak_gauge_cfs", "peak_model_cfs", "peak_diff_cfs", "peak_ratio", "dt_peak_h", "dt10_h", "dt50_h", "dt90_h", "arr20_diff_h", "arr50_diff_h", "lag_h", "tail_vol_ratio"]
    for st in [s for s in settings if s.replace("_old", "") in ("Retuned", "Original")]:
        d = df[df["setting"] == st]
        L += ["=" * 110, "SETTING {}   (runs: {})".format(st, d["runs"].iloc[0] if len(d) else "?"), "",
              "Volume, PBIAS and KGE' by window", fmt_df(d[d["kind"] != "off"], main_cols, nd=3), "",
              "OFF windows (outside the 12 storms)", fmt_df(d[d["kind"] == "off"], ["window", "gauge_class", "vol_gauge_Mm3", "vol_model_Mm3", "vol_diff_Mm3", "pbias_pct", "days_model_gt1cfs", "days_gauge_gt1cfs",
                                                                                      "days_model_gt10cfs", "days_gauge_gt10cfs", "model_max_cfs", "gauge_max_cfs"], nd=3), "",
              "Peak flow and timing (storm windows; the AllStorms rows are medians)", fmt_df(d[d["kind"].isin(["storm", "all_storms"])], pt_cols, nd=2), ""]
        m = mech[mech["setting"] == st]
        L += ["Mechanisms: runoff the model MADE on the land (M m3; shares in % of the four types)",
              fmt_df(m, ["window", "rain_mm", "Hsrf_Mm3", "Sbsrf_Mm3", "Psrf_Mm3", "Satsrf_Mm3", "made_total_Mm3", "Hsrf_pct", "Sbsrf_pct", "Psrf_pct", "Satsrf_pct"], nd=2), ""]
    L += ["=" * 110, "ALL SETTINGS, KEY WINDOWS (PBIAS % model minus gauge | KGE')", ""]
    key = ["OctMay2017", "OctMay2020", "Both", "AllStorms", "OFF_both"]
    tab = []
    for st in settings:
        r = {"setting": st}
        for w in key:
            x = df[(df["setting"] == st) & (df["window"] == w)]
            if len(x):
                r[w + " PBIAS"] = x.iloc[0]["pbias_pct"]
                r[w + " KGE'"] = x.iloc[0]["kge_prime"]
        tab.append(r)
    L += [fmt_df(pd.DataFrame(tab), ["setting"] + [c for w in key for c in (w + " PBIAS", w + " KGE'")], nd=2), ""]
    if len(spin):
        L += ["=" * 110, "SPIN-UP CHECK: new Sep-1-start WY2017 runs vs old Oct-1-start WY2017 runs (same 2017 windows, same clock; new minus old)", ""]
        L += [fmt_df(spin, ["setting", "window", "metric", "new_Sep1_start", "old_Oct1_start", "new_minus_old"], nd=3), ""]
    else:
        L += ["SPIN-UP CHECK: not possible (the old W5 runs or the new WY2017sep Re-tuned/Original runs were not found).", ""]
    blanks = df[(df["blank_reasons"].astype(str).str.len() > 0) & (df["kind"] == "storm")]
    L += ["BLANK VALUES AND NOTES (storm windows of Re-tuned and Original; every blank has its reason in the CSV)"]
    for _, r in blanks[blanks["setting"].isin(["Retuned", "Original"])].iterrows():
        L.append("  {} {}: {}".format(r["setting"], r["window"], r["blank_reasons"]))
    L += [""]
    report = "\n".join(L) + "\n"

    hdr = ["window scores, one row per setting and window. " + DEF_LINES[0]] + DEF_LINES[1:] + ["KGE' basis: {}".format(a.kge_basis)]
    write_csv(out / "window_scores.csv", hdr, df)
    write_csv(out / "window_mechanisms.csv", ["runoff the model MADE on the land, by type, from the .rft file (M m3; all cells added; before routing); rain_mm = basin-mean rain from the .mrf file",
                                              "pct columns = share of the four types added; made_over_rain_pct = runoff made / rain volume over {:g} km2".format(BASIN_KM2)], mech)
    write_csv(out / "spinup_check.csv", ["new = Sycamore_WY2017sep_* (Time_hr 0 = Sep 1 2016); old = the Oct-1-start W5 run; same 2017 windows; new_minus_old"], spin if len(spin) else pd.DataFrame(columns=["setting", "window", "metric", "new_Sep1_start", "old_Oct1_start", "new_minus_old"]))
    write_text(out / "score_report.txt", report)
    print(report)
    log("WROTE {}".format(", ".join(str(out / f) for f in ("window_scores.csv", "window_mechanisms.csv", "spinup_check.csv", "score_report.txt"))))
    return 0


# ================================================================================================ self-test
def _pulse(t_h, t0_h, peak, rise_h, k):
    x = t_h - t0_h
    return np.where(x < 0, 0.0, np.where(x <= rise_h, peak * x / rise_h, peak * np.exp(-k * (x - rise_h))))


def _build_world(root: Path):
    """A fake bundle: gauge files, rain, outlet/.mrf/.rft/.in files for several settings with known relationships to the gauge."""
    (root / "gauge").mkdir(parents=True)
    (root / "results").mkdir()
    (root / "inputs").mkdir()
    ref0 = pd.Timestamp("2016-09-01")
    named = [s for s in STORM_TABLE if s[5]]
    flows = []          # (year, rain_start_ts, peak, rise, k)
    rain_events = []    # (year, start_ts, end_ts, mm_per_hour)
    for (y, name, f, l, mm, nm) in STORM_TABLE:
        f0, l0 = pd.Timestamp(f), pd.Timestamp(l)
        if nm:
            rs, re_ = f0 + 3 * HOUR, l0 + 15 * HOUR
            rain_events.append((y, rs, re_, 2.0, name))
            idx = [s for s in named if s[0] == y].index([s for s in named if s[1] == name and s[0] == y][0])
            peak = 80.0 + 10.0 * idx
            if (y, name) == (2017, "Dec 21-25"):
                peak = 0.0
            if (y, name) == (2020, "Nov 19-22"):
                peak = 0.085
            flows.append((y, rs + 4 * HOUR, peak, 5.0, 0.05, name))
        else:
            rain_events.append((y, f0 + 5 * HOUR, f0 + 13 * HOUR, 3.0, name))
    rain_events.append((2017, pd.Timestamp("2017-02-21 10:00"), pd.Timestamp("2017-02-21 12:00"), 0.2, "tiny"))      # 0.4 mm: must not cut a window
    rain_events.append((2017, pd.Timestamp("2016-11-15 00:00"), pd.Timestamp("2016-11-15 06:00"), 0.05, "drizzle"))  # below THR

    def rain_series(start, end):
        idx = pd.date_range(start, end, freq="60min", inclusive="left")
        r = pd.Series(0.0, index=idx)
        for (y, a, b, rate, name) in rain_events:
            seg = r.index[(r.index >= a) & (r.index < b)]
            r[seg] = rate
        # a dry break inside Dec 23-28 2019 (10 hours, shorter than GAP 12)
        brk = r.index[(r.index >= pd.Timestamp("2019-12-25 10:00")) & (r.index < pd.Timestamp("2019-12-25 20:00"))]
        r[brk] = 0.0
        return r

    def truth_cms(ts_index: pd.DatetimeIndex, year: int):
        h = ((ts_index - ref0) / HOUR).to_numpy(float) if year == 2017 else ((ts_index - pd.Timestamp("2019-09-01")) / HOUR).to_numpy(float)
        q = np.zeros(len(ts_index))
        base = ref0 if year == 2017 else pd.Timestamp("2019-09-01")
        for (y, t0, peak, rise, k, name) in flows:
            if y == year and peak > 0:
                q += _pulse(h, float((t0 - base) / HOUR), peak, rise, k)
        return q

    spec = {2017: (pd.Timestamp("2016-09-01"), pd.Timestamp("2017-06-01")), 2020: (pd.Timestamp("2019-09-01"), pd.Timestamp("2020-06-01"))}
    # gauge files
    for year, fname in ((2017, "sycamore_09510200_obs_flow_WY2017.csv"), (2020, "sycamore_09510200_obs_flow_Sep2019_to_May2020.csv")):
        lo = PERIODS[year][0]
        idx = pd.date_range(lo, PERIODS[year][1], freq="15min", inclusive="left")
        q = truth_cms(idx, year)
        q = np.where(q < 1e-4, 0.0, q)
        g = pd.DataFrame({"datetime": idx, "q_cms": q, "q_cfs": q / CFS, "tier": "A"})
        if year == 2017:
            g.loc[(idx >= "2017-01-01") & (idx < "2017-01-02"), "q_cms"] = np.nan
        g.to_csv(root / "gauge" / fname, index=False)

    def model_variant(kind, ts_index, year):
        q = truth_cms(ts_index, year)
        if kind == "perfect":
            return q
        if kind == "x1.2":
            return 1.2 * q
        if kind == "x1.25":
            return 1.25 * q
        if kind == "late2h":
            return truth_cms(ts_index - 2 * HOUR, year)
        if kind == "zero":
            return np.zeros(len(ts_index))
        if kind == "false":
            q2 = q.copy()
            if year == 2017:
                hh = ((ts_index - ref0) / HOUR).to_numpy(float)
                z0 = float((pd.Timestamp("2016-12-21 07:00") - ref0) / HOUR)
                q2 += _pulse(hh, z0, 30.0, 5.0, 0.05)                                            # flow where the gauge recorded none (zero storm)
                q2[(ts_index >= "2016-10-01") & (ts_index < "2016-11-01")] += 0.05               # false flow in October (OFF window)
            return q2
        raise ValueError(kind)

    def write_run(stem, start, end, kind, year, with_in=True):
        runtime = int((end - start) / HOUR) - 1
        th = np.arange(0.0, runtime + 1e-9, 0.0625)
        idx = start + pd.to_timedelta(np.round(th * 3600).astype("int64"), unit="s")
        pd.DataFrame({"Time_hr": th, "Q": model_variant(kind, idx, year)}).to_csv(root / "results" / "{}_Outlet.qout".format(stem), index=False)
        hr = np.arange(1.0, runtime + 1e-9, 1.0)
        hidx = start + pd.to_timedelta(hr.astype(int), unit="h") - HOUR
        rs = rain_series(start, end)
        rain = rs.reindex(hidx).fillna(0.0).to_numpy()
        mrf = pd.DataFrame({"Time_hr": hr, "RainMin_mm_hr": 0.0, "MAP_mm_hr": rain, "Srf_m3_s": 0.0, "ChannelPerc_m3": 0.0, "MDGW_mm": 1300.0, "AvSWE_cm": 0.0})
        mrf = pd.concat([mrf, pd.DataFrame({"Time_hr": [hr[-1] + 1], "RainMin_mm_hr": 9999.0, "MAP_mm_hr": 0.0, "Srf_m3_s": 0.0, "ChannelPerc_m3": 0.0, "MDGW_mm": 0.0, "AvSWE_cm": 0.0})])
        mrf.to_csv(root / "results" / "{}{}_00.mrf".format(stem, runtime), index=False)
        pd.DataFrame({"Time_hr": hr, "Hsrf_m3_s": 3.0, "Sbsrf_m3_s": 1.0, "Psrf_m3_s": 0.0, "Satsrf_m3_s": 0.0}).to_csv(root / "results" / "{}{}_00.rft".format(stem, runtime), index=False)
        if with_in:
            (root / "inputs" / "{}.in".format(stem)).write_text("STARTDATE:  # c\n{:%m/%d/%Y/%H/%M}\nRUNTIME:  # c\n{}\nOUTFILENAME:\nresults/{}\n".format(start, runtime, stem))

    plan = [("Retuned", "perfect"), ("Original", "x1.2"), ("S1_Ks1p0", "late2h"), ("S5_Ks1p0_Bed2p0", "zero"), ("S6_Ks0p5_Bed2p0", "false")]
    for tail, kind in plan:
        write_run("Sycamore_WY2017sep_" + tail, spec[2017][0], spec[2017][1], kind, 2017)
        write_run("Sycamore_WY2020_" + tail, spec[2020][0], spec[2020][1], kind, 2020)
    write_run(OLD_STEMS["Retuned"], pd.Timestamp("2016-10-01"), spec[2017][1], "perfect", 2017)
    write_run(OLD_STEMS["Original"], pd.Timestamp("2016-10-01"), spec[2017][1], "x1.25", 2017)
    return rain_events


def do_self_test() -> int:
    ok = True

    def check(cond, text):
        nonlocal ok
        print(("  PASS  " if cond else "  FAIL  ") + text)
        ok = ok and bool(cond)

    def near(a, b, tol):
        return bool(isnum(a) and abs(a - b) <= tol)

    print("SELF-TEST: a fake bundle with known answers (nothing here is real data)")
    tmp = Path(tempfile.mkdtemp(prefix="score_windows_selftest_"))
    try:
        bundle = tmp / "bundle"
        _build_world(bundle)
        out = tmp / "out"
        wargs = argparse.Namespace(data=str(bundle), out_dir=str(out), thr=0.1, gap=12, min_next=5.0, R=4, gauge_2017=None, gauge_2020=None)
        import io
        import contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = cmd_windows(wargs)
        wtxt = buf.getvalue()
        sw = pd.read_csv(out / "storm_windows.csv", comment="#", parse_dates=["window_start", "window_end", "rain_start", "rain_last_end"])
        check(len(sw) == 12, "windows: 12 storm windows found ({})".format(len(sw)))
        check("MISMATCH" in wtxt, "windows: the gauge check reports MISMATCH on the fake gauge (it must not pass by accident)")
        s1 = sw[(sw["year"] == 2017) & (sw["name"] == "Dec 30-Jan 2")].iloc[0]
        check(s1["window_start"] == pd.Timestamp("2016-12-30 03:00"), "window start = first rain hour (Dec 30 2016 03:00; got {})".format(s1["window_start"]))
        last_end = pd.Timestamp("2017-01-02 15:00") + HOUR - HOUR
        check(s1["rain_last_end"] == pd.Timestamp("2017-01-02 15:00"), "last rain end is Jan 2 15:00 (got {})".format(s1["rain_last_end"]))
        check(s1["window_end"] == pd.Timestamp("2017-01-06 15:00"), "window end = last rain end + 4 days (got {})".format(s1["window_end"]))
        j = sw[(sw["year"] == 2017) & (sw["name"] == "Jan 19-21")].iloc[0]
        check(j["window_end"] == pd.Timestamp("2017-01-23 05:00") and bool(j["cut_by_next_storm"]), "Jan 19-21 is cut at the start of the Jan 23 storm, 05:00 (got {})".format(j["window_end"]))
        f = sw[(sw["year"] == 2017) & (sw["name"] == "Feb 18-19")].iloc[0]
        check(f["window_end"] == pd.Timestamp("2017-02-23 15:00") and not bool(f["cut_by_next_storm"]), "a 0.4 mm rain cluster after Feb 18-19 does not cut the window (got {})".format(f["window_end"]))
        d = sw[(sw["year"] == 2020) & (sw["name"] == "Dec 23-28")].iloc[0]
        check(int(d["n_clusters_merged"]) == 1 and d["rain_last_end"] == pd.Timestamp("2019-12-28 15:00"), "a 10-hour dry break inside Dec 23-28 does not split the storm at GAP 12")
        sw2 = pd.read_csv(out / "windows_grid_test.csv", comment="#")
        row = sw2[(sw2["thr_mm_h"] == 0.1) & (sw2["gap_h"] == 12)].iloc[0]
        check(int(row["named_one_cluster"]) == 12, "grid test: all 12 named storms are one cluster at THR 0.1, GAP 12")
        row6 = sw2[(sw2["thr_mm_h"] == 0.1) & (sw2["gap_h"] == 6)].iloc[0]
        check(int(row6["named_split"]) >= 1, "grid test: GAP 6 splits the storm with the 10-hour break ({} split)".format(int(row6["named_split"])))

        sargs = argparse.Namespace(data=str(bundle), windows=str(out / "storm_windows.csv"), out_dir=str(out), kge_basis="hourly", settings=None,
                                   old_retuned=None, old_original=None, gauge_2017=None, gauge_2020=None)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = cmd_score(sargs)
        check(rc == 0, "score command ran")
        sc = pd.read_csv(out / "window_scores.csv", comment="#")
        me = pd.read_csv(out / "window_mechanisms.csv", comment="#")
        sp = pd.read_csv(out / "spinup_check.csv", comment="#")

        def g(setting, window):
            x = sc[(sc["setting"] == setting) & (sc["window"] == window)]
            return x.iloc[0] if len(x) else None

        r = g("Retuned", "OctMay2017")
        check(near(r["pbias_pct"], 0.0, 0.5) and near(r["kge_prime"], 1.0, 0.01), "perfect run: OctMay2017 PBIAS {:.3f} % and KGE' {:.4f}".format(r["pbias_pct"], r["kge_prime"]))
        r = g("Retuned", "Both")
        a17, a20 = g("Retuned", "OctMay2017"), g("Retuned", "OctMay2020")
        check(near(r["vol_gauge_Mm3"], a17["vol_gauge_Mm3"] + a20["vol_gauge_Mm3"], 2e-3) and near(r["vol_model_Mm3"], a17["vol_model_Mm3"] + a20["vol_model_Mm3"], 2e-3),
              "Both = the two years added (gauge {:.2f} M m3)".format(r["vol_gauge_Mm3"]))
        o = g("Original", "Both")
        check(near(o["pbias_pct"], 20.0, 0.5) and near(o["kge_beta"], 1.2, 0.01), "x1.2 run: PBIAS is +20 ({:.2f}) and beta 1.2 ({:.3f}): positive = model too high".format(o["pbias_pct"], o["kge_beta"]))
        check(o["kge_prime"] < 1.0 and o["kge_prime"] > 0.5, "x1.2 run: KGE' below 1 ({:.3f})".format(o["kge_prime"]))
        z = sc[(sc["setting"] == "S5_Ks1p0_Bed2p0".split("_")[0]) & (sc["window"] == "OctMay2017")]
        check(len(z) == 1 and near(z.iloc[0]["pbias_pct"], -100.0, 1e-6), "zero-flow model: PBIAS -100")
        check("model flow is zero" in str(z.iloc[0]["blank_reasons"]), "zero-flow model: KGE' blank with the reason")
        # storms
        st = sc[(sc["setting"] == "Retuned") & (sc["kind"] == "storm")]
        zs = st[st["name"] == "Dec 21-25"].iloc[0]
        check(zs["gauge_class"] == "zero" and not isnum(zs["pbias_pct"]) and "gauge class is zero" in zs["blank_reasons"], "zero-gauge storm (Dec 21-25 2016): PBIAS and timing blank with the reason")
        ns = st[st["name"] == "Nov 19-22"].iloc[0]
        check(ns["gauge_class"] == "near-zero" and not isnum(ns["kge_prime"]), "near-zero storm (Nov 19-22 2019): ratios blank (class {})".format(ns["gauge_class"]))
        fs = sc[(sc["setting"] == "S6") & (sc["name"] == "Dec 21-25")].iloc[0]
        check(isnum(fs["peak_diff_cfs"]) and fs["peak_diff_cfs"] > 100 and fs["vol_diff_Mm3"] > 0.5 and not isnum(fs["peak_ratio"]),
              "false-flow run on the zero storm: plain differences shown (peak +{:.0f} cfs, volume +{:.2f} M m3), ratio blank".format(fs["peak_diff_cfs"], fs["vol_diff_Mm3"]))
        big = st[st["name"] == "Jan 19-21"].iloc[0]
        check(near(big["pbias_pct"], 0.0, 1.0) and near(big["peak_ratio"], 1.0, 0.02) and near(big["dt_peak_h"], 0.0, 0.1) and near(big["dt50_h"], 0.0, 0.1) and near(big["lag_h"], 0.0, 0.1),
              "perfect run, Jan 19-21: PBIAS {:.2f}, peak ratio {:.3f}, dt_peak {:.2f}, dt50 {:.2f}, lag {}".format(big["pbias_pct"], big["peak_ratio"], big["dt_peak_h"], big["dt50_h"], big["lag_h"]))
        late = sc[(sc["setting"] == "S1") & (sc["name"] == "Feb 22-23")].iloc[0]
        check(near(late["dt_peak_h"], 2.0, 0.1) and near(late["dt10_h"], 2.0, 0.15) and near(late["dt50_h"], 2.0, 0.15) and near(late["dt90_h"], 2.0, 0.3) and near(late["lag_h"], 2.0, 0.1),
              "2-hour-late run, Feb 22-23: dt_peak {:.2f}, dt10 {:.2f}, dt50 {:.2f}, dt90 {:.2f}, lag {} (all about +2: model later)".format(late["dt_peak_h"], late["dt10_h"], late["dt50_h"], late["dt90_h"], late["lag_h"]))
        check(near(late["arr20_diff_h"], 2.0, 0.2) and near(late["arr50_diff_h"], 2.0, 0.2), "late run: arrival through 20% and 50% of the gauge peak about +2 h ({:.2f}, {:.2f})".format(late["arr20_diff_h"], late["arr50_diff_h"]))
        sc_big = sc[(sc["setting"] == "Original") & (sc["name"] == "Jan 19-21")].iloc[0]
        check(near(sc_big["peak_ratio"], 1.2, 0.02) and near(sc_big["peak_diff_cfs"], 0.2 * sc_big["peak_gauge_cfs"], 0.5), "x1.2 run: peak ratio 1.2 ({:.3f})".format(sc_big["peak_ratio"]))
        # AllStorms
        al = g("S1", "AllStorms")
        sts = sc[(sc["setting"] == "S1") & (sc["kind"] == "storm")]
        okst = sts[sts["gauge_class"] == "ok"]
        check(int(al["n_ok_storms"]) == 10 and near(al["dt_peak_h"], float(np.median(okst["dt_peak_h"])), 1e-9), "AllStorms: median of the 10 storms with a meaningful gauge hydrograph (dt_peak {:.2f})".format(al["dt_peak_h"]))
        allv = g("Retuned", "AllStorms")
        check(near(allv["vol_gauge_Mm3"], float(st["vol_gauge_Mm3"].sum()), 2e-3), "AllStorms volume = the 12 storm windows added ({:.2f} M m3)".format(allv["vol_gauge_Mm3"]))
        # OFF complement
        for y in (2017, 2020):
            parts = sc[(sc["setting"] == "Retuned") & (sc["kind"] == "storm") & (sc["years"] == str(y))]
            offr = g("Retuned", "OFF_{}".format(y))
            season = g("Retuned", "OctMay{}".format(y))
            check(int(parts["n_hours"].sum() + offr["n_hours"]) == int(season["n_hours"]), "OFF_{} + the {} storm windows = all Oct-May hours ({})".format(y, len(parts), int(season["n_hours"])))
        of = g("S6", "OFF_2017")
        check(of["days_model_gt1cfs"] > 25 and of["days_gauge_gt1cfs"] < of["days_model_gt1cfs"], "false-flow run: OFF_2017 days above 1 cfs: model {} vs gauge {}".format(of["days_model_gt1cfs"], of["days_gauge_gt1cfs"]))
        # Jan 1 outage
        jj = st[st["name"] == "Dec 30-Jan 2"].iloc[0]
        check(int(jj["gauge_hours"]) == int(jj["n_hours"]) - 24 and near(jj["pbias_pct"], 0.0, 1.0), "Jan 1 outage: gauge has {} of {} hours; matched-hour PBIAS {:.2f} still about 0".format(int(jj["gauge_hours"]), int(jj["n_hours"]), jj["pbias_pct"]))
        check(jj["vol_model_allhours_Mm3"] > jj["vol_model_Mm3"], "Jan 1 outage: the all-hours model volume is larger than the matched-hours volume")
        # mechanisms
        m = me[(me["setting"] == "Retuned") & (me["window"] == "OctMay2017")].iloc[0]
        check(near(m["Hsrf_pct"], 75.0, 1e-6) and near(m["Sbsrf_pct"], 25.0, 1e-6) and near(m["made_total_Mm3"], 4 * 3600 * m_hours(sc) / 1e6, 1e-3), "mechanisms: Hsrf 75 %, Sbsrf 25 %, total {:.3f} M m3 = 4 m3/s x hours".format(m["made_total_Mm3"]))
        mr = me[(me["setting"] == "Retuned") & (me["window"] == "storm01_2017_Dec21-25")].iloc[0]
        expect_rain = 2.0 * float((pd.Timestamp("2016-12-25 15:00") - pd.Timestamp("2016-12-21 03:00")) / HOUR)
        check(near(mr["rain_mm"], expect_rain, 1e-6), "rain depth in the first storm window = {:.1f} mm (got {:.1f}); the hour-end labels were moved back".format(expect_rain, mr["rain_mm"]))
        # spin-up
        sr = sp[(sp["setting"] == "Retuned") & (sp["metric"] == "model volume, every hour (M m3)")]
        check(len(sr) > 5 and sr["new_minus_old"].abs().max() < 1e-6, "spin-up: Retuned new vs old identical in all {} windows ({} max difference)".format(len(sr), "0" if not len(sr) else "{:.2g}".format(sr["new_minus_old"].abs().max())))
        so = sp[(sp["setting"] == "Original") & (sp["metric"] == "model volume, every hour (M m3)") & (sp["window"] == "OctMay2017")].iloc[0]
        check(so["new_minus_old"] < 0 and near(so["new_Sep1_start"] / so["old_Oct1_start"], 1.2 / 1.25, 1e-3), "spin-up: Original new (x1.2) / old (x1.25) = {:.4f}, new minus old is negative".format(so["new_Sep1_start"] / so["old_Oct1_start"]))
        check("PBIAS (%)" in set(sp["metric"]) and "old_Oct1_start" in sp.columns, "spin-up table has PBIAS and the old-start column")
        # the old run's clock: Time_hr 0 = Oct 1, so its Oct 1 flow must line up with the new run's
        # 15-minute basis
        sargs.kge_basis = "15min"
        sargs.out_dir = str(tmp / "out15")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            cmd_score(sargs)
        s15 = pd.read_csv(tmp / "out15" / "window_scores.csv", comment="#")
        k15 = s15[(s15["setting"] == "Original") & (s15["window"] == "Both")].iloc[0]
        check(near(k15["kge_beta"], 1.2, 0.01), "15-minute KGE' basis runs and gives beta 1.2 ({:.3f})".format(k15["kge_beta"]))
        # never overwrite a file that is not ours
        (tmp / "out2").mkdir()
        (tmp / "out2" / "storm_windows.csv").write_text("my own file\n")
        try:
            wargs.out_dir = str(tmp / "out2")
            with contextlib.redirect_stdout(io.StringIO()):
                cmd_windows(wargs)
            check(False, "refuses to overwrite a file it did not write")
        except SystemExit:
            check((tmp / "out2" / "storm_windows.csv").read_text() == "my own file\n", "refuses to overwrite a file it did not write (the file is unchanged)")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\nSELF-TEST " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


def m_hours(sc):
    r = sc[(sc["setting"] == "Retuned") & (sc["window"] == "OctMay2017")].iloc[0]
    return float(r["n_hours"])


# ================================================================================================ command line
def build_parser():
    p = argparse.ArgumentParser(description="Window scorer for the Lower Sycamore runs (see the text at the top of this file).",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--data", required=False, help="folder with the unzipped bundle (or the lab folder); searched below it")
    common.add_argument("--out-dir", default="score_windows_out")
    common.add_argument("--gauge-2017", default=None, help="WY2017 gauge file (default: found under --data)")
    common.add_argument("--gauge-2020", default=None, help="WY2020 gauge file (default: found under --data)")
    w = sub.add_parser("windows", parents=[common], help="find the storm windows from the rain and test THR, GAP and R")
    w.add_argument("--thr", type=float, default=0.1, help="basin-mean rain (mm/h) at or above which an hour counts as rain (default 0.1)")
    w.add_argument("--gap", type=int, default=12, help="dry hours that separate two storms (default 12)")
    w.add_argument("--min-next", type=float, default=5.0, help="smallest rain cluster (mm) that counts as a next storm (default 5)")
    w.add_argument("--R", type=int, default=4, help="days after the last rain that the window runs (default 4; test 3, 4, 5 first)")
    s = sub.add_parser("score", parents=[common], help="score every run in every window")
    s.add_argument("--windows", help="storm_windows.csv from the windows command")
    s.add_argument("--kge-basis", choices=["hourly", "15min"], default="hourly", help="KGE' on hourly means (default) or on the gauge's 15-minute times")
    s.add_argument("--settings", help="only these settings, comma separated, for example Retuned,Original")
    s.add_argument("--old-retuned", help="outlet run name of the old Oct-1-start Re-tuned run (default {})".format(OLD_STEMS["Retuned"]))
    s.add_argument("--old-original", help="outlet run name of the old Oct-1-start Original run (default {})".format(OLD_STEMS["Original"]))
    sub.add_parser("self-test", help="check every number on a fake bundle with known answers")
    return p


def main(argv=None) -> int:
    warnings.filterwarnings("ignore", category=FutureWarning)
    a = build_parser().parse_args(argv)
    if a.cmd == "self-test":
        return do_self_test()
    if not a.data:
        die("give --data (the folder with the unzipped bundle).")
    return cmd_windows(a) if a.cmd == "windows" else cmd_score(a)


if __name__ == "__main__":
    sys.exit(main())
