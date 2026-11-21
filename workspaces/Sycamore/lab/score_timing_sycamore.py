#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
score_timing_sycamore.py   (v1.3, 2026-11-17)
=============================================

Timing and shape scores for the Lower Sycamore Wash storms: WY2017 (the calibration year, 5 storms, 4 whole-year runs)
and WY2020 (the validation year, 7 storms, 12 runs).

NEW IN v1.3: the pass 2 CSV has the `role` of every RUN (not every year) and a second column, `fixed_before_wy2020_results`.
  Nothing else changed: every number, flag, blank and reason is identical to v1.2. See RUN_DESIGN below for the labels and why.

READ-ONLY. It reads outlet files, the two gauge files and storm_table.csv. It runs no model, changes none of your
files, and writes only into its own output folder (timing_out/ by default). Every definition it uses is written at the
top of each CSV it makes (lines that start with "#"); `--definitions` prints them without running anything.

WHAT IT DOES, IN ORDER
  1. VOLUME CHECK (always first, in both passes). It recomputes the storm volumes you already have on record
     (the WY2020 write-up, the volumes in storm_table.csv for WY2017, and the running-log volumes for S1, S9, S10, S11
     and the seven-storm totals for S5 to S8) and compares them with what the outlet files and gauge files give now.
     If any value differs by more than 0.005 M m3 it STOPS and writes nothing.
  2. PASS 1  (--pass 1)   GAUGE ONLY. One row per storm, both years: window, gauge volume, gauge peak and its time,
     the gauge times by which 10, 50 and 90 percent of the window volume had passed, flow at both window edges,
     gauge quality, outage and gap flags, and the rain columns. It stops there so you can choose the "too small to
     time" floor.
  3. PASS 2  (--pass 2 --floor-mm3 X)   every run against the gauge, per storm (storms whose GAUGE volume is under X get no timing): volume-time differences, peak ratio and
     peak-time difference with flags, hourly Pearson r, cross-correlation lag within +/-12 h, tail volume ratio.
     Writes ONE long-format CSV. Then it prints CHECK 5 (the Original run on the WY2017 Feb 27-28 flood).

COMMANDS (from workspaces/Sycamore/lab/; one at a time, foreground)
  python score_timing_sycamore.py --show-arrival              # step 0: prints the first-arrival rule in YOUR score_run_sycamore.py
  python score_timing_sycamore.py --pass 1 --dry-run          # lists the files and windows it found; computes nothing
  python score_timing_sycamore.py --pass 1                    # volume check, then the gauge-only table
  python score_timing_sycamore.py --pass 2 --floor-mm3 0.05   # after you have chosen the floor (0.05 is only an example)
  python score_timing_sycamore.py --definitions               # prints every definition

OPTIONS YOU MAY NEED
  --windows as-scored|cut|uncut   as-scored (default): each year's existing windows (WY2020 uncut, WY2017 cut at the next storm).
                                  cut / uncut: one rule for every storm (storms with no reference volume then say NOT CHECKED).
  --stem 2017:original=NAME       use another outlet-file name (without _Outlet.qout) for one run. Repeat the option as needed.
                                  Years are 2017 or 2020; run keys are retuned, original, kchan05, kchan02, S1, S5 ... S11.
  --obs2017 FILE, --obs2020 FILE  other gauge files.   --results DIR   other outlet folder.   --out-dir DIR   other output folder.
  --print-all                     in pass 2, print every run (default prints Re-tuned and Original only).
  --ignore-check-failure          DIAGNOSIS ONLY: carry on after a volume-check mismatch; every row then says MISMATCH.

WHAT IT READS (never changed)
  results/<run>_Outlet.qout   columns Time_hr, Q (m3/s). Time_hr 0 = 2019-09-01 00:00 (WY2020) or 2016-10-01 00:00 (WY2017).
  ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv                  datetime, q_cms [, tier]
  ../init_data/obs/sycamore_09510200_obs_flow_Sep2019_to_May2020.csv      (or ..._WY2020.csv if that one is missing)
  storm_table.csv             only for the rain columns (matched on water year and first rain day)
  score_run_sycamore.py       only with --show-arrival

WORDS
  hourly mean      the average of all readings in one clock hour. Hours with no reading are left out, never filled, never zero.
  M m3             million cubic meters.      cfs   cubic feet per second (m3/s divided by 0.0283168).
  volume time t50  hours from the window start by which 50 percent of the window's volume had passed.
  window           00:00 on the first rain day to 00:00 four days after the last rain day (see --windows for the cut rule).
  Pearson r        the usual correlation between two hourly series (1 = same shape, 0 = unrelated).
  lag              the shift in hours that makes the model hourly series best match the gauge; positive = model later.
  floor            gauge volume below which a storm is "too small to time"; you choose it after pass 1. It is applied to the GAUGE only.
                   --model-floor-mm3 Y (optional) also blanks the timing of any run whose own volume is under Y.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = "1.3"
SELF_MARK = "score_timing_sycamore.py"
HERE = Path(__file__).resolve().parent
CFS_PER_CMS = 1.0 / 0.028316846592
HOUR = pd.Timedelta(hours=1)
TOL_M3 = 0.005                 # volume-check tolerance, M m3
LOW_VOL_RATIO = 0.5             # console note only: model volume under this fraction of the gauge's
EDGE_H = 6.0                   # peak within this many hours of a window edge is flagged
BUMP_FRAC = 0.8                # a second bump at least this fraction of the highest hourly flow is flagged
DIP_FRAC = 0.8                 # ... and only if the flow between the two bumps falls to this fraction of the lower bump or less
FLAT_SPAN_H = 1.0              # a peak whose readings within 1 percent of it span more than this is called flat-topped
MAX_LAG_H = 12                 # cross-correlation search range, hours
MIN_COMMON_H = 12              # fewest hours two hourly series need in common for a lag to be scored
MAX_STEP_H = 3.0               # a rise through a level is not trusted across a reading gap longer than this (as in the scorer)
OLD_T50_PEAK_H = -0.58         # the old time-to-50%-of-gauge-peak figure for CHECK 5, hours (your record)
OUTAGE = (pd.Timestamp("2017-01-01 00:00"), pd.Timestamp("2017-01-03 00:00"))   # the Jan 1-2 2017 gauge outage

YEARS = {
    "WY2017": dict(role="calibration", start=pd.Timestamp("2016-10-01 00:00"),
                   obs=["sycamore_09510200_obs_flow_WY2017.csv"]),
    "WY2020": dict(role="validation", start=pd.Timestamp("2019-09-01 00:00"),
                   obs=["sycamore_09510200_obs_flow_Sep2019_to_May2020.csv", "sycamore_09510200_obs_flow_WY2020.csv"]),
}

# (key, label, outlet-file stem without _Outlet.qout)
RUNS = {
    "WY2020": [
        ("retuned", "Re-tuned", "Sycamore_WY2020_Retuned"),
        ("original", "Original", "Sycamore_WY2020_Original"),
        ("kchan05", "Original + Kchan 0.5", "Sycamore_WY2020_Original_kchan0p5"),
        ("kchan02", "Original + Kchan 0.2", "Sycamore_WY2020_Original_kchan0p2"),
        ("S1", "S1", "Sycamore_WY2020_S1_Ks1p0"),
        ("S5", "S5", "Sycamore_WY2020_S5_Ks1p0_Bed2p0"),
        ("S6", "S6", "Sycamore_WY2020_S6_Ks0p5_Bed2p0"),
        ("S7", "S7", "Sycamore_WY2020_S7_Ks1p0_Bed2p0_As500"),
        ("S8", "S8", "Sycamore_WY2020_S8_Ks1p0_Bed2p0_f0p02"),
        ("S9", "S9", "Sycamore_WY2020_S9_Ks1p0_Bed0p5"),
        ("S10", "S10", "Sycamore_WY2020_S10_Ks1p0_GW0p5"),
        ("S11", "S11", "Sycamore_WY2020_S11_Ks1p0_PsiB1p0"),
    ],
    "WY2017": [
        ("retuned", "Re-tuned", "W5_d012_k10_f0p15_as300_ksx1p03_nocl"),
        ("original", "Original", "W5_d012_k10_base_nocl"),
        ("kchan05", "Original + Kchan 0.5", "W5_d012_k10_base_kchan0p5"),
        ("kchan02", "Original + Kchan 0.2", "W5_d012_k10_base_kchan0p2"),
    ],
}
PRIORITY = ("retuned", "original")

# Design role of each RUN (pass 2 only) and whether its settings were written down BEFORE any WY2020 gauge volume or model
# result was opened.  (role, fixed_before_wy2020_results).  Sources: Tail_Screen_W4_Running_Log.md (Re-tuned frozen Nov 3 20:15,
# Original named in the same entry; Kchan 0.5 and 0.2 written down Nov 4 02:10, after the gauge volumes and the Re-tuned and
# Original results had been opened; their values came from the WY2017 probes) and the screen scripts (S1, S5 to S11 built from Nov 7 on).
RUN_DESIGN = {
    ("WY2017", "retuned"): ("calibration", "n/a"),
    ("WY2017", "original"): ("calibration", "n/a"),
    ("WY2017", "kchan05"): ("calibration", "n/a"),
    ("WY2017", "kchan02"): ("calibration", "n/a"),
    ("WY2020", "retuned"): ("validation", "yes"),
    ("WY2020", "original"): ("validation", "yes"),
    ("WY2020", "kchan05"): ("validation", "no"),
    ("WY2020", "kchan02"): ("validation", "no"),
    **{("WY2020", k): ("post_validation_screen", "no") for k in ("S1", "S5", "S6", "S7", "S8", "S9", "S10", "S11")},
}
assert all((y, k) in RUN_DESIGN for y, runs_ in RUNS.items() for k, _, _ in runs_), "RUN_DESIGN must cover every run in RUNS"


def run_design(year: str, key: str):
    """(role, fixed_before_wy2020_results) for one run; a key that is not in RUNS (only possible in a unit test) is 'unassigned'."""
    return RUN_DESIGN.get((year, key), ("unassigned", "n/a"))

# (year, name, first rain day, last rain day, first day of the next storm in storm_table.csv or None)
STORMS = [
    ("WY2017", "Dec 21-25", "2016-12-21", "2016-12-25", "2016-12-30"),
    ("WY2017", "Dec 30-Jan 2", "2016-12-30", "2017-01-02", "2017-01-13"),
    ("WY2017", "Jan 19-21", "2017-01-19", "2017-01-21", "2017-01-23"),
    ("WY2017", "Feb 18-19", "2017-02-18", "2017-02-19", "2017-02-27"),
    ("WY2017", "Feb 27-28", "2017-02-27", "2017-02-28", "2017-03-23"),
    ("WY2020", "Nov 19-22", "2019-11-19", "2019-11-22", "2019-11-28"),
    ("WY2020", "Nov 28-30", "2019-11-28", "2019-11-30", "2019-12-03"),
    ("WY2020", "Dec 23-28", "2019-12-23", "2019-12-28", "2020-01-21"),
    ("WY2020", "Feb 10-11", "2020-02-10", "2020-02-11", "2020-02-22"),
    ("WY2020", "Feb 22-23", "2020-02-22", "2020-02-23", "2020-03-10"),
    ("WY2020", "Mar 10-13", "2020-03-10", "2020-03-13", "2020-03-18"),
    ("WY2020", "Mar 18-19", "2020-03-18", "2020-03-19", None),
]

RAIN_COLS = ["basin_mm", "peak_1h_mm", "peak_3h_mm", "peak_6h_mm", "top6h_share", "t50_frac", "dry_days_before", "prev7_mm", "prev30_mm"]
# the same nine columns under the readable headers of the lab storm_table.csv (short name first, then the readable one)
RAIN_ALIASES = {
    "basin_mm": ["basin_mm", "Basin rain (mm)"],
    "peak_1h_mm": ["peak_1h_mm", "Peak 1 h (mm)"],
    "peak_3h_mm": ["peak_3h_mm", "Peak 3 h (mm)"],
    "peak_6h_mm": ["peak_6h_mm", "Peak 6 h (mm)"],
    "top6h_share": ["top6h_share", "Top-6 h share (fraction)"],
    "t50_frac": ["t50_frac", "Half-rain point (fraction)"],
    "dry_days_before": ["dry_days_before", "Dry days before (d)"],
    "prev7_mm": ["prev7_mm", "Prior 7 days (mm)"],
    "prev30_mm": ["prev30_mm", "Prior 30 days (mm)"],
}

# ------------------------------------------------------------------------------------------------ reference volumes (M m3)
# Source labels: W = WY2020 write-up (as coded in add_wy2020_results_to_storm_table.py); T = storm_table.csv (WY2017 columns);
# L = running log (Soil_Storage_Running_Log / part 2).  Order of the WY2020 tuples: gauge, retuned, original, kchan05, kchan02.
WY2020_EXPECTED = {
    "Nov 19-22": (0.015, 2.800, 3.165, 3.019, 3.091), "Nov 28-30": (1.473, 3.235, 3.538, 3.341, 3.433),
    "Dec 23-28": (4.069, 0.218, 0.154, 0.044, 0.083), "Feb 10-11": (0.386, 1.767, 1.765, 1.641, 1.700),
    "Feb 22-23": (4.131, 4.839, 4.981, 4.825, 4.899), "Mar 10-13": (5.626, 0.875, 0.493, 0.340, 0.413),
    "Mar 18-19": (3.317, 1.031, 0.841, 0.683, 0.775),
}
WY2020_ORDER = ["Nov 19-22", "Nov 28-30", "Dec 23-28", "Feb 10-11", "Feb 22-23", "Mar 10-13", "Mar 18-19"]
WY2020_LOG = {   # per-storm volumes in WY2020_ORDER order
    "S1": (0.024, 0.164, 2.487, 0.634, 1.222, 1.966, 2.055),
    "S9": (0.005, 0.602, 3.777, 0.396, 1.085, 1.815, 2.484),
    "S10": (0.399, 1.236, 4.258, 0.796, 1.581, 2.239, 2.196),
    "S11": (0.0105, 0.4414, 3.1876, 0.2578, 0.7194, 1.0855, 1.9640),
}
WY2020_TOTALS = {"S5": 6.29, "S6": 4.11, "S7": 7.97, "S8": 9.06}   # seven-storm totals (log); two decimals, so tolerance 0.006
WY2017_REF = {   # storm_table.csv, event (cut) windows: gauge_hours, gauge_Mm3, model_Mm3, asdrawn_nocl, asdrawn_k05, asdrawn_k02
    "Dec 21-25": ("192/192", 0.0, 1.1794, 1.3767, 1.2078, 1.2861),
    "Dec 30-Jan 2": ("144/168", 0.8158, 0.7723, 1.01, 0.8435, 0.9256),
    "Jan 19-21": ("96/96", 5.5667, 7.3129, 7.8355, 7.7118, 7.761),
    "Feb 18-19": ("120/120", 1.2193, 2.4445, 2.1989, 2.0556, 2.1292),
    "Feb 27-28": ("120/120", 3.8067, 3.5866, 3.2561, 3.0908, 3.1617),
}
WY2017_JAN_UNCUT = (7.1497, 7.6066)          # gauge_storm_Mm3, model_storm_Mm3 (Re-tuned) for Jan 19-21 to Jan 25 00:00
WY2017_MATCHED_RATIO = ("Dec 30-Jan 2", 0.58, 0.01)   # Re-tuned model over the gauge's 144 hours / gauge volume, from the table's note


def build_refs():
    """(year, storm) -> list of (window_end, {series_key: (value, source)}, gauge_hours_text or None)."""
    refs = {}
    for name, tup in WY2020_EXPECTED.items():
        first, last = [(f, l) for y, n, f, l, _ in STORMS if y == "WY2020" and n == name][0]
        hi = pd.Timestamp(last) + pd.Timedelta(days=4)
        d = {k: (v, "write-up") for k, v in zip(["gauge", "retuned", "original", "kchan05", "kchan02"], tup)}
        for sk, vals in WY2020_LOG.items():
            d[sk] = (vals[WY2020_ORDER.index(name)], "running log")
        refs[("WY2020", name)] = [(hi, d, None)]
    for name, (gh, g, m, o, k5, k2) in WY2017_REF.items():
        first, last, nxt = [(f, l, n) for y, nm, f, l, n in STORMS if y == "WY2017" and nm == name][0]
        hi_unc = pd.Timestamp(last) + pd.Timedelta(days=4)
        hi_cut = min(hi_unc, pd.Timestamp(nxt)) if nxt else hi_unc
        d = {"gauge": (g, "storm_table.csv"), "retuned": (m, "storm_table.csv"), "original": (o, "storm_table.csv"),
             "kchan05": (k5, "storm_table.csv"), "kchan02": (k2, "storm_table.csv")}
        refs[("WY2017", name)] = [(hi_cut, d, gh)]
        if name == "Jan 19-21":
            refs[("WY2017", name)].append((hi_unc, {"gauge": (WY2017_JAN_UNCUT[0], "storm_table.csv"),
                                                    "retuned": (WY2017_JAN_UNCUT[1], "storm_table.csv")}, None))
    return refs


REFS = build_refs()


# ------------------------------------------------------------------------------------------------ small helpers
def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    sys.exit("\nSTOPPED: " + msg + "\n")


def isnum(x) -> bool:
    return x is not None and isinstance(x, (int, float, np.floating, np.integer)) and math.isfinite(float(x))


def fnum(x, nd: int = 3, dash: str = "-") -> str:
    return dash if not isnum(x) else f"{float(x):.{nd}f}"


def cfs(q_cms) -> float:
    return float(q_cms) * CFS_PER_CMS


class Win:
    """One scored window."""

    def __init__(self, year, name, first, last, nxt, mode):
        self.year, self.name = year, name
        self.first, self.last = pd.Timestamp(first), pd.Timestamp(last)
        self.next_first = pd.Timestamp(nxt) if nxt else None
        self.lo = self.first
        self.hi_uncut = self.last + pd.Timedelta(days=4)
        self.hi_cut = min(self.hi_uncut, self.next_first) if self.next_first is not None else self.hi_uncut
        if mode == "uncut":
            self.hi = self.hi_uncut
        elif mode == "cut":
            self.hi = self.hi_cut
        else:
            self.hi = self.hi_uncut if year == "WY2020" else self.hi_cut
        self.hours = int(round((self.hi - self.lo) / HOUR))
        self.rule = ("cut at the next storm's first day" if self.hi < self.hi_uncut else "last rain day + 4 d")
        self.role = YEARS[year]["role"]

    @property
    def key(self):
        return (self.year, self.name)


def make_windows(mode: str):
    return [Win(y, n, f, l, nx, mode) for y, n, f, l, nx in STORMS]


# ------------------------------------------------------------------------------------------------ readers
def read_qout(path: Path, start: pd.Timestamp) -> pd.Series:
    with open(path) as fh:
        fh.readline()
        second = fh.readline()
    if "," in second:
        raw = pd.read_csv(path, sep=",", header=0, float_precision="round_trip")
    else:
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None)
    raw = raw.iloc[:, :2].copy()
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    s = pd.Series(raw["Q"].to_numpy(float), index=start + pd.to_timedelta(raw["Time_hr"].to_numpy(float), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_gauge(path: Path) -> pd.DataFrame:
    o = pd.read_csv(path, float_precision="round_trip")
    for c in ("datetime", "q_cms"):
        if c not in o.columns:
            die(f"{path} has no column '{c}'. Columns found: {list(o.columns)[:10]}")
    idx = pd.to_datetime(o["datetime"])
    df = pd.DataFrame({"q_cms": pd.to_numeric(o["q_cms"], errors="coerce").to_numpy(float)}, index=idx)
    if "tier" in o.columns:
        df["tier"] = o["tier"].astype(str).str.strip().str.lower().to_numpy()
    df = df[~df.index.duplicated(keep="last")].sort_index()
    return df


def read_rain_table(path: Path):
    """Returns (dict (year, first_day) -> {col: text}, notes list). Soft: problems become notes, not stops."""
    notes = []
    if not path.is_file():
        return {}, [f"storm table not found: {path}. The rain columns are blank."]
    raw = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    need = {"year": ["year", "Water year"], "first_day": ["first_day", "First rain day"], "last_day": ["last_day", "Last rain day"]}
    hdr = {}
    for k, opts in need.items():
        hit = [c for c in opts if c in raw.columns]
        if not hit:
            return {}, [f"storm table {path.name} has no '{opts[0]}' column; the rain columns are blank. Columns found: {list(raw.columns)[:12]}"]
        hdr[k] = hit[0]
    out = {}
    src = {}                                   # short name -> the header actually found in this table
    for c in RAIN_COLS:
        hit = [h for h in RAIN_ALIASES[c] if h in raw.columns]
        if hit:
            src[c] = hit[0]
    missing = [c for c in RAIN_COLS if c not in src]
    if missing:
        notes.append(f"storm table {path.name} has no column for {missing} (looked for the short headers and the readable ones); "
                     f"those rain columns are blank. Columns found: {list(raw.columns)[:14]} ...")
    for _, r in raw.iterrows():
        out[(("WY" + re.sub(r"\D", "", r[hdr["year"]])[-4:]), str(r[hdr["first_day"]]).strip())] = {
            **{c: (r[src[c]].strip() if c in src else "") for c in RAIN_COLS}, "_last": str(r[hdr["last_day"]]).strip()}
    return out, notes


# ------------------------------------------------------------------------------------------------ series numbers
def hourly(s: pd.Series, lo, hi) -> pd.Series:
    x = s[(s.index >= lo) & (s.index < hi)]
    return x.resample("60min").mean().dropna()


def vol_mm3(h: pd.Series) -> float:
    return float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")


def volume_times(h: pd.Series, lo, fracs=(0.1, 0.5, 0.9)):
    """Hours from lo by which each fraction of the window volume had passed. Returns (total M m3, {frac: hours} or None, reason)."""
    if len(h) == 0:
        return float("nan"), None, "no hourly data in the window"
    v = h.to_numpy(float) * 3600.0
    t0 = ((h.index - lo) / HOUR).to_numpy(float)
    c = np.cumsum(v)
    tot = float(c[-1])
    if not tot > 0:
        return tot / 1e6, None, "window volume is zero"
    out = {}
    for f in fracs:
        target = f * tot
        k = int(np.searchsorted(c, target, side="left"))
        prev = float(c[k - 1]) if k > 0 else 0.0
        out[f] = float(t0[k] + (target - prev) / v[k])
    return tot / 1e6, out, ""


def volume_after(h: pd.Series, lo, hi, t_clock) -> float:
    """Volume (M m3) of the hourly series after the clock time t_clock (each hour spread evenly over itself)."""
    if len(h) == 0:
        return float("nan")
    ends = h.index + HOUR
    frac = np.clip(((ends - t_clock) / HOUR).to_numpy(float), 0.0, 1.0)     # share of each hour that lies after t_clock
    return float(np.sum(h.to_numpy(float) * 3600.0 * frac) / 1e6)


def second_bump(h: pd.Series) -> bool:
    v = h.to_numpy(float)
    if len(v) < 3:
        return False
    i = int(np.argmax(v))
    p = v[i]
    if not p > 0:
        return False
    n = len(v)
    for j in range(n):
        if j == i or v[j] < BUMP_FRAC * p:
            continue
        if (j > 0 and v[j] < v[j - 1]) or (j < n - 1 and v[j] < v[j + 1]):
            continue                                              # not a local maximum
        a, b = (i, j) if i < j else (j, i)
        if float(v[a:b + 1].min()) <= DIP_FRAC * v[j]:
            return True
    return False


def peak_info(raw: pd.Series, h: pd.Series, lo, hi):
    """Highest reading in the window (raw), its time, the highest hourly mean, and flags. Returns (dict or None, reason)."""
    x = raw[(raw.index >= lo) & (raw.index < hi)].dropna()
    if len(x) == 0:
        return None, "no readings in the window"
    v = x.to_numpy(float)
    i = int(np.argmax(v))
    pk = float(v[i])
    if not pk > 0:
        return None, "highest flow in the window is zero"
    thr = 0.99 * pk
    a = i
    while a > 0 and v[a - 1] >= thr:
        a -= 1
    b = i
    while b < len(v) - 1 and v[b + 1] >= thr:
        b += 1
    flat_h = float((x.index[b] - x.index[a]) / HOUR)
    t = x.index[i]
    flags = []
    if min((t - lo) / HOUR, (hi - t) / HOUR) < EDGE_H:
        flags.append(f"peak within {EDGE_H:g} h of a window edge")
    if second_bump(h):
        flags.append(f"second bump >= {int(BUMP_FRAC * 100)}% as high")
    if flat_h > FLAT_SPAN_H:
        flags.append(f"flat top ({flat_h:.1f} h within 1% of the peak)")
    return {"cms": pk, "cfs": cfs(pk), "time": t, "hourly_cms": float(h.max()) if len(h) else float("nan"),
            "flags": "; ".join(flags)}, ""


def pearson(g: pd.Series, m: pd.Series):
    """Pearson r on the hours both series have. Returns (r or None, n, reason)."""
    idx = g.index.intersection(m.index)
    if len(idx) < 3:
        return None, len(idx), "fewer than 3 hours in common"
    a, b = g.loc[idx].to_numpy(float), m.loc[idx].to_numpy(float)
    if a.std() == 0:
        return None, len(idx), "gauge hourly flow is constant in the window"
    if b.std() == 0:
        return None, len(idx), "model hourly flow is constant in the window"
    return float(np.corrcoef(a, b)[0, 1]), len(idx), ""


def best_lag(g: pd.Series, m: pd.Series):
    """Lag in hours (within +/-MAX_LAG_H) with the highest correlation. Positive = model later than gauge.
    Returns (lag, r, at_limit, reason)."""
    best = None
    for L in range(-MAX_LAG_H, MAX_LAG_H + 1):
        m2 = m.copy()
        m2.index = m2.index - pd.Timedelta(hours=L)               # model value at t+L set against the gauge value at t
        idx = g.index.intersection(m2.index)
        if len(idx) < MIN_COMMON_H:
            continue
        a, b = g.loc[idx].to_numpy(float), m2.loc[idx].to_numpy(float)
        if a.std() == 0 or b.std() == 0:
            continue
        r = float(np.corrcoef(a, b)[0, 1])
        if best is None or r > best[1] + 1e-12 or (abs(r - best[1]) <= 1e-12 and abs(L) < abs(best[0])):
            best = (L, r)
    if best is None:
        return None, None, False, "no lag had enough hours in common with varying flow in both series"
    return best[0], best[1], abs(best[0]) == MAX_LAG_H, ""


def arrival(raw: pd.Series, level: float, lo, hi):
    """Hours from lo at which the flow first RISES THROUGH level (linear between readings). Returns (hours or None, reason)."""
    x = raw[(raw.index >= lo) & (raw.index < hi)].dropna()
    if len(x) == 0:
        return None, "no readings in the window"
    v = x.to_numpy(float)
    t = ((x.index - lo) / HOUR).to_numpy(float)
    if v[0] >= level:
        return None, "already at or above the level at the window start"
    idx = np.nonzero(v >= level)[0]
    if len(idx) == 0:
        return None, "never reaches the level"
    j = int(idx[0])
    if t[j] - t[j - 1] > MAX_STEP_H:
        return None, "the first rise through the level falls inside a reading gap"
    return float(t[j - 1] + (level - v[j - 1]) / (v[j] - v[j - 1]) * (t[j] - t[j - 1])), ""


def missing_runs(h: pd.Series, lo, hi):
    grid = pd.date_range(lo, hi - HOUR, freq="60min")
    miss = grid.difference(h.index)
    runs = []
    if len(miss) == 0:
        return runs
    s = p = miss[0]
    for t in miss[1:]:
        if t - p == HOUR:
            p = t
            continue
        runs.append((s, p + HOUR))
        s = p = t
    runs.append((s, p + HOUR))
    return runs


def fmt_runs(runs, limit=4):
    parts = [f"{a:%b %d %H:%M}-{b:%b %d %H:%M} ({int(round((b - a) / HOUR))} h)" for a, b in runs[:limit]]
    if len(runs) > limit:
        parts.append(f"+{len(runs) - limit} more")
    return "; ".join(parts)


# ------------------------------------------------------------------------------------------------ gauge bundle (pass 1)
def gauge_bundle(w: Win, gdf: pd.DataFrame) -> dict:
    q = gdf["q_cms"]
    h = hourly(q, w.lo, w.hi)
    why = []
    tot, times, wt = volume_times(h, w.lo)
    if times is None:
        why.append(f"gauge volume times: {wt}")
    pk, wp = peak_info(q, h, w.lo, w.hi)
    if pk is None:
        why.append(f"gauge peak: {wp}")
    runs = missing_runs(h, w.lo, w.hi)
    left, right = h.get(w.lo), h.get(w.hi - HOUR)
    if left is None:
        why.append("gauge left edge flow: no gauge data in the first hour")
    if right is None:
        why.append("gauge right edge flow: no gauge data in the last hour")
    quality = {"pct_interp_or_est": None, "pct_filled_zero": None, "pct_provisional": None}
    if "tier" in gdf.columns:
        sub = gdf[(gdf.index >= w.lo) & (gdf.index < w.hi)]
        t = sub.loc[sub["q_cms"].notna(), "tier"]
        n = len(t)
        if n:
            quality["pct_interp_or_est"] = 100.0 * float(t.isin(["estimated", "filled_interp"]).sum()) / n
            quality["pct_filled_zero"] = 100.0 * float((t == "filled_zero").sum()) / n
            quality["pct_provisional"] = 100.0 * float((t == "provisional").sum()) / n
        else:
            why.append("gauge quality percentages: no readings in the window")
    else:
        why.append("gauge quality percentages: the gauge file has no 'tier' column")
    gap_before = ""
    if times is not None and runs:
        t90 = w.lo + pd.Timedelta(hours=times[0.9])
        gap_before = "yes" if any(a < t90 for a, _ in runs) else "no"
    elif times is not None:
        gap_before = "no"
    touches = (w.lo < OUTAGE[1]) and (w.hi > OUTAGE[0])
    return {"h": h, "raw": q, "tot": tot, "times": times, "times_why": wt, "peak": pk, "peak_why": wp, "runs": runs,
            "left": left, "right": right, "quality": quality, "gap_before_t90": gap_before, "touches_outage": touches, "why": why}


def rain_cells(rain, w: Win):
    r = rain.get((w.year, w.first.strftime("%Y-%m-%d")))
    return {f"rain_{c}": (r[c] if r else "") for c in RAIN_COLS}


# ------------------------------------------------------------------------------------------------ volume check
def volume_check(wins, gauges, models, tol=TOL_M3):
    """Returns (status dict {(year, storm, series_key): 'MATCH'|'MISMATCH'|'NOT CHECKED'}, mismatches list, printable lines)."""
    status, bad, lines = {}, [], []
    here = {}
    for w in wins:
        gh = hourly(gauges[w.year]["q_cms"], w.lo, w.hi)
        here[(w.year, w.name, "gauge")] = vol_mm3(gh)
        here[(w.year, w.name, "gauge_hours")] = f"{len(gh)}/{w.hours}"
        for key, q in models[w.year].items():
            here[(w.year, w.name, key)] = vol_mm3(hourly(q, w.lo, w.hi))
    labels = {k: l for yr in RUNS for k, l, _ in RUNS[yr]}
    labels["gauge"] = "gauge"
    for yr in ("WY2017", "WY2020"):
        lines.append(f"  {yr}")
        lines.append("    {:<13}{:>7} {:>16} {:>16}   {}".format("storm", "gauge", "Re-tuned (ratio)", "Original (ratio)", "check (S5-S8: by their totals, below)"))
        for w in [x for x in wins if x.year == yr]:
            refs = [r for r in REFS.get(w.key, []) if r[0] == w.hi]
            ref = refs[0] if refs else None
            n_ok = n_bad = n_nc = 0
            bad_here = []
            keys = ["gauge"] + list(models[yr].keys())
            for key in keys:
                got = here[(yr, w.name, key)]
                if ref is not None and key in ref[1]:
                    exp, src = ref[1][key]
                    ok = isnum(got) and abs(got - exp) <= tol
                    status[(yr, w.name, key)] = "MATCH" if ok else "MISMATCH"
                    if ok:
                        n_ok += 1
                    else:
                        n_bad += 1
                        bad_here.append(f"{yr} {w.name} {labels.get(key, key)}: here {fnum(got, 4)} vs {fnum(exp, 4)} ({src})")
                else:
                    status[(yr, w.name, key)] = "NOT CHECKED"
                    n_nc += 1
            if ref is not None and ref[2] is not None:                       # gauge hours (WY2017 table)
                got_h = here[(yr, w.name, "gauge_hours")]
                if got_h == ref[2]:
                    n_ok += 1
                else:
                    n_bad += 1
                    bad_here.append(f"{yr} {w.name} gauge hours with data: here {got_h} vs {ref[2]} (storm_table.csv)")
            bad.extend(bad_here)
            g = here[(yr, w.name, "gauge")]

            def vr(k):
                if k not in models[yr]:
                    return "      not found"
                v = here[(yr, w.name, k)]
                return f"{v:7.3f} ({v / g:5.2f})" if isnum(g) and g > 0 else f"{v:7.3f} (   - )"

            lines.append("    {:<13}{:>7.3f} {:>16} {:>16}   {} ok, {} not checked{}".format(
                w.name, g, vr("retuned"), vr("original"), n_ok, n_nc, f", {n_bad} MISMATCH" if n_bad else ""))
    # like-for-like ratio for the outage storm (Re-tuned), from the table's note
    nm, exp, tl = WY2017_MATCHED_RATIO
    w = [x for x in wins if x.year == "WY2017" and x.name == nm][0]
    if "retuned" in models["WY2017"] and any(r[0] == w.hi for r in REFS[w.key]):
        gh = hourly(gauges["WY2017"]["q_cms"], w.lo, w.hi)
        mh = hourly(models["WY2017"]["retuned"], w.lo, w.hi).reindex(gh.index).dropna()
        gv = vol_mm3(gh)
        ratio = vol_mm3(mh) / gv if isnum(gv) and gv > 0 else float("nan")
        ok = isnum(ratio) and abs(ratio - exp) <= tl
        lines.append(f"    like-for-like ratio, Re-tuned, {nm} (model over the {len(gh)} gauge hours / gauge): here {fnum(ratio, 3)} vs {exp} (storm_table.csv note)  {'ok' if ok else 'MISMATCH'}")
        if not ok:
            bad.append(f"WY2017 {nm} like-for-like Re-tuned/gauge ratio: here {fnum(ratio, 3)} vs {exp}")
    # seven-storm totals for S5..S8 (only when the WY2020 windows are the scorer's)
    wy20 = [x for x in wins if x.year == "WY2020"]
    if all(x.hi == x.hi_uncut for x in wy20):
        for key, exp in WY2020_TOTALS.items():
            if key in models["WY2020"]:
                tot = sum(here[("WY2020", n, key)] for n in WY2020_ORDER)
                ok = isnum(tot) and abs(tot - exp) <= 0.006
                lines.append(f"    WY2020 {key} seven-storm total: here {fnum(tot, 3)} vs {exp} (running log)  {'ok' if ok else 'MISMATCH'}")
                for n in WY2020_ORDER:
                    status[("WY2020", n, key)] = "MATCH" if ok else "MISMATCH"
                if not ok:
                    bad.append(f"WY2020 {key} seven-storm total: here {fnum(tot, 3)} vs {exp} (running log)")
    else:
        lines.append("    WY2020 S5-S8 seven-storm totals not checked (the WY2020 windows are not the scorer's).")
    return status, bad, lines


# ------------------------------------------------------------------------------------------------ definitions
def definition_lines(mode: str, floor, model_floor=None):
    m = {"as-scored": "WY2020 uncut (exactly as score_wy2020_run_sycamore.py: Nov 28-30 runs to Dec 4); WY2017 cut at the next storm's first day when that comes sooner (only Jan 19-21 is affected: it ends Jan 23, not Jan 25)",
         "cut": "every storm cut at the next storm's first day when that comes sooner (WY2020 Nov 28-30 ends Dec 3; WY2017 Jan 19-21 ends Jan 23)",
         "uncut": "no storm is cut (WY2020 Nov 28-30 ends Dec 4; WY2017 Jan 19-21 ends Jan 25, which includes the Jan 23-24 storm)"}[mode]
    L = [
        f"{SELF_MARK} v{VERSION} (2026-11-17). READ-ONLY: reads outlet files, gauge files and storm_table.csv; changes none of them.",
        "ROLES: in the pass 2 CSV the role belongs to the RUN: calibration = the four WY2017 whole-year runs; validation = the WY2020 frozen setting (Re-tuned) and its three comparison runs (Original, Original + Kchan 0.5, Original + Kchan 0.2); post_validation_screen = the eight WY2020 one-change screens S1 and S5 to S11, built after the WY2020 gauge and the Re-tuned and Original results had been opened, so they are NOT out of sample. In the pass 1 CSV (gauge only, one row per storm) the role is the storm's year: WY2017 calibration, WY2020 validation. Times are Arizona local time, as in the gauge and model files.",
        "fixed_before_wy2020_results (pass 2 only): yes = the run's settings were written down before any WY2020 gauge volume or model result was opened (Re-tuned, frozen Nov 3 20:15; Original, named as a comparison in the same entry); no = written down or built after results had been opened (Kchan 0.5 and 0.2, Nov 4 02:10, with values taken from the WY2017 probes; the eight screens, Nov 7 onward); n/a = WY2017 runs (calibration year). Source: Tail_Screen_W4_Running_Log.md and the screen scripts; this script does not read them, the labels are fixed in RUN_DESIGN.",
        "HOURLY MEAN: average of all readings in one clock hour [hh:00, hh+1:00), m3/s. Hours with no reading are left out (never filled, never zero). Same as score_wy2020_run_sycamore.py.",
        "VOLUME: sum of hourly means x 3600 s over the hours that have data, M m3 (million cubic meters). A gauge window with missing hours gives a minimum volume. The model volume counts all its hours.",
        "WINDOW: starts 00:00 on the first rain day; ends 00:00 four days after the last rain day. Rule used in this run (--windows " + mode + "): " + m + ".",
        "VOLUME TIME tNN (NN = 10, 50, 90): hours from the window start by which NN percent of the window volume had passed. Cumulative sum of the hourly volumes, each hour's volume spread evenly over its own hour, linear inside the hour where the target is crossed. Gauge hours with no data add no volume.",
        "TIME DIFFERENCES are always model minus gauge, in hours: positive = the model is later than the gauge.",
        "RESOLUTION: the gauge is read every 15 min and the model written every 3.75 min, so an hourly mean of gauge readings is centered about 5.6 min (0.09 h) earlier in true time than the model's. In the volume-time columns the model therefore looks about 0.09 h earlier than it is (checked on synthetic data). Differences smaller than about 0.1 h are inside this convention (the one the existing scorers use). The peak, lag-free arrival and peak-time columns use the raw readings and do not have this offset, but are limited by the 15-min gauge spacing.",
        "PEAK: the highest single reading in the window (raw: gauge every 15 min, model every output step), shown in cfs (m3/s / 0.0283168); the first one if several are equal. PEAK (HOURLY) = highest hourly mean, shown only so it can be compared with earlier hourly-based figures. Peak ratio = model peak / gauge peak.",
        f"PEAK FLAGS: (1) the highest reading is within {EDGE_H:g} h of a window edge; (2) a second bump at least {int(BUMP_FRAC * 100)}% as high as the highest hourly mean, with the flow between the two bumps falling to {int(DIP_FRAC * 100)}% of the lower bump or less (checked on hourly means; hours with no data are skipped); (3) flat top: readings within 1% of the peak span more than {FLAT_SPAN_H:g} h. A flagged peak time is still shown, not blanked.",
        "PEARSON r: ordinary correlation of the gauge and model hourly means over the hours both have, in the window. n_hours_r = how many hours.",
        f"LAG: shift in whole hours, within +/-{MAX_LAG_H} h, at which the model hourly series correlates best with the gauge hourly series (positive = model later; ties go to the smaller shift). Only shifts with at least {MIN_COMMON_H} common hours count. xcorr_at_limit = the best lag is at +/-{MAX_LAG_H}.",
        "TAIL: the tail starts at the GAUGE's t90 (so both series use the same clock interval) and runs to the window end. Tail volume ratio = model volume in the tail / gauge volume in the tail. model_tail_fraction = model volume in the tail / the model's own window volume (the gauge's is 0.10 by construction). Hours with no gauge data add no tail volume.",
        "BLANKS: an undefined value is left blank and the reason is in blank_reasons; never zero.",
        ("FLOOR: " + (f"{floor:g} M m3 (your choice), applied to the GAUGE window volume only. If the gauge window volume is below it, the timing columns (t10, t50, t90 differences, peak-time difference, lag, tail, arrival extras) are blank with that reason; volumes, peak ratio and Pearson r are still shown. A run whose own volume is small still gets timing numbers: read them next to vol_ratio_model_over_gauge (the printed table adds a note when the model has under " + f"{LOW_VOL_RATIO:.0%} of the gauge volume; that note is display only)." + ((" MODEL FLOOR " + f"{model_floor:g} M m3 (your choice, --model-floor-mm3): a run whose own window volume is below it also gets blank timing columns; below_floor says which side.") if model_floor is not None else " No model floor was set.") if floor is not None else "not applied in pass 1 (you choose it after seeing this table).")),
        "EXTRA COLUMNS ON A FRACTION OF THE GAUGE PEAK (labeled that way, not used in anything else): extra_arr20_of_gauge_peak_diff_h and extra_arr50_of_gauge_peak_diff_h = model minus gauge time of the first rise through 20% / 50% of the GAUGE's highest reading, linear between readings, not trusted across a reading gap longer than 3 h.",
        "EDGE FLOWS: hourly mean of the first hour and of the last hour of the window, cfs.",
        "GAUGE QUALITY: from the 'tier' column of the gauge file. pct_interp_or_est = readings tagged estimated or filled_interp as a percent of the window's readings that have a value; pct_filled_zero and pct_provisional likewise; pct_missing_hours = hours with no gauge data / hours in the window. Blank if the file has no 'tier' column.",
        "OUTAGE: jan1_2_outage_flag is set when the window overlaps 2017-01-01 00:00 to 2017-01-03 00:00. gauge_missing_intervals lists every run of hours with no gauge data in the window (this also catches the WY2020 dropouts). gap_before_t90 = yes when any missing gauge hour starts before the gauge's t90, so the gauge volume times are then biased.",
        "RAIN COLUMNS (rain_*): copied unchanged from storm_table.csv for the storm with the same water year and first rain day: basin rain (mm), highest 1, 3 and 6 h rain (mm), share of the rain in the six wettest hours, timing of half the rain (t50_frac, as defined there), dry days before, rain in the previous 7 and 30 days (mm). Not recomputed here.",
        "VOLUME CHECK: before anything else the volumes already on record (write-up, storm_table.csv, running log) are recomputed; volume_check = MATCH / MISMATCH / NOT CHECKED per run and storm (tolerance 0.005 M m3).",
    ]
    return L


def write_csv(path: Path, header_lines, df: pd.DataFrame):
    if path.exists():
        with open(path, encoding="utf-8") as fh:
            first = fh.readline()
        if SELF_MARK not in first:
            die(f"{path} already exists and was not written by this script. Move it or use --out-dir.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        for ln in header_lines:
            fh.write("# " + ln + "\n")
        df.to_csv(fh, index=False)


def r4(x):
    return round(float(x), 4) if isnum(x) else None


# ------------------------------------------------------------------------------------------------ pass 1
P1_COLS = ["water_year", "role", "storm", "window_rule", "window_start", "window_end", "window_hours",
           "gauge_vol_Mm3", "gauge_volume_check", "gauge_hours_with_data", "gauge_peak_cfs", "gauge_peak_time",
           "gauge_peak_hourly_cfs", "gauge_t10_h", "gauge_t50_h", "gauge_t90_h", "gauge_left_edge_cfs", "gauge_right_edge_cfs",
           "pct_interp_or_est", "pct_filled_zero", "pct_provisional", "pct_missing_hours", "jan1_2_outage_flag",
           "gauge_missing_intervals", "gap_before_t90", "gauge_peak_flag"] + [f"rain_{c}" for c in RAIN_COLS] + ["blank_reasons"]


def pass1_row(w: Win, G: dict, status, rain):
    q = G["quality"]
    why = list(G["why"])
    t = G["times"]
    row = {
        "water_year": w.year, "role": w.role, "storm": w.name, "window_rule": w.rule,
        "window_start": f"{w.lo:%Y-%m-%d %H:%M}", "window_end": f"{w.hi:%Y-%m-%d %H:%M}", "window_hours": w.hours,
        "gauge_vol_Mm3": r4(G["tot"]), "gauge_volume_check": status.get((w.year, w.name, "gauge"), "NOT CHECKED"),
        "gauge_hours_with_data": len(G["h"]),
        "gauge_peak_cfs": r4(G["peak"]["cfs"]) if G["peak"] else None,
        "gauge_peak_time": f"{G['peak']['time']:%Y-%m-%d %H:%M}" if G["peak"] else None,
        "gauge_peak_hourly_cfs": r4(cfs(G["peak"]["hourly_cms"])) if G["peak"] else None,
        "gauge_t10_h": r4(t[0.1]) if t else None, "gauge_t50_h": r4(t[0.5]) if t else None, "gauge_t90_h": r4(t[0.9]) if t else None,
        "gauge_left_edge_cfs": r4(cfs(G["left"])) if G["left"] is not None else None,
        "gauge_right_edge_cfs": r4(cfs(G["right"])) if G["right"] is not None else None,
        "pct_interp_or_est": r4(q["pct_interp_or_est"]), "pct_filled_zero": r4(q["pct_filled_zero"]),
        "pct_provisional": r4(q["pct_provisional"]),
        "pct_missing_hours": r4(100.0 * (w.hours - len(G["h"])) / w.hours),
        "jan1_2_outage_flag": "touches the Jan 1-2 2017 outage" if G["touches_outage"] else "",
        "gauge_missing_intervals": fmt_runs(G["runs"]), "gap_before_t90": G["gap_before_t90"],
        "gauge_peak_flag": G["peak"]["flags"] if G["peak"] else "",
        "blank_reasons": " | ".join(why),
    }
    row.update(rain_cells(rain, w))
    return row


def print_pass1(rows):
    log("\nPASS 1 TABLE (gauge only; volumes M m3, times hours from the window start, flows cfs)")
    hdr = "{:<7}{:<13}{:<26}{:>7}{:>8}{:<12}{:>8}{:>8}{:>8}{:>9}{:>9}{:>8}{:>8}  {}".format(
        "year", "storm", "window", "vol", "peak", "  peak at", "t10", "t50", "t90", "left", "right", "%i/e", "%miss", "flags")
    log(hdr)
    for r in rows:
        fl = []
        if r["jan1_2_outage_flag"]:
            fl.append("JAN1-2 OUTAGE")
        if r["gauge_missing_intervals"]:
            fl.append("gap " + ("before t90" if r["gap_before_t90"] == "yes" else "after t90"))
        if r["gauge_peak_flag"]:
            fl.append("peak: " + r["gauge_peak_flag"].split(";")[0])
        pk_at = (r["gauge_peak_time"] or "")[5:16] if r["gauge_peak_time"] else "-"
        log("{:<7}{:<13}{:<26}{:>7}{:>8}  {:<10}{:>8}{:>8}{:>8}{:>9}{:>9}{:>8}{:>8}  {}".format(
            r["water_year"][2:], r["storm"], f"{r['window_start'][5:10]} to {r['window_end'][5:10]} ({r['window_hours']} h)",
            fnum(r["gauge_vol_Mm3"], 3), fnum(r["gauge_peak_cfs"], 0), pk_at, fnum(r["gauge_t10_h"], 1), fnum(r["gauge_t50_h"], 1),
            fnum(r["gauge_t90_h"], 1), fnum(r["gauge_left_edge_cfs"], 1), fnum(r["gauge_right_edge_cfs"], 1),
            fnum(r["pct_interp_or_est"], 1), fnum(r["pct_missing_hours"], 1), "; ".join(fl)))
    log("  (%i/e = percent of readings interpolated or estimated; %miss = percent of window hours with no gauge data)")


# ------------------------------------------------------------------------------------------------ pass 2
P2_COLS = ["water_year", "role", "fixed_before_wy2020_results", "run_key", "run_label", "run_stem", "storm", "window_rule", "window_start", "window_end",
           "volume_check", "gauge_vol_Mm3", "model_vol_Mm3", "vol_ratio_model_over_gauge", "below_floor",
           "gauge_t10_h", "model_t10_h", "t10_model_minus_gauge_h",
           "gauge_t50_h", "model_t50_h", "t50_model_minus_gauge_h",
           "gauge_t90_h", "model_t90_h", "t90_model_minus_gauge_h",
           "gauge_peak_cfs", "model_peak_cfs", "peak_ratio_model_over_gauge", "peak_hourly_ratio_model_over_gauge",
           "gauge_peak_time", "model_peak_time", "peak_time_model_minus_gauge_h", "gauge_peak_flag", "model_peak_flag",
           "pearson_r_hourly", "n_hours_r", "xcorr_best_lag_h", "xcorr_best_r", "xcorr_at_limit",
           "tail_gauge_Mm3", "tail_model_Mm3", "tail_volume_ratio", "model_tail_fraction",
           "extra_arr20_of_gauge_peak_diff_h", "extra_arr50_of_gauge_peak_diff_h",
           "jan1_2_outage_flag", "gauge_missing_intervals", "gap_before_t90"] + [f"rain_{c}" for c in RAIN_COLS] + ["blank_reasons"]


def pass2_row(w: Win, G: dict, key, label, stem, mq: pd.Series, floor: float, status, rain, model_floor=None):
    R = []
    row = {c: None for c in P2_COLS}
    role, fixed = run_design(w.year, key)
    row.update({"water_year": w.year, "role": role, "fixed_before_wy2020_results": fixed, "run_key": key, "run_label": label, "run_stem": stem, "storm": w.name,
                "window_rule": w.rule, "window_start": f"{w.lo:%Y-%m-%d %H:%M}", "window_end": f"{w.hi:%Y-%m-%d %H:%M}",
                "volume_check": status.get((w.year, w.name, key), "NOT CHECKED"),
                "jan1_2_outage_flag": "touches the Jan 1-2 2017 outage" if G["touches_outage"] else "",
                "gauge_missing_intervals": fmt_runs(G["runs"]), "gap_before_t90": G["gap_before_t90"]})
    row.update(rain_cells(rain, w))
    gh = G["h"]
    mh = hourly(mq, w.lo, w.hi)
    gtot, mtot = G["tot"], vol_mm3(mh)
    row["gauge_vol_Mm3"], row["model_vol_Mm3"] = r4(gtot), r4(mtot)
    if len(mh) < w.hours:
        R.append(f"model: the outlet file has hourly means for only {len(mh)} of the window's {w.hours} h; every model-based value is blank")
        row["blank_reasons"] = " | ".join(R)
        return row
    if isnum(gtot) and gtot > 0:
        row["vol_ratio_model_over_gauge"] = r4(mtot / gtot)
    else:
        R.append("vol_ratio: gauge volume is zero or missing")
    gt, gpk = G["times"], G["peak"]
    mtot_chk, mt, mwhy = volume_times(mh, w.lo)
    mpk, mpwhy = peak_info(mq, mh, w.lo, w.hi)
    if gt:
        row["gauge_t10_h"], row["gauge_t50_h"], row["gauge_t90_h"] = r4(gt[0.1]), r4(gt[0.5]), r4(gt[0.9])
    if mt:
        row["model_t10_h"], row["model_t50_h"], row["model_t90_h"] = r4(mt[0.1]), r4(mt[0.5]), r4(mt[0.9])
    if gpk:
        row["gauge_peak_cfs"], row["gauge_peak_time"], row["gauge_peak_flag"] = r4(gpk["cfs"]), f"{gpk['time']:%Y-%m-%d %H:%M}", gpk["flags"]
    if mpk:
        row["model_peak_cfs"], row["model_peak_time"], row["model_peak_flag"] = r4(mpk["cfs"]), f"{mpk['time']:%Y-%m-%d %H:%M}", mpk["flags"]
    # peak ratio (not gated by the floor)
    if gpk and mpk:
        row["peak_ratio_model_over_gauge"] = r4(mpk["cms"] / gpk["cms"])
        if gpk["hourly_cms"] > 0:
            row["peak_hourly_ratio_model_over_gauge"] = r4(mpk["hourly_cms"] / gpk["hourly_cms"])
    else:
        R.append("peak ratio: " + (("gauge " + G["peak_why"]) if not gpk else ("model " + mpwhy)))
    # Pearson r (not gated by the floor)
    r, n, rwhy = pearson(gh, mh)
    row["n_hours_r"] = n
    if r is None:
        R.append("pearson_r_hourly: " + rwhy)
    else:
        row["pearson_r_hourly"] = r4(r)
    # floor gate for the timing block
    gate = None
    below = []
    if floor is not None and (not isnum(gtot) or gtot < floor):
        below.append("gauge")
    if model_floor is not None and (not isnum(mtot) or mtot < model_floor):
        below.append("model")
    if below:
        row["below_floor"] = " and ".join(below)
        gate = "; ".join(f"{s} volume {fnum(gtot if s == 'gauge' else mtot, 4)} M m3 is below the {'floor' if s == 'gauge' else 'model floor'} {(floor if s == 'gauge' else model_floor):g}" for s in below)
    timing_keys = ["t10_model_minus_gauge_h", "t50_model_minus_gauge_h", "t90_model_minus_gauge_h", "peak_time_model_minus_gauge_h",
                   "xcorr_best_lag_h", "tail_volume_ratio", "extra_arr20_of_gauge_peak_diff_h"]
    if gate:
        R.append("timing columns (t10/t50/t90 difference, peak time, lag, tail, arrival extras): " + gate)
        row["blank_reasons"] = " | ".join(R)
        return row
    # volume-time differences
    if gt is None:
        R.append("t10/t50/t90 differences: gauge " + G["times_why"])
    elif mt is None:
        R.append("t10/t50/t90 differences: model " + mwhy)
    else:
        for f, c in ((0.1, "t10"), (0.5, "t50"), (0.9, "t90")):
            row[f"{c}_model_minus_gauge_h"] = r4(mt[f] - gt[f])
    # peak time difference
    if gpk and mpk:
        row["peak_time_model_minus_gauge_h"] = r4((mpk["time"] - gpk["time"]) / HOUR)
    else:
        R.append("peak_time_model_minus_gauge_h: gauge or model peak undefined")
    # lag
    lag, lr, atl, lwhy = best_lag(gh, mh)
    if lag is None:
        R.append("xcorr_best_lag_h: " + lwhy)
    else:
        row["xcorr_best_lag_h"], row["xcorr_best_r"], row["xcorr_at_limit"] = lag, r4(lr), ("yes" if atl else "no")
    # tail
    if gt is None:
        R.append("tail: gauge t90 undefined")
    else:
        t0 = w.lo + pd.Timedelta(hours=gt[0.9])
        tg, tm = volume_after(gh, w.lo, w.hi, t0), volume_after(mh, w.lo, w.hi, t0)
        row["tail_gauge_Mm3"], row["tail_model_Mm3"] = r4(tg), r4(tm)
        if isnum(tg) and tg > 0:
            row["tail_volume_ratio"] = r4(tm / tg)
        else:
            R.append("tail_volume_ratio: gauge tail volume is zero")
        if isnum(mtot) and mtot > 0:
            row["model_tail_fraction"] = r4(tm / mtot)
    # extras on a fraction of the gauge peak
    if gpk:
        for frac, c in ((0.2, "extra_arr20_of_gauge_peak_diff_h"), (0.5, "extra_arr50_of_gauge_peak_diff_h")):
            tg_, wg = arrival(G["raw"], frac * gpk["cms"], w.lo, w.hi)
            tm_, wm = arrival(mq, frac * gpk["cms"], w.lo, w.hi)
            if tg_ is None or tm_ is None:
                R.append(f"{c}: " + (("gauge " + wg) if tg_ is None else ("model " + wm)))
            else:
                row[c] = r4(tm_ - tg_)
    else:
        R.append("extra arrival columns: gauge peak undefined")
    row["blank_reasons"] = " | ".join(R)
    return row


def print_pass2(rows, print_all):
    log("\nPASS 2 SUMMARY (model minus gauge; hours; blanks shown as -)")
    log("{:<7}{:<13}{:<22}{:>7}{:>7}{:>7}{:>10}{:>9}{:>7}{:>6}{:>7}   {}".format(
        "year", "storm", "run", "dt10", "dt50", "dt90", "pk ratio", "dt peak", "r", "lag", "tail", "flags / floor"))
    for r in rows:
        if not print_all and r["run_key"] not in PRIORITY:
            continue
        fl = []
        if r["below_floor"]:
            fl.append("BELOW FLOOR: " + r["below_floor"])
        vr = r["vol_ratio_model_over_gauge"]
        if isnum(vr) and vr < LOW_VOL_RATIO and not r["below_floor"]:
            fl.append(f"model volume only {vr:.0%} of gauge")
        if r["model_peak_flag"]:
            fl.append("model " + r["model_peak_flag"].split(";")[0])
        if r["gauge_peak_flag"]:
            fl.append("gauge " + r["gauge_peak_flag"].split(";")[0])
        if r["jan1_2_outage_flag"]:
            fl.append("JAN1-2 OUTAGE")
        if r["gap_before_t90"] == "yes":
            fl.append("gauge gap before t90")
        log("{:<7}{:<13}{:<22}{:>7}{:>7}{:>7}{:>10}{:>9}{:>7}{:>6}{:>7}   {}".format(
            r["water_year"][2:], r["storm"], r["run_label"], fnum(r["t10_model_minus_gauge_h"], 1), fnum(r["t50_model_minus_gauge_h"], 1),
            fnum(r["t90_model_minus_gauge_h"], 1), fnum(r["peak_ratio_model_over_gauge"], 2), fnum(r["peak_time_model_minus_gauge_h"], 1),
            fnum(r["pearson_r_hourly"], 2), fnum(r["xcorr_best_lag_h"], 0), fnum(r["tail_volume_ratio"], 2), "; ".join(fl)))


def check5(rows):
    lines = ["CHECK 5: Original (loss off), WY2017 Feb 27-28 flood: 50% volume time next to the old time-to-50%-of-gauge-peak figure"]
    sel = [r for r in rows if r["water_year"] == "WY2017" and r["storm"] == "Feb 27-28" and r["run_key"] in ("original", "retuned")]
    if not [r for r in sel if r["run_key"] == "original"]:
        lines.append("  The WY2017 Original outlet file was not found, so CHECK 5 cannot be reported (see the missing-file list above).")
        return lines
    for r in sorted(sel, key=lambda x: x["run_key"] != "original"):
        lines.append(f"  {r['run_label']:<10} ({r['run_stem']})")
        lines.append(f"     volume t50: gauge {fnum(r['gauge_t50_h'], 2)} h, model {fnum(r['model_t50_h'], 2)} h, model minus gauge {fnum(r['t50_model_minus_gauge_h'], 2)} h")
        lines.append(f"     extra (fraction of the GAUGE peak): first rise through 50%: {fnum(r['extra_arr50_of_gauge_peak_diff_h'], 2)} h;  through 20%: {fnum(r['extra_arr20_of_gauge_peak_diff_h'], 2)} h")
    lines.append(f"  old figure on record, time to 50% of the gauge peak: {OLD_T50_PEAK_H:+.2f} h")
    lines.append("  (My 50%/20% arrival columns use linear interpolation between readings and this script's window; the old scorer's exact interpolation was not compared.)")
    return lines


# ------------------------------------------------------------------------------------------------ --show-arrival
def show_arrival(path: Path) -> int:
    if not path.is_file():
        die(f"cannot find {path}. Run this from workspaces/Sycamore/lab/, or give the path: --show-arrival PATH")
    text = path.read_text(encoding="utf-8", errors="replace")
    md5 = hashlib.md5(text.encode("utf-8")).hexdigest()
    log(f"\nFile: {path}\n  md5 {md5}   {len(text.splitlines())} lines\n")
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        log(f"  (could not parse the file: {exc})")
        tree = None
    found = False
    if tree is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and node.name == "arrival":
                found = True
                log(f"The function `arrival` (lines {node.lineno} to {node.end_lineno}):\n")
                log(ast.get_source_segment(text, node))
    if not found:
        log("No function named `arrival` was found.")
    log("\nEvery line that mentions `level` or calls arrival( (with line numbers):")
    for i, ln in enumerate(text.splitlines(), 1):
        if re.search(r"\blevel\b|arrival\(|0\.[25]\s*\*", ln):
            log(f"  {i:>5}: {ln.rstrip()}")
    log("\nWhat to look for: a line like `level = 0.2 * pg` means 20 percent of the GAUGE peak (pg); `0.5 * pg` would mean 50 percent.")
    return 0


# ------------------------------------------------------------------------------------------------ main
def build_parser():
    ap = argparse.ArgumentParser(description="Timing and shape scores for the Lower Sycamore storms (read-only; see the top of this file).")
    ap.add_argument("--pass", dest="pass_", type=int, choices=[1, 2], help="1 = gauge only, 2 = every run against the gauge")
    ap.add_argument("--floor-mm3", type=float, default=None, help="pass 2: volume (M m3) below which a storm is too small to time")
    ap.add_argument("--model-floor-mm3", type=float, default=None, help="pass 2, optional: also blank the timing of any run whose own volume (M m3) is below this")
    ap.add_argument("--root", default=str(HERE), help="the lab folder (default: the folder this script is in)")
    ap.add_argument("--results", default=None, help="folder with the outlet files (default <root>/results)")
    ap.add_argument("--table", default=None, help="storm_table.csv (default <root>/storm_table.csv)")
    ap.add_argument("--obs2017", default=None)
    ap.add_argument("--obs2020", default=None)
    ap.add_argument("--out-dir", default=None, help="output folder (default <root>/timing_out)")
    ap.add_argument("--windows", choices=["as-scored", "cut", "uncut"], default="as-scored")
    ap.add_argument("--stem", action="append", default=[], metavar="YEAR:KEY=NAME", help="other outlet-file name for one run, e.g. 2017:original=W5_xxx")
    ap.add_argument("--print-all", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--ignore-check-failure", action="store_true", help="diagnosis only")
    ap.add_argument("--definitions", action="store_true", help="print every definition and exit")
    ap.add_argument("--show-arrival", nargs="?", const="", default=None, metavar="PATH", help="print the first-arrival rule from score_run_sycamore.py and exit")
    return ap


def apply_stem_overrides(specs):
    runs = {y: list(v) for y, v in RUNS.items()}
    for s in specs:
        m = re.fullmatch(r"(20(?:17|20)):([A-Za-z0-9_]+)=(\S+)", s.strip())
        if not m:
            die(f"--stem '{s}' is not in the form 2017:original=NAME")
        yr, key, name = "WY" + m.group(1), m.group(2), m.group(3)
        if name.endswith("_Outlet.qout"):
            name = name[:-len("_Outlet.qout")]
        hit = [i for i, (k, _, _) in enumerate(runs[yr]) if k == key]
        if not hit:
            die(f"--stem: no run '{key}' for {yr}. Keys: {[k for k, _, _ in runs[yr]]}")
        k, l, _ = runs[yr][hit[0]]
        runs[yr][hit[0]] = (k, l, name)
    return runs


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    root = Path(a.root)
    if a.definitions:
        log("\n".join(definition_lines(a.windows, a.floor_mm3, a.model_floor_mm3)))
        return 0
    if a.show_arrival is not None:
        return show_arrival(Path(a.show_arrival) if a.show_arrival else root / "score_run_sycamore.py")
    if a.pass_ is None:
        die("choose --pass 1 or --pass 2 (or --show-arrival, or --definitions). See the top of this file.")
    if a.pass_ == 2 and a.floor_mm3 is None and not a.dry_run:
        die("pass 2 needs --floor-mm3 X: the volume (M m3) below which a storm is too small to time. Look at the pass 1 table and choose it; "
            "for example  python score_timing_sycamore.py --pass 2 --floor-mm3 0.05")
    if a.floor_mm3 is not None and not (math.isfinite(a.floor_mm3) and a.floor_mm3 >= 0):
        die("--floor-mm3 must be a number of 0 or more.")
    if a.model_floor_mm3 is not None:
        if not (math.isfinite(a.model_floor_mm3) and a.model_floor_mm3 >= 0):
            die("--model-floor-mm3 must be a number of 0 or more.")
        if a.floor_mm3 is None:
            die("--model-floor-mm3 is an extra on top of --floor-mm3 X; give both.")

    results = Path(a.results) if a.results else root / "results"
    table = Path(a.table) if a.table else root / "storm_table.csv"
    out_dir = Path(a.out_dir) if a.out_dir else root / "timing_out"
    if out_dir.resolve() == results.resolve():
        die("--out-dir is the results folder. This script never writes there.")
    runs = apply_stem_overrides(a.stem)
    if not results.is_dir():
        die(f"cannot find the results folder: {results}\n  Run this from workspaces/Sycamore/lab/, or point to it with --results.")

    log(f"\nscore_timing_sycamore v{VERSION}   pass {a.pass_}   windows: {a.windows}" + ("   [DRY RUN]" if a.dry_run else ""))
    # gauges
    gauges, gpaths = {}, {}
    for yr, spec in YEARS.items():
        given = a.obs2017 if yr == "WY2017" else a.obs2020
        cands = [Path(given)] if given else [root.parent / "init_data" / "obs" / n for n in spec["obs"]]
        hit = [p for p in cands if p.is_file()]
        if not hit:
            die(f"cannot find the {yr} gauge file. Looked for:\n    " + "\n    ".join(str(p) for p in cands)
                + "\n  Run this from workspaces/Sycamore/lab/, or give --obs2017 / --obs2020.")
        gpaths[yr] = hit[0]
        gauges[yr] = read_gauge(hit[0])
        g = gauges[yr]
        log(f"  {yr} gauge: {hit[0].name}: {int(g['q_cms'].notna().sum()):,} readings, {g.index.min():%Y-%m-%d} to {g.index.max():%Y-%m-%d}"
            + ("" if "tier" in g.columns else "   (no 'tier' column: gauge quality percentages will be blank)"))
    # outlet files
    models, missing = {"WY2017": {}, "WY2020": {}}, []
    for yr in ("WY2017", "WY2020"):
        for key, label, stem in runs[yr]:
            p = results / f"{stem}_Outlet.qout"
            if not p.is_file():
                missing.append((yr, key, label, p))
                continue
            q = read_qout(p, YEARS[yr]["start"])
            models[yr][key] = q
            log(f"  {yr} {label:<22} {stem}_Outlet.qout   {len(q):,} rows, last hour {(q.index.max() - YEARS[yr]['start']) / HOUR:.1f}")
    if missing:
        log("\nOUTLET FILES NOT FOUND (these runs are skipped):")
        for yr, key, label, p in missing:
            log(f"  {yr} {label:<22} {p.name}   (key {key})")
        have = sorted(f.name[:-len('_Outlet.qout')] for f in results.glob("*_Outlet.qout")) if results.is_dir() else []
        w5 = [h for h in have if h.startswith("W5_")]
        if w5 and any(y == "WY2017" for y, *_ in missing):
            log("  Outlet files in results/ that start with W5_ (use --stem 2017:original=NAME if one of these is the run):")
            for h in w5[:30]:
                log("    " + h)
    for yr in ("WY2017", "WY2020"):
        if not any(k in models[yr] for k in PRIORITY):
            die(f"neither Re-tuned nor Original was found for {yr}, and they are the priority runs. See the list above; "
                f"check the names or use --stem {yr[2:]}:retuned=NAME.")
    rain, rnotes = read_rain_table(table)
    for n in rnotes:
        log("  NOTE: " + n)
    wins = make_windows(a.windows)
    # sanity: storm_table agrees with the storms and last rain days used here
    for w in wins:
        r = rain.get((w.year, w.first.strftime("%Y-%m-%d")))
        if rain and r is None:
            log(f"  NOTE: storm_table.csv has no row for {w.year} {w.name} (first day {w.first:%Y-%m-%d}); its rain columns are blank.")
        elif r is not None and r["_last"] and r["_last"][:10] != w.last.strftime("%Y-%m-%d"):
            log(f"  NOTE: storm_table.csv gives the last rain day of {w.year} {w.name} as {r['_last']}, but this script uses {w.last:%Y-%m-%d}.")

    log("\nWINDOWS (start to end, hours):")
    for w in wins:
        log(f"  {w.year}  {w.name:<13} {w.lo:%Y-%m-%d %H:%M} to {w.hi:%Y-%m-%d %H:%M}  {w.hours:>4} h   {w.rule}")
    if a.dry_run:
        log("\nDRY RUN: every input was read and the windows listed. Nothing was computed or written.")
        return 0

    # 1. volume check
    log("\nVOLUME CHECK (recomputed volumes against the volumes on record; M m3; ratio = model / gauge)")
    status, bad, lines = volume_check(wins, gauges, models)
    for ln in lines:
        log(ln)
    if bad:
        log(f"\nVOLUME CHECK: MISMATCH ({len(bad)} value(s) differ by more than {TOL_M3} M m3):")
        for b in bad:
            log("  " + b)
        if not a.ignore_check_failure:
            die("the volumes do not reproduce, so nothing was written. Do not use any timing number until this is explained.")
        log("  --ignore-check-failure was given: carrying on. Rows are labeled MISMATCH.")
    else:
        n_m = sum(1 for v in status.values() if v == "MATCH")
        log(f"\nVOLUME CHECK: MATCH. {n_m} run-and-storm volumes (plus the gauge) reproduce the record within {TOL_M3} M m3; the rest have no reference (NOT CHECKED).")

    bundles = {w.key: gauge_bundle(w, gauges[w.year]) for w in wins}
    header = definition_lines(a.windows, a.floor_mm3 if a.pass_ == 2 else None, a.model_floor_mm3 if a.pass_ == 2 else None)
    header.append("VOLUME CHECK RESULT: " + ("MISMATCH, carried on with --ignore-check-failure" if bad else "MATCH") + ". Gauge files: "
                  + "; ".join(f"{y} {p.name}" for y, p in gpaths.items()) + ".")
    if missing:
        header.append("RUNS NOT FOUND (skipped): " + "; ".join(f"{y} {l}" for y, _, l, _ in missing))

    # 2. pass 1
    if a.pass_ == 1:
        rows = [pass1_row(w, bundles[w.key], status, rain) for w in wins]
        print_pass1(rows)
        path = out_dir / "pass1_gauge_table.csv"
        write_csv(path, header, pd.DataFrame(rows, columns=P1_COLS, dtype=object))
        log(f"\nWrote {len(rows)} rows to {path}")
        log("\nPASS 1 DONE. Stopping here as agreed: look at the table, choose the 'too small to time' floor, then run pass 2.")
        return 0

    # 3. pass 2
    rows = []
    for yr in ("WY2017", "WY2020"):
        for key, label, stem in runs[yr]:
            if key not in models[yr]:
                continue
            for w in [x for x in wins if x.year == yr]:
                rows.append(pass2_row(w, bundles[w.key], key, label, stem, models[yr][key], a.floor_mm3, status, rain, a.model_floor_mm3))
    print_pass2(rows, a.print_all)
    c5 = check5(rows)
    log("")
    for ln in c5:
        log(ln)
    n_blank = sum(1 for r in rows if r["blank_reasons"])
    log(f"\n{len(rows)} run-and-storm rows; {n_blank} have at least one blank value (reason in blank_reasons).")
    header.extend(c5)
    path = out_dir / "pass2_timing_long.csv"
    write_csv(path, header, pd.DataFrame(rows, columns=P2_COLS, dtype=object))
    log(f"Wrote {len(rows)} rows to {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
