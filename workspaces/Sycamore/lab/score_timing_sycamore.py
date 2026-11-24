#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
score_timing_sycamore.py   (v1.4, 2026-11-22)
=============================================

Timing and shape scores for the Lower Sycamore Wash storms: WY2017 (5 storms) and WY2020 (7 storms), 12 storms in all.

NEW IN v1.4.  Every definition (t10/t50/t90, peak and its flags, the 80 percent second-bump rule, tail, lag within +/-12 h,
the resolution note) is exactly as in v1.3.  What changed:
  1. WINDOWS.  --windows step0 is the default: the 12 storm windows are read from step0_windows.csv (water_year, storm, start, end,
     gauge_response), not derived from the first and last rain days.  Arizona local time, half-open [start, end), each hour named by its
     start.  The script STOPS if any window's length in hours differs from the required list (see STEP0_HOURS).  The old rule is still
     there: --windows as-scored (v1.3 behavior), cut, uncut.
  2. RUNS.  Re-tuned, Original and the eight Phase 3 screens, both years, all started 1 September (WY2017 from 2016-09-01 00:00,
     6,551 h; WY2020 from 2019-09-01 00:00, 6,575 h).  The start time and the expected length belong to each RUN, not to the year.
     The script STOPS (after printing the inventory) if a file's last Time_hr does not match its run's expected length.  The four
     Kchan runs are left out unless --include-kchan (the two WY2017 Kchan runs are the old 1 October whole-year runs, 8,759 h).
     --wy2017-runs oct1 uses the four old WY2017 whole-year runs (v1.3 files) in place of the 1 September runs.
  3. NO FLOOR.  --floor-mm3 is optional and off by default: every storm is scored and nothing is blanked for small volume.  The new
     column no_gauge_response (from the windows table) lets you take medians with and without the two no-response storms.
  4. VOLUME RATIO.  vol_ratio_gauge_hours = model volume over the gauge's hours with data / gauge volume (the thesis's volume ratio);
     vol_ratio_model_over_gauge (all model hours) is still there.  With --windows step0 the old reference volumes (tied to the old
     windows) are replaced by one check on Re-tuned and Original: vol_ratio_gauge_hours of the 10 storms with a gauge response must
     match the per-storm ratios on record within 0.01, or the script stops.  The record is read from window_scores_R4.csv (the
     window scorer's matched-hours volumes) when --ratio-source (or <root>/window_scores_R4.csv) is there; otherwise the typed values
     in RATIO_TYPED are used and the printout says so.
  5. The rain_* columns are DROPPED (storm_table.csv is matched on the old first rain days; it is no longer read).
  6. pass2_medians.csv: one row per run, measure and storm set (median, min, max, n); both years pooled, plus each year alone.
  7. --compare-old OLD_pass2_timing_long.csv: after a pass 2 run, compares the new rows with an old file line by line.

READ-ONLY.  It reads outlet files, the two gauge files and step0_windows.csv.  It runs no model, changes none of your files, and
writes only into its own output folder (timing_out_step0/ by default; timing_out_check_<mode>/ for the old window rules, so the old
folders are never touched).  It will not overwrite a CSV that a different version of this script wrote.  Every definition it uses
is written at the top of each CSV it makes (lines that start with "#"); --definitions prints them without running anything.

WHAT IT DOES, IN ORDER
  1. RUN INVENTORY AND WINDOW TABLE (always first): every run's file, its modified date, its start clock, first and last Time_hr and
     its expected length; then the window table (water year, storm, start, end, hours, gauge response).
  2. VOLUME CHECK (both passes): see item 4 above (step0), or the v1.3 check against the volumes on record (the old window rules).
     A mismatch STOPS the run and writes nothing.
  3. PASS 1  (--pass 1)   GAUGE ONLY.  One row per storm: window, gauge volume, gauge peak and its time, the gauge times by which 10,
     50 and 90 percent of the window volume had passed, flow at both window edges, gauge quality, outage and gap flags.
  4. PASS 2  (--pass 2)   every run against the gauge, per storm: volume-time differences, peak ratio and peak-time difference with
     flags, hourly Pearson r, cross-correlation lag within +/-12 h, tail volume ratio.  Writes pass2_timing_long.csv and
     pass2_medians.csv, prints the Re-tuned and Original medians, then CHECK 5.

COMMANDS (from workspaces/Sycamore/lab/; one at a time, foreground)
  python score_timing_sycamore.py --dry-run                   # inventory + window table + volume check; computes no scores, writes nothing
  python score_timing_sycamore.py --pass 2                    # the full run (no floor)
  python score_timing_sycamore.py --pass 1                    # gauge-only table
  python score_timing_sycamore.py --definitions               # prints every definition
  # reproduce the 20 Nov files (v1.3 behavior; the old folders are not touched):
  python score_timing_sycamore.py --pass 2 --windows as-scored --wy2017-runs oct1 --include-kchan --floor-mm3 0 \\
         --compare-old timing_out_all_storms/pass2_timing_long.csv
  (the 2.0 M m3 file: the same with --floor-mm3 2.0 and --compare-old timing_out/pass2_timing_long.csv)

OPTIONS YOU MAY NEED
  --windows step0|as-scored|cut|uncut   step0 (default): the table in step0_windows.csv.  as-scored: the v1.3 windows (WY2020 uncut,
                                  WY2017 cut at the next storm).  cut / uncut: one rule for every storm.
  --windows-file FILE             another windows table (default <root>/step0_windows.csv).
  --ratio-source FILE             window_scores CSV from score_windows_sycamore.py for the ratio check (default <root>/window_scores_R4.csv).
  --stem 2017:original=NAME       use another outlet-file name (without _Outlet.qout) for one run.  Repeat as needed.  The run keeps its
                                  start and expected length.  Years are 2017 or 2020; run keys are retuned, original, S1, S5 ... S11,
                                  kchan05, kchan02.
  --obs2017 FILE, --obs2020 FILE  other gauge files.   --results DIR   other outlet folder.   --out-dir DIR   other output folder.
  --print-all                     in pass 2, print every run (default prints Re-tuned and Original only).
  --ignore-check-failure          DIAGNOSIS ONLY: carry on after a volume-check mismatch or a length mismatch; every row then says MISMATCH.

WHAT IT READS (never changed)
  results/<run>_Outlet.qout   columns Time_hr, Q (m3/s).  Time_hr 0 = the run's own start (2016-09-01, 2016-10-01 or 2019-09-01 00:00).
  ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv                  datetime, q_cms [, tier]
  ../init_data/obs/sycamore_09510200_obs_flow_Sep2019_to_May2020.csv      (or ..._WY2020.csv if that one is missing)
  step0_windows.csv           the 12 windows (step0 only)
  window_scores_R4.csv        optional, the ratio check's record (step0 only)
  score_run_sycamore.py       only with --show-arrival

WORDS
  hourly mean      the average of all readings in one clock hour. Hours with no reading are left out, never filled, never zero.
  M m3             million cubic meters.      cfs   cubic feet per second (m3/s divided by 0.0283168).
  volume time t50  hours from the window start by which 50 percent of the window's volume had passed.
  Pearson r        the usual correlation between two hourly series (1 = same shape, 0 = unrelated).
  lag              the shift in hours that makes the model hourly series best match the gauge; positive = model later.
  vol_ratio_gauge_hours   model volume over the gauge's hours with data, divided by the gauge volume.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import io
import math
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = "1.4"
SELF_MARK = "score_timing_sycamore.py"
SELF_VERSION_MARK = f"{SELF_MARK} v{VERSION}"
HERE = Path(__file__).resolve().parent
CFS_PER_CMS = 1.0 / 0.028316846592
HOUR = pd.Timedelta(hours=1)
TOL_M3 = 0.005                 # volume-check tolerance of the old window rules, M m3
RATIO_TOL = 0.01               # volume-ratio check tolerance (step0)
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

SEP1_2016 = pd.Timestamp("2016-09-01 00:00")
OCT1_2016 = pd.Timestamp("2016-10-01 00:00")
SEP1_2019 = pd.Timestamp("2019-09-01 00:00")
HOURS_SEP1_2017, HOURS_OCT1_2017, HOURS_SEP1_2020 = 6551, 8759, 6575

YEARS = {
    "WY2017": dict(role="calibration", obs=["sycamore_09510200_obs_flow_WY2017.csv"]),
    "WY2020": dict(role="validation", obs=["sycamore_09510200_obs_flow_Sep2019_to_May2020.csv", "sycamore_09510200_obs_flow_WY2020.csv"]),
}

# ------------------------------------------------------------------------------------------------ runs
# Outlet-file names (without _Outlet.qout).  WY2020 and WY2017 September-start names are the run names in window_scores_R4.csv
# (score_windows_sycamore.py v1.2 bundle); the old WY2017 names are the 1 October whole-year runs of v1.3.
STEMS_2020 = {
    "retuned": "Sycamore_WY2020_Retuned", "original": "Sycamore_WY2020_Original",
    "kchan05": "Sycamore_WY2020_Original_kchan0p5", "kchan02": "Sycamore_WY2020_Original_kchan0p2",
    "S1": "Sycamore_WY2020_S1_Ks1p0", "S5": "Sycamore_WY2020_S5_Ks1p0_Bed2p0", "S6": "Sycamore_WY2020_S6_Ks0p5_Bed2p0",
    "S7": "Sycamore_WY2020_S7_Ks1p0_Bed2p0_As500", "S8": "Sycamore_WY2020_S8_Ks1p0_Bed2p0_f0p02", "S9": "Sycamore_WY2020_S9_Ks1p0_Bed0p5",
    "S10": "Sycamore_WY2020_S10_Ks1p0_GW0p5", "S11": "Sycamore_WY2020_S11_Ks1p0_PsiB1p0",
}
STEMS_2017_SEP1 = {
    "retuned": "Sycamore_WY2017sep_Retuned", "original": "Sycamore_WY2017sep_Original",
    "S1": "Sycamore_WY2017sep_S1_Ks1p0", "S5": "Sycamore_WY2017sep_S5_Ks1p0_Bed2p0", "S6": "Sycamore_WY2017sep_S6_Ks0p5_Bed2p0",
    "S7": "Sycamore_WY2017sep_S7_Ks1p0_Bed2p0_As500", "S8": "Sycamore_WY2017sep_S8_Ks1p0_Bed2p0_f0p02",
    "S9": "Sycamore_WY2017sep_S9_Ks1p0_Bed0p5", "S10": "Sycamore_WY2017sep_S10_Ks1p0_GW0p5", "S11": "Sycamore_WY2017sep_S11_Ks1p0_PsiB1p0",
}
STEMS_2017_OCT1 = {
    "retuned": "W5_d012_k10_f0p15_as300_ksx1p03_nocl", "original": "W5_d012_k10_base_nocl",
    "kchan05": "W5_d012_k10_base_kchan0p5", "kchan02": "W5_d012_k10_base_kchan0p2",
}
CORE_RUNS = [("retuned", "Re-tuned"), ("original", "Original")]
KCHAN_RUNS = [("kchan05", "Original + Kchan 0.5"), ("kchan02", "Original + Kchan 0.2")]
SCREEN_RUNS = [(1, "S1"), (2, "S5"), (3, "S6"), (4, "S7"), (5, "S8"), (6, "S9"), (7, "S10"), (8, "S11")]    # (Screen number, run key)
ALL_KEYS = [k for k, _ in CORE_RUNS] + [k for k, _ in KCHAN_RUNS] + [k for _, k in SCREEN_RUNS]
PRIORITY = ("retuned", "original")

# Design role of each RUN (pass 2 only) and whether its settings were written down BEFORE any WY2020 gauge volume or model
# result was opened.  (role, fixed_before_wy2020_results).  Sources: Tail_Screen_W4_Running_Log.md (Re-tuned frozen Nov 3 20:15,
# Original named in the same entry; Kchan 0.5 and 0.2 written down Nov 4 02:10, after the gauge volumes and the Re-tuned and
# Original results had been opened; their values came from the WY2017 probes) and the screen scripts (S1, S5 to S11 built from Nov 7 on).
# v1.4: the eight screens also exist for WY2017 (September-start runs built Nov 20-21); their settings are the same ones, written down
# after the WY2020 results had been opened, so they carry the same label ("post_validation_screen", "no").
RUN_DESIGN = {
    ("WY2017", "retuned"): ("calibration", "n/a"),
    ("WY2017", "original"): ("calibration", "n/a"),
    ("WY2017", "kchan05"): ("calibration", "n/a"),
    ("WY2017", "kchan02"): ("calibration", "n/a"),
    ("WY2020", "retuned"): ("validation", "yes"),
    ("WY2020", "original"): ("validation", "yes"),
    ("WY2020", "kchan05"): ("validation", "no"),
    ("WY2020", "kchan02"): ("validation", "no"),
    **{(y, k): ("post_validation_screen", "no") for y in ("WY2017", "WY2020") for k in ("S1", "S5", "S6", "S7", "S8", "S9", "S10", "S11")},
}
assert all((y, k) in RUN_DESIGN for y in ("WY2017", "WY2020") for k in ALL_KEYS), "RUN_DESIGN must cover every run key in both years"


def run_design(year: str, key: str):
    """(role, fixed_before_wy2020_results) for one run; a key that is not in RUN_DESIGN (only possible in a unit test) is 'unassigned'."""
    return RUN_DESIGN.get((year, key), ("unassigned", "n/a"))


class Run:
    """One outlet file: its key, label, stem, the clock its Time_hr 0 means, and the length it must have."""

    def __init__(self, year, key, label, stem, start, hours):
        self.year, self.key, self.label, self.stem = year, key, label, stem
        self.start, self.hours = pd.Timestamp(start), int(hours)

    @property
    def fname(self):
        return f"{self.stem}_Outlet.qout"


def build_runs(wy2017_runs: str = "sep1", include_kchan: bool = False, overrides=()):
    """The run list in v1.3 order (Re-tuned, Original, Kchan if asked, screens).  Start and expected length are set PER RUN."""
    runs = []
    for yr in ("WY2017", "WY2020"):
        wanted = list(CORE_RUNS) + (list(KCHAN_RUNS) if include_kchan else []) + [(k, f"Screen {n} ({k})") for n, k in SCREEN_RUNS]
        for key, label in wanted:
            if yr == "WY2020":
                runs.append(Run(yr, key, label, STEMS_2020[key], SEP1_2019, HOURS_SEP1_2020))
            elif key in ("kchan05", "kchan02"):
                runs.append(Run(yr, key, label, STEMS_2017_OCT1[key], OCT1_2016, HOURS_OCT1_2017))
            elif wy2017_runs == "oct1":
                if key in STEMS_2017_OCT1:
                    runs.append(Run(yr, key, label, STEMS_2017_OCT1[key], OCT1_2016, HOURS_OCT1_2017))
                # the old whole-year set has no screens: nothing to add
            else:
                runs.append(Run(yr, key, label, STEMS_2017_SEP1[key], SEP1_2016, HOURS_SEP1_2017))
    for s in overrides:
        m = re.fullmatch(r"(20(?:17|20)):([A-Za-z0-9_]+)=(\S+)", s.strip())
        if not m:
            die(f"--stem '{s}' is not in the form 2017:original=NAME")
        yr, key, name = "WY" + m.group(1), m.group(2), m.group(3)
        if name.endswith("_Outlet.qout"):
            name = name[:-len("_Outlet.qout")]
        hit = [r for r in runs if r.year == yr and r.key == key]
        if not hit:
            die(f"--stem: no run '{key}' for {yr} in this run list. Keys: {[r.key for r in runs if r.year == yr]}"
                + ("  (the Kchan runs need --include-kchan)" if key.startswith("kchan") else ""))
        hit[0].stem = name
    return runs


# ------------------------------------------------------------------------------------------------ windows
# Required window lengths (hours) of the Step 0 windows.  The script stops if step0_windows.csv gives anything else.
STEP0_HOURS = {
    ("WY2017", "Dec 21-25"): 188, ("WY2017", "Dec 30-Jan 2"): 153, ("WY2017", "Jan 19-21"): 85, ("WY2017", "Feb 18-19"): 135,
    ("WY2017", "Feb 27-28"): 121, ("WY2020", "Nov 19-22"): 160, ("WY2020", "Nov 28-30"): 141, ("WY2020", "Dec 23-28"): 210,
    ("WY2020", "Feb 10-11"): 132, ("WY2020", "Feb 22-23"): 121, ("WY2020", "Mar 10-13"): 187, ("WY2020", "Mar 18-19"): 120,
}
NO_RESPONSE = {("WY2017", "Dec 21-25"), ("WY2020", "Nov 19-22")}        # gauge_response = none
OUTAGE_STORM = ("WY2017", "Dec 30-Jan 2")

# (year, name, first rain day, last rain day, first day of the next storm in storm_table.csv or None) -- the v1.3 window rule
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

WINDOW_SET_NAMES = {
    "step0": "Step 0 windows (the 12 explicit windows of step0_windows.csv; not derived from rain days)",
    "as-scored": "old rule, as-scored (v1.3: WY2020 uncut, WY2017 cut at the next storm; 00:00 first rain day to 00:00 four days after the last rain day)",
    "cut": "old rule, cut (every storm cut at the next storm's first day)",
    "uncut": "old rule, uncut (no storm cut)",
}


class Win:
    """One scored window [lo, hi)."""

    def __init__(self, year, name, lo, hi, rule, gauge_response="yes", hi_uncut=None):
        self.year, self.name = year, name
        self.lo, self.hi = pd.Timestamp(lo), pd.Timestamp(hi)
        self.hours = int(round((self.hi - self.lo) / HOUR))
        self.rule = rule
        self.gauge_response = gauge_response
        self.hi_uncut = pd.Timestamp(hi_uncut) if hi_uncut is not None else self.hi
        self.role = YEARS[year]["role"]

    @property
    def key(self):
        return (self.year, self.name)

    @property
    def no_gauge_response(self):
        return "yes" if self.gauge_response == "none" else "no"


def legacy_windows(mode: str):
    """The v1.3 window rules: as-scored, cut, uncut."""
    out = []
    for y, n, f, l, nx in STORMS:
        first, last = pd.Timestamp(f), pd.Timestamp(l)
        hi_uncut = last + pd.Timedelta(days=4)
        hi_cut = min(hi_uncut, pd.Timestamp(nx)) if nx else hi_uncut
        hi = hi_uncut if mode == "uncut" else hi_cut if mode == "cut" else (hi_uncut if y == "WY2020" else hi_cut)
        rule = "cut at the next storm's first day" if hi < hi_uncut else "last rain day + 4 d"
        out.append(Win(y, n, first, hi, rule, "none" if (y, n) in NO_RESPONSE else "yes", hi_uncut=hi_uncut))
    return out


def load_step0_windows(path: Path):
    """Read step0_windows.csv; STOP on anything that is not exactly the 12 windows with the required lengths."""
    if not path.is_file():
        die(f"cannot find the windows table: {path}\n  It belongs next to this script (step0_windows.csv), or give --windows-file.")
    raw = pd.read_csv(path, comment="#", dtype=str, keep_default_na=False)
    raw.columns = [c.strip() for c in raw.columns]
    need = ["water_year", "storm", "start", "end", "gauge_response"]
    miss = [c for c in need if c not in raw.columns]
    if miss:
        die(f"{path.name} has no column {miss}. Columns found: {list(raw.columns)}")
    wins, problems, seen = [], [], set()
    for i, r in raw.iterrows():
        yr = "WY" + re.sub(r"\D", "", r["water_year"])[-4:]
        name = r["storm"].strip()
        key = (yr, name)
        if key in seen:
            problems.append(f"{yr} {name}: listed twice")
            continue
        seen.add(key)
        if key not in STEP0_HOURS:
            problems.append(f"{yr} {name}: not one of the 12 Step 0 storms")
            continue
        try:
            lo, hi = pd.Timestamp(r["start"].strip()), pd.Timestamp(r["end"].strip())
        except Exception:
            problems.append(f"{yr} {name}: start or end is not a date and time ({r['start']!r}, {r['end']!r})")
            continue
        if lo.minute or lo.second or hi.minute or hi.second:
            problems.append(f"{yr} {name}: start and end must be whole hours (each hour is named by its start); got {lo} and {hi}")
            continue
        gr = r["gauge_response"].strip().lower()
        if gr not in ("yes", "none"):
            problems.append(f"{yr} {name}: gauge_response must be yes or none, got {r['gauge_response']!r}")
            continue
        hours = int(round((hi - lo) / HOUR))
        if hours != STEP0_HOURS[key]:
            problems.append(f"{yr} {name}: {hours} h from {lo:%Y-%m-%d %H:%M} to {hi:%Y-%m-%d %H:%M}, but {STEP0_HOURS[key]} h is required")
            continue
        wins.append(Win(yr, name, lo, hi, "Step 0 window (step0_windows.csv)", gr))
    missing = [f"{y} {n}" for (y, n) in STEP0_HOURS if (y, n) not in seen]
    if missing:
        problems.append("missing windows: " + "; ".join(missing))
    if problems:
        die(f"{path.name} does not match the required Step 0 windows:\n  " + "\n  ".join(problems))
    order = {k: i for i, k in enumerate(STEP0_HOURS)}
    wins.sort(key=lambda w: order[w.key])
    return wins


def window_table_lines(wins):
    L = ["  {:<7}{:<14}{:<18}{:<18}{:>6}   {}".format("year", "storm", "start", "end", "hours", "gauge response")]
    for w in wins:
        L.append("  {:<7}{:<14}{:<18}{:<18}{:>6}   {}".format(w.year, w.name, f"{w.lo:%Y-%m-%d %H:%M}", f"{w.hi:%Y-%m-%d %H:%M}", w.hours, w.gauge_response))
    return L


# ------------------------------------------------------------------------------------------------ reference volumes of the OLD window rules (M m3)
# Source labels: W = WY2020 write-up (as coded in add_wy2020_results_to_storm_table.py); T = storm_table.csv (WY2017 columns);
# L = running log (Soil_Storage_Running_Log / part 2).  Order of the WY2020 tuples: gauge, retuned, original, kchan05, kchan02.
# Used ONLY with --windows as-scored|cut|uncut.  (They are tied to the old windows and the old Oct-1 WY2017 runs.)
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
    """(year, storm) -> list of (window_end, {series_key: (value, source)}, gauge_hours_text or None).  Old window rules only."""
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

# ------------------------------------------------------------------------------------------------ volume-ratio record (step0)
# vol_ratio_gauge_hours per storm on record, Re-tuned and Original, for the 10 storms with a gauge response.  Typed by Tracy
# (rounded to 2 decimals); they equal 1 + PBIAS/100 of window_scores_R4.csv, i.e. vol_model_Mm3 / vol_gauge_Mm3 on matched hours.
RATIO_STORMS = [("WY2017", "Dec 30-Jan 2"), ("WY2017", "Jan 19-21"), ("WY2017", "Feb 18-19"), ("WY2017", "Feb 27-28"),
                ("WY2020", "Nov 28-30"), ("WY2020", "Dec 23-28"), ("WY2020", "Feb 10-11"), ("WY2020", "Feb 22-23"),
                ("WY2020", "Mar 10-13"), ("WY2020", "Mar 18-19")]
RATIO_TYPED = {
    "retuned": (0.58, 1.32, 1.90, 0.92, 2.22, 0.05, 4.28, 1.16, 0.15, 0.31),
    "original": (0.78, 1.42, 1.70, 0.84, 2.43, 0.04, 4.27, 1.20, 0.08, 0.25),
}
RATIO_SETTING_NAMES = {"Retuned": "retuned", "Original": "original"}      # the setting names in window_scores_R4.csv


def read_ratio_source(path: Path):
    """{(run key, year, storm): vol_model_Mm3 / vol_gauge_Mm3} from a window_scores CSV (score_windows_sycamore.py), storm rows only."""
    try:
        df = pd.read_csv(path, comment="#", dtype=str, keep_default_na=False)
    except Exception as exc:
        die(f"cannot read the ratio source {path}: {exc}")
    need = ["setting", "kind", "name", "years", "vol_gauge_Mm3", "vol_model_Mm3"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        die(f"the ratio source {path.name} has no column {miss}. Columns found: {list(df.columns)[:12]}")
    out = {}
    for _, r in df[df["kind"] == "storm"].iterrows():
        key = RATIO_SETTING_NAMES.get(r["setting"])
        if key is None:
            continue
        try:
            g, m = float(r["vol_gauge_Mm3"]), float(r["vol_model_Mm3"])
        except ValueError:
            continue
        if g > 0:
            out[(key, "WY" + r["years"].strip(), r["name"].strip())] = m / g
    return out


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


def make_windows(mode: str, windows_file=None):
    return legacy_windows(mode) if mode != "step0" else load_step0_windows(Path(windows_file))


# ------------------------------------------------------------------------------------------------ readers
def read_qout_frame(path: Path) -> pd.DataFrame:
    """Time_hr and Q of an outlet file (comma or whitespace separated), numeric rows only."""
    with open(path) as fh:
        fh.readline()
        second = fh.readline()
    if "," in second:
        raw = pd.read_csv(path, sep=",", header=0, float_precision="round_trip")
    else:
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None)
    raw = raw.iloc[:, :2].copy()
    raw.columns = ["Time_hr", "Q"]
    return raw.apply(pd.to_numeric, errors="coerce").dropna()


def qout_series(raw: pd.DataFrame, start: pd.Timestamp) -> pd.Series:
    s = pd.Series(raw["Q"].to_numpy(float), index=start + pd.to_timedelta(raw["Time_hr"].to_numpy(float), unit="h"))
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_qout(path: Path, start: pd.Timestamp) -> pd.Series:
    return qout_series(read_qout_frame(path), start)


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


def fmt_mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")


def run_file_info(run: Run, path: Path, raw: pd.DataFrame) -> dict:
    """First and last Time_hr of a run's file, its row step, and whether its length matches the run's expected hours."""
    t = raw["Time_hr"].to_numpy(float)
    first, last = float(t.min()), float(t.max())
    step = float(np.median(np.diff(np.sort(t)))) if len(t) > 1 else float("nan")
    tol = max(0.1, 1.5 * step) if math.isfinite(step) else 0.1
    ok = (first <= tol) and abs(last - run.hours) <= tol
    why = ""
    if not ok:
        why = f"last Time_hr {last:g}, first {first:g}; {run.hours} h expected"
        if abs(last - HOURS_OCT1_2017) <= tol and run.hours != HOURS_OCT1_2017:
            why += " (8,759 h is the old 1 October whole-year run: wrong file for a 1 September run)"
        elif abs(last - HOURS_SEP1_2017) <= tol and run.hours != HOURS_SEP1_2017:
            why += " (6,551 h is a 1 September WY2017 run)"
        elif abs(last - HOURS_SEP1_2020) <= tol and run.hours != HOURS_SEP1_2020:
            why += " (6,575 h is a 1 September WY2020 run)"
    return {"path": path, "modified": fmt_mtime(path), "first": first, "last": last, "rows": len(t), "step": step, "ok": ok, "why": why}


def load_models(runs, results: Path):
    """Read every run's outlet file. Returns (models {year: {key: Series}}, infos {(year, key): dict}, missing list)."""
    models, infos, missing = {"WY2017": {}, "WY2020": {}}, {}, []
    for r in runs:
        p = results / r.fname
        if not p.is_file():
            missing.append(r)
            continue
        raw = read_qout_frame(p)
        if len(raw) == 0:
            infos[(r.year, r.key)] = {"path": p, "modified": fmt_mtime(p), "first": float("nan"), "last": float("nan"), "rows": 0,
                                      "step": float("nan"), "ok": False, "why": "no numeric rows"}
            continue
        infos[(r.year, r.key)] = run_file_info(r, p, raw)
        models[r.year][r.key] = qout_series(raw, r.start)
    return models, infos, missing


def inventory_lines(runs, infos):
    L = ["  {:<7}{:<9}{:<22}{:<42}{:<17}{:<17}{:>10}{:>11}{:>9}  {}".format(
        "year", "key", "label", "stem", "modified", "Time_hr 0 =", "first hr", "last hr", "expect", "status")]
    for r in runs:
        i = infos.get((r.year, r.key))
        if i is None:
            L.append("  {:<7}{:<9}{:<22}{:<42}{:<17}{:<17}{:>10}{:>11}{:>9}  {}".format(
                r.year, r.key, r.label, r.stem, "-", f"{r.start:%Y-%m-%d %H:%M}", "-", "-", r.hours, "NOT FOUND"))
        else:
            L.append("  {:<7}{:<9}{:<22}{:<42}{:<17}{:<17}{:>10}{:>11}{:>9}  {}".format(
                r.year, r.key, r.label, r.stem, i["modified"], f"{r.start:%Y-%m-%d %H:%M}", fnum(i["first"], 4), fnum(i["last"], 4),
                r.hours, "ok" if i["ok"] else "LENGTH DOES NOT MATCH: " + i["why"]))
    return L


def candidate_lines(results: Path, runs):
    """Outlet files in results/ that could be a WY2017 run (names with WY2017, or the old W5_ runs), with what they are used as."""
    used = {r.stem: f"{r.year} {r.key}" for r in runs}
    cands = sorted(f for f in results.glob("*_Outlet.qout") if "WY2017" in f.name or f.name.startswith("W5_") or used.get(f.name[:-len("_Outlet.qout")], "").startswith("WY2017"))
    L = []
    for f in cands:
        stem = f.name[:-len("_Outlet.qout")]
        try:
            raw = read_qout_frame(f)
            fl = f"first {raw['Time_hr'].min():g}, last {raw['Time_hr'].max():g}, {len(raw):,} rows"
        except Exception as exc:
            fl = f"could not be read ({exc})"
        L.append(f"    {stem:<46} {fmt_mtime(f)}   {fl}   " + (f"<- used as {used[stem]}" if stem in used else "(not in the run list)"))
    return L


# ------------------------------------------------------------------------------------------------ series numbers
def hourly(s: pd.Series, lo, hi) -> pd.Series:
    x = s[(s.index >= lo) & (s.index < hi)]
    return x.resample("60min").mean().dropna()


def vol_mm3(h: pd.Series) -> float:
    return float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")


def vol_over_gauge_hours(gh: pd.Series, mh: pd.Series) -> float:
    """Model volume (M m3) over the hours the gauge has data."""
    return vol_mm3(mh.reindex(gh.index).dropna())


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


# ------------------------------------------------------------------------------------------------ volume check, old window rules (v1.3)
def volume_check(wins, gauges, models, tol=TOL_M3, skip_models=frozenset()):
    """Old window rules only.  Returns (status dict {(year, storm, series_key): 'MATCH'|'MISMATCH'|'NOT CHECKED'}, mismatches list, printable lines).
    skip_models: {(year, key)} whose model volumes are not on record (e.g. the 1 September WY2017 runs): their volumes are NOT CHECKED."""
    status, bad, lines = {}, [], []
    here = {}
    for w in wins:
        gh = hourly(gauges[w.year]["q_cms"], w.lo, w.hi)
        here[(w.year, w.name, "gauge")] = vol_mm3(gh)
        here[(w.year, w.name, "gauge_hours")] = f"{len(gh)}/{w.hours}"
        for key, q in models[w.year].items():
            here[(w.year, w.name, key)] = vol_mm3(hourly(q, w.lo, w.hi))
    labels = {k: l for k, l in CORE_RUNS + KCHAN_RUNS}
    labels.update({k: f"Screen {n} ({k})" for n, k in SCREEN_RUNS})
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
                if ref is not None and key in ref[1] and (yr, key) not in skip_models:
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
    if "retuned" in models["WY2017"] and ("WY2017", "retuned") not in skip_models and any(r[0] == w.hi for r in REFS[w.key]):
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
    skipped = sorted({f"{y} {k}" for (y, k) in skip_models if k in models.get(y, {})})
    if skipped:
        lines.append("    Model volumes NOT CHECKED for " + ", ".join(skipped) + ": the volumes on record belong to the old 1 October runs (--wy2017-runs oct1).")
    return status, bad, lines


# ------------------------------------------------------------------------------------------------ volume-ratio check (step0)
def ratio_check(wins, gauges, models, source=None, tol=RATIO_TOL):
    """Step 0 check: vol_ratio_gauge_hours of Re-tuned and Original on the 10 storms with a gauge response against the record.
    source = {(key, year, storm): ratio} from window_scores_R4.csv, or None to use RATIO_TYPED.
    Returns (status dict, mismatches list, printable lines)."""
    status, bad, lines = {}, [], []
    by = {w.key: w for w in wins}
    lines.append("    {:<7}{:<14}{:>9}{:>9}{:>9}   {:>9}{:>9}{:>9}   {}".format(
        "year", "storm", "Re-tuned", "record", "diff", "Original", "record", "diff", "result"))
    for i, (yr, nm) in enumerate(RATIO_STORMS):
        w = by[(yr, nm)]
        gh = hourly(gauges[yr]["q_cms"], w.lo, w.hi)
        gv = vol_mm3(gh)
        cells, res = [], []
        for key in PRIORITY:
            mq = models[yr][key]
            mh = hourly(mq, w.lo, w.hi)
            got = vol_over_gauge_hours(gh, mh) / gv if (isnum(gv) and gv > 0 and len(mh) == w.hours) else float("nan")
            exp = source[(key, yr, nm)] if source is not None and (key, yr, nm) in source else RATIO_TYPED[key][i] if source is None else float("nan")
            ok = isnum(got) and isnum(exp) and abs(got - exp) <= tol
            status[(yr, nm, key)] = "MATCH" if ok else "MISMATCH"
            res.append(ok)
            cells.append((got, exp))
            if not ok:
                bad.append(f"{yr} {nm} {key}: vol_ratio_gauge_hours here {fnum(got, 3)} vs record {fnum(exp, 3)} (tolerance {tol:g})")
        lines.append("    {:<7}{:<14}{:>9}{:>9}{:>9}   {:>9}{:>9}{:>9}   {}".format(
            yr, nm, fnum(cells[0][0], 3), fnum(cells[0][1], 3), fnum(cells[0][0] - cells[0][1], 3) if isnum(cells[0][0]) and isnum(cells[0][1]) else "-",
            fnum(cells[1][0], 3), fnum(cells[1][1], 3), fnum(cells[1][0] - cells[1][1], 3) if isnum(cells[1][0]) and isnum(cells[1][1]) else "-",
            "ok" if all(res) else "MISMATCH"))
    return status, bad, lines


# ------------------------------------------------------------------------------------------------ definitions
def definition_lines(mode: str, floor, model_floor=None, ratio_source_text=""):
    m = {"step0": "the 12 storm windows of the Step 0 window finder, read from step0_windows.csv: Arizona local time, half-open [start, end), each hour named by its start; the two storms that the named storm days join (21-25 Dec 2016 and 23-28 Dec 2019) are as the finder joined them; nothing is derived from rain days here",
         "as-scored": "WY2020 uncut (exactly as score_wy2020_run_sycamore.py: Nov 28-30 runs to Dec 4); WY2017 cut at the next storm's first day when that comes sooner (only Jan 19-21 is affected: it ends Jan 23, not Jan 25)",
         "cut": "every storm cut at the next storm's first day when that comes sooner (WY2020 Nov 28-30 ends Dec 3; WY2017 Jan 19-21 ends Jan 23)",
         "uncut": "no storm is cut (WY2020 Nov 28-30 ends Dec 4; WY2017 Jan 19-21 ends Jan 25, which includes the Jan 23-24 storm)"}[mode]
    window_line = (f"WINDOW: Step 0 window of each storm (--windows step0): start and end are read from step0_windows.csv, [start, end), the hour named by its start. Rule used in this run: {m}."
                   if mode == "step0" else
                   "WINDOW: starts 00:00 on the first rain day; ends 00:00 four days after the last rain day. Rule used in this run (--windows " + mode + "): " + m + ".")
    L = [
        f"{SELF_MARK} v{VERSION} (2026-11-22). READ-ONLY: reads outlet files, gauge files" + (" and step0_windows.csv" if mode == "step0" else "") + "; changes none of them.",
        "WINDOW SET USED IN THIS FILE: " + WINDOW_SET_NAMES[mode] + ".",
        "RAIN COLUMNS: the rain_* columns of v1.3 are DROPPED in v1.4. storm_table.csv is matched on the old first rain days, which are not these windows, so it is no longer read.",
        "RUN START TIMES: set per RUN, not per year. Time_hr 0 = 2016-09-01 00:00 (WY2017 September-start runs, 6,551 h), 2016-10-01 00:00 (the old WY2017 whole-year runs, 8,759 h) or 2019-09-01 00:00 (WY2020 runs, 6,575 h). The script stops if a file's last Time_hr does not match its run's length.",
        "ROLES: in the pass 2 CSV the role belongs to the RUN: calibration = the WY2017 Re-tuned and Original runs (and the two WY2017 Kchan runs when included); validation = the WY2020 frozen setting (Re-tuned) and its three comparison runs (Original, Original + Kchan 0.5, Original + Kchan 0.2); post_validation_screen = the eight one-change screens S1 and S5 to S11, built after the WY2020 gauge and the Re-tuned and Original results had been opened, so they are NOT out of sample (the same label is given to their WY2017 runs). In the pass 1 CSV (gauge only, one row per storm) the role is the storm's year: WY2017 calibration, WY2020 validation. Times are Arizona local time, as in the gauge and model files.",
        "fixed_before_wy2020_results (pass 2 only): yes = the run's settings were written down before any WY2020 gauge volume or model result was opened (Re-tuned, frozen Nov 3 20:15; Original, named as a comparison in the same entry); no = written down or built after results had been opened (Kchan 0.5 and 0.2, Nov 4 02:10, with values taken from the WY2017 probes; the eight screens, Nov 7 onward); n/a = WY2017 Re-tuned, Original and Kchan runs (calibration year). Source: Tail_Screen_W4_Running_Log.md and the screen scripts; this script does not read them, the labels are fixed in RUN_DESIGN.",
        "HOURLY MEAN: average of all readings in one clock hour [hh:00, hh+1:00), m3/s. Hours with no reading are left out (never filled, never zero). Same as score_wy2020_run_sycamore.py.",
        "VOLUME: sum of hourly means x 3600 s over the hours that have data, M m3 (million cubic meters). A gauge window with missing hours gives a minimum volume. The model volume counts all its hours.",
        "VOLUME RATIO (v1.4): vol_ratio_gauge_hours = model volume over the gauge's hours with data (model_vol_gauge_hours_Mm3) / gauge volume (the thesis's volume ratio). vol_ratio_model_over_gauge = model volume over ALL its hours in the window / gauge volume (as in v1.3). They differ only where the gauge has missing hours (Dec 30 2016-Jan 2 2017).",
        window_line,
        "VOLUME TIME tNN (NN = 10, 50, 90): hours from the window start by which NN percent of the window volume had passed. Cumulative sum of the hourly volumes, each hour's volume spread evenly over its own hour, linear inside the hour where the target is crossed. Gauge hours with no data add no volume.",
        "TIME DIFFERENCES are always model minus gauge, in hours: positive = the model is later than the gauge.",
        "RESOLUTION: the gauge is read every 15 min and the model written every 3.75 min, so an hourly mean of gauge readings is centered about 5.6 min (0.09 h) earlier in true time than the model's. In the volume-time columns the model therefore looks about 0.09 h earlier than it is (checked on synthetic data). Differences smaller than about 0.1 h are inside this convention (the one the existing scorers use). The peak, lag-free arrival and peak-time columns use the raw readings and do not have this offset, but are limited by the 15-min gauge spacing.",
        "PEAK: the highest single reading in the window (raw: gauge every 15 min, model every output step), shown in cfs (m3/s / 0.0283168); the first one if several are equal. PEAK (HOURLY) = highest hourly mean, shown only so it can be compared with earlier hourly-based figures. Peak ratio = model peak / gauge peak.",
        f"PEAK FLAGS: (1) the highest reading is within {EDGE_H:g} h of a window edge; (2) a second bump at least {int(BUMP_FRAC * 100)}% as high as the highest hourly mean, with the flow between the two bumps falling to {int(DIP_FRAC * 100)}% of the lower bump or less (checked on hourly means; hours with no data are skipped); (3) flat top: readings within 1% of the peak span more than {FLAT_SPAN_H:g} h. A flagged peak time is still shown, not blanked.",
        "PEARSON r: ordinary correlation of the gauge and model hourly means over the hours both have, in the window. n_hours_r = how many hours.",
        f"LAG: shift in whole hours, within +/-{MAX_LAG_H} h, at which the model hourly series correlates best with the gauge hourly series (positive = model later; ties go to the smaller shift). Only shifts with at least {MIN_COMMON_H} common hours count. xcorr_at_limit = the best lag is at +/-{MAX_LAG_H}.",
        "TAIL: the tail starts at the GAUGE's t90 (so both series use the same clock interval) and runs to the window end. Tail volume ratio = model volume in the tail / gauge volume in the tail. model_tail_fraction = model volume in the tail / the model's own window volume (the gauge's is 0.10 by construction). Hours with no gauge data add no tail volume.",
        "BLANKS: an undefined value is left blank and the reason is in blank_reasons; never zero.",
        ("FLOOR: " + (f"{floor:g} M m3 (your choice), applied to the GAUGE window volume only. If the gauge window volume is below it, the timing columns (t10, t50, t90 differences, peak-time difference, lag, tail, arrival extras) are blank with that reason; volumes, peak ratio and Pearson r are still shown. A run whose own volume is small still gets timing numbers: read them next to vol_ratio_model_over_gauge (the printed table adds a note when the model has under " + f"{LOW_VOL_RATIO:.0%} of the gauge volume; that note is display only)." + ((" MODEL FLOOR " + f"{model_floor:g} M m3 (your choice, --model-floor-mm3): a run whose own window volume is below it also gets blank timing columns; below_floor says which side.") if model_floor is not None else " No model floor was set.") if floor is not None else "none. Every storm is scored and nothing is blanked for small volume (v1.4 default; --floor-mm3 is optional). Read the timing columns next to vol_ratio_gauge_hours and no_gauge_response (the printed table adds a note when the model has under " + f"{LOW_VOL_RATIO:.0%} of the gauge volume; that note is display only).")),
        "no_gauge_response: yes for the storms the windows table marks gauge_response = none (21-25 Dec 2016 and 19-22 Nov 2019), otherwise no; use it to take medians with and without those two storms.",
        "EXTRA COLUMNS ON A FRACTION OF THE GAUGE PEAK (labeled that way, not used in anything else): extra_arr20_of_gauge_peak_diff_h and extra_arr50_of_gauge_peak_diff_h = model minus gauge time of the first rise through 20% / 50% of the GAUGE's highest reading, linear between readings, not trusted across a reading gap longer than 3 h.",
        "EDGE FLOWS: hourly mean of the first hour and of the last hour of the window, cfs.",
        "GAUGE QUALITY: from the 'tier' column of the gauge file. pct_interp_or_est = readings tagged estimated or filled_interp as a percent of the window's readings that have a value; pct_filled_zero and pct_provisional likewise; pct_missing_hours = hours with no gauge data / hours in the window. Blank if the file has no 'tier' column.",
        "OUTAGE: jan1_2_outage_flag is set when the window overlaps 2017-01-01 00:00 to 2017-01-03 00:00. gauge_missing_intervals lists every run of hours with no gauge data in the window (this also catches the WY2020 dropouts). gap_before_t90 = yes when any missing gauge hour starts before the gauge's t90, so the gauge volume times are then biased. The outage storm (30 Dec 2016-2 Jan 2017) is scored.",
        ("VOLUME CHECK (Step 0): before anything else, for Re-tuned and Original in both years, vol_ratio_gauge_hours of the 10 storms with a gauge response is recomputed and compared with the record, tolerance " + f"{RATIO_TOL:g}" + "; MATCH or the script stops. Record used: " + (ratio_source_text or "the typed per-storm ratios (RATIO_TYPED)") + ". volume_check = MATCH / MISMATCH / NOT CHECKED per run and storm.")
        if mode == "step0" else
        "VOLUME CHECK: before anything else the volumes already on record (write-up, storm_table.csv, running log) are recomputed; volume_check = MATCH / MISMATCH / NOT CHECKED per run and storm (tolerance " + f"{TOL_M3}" + " M m3). They belong to the old windows and the old Oct-1 WY2017 runs.",
    ]
    return L


def guard_overwrite(path: Path):
    """Stop if path exists and was not written by this version of the script."""
    if path.exists():
        with open(path, encoding="utf-8") as fh:
            first = fh.readline()
        if SELF_VERSION_MARK not in first:
            die(f"{path} already exists and was not written by {SELF_VERSION_MARK}. Move it or use --out-dir.")


def write_csv(path: Path, header_lines, df: pd.DataFrame):
    guard_overwrite(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        for ln in header_lines:
            fh.write("# " + ln + "\n")
        df.to_csv(fh, index=False)


def r4(x):
    return round(float(x), 4) if isnum(x) else None


# ------------------------------------------------------------------------------------------------ pass 1
P1_COLS = ["water_year", "role", "storm", "no_gauge_response", "window_rule", "window_start", "window_end", "window_hours",
           "gauge_vol_Mm3", "gauge_volume_check", "gauge_hours_with_data", "gauge_peak_cfs", "gauge_peak_time",
           "gauge_peak_hourly_cfs", "gauge_t10_h", "gauge_t50_h", "gauge_t90_h", "gauge_left_edge_cfs", "gauge_right_edge_cfs",
           "pct_interp_or_est", "pct_filled_zero", "pct_provisional", "pct_missing_hours", "jan1_2_outage_flag",
           "gauge_missing_intervals", "gap_before_t90", "gauge_peak_flag", "blank_reasons"]


def pass1_row(w: Win, G: dict, status):
    q = G["quality"]
    why = list(G["why"])
    t = G["times"]
    row = {
        "water_year": w.year, "role": w.role, "storm": w.name, "no_gauge_response": w.no_gauge_response, "window_rule": w.rule,
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
P2_COLS = ["water_year", "role", "fixed_before_wy2020_results", "run_key", "run_label", "run_stem", "storm", "no_gauge_response",
           "window_rule", "window_start", "window_end",
           "volume_check", "gauge_vol_Mm3", "model_vol_Mm3", "model_vol_gauge_hours_Mm3", "vol_ratio_model_over_gauge", "vol_ratio_gauge_hours", "below_floor",
           "gauge_t10_h", "model_t10_h", "t10_model_minus_gauge_h",
           "gauge_t50_h", "model_t50_h", "t50_model_minus_gauge_h",
           "gauge_t90_h", "model_t90_h", "t90_model_minus_gauge_h",
           "gauge_peak_cfs", "model_peak_cfs", "peak_ratio_model_over_gauge", "peak_hourly_ratio_model_over_gauge",
           "gauge_peak_time", "model_peak_time", "peak_time_model_minus_gauge_h", "gauge_peak_flag", "model_peak_flag",
           "pearson_r_hourly", "n_hours_r", "xcorr_best_lag_h", "xcorr_best_r", "xcorr_at_limit",
           "tail_gauge_Mm3", "tail_model_Mm3", "tail_volume_ratio", "model_tail_fraction",
           "extra_arr20_of_gauge_peak_diff_h", "extra_arr50_of_gauge_peak_diff_h",
           "jan1_2_outage_flag", "gauge_missing_intervals", "gap_before_t90", "blank_reasons"]


def pass2_row(w: Win, G: dict, key, label, stem, mq: pd.Series, floor, status, model_floor=None):
    R = []
    row = {c: None for c in P2_COLS}
    role, fixed = run_design(w.year, key)
    row.update({"water_year": w.year, "role": role, "fixed_before_wy2020_results": fixed, "run_key": key, "run_label": label, "run_stem": stem, "storm": w.name,
                "no_gauge_response": w.no_gauge_response,
                "window_rule": w.rule, "window_start": f"{w.lo:%Y-%m-%d %H:%M}", "window_end": f"{w.hi:%Y-%m-%d %H:%M}",
                "volume_check": status.get((w.year, w.name, key), "NOT CHECKED"),
                "jan1_2_outage_flag": "touches the Jan 1-2 2017 outage" if G["touches_outage"] else "",
                "gauge_missing_intervals": fmt_runs(G["runs"]), "gap_before_t90": G["gap_before_t90"]})
    gh = G["h"]
    mh = hourly(mq, w.lo, w.hi)
    gtot, mtot = G["tot"], vol_mm3(mh)
    row["gauge_vol_Mm3"], row["model_vol_Mm3"] = r4(gtot), r4(mtot)
    if len(mh) < w.hours:
        R.append(f"model: the outlet file has hourly means for only {len(mh)} of the window's {w.hours} h; every model-based value is blank")
        row["blank_reasons"] = " | ".join(R)
        return row
    mgh = vol_over_gauge_hours(gh, mh)
    row["model_vol_gauge_hours_Mm3"] = r4(mgh)
    if isnum(gtot) and gtot > 0:
        row["vol_ratio_model_over_gauge"] = r4(mtot / gtot)
        row["vol_ratio_gauge_hours"] = r4(mgh / gtot)
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
    # floor gate for the timing block (only when a floor was asked for)
    gate = None
    below = []
    if floor is not None and (not isnum(gtot) or gtot < floor):
        below.append("gauge")
    if model_floor is not None and (not isnum(mtot) or mtot < model_floor):
        below.append("model")
    if below:
        row["below_floor"] = " and ".join(below)
        gate = "; ".join(f"{s} volume {fnum(gtot if s == 'gauge' else mtot, 4)} M m3 is below the {'floor' if s == 'gauge' else 'model floor'} {(floor if s == 'gauge' else model_floor):g}" for s in below)
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
        if r["no_gauge_response"] == "yes":
            fl.append("NO GAUGE RESPONSE")
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


def check5(rows, mode):
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
    if mode == "step0":
        lines.append("  (Step 0 window: Feb 27 13:00 to Mar 4 14:00; the old figure used the old window, which started Feb 27 00:00.)")
    return lines


# ------------------------------------------------------------------------------------------------ medians file
MEASURES = [("peak ratio", "peak_ratio_model_over_gauge"), ("peak time (h)", "peak_time_model_minus_gauge_h"),
            ("t10 (h)", "t10_model_minus_gauge_h"), ("t50 (h)", "t50_model_minus_gauge_h"), ("t90 (h)", "t90_model_minus_gauge_h"),
            ("lag (h)", "xcorr_best_lag_h"), ("tail ratio", "tail_volume_ratio"), ("Pearson r", "pearson_r_hourly")]
STORM_SETS = [("a", "storms with a gauge response"),
              ("b", "storms with a gauge response, without 30 Dec 2016-2 Jan 2017 (outage)"),
              ("c", "all 12 storms, where defined")]
SCOPES = [("both years", None), ("WY2017", "WY2017"), ("WY2020", "WY2020")]
P3_COLS = ["run_key", "run_label", "years", "years_scored", "storm_set", "storm_set_description", "measure", "median", "min", "max", "n", "n_storms_in_set"]


def storm_set_mask(df: pd.DataFrame, sid: str):
    if sid == "a":
        return df["no_gauge_response"] == "no"
    if sid == "b":
        return (df["no_gauge_response"] == "no") & ~((df["water_year"] == OUTAGE_STORM[0]) & (df["storm"] == OUTAGE_STORM[1]))
    return pd.Series(True, index=df.index)


def build_medians(rows):
    """One row per run, years, storm set and measure: median, min, max over the storms where the measure is defined; n = how many."""
    df = pd.DataFrame(rows)
    out = []
    for key in list(dict.fromkeys(df["run_key"])):
        sub = df[df["run_key"] == key]
        label = sub["run_label"].iloc[0]
        for scope, yr in SCOPES:
            s2 = sub if yr is None else sub[sub["water_year"] == yr]
            if s2.empty:
                continue
            scored = "+".join(sorted(set(sub["water_year"])))
            for sid, desc in STORM_SETS:
                s3 = s2[storm_set_mask(s2, sid)]
                for mlabel, col in MEASURES:
                    v = pd.to_numeric(s3[col], errors="coerce").dropna()
                    out.append({"run_key": key, "run_label": label, "years": scope, "years_scored": scored, "storm_set": sid,
                                "storm_set_description": desc, "measure": mlabel,
                                "median": r4(v.median()) if len(v) else None, "min": r4(v.min()) if len(v) else None,
                                "max": r4(v.max()) if len(v) else None, "n": int(len(v)), "n_storms_in_set": int(len(s3))})
    return out


def print_medians(med, keys=PRIORITY):
    log("\nMEDIANS (pass2_medians.csv rows for Re-tuned and Original; model minus gauge in hours; * = fewer storms than in the set had the value defined, n in brackets)")
    mlabels = [m for m, _ in MEASURES]
    for key in keys:
        rows = [m for m in med if m["run_key"] == key]
        if not rows:
            continue
        log(f"\n  {rows[0]['run_label']}   (years scored: {rows[0]['years_scored']})")
        for what in ("median", "min .. max"):
            log("  " + what)
            wd = 16 if what == "median" else 22
            log("    {:<11}{:<5}{:>6}".format("years", "set", "storms") + "".join(f"{m:>{wd}}" for m in mlabels))
            for scope, _ in SCOPES:
                for sid, _ in STORM_SETS:
                    cells = [m for m in rows if m["years"] == scope and m["storm_set"] == sid]
                    if not cells:
                        continue
                    n_set = cells[0]["n_storms_in_set"]
                    txt = ""
                    for ml in mlabels:
                        c = [x for x in cells if x["measure"] == ml][0]
                        if what == "median":
                            s = fnum(c["median"], 2)
                        else:
                            s = "-" if c["min"] is None else f"{c['min']:.2f}..{c['max']:.2f}"
                        if c["n"] < n_set:
                            s += f"*({c['n']})"
                        txt += f"{s:>{wd}}"
                    log("    {:<11}({}){:>6}{}".format(scope, sid, n_set, txt))
    log("\n  storm sets: " + "; ".join(f"({s}) {d}" for s, d in STORM_SETS))


# ------------------------------------------------------------------------------------------------ compare with an old pass 2 file
def compare_old(rows, old_path: Path):
    """Line-by-line comparison of the new pass 2 rows with an old pass2_timing_long.csv (any version), keyed on water_year, run_key, storm.
    Returns (lines, ok)."""
    L = [f"COMPARE WITH OLD FILE: {old_path}"]
    if not old_path.is_file():
        return L + ["  The old file was not found, nothing compared."], False
    old = pd.read_csv(old_path, comment="#", dtype=str, keep_default_na=False)
    buf = io.StringIO()
    pd.DataFrame(rows, columns=P2_COLS, dtype=object).to_csv(buf, index=False)
    buf.seek(0)
    new = pd.read_csv(buf, dtype=str, keep_default_na=False)
    keys = ["water_year", "run_key", "storm"]
    if any(k not in old.columns for k in keys):
        return L + [f"  The old file has no columns {keys}; nothing compared."], False
    skip = {"run_label"}                         # screens are labeled 'Screen n (Sx)' from v1.4 on
    shared = [c for c in new.columns if c in old.columns and c not in skip]
    new_only = [c for c in new.columns if c not in old.columns]
    old_only = [c for c in old.columns if c not in new.columns]
    o = old.set_index(keys, drop=False)
    n = new.set_index(keys, drop=False)
    common = [k for k in n.index if k in o.index]
    only_new = sorted({k[1] for k in n.index if k not in o.index})
    only_old = sorted({k[1] for k in o.index if k not in n.index})
    diffs = []
    for k in common:
        for c in shared:
            if str(o.loc[k, c]) != str(n.loc[k, c]):
                diffs.append((k, c, o.loc[k, c], n.loc[k, c]))
    L.append(f"  rows: old {len(old)}, new {len(new)}, in both {len(common)}; columns compared {len(shared)}")
    if new_only:
        L.append("  columns only in the new file (not compared): " + ", ".join(new_only))
    if old_only:
        L.append("  columns only in the old file (not compared): " + ", ".join(old_only))
    if only_new:
        L.append("  runs only in the new file: " + ", ".join(only_new))
    if only_old:
        L.append("  runs only in the old file: " + ", ".join(only_old) + "  (add --include-kchan, or --wy2017-runs oct1, to reproduce them)")
    L.append(f"  cells that differ: {len(diffs)}")
    for k, c, a, b in diffs[:20]:
        L.append(f"    {k[0]} {k[1]} {k[2]}  {c}: old {a!r}  new {b!r}")
    if len(diffs) > 20:
        L.append(f"    ... and {len(diffs) - 20} more")
    ok = (not diffs) and not only_new and not only_old and len(common) > 0
    L.append("  RESULT: " + ("IDENTICAL on every compared cell, and the same rows." if ok else "NOT IDENTICAL (see above)."))
    return L, ok


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
    ap.add_argument("--pass", dest="pass_", type=int, choices=[1, 2], help="1 = gauge only, 2 = every run against the gauge (--dry-run alone means 2)")
    ap.add_argument("--floor-mm3", type=float, default=None, help="OPTIONAL (off by default): gauge volume (M m3) below which a storm's timing is blanked")
    ap.add_argument("--model-floor-mm3", type=float, default=None, help="optional, with --floor-mm3: also blank the timing of any run whose own volume (M m3) is below this")
    ap.add_argument("--root", default=str(HERE), help="the lab folder (default: the folder this script is in)")
    ap.add_argument("--results", default=None, help="folder with the outlet files (default <root>/results)")
    ap.add_argument("--obs2017", default=None)
    ap.add_argument("--obs2020", default=None)
    ap.add_argument("--out-dir", default=None, help="output folder (default <root>/timing_out_step0; timing_out_check_<mode> for the old window rules)")
    ap.add_argument("--windows", choices=["step0", "as-scored", "cut", "uncut"], default="step0")
    ap.add_argument("--windows-file", default=None, help="the Step 0 windows table (default <root>/step0_windows.csv)")
    ap.add_argument("--wy2017-runs", choices=["sep1", "oct1"], default="sep1", help="sep1 (default): the 1 September 2016 runs; oct1: the four old 1 October whole-year runs")
    ap.add_argument("--include-kchan", action="store_true", help="also score the four Kchan runs (left out by default)")
    ap.add_argument("--ratio-source", default=None, help="window_scores CSV for the Step 0 ratio check (default <root>/window_scores_R4.csv; typed values if absent)")
    ap.add_argument("--compare-old", default=None, metavar="OLD_CSV", help="pass 2: compare the new rows with an old pass2_timing_long.csv, line by line")
    ap.add_argument("--stem", action="append", default=[], metavar="YEAR:KEY=NAME", help="other outlet-file name for one run, e.g. 2017:original=NAME")
    ap.add_argument("--print-all", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--ignore-check-failure", action="store_true", help="diagnosis only")
    ap.add_argument("--definitions", action="store_true", help="print every definition and exit")
    ap.add_argument("--show-arrival", nargs="?", const="", default=None, metavar="PATH", help="print the first-arrival rule from score_run_sycamore.py and exit")
    return ap


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    root = Path(a.root)
    mode = a.windows
    if a.definitions:
        log("\n".join(definition_lines(mode, a.floor_mm3, a.model_floor_mm3)))
        return 0
    if a.show_arrival is not None:
        return show_arrival(Path(a.show_arrival) if a.show_arrival else root / "score_run_sycamore.py")
    if a.dry_run and a.pass_ is None:
        a.pass_ = 2
    if a.pass_ is None:
        die("choose --pass 1 or --pass 2 (or --dry-run, --show-arrival, --definitions). See the top of this file.")
    if a.floor_mm3 is not None and not (math.isfinite(a.floor_mm3) and a.floor_mm3 >= 0):
        die("--floor-mm3 must be a number of 0 or more.")
    if a.model_floor_mm3 is not None:
        if not (math.isfinite(a.model_floor_mm3) and a.model_floor_mm3 >= 0):
            die("--model-floor-mm3 must be a number of 0 or more.")
        if a.floor_mm3 is None:
            die("--model-floor-mm3 is an extra on top of --floor-mm3 X; give both.")
    if a.compare_old and a.pass_ != 2:
        die("--compare-old compares pass 2 rows; use --pass 2.")

    results = Path(a.results) if a.results else root / "results"
    out_dir = Path(a.out_dir) if a.out_dir else root / ("timing_out_step0" if mode == "step0" else "timing_out_check_" + mode.replace("-", "_"))
    if out_dir.resolve() == results.resolve():
        die("--out-dir is the results folder. This script never writes there.")
    runs = build_runs(a.wy2017_runs, a.include_kchan, a.stem)
    if not results.is_dir():
        die(f"cannot find the results folder: {results}\n  Run this from workspaces/Sycamore/lab/, or point to it with --results.")

    log(f"\nscore_timing_sycamore v{VERSION}   pass {a.pass_}   windows: {mode}   WY2017 runs: {a.wy2017_runs}   Kchan: {'included' if a.include_kchan else 'left out'}"
        + ("   [DRY RUN]" if a.dry_run else ""))
    log(f"  floor: " + ("none (every storm is scored)" if a.floor_mm3 is None else f"{a.floor_mm3:g} M m3 (asked for)") + f"   output folder: {out_dir}")
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
    # windows
    wfile = Path(a.windows_file) if a.windows_file else root / "step0_windows.csv"
    wins = make_windows(mode, wfile)
    log(f"\nWINDOW TABLE ({WINDOW_SET_NAMES[mode]}; Arizona local time, [start, end), hours named by their start)")
    for ln in window_table_lines(wins):
        log(ln)
    if mode == "step0":
        log(f"  window lengths: all 12 match the required list ({', '.join(str(STEP0_HOURS[k]) for k in STEP0_HOURS)}).")
    # outlet files
    models, infos, missing = load_models(runs, results)
    log("\nRUN INVENTORY (Time_hr is read as hours after the run's own start; file dates are this machine's clock, a Codespace shows UTC)")
    for ln in inventory_lines(runs, infos):
        log(ln)
    bad_len = [(r, infos[(r.year, r.key)]) for r in runs if (r.year, r.key) in infos and not infos[(r.year, r.key)]["ok"]]
    if missing or a.dry_run:
        w17 = [r for r in runs if r.year == "WY2017"]
        if any(r in missing for r in w17) or a.dry_run:
            log("\nWY2017 outlet files in results/ (candidates; compare with the run list above and confirm the screen names):")
            cl = candidate_lines(results, runs)
            for ln in cl:
                log(ln)
            if not cl:
                log("    (none: no outlet file in results/ has WY2017 or W5_ in its name)")
    if missing:
        log("\nOUTLET FILES NOT FOUND (these runs are skipped):")
        for r in missing:
            log(f"  {r.year} {r.label:<22} {r.fname}   (key {r.key})")
    if bad_len and not (a.ignore_check_failure and not a.dry_run):
        msg = "a run file's length does not match its expected hours, so reading it with its start would shift every time:\n  " + "\n  ".join(
            f"{r.year} {r.label} ({r.fname}): {i['why']}" for r, i in bad_len)
        if a.dry_run:
            log("\nLENGTH CHECK FAILED (a real run would stop here): " + msg)
        else:
            die(msg + "\n  Use --stem YEAR:KEY=NAME to point the run at the right file.")
    elif bad_len:
        log("\n  --ignore-check-failure was given: carrying on with files whose length does not match.")
    # priority runs
    for yr in ("WY2017", "WY2020"):
        have = [k for k in PRIORITY if k in models[yr]]
        need = len(PRIORITY) if mode == "step0" else 1
        if len(have) < need:
            msg = (f"{'both Re-tuned and Original are' if need == 2 else 'neither Re-tuned nor Original was found; they are'} needed for {yr}"
                   f" (the Step 0 ratio check covers both). See the inventory above; check the names or use --stem {yr[2:]}:retuned=NAME.")
            if a.dry_run:
                log("\nWOULD STOP: " + msg)
            else:
                die(msg)
    # ratio source
    ratio_source, ratio_text = None, ""
    if mode == "step0":
        rs = Path(a.ratio_source) if a.ratio_source else root / "window_scores_R4.csv"
        if a.ratio_source and not rs.is_file():
            die(f"cannot find the ratio source {rs}")
        if rs.is_file():
            ratio_source = read_ratio_source(rs)
            ratio_text = f"{rs.name} (vol_model_Mm3 / vol_gauge_Mm3, matched hours, from score_windows_sycamore.py)"
            need_keys = [(k, y, n) for k in PRIORITY for (y, n) in RATIO_STORMS]
            nf = [k for k in need_keys if k not in ratio_source]
            if nf:
                die(f"the ratio source {rs.name} has no storm row for: " + "; ".join(f"{k[0]} {k[1]} {k[2]}" for k in nf[:6])
                    + (" ..." if len(nf) > 6 else ""))
        else:
            ratio_text = "the typed per-storm ratios (RATIO_TYPED; window_scores_R4.csv was not found next to this script)"
    # 1. volume check
    can_check = all(k in models[y] for y in ("WY2017", "WY2020") for k in (PRIORITY if mode == "step0" else ()))
    if mode == "step0":
        log(f"\nVOLUME CHECK (Step 0): vol_ratio_gauge_hours of Re-tuned and Original, 10 storms with a gauge response, tolerance {RATIO_TOL:g}\n  record: {ratio_text}")
        if can_check:
            status, bad, lines = ratio_check(wins, gauges, models, ratio_source)
        else:
            status, bad, lines = {}, ["Re-tuned or Original is missing for a year, so the check cannot be made"], ["    (cannot be made)"]
    else:
        log("\nVOLUME CHECK (old window rules; recomputed volumes against the volumes on record; M m3; ratio = model / gauge)")
        skip = frozenset((r.year, r.key) for r in runs if r.year == "WY2017" and r.start != OCT1_2016)
        status, bad, lines = volume_check(wins, gauges, models, skip_models=skip)
    for ln in lines:
        log(ln)
    if bad:
        log(f"\nVOLUME CHECK: MISMATCH ({len(bad)} value(s)):")
        for b in bad:
            log("  " + b)
        if a.dry_run:
            log("  (a real run would stop here; nothing is written in a dry run)")
        elif not a.ignore_check_failure:
            die("the volumes do not reproduce, so nothing was written. Do not use any timing number until this is explained.")
        else:
            log("  --ignore-check-failure was given: carrying on. Rows are labeled MISMATCH.")
    else:
        n_m = sum(1 for v in status.values() if v == "MATCH")
        log(f"\nVOLUME CHECK: MATCH. {n_m} run-and-storm values reproduce the record" + (f" within {RATIO_TOL:g}." if mode == "step0" else f" within {TOL_M3} M m3; the rest have no reference (NOT CHECKED)."))

    if a.dry_run:
        n_found = sum(1 for r in runs if (r.year, r.key) in infos)
        log(f"\nDRY RUN: {len(wins)} windows, {n_found} of {len(runs)} run files found"
            + (", every length matches" if not bad_len else ", LENGTH PROBLEMS ABOVE") + ". No scores were computed and nothing was written.")
        return 0 if not (bad or bad_len or missing) else 1

    # headers
    bundles = {w.key: gauge_bundle(w, gauges[w.year]) for w in wins}
    header = definition_lines(mode, a.floor_mm3 if a.pass_ == 2 else None, a.model_floor_mm3 if a.pass_ == 2 else None, ratio_text)
    header.append("VOLUME CHECK RESULT: " + ("MISMATCH, carried on with --ignore-check-failure" if bad else "MATCH") + ". Gauge files: "
                  + "; ".join(f"{y} {p.name}" for y, p in gpaths.items()) + ".")
    header.append("WINDOW TABLE: " + "; ".join(f"{w.year} {w.name} {w.lo:%Y-%m-%d %H:%M} to {w.hi:%Y-%m-%d %H:%M} ({w.hours} h, gauge response {w.gauge_response})" for w in wins) + ".")
    for r in runs:
        i = infos.get((r.year, r.key))
        if i:
            header.append(f"RUN FILE: {r.year} {r.key} '{r.label}' {r.fname}; modified {i['modified']} (machine clock); Time_hr 0 = {r.start:%Y-%m-%d %H:%M}; first Time_hr {i['first']:g}, last Time_hr {i['last']:g}; expected {r.hours} h.")
    if missing:
        header.append("RUNS NOT FOUND (skipped): " + "; ".join(f"{r.year} {r.label}" for r in missing))

    for fn in (("pass1_gauge_table.csv",) if a.pass_ == 1 else ("pass2_timing_long.csv", "pass2_medians.csv")) + (("step0_windows.csv",) if mode == "step0" else ()):
        guard_overwrite(out_dir / fn)
    # 2. pass 1
    if a.pass_ == 1:
        rows = [pass1_row(w, bundles[w.key], status) for w in wins]
        print_pass1(rows)
        path = out_dir / "pass1_gauge_table.csv"
        write_csv(path, header, pd.DataFrame(rows, columns=P1_COLS, dtype=object))
        log(f"\nWrote {len(rows)} rows to {path}")
        save_windows_copy(out_dir, wins, header, mode, wfile)
        log("\nPASS 1 DONE.")
        return 0

    # 3. pass 2
    rows = []
    for r in runs:
        if r.key not in models[r.year]:
            continue
        for w in [x for x in wins if x.year == r.year]:
            rows.append(pass2_row(w, bundles[w.key], r.key, r.label, r.stem, models[r.year][r.key], a.floor_mm3, status, a.model_floor_mm3))
    print_pass2(rows, a.print_all)
    med = build_medians(rows)
    print_medians(med)
    c5 = check5(rows, mode)
    log("")
    for ln in c5:
        log(ln)
    n_blank = sum(1 for r in rows if r["blank_reasons"])
    log(f"\n{len(rows)} run-and-storm rows; {n_blank} have at least one blank value (reason in blank_reasons).")
    header.extend(c5)
    path = out_dir / "pass2_timing_long.csv"
    write_csv(path, header, pd.DataFrame(rows, columns=P2_COLS, dtype=object))
    log(f"Wrote {len(rows)} rows to {path}")
    mh_ = list(header)
    mh_.insert(2, "THIS FILE (pass2_medians.csv): one row per run, years, storm set and measure; median, min, max over the storms where the measure is defined; n = how many storms had it defined; n_storms_in_set = storms in the set. Measures: peak ratio, peak time (model minus gauge, h), t10, t50, t90 (model minus gauge, h), lag (h), tail volume ratio, Pearson r (hourly). Storm sets: " + "; ".join(f"({s}) {d}" for s, d in STORM_SETS) + ". 'both years' pools WY2017 and WY2020 storms of the same run key; years_scored says which years that run was scored for.")
    mpath = out_dir / "pass2_medians.csv"
    write_csv(mpath, mh_, pd.DataFrame(med, columns=P3_COLS, dtype=object))
    log(f"Wrote {len(med)} rows to {mpath}")
    save_windows_copy(out_dir, wins, header, mode, wfile)
    if a.compare_old:
        cl, ok = compare_old(rows, Path(a.compare_old))
        log("")
        for ln in cl:
            log(ln)
    return 0


def save_windows_copy(out_dir: Path, wins, header, mode: str, src: Path):
    """step0_windows.csv in the output folder: the table the run used (with hours), headed with the window set."""
    if mode != "step0":
        return
    dst = out_dir / "step0_windows.csv"
    if src.exists() and dst.exists() and src.resolve() == dst.resolve():
        return
    df = pd.DataFrame([{"water_year": w.year, "storm": w.name, "start": f"{w.lo:%Y-%m-%d %H:%M}", "end": f"{w.hi:%Y-%m-%d %H:%M}",
                        "hours": w.hours, "gauge_response": w.gauge_response} for w in wins], dtype=object)
    hl = [header[0], header[1], f"THIS FILE (step0_windows.csv): the table the run read from {src.name}, with the hours. Required hours: " + ", ".join(f"{k[1]} {k[0][2:]} {v}" for k, v in STEP0_HOURS.items()) + "."]
    write_csv(dst, hl, df)
    log(f"Wrote {len(df)} rows to {dst}")


if __name__ == "__main__":
    sys.exit(main())
