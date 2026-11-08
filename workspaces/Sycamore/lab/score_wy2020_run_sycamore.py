#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
score_wy2020_run_sycamore.py   (v1.0, 2026-11-07)
=================================================

Scores ONE WY2020 model run (default: S1, the Ks x1.0 screen) on the seven scored storms, against the USGS gauge and against the
Re-tuned run. READ-ONLY: it opens result files and writes two small CSV files. It does not run a model, change a file in results/,
or fit anything.

THE READ-OUT (agreed with Tracy on 2026-11-07 at 12:59, BEFORE any S1 number was opened)
  R1  Nov 19-22   S1 volume at most  1.40 M m3   (Re-tuned 2.800 cut by half or more)
  R2  Dec 23-28   S1 volume at least 0.44 M m3   (Re-tuned 0.218 at least doubled)
  R3  Mar 10-13   S1 volume at least 1.75 M m3   (Re-tuned 0.875 at least doubled)
  The script prints MET or NOT MET for each. It does not combine them into a verdict. Those three are for S1 only: if you score another
  run, the same three lines are printed as description, and the script says so.
  Also printed, with no pass or fail (they only describe):
    - S1 / Re-tuned volume for each storm, and how many storms went up and how many went down ("lockstep": did every storm move the same way?)
    - Spearman rank correlation against the gauge (Re-tuned: -0.286) and mean absolute error (Re-tuned: 2.50 M m3)
    - the seven-storm totals (gauge 19.02, Re-tuned 14.76 M m3)
    - the runoff the model MADE on the land, by type, for Re-tuned and for the run, and the share that is saturation excess (Satsrf)

WHAT IT CAN AND CANNOT SAY
  Seven storms is a very small sample, and S1 is a one-change screen (only Ks differs from Re-tuned), not a calibrated setting.
  Nothing here is re-fit on WY2020. A storm-volume match is not proof that the model has the right runoff process, and a miss on one storm
  can come from the rain data (the radar) as well as from the model.

THE SEVEN STORMS AND THE WINDOW (same rule as the earlier scoring)
  Nov 19-22, Nov 28-30, Dec 23-28, Feb 10-11, Feb 22-23, Mar 10-13, Mar 18-19.
  Window = 00:00 on the first rain day to 00:00 four days after the last rain day.
  Volume = hourly mean flow x 3,600 s, added up over the window, in M m3 (million cubic meters).
  The gauge and Re-tuned volumes are checked against the numbers in the WY2020 write-up first (MATCH or MISMATCH). If they do not match,
  a wrong file or window is being used: do not use the output.

WHAT IT READS (read-only)
  results/Sycamore_WY2020_Retuned_Outlet.qout          the reference run (Time_hr 0 = Sep 1 2019 00:00)
  results/Sycamore_WY2020_<run>_Outlet.qout            the run being scored
  ../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv                columns datetime, q_cms
  results/Sycamore_WY2020_Retuned*.rft and ..._<run>*.rft     runoff made, by type (columns Time_hr, Hsrf_m3_s, Sbsrf_m3_s, Psrf_m3_s, Satsrf_m3_s).
                                                       The file name has the run hours after the run name (for example ...4368_00.rft for WY2017),
                                                       so it is found by name pattern. If it is not found, or more than one file fits, that
                                                       section is skipped and the script says why (use --rft / --rft-ref to name a file).
WHAT IT WRITES (in score_out/ unless you set --out-dir)
  score_wy2020_<run>_per_storm.csv     one row per storm: window, gauge, Re-tuned and run volumes, ratios, and the runoff made by type
  score_wy2020_<run>_summary.csv       one row per number: the read-out lines, Spearman, mean absolute error, totals, runoff by type

COMMANDS (from workspaces/Sycamore/lab/)
  python score_wy2020_run_sycamore.py --dry-run      # finds and reads the files, lists the windows; prints no volumes, writes nothing
  python score_wy2020_run_sycamore.py                # scores S1 (run name S1_Ks1p0)
  python score_wy2020_run_sycamore.py --run S2_Ks0p5          # another run; the prefix Sycamore_WY2020_ may be left off
  python score_wy2020_run_sycamore.py --no-types              # skip the runoff-by-type section (volumes only)
  python score_wy2020_run_sycamore.py --rft results/<file>.rft --rft-ref results/<file>.rft     # name the .rft files yourself

THE LINE "NOTE: ... within 5% of Re-tuned on all seven storms"
  If the run is within 5% of Re-tuned on every storm, that is what a Ks change that never reached the model would look like. The script
  prints a reminder to check the run's .in file and the notebook's Ks printout before reading anything into the result. It is a setup
  check, not a result.

WORDS
  Spearman rank correlation   asks only whether the model puts the storms in the same ORDER of size as the gauge. +1 = same order,
                              -1 = reversed, 0 = no relation. With seven storms, treat values between about -0.5 and +0.5 as "no clear order".
  mean absolute error (MAE)   the average size of the miss, ignoring whether it was too high or too low (here in M m3 per storm).
  runoff made, by type        what the model generates on the land, all cells added, BEFORE it is routed to the outlet (so it is not outlet flow):
        Hsrf    infiltration excess  (rain falls faster than the soil takes it)
        Sbsrf   subsurface stormflow (water coming out of the soil along the slope)
        Psrf    perched-water runoff
        Satsrf  saturation excess    (the soil is already full, so more rain runs off)
  lockstep                    every storm moving the same way (all up or all down) when one setting changes.
  M m3                        million cubic meters.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
T0 = pd.Timestamp("2019-09-01 00:00")          # Time_hr 0 of every WY2020 run
RUN_HOURS = 6575
AFTER_DAYS = 3                                  # window ends at 00:00, (1 + 3) days after the last rain day's 00:00
MATCH_TOL = 0.005                               # M m3
MIN_GAUGE_VOL = 0.05                            # M m3; below this a ratio to the gauge is not shown
SAME_BAND = 0.05                                # all storms within +-5% of Re-tuned: print the setup reminder
WRITEUP_SPEARMAN, WRITEUP_MAE = -0.286, 2.50    # Re-tuned, from the WY2020 write-up
PREFIX = "Sycamore_WY2020_"
REF_KEY, REF_STEM, REF_LABEL = "retuned", PREFIX + "Retuned", "Re-tuned"
DEFAULT_RUN = "S1_Ks1p0"
LABELS = {"S1_Ks1p0": "S1", "S2_Ks0p5": "S2", "S3_Bed0p5": "S3", "S4_GW0p5": "S4"}
READOUT_AGREED = "2026-11-07 12:59"
TYPES = ["Hsrf", "Sbsrf", "Psrf", "Satsrf"]

# (storm, first rain day, last rain day)
STORMS = [
    ("Nov 19-22", "2019-11-19", "2019-11-22"), ("Nov 28-30", "2019-11-28", "2019-11-30"), ("Dec 23-28", "2019-12-23", "2019-12-28"),
    ("Feb 10-11", "2020-02-10", "2020-02-11"), ("Feb 22-23", "2020-02-22", "2020-02-23"),
    ("Mar 10-13", "2020-03-10", "2020-03-13"), ("Mar 18-19", "2020-03-18", "2020-03-19"),
]
# WY2020 write-up volumes (M m3): gauge and Re-tuned, scored storm windows
EXPECTED = {
    "Nov 19-22": (0.015, 2.800), "Nov 28-30": (1.473, 3.235), "Dec 23-28": (4.069, 0.218), "Feb 10-11": (0.386, 1.767),
    "Feb 22-23": (4.131, 4.839), "Mar 10-13": (5.626, 0.875), "Mar 18-19": (3.317, 1.031),
}
# the agreed read-out: (id, storm, "max" or "min", limit, what it means)
READOUT = [
    ("R1", "Nov 19-22", "max", 1.40, "Re-tuned 2.800 cut by half or more"),
    ("R2", "Dec 23-28", "min", 0.44, "Re-tuned 0.218 at least doubled"),
    ("R3", "Mar 10-13", "min", 1.75, "Re-tuned 0.875 at least doubled"),
]


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    sys.exit("\nSTOPPED: " + msg + "\n")


def need(path: Path, what: str) -> Path:
    if not Path(path).is_file():
        die("cannot find the {} file:\n    {}\n  Run this from workspaces/Sycamore/lab/, or point to it with the matching option "
            "(see the top of this script).".format(what, path))
    return Path(path)


def tag_of(key: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", key).strip("_")


def fmt(v, nd: int = 2) -> str:
    return "n/a" if v is None or not np.isfinite(v) else "{:.{n}f}".format(v, n=nd)


def fmt_signed(v, nd: int = 2) -> str:
    return "n/a" if v is None or not np.isfinite(v) else "{:+.{n}f}".format(v, n=nd)


def fmt_ratio(v) -> str:
    if v is None or not np.isfinite(v):
        return "n/a"
    if abs(v) < 10:
        return "{:.2f}".format(v)
    if abs(v) < 100:
        return "{:.1f}".format(v)
    return "{:.0f}".format(v)


# ------------------------------------------------------------------------------------------------ reading
def read_qout(path: Path) -> pd.Series:
    raw = pd.read_csv(path).iloc[:, :2]
    raw.columns = ["Time_hr", "Q"]
    raw = raw.apply(pd.to_numeric, errors="coerce").dropna()
    idx = T0 + pd.to_timedelta(raw["Time_hr"].to_numpy(float), unit="h")
    s = pd.Series(raw["Q"].to_numpy(float), index=idx)
    return s[~s.index.duplicated(keep="last")].sort_index()


def read_obs(path: Path) -> pd.Series:
    o = pd.read_csv(path)
    for c in ("datetime", "q_cms"):
        if c not in o.columns:
            die("{} has no column '{}'. Columns found: {}".format(path, c, list(o.columns)[:10]))
    o["datetime"] = pd.to_datetime(o["datetime"])
    s = pd.Series(pd.to_numeric(o["q_cms"], errors="coerce").to_numpy(float), index=o["datetime"])
    return s[~s.index.duplicated(keep="last")].sort_index()


def find_rft(results: Path, stem: str):
    """The .rft file(s) of a run: the run name, then the run hours and a suffix (for example Sycamore_WY2020_Retuned6575_00.rft)."""
    pat = re.compile("^" + re.escape(stem) + r"\d.*\.rft$")
    return sorted(p for p in results.glob(stem + "*.rft") if pat.match(p.name))


def read_rft(path: Path):
    """Hourly runoff made by type (m3/s, all cells added). Index = hour END time, as the files label it. Returns (frame, note)."""
    d = pd.read_csv(path)
    want = ["Time_hr"] + ["{}_m3_s".format(t) for t in TYPES]
    miss = [c for c in want if c not in d.columns]
    if miss:
        die("{} has no column(s) {}. Columns found: {}".format(path, miss, list(d.columns)))
    d = d[want].apply(pd.to_numeric, errors="coerce").dropna()
    big = (d[want[1:]] >= 9000).any(axis=1)            # an unfilled last row is flagged 9999 in the .mrf file; the same guard here
    note = "{} row(s) with a value of 9000 or more were dropped (unfilled)".format(int(big.sum())) if big.any() else ""
    d = d[~big]
    d.index = T0 + pd.to_timedelta(d["Time_hr"].to_numpy(float), unit="h")
    return d[~d.index.duplicated(keep="last")].sort_index(), note


# ------------------------------------------------------------------------------------------------ numbers
def volume(q: pd.Series, lo, hi):
    """Million m3 from the hourly mean flow x 3,600 s over [lo, hi), counting only hours with data. Returns (volume, hours)."""
    s = q[(q.index >= lo) & (q.index < hi)]
    h = s.resample("60min").mean().dropna()
    return (float(h.sum() * 3600.0 / 1e6) if len(h) else float("nan")), int(len(h))


def spearman(x, y) -> float:
    """Spearman rank correlation (ties share the average rank), storms with a missing value left out."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    rx, ry = pd.Series(x[ok]).rank().to_numpy(), pd.Series(y[ok]).rank().to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def mae(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.abs(x[ok] - y[ok]).mean()) if ok.any() else float("nan")


def default_label(key: str) -> str:
    return LABELS.get(key, key[:8])


def parse_run(spec: str):
    nm = spec.strip()
    if nm.lower() in ("retuned", "re-tuned"):
        return "retuned", REF_STEM
    stem = nm if nm.startswith("Sycamore_") else PREFIX + nm
    if stem.endswith("_Outlet.qout"):
        stem = stem[:-len("_Outlet.qout")]
    key = stem[len(PREFIX):] if stem.startswith(PREFIX) else stem
    return key, stem


def need_run(results: Path, stem: str, what: str) -> Path:
    p = results / (stem + "_Outlet.qout")
    if p.is_file():
        return p
    have = sorted(f.name[:-len("_Outlet.qout")] for f in results.glob("*_Outlet.qout")) if results.is_dir() else []
    msg = "cannot find the outlet file for {}:\n    {}".format(what, p)
    if have:
        msg += "\n  Outlet files found in {}:\n    ".format(results) + "\n    ".join(have[:25])
        msg += "\n  Give the run name without _Outlet.qout (the prefix {} may be left off).".format(PREFIX)
    else:
        msg += "\n  Run this from workspaces/Sycamore/lab/, or point to the folder with --results."
    die(msg)


# ------------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Score one WY2020 run on the seven scored storms against the gauge and Re-tuned (see the top of this file).")
    ap.add_argument("--root", default=str(HERE), help="the lab folder (default: the folder this script is in)")
    ap.add_argument("--results", default=None, help="folder with the result files (default <root>/results)")
    ap.add_argument("--obs", default=None, help="USGS WY2020 gauge file (default <root>/../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv)")
    ap.add_argument("--out-dir", default=None, help="where the two CSV files go (default <root>/score_out)")
    ap.add_argument("--run", default=DEFAULT_RUN, help="the run to score, with or without the prefix {} (default {})".format(PREFIX, DEFAULT_RUN))
    ap.add_argument("--label", default=None, help="a short name for the run in the tables (default: S1 for S1_Ks1p0)")
    ap.add_argument("--rft", default=None, help="the run's .rft file, if its name is not found automatically")
    ap.add_argument("--rft-ref", default=None, help="the Re-tuned run's .rft file, if its name is not found automatically")
    ap.add_argument("--no-types", action="store_true", help="skip the runoff-by-type section")
    ap.add_argument("--dry-run", action="store_true", help="find and read the files and list the windows; print no volumes, write nothing")
    a = ap.parse_args(argv)

    root = Path(a.root)
    results = Path(a.results) if a.results else root / "results"
    obs_path = Path(a.obs) if a.obs else root.parent / "init_data" / "obs" / "sycamore_09510200_obs_flow_WY2020.csv"
    out_dir = Path(a.out_dir) if a.out_dir else root / "score_out"
    key, stem = parse_run(a.run)
    if key == REF_KEY:
        die("--run is the Re-tuned run itself. Give the run you want to compare with it, for example --run S1_Ks1p0.")
    label = (a.label or default_label(key)).strip()[:8] or key[:8]
    if label.lower().replace("-", "") == "retuned":
        die("the run's short name (--label) cannot be 'Re-tuned', because that name is used for the reference run. Pick another, for example --label S1.")
    is_s1 = key == DEFAULT_RUN

    log("\nscore_wy2020_run_sycamore.py v1.0: reading files (read-only)")
    gauge = read_obs(need(obs_path, "WY2020 gauge"))
    ref_q = read_qout(need_run(results, REF_STEM, "the Re-tuned run"))
    run_q = read_qout(need_run(results, stem, "the run " + key))
    log("  gauge: {} readings, {:%b %d %Y} to {:%b %d %Y}".format(gauge.notna().sum(), gauge.index.min(), gauge.index.max()))
    log("  Re-tuned outlet: {} rows".format(len(ref_q)))
    log("  {} outlet [{}]: {} rows".format(label, stem, len(run_q)))
    for nm, q in (("Re-tuned", ref_q), (label, run_q)):
        last_h = (q.index.max() - T0) / pd.Timedelta(hours=1)
        if last_h < RUN_HOURS - 0.01:
            log("  WARNING: the {} outlet file ends at hour {:.1f}, not {} (the run did not finish, or the file is cut). "
                "Storms after that hour will show n/a.".format(nm, last_h, RUN_HOURS))

    # runoff by type
    rft = {}
    if a.no_types:
        log("  .rft files: skipped (--no-types)")
    else:
        for nm, spec, st in (("Re-tuned", a.rft_ref, REF_STEM), (label, a.rft, stem)):
            if spec:
                found = [need(Path(spec), nm + " .rft")]
            else:
                found = find_rft(results, st)
            if len(found) == 1:
                df, note = read_rft(found[0])
                rft[nm] = df
                log("  {} .rft: {} ({} rows){}".format(nm, found[0].name, len(df), "; " + note if note else ""))
            elif not found:
                log("  {} .rft: NOT FOUND (looked for {}*.rft with the run hours after the name in {}). The runoff-by-type section will be skipped; "
                    "use --rft / --rft-ref to name the file.".format(nm, st, results))
            else:
                log("  {} .rft: MORE THAN ONE file fits: {}. The runoff-by-type section will be skipped; use --rft / --rft-ref to name the one you want.".format(
                    nm, ", ".join(p.name for p in found)))
    types_ok = (not a.no_types) and ("Re-tuned" in rft) and (label in rft)

    windows = []
    for name, f, l in STORMS:
        lo = pd.Timestamp(f)
        hi = pd.Timestamp(l) + pd.Timedelta(days=1 + AFTER_DAYS)
        windows.append((name, lo, hi))
    log("\nWindows (seven scored storms):")
    for name, lo, hi in windows:
        log("  {:10s} {:%Y-%m-%d %H:%M} to {:%Y-%m-%d %H:%M}".format(name, lo, hi))
    if a.dry_run:
        log("\nDRY RUN: all inputs were found and read. No volumes were printed and nothing was written.")
        return 0

    # ---- volumes
    rows = []
    for name, lo, hi in windows:
        gv, gh = volume(gauge, lo, hi)
        rv, _ = volume(ref_q, lo, hi)
        sv, sh = volume(run_q, lo, hi)
        n_hours = int((hi - lo) / pd.Timedelta(hours=1))
        r = {"storm": name, "lo": lo, "hi": hi, "gauge": gv, "gauge_hours": gh, "ref": rv, "run": sv, "run_hours": sh, "n_hours": n_hours}
        r["run_minus_gauge"] = sv - gv if np.isfinite(sv) and np.isfinite(gv) else float("nan")
        r["run_over_gauge"] = sv / gv if (np.isfinite(sv) and np.isfinite(gv) and gv >= MIN_GAUGE_VOL) else float("nan")
        r["run_over_ref"] = sv / rv if (np.isfinite(sv) and np.isfinite(rv) and rv > 0) else float("nan")
        if np.isfinite(sv) and np.isfinite(rv):
            r["move"] = "up" if sv - rv > 0.0005 else ("down" if rv - sv > 0.0005 else "same")
        else:
            r["move"] = "n/a"
        if types_ok:
            for nm, tag in (("Re-tuned", "ref"), (label, "run")):
                for t in TYPES:
                    r["{}_{}".format(tag, t)] = volume(rft[nm]["{}_m3_s".format(t)], lo, hi)[0]
        rows.append(r)
    T = pd.DataFrame(rows)

    # ---- 0. check against the write-up
    bad = []
    for r in rows:
        eg, er = EXPECTED[r["storm"]]
        if not (np.isfinite(r["gauge"]) and abs(r["gauge"] - eg) <= MATCH_TOL):
            bad.append("{} gauge: here {} vs write-up {:.3f}".format(r["storm"], fmt(r["gauge"], 3), eg))
        if not (np.isfinite(r["ref"]) and abs(r["ref"] - er) <= MATCH_TOL):
            bad.append("{} Re-tuned: here {} vs write-up {:.3f}".format(r["storm"], fmt(r["ref"], 3), er))
        if r["gauge_hours"] < 0.97 * r["n_hours"]:
            log("  note: the gauge has data for {} of {} hours in the {} window; its volume counts only those hours.".format(
                r["gauge_hours"], r["n_hours"], r["storm"]))
    log("\n0. Check of the gauge and Re-tuned volumes against the WY2020 write-up")
    if bad:
        log("   MISMATCH ({} value(s) differ by more than {:g} M m3):".format(len(bad), MATCH_TOL))
        for x in bad:
            log("     " + x)
        log("   Wrong file, window or run. Do not use the numbers below until this is explained.")
    else:
        log("   MATCH. Gauge and Re-tuned volumes equal the write-up values for all seven storms (within {:g} M m3).".format(MATCH_TOL))

    # ---- 1. volumes
    log("\n1. Volumes over the scored window (M m3); run = {} ({})".format(label, key))
    L = label
    log("   {:10s} {:>7s} {:>9s} {:>8s} {:>10s} {:>9s} {:>9s}  {}".format(
        "storm", "gauge", "Re-tuned", L, L + "-gauge", L + "/gauge", L + "/Re-t", "vs Re-tuned"))
    for r in rows:
        log("   {:10s} {:>7s} {:>9s} {:>8s} {:>10s} {:>9s} {:>9s}  {}".format(
            r["storm"], fmt(r["gauge"], 3), fmt(r["ref"], 3), fmt(r["run"], 3), fmt_signed(r["run_minus_gauge"], 2),
            fmt_ratio(r["run_over_gauge"]), fmt_ratio(r["run_over_ref"]), r["move"]))
    log("   ({}/gauge is n/a when the gauge volume is under {:g} M m3.)".format(L, MIN_GAUGE_VOL))

    # ---- 2. read-out
    log("\n2. The read-out agreed {} (before S1 was opened)".format(READOUT_AGREED))
    if not is_s1:
        log("   NOTE: this read-out was agreed for S1 only. For {} the three lines below are description, not a test.".format(label))
    by_name = {r["storm"]: r for r in rows}
    met = []
    for rid, storm, kind, limit, why in READOUT:
        v = by_name[storm]["run"]
        if not np.isfinite(v):
            res = "CANNOT BE JUDGED (no {} volume for this storm)".format(label)
            ok = None
        elif kind == "max":
            ok = v <= limit
            res = "MET" if ok else "NOT MET"
        else:
            ok = v >= limit
            res = "MET" if ok else "NOT MET"
        met.append(ok)
        log("   {}  {:10s} {} {} at {} {:.2f} M m3 ({})  ->  {}".format(
            rid, storm, label, fmt(v, 3), "most" if kind == "max" else "least", limit, why, res))
    log("   {} of 3 met. (Each line stands on its own; this script gives no overall verdict.)".format(sum(1 for m in met if m is True)))

    # ---- 3. descriptive numbers
    log("\n3. Other numbers (description only, no pass or fail)")
    g, rf, rn = T["gauge"].to_numpy(float), T["ref"].to_numpy(float), T["run"].to_numpy(float)
    sp_ref, sp_run = spearman(rf, g), spearman(rn, g)
    mae_ref, mae_run = mae(rf, g), mae(rn, g)
    log("   Spearman against the gauge, seven storms:   Re-tuned {}   {} {}".format(fmt_signed(sp_ref, 3), L, fmt_signed(sp_run, 3)))
    log("   Mean absolute error per storm (M m3):       Re-tuned {}    {} {}".format(fmt(mae_ref, 2), L, fmt(mae_run, 2)))
    if np.isfinite(sp_ref) and np.isfinite(mae_ref):
        same = abs(sp_ref - WRITEUP_SPEARMAN) <= 0.001 and abs(mae_ref - WRITEUP_MAE) <= 0.01
        log("   (Re-tuned here {} the write-up's {:+.3f} and {:.2f}.)".format("equals" if same else "DOES NOT EQUAL", WRITEUP_SPEARMAN, WRITEUP_MAE))
    both = np.isfinite(g) & np.isfinite(rf) & np.isfinite(rn)
    if both.all():
        log("   Seven-storm totals (M m3):                  gauge {}   Re-tuned {}   {} {}".format(fmt(g.sum(), 2), fmt(rf.sum(), 2), L, fmt(rn.sum(), 2)))
    else:
        log("   Totals over the {} storms that have all three volumes (M m3):  gauge {}   Re-tuned {}   {} {}".format(
            int(both.sum()), fmt(g[both].sum(), 2), fmt(rf[both].sum(), 2), L, fmt(rn[both].sum(), 2)))
    ratios = T["run_over_ref"].to_numpy(float)
    n_up, n_down = int((T["move"] == "up").sum()), int((T["move"] == "down").sum())
    if np.isfinite(ratios).any():
        lock = ("all {} storms moved the same way (lockstep)".format(n_up + n_down) if (n_up == 0 or n_down == 0) and (n_up + n_down) > 1
                else "the storms moved in different directions")
        log("   {}/Re-tuned volume ratio: smallest {}, largest {}; {} storm(s) up, {} down: {}.".format(
            L, fmt_ratio(np.nanmin(ratios)), fmt_ratio(np.nanmax(ratios)), n_up, n_down, lock))
        if np.isfinite(ratios).all() and np.all(np.abs(ratios - 1.0) <= SAME_BAND):
            log("   NOTE: {} is within {:.0f}% of Re-tuned on all seven storms. That is what a Ks change that never reached the model would look "
                "like. Check the run's .in file and the notebook's Ks printout (step 8) before reading anything into this.".format(L, 100 * SAME_BAND))

    # ---- 4. runoff made by type
    type_rows = {}
    if types_ok:
        log("\n4. Runoff the model MADE on the land, by type (M m3; all cells added; before routing, so not outlet flow)")
        for nm, tag in (("Re-tuned", "ref"), (label, "run")):
            df = rft[nm]
            log("   {}".format(nm))
            log("   {:10s} {:>7s} {:>7s} {:>7s} {:>7s} {:>8s} {:>9s}".format("storm", "Hsrf", "Sbsrf", "Psrf", "Satsrf", "total", "Satsrf %"))
            sums = {t: 0.0 for t in TYPES}
            for r in rows:
                v = {t: r["{}_{}".format(tag, t)] for t in TYPES}
                for t in TYPES:
                    sums[t] += v[t] if np.isfinite(v[t]) else 0.0
                tot = sum(x for x in v.values() if np.isfinite(x))
                log("   {:10s} {:>7s} {:>7s} {:>7s} {:>7s} {:>8s} {:>9s}".format(
                    r["storm"], fmt(v["Hsrf"], 3), fmt(v["Sbsrf"], 3), fmt(v["Psrf"], 3), fmt(v["Satsrf"], 3), fmt(tot, 3),
                    "{:.1f}".format(100 * v["Satsrf"] / tot) if tot > 0 else "n/a"))
            tot7 = sum(sums.values())
            log("   {:10s} {:>7s} {:>7s} {:>7s} {:>7s} {:>8s} {:>9s}".format(
                "7 storms", fmt(sums["Hsrf"], 3), fmt(sums["Sbsrf"], 3), fmt(sums["Psrf"], 3), fmt(sums["Satsrf"], 3), fmt(tot7, 3),
                "{:.1f}".format(100 * sums["Satsrf"] / tot7) if tot7 > 0 else "n/a"))
            whole = {t: volume(df["{}_m3_s".format(t)], df.index.min(), df.index.max() + pd.Timedelta(hours=1))[0] for t in TYPES}
            totw = sum(whole.values())
            log("   {:10s} {:>7s} {:>7s} {:>7s} {:>7s} {:>8s} {:>9s}".format(
                "whole run", fmt(whole["Hsrf"], 3), fmt(whole["Sbsrf"], 3), fmt(whole["Psrf"], 3), fmt(whole["Satsrf"], 3), fmt(totw, 3),
                "{:.1f}".format(100 * whole["Satsrf"] / totw) if totw > 0 else "n/a"))
            type_rows[nm] = {"seven": dict(sums, total=tot7), "whole": dict(whole, total=totw)}
        a_, b_ = type_rows["Re-tuned"]["seven"], type_rows[label]["seven"]
        if a_["total"] > 0 and b_["total"] > 0:
            log("   Satsrf (saturation excess) share of the runoff made, seven storms together: Re-tuned {:.1f}%  ->  {} {:.1f}%".format(
                100 * a_["Satsrf"] / a_["total"], L, 100 * b_["Satsrf"] / b_["total"]))
            log("   Hsrf (infiltration excess) share, seven storms together:                    Re-tuned {:.1f}%  ->  {} {:.1f}%".format(
                100 * a_["Hsrf"] / a_["total"], L, 100 * b_["Hsrf"] / b_["total"]))
    elif not a.no_types:
        log("\n4. Runoff made, by type: skipped (the .rft file(s) were not found; see the messages near the top).")

    # ---- write the CSV files
    out_dir.mkdir(parents=True, exist_ok=True)
    tg = tag_of(key)
    per = pd.DataFrame({
        "storm": T["storm"], "window_start": T["lo"], "window_end": T["hi"], "gauge_Mm3": T["gauge"], "retuned_Mm3": T["ref"],
        "run": key, "run_Mm3": T["run"], "run_minus_gauge_Mm3": T["run_minus_gauge"], "run_over_gauge": T["run_over_gauge"],
        "run_over_retuned": T["run_over_ref"], "run_vs_retuned": T["move"]})
    if types_ok:
        for tag, nm in (("retuned", "ref"), ("run", "run")):
            for t in TYPES:
                per["{}_{}_made_Mm3".format(tag, t)] = T["{}_{}".format(nm, t)]
    per_path = out_dir / "score_wy2020_{}_per_storm.csv".format(tg)
    per.round(4).to_csv(per_path, index=False)
    summ = [("writeup_check", "MATCH" if not bad else "MISMATCH", "")]
    for (rid, storm, kind, limit, why), ok in zip(READOUT, met):
        summ.append((rid + "_" + storm.replace(" ", ""), "MET" if ok is True else ("NOT MET" if ok is False else "CANNOT BE JUDGED"),
                     "{} {}".format("at most" if kind == "max" else "at least", limit)))
    summ += [("spearman_retuned", fmt(sp_ref, 4), ""), ("spearman_run", fmt(sp_run, 4), ""),
             ("mae_retuned_Mm3", fmt(mae_ref, 4), ""), ("mae_run_Mm3", fmt(mae_run, 4), ""),
             ("total_gauge_Mm3", fmt(g[both].sum(), 4), "storms with all three volumes"), ("total_retuned_Mm3", fmt(rf[both].sum(), 4), ""),
             ("total_run_Mm3", fmt(rn[both].sum(), 4), ""), ("storms_up_vs_retuned", str(n_up), ""), ("storms_down_vs_retuned", str(n_down), "")]
    for nm, blk in type_rows.items():
        for span in ("seven", "whole"):
            for t in TYPES + ["total"]:
                summ.append(("made_{}_{}_{}_Mm3".format("retuned" if nm == "Re-tuned" else "run", span, t), fmt(blk[span][t], 4), ""))
    summ_path = out_dir / "score_wy2020_{}_summary.csv".format(tg)
    pd.DataFrame(summ, columns=["item", "value", "note"]).to_csv(summ_path, index=False)
    log("\nWrote:\n  {}\n  {}".format(per_path, summ_path))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
