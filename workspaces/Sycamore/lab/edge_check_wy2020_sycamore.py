#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
edge_check_wy2020_sycamore.py   (v1.0, 2026-11-07)
==================================================

READ-ONLY check of what the model and the gauge do BETWEEN and AFTER the seven scored WY2020 storms. It answers three questions that the
storm-by-storm figure (all_hydrographs_log1_v1_2.py) can only show by eye:

  1. LEFT EDGE   How much flow is there at the left edge of each storm panel (the start of the gray strip), for the gauge, Re-tuned and the run?
                 That is the leftover and slow flow from earlier storms, NOT an initial condition: the model runs continuously from Sep 1 2019.
  2. RIGHT EDGE  How much flow is there at the right edge (the end of the storm window, four days after the last rain day)? This is the
                 recession tail: the slow decline after the peak. Short model tails show up here.
  3. BY MONTH    How much water leaves the outlet each month (gauge, Re-tuned, run), how much of it falls OUTSIDE the seven scored windows, and how
                 much runoff the model MADE on the land outside those windows, by type. (This is where the "19.6 M m3 made outside the scored
                 storms" in the S1 scoring comes from.)

NOTHING HERE IS A TEST. No pass or fail was agreed for this check, and nothing is re-fit on WY2020. The script only prints and saves numbers, so the
eye reading of the figure can be checked against them. It runs no model and changes no file in results/.

WHAT IT READS (read-only; same files and same rules as score_wy2020_run_sycamore.py)
  results/Sycamore_WY2020_Retuned_Outlet.qout       the reference run (Time_hr 0 = Sep 1 2019 00:00)
  results/Sycamore_WY2020_<run>_Outlet.qout         the run (default S1_Ks1p0)
  ../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv       columns datetime, q_cms   (the gauge file starts Oct 1 2019)
  results/Sycamore_WY2020_Retuned*.rft and ..._<run>*.rft      runoff made by type; found by name pattern; skipped (with a message) if not found
WHAT IT WRITES (in edge_out/ unless you set --out-dir)
  edge_check_<run>_panel_edges.csv        one row per storm: left-edge and right-edge flows (CFS), ratios, gaps
  edge_check_<run>_monthly_outlet.csv     one row per month: outlet volumes (all hours, and outside the scored windows), with hours of data
  edge_check_<run>_monthly_made.csv       one row per month: runoff made OUTSIDE the scored windows, by type (only if the .rft files were found)

COMMANDS (from workspaces/Sycamore/lab/)
  python edge_check_wy2020_sycamore.py --dry-run      # finds and reads the files, lists the panel edges; prints no flows, writes nothing
  python edge_check_wy2020_sycamore.py                # S1 (run name S1_Ks1p0) against the gauge and Re-tuned
  python edge_check_wy2020_sycamore.py --run S2_Ks0p5          # another run; the prefix Sycamore_WY2020_ may be left off
  python edge_check_wy2020_sycamore.py --pad-days 0            # left edge = start of the storm window instead of 1 day before (figure default 1)
  python edge_check_wy2020_sycamore.py --no-types              # skip the runoff-made section
  python edge_check_wy2020_sycamore.py --rft results/<file>.rft --rft-ref results/<file>.rft     # name the .rft files yourself

HOW THE NUMBERS ARE MADE
  Flow at an edge = the HOURLY MEAN flow of the hour that starts at the left edge (left) or ends at the right edge (right), in CFS
  (cubic feet per second = m3/s x 35.3147). The figure draws every reading, so a line in the figure can differ a little from this hourly mean.
  Panel left edge = 00:00 on the first rain day minus --pad-days days (default 1, as in the figure). Panel right edge = the end of the storm window
  (00:00 four days after the last rain day).
  Monthly volume = hourly mean flow x 3,600 s added up over the hours of a calendar month, in M m3 (million cubic meters). The model starts Sep 1 2019
  and the gauge file Oct 1 2019, so September has no gauge. The run ends at hour 6,575 (the last hour of May 2020).
  Runoff made = what the model generates on the land (all cells added) BEFORE it is routed to the outlet, from the .rft file, by hour-END label with
  no shift (as in the earlier scripts). It is not outlet flow, so it is not compared with the gauge.
  The gauge and Re-tuned storm volumes are first checked against the WY2020 write-up (MATCH or MISMATCH). A MISMATCH means a wrong file or window: do
  not use the output.

WORDS
  left edge, right edge     the first and last moments drawn in a storm panel of the figure.
  recession tail            the slow decline of flow after a storm's peak.
  outside the scored windows   all hours that are not inside one of the seven storm windows.
  Hsrf infiltration excess (rain falls faster than the soil takes it)      Sbsrf subsurface stormflow (water coming out of the soil along the slope)
  Psrf perched-water runoff        Satsrf saturation excess (the soil is already full, so more rain runs off)
  M m3 million cubic meters        CFS cubic feet per second        M m3/day  a flow of 1 m3/s for one day is 0.0864 M m3
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
T0 = pd.Timestamp("2019-09-01 00:00")          # Time_hr 0 of every WY2020 run
RUN_HOURS = 6575
AFTER_DAYS = 3                                  # window ends at 00:00, (1 + 3) days after the last rain day's 00:00
CFS_PER_CMS = 35.3146667
MATCH_TOL = 0.005                               # M m3
WITHIN = 0.25                                   # same 25% band as the flags in the figure
MIN_GAUGE_CFS = 0.01                            # below this the gauge flow is treated as zero and no ratio is shown
MIN_GAUGE_VOL = 0.05                            # M m3; a monthly ratio to the gauge is not shown below this
COVER = 0.97                                    # share of the model's hours the gauge must have in a month for a ratio to be shown
PREFIX = "Sycamore_WY2020_"
REF_KEY, REF_STEM = "retuned", PREFIX + "Retuned"
DEFAULT_RUN = "S1_Ks1p0"
LABELS = {"S1_Ks1p0": "S1", "S2_Ks0p5": "S2", "S3_Bed0p5": "S3", "S4_GW0p5": "S4"}
TYPES = ["Hsrf", "Sbsrf", "Psrf", "Satsrf"]
# the S1 scoring (2026-11-07 15:57): whole run minus the seven storms, runoff made on the land, M m3
EARLIER_OUTSIDE = {"Re-tuned": 31.07, DEFAULT_RUN: 19.60}

STORMS = [
    ("Nov 19-22", "2019-11-19", "2019-11-22"), ("Nov 28-30", "2019-11-28", "2019-11-30"), ("Dec 23-28", "2019-12-23", "2019-12-28"),
    ("Feb 10-11", "2020-02-10", "2020-02-11"), ("Feb 22-23", "2020-02-22", "2020-02-23"),
    ("Mar 10-13", "2020-03-10", "2020-03-13"), ("Mar 18-19", "2020-03-18", "2020-03-19"),
]
EXPECTED = {   # WY2020 write-up volumes (M m3): gauge, Re-tuned
    "Nov 19-22": (0.015, 2.800), "Nov 28-30": (1.473, 3.235), "Dec 23-28": (4.069, 0.218), "Feb 10-11": (0.386, 1.767),
    "Feb 22-23": (4.131, 4.839), "Mar 10-13": (5.626, 0.875), "Mar 18-19": (3.317, 1.031),
}


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


def fmt_flow(v) -> str:
    """CFS, with fewer decimals as the flow gets bigger."""
    if v is None or not np.isfinite(v):
        return "n/a"
    if abs(v) < 10:
        return "{:.2f}".format(v)
    if abs(v) < 100:
        return "{:.1f}".format(v)
    return "{:.0f}".format(v)


def fmt_ratio(v) -> str:
    if v is None or not np.isfinite(v):
        return "n/a"
    if abs(v) < 10:
        return "{:.2f}".format(v)
    if abs(v) < 100:
        return "{:.1f}".format(v)
    return "{:.0f}".format(v)


# ------------------------------------------------------------------------------------------------ reading (same rules as the scoring script)
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
    pat = re.compile("^" + re.escape(stem) + r"\d.*\.rft$")
    return sorted(p for p in results.glob(stem + "*.rft") if pat.match(p.name))


def read_rft(path: Path):
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


def pad_type(v: str) -> int:
    try:
        n = int(v)
    except ValueError:
        raise argparse.ArgumentTypeError("--pad-days must be a whole number from 0 to 7")
    if not 0 <= n <= 7:
        raise argparse.ArgumentTypeError("--pad-days must be from 0 to 7 (you gave {})".format(n))
    return n


# ------------------------------------------------------------------------------------------------ numbers
def hourly(q: pd.Series) -> pd.Series:
    """Hourly mean flow, indexed by the start of each hour; only hours that have data."""
    return q.resample("60min").mean().dropna()


def volume(h: pd.Series, lo, hi) -> float:
    """Million m3 from hourly mean flow x 3,600 s over [lo, hi). NaN if there is no hour of data."""
    s = h[(h.index >= lo) & (h.index < hi)]
    return float(s.sum() * 3600.0 / 1e6) if len(s) else float("nan")


def hour_at(h: pd.Series, t) -> float:
    """Hourly mean flow (m3/s) of the hour that starts at t; NaN if that hour has no data."""
    return float(h.get(t, np.nan))


def in_windows(idx, windows) -> np.ndarray:
    m = np.zeros(len(idx), dtype=bool)
    for _, lo, hi in windows:
        m |= (np.asarray(idx >= lo)) & (np.asarray(idx < hi))
    return m


def by_month(h: pd.Series, mask=None):
    """(volume M m3, hours of data) for each calendar month, by the start of the hour. mask = True for the hours to keep."""
    s = h if mask is None else h[mask]
    per = s.index.to_period("M")
    vol = (s * 3600.0 / 1e6).groupby(per).sum()
    n = s.groupby(per).size()
    return vol, n


def word(ratio: float) -> str:
    if not np.isfinite(ratio):
        return "n/a"
    if ratio > 1 + WITHIN:
        return "above"
    if ratio < 1 - WITHIN:
        return "below"
    return "within 25%"


# ------------------------------------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Read-only check of the flow at the edges of the storm panels and by month (see the top of this file).")
    ap.add_argument("--root", default=str(HERE), help="the lab folder (default: the folder this script is in)")
    ap.add_argument("--results", default=None, help="folder with the result files (default <root>/results)")
    ap.add_argument("--obs", default=None, help="USGS WY2020 gauge file (default <root>/../init_data/obs/sycamore_09510200_obs_flow_WY2020.csv)")
    ap.add_argument("--out-dir", default=None, help="where the CSV files go (default <root>/edge_out)")
    ap.add_argument("--run", default=DEFAULT_RUN, help="the run to check, with or without the prefix {} (default {})".format(PREFIX, DEFAULT_RUN))
    ap.add_argument("--label", default=None, help="a short name for the run in the tables (default: S1 for S1_Ks1p0)")
    ap.add_argument("--pad-days", type=pad_type, default=1, help="days before the storm window where the panel's left edge is, 0 to 7 (default 1, as in the figure)")
    ap.add_argument("--rft", default=None, help="the run's .rft file, if its name is not found automatically")
    ap.add_argument("--rft-ref", default=None, help="the Re-tuned run's .rft file, if its name is not found automatically")
    ap.add_argument("--no-types", action="store_true", help="skip the runoff-made section")
    ap.add_argument("--dry-run", action="store_true", help="find and read the files and list the panel edges; print no flows, write nothing")
    a = ap.parse_args(argv)

    root = Path(a.root)
    results = Path(a.results) if a.results else root / "results"
    obs_path = Path(a.obs) if a.obs else root.parent / "init_data" / "obs" / "sycamore_09510200_obs_flow_WY2020.csv"
    out_dir = Path(a.out_dir) if a.out_dir else root / "edge_out"
    key, stem = parse_run(a.run)
    if key == REF_KEY:
        die("--run is the Re-tuned run itself. Give the run you want to compare with it, for example --run S1_Ks1p0.")
    label = (a.label or LABELS.get(key, key[:8])).strip()[:8] or key[:8]
    if label.lower().replace("-", "") == "retuned":
        die("the run's short name (--label) cannot be 'Re-tuned', because that name is used for the reference run. Pick another, for example --label S1.")
    L = label

    log("\nedge_check_wy2020_sycamore.py v1.0: reading files (read-only)")
    gauge = read_obs(need(obs_path, "WY2020 gauge"))
    ref_q = read_qout(need_run(results, REF_STEM, "the Re-tuned run"))
    run_q = read_qout(need_run(results, stem, "the run " + key))
    log("  gauge: {} readings, {:%b %d %Y} to {:%b %d %Y}".format(gauge.notna().sum(), gauge.index.min(), gauge.index.max()))
    log("  Re-tuned outlet: {} rows".format(len(ref_q)))
    log("  {} outlet [{}]: {} rows".format(L, stem, len(run_q)))
    for nm, q in (("Re-tuned", ref_q), (L, run_q)):
        last_h = (q.index.max() - T0) / pd.Timedelta(hours=1)
        if last_h < RUN_HOURS - 0.01:
            log("  WARNING: the {} outlet file ends at hour {:.1f}, not {} (the run did not finish, or the file is cut). "
                "Anything after that hour will show n/a.".format(nm, last_h, RUN_HOURS))

    rft = {}
    if a.no_types:
        log("  .rft files: skipped (--no-types)")
    else:
        for nm, spec, st in (("Re-tuned", a.rft_ref, REF_STEM), (L, a.rft, stem)):
            found = [need(Path(spec), nm + " .rft")] if spec else find_rft(results, st)
            if len(found) == 1:
                df, note = read_rft(found[0])
                rft[nm] = df
                log("  {} .rft: {} ({} rows){}".format(nm, found[0].name, len(df), "; " + note if note else ""))
            elif not found:
                log("  {} .rft: NOT FOUND (looked for {}*.rft with the run hours after the name in {}). The runoff-made section will be skipped; "
                    "use --rft / --rft-ref to name the file.".format(nm, st, results))
            else:
                log("  {} .rft: MORE THAN ONE file fits: {}. The runoff-made section will be skipped; use --rft / --rft-ref to name the one you want.".format(
                    nm, ", ".join(p.name for p in found)))
    types_ok = (not a.no_types) and ("Re-tuned" in rft) and (L in rft)

    windows = []
    for name, f, l in STORMS:
        lo = pd.Timestamp(f)
        hi = pd.Timestamp(l) + pd.Timedelta(days=1 + AFTER_DAYS)
        windows.append((name, lo, hi))
    log("\nStorm windows and panel edges (left edge = window start minus {} day(s); right edge = window end):".format(a.pad_days))
    for name, lo, hi in windows:
        log("  {:10s} left edge {:%Y-%m-%d %H:%M}   window {:%Y-%m-%d %H:%M} to {:%Y-%m-%d %H:%M}".format(
            name, lo - pd.Timedelta(days=a.pad_days), lo, hi))
    if a.dry_run:
        log("\nDRY RUN: all inputs were found and read. No flows were printed and nothing was written.")
        return 0

    hg, hr, hn = hourly(gauge), hourly(ref_q), hourly(run_q)

    # ---- 0. check against the write-up
    bad = []
    for name, lo, hi in windows:
        eg, er = EXPECTED[name]
        gv, rv = volume(hg, lo, hi), volume(hr, lo, hi)
        if not (np.isfinite(gv) and abs(gv - eg) <= MATCH_TOL):
            bad.append("{} gauge: here {} vs write-up {:.3f}".format(name, fmt(gv, 3), eg))
        if not (np.isfinite(rv) and abs(rv - er) <= MATCH_TOL):
            bad.append("{} Re-tuned: here {} vs write-up {:.3f}".format(name, fmt(rv, 3), er))
    log("\n0. Check of the gauge and Re-tuned storm volumes against the WY2020 write-up")
    if bad:
        log("   MISMATCH ({} value(s) differ by more than {:g} M m3):".format(len(bad), MATCH_TOL))
        for x in bad:
            log("     " + x)
        log("   Wrong file, window or run. Do not use the numbers below until this is explained.")
    else:
        log("   MATCH. Gauge and Re-tuned volumes equal the write-up values for all seven storms (within {:g} M m3).".format(MATCH_TOL))

    # ---- 1 and 2. panel edges
    rows = []
    for name, lo, hi in windows:
        t_left = lo - pd.Timedelta(days=a.pad_days)
        t_right = hi - pd.Timedelta(hours=1)                    # the hour that ENDS at the window end
        r = {"storm": name, "panel_left": t_left, "window_start": lo, "window_end": hi}
        for side, t in (("left", t_left), ("right", t_right)):
            g, rf, rn = (hour_at(x, t) * CFS_PER_CMS for x in (hg, hr, hn))
            r[side] = {"gauge": g, "ref": rf, "run": rn,
                       "run_over_gauge": rn / g if (np.isfinite(rn) and np.isfinite(g) and g >= MIN_GAUGE_CFS) else float("nan"),
                       "ref_over_gauge": rf / g if (np.isfinite(rf) and np.isfinite(g) and g >= MIN_GAUGE_CFS) else float("nan"),
                       "gap_Mm3_per_day": (rn - g) / CFS_PER_CMS * 86400.0 / 1e6 if (np.isfinite(rn) and np.isfinite(g)) else float("nan")}
        rows.append(r)

    def edge_table(side: str, title: str, when) -> None:
        log("\n{}".format(title))
        log("   Flow in CFS: hourly mean of the hour {}. 'gap' = ({} - gauge) in M m3 per day (positive = the model has more water going out).".format(
            "that starts at the panel's left edge" if side == "left" else "that ends at the window end (the panel's right edge)", L))
        log("   {:10s} {:16s} {:>8s} {:>9s} {:>8s} {:>8s} {:>9s}  {:11s} {:>9s}".format(
            "storm", "hour starting", "gauge", "Re-tuned", L, L + "/gauge", "Re-t/gauge", L + " is", "gap"))
        for r in rows:
            e = r[side]
            log("   {:10s} {:%Y-%m-%d %H:%M}   {:>8s} {:>9s} {:>8s} {:>8s} {:>9s}  {:11s} {:>9s}".format(
                r["storm"], when(r), fmt_flow(e["gauge"]), fmt_flow(e["ref"]), fmt_flow(e["run"]), fmt_ratio(e["run_over_gauge"]),
                fmt_ratio(e["ref_over_gauge"]), word(e["run_over_gauge"]), fmt_signed(e["gap_Mm3_per_day"], 3)))
        w = [word(r[side]["run_over_gauge"]) for r in rows]
        log("   {} is above the gauge in {} panel(s), within 25% in {}, below in {} (a ratio is n/a when the gauge flow is under {:g} CFS).".format(
            L, w.count("above"), w.count("within 25%"), w.count("below"), MIN_GAUGE_CFS))

    edge_table("left", "1. LEFT EDGE of each storm panel (leftover and slow flow; not an initial condition)", lambda r: r["panel_left"])
    edge_table("right", "2. RIGHT EDGE of each storm panel (the recession tail at the end of the storm window)", lambda r: r["window_end"] - pd.Timedelta(hours=1))
    for r1 in rows:
        for r2 in rows:
            if r1 is not r2 and r1["window_end"] == r2["panel_left"]:
                log("   Note: the right edge of {} and the left edge of {} are the same moment ({:%Y-%m-%d %H:%M}); the flows quoted are for the hour just before "
                    "and the hour just after it, so the two panels should show nearly the same values there.".format(r1["storm"], r2["storm"], r2["panel_left"]))

    # ---- 3. by month
    first = min(hn.index.min(), hr.index.min())
    last = max(hn.index.max(), hr.index.max())
    months = pd.period_range(first.to_period("M"), last.to_period("M"), freq="M")      # months of Re-tuned and the run together
    out_g, out_r, out_n = (~in_windows(hg.index, windows)), (~in_windows(hr.index, windows)), (~in_windows(hn.index, windows))
    vg, ng = by_month(hg)
    vr, nr = by_month(hr)
    vn, nn = by_month(hn)
    vgo, ngo = by_month(hg, out_g)
    vro, _ = by_month(hr, out_r)
    vno, _ = by_month(hn, out_n)
    M = pd.DataFrame({"month": months.astype(str)}, index=months)
    for nm, s in (("gauge", vg), ("ref", vr), ("run", vn), ("gauge_out", vgo), ("ref_out", vro), ("run_out", vno), ("gauge_hours", ng), ("ref_hours", nr), ("run_hours", nn)):
        M[nm] = s.reindex(months)
    for c in ("gauge_hours", "ref_hours", "run_hours"):
        M[c] = M[c].fillna(0).astype(int)
    M["model_hours"] = np.maximum(M["ref_hours"], M["run_hours"])         # hours the model has in the month (Re-tuned or the run, whichever has more)
    M["gauge_complete"] = (M["gauge_hours"] >= COVER * M["model_hours"]) & (M["gauge_hours"] > 0)

    def month_table(title: str, cg: str, cr: str, cn: str) -> None:
        log("\n{}".format(title))
        log("   {:8s} {:>6s} {:>9s} {:>9s} {:>9s} {:>9s} {:>10s}".format("month", "hours", "gauge", "Re-tuned", L, L + "/gauge", L + "-gauge"))
        tg = tr = tn = 0.0
        skipped = []
        for p, m in M.iterrows():
            gv = m[cg] if m["gauge_hours"] > 0 else float("nan")
            ratio = m[cn] / gv if (m["gauge_complete"] and np.isfinite(gv) and gv >= MIN_GAUGE_VOL and np.isfinite(m[cn])) else float("nan")
            diff = m[cn] - gv if (m["gauge_complete"] and np.isfinite(gv)) else float("nan")
            log("   {:8s} {:>6d} {:>9s} {:>9s} {:>9s} {:>9s} {:>10s}{}".format(
                m["month"], int(m["model_hours"]), fmt(gv, 3), fmt(m[cr], 3), fmt(m[cn], 3), fmt_ratio(ratio), fmt_signed(diff, 3),
                "" if m["gauge_complete"] else "   (gauge has {} of {} hours)".format(int(m["gauge_hours"]), int(m["model_hours"]))))
            if m["gauge_complete"] and np.isfinite(gv) and np.isfinite(m[cr]) and np.isfinite(m[cn]):
                tg += gv
                tr += m[cr]
                tn += m[cn]
            else:
                skipped.append(m["month"])
        log("   {:8s} {:>6s} {:>9s} {:>9s} {:>9s}".format("total*", "", fmt(tg, 3), fmt(tr, 3), fmt(tn, 3)))
        log("   * Only months where the gauge has data for nearly the whole month, so the three totals cover the same months{}.".format(
            "; left out: " + ", ".join(skipped) if skipped else ""))

    log("\n3. BY MONTH (calendar months by the start of each hour; M m3; ratio shown only where the gauge has data for nearly the whole month "
        "and at least {:g} M m3)".format(MIN_GAUGE_VOL))
    month_table("3a. Outlet volume, ALL hours", "gauge", "ref", "run")
    month_table("3b. Outlet volume OUTSIDE the seven scored windows (hours inside a storm window left out of all three)", "gauge_out", "ref_out", "run_out")

    made = None
    if types_ok:
        made = {}
        for nm in ("Re-tuned", L):
            d = rft[nm]
            hh = {t: hourly(d["{}_m3_s".format(t)]) for t in TYPES}
            idx = hh[TYPES[0]].index
            mask = ~in_windows(idx, windows)
            made[nm] = {t: by_month(hh[t], mask)[0].reindex(months) for t in TYPES}          # a month with no rows stays n/a, not 0
            made[nm]["total"] = pd.DataFrame({t: made[nm][t] for t in TYPES}).sum(axis=1, min_count=1)
        log("\n3c. Runoff the model MADE on the land OUTSIDE the seven scored windows, by type (M m3; all cells added; before routing, so not outlet flow)")
        for nm in ("Re-tuned", L):
            log("   {}".format(nm))
            log("   {:8s} {:>7s} {:>7s} {:>7s} {:>7s} {:>8s} {:>9s}".format("month", "Hsrf", "Sbsrf", "Psrf", "Satsrf", "total", "Satsrf %"))
            for p in months:
                v = {t: made[nm][t][p] for t in TYPES + ["total"]}
                log("   {:8s} {:>7s} {:>7s} {:>7s} {:>7s} {:>8s} {:>9s}".format(
                    str(p), fmt(v["Hsrf"], 3), fmt(v["Sbsrf"], 3), fmt(v["Psrf"], 3), fmt(v["Satsrf"], 3), fmt(v["total"], 3),
                    "{:.1f}".format(100 * v["Satsrf"] / v["total"]) if (np.isfinite(v["total"]) and v["total"] > 0) else "n/a"))
            s = {t: float(made[nm][t].sum()) for t in TYPES + ["total"]}
            log("   {:8s} {:>7s} {:>7s} {:>7s} {:>7s} {:>8s} {:>9s}".format(
                "total", fmt(s["Hsrf"], 3), fmt(s["Sbsrf"], 3), fmt(s["Psrf"], 3), fmt(s["Satsrf"], 3), fmt(s["total"], 3),
                "{:.1f}".format(100 * s["Satsrf"] / s["total"]) if s["total"] > 0 else "n/a"))
            earlier = EARLIER_OUTSIDE.get("Re-tuned" if nm == "Re-tuned" else key)
            if earlier is not None:
                same = abs(s["total"] - earlier) <= 0.01
                log("   (The scoring script's S1 run on 2026-11-07 gave {:.2f} M m3 for {} over the whole run minus the seven storms; this total {} it.)".format(
                    earlier, nm, "equals" if same else "DOES NOT EQUAL"))
    elif not a.no_types:
        log("\n3c. Runoff made outside the scored windows, by type: skipped (the .rft file(s) were not found; see the messages near the top).")

    # ---- write the CSV files
    out_dir.mkdir(parents=True, exist_ok=True)
    tg_ = tag_of(key)
    erows = []
    for r in rows:
        row = {"storm": r["storm"], "panel_left": r["panel_left"], "window_start": r["window_start"], "window_end": r["window_end"], "pad_days": a.pad_days}
        for side in ("left", "right"):
            e = r[side]
            row.update({"gauge_{}_cfs".format(side): e["gauge"], "retuned_{}_cfs".format(side): e["ref"], "run_{}_cfs".format(side): e["run"],
                        "run_over_gauge_{}".format(side): e["run_over_gauge"], "retuned_over_gauge_{}".format(side): e["ref_over_gauge"],
                        "run_is_{}".format(side): word(e["run_over_gauge"]), "gap_{}_Mm3_per_day".format(side): e["gap_Mm3_per_day"]})
        row["run"] = key
        erows.append(row)
    p1 = out_dir / "edge_check_{}_panel_edges.csv".format(tg_)
    pe = pd.DataFrame(erows)
    for c in ("panel_left", "window_start", "window_end"):
        pe[c] = pe[c].dt.strftime("%Y-%m-%d %H:%M")
    pe.round(4).to_csv(p1, index=False)
    mo = pd.DataFrame({
        "month": M["month"], "model_hours": M["model_hours"], "run_hours": M["run_hours"], "gauge_hours": M["gauge_hours"],
        "gauge_Mm3": M["gauge"].where(M["gauge_hours"] > 0), "retuned_Mm3": M["ref"], "run_Mm3": M["run"],
        "gauge_outside_windows_Mm3": M["gauge_out"].where(M["gauge_hours"] > 0), "retuned_outside_windows_Mm3": M["ref_out"], "run_outside_windows_Mm3": M["run_out"],
        "run": key})
    p2 = out_dir / "edge_check_{}_monthly_outlet.csv".format(tg_)
    mo.round(4).to_csv(p2, index=False)
    wrote = [p1, p2]
    if made is not None:
        mm = pd.DataFrame({"month": months.astype(str)}, index=months)
        for nm, tag in (("Re-tuned", "retuned"), (L, "run")):
            for t in TYPES + ["total"]:
                mm["{}_{}_made_outside_Mm3".format(tag, t)] = made[nm][t]
        mm["run"] = key
        p3 = out_dir / "edge_check_{}_monthly_made.csv".format(tg_)
        mm.round(4).to_csv(p3, index=False)
        wrote.append(p3)
    log("\nWrote:\n  " + "\n  ".join(str(p) for p in wrote))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
