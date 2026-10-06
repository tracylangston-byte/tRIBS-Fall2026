"""
score_storms_110.py   (version 1)
=================================
Series 110 -- STORM-BY-STORM scoring of the long run (label ms1: 240 runs, cc ON, 1 Aug 2014 00:00 for 1416 h)
against the REAL SMF gauge, and a reading of where Ks, f and cc sit. No tRIBS runs; it only reads files.

WHAT IT DOES, IN PLAIN WORDS
-----------------------------
Every ms1 run left one long flow record (its Outlet .qout). This script cuts four storms out of each record
(12 Aug, 19 Aug, 8 Sep, 27 Sep), scores each cut against the real gauge, combines the scores, and reports
which values of Ks, f and cc the best runs share. The windows come from storm_windows_110.json, which
define_storm_windows_110.py wrote from the REAL gauge alone before any ms1 result was looked at.

THE RULE (approved by Tracy 4 Oct 2026; every number below is a constant you can read and change BEFORE the run)
------------------------------------------------------------------------------------------------------------
1. Windows: 12 Aug keeps the Stage B window (16:00 12 Aug to 12:00 13 Aug). The others start 2 h before and end
   6 h after the real flow at or above 1 cfs (see define_storm_windows_110.py).
2. Score per run per storm.
   * Big storms (12 Aug, 8 Sep, 27 Sep): KGE_2012 on the 5-minute grid, the same formula and the same reading of
     the .qout as the scorer. Two ways to treat the gauge's gaps ("arms"), as in Stage B: FILLED (primary; all
     5-minute bins, gaps filled from the record) and AS-IS (only bins that hold a reading). Same code as
     rescore_real_gauge_110.py.
   * 19 Aug (13 cfs, under 1 acre-foot) is too small for KGE. Three plain numbers instead: did the run make
     at least 1 cfs, model peak / real peak, model volume / real volume. A run PASSES if it made flow and both
     ratios are between 1/3 and 3. Pass/fail is shown beside the combined score, not folded into it.
   * A run that makes under 1 cfs in a storm that really flowed gets the WORST score on that storm (0 after
     rescaling). It is not left out: leaving it out would reward runs for missing storms.
3. Combine (Ty and Josh's rule): rescale each big storm's KGE across the runs, worst run 0 and best run 1
   (min-max), then average over the big storms. Each run's worst single-storm score is shown beside it, for
   display only. EDGE CASE I had to settle (NOT in the approved text): if fewer than MIN_FLOW_RUNS runs make flow
   in a storm, or all scored runs get the same score, rescaling is meaningless, so that storm is left OUT of the
   mean and reported on its own (it is never silently dropped; the report says so). If a storm's best raw KGE is
   under WEAK_KGE it is flagged; leave-one-out (below) then shows the result without it.
4. What counts as CLEAR. Take the best TOP_FRAC of the runs (10 per cent: 24 of 240) by combined score.
   For each of Ks (linear), f and cc (log), measure position in the sampled box from 0 to 1.
   * NARROWED: the middle 80 per cent of the best runs (10th to 90th percentile) cover less than NARROW_SPAN
     (half) of the box.
   * EDGE: the median position of the best runs is within EDGE_BAND (0.15) of either end of the box. The box was
     too small on that side: no claim is made.
   * BRACKETED: cut the box into N_SLICES (5) equal slices; take the mean combined score of ALL runs in each.
     The best slice must be an interior one (not the first or last) and beat BOTH end slices by at least
     BRACKET_GAP (0.10 on the 0-to-1 scale). A slice with fewer than MIN_SLICE_RUNS runs does not count.
   * Labels: PINNED (narrowed and bracketed), NARROW (narrowed, not bracketed), WIDE-PEAK (bracketed, not
     narrowed), FLAT (neither), EDGE-LO / EDGE-HI (edge wins over everything else).
   * The reading is repeated leaving each big storm out in turn. A parameter is CLEAR only if it is PINNED in the
     all-storms version and in every leave-one-out version. 12 Aug ALONE is shown too, as the reference the extra
     storms are compared with.
   * Also shown (information only, not part of the verdict): the Spearman correlation between parameters among
     the best runs (|rho| of RIDGE_RHO or more means they trade off and should be read as a pair), and whether
     the 19 Aug passers overlap the best runs, with the overlap chance would give.

CAVEATS (printed again at the end of every run)
------------------------------------------------
All four storms use the same two rain gauges, so the SMPHQ timing and intensity offset is in all of them. One
monsoon season. Routing is pinned. Gauge error is not quantified (no uncertainty statement is made; that needs a
cited measurement-error value). The runs are cc ON only, so this shows where cc, Ks and f sit, not whether cc is
needed. In the filled arm each 5-minute real bin is labelled by its start while the model value is AT that time
(about 2.5 minutes of offset, as in all earlier scoring). The scorer puts the 3.75-minute .qout onto a 5-minute grid
by linear interpolation in time between the neighbouring samples (resample('5min').interpolate('time'), tested:
it equals an interpolation from all the samples); this script makes the very same call so its 12 Aug numbers equal
the stored ones.

SELF-CHECKS ("gates"; the script stops if one fails)
-----------------------------------------------------
  windows file     the right four storms, hours consistent with 1 Aug 00:00, same run length as the design,
                   and no real flow above 0 cfs outside the windows
  results/design   the ms1 design and results agree (n, ON series, run length); partial sets stop unless
                   --allow_partial (they are then labelled PARTIAL everywhere)
  .qout files      every run has its kept Outlet .qout, and it reaches the end of the run
  GATE 2           the 12 Aug window of every run is re-scored against the SYNTHETIC truth (the single .qout in
                   synth_truth/) and must reproduce the metrics stored in the results CSV (KGE_2012, r, beta,
                   gamma, PBIAS, volumes, peak, RMSE, NSE). This proves the .qout reading, the window cut and the
                   metric code, on real ms1 files, before any real-gauge score is trusted.
  GATE 1           (only if present) the 12 Aug real series built here equals real_gauge_observed_110.csv from
                   the Stage 0 script, bin for bin.

USAGE (from the lab/ folder; the ms1 sweep must have finished, and nothing may be rewriting its files)
-------------------------------------------------------------------------------------------------------
    python score_storms_110.py --check_only                 # gates only; writes nothing
    python score_storms_110.py                              # gates, then the full scoring and reading
    python score_storms_110.py --allow_partial              # a partial sweep (labelled PARTIAL)
    python score_storms_110.py --label ms1 --primary_arm filled --no_plots

WHAT YOU GET   (calibration_work/03_comparisons/summary_tables/storm_scores_<label>_110/)
------------------------------------------------------------------------------------------
    storm_scores_filled_110.csv / storm_scores_asis_110.csv   one row per run: parameters, per-storm KGE, r, beta,
                                gamma, PBIAS, rescaled score, no-flow flag, peak/volume ratios, 19 Aug pass,
                                combined scores for every version
    storm_ranges_110.csv        per arm and storm: runs making flow, best/median/worst raw KGE, flag
    storm_readings_110.csv      per arm, version and parameter: label, span, median, 10th-90th values, slice means
    storm_top_runs_110.csv      the best runs of each version
    fig_storm_hydrographs_110.png  real gauge with the best runs, four storms
    fig_storm_param_scores_110.png rescaled score against Ks, f and cc, storm by storm
    fig_storm_19aug_110.png     the 19 Aug pass/fail against Ks, f and cc
    fig_storm_versions_110.png  where the best runs sit, version by version
    PROVENANCE_storm_scores_110.json  checksums, versions, gate results, settings
It never writes outside that folder and never changes an input file.
"""

import argparse
import hashlib
import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

VERSION = "1"

# ------------------------------------------------------------------
# Constants. The first block must agree with the scorer; the second block is the approved rule.
# ------------------------------------------------------------------
GAUGE_REL = ("init_data", "met", "SMF_Observations_1993-2025.xlsx")      # relative to the project root (lab/'s parent)
GAUGE_SHEET = "Discharge"
GAUGE_SKIPROWS = 6
CFS_TO_CMS = 0.0283168
M3_PER_ACFT = 1233.4818
GRID_FREQ = "5min"
GRID_SEC = 300.0
ORIGIN = pd.Timestamp("2014-08-01")                 # .qout time is hours from here (scorer)
QOUT_NAMES = ["Time_hr", "Qstrm_m3s", "Hlev_m"]
STAGE_B_START, STAGE_B_END = pd.Timestamp("2014-08-12 16:00"), pd.Timestamp("2014-08-13 12:00")
RTOL, ATOL = 1e-7, 1e-9                             # gate 2: stored against recomputed metrics
EXPECT_PANDAS_MINOR = "3.0"
COVER_TOL_H = 2.0                                   # a .qout must reach within this many hours of the run length
WINDOWS_NAME = "storm_windows_110.json"
LABELS_EXPECTED = ["12 Aug", "19 Aug", "8 Sep", "27 Sep"]
BIG = ["12 Aug", "8 Sep", "27 Sep"]                 # scored by KGE and combined
SMALL = "19 Aug"                                    # scored by flow / peak / volume checks
DEFAULT_MAX_GAP_MIN = 60.0                          # filled arm: straight line between readings at most this far apart
DEFAULT_MAX_HOLD_HR = 24.0                          # filled arm: hold the last reading at most this long
MIN_ASIS_POINTS = 5                                 # fewer scorable as-is bins than this: that storm is not scored in the as-is arm
GATE2_COLS = ["kge_2012", "kge_r", "kge_beta", "kge_gamma", "pbias_pct", "sim_volume_m3", "obs_volume_m3",
              "sim_peak_m3s", "rmse_m3s", "nse"]

# ---- the approved rule ----------------------------------------------------------------------
FLOW_CFS = 1.0                  # a run "makes flow" in a storm if its highest flow in the window reaches this
FLOW_CMS = FLOW_CFS * CFS_TO_CMS
PASS_LO, PASS_HI = 1.0 / 3.0, 3.0   # 19 Aug: peak and volume ratios must both fall between these
WEAK_KGE = 0.5                  # a storm whose best raw KGE_2012 is below this is flagged
MIN_FLOW_RUNS = 10              # fewer runs than this making flow in a storm: the storm is left out of the mean
TOP_FRAC = 0.10                 # share of the runs counted as "the best runs"
TOP_MIN = 5
NARROW_SPAN = 0.5               # best runs' 10th-90th percentile span, as a share of the box, below which = narrowed
EDGE_BAND = 0.15                # median position within this of either end of the box = piled at an edge
N_SLICES = 5                    # equal slices of each box for the bracket test
BRACKET_GAP = 0.10              # the best slice must beat BOTH end slices by this much (0-to-1 score scale)
MIN_SLICE_RUNS = 10
RIDGE_RHO = 0.6                 # information only
ARM_OVERLAP = 0.5               # arms "agree" if at least this share of the best runs is common (and labels match)

PARAMS = [("Ks_mult", "Ks", "linear"), ("f_RS_abs", "f", "log"), ("channelconductivity_mmhr", "cc", "log")]
PNAME = {p[0]: p[1] for p in PARAMS}
OUT_PREFIX = "storm_scores_"

# palette (same as the other Series 110 scripts); one fixed hue per parameter
C_KS, C_F, C_CC = "#2a78d6", "#eb6834", "#1baf7a"
PCOLOR = {"Ks_mult": C_KS, "f_RS_abs": C_F, "channelconductivity_mmhr": C_CC}
INK, INK2, GRID, SURFACE, NULLC = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb", "#8a8984"


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
class Checks(object):
    def __init__(self):
        self.rows = []

    def add(self, status, name, text):
        self.rows.append((status, name, text))

    def failures(self):
        return [r for r in self.rows if r[0] == "FAIL"]

    def show(self):
        for s, n, t in self.rows:
            print("  [%-4s] %s: %s" % (s, n, t))


def tag_of(label):
    return label.replace(" ", "")


def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def write_csv_atomic(df, path):
    tmp = path.with_name(path.name + ".part")
    df.to_csv(tmp, index=False)
    os.replace(str(tmp), str(path))


def write_json_atomic(path, obj):
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=2, default=str))
    os.replace(str(tmp), str(path))


def fnum(v, spec):
    try:
        if v is None or not np.isfinite(v):
            return "n/a"
    except TypeError:
        return "n/a"
    return format(v, spec)


def hours_from_origin(t):
    return (t - ORIGIN).total_seconds() / 3600.0


def fmt_t(t):
    return t.strftime("%d %b %H:%M")


def unit_pos(values, lo, hi, scale):
    v = np.asarray(values, dtype=float)
    if scale == "log":
        return (np.log10(v) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
    return (v - lo) / (hi - lo)


def from_unit(u, lo, hi, scale):
    u = np.asarray(u, dtype=float)
    if scale == "log":
        return 10.0 ** (np.log10(lo) + u * (np.log10(hi) - np.log10(lo)))
    return lo + u * (hi - lo)


def fval(name, v):
    """A parameter value for printing."""
    if v is None or not np.isfinite(v):
        return "n/a"
    if name == "Ks_mult":
        return "%.2f" % v
    if name == "f_RS_abs":
        return "%.5f" % v
    return "%.0f" % v


# ------------------------------------------------------------------
# Reading a .qout the way the scorer does, and the scorer's metrics
# ------------------------------------------------------------------
def read_qout_5min(path):
    """Same calls as rescore_cc_truths_110.read_truth_5min / the scorer: sniff the format (tRIBS 6.0.0 is
    comma-delimited, older tRIBS whitespace), hours from 1 Aug 2014 00:00, resample('5min').interpolate('time')."""
    with open(path) as fh:
        fh.readline()
        first_data = fh.readline()
    if "," in first_data:
        raw = pd.read_csv(path, sep=",", header=0, names=QOUT_NAMES)
    else:
        raw = pd.read_csv(path, sep=r"\s+", skiprows=1, names=QOUT_NAMES)
    raw["datetime"] = pd.to_datetime(raw["Time_hr"] * 3600, unit="s", origin=ORIGIN)
    raw = raw.set_index("datetime")
    return raw["Qstrm_m3s"].resample(GRID_FREQ).interpolate(method="time")


def compute_metrics(obs, sim):
    """obs, sim: Series on the same DatetimeIndex. The scorer's formulas (as in rescore_cc_truths_110.compute_metrics),
    without the phase metrics."""
    with np.errstate(all="ignore"):
        obs_peak, sim_peak = obs.max(), sim.max()
        dt_seconds = (obs.index[1] - obs.index[0]).total_seconds() if len(obs) > 1 else GRID_SEC
        obs_vol, sim_vol = obs.sum() * dt_seconds, sim.sum() * dt_seconds
        rmse = np.sqrt(np.mean((sim - obs) ** 2))
        nse = 1 - (np.sum((sim - obs) ** 2) / np.sum((obs - obs.mean()) ** 2))
        pbias = 100 * (np.sum(sim - obs) / np.sum(obs))
        r = np.corrcoef(sim, obs)[0, 1]
        alpha = np.std(sim) / np.std(obs)
        beta = np.mean(sim) / np.mean(obs)
        gamma = alpha / beta
        kge_2012 = 1 - np.sqrt((r - 1) ** 2 + (gamma - 1) ** 2 + (beta - 1) ** 2)
    return {"obs_peak_m3s": float(obs_peak), "sim_peak_m3s": float(sim_peak), "obs_volume_m3": float(obs_vol),
            "sim_volume_m3": float(sim_vol), "rmse_m3s": float(rmse), "nse": float(nse), "pbias_pct": float(pbias),
            "kge_r": float(r), "kge_beta": float(beta), "kge_gamma": float(gamma), "kge_2012": float(kge_2012)}


# ------------------------------------------------------------------
# The real gauge, per storm window (same recipe as rescore_real_gauge_110.py, applied to any window)
# ------------------------------------------------------------------
def read_gauge(path):
    """Readings in m3/s indexed by time (same calls as the scorer's real-gauge mode)."""
    obs = pd.read_excel(path, sheet_name=GAUGE_SHEET, skiprows=GAUGE_SKIPROWS)
    obs["datetime"] = pd.to_datetime(obs["Date"].astype(str) + " " + obs["Time"].astype(str))
    obs["cfs"] = pd.to_numeric(obs["cfs"], errors="coerce")
    obs = obs.dropna(subset=["datetime", "cfs"]).sort_values("datetime")
    return obs.set_index("datetime")["cfs"] * CFS_TO_CMS


def build_real_window(raw, ws, we, max_gap_min, max_hold_hr):
    """Returns (asis, filled, ftype) on the 5-minute grid ws..we inclusive.
    asis   the scorer's bin means (cfs x 0.0283168, .resample('5min').mean()); NaN where a bin holds no reading
    filled the same values plus the empty bins filled: a straight line between readings at most max_gap_min apart,
           otherwise the last reading held for at most max_hold_hr (a logger that writes only on change writes
           nothing when nothing changes); NaN ('uncovered') if there is no recent enough reading
    """
    grid = pd.date_range(ws, we, freq=GRID_FREQ)
    win = raw[(raw.index >= ws) & (raw.index < we + pd.Timedelta(minutes=5))]
    asis = win.resample(GRID_FREQ).mean().reindex(grid)
    keep = raw[(raw.index >= ws - pd.Timedelta(hours=max_hold_hr + 2)) & (raw.index <= we + pd.Timedelta(hours=2))]
    keep = keep.groupby(level=0).mean().sort_index()
    epoch = pd.Timestamp("1970-01-01")
    sec = np.asarray((keep.index - epoch) / pd.Timedelta(seconds=1), dtype=float)
    val = keep.to_numpy(dtype=float)
    asis_v = asis.to_numpy(dtype=float)
    gsec = np.asarray((grid - epoch) / pd.Timedelta(seconds=1), dtype=float)
    half = 0.5 * GRID_SEC
    max_gap, max_hold = float(max_gap_min) * 60.0, float(max_hold_hr) * 3600.0
    fv = asis_v.copy()
    ft = np.empty(len(grid), dtype=object)
    for k in range(len(grid)):
        if np.isfinite(asis_v[k]):
            ft[k] = "record"
            continue
        c = gsec[k] + half
        i = int(np.searchsorted(sec, c))
        if i == 0:
            ft[k] = "uncovered"
            continue
        t_prev, v_prev = sec[i - 1], val[i - 1]
        if i < len(sec) and (sec[i] - t_prev) <= max_gap:
            w = (c - t_prev) / (sec[i] - t_prev)
            fv[k] = v_prev + w * (val[i] - v_prev)
            ft[k] = "linear"
        elif (c - t_prev) <= max_hold:
            fv[k] = v_prev
            ft[k] = "hold"
        else:
            ft[k] = "uncovered"
    return asis, pd.Series(fv, index=grid, name="filled"), pd.Series(ft, index=grid, name="fill_type")


# ------------------------------------------------------------------
# Loading the inputs
# ------------------------------------------------------------------
def load_windows(path, ck):
    """Returns the list of storm dicts, or None (a FAIL is recorded)."""
    if not path.exists():
        ck.add("FAIL", "windows file", "%s not found here. Run define_storm_windows_110.py first." % path.name)
        return None, None
    try:
        d = json.loads(path.read_text())
        storms = []
        for s in d["storms"]:
            ws, we = pd.Timestamp(s["window_start"]), pd.Timestamp(s["window_end"])
            storms.append({"label": s["label"], "tag": tag_of(s["label"]), "date": s["date"], "ws": ws, "we": we,
                           "start_hour": float(s["start_hour"]), "end_hour": float(s["end_hour"]),
                           "grid": pd.date_range(ws, we, freq=GRID_FREQ)})
        runtime = int(d["runtime_hours"])
    except Exception as exc:
        ck.add("FAIL", "windows file", "cannot read %s: %s" % (path.name, exc))
        return None, None
    labels = [s["label"] for s in storms]
    if labels != LABELS_EXPECTED:
        ck.add("FAIL", "windows file", "storms are %s, expected %s" % (labels, LABELS_EXPECTED))
        return None, None
    bad = []
    for s in storms:
        if abs(hours_from_origin(s["ws"]) - s["start_hour"]) > 1e-6 or abs(hours_from_origin(s["we"]) - s["end_hour"]) > 1e-6:
            bad.append("%s: stated hours do not match the window times" % s["label"])
        if s["we"] <= s["ws"] or s["end_hour"] > runtime:
            bad.append("%s: window is empty or runs past the %d-hour run" % (s["label"], runtime))
    for a, b in zip(storms, storms[1:]):
        if a["we"] > b["ws"]:
            bad.append("%s and %s overlap" % (a["label"], b["label"]))
    if (storms[0]["ws"], storms[0]["we"]) != (STAGE_B_START, STAGE_B_END):
        ck.add("NOTE", "windows file", "the 12 Aug window is not the Stage B window; gate 2 (the stored-metric check) "
               "is skipped because the stored metrics were cut on the Stage B window")
    if bad:
        ck.add("FAIL", "windows file", "; ".join(bad))
        return None, None
    ck.add("PASS", "windows file", "%s: four storms, %s; hours consistent with 1 Aug 00:00; run length %d h"
           % (path.name, ", ".join("%s %s to %s (h %.0f-%.0f)" % (s["label"], fmt_t(s["ws"]), fmt_t(s["we"]),
                                                                  s["start_hour"], s["end_hour"]) for s in storms), runtime))
    return storms, runtime


def read_design(path, label, ck):
    if not path.exists():
        ck.add("FAIL", "design file", "%s not found (has the ms1 sweep been started?)" % path.name)
        return None
    try:
        d = json.loads(path.read_text())
        box = d["box"]
        out = {"path": path, "n": int(d["n"]), "seed": int(d["seed"]), "runtime_hours": d.get("runtime_hours"),
               "f_fixed": d.get("f_fixed"), "routing": d.get("routing_pinned"), "version": d.get("version"),
               "box": {k: (float(box[k]["lo"]), float(box[k]["hi"]), box[k]["scale"]) for k, _, _ in PARAMS}}
    except Exception as exc:
        ck.add("FAIL", "design file", "cannot read %s: %s" % (path.name, exc))
        return None
    if d.get("series") != "ON" or d.get("label") != label:
        ck.add("FAIL", "design file", "%s is for series %s label %s, not ON / %s" % (path.name, d.get("series"), d.get("label"), label))
        return None
    if out["f_fixed"] is not None:
        ck.add("FAIL", "design file", "this design holds f fixed at %s; this scorer reads f, so it needs a sweep that samples f" % out["f_fixed"])
        return None
    for k, nm, sc in PARAMS:
        if out["box"][k][2] != sc:
            ck.add("FAIL", "design file", "box scale of %s is %s, expected %s" % (k, out["box"][k][2], sc))
            return None
    ck.add("PASS", "design file", "%s: label %s, ON, n %d, seed %d, run length %s h, Ks %g-%g, f %g-%g, cc %g-%g"
           % (path.name, label, out["n"], out["seed"], out["runtime_hours"], out["box"]["Ks_mult"][0], out["box"]["Ks_mult"][1],
              out["box"]["f_RS_abs"][0], out["box"]["f_RS_abs"][1], out["box"]["channelconductivity_mmhr"][0],
              out["box"]["channelconductivity_mmhr"][1]))
    return out


def load_results(path, design, ck, allow_partial):
    if not path.exists():
        ck.add("FAIL", "results file", "%s not found" % path.name)
        return None, False
    df = pd.read_csv(path)
    need = ["run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "optpercolation"] + GATE2_COLS
    missing = [c for c in need if c not in df.columns]
    if missing:
        ck.add("FAIL", "results file", "%s is missing columns %s" % (path.name, missing))
        return None, False
    if df["run_id"].duplicated().any():
        ck.add("FAIL", "results file", "duplicate run_id values in %s" % path.name)
        return None, False
    for c in ["Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "optpercolation"] + GATE2_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    if sorted(df["optpercolation"].dropna().unique().tolist()) != [1.0]:
        ck.add("FAIL", "results file", "optpercolation values are not all 1 (not a cc ON file)")
        return None, False
    lo_hi = {k: design["box"][k] for k, _, _ in PARAMS}
    outside = []
    for k, _, _ in PARAMS:
        lo, hi, _ = lo_hi[k]
        if ((df[k] < lo * (1 - 1e-9)) | (df[k] > hi * (1 + 1e-9))).any():
            outside.append(k)
    if outside:
        ck.add("FAIL", "results file", "some runs lie outside the design box in %s" % outside)
        return None, False
    partial = len(df) < design["n"]
    if len(df) > design["n"]:
        ck.add("FAIL", "results file", "%d rows but the design has n = %d" % (len(df), design["n"]))
        return None, False
    if partial and not allow_partial:
        ck.add("FAIL", "results file", "only %d of %d runs are in %s. The reading rule is meant for the finished sweep; "
               "finish the sweep (run_lowf_sweep_110.py ... --skip_existing) or add --allow_partial" % (len(df), design["n"], path.name))
        return None, True
    ck.add("WARN" if partial else "PASS", "results file", "%s: %d of %d runs%s" % (path.name, len(df), design["n"],
           "   <-- PARTIAL; every table is labelled so" if partial else ""))
    return df.reset_index(drop=True), partial


# ------------------------------------------------------------------
# Scoring every run: one pass over the kept .qout files
# ------------------------------------------------------------------
def score_runs(res, keep_dir, storms, real, runtime_hours, synth5, ck, drop_bad):
    """real: {tag: {'asis','filled','ftype'}}. Returns (records list, sim_cfs dict, bad dict)."""
    s12 = storms[0]
    do_gate2 = synth5 is not None and (s12["ws"], s12["we"]) == (STAGE_B_START, STAGE_B_END)
    obs12 = synth5.reindex(s12["grid"]) if do_gate2 else None
    if do_gate2 and obs12.isna().any():
        ck.add("FAIL", "GATE 2", "the synthetic truth has no value at some 12 Aug window times")
        do_gate2 = False
    recs, sim_cfs, bad = [], {}, {}
    worst2, n2bad, n2ok = {}, 0, 0
    total = len(res)
    for i, row in res.iterrows():
        rid = row["run_id"]
        p = keep_dir / ("%s_Outlet.qout" % rid)
        if not p.exists() or p.stat().st_size == 0:
            bad[rid] = "the kept Outlet .qout is missing"
            continue
        try:
            sim5 = read_qout_5min(p)
        except Exception as exc:
            bad[rid] = "the .qout is unreadable: %s" % str(exc)[:80]
            continue
        last_h = hours_from_origin(sim5.index[-1])
        if last_h < runtime_hours - COVER_TOL_H:
            bad[rid] = "the .qout stops at hour %.1f of the %d-hour run" % (last_h, runtime_hours)
            continue
        rec = {"run_id": rid, "Ks_mult": float(row["Ks_mult"]), "f_RS_abs": float(row["f_RS_abs"]),
               "channelconductivity_mmhr": float(row["channelconductivity_mmhr"])}
        sims = {}
        ok = True
        for s in storms:
            w = sim5.reindex(s["grid"])
            if w.isna().any():
                bad[rid] = "the .qout has no value at some %s window times" % s["label"]
                ok = False
                break
            sims[s["tag"]] = w
        if not ok:
            continue
        if do_gate2:
            m = compute_metrics(obs12, sims[s12["tag"]])
            off = []
            for c in GATE2_COLS:
                stored, new = float(row[c]), m[c]
                if not np.isclose(new, stored, rtol=RTOL, atol=ATOL, equal_nan=True):
                    off.append("%s (stored %.6g, recomputed %.6g)" % (c, stored, new))
                elif np.isfinite(stored) and np.isfinite(new):
                    worst2[c] = max(worst2.get(c, 0.0), abs(new - stored))
            if off:
                bad[rid] = "GATE 2: " + "; ".join(off[:3])
                n2bad += 1
                continue
            n2ok += 1
        for s in storms:
            tg = s["tag"]
            w = sims[tg]
            fl = real[tg]["filled"]
            sim_peak_cms = float(w.max())
            rec[tg + "_sim_peak_cfs"] = sim_peak_cms / CFS_TO_CMS
            rec[tg + "_sim_vol_acft"] = float(w.sum()) * GRID_SEC / M3_PER_ACFT
            rec[tg + "_made_flow"] = bool(sim_peak_cms >= FLOW_CMS)
            rec[tg + "_peak_ratio"] = rec[tg + "_sim_peak_cfs"] / real[tg]["real_peak_cfs"]
            rec[tg + "_vol_ratio"] = rec[tg + "_sim_vol_acft"] / real[tg]["real_vol_acft"]
            if rec[tg + "_made_flow"]:
                rec[tg + "_peak_dt_min"] = (w.idxmax() - fl.idxmax()).total_seconds() / 60.0
            else:
                rec[tg + "_peak_dt_min"] = np.nan
            for arm in ("filled", "asis"):
                obs = real[tg][arm]
                if obs is None:
                    continue
                o = obs.dropna()
                mm = compute_metrics(o, w.reindex(o.index))
                for k in ("kge_2012", "kge_r", "kge_beta", "kge_gamma", "pbias_pct"):
                    rec["%s_%s_%s" % (arm, tg, k)] = mm[k]
        rec["pass_19Aug"] = bool(rec[tag_of(SMALL) + "_made_flow"] and PASS_LO <= rec[tag_of(SMALL) + "_peak_ratio"] <= PASS_HI
                                 and PASS_LO <= rec[tag_of(SMALL) + "_vol_ratio"] <= PASS_HI)
        recs.append(rec)
        sim_cfs[rid] = {tg: sims[tg].to_numpy(dtype=float) / CFS_TO_CMS for tg in sims}
        if (i + 1) % 40 == 0 or i + 1 == total:
            print("    scored %d of %d runs ..." % (i + 1, total))
    if do_gate2:
        worst_txt = ""
        if worst2:
            wk = max(worst2, key=worst2.get)
            worst_txt = " (worst absolute deviation %.2e in %s)" % (worst2[wk], wk)
        ck.add("FAIL" if n2bad and not drop_bad else ("WARN" if n2bad else "PASS"), "GATE 2",
               "%d of %d runs reproduce their stored metrics when the 12 Aug window is re-scored against the synthetic "
               "truth%s%s" % (n2ok, n2ok + n2bad, worst_txt, "" if not n2bad else "; %d run(s) DO NOT (listed below)" % n2bad))
    else:
        ck.add("WARN", "GATE 2", "not run (no synthetic truth read, or the 12 Aug window is not the Stage B window): "
               "the .qout reading and metric code were NOT checked against the stored metrics")
    return recs, sim_cfs, bad


# ------------------------------------------------------------------
# Tables, rescaling, versions
# ------------------------------------------------------------------
SHARED_KEYS = ["sim_peak_cfs", "sim_vol_acft", "made_flow", "peak_ratio", "vol_ratio", "peak_dt_min"]
ARM_KEYS = ["kge_2012", "kge_r", "kge_beta", "kge_gamma", "pbias_pct"]


def arm_table(df_all, arm, storms):
    """One row per run for one arm: parameters, per-storm raw metrics (renamed without the arm prefix)."""
    t = df_all[["run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr"]].copy()
    for s in storms:
        tg = s["tag"]
        for k in SHARED_KEYS:
            t["%s_%s" % (tg, k)] = df_all["%s_%s" % (tg, k)].to_numpy()
        for k in ARM_KEYS:
            col = "%s_%s_%s" % (arm, tg, k)
            t["%s_%s" % (tg, k)] = df_all[col].to_numpy() if col in df_all.columns else np.nan
        t["%s_noflow" % tg] = ~t["%s_made_flow" % tg].to_numpy(dtype=bool)
    t["pass_19Aug"] = df_all["pass_19Aug"].to_numpy()
    return t.reset_index(drop=True)


def rescale_storms(t):
    """Min-max each big storm's KGE_2012 across the runs that made flow (worst 0, best 1); runs that made no flow
    (or whose KGE is undefined) get 0. Returns (list of per-storm info dicts, list of labels that are in the mean)."""
    infos, included = [], []
    for lab in BIG:
        tg = tag_of(lab)
        raw = t[tg + "_kge_2012"].to_numpy(dtype=float)
        mf = t[tg + "_made_flow"].to_numpy(dtype=bool)
        ok = mf & np.isfinite(raw)
        info = {"storm": lab, "n_runs": int(len(t)), "n_flow": int(mf.sum()), "n_scored": int(ok.sum()),
                "n_noflow": int((~mf).sum()), "best_raw": np.nan, "median_raw": np.nan, "worst_raw": np.nan,
                "status": "", "flag": ""}
        norm = np.full(len(t), np.nan)
        if ok.sum():
            info["best_raw"], info["worst_raw"] = float(raw[ok].max()), float(raw[ok].min())
            info["median_raw"] = float(np.median(raw[ok]))
        if ok.sum() < MIN_FLOW_RUNS:
            info["status"] = "LEFT OUT: fewer than %d runs make flow in this storm" % MIN_FLOW_RUNS
        elif info["best_raw"] - info["worst_raw"] < 1e-9:
            info["status"] = "LEFT OUT: all scored runs have the same KGE"
        else:
            lo, hi = info["worst_raw"], info["best_raw"]
            norm = np.where(ok, (raw - lo) / (hi - lo), 0.0)
            info["status"] = "in the mean"
            included.append(lab)
            if info["best_raw"] < WEAK_KGE:
                info["flag"] = "WEAK: best raw KGE_2012 is under %g" % WEAK_KGE
        t[tg + "_norm"] = norm
        infos.append(info)
    return infos, included


def make_versions(included):
    v = {}
    if "12 Aug" in included:
        v["ONLY 12 Aug"] = ["12 Aug"]
    # A combined score needs at least two storms. With exactly two, leaving one out leaves the other alone, so each
    # leave-one-out version is then a single storm (one of them is the same as ONLY 12 Aug; it is kept so that
    # "every leave-one-out version" always means every one).
    if len(included) >= 2:
        v["ALL"] = list(included)
        for s in included:
            v["LOO " + s] = [x for x in included if x != s]
    return v


def add_combined(t, versions):
    for key, stl in versions.items():
        cols = [tag_of(s) + "_norm" for s in stl]
        t["combined_" + key.replace(" ", "")] = t[cols].mean(axis=1).to_numpy()
    if "ALL" in versions:
        cols = [tag_of(s) + "_norm" for s in versions["ALL"]]
        t["worst_storm_norm"] = t[cols].min(axis=1).to_numpy()


def ckey(key):
    return "combined_" + key.replace(" ", "")


def top_positions(t, key, n_top):
    c = t[ckey(key)]
    order = t.assign(_c=c).sort_values(["_c", "run_id"], ascending=[False, True]).index.to_numpy()
    return order[:n_top]


# ------------------------------------------------------------------
# The reading rule
# ------------------------------------------------------------------
def read_parameter(u_all, comb_all, top_pos):
    """u_all: position of every run in the box (0-1); comb_all: combined score of every run; top_pos: positions of the
    best runs. Returns the label and the numbers behind it."""
    u_top = u_all[top_pos]
    p10, p50, p90 = [float(x) for x in np.percentile(u_top, [10, 50, 90])]
    span = p90 - p10
    narrowed = bool(span < NARROW_SPAN)
    edge = "LO" if p50 <= EDGE_BAND else ("HI" if p50 >= 1.0 - EDGE_BAND else "")
    k_idx = np.clip((u_all * N_SLICES).astype(int), 0, N_SLICES - 1)
    means = np.full(N_SLICES, np.nan)
    counts = np.zeros(N_SLICES, dtype=int)
    for k in range(N_SLICES):
        sel = k_idx == k
        counts[k] = int(sel.sum())
        if counts[k] >= MIN_SLICE_RUNS:
            means[k] = float(np.mean(comb_all[sel]))
    bracketed, peak = False, -1
    if np.isfinite(means).any():
        peak = int(np.nanargmax(means))
        if (1 <= peak <= N_SLICES - 2 and np.isfinite(means[0]) and np.isfinite(means[-1])
                and means[peak] - means[0] >= BRACKET_GAP and means[peak] - means[-1] >= BRACKET_GAP):
            bracketed = True
    if edge:
        label = "EDGE-" + edge
    elif narrowed and bracketed:
        label = "PINNED"
    elif narrowed:
        label = "NARROW"
    elif bracketed:
        label = "WIDE-PEAK"
    else:
        label = "FLAT"
    return {"label": label, "span": span, "p10": p10, "p50": p50, "p90": p90, "narrowed": narrowed, "edge": edge,
            "bracketed": bracketed, "peak_slice": peak, "slice_means": means, "slice_counts": counts}


def analyze_arm(t, versions, boxes, n_top):
    """readings[version][param_key] -> dict from read_parameter (plus values), tops[version] -> positions."""
    readings, tops = {}, {}
    for key in versions:
        top = top_positions(t, key, n_top)
        tops[key] = top
        comb = t[ckey(key)].to_numpy(dtype=float)
        rd = {}
        for pk, nm, sc in PARAMS:
            lo, hi, scale = boxes[pk]
            u = unit_pos(t[pk].to_numpy(dtype=float), lo, hi, scale)
            r = read_parameter(u, comb, top)
            for q in ("p10", "p50", "p90"):
                r["v" + q[1:]] = float(from_unit(r[q], lo, hi, scale))
            rd[pk] = r
        readings[key] = rd
    return readings, tops


def clear_verdicts(readings):
    """CLEAR only if PINNED in ALL and in every leave-one-out version."""
    out = {}
    loo = [k for k in readings if k.startswith("LOO ")]
    for pk, nm, sc in PARAMS:
        if "ALL" not in readings:
            out[pk] = (False, "no combined score (fewer than two storms could be rescaled)")
            continue
        labs = {"ALL": readings["ALL"][pk]["label"]}
        for k in loo:
            labs[k] = readings[k][pk]["label"]
        if all(v == "PINNED" for v in labs.values()):
            out[pk] = (True, "PINNED in all versions" + ("" if loo else " (no leave-one-out: fewer than three storms in the mean)"))
        else:
            edge = [k for k, v in labs.items() if v.startswith("EDGE")]
            if edge:
                why = "at a box edge in %s: the box was too small on that side, no claim" % ", ".join(edge)
            else:
                why = "not PINNED in " + ", ".join("%s (%s)" % (k, v) for k, v in labs.items() if v != "PINNED")
            out[pk] = (False, why)
    return out


def pair_rhos(t, top, boxes):
    u = {}
    for pk, nm, sc in PARAMS:
        lo, hi, scale = boxes[pk]
        u[pk] = unit_pos(t[pk].to_numpy(dtype=float)[top], lo, hi, scale)
    out = []
    for i in range(len(PARAMS)):
        for j in range(i + 1, len(PARAMS)):
            a, b = PARAMS[i][0], PARAMS[j][0]
            if np.ptp(u[a]) > 0 and np.ptp(u[b]) > 0:
                rho = float(stats.spearmanr(u[a], u[b])[0])
            else:
                rho = np.nan
            out.append((PNAME[a], PNAME[b], rho))
    return out


def passers_view(t, top, boxes):
    p = t["pass_19Aug"].to_numpy(dtype=bool)
    n_pass = int(p.sum())
    in_top = int(p[top].sum())
    chance = len(top) * n_pass / float(len(t))
    rng = {}
    if n_pass >= 5:
        for pk, nm, sc in PARAMS:
            lo, hi, scale = boxes[pk]
            u = unit_pos(t[pk].to_numpy(dtype=float)[p], lo, hi, scale)
            q = np.percentile(u, [10, 50, 90])
            rng[pk] = [float(from_unit(x, lo, hi, scale)) for x in q]
    return {"n_pass": n_pass, "in_top": in_top, "chance": chance, "ranges": rng}


# ------------------------------------------------------------------
# Printing
# ------------------------------------------------------------------
LABEL_HELP = [
    "PINNED     best runs share a narrow range (under half the box) AND runs on both sides of it score worse",
    "NARROW     best runs share a narrow range, but the box shows no clear peak (the range may be a slope, not a peak)",
    "WIDE-PEAK  the box shows an interior peak, but the best runs are spread over more than half of it",
    "FLAT       neither: this storm set does not pin the parameter",
    "EDGE-LO/HI the best runs pile against one end of the box: the box was too small there, so no claim is made",
]


def print_rules(args, n_runs, n_top):
    print("\nRULES IN FORCE (fixed before any ms1 result was looked at)")
    print("  storms in the mean: %s;  %s is read by pass/fail (flow of at least %g cfs, peak and volume within %.2f to %.0f times the real)"
          % (", ".join(BIG), SMALL, FLOW_CFS, PASS_LO, PASS_HI))
    print("  no flow (under %g cfs) in a storm that flowed = worst score (0); a storm is left out of the mean if fewer than %d runs make flow in it"
          % (FLOW_CFS, MIN_FLOW_RUNS))
    print("  each storm's KGE_2012 rescaled min-max over the runs (worst 0, best 1), then averaged; a storm whose best KGE_2012 is under %g is flagged"
          % WEAK_KGE)
    print("  best runs = top %d of %d (%.0f per cent).  NARROWED: middle 80 per cent of them cover under %g of the box.  EDGE: their median is within %g of an end."
          % (n_top, n_runs, 100 * TOP_FRAC, NARROW_SPAN, EDGE_BAND))
    print("  BRACKETED: the best of %d equal slices is an interior one and beats BOTH end slices by %g or more (slices with under %d runs do not count)."
          % (N_SLICES, BRACKET_GAP, MIN_SLICE_RUNS))
    print("  gauge arms: %s (primary) and %s; gap fill: line up to %g min, hold up to %g h"
          % (args.primary_arm, "asis" if args.primary_arm == "filled" else "filled", args.max_gap_min, args.max_hold_hr))


def print_storm_table(infos, storms, real):
    print("\nSTORMS (rescaling uses only runs that made flow; no-flow runs score 0)")
    print("  %-7s %-26s %9s %9s  %-15s %7s %7s %7s  %s" % ("storm", "window", "real peak", "real vol", "runs with flow", "best", "median", "worst", "status"))
    byl = {s["label"]: s for s in storms}
    for i in infos:
        s = byl[i["storm"]]
        tg = s["tag"]
        print("  %-7s %-26s %6.0f cfs %6.1f af  %4d of %-8d %7s %7s %7s  %s%s"
              % (i["storm"], "%s - %s" % (fmt_t(s["ws"]), fmt_t(s["we"])), real[tg]["real_peak_cfs"], real[tg]["real_vol_acft"],
                 i["n_flow"], i["n_runs"], fnum(i["best_raw"], ".3f"), fnum(i["median_raw"], ".3f"), fnum(i["worst_raw"], ".3f"),
                 i["status"], ("   *** " + i["flag"]) if i["flag"] else ""))
    print("  (best / median / worst are raw KGE_2012 among the runs that made flow; af = acre-feet; real volume = the filled gauge series summed over the window)")


def print_19aug(t, storms, real, view):
    tg = tag_of(SMALL)
    s = [x for x in storms if x["label"] == SMALL][0]
    mf = t[tg + "_made_flow"].to_numpy(dtype=bool)
    print("\n%s (real peak %.0f cfs, real volume %.2f acre-ft; read by pass/fail, not by KGE)"
          % (SMALL, real[tg]["real_peak_cfs"], real[tg]["real_vol_acft"]))
    print("  %d of %d runs made flow of at least %g cfs; %d PASS (peak and volume both within %.2f to %.0f times the real)"
          % (int(mf.sum()), len(t), FLOW_CFS, view["n_pass"], PASS_LO, PASS_HI))
    if mf.sum():
        print("  among runs that made flow: median peak ratio %.2f, median volume ratio %.2f"
              % (float(np.median(t[tg + "_peak_ratio"][mf])), float(np.median(t[tg + "_vol_ratio"][mf]))))


def print_reading(arm, t, versions, readings, tops, boxes, n_top, included, view, rhos, clear):
    n = len(t)
    print("\n" + "=" * 100)
    print("READING, %s ARM   (best %d of %d runs by combined score; storms in the mean: %s)"
          % (arm.upper(), n_top, n, ", ".join(included) if included else "none"))
    print("=" * 100)
    if "ALL" not in readings:
        print("  Fewer than two big storms could be rescaled: there is no combined score and no reading.")
        return
    top = tops["ALL"]
    best = t.iloc[top[0]]
    print("  best combined run: %s   score %.3f   Ks %s  f %s  cc %s"
          % (best["run_id"], best[ckey("ALL")], fval("Ks_mult", best["Ks_mult"]), fval("f_RS_abs", best["f_RS_abs"]),
             fval("channelconductivity_mmhr", best["channelconductivity_mmhr"])))
    if n_top > int((t[ckey("ALL")] > 0).sum()):
        print("  WARNING: fewer than %d runs score above 0, so the 'best runs' include zero scores; read nothing from this." % n_top)
    print("\n  ALL STORMS      %-9s  %-26s %9s  %6s   %s" % ("verdict", "middle 80% of best runs", "median", "span", "in every version?"))
    for pk, nm, sc in PARAMS:
        r = readings["ALL"][pk]
        ok, why = clear[pk]
        print("  %-14s  %-9s  %-26s %9s  %6.2f   %s" % (nm, r["label"], "%s to %s" % (fval(pk, r["v10"]), fval(pk, r["v90"])),
                                                       fval(pk, r["v50"]), r["span"], "CLEAR" if ok else "not clear: " + why))
    print("  (span = share of the box the middle 80% of the best runs cover; Ks linear, f and cc log scale)")
    print("\n  VERSIONS (verdict, with the span in brackets)")
    keys = [k for k in versions]
    print("  %-6s " % "" + " ".join("%-18s" % k for k in keys))
    for pk, nm, sc in PARAMS:
        print("  %-6s " % nm + " ".join("%-18s" % ("%s (%.2f)" % (readings[k][pk]["label"], readings[k][pk]["span"])) for k in keys))
    if "ONLY 12 Aug" in readings:
        print("\n  Compared with 12 Aug alone:")
        for pk, nm, sc in PARAMS:
            a, b = readings["ONLY 12 Aug"][pk], readings["ALL"][pk]
            ch = "narrower" if b["span"] < a["span"] - 0.05 else ("wider" if b["span"] > a["span"] + 0.05 else "about the same")
            print("    %-3s span %.2f (%s) -> %.2f (%s): %s with the other storms" % (nm, a["span"], a["label"], b["span"], b["label"], ch))
    print("\n  Parameter trade-offs among the best runs (Spearman rho; |rho| >= %g means they trade off, so read them as a pair):" % RIDGE_RHO)
    print("    " + "   ".join("%s-%s %s%s" % (a, b, fnum(r, "+.2f"), " (trade-off)" if np.isfinite(r) and abs(r) >= RIDGE_RHO else "")
                              for a, b, r in rhos))
    print("\n  19 Aug passers: %d runs pass; %d of them are among the best %d (chance alone would put about %.1f there)"
          % (view["n_pass"], view["in_top"], n_top, view["chance"]))
    if view["ranges"]:
        print("    where the passers sit (10th, median, 90th): " + ";  ".join(
            "%s %s / %s / %s" % (PNAME[pk], fval(pk, v[0]), fval(pk, v[1]), fval(pk, v[2])) for pk, v in view["ranges"].items()))
    elif view["n_pass"]:
        print("    (fewer than 5 passers: no ranges given)")
    print("\n  What the labels mean:")
    for line in LABEL_HELP:
        print("    " + line)


def print_arms(primary, other, ra, rb, ta, tb, tops_a, tops_b):
    print("\nARM AGREEMENT (primary %s against %s)" % (primary, other))
    if "ALL" not in ra or "ALL" not in rb:
        print("  not available: one arm has no combined score.")
        return None
    sa, sb = set(ta.iloc[tops_a["ALL"]]["run_id"]), set(tb.iloc[tops_b["ALL"]]["run_id"])
    share = len(sa & sb) / float(max(len(sa), 1))
    same = all(ra["ALL"][pk]["label"] == rb["ALL"][pk]["label"] for pk, _, _ in PARAMS)
    for pk, nm, sc in PARAMS:
        print("  %-3s %s: %-9s   %s: %-9s" % (nm, primary, ra["ALL"][pk]["label"], other, rb["ALL"][pk]["label"]))
    print("  best runs in common: %d of %d (%.0f per cent)" % (len(sa & sb), len(sa), 100 * share))
    agree = bool(same and share >= ARM_OVERLAP)
    print("  => the arms %s (agreement needs the same verdict for all three parameters and at least %.0f per cent of the best runs in common)"
          % ("AGREE" if agree else "DISAGREE: treat the primary arm's reading as unsettled", 100 * ARM_OVERLAP))
    return agree


CAVEATS = [
    "All four storms use the same two rain gauges, so the SMPHQ timing and intensity offset is in every one of them.",
    "One monsoon season (2014). Routing (cv, flowexp, n) is pinned. The runs are cc ON only: this shows where cc, Ks and f sit, not whether cc is needed.",
    "Gauge error is not quantified, so no uncertainty is stated. A measurement-error value from the literature is needed before any claim of that kind.",
    "min-max rescaling stretches a storm where nothing fits well; the WEAK flag and the leave-one-out versions are there for that.",
    "Marginal verdicts can understate a trade-off between parameters (a ridge); read the trade-off line.",
]


# ------------------------------------------------------------------
# Figures (matplotlib, light surface; one fixed hue per parameter; text in ink, never in the series colour)
# ------------------------------------------------------------------
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _style(ax):
    ax.set_facecolor(SURFACE)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8, length=3)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def _axis_for(ax, pk, boxes):
    """x limits with a little padding, and readable ticks (plain numbers, no minor-tick clutter on log axes)."""
    from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter
    lo, hi, scale = boxes[pk]
    if scale == "log":
        ax.set_xscale("log")
        pad = 0.03 * np.log10(hi / lo)
        ax.set_xlim(lo * 10 ** -pad, hi * 10 ** pad)
        ticks = [v for v in ([0.0003, 0.0005, 0.0008, 0.0014] if pk == "f_RS_abs" else [30, 100, 300, 1000]) if lo * 0.999 <= v <= hi * 1.001]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FixedFormatter([("%.4f" % v) if pk == "f_RS_abs" else ("%g" % v) for v in ticks]))
        ax.xaxis.set_minor_formatter(NullFormatter())
    else:
        pad = 0.03 * (hi - lo)
        ax.set_xlim(lo - pad, hi + pad)


XLAB = {"Ks_mult": "Ks multiplier", "f_RS_abs": "f (1/mm, log scale)", "channelconductivity_mmhr": "cc (mm/hr, log scale)"}


def fig_hydrographs(path, storms, real, t, tops, sim_cfs, arm, partial_note):
    plt = _plt()
    import matplotlib.dates as mdates
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6))
    fig.patch.set_facecolor(SURFACE)
    top = tops["ALL"]
    best_rid = t.iloc[top[0]]["run_id"]
    for ax, s in zip(axes.ravel(), storms):
        tg, x = s["tag"], s["grid"]
        _style(ax)
        for pos in top[1:]:
            ax.plot(x, sim_cfs[t.iloc[pos]["run_id"]][tg], color=C_KS, alpha=0.2, lw=0.8)
        ax.plot(x, sim_cfs[best_rid][tg], color=C_KS, lw=1.8)
        ax.plot(x, real[tg]["filled"].to_numpy() / CFS_TO_CMS, color=INK, lw=1.1)
        a = real[tg]["asis"]
        if a is not None:
            d = a.dropna()
            ax.scatter(d.index, d.to_numpy() / CFS_TO_CMS, s=9, color=INK, zorder=4)
        ax.set_title("%s   window %s to %s" % (s["label"], fmt_t(s["ws"]), fmt_t(s["we"])), fontsize=9, color=INK, loc="left")
        ax.set_ylabel("flow (cfs)", fontsize=8, color=INK2)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.set_ylim(bottom=0)
    handles = [Line2D([0], [0], color=INK, lw=1.1, label="real gauge (%s)" % ("filled arm" if arm == "filled" else "gaps filled")),
               Line2D([0], [0], color=INK, marker="o", lw=0, ms=4, label="real gauge bins that hold a reading"),
               Line2D([0], [0], color=C_KS, alpha=0.4, lw=1, label="the other best runs (all-storms score)"),
               Line2D([0], [0], color=C_KS, lw=1.8, label="best run (all-storms score)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8, frameon=False, labelcolor=INK2)
    fig.suptitle("Real gauge and the best runs, storm by storm%s" % partial_note, fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(path, dpi=140, facecolor=SURFACE)
    plt.close(fig)


def fig_param_scores(path, t, tops, boxes, included, arm, partial_note):
    plt = _plt()
    cols = [(lab, tag_of(lab) + "_norm") for lab in BIG] + [("all storms combined", ckey("ALL"))]
    fig, axes = plt.subplots(3, 4, figsize=(13, 9.4), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    top = set(int(i) for i in tops["ALL"])
    for r, (pk, nm, sc) in enumerate(PARAMS):
        lo, hi, scale = boxes[pk]
        xv = t[pk].to_numpy(dtype=float)
        u = unit_pos(xv, lo, hi, scale)
        for c, (title, col) in enumerate(cols):
            ax = axes[r, c]
            _style(ax)
            y = t[col].to_numpy(dtype=float)
            if not np.isfinite(y).any():
                ax.text(0.5, 0.5, "left out of the mean\n(see the storm table)", ha="center", va="center", transform=ax.transAxes,
                        fontsize=9, color=INK2)
            else:
                is_top = np.array([i in top for i in range(len(t))])
                ax.scatter(xv[~is_top], y[~is_top], s=14, color=PCOLOR[pk], alpha=0.45, lw=0)
                ax.scatter(xv[is_top], y[is_top], s=34, color=PCOLOR[pk], edgecolor=INK, lw=0.8, zorder=4)
                k_idx = np.clip((u * N_SLICES).astype(int), 0, N_SLICES - 1)
                cx, cy = [], []
                for k in range(N_SLICES):
                    sel = (k_idx == k) & np.isfinite(y)
                    if sel.sum() >= MIN_SLICE_RUNS:
                        cx.append(float(from_unit((k + 0.5) / N_SLICES, lo, hi, scale)))
                        cy.append(float(np.mean(y[sel])))
                ax.plot(cx, cy, color=INK, lw=1.4, marker="s", ms=3.5, zorder=3)
            _axis_for(ax, pk, boxes)
            ax.set_ylim(-0.04, 1.06)
            if r == 0:
                ax.set_title(title, fontsize=9.5, color=INK, loc="left")
            ax.set_xlabel(XLAB[pk], fontsize=8, color=INK2)
            if c == 0:
                ax.set_ylabel("%s\nrescaled score" % nm, fontsize=8.5, color=INK2)
    fig.suptitle("Rescaled storm scores against each parameter (0 = worst run, 1 = best). Ringed = the best %d runs by the all-storms score; line = mean of each of %d slices%s"
                 % (len(top), N_SLICES, partial_note), fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=1.6)
    fig.savefig(path, dpi=140, facecolor=SURFACE)
    plt.close(fig)


def fig_19aug(path, t, boxes, partial_note):
    plt = _plt()
    from matplotlib.lines import Line2D
    tg = tag_of(SMALL)
    mf = t[tg + "_made_flow"].to_numpy(dtype=bool)
    pas = t["pass_19Aug"].to_numpy(dtype=bool)
    vr = t[tg + "_vol_ratio"].to_numpy(dtype=float)
    pos = vr[mf & (vr > 0)]
    ymin = (pos.min() / 3.0) if pos.size else 0.01
    ymax = (pos.max() * 3.0) if pos.size else 10.0
    ymin, ymax = min(ymin, PASS_LO / 3.0), max(ymax, PASS_HI * 3.0)
    floor = ymin
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, (pk, nm, sc) in zip(axes, PARAMS):
        _style(ax)
        xv = t[pk].to_numpy(dtype=float)
        ax.axhspan(PASS_LO, PASS_HI, color=GRID, alpha=0.7, lw=0)
        ax.scatter(xv[~mf], np.full((~mf).sum(), floor * 1.4), marker="x", s=16, color=NULLC, lw=0.9)
        sel = mf & ~pas
        ax.scatter(xv[sel], np.maximum(vr[sel], floor * 1.4), s=22, facecolors="none", edgecolors=NULLC, lw=0.9)
        ax.scatter(xv[pas], vr[pas], s=30, color=C_KS, edgecolor=INK, lw=0.6, zorder=4)
        ax.set_yscale("log")
        ax.set_ylim(floor, ymax)
        _axis_for(ax, pk, boxes)
        ax.set_xlabel(XLAB[pk], fontsize=8, color=INK2)
    axes[0].set_ylabel("model volume / real volume (19 Aug)", fontsize=8.5, color=INK2)
    handles = [Line2D([0], [0], color=C_KS, marker="o", lw=0, ms=6, markeredgecolor=INK, label="PASS: flow, peak and volume within 1/3 to 3 times the real"),
               Line2D([0], [0], color=NULLC, marker="o", lw=0, ms=6, markerfacecolor="none", label="made flow but fails a ratio"),
               Line2D([0], [0], color=NULLC, marker="x", lw=0, ms=6, label="no flow (under %g cfs; drawn at the bottom)" % FLOW_CFS)]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8, frameon=False, labelcolor=INK2)
    fig.suptitle("19 Aug: which runs reproduce the small storm%s   (grey band = volume ratio 1/3 to 3)" % partial_note,
                 fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    fig.savefig(path, dpi=140, facecolor=SURFACE)
    plt.close(fig)


def fig_versions(path, readings, versions, boxes, partial_note):
    plt = _plt()
    keys = list(versions.keys())
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8.0))
    fig.patch.set_facecolor(SURFACE)
    for ax, (pk, nm, sc) in zip(axes, PARAMS):
        _style(ax)
        ax.axvspan(0, EDGE_BAND, color=GRID, alpha=0.8, lw=0)
        ax.axvspan(1 - EDGE_BAND, 1, color=GRID, alpha=0.8, lw=0)
        for i, k in enumerate(keys):
            r = readings[k][pk]
            y = len(keys) - 1 - i
            ax.plot([r["p10"], r["p90"]], [y, y], color=PCOLOR[pk], lw=6, solid_capstyle="round", alpha=0.85)
            ax.plot([r["p50"]], [y], marker="|", color=INK, ms=14, mew=2)
            ax.text(1.015, y, r["label"], va="center", fontsize=8.5, color=INK, transform=ax.get_yaxis_transform())
        ax.set_yticks(range(len(keys)))
        ax.set_yticklabels(list(reversed(keys)), fontsize=8.5, color=INK2)
        ax.set_ylim(-0.6, len(keys) - 0.4)
        ax.set_xlim(0, 1)
        lo, hi, scale = boxes[pk]
        ax.set_ylabel(nm, fontsize=10, color=INK, rotation=0, labelpad=22, va="center")
        ticks = [0, 0.25, 0.5, 0.75, 1.0]
        ax.set_xticks(ticks)
        ax.set_xticklabels([fval(pk, float(from_unit(u, lo, hi, scale))) for u in ticks], fontsize=8, color=INK2)
    axes[-1].set_xlabel("position in the sampled box, labelled with each parameter's own values (Ks linear, f and cc log); grey = edge bands", fontsize=8.5, color=INK2)
    fig.suptitle("Where the best runs sit, version by version (bar = middle 80 per cent, tick = median)%s" % partial_note,
                 fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 0.93, 0.96), h_pad=1.4)
    fig.savefig(path, dpi=140, facecolor=SURFACE)
    plt.close(fig)


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110: storm-by-storm scoring of the long ms1 run against the real gauge (no tRIBS).")
    ap.add_argument("--label", default="ms1", help="the --label of the long sweep (default ms1)")
    ap.add_argument("--primary_arm", choices=["filled", "asis"], default="filled", help="the arm whose reading is the verdict (default filled)")
    ap.add_argument("--gauge_mode", choices=["both", "filled", "asis"], default="both")
    ap.add_argument("--max_gap_min", type=float, default=DEFAULT_MAX_GAP_MIN)
    ap.add_argument("--max_hold_hr", type=float, default=DEFAULT_MAX_HOLD_HR)
    ap.add_argument("--top_frac", type=float, default=TOP_FRAC, help="share of the runs counted as the best runs (default %g)" % TOP_FRAC)
    ap.add_argument("--check_only", action="store_true", help="run the gates and stop; writes nothing")
    ap.add_argument("--allow_partial", action="store_true", help="accept a sweep that has not finished (labelled PARTIAL)")
    ap.add_argument("--drop_bad_runs", action="store_true", help="exclude runs that fail a gate and continue")
    ap.add_argument("--skip_gate2", action="store_true", help="skip the stored-metric check (not advised)")
    ap.add_argument("--no_plots", action="store_true")
    ap.add_argument("--windows", type=Path, default=None)
    ap.add_argument("--gauge_xlsx", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    args = ap.parse_args()
    if not (0.02 <= args.top_frac <= 0.5):
        ap.error("--top_frac must be between 0.02 and 0.5.")
    if args.gauge_mode != "both":
        args.primary_arm = args.gauge_mode
    arms = [args.primary_arm] + ([("asis" if args.primary_arm == "filled" else "filled")] if args.gauge_mode == "both" else [])

    script_dir = Path.cwd()
    project_root = script_dir.parent
    calib_dir = project_root / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    csv_dir = calib_dir / "03_comparisons" / "csv_exports"
    synth_dir = calib_dir / "synth_truth"
    keep_dir = calib_dir / "kept_qout_110"
    win_path = args.windows or (script_dir / WINDOWS_NAME)
    xlsx_path = args.gauge_xlsx or project_root.joinpath(*GAUGE_REL)
    out_dir = args.out_dir or (summary_dir / ("%s%s_110" % (OUT_PREFIX, args.label)))

    print("\n" + "=" * 100)
    print("Series 110 -- storm-by-storm scoring of the long run, real gauge (version %s)" % VERSION)
    print("label %s   primary arm %s   python %s  pandas %s  numpy %s" % (args.label, args.primary_arm, sys.version.split()[0], pd.__version__, np.__version__))
    print("=" * 100)
    if ".".join(pd.__version__.split(".")[:2]) != EXPECT_PANDAS_MINOR:
        print("  *** NOTE: pandas here is %s; the runs were scored with pandas %s.x. Gate 2 will tell. ***" % (pd.__version__, EXPECT_PANDAS_MINOR))

    out_res = out_dir.resolve()
    for guarded in (synth_dir, csv_dir, keep_dir):
        g = guarded.resolve()
        if out_res == g or g in out_res.parents or out_res in g.parents:
            sys.exit("Output folder %s overlaps an input folder (%s); refusing." % (out_dir, guarded))

    ck = Checks()
    storms, runtime = load_windows(win_path, ck)
    design = read_design(summary_dir / ("lhs_design_lowf_%s_ON_110.json" % args.label), args.label, ck)
    if storms and design:
        if design["runtime_hours"] is None or int(design["runtime_hours"]) != runtime:
            ck.add("FAIL", "run length", "the design says %s h, the windows file says %d h" % (design["runtime_hours"], runtime))
        else:
            ck.add("PASS", "run length", "design and windows file both say %d h" % runtime)
    res, partial = (None, False)
    if design:
        res, partial = load_results(summary_dir / ("lhs_results_lowf_%s_ON_110.csv" % args.label), design, ck, args.allow_partial)
    if not xlsx_path.exists():
        ck.add("FAIL", "gauge workbook", "not found: %s (--gauge_xlsx points elsewhere)" % xlsx_path)
    if not keep_dir.exists():
        ck.add("FAIL", "kept .qout folder", "not found: %s" % keep_dir)
    if ck.failures():
        print("\nGATES")
        ck.show()
        sys.exit("\nSTOPPED: nothing was scored. Fix the cause above.")
    boxes = design["box"]
    partial_note = "   [PARTIAL: %d of %d runs]" % (len(res), design["n"]) if partial else ""

    # ---- the real gauge, cut into the four windows ---------------------------------------------
    print("\nReading the real gauge ...")
    raw = read_gauge(xlsx_path)
    t_end = ORIGIN + pd.Timedelta(hours=runtime)
    in_run = raw[(raw.index >= ORIGIN) & (raw.index <= t_end)]
    covered = pd.Series(False, index=in_run.index)
    for s in storms:
        covered |= (in_run.index >= s["ws"]) & (in_run.index <= s["we"])
    outside = in_run[(in_run > 0) & ~covered]
    if len(outside):
        ck.add("FAIL", "real flow outside the windows", "%d reading(s) above 0 cfs fall outside all four windows (first: %s, %.0f cfs)"
               % (len(outside), outside.index[0], outside.iloc[0] / CFS_TO_CMS))
    else:
        ck.add("PASS", "real flow outside the windows", "no reading above 0 cfs lies outside the four windows")
    real = {}
    for s in storms:
        asis, filled, ftype = build_real_window(raw, s["ws"], s["we"], args.max_gap_min, args.max_hold_hr)
        n_asis, n_unc = int(asis.notna().sum()), int((ftype == "uncovered").sum())
        counts = {k: int((ftype == k).sum()) for k in ("record", "linear", "hold", "uncovered")}
        real[s["tag"]] = {"asis": asis if n_asis >= MIN_ASIS_POINTS else None, "asis_full": asis, "filled": filled, "ftype": ftype,
                          "real_peak_cfs": float(np.nanmax(filled.to_numpy())) / CFS_TO_CMS,
                          "real_vol_acft": float(np.nansum(filled.to_numpy())) * GRID_SEC / M3_PER_ACFT}
        ck.add("WARN" if (n_asis < MIN_ASIS_POINTS or n_unc) else "NOTE", "gauge, " + s["label"],
               "%d bins: %d hold a reading, %d line, %d hold, %d uncovered%s" % (len(filled), counts["record"], counts["linear"], counts["hold"], n_unc,
               "" if n_asis >= MIN_ASIS_POINTS else "; only %d as-is bins (under %d): the as-is arm cannot score this storm" % (n_asis, MIN_ASIS_POINTS)))
    stage0 = summary_dir / "real_gauge_110" / "real_gauge_observed_110.csv"
    if stage0.exists() and (storms[0]["ws"], storms[0]["we"]) == (STAGE_B_START, STAGE_B_END) \
            and args.max_gap_min == DEFAULT_MAX_GAP_MIN and args.max_hold_hr == DEFAULT_MAX_HOLD_HR:
        try:
            o = pd.read_csv(stage0, parse_dates=["datetime"]).set_index("datetime")
            same_idx = o.index.equals(real["12Aug"]["filled"].index)
            a_ok = np.allclose(o["Observed_asis_m3s"].to_numpy(float), real["12Aug"]["asis_full"].to_numpy(float), rtol=1e-9, atol=1e-12, equal_nan=True)
            f_ok = np.allclose(o["Observed_filled_m3s"].to_numpy(float), real["12Aug"]["filled"].to_numpy(float), rtol=1e-9, atol=1e-12, equal_nan=True)
            ck.add("PASS" if (same_idx and a_ok and f_ok) else "FAIL", "GATE 1",
                   "the 12 Aug real series (as-is and filled) built here %s %s from Stage 0" %
                   ("equals" if (same_idx and a_ok and f_ok) else "DIFFERS from", stage0.name))
        except Exception as exc:
            ck.add("WARN", "GATE 1", "could not compare with %s: %s" % (stage0.name, str(exc)[:80]))
    else:
        ck.add("NOTE", "GATE 1", "%s not found (or non-default gap settings): the 12 Aug real series was not compared with Stage 0" % stage0.name)

    synth5 = None
    if not args.skip_gate2:
        qs = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
        if len(qs) != 1:
            ck.add("FAIL", "GATE 2", "expected exactly one *.qout in %s (the synthetic truth the runs were scored against), found %d"
                   % (synth_dir, len(qs)))
        else:
            synth5 = read_qout_5min(qs[0])
            ck.add("NOTE", "GATE 2", "synthetic truth (read-only): %s" % qs[0].name)
    else:
        ck.add("WARN", "GATE 2", "skipped by --skip_gate2")
    if ck.failures():
        print("\nGATES")
        ck.show()
        sys.exit("\nSTOPPED: nothing was scored. Fix the cause above.")

    print("\nScoring %d runs (reading each kept .qout once) ..." % len(res))
    recs, sim_cfs, bad = score_runs(res, keep_dir, storms, real, runtime, synth5, ck, args.drop_bad_runs)
    print("\nGATES")
    ck.show()
    if bad:
        print("\n  %d run(s) failed a gate or could not be read. First few:" % len(bad))
        for rid, why in list(bad.items())[:8]:
            print("    %s: %s" % (rid, why))
        if not args.drop_bad_runs:
            sys.exit("\nSTOPPED: nothing was scored. Fix the cause (or add --drop_bad_runs to leave these runs out and continue).")
        print("  --drop_bad_runs set: these runs are excluded.")
    if ck.failures():
        sys.exit("\nSTOPPED: a gate failed. Nothing was scored.")
    if len(recs) < 20:
        sys.exit("\nSTOPPED: only %d usable runs; too few to read anything." % len(recs))
    if args.check_only:
        print("\n--check_only: gates finished; %d runs readable; nothing was written." % len(recs))
        return

    # ---- scoring tables, rescaling, versions, reading ---------------------------------------------
    df_all = pd.DataFrame(recs)
    n_runs = len(df_all)
    n_top = int(min(max(TOP_MIN, int(round(args.top_frac * n_runs))), n_runs))
    print_rules(args, n_runs, n_top)
    per_arm = {}
    for arm in arms:
        t = arm_table(df_all, arm, storms)
        infos, included = rescale_storms(t)
        versions = make_versions(included)
        add_combined(t, versions)
        readings, tops = analyze_arm(t, versions, boxes, n_top) if versions else ({}, {})
        per_arm[arm] = {"t": t, "infos": infos, "included": included, "versions": versions, "readings": readings, "tops": tops}

    P = per_arm[args.primary_arm]
    t = P["t"]
    print("\n" + "=" * 100)
    print("RESULTS, %s ARM (primary)%s" % (args.primary_arm.upper(), partial_note))
    print("=" * 100)
    print_storm_table(P["infos"], storms, real)
    tops_all = P["tops"].get("ALL")
    top_for_view = tops_all if tops_all is not None else np.array([], dtype=int)
    view = passers_view(t, top_for_view, boxes) if len(top_for_view) else passers_view(t, np.array([0]), boxes)
    print_19aug(t, storms, real, view)
    rhos = pair_rhos(t, tops_all, boxes) if tops_all is not None else []
    clear = clear_verdicts(P["readings"])
    print_reading(args.primary_arm, t, P["versions"], P["readings"], P["tops"], boxes, n_top, P["included"], view, rhos, clear)
    agree = None
    if len(arms) == 2:
        Q = per_arm[arms[1]]
        agree = print_arms(arms[0], arms[1], P["readings"], Q["readings"], t, Q["t"], P["tops"], Q["tops"])
        if [i["status"] for i in Q["infos"]] != [i["status"] for i in P["infos"]]:
            print("  note: the arms differ in which storms could be rescaled: %s" % "; ".join(
                "%s %s: %s" % (arms[1], i["storm"], i["status"]) for i in Q["infos"]))
    print("\nCAVEATS")
    for c in CAVEATS:
        print("  - " + c)

    # ---- files ---------------------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    for arm in arms:
        tt = per_arm[arm]["t"].copy()
        write_csv_atomic(tt, out_dir / ("storm_scores_%s_110.csv" % arm))
    rng_rows = []
    for arm in arms:
        for i in per_arm[arm]["infos"]:
            rng_rows.append(dict(arm=arm, **i))
    write_csv_atomic(pd.DataFrame(rng_rows), out_dir / "storm_ranges_110.csv")
    rd_rows, top_rows = [], []
    for arm in arms:
        A = per_arm[arm]
        for key in A["versions"]:
            for pk, nm, sc in PARAMS:
                r = A["readings"][key][pk]
                rd_rows.append({"arm": arm, "version": key, "storms": "+".join(A["versions"][key]), "parameter": nm, "verdict": r["label"],
                                "span": r["span"], "median": r["v50"], "value_p10": r["v10"], "value_p90": r["v90"],
                                "narrowed": r["narrowed"], "edge": r["edge"], "bracketed": r["bracketed"], "peak_slice": r["peak_slice"],
                                **{"slice%d_mean" % (k + 1): r["slice_means"][k] for k in range(N_SLICES)},
                                **{"slice%d_runs" % (k + 1): int(r["slice_counts"][k]) for k in range(N_SLICES)}})
            tt = A["t"]
            for rank, pos in enumerate(A["tops"][key], 1):
                row = tt.iloc[pos]
                top_rows.append({"arm": arm, "version": key, "rank": rank, "run_id": row["run_id"], "combined": row[ckey(key)],
                                 "Ks_mult": row["Ks_mult"], "f_RS_abs": row["f_RS_abs"], "channelconductivity_mmhr": row["channelconductivity_mmhr"],
                                 **{tag_of(b) + "_norm": row[tag_of(b) + "_norm"] for b in BIG}})
    write_csv_atomic(pd.DataFrame(rd_rows), out_dir / "storm_readings_110.csv")
    write_csv_atomic(pd.DataFrame(top_rows), out_dir / "storm_top_runs_110.csv")
    made = ["storm_scores_%s_110.csv" % a for a in arms] + ["storm_ranges_110.csv", "storm_readings_110.csv", "storm_top_runs_110.csv"]

    if not args.no_plots:
        # three of the figures show the combined reading, so they are skipped when there is none
        has_all = "ALL" in P["readings"]
        figs = [("fig_storm_hydrographs_110.png", lambda p: fig_hydrographs(p, storms, real, t, P["tops"], sim_cfs, args.primary_arm, partial_note), True),
                ("fig_storm_param_scores_110.png", lambda p: fig_param_scores(p, t, P["tops"], boxes, P["included"], args.primary_arm, partial_note), True),
                ("fig_storm_19aug_110.png", lambda p: fig_19aug(p, t, boxes, partial_note), False),
                ("fig_storm_versions_110.png", lambda p: fig_versions(p, P["readings"], P["versions"], boxes, partial_note), True)]
        if not has_all:
            print("\n  (no combined score, so the hydrograph, parameter-score and versions figures are not drawn)")
        for name, fn, needs_all in figs:
            if needs_all and not has_all:
                continue
            try:
                fn(out_dir / name)
                made.append(name)
            except Exception as exc:
                print("\n  (figure %s skipped: %s -- the tables are saved regardless)" % (name, exc))

    prov = {"script": Path(__file__).name, "version": VERSION, "script_md5": md5_of(Path(__file__)),
            "created_local": datetime.now().isoformat(timespec="seconds"), "python": sys.version.split()[0],
            "pandas": pd.__version__, "numpy": np.__version__,
            "label": args.label, "partial": bool(partial), "runs_scored": int(n_runs), "runs_excluded": bad,
            "inputs": {"windows": {"file": win_path.name, "md5": md5_of(win_path)},
                       "design": {"file": design["path"].name, "md5": md5_of(design["path"])},
                       "results": {"file": "lhs_results_lowf_%s_ON_110.csv" % args.label,
                                   "md5": md5_of(summary_dir / ("lhs_results_lowf_%s_ON_110.csv" % args.label))},
                       "gauge": {"file": xlsx_path.name, "md5": md5_of(xlsx_path)}},
            "gates": [{"status": s, "name": n, "text": x} for s, n, x in ck.rows],
            "rule": {"flow_cfs": FLOW_CFS, "pass_ratio": [PASS_LO, PASS_HI], "weak_kge": WEAK_KGE, "min_flow_runs": MIN_FLOW_RUNS,
                     "top_frac": args.top_frac, "n_top": n_top, "narrow_span": NARROW_SPAN, "edge_band": EDGE_BAND,
                     "n_slices": N_SLICES, "bracket_gap": BRACKET_GAP, "min_slice_runs": MIN_SLICE_RUNS,
                     "gauge_gap_min": args.max_gap_min, "gauge_hold_hr": args.max_hold_hr},
            "storms": [{"label": s["label"], "window_start": str(s["ws"]), "window_end": str(s["we"]),
                        "real_peak_cfs": real[s["tag"]]["real_peak_cfs"], "real_volume_acft": real[s["tag"]]["real_vol_acft"]} for s in storms],
            "primary_arm": args.primary_arm, "arms_agree": agree,
            "clear": {PNAME[pk]: {"clear": bool(v[0]), "why": v[1]} for pk, v in clear.items()}}
    write_json_atomic(out_dir / "PROVENANCE_storm_scores_110.json", prov)
    made.append("PROVENANCE_storm_scores_110.json")
    print("\nSaved to: %s" % out_dir)
    print("  " + "   ".join(made))


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        main()
