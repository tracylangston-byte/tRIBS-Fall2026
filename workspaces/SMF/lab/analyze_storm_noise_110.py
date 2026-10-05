"""
analyze_storm_noise_110.py
==========================
Series 110 -- THREE STORMS + OBSERVATION ERROR: does "the best runs follow the true cc"
survive an imperfect observed hydrograph, in each storm and when the storms are used together?

THE QUESTION
------------
analyze_cc_recovery_noise_110.py asked this for ONE storm with Ks and f pinned at the truth.
analyze_storm_compare_110.py then repeated the noise-free analysis for three storms (rain x0.8,
x1.0, x1.25) with 27 truths per storm (3 Ks x 3 f x 3 cc), where Ks and f are also unknown.
This script puts the two together: it adds observation error to the "observed" hydrograph of
every truth in every storm and measures how well the best runs still follow the true cc.

Pure post-hoc analysis. Runs no tRIBS. It re-reads the stored hydrograph of every run
(<run_id>_compare_obs_sim.csv) and the truth files, scores every run against a noisy copy of
every truth, picks the best runs, and measures the tracking slope again. It writes only into its
own new folder (calibration_work/03_comparisons/summary_tables/storm_noise_110/).

Keep this file in the same folder as analyze_storm_compare_110.py, analyze_truth_location_110.py,
analyze_cc_recovery_110.py, analyze_cc_tests_110.py, analyze_cc_pca_110.py and
rescore_cc_truths_110.py (it reuses their code). Run it from the lab/ directory.

WHAT IT DOES
------------
For each noise level (none, then 5%, 10%, 20% by default) and for each of 100 noisy replicates:
  1. For every truth in every storm, build a noisy "observed" hydrograph (noise model below).
     Every (storm, truth) gets its own independent noise draw.
  2. Score all 250 cc-ON runs and all 250 control runs (channel loss OFF) against that noisy
     hydrograph with KGE_2012, exactly as the main analysis does. The control sees the same noisy
     hydrograph as the cc-ON runs.
  3. For each storm: take the best 20% / 10% / 5% of runs per truth and fit the pooled tracking
     slope of (median best-run log10 cc) against (true log10 cc) over all 27 truths.
  4. "Combined": average each run's noisy KGE_2012 over the three storms (the same parameter
     point is a run in every storm), then pick the best runs on that average. Two versions:
       combined            each storm has its own independent whole-hydrograph error b
       combined (shared)   the SAME b is used in all three storms for a given truth. That is
                           how a rating-curve error at one gauge would behave, so this is the
                           honest test of whether combining storms still helps. (Skipped with
                           --noise_parts shape, where there is no b.)
  5. Shuffle cc among the runs (1000 shuffles, one shuffle shared by all truths, storms and
     versions) to get the p-value for "could a fit that ignores cc give a slope this large?".
Then it summarizes over the replicates: the typical slope, how much it varies with the luck of
the noise draw, the share of replicates in which it is still significant, ON minus control, and
how many of the 9 (Ks, f) cells still show cc tracking.

NOISE MODEL (the same as the single-storm noise test)
-----------------------------------------------------
    observed = truth x max(1 + b, 0.05) x exp(e)
    b ~ Normal(0, sd)        one error for the whole hydrograph (a volume or rating-curve error)
    e ~ Normal(0, sd / 2)    independent at every 5-minute step (measurement scatter)
sd is 5%, 10% and 20% by default. --noise_parts volume uses only b; --noise_parts shape only e.
The draws depend only on (seed, noise level, replicate, storm, truth), so selecting different
--storms or --noise_levels does not change the numbers of the ones you keep.

VOCABULARY
-----------
  Storm            "x0.8" = the original 12 Aug 2014 rain times 0.8; "x1.0" = the original;
                   "x1.25" = times 1.25.
  Replicate        One noisy "observed" data set (for the combined case, one noisy data set per
                   storm). 100 replicates show how much the result changes with the luck of the
                   noise draw.
  Middle 80% of replicates   The range holding the central 80% of the replicate values (10th to
                   90th percentile). It measures the effect of noise, not of the choice of runs.
  Tracking slope   Slope of (median best-run log10 cc) against (true log10 cc), over all 27
                   truths. 1 = follows the truth one-for-one; 0 = ignores it. The edges of the
                   sampled range squash it below 1 even for a perfect method.
  Significant      The share of replicates in which shuffling cc almost never (p < 0.05) gives a
                   slope this large. With no noise it is a single yes or no.
  Typical error    The factor by which the best runs' median cc typically misses the true cc
                   (x1.00 = perfect), next to what a no-information fit gives (in brackets).
  Control          The same 250 parameter points with channel loss OFF, so cc is inert. Its
                   slope must stay near 0 at every noise level; if not, noise alone creates a
                   pattern.
  ON minus control The ON slope minus the control slope in the same replicate (both see the same
                   noisy hydrograph): how much of the tracking is more than chance.
  Cell             One (Ks, f) combination, with its three truths at three cc values. A cell
                   "tracks" in a replicate when its cc slope is above 0 and shuffling cc rarely
                   (p < 0.05) gives a slope that large. This is a permutation rule, not the
                   bootstrap rule of analyze_storm_compare_110.py, so the counts are similar in
                   spirit but not identical.
  Combined         Each run's KGE_2012 averaged over the storms, then the best runs are picked.
  Shared b         In "combined (shared)", one volume error per truth and replicate applies to
                   all storms (a rating-curve error at one gauge). The scatter e stays
                   independent in each storm.

HOW TO READ THE RESULTS
------------------------
  * Slope stays near its no-noise value and the share significant stays high: the tracking
    survives that noise level.
  * Slope shrinks toward 0 as noise grows: the tracking is fragile. The noise level where its
    middle 80% first reaches 0 is the limit of what these data can tell you about cc.
  * Combined above the single storms: using several storms together helps. If "combined
    (shared)" is much lower than "combined", the gain depends on the storms' errors being
    independent, which a real gauge would not guarantee.
  * Slope stays positive but the share significant falls: the average tracking survives but a
    single real data set could easily fail to show it.
  * The numbers say what these storms, these 250 runs and these 27 truths can do. They do not
    make a real basin's gauge error equal to 5%, 10% or 20%.

OUTPUTS   (calibration_work/03_comparisons/summary_tables/storm_noise_110/)
---------------------------------------------------------------------------
  storm_noise_summary_110.csv          cut x entry x series x noise level: slope and spread over
                                       replicates, share significant, typical error, cells tracking
  storm_noise_on_minus_control_110.csv cut x entry x noise level: ON minus control slope
  storm_noise_per_truth_110.csv        cut x entry x series x noise x truth: median best-run cc
  storm_noise_cells_110.csv            cut x entry x series x noise x cell: cc slope in that cell
  storm_noise_realizations_110.csv     every replicate (slope, p-value, typical error, cells tracking)
  storm_noise_noisefree_check_110.csv  the noise-free re-computation check, per storm and series
  storm_noise_input_checks_110.csv
  fig_stormnoise_slope_110.png         slope against noise level, one panel per cut and entry
  fig_stormnoise_share_110.png         share of replicates still significant, one panel per cut
  fig_stormnoise_recovery_110.png      best-run cc against true cc, per noise level
  PROVENANCE_storm_noise_110.json
(With --noise_parts volume or shape the file names end in _volume or _shape.)

USAGE (run from the lab/ directory)
------------------------------------
Allow several minutes for a default run (about 4 minutes on a test machine; progress is printed).
Start with --check_only: it loads everything and runs all checks in about a minute and writes nothing.
    python -u analyze_storm_noise_110.py --check_only            # loads and checks, writes nothing
    python -u analyze_storm_noise_110.py 2>&1 | tee storm_noise_110.log
    python -u analyze_storm_noise_110.py --noise_reps 200 --perm 2000
    python -u analyze_storm_noise_110.py --noise_levels 0.02 0.05 0.10
    python -u analyze_storm_noise_110.py --noise_parts volume    # only the whole-hydrograph error
    python -u analyze_storm_noise_110.py --no_plots

CAVEATS
--------
The storms share the same 250 sweep points, routing and rainfall pattern (only scaled), so they
are not independent experiments. Replicates reuse the same 250 runs, so the control column is not
the 5% calibration of the permutation test. The noise is multiplicative and has no memory in time;
real gauge errors can be larger at high flows and correlated in time. The 27 truths share the same
runs, so p-values describe the runs, not 27 separate experiments.
"""

import argparse
import hashlib
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import analyze_storm_compare_110 as SC
except ImportError:
    sys.exit("This script needs analyze_storm_compare_110.py, analyze_truth_location_110.py, analyze_cc_recovery_110.py, "
             "analyze_cc_tests_110.py and analyze_cc_pca_110.py in the same folder (it reuses their code). "
             "Put all of them in lab/ and run from there.")
try:
    import rescore_cc_truths_110 as RS
except ImportError:
    sys.exit("This script needs rescore_cc_truths_110.py in the same folder (it re-reads the stored hydrographs and truths).")

L, R, T, P = SC.L, SC.R, SC.T, SC.P

_NEED = {
    "analyze_storm_compare_110": (SC, ["entry_paths", "load_long", "check_entry", "check_cross", "build_entry_data", "make_combined",
                                       "scale_name", "storm_label", "md5_file", "REF_SCALE", "REF_KEY", "COMB_KEY", "SER_ON", "SER_OFF"]),
    "analyze_truth_location_110": (L, ["truth_coords", "find_cells", "cell_name", "LOC_DIRNAME"]),
    "analyze_cc_recovery_110": (R, ["top_index", "slope_of", "LOG_LO", "LOG_HI", "CENTER_LOG"]),
    "analyze_cc_tests_110": (T, ["kge_matrix", "make_perm_matrix", "p_high", "p_low", "fmt_p", "_style", "n_top_for",
                                 "COL_CC", "COL_NULL", "INK", "INK2", "GRID", "SURFACE"]),
    "analyze_cc_pca_110": (P, ["Checks", "print_checks", "generate_lhs_samples", "LHS_PARAMS", "DESIGN_SEED", "PARAMS", "STAGE2_NAME",
                               "CONTROL_NAME"]),
    "rescore_cc_truths_110": (RS, ["load_all_compares", "read_truth_5min", "md5_of"]),
}
for _label, (_mod, _names) in _NEED.items():
    _miss = [n for n in _names if not hasattr(_mod, n)]
    if _miss:
        sys.exit("%s.py in this folder is an older version (missing: %s). Replace it with the current copy."
                 % (_label, ", ".join(_miss)))

SER_ON, SER_OFF = SC.SER_ON, SC.SER_OFF
REF_SCALE, REF_KEY, COMB_KEY = SC.REF_SCALE, SC.REF_KEY, SC.COMB_KEY
SHARED_KEY = "combined_shared"
OUT_DIRNAME = "storm_noise_110"
EARLIER_DIR = SC.OUT_DIRNAME
ALPHA = 0.05
NOISE_FREE_TOL = 1e-6
COL_ON, COL_NULL = T.COL_CC, T.COL_NULL
INK, INK2, GRID, SURFACE = T.INK, T.INK2, T.GRID, T.SURFACE
ENTRY_MARKERS = {"single0": "v", "single1": "o", "single2": "^", "single3": "P", "single4": "X", COMB_KEY: "D", SHARED_KEY: "s"}


def level_label(lv):
    return "none" if lv == 0 else "%g%%" % round(100 * lv, 6)


def cut_label(frac):
    return "%d%%" % round(100 * frac)


def fnum(v, spec):
    return "n/a" if (v is None or not np.isfinite(v)) else format(v, spec)


# ------------------------------------------------------------------
# Stored hydrographs and truth hydrographs
# ------------------------------------------------------------------
def truth_dir(calib, scale):
    if abs(scale - REF_SCALE) < 1e-9:
        return calib / L.LOC_DIRNAME
    return calib / ("synth_truth_location_%s" % SC.storm_label(scale))


def load_storm_hydrographs(ent, csv_dir):
    """Simulated hydrographs (runs x time) for each series of this storm, in the order of ent["data"][s]["run_id"].
    Returns (sims, time index)."""
    sims, idx0 = {}, None
    for s in ent["data"]:
        ids = ent["data"][s]["run_id"]
        cache, prob = RS.load_all_compares(pd.DataFrame({"run_id": ids}), csv_dir)
        if prob:
            sys.exit("STOPPED: %d compare CSV(s) missing or unreadable for %s (%s), e.g. %s\n  Looked in: %s\n  This script re-reads "
                     "every run's stored hydrograph (<run_id>_compare_obs_sim.csv). Nothing was written."
                     % (len(prob), ent["name"], s, list(prob.items())[0], csv_dir))
        mats = []
        for rid in ids:
            d = cache[rid]
            if idx0 is None:
                idx0 = d.index
            if not d.index.equals(idx0):
                sys.exit("STOPPED: compare CSV %s (%s, %s) has a different time index from the others of this storm; refusing."
                         % (rid, ent["name"], s))
            mats.append(d["Simulated"].to_numpy(float))
        sims[s] = np.vstack(mats)
    return sims, idx0


def load_truth_hydrographs(ent, calib, idx, truth_ids):
    """Each truth's hydrograph on this storm's compare-CSV time index, verified against the checksum recorded at re-scoring."""
    prov_path = ent["paths"]["prov"]
    if not prov_path.exists():
        sys.exit("STOPPED: %s is missing, so there is no record of which truth file belongs to which truth for %s. "
                 "Re-run the re-scoring script for this storm." % (prov_path.name, ent["name"]))
    prov = json.loads(prov_path.read_text())
    verified = prov.get("truths", {})
    tdir = truth_dir(calib, ent["scale"])
    out, md5s = {}, {}
    for tid in truth_ids:
        info = verified.get(tid)
        if info is None:
            sys.exit("STOPPED: truth %s is not in the verified list of %s (%s)." % (tid, prov_path.name, ent["name"]))
        p = tdir / info["file"]
        if not p.exists():
            sys.exit("STOPPED: truth file %s is missing from %s (%s)." % (info["file"], tdir, ent["name"]))
        m = RS.md5_of(p)
        if m != info.get("md5"):
            sys.exit("STOPPED: truth file %s no longer matches the checksum recorded at re-scoring (%s); refusing. NEVER move or "
                     "replace files in the truth folders." % (info["file"], ent["name"]))
        s5 = RS.read_truth_5min(p).reindex(idx)
        if s5.isna().any():
            sys.exit("STOPPED: truth file %s has no value at some timestamps of the compare CSVs (%s)." % (info["file"], ent["name"]))
        out[tid] = s5.to_numpy(float)
        md5s[tid] = m
    return out, md5s


def noise_free_check(ent, sims, truth_obs, truth_ids):
    """Recompute KGE_2012 of every (truth, run) without noise, with the function the noisy runs use, and compare with the stored
    values of the main analysis. Returns {series: dict}."""
    res = {}
    for s in ent["data"]:
        stored = ent["data"][s]["K"]
        worst, mism, npairs = 0.0, 0, 0
        for ti, tid in enumerate(truth_ids):
            k = T.kge_matrix(sims[s], truth_obs[tid])
            a, b = np.isfinite(k), np.isfinite(stored[ti])
            mism += int(np.sum(a != b))
            ok = a & b
            npairs += int(ok.sum())
            if ok.any():
                worst = max(worst, float(np.max(np.abs(k[ok] - stored[ti][ok]))))
        res[s] = {"worst_abs_diff": worst, "n_pairs": npairs, "n_finite_mismatch": mism}
    return res


# ------------------------------------------------------------------
# Noise
# ------------------------------------------------------------------
def truth_key_ints(tr):
    """Non-negative integers identifying a truth, for seeding."""
    return [int(round(tr["Ks"] * 1000)), int(round(tr["f"] * 1e6)), int(round(tr["cc"] * 1000))]


def draw_noise(lv, parts, seed, lcode, rep, scode, key, n_t):
    """(b_independent, b_shared, e) for one (level, replicate, storm, truth). b_independent and e depend on the storm;
    b_shared does not (the same value for all storms of this truth and replicate)."""
    rng = np.random.default_rng([seed, 88, lcode, rep, scode] + key)
    b_i = float(rng.normal(0.0, lv))
    e = rng.normal(0.0, lv / 2.0, n_t)
    b_s = float(np.random.default_rng([seed, 89, lcode, rep] + key).normal(0.0, lv))
    if parts == "shape":
        b_i = b_s = 0.0
    elif parts == "volume":
        e = np.zeros(n_t)
    return b_i, b_s, e


def make_obs(base, b, e):
    return base * max(1.0 + b, 0.05) * np.exp(e)


# ------------------------------------------------------------------
# One replicate
# ------------------------------------------------------------------
def rep_stats(Kd, ccs, xcc, perm, fracs, keys, cells):
    """Tracking statistics of one entry for one replicate. Kd[s] has shape (truths, runs); ccs[s] is the log10 cc of every run;
    xcc is the log10 true cc of every truth; perm[s] is a (shuffles, runs) matrix shared by all entries.
    Returns rows (one per series and cut) and, per (series, cut), the best-run median log10 cc per truth and the cell results."""
    rows, pts, cell_out = [], {}, {}
    for s, K in Kd.items():
        cc, Pm = ccs[s], perm[s]
        for frac in fracs:
            idxs = [R.top_index(K[ti], frac) for ti in range(K.shape[0])]
            pt = np.array([np.median(cc[ix]) for ix in idxs])
            if len({len(ix) for ix in idxs}) == 1:
                nl = np.median(cc[Pm[:, np.array(idxs)]], axis=2)
            else:
                nl = np.column_stack([np.median(cc[Pm[:, ix]], axis=1) for ix in idxs])
            sl = float(R.slope_of(xcc, pt))
            sl_null = R.slope_of(xcc, nl)
            err = float(np.mean(np.abs(pt - xcc)))
            err_null = np.mean(np.abs(nl - xcc), axis=1)
            csl, cp = np.empty(len(keys)), np.empty(len(keys))
            for ci, k in enumerate(keys):
                ix = cells[k]
                csl[ci] = float(R.slope_of(xcc[ix], pt[ix]))
                cp[ci] = float(T.p_high(R.slope_of(xcc[ix], nl[:, ix]), csl[ci]))
            ntrack = int(np.sum((csl > 0) & (cp < ALPHA)))
            rows.append({"series": s, "top_frac": frac, "slope": sl, "slope_null_mean": float(sl_null.mean()),
                         "p_slope": float(T.p_high(sl_null, sl)), "typical_error_factor": float(10 ** err),
                         "null_typical_error_factor": float(10 ** err_null.mean()),
                         "p_error_smaller": float(T.p_low(err_null, err)), "n_cells": len(keys), "n_cells_tracking": ntrack})
            pts[(s, frac)] = pt
            cell_out[(s, frac)] = (csl, cp)
    return rows, pts, cell_out


def run_replicates(args, singles, extra, truths, xcc, keys, cells, sims, truth_obs, perm, fracs, n_t_by_storm):
    """Loops over noise levels and replicates. Returns (realizations DataFrame, per-truth arrays, cell arrays)."""
    levels = [0.0] + list(args.noise_levels)
    need_shared = SHARED_KEY in [e["key"] for e in extra]
    comb = [e for e in extra if e["key"] == COMB_KEY]
    all_entries = list(singles) + list(extra)
    ccs = {e["key"]: {s: e["data"][s]["X"][:, 2] for s in e["data"]} for e in all_entries}
    keyints = [truth_key_ints(tr) for tr in truths]
    rz = []
    ptz, clz = {}, {}                                   # (level, entry, series, frac) -> list over replicates
    t0 = time.time()
    for lv in levels:
        lcode = int(round(lv * 100000))
        reps = 1 if lv == 0 else args.noise_reps
        for rep in range(reps):
            K_ind, K_sh = {}, {}
            for e in singles:
                k = e["key"]
                if lv == 0:
                    K_ind[k] = {s: e["data"][s]["K"] for s in e["data"]}
                    K_sh[k] = K_ind[k]
                    continue
                scode = int(round(e["scale"] * 100))
                K_ind[k] = {s: np.full(e["data"][s]["K"].shape, np.nan) for s in e["data"]}
                K_sh[k] = {s: np.full(e["data"][s]["K"].shape, np.nan) for s in e["data"]} if need_shared else None
                for ti, tr in enumerate(truths):
                    base = truth_obs[k][tr["id"]]
                    b_i, b_s, ez = draw_noise(lv, args.noise_parts, args.seed, lcode, rep, scode, keyints[ti], n_t_by_storm[k])
                    o_i = make_obs(base, b_i, ez)
                    o_s = make_obs(base, b_s, ez) if need_shared else None
                    for s in e["data"]:
                        K_ind[k][s][ti] = T.kge_matrix(sims[k][s], o_i)
                        if need_shared:
                            K_sh[k][s][ti] = T.kge_matrix(sims[k][s], o_s)
            Kd_by_entry = {e["key"]: K_ind[e["key"]] for e in singles}
            if comb:
                sers = list(comb[0]["data"])
                Kd_by_entry[COMB_KEY] = {s: np.mean([K_ind[e["key"]][s] for e in singles], axis=0) for s in sers}
                if need_shared:
                    Kd_by_entry[SHARED_KEY] = {s: np.mean([K_sh[e["key"]][s] for e in singles], axis=0) for s in sers}
            for e in all_entries:
                rows, pts, cell_out = rep_stats(Kd_by_entry[e["key"]], ccs[e["key"]], xcc, perm, fracs, keys, cells)
                for row in rows:
                    row.update({"noise_level": lv, "rep": rep, "entry": e["key"]})
                    rz.append(row)
                for sk, pt in pts.items():
                    ptz.setdefault((lv, e["key"], sk[0], sk[1]), []).append(pt)
                    clz.setdefault((lv, e["key"], sk[0], sk[1]), []).append(cell_out[sk])
            if lv > 0 and ((rep + 1) % 10 == 0 or rep + 1 == reps):
                print("    noise %-5s replicate %d/%d   (%.0f s so far)" % (level_label(lv), rep + 1, reps, time.time() - t0), flush=True)
        print("  noise %-5s done (%d replicate%s)" % (level_label(lv), reps, "" if reps == 1 else "s"), flush=True)
    return pd.DataFrame(rz), ptz, clz


# ------------------------------------------------------------------
# Summaries
# ------------------------------------------------------------------
def summarize(rz, names):
    rows = []
    for (lv, frac, ek, s), g in rz.groupby(["noise_level", "top_frac", "entry", "series"], sort=False):
        sl = g["slope"].to_numpy(float)
        nt = g["n_cells_tracking"].to_numpy(float)
        rows.append({
            "noise_level": lv, "top_frac": frac, "entry": ek, "entry_name": names[ek], "series": s, "n_reps": len(g),
            "slope_median": float(np.median(sl)), "slope_p10": float(np.percentile(sl, 10)), "slope_p90": float(np.percentile(sl, 90)),
            "share_significant": float(np.mean(g["p_slope"].to_numpy(float) < ALPHA)),
            "share_slope_positive": float(np.mean(sl > 0)),
            "typical_error_median": float(np.median(g["typical_error_factor"])),
            "null_typical_error_median": float(np.median(g["null_typical_error_factor"])),
            "share_error_smaller": float(np.mean(g["p_error_smaller"].to_numpy(float) < ALPHA)),
            "n_cells": int(g["n_cells"].iloc[0]), "cells_tracking_median": float(np.median(nt)),
            "cells_tracking_p10": float(np.percentile(nt, 10)), "cells_tracking_p90": float(np.percentile(nt, 90))})
    return pd.DataFrame(rows)


def on_minus_control(rz, names):
    on = rz[rz.series == SER_ON][["noise_level", "rep", "entry", "top_frac", "slope"]]
    off = rz[rz.series == SER_OFF][["noise_level", "rep", "entry", "top_frac", "slope"]]
    m = on.merge(off, on=["noise_level", "rep", "entry", "top_frac"], suffixes=("_on", "_off"))
    rows = []
    for (lv, frac, ek), g in m.groupby(["noise_level", "top_frac", "entry"], sort=False):
        d = (g["slope_on"] - g["slope_off"]).to_numpy(float)
        rows.append({"noise_level": lv, "top_frac": frac, "entry": ek, "entry_name": names[ek], "n_reps": len(d),
                     "diff_median": float(np.median(d)), "diff_p10": float(np.percentile(d, 10)),
                     "diff_p90": float(np.percentile(d, 90)), "share_diff_positive": float(np.mean(d > 0))})
    return pd.DataFrame(rows)


def per_truth_table(ptz, truths, names):
    rows = []
    for (lv, ek, s, frac), lst in ptz.items():
        A = np.vstack(lst)
        for ti, tr in enumerate(truths):
            v = A[:, ti]
            rows.append({"noise_level": lv, "top_frac": frac, "entry": ek, "entry_name": names[ek], "series": s, "truth_id": tr["id"],
                         "truth_Ks": tr["Ks"], "truth_f": tr["f"], "truth_cc_mmhr": tr["cc"], "n_reps": len(v),
                         "median_cc": float(10 ** np.median(v)), "cc_p10": float(10 ** np.percentile(v, 10)),
                         "cc_p90": float(10 ** np.percentile(v, 90)), "factor_vs_truth": float(10 ** (np.median(v) - np.log10(tr["cc"])))})
    return pd.DataFrame(rows)


def cells_table(clz, keys, names):
    rows = []
    for (lv, ek, s, frac), lst in clz.items():
        S = np.vstack([a for a, _ in lst])
        Pv = np.vstack([b for _, b in lst])
        for ci, k in enumerate(keys):
            sl, pv = S[:, ci], Pv[:, ci]
            rows.append({"noise_level": lv, "top_frac": frac, "entry": ek, "entry_name": names[ek], "series": s, "cell": L.cell_name(k),
                         "cell_Ks": k[0], "cell_f": k[1], "n_reps": len(sl), "slope_median": float(np.median(sl)),
                         "slope_p10": float(np.percentile(sl, 10)), "slope_p90": float(np.percentile(sl, 90)),
                         "share_tracking": float(np.mean((sl > 0) & (pv < ALPHA)))})
    return pd.DataFrame(rows)


def compare_with_saved(summ, summary_dir, fracs):
    """The no-noise rows should reproduce the pooled cc slopes saved by analyze_storm_compare_110.py. Returns (worst difference
    or None, message)."""
    p = summary_dir / EARLIER_DIR / "storm_compare_pooled_110.csv"
    if not p.exists():
        return None, ("no saved %s/%s found to compare with (run analyze_storm_compare_110.py first if you want this check)."
                      % (EARLIER_DIR, p.name))
    old = pd.read_csv(p)
    old = old[old["param"] == "cc"]
    new = summ[summ.noise_level == 0.0][["entry", "series", "top_frac", "slope_median"]]
    m = new.merge(old[["entry", "series", "top_frac", "slope"]], on=["entry", "series", "top_frac"], how="inner")
    if m.empty:
        return None, "the saved storm_compare_pooled_110.csv has no matching rows (different storms or cuts), so no comparison was made."
    d = float(np.max(np.abs(m["slope_median"].to_numpy(float) - m["slope"].to_numpy(float))))
    return d, ("the no-noise slopes (cc; ON and control; storms and combined; %d values) compared with the saved storm_compare_pooled_110.csv: "
               "largest difference %.1e%s" % (len(m), d, " (identical to rounding)" if d <= 1e-9 else ""))


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------
def _lookup(df, **kw):
    q = df
    for k, v in kw.items():
        q = q[q[k] == v]
    return None if q.empty else q.iloc[0]


def print_summary(summ, omc, entries, fracs, levels, n_by_frac, n_runs, n_truths):
    for frac in fracs:
        print("\n" + "=" * 156)
        print("CC RECOVERY WITH NOISY OBSERVATIONS, THREE STORMS -- best %s of runs by KGE_2012 (%d of %d runs), pooled over %d truths"
              % (cut_label(frac), n_by_frac[frac], n_runs, n_truths))
        print("=" * 156)
        print("%-18s %-6s | %-24s %-6s| %-15s | %-24s %-6s| %-24s | %-16s"
              % ("entry", "noise", "ON: tracking slope", "signif", "typ.error", "control: slope", "signif", "ON minus control", "cells tracking"))
        print("%-18s %-6s | %-24s %-6s| %-15s | %-24s %-6s| %-24s | %-16s"
              % ("", "", "[middle 80% of repl.]", "", "(no-info)", "[middle 80% of repl.]", "", "[middle 80% of repl.]", "of %d: ON / ctl" % summ["n_cells"].iloc[0]))
        for e in entries:
            for lv in levels:
                a = _lookup(summ, top_frac=frac, entry=e["key"], series=SER_ON, noise_level=lv)
                c = _lookup(summ, top_frac=frac, entry=e["key"], series=SER_OFF, noise_level=lv)
                d = _lookup(omc, top_frac=frac, entry=e["key"], noise_level=lv)
                if a is None:
                    continue

                def cell(r):
                    if r is None:
                        return "n/a", "n/a"
                    if r.n_reps == 1:
                        return "%+5.2f" % r.slope_median, "yes" if r.share_significant > 0 else "no"
                    return ("%+5.2f [%+5.2f, %+5.2f]" % (r.slope_median, r.slope_p10, r.slope_p90), "%3.0f%%" % (100 * r.share_significant))
                ca, sa = cell(a)
                cc_, sc_ = cell(c)
                terr = "x%.2f (x%.2f)" % (a.typical_error_median, a.null_typical_error_median)
                if d is None:
                    dtxt = "n/a"
                elif d.n_reps == 1:
                    dtxt = "%+5.2f" % d.diff_median
                else:
                    dtxt = "%+5.2f [%+5.2f, %+5.2f]" % (d.diff_median, d.diff_p10, d.diff_p90)
                if a.n_reps == 1:
                    ctxt = "%d / %s" % (round(a.cells_tracking_median), "n/a" if c is None else "%d" % round(c.cells_tracking_median))
                else:
                    ctxt = "%d [%d-%d] / %s" % (round(a.cells_tracking_median), round(a.cells_tracking_p10), round(a.cells_tracking_p90),
                                                "n/a" if c is None else "%d" % round(c.cells_tracking_median))
                print("%-18s %-6s | %-24s %-6s| %-15s | %-24s %-6s| %-24s | %-16s"
                      % (e["name"] if lv == levels[0] else "", level_label(lv), ca, sa, terr, cc_, sc_, dtxt, ctxt))
            print()
    print("  'none' = the exact truth, as in analyze_storm_compare_110.py. Other rows: the median over the noisy replicates, with the middle 80% "
          "of the replicates in brackets.")
    print("  signif: with no noise a yes or no; with noise, the share of replicates in which shuffling cc rarely (p < 0.05) gives a slope this large.")
    print("  typ.error: the typical factor by which the best runs' median misses the true cc (x1.00 = perfect); in brackets what a no-information "
          "fit gives.")
    print("  cells tracking: of the 9 (Ks, f) cells, how many show cc tracking (slope > 0 and p < 0.05), median over replicates [middle 80%].")
    if any(e["key"] == COMB_KEY for e in entries):
        print("  'combined' = each run's KGE averaged over the storms, each storm with its own independent volume error b"
              + ("; 'combined (shared)' = the same b in all storms (a rating-curve error at one gauge)." if any(e["key"] == SHARED_KEY for e in entries) else "."))
    print("  The control should keep a slope near 0 at every noise level and a low 'signif' share; if not, noise alone creates a pattern. "
          "Replicates reuse the same 250 runs, so the control share is not a 5% calibration.")


def summary_lines(summ, omc, entries, levels, frac):
    out = []
    keys = [e["key"] for e in entries]
    nm = {e["key"]: e["name"] for e in entries}
    singles = [e["key"] for e in entries if e["key"] not in (COMB_KEY, SHARED_KEY)]

    def sl(ek, s, lv):
        r = _lookup(summ, top_frac=frac, entry=ek, series=s, noise_level=lv)
        return None if r is None else r
    parts = []
    for ek in keys:
        r = sl(ek, SER_ON, 0.0)
        if r is not None and ek != SHARED_KEY:
            c = sl(ek, SER_OFF, 0.0)
            parts.append("%s %+.2f%s" % (nm[ek], r.slope_median, "" if c is None else " (control %+.2f)" % c.slope_median))
    out.append("No noise, best %s of runs, pooled cc slope: %s.  (0 = ignores the truth, 1 = follows it one-for-one.)" % (cut_label(frac), "; ".join(parts)))
    for lv in [x for x in levels if x > 0]:
        parts = []
        for ek in keys:
            r = sl(ek, SER_ON, lv)
            if r is not None:
                parts.append("%s %+.2f [%+.2f, %+.2f], significant in %.0f%%" % (nm[ek], r.slope_median, r.slope_p10, r.slope_p90, 100 * r.share_significant))
        out.append("Noise %s, ON slope (middle 80%% of replicates): %s." % (level_label(lv), "; ".join(parts)))
    lim = []
    for ek in keys:
        first = None
        for lv in [x for x in levels if x > 0]:
            r = sl(ek, SER_ON, lv)
            if r is not None and r.slope_p10 <= 0:
                first = lv
                break
        lim.append("%s: %s" % (nm[ek], "reaches 0 at %s" % level_label(first) if first is not None else "stays above 0 at every level tried"))
    out.append("Lowest noise level at which the middle 80%% of the ON slope reaches 0 or below -- %s." % "; ".join(lim))
    if COMB_KEY in keys:
        for lv in [x for x in levels if x > 0]:
            best = max(singles, key=lambda k: -9 if sl(k, SER_ON, lv) is None else sl(k, SER_ON, lv).slope_median)
            rc, rb = sl(COMB_KEY, SER_ON, lv), sl(best, SER_ON, lv)
            txt = "Noise %s: combined %+.2f against the best single storm (%s) %+.2f" % (level_label(lv), rc.slope_median, nm[best], rb.slope_median)
            if SHARED_KEY in keys:
                rs = sl(SHARED_KEY, SER_ON, lv)
                txt += "; combined with a shared volume error %+.2f" % rs.slope_median
            out.append(txt + ".")
    cmax, cshare, cexcl = 0.0, 0.0, []
    for lv in levels:
        for ek in keys:
            r = sl(ek, SER_OFF, lv)
            if r is not None:
                cmax = max(cmax, abs(r.slope_median))
                if lv > 0:
                    cshare = max(cshare, r.share_significant)
    out.append("Control (cc inert): the largest |median slope| over all entries and noise levels is %.2f, and the largest share of replicates called "
               "significant is %.0f%%. A control far from 0 would mean noise alone creates a pattern." % (cmax, 100 * cshare))
    low_by_level = []
    for lv in [x for x in levels if x > 0]:
        low = []
        for ek in keys:
            r = _lookup(omc, top_frac=frac, entry=ek, noise_level=lv)
            if r is not None and r.diff_p10 <= 0:
                low.append(nm[ek])
        if low:
            low_by_level.append("noise %s: %s" % (level_label(lv), ", ".join(low)))
    if low_by_level:
        out.append("ON minus control: the middle 80%% of replicates reaches 0 or below at -- %s. Everywhere else it stays above 0." % "; ".join(low_by_level))
    else:
        out.append("ON minus control: for every entry and every noise level tried, the middle 80% of replicates stays above 0.")
    return out


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def fig_slope(summ, entries, fracs, levels, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    nr, nc = len(fracs), len(entries)
    fig, axes = plt.subplots(nr, nc, figsize=(2.9 * nc + 0.7, 2.7 * nr + 1.5), facecolor=SURFACE, sharey=True, squeeze=False)
    xs = np.arange(len(levels))
    ymin = min(-0.35, float(np.nanmin(summ.slope_p10)) - 0.1)
    ymax = max(1.1, float(np.nanmax(summ.slope_p90)) + 0.1)
    for ri, frac in enumerate(fracs):
        for ci, e in enumerate(entries):
            ax = axes[ri][ci]
            T._style(ax)
            ax.axhline(1.0, color=INK, linestyle="--", linewidth=1.0, zorder=2)
            ax.axhline(0.0, color=COL_NULL, linestyle=":", linewidth=1.4, zorder=2)
            for s, col, off, mk, filled in ((SER_OFF, COL_NULL, +0.08, "s", False), (SER_ON, COL_ON, -0.08, "o", True)):
                g = summ[(summ.top_frac == frac) & (summ.entry == e["key"]) & (summ.series == s)].set_index("noise_level")
                lv_ok = [i for i, lv in enumerate(levels) if lv in g.index]
                if not lv_ok:
                    continue
                x = xs[lv_ok] + off
                med = np.array([g.loc[levels[i], "slope_median"] for i in lv_ok])
                lo = np.array([g.loc[levels[i], "slope_p10"] for i in lv_ok])
                hi = np.array([g.loc[levels[i], "slope_p90"] for i in lv_ok])
                ax.vlines(x, lo, hi, color=col, linewidth=2.2 if s == SER_ON else 1.5, zorder=3)
                ax.plot(x, med, linestyle="-", color=col, linewidth=1.0, alpha=0.6, zorder=3)
                ax.plot(x, med, linestyle="none", marker=mk, markersize=5.5, markerfacecolor=col if filled else SURFACE,
                        markeredgecolor=col, markeredgewidth=1.4, zorder=4)
            ax.set_xticks(xs)
            ax.set_xticklabels([level_label(lv) for lv in levels], fontsize=8)
            ax.set_xlim(-0.5, len(levels) - 0.5)
            ax.set_ylim(ymin, ymax)
            if ri == 0:
                ax.set_title(e["name"], color=INK, fontsize=9.5, loc="left")
            if ri == nr - 1:
                ax.set_xlabel("observation error (sd)", color=INK2, fontsize=8.5)
            if ci == 0:
                ax.set_ylabel("best %s of runs\ntracking slope" % cut_label(frac), color=INK2, fontsize=8.5)
    handles = [
        Line2D([0], [0], color=COL_ON, marker="o", markersize=5.5, linewidth=2.2, label="channel loss ON: median over replicates, bar = middle 80% of replicates"),
        Line2D([0], [0], color=COL_NULL, marker="s", markersize=5.5, markerfacecolor=SURFACE, markeredgewidth=1.4, linewidth=1.5,
               label="control (cc inert): the same"),
        Line2D([0], [0], color=INK, linestyle="--", linewidth=1.0, label="1 = best runs follow the true cc one-for-one"),
        Line2D([0], [0], color=COL_NULL, linestyle=":", linewidth=1.4, label="0 = best runs ignore cc"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=INK2)
    fig.suptitle("Do the best runs still follow the true cc when the observed hydrograph has error?  (27 truths per storm)",
                 color=INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.1 if nr <= 2 else 0.075, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True


def fig_share(summ, entries, fracs, levels, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    nc = len(fracs)
    fig, axes = plt.subplots(1, nc, figsize=(4.4 * nc + 0.6, 4.6), facecolor=SURFACE, sharey=True, squeeze=False)
    xs = np.arange(len(levels))
    singles = [e for e in entries if e["key"] not in (COMB_KEY, SHARED_KEY)]
    mk = {}
    for i, e in enumerate(singles):
        mk[e["key"]] = ENTRY_MARKERS.get("single%d" % i, "o")
    mk[COMB_KEY], mk[SHARED_KEY] = ENTRY_MARKERS[COMB_KEY], ENTRY_MARKERS[SHARED_KEY]
    for ci, frac in enumerate(fracs):
        ax = axes[0][ci]
        T._style(ax)
        ax.axhline(ALPHA, color=COL_NULL, linestyle=":", linewidth=1.4, zorder=2)
        n = len(entries)
        for k, e in enumerate(entries):
            g = summ[(summ.top_frac == frac) & (summ.entry == e["key"]) & (summ.series == SER_ON)].set_index("noise_level")
            lv_ok = [i for i, lv in enumerate(levels) if lv in g.index]
            if not lv_ok:
                continue
            x = xs[lv_ok] + 0.10 * (k - (n - 1) / 2.0)
            y = np.array([g.loc[levels[i], "share_significant"] for i in lv_ok])
            dashed = e["key"] in (COMB_KEY, SHARED_KEY)
            ax.plot(x, y, linestyle="--" if dashed else "-", color=COL_ON, linewidth=1.1, alpha=0.7, zorder=3)
            ax.plot(x, y, linestyle="none", marker=mk[e["key"]], markersize=6, markerfacecolor=COL_ON if e["key"] != SHARED_KEY else SURFACE,
                    markeredgecolor=COL_ON, markeredgewidth=1.4, zorder=4)
        cmax = []
        for lv in levels:
            v = summ[(summ.top_frac == frac) & (summ.series == SER_OFF) & (summ.noise_level == lv)]["share_significant"]
            cmax.append(float(v.max()) if len(v) else np.nan)
        ax.plot(xs, cmax, linestyle="-", color=COL_NULL, linewidth=1.6, marker="s", markersize=5, markerfacecolor=SURFACE,
                markeredgecolor=COL_NULL, markeredgewidth=1.3, zorder=3)
        ax.set_xticks(xs)
        ax.set_xticklabels([level_label(lv) for lv in levels], fontsize=8.5)
        ax.set_xlim(-0.5, len(levels) - 0.5)
        ax.set_ylim(-0.04, 1.05)
        ax.set_title("best %s of runs" % cut_label(frac), color=INK, fontsize=9.5, loc="left")
        ax.set_xlabel("observation error (sd)", color=INK2, fontsize=8.5)
        if ci == 0:
            ax.set_ylabel("share of replicates with a significant cc slope", color=INK2, fontsize=8.5)
    handles = []
    for e in entries:
        dashed = e["key"] in (COMB_KEY, SHARED_KEY)
        handles.append(Line2D([0], [0], color=COL_ON, marker=mk[e["key"]], markersize=6, linestyle="--" if dashed else "-", linewidth=1.1,
                              markerfacecolor=COL_ON if e["key"] != SHARED_KEY else SURFACE, markeredgewidth=1.4, label="ON: " + e["name"]))
    handles.append(Line2D([0], [0], color=COL_NULL, marker="s", markersize=5, markerfacecolor=SURFACE, markeredgewidth=1.3, linewidth=1.6,
                          label="control (cc inert): the highest of the entries"))
    handles.append(Line2D([0], [0], color=COL_NULL, linestyle=":", linewidth=1.4, label="5% = what chance should give"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8, labelcolor=INK2)
    fig.suptitle("How often is the cc tracking still detected? (none = one value, so 0 or 1)", color=INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.15, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True


def fig_recovery(ptab, entries, show_keys, frac, levels, truths, tx, keys, cells, path):
    import matplotlib
    import matplotlib.ticker
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    names = {e["key"]: e["name"] for e in entries}
    show = [k for k in show_keys if k in names]
    if not show:
        return False
    nr, nc = len(show), len(levels)
    fig, axes = plt.subplots(nr, nc, figsize=(3.1 * nc + 0.7, 3.0 * nr + 1.4), facecolor=SURFACE, sharex=True, sharey=True, squeeze=False)
    lo, hi = 10 ** R.LOG_LO, 10 ** R.LOG_HI
    lim = (lo * 0.75, hi * 1.3)
    cc_all = sorted({t["cc"] for t in truths})
    rank = {}
    for ri_, k in enumerate(keys):
        for ti in cells[k]:
            rank[ti] = ri_
    nk = max(1, len(keys) - 1)
    d = ptab[ptab.top_frac == frac]
    for ri, ek in enumerate(show):
        for ci, lv in enumerate(levels):
            ax = axes[ri][ci]
            T._style(ax)
            ax.axhspan(10 ** (R.LOG_LO + 0.1 * (R.LOG_HI - R.LOG_LO)), 10 ** (R.LOG_LO + 0.9 * (R.LOG_HI - R.LOG_LO)),
                       color=COL_NULL, alpha=0.16, linewidth=0, zorder=1)
            ax.axhline(10 ** R.CENTER_LOG, color=COL_NULL, linestyle=":", linewidth=1.3, zorder=2)
            ax.plot(lim, lim, color=INK, linestyle="--", linewidth=1.0, zorder=2)
            for s, col, sh, mk, filled in ((SER_OFF, COL_NULL, +0.012, "s", False), (SER_ON, COL_ON, -0.012, "o", True)):
                g = d[(d.entry == ek) & (d.series == s) & (d.noise_level == lv)].set_index("truth_id")
                if g.empty:
                    continue
                ti_list = [ti for ti, tr in enumerate(truths) if tr["id"] in g.index]
                xx = np.array([10 ** (tx[ti, 2] + sh + 0.03 * (rank[ti] / nk - 0.5)) for ti in ti_list])
                med = np.array([g.loc[truths[ti]["id"], "median_cc"] for ti in ti_list])
                lo_ = np.array([g.loc[truths[ti]["id"], "cc_p10"] for ti in ti_list])
                hi_ = np.array([g.loc[truths[ti]["id"], "cc_p90"] for ti in ti_list])
                if lv > 0:
                    ax.vlines(xx, lo_, hi_, color=col, linewidth=1.5 if s == SER_ON else 1.1, alpha=0.8, zorder=3)
                ax.plot(xx, med, linestyle="none", marker=mk, markersize=4.6, markerfacecolor=col if filled else SURFACE,
                        markeredgecolor=col, markeredgewidth=1.2, alpha=0.95, zorder=4 if s == SER_ON else 3)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(lim)
            ax.set_ylim(lim)
            ax.set_xticks(cc_all)
            ax.set_xticklabels(["%g" % v for v in cc_all], fontsize=7.5)
            ax.set_yticks([30, 100, 300, 1000])
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
            ax.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
            if ri == 0:
                ax.set_title("noise: none" if lv == 0 else "noise sd %s" % level_label(lv), color=INK, fontsize=9.5, loc="left")
            if ri == nr - 1:
                ax.set_xlabel("true cc (mm/hr)", color=INK2, fontsize=8.5)
            if ci == 0:
                ax.set_ylabel("%s\nmedian cc of best runs (mm/hr)" % names[ek], color=INK2, fontsize=8.5)
    handles = [
        Line2D([0], [0], color=COL_ON, marker="o", markersize=4.6, linewidth=1.5, label="channel loss ON: median over replicates (one dot per truth), bar = middle 80%"),
        Line2D([0], [0], color=COL_NULL, marker="s", markersize=4.6, markerfacecolor=SURFACE, markeredgewidth=1.2, linewidth=1.1, label="control (cc inert): the same"),
        Line2D([0], [0], color=INK, linestyle="--", linewidth=1.0, label="perfect recovery"),
        Patch(facecolor=COL_NULL, alpha=0.16, label="no information: middle 80% of the sampled range (dotted = its middle)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=INK2)
    fig.suptitle("Best %s of runs: where does the median cc land when the observed hydrograph has error?  (27 truths per panel, 9 per true-cc column)"
                 % cut_label(frac), color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.1, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110 -- three storms plus observation error: does the tracking of the true cc survive?")
    ap.add_argument("--storms", type=float, nargs="+", default=[0.8, 1.25],
                    help="Rain multipliers to use besides the original storm (default 0.8 1.25). Storms whose re-scored files do not exist "
                         "are skipped with a note.")
    ap.add_argument("--top_fracs", type=float, nargs="+", default=[0.20, 0.10, 0.05],
                    help="Top fractions by KGE_2012; the first is the main one (default 0.20 0.10 0.05)")
    ap.add_argument("--noise_levels", type=float, nargs="+", default=[0.05, 0.10, 0.20], help="Noise sd levels as fractions (default 0.05 0.10 0.20)")
    ap.add_argument("--noise_reps", type=int, default=100, help="Noisy replicates per level (default 100)")
    ap.add_argument("--noise_parts", choices=["both", "volume", "shape"], default="both",
                    help="both = whole-hydrograph error b and per-step scatter e (default); volume = only b; shape = only e")
    ap.add_argument("--perm", type=int, default=1000, help="Permutation shuffles per replicate (default 1000)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--expect_n", type=int, default=250)
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--check_only", action="store_true", help="Load everything, run the input checks and the noise-free check, then stop. Writes nothing.")
    ap.add_argument("--continue_on_fail", action="store_true")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()
    if any(lv <= 0 for lv in args.noise_levels):
        sys.exit("--noise_levels must be positive fractions, e.g. 0.05 0.10 0.20 (the no-noise case is always included).")
    if len(set(args.noise_levels)) != len(args.noise_levels):
        sys.exit("--noise_levels lists the same level twice.")
    if args.noise_reps < 2:
        sys.exit("--noise_reps must be at least 2.")
    if args.perm < 200:
        sys.exit("--perm must be at least 200.")
    if args.seed < 0:
        sys.exit("--seed must be 0 or more.")
    if not args.top_fracs or any(not (0 < f <= 1) for f in args.top_fracs):
        sys.exit("--top_fracs must be fractions between 0 and 1, e.g. 0.20 0.10 0.05.")
    fracs = list(args.top_fracs)

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    calib = summary_dir.parent.parent
    csv_dir = calib / "03_comparisons" / "csv_exports"
    out_dir = summary_dir / OUT_DIRNAME
    tag = "" if args.noise_parts == "both" else "_" + args.noise_parts
    print("\n" + "=" * 78 + "\nSeries 110 -- three storms plus observation error: does the tracking of the true cc survive?\n" + "=" * 78)
    print("python %s, pandas %s, numpy %s" % (sys.version.split()[0], pd.__version__, np.__version__))
    print("Reading from: %s" % summary_dir)

    # ---- load every storm that exists (same loading and checks as analyze_storm_compare_110.py) ----------------------
    scales = [REF_SCALE] + sorted(set(float(s) for s in args.storms if abs(float(s) - REF_SCALE) > 1e-9))
    ents, notes = [], []
    for sc in scales:
        pth = SC.entry_paths(summary_dir, sc)
        name = SC.scale_name(sc)
        if not pth["long"].exists():
            if abs(sc - REF_SCALE) < 1e-9:
                sys.exit("Required file not found: %s\nThe original storm's re-scored file is needed as the reference." % pth["long"])
            notes.append("%s (%s): %s not found, so this storm is skipped. Run: python rescore_storm_110.py --rain_scale %g"
                         % (name, pth["label"], pth["long"].name, sc))
            continue
        ents.append({"key": pth["label"], "name": name, "scale": sc, "paths": pth, "df": SC.load_long(pth["long"], pth["label"], sc)})
    for n in notes:
        print("NOTE: " + n)
    D = P.generate_lhs_samples(args.expect_n, P.LHS_PARAMS, P.DESIGN_SEED)[P.PARAMS].to_numpy(float)
    ck = P.Checks()
    for e in ents:
        SC.check_entry(ck, e, D, args.expect_n)
    if len(ents) >= 2:
        truth_ids = SC.check_cross(ck, ents)
    else:
        td0 = ents[0]["df"].drop_duplicates("truth_id").set_index("truth_id")[["truth_Ks", "truth_f", "truth_cc"]]
        truth_ids = list(td0.reset_index().sort_values(["truth_Ks", "truth_f", "truth_cc"])["truth_id"])
    P.print_checks(ck)
    if ck.n("FAIL") and not args.continue_on_fail:
        sys.exit("\nSTOPPED: an input check failed. Nothing was analysed or written. (--continue_on_fail to override.)")
    if len(truth_ids) < 3:
        sys.exit("\nFewer than 3 truths are common to all storms; cannot compute tracking slopes.")
    for e in ents:
        e["data"] = SC.build_entry_data(e["df"], truth_ids, D)
    singles = sorted(ents, key=lambda e: e["scale"])
    ref = [e for e in singles if e["key"] == REF_KEY][0]
    td = ref["df"].drop_duplicates("truth_id").set_index("truth_id")
    truths = [{"id": t, "Ks": float(td.loc[t, "truth_Ks"]), "f": float(td.loc[t, "truth_f"]), "cc": float(td.loc[t, "truth_cc"])} for t in truth_ids]
    tx = L.truth_coords(truths)
    xcc = tx[:, 2]
    keys, cells, centre = L.find_cells(truths)
    n_runs = ref["data"][SER_ON]["X"].shape[0]

    extra = []
    if len(singles) >= 2:
        comb = SC.make_combined(singles, REF_KEY)
        extra.append(comb)
        if args.noise_parts != "shape":
            extra.append({"key": SHARED_KEY, "name": "combined (shared)", "scale": None, "data": comb["data"]})
        if SER_OFF not in comb["data"]:
            print("NOTE: not every storm has a control, so there is no control for the combined entries.")
    entries = list(singles) + extra
    names = {e["key"]: e["name"] for e in entries}
    levels = [0.0] + list(args.noise_levels)
    n_by_frac = {f: T.n_top_for(n_runs, f) for f in fracs}
    print("\nStorms: %s.  Entries analysed: %s." % (", ".join(e["name"] for e in singles), ", ".join(e["name"] for e in entries)))
    print("%d truth(s) in %d cell(s), %d runs per series, cuts %s, noise sd %s, %d replicates per level, %d shuffles per replicate, noise = %s"
          % (len(truths), len(keys), n_runs, [cut_label(f) for f in fracs], ["%g%%" % round(100 * lv, 6) for lv in args.noise_levels],
             args.noise_reps, args.perm, args.noise_parts))

    # ---- stored hydrographs, verified truths, noise-free check -------------------------------------------------------
    print("\nReading the stored hydrographs and truth files ...", flush=True)
    sims, truth_obs, md5s, nf, n_t_by_storm = {}, {}, {}, {}, {}
    for e in singles:
        sims[e["key"]], idx = load_storm_hydrographs(e, csv_dir)
        n_t_by_storm[e["key"]] = len(idx)
        truth_obs[e["key"]], md5s[e["key"]] = load_truth_hydrographs(e, calib, idx, truth_ids)
        nf[e["key"]] = noise_free_check(e, sims[e["key"]], truth_obs[e["key"]], truth_ids)
        for s, r in nf[e["key"]].items():
            bad = r["worst_abs_diff"] > NOISE_FREE_TOL or r["n_finite_mismatch"] > 0
            print("  %-6s %-6s: %d runs x %d five-minute steps; noise-free KGE_2012 vs stored: worst difference %.1e over %d (truth, run) pairs%s"
                  % (e["name"], s, sims[e["key"]][s].shape[0], n_t_by_storm[e["key"]], r["worst_abs_diff"], r["n_pairs"],
                     "" if not bad else "   <-- MISMATCH (%d finite/non-finite mismatches)" % r["n_finite_mismatch"]))
            if bad:
                sys.exit("STOPPED: re-computing KGE_2012 without noise does not reproduce the stored values for %s (%s): worst difference %.3g, "
                         "%d finite/non-finite mismatches. The noisy results would not be comparable to the main analysis."
                         % (e["name"], s, r["worst_abs_diff"], r["n_finite_mismatch"]))
    print("  Truth files: every one matches the checksum recorded at re-scoring.")
    nf_rows = [{"storm": e["name"], "series": s, **r} for e in singles for s, r in nf[e["key"]].items()]
    n_draws = (len(levels) - 1) * args.noise_reps * len(singles) * len(truths)
    if args.check_only:
        print("\n--check_only: loaded and checked, nothing written. A full run would draw %d noisy hydrographs (%d storm(s) x %d truths x %d "
              "replicates x %d noise level(s)), each scored against the cc-ON and the control runs."
              % (n_draws, len(singles), len(truths), args.noise_reps, len(args.noise_levels)))
        return

    # ---- noise ---------------------------------------------------------------------------------------------------
    perm = {s: T.make_perm_matrix(n_runs, args.perm, np.random.default_rng([args.seed, 99, si])) for si, s in enumerate((SER_ON, SER_OFF))}
    print("\nAdding noise and re-picking the best runs (several minutes at the default settings; progress is printed) ...", flush=True)
    t0 = time.time()
    rz, ptz, clz = run_replicates(args, singles, extra, truths, xcc, keys, cells, sims, truth_obs, perm, fracs, n_t_by_storm)
    print("  finished in %.0f s" % (time.time() - t0))
    summ = summarize(rz, names)
    omc = on_minus_control(rz, names)
    ptab = per_truth_table(ptz, truths, names)
    ctab = cells_table(clz, keys, names)
    worst, msg = compare_with_saved(summ, summary_dir, fracs)
    print("\nCross-check: " + msg)
    if worst is not None and worst > 1e-9:
        print("  WARN: not identical. The no-noise numbers here should equal the saved ones; a difference means the runs were lined up "
              "differently, different cuts were used, or best runs were tied. Treat the results with caution until this is explained.")

    print_summary(summ, omc, entries, fracs, levels, n_by_frac, n_runs, len(truths))
    print("\n" + "=" * 156)
    print("WHAT THE NUMBERS SAY   (best %s of runs; descriptive, assembled from the tables above)" % cut_label(fracs[0]))
    print("=" * 156)
    for line in summary_lines(summ, omc, entries, levels, fracs[0]):
        print("  - " + line)
    print("  Caveat: the storms share the same 250 sweep points, routing and rainfall pattern (only scaled), the noise has no memory in time, and "
          "real gauge error can be larger at high flows. This shows how much observation error these runs can absorb, not what a real gauge's "
          "error is.")

    # ---- save --------------------------------------------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ck.rows).to_csv(out_dir / "storm_noise_input_checks_110.csv", index=False)
    pd.DataFrame(nf_rows).to_csv(out_dir / "storm_noise_noisefree_check_110.csv", index=False)
    rz_out = rz.copy()
    rz_out.insert(2, "entry_name", rz_out["entry"].map(names))
    rz_out.to_csv(out_dir / ("storm_noise_realizations_110%s.csv" % tag), index=False)
    summ.to_csv(out_dir / ("storm_noise_summary_110%s.csv" % tag), index=False)
    omc.to_csv(out_dir / ("storm_noise_on_minus_control_110%s.csv" % tag), index=False)
    ptab.to_csv(out_dir / ("storm_noise_per_truth_110%s.csv" % tag), index=False)
    ctab.to_csv(out_dir / ("storm_noise_cells_110%s.csv" % tag), index=False)
    if not args.no_plots:
        show_keys = [REF_KEY, COMB_KEY, SHARED_KEY]
        for name, fn in (
                ("fig_stormnoise_slope_110%s.png" % tag, lambda p: fig_slope(summ, entries, fracs, levels, p)),
                ("fig_stormnoise_share_110%s.png" % tag, lambda p: fig_share(summ, entries, fracs, levels, p)),
                ("fig_stormnoise_recovery_110%s.png" % tag,
                 lambda p: fig_recovery(ptab, entries, show_keys, fracs[0], levels, truths, tx, keys, cells, p))):
            try:
                fn(out_dir / name)
            except Exception as ex:
                print("\n  (figure %s skipped: %s -- the CSVs are saved regardless)" % (name, ex))

    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
            "settings": {k: (v if not isinstance(v, Path) else str(v)) for k, v in vars(args).items()},
            "storms": [{"key": e["key"], "name": e["name"], "rain_scale": e["scale"], "has_control": SER_OFF in e["data"],
                        "long_file": str(e["paths"]["long"]), "long_md5": SC.md5_file(e["paths"]["long"]),
                        "sweep_on_md5": SC.md5_file(e["paths"]["on_csv"]),
                        "sweep_control_md5": SC.md5_file(e["paths"]["off_csv"]) if e["paths"]["off_csv"].exists() else None,
                        "truth_md5": md5s[e["key"]], "rescoring_pandas": e.get("pandas"),
                        "noise_free_check": nf[e["key"]]} for e in singles],
            "storms_skipped": notes, "entries": [e["key"] for e in entries], "truths": truths, "cells": [list(k) for k in keys],
            "noise_free_tolerance": NOISE_FREE_TOL, "checks": ck.rows,
            "saved_slopes_crosscheck": {"worst_abs_difference": worst, "message": msg}}
    (out_dir / ("PROVENANCE_storm_noise_110%s.json" % tag)).write_text(json.dumps(prov, indent=2, default=str))

    print("\nSaved to: %s" % out_dir)
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results." % ck.n("FAIL"))


if __name__ == "__main__":
    main()
