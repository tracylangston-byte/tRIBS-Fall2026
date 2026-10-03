"""
analyze_storm_compare_110.py
============================
Series 110 -- STORM COMPARISON: does the answer about cc (and Ks, f) change when the
storm is bigger or smaller?

THE QUESTION
------------
Everything so far used one storm (the 12 Aug 2014 event). A skeptic's worry is that
cc was only recoverable because of that one storm. Two more storms were made by
multiplying the rainfall (both gauges) by 0.8 and by 1.25, and every one of the 250
sweep runs (and the 250 control runs with channel loss OFF) was re-scored against 27
cc-ON truths per storm (rescore_storm_110.py). This script puts the storms side by side:

  1. Per storm: do the best runs (by KGE_2012) put cc, Ks and f where the truth is?
     ("tracking slope" -- same definition and same cuts as analyze_truth_location_110.py)
  2. Does that change with storm size? (paired differences from the original storm)
  3. Does cc tracking hold in every (Ks, f) cell, in every storm?
  4. What if the storms are used TOGETHER? Every run is the same parameter point in
     every storm, so its scores can be averaged across storms ("combined"). That is
     what calibrating to several events at once would do. It needs no new runs.

Pure post-hoc analysis. Runs no tRIBS. Writes only into its own new folder
(calibration_work/03_comparisons/summary_tables/storm_compare_110/).

Keep this file in the same folder as analyze_truth_location_110.py,
analyze_cc_recovery_110.py, analyze_cc_tests_110.py and analyze_cc_pca_110.py
(it reuses their code). Run it from the lab/ directory.

VOCABULARY
-----------
  Storm            A rain event. "x0.8" = the original rain times 0.8, "x1.0" = the
                   original (reference), "x1.25" = times 1.25.
  Best runs        The top 20% / 10% / 5% of the 250 runs by KGE_2012 against a truth
                   (at least 10 runs).
  Tracking slope   Plot the best runs' median parameter against the true parameter (cc and
                   f on log scales, Ks linear) and fit a line. 1 = follows the truth,
                   0 = ignores it. The edges of the sampled box squash it below 1.
  Control          The same 250 parameter points with channel loss OFF. cc has no effect
                   there, so its "cc slope" shows what chance gives.
  ON minus control The ON slope minus the control slope: how much of the tracking is more
                   than chance.
  Paired bootstrap Re-draw the 250 runs with replacement and redo everything. The SAME
                   re-draw is used for every storm (the runs are the same parameter points
                   in every storm), so a difference between storms is measured on the same
                   wobble. If the 95% interval of a difference excludes 0 the storms
                   really differ in that statistic.
  Permutation test Shuffle a parameter among the runs to show what "no information" looks like.
  Combined         Each run's KGE_2012 averaged over the storms (same truth parameters in
                   every storm), then the best runs are chosen on that average.
  Median shift     How far the best runs' median moves between a storm and the reference
                   storm, averaged over the truths (a factor for cc and f, a difference for
                   Ks). The control's shift shows how much the median wanders by chance.

HOW TO READ THE RESULTS
------------------------
  * Slopes similar across storms, intervals of the differences including 0: the answer does
    not depend detectably on storm size (within what 250 runs can resolve).
  * ON slope above the control's in every storm: cc tracking is more than chance in each.
  * A storm whose slope is clearly lower: cc identifiability depends on the storm.
  * Combined higher than any single storm: using several storms together helps.
  * The control is the yardstick for "chance". A difference between storms that is no bigger
    than the control's differences is not worth reading anything into.

OUTPUTS   (calibration_work/03_comparisons/summary_tables/storm_compare_110/)
-----------------------------------------------------------------------------
  storm_compare_pooled_110.csv      entry x series x cut x parameter: tracking slope, interval, p
  storm_compare_cells_110.csv       entry x series x cut x cell: cc tracking slope
  storm_compare_per_truth_110.csv   entry x series x cut x truth: best-run medians
  storm_compare_diffs_110.csv       paired differences (from the reference storm; ON minus control)
  storm_compare_input_checks_110.csv
  fig_stormcmp_cell_slopes_110.png  cc slope per cell, per storm; ON beside control
  fig_stormcmp_pooled_slopes_110.png  pooled slope of cc, Ks and f against storm size
  fig_stormcmp_cc_recovery_110.png  best-run cc against true cc, one panel per storm
  PROVENANCE_storm_compare_110.json

USAGE (run from the lab/ directory)
------------------------------------
    python -u analyze_storm_compare_110.py 2>&1 | tee storm_compare_110.log
    python -u analyze_storm_compare_110.py --storms 0.8          # compare with one storm only
    python -u analyze_storm_compare_110.py --boot 2000 --perm 5000 --no_plots

CAVEATS
-------
The storms are NOT independent experiments: they share the same 250 sweep points, the same
routing parameters, the same rainfall pattern (only scaled) and noise-free truths, so
agreement between them is partly built in. This tests how much the answer depends on storm
SIZE, not whether it would survive a different storm shape or real noise. The 27 truths
share the same runs, so p-values describe the runs, not 27 separate experiments. Many
intervals are shown (cells x storms x cuts); a few will exclude 0 by chance alone.
"""

import argparse
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

try:
    import analyze_truth_location_110 as L
except ImportError:
    sys.exit("This script needs analyze_truth_location_110.py, analyze_cc_recovery_110.py, analyze_cc_tests_110.py and "
             "analyze_cc_pca_110.py in the same folder (it reuses their code). Put all of them in lab/ and run from there.")

R = L.R
T = L.T
P = L.P

_NEED = {
    "analyze_truth_location_110": (L, ["KEY_COLS", "RS_DIRNAME", "LONG_NAME", "RS_PROV", "truth_list", "truth_coords", "find_cells",
                                       "cell_name", "track", "ci", "PCOL", "PAR", "cut_label", "iv", "param_err_fmt", "box_pos",
                                       "_setup_log_axis"]),
    "analyze_cc_recovery_110": (R, ["top_index", "slope_of", "fcc", "ffactor", "CENTER_LOG", "LOG_LO", "LOG_HI"]),
    "analyze_cc_tests_110": (T, ["n_top_for", "make_perm_matrix", "fmt_p", "_style", "COL_NULL", "INK", "INK2", "GRID", "SURFACE"]),
    "analyze_cc_pca_110": (P, ["Checks", "print_checks", "generate_lhs_samples", "LHS_PARAMS", "DESIGN_SEED", "PARAMS", "md5_of",
                               "SER_ON", "SER_OFF", "STAGE2_NAME", "CONTROL_NAME", "NULL_RHO", "KEY_METRICS", "C_KS", "C_F", "C_CC"]),
}
for _label, (_mod, _names) in _NEED.items():
    _miss = [n for n in _names if not hasattr(_mod, n)]
    if _miss:
        sys.exit("%s.py in this folder is an older version (missing: %s). Replace it with the current copy."
                 % (_label, ", ".join(_miss)))

SER_ON, SER_OFF = P.SER_ON, P.SER_OFF
SERIES = (SER_ON, SER_OFF)
PAR = L.PAR
PCOL = L.PCOL
REF_SCALE = 1.0
REF_KEY = "storm100"
COMB_KEY = "combined"
LHS_SERIES = "110"
OUT_DIRNAME = "storm_compare_110"
EARLIER_DIR = "location_cc_110"
INK, INK2, GRID, SURFACE, NULLC = T.INK, T.INK2, T.GRID, T.SURFACE, T.COL_NULL
PHUE = {"Ks": P.C_KS, "f": P.C_F, "cc": P.C_CC}
SINGLE_MARKERS = ["v", "o", "^", "s", "P", "X"]


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def storm_label(scale):
    pct = int(round(scale * 100))
    if not np.isfinite(scale) or scale <= 0 or abs(pct / 100.0 - scale) > 1e-9:
        sys.exit("--storms needs positive rain multipliers with at most 2 decimal places (0.8, 1.25), got %r" % scale)
    return "storm%03d" % pct


def scale_name(scale):
    t = "%g" % scale
    if "." not in t:
        t += ".0"
    return "x" + t


def md5_file(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def minor(v):
    return ".".join(str(v).split(".")[:2])


def entry_paths(summary_dir, scale):
    if abs(scale - REF_SCALE) < 1e-9:
        d = summary_dir / L.RS_DIRNAME
        return {"label": REF_KEY, "long": d / L.LONG_NAME, "prov": d / L.RS_PROV,
                "on_csv": summary_dir / P.STAGE2_NAME, "off_csv": summary_dir / P.CONTROL_NAME}
    label = storm_label(scale)
    d = summary_dir / ("rescored_storm_110_%s" % label)
    return {"label": label, "long": d / ("rescored_location_long_110_%s.csv" % label),
            "prov": d / ("PROVENANCE_rescored_storm_110_%s.json" % label),
            "on_csv": summary_dir / ("lhs_results_joint_Ks_f_cc_%s_%s.csv" % (label, LHS_SERIES)),
            "off_csv": summary_dir / ("lhs_results_joint_Ks_f_cc_CONTROL_%s_%s.csv" % (label, LHS_SERIES))}


def load_long(path, label, scale):
    df = pd.read_csv(path)
    miss = [c for c in L.KEY_COLS if c not in df.columns]
    if miss:
        sys.exit("%s is missing columns: %s. Re-run the re-scoring script for this storm." % (path.name, miss))
    if "storm_label" in df.columns:
        other = sorted(set(df["storm_label"].dropna().astype(str)) - {label})
        if other:
            sys.exit("%s says it belongs to %s, not %s. Wrong file?" % (path.name, ", ".join(other), label))
    if "rain_scale" in df.columns:
        vals = pd.to_numeric(df["rain_scale"], errors="coerce").dropna()
        if len(vals) and not np.allclose(vals.to_numpy(float), scale, rtol=0, atol=1e-9):
            sys.exit("%s has rain_scale values other than %g. Wrong file?" % (path.name, scale))
    for c in L.KEY_COLS:
        if c not in ("truth_id", "series", "run_id"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def design_match(Rm, D):
    """Match each run's (Ks, f, cc) to a point of the intended design. Returns (design index per run or None,
    n unmatched, n matched more than once, n design points used)."""
    if len(Rm) == 0 or np.isnan(Rm).any():
        return None, len(Rm), 0, 0
    close = np.all(np.isclose(Rm[:, None, :], D[None, :, :], rtol=1e-3, atol=0.0), axis=2)
    cnt = close.sum(axis=1)
    unmatched, multi = int((cnt == 0).sum()), int((cnt > 1).sum())
    used = int(close.any(axis=0).sum())
    ok = unmatched == 0 and multi == 0 and used == len(D) == len(Rm)
    return (close.argmax(axis=1) if ok else None), unmatched, multi, used


# ------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------
def check_entry(ck, ent, D, expect_n):
    nm, df, pth = ent["name"], ent["df"], ent["paths"]
    truth_ids = sorted(set(df["truth_id"]))
    has_off = SER_OFF in set(df["series"])
    ent["has_control"] = has_off

    # 1. provenance of the re-scoring ------------------------------------------------
    if not pth["prov"].exists():
        ck.add("WARN", "%s provenance" % nm, "%s missing; cannot confirm the re-scoring gates passed" % pth["prov"].name)
    else:
        prov = json.loads(pth["prov"].read_text())
        ent["pandas"] = prov.get("pandas")
        names = [pth["on_csv"]] + ([pth["off_csv"]] if has_off else [])
        bad = []
        for p in names:
            rec = prov.get("inputs", {}).get(p.name, {}).get("md5")
            if rec != (md5_file(p) if p.exists() else None):
                bad.append(p.name)
        if bad:
            ck.add("FAIL", "%s inputs unchanged" % nm,
                   "checksum differs from the one recorded when this storm was re-scored (or none recorded): %s. "
                   "Re-run the re-scoring script for this storm." % ", ".join(bad))
        else:
            ck.add("PASS", "%s inputs unchanged" % nm, "%d sweep CSV checksum(s) match the re-scoring provenance" % len(names))
        gates = prov.get("gates", {})
        if gates.get("runs_excluded"):
            ck.add("WARN", "%s excluded runs" % nm, "%d run(s) were excluded at re-scoring (--drop_bad_runs)"
                   % len(gates["runs_excluded"]))
        else:
            ck.add("PASS", "%s gates A and B" % nm, "all %s runs passed gates A and B at re-scoring" % gates.get("runs_checked", "?"))
        if gates.get("D_overridden"):
            ck.add("FAIL", "%s gate D" % nm, "gate D failed at re-scoring and was overridden: the truths may not match their runs")
        dry = prov.get("dry_truths") or []
        skipped = prov.get("skipped_truths") or []
        if dry or skipped:
            ck.add("WARN", "%s truths set aside" % nm, "%d dry and %d skipped truth(s) were left out at re-scoring" % (len(dry), len(skipped)))
        verified = prov.get("truths", {})
        td = df.drop_duplicates("truth_id").set_index("truth_id")
        wrong = [t for t in truth_ids if t not in verified or not (
            np.isclose(verified[t]["Ks_mult"], td.loc[t, "truth_Ks"]) and np.isclose(verified[t]["f_RS_abs"], td.loc[t, "truth_f"])
            and np.isclose(verified[t]["cc_mmhr"], td.loc[t, "truth_cc"]))]
        if wrong:
            ck.add("FAIL", "%s truths verified" % nm, "not in the verified list or with different parameters: %s" % ", ".join(wrong[:4]))
        else:
            ck.add("PASS", "%s truths verified" % nm, "%d truth(s), all checksum-verified at re-scoring" % len(truth_ids))

    # 2. control status -------------------------------------------------------------------
    if has_off:
        ck.add("PASS", "%s control" % nm, "cc-OFF control runs are included")
    elif pth["off_csv"].exists():
        ck.add("WARN", "%s control" % nm, "a control sweep file exists but this re-scored file does not include it: re-run the "
               "re-scoring script for this storm to pick it up. Until then there is no chance baseline for this storm.")
    else:
        ck.add("WARN", "%s control" % nm, "no control (cc OFF) sweep for this storm: no chance baseline for it")

    # 3. counts and run identity ----------------------------------------------------------
    problems, sets = [], {}
    for (tid, s), g in df.groupby(["truth_id", "series"]):
        if len(g) != expect_n:
            problems.append("truth %s %s: %d runs (expected %d)" % (tid, s, len(g), expect_n))
        if g["run_id"].duplicated().any():
            problems.append("truth %s %s: duplicate run_id" % (tid, s))
        sets.setdefault(s, []).append(frozenset(g["run_id"]))
    for s, ss in sets.items():
        if len(set(ss)) > 1:
            problems.append("%s: run_id sets differ between truths" % s)
    if SER_ON not in sets:
        problems.append("no rows for series %s" % SER_ON)
    if problems:
        ck.add("FAIL", "%s run counts" % nm, "; ".join(problems[:4]))
    else:
        ck.add("PASS", "%s run counts" % nm, "%d truth(s) x %d series, %d runs each, same runs for every truth"
               % (len(truth_ids), len(sets), expect_n))

    # 4. series flags ---------------------------------------------------------------------------
    on_bad = int((df.loc[df.series == SER_ON, "optpercolation"] != 1).sum())
    off_bad = int((df.loc[df.series == SER_OFF, "optpercolation"] != 0).sum())
    if on_bad or off_bad:
        ck.add("FAIL", "%s optpercolation" % nm, "%d cc-ON rows not 1, %d cc-OFF rows not 0" % (on_bad, off_bad))
    else:
        ck.add("PASS", "%s optpercolation" % nm, "cc-ON rows all 1%s" % (", cc-OFF rows all 0" if has_off else ""))

    # 5. design: the parameters equal the intended seed-42 sweep ------------------------------------
    for s in SERIES:
        g = df[(df.truth_id == truth_ids[0]) & (df.series == s)]
        if g.empty:
            continue
        di, unmatched, multi, used = design_match(g[P.PARAMS].to_numpy(float), D)
        if di is not None:
            ck.add("PASS", "%s design (%s)" % (nm, s), "all %d (Ks, f, cc) points equal the seed-%d, n=%d design"
                   % (len(g), P.DESIGN_SEED, expect_n))
        else:
            ck.add("FAIL", "%s design (%s)" % (nm, s), "%d point(s) unmatched, %d matched twice; %d of %d design points used (seed %d)"
                   % (unmatched, multi, used, len(D), P.DESIGN_SEED))

    # 6. finite metrics ----------------------------------------------------------------------------------
    worst = (1.0, "")
    for (tid, s), gg in df.groupby(["truth_id", "series"]):
        for m in P.KEY_METRICS:
            fr = float(np.isfinite(gg[m]).mean())
            if fr < worst[0]:
                worst = (fr, "truth %s %s %s" % (tid, s, m))
    if worst[0] < 0.99:
        ck.add("FAIL", "%s finite metrics" % nm, "worst finite fraction %.3f (%s)" % worst)
    else:
        ck.add("PASS", "%s finite metrics" % nm, "all key metrics >= 99%% finite (worst %.3f)" % worst[0])

    # 7. null check: cc is inert in the control -------------------------------------------------------------
    if has_off:
        worst_rho, worst_t = 0.0, None
        for tid in truth_ids:
            g = df[(df.truth_id == tid) & (df.series == SER_OFF)]
            ok = np.isfinite(g["kge_2012"]) & np.isfinite(g["channelconductivity_mmhr"])
            if ok.sum() > 10:
                rho = stats.spearmanr(np.log10(g.loc[ok, "channelconductivity_mmhr"]), g.loc[ok, "kge_2012"])[0]
                if np.isfinite(rho) and abs(rho) > abs(worst_rho):
                    worst_rho, worst_t = rho, tid
        if worst_t is not None and abs(worst_rho) > P.NULL_RHO:
            ck.add("WARN", "%s null check" % nm, "rho(log cc, KGE_2012) = %+.2f at truth %s in the cc-OFF series (limit %.2f)"
                   % (worst_rho, worst_t, P.NULL_RHO))
        else:
            ck.add("PASS", "%s null check" % nm, "largest |rho(log cc, KGE_2012)| in the control is %.2f (limit %.2f)"
                   % (abs(worst_rho), P.NULL_RHO))


def check_cross(ck, ents):
    """Checks between storms. Returns the truth ids present in every storm, in (Ks, f, cc) order."""
    per = {}
    for e in ents:
        td = e["df"].drop_duplicates("truth_id").set_index("truth_id")[["truth_Ks", "truth_f", "truth_cc"]]
        per[e["key"]] = td
    common = set.intersection(*[set(td.index) for td in per.values()])
    ref_td = per[ents[0]["key"]]
    union = set.union(*[set(td.index) for td in per.values()])
    if len(common) < len(union):
        miss = []
        for e in ents:
            gone = sorted(union - set(per[e["key"]].index))
            if gone:
                miss.append("%s lacks %d (%s)" % (e["name"], len(gone), ", ".join(gone[:3]) + (" ..." if len(gone) > 3 else "")))
        ck.add("WARN", "common truths", "only %d of %d truths exist in every storm and are used for all storms: %s"
               % (len(common), len(union), "; ".join(miss)))
    else:
        ck.add("PASS", "common truths", "all %d truths exist in every storm" % len(common))
    bad = []
    for tid in sorted(common):
        base = ref_td.loc[tid]
        for e in ents[1:]:
            o = per[e["key"]].loc[tid]
            if not (np.isclose(o["truth_Ks"], base["truth_Ks"]) and np.isclose(o["truth_f"], base["truth_f"])
                    and np.isclose(o["truth_cc"], base["truth_cc"])):
                bad.append("%s (%s)" % (tid, e["name"]))
    if bad:
        ck.add("FAIL", "truth parameters", "the same truth id has different (Ks, f, cc) in different storms: %s" % ", ".join(bad[:4]))
    else:
        ck.add("PASS", "truth parameters", "every truth has the same (Ks, f, cc) in every storm")
    vers = {e["name"]: minor(e["pandas"]) for e in ents if e.get("pandas")}
    if len(set(vers.values())) > 1:
        ck.add("WARN", "pandas versions", "storms were re-scored with different pandas versions (%s). The 5-minute interpolation "
               "differs slightly between pandas 2 and 3, so storm-to-storm differences smaller than that are not meaningful."
               % ", ".join("%s: %s" % kv for kv in vers.items()))
    elif vers:
        ck.add("PASS", "pandas versions", "all storms re-scored with pandas %s" % list(vers.values())[0])
    t_sorted = (ref_td.loc[sorted(common)].reset_index().sort_values(["truth_Ks", "truth_f", "truth_cc"]))
    return list(t_sorted["truth_id"])


# ------------------------------------------------------------------
# Data tensors, aligned to the design order so runs pair across storms
# ------------------------------------------------------------------
def build_entry_data(df, truth_ids, D):
    out = {}
    for s in SERIES:
        sub = df[(df.series == s) & df.truth_id.isin(truth_ids)]
        if sub.empty:
            continue
        base = sub[sub.truth_id == truth_ids[0]].reset_index(drop=True)
        di = design_match(base[P.PARAMS].to_numpy(float), D)[0]
        if di is None:
            sys.exit("Cannot line the runs up with the design (%s). The input checks above say why." % s)
        order = np.argsort(di, kind="stable")
        run_order = [str(r) for r in base["run_id"].to_numpy()[order]]
        X = np.column_stack([base["Ks_mult"].to_numpy(float), np.log10(base["f_RS_abs"].to_numpy(float)),
                             np.log10(base["channelconductivity_mmhr"].to_numpy(float))])[order]
        piv = sub.pivot(index="truth_id", columns="run_id", values="kge_2012")
        K = piv.reindex(index=truth_ids, columns=run_order).to_numpy(float)
        out[s] = {"X": X, "K": K, "run_id": run_order}
    return out


def make_combined(ents, ref_key):
    ref = [e for e in ents if e["key"] == ref_key][0]
    data = {}
    for s in SERIES:
        if all(s in e["data"] for e in ents):
            K = np.mean([e["data"][s]["K"] for e in ents], axis=0)
            data[s] = {"X": ref["data"][s]["X"], "K": K, "run_id": ref["data"][s]["run_id"]}
    return {"key": COMB_KEY, "name": "combined", "scale": None, "data": data}


# ------------------------------------------------------------------
# Resampling machinery. One re-draw / one shuffle is shared by every storm.
# ------------------------------------------------------------------
def point_sets(entries, fracs):
    sets, pt = {}, {}
    for e in entries:
        for s in e["data"]:
            X, K = e["data"][s]["X"], e["data"][s]["K"]
            for frac in fracs:
                idxs = [R.top_index(K[ti], frac) for ti in range(K.shape[0])]
                sets[(e["key"], s, frac)] = idxs
                med = np.array([np.median(X[ix], axis=0) for ix in idxs])
                q = np.array([np.percentile(X[ix], [10, 90], axis=0) for ix in idxs])
                pt[(e["key"], s, frac)] = {"med": med, "q10": q[:, 0, :], "q90": q[:, 1, :]}
    return sets, pt


def run_boot(entries, fracs, B, seed):
    boot = {}
    for si, s in enumerate(SERIES):
        ents = [e for e in entries if s in e["data"]]
        if not ents:
            continue
        n = ents[0]["data"][s]["X"].shape[0]
        idxb = np.random.default_rng([seed, 77, si]).integers(0, n, size=(B, n))
        for e in ents:
            X, K = e["data"][s]["X"], e["data"][s]["K"]
            Xb = X[idxb]
            for frac in fracs:
                arr = np.full((B, K.shape[0], 3), np.nan)
                for ti in range(K.shape[0]):
                    kb = K[ti][idxb]
                    kb = np.where(np.isfinite(kb), kb, -np.inf)
                    m = T.n_top_for(int(np.isfinite(K[ti]).sum()), frac)
                    order = np.argsort(-kb, axis=1, kind="stable")[:, :m]
                    v = np.take_along_axis(Xb, order[:, :, None], axis=1)
                    arr[:, ti, :] = np.median(v, axis=1)
                boot[(e["key"], s, frac)] = arr
    return boot


def run_perm(entries, sets, fracs, Bp, seed):
    null = {}
    for si, s in enumerate(SERIES):
        ents = [e for e in entries if s in e["data"]]
        if not ents:
            continue
        n = ents[0]["data"][s]["X"].shape[0]
        Pm = T.make_perm_matrix(n, Bp, np.random.default_rng([seed, 66, si]))
        for e in ents:
            X = e["data"][s]["X"]
            J = e["data"][s]["K"].shape[0]
            for frac in fracs:
                med = np.full((Bp, J, 3), np.nan)
                for ti in range(J):
                    med[:, ti, :] = np.median(X[Pm[:, sets[(e["key"], s, frac)][ti]]], axis=1)
                null[(e["key"], s, frac)] = med
    return null


# ------------------------------------------------------------------
# Tables
# ------------------------------------------------------------------
def build_tables(entries, ref_key, truths, tx, fracs, sets, pt, boot, null):
    keys, cells, centre = L.find_cells(truths)
    all_ix = list(range(len(truths)))
    SL, SB = {}, {}                        # (entry, series, frac, scope, param) -> (slope, lo, hi) / bootstrap slopes
    pool_rows, cell_rows, truth_rows, diff_rows = [], [], [], []
    for e in entries:
        for s in e["data"]:
            K = e["data"][s]["K"]
            for frac in fracs:
                P_, bt, nl = pt[(e["key"], s, frac)], boot[(e["key"], s, frac)], null[(e["key"], s, frac)]
                base = {"entry": e["key"], "entry_name": e["name"], "rain_scale": e["scale"], "series": s, "top_frac": frac}
                for ti, t in enumerate(truths):
                    idx = sets[(e["key"], s, frac)][ti]
                    row = dict(base)
                    row.update({"truth_id": t["id"], "truth_Ks": t["Ks"], "truth_f": t["f"], "truth_cc": t["cc"],
                                "n_top": int(len(idx)), "best_kge_2012": float(np.nanmax(K[ti]))})
                    for p in PAR:
                        col = PCOL[p]
                        med, q10, q90 = P_["med"][ti, col], P_["q10"][ti, col], P_["q90"][ti, col]
                        lo, hi = L.ci(bt[:, ti, col])
                        if p == "Ks":
                            row.update({"median_Ks": med, "q10_Ks": q10, "q90_Ks": q90, "median_Ks_ci_lo": lo,
                                        "median_Ks_ci_hi": hi, "Ks_bias": med - tx[ti, col]})
                        else:
                            row.update({"median_%s" % p: 10 ** med, "q10_%s" % p: 10 ** q10, "q90_%s" % p: 10 ** q90,
                                        "median_%s_ci_lo" % p: 10 ** lo, "median_%s_ci_hi" % p: 10 ** hi,
                                        "%s_factor" % p: 10 ** (med - tx[ti, col])})
                    truth_rows.append(row)
                for p in PAR:
                    st, sl_b = L.track(all_ix, p, tx, P_["med"], bt, nl)
                    if st is None:
                        continue
                    r = dict(base)
                    r.update({"param": p})
                    r.update(st)
                    pool_rows.append(r)
                    SL[(e["key"], s, frac, "ALL", p)] = (st["slope"], st["slope_ci_lo"], st["slope_ci_hi"])
                    SB[(e["key"], s, frac, "ALL", p)] = sl_b
                for k in keys:
                    st, sl_b = L.track(cells[k], "cc", tx, P_["med"], bt, nl)
                    if st is None:
                        continue
                    r = dict(base)
                    r.update({"cell_Ks": k[0], "cell_f": k[1], "cell": L.cell_name(k)})
                    r.update(st)
                    cell_rows.append(r)
                    SL[(e["key"], s, frac, L.cell_name(k), "cc")] = (st["slope"], st["slope_ci_lo"], st["slope_ci_hi"])
                    SB[(e["key"], s, frac, L.cell_name(k), "cc")] = sl_b

    def add(kind, ent, ref, s, frac, scope, p, a, b, d_pt, d_b):
        lo, hi = L.ci(d_b)
        diff_rows.append({"kind": kind, "entry": ent, "reference": ref, "series": s, "top_frac": frac, "scope": scope,
                          "param": p, "value_entry": a, "value_reference": b, "diff": d_pt, "ci_lo": lo, "ci_hi": hi,
                          "excludes_0": bool(lo > 0 or hi < 0)})

    scopes = [("ALL", p) for p in PAR] + [(L.cell_name(k), "cc") for k in keys]
    for e in entries:
        for s in e["data"]:
            for frac in fracs:
                for scope, p in scopes:
                    ka = (e["key"], s, frac, scope, p)
                    if ka not in SL:
                        continue
                    kr = (ref_key, s, frac, scope, p)
                    if e["key"] != ref_key and kr in SL:                    # paired difference from the reference storm
                        add("slope_vs_ref", e["key"], ref_key, s, frac, scope, p, SL[ka][0], SL[kr][0],
                            SL[ka][0] - SL[kr][0], SB[ka] - SB[kr])
                    if s == SER_ON:                                         # ON minus control (independent re-draws)
                        ko = (e["key"], SER_OFF, frac, scope, p)
                        if ko in SL:
                            add("on_minus_control", e["key"], "control", "ON-OFF", frac, scope, p, SL[ka][0], SL[ko][0],
                                SL[ka][0] - SL[ko][0], SB[ka] - SB[ko])
                if e["key"] != ref_key and (ref_key, s, frac) in pt:        # how far the best-run median moves between storms
                    for p in PAR:
                        col = PCOL[p]
                        d_pt = float(np.mean(np.abs(pt[(e["key"], s, frac)]["med"][:, col] - pt[(ref_key, s, frac)]["med"][:, col])))
                        d_b = np.mean(np.abs(boot[(e["key"], s, frac)][:, :, col] - boot[(ref_key, s, frac)][:, :, col]), axis=1)
                        add("shift_vs_ref", e["key"], ref_key, s, frac, "ALL", p, np.nan, np.nan, d_pt, d_b)
    return (pd.DataFrame(pool_rows), pd.DataFrame(cell_rows), pd.DataFrame(truth_rows), pd.DataFrame(diff_rows), SL, SB, keys)


# ------------------------------------------------------------------
# Printing
# ------------------------------------------------------------------
def fnum(v, spec):
    return "n/a" if (v is None or not np.isfinite(v)) else format(v, spec)


def slope_txt(SL, key, s, frac, scope, p):
    v = SL.get((key, s, frac, scope, p))
    return "n/a" if v is None else L.iv(v[0], v[1], v[2])


def dtxt(D, kind, key, s, frac, scope, p):
    r = D.get((kind, key, s, frac, scope, p))
    if r is None:
        return "n/a"
    return "%+5.2f [%+5.2f, %+5.2f]%s" % (r["diff"], r["ci_lo"], r["ci_hi"], " *" if r["excludes_0"] else "  ")


def shift_txt(D, key, s, frac, p):
    r = D.get(("shift_vs_ref", key, s, frac, "ALL", p))
    if r is None:
        return "n/a"
    return "%.2f" % r["diff"] if p == "Ks" else "x%.2f" % (10 ** r["diff"])


def print_inputs(entries, ref_key):
    print("\nStorms compared (rain multiplier; control = same 250 points with channel loss OFF):")
    for e in entries:
        if e["key"] == COMB_KEY:
            continue
        print("  %-6s %-9s %s%s" % (e["name"], e["key"], "cc ON + control" if SER_OFF in e["data"] else "cc ON only (no control)",
                                    "   <- reference for the differences" if e["key"] == ref_key else ""))


def print_pooled(pool_df, SL, D, entries, fracs, main_frac, truths):
    print("\n" + "=" * 150)
    print("1. WHAT EACH STORM GIVES: pooled tracking slope over all %d truths  (best runs by KGE_2012; slope of best-run parameter "
          "against true parameter)" % len(truths))
    print("=" * 150)
    print("  1 = the best runs follow the truth one-for-one; 0 = they ignore it. The control (cc inert) shows what 0 looks like by "
          "chance. 'combined' = each run's KGE averaged over the storms.")
    print("\n  cc, all cuts")
    print("  %-4s %-9s | %-24s %7s | %-9s %9s | %-24s | %-26s"
          % ("cut", "storm", "ON slope [95% interval]", "p", "typ.error", "no-info", "control slope [95%]", "ON minus control [95%]"))
    for frac in fracs:
        for e in entries:
            g = pool_df[(pool_df.entry == e["key"]) & (pool_df.series == SER_ON) & (pool_df.top_frac == frac) & (pool_df.param == "cc")]
            if g.empty:
                continue
            r = g.iloc[0]
            print("  %-4s %-9s | %-24s %7s | %-9s %9s | %-24s | %-26s"
                  % (L.cut_label(frac), e["name"], slope_txt(SL, e["key"], SER_ON, frac, "ALL", "cc"), T.fmt_p(r.p_slope),
                     L.param_err_fmt("cc", r.typical_error), L.param_err_fmt("cc", r.null_typical_error),
                     slope_txt(SL, e["key"], SER_OFF, frac, "ALL", "cc"),
                     dtxt(D, "on_minus_control", e["key"], "ON-OFF", frac, "ALL", "cc")))
        print("")
    for p in ("Ks", "f"):
        print("  %s, best %s of runs" % (p, L.cut_label(main_frac)))
        print("  %-9s | %-24s %7s | %-9s %9s | %-24s" % ("storm", "ON slope [95% interval]", "p", "typ.error", "no-info",
                                                         "control slope [95%]"))
        for e in entries:
            g = pool_df[(pool_df.entry == e["key"]) & (pool_df.series == SER_ON) & (pool_df.top_frac == main_frac) & (pool_df.param == p)]
            if g.empty:
                continue
            r = g.iloc[0]
            print("  %-9s | %-24s %7s | %-9s %9s | %-24s"
                  % (e["name"], slope_txt(SL, e["key"], SER_ON, main_frac, "ALL", p), T.fmt_p(r.p_slope),
                     L.param_err_fmt(p, r.typical_error), L.param_err_fmt(p, r.null_typical_error),
                     slope_txt(SL, e["key"], SER_OFF, main_frac, "ALL", p)))
        print("")
    print("  typ.error: typical distance between the best runs' median and the true value (a factor for cc and f, a difference for Ks); "
          "no-info: the same for shuffled runs. p: could that slope arise if the parameter were shuffled among the runs?")
    print("  ON minus control: how much of the cc tracking is more than chance.  * = its 95% interval excludes 0.")
    print("  The control's Ks and f slopes show what Ks and f alone can do when the model has no channel loss to explain.")


def print_diffs(D, entries, ref_key, fracs, main_frac, ref_name):
    others = [e for e in entries if e["key"] != ref_key]
    print("\n" + "=" * 150)
    print("2. DOES IT CHANGE WITH STORM SIZE?  (slope in the storm MINUS slope in the %s storm; paired bootstrap; * = the 95%% interval "
          "excludes 0)" % ref_name)
    print("=" * 150)
    print("  The control columns are the yardstick: they show how much a slope moves between storms when there is nothing to find (cc) "
          "or when channel loss is absent (Ks, f).")
    for p, cuts in (("cc", fracs), ("Ks", [main_frac]), ("f", [main_frac])):
        print("\n  %s" % p)
        print("  %-4s %-9s | %-26s | %-26s | %-10s %-10s"
              % ("cut", "storm", "ON: slope difference", "control: slope difference", "ON shift", "ctl shift"))
        for frac in cuts:
            for e in others:
                print("  %-4s %-9s | %-26s | %-26s | %-10s %-10s"
                      % (L.cut_label(frac), e["name"], dtxt(D, "slope_vs_ref", e["key"], SER_ON, frac, "ALL", p),
                         dtxt(D, "slope_vs_ref", e["key"], SER_OFF, frac, "ALL", p),
                         shift_txt(D, e["key"], SER_ON, frac, p), shift_txt(D, e["key"], SER_OFF, frac, p)))
    print("\n  shift: how far the best runs' median moves between the storm and the %s storm, averaged over the truths (x1.00 = not at "
          "all; for Ks a plain difference). Shifts are always >= 0, so compare ON with the control, not with 0." % ref_name)
    print("  'combined' minus %s: a positive slope difference means using all the storms together tracks the truth better than the "
          "%s storm alone." % (ref_name, ref_name))


def print_cells(SL, D, entries, ref_key, keys, frac, ref_name):
    print("\n" + "=" * 150)
    print("3. DOES cc TRACK THE TRUTH IN EVERY (Ks, f) CELL, IN EVERY STORM?  cc slope over each cell's cc values, best %s of runs, "
          "channel loss ON" % L.cut_label(frac))
    print("=" * 150)
    head = "  %-16s |" % "cell"
    for e in entries:
        head += " %-24s|" % e["name"]
    print(head)
    for k in keys:
        nm = L.cell_name(k)
        line = "  %-16s |" % nm
        for e in entries:
            line += " %-24s|" % slope_txt(SL, e["key"], SER_ON, frac, nm, "cc")
        print(line)
    for lab, s in (("ALL TRUTHS, ON", SER_ON), ("ALL TRUTHS, control", SER_OFF)):
        line = "  %-16s |" % lab[:16]
        for e in entries:
            line += " %-24s|" % slope_txt(SL, e["key"], s, frac, "ALL", "cc")
        print(line)
    print("")
    others = [e for e in entries if e["key"] != ref_key]
    for e in others:
        flagged = []
        for k in keys:
            r = D.get(("slope_vs_ref", e["key"], SER_ON, frac, L.cell_name(k), "cc"))
            if r is not None and r["excludes_0"]:
                flagged.append("%s (%+.2f [%+.2f, %+.2f])" % (L.cell_name(k), r["diff"], r["ci_lo"], r["ci_hi"]))
        n_have = sum(1 for k in keys if ("slope_vs_ref", e["key"], SER_ON, frac, L.cell_name(k), "cc") in D)
        print("  %-9s cells whose cc slope differs detectably from the %s storm: %d of %d%s"
              % (e["name"], ref_name, len(flagged), n_have, (":  " + "; ".join(flagged)) if flagged else ""))
    print("  A cell's slope rests on only a few cc values, so it is rough. ! = the slope lies outside its own bootstrap interval "
          "(few, nearly tied best runs). With about %d such comparisons per storm, one or two intervals excluding 0 can occur by "
          "chance." % len(keys))


def summary_lines(SL, D, pool_df, entries, ref_key, keys, frac, ref_name):
    out = []
    singles = [e for e in entries if e["key"] not in (COMB_KEY,)]
    parts = []
    for e in singles:
        t = "%s %s" % (e["name"], slope_txt(SL, e["key"], SER_ON, frac, "ALL", "cc").strip().rstrip("!"))
        if (e["key"], SER_OFF, frac, "ALL", "cc") in SL:
            c = SL[(e["key"], SER_OFF, frac, "ALL", "cc")]
            t += " (control %+.2f)" % c[0]
        parts.append(t)
    out.append("Pooled cc tracking slope, best %s of runs, channel loss ON: %s.  (0 = ignores the truth, 1 = follows it one-for-one.)"
               % (L.cut_label(frac), "; ".join(parts)))
    nm_missing = [e["name"] for e in singles if (e["key"], SER_OFF, frac, "ALL", "cc") not in SL]
    if nm_missing:
        out.append("No control for %s, so there is no chance baseline for that storm." % ", ".join(nm_missing))
    ex = [e for e in singles if D.get(("on_minus_control", e["key"], "ON-OFF", frac, "ALL", "cc"))]
    if ex:
        pos = [e["name"] for e in ex if D[("on_minus_control", e["key"], "ON-OFF", frac, "ALL", "cc")]["ci_lo"] > 0]
        out.append("ON minus control (cc): the 95%% interval is entirely above 0 in %d of %d storms with a control (%s)."
                   % (len(pos), len(ex), ", ".join(pos) if pos else "none"))
    for e in [x for x in entries if x["key"] != ref_key]:
        r = D.get(("slope_vs_ref", e["key"], SER_ON, frac, "ALL", "cc"))
        if r is None:
            continue
        rc = D.get(("slope_vs_ref", e["key"], SER_OFF, frac, "ALL", "cc"))
        ctl = (" (the control moves by %+.2f [%+.2f, %+.2f] for scale)" % (rc["diff"], rc["ci_lo"], rc["ci_hi"])) if rc else ""
        what = "combined" if e["key"] == COMB_KEY else "the %s storm" % e["name"]
        out.append("cc slope, %s minus the %s storm: %+.2f [%+.2f, %+.2f]%s -- %s."
                   % (what, ref_name, r["diff"], r["ci_lo"], r["ci_hi"], ctl,
                      "the interval excludes 0" if r["excludes_0"] else "the interval includes 0: no detectable difference"))
    out.append("'No detectable difference' means the data cannot tell the storms apart at this sample size (250 runs); it does not prove "
               "they are equal.")
    for p in ("Ks", "f"):
        parts = []
        for e in entries:
            v = SL.get((e["key"], SER_ON, frac, "ALL", p))
            if v is not None:
                parts.append("%s %+.2f [%+.2f, %+.2f]" % (e["name"], v[0], v[1], v[2]))
        if parts:
            out.append("Pooled %s tracking slope: %s." % (p, "; ".join(parts)))
    for e in singles:
        pos = [k for k in keys if (e["key"], SER_ON, frac, L.cell_name(k), "cc") in SL
               and SL[(e["key"], SER_ON, frac, L.cell_name(k), "cc")][1] > 0]
        have = [k for k in keys if (e["key"], SER_ON, frac, L.cell_name(k), "cc") in SL]
        if have:
            out.append("%s: %d of %d cells have a cc slope whose 95%% interval is entirely above 0." % (e["name"], len(pos), len(have)))
    return out


def cross_check_earlier(summary_dir, pool_df, fracs):
    """The reference storm here should give the same point estimates as the earlier analysis of that storm."""
    path = summary_dir / EARLIER_DIR / "location_pooled_110.csv"
    if not path.exists():
        return None, "earlier analysis not found (%s/%s): nothing to compare the %s numbers with" % (EARLIER_DIR, path.name, REF_KEY)
    old = pd.read_csv(path)
    worst, n = 0.0, 0
    for _, o in old.iterrows():
        frac = None
        for f_ in fracs:
            if abs(f_ - float(o["top_frac"])) < 1e-9:
                frac = f_
        if frac is None:
            continue
        g = pool_df[(pool_df.entry == REF_KEY) & (pool_df.series == o["series"]) & (pool_df.top_frac == frac) & (pool_df.param == o["param"])]
        if g.empty:
            continue
        worst = max(worst, abs(float(g.iloc[0]["slope"]) - float(o["slope"])))
        n += 1
    if n == 0:
        return None, "earlier analysis found but no overlapping cuts to compare"
    return worst, "the x1.0 slopes (cc, Ks, f; ON and control; %d values) compared with the earlier analyze_truth_location_110.py " \
                  "results: largest difference %.1e%s" % (n, worst, " (identical to rounding)" if worst <= 1e-9 else "")


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def _mix(c1, c2, t):
    import matplotlib.colors as mc
    a, b = np.array(mc.to_rgb(c1)), np.array(mc.to_rgb(c2))
    return mc.to_hex(a * (1 - t) + b * t)


def _ramp(anchors, n):
    if n == 1:
        return [anchors[1]]
    pos = np.linspace(0, 2, n)
    out = []
    for v in pos:
        i = min(int(np.floor(v)), 1)
        out.append(_mix(anchors[i], anchors[i + 1], v - i))
    return out


def storm_styles(entries):
    singles = [e for e in entries if e["key"] != COMB_KEY]
    on = _ramp([_mix(P.C_CC, "#ffffff", 0.35), P.C_CC, _mix(P.C_CC, "#000000", 0.40)], len(singles))
    off = _ramp(["#b9b8b3", "#8a8984", "#55544f"], len(singles))
    sty = {}
    for i, e in enumerate(singles):
        sty[e["key"]] = {SER_ON: on[i], SER_OFF: off[i], "marker": SINGLE_MARKERS[i % len(SINGLE_MARKERS)]}
    sty[COMB_KEY] = {SER_ON: INK, SER_OFF: INK, "marker": "D"}
    return sty


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def fig_cell_slopes(SL, entries, keys, frac, path):
    plt = _plt()
    from matplotlib.lines import Line2D
    sty = storm_styles(entries)
    labels = [L.cell_name(k) for k in keys] + ["all truths pooled"]
    scopes = [L.cell_name(k) for k in keys] + ["ALL"]
    ny, ne = len(labels), len(entries)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 0.5 * ny + 2.9), facecolor=SURFACE, sharey=True, sharex=True)
    offs = np.linspace(-0.3, 0.3, ne) if ne > 1 else [0.0]
    allv = []
    for ax, s in zip(axes, SERIES):
        T._style(ax)
        ax.axvline(0.0, color=INK, linestyle="--", linewidth=0.9, zorder=1)
        ax.axvline(1.0, color=NULLC, linestyle=":", linewidth=1.1, zorder=1)
        for ei, e in enumerate(entries):
            c = sty[e["key"]][s]
            mk = sty[e["key"]]["marker"]
            for yi, sc in enumerate(scopes):
                v = SL.get((e["key"], s, frac, sc, "cc"))
                if v is None:
                    continue
                y = yi + offs[ei]
                allv.extend([v[1], v[2]])
                ax.hlines(y, v[1], v[2], color=c, linewidth=1.7, zorder=3)
                ax.plot([v[0]], [y], linestyle="none", marker=mk, markersize=5.5, markerfacecolor=c if s == SER_ON else SURFACE,
                        markeredgecolor=c, markeredgewidth=1.4, zorder=4)
        ax.set_yticks(range(ny))
        ax.set_yticklabels(labels, fontsize=8)
        ax.axhline(ny - 1.5, color=GRID, linewidth=1.0, zorder=1)
        ax.set_title("channel loss ON" if s == SER_ON else "control (cc inert: what chance gives)", color=INK, fontsize=9.5, loc="left")
        ax.set_xlabel("cc tracking slope (95% interval)", color=INK2, fontsize=9)
    axes[0].invert_yaxis()                      # the two panels share their y axis: invert it once, not once per panel
    if allv:
        lo, hi = min(-0.3, min(allv) - 0.05), max(1.3, max(allv) + 0.05)
        axes[0].set_xlim(lo, hi)
    handles = [Line2D([0], [0], color=sty[e["key"]][SER_ON] if e["key"] != COMB_KEY else INK, marker=sty[e["key"]]["marker"],
                      markersize=5.5, linewidth=1.7, label=e["name"] + (" (all storms averaged)" if e["key"] == COMB_KEY else ""))
               for e in entries]
    handles += [Line2D([0], [0], color=INK, linestyle="--", linewidth=0.9, label="0 = ignores the true cc"),
                Line2D([0], [0], color=NULLC, linestyle=":", linewidth=1.1, label="1 = follows it one-for-one")]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 6), frameon=False, fontsize=8, labelcolor=INK2)
    fig.suptitle("Does cc tracking hold in every cell and every storm?  (best %s of runs; storms in rain order, top to bottom within a row; "
                 "control panel in greys)" % L.cut_label(frac), color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.1, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_pooled(SL, entries, fracs, path):
    plt = _plt()
    from matplotlib.lines import Line2D
    singles = [e for e in entries if e["key"] != COMB_KEY]
    comb = [e for e in entries if e["key"] == COMB_KEY]
    xs = {e["key"]: float(i) for i, e in enumerate(singles)}
    if comb:
        xs[COMB_KEY] = len(singles) + 0.7
    nr = len(fracs)
    fig, axes = plt.subplots(nr, 3, figsize=(10.8, 2.75 * nr + 1.4), facecolor=SURFACE, sharex=True, sharey=True, squeeze=False)
    for ri, frac in enumerate(fracs):
        for ci_, p in enumerate(("cc", "Ks", "f")):
            ax = axes[ri][ci_]
            T._style(ax)
            ax.axhline(0.0, color=INK, linestyle="--", linewidth=0.9, zorder=1)
            ax.axhline(1.0, color=NULLC, linestyle=":", linewidth=1.1, zorder=1)
            for s, col, off_, mk, filled in ((SER_OFF, NULLC, 0.13, "s", False), (SER_ON, PHUE[p], -0.13, "o", True)):
                xx, yy = [], []
                for e in entries:
                    v = SL.get((e["key"], s, frac, "ALL", p))
                    if v is None:
                        continue
                    x = xs[e["key"]] + off_
                    ax.vlines(x, v[1], v[2], color=col, linewidth=2.0 if s == SER_ON else 1.4, zorder=3)
                    ax.plot([x], [v[0]], linestyle="none", marker=mk, markersize=6 if s == SER_ON else 5.5,
                            markerfacecolor=col if filled else SURFACE, markeredgecolor=col, markeredgewidth=1.5, zorder=4)
                    if e["key"] != COMB_KEY:
                        xx.append(x)
                        yy.append(v[0])
                if len(xx) > 1:
                    ax.plot(xx, yy, color=col, linewidth=0.9, alpha=0.6, linestyle="-" if s == SER_ON else "--", zorder=2)
            ax.set_xticks([xs[e["key"]] for e in entries])
            ax.set_xticklabels([e["name"] for e in entries], fontsize=8)
            ax.set_xlim(-0.6, max(xs.values()) + 0.6)
            if ri == 0:
                ax.set_title(p, color=INK, fontsize=10, loc="left", fontweight="bold")
            if ci_ == 0:
                ax.set_ylabel("best %s of runs\npooled tracking slope" % L.cut_label(frac), color=INK2, fontsize=8.5)
            if ri == nr - 1:
                ax.set_xlabel("rain multiplier", color=INK2, fontsize=8.5)
    handles = [Line2D([0], [0], color=INK, marker="o", markersize=6, linestyle="none", label="channel loss ON (colour = parameter)"),
               Line2D([0], [0], color=NULLC, marker="s", markersize=5.5, markerfacecolor=SURFACE, markeredgewidth=1.5,
                      linestyle="none", label="control (cc inert)"),
               Line2D([0], [0], color=INK, linestyle="--", linewidth=0.9, label="0 = ignores the truth"),
               Line2D([0], [0], color=NULLC, linestyle=":", linewidth=1.1, label="1 = follows it one-for-one")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8, labelcolor=INK2)
    fig.suptitle("Does recovery of cc, Ks and f depend on storm size?  (pooled over all truths; bars are 95% intervals; "
                 "'combined' = runs scored on all storms together)", color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_recovery(pt, SL, entries, truths, tx, keys_cells, cells, frac, path):
    plt = _plt()
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ne = len(entries)
    fig, axes = plt.subplots(1, ne, figsize=(3.1 * ne + 0.9, 4.7), facecolor=SURFACE, sharex=True, sharey=True, squeeze=False)
    lo, hi = 10 ** R.LOG_LO, 10 ** R.LOG_HI
    lim = (lo * 0.75, hi * 1.3)
    cc_all = sorted({t["cc"] for t in truths})
    rank = {}
    for ri, k in enumerate(keys_cells):
        for ti in cells[k]:
            rank[ti] = ri
    nk = max(1, len(keys_cells) - 1)
    sty = storm_styles(entries)
    for ax, e in zip(axes[0], entries):
        T._style(ax)
        ax.axhspan(10 ** (R.LOG_LO + 0.1 * (R.LOG_HI - R.LOG_LO)), 10 ** (R.LOG_LO + 0.9 * (R.LOG_HI - R.LOG_LO)),
                   color=NULLC, alpha=0.16, linewidth=0, zorder=1)
        ax.axhline(10 ** R.CENTER_LOG, color=NULLC, linestyle=":", linewidth=1.1, zorder=2)
        ax.plot(lim, lim, color=INK, linestyle="--", linewidth=1.0, zorder=2)
        for s, col, sh, mk, filled in ((SER_OFF, NULLC, +0.012, "s", False), (SER_ON, sty[e["key"]][SER_ON] if e["key"] != COMB_KEY
                                                                              else INK, -0.012, "o", True)):
            if (e["key"], s, frac) not in pt:
                continue
            med = pt[(e["key"], s, frac)]["med"][:, 2]
            xx = np.array([10 ** (tx[ti, 2] + sh + 0.03 * (rank[ti] / nk - 0.5)) for ti in range(len(truths))])
            ax.plot(xx, 10 ** med, linestyle="none", marker=mk, markersize=5.2, markerfacecolor=col if filled else SURFACE,
                    markeredgecolor=col, markeredgewidth=1.3, alpha=0.95, zorder=3 if s == SER_ON else 2)
        L._setup_log_axis(ax, lim)
        ax.set_xticks(cc_all)
        ax.set_xticklabels(["%g" % v for v in cc_all], fontsize=7.5)
        ax.set_yticks([30, 100, 300, 1000])
        ax.set_xlabel("true cc (mm/hr)", color=INK2, fontsize=8.5)
        v = SL.get((e["key"], SER_ON, frac, "ALL", "cc"))
        c = SL.get((e["key"], SER_OFF, frac, "ALL", "cc"))
        txt = "ON slope %+.2f" % v[0] if v else ""
        if c:
            txt += "\ncontrol %+.2f" % c[0]
        ax.set_title(e["name"] + ("  (all storms averaged)" if e["key"] == COMB_KEY else ""), color=INK, fontsize=9.5, loc="left")
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, ha="left", va="top", fontsize=8, color=INK2)
    axes[0][0].set_ylabel("cc of the best runs, median (mm/hr)", color=INK2, fontsize=8.5)
    handles = [Line2D([0], [0], color=P.C_CC, marker="o", markersize=5.2, linestyle="none", label="channel loss ON (one dot per truth)"),
               Line2D([0], [0], color=NULLC, marker="s", markersize=5.2, markerfacecolor=SURFACE, markeredgewidth=1.3,
                      linestyle="none", label="control (cc inert)"),
               Line2D([0], [0], color=INK, linestyle="--", linewidth=1.0, label="perfect recovery"),
               Patch(facecolor=NULLC, alpha=0.16, label="no information: middle 80% of the sampled range")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8, labelcolor=INK2)
    fig.suptitle("Where do the best runs put cc, storm by storm?  (best %s of runs; 27 truths per panel, 9 per true-cc column)"
                 % L.cut_label(frac), color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.1, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110 -- compare the storms: does the cc / Ks / f answer depend on storm size?")
    ap.add_argument("--storms", type=float, nargs="+", default=[0.8, 1.25],
                    help="Rain multipliers to compare with the original storm (default 0.8 1.25). Storms whose re-scored files "
                         "do not exist yet are skipped with a note.")
    ap.add_argument("--top_fracs", type=float, nargs="+", default=[0.20, 0.10, 0.05],
                    help="Top fractions by KGE_2012; the first is the main one (default 0.20 0.10 0.05)")
    ap.add_argument("--boot", type=int, default=1000, help="Paired-bootstrap re-draws (default 1000)")
    ap.add_argument("--perm", type=int, default=2000, help="Permutation shuffles (default 2000)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--expect_n", type=int, default=250)
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--continue_on_fail", action="store_true")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()
    fracs = args.top_fracs
    main_frac = fracs[0]
    if args.boot < 100 or args.perm < 100:
        sys.exit("--boot and --perm should be at least 100.")

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    out_dir = summary_dir / OUT_DIRNAME
    print("\n" + "=" * 78 + "\nSeries 110 -- storm comparison: does the cc / Ks / f answer depend on storm size?\n" + "=" * 78)
    print("python %s, pandas %s, numpy %s" % (sys.version.split()[0], pd.__version__, np.__version__))
    print("Reading from: %s" % summary_dir)

    # ---- load every storm that exists ---------------------------------------------------------------
    scales = [REF_SCALE] + sorted(set(float(s) for s in args.storms if abs(float(s) - REF_SCALE) > 1e-9))
    ents, notes = [], []
    for sc in scales:
        pth = entry_paths(summary_dir, sc)
        name = scale_name(sc)
        if not pth["long"].exists():
            if abs(sc - REF_SCALE) < 1e-9:
                sys.exit("Required file not found: %s\nThe original storm's re-scored file is needed as the reference." % pth["long"])
            notes.append("%s (%s): %s not found, so this storm is skipped. Run: python rescore_storm_110.py --rain_scale %g"
                         % (name, pth["label"], pth["long"].name, sc))
            continue
        ents.append({"key": pth["label"], "name": name, "scale": sc, "paths": pth, "df": load_long(pth["long"], pth["label"], sc)})
    for n in notes:
        print("NOTE: " + n)
    if len(ents) < 2:
        sys.exit("\nNo second storm to compare with. Run rescore_storm_110.py for at least one of the new storms first.")

    D = P.generate_lhs_samples(args.expect_n, P.LHS_PARAMS, P.DESIGN_SEED)[P.PARAMS].to_numpy(float)
    ck = P.Checks()
    for e in ents:
        check_entry(ck, e, D, args.expect_n)
    truth_ids = check_cross(ck, ents)
    P.print_checks(ck)
    if ck.n("FAIL") and not args.continue_on_fail:
        sys.exit("\nSTOPPED: an input check failed. Nothing was analysed or written. (--continue_on_fail to override.)")
    if len(truth_ids) < 3:
        sys.exit("\nFewer than 3 truths are common to all storms; cannot compute tracking slopes.")

    # ---- tensors --------------------------------------------------------------------------------------
    for e in ents:
        e["data"] = build_entry_data(e["df"], truth_ids, D)
    ref_key = REF_KEY
    ref_name = scale_name(REF_SCALE)
    entries = list(sorted(ents, key=lambda e: e["scale"]))
    if len(entries) >= 2:
        entries.append(make_combined(entries, ref_key))
        if SER_OFF not in entries[-1]["data"]:
            print("NOTE: not every storm has a control, so there is no control for the 'combined' score.")
    td = ents[0]["df"].drop_duplicates("truth_id").set_index("truth_id")
    truths = [{"id": t, "Ks": float(td.loc[t, "truth_Ks"]), "f": float(td.loc[t, "truth_f"]), "cc": float(td.loc[t, "truth_cc"])}
              for t in truth_ids]
    tx = L.truth_coords(truths)
    keys, cells, centre = L.find_cells(truths)
    print_inputs(entries, ref_key)
    print("\n%d truth(s) in %d cell(s), %d runs per series, cuts %s, %d shuffles, %d paired-bootstrap re-draws"
          % (len(truths), len(keys), args.expect_n, [L.cut_label(f) for f in fracs], args.perm, args.boot))

    print("\nRunning the shuffles and bootstrap ...", flush=True)
    sets, pt = point_sets(entries, fracs)
    boot = run_boot(entries, fracs, args.boot, args.seed)
    null = run_perm(entries, sets, fracs, args.perm, args.seed)
    pool_df, cell_df, truth_df, diff_df, SL, SB, keys = build_tables(entries, ref_key, truths, tx, fracs, sets, pt, boot, null)
    D_ = {}
    for r in diff_df.to_dict("records"):
        D_[(r["kind"], r["entry"], r["series"], r["top_frac"], r["scope"], r["param"])] = r

    worst, msg = cross_check_earlier(summary_dir, pool_df, fracs)
    print("\nCross-check: " + msg)
    if worst is not None and worst > 1e-9:
        print("  WARN: not identical. Point estimates should be the same; a difference means the runs were lined up differently "
              "or best runs were tied. Treat the x1.0 numbers here with caution until this is explained.")

    print_pooled(pool_df, SL, D_, entries, fracs, main_frac, truths)
    print_diffs(D_, entries, ref_key, fracs, main_frac, ref_name)
    print_cells(SL, D_, entries, ref_key, keys, main_frac, ref_name)
    print("\n" + "=" * 150)
    print("4. WHAT THE NUMBERS SAY   (best %s of runs; descriptive, assembled from the tables above)" % L.cut_label(main_frac))
    print("=" * 150)
    for line in summary_lines(SL, D_, pool_df, entries, ref_key, keys, main_frac, ref_name):
        print("  - " + line)
    print("  Caveat: the storms share the same 250 sweep points, the same routing and the same rainfall pattern (only scaled), and the "
          "truths are noise-free. This shows how much the answer depends on storm SIZE, not whether it would hold for a different "
          "storm shape or for real, noisy data.")

    # ---- save ----------------------------------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ck.rows).to_csv(out_dir / "storm_compare_input_checks_110.csv", index=False)
    pool_df.to_csv(out_dir / "storm_compare_pooled_110.csv", index=False)
    cell_df.to_csv(out_dir / "storm_compare_cells_110.csv", index=False)
    truth_df.to_csv(out_dir / "storm_compare_per_truth_110.csv", index=False)
    diff_df.to_csv(out_dir / "storm_compare_diffs_110.csv", index=False)
    if not args.no_plots:
        for name, fn in (
                ("fig_stormcmp_cell_slopes_110.png", lambda p: fig_cell_slopes(SL, entries, keys, main_frac, p)),
                ("fig_stormcmp_pooled_slopes_110.png", lambda p: fig_pooled(SL, entries, fracs, p)),
                ("fig_stormcmp_cc_recovery_110.png", lambda p: fig_recovery(pt, SL, entries, truths, tx, keys, cells, main_frac, p))):
            try:
                fn(out_dir / name)
            except Exception as ex:
                print("\n  (figure %s skipped: %s -- the CSVs are saved regardless)" % (name, ex))

    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__,
            "settings": {k: (v if not isinstance(v, Path) else str(v)) for k, v in vars(args).items()},
            "storms": [{"key": e["key"], "name": e["name"], "rain_scale": e["scale"], "has_control": SER_OFF in e["data"],
                        "long_file": str(e["paths"]["long"]), "long_md5": md5_file(e["paths"]["long"]),
                        "rescoring_pandas": e.get("pandas")} for e in entries if e["key"] != COMB_KEY],
            "storms_skipped": notes, "truths": truths, "cells": [list(k) for k in keys], "checks": ck.rows,
            "earlier_analysis_crosscheck": {"worst_abs_slope_difference": worst, "message": msg}}
    (out_dir / "PROVENANCE_storm_compare_110.json").write_text(json.dumps(prov, indent=2, default=str))
    print("\nSaved to: %s" % out_dir)
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results." % ck.n("FAIL"))


if __name__ == "__main__":
    main()
