"""
diagnose_truth_location_110.py
==============================
Series 110 -- truth-location test, part 4 (a closer look): WHY does cc get recovered
well in some (Ks, f) cells and hardly at all in others?

THE QUESTION
------------
analyze_truth_location_110.py found that the best runs follow the true cc in most
cells, but not in every one (for example the cells with true Ks 5 and f 0.012 or 0.02).
Two very different explanations are possible, and they lead to different conclusions:

  (1) cc has little pull on the hydrograph in those cells. The truth barely changes when
      its true cc changes, so no method could tell the cc values apart there.
  (2) cc does matter, but the best runs can hide the difference by trading cc off against
      Ks and f (a different Ks and f paired with a different cc fits about as well). In
      that case cc is only hard to identify while Ks and f are also uncertain.

This script puts numbers on both. It runs NO tRIBS and changes nothing. It reads the
table written by rescore_truth_location_110.py and writes only into its own files in
calibration_work/03_comparisons/summary_tables/location_cc_110/.

THE FOUR LOOKS
--------------
  1. LEVERAGE        How much does the TRUTH itself change when the true cc goes from the
                     lowest to the highest value? (peak and volume of the truth hydrograph)
                     A cell where it barely changes cannot identify cc, whatever the method.
  2. SLICE           Keep only the runs whose Ks and f are near the truth's (within
                     --slice_radius, in units of the sampled box). Among those, do the best
                     runs' cc values follow the true cc? Compared with the all-runs slope, this
                     says whether uncertainty in Ks and f was hiding cc. A paired bootstrap
                     gives an interval for the difference.
  3. RIDGE           Among the best runs, how strongly are Ks, f and cc correlated with each
                     other? A strong correlation means the best runs lie along a trade-off line
                     (a "ridge"). Compared with randomly chosen sets of the same size.
  4. SIDE BY SIDE    The actual best runs (Ks, f, cc, KGE, PBIAS) for the weakest- and the
                     strongest-tracked cell, so you can see the trade-off yourself.

VOCABULARY
-----------
  Cell              One (true Ks, true f) pair; it holds the truths at 3 true cc values.
  Weak / strong     Weak = the cell's all-runs cc slope (best 20% by default) has a 95% bootstrap
                    interval that includes 0. Strong = the interval is entirely above 0. This is
                    just a label for this report, not a test.
  Slice             The runs near the truth in Ks and f (see above). A slice holds roughly 30-50
                    of the 250 runs, so slice results are rougher than all-runs results.
  Spearman rho      A rank correlation between two lists: +1 = they rise together, -1 = one rises
                    as the other falls, 0 = unrelated. It does not assume a straight line.
  Ridge             A line (or narrow band) of parameter combinations that fit about equally
                    well. Along a ridge, a change in one parameter is cancelled by a change in
                    another, so the data cannot tell the combinations apart.
  Paired bootstrap  Re-draw the 250 runs with replacement and redo everything; the SAME re-draw
                    is used for the all-runs slope and the slice slope, so their difference is
                    measured on the same wobble. An interval for the difference that excludes 0
                    means conditioning on Ks and f really changed the answer.
  Control           The same 250 parameter points run with channel loss OFF. cc has no effect
                    there, so its slope shows what chance gives.

HOW TO READ THE RESULTS (rules of thumb, mine)
-----------------------------------------------
  * Weak cells with a SMALL change in the truth (look 1): explanation (1). cc has little leverage
    there; the weakness is physical, not statistical.
  * Weak cells with a LARGE change in the truth, a slice slope clearly above the all-runs slope,
    and a strong Ks-cc or f-cc correlation among the best runs (looks 1-3): explanation (2). cc
    is identifiable only if Ks and f are pinned down as well.
  * Weak cells where everything is weak: cc matters little, or it matters in a way this storm
    cannot reveal. The side-by-side listing (look 4) shows the actual runs.
  * Each slice slope uses about 8-12 runs per truth and 3 truths per cell: treat a single cell
    as a hint, and the pattern across weak versus strong cells as the evidence.

OUTPUTS   (calibration_work/03_comparisons/summary_tables/location_cc_110/)
-----------------------------------------------------------------------------
  diag_cells_110.csv        one row per cell: leverage, all-runs and slice slopes (with intervals),
                            ridge correlations
  diag_per_truth_110.csv    one row per truth: its peak and volume, slice size, slice statistics,
                            ridge correlations and null limits
  diag_best_runs_110.csv    the --list_runs best cc-ON runs for every truth
  fig_diag_slice_leverage_110.png   slopes (all runs vs slice vs control) next to the leverage
  fig_diag_ridge_cells_110.png      the best runs in the Ks-cc plane for the weakest and strongest cells
  PROVENANCE_diag_110.json

USAGE (run from the lab/ directory; keep analyze_truth_location_110.py and its helpers there)
------------------------------------------------------------------------------------------------
    python diagnose_truth_location_110.py
    python diagnose_truth_location_110.py --slice_radius 0.3 --boot 2000
    python diagnose_truth_location_110.py --show_cells 5,0.02 7,0.012      # choose the cells for looks 4

CAVEATS
-------
One storm, 250 runs, routing pinned at truth, noise-free truths: a best case. The slice is a
rough way of holding Ks and f roughly fixed (the runs inside it still differ by up to the
radius), and 3 truths per cell give 3-point slopes. The 27 truths share the same 250 runs, so
p-values and intervals describe the runs, not 27 separate experiments. Everything here is
descriptive; it shows which explanation the numbers favour, it does not prove a mechanism.
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
    import analyze_truth_location_110 as A
except ImportError:
    sys.exit("This script needs analyze_truth_location_110.py (and the helpers it needs: analyze_cc_recovery_110.py, "
             "analyze_cc_tests_110.py, analyze_cc_pca_110.py) in the same folder. Put them all in lab/ and run from there.")

_NEED_A = ["load_long", "truth_list", "truth_coords", "build_tensors", "find_cells", "cell_name", "run_checks", "point_sets",
           "point_stats", "run_boot", "cut_label", "SER_ON", "SER_OFF", "OUT_DIRNAME", "LOC_DIRNAME", "RS_DIRNAME",
           "LONG_NAME", "MIN_TRACK", "STANDARD_KS", "STANDARD_F", "PCOL", "R", "T", "P"]
_miss = [n for n in _NEED_A if not hasattr(A, n)]
if _miss:
    sys.exit("analyze_truth_location_110.py in this folder is an older version (missing: %s). Replace it with the "
             "current copy." % ", ".join(_miss))

R, T, P = A.R, A.T, A.P
SER_ON, SER_OFF = A.SER_ON, A.SER_OFF
MIN_SLICE_N = 12            # fewest runs in a slice for slice statistics
MIN_SLICE_TOP = 5           # fewest best runs kept from a slice
PAIRS = [("Ks", "cc", 0, 2), ("f", "cc", 1, 2), ("Ks", "f", 0, 1)]   # columns of X: Ks, log10 f, log10 cc
PAIR_LABEL = {("Ks", "cc"): "Ks-cc", ("f", "cc"): "f-cc", ("Ks", "f"): "Ks-f"}
EXTRA_COLS = ["obs_volume_m3", "obs_peak_m3s", "pbias_pct"]


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def safe_spearman(a, b, min_n=8):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < min_n or np.ptp(a[ok]) <= 0 or np.ptp(b[ok]) <= 0:
        return np.nan
    return float(stats.spearmanr(a[ok], b[ok])[0])


def rank_corr(a, b):
    """Spearman correlation along the last axis (works on stacks of equal-length lists)."""
    ra = stats.rankdata(a, axis=-1)
    rb = stats.rankdata(b, axis=-1)
    ra = ra - ra.mean(axis=-1, keepdims=True)
    rb = rb - rb.mean(axis=-1, keepdims=True)
    den = np.sqrt((ra ** 2).sum(axis=-1) * (rb ** 2).sum(axis=-1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return (ra * rb).sum(axis=-1) / den


def pct_ci(v, min_frac=0.9):
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    if ok.mean() < min_frac:
        return np.nan, np.nan
    return float(np.percentile(v[ok], 2.5)), float(np.percentile(v[ok], 97.5))


def nz(v):
    return v is not None and np.isfinite(v)


def sgn(v, spec="%+5.2f", blank="   n/a"):
    return spec % v if nz(v) else blank


def interval_str(sl, lo, hi):
    if not nz(sl):
        return "n/a"
    if not (nz(lo) and nz(hi)):
        return "%+5.2f [  n/a  ]" % sl
    return "%+5.2f [%+5.2f, %+5.2f]" % (sl, lo, hi)


def box_units(X2):
    """(n, 2) array of (Ks, log10 f) -> position in the sampled box, 0..1 on each axis."""
    kb, fb = P.BOX["Ks_mult"], P.BOX["f_RS_abs"]
    return np.column_stack([(X2[:, 0] - kb[0]) / (kb[1] - kb[0]),
                            (X2[:, 1] - np.log10(fb[0])) / (np.log10(fb[1]) - np.log10(fb[0]))])


def metric_tensor(df, truths, order, series, col):
    M = np.full((len(truths), len(order)), np.nan)
    for ti, t in enumerate(truths):
        g = df[(df.truth_id == t["id"]) & (df.series == series)].set_index("run_id")[col]
        M[ti] = g.reindex(order).to_numpy(float)
    return M


# ------------------------------------------------------------------
# Look 1: leverage (how much does the truth itself change with cc?)
# ------------------------------------------------------------------
def truth_facts(df, truths):
    facts, bad = {}, []
    for t in truths:
        g = df[df.truth_id == t["id"]]
        v, p = g["obs_volume_m3"].to_numpy(float), g["obs_peak_m3s"].to_numpy(float)
        fine = bool(np.isfinite(v).all() and np.isfinite(p).all() and v.size > 0
                    and np.ptp(v) <= 1e-6 * max(1.0, abs(v).max()) and np.ptp(p) <= 1e-6 * max(1.0, abs(p).max()))
        if not fine:
            bad.append(t["id"])
        facts[t["id"]] = {"vol": float(np.nanmedian(v)) if v.size else np.nan,
                          "peak": float(np.nanmedian(p)) if p.size else np.nan}
    return facts, bad


def leverage(ix, truths, facts):
    v = np.array([facts[truths[i]["id"]]["vol"] for i in ix], float)
    p = np.array([facts[truths[i]["id"]]["peak"] for i in ix], float)
    cc = np.array([truths[i]["cc"] for i in ix], float)
    out = {"vols": v, "peaks": p, "ccs": cc, "vol_change_pct": np.nan, "peak_change_pct": np.nan}
    if len(ix) >= 2 and np.isfinite(v).all() and v[0] > 0 and p[0] > 0:
        out["vol_change_pct"] = float((v[-1] / v[0] - 1.0) * 100.0)
        out["peak_change_pct"] = float((p[-1] / p[0] - 1.0) * 100.0)
    return out


# ------------------------------------------------------------------
# Look 2: the slice (Ks and f held near the truth)
# ------------------------------------------------------------------
def slice_masks(X, tx, radius):
    uX, uT = box_units(X[:, :2]), box_units(tx[:, :2])
    d = np.sqrt(((uX[None, :, :] - uT[:, None, :]) ** 2).sum(axis=2))          # (J, n)
    return d <= radius + 1e-12


def slice_point(data, masks, q, PB):
    """Per series: slice size, median log10 cc of the slice's best runs, Spearman of cc with KGE and with PBIAS."""
    out = {}
    for s in (SER_ON, SER_OFF):
        X, K = data[s]["X"], data[s]["K"]
        J = K.shape[0]
        d = {"n": np.zeros(J, int), "med": np.full(J, np.nan), "best_cc": np.full(J, np.nan),
             "rho_kge": np.full(J, np.nan), "rho_pbias": np.full(J, np.nan)}
        for ti in range(J):
            sel = np.flatnonzero(masks[ti] & np.isfinite(K[ti]))
            d["n"][ti] = len(sel)
            if len(sel) < MIN_SLICE_N:
                continue
            m = max(MIN_SLICE_TOP, int(round(q * len(sel))))
            top = sel[np.argsort(-K[ti][sel], kind="stable")][:m]
            d["med"][ti] = float(np.median(X[top, 2]))
            d["best_cc"][ti] = float(X[top[0], 2])
            d["rho_kge"][ti] = safe_spearman(X[sel, 2], K[ti][sel])
            d["rho_pbias"][ti] = safe_spearman(X[sel, 2], PB[s][ti][sel])
        out[s] = d
    return out


def boot_slice(data, masks, q, B, seed):
    """(B, J) median log10 cc of each slice's best runs, per series. Uses exactly the same re-drawn runs as
    analyze_truth_location_110.run_boot (same seed spec), so slice and all-runs results are paired."""
    out = {}
    for si, s in enumerate((SER_ON, SER_OFF)):
        X, K = data[s]["X"], data[s]["K"]
        n, J = X.shape[0], K.shape[0]
        idxb = np.random.default_rng([0 if seed is None else seed, 77, si]).integers(0, n, size=(B, n))
        lcc = X[:, 2][idxb]
        arr = np.full((B, J), np.nan)
        for ti in range(J):
            kb = K[ti][idxb]
            mb = masks[ti][idxb] & np.isfinite(kb)
            for b in range(B):
                sel = np.flatnonzero(mb[b])
                if len(sel) < MIN_SLICE_N:
                    continue
                m = max(MIN_SLICE_TOP, int(round(q * len(sel))))
                top = sel[np.argsort(-kb[b][sel], kind="stable")[:m]]
                arr[b, ti] = np.median(lcc[b][top])
        out[s] = arr
    return out


def cell_slope(ix, tx, med_pt, med_boot):
    """Tracking slope of log10 cc (best-run median against truth) over a cell's truths, with its bootstrap draws."""
    x = tx[ix, 2]
    if len(ix) < A.MIN_TRACK or np.ptp(x) <= 1e-12 or not np.isfinite(med_pt[ix]).all():
        return None, None
    sl = float(R.slope_of(x, med_pt[ix]))
    return sl, R.slope_of(x, med_boot[:, ix])


# ------------------------------------------------------------------
# Look 3: the ridge (correlations among Ks, f and cc within the best runs)
# ------------------------------------------------------------------
def ridge_null(X, sizes, Bn, seed):
    """Null limit for |rho| in a set of m runs chosen at random: the 95th percentile of |rho| over Bn random sets."""
    n = X.shape[0]
    lim = {}
    rng = np.random.default_rng([0 if seed is None else seed, 55])
    for m in sorted(set(int(v) for v in sizes)):
        sub = np.argsort(rng.random((Bn, n)), axis=1)[:, :m]
        for (a, b, ca, cb) in PAIRS:
            lim[(m, a, b)] = float(np.percentile(np.abs(rank_corr(X[sub, ca], X[sub, cb])), 95))
    return lim


def ridge_point(data, sets_on, truths, Bn, seed):
    X = data[SER_ON]["X"]
    sizes = [len(ix) for ix in sets_on]
    lim = ridge_null(X, sizes, Bn, seed)
    rows = []
    for ti in range(len(truths)):
        ix = sets_on[ti]
        row = {"n_best": len(ix)}
        for (a, b, ca, cb) in PAIRS:
            rho = float(rank_corr(X[ix, ca], X[ix, cb]))
            row["rho_%s_%s" % (a, b)] = rho
            row["null95_%s_%s" % (a, b)] = lim[(len(ix), a, b)]
        rows.append(row)
    return rows


# ------------------------------------------------------------------
# Assembling tables
# ------------------------------------------------------------------
def build_all(df, data, truths, tx, keys, cells, centre, frac, radius, q, B, Bn, seed):
    J = len(truths)
    PB = {s: metric_tensor(df, truths, data[s]["run_id"], s, "pbias_pct") for s in (SER_ON, SER_OFF)}
    facts, bad_facts = truth_facts(df, truths)
    masks = slice_masks(data[SER_ON]["X"], tx, radius)
    sp = slice_point(data, masks, q, PB)
    sb = boot_slice(data, masks, q, B, seed)
    sets = A.point_sets(data, [frac])
    pt = A.point_stats(data, sets, [frac])
    boot = A.run_boot(data, [frac], B, seed)
    ridge = ridge_point(data, sets[(SER_ON, frac)], truths, Bn, seed)

    per_truth = []
    for ti, t in enumerate(truths):
        row = {"truth_id": t["id"], "truth_Ks": t["Ks"], "truth_f": t["f"], "truth_cc": t["cc"],
               "truth_volume_m3": facts[t["id"]]["vol"], "truth_peak_m3s": facts[t["id"]]["peak"],
               "n_slice": int(sp[SER_ON]["n"][ti])}
        for s, tag in ((SER_ON, "on"), (SER_OFF, "ctl")):
            row["%s_slice_best_median_cc" % tag] = 10 ** sp[s]["med"][ti] if nz(sp[s]["med"][ti]) else np.nan
            row["%s_slice_best_run_cc" % tag] = 10 ** sp[s]["best_cc"][ti] if nz(sp[s]["best_cc"][ti]) else np.nan
            row["%s_slice_rho_cc_kge" % tag] = sp[s]["rho_kge"][ti]
            row["%s_slice_rho_cc_pbias" % tag] = sp[s]["rho_pbias"][ti]
        row["on_all_median_cc"] = 10 ** pt[(SER_ON, frac)]["med"][ti, 2]
        row.update(ridge[ti])
        per_truth.append(row)

    cell_rows, slope_draws = [], {}
    for k in keys:
        ix = cells[k]
        lev = leverage(ix, truths, facts)
        row = {"cell_Ks": k[0], "cell_f": k[1], "is_centre": bool(k == centre), "n_truths": len(ix),
               "n_slice": int(np.mean([sp[SER_ON]["n"][i] for i in ix])),
               "vol_change_pct": lev["vol_change_pct"], "peak_change_pct": lev["peak_change_pct"],
               "true_cc_low": lev["ccs"][0], "true_cc_high": lev["ccs"][-1]}
        for j, i in enumerate(ix):
            row["volume_m3_at_cc_%d" % j] = lev["vols"][j]
            row["peak_m3s_at_cc_%d" % j] = lev["peaks"][j]
        sl_all, b_all = cell_slope(ix, tx, pt[(SER_ON, frac)]["med"][:, 2], boot[(SER_ON, frac)][:, :, 2])
        sl_on, b_on = cell_slope(ix, tx, sp[SER_ON]["med"], sb[SER_ON])
        sl_ct, b_ct = cell_slope(ix, tx, sp[SER_OFF]["med"], sb[SER_OFF])
        for tag, sl, b in (("on_all", sl_all, b_all), ("on_slice", sl_on, b_on), ("ctl_slice", sl_ct, b_ct)):
            row["%s_slope" % tag] = sl if sl is not None else np.nan
            lo, hi = pct_ci(b) if b is not None else (np.nan, np.nan)
            row["%s_slope_ci_lo" % tag], row["%s_slope_ci_hi" % tag] = lo, hi
        if sl_all is not None and sl_on is not None:
            d = b_on - b_all
            row["slice_minus_all"] = sl_on - sl_all
            row["slice_minus_all_ci_lo"], row["slice_minus_all_ci_hi"] = pct_ci(d)
        else:
            row["slice_minus_all"] = row["slice_minus_all_ci_lo"] = row["slice_minus_all_ci_hi"] = np.nan
        row["weak"] = bool(nz(row["on_all_slope_ci_lo"]) and row["on_all_slope_ci_lo"] <= 0)
        row["strong"] = bool(nz(row["on_all_slope_ci_lo"]) and row["on_all_slope_ci_lo"] > 0)
        for (a, b, ca, cb) in PAIRS:
            key = "rho_%s_%s" % (a, b)
            vals = np.array([ridge[i][key] for i in ix])
            lim = np.array([ridge[i]["null95_%s_%s" % (a, b)] for i in ix])
            m = float(np.mean(vals))
            row[key + "_mean"] = m
            row[key + "_n_flagged"] = int(np.sum((np.abs(vals) > lim) & (np.sign(vals) == np.sign(m))))
        for s, tag in ((SER_ON, "on"), (SER_OFF, "ctl")):
            row["%s_slice_rho_cc_kge" % tag] = float(np.nanmean(sp[s]["rho_kge"][ix])) if np.isfinite(sp[s]["rho_kge"][ix]).any() else np.nan
            row["%s_slice_rho_cc_pbias" % tag] = float(np.nanmean(sp[s]["rho_pbias"][ix])) if np.isfinite(sp[s]["rho_pbias"][ix]).any() else np.nan
        cell_rows.append(row)
        slope_draws[k] = (b_all, b_on)
    return {"per_truth": pd.DataFrame(per_truth), "cells": pd.DataFrame(cell_rows), "facts": facts, "bad_facts": bad_facts,
            "sets": sets, "pt": pt, "masks": masks, "PB": PB, "sp": sp, "ridge": ridge}


def best_runs_table(df, data, truths, sets_on, PB_on, n_list):
    X, K = data[SER_ON]["X"], data[SER_ON]["K"]
    order = data[SER_ON]["run_id"]
    rows = []
    for ti, t in enumerate(truths):
        ok = np.flatnonzero(np.isfinite(K[ti]))
        top = ok[np.argsort(-K[ti][ok], kind="stable")][:n_list]
        for r, i in enumerate(top):
            rows.append({"truth_id": t["id"], "truth_Ks": t["Ks"], "truth_f": t["f"], "truth_cc": t["cc"], "rank": r + 1,
                         "run_id": order[i], "Ks_mult": X[i, 0], "f_RS_abs": 10 ** X[i, 1], "cc_mmhr": 10 ** X[i, 2],
                         "kge_2012": K[ti][i], "pbias_pct": PB_on[ti][i]})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------
def kfmt(v):
    return "%7.0f" % (v / 1000.0) if nz(v) else "    n/a"


def print_leverage(cells_df, keys, cells, truths, facts, tx):
    print("\n" + "=" * 150)
    print("1. LEVERAGE: how much does the TRUTH itself change when the true cc goes from its lowest to its highest value?")
    print("=" * 150)
    print("  If a cell's truth hardly changes, cc cannot be identified there by ANY method. Volumes are the truth hydrograph's total "
          "runoff (thousand m3); peaks in m3/s.")
    cc_vals = sorted({t["cc"] for t in truths})
    head = " / ".join(("%g" % v) for v in cc_vals)
    print("  %-16s | %-30s | %-9s | %-26s | %-9s" % ("cell", "volume at true cc " + head, "change", "peak at true cc " + head, "change"))
    for _, r in cells_df.iterrows():
        ix = cells[(r.cell_Ks, r.cell_f)]
        vols = " / ".join(kfmt(facts[truths[i]["id"]]["vol"]).strip().rjust(7) for i in ix)
        pks = " / ".join(("%7.1f" % facts[truths[i]["id"]]["peak"]) for i in ix)
        name = A.cell_name((r.cell_Ks, r.cell_f)) + (" *" if r.is_centre else "")
        print("  %-16s | %-30s | %+8.1f%% | %-26s | %+8.1f%%"
              % (name, vols, r.vol_change_pct, pks, r.peak_change_pct))
    print("  * = original centre cell. 'change' = the truth at the highest true cc compared with the truth at the lowest (a negative number "
          "means more channel loss gives less runoff).")


def print_slice(cells_df, frac, radius, q):
    print("\n" + "=" * 150)
    print("2. SLICE: with Ks and f held near the truth, does the best cc follow the true cc?   (all runs vs the slice; best %s of runs)"
          % A.cut_label(frac))
    print("=" * 150)
    print("  Slice = runs within %.2f (box units) of the truth's Ks and f; its best %d%% of runs (at least %d) are used. 1 = the median best cc "
          "follows the true cc one-for-one; 0 = it ignores it." % (radius, round(100 * q), MIN_SLICE_TOP))
    print("  'weak' = the all-runs slope's 95% interval includes 0; 'strong' = it lies entirely above 0 (labels for this report only).")
    print("  %-16s %-6s %7s | %-22s | %-22s | %-24s | %-22s | %7s %7s"
          % ("cell", "label", "n slice", "all runs [95%]", "slice, cc ON [95%]", "slice minus all [95%]", "slice, control [95%]",
             "rho cc-", "rho cc-"))
    print("  %-16s %-6s %7s | %-22s | %-22s | %-24s | %-22s | %7s %7s" % ("", "", "(runs)", "", "", "", "", "KGE", "PBIAS"))
    for _, r in cells_df.iterrows():
        name = A.cell_name((r.cell_Ks, r.cell_f)) + (" *" if r.is_centre else "")
        label = "weak" if r.weak else ("strong" if r.strong else "n/a")
        d = "n/a"
        if nz(r.slice_minus_all):
            d = interval_str(r.slice_minus_all, r.slice_minus_all_ci_lo, r.slice_minus_all_ci_hi)
            if nz(r.slice_minus_all_ci_lo) and (r.slice_minus_all_ci_lo > 0 or r.slice_minus_all_ci_hi < 0):
                d += " *"
        print("  %-16s %-6s %7d | %-22s | %-22s | %-24s | %-22s | %7s %7s"
              % (name, label, r.n_slice, interval_str(r.on_all_slope, r.on_all_slope_ci_lo, r.on_all_slope_ci_hi),
                 interval_str(r.on_slice_slope, r.on_slice_slope_ci_lo, r.on_slice_slope_ci_hi), d,
                 interval_str(r.ctl_slice_slope, r.ctl_slice_slope_ci_lo, r.ctl_slice_slope_ci_hi),
                 sgn(r.on_slice_rho_cc_kge, "%+5.2f", "  n/a"), sgn(r.on_slice_rho_cc_pbias, "%+5.2f", "  n/a")))
    print("  * in the 'slice minus all' column = that interval excludes 0 (holding Ks and f near the truth really changed the answer).")
    print("  rho cc-KGE / rho cc-PBIAS: Spearman correlation across the slice's runs between cc and KGE_2012 / PBIAS, averaged over the cell's "
          "truths. A cc that matters pushes PBIAS (volume error) down as it rises; a cc near the truth's pushes KGE up.")
    print("  Slice slopes use only about 8-12 best runs per truth and 3 truths per cell. Read single cells as hints.")


def print_ridge(cells_df, per_truth, frac):
    print("\n" + "=" * 150)
    print("3. RIDGE: are the best %s of runs strung along a trade-off between Ks, f and cc?" % A.cut_label(frac))
    print("=" * 150)
    n_best = int(per_truth.n_best.iloc[0])
    lims = {p: float(per_truth["null95_%s_%s" % p].mean()) for p in (("Ks", "cc"), ("f", "cc"), ("Ks", "f"))}
    print("  Spearman correlation among the %d best runs' parameters, averaged over the cell's 3 truths (in brackets: how many of the 3 exceed "
          "the chance limit with the same sign)." % n_best)
    print("  Chance limit: |rho| larger than about %.2f (Ks-cc), %.2f (f-cc), %.2f (Ks-f) happens in only 5%% of randomly chosen sets of the "
          "same size." % (lims[("Ks", "cc")], lims[("f", "cc")], lims[("Ks", "f")]))
    print("  A negative Ks-cc value means higher Ks pairs with lower cc: both remove water, so one can offset the other.")
    print("  %-16s %-6s | %-14s | %-14s | %-14s" % ("cell", "label", "Ks vs cc", "f vs cc", "Ks vs f"))
    for _, r in cells_df.iterrows():
        name = A.cell_name((r.cell_Ks, r.cell_f)) + (" *" if r.is_centre else "")
        label = "weak" if r.weak else ("strong" if r.strong else "n/a")
        cols = []
        for key in ("Ks_cc", "f_cc", "Ks_f"):
            m, nf = r["rho_%s_mean" % key], int(r["rho_%s_n_flagged" % key])
            cols.append("%+5.2f  (%d/3)%s" % (m, nf, " *" if nf >= 2 else "  "))
        print("  %-16s %-6s | %-14s | %-14s | %-14s" % (name, label, cols[0], cols[1], cols[2]))
    print("  * = at least 2 of the cell's 3 truths exceed the chance limit. The 3 truths in a cell share most of their best runs, so these "
          "counts are not independent, and with 9 cells x 3 pairs a few stray flags are expected by chance. Look for a consistent pattern "
          "(especially a negative Ks-cc), not a single flag.")


def print_side_by_side(df, data, truths, tx, cells, show, n_show, PB_on):
    X, K = data[SER_ON]["X"], data[SER_ON]["K"]
    order = data[SER_ON]["run_id"]
    print("\n" + "=" * 150)
    print("4. SIDE BY SIDE: the %d best cc-ON runs for each true cc, in the cells chosen below" % n_show)
    print("=" * 150)
    print("  Compare each run's Ks, f and cc with the truth's. If the best runs differ from the truth in Ks and f but compensate with cc "
          "(or the other way round), that is a trade-off.")
    for label, k in show:
        print("\n  Cell: %s   (%s)" % (A.cell_name(k), label))
        for i in cells[k]:
            t = truths[i]
            ok = np.flatnonzero(np.isfinite(K[i]))
            top = ok[np.argsort(-K[i][ok], kind="stable")][:n_show]
            print("    truth: Ks %g, f %g, cc %g mm/hr" % (t["Ks"], t["f"], t["cc"]))
            print("      %4s %7s %9s %9s | %9s %9s" % ("rank", "Ks", "f", "cc", "KGE_2012", "PBIAS %"))
            for r, j in enumerate(top):
                print("      %4d %7.2f %9.4f %9.1f | %9.3f %9.1f" % (r + 1, X[j, 0], 10 ** X[j, 1], 10 ** X[j, 2], K[i][j], PB_on[i][j]))
            print("      median of these: Ks %.2f, f %.4f, cc %.1f" % (np.median(X[top, 0]), 10 ** np.median(X[top, 1]), 10 ** np.median(X[top, 2])))


def summary_lines(cells_df, frac):
    out = []
    c = cells_df
    weak, strong = c[c.weak], c[c.strong]
    out.append("Cells by cc tracking (all-runs slope, best %s): %d weak (interval includes 0), %d strong (interval above 0), %d undefined."
               % (A.cut_label(frac), len(weak), len(strong), len(c) - len(weak) - len(strong)))
    if len(weak) + len(strong) == 0:
        out.append("No cell has %d or more cc values, so per-cell slopes and the weak/strong labels are undefined; the leverage and ridge "
                   "tables above still apply." % A.MIN_TRACK)
    if len(weak):
        out.append("Weak cells: " + "; ".join(A.cell_name((r.cell_Ks, r.cell_f)) for _, r in weak.iterrows()) + ".")
    for nm, g in (("weak", weak), ("strong", strong)):
        if not len(g):
            continue
        out.append("In the %s cells the truth's volume changes by %+.1f%% on median (range %+.1f%% to %+.1f%%) from the lowest to the highest "
                   "true cc, and its peak by %+.1f%% (range %+.1f%% to %+.1f%%)."
                   % (nm, g.vol_change_pct.median(), g.vol_change_pct.min(), g.vol_change_pct.max(), g.peak_change_pct.median(),
                      g.peak_change_pct.min(), g.peak_change_pct.max()))
    for nm, g in (("weak", weak), ("strong", strong)):
        g = g[g.on_slice_slope.notna()] if len(g) else g
        if not len(g):
            continue
        sig = int(((g.slice_minus_all_ci_lo > 0)).sum())
        out.append("In the %s cells the slice slope averages %+.2f (all-runs slope %+.2f; control slice %+.2f); the slice slope is detectably "
                   "above the all-runs slope in %d of %d."
                   % (nm, g.on_slice_slope.mean(), g.on_all_slope.mean(), g.ctl_slice_slope.mean(), sig, len(g)))
    for nm, g in (("weak", weak), ("strong", strong)):
        if not len(g):
            continue
        out.append("In the %s cells the mean rank correlation among the best runs is Ks-cc %+.2f, f-cc %+.2f, Ks-f %+.2f."
                   % (nm, g.rho_Ks_cc_mean.mean(), g.rho_f_cc_mean.mean(), g.rho_Ks_f_mean.mean()))
    return out


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def fig_slice_leverage(cells_df, centre, frac, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    ny = len(cells_df)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.5, 0.5 * ny + 2.8), facecolor=T.SURFACE, sharey=True,
                                 gridspec_kw={"width_ratios": [3, 2]})
    labels = [A.cell_name((r.cell_Ks, r.cell_f)) + ("  (original centre)" if r.is_centre else "") for _, r in cells_df.iterrows()]
    for a in (ax, bx):
        T._style(a)
        for yi, (_, r) in enumerate(cells_df.iterrows()):
            if r.is_centre:
                a.axhspan(yi - 0.5, yi + 0.5, color=T.COL_CC, alpha=0.08, linewidth=0, zorder=0)
    ax.axvline(0.0, color=T.INK, linestyle="--", linewidth=0.9, zorder=1)
    ax.axvline(1.0, color=T.COL_NULL, linestyle=":", linewidth=1.1, zorder=1)
    series = (("on_all", -0.2, "o", False, T.COL_CC, 1.6), ("on_slice", 0.0, "o", True, T.COL_CC, 2.0),
              ("ctl_slice", 0.2, "s", False, T.COL_NULL, 1.4))
    for yi, (_, r) in enumerate(cells_df.iterrows()):
        for tag, off, mk, filled, col, lw in series:
            sl, lo, hi = r["%s_slope" % tag], r["%s_slope_ci_lo" % tag], r["%s_slope_ci_hi" % tag]
            if not nz(sl):
                continue
            y = yi + off
            if nz(lo) and nz(hi):
                ax.hlines(y, lo, hi, color=col, linewidth=lw, zorder=3)
            ax.plot([sl], [y], linestyle="none", marker=mk, markersize=5.5, markerfacecolor=col if filled else T.SURFACE,
                    markeredgecolor=col, markeredgewidth=1.5, zorder=4)
    ax.set_yticks(range(ny))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("cc tracking slope (95% interval)", color=T.INK2, fontsize=9)
    ax.set_title("Does the best cc follow the true cc?", color=T.INK, fontsize=9.5, loc="left")
    for yi, (_, r) in enumerate(cells_df.iterrows()):
        v = r.vol_change_pct
        if not nz(v):
            continue
        bx.hlines(yi, 0, v, color=T.COL_CC, linewidth=2.0, zorder=3)
        bx.plot([v], [yi], linestyle="none", marker="o", markersize=6, markerfacecolor=T.COL_CC, markeredgecolor=T.COL_CC, zorder=4)
    bx.axvline(0.0, color=T.INK, linestyle="--", linewidth=0.9, zorder=1)
    bx.set_xlabel("change in the truth's volume, lowest to highest true cc (%)", color=T.INK2, fontsize=9)
    bx.set_title("How much does cc change the truth itself?", color=T.INK, fontsize=9.5, loc="left")
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5.5, markerfacecolor=T.SURFACE, markeredgewidth=1.5, linewidth=1.6,
               label="all runs (best %s)" % A.cut_label(frac)),
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5.5, linewidth=2.0, label="only runs with Ks and f near the truth"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=5.5, markerfacecolor=T.SURFACE, markeredgewidth=1.5, linewidth=1.4,
               label="control (cc inert), same slice"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=0.9, label="0 = ignores the truth"),
        Line2D([0], [0], color=T.COL_NULL, linestyle=":", linewidth=1.1, label="1 = follows it"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Why is cc recovered in some cells and not others?", color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.1, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


def fig_ridge_cells(data, truths, cells, show, frac, sets_on, ridge, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    X = data[SER_ON]["X"]
    nr = len(show)
    nc = max(len(cells[k]) for _, k in show)
    fig, axes = plt.subplots(nr, nc, figsize=(3.5 * nc + 0.8, 3.2 * nr + 1.4), facecolor=T.SURFACE, sharex=True, sharey=True, squeeze=False)
    kb, cb = P.BOX["Ks_mult"], P.BOX["channelconductivity_mmhr"]
    for ri, (label, k) in enumerate(show):
        for ci_ in range(nc):
            ax = axes[ri][ci_]
            if ci_ >= len(cells[k]):
                ax.set_visible(False)
                continue
            T._style(ax)
            ti = cells[k][ci_]
            t = truths[ti]
            best = sets_on[ti]
            ax.scatter(X[:, 0], 10 ** X[:, 2], s=10, color=T.COL_NULL, alpha=0.45, linewidths=0, zorder=2)
            ax.scatter(X[best, 0], 10 ** X[best, 2], s=26, color=T.COL_CC, edgecolors=T.SURFACE, linewidths=0.5, zorder=3)
            ax.plot([t["Ks"]], [t["cc"]], linestyle="none", marker="*", markersize=15, markerfacecolor=T.INK,
                    markeredgecolor=T.SURFACE, markeredgewidth=0.8, zorder=5)
            ax.set_yscale("log")
            ax.set_xlim(kb[0] - 0.2, kb[1] + 0.2)
            ax.set_ylim(cb[0] * 0.8, cb[1] * 1.25)
            ax.set_yticks([30, 100, 300, 1000])
            ax.set_yticklabels(["30", "100", "300", "1000"], fontsize=8)
            ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            rho = ridge[ti]["rho_Ks_cc"]
            ax.set_title("%s, true cc %g" % (A.cell_name(k), t["cc"]), color=T.INK, fontsize=9, loc="left")
            ax.text(0.97, 0.04, "best runs: rho(Ks, cc) = %+.2f" % rho, transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
                    color=T.INK2)
            if ci_ == 0:
                ax.set_ylabel(("%s\n" % label if label != "chosen" else "") + "cc (mm/hr)", color=T.INK2, fontsize=8.5)
            if ri == nr - 1:
                ax.set_xlabel("Ks multiplier", color=T.INK2, fontsize=8.5)
    handles = [
        Line2D([0], [0], color=T.COL_NULL, marker="o", markersize=4, linestyle="none", alpha=0.7, label="all 250 cc-ON runs"),
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5.5, linestyle="none", label="best %s of runs for this truth" % A.cut_label(frac)),
        Line2D([0], [0], color=T.INK, marker="*", markersize=11, linestyle="none", label="the truth (its Ks and cc)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Do the best runs lie along a trade-off between Ks and cc?  (a 2-D view of a 3-parameter set; f is not shown)",
                 color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
def parse_cells(items, keys):
    out = []
    for s in items:
        try:
            a, b = s.split(",")
            ks, f = float(a), float(b)
        except Exception:
            sys.exit("--show_cells needs entries like 5,0.02 (true Ks, true f). Could not read: %r" % s)
        hit = [k for k in keys if np.isclose(k[0], ks) and np.isclose(k[1], f)]
        if not hit:
            sys.exit("--show_cells: no cell with Ks %g, f %g. Cells present: %s" % (ks, f, "; ".join(A.cell_name(k) for k in keys)))
        out.append(("chosen", hit[0]))
    return out


def main():
    ap = argparse.ArgumentParser(description="Series 110 -- why is cc recovered in some (Ks, f) cells and not others?")
    ap.add_argument("--frac", type=float, default=0.20, help="Best-run fraction (by KGE_2012) for the all-runs slope and the ridge look (default 0.20)")
    ap.add_argument("--slice_radius", type=float, default=0.25,
                    help="Slice = runs within this distance (box units, Ks and log f) of the truth's Ks and f (default 0.25)")
    ap.add_argument("--slice_top", type=float, default=0.25, help="Share of a slice's runs kept as 'best' (default 0.25)")
    ap.add_argument("--boot", type=int, default=1000, help="Paired bootstrap re-draws (default 1000)")
    ap.add_argument("--null_sets", type=int, default=3000, help="Random sets for the ridge chance limit (default 3000)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--expect_n", type=int, default=250)
    ap.add_argument("--list_runs", type=int, default=5, help="Best runs listed per truth in look 4 (default 5)")
    ap.add_argument("--show_cells", nargs="*", default=None,
                    help="Cells for look 4 and the ridge figure, like 5,0.02 7,0.012 (default: the weakest- and strongest-tracked cell)")
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--continue_on_fail", action="store_true")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()
    if not 0 < args.slice_radius <= 1.5:
        sys.exit("--slice_radius must be between 0 and 1.5 (box units).")
    if not 0 < args.slice_top <= 1 or not 0 < args.frac <= 1:
        sys.exit("--frac and --slice_top must be between 0 and 1.")

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    loc_dir = summary_dir.parent.parent / A.LOC_DIRNAME
    out_dir = summary_dir / A.OUT_DIRNAME
    print("\n" + "=" * 78 + "\nSeries 110 -- closer look: why is cc recovered in some cells and not others?\n" + "=" * 78)
    print("Reading from: %s" % summary_dir)

    df = A.load_long(summary_dir)
    miss = [c for c in EXTRA_COLS if c not in df.columns]
    if miss:
        sys.exit("The re-scored table is missing columns %s. Re-run rescore_truth_location_110.py." % miss)
    for c in EXTRA_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    truths = A.truth_list(df)
    ck = A.run_checks(df, truths, summary_dir, loc_dir, args.expect_n)
    facts, bad = truth_facts(df, truths)
    if bad:
        ck.add("FAIL", "truth hydrograph facts", "the truth's own peak/volume is missing or differs between runs for: %s" % ", ".join(bad[:4]))
    else:
        ck.add("PASS", "truth hydrograph facts", "each truth has one consistent peak and volume across all its rows")
    P.print_checks(ck)
    if ck.n("FAIL") and not args.continue_on_fail:
        sys.exit("\nSTOPPED: an input check failed. Nothing was analysed or written. (--continue_on_fail to override.)")

    data = A.build_tensors(df, truths)
    tx = A.truth_coords(truths)
    keys, cells, centre = A.find_cells(truths)
    n = data[SER_ON]["X"].shape[0]
    show_user = parse_cells(args.show_cells, keys) if args.show_cells else None      # fail early on a bad cell name
    print("\n%d truth(s) in %d cell(s), %d runs per series; best %s for the all-runs slope; slice radius %.2f (about %d runs per slice); "
          "%d bootstrap re-draws" % (len(truths), len(keys), n, A.cut_label(args.frac), args.slice_radius,
                                     int(np.round(np.mean(slice_masks(data[SER_ON]["X"], tx, args.slice_radius).sum(axis=1)))), args.boot))
    print("\nRunning the bootstrap and the chance limits ...", flush=True)
    res = build_all(df, data, truths, tx, keys, cells, centre, args.frac, args.slice_radius, args.slice_top, args.boot, args.null_sets, args.seed)
    cells_df, per_truth = res["cells"], res["per_truth"]
    PB_on = res["PB"][SER_ON]

    # cells for look 4 and the ridge figure
    if show_user:
        show = show_user
    else:
        have = cells_df[cells_df.on_all_slope.notna()]
        if len(have) >= 2:
            lo, hi = have.loc[have.on_all_slope.idxmin()], have.loc[have.on_all_slope.idxmax()]
            show = [("weakest all-runs cc slope", (lo.cell_Ks, lo.cell_f)), ("strongest all-runs cc slope", (hi.cell_Ks, hi.cell_f))]
        else:
            show = [("first cell", keys[0])]

    print_leverage(cells_df, keys, cells, truths, facts, tx)
    print_slice(cells_df, args.frac, args.slice_radius, args.slice_top)
    print_ridge(cells_df, per_truth, args.frac)
    print_side_by_side(df, data, truths, tx, cells, show, args.list_runs, PB_on)
    print("\n" + "=" * 150)
    print("5. WHAT THE NUMBERS SAY   (descriptive; assembled from the tables above)")
    print("=" * 150)
    for line in summary_lines(cells_df, args.frac):
        print("  - " + line)
    print("  Caveat: one storm, %d runs, routing pinned at truth, noise-free truths, 3 truths per cell and a rough slice. The pattern across "
          "weak and strong cells is the evidence; a single cell is a hint." % n)

    out_dir.mkdir(parents=True, exist_ok=True)
    cells_df.to_csv(out_dir / "diag_cells_110.csv", index=False)
    per_truth.to_csv(out_dir / "diag_per_truth_110.csv", index=False)
    best_runs_table(df, data, truths, res["sets"][(SER_ON, args.frac)], PB_on, 10).to_csv(out_dir / "diag_best_runs_110.csv", index=False)
    if not args.no_plots:
        for name, fn in (("fig_diag_slice_leverage_110.png", lambda p: fig_slice_leverage(cells_df, centre, args.frac, p)),
                         ("fig_diag_ridge_cells_110.png",
                          lambda p: fig_ridge_cells(data, truths, cells, show, args.frac, res["sets"][(SER_ON, args.frac)], res["ridge"], p))):
            try:
                fn(out_dir / name)
            except Exception as e:
                print("\n  (figure %s skipped: %s -- the CSVs are saved regardless)" % (name, e))

    def md5_of(p):
        return hashlib.md5(Path(p).read_bytes()).hexdigest()
    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "pandas": pd.__version__, "numpy": np.__version__,
            "rescored_long": md5_of(summary_dir / A.RS_DIRNAME / A.LONG_NAME),
            "settings": {k: (v if not isinstance(v, Path) else str(v)) for k, v in vars(args).items()},
            "cells_shown": [[lab, list(k)] for lab, k in show], "checks": ck.rows}
    (out_dir / "PROVENANCE_diag_110.json").write_text(json.dumps(prov, indent=2, default=str))
    print("\nSaved to: %s" % out_dir)
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results." % ck.n("FAIL"))


if __name__ == "__main__":
    main()
