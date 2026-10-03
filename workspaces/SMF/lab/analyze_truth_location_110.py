"""
analyze_truth_location_110.py
=============================
Series 110 -- the TRUTH-LOCATION test, part 3: does cc (and Ks, f) still get
recovered when the TRUE Ks and f are not the centre of the sampled box?

THE QUESTION
------------
Everything so far used one truth location: Ks 7.0, f 0.012, in the middle of the
box the 250 runs were drawn from. A skeptic's worry is that the answer was easy
because it sat in the middle. run_truth_location_110.py made new truths with
different TRUE Ks and f (and three true cc values at each), and
rescore_truth_location_110.py scored the 250 + 250 saved runs against all of
them. This script reads those scores and asks, truth by truth:

  1. Where do the best runs (by KGE_2012) put cc, Ks and f, compared with the
     values that made the truth?
  2. As the true cc changes (95 -> 201 -> 425), do the best runs' cc follow it?
     ("tracking slope", computed separately for each (Ks, f) cell and pooled.)
  3. Does a cell with a different true Ks/f track cc as well as the ORIGINAL
     centre cell does? (the same three cc values, so it is a like-for-like
     comparison; paired bootstrap)
  4. Do the best runs recover Ks and f as the truth moves?
  5. Is any weakness just because the sweep had few samples near that truth?
     ("nearest sample" distance, in box units)

Pure post-hoc analysis. Runs no tRIBS. Writes only into its own new folder
(calibration_work/03_comparisons/summary_tables/location_cc_110/).

Keep this file in the same folder as analyze_cc_recovery_110.py,
analyze_cc_tests_110.py and analyze_cc_pca_110.py (it reuses their helpers).
Run it from the lab/ directory.

VOCABULARY
-----------
  Truth location  The true (Ks, f, cc) that made a synthetic truth hydrograph.
  Cell            One (true Ks, true f) pair. Each cell holds the truths at its
                  different true cc values (default 3).
  Best runs       The top 20% / 10% / 5% of the 250 runs by KGE_2012 against that
                  truth (at least 10 runs).
  Median          The middle value of a parameter among the best runs.
  Tracking slope  Plot the best runs' median parameter against the true parameter,
                  on the scale the sweep was drawn on (cc and f logarithmic, Ks
                  linear), and fit a line. Slope 1 = the median follows the truth
                  one-for-one; slope 0 = it ignores the truth. The edges of the
                  sampled box squash the slope below 1 even for a perfect method,
                  so compare it with 0 and with the control, not only with 1.
  Control         The same 250 parameter points run with channel loss OFF. In that
                  series cc has no effect, so its "cc tracking" shows what pure
                  chance gives. (Its Ks and f tracking show what Ks and f alone can
                  do when the model cannot represent channel loss.)
  Paired bootstrap  Re-draw the 250 runs with replacement and redo everything; the
                  SAME re-draw is used for every truth, so differences between
                  cells are measured on the same wobble. The 95% interval of a
                  difference that excludes 0 means the cells really differ.
  Permutation test  Shuffle a parameter among the runs to show what "no
                  information" looks like.
  Nearest sample  Distance from the truth point to the closest of the 250 sampled
                  points, in units of the box (0 = same point; 1 = a full side).
                  A uniform sweep puts about 3-4 samples within 0.15 of a point in
                  the middle. A large distance means the sweep barely looked near
                  that truth, so a poor recovery there says little about cc.

HOW TO READ THE RESULTS
------------------------
  * ON slope clearly above the control's, in every cell, and about as high as the
    centre cell's: cc is recovered wherever Ks and f are, not just at the centre.
  * Centre cell good, off-centre cells low: the earlier success depended on the truth
    sitting in the middle of the box. That is the skeptic's case.
  * Low slopes in cells whose truth has few samples nearby (large "nearest"): treat as
    inconclusive rather than as failure.
  * The control's cc slope should stay near 0. If it does not, something other than
    cc is driving the pattern.

OUTPUTS   (calibration_work/03_comparisons/summary_tables/location_cc_110/)
-----------------------------------------------------------------------------
  location_per_truth_110.csv    series x cut x truth: median / middle 80% / bias for cc, Ks
                                and f, bootstrap intervals, p-values, nearest sample
  location_cells_110.csv        series x cut x cell: cc tracking slope, mean biases, and the
                                paired difference from the centre cell
  location_pooled_110.csv       series x cut: tracking slopes pooled over all truths
                                (cc, Ks, f), centre vs off-centre
  location_input_checks_110.csv
  fig_location_recovery_cc_110.png       one small panel per cell: best-run cc against true cc
  fig_location_cell_slopes_110.png       cc tracking slope per cell, with the control
  fig_location_Ks_f_recovery_110.png     best-run Ks and f against the true Ks and f
  PROVENANCE_location_cc_110.json

USAGE (run from the lab/ directory)
------------------------------------
    python analyze_truth_location_110.py                    # about a minute
    python analyze_truth_location_110.py --top_fracs 0.2 0.1 0.05 --boot 2000 --perm 10000
    python analyze_truth_location_110.py --within 3 --no_plots

CAVEATS
-------
One storm, 250 runs, routing pinned at truth, noise-free truths: a best case. Routing
at the truth's value means a truth with different routing is a different test (not done
here). Each cell has only a few cc values, so a per-cell slope is rough; trust the pooled
numbers and the pattern across cells more than any single cell. The 27 truths are not
independent (they share the same 250 runs), so p-values describe the runs, not 27
separate experiments.
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
    import analyze_cc_recovery_110 as R
except ImportError:
    sys.exit("This script needs analyze_cc_recovery_110.py, analyze_cc_tests_110.py and analyze_cc_pca_110.py in the "
             "same folder (it reuses their helpers). Put all of them in lab/ and run from there.")

T = R.T
P = T.P

_NEED = {
    "R": (R, ["top_index", "slope_of", "fcc", "ffactor", "CENTER_LOG", "LOG_LO", "LOG_HI"]),
    "T": (T, ["n_top_for", "make_perm_matrix", "p_low", "p_high", "fmt_p", "_style", "COL_CC", "COL_NULL", "INK",
              "INK2", "GRID", "SURFACE"]),
    "P": (P, ["Checks", "print_checks", "to_X", "generate_lhs_samples", "LHS_PARAMS", "DESIGN_SEED", "PARAMS", "BOX",
              "SER_ON", "SER_OFF", "md5_of", "STAGE2_NAME", "CONTROL_NAME", "NULL_RHO", "KEY_METRICS"]),
}
for _label, (_mod, _names) in _NEED.items():
    _miss = [n for n in _names if not hasattr(_mod, n)]
    if _miss:
        sys.exit("%s in this folder is an older version (missing: %s). Replace it with the current copy."
                 % (_mod.__name__ + ".py", ", ".join(_miss)))

SER_ON, SER_OFF = P.SER_ON, P.SER_OFF
RS_DIRNAME = "rescored_truth_location_110"
LONG_NAME = "rescored_location_long_110.csv"
RS_PROV = "PROVENANCE_rescored_truth_location_110.json"
LOC_DIRNAME = "synth_truth_location"
LOC_PROV = "PROVENANCE_truth_location_110.json"
OUT_DIRNAME = "location_cc_110"

STANDARD_KS, STANDARD_F = 7.0, 0.012
STAGE2_ROUTING = {"kinemvelcoef": 4.5, "flowexp": 0.24, "channelroughness": 0.026}
MIN_TRACK = 3                       # fewest truths for a tracking slope
NEAR_RADIUS = 0.15
PAR = ["cc", "Ks", "f"]             # reporting order
PCOL = {"Ks": 0, "f": 1, "cc": 2}   # column in the X matrix (Ks, log10 f, log10 cc)
PNAME = {"cc": "cc", "Ks": "Ks", "f": "f"}
KEY_COLS = ["truth_id", "truth_Ks", "truth_f", "truth_cc", "series", "run_id", "Ks_mult", "f_RS_abs",
            "channelconductivity_mmhr", "optpercolation", "kge_2012", "pbias_pct", "kge_r", "kge_beta", "kge_gamma"]


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def cut_label(frac):
    return "%d%%" % round(100 * frac)


def gfmt(v):
    return "%g" % v


def fnum(v, spec):
    return "   n/a" if (v is None or not np.isfinite(v)) else format(v, spec)


def box_pos(v, key):
    lo, hi = P.BOX[key]
    if key == "Ks_mult":
        return (v - lo) / (hi - lo)
    return (np.log10(v) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))


def iv(sl, lo, hi):
    """'+0.67 [+0.49, +0.94]' with a dagger when the point estimate falls outside its own bootstrap interval."""
    odd = not (lo - 1e-9 <= sl <= hi + 1e-9)
    return "%+5.2f [%+5.2f, %+5.2f]%s" % (sl, lo, hi, "!" if odd else " ")


def param_err_fmt(p, e):
    """A typical error for printing: factor for cc and f (log scale), absolute for Ks."""
    if not np.isfinite(e):
        return "n/a"
    return "%.2f" % e if p == "Ks" else "x%.2f" % (10 ** e)


# ------------------------------------------------------------------
# Loading
# ------------------------------------------------------------------
def load_long(summary_dir):
    path = summary_dir / RS_DIRNAME / LONG_NAME
    if not path.exists():
        sys.exit("Required file not found: %s\nRun rescore_truth_location_110.py first (after "
                 "run_truth_location_110.py has finished)." % path)
    df = pd.read_csv(path)
    miss = [c for c in KEY_COLS if c not in df.columns]
    if miss:
        sys.exit("%s is missing columns: %s. Re-run rescore_truth_location_110.py." % (path.name, miss))
    for c in KEY_COLS:
        if c not in ("truth_id", "series", "run_id"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def truth_list(df):
    g = (df.drop_duplicates("truth_id")[["truth_id", "truth_Ks", "truth_f", "truth_cc"]]
         .sort_values(["truth_Ks", "truth_f", "truth_cc"]).reset_index(drop=True))
    out = []
    for r in g.itertuples(index=False):
        out.append({"id": r.truth_id, "Ks": float(r.truth_Ks), "f": float(r.truth_f), "cc": float(r.truth_cc)})
    return out


def truth_coords(truths):
    """(J, 3) true parameters in X-space: Ks, log10 f, log10 cc."""
    return np.array([[t["Ks"], np.log10(t["f"]), np.log10(t["cc"])] for t in truths], dtype=float)


def build_tensors(df, truths):
    """Align runs across truths by run_id (the input checks guarantee identical run sets)."""
    out = {}
    for s in (SER_ON, SER_OFF):
        base = df[(df.truth_id == truths[0]["id"]) & (df.series == s)].reset_index(drop=True)
        order = base["run_id"].tolist()
        K = np.full((len(truths), len(order)), np.nan)
        for ti, t in enumerate(truths):
            g = df[(df.truth_id == t["id"]) & (df.series == s)].set_index("run_id")["kge_2012"]
            K[ti] = g.reindex(order).to_numpy(float)
        X = np.column_stack([base["Ks_mult"].to_numpy(float), np.log10(base["f_RS_abs"].to_numpy(float)),
                             np.log10(base["channelconductivity_mmhr"].to_numpy(float))])
        out[s] = {"run_id": order, "X": X, "K": K}
    return out


def find_cells(truths):
    """cells[(Ks, f)] = truth indices sorted by cc. Returns (ordered list of keys, dict, centre key or None)."""
    cells = {}
    for i, t in enumerate(truths):
        cells.setdefault((t["Ks"], t["f"]), []).append(i)
    for k in cells:
        cells[k].sort(key=lambda i: truths[i]["cc"])
    keys = sorted(cells)
    centre = None
    for k in keys:
        if np.isclose(k[0], STANDARD_KS) and np.isclose(k[1], STANDARD_F):
            centre = k
    return keys, cells, centre


def cell_name(k):
    return "Ks %g, f %g" % (k[0], k[1])


# ------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------
def run_checks(df, truths, summary_dir, loc_dir, expect_n):
    ck = P.Checks()
    ids = [t["id"] for t in truths]

    # 1. provenance of the re-scoring ---------------------------------------
    prov_path = summary_dir / RS_DIRNAME / RS_PROV
    if not prov_path.exists():
        ck.add("WARN", "re-scoring provenance", "%s missing; cannot confirm the re-scoring gates passed" % RS_PROV)
    else:
        prov = json.loads(prov_path.read_text())
        bad = []
        for name in (P.STAGE2_NAME, P.CONTROL_NAME):
            rec = prov.get("inputs", {}).get(name, {}).get("md5")
            p = summary_dir / name
            if rec != (P.md5_of(p) if p.exists() else None):
                bad.append(name)
        if bad:
            ck.add("FAIL", "inputs unchanged since re-scoring",
                   "checksum differs from the one recorded by the re-scoring run: %s. Re-run "
                   "rescore_truth_location_110.py." % ", ".join(bad))
        else:
            ck.add("PASS", "inputs unchanged since re-scoring",
                   "Stage 2 and control CSV checksums match the re-scoring provenance")
        gates = prov.get("gates", {})
        excl = gates.get("runs_excluded", {})
        if excl:
            ck.add("WARN", "re-scoring excluded runs", "%d run(s) were excluded at re-scoring (--drop_bad_runs)" % len(excl))
        else:
            ck.add("PASS", "re-scoring gates A and B", "all %s runs passed gates A and B at re-scoring"
                   % gates.get("runs_checked", "?"))
        rep = gates.get("D_reproduction", [])
        if gates.get("D_overridden"):
            ck.add("FAIL", "gate D (reproduces the earlier cc candidates)",
                   "gate D FAILED at re-scoring and was overridden (--ignore_reproduction): the new truths may not be "
                   "comparable with the 500 stored runs")
        elif not rep:
            ck.add("WARN", "gate D (reproduces the earlier cc candidates)",
                   "not checked at re-scoring (no centre truth matched an earlier candidate)")
        else:
            ck.add("PASS", "gate D (reproduces the earlier cc candidates)",
                   "%d centre truth(s) reproduce the earlier cc candidates" % len(rep))
        verified = prov.get("truths", {})
        extra = [t["id"] for t in truths if t["id"] not in verified]
        wrong = [t["id"] for t in truths if t["id"] in verified and not (
            np.isclose(verified[t["id"]]["Ks_mult"], t["Ks"]) and np.isclose(verified[t["id"]]["f_RS_abs"], t["f"])
            and np.isclose(verified[t["id"]]["cc_mmhr"], t["cc"]))]
        if extra or wrong:
            ck.add("FAIL", "truths verified at re-scoring",
                   "not in the verified list or with different parameters: %s" % ", ".join((extra + wrong)[:4]))
        else:
            ck.add("PASS", "truths verified at re-scoring", "%d truth(s), all checksum-verified at re-scoring" % len(truths))

    # 2. counts and run identity ----------------------------------------------
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
    for s in (SER_ON, SER_OFF):
        if s not in sets:
            problems.append("no rows for series %s" % s)
    if problems:
        ck.add("FAIL", "run counts and identity", "; ".join(problems[:4]))
    else:
        ck.add("PASS", "run counts and identity",
               "%d truth(s) x 2 series, %d unique runs each, same runs for every truth" % (len(truths), expect_n))

    # 3. series flags ----------------------------------------------------------
    on_bad = int((df.loc[df.series == SER_ON, "optpercolation"] != 1).sum())
    off_bad = int((df.loc[df.series == SER_OFF, "optpercolation"] != 0).sum())
    if on_bad or off_bad:
        ck.add("FAIL", "optpercolation flags", "%d cc-ON rows not 1, %d cc-OFF rows not 0" % (on_bad, off_bad))
    else:
        ck.add("PASS", "optpercolation flags", "cc-ON rows all 1, cc-OFF rows all 0")

    # 4. design: parameters equal the intended LHS --------------------------------
    design = P.generate_lhs_samples(expect_n, P.LHS_PARAMS, P.DESIGN_SEED)
    D = design[P.PARAMS].to_numpy(float)
    for s in (SER_ON, SER_OFF):
        g = df[(df.truth_id == ids[0]) & (df.series == s)]
        Rm = g[P.PARAMS].to_numpy(float)
        if len(Rm) == 0 or np.isnan(Rm).any():
            ck.add("FAIL", "design match (%s)" % s, "missing or non-numeric parameter values")
            continue
        close = np.all(np.isclose(Rm[:, None, :], D[None, :, :], rtol=1e-3, atol=0.0), axis=2)
        unmatched = int((close.sum(axis=1) == 0).sum())
        used = close.any(axis=0)
        if unmatched == 0 and used.sum() == min(len(Rm), len(D)):
            ck.add("PASS", "design match (%s)" % s,
                   "all %d (Ks, f, cc) points equal the seed-%d, n=%d design (rtol 1e-3)" % (len(Rm), P.DESIGN_SEED, expect_n))
        else:
            ck.add("FAIL", "design match (%s)" % s,
                   "%d point(s) do not match the seed-%d design; %d of %d design points are used"
                   % (unmatched, P.DESIGN_SEED, int(used.sum()), len(D)))

    # 5. truths inside the sampled box ---------------------------------------------
    outside = [t["id"] for t in truths
               if not (0.0 <= box_pos(t["Ks"], "Ks_mult") <= 1.0 and 0.0 <= box_pos(t["f"], "f_RS_abs") <= 1.0
                       and 0.0 <= box_pos(t["cc"], "channelconductivity_mmhr") <= 1.0)]
    if outside:
        ck.add("FAIL", "truths inside the sampled box", "outside the box: %s" % ", ".join(outside[:4]))
    else:
        lows = min(min(box_pos(t["Ks"], "Ks_mult"), box_pos(t["f"], "f_RS_abs"),
                       box_pos(t["cc"], "channelconductivity_mmhr")) for t in truths)
        highs = max(max(box_pos(t["Ks"], "Ks_mult"), box_pos(t["f"], "f_RS_abs"),
                        box_pos(t["cc"], "channelconductivity_mmhr")) for t in truths)
        ck.add("PASS", "truths inside the sampled box",
               "all %d truths inside the box (positions %.2f to %.2f of each range)" % (len(truths), lows, highs))

    # 6. routing of the truths equals the routing of every Stage 2 run -----------------
    lp = loc_dir / LOC_PROV
    if not lp.exists():
        ck.add("WARN", "truth routing", "%s not found; cannot confirm the truths use the Stage 2 routing" % LOC_PROV)
    else:
        try:
            rt = json.loads(lp.read_text()).get("routing_truth", {})
        except Exception:
            rt = {}
        diff = [k for k, v in STAGE2_ROUTING.items() if k in rt and not np.isclose(float(rt[k]), v)]
        if not rt or any(k not in rt for k in STAGE2_ROUTING):
            ck.add("WARN", "truth routing", "%s does not record the routing values, so I cannot confirm the truths use "
                   "the Stage 2 routing (cv %g, r %g, n %g)" % (LOC_PROV, STAGE2_ROUTING["kinemvelcoef"],
                                                                 STAGE2_ROUTING["flowexp"], STAGE2_ROUTING["channelroughness"]))
        elif diff:
            ck.add("FAIL", "truth routing", "routing recorded for the new truths differs from Stage 2's "
                   "(cv %g, r %g, n %g): %s" % (STAGE2_ROUTING["kinemvelcoef"], STAGE2_ROUTING["flowexp"],
                                               STAGE2_ROUTING["channelroughness"], ", ".join(diff)))
        else:
            ck.add("PASS", "truth routing", "truths use the Stage 2 routing (cv %g, r %g, n %g)"
                   % (STAGE2_ROUTING["kinemvelcoef"], STAGE2_ROUTING["flowexp"], STAGE2_ROUTING["channelroughness"]))

    # 7. finite metrics ----------------------------------------------------------------
    worst = (1.0, "")
    for (tid, s), gg in df.groupby(["truth_id", "series"]):
        for m in P.KEY_METRICS:
            fr = float(np.isfinite(gg[m]).mean())
            if fr < worst[0]:
                worst = (fr, "truth %s %s %s" % (tid, s, m))
    if worst[0] < 0.99:
        ck.add("FAIL", "finite metrics", "worst finite fraction %.3f (%s)" % worst)
    else:
        ck.add("PASS", "finite metrics", "all key metrics >= 99%% finite (worst %.3f)" % worst[0])

    # 8. null check: cc is inert in the control -------------------------------------------
    worst_rho, worst_t = 0.0, None
    for tid in ids:
        g = df[(df.truth_id == tid) & (df.series == SER_OFF)]
        ok = np.isfinite(g["kge_2012"]) & np.isfinite(g["channelconductivity_mmhr"])
        if ok.sum() > 10:
            rho = stats.spearmanr(np.log10(g.loc[ok, "channelconductivity_mmhr"]), g.loc[ok, "kge_2012"])[0]
            if np.isfinite(rho) and abs(rho) > abs(worst_rho):
                worst_rho, worst_t = rho, tid
    if worst_t is not None and abs(worst_rho) > P.NULL_RHO:
        ck.add("WARN", "null check (cc inert in control)",
               "rho(log cc, KGE_2012) = %+.2f at truth %s in the cc-OFF series (limit %.2f)" % (worst_rho, worst_t, P.NULL_RHO))
    else:
        ck.add("PASS", "null check (cc inert in control)",
               "largest |rho(log cc, KGE_2012)| in the cc-OFF series is %.2f (limit %.2f)" % (abs(worst_rho), P.NULL_RHO))

    # 9. shape of the grid ---------------------------------------------------------------------
    keys, cells, centre = find_cells(truths)
    small = [cell_name(k) for k in keys if len(cells[k]) < MIN_TRACK]
    if small:
        ck.add("WARN", "grid shape", "%d cell(s) have fewer than %d truths, so no per-cell slope: %s"
               % (len(small), MIN_TRACK, "; ".join(small[:3])))
    elif centre is None:
        ck.add("WARN", "grid shape", "%d cells, but none is the original centre (Ks %g, f %g): no comparison with the "
               "centre cell" % (len(keys), STANDARD_KS, STANDARD_F))
    else:
        ck.add("PASS", "grid shape", "%d cells of %d-%d truths each, including the original centre"
               % (len(keys), min(len(cells[k]) for k in keys), max(len(cells[k]) for k in keys)))
    return ck


# ------------------------------------------------------------------
# Fairness: how close is the nearest sample?
# ------------------------------------------------------------------
def nearest_sample(X, tx):
    """Distance (box units) from each truth to its nearest sample, and the number within NEAR_RADIUS.
    X is (n, 3) in X-space; tx is (J, 3)."""
    lo = np.array([P.BOX["Ks_mult"][0], np.log10(P.BOX["f_RS_abs"][0]), np.log10(P.BOX["channelconductivity_mmhr"][0])])
    hi = np.array([P.BOX["Ks_mult"][1], np.log10(P.BOX["f_RS_abs"][1]), np.log10(P.BOX["channelconductivity_mmhr"][1])])
    span = hi - lo
    d = np.sqrt((((X[None, :, :] - tx[:, None, :]) / span) ** 2).sum(axis=2))      # (J, n)
    return d.min(axis=1), (d <= NEAR_RADIUS).sum(axis=1)


# ------------------------------------------------------------------
# Resampling machinery
# ------------------------------------------------------------------
def point_sets(data, fracs):
    """sets[(series, frac)] = list over truths of run-index arrays (the best runs)."""
    sets = {}
    for s in (SER_ON, SER_OFF):
        K = data[s]["K"]
        for frac in fracs:
            sets[(s, frac)] = [R.top_index(K[ti], frac) for ti in range(K.shape[0])]
    return sets


def point_stats(data, sets, fracs):
    """pt[(s, frac)] = dict(med, q10, q90) each (J, 3)."""
    pt = {}
    for s in (SER_ON, SER_OFF):
        X = data[s]["X"]
        for frac in fracs:
            idxs = sets[(s, frac)]
            med = np.array([np.median(X[ix], axis=0) for ix in idxs])
            q = np.array([np.percentile(X[ix], [10, 90], axis=0) for ix in idxs])
            pt[(s, frac)] = {"med": med, "q10": q[:, 0, :], "q90": q[:, 1, :]}
    return pt


def run_boot(data, fracs, B, seed):
    """boot[(s, frac)] = (B, J, 3) medians of the best runs. The SAME resampled run indices are used for every truth
    and every cut within a series; ON and control are drawn independently."""
    boot = {}
    for si, s in enumerate((SER_ON, SER_OFF)):
        X, K = data[s]["X"], data[s]["K"]
        n = X.shape[0]
        idxb = np.random.default_rng([0 if seed is None else seed, 77, si]).integers(0, n, size=(B, n))
        Xb = X[idxb]                                       # (B, n, 3)
        for frac in fracs:
            arr = np.full((B, K.shape[0], 3), np.nan)
            for ti in range(K.shape[0]):
                kb = K[ti][idxb]
                kb = np.where(np.isfinite(kb), kb, -np.inf)
                m = T.n_top_for(int(np.isfinite(K[ti]).sum()), frac)
                order = np.argsort(-kb, axis=1, kind="stable")[:, :m]
                v = np.take_along_axis(Xb, order[:, :, None], axis=1)          # (B, m, 3)
                arr[:, ti, :] = np.median(v, axis=1)
            boot[(s, frac)] = arr
    return boot


def run_perm(data, sets, truths, fracs, Bp, seed, half_width):
    """null[(s, frac)] = dict(med (Bp, J, 3), within (Bp, J)): the best-runs' medians, and the share of them within
    half_width of the true cc, when run parameters are shuffled among runs and the real top sets are kept. ONE shuffle is
    shared by all truths."""
    null = {}
    tl = np.log10([t["cc"] for t in truths])
    for si, s in enumerate((SER_ON, SER_OFF)):
        X = data[s]["X"]
        Pm = T.make_perm_matrix(X.shape[0], Bp, np.random.default_rng([0 if seed is None else seed, 66, si]))
        for frac in fracs:
            med = np.full((Bp, len(truths), 3), np.nan)
            within = np.full((Bp, len(truths)), np.nan)
            for ti in range(len(truths)):
                v = X[Pm[:, sets[(s, frac)][ti]]]                   # (Bp, m, 3)
                med[:, ti, :] = np.median(v, axis=1)
                within[:, ti] = np.mean(np.abs(v[:, :, 2] - tl[ti]) <= half_width + 1e-12, axis=1)
            null[(s, frac)] = {"med": med, "within": within}
    return null


# ------------------------------------------------------------------
# Tracking
# ------------------------------------------------------------------
def track(ix, p, tx, pt_med, bmed, nmed):
    """Tracking statistics of parameter p over the truths ix. Returns (dict, bootstrap slopes) or (None, None)."""
    col = PCOL[p]
    x = tx[ix, col]
    if len(ix) < MIN_TRACK or np.ptp(x) <= 1e-12:
        return None, None
    obs = pt_med[ix, col]
    sl, err = float(R.slope_of(x, obs)), float(np.mean(np.abs(obs - x)))
    nl, bt = nmed[:, ix, col], bmed[:, ix, col]
    sl_null, err_null = R.slope_of(x, nl), np.mean(np.abs(nl - x), axis=1)
    sl_b, err_b = R.slope_of(x, bt), np.mean(np.abs(bt - x), axis=1)
    out = {
        "n_truths": int(len(ix)), "slope": sl,
        "slope_ci_lo": float(np.percentile(sl_b, 2.5)), "slope_ci_hi": float(np.percentile(sl_b, 97.5)),
        "slope_null_mean": float(sl_null.mean()), "p_slope": float(T.p_high(sl_null, sl)),
        "typical_error": err,
        "typical_error_ci_lo": float(np.percentile(err_b, 2.5)), "typical_error_ci_hi": float(np.percentile(err_b, 97.5)),
        "null_typical_error": float(err_null.mean()), "p_error_smaller": float(T.p_low(err_null, err)),
    }
    return out, sl_b


def mean_bias(ix, p, tx, pt_med, bmed):
    """Mean over the truths ix of (median best-run parameter - true parameter), point estimate and bootstrap draws.
    Units: Ks as is; f and cc as log10."""
    col = PCOL[p]
    x = tx[ix, col]
    return float(np.mean(pt_med[ix, col] - x)), np.mean(bmed[:, ix, col] - x, axis=1)


def ci(v):
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def build_tables(data, truths, tx, fracs, sets, pt, boot, null, half_width, dist, nwithin):
    keys, cells, centre = find_cells(truths)
    per_rows, cell_rows, pool_rows = [], [], []
    all_ix = list(range(len(truths)))
    off_keys = [k for k in keys if k != centre]
    for s in (SER_ON, SER_OFF):
        X, K = data[s]["X"], data[s]["K"]
        for frac in fracs:
            P_ = pt[(s, frac)]
            bt, nl = boot[(s, frac)], null[(s, frac)]
            # ---- per truth ----------------------------------------------------------
            for ti, t in enumerate(truths):
                idx = sets[(s, frac)][ti]
                tl = tx[ti, 2]
                share_all = float(np.mean(np.abs(X[:, 2] - tl) <= half_width + 1e-12))
                within = float(np.mean(np.abs(X[idx, 2] - tl) <= half_width + 1e-12))
                row = {"truth_id": t["id"], "truth_Ks": t["Ks"], "truth_f": t["f"], "truth_cc": t["cc"],
                       "series": s, "top_frac": frac, "n_top": int(len(idx)),
                       "best_kge_2012": float(np.nanmax(K[ti])),
                       "box_pos_Ks": float(box_pos(t["Ks"], "Ks_mult")), "box_pos_f": float(box_pos(t["f"], "f_RS_abs")),
                       "box_pos_cc": float(box_pos(t["cc"], "channelconductivity_mmhr")),
                       "nearest_sample_dist": float(dist[ti]), "n_samples_within_0p15": int(nwithin[ti])}
                for p in PAR:
                    col = PCOL[p]
                    med, q10, q90 = P_["med"][ti, col], P_["q10"][ti, col], P_["q90"][ti, col]
                    lo, hi = ci(bt[:, ti, col])
                    nmean = float(np.mean(nl["med"][:, ti, col]))
                    if p == "Ks":
                        row.update({"median_Ks": med, "q10_Ks": q10, "q90_Ks": q90, "median_Ks_ci_lo": lo,
                                    "median_Ks_ci_hi": hi, "null_median_Ks": nmean, "Ks_bias": med - tx[ti, col],
                                    "p_Ks_closer": float(T.p_low(np.abs(nl["med"][:, ti, col] - tx[ti, col]),
                                                                 abs(med - tx[ti, col])))})
                    else:
                        row.update({"median_%s" % p: 10 ** med, "q10_%s" % p: 10 ** q10, "q90_%s" % p: 10 ** q90,
                                    "median_%s_ci_lo" % p: 10 ** lo, "median_%s_ci_hi" % p: 10 ** hi,
                                    "null_median_%s" % p: 10 ** nmean, "%s_log10_bias" % p: med - tx[ti, col],
                                    "%s_factor" % p: 10 ** (med - tx[ti, col]),
                                    "p_%s_closer" % p: float(T.p_low(np.abs(nl["med"][:, ti, col] - tx[ti, col]),
                                                                     abs(med - tx[ti, col])))})
                row.update({"cc_truth_in_middle80": float(P_["q10"][ti, 2] <= tl <= P_["q90"][ti, 2]),
                            "cc_share_near": within, "cc_share_near_all_runs": share_all,
                            "cc_enrichment": within / share_all if share_all > 0 else np.nan,
                            "p_cc_near": float(T.p_high(nl["within"][:, ti], within))})
                per_rows.append(row)

            # ---- per cell (cc tracking) and pooled -----------------------------------------
            cell_stats, cell_slope_b, cell_bias_b = {}, {}, {}
            for k in keys:
                ix = cells[k]
                st, sl_b = track(ix, "cc", tx, P_["med"], bt["med"] if False else bt, nl["med"]) if False else \
                    track(ix, "cc", tx, P_["med"], bt, nl["med"])
                row = {"series": s, "top_frac": frac, "cell_Ks": k[0], "cell_f": k[1], "is_centre": bool(k == centre),
                       "n_truths": len(ix),
                       "nearest_sample_dist_mean": float(np.mean(dist[ix])), "nearest_sample_dist_max": float(np.max(dist[ix]))}
                if st is not None:
                    row.update({"cc_" + a: b for a, b in st.items() if a != "n_truths"})
                    cell_slope_b[k] = sl_b
                for p in PAR:
                    b, bb = mean_bias(ix, p, tx, P_["med"], bt)
                    lo, hi = ci(bb)
                    row.update({"%s_mean_bias" % p: b, "%s_mean_bias_ci_lo" % p: lo, "%s_mean_bias_ci_hi" % p: hi})
                    cell_bias_b[(k, p)] = bb
                cell_stats[k] = row
            if centre is not None:
                for k in keys:
                    row = cell_stats[k]
                    if k == centre:
                        continue
                    if "cc_slope" in row and "cc_slope" in cell_stats[centre]:
                        dd = cell_slope_b[k] - cell_slope_b[centre]
                        row["cc_slope_minus_centre"] = row["cc_slope"] - cell_stats[centre]["cc_slope"]
                        row["cc_slope_minus_centre_ci_lo"], row["cc_slope_minus_centre_ci_hi"] = ci(dd)
                    for p in PAR:
                        dd = cell_bias_b[(k, p)] - cell_bias_b[(centre, p)]
                        row["%s_bias_minus_centre" % p] = row["%s_mean_bias" % p] - cell_stats[centre]["%s_mean_bias" % p]
                        row["%s_bias_minus_centre_ci_lo" % p], row["%s_bias_minus_centre_ci_hi" % p] = ci(dd)
            cell_rows.extend(cell_stats[k] for k in keys)

            for p in PAR:
                st, sl_b = track(all_ix, p, tx, P_["med"], bt, nl["med"])
                if st is None:
                    continue
                row = {"series": s, "top_frac": frac, "param": p}
                row.update(st)
                if p == "cc":
                    cs = [cell_slope_b[k] for k in off_keys if k in cell_slope_b]
                    if centre is not None and centre in cell_slope_b:
                        row["centre_slope"] = cell_stats[centre]["cc_slope"]
                        row["centre_slope_ci_lo"], row["centre_slope_ci_hi"] = ci(cell_slope_b[centre])
                        if cs:
                            off_b = np.mean(np.vstack(cs), axis=0)
                            off_pt = float(np.mean([cell_stats[k]["cc_slope"] for k in off_keys if k in cell_slope_b]))
                            row["offcentre_mean_slope"] = off_pt
                            row["offcentre_mean_slope_ci_lo"], row["offcentre_mean_slope_ci_hi"] = ci(off_b)
                            dd = off_b - cell_slope_b[centre]
                            row["offcentre_minus_centre"] = off_pt - cell_stats[centre]["cc_slope"]
                            row["offcentre_minus_centre_ci_lo"], row["offcentre_minus_centre_ci_hi"] = ci(dd)
                pool_rows.append(row)
    return pd.DataFrame(per_rows), pd.DataFrame(cell_rows), pd.DataFrame(pool_rows)


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------
def print_truth_table(pt_df, main_frac, factor):
    on = pt_df[(pt_df.series == SER_ON) & (pt_df.top_frac == main_frac)].sort_values(["truth_Ks", "truth_f", "truth_cc"])
    off = pt_df[(pt_df.series == SER_OFF) & (pt_df.top_frac == main_frac)].set_index("truth_id")
    n = int(on.n_top.iloc[0])
    lo10, hi10 = 10 ** (R.LOG_LO + 0.1 * (R.LOG_HI - R.LOG_LO)), 10 ** (R.LOG_LO + 0.9 * (R.LOG_HI - R.LOG_LO))
    print("\n" + "=" * 150)
    print("1. WHERE DO THE BEST RUNS PUT cc, Ks AND f?  (cc-ON runs, best %s = %d runs by KGE_2012; 'near' = within a factor of %g)"
          % (cut_label(main_frac), n, factor))
    print("=" * 150)
    print("  No information (runs that ignore cc) would give a cc median near %s mm/hr (middle 80%% about %s-%s) whatever the true cc."
          % (R.fcc(10 ** R.CENTER_LOG), R.fcc(lo10), R.fcc(hi10)))
    print("%5s %7s %6s | %6s %7s %5s | %7s [%6s-%6s] %8s %5s %-11s %6s | %6s | %6s %6s | %7s %6s"
          % ("true", "true", "true", "best", "nearest", "n<", "cc", "10%", "90%", "vs truth", "in80%", "near/chance", "p",
             "ctl cc", "Ks", "vs", "f", "vs"))
    print("%5s %7s %6s | %6s %7s %5s | %7s %6s %6s %8s %5s %-11s %6s | %6s | %6s %6s | %7s %6s"
          % ("Ks", "f", "cc", "KGE", "sample", ".15", "median", "", "", "(x)", "", "(best/all)", "", "median", "med", "truth",
             "med", "(x)"))
    last = None
    for _, a in on.iterrows():
        cell = (a.truth_Ks, a.truth_f)
        if last is not None and cell != last:
            print("")
        last = cell
        c = off.loc[a.truth_id] if a.truth_id in off.index else None
        print("%5g %7g %6g | %6.3f %7.3f %5d | %7s [%6s-%6s] %8s %5s %-11s %6s | %6s | %6.2f %+6.2f | %7.4f %6s"
              % (a.truth_Ks, a.truth_f, a.truth_cc, a.best_kge_2012, a.nearest_sample_dist, a.n_samples_within_0p15,
                 R.fcc(a.median_cc), R.fcc(a.q10_cc), R.fcc(a.q90_cc), R.ffactor(a.cc_factor),
                 "yes" if a.cc_truth_in_middle80 == 1 else "NO",
                 "%d%% / %d%%" % (round(100 * a.cc_share_near), round(100 * a.cc_share_near_all_runs)),
                 T.fmt_p(a.p_cc_near), R.fcc(c.median_cc) if c is not None else "n/a", a.median_Ks, a.Ks_bias,
                 a.median_f, R.ffactor(a.f_factor)))
    print("  nearest sample / n<.15: distance (box units) from this truth to the closest of the sampled runs / how many runs lie within 0.15 of it.")
    print("  vs truth (x): median of the best runs / true value (x1.00 = centred on the truth).  in80%: is the true cc inside the "
          "middle 80% of the best runs' cc?")
    print("  near/chance: share of the best runs within a factor of %g of the true cc / the share of ALL runs that are (what chance gives)." % factor)
    print("  p: could that much enrichment arise if cc were shuffled among the runs?   ctl cc: median cc of the best control runs "
          "(it should wander near %s)." % R.fcc(10 ** R.CENTER_LOG))
    print("  Ks med / f med: medians of the best runs; Ks 'vs truth' is the difference (best-run median minus true Ks), f is a factor.")


def print_cells(cell_df, pool_df, fracs):
    print("\n" + "=" * 150)
    print("2. DOES cc TRACK THE TRUTH IN EVERY CELL?  slope of (median best-run cc) against (true cc), both on log scales, over each "
          "cell's cc values")
    print("=" * 150)
    print("  1 = follows the truth one-for-one;  0 = ignores it.  The control (cc inert) is the reference for 0.  'vs centre' compares "
          "each cell with the original centre cell on the same re-drawn runs.")
    for frac in fracs:
        print("\n  -- best %s of runs ------------------------------------------------------------------------------------------------"
              % cut_label(frac))
        print("  %-16s %2s | %-22s %7s | %-22s | %-30s | %6s %6s"
              % ("cell", "n", "ON slope [95% interval]", "p", "control slope [95%]", "ON minus centre [95%]", "nearest",
                 "cc x"))
        on = cell_df[(cell_df.series == SER_ON) & (cell_df.top_frac == frac)]
        off = cell_df[(cell_df.series == SER_OFF) & (cell_df.top_frac == frac)].set_index(["cell_Ks", "cell_f"])
        for _, a in on.iterrows():
            c = off.loc[(a.cell_Ks, a.cell_f)] if (a.cell_Ks, a.cell_f) in off.index else None
            name = cell_name((a.cell_Ks, a.cell_f)) + (" *" if a.is_centre else "")
            if "cc_slope" in a and np.isfinite(a.get("cc_slope", np.nan)):
                s_on = iv(a.cc_slope, a.cc_slope_ci_lo, a.cc_slope_ci_hi)
                p_on = T.fmt_p(a.cc_p_slope)
            else:
                s_on, p_on = "n/a (needs >= %d cc values)" % MIN_TRACK, ""
            if c is not None and "cc_slope" in c and np.isfinite(c.get("cc_slope", np.nan)):
                s_off = iv(c.cc_slope, c.cc_slope_ci_lo, c.cc_slope_ci_hi)
            else:
                s_off = "n/a"
            if a.is_centre:
                d = "(reference)"
            elif "cc_slope_minus_centre" in a and np.isfinite(a.get("cc_slope_minus_centre", np.nan)):
                d = "%+5.2f [%+5.2f, %+5.2f]%s" % (a.cc_slope_minus_centre, a.cc_slope_minus_centre_ci_lo,
                                                   a.cc_slope_minus_centre_ci_hi,
                                                   "  *" if (a.cc_slope_minus_centre_ci_lo > 0 or a.cc_slope_minus_centre_ci_hi < 0) else "")
            else:
                d = "n/a"
            print("  %-16s %2d | %-22s %7s | %-22s | %-30s | %6.3f %6s"
                  % (name, a.n_truths, s_on, p_on, s_off, d, a.nearest_sample_dist_mean,
                     R.ffactor(10 ** a.cc_mean_bias)))
        print("")
        for s in (SER_ON, SER_OFF):
            g = pool_df[(pool_df.series == s) & (pool_df.top_frac == frac) & (pool_df.param == "cc")]
            if g.empty:
                continue
            r = g.iloc[0]
            line = ("  %-16s %2d | %+5.2f [%+5.2f, %+5.2f]  p %-6s  no-info %+5.2f | typical cc error x%.2f (no-info x%.2f)"
                    % ("ALL TRUTHS, " + s.replace("cc ", ""), r.n_truths, r.slope, r.slope_ci_lo, r.slope_ci_hi,
                       T.fmt_p(r.p_slope), r.slope_null_mean, 10 ** r.typical_error, 10 ** r.null_typical_error))
            print(line)
            if s == SER_ON and "centre_slope" in r.index and np.isfinite(r.get("centre_slope", np.nan)):
                print("  %-16s    | centre cell %+5.2f [%+5.2f, %+5.2f];  mean of the %d off-centre cells %+5.2f [%+5.2f, %+5.2f];  "
                      "difference %+5.2f [%+5.2f, %+5.2f]"
                      % ("", r.centre_slope, r.centre_slope_ci_lo, r.centre_slope_ci_hi,
                         int((cell_df[(cell_df.series == SER_ON) & (cell_df.top_frac == frac)].is_centre == False).sum()),
                         r.offcentre_mean_slope, r.offcentre_mean_slope_ci_lo, r.offcentre_mean_slope_ci_hi,
                         r.offcentre_minus_centre, r.offcentre_minus_centre_ci_lo, r.offcentre_minus_centre_ci_hi))
    print("\n  * = original centre cell (Ks %g, f %g).  In the 'ON minus centre' column, a trailing * marks an interval that excludes 0 "
          "(the cell really differs from the centre)." % (STANDARD_KS, STANDARD_F))
    print("  cc x: mean over the cell's truths of (median best-run cc / true cc); x1.00 = unbiased.  nearest: mean distance from the "
          "cell's truths to their nearest sample (box units).")
    print("  ! = the slope lies outside its own bootstrap interval. That happens when the best runs are few and nearly tied (typical "
          "for the control); trust the p-value and the pattern across cells over that interval.")
    print("  Each cell has only a few cc values, so a per-cell slope is rough. The 'ALL TRUTHS' line pools every truth.")


def print_params(cell_df, pool_df, main_frac):
    print("\n" + "=" * 150)
    print("3. DO THE BEST RUNS RECOVER Ks AND f AS THE TRUTH MOVES?   (best %s of runs)" % cut_label(main_frac))
    print("=" * 150)
    print("  Pooled tracking slope across all truths (Ks on a linear scale, f on a log scale). 1 = follows; 0 = ignores.")
    print("  %-6s %-8s %2s | %-26s %7s | %-9s %9s %7s"
          % ("param", "series", "n", "slope [95% interval]", "p", "typ. err", "no-info", "p"))
    for p in ("Ks", "f"):
        for s in (SER_ON, SER_OFF):
            g = pool_df[(pool_df.series == s) & (pool_df.top_frac == main_frac) & (pool_df.param == p)]
            if g.empty:
                print("  %-6s %-8s    | n/a (needs at least %d truths and 2 distinct true values)" % (p, s.replace("cc ", ""), MIN_TRACK))
                continue
            r = g.iloc[0]
            print("  %-6s %-8s %2d | %+5.2f [%+5.2f, %+5.2f]        %7s | %-9s %9s %7s"
                  % (p, s.replace("cc ", ""), r.n_truths, r.slope, r.slope_ci_lo, r.slope_ci_hi, T.fmt_p(r.p_slope),
                     param_err_fmt(p, r.typical_error), param_err_fmt(p, r.null_typical_error), T.fmt_p(r.p_error_smaller)))
    print("\n  Mean bias of the best runs' median in each cell (average over the cell's truths): Ks as a difference, f and cc as factors.")
    print("  %-16s | %-26s | %-26s | %-26s" % ("cell", "Ks (ON / control)", "f (ON / control)", "cc (ON / control)"))
    on = cell_df[(cell_df.series == SER_ON) & (cell_df.top_frac == main_frac)]
    off = cell_df[(cell_df.series == SER_OFF) & (cell_df.top_frac == main_frac)].set_index(["cell_Ks", "cell_f"])
    for _, a in on.iterrows():
        c = off.loc[(a.cell_Ks, a.cell_f)]
        name = cell_name((a.cell_Ks, a.cell_f)) + (" *" if a.is_centre else "")
        print("  %-16s | %+7.2f / %+-7.2f         | %-9s / %-9s      | %-9s / %-9s"
              % (name, a.Ks_mean_bias, c.Ks_mean_bias, R.ffactor(10 ** a.f_mean_bias), R.ffactor(10 ** c.f_mean_bias),
                 R.ffactor(10 ** a.cc_mean_bias), R.ffactor(10 ** c.cc_mean_bias)))
    print("  The control cannot represent channel loss, so its Ks and f absorb the missing loss: its biases show how far Ks and f "
          "alone can drift to imitate a truth that has loss.")


def summary_lines(pt_df, cell_df, pool_df, main_frac):
    """Plain statements assembled from the numbers above. Descriptive only."""
    out = []
    on_c = pool_df[(pool_df.series == SER_ON) & (pool_df.top_frac == main_frac) & (pool_df.param == "cc")]
    off_c = pool_df[(pool_df.series == SER_OFF) & (pool_df.top_frac == main_frac) & (pool_df.param == "cc")]
    if not on_c.empty:
        r = on_c.iloc[0]
        ctl = ("; control %+.2f [%+.2f, %+.2f]" % (off_c.iloc[0].slope, off_c.iloc[0].slope_ci_lo, off_c.iloc[0].slope_ci_hi)
               if not off_c.empty else "")
        out.append("Pooled over all %d truths, the best runs' cc follows the true cc with slope %+.2f [%+.2f, %+.2f]%s "
                   "(0 = ignores the truth, 1 = follows it one-for-one). Typical cc error x%.2f against x%.2f for no information."
                   % (r.n_truths, r.slope, r.slope_ci_lo, r.slope_ci_hi, ctl, 10 ** r.typical_error, 10 ** r.null_typical_error))
        if "centre_slope" in on_c.columns and np.isfinite(on_c.iloc[0].get("centre_slope", np.nan)):
            out.append("The original centre cell gives %+.2f [%+.2f, %+.2f]; the off-centre cells average %+.2f [%+.2f, %+.2f]; "
                       "difference %+.2f [%+.2f, %+.2f]%s."
                       % (r.centre_slope, r.centre_slope_ci_lo, r.centre_slope_ci_hi, r.offcentre_mean_slope,
                          r.offcentre_mean_slope_ci_lo, r.offcentre_mean_slope_ci_hi, r.offcentre_minus_centre,
                          r.offcentre_minus_centre_ci_lo, r.offcentre_minus_centre_ci_hi,
                          " (the interval excludes 0)" if (r.offcentre_minus_centre_ci_lo > 0 or r.offcentre_minus_centre_ci_hi < 0)
                          else " (the interval includes 0: no detectable difference)"))
    on = cell_df[(cell_df.series == SER_ON) & (cell_df.top_frac == main_frac)]
    have = on[on.get("cc_slope", pd.Series(dtype=float)).notna()] if "cc_slope" in on.columns else on.iloc[0:0]
    off = have[~have.is_centre]
    has_centre = bool(on.is_centre.any()) if len(on) else False
    if len(off):
        pos = int((off.cc_slope_ci_lo > 0).sum())
        out.append("%s cells whose cc slope has a 95%% interval entirely above 0: %d of %d (slope range %+.2f to %+.2f)."
                   % ("Off-centre" if has_centre else "Of the", pos, len(off), off.cc_slope.min(), off.cc_slope.max()))
        if "cc_slope_minus_centre_ci_lo" in off.columns:
            lower = int((off.cc_slope_minus_centre_ci_hi < 0).sum())
            higher = int((off.cc_slope_minus_centre_ci_lo > 0).sum())
            out.append("Compared with the centre cell, %d off-centre cell(s) track detectably WORSE, %d detectably BETTER, and %d "
                       "are not distinguishable." % (lower, higher, len(off) - lower - higher))
    if not len(have):
        out.append("No cell has %d or more cc values, so no per-cell slopes (and no centre-versus-off-centre comparison) could be "
                   "computed; only the pooled lines are reported." % MIN_TRACK)
    for p in ("Ks", "f"):
        g = pool_df[(pool_df.series == SER_ON) & (pool_df.top_frac == main_frac) & (pool_df.param == p)]
        if not g.empty:
            r = g.iloc[0]
            out.append("%s pooled tracking slope %+.2f [%+.2f, %+.2f] (typical error %s against %s for no information)."
                       % (p, r.slope, r.slope_ci_lo, r.slope_ci_hi, param_err_fmt(p, r.typical_error),
                          param_err_fmt(p, r.null_typical_error)))
    d = pt_df[(pt_df.series == SER_ON) & (pt_df.top_frac == main_frac)]
    if len(d):
        cen = d[np.isclose(d.truth_Ks, STANDARD_KS) & np.isclose(d.truth_f, STANDARD_F)]
        oth = d[~(np.isclose(d.truth_Ks, STANDARD_KS) & np.isclose(d.truth_f, STANDARD_F))]
        txt = "Coverage: the nearest sample to a truth is %.3f box units away on median (range %.3f to %.3f)" % (
            d.nearest_sample_dist.median(), d.nearest_sample_dist.min(), d.nearest_sample_dist.max())
        if len(cen) and len(oth):
            txt += "; centre truths %.3f, off-centre %.3f" % (cen.nearest_sample_dist.median(), oth.nearest_sample_dist.median())
        if len(d) >= 6 and d.nearest_sample_dist.nunique() > 2 and d.cc_log10_bias.abs().nunique() > 2:
            rho = stats.spearmanr(d.nearest_sample_dist, d.cc_log10_bias.abs())[0]
            txt += "; Spearman rho(nearest-sample distance, |cc error|) = %+.2f" % rho
        out.append(txt + ".")
    return out


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def _setup_log_axis(ax, lim):
    import matplotlib
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    for axis in (ax.get_xaxis(), ax.get_yaxis()):
        axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
        axis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def fig_recovery(pt_df, keys, cells, truths, frac, factor, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    Ks_vals = sorted({k[0] for k in keys})
    f_vals = sorted({k[1] for k in keys})
    nr, nc = len(Ks_vals), len(f_vals)
    fig, axes = plt.subplots(nr, nc, figsize=(3.3 * nc + 0.6, 3.1 * nr + 1.3), facecolor=T.SURFACE, sharex=True, sharey=True,
                             squeeze=False)
    lo, hi = 10 ** R.LOG_LO, 10 ** R.LOG_HI
    lim = (lo * 0.75, hi * 1.3)
    d = pt_df[pt_df.top_frac == frac]
    cc_all = sorted({t["cc"] for t in truths})
    for ri, ks in enumerate(Ks_vals):
        for ci_, fv in enumerate(f_vals):
            ax = axes[ri][ci_]
            T._style(ax)
            if (ks, fv) not in cells:
                ax.set_visible(False)
                continue
            ax.axhspan(10 ** (R.LOG_LO + 0.1 * (R.LOG_HI - R.LOG_LO)), 10 ** (R.LOG_LO + 0.9 * (R.LOG_HI - R.LOG_LO)),
                       color=T.COL_NULL, alpha=0.16, linewidth=0, zorder=1)
            ax.axhline(10 ** R.CENTER_LOG, color=T.COL_NULL, linestyle=":", linewidth=1.2, zorder=2)
            ax.plot(lim, lim, color=T.INK, linestyle="--", linewidth=1.0, zorder=2)
            ids = [truths[i]["id"] for i in cells[(ks, fv)]]
            for s, col, off_, mk, filled in ((SER_OFF, T.COL_NULL, +1, "s", False), (SER_ON, T.COL_CC, -1, "o", True)):
                g = d[(d.series == s) & d.truth_id.isin(ids)].sort_values("truth_cc")
                x = g.truth_cc.values * 10 ** (off_ * 0.012)
                ax.vlines(x, g.q10_cc, g.q90_cc, color=col, linewidth=2.0 if s == SER_ON else 1.4, zorder=3)
                ax.plot(x, g.median_cc.values, linestyle="none", marker=mk, markersize=5.5,
                        markerfacecolor=col if filled else T.SURFACE, markeredgecolor=col, markeredgewidth=1.5, zorder=4)
            _setup_log_axis(ax, lim)
            ax.set_xticks(cc_all)
            ax.set_xticklabels(["%g" % v for v in cc_all], fontsize=7.5)
            ax.set_yticks([30, 100, 300, 1000])
            centre = bool(np.isclose(ks, STANDARD_KS) and np.isclose(fv, STANDARD_F))
            ax.set_title("true Ks %g, f %g%s" % (ks, fv, "  (original centre)" if centre else ""), color=T.INK,
                         fontsize=9, loc="left", fontweight="bold" if centre else "normal")
            dd = d[(d.series == SER_ON) & d.truth_id.isin(ids)].sort_values("truth_cc")
            ax.text(0.97, 0.04, "nearest sample: " + " / ".join("%.2f" % v for v in dd.nearest_sample_dist.values),
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, color=T.INK2)
            if ri == nr - 1:
                ax.set_xlabel("true cc (mm/hr)", color=T.INK2, fontsize=8.5)
            if ci_ == 0:
                ax.set_ylabel("cc of the best runs (mm/hr)", color=T.INK2, fontsize=8.5)
    n_top = int(d[d.series == SER_ON].n_top.iloc[0])
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5.5, linewidth=2.0,
               label="channel loss ON: median and middle 80% of the best runs"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=5.5, markerfacecolor=T.SURFACE, markeredgewidth=1.5,
               linewidth=1.4, label="control (cc inert): the same, for runs that cannot feel cc"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=1.0, label="perfect recovery"),
        Patch(facecolor=T.COL_NULL, alpha=0.16, label="no information: middle 80% of the sampled range (dotted = its middle)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Do the best runs choose the cc that made the truth, wherever the true Ks and f are?  "
                 "(best %s of runs, n = %d)" % (cut_label(frac), n_top), color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.09, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


def fig_cell_slopes(cell_df, pool_df, keys, centre, fracs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    labels = [cell_name(k) + ("  (original centre)" if k == centre else "") for k in keys] + ["all truths pooled"]
    ny = len(labels)
    nf = len(fracs)
    fig, axes = plt.subplots(1, nf, figsize=(4.9 * nf + 2.2, 0.42 * ny + 2.4), facecolor=T.SURFACE, sharey=True, sharex=True, squeeze=False)
    for ci_, frac in enumerate(fracs):
        ax = axes[0][ci_]
        T._style(ax)
        ax.axvline(0.0, color=T.INK, linestyle="--", linewidth=0.9, zorder=1)
        ax.axvline(1.0, color=T.COL_NULL, linestyle=":", linewidth=1.1, zorder=1)
        if centre is not None:
            ax.axhspan(keys.index(centre) - 0.5, keys.index(centre) + 0.5, color=T.COL_CC, alpha=0.08, linewidth=0, zorder=0)
        for s, col, off_, mk, filled in ((SER_ON, T.COL_CC, -0.14, "o", True), (SER_OFF, T.COL_NULL, +0.14, "s", False)):
            for yi, k in enumerate(keys + ["ALL"]):
                if k == "ALL":
                    g = pool_df[(pool_df.series == s) & (pool_df.top_frac == frac) & (pool_df.param == "cc")]
                    if g.empty:
                        continue
                    r = g.iloc[0]
                    sl, lo, hi = r.slope, r.slope_ci_lo, r.slope_ci_hi
                else:
                    g = cell_df[(cell_df.series == s) & (cell_df.top_frac == frac) & (cell_df.cell_Ks == k[0]) & (cell_df.cell_f == k[1])]
                    if g.empty or "cc_slope" not in g.columns or not np.isfinite(g.iloc[0].cc_slope):
                        continue
                    r = g.iloc[0]
                    sl, lo, hi = r.cc_slope, r.cc_slope_ci_lo, r.cc_slope_ci_hi
                y = yi + off_
                ax.hlines(y, lo, hi, color=col, linewidth=2.0 if s == SER_ON else 1.4, zorder=3)
                ax.plot([sl], [y], linestyle="none", marker=mk, markersize=5.5, markerfacecolor=col if filled else T.SURFACE,
                        markeredgecolor=col, markeredgewidth=1.5, zorder=4)
        ax.set_yticks(range(ny))
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.axhline(ny - 1.5, color=T.GRID, linewidth=1.0, zorder=1)
        ax.set_title("Best %s of runs" % cut_label(frac), color=T.INK, fontsize=9.5, loc="left")
        ax.set_xlabel("cc tracking slope (95% interval)", color=T.INK2, fontsize=9)
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5.5, linewidth=2.0, label="channel loss ON"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=5.5, markerfacecolor=T.SURFACE, markeredgewidth=1.5,
               linewidth=1.4, label="control (cc inert)"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=0.9, label="0 = the best runs ignore the true cc"),
        Line2D([0], [0], color=T.COL_NULL, linestyle=":", linewidth=1.1, label="1 = they follow it one-for-one"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Does cc tracking hold away from the centre?  (slope of best-run cc against true cc, per (Ks, f) cell)",
                 color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


def fig_params(pt_df, frac, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    d = pt_df[pt_df.top_frac == frac]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 5.3), facecolor=T.SURFACE)
    spec = [("Ks", "truth_Ks", "median_Ks", "true Ks multiplier", "Ks multiplier of the best runs", False,
             P.BOX["Ks_mult"], [4, 5, 6, 7, 8, 9, 10]),
            ("f", "truth_f", "median_f", "true f_RS_abs", "f_RS_abs of the best runs", True,
             P.BOX["f_RS_abs"], [0.005, 0.01, 0.02, 0.03])]
    for ax, (nm, xc, yc, xl, yl, log, box, yticks) in zip(axes, spec):
        T._style(ax)
        lo, hi = box[0], box[1]
        if log:
            pad = (hi / lo) ** 0.06
            lim = (lo / pad, hi * pad)
            ax.set_xscale("log")
            ax.set_yscale("log")
        else:
            pad = 0.06 * (hi - lo)
            lim = (lo - pad, hi + pad)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.axhspan(lim[0], lo, color=T.COL_NULL, alpha=0.16, linewidth=0, zorder=1)
        ax.axhspan(hi, lim[1], color=T.COL_NULL, alpha=0.16, linewidth=0, zorder=1)
        ax.plot(lim, lim, color=T.INK, linestyle="--", linewidth=1.0, zorder=2)
        xs = sorted(d[xc].unique())
        for s_, col, sgn, mk, filled in ((SER_OFF, T.COL_NULL, +1, "s", False), (SER_ON, T.COL_CC, -1, "o", True)):
            g = d[d.series == s_]
            rng = np.random.default_rng(3)
            for xv in xs:
                gg = g[g[xc] == xv]
                if log:
                    xx = xv * (1 + 0.02 * sgn + rng.uniform(-0.012, 0.012, len(gg)))
                else:
                    xx = xv + 0.08 * sgn + rng.uniform(-0.05, 0.05, len(gg))
                ax.plot(xx, gg[yc].values, linestyle="none", marker=mk, markersize=5,
                        markerfacecolor=col if filled else T.SURFACE, markeredgecolor=col, markeredgewidth=1.3,
                        alpha=0.9, zorder=3)
        ax.set_xticks(xs)
        ax.set_xticklabels(["%g" % v for v in xs], fontsize=8)
        ax.set_yticks(yticks)
        ax.set_yticklabels(["%g" % v for v in yticks], fontsize=8)
        if log:
            for axis in (ax.xaxis, ax.yaxis):
                axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel(xl, color=T.INK2, fontsize=9)
        ax.set_ylabel(yl, color=T.INK2, fontsize=9)
        ax.set_title("%s (one dot per truth)" % nm, color=T.INK, fontsize=9.5, loc="left")
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5, linestyle="none",
               label="channel loss ON: median of the best runs"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=5, markerfacecolor=T.SURFACE, markeredgewidth=1.3,
               linestyle="none", label="control (cc inert)"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=1.0, label="perfect recovery"),
        Patch(facecolor=T.COL_NULL, alpha=0.16, label="outside the sampled range"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Do the best runs recover the true Ks and f?  (best %s of runs)" % cut_label(frac),
                 color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110 -- is cc (and Ks, f) recovered when the true Ks and f are off-centre?")
    ap.add_argument("--top_fracs", type=float, nargs="+", default=[0.20, 0.10, 0.05],
                    help="Top fractions by KGE_2012; the first is the main one (default 0.20 0.10 0.05)")
    ap.add_argument("--within", type=float, default=2.0,
                    help="'Near the truth' means within this FACTOR of the true cc (default 2)")
    ap.add_argument("--boot", type=int, default=1000, help="Joint-bootstrap re-draws (default 1000)")
    ap.add_argument("--perm", type=int, default=5000, help="Permutation shuffles (default 5000)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--expect_n", type=int, default=250)
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--continue_on_fail", action="store_true")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()
    if args.within <= 1.0:
        sys.exit("--within must be greater than 1 (it is a factor, e.g. 2 = between half and double).")
    half_width = float(np.log10(args.within))
    fracs = args.top_fracs
    main_frac = fracs[0]

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    loc_dir = summary_dir.parent.parent / LOC_DIRNAME
    out_dir = summary_dir / OUT_DIRNAME
    print("\n" + "=" * 78 + "\nSeries 110 -- truth-location test: is cc recovered away from the centre?\n" + "=" * 78)
    print("Reading from: %s" % summary_dir)

    df = load_long(summary_dir)
    truths = truth_list(df)
    ck = run_checks(df, truths, summary_dir, loc_dir, args.expect_n)
    P.print_checks(ck)
    if ck.n("FAIL") and not args.continue_on_fail:
        sys.exit("\nSTOPPED: an input check failed. Nothing was analysed or written. (--continue_on_fail to override.)")

    data = build_tensors(df, truths)
    tx = truth_coords(truths)
    keys, cells, centre = find_cells(truths)
    dist, nwithin = nearest_sample(data[SER_ON]["X"], tx)
    print("\n%d truth(s) in %d cell(s)%s, %d runs per series, cuts %s, %d shuffles, %d bootstrap re-draws, near = within x%g"
          % (len(truths), len(keys), " (original centre present)" if centre is not None else " (no original centre)",
             data[SER_ON]["X"].shape[0], [cut_label(f) for f in fracs], args.perm, args.boot, args.within))

    print("\nRunning the shuffles and bootstrap ...", flush=True)
    sets = point_sets(data, fracs)
    pt = point_stats(data, sets, fracs)
    boot = run_boot(data, fracs, args.boot, args.seed)
    null = run_perm(data, sets, truths, fracs, args.perm, args.seed, half_width)
    pt_df, cell_df, pool_df = build_tables(data, truths, tx, fracs, sets, pt, boot, null, half_width, dist, nwithin)

    print_truth_table(pt_df, main_frac, args.within)
    print_cells(cell_df, pool_df, fracs)
    print_params(cell_df, pool_df, main_frac)
    print("\n" + "=" * 150)
    print("4. WHAT THE NUMBERS SAY   (best %s of runs; descriptive, assembled from the tables above)" % cut_label(main_frac))
    print("=" * 150)
    for line in summary_lines(pt_df, cell_df, pool_df, main_frac):
        print("  - " + line)
    print("  Caveat: one storm, %d runs, routing pinned at truth, noise-free truths, and only a few cc values per cell. A cell's slope "
          "is rough; the pattern across cells and the pooled lines are the better guide." % data[SER_ON]["X"].shape[0])

    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ck.rows).to_csv(out_dir / "location_input_checks_110.csv", index=False)
    pt_df.to_csv(out_dir / "location_per_truth_110.csv", index=False)
    cell_df.to_csv(out_dir / "location_cells_110.csv", index=False)
    pool_df.to_csv(out_dir / "location_pooled_110.csv", index=False)
    if not args.no_plots:
        for name, fn in (
                ("fig_location_recovery_cc_110.png", lambda p: fig_recovery(pt_df, keys, cells, truths, main_frac, args.within, p)),
                ("fig_location_cell_slopes_110.png", lambda p: fig_cell_slopes(cell_df, pool_df, keys, centre, fracs, p)),
                ("fig_location_Ks_f_recovery_110.png", lambda p: fig_params(pt_df, main_frac, p))):
            try:
                fn(out_dir / name)
            except Exception as e:
                print("\n  (figure %s skipped: %s -- the CSVs are saved regardless)" % (name, e))

    def md5_of(p):
        return hashlib.md5(Path(p).read_bytes()).hexdigest()
    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "pandas": pd.__version__, "numpy": np.__version__,
            "inputs": {n: md5_of(summary_dir / n) for n in (P.STAGE2_NAME, P.CONTROL_NAME)},
            "rescored_long": md5_of(summary_dir / RS_DIRNAME / LONG_NAME),
            "settings": {k: (v if not isinstance(v, Path) else str(v)) for k, v in vars(args).items()},
            "truths": [{"id": t["id"], "Ks": t["Ks"], "f": t["f"], "cc": t["cc"]} for t in truths],
            "cells": [list(k) for k in keys], "centre_present": centre is not None, "checks": ck.rows}
    (out_dir / "PROVENANCE_location_cc_110.json").write_text(json.dumps(prov, indent=2, default=str))

    print("\nSaved to: %s" % out_dir)
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results."
              % ck.n("FAIL"))


if __name__ == "__main__":
    main()
