"""
analyze_cc_stability_110.py
===========================
Series 110 -- is the cc tracking result stable, or does it depend on which 250 runs we happened to draw?

The earlier scripts report one answer from one set of 250 runs, and the "best 5%" is only 12 runs. This script
asks, using ONLY the runs you already have (no new tRIBS runs):

  1. SPLIT-HALF.  Split the 250 runs at random into two halves of 125. Do the two halves, each analysed on its own,
     give the same tracking slope? Repeated for many random splits. Two independent halves disagreeing by a lot means
     the result depends on the sample; agreeing means it is stable.
  2. LEARNING CURVE, same cut share.  Draw random subsets of 75, 125, 175 and 225 runs and keep the best 20%, 10% and 5% of
     each. Does the slope depend on how many runs we have? (If it does not, the headline numbers are not an accident of n.)
  3. LEARNING CURVE, same number of best runs.  Same subsets, but keep the SAME NUMBER of best runs (50, 25 and 12, the sizes
     of today's cuts). With more runs to choose from, the same number of best runs is a tighter, more selective cut. If the slope
     rises with the number of runs, a bigger set of runs would let the tight cuts do better than they do today.
  4. PROJECTION.  From how much the two halves disagree, estimate how wide the 95% interval of the slope would be with
     250, 500 or 1,000 runs, and compare the 250-run estimate with the bootstrap interval used in the other scripts.
     If the two agree, the bootstrap intervals you have been quoting are about the right size.

Everything uses KGE_2012 to pick the best runs, exactly as analyze_cc_recovery_110.py does, and the same cuts (best 20%,
10%, 5% of the runs used). It reads the stored tables only and writes into its own folder
(calibration_work/03_comparisons/summary_tables/stability_cc_110/).

Keep this file in the same folder as analyze_cc_recovery_110.py, analyze_cc_tests_110.py and analyze_cc_pca_110.py.
Run it from the lab/ directory.

VOCABULARY
-----------
  Subset            A random selection of the 250 runs (no run used twice in the same subset).
  Half              One of the two halves in a split-half test. The two halves of a split share no runs.
  Tracking slope    As before: slope of (median log10 cc of the best runs) against (log10 true cc) over the truths.
                    1 = follows the true cc one-for-one, 0 = ignores it.
  Spread (sd)       Standard deviation: how far a number typically wanders from its average when the sample changes.
  Learning curve    A result plotted against how many runs it used.
  95% half-width    About 2 spreads. The slope is "slope +/- half-width" with roughly 95% confidence.

HOW TO READ THE RESULTS
------------------------
  * Split-half: compare "halves differ by" with the slope. If the halves typically differ by less than half the slope and
    "both halves > 0" is near 100%, the result does not hinge on the sample. If they differ by as much as the slope itself, it does.
  * Learning curve, same cut share: a flat curve means the slope of a "best 20%" does not depend on how many runs there are.
  * Learning curve, same number of best runs: if the mean slope rises with the runs used, the same 12 (or 25, or 50) best runs
    become a tighter cut and track cc better, so more runs would push the tight-cut slopes UP. A flat curve means more runs would
    tighten the interval but leave the slope where it is. This is the curve that answers "would more runs help the headline?".
  * Small subsets: with few runs the best-5% and best-10% cuts are raised to the 10-run minimum (marked *). Do not read those
    rows as "5%" or "10%".
  * Projection: assumes the spread falls with the square root of the number of runs (4 times the runs, half the width). It is
    checked against the bootstrap at 250 runs. In simulations of designs like yours it slightly overstated the true width (by up to
    about 1.4 times) and got the shrinkage with more runs about right. It says nothing about the single storm, observation error,
    or pinned routing.

OUTPUTS   (calibration_work/03_comparisons/summary_tables/stability_cc_110/)
-------------------------------------------------------------------------------
  stability_learning_curve_110.csv   slope and typical error by number of runs used (both kinds of cut)
  stability_split_half_110.csv       agreement between the two halves
  stability_projection_110.csv       expected 95% half-width of the slope for 125 ... 1000 runs
  stability_input_checks_110.csv
  fig_stability_learning_110.png     slope against number of runs used, same cut share
  fig_stability_learning_count_110.png   slope against number of runs used, same number of best runs
  fig_stability_projection_110.png   expected interval width against number of runs
  PROVENANCE_stability_cc_110.json

USAGE (run from the lab/ directory)
------------------------------------
    python analyze_cc_stability_110.py                       # well under a minute
    python analyze_cc_stability_110.py --splits 2000 --subsets 500
    python analyze_cc_stability_110.py --no_plots

CAVEATS
--------
Subsets of one 250-run design are a stand-in for fresh independent designs: a real second set of runs would also vary the
space-filling pattern of the design. Slopes from small subsets are noisy and the best-5% cut is only a handful of runs there.
One storm, routing pinned at truth, exact observations (a best case).
"""

import argparse
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import analyze_cc_recovery_110 as R
except ImportError:
    sys.exit("This script needs analyze_cc_recovery_110.py, analyze_cc_tests_110.py and analyze_cc_pca_110.py in the same "
             "folder (it reuses their loader, input checks and tests). Put the files in lab/ and run from there.")

_NEEDED_R = ["T", "P", "slope_of", "tracking_sets", "run_boot", "MIN_TRACK_TRUTHS", "OUT_DIRNAME"]
_missing = [n for n in _NEEDED_R if not hasattr(R, n)]
if _missing:
    sys.exit("analyze_cc_recovery_110.py in this folder is an older version (missing: %s). Replace it with the current copy."
             % ", ".join(_missing))
T, P = R.T, R.P
_NEEDED_T = ["build_tensors", "n_top_for", "_style", "COL_CC", "COL_NULL", "INK", "INK2", "GRID", "SURFACE", "SER_ON", "SER_OFF", "OFF_TRUTH"]
_missing = [n for n in _NEEDED_T if not hasattr(T, n)]
if _missing:
    sys.exit("analyze_cc_tests_110.py in this folder is an older version (missing: %s). Replace it with the current copy."
             % ", ".join(_missing))

SER_ON, SER_OFF, OFF_TRUTH = T.SER_ON, T.SER_OFF, T.OFF_TRUTH
SERIES = (SER_ON, SER_OFF)
OUT_DIRNAME = "stability_cc_110"
MIN_SUBSET = 40


# ------------------------------------------------------------------
# One subset of runs -> tracking slope and typical error
# ------------------------------------------------------------------
# A "cut" says how many of the runs count as the best ones. cutfns maps a cut key to a function of the number of usable runs:
#   share cut  (key = 0.2):  best 20% of the runs used (never fewer than the 10-run minimum)
#   count cut  (key = 25):   always the best 25 runs
def share_cuts(fracs):
    return {f: (lambda n, f=f: T.n_top_for(n, f)) for f in fracs}


def count_cuts(counts):
    return {c: (lambda n, c=c: min(int(c), n)) for c in counts}


def medians_for(cc, K, cutfns, sub=None):
    """Median log10 cc of the best runs for every truth and cut. K has shape (truths, runs); sub = run positions or None (all).
    Same rule as analyze_cc_recovery_110.py: finite runs, best KGE first (stable), first m of them."""
    c = cc if sub is None else cc[sub]
    Ks = K if sub is None else K[:, sub]
    out = {key: np.empty(K.shape[0]) for key in cutfns}
    for ti in range(K.shape[0]):
        k = Ks[ti]
        ok = np.flatnonzero(np.isfinite(k))
        order = ok[np.argsort(-k[ok], kind="stable")]
        for key, fn in cutfns.items():
            out[key][ti] = np.median(c[order[:fn(len(ok))]])
    return out


def set_axes(truths, sets):
    """log10 true cc for each tracking set: list of (label, idx, x)."""
    return [(lab, np.array(idx), np.log10(np.array([truths[i] for i in idx], float))) for lab, idx in sets]


def stats_from_medians(med, axes):
    """(slope, mean |log10 median - log10 truth|) for each cut and truth set."""
    out = {}
    for key, y in med.items():
        for lab, idx, x in axes:
            out[(key, lab)] = (float(R.slope_of(x, y[idx])), float(np.mean(np.abs(y[idx] - x))))
    return out


def collect(data, cutfns, axes, subs):
    """subs = {series: run positions or None}. Returns {(series, cut, label): (slope, typical error in log10)}."""
    res = {}
    for s in SERIES:
        med = medians_for(data[s]["X"][:, 2], data[s]["K"], cutfns, subs[s])
        for (key, lab), v in stats_from_medians(med, axes).items():
            res[(s, key, lab)] = v
    return res


# ------------------------------------------------------------------
# The analyses
# ------------------------------------------------------------------
def run_learning(data, cutfns, axes, sizes, reps, seed, n_runs, tag):
    """Random subsets (without replacement) of each size. Returns {(size, series, cut, label): (slopes[reps], typical errors[reps])}."""
    out = {}
    for m in sizes:
        acc = {}
        for si, s in enumerate(SERIES):
            rng = np.random.default_rng([seed, 31 + tag, m, si])
            for r in range(reps):
                sub = rng.choice(n_runs, size=m, replace=False)
                med = medians_for(data[s]["X"][:, 2], data[s]["K"], cutfns, sub)
                for (key, lab), v in stats_from_medians(med, axes).items():
                    acc.setdefault((s, key, lab), []).append(v)
        for key, vals in acc.items():
            a = np.array(vals)
            out[(m,) + key] = (a[:, 0], a[:, 1])
    return out


def run_split_half(data, cutfns, axes, splits, seed, n_runs):
    """Random splits into two disjoint halves. Returns {(series, cut, label): (slopes A[splits], slopes B[splits])} and the half size."""
    h = n_runs // 2
    out = {}
    for si, s in enumerate(SERIES):
        rng = np.random.default_rng([seed, 41, si])
        for r in range(splits):
            perm = rng.permutation(n_runs)
            ra = stats_from_medians(medians_for(data[s]["X"][:, 2], data[s]["K"], cutfns, perm[:h]), axes)
            rb = stats_from_medians(medians_for(data[s]["X"][:, 2], data[s]["K"], cutfns, perm[h:2 * h]), axes)
            for (key, lab), va in ra.items():
                out.setdefault((s, key, lab), ([], []))
                out[(s, key, lab)][0].append(va[0])
                out[(s, key, lab)][1].append(rb[(key, lab)][0])
    return {k: (np.array(v[0]), np.array(v[1])) for k, v in out.items()}, h


def bootstrap_sd(boot, axes):
    """Spread of the tracking slope over the joint bootstrap draws of analyze_cc_recovery_110.py."""
    out = {}
    for (s, f), arr in boot.items():
        for lab, idx, x in axes:
            out[(s, f, lab)] = float(np.std(R.slope_of(x, arr["med"][:, idx]), ddof=1))
    return out


# ------------------------------------------------------------------
# Tables
# ------------------------------------------------------------------
def build_learning_table(lc, full, sizes, kind, cutfns, axes, n_runs):
    rows = []
    for s in SERIES:
        for key, fn in cutfns.items():
            for lab, idx, x in axes:
                for m in list(sizes) + [n_runs]:
                    nb = int(fn(m))
                    raised = bool(kind == "share" and nb > int(round(key * m)))
                    base = {"series": s, "truth_set": lab, "n_truths": len(idx), "cut_kind": kind, "cut": key,
                            "subset_size": m, "n_best": nb, "share_of_runs": nb / float(m), "cut_raised_to_minimum": raised}
                    if m == n_runs:
                        sl, te = full[(s, key, lab)]
                        base.update({"subsets": 1, "slope_mean": sl, "slope_median": sl, "slope_p10": sl, "slope_p90": sl,
                                     "slope_sd": np.nan, "typ_error_median": float(10 ** te)})
                    else:
                        sl, te = lc[(m, s, key, lab)]
                        base.update({"subsets": len(sl), "slope_mean": float(np.mean(sl)), "slope_median": float(np.median(sl)),
                                     "slope_p10": float(np.percentile(sl, 10)), "slope_p90": float(np.percentile(sl, 90)),
                                     "slope_sd": float(np.std(sl, ddof=1)), "typ_error_median": float(10 ** np.median(te))})
                    rows.append(base)
    return pd.DataFrame(rows)


def build_split_table(sh, h, full, fracs, axes):
    rows = []
    for s in SERIES:
        for f in fracs:
            for lab, idx, x in axes:
                a, b = sh[(s, f, lab)]
                both = np.concatenate([a, b])
                rows.append({"series": s, "truth_set": lab, "n_truths": len(idx), "top_frac": f, "n_half": h,
                             "n_best_half": T.n_top_for(h, f), "cut_raised_to_minimum": T.n_top_for(h, f) > int(round(f * h)),
                             "splits": len(a), "slope_all_runs": full[(s, f, lab)][0],
                             "half_slope_mean": float(np.mean(both)), "half_slope_sd": float(np.std(both, ddof=1)),
                             "diff_sd": float(np.std(a - b, ddof=1)), "diff_median_abs": float(np.median(np.abs(a - b))),
                             "share_both_positive": float(np.mean((a > 0) & (b > 0))),
                             "share_both_below_zero": float(np.mean((a < 0) & (b < 0)))})
    return pd.DataFrame(rows)


def build_projection(split_tab, bsd, n_runs, h, proj_n):
    rows = []
    for _, r in split_tab.iterrows():
        sd_h = r.diff_sd / np.sqrt(2.0)                       # spread of ONE half's slope (two independent halves differ by sqrt(2) times it)
        for n in sorted(set([h, n_runs] + list(proj_n))):
            sd = sd_h * np.sqrt(h / float(n))
            rows.append({"series": r.series, "truth_set": r.truth_set, "n_truths": r.n_truths, "top_frac": r.top_frac, "n_runs": int(n),
                         "source": "measured on halves" if n == h else "projected from halves",
                         "slope_sd": float(sd), "half_width_95": float(1.96 * sd)})
        b = bsd.get((r.series, r.top_frac, r.truth_set))
        if b is not None:
            rows.append({"series": r.series, "truth_set": r.truth_set, "n_truths": r.n_truths, "top_frac": r.top_frac, "n_runs": int(n_runs),
                         "source": "bootstrap (other scripts)", "slope_sd": float(b), "half_width_95": float(1.96 * b)})
    return pd.DataFrame(rows)


def compare_with_saved(full, summary_dir):
    p = summary_dir / R.OUT_DIRNAME / "recovery_tracking_110.csv"
    if not p.exists():
        return "no saved recovery_tracking_110.csv found to compare with (run analyze_cc_recovery_110.py first if you want this check)."
    old = pd.read_csv(p)
    rows = []
    for (s, f, lab), (sl, te) in full.items():
        g = old[(old.series == s) & np.isclose(old.top_frac, f) & (old.truth_set == lab)]
        if len(g):
            rows.append(abs(sl - float(g.iloc[0]["slope"])))
    if not rows:
        return "the saved recovery_tracking_110.csv has no matching rows (different cuts or truths), so no comparison was made."
    d = max(rows)
    if d < 1e-9:
        return "the slopes on all runs reproduce analyze_cc_recovery_110.py (largest difference %.1e, %d rows)." % (d, len(rows))
    return ("NOTE: the slopes on all runs differ from the saved recovery_tracking_110.csv by up to %.3g (%d rows). That file may come from "
            "different settings (cuts, --truths) or different inputs." % (d, len(rows)))


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------
def print_split(tab, n_runs, h, splits):
    print("\n" + "=" * 150)
    print("1. SPLIT-HALF: %d random splits of the %d runs into two independent halves of %d, each half analysed on its own (cc ON)"
          % (splits, n_runs, h))
    print("=" * 150)
    print("%-17s %-5s %-6s | %-12s | %-24s | %-22s | %-11s | %-26s"
          % ("truths", "cut", "best", "slope on all", "slope in a half", "halves differ by", "both halves", "control: slope in a half"))
    print("%-17s %-5s %-6s | %-12s | %-24s | %-22s | %-11s | %-26s"
          % ("", "", "(half)", "%d runs" % n_runs, "mean   (spread)", "(median, |A - B|)", "> 0", "mean   (spread)"))
    for lab in tab.truth_set.unique():
        for f in sorted(tab.top_frac.unique(), reverse=True):
            a = tab[(tab.series == SER_ON) & (tab.truth_set == lab) & np.isclose(tab.top_frac, f)].iloc[0]
            c = tab[(tab.series == SER_OFF) & (tab.truth_set == lab) & np.isclose(tab.top_frac, f)].iloc[0]
            best = "%d%s" % (a.n_best_half, "*" if a.cut_raised_to_minimum else "")
            print("%-17s %-5s %-6s | %+8.2f     | %+6.2f   (%5.2f)        | %-22s | %-11s | %+6.2f   (%5.2f)"
                  % (lab if np.isclose(f, max(tab.top_frac)) else "", "%d%%" % round(100 * f), best, a.slope_all_runs, a.half_slope_mean,
                     a.half_slope_sd, "%.2f" % a.diff_median_abs, "%3.0f%%" % (100 * a.share_both_positive), c.half_slope_mean, c.half_slope_sd))
        print()
    print("  slope in a half: the tracking slope computed from only that half of the runs.  spread: how much it wanders between random halves.")
    print("  halves differ by: the typical gap between the slopes of the two halves of one split (median over splits).  both halves > 0: share of")
    print("  splits in which both halves show tracking.  * = the cut was raised to the 10-run minimum for a half.")
    print("  Both halves come from the same 250 runs, so any chance pattern in those runs shows up in both: the halves can agree more than two")
    print("  truly new samples would. That is why the learning curves and the projection below matter too.")


def print_learning(tab, kind, n_runs):
    d0 = tab[tab.cut_kind == kind]
    cuts = sorted(d0.cut.unique(), reverse=True)
    print("\n" + "=" * 150)
    if kind == "share":
        print("2. LEARNING CURVE, SAME CUT SHARE: does the slope of the best 20%% / 10%% / 5%% depend on how many runs there are?  (random subsets of the %d runs)" % n_runs)
    else:
        print("3. LEARNING CURVE, SAME NUMBER OF BEST RUNS: with more runs to choose from, does the same number of best runs track cc better?  (random subsets of the %d)" % n_runs)
    print("=" * 150)
    for c in cuts:
        if kind == "share":
            print("\nbest %d%% of the runs used" % round(100 * c))
        else:
            print("\nalways the best %d runs" % int(c))
        print("%-17s %9s %6s %7s | %-34s %-13s | %-14s"
              % ("truths", "runs used", "best", "share", "ON slope: mean [10th, 90th pct]", "typ. error", "control: mean"))
        for lab in d0.truth_set.unique():
            d = d0[(d0.truth_set == lab) & np.isclose(d0.cut, c)]
            on = d[d.series == SER_ON].sort_values("subset_size")
            off = d[d.series == SER_OFF].set_index("subset_size")
            for i, (_, r) in enumerate(on.iterrows()):
                full = int(r.subset_size) == n_runs
                best = "%d%s" % (r.n_best, "*" if r.cut_raised_to_minimum else "")
                if full:
                    sl = "%+6.2f   (the actual result)" % r.slope_mean
                else:
                    sl = "%+6.2f [%+5.2f, %+5.2f]" % (r.slope_mean, r.slope_p10, r.slope_p90)
                print("%-17s %9s %6s %6.0f%% | %-34s %-13s | %+6.2f"
                      % (lab if i == 0 else "", "%d%s" % (r.subset_size, " (all)" if full else ""), best, 100 * r.share_of_runs, sl,
                         "x%.2f" % r.typ_error_median, off.loc[r.subset_size, "slope_mean"]))
            print()
    print("  mean [10th, 90th pct]: average slope over the random subsets and the range containing the middle 80% of them.  typ. error: typical factor")
    print("  by which the best runs' median misses the true cc (median over subsets).  share: best runs as a share of the runs used.")
    if kind == "share":
        print("  * = cut raised to the 10-run minimum.")
    else:
        print("  Read down each block: the number of best runs is fixed, so more runs used means a smaller share, i.e. a tighter cut.")
    print("  Subsets of 225 out of 250 runs overlap heavily, so their range is narrower than a genuinely new sample would give.")


def print_projection(proj, n_runs, h):
    print("\n" + "=" * 150)
    print("4. PROJECTION: expected 95% half-width of the tracking slope (slope +/- this) with more runs  (cc ON; spread falls as 1/sqrt(runs); share cuts)")
    print("=" * 150)
    ns = sorted(proj.n_runs.unique())
    head = "%-17s %-5s | " % ("truths", "cut") + " ".join("%9s" % ("%d" % n) for n in ns) + " | %-22s" % ("bootstrap at %d" % n_runs)
    print(head)
    print("%-17s %-5s | " % ("", "") + " ".join("%9s" % ("(halves)" if n == h else "") for n in ns) + " | %-22s" % "(other scripts)")
    for lab in proj.truth_set.unique():
        for f in sorted(proj.top_frac.unique(), reverse=True):
            d = proj[(proj.series == SER_ON) & (proj.truth_set == lab) & np.isclose(proj.top_frac, f)]
            hw = d[d.source != "bootstrap (other scripts)"].set_index("n_runs")["half_width_95"]
            b = d[d.source == "bootstrap (other scripts)"]
            btxt = ""
            if len(b):
                bw = float(b.iloc[0].half_width_95)
                btxt = "%.2f   (halves give %.2f)" % (bw, float(hw.loc[n_runs]))
            print("%-17s %-5s | " % (lab if np.isclose(f, max(proj.top_frac)) else "", "%d%%" % round(100 * f))
                  + " ".join("%9.2f" % float(hw.loc[n]) for n in ns) + " | %-22s" % btxt)
        print()
    print("  Columns are the number of runs. The column marked (halves) is measured; the others scale it by sqrt(runs).  If the bootstrap value is close")
    print("  to the 'halves give' value, the bootstrap intervals in the other scripts are about the right size.")
    print("  Both treat the runs as a plain random sample. Your runs are a Latin hypercube, which covers the space more evenly, so in my simulations of such")
    print("  designs both OVERSTATED the true spread, by up to about 1.4 times: read these widths as slightly conservative, not as too optimistic.")
    print("  This projects the INTERVAL WIDTH at a fixed cut share. It does not say how the slope itself would change: that is what table 3 is for.")


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def fig_learning(tab, kind, n_runs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    d0 = tab[tab.cut_kind == kind]
    labels = list(d0.truth_set.unique())
    cuts = sorted(d0.cut.unique(), reverse=True)
    nr, nc = len(labels), len(cuts)
    if not nr:
        return False
    fig, axes = plt.subplots(nr, nc, figsize=(4.4 * nc + 0.5, 3.5 * nr + 1.2), facecolor=T.SURFACE, sharex=True, sharey=True, squeeze=False)
    ymin = min(-0.3, float(np.nanmin(d0.slope_p10)) - 0.1)
    ymax = max(1.1, float(np.nanmax(d0.slope_p90)) + 0.1)
    sizes = sorted(d0.subset_size.unique())
    for ri, lab in enumerate(labels):
        for ci, c in enumerate(cuts):
            ax = axes[ri][ci]
            T._style(ax)
            ax.axhline(1.0, color=T.INK, linestyle="--", linewidth=1.0, zorder=2)
            ax.axhline(0.0, color=T.COL_NULL, linestyle=":", linewidth=1.4, zorder=2)
            d = d0[(d0.truth_set == lab) & np.isclose(d0.cut, c)]
            off = d[d.series == SER_OFF].sort_values("subset_size")
            on = d[d.series == SER_ON].sort_values("subset_size")
            ax.plot(off.subset_size, off.slope_mean, color=T.COL_NULL, linewidth=1.3, zorder=3)
            ax.plot(off.subset_size, off.slope_mean, linestyle="none", marker="s", markersize=5, markerfacecolor=T.SURFACE,
                    markeredgecolor=T.COL_NULL, markeredgewidth=1.4, zorder=4)
            sub = on[on.subset_size < n_runs]
            ax.fill_between(sub.subset_size, sub.slope_p10, sub.slope_p90, color=T.COL_CC, alpha=0.18, linewidth=0, zorder=2)
            ax.plot(on.subset_size, on.slope_mean, color=T.COL_CC, linewidth=2.0, zorder=5)
            ax.plot(sub.subset_size, sub.slope_mean, linestyle="none", marker="o", markersize=5.5, markerfacecolor=T.COL_CC,
                    markeredgecolor=T.COL_CC, zorder=6)
            full = on[on.subset_size == n_runs]
            ax.plot(full.subset_size, full.slope_mean, linestyle="none", marker="D", markersize=7, markerfacecolor=T.COL_CC,
                    markeredgecolor=T.INK, markeredgewidth=1.2, zorder=7)
            raised = on[on.cut_raised_to_minimum & (on.subset_size < n_runs)]
            if len(raised):
                ax.plot(raised.subset_size, raised.slope_mean, linestyle="none", marker="o", markersize=9, markerfacecolor="none",
                        markeredgecolor=T.INK2, markeredgewidth=0.9, zorder=6)
            ax.set_xticks(sizes)
            ax.set_xticklabels(["%d" % m for m in sizes], fontsize=8.5)
            ax.set_ylim(ymin, ymax)
            nt = int(d.n_truths.iloc[0])
            what = "Best %d%% of runs used" % round(100 * c) if kind == "share" else "Always the best %d runs" % int(c)
            ax.set_title("%s, %s (%d)" % (what, lab, nt), color=T.INK, fontsize=9, loc="left")
            if ri == nr - 1:
                ax.set_xlabel("runs used (random subsets of the %d)" % n_runs, color=T.INK2, fontsize=9)
            if ci == 0:
                ax.set_ylabel("tracking slope", color=T.INK2, fontsize=9)
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5.5, linewidth=2.0, label="channel loss ON: mean over random subsets, band = middle 80%"),
        Line2D([0], [0], color=T.COL_CC, marker="D", markersize=7, markeredgecolor=T.INK, linestyle="none", label="the actual result on all runs"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=5, markerfacecolor=T.SURFACE, markeredgewidth=1.4, linewidth=1.3,
               label="control (cc inert): mean over random subsets"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=1.0, label="1 = follows the true cc one-for-one"),
        Line2D([0], [0], color=T.COL_NULL, linestyle=":", linewidth=1.4, label="0 = ignores cc"),
    ]
    if kind == "share":
        handles.insert(3, Line2D([0], [0], color="none", marker="o", markersize=9, markerfacecolor="none", markeredgecolor=T.INK2,
                                 markeredgewidth=0.9, label="ring = cut raised to the 10-run minimum"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8, labelcolor=T.INK2)
    title = ("Does the cc tracking slope depend on the number of runs?  (same cut share)" if kind == "share"
             else "Would more runs help? The same number of best runs, chosen from more runs  (a tighter cut as runs increase)")
    fig.suptitle(title, color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.10 if nr == 1 else 0.09, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


def fig_projection(proj, fracs, n_runs, h, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    d0 = proj[proj.series == SER_ON]
    labels = list(d0.truth_set.unique())
    if not labels:
        return False
    fig, axes = plt.subplots(1, len(labels), figsize=(5.4 * len(labels) + 0.5, 4.8), facecolor=T.SURFACE, sharey=True, squeeze=False)
    styles = {}
    for f, ls, mk in zip(sorted(fracs, reverse=True), ("-", "--", ":"), ("o", "s", "^")):
        styles[f] = (ls, mk)
    ns = sorted(d0.n_runs.unique())
    ymax = 1.08 * float(d0.half_width_95.max())
    for ci, lab in enumerate(labels):
        ax = axes[0][ci]
        T._style(ax)
        for f in sorted(fracs, reverse=True):
            ls, mk = styles.get(f, ("-", "o"))
            d = d0[(d0.truth_set == lab) & np.isclose(d0.top_frac, f)]
            line = d[d.source != "bootstrap (other scripts)"].sort_values("n_runs")
            ax.plot(line.n_runs, line.half_width_95, color=T.COL_CC, linestyle=ls, linewidth=2.0 if ls == "-" else 1.8, zorder=3)
            m = line[line.n_runs == h]
            ax.plot(m.n_runs, m.half_width_95, linestyle="none", marker=mk, markersize=6.5, markerfacecolor=T.COL_CC,
                    markeredgecolor=T.COL_CC, zorder=5)
            b = d[d.source == "bootstrap (other scripts)"]
            ax.plot(b.n_runs, b.half_width_95, linestyle="none", marker="D", markersize=7, markerfacecolor=T.SURFACE,
                    markeredgecolor=T.INK, markeredgewidth=1.4, zorder=6)
        ax.set_xscale("log")
        ax.set_xticks(ns)
        ax.set_xticklabels(["%d" % n for n in ns], fontsize=8.5)
        ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlim(min(ns) * 0.85, max(ns) * 1.2)
        ax.set_ylim(0, ymax)
        ax.set_title("%s (%d)" % (lab, int(d0[d0.truth_set == lab].n_truths.iloc[0])), color=T.INK, fontsize=9, loc="left")
        ax.set_xlabel("number of runs", color=T.INK2, fontsize=9)
        if ci == 0:
            ax.set_ylabel("95% half-width of the tracking slope", color=T.INK2, fontsize=9)
    handles = [Line2D([0], [0], color=T.COL_CC, linestyle=styles[f][0], linewidth=1.9, label="best %d%% of runs used" % round(100 * f))
               for f in sorted(fracs, reverse=True)]
    handles += [
        Line2D([0], [0], color=T.COL_CC, marker="o", linestyle="none", markersize=6, label="measured on halves (%d runs)" % h),
        Line2D([0], [0], color=T.INK, marker="D", markerfacecolor=T.SURFACE, markeredgewidth=1.4, linestyle="none", markersize=7,
               label="bootstrap on all %d runs (other scripts)" % n_runs),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("How wide would the slope's 95% interval be with more runs?  (projected as 1/sqrt(runs) from the two halves)",
                 color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.13, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110 -- is the cc tracking result stable across subsets of the runs?")
    ap.add_argument("--top_fracs", type=float, nargs="+", default=[0.20, 0.10, 0.05],
                    help="Top fractions (default 0.20 0.10 0.05)")
    ap.add_argument("--counts", type=int, nargs="+", default=None,
                    help="Numbers of best runs for the same-number-of-best-runs learning curve (default: the size of each cut on all runs, e.g. 50 25 12)")
    ap.add_argument("--sizes", type=int, nargs="+", default=None,
                    help="Subset sizes for the learning curves (default: about 30%%, 50%%, 70%% and 90%% of the runs)")
    ap.add_argument("--subsets", type=int, default=300, help="Random subsets per size (default 300; minimum 50)")
    ap.add_argument("--splits", type=int, default=1000, help="Random split-half splits (default 1000; minimum 100)")
    ap.add_argument("--boot", type=int, default=1000, help="Bootstrap re-draws for the cross-check (default 1000)")
    ap.add_argument("--project_n", type=int, nargs="+", default=[500, 1000], help="Run counts to project to (default 500 1000)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--expect_n", type=int, default=250)
    ap.add_argument("--truths", type=float, nargs="+", default=None,
                    help="Only these cc-ON truths (the cc-OFF reference is always kept)")
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--continue_on_fail", action="store_true")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()
    if args.subsets < 50 or args.splits < 100 or args.boot < 200:
        sys.exit("Too few repeats: --subsets must be at least 50, --splits at least 100, --boot at least 200.")
    if any(f <= 0 or f >= 1 for f in args.top_fracs):
        sys.exit("--top_fracs must be fractions between 0 and 1 (for example 0.2 0.1 0.05).")
    if any(n <= 0 for n in args.project_n):
        sys.exit("--project_n must be positive run counts.")

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    out_dir = summary_dir / OUT_DIRNAME
    print("\n" + "=" * 78 + "\nSeries 110 -- is the cc tracking result stable, or does it depend on which runs we drew?\n" + "=" * 78)
    print("Reading from: %s" % summary_dir)

    df, notes = P.load_all(summary_dir)
    for nt in notes:
        print("  NOTE: %s" % nt)
    ck = P.run_checks(df, summary_dir, args.expect_n)
    if args.truths:
        keep = [OFF_TRUTH] + [t for t in sorted(df.truth_cc_mmhr.unique())
                               if t != OFF_TRUTH and any(np.isclose(t, v, rtol=1e-3) for v in args.truths)]
        df = df[df.truth_cc_mmhr.isin(keep)]
    truths = sorted(df["truth_cc_mmhr"].unique())
    fracs = sorted(args.top_fracs, reverse=True)
    data = T.build_tensors(df, truths)
    n_runs = data[SER_ON]["X"].shape[0]
    h = n_runs // 2
    sizes = args.sizes if args.sizes else sorted(set(int(round(n_runs * q)) for q in (0.3, 0.5, 0.7, 0.9)))
    if any(m < MIN_SUBSET or m >= n_runs for m in sizes):
        sys.exit("STOPPED: every --sizes value must be at least %d and smaller than the %d runs available (got %s)." % (MIN_SUBSET, n_runs, sizes))
    if h < MIN_SUBSET:
        sys.exit("STOPPED: only %d runs available; the split-half test needs at least %d runs per half." % (n_runs, MIN_SUBSET))
    sizes = sorted(set(sizes))
    counts = args.counts if args.counts else sorted(set(T.n_top_for(n_runs, f) for f in fracs), reverse=True)
    if any(c < P.MIN_TOP or c >= n_runs for c in counts):
        sys.exit("STOPPED: every --counts value must be at least %d and smaller than the %d runs available (got %s)." % (P.MIN_TOP, n_runs, counts))
    counts = sorted(set(counts), reverse=True)
    P.print_checks(ck)
    if ck.n("FAIL") and not args.continue_on_fail:
        sys.exit("\nSTOPPED: an input check failed. Nothing was analysed or written. (--continue_on_fail to override.)")

    sets = R.tracking_sets(truths)
    if not sets:
        sys.exit("STOPPED: the tracking slope needs at least %d cc-ON truths." % R.MIN_TRACK_TRUTHS)
    axes = set_axes(truths, sets)
    cf_share, cf_count = share_cuts(fracs), count_cuts(counts)
    print("\n%d cc-ON truth(s), %d runs per series, cuts %s (or the best %s runs), learning-curve sizes %s x %d subsets, %d split-half splits, %d bootstrap re-draws"
          % (sum(1 for t in truths if t > 0), n_runs, ["%d%%" % round(100 * f) for f in fracs], counts, sizes, args.subsets, args.splits, args.boot))

    print("\nSlopes on all runs ...")
    full_share = collect(data, cf_share, axes, {s: None for s in SERIES})
    full_count = collect(data, cf_count, axes, {s: None for s in SERIES})
    print("Split-half test ...")
    sh, h = run_split_half(data, cf_share, axes, args.splits, args.seed, n_runs)
    print("Learning curve, same cut share ...")
    lc_share = run_learning(data, cf_share, axes, sizes, args.subsets, args.seed, n_runs, 0)
    print("Learning curve, same number of best runs ...")
    lc_count = run_learning(data, cf_count, axes, sizes, args.subsets, args.seed, n_runs, 1)
    print("Bootstrap cross-check ...")
    boot = R.run_boot(data, truths, fracs, args.boot, args.seed)
    bsd = bootstrap_sd(boot, axes)

    split_tab = build_split_table(sh, h, full_share, fracs, axes)
    lc_tab = pd.concat([build_learning_table(lc_share, full_share, sizes, "share", cf_share, axes, n_runs),
                        build_learning_table(lc_count, full_count, sizes, "count", cf_count, axes, n_runs)], ignore_index=True)
    proj = build_projection(split_tab, bsd, n_runs, h, args.project_n)

    print_split(split_tab, n_runs, h, args.splits)
    print_learning(lc_tab, "share", n_runs)
    print_learning(lc_tab, "count", n_runs)
    print_projection(proj, n_runs, h)
    print("\n  " + compare_with_saved(full_share, summary_dir))

    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ck.rows).to_csv(out_dir / "stability_input_checks_110.csv", index=False)
    lc_tab.to_csv(out_dir / "stability_learning_curve_110.csv", index=False)
    split_tab.to_csv(out_dir / "stability_split_half_110.csv", index=False)
    proj.to_csv(out_dir / "stability_projection_110.csv", index=False)
    if not args.no_plots:
        try:
            fig_learning(lc_tab, "share", n_runs, out_dir / "fig_stability_learning_110.png")
            fig_learning(lc_tab, "count", n_runs, out_dir / "fig_stability_learning_count_110.png")
            fig_projection(proj, fracs, n_runs, h, out_dir / "fig_stability_projection_110.png")
        except Exception as e:
            print("\n  (figures skipped: %s -- the CSVs are saved regardless)" % e)

    def md5_of(p):
        return hashlib.md5(Path(p).read_bytes()).hexdigest()
    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "pandas": pd.__version__, "numpy": np.__version__,
            "inputs": {n: md5_of(summary_dir / n) for n in (P.STAGE2_NAME, P.CONTROL_NAME)},
            "settings": {k: (v if not isinstance(v, Path) else str(v)) for k, v in vars(args).items()},
            "sizes": sizes, "counts": counts, "n_runs": n_runs, "half": h, "truths": [float(t) for t in truths], "checks": ck.rows}
    (out_dir / "PROVENANCE_stability_cc_110.json").write_text(json.dumps(prov, indent=2, default=str))

    print("\nSaved to: %s" % out_dir)
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results." % ck.n("FAIL"))
    print("  Caveat: random subsets of ONE %d-run design stand in for new independent designs. One storm, routing pinned at truth, exact "
          "observations: a best case." % n_runs)


if __name__ == "__main__":
    main()
