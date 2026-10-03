"""
analyze_cc_recovery_110.py
==========================
Series 110 -- do the best runs choose the cc that made the truth?

The PCA and test scripts showed that cc is detectably informative but only weakly
constraining. They did not show WHERE the best runs put cc. Low "shrinkage" means the
best runs cluster; it does not say they cluster on the true value. In a synthetic
experiment the true value is known, so we can check. That is what this script does.

For every cc-ON truth (a truth made with a known channel conductivity) it asks:

  1. Where does cc sit among the best runs?   median, and the middle 80% of values
  2. Is that centred on the true cc?          best-run median divided by the true cc
  3. How many best runs are near the truth?   share within a factor N (default 2) of
                                              the true cc, versus the share of ALL runs
                                              that are (that share is what chance gives)
  4. Do the best runs FOLLOW the true cc as the truth changes?   "tracking slope"

Pure post-hoc analysis. Runs no tRIBS. Writes only into its own new folder
(calibration_work/03_comparisons/summary_tables/recovery_cc_110/).

Keep this file in the same folder as analyze_cc_tests_110.py and analyze_cc_pca_110.py
(it reuses their loader, input checks and plot style). Run it from the lab/ directory.

VOCABULARY
-----------
  Median      The middle value: half the best runs have a lower cc, half a higher one.
  Middle 80%  The range holding the central 80% of the best runs' cc values (from the
              10th to the 90th percentile). A "percentile" is the value below which
              that percent of the runs fall.
  Recovery    Getting back the value that made the truth. Perfect recovery = the best
              runs sit on the true cc.
  Factor      cc differences are multiplicative (cc spans 30-1000), so errors are given
              as factors: "x2" means twice the true cc, "/2" means half of it.
  No information   What you would see if the fit ignored cc completely: the best runs'
              cc would be a random draw from the whole sampled range, so their median
              would land near the middle of that range (about 173 mm/hr, the geometric
              middle of 30-1000) WHATEVER the true cc is.
  Tracking slope   Plot the best runs' median cc against the true cc on log scales and
              fit a line across all the truths. Slope 1 = the median follows the truth
              one-for-one; slope 0 = it does not move (no information). Edges of the
              sampled range squash the slope below 1 even for a perfect method, so
              compare it with 0 and with the control, not only with 1.
  Permutation test / bootstrap   As in analyze_cc_tests_110.py. The permutation test
              shuffles cc among the runs to show what "no information" looks like; the
              joint bootstrap re-draws the runs (the same re-draw for every truth) to
              show how much a number wobbles.

HOW TO READ THE RESULTS
------------------------
  * Median centred on the truth, bar (middle 80%) much narrower than the grey
    no-information band, slope near 1, control slope near 0: cc is recovered.
  * Median stuck near 173 whatever the truth, bar as wide as the grey band, slope near 0:
    the best runs ignore cc; any shrinkage seen earlier was not tracking the truth.
  * Median moves the right way but not far enough (slope between 0 and 1): cc is only
    partly recovered; the data pull the estimate toward the truth but weakly.
  * Truths near 173 cannot be told from "no information" by their own row, because runs
    that ignore cc would land there too. Use the slope, which uses all the truths.
  * The control series (cc inert) must show a slope near 0. If it does not, something
    other than cc is driving the pattern.

OUTPUTS   (calibration_work/03_comparisons/summary_tables/recovery_cc_110/)
----------------------------------------------------------------------------
  recovery_per_truth_110.csv   per series x top_frac x truth: median, middle 80%, factor,
                               share near truth, bootstrap interval for the median, p-values
  recovery_tracking_110.csv    tracking slope and typical error across truths, with
                               bootstrap interval and permutation p-value
  recovery_input_checks_110.csv
  fig_recovery_cc_110.png      best-run cc against true cc, one panel per cut
  PROVENANCE_recovery_cc_110.json

USAGE (run from the lab/ directory)
------------------------------------
    python analyze_cc_recovery_110.py                  # about a minute
    python analyze_cc_recovery_110.py --within 3       # "near" = within a factor of 3
    python analyze_cc_recovery_110.py --top_fracs 0.2 0.1 0.05 --boot 2000 --perm 10000
    python analyze_cc_recovery_110.py --truths 157 425 --no_plots

CAVEATS
--------
One storm, 250 runs, routing pinned at truth: a best case. With only 12-50 runs in a
top set, a median is noisy; the bootstrap interval says how noisy. The sampled cc range
(30-1000) is finite, so truths near its edges are pulled toward the middle by the box
itself. A truth is flagged "edge" when it lies within 25% of the log range of an edge.
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
    import analyze_cc_tests_110 as T
except ImportError:
    sys.exit("This script needs analyze_cc_tests_110.py and analyze_cc_pca_110.py in the same folder "
             "(it reuses their loader and input checks). Put all three files in lab/ and run from there.")

_NEEDED = ["P", "build_tensors", "n_top_for", "make_perm_matrix", "p_low", "p_high", "fmt_p", "_style", "box_pos",
           "COL_CC", "COL_NULL", "INK", "INK2", "GRID", "SURFACE", "PARAM_BOX", "SER_ON", "SER_OFF", "OFF_TRUTH"]
_missing = [n for n in _NEEDED if not hasattr(T, n)]
if _missing:
    sys.exit("analyze_cc_tests_110.py in this folder is an older version (missing: %s). Replace it with the current copy."
             % ", ".join(_missing))

P = T.P
SER_ON, SER_OFF, OFF_TRUTH = T.SER_ON, T.SER_OFF, T.OFF_TRUTH
OUT_DIRNAME = "recovery_cc_110"
LOG_LO, LOG_HI = np.log10(T.PARAM_BOX[2][0]), np.log10(T.PARAM_BOX[2][1])
CENTER_LOG = 0.5 * (LOG_LO + LOG_HI)
MIN_TRACK_TRUTHS = 3


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def fcc(v):
    """cc in mm/hr for tables."""
    if not np.isfinite(v):
        return "n/a"
    return "%.0f" % v if v >= 100 else "%.1f" % v


def ffactor(f):
    """x1.20 = 20% above the truth, /1.50 = a third below it."""
    if not np.isfinite(f):
        return "n/a"
    return "x%.2f" % f if f >= 1 else "/%.2f" % (1.0 / f)


def top_index(k, frac):
    """Run indices of the best top_frac of runs by KGE_2012 (same rule as the other scripts)."""
    ok = np.flatnonzero(np.isfinite(k))
    m = T.n_top_for(len(ok), frac)
    return ok[np.argsort(-k[ok], kind="stable")][:m]


def set_stats(v, tl, half_width):
    """Statistics of the log10 cc values of a top set. v has shape (..., m).
    tl = log10 of the true cc (nan for the cc-OFF reference truth, which has no cc)."""
    med = np.median(v, axis=-1)
    q = np.percentile(v, [10, 90], axis=-1)
    if np.isfinite(tl):
        within = np.mean(np.abs(v - tl) <= half_width + 1e-12, axis=-1)
    else:
        within = np.full(np.shape(med), np.nan)
    return med, q[0], q[1], within


def slope_of(x, Y):
    """Least-squares slope of Y (..., k) on x (k,)."""
    xc = x - x.mean()
    return (Y * xc).sum(axis=-1) / (xc ** 2).sum()


def truth_log(t):
    return float(np.log10(t)) if t > 0 else np.nan


def tracking_sets(truths):
    Ton = [i for i, t in enumerate(truths) if t > 0]
    interior = [i for i in Ton if P.NEAR_EDGE <= T.box_pos(truths[i]) <= 1 - P.NEAR_EDGE]
    sets = []
    if len(Ton) >= MIN_TRACK_TRUTHS:
        sets.append(("all cc-ON truths", Ton))
        if len(interior) >= MIN_TRACK_TRUTHS and len(interior) != len(Ton):
            sets.append(("interior truths", interior))
    return sets


# ------------------------------------------------------------------
# Permutation null: what "the fit ignores cc" looks like
# ------------------------------------------------------------------
def run_perm(data, truths, fracs, Bp, seed, half_width):
    """null[(series, frac)] = dict of (Bp, T) arrays: median, q10, q90, share near truth when cc is shuffled
    among all runs and the real top set is kept. ONE shuffle is shared by all truths."""
    null = {}
    for si, s in enumerate((SER_ON, SER_OFF)):
        X, K = data[s]["X"], data[s]["K"]
        cc = X[:, 2]
        Pm = T.make_perm_matrix(len(cc), Bp, np.random.default_rng([seed, 66, si]))
        for frac in fracs:
            arr = {k: np.full((Bp, len(truths)), np.nan) for k in ("med", "q10", "q90", "within")}
            for ti, t in enumerate(truths):
                idx = top_index(K[ti], frac)
                med, q10, q90, within = set_stats(cc[Pm[:, idx]], truth_log(t), half_width)
                arr["med"][:, ti], arr["q10"][:, ti], arr["q90"][:, ti], arr["within"][:, ti] = med, q10, q90, within
            null[(s, frac)] = arr
    return null


# ------------------------------------------------------------------
# Joint bootstrap: how much do the medians wobble?
# ------------------------------------------------------------------
def run_boot(data, truths, fracs, B, seed):
    """boot[(series, frac)] = dict of (B, T) arrays. The same resampled run indices are used for every truth and
    every cut within a series; ON and control are drawn independently."""
    boot = {}
    for si, s in enumerate((SER_ON, SER_OFF)):
        X, K = data[s]["X"], data[s]["K"]
        cc = X[:, 2]
        n = len(cc)
        idxb = np.random.default_rng([seed, 77, si]).integers(0, n, size=(B, n))
        ccb = cc[idxb]
        for frac in fracs:
            arr = {k: np.full((B, len(truths)), np.nan) for k in ("med", "q10", "q90")}
            for ti in range(len(truths)):
                kb = K[ti][idxb]
                kb = np.where(np.isfinite(kb), kb, -np.inf)
                m = T.n_top_for(int(np.isfinite(K[ti]).sum()), frac)
                order = np.argsort(-kb, axis=1, kind="stable")[:, :m]
                v = np.take_along_axis(ccb, order, axis=1)
                med, q10, q90, _ = set_stats(v, np.nan, 0.0)
                arr["med"][:, ti], arr["q10"][:, ti], arr["q90"][:, ti] = med, q10, q90
            boot[(s, frac)] = arr
    return boot


# ------------------------------------------------------------------
# Tables
# ------------------------------------------------------------------
def build_per_truth(data, truths, fracs, half_width, null, boot):
    rows, point_med = [], {}
    for s in (SER_ON, SER_OFF):
        X, K = data[s]["X"], data[s]["K"]
        cc = X[:, 2]
        for frac in fracs:
            pm = np.full(len(truths), np.nan)
            for ti, t in enumerate(truths):
                tl = truth_log(t)
                idx = top_index(K[ti], frac)
                med, q10, q90, within = (float(a) for a in set_stats(cc[idx], tl, half_width))
                pm[ti] = med
                share_all = float(np.mean(np.abs(cc - tl) <= half_width + 1e-12)) if np.isfinite(tl) else np.nan
                nl, bt = null[(s, frac)], boot[(s, frac)]
                row = {
                    "truth_cc_mmhr": t, "series": s, "top_frac": frac, "n_top": int(len(idx)),
                    "median_cc": 10 ** med, "q10_cc": 10 ** q10, "q90_cc": 10 ** q90,
                    "median_cc_ci_lo": float(10 ** np.percentile(bt["med"][:, ti], 2.5)),
                    "median_cc_ci_hi": float(10 ** np.percentile(bt["med"][:, ti], 97.5)),
                    "null_median_cc": float(10 ** np.mean(nl["med"][:, ti])),
                    "median_Ks": float(np.median(X[idx, 0])), "median_f": float(10 ** np.median(X[idx, 1])),
                    "box_pos": T.box_pos(t) if t > 0 else np.nan,
                }
                row["near_box_edge"] = bool(t > 0 and not (P.NEAR_EDGE <= row["box_pos"] <= 1 - P.NEAR_EDGE))
                if np.isfinite(tl):
                    row.update({
                        "log10_bias": med - tl, "factor": 10 ** (med - tl),
                        "truth_in_middle80": float(q10 <= tl <= q90),
                        "share_near": within, "share_near_all_runs": share_all,
                        "enrichment": within / share_all if share_all > 0 else np.nan,
                        "p_near": T.p_high(nl["within"][:, ti], within),
                        "p_median_closer": T.p_low(np.abs(nl["med"][:, ti] - tl), abs(med - tl)),
                    })
                else:
                    row.update({"log10_bias": np.nan, "factor": np.nan, "truth_in_middle80": np.nan, "share_near": np.nan,
                                "share_near_all_runs": np.nan, "enrichment": np.nan, "p_near": np.nan, "p_median_closer": np.nan})
                rows.append(row)
            point_med[(s, frac)] = pm
    return pd.DataFrame(rows), point_med


def build_tracking(truths, fracs, half_width, null, boot, point_med):
    rows = []
    for s in (SER_ON, SER_OFF):
        for frac in fracs:
            for label, idx in tracking_sets(truths):
                x = np.log10(np.array([truths[i] for i in idx]))
                obs = point_med[(s, frac)][idx]
                sl, err = float(slope_of(x, obs)), float(np.mean(np.abs(obs - x)))
                nl = null[(s, frac)]["med"][:, idx]
                sl_null, err_null = slope_of(x, nl), np.mean(np.abs(nl - x), axis=1)
                bt = boot[(s, frac)]["med"][:, idx]
                sl_b, err_b = slope_of(x, bt), np.mean(np.abs(bt - x), axis=1)
                rows.append({
                    "series": s, "top_frac": frac, "truth_set": label, "n_truths": len(idx),
                    "slope": sl, "slope_ci_lo": float(np.percentile(sl_b, 2.5)), "slope_ci_hi": float(np.percentile(sl_b, 97.5)),
                    "slope_null_mean": float(sl_null.mean()), "p_slope": T.p_high(sl_null, sl),
                    "typical_error_factor": float(10 ** err),
                    "typical_error_ci_lo": float(10 ** np.percentile(err_b, 2.5)),
                    "typical_error_ci_hi": float(10 ** np.percentile(err_b, 97.5)),
                    "null_typical_error_factor": float(10 ** err_null.mean()),
                    "p_error_smaller": T.p_low(err_null, err),
                })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------
def print_per_truth(pt, main_frac, factor):
    on = pt[(pt.series == SER_ON) & (pt.top_frac == main_frac)].sort_values("truth_cc_mmhr")
    off = pt[(pt.series == SER_OFF) & (pt.top_frac == main_frac)].set_index("truth_cc_mmhr")
    n = int(on.n_top.iloc[0])
    lo10, hi10 = 10 ** (LOG_LO + 0.1 * (LOG_HI - LOG_LO)), 10 ** (LOG_LO + 0.9 * (LOG_HI - LOG_LO))
    print("\n" + "=" * 126)
    print("1. WHERE DO THE BEST RUNS PUT cc?  (cc-ON runs, best %d%% = %d runs by KGE_2012; 'near' = within a factor of %g)"
          % (round(100 * main_frac), n, factor))
    print("=" * 126)
    print("  No information (runs that ignore cc) would give a median near %s mm/hr and a middle 80%% of about %s-%s, "
          "whatever the true cc." % (fcc(10 ** CENTER_LOG), fcc(lo10), fcc(hi10)))
    print("%12s | %7s [%7s - %7s] | %8s | %5s | %-12s %7s | %7s | %7s | %6s"
          % ("true cc", "median", "10%", "90%", "vs truth", "in80%", "near/chance", "p", "err p", "ctl med", "Ks med"))
    for _, a in on.iterrows():
        t = a.truth_cc_mmhr
        c = off.loc[t] if t in off.index else None
        ctl = fcc(c.median_cc) if c is not None else "n/a"
        if t > 0:
            print("%12s | %7s [%7s - %7s] | %8s | %5s | %-12s %7s | %7s | %7s | %6.1f%s"
                  % (fcc(t), fcc(a.median_cc), fcc(a.q10_cc), fcc(a.q90_cc), ffactor(a.factor),
                     "yes" if a.truth_in_middle80 == 1 else "NO",
                     "%d%% / %d%%" % (round(100 * a.share_near), round(100 * a.share_near_all_runs)),
                     T.fmt_p(a.p_near), T.fmt_p(a.p_median_closer), ctl, a.median_Ks, "  edge" if a.near_box_edge else ""))
        else:
            print("%12s | %7s [%7s - %7s] | %8s | %5s | %-12s %7s | %7s | %7s | %6.1f"
                  % ("cc-OFF truth", fcc(a.median_cc), fcc(a.q10_cc), fcc(a.q90_cc), "no cc", "", "", "", "", ctl, a.median_Ks))
    print("  vs truth: median of the best runs / true cc (x1.00 = centred on the truth).  in80%: is the true cc inside the middle 80%?")
    print("  near/chance: share of the best runs within a factor of %g of the true cc / the share of ALL runs that are (what chance gives)." % factor)
    print("  p: could that much enrichment arise if cc were shuffled?   err p: is the median closer to the truth than a shuffled one?")
    print("  ctl med: median cc of the best runs in the control series (cc inert); it should wander around %s whatever the truth."
          % fcc(10 ** CENTER_LOG))
    print("  Ks med: median Ks multiplier of the best runs (the true Ks multiplier is 7.0).   'edge': truth within 25% of an edge "
          "of the sampled range, where the box squeezes the answer toward the middle.")
    print("  Rows near %s mm/hr cannot be told from 'no information' on their own; read section 2." % fcc(10 ** CENTER_LOG))


def print_tracking(tr, fracs):
    print("\n" + "=" * 126)
    print("2. DO THE BEST RUNS FOLLOW THE TRUE cc?  slope of (median best-run cc) against (true cc), both on log scales")
    print("=" * 126)
    print("  1 = follows the truth one-for-one;  0 = ignores it.  The control (cc inert) is the reference for 0.")
    print("  %-5s %-8s %-17s %-24s %9s %7s | %-9s %9s %7s"
          % ("cut", "series", "truths", "slope [95% interval]", "no-info", "p", "typ. error", "no-info", "p"))
    for frac in fracs:
        for label in tr.truth_set.unique():
            for s in (SER_ON, SER_OFF):
                g = tr[(tr.top_frac == frac) & (tr.truth_set == label) & (tr.series == s)]
                if g.empty:
                    continue
                r = g.iloc[0]
                print("  %-5s %-8s %-17s %+6.2f [%+5.2f, %+5.2f]       %+8.2f %7s | x%-8.2f x%-8.2f %7s"
                      % ("%d%%" % round(100 * frac), s.replace("cc ", ""), label, r.slope, r.slope_ci_lo, r.slope_ci_hi,
                         r.slope_null_mean, T.fmt_p(r.p_slope), r.typical_error_factor, r.null_typical_error_factor,
                         T.fmt_p(r.p_error_smaller)))
        print()
    print("  typ. error: typical factor by which the best runs' median misses the true cc (x1.00 = perfect), versus what a "
          "shuffled (no-information) data set gives.")
    print("  Slopes are squashed below 1 by the edges of the sampled range even for a perfect method; the 'interior truths' rows "
          "reduce that.")
    print("  The control should show a slope near 0 and a large p. If it does not, something other than cc drives the pattern.")


# ------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------
def fig_recovery(pt, fracs, factor, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    d = pt[pt.truth_cc_mmhr > 0]
    if d.empty:
        return False
    tv = sorted(d.truth_cc_mmhr.unique())
    nf = len(fracs)
    fig, axes = plt.subplots(1, nf, figsize=(4.7 * nf + 0.5, 5.6), facecolor=T.SURFACE, sharey=True, squeeze=False)
    lo, hi = 10 ** LOG_LO, 10 ** LOG_HI
    lim = (lo * 0.75, hi * 1.3)
    for ci, frac in enumerate(fracs):
        ax = axes[0][ci]
        T._style(ax)
        ax.axhspan(10 ** (LOG_LO + 0.1 * (LOG_HI - LOG_LO)), 10 ** (LOG_LO + 0.9 * (LOG_HI - LOG_LO)),
                   color=T.COL_NULL, alpha=0.16, linewidth=0, zorder=1)
        ax.axhline(10 ** CENTER_LOG, color=T.COL_NULL, linestyle=":", linewidth=1.3, zorder=2)
        ax.plot(lim, lim, color=T.INK, linestyle="--", linewidth=1.0, zorder=2)
        for s, col, off, mk, filled in ((SER_OFF, T.COL_NULL, +1, "s", False), (SER_ON, T.COL_CC, -1, "o", True)):
            g = d[(d.series == s) & (d.top_frac == frac)].sort_values("truth_cc_mmhr")
            x = g.truth_cc_mmhr.values * 10 ** (off * 0.012)
            ax.vlines(x, g.q10_cc, g.q90_cc, color=col, linewidth=2.0 if s == SER_ON else 1.4, zorder=3)
            ax.plot(x, g.median_cc.values, linestyle="none", marker=mk, markersize=5.5,
                    markerfacecolor=col if filled else T.SURFACE, markeredgecolor=col, markeredgewidth=1.5, zorder=4)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xticks(tv)
        ax.set_xticklabels(["%g" % v for v in tv], rotation=60, fontsize=7.5)
        ax.set_yticks([30, 100, 300, 1000])
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
        ax.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
        n = int(d[(d.series == SER_ON) & (d.top_frac == frac)].n_top.iloc[0])
        ax.set_title("Best %d%% of runs (n = %d)" % (round(100 * frac), n), color=T.INK, fontsize=9.5, loc="left")
        ax.set_xlabel("true cc of the truth (mm/hr)", color=T.INK2, fontsize=9)
        if ci == 0:
            ax.set_ylabel("cc of the best runs (mm/hr)", color=T.INK2, fontsize=9)
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5.5, linewidth=2.0, label="channel loss ON: median and middle 80% of the best runs"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=5.5, markerfacecolor=T.SURFACE, markeredgewidth=1.5, linewidth=1.4,
               label="control (cc inert): the same, for runs that cannot feel cc"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=1.0, label="perfect recovery (best runs sit on the true cc)"),
        Patch(facecolor=T.COL_NULL, alpha=0.16, label="no information: middle 80% of the sampled range (dotted line = its middle)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Do the best runs choose the cc that made the truth?  (dot = median, bar = middle 80%)",
                 color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.11, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110 -- where do the best runs put cc, compared with the true cc?")
    ap.add_argument("--top_fracs", type=float, nargs="+", default=[0.20, 0.10, 0.05],
                    help="Top fractions by KGE_2012; the first is the main one (default 0.20 0.10 0.05)")
    ap.add_argument("--within", type=float, default=2.0,
                    help="'Near the truth' means within this FACTOR of the true cc (default 2 = between half and double)")
    ap.add_argument("--boot", type=int, default=1000, help="Joint-bootstrap re-draws (default 1000)")
    ap.add_argument("--perm", type=int, default=5000, help="Permutation shuffles (default 5000)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--expect_n", type=int, default=250)
    ap.add_argument("--truths", type=float, nargs="+", default=None,
                    help="Only these cc-ON truths (the cc-OFF reference is always kept)")
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--continue_on_fail", action="store_true")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()
    if args.within <= 1.0:
        sys.exit("--within must be greater than 1 (it is a factor, e.g. 2 = between half and double).")
    half_width = float(np.log10(args.within))

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    out_dir = summary_dir / OUT_DIRNAME
    print("\n" + "=" * 78 + "\nSeries 110 -- where do the best runs put cc, compared with the true cc?\n" + "=" * 78)
    print("Reading from: %s" % summary_dir)

    df, notes = P.load_all(summary_dir)
    for nt in notes:
        print("  NOTE: %s" % nt)
    ck = P.run_checks(df, summary_dir, args.expect_n)
    P.print_checks(ck)
    if ck.n("FAIL") and not args.continue_on_fail:
        sys.exit("\nSTOPPED: an input check failed. Nothing was analysed or written. (--continue_on_fail to override.)")

    if args.truths:
        keep = [OFF_TRUTH] + [t for t in sorted(df.truth_cc_mmhr.unique())
                               if t != OFF_TRUTH and any(np.isclose(t, v, rtol=1e-3) for v in args.truths)]
        df = df[df.truth_cc_mmhr.isin(keep)]
    truths = sorted(df["truth_cc_mmhr"].unique())
    fracs = args.top_fracs
    main_frac = fracs[0]
    data = T.build_tensors(df, truths)
    Ton = [i for i, t in enumerate(truths) if t > 0]
    print("\n%d truth(s) (%d cc-ON + %s), %d runs per series, cuts %s, %d shuffles, %d bootstrap re-draws, near = within x%g"
          % (len(truths), len(Ton), "cc-OFF reference" if OFF_TRUTH in truths else "no reference",
             data[SER_ON]["X"].shape[0], ["%d%%" % round(100 * f) for f in fracs], args.perm, args.boot, args.within))

    print("\nRunning the shuffles and bootstrap ...")
    null = run_perm(data, truths, fracs, args.perm, args.seed, half_width)
    boot = run_boot(data, truths, fracs, args.boot, args.seed)
    pt, point_med = build_per_truth(data, truths, fracs, half_width, null, boot)
    tr = build_tracking(truths, fracs, half_width, null, boot, point_med)

    print_per_truth(pt, main_frac, args.within)
    if len(tr):
        print_tracking(tr, fracs)
    else:
        print("\n  (tracking slope skipped: it needs at least %d cc-ON truths)" % MIN_TRACK_TRUTHS)

    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ck.rows).to_csv(out_dir / "recovery_input_checks_110.csv", index=False)
    pt.to_csv(out_dir / "recovery_per_truth_110.csv", index=False)
    tr.to_csv(out_dir / "recovery_tracking_110.csv", index=False)
    if not args.no_plots:
        try:
            if not fig_recovery(pt, fracs, args.within, out_dir / "fig_recovery_cc_110.png"):
                print("\n  (figure skipped: no cc-ON truths)")
        except Exception as e:
            print("\n  (figure skipped: %s -- the CSVs are saved regardless)" % e)

    def md5_of(p):
        return hashlib.md5(Path(p).read_bytes()).hexdigest()
    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "pandas": pd.__version__, "numpy": np.__version__,
            "inputs": {n: md5_of(summary_dir / n) for n in (P.STAGE2_NAME, P.CONTROL_NAME)},
            "settings": {k: (v if not isinstance(v, Path) else str(v)) for k, v in vars(args).items()},
            "truths": [float(t) for t in truths], "checks": ck.rows}
    (out_dir / "PROVENANCE_recovery_cc_110.json").write_text(json.dumps(prov, indent=2, default=str))

    print("\nSaved to: %s" % out_dir)
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results."
              % ck.n("FAIL"))
    edge = [fcc(t) for t in truths if t > 0 and not (P.NEAR_EDGE <= T.box_pos(t) <= 1 - P.NEAR_EDGE)]
    print("  Caveat: one storm, %d runs, routing pinned at truth. Medians of 12-50 runs are noisy (see the bootstrap interval in the CSV). "
          "Edge truths (%s) are pulled toward the middle by the sampled range itself."
          % (data[SER_ON]["X"].shape[0], ", ".join(edge) if edge else "none"))


if __name__ == "__main__":
    main()
