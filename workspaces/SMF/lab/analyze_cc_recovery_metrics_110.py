"""
analyze_cc_recovery_metrics_110.py
==================================
Series 110 -- does ranking the runs by more than KGE help the best runs find the true cc?

Every "best run" in the earlier scripts was chosen by KGE_2012 alone. The tables you saw show
that cc leaves fingerprints in metrics KGE does not contain (how long the flow stays above a
threshold, when it first arrives, how it recedes). This script asks the direct question: if the
best runs are picked using those other metrics too, do they land closer to the true cc?

It repeats the tracking-slope test of analyze_cc_recovery_110.py for five different ways of
choosing the best runs:

  1. KGE alone                 KGE_2012, exactly as before (reproduces the earlier output)
  2. KGE + timing metrics      half the vote to KGE, half shared by 4 timing/duration metrics
  3. timing metrics alone      only the 4 timing/duration metrics
  4. KGE + all 7 other metrics half the vote to KGE, half shared by the 4 timing/duration
                               metrics and 3 shape metrics
  5. all 7 other metrics alone only those 7 metrics

  timing/duration metrics:  peak timing error, first-arrival error, time-to-peak error,
                            duration-above-threshold error
  shape metrics:            peak magnitude error, rising-limb steepness ratio, recession-rate ratio

The weights are fixed in advance (equal weights, nothing tuned). The only choice is the 50/50 split
between KGE and the extras (--weight, default 0.5), so rows 2 and 4 are one point on a dial, and
rows 3 and 5 are the end of the dial where KGE has no vote at all.

Pure post-hoc analysis. Runs no tRIBS and does not read the hydrographs: it uses the metrics
already stored by the re-scoring (rescored_long_110.csv and the Stage 2 / control tables).
Writes only into its own new folder (calibration_work/03_comparisons/summary_tables/recovery_metrics_cc_110/).

Keep this file in the same folder as analyze_cc_recovery_110.py, analyze_cc_tests_110.py and
analyze_cc_pca_110.py. Run it from the lab/ directory.

VOCABULARY
-----------
  Percentile rank   Order the 250 runs from best to worst on one metric. Each run gets a score from
                    0 (best) to 1 (worst); ties share their average. The metric's units do not matter.
  Error             How far a metric is from a perfect fit: |value| for metrics that are 0 when perfect
                    (timing and magnitude errors), |value - 1| for ratios that are 1 when perfect.
  Ranking score     The weighted average of the percentile ranks of the metrics in that ranking. The
                    best runs are those with the lowest score (good on everything, not on one thing).
                    A run with a missing metric (for example, it never crosses the flow threshold) counts
                    as worst on that metric.
  Tracking slope, typical error, no-information, control, permutation test, bootstrap
                    As in analyze_cc_recovery_110.py. The same shuffles and the same bootstrap re-draws
                    of the runs are used for every ranking, so rankings can be compared directly.
  Paired difference The tracking slope of a ranking minus the slope of KGE alone, computed inside every
                    bootstrap re-draw. Its 95% interval is the honest way to say whether the new ranking
                    is better; the two separate slope intervals overlap even when the difference is clear,
                    because they share the same runs.

HOW TO READ THE RESULTS
------------------------
  * A ranking helps if its paired difference is positive with a 95% interval that stays above 0 at
    several cuts, and the control difference stays near 0. One interval above 0 among many is not enough:
    there are 4 rankings x 3 cuts x 2 truth sets = 24 intervals and none is corrected for that.
  * "timing alone" near 0 means the timing metrics do not carry cc by themselves.
  * A ranking that is worse than KGE alone usually means the extra metrics add noise or dilute the volume
    information that KGE is built on.
  * Best case only: no observation error, one storm, routing pinned at truth. The extra metrics are the
    ones the routing parameters own (recession rate: flowexp; first arrival and time to peak: channel
    roughness and hillslope velocity). A gain here is a ceiling for what these metrics could add when the
    routing parameters are free. No gain here would be a clear sign to stop.

OUTPUTS   (calibration_work/03_comparisons/summary_tables/recovery_metrics_cc_110/)
-----------------------------------------------------------------------------------
  recovery_metrics_tracking_110.csv    per ranking x cut x truth set x series: slope, interval, p, typical error
  recovery_metrics_difference_110.csv  per ranking vs KGE alone: paired slope and error differences with intervals
  recovery_metrics_per_truth_110.csv   per ranking x cut x series x truth: median cc of the best runs, middle 80%
  recovery_metrics_input_checks_110.csv
  fig_recovery_metrics_slope_110.png   slope for each ranking, one panel per cut and truth set
  fig_recovery_metrics_truth_110.png   best-run cc against true cc, one panel per ranking
  PROVENANCE_recovery_metrics_cc_110.json

USAGE (run from the lab/ directory)
------------------------------------
    python analyze_cc_recovery_metrics_110.py                 # about a minute
    python analyze_cc_recovery_metrics_110.py --weight 0.25   # give KGE 75% of the vote in rows 2 and 4
    python analyze_cc_recovery_metrics_110.py --boot 500 --perm 2000   # faster
    python analyze_cc_recovery_metrics_110.py --no_plots

CAVEATS
--------
One storm, 250 runs, routing pinned at truth, exact (noise-free) observations. The metrics were chosen
for this test after looking at how cc correlates with them, so a gain should be re-checked on new runs
before it is trusted. Timing and duration metrics move in 5-minute steps, so many runs tie.
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

_NEEDED_R = ["T", "P", "top_index", "slope_of", "tracking_sets", "fcc", "ffactor", "LOG_LO", "LOG_HI", "CENTER_LOG",
             "run_perm", "run_boot", "build_per_truth", "build_tracking", "MIN_TRACK_TRUTHS", "OUT_DIRNAME"]
_missing = [n for n in _NEEDED_R if not hasattr(R, n)]
if _missing:
    sys.exit("analyze_cc_recovery_110.py in this folder is an older version (missing: %s). Replace it with the current copy."
             % ", ".join(_missing))
T, P = R.T, R.P
_NEEDED_T = ["build_tensors", "p_low", "p_high", "fmt_p", "_style", "_ranks", "box_pos", "COL_CC", "COL_NULL", "INK", "INK2",
             "GRID", "SURFACE", "SER_ON", "SER_OFF", "OFF_TRUTH"]
_missing = [n for n in _NEEDED_T if not hasattr(T, n)]
if _missing:
    sys.exit("analyze_cc_tests_110.py in this folder is an older version (missing: %s). Replace it with the current copy."
             % ", ".join(_missing))

SER_ON, SER_OFF, OFF_TRUTH = T.SER_ON, T.SER_OFF, T.OFF_TRUTH
OUT_DIRNAME = "recovery_metrics_cc_110"
KGE_COL = "kge_2012"

TIMING = ["peak_timing_error_hr", "first_arrival_error_min", "time_to_peak_from_exc_min", "duration_above_thresh_error_min"]
SHAPE = ["peak_error_pct", "rising_limb_steepness_ratio", "recession_rate_ratio"]
EXTRA = TIMING + SHAPE
KIND = {"peak_timing_error_hr": "abs", "first_arrival_error_min": "abs", "time_to_peak_from_exc_min": "abs",
        "duration_above_thresh_error_min": "abs", "peak_error_pct": "abs",
        "rising_limb_steepness_ratio": "ratio", "recession_rate_ratio": "ratio"}
# key, label, short label, metrics, mode ("kge" = KGE_2012 as before, "mixed" = KGE gets (1 - weight), "alone" = no KGE vote)
RANKINGS = [
    ("kge", "KGE alone (as before)", "KGE alone", None, "kge"),
    ("kge_timing", "KGE + timing metrics", "KGE + timing", TIMING, "mixed"),
    ("timing", "timing metrics alone", "timing alone", TIMING, "alone"),
    ("kge_extra", "KGE + all 7 other metrics", "KGE + all 7", EXTRA, "mixed"),
    ("extra", "all 7 other metrics alone", "all 7 alone", EXTRA, "alone"),
]
MIN_FINITE_WARN = 0.90
MIN_FINITE_STOP = 0.50


# ------------------------------------------------------------------
# Ranking scores
# ------------------------------------------------------------------
def error_of(name, v):
    """Distance from a perfect fit: |v| for errors, |v - 1| for ratios. NaN stays NaN."""
    v = np.asarray(v, float)
    return np.abs(v - 1.0) if KIND[name] == "ratio" else np.abs(v)


def pct_rank(err):
    """Percentile rank of each run within its row (0 = best = smallest error, 1 = worst). Missing values count as worst.
    err has shape (truths, runs)."""
    e = np.where(np.isfinite(err), err, np.inf)
    return T._ranks(e) / e.shape[1]


def ranking_score(kge, metrics, group, weight):
    """Score for one ranking, shape (truths, runs), HIGHER = BETTER (so the existing top-set code can use it as it uses KGE).
    kge: (truths, runs) KGE_2012.  metrics: dict name -> (truths, runs).  group: list of metric names.
    weight: share of the vote given to the group (1 = the group alone, 0.5 = half to KGE and half to the group)."""
    g = np.mean([pct_rank(error_of(m, metrics[m])) for m in group], axis=0)
    if weight >= 1.0:
        return -g
    return -((1.0 - weight) * pct_rank(-kge) + weight * g)


def build_metric_tensors(df, truths, data):
    """Align the extra metrics with the runs of build_tensors(): out[series][metric] has shape (truths, runs)."""
    out = {}
    for s in (SER_ON, SER_OFF):
        order = data[s]["run_id"]
        out[s] = {}
        for m in EXTRA:
            M = np.full((len(truths), len(order)), np.nan)
            for ti, t in enumerate(truths):
                g = df[(df.truth_cc_mmhr == t) & (df.series == s)].set_index("run_id")[m]
                M[ti] = g.reindex(order).to_numpy(float)
            out[s][m] = M
    return out


def check_metrics(metrics, truths):
    """Availability of the extra metrics. Returns a list of (status, name, detail)."""
    rows = []
    for m in EXTRA:
        shares = []
        for s in (SER_ON, SER_OFF):
            for ti, t in enumerate(truths):
                shares.append(float(np.mean(np.isfinite(metrics[s][m][ti]))))
        lo, mean = min(shares), float(np.mean(shares))
        if mean < MIN_FINITE_STOP:
            rows.append(("FAIL", "metric available: %s" % m, "only %.0f%% of values are finite (column missing or blank: re-run rescore_cc_truths_110.py)"
                         % (100 * mean)))
        elif lo < MIN_FINITE_WARN:
            rows.append(("WARN", "metric available: %s" % m, "worst (truth, series) has %.0f%% finite; missing runs count as worst" % (100 * lo)))
        else:
            rows.append(("PASS", "metric available: %s" % m, "finite for %.0f%% or more of the runs in every truth and series" % (100 * lo)))
    return rows


# ------------------------------------------------------------------
# Running the five rankings
# ------------------------------------------------------------------
def run_rankings(data, metrics, truths, fracs, weight, half_width, Bp, B, seed, say=None):
    """One entry per ranking: the permutation null, the joint bootstrap, the per-truth table and the tracking table.
    The same shuffles and the same bootstrap re-draws are used for every ranking (same seed)."""
    out = {}
    for key, label, short, group, mode in RANKINGS:
        dr = {}
        for s in (SER_ON, SER_OFF):
            if mode == "kge":
                K = data[s]["K"]
            else:
                K = ranking_score(data[s]["K"], metrics[s], group, 1.0 if mode == "alone" else weight)
            dr[s] = {"run_id": data[s]["run_id"], "X": data[s]["X"], "K": K}
        null = R.run_perm(dr, truths, fracs, Bp, seed, half_width)
        boot = R.run_boot(dr, truths, fracs, B, seed)
        pt, pm = R.build_per_truth(dr, truths, fracs, half_width, null, boot)
        tr = R.build_tracking(truths, fracs, half_width, null, boot, pm)
        pt.insert(0, "ranking", key)
        tr.insert(0, "ranking", key)
        out[key] = {"label": label, "short": short, "boot": boot, "pt": pt, "tr": tr, "pm": pm}
        if say:
            say("  done: %s" % label)
    return out


def paired_differences(out, truths, fracs, base="kge"):
    """Ranking minus KGE alone, inside every bootstrap re-draw (the re-draws are shared)."""
    rows = []
    for key, o in out.items():
        if key == base:
            continue
        for s in (SER_ON, SER_OFF):
            for frac in fracs:
                for label, idx in R.tracking_sets(truths):
                    x = np.log10(np.array([truths[i] for i in idx]))
                    bb = out[base]["boot"][(s, frac)]["med"][:, idx]
                    br = o["boot"][(s, frac)]["med"][:, idx]
                    d_slope = R.slope_of(x, br) - R.slope_of(x, bb)
                    d_err = np.mean(np.abs(br - x), axis=1) - np.mean(np.abs(bb - x), axis=1)
                    pb, pr = out[base]["pm"][(s, frac)][idx], o["pm"][(s, frac)][idx]
                    rows.append({
                        "ranking": key, "series": s, "top_frac": frac, "truth_set": label, "n_truths": len(idx),
                        "slope_diff": float(R.slope_of(x, pr) - R.slope_of(x, pb)),
                        "slope_diff_ci_lo": float(np.percentile(d_slope, 2.5)), "slope_diff_ci_hi": float(np.percentile(d_slope, 97.5)),
                        "share_diff_positive": float(np.mean(d_slope > 0)),
                        "error_ratio": float(10 ** (np.mean(np.abs(pr - x)) - np.mean(np.abs(pb - x)))),
                        "error_ratio_ci_lo": float(10 ** np.percentile(d_err, 2.5)),
                        "error_ratio_ci_hi": float(10 ** np.percentile(d_err, 97.5)),
                    })
    return pd.DataFrame(rows)


def compare_with_saved(out, summary_dir):
    """The 'KGE alone' slopes must equal the ones saved by analyze_cc_recovery_110.py."""
    p = summary_dir / R.OUT_DIRNAME / "recovery_tracking_110.csv"
    if not p.exists():
        return "no saved recovery_tracking_110.csv found to compare with (run analyze_cc_recovery_110.py first if you want this check)."
    old = pd.read_csv(p)
    new = out["kge"]["tr"][["series", "top_frac", "truth_set", "slope"]]
    m = new.merge(old[["series", "top_frac", "truth_set", "slope"]], on=["series", "top_frac", "truth_set"], how="inner", suffixes=("", "_old"))
    if m.empty:
        return "the saved recovery_tracking_110.csv has no matching rows (different cuts or truths), so no comparison was made."
    d = float(np.max(np.abs(m["slope"].to_numpy(float) - m["slope_old"].to_numpy(float))))
    if d < 1e-9:
        return "'KGE alone' reproduces the slopes saved by analyze_cc_recovery_110.py (largest difference %.1e, %d rows)." % (d, len(m))
    return ("NOTE: 'KGE alone' differs from the saved recovery_tracking_110.csv by up to %.3g (%d rows). That file may come from "
            "different settings (cuts, --truths) or different inputs." % (d, len(m)))


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------
def print_tracking(out, truths, fracs, n_by_frac, weight):
    sets = [lab for lab, _ in R.tracking_sets(truths)]
    for frac in fracs:
        print("\n" + "=" * 140)
        print("DO THE BEST RUNS FOLLOW THE TRUE cc?  best %d%% of runs (%d runs), by each ranking" % (round(100 * frac), n_by_frac[frac]))
        print("=" * 140)
        print("%-17s %-27s | %-26s %7s | %-18s | %-26s"
              % ("truths", "ranking", "ON: slope [95% interval]", "p", "typ. error (no-info)", "control: slope [95% interval]"))
        for lab in sets:
            for key, label, short, group, mode in RANKINGS:
                tr = out[key]["tr"]
                a = tr[(tr.top_frac == frac) & (tr.truth_set == lab) & (tr.series == SER_ON)]
                c = tr[(tr.top_frac == frac) & (tr.truth_set == lab) & (tr.series == SER_OFF)]
                if a.empty or c.empty:
                    continue
                a, c = a.iloc[0], c.iloc[0]
                terr = "x%.2f (x%.2f)" % (a.typical_error_factor, a.null_typical_error_factor)
                print("%-17s %-27s | %+6.2f [%+5.2f, %+5.2f]      %7s | %-20s | %+6.2f [%+5.2f, %+5.2f]"
                      % (lab if key == "kge" else "", short, a.slope, a.slope_ci_lo, a.slope_ci_hi, T.fmt_p(a.p_slope),
                         terr, c.slope, c.slope_ci_lo, c.slope_ci_hi))
            print()
    print("  slope: 1 = the best runs' median cc follows the true cc one-for-one, 0 = it ignores it.  typ. error: typical factor by which the")
    print("  median misses the true cc (x1.00 = perfect), and in brackets what a no-information fit gives.  The control (cc inert) must stay near 0.")
    print("  Mixed rankings give KGE %d%% of the vote and share the other %d%% equally among the extra metrics; 'alone' rankings give KGE no vote."
          % (round(100 * (1 - weight)), round(100 * weight)))


def print_differences(diff, truths, fracs):
    sets = [lab for lab, _ in R.tracking_sets(truths)]
    print("\n" + "=" * 140)
    print("IS A RANKING BETTER THAN KGE ALONE?  paired difference in tracking slope (ranking minus KGE alone), 95% interval from the shared bootstrap")
    print("=" * 140)
    print("%-5s %-17s %-17s | %-30s %6s | %-22s | %-30s"
          % ("cut", "truths", "ranking", "ON: slope difference", ">0 in", "error ratio vs KGE", "control: slope difference"))
    for frac in fracs:
        for lab in sets:
            for key, label, short, group, mode in RANKINGS:
                if key == "kge":
                    continue
                a = diff[(diff.ranking == key) & (diff.top_frac == frac) & (diff.truth_set == lab) & (diff.series == SER_ON)]
                c = diff[(diff.ranking == key) & (diff.top_frac == frac) & (diff.truth_set == lab) & (diff.series == SER_OFF)]
                if a.empty or c.empty:
                    continue
                a, c = a.iloc[0], c.iloc[0]
                flag = ""
                if a.slope_diff_ci_lo > 0:
                    flag = "  better"
                elif a.slope_diff_ci_hi < 0:
                    flag = "  worse"
                er = "x%.2f [%.2f, %.2f]" % (a.error_ratio, a.error_ratio_ci_lo, a.error_ratio_ci_hi)
                print("%-5s %-17s %-17s | %+6.2f [%+5.2f, %+5.2f]        %4.0f%% | %-20s | %+6.2f [%+5.2f, %+5.2f]%s"
                      % ("%d%%" % round(100 * frac), lab if key == "kge_timing" else "", short, a.slope_diff, a.slope_diff_ci_lo,
                         a.slope_diff_ci_hi, 100 * a.share_diff_positive, er,
                         c.slope_diff, c.slope_diff_ci_lo, c.slope_diff_ci_hi, flag))
        print()
    print("  'better' / 'worse' = the 95% interval of the slope difference lies entirely above / below 0.  >0 in = share of bootstrap re-draws in which")
    print("  the ranking's slope beats KGE alone.  error ratio: typical error of the ranking divided by that of KGE alone (below 1 = closer to the truth).")
    print("  24 intervals are shown and none is corrected for that, so look for the same sign at several cuts, not for one 'better'.")
    print("  The control difference should straddle 0: the control cannot feel cc, so a ranking cannot track it.")


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def fig_slopes(out, truths, fracs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sets = [lab for lab, _ in R.tracking_sets(truths)]
    if not sets:
        return False
    nr, nc = len(sets), len(fracs)
    order = [k for k, *_ in RANKINGS]
    allv = pd.concat([out[k]["tr"] for k in order])
    xmin = min(-0.4, float(allv.slope_ci_lo.min()) - 0.05)
    xmax = max(1.2, float(allv.slope_ci_hi.max()) + 0.05)
    fig, axes = plt.subplots(nr, nc, figsize=(4.4 * nc + 1.3, 3.1 * nr + 1.3), facecolor=T.SURFACE, sharex=True, sharey=True, squeeze=False)
    ys = np.arange(len(order))[::-1]
    for ri, lab in enumerate(sets):
        for ci, frac in enumerate(fracs):
            ax = axes[ri][ci]
            T._style(ax)
            ax.axvline(1.0, color=T.INK, linestyle="--", linewidth=1.0, zorder=2)
            ax.axvline(0.0, color=T.COL_NULL, linestyle=":", linewidth=1.4, zorder=2)
            for s, col, off, mk, filled in ((SER_OFF, T.COL_NULL, -0.17, "s", False), (SER_ON, T.COL_CC, +0.17, "o", True)):
                for k, y in zip(order, ys):
                    tr = out[k]["tr"]
                    g = tr[(tr.top_frac == frac) & (tr.truth_set == lab) & (tr.series == s)]
                    if g.empty:
                        continue
                    g = g.iloc[0]
                    ax.hlines(y + off, g.slope_ci_lo, g.slope_ci_hi, color=col, linewidth=2.2 if s == SER_ON else 1.5, zorder=3)
                    ax.plot([g.slope], [y + off], linestyle="none", marker=mk, markersize=6, markerfacecolor=col if filled else T.SURFACE,
                            markeredgecolor=col, markeredgewidth=1.5, zorder=4)
            ax.set_yticks(ys)
            ax.set_yticklabels([out[k]["short"] for k in order], fontsize=8.5)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(-0.6, len(order) - 0.4)
            nt = int(out["kge"]["tr"][(out["kge"]["tr"].truth_set == lab)].n_truths.iloc[0])
            ax.set_title("Best %d%% of runs, %s (%d)" % (round(100 * frac), lab, nt), color=T.INK, fontsize=9, loc="left")
            if ri == nr - 1:
                ax.set_xlabel("tracking slope (95% bootstrap interval)", color=T.INK2, fontsize=9)
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=6, linewidth=2.2, label="channel loss ON"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=6, markerfacecolor=T.SURFACE, markeredgewidth=1.5, linewidth=1.5,
               label="control (cc inert)"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=1.0, label="1 = follows the true cc one-for-one"),
        Line2D([0], [0], color=T.COL_NULL, linestyle=":", linewidth=1.4, label="0 = ignores cc"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Does ranking the runs by more than KGE help the best runs find the true cc?", color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.07 if nr > 1 else 0.10, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


def fig_truth(out, frac, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    order = [k for k, *_ in RANKINGS]
    d0 = out["kge"]["pt"]
    tv = sorted(d0[d0.truth_cc_mmhr > 0].truth_cc_mmhr.unique())
    if not tv:
        return False
    lo, hi = 10 ** R.LOG_LO, 10 ** R.LOG_HI
    lim = (lo * 0.75, hi * 1.3)
    n = len(order)
    fig, axes = plt.subplots(1, n, figsize=(3.7 * n + 0.5, 5.6), facecolor=T.SURFACE, sharey=True, squeeze=False)
    for ci, k in enumerate(order):
        ax = axes[0][ci]
        T._style(ax)
        ax.axhspan(10 ** (R.LOG_LO + 0.1 * (R.LOG_HI - R.LOG_LO)), 10 ** (R.LOG_LO + 0.9 * (R.LOG_HI - R.LOG_LO)),
                   color=T.COL_NULL, alpha=0.16, linewidth=0, zorder=1)
        ax.axhline(10 ** R.CENTER_LOG, color=T.COL_NULL, linestyle=":", linewidth=1.3, zorder=2)
        ax.plot(lim, lim, color=T.INK, linestyle="--", linewidth=1.0, zorder=2)
        pt = out[k]["pt"]
        for s, col, off, mk, filled in ((SER_OFF, T.COL_NULL, +1, "s", False), (SER_ON, T.COL_CC, -1, "o", True)):
            g = pt[(pt.series == s) & (pt.top_frac == frac) & (pt.truth_cc_mmhr > 0)].sort_values("truth_cc_mmhr")
            if g.empty:
                continue
            x = g.truth_cc_mmhr.values * 10 ** (off * 0.012)
            ax.vlines(x, g.q10_cc, g.q90_cc, color=col, linewidth=2.0 if s == SER_ON else 1.4, zorder=3)
            ax.plot(x, g.median_cc.values, linestyle="none", marker=mk, markersize=5.2, markerfacecolor=col if filled else T.SURFACE,
                    markeredgecolor=col, markeredgewidth=1.5, zorder=4)
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
        ax.set_title(out[k]["short"], color=T.INK, fontsize=9.5, loc="left")
        ax.set_xlabel("true cc of the truth (mm/hr)", color=T.INK2, fontsize=9)
        if ci == 0:
            ax.set_ylabel("cc of the best runs (mm/hr)", color=T.INK2, fontsize=9)
    nrun = int(d0[(d0.series == SER_ON) & (d0.top_frac == frac)].n_top.iloc[0])
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5.5, linewidth=2.0, label="channel loss ON: median and middle 80% of the best runs"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=5.5, markerfacecolor=T.SURFACE, markeredgewidth=1.5, linewidth=1.4,
               label="control (cc inert): the same"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=1.0, label="perfect recovery (median sits on the true cc)"),
        Patch(facecolor=T.COL_NULL, alpha=0.16, label="no information: middle 80% of the sampled range (dotted = its middle)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Best %d%% of runs (n = %d) by each ranking: where do they put cc?  (dot = median, bar = middle 80%%)" % (round(100 * frac), nrun),
                 color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.11, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110 -- does ranking the runs by more than KGE help the best runs find the true cc?")
    ap.add_argument("--top_fracs", type=float, nargs="+", default=[0.20, 0.10, 0.05],
                    help="Top fractions; the first is the main one (default 0.20 0.10 0.05)")
    ap.add_argument("--weight", type=float, default=0.5,
                    help="Share of the vote given to the extra metrics in the 'KGE + ...' rankings (default 0.5; must be between 0 and 1, exclusive)")
    ap.add_argument("--within", type=float, default=2.0,
                    help="'Near the truth' means within this FACTOR of the true cc (default 2); only used in the per-truth CSV")
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
    if not 0.0 < args.weight < 1.0:
        sys.exit("--weight must be between 0 and 1 (exclusive), e.g. 0.5 = half the vote to KGE and half to the extra metrics.")
    if args.within <= 1.0:
        sys.exit("--within must be greater than 1 (it is a factor, e.g. 2 = between half and double).")
    half_width = float(np.log10(args.within))

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    out_dir = summary_dir / OUT_DIRNAME
    print("\n" + "=" * 78 + "\nSeries 110 -- does ranking the runs by more than KGE help the best runs find the true cc?\n" + "=" * 78)
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
    fracs = args.top_fracs
    data = T.build_tensors(df, truths)
    metrics = build_metric_tensors(df, truths, data)
    for status, name, detail in check_metrics(metrics, truths):
        ck.add(status, name, detail)
    P.print_checks(ck)
    if ck.n("FAIL") and not args.continue_on_fail:
        sys.exit("\nSTOPPED: an input check failed. Nothing was analysed or written. (--continue_on_fail to override.)")

    sets = R.tracking_sets(truths)
    if not sets:
        sys.exit("STOPPED: the tracking slope needs at least %d cc-ON truths." % R.MIN_TRACK_TRUTHS)
    n_runs = data[SER_ON]["X"].shape[0]
    n_by_frac = {f: T.n_top_for(n_runs, f) for f in fracs}
    print("\n%d cc-ON truth(s), %d runs per series, cuts %s, %d shuffles, %d bootstrap re-draws, KGE vote in mixed rankings %d%%"
          % (sum(1 for t in truths if t > 0), n_runs, ["%d%%" % round(100 * f) for f in fracs], args.perm, args.boot,
             round(100 * (1 - args.weight))))
    print("\nRankings compared:")
    for key, label, short, group, mode in RANKINGS:
        what = "KGE_2012" if group is None else "%d metric%s" % (len(group), "" if len(group) == 1 else "s")
        print("  %-28s %s" % (label, what))
    print("\nRanking the runs and re-running the tracking test for each ranking ...")
    out = run_rankings(data, metrics, truths, fracs, args.weight, half_width, args.perm, args.boot, args.seed, say=print)
    diff = paired_differences(out, truths, fracs)

    print_tracking(out, truths, fracs, n_by_frac, args.weight)
    print_differences(diff, truths, fracs)
    print("\n  " + compare_with_saved(out, summary_dir))

    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ck.rows).to_csv(out_dir / "recovery_metrics_input_checks_110.csv", index=False)
    pd.concat([out[k]["tr"] for k, *_ in RANKINGS], ignore_index=True).to_csv(out_dir / "recovery_metrics_tracking_110.csv", index=False)
    diff.to_csv(out_dir / "recovery_metrics_difference_110.csv", index=False)
    pd.concat([out[k]["pt"] for k, *_ in RANKINGS], ignore_index=True).to_csv(out_dir / "recovery_metrics_per_truth_110.csv", index=False)
    if not args.no_plots:
        try:
            fig_slopes(out, truths, fracs, out_dir / "fig_recovery_metrics_slope_110.png")
            fig_truth(out, fracs[0], out_dir / "fig_recovery_metrics_truth_110.png")
        except Exception as e:
            print("\n  (figures skipped: %s -- the CSVs are saved regardless)" % e)

    def md5_of(p):
        return hashlib.md5(Path(p).read_bytes()).hexdigest()
    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "pandas": pd.__version__, "numpy": np.__version__,
            "inputs": {n: md5_of(summary_dir / n) for n in (P.STAGE2_NAME, P.CONTROL_NAME)},
            "rankings": {k: {"label": lab, "metrics": g, "mode": mode} for k, lab, short, g, mode in RANKINGS},
            "settings": {k: (v if not isinstance(v, Path) else str(v)) for k, v in vars(args).items()},
            "truths": [float(t) for t in truths], "checks": ck.rows}
    (out_dir / "PROVENANCE_recovery_metrics_cc_110.json").write_text(json.dumps(prov, indent=2, default=str))

    print("\nSaved to: %s" % out_dir)
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results." % ck.n("FAIL"))
    print("  Caveat: one storm, %d runs, routing pinned at truth, exact observations: a best case. The extra metrics are the ones the routing "
          "parameters own, so a gain here is a ceiling; no gain here means stop." % n_runs)


if __name__ == "__main__":
    main()
