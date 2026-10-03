"""
analyze_cc_recovery_noise_110.py
================================
Series 110 -- does "the best runs follow the true cc" survive an imperfect observed hydrograph?

analyze_cc_recovery_110.py showed that, when the "observed" hydrograph is the exact model output of
a known cc, the best runs' cc follows the true cc (tracking slope 0.49 / 0.65 / 0.83 for the best
20% / 10% / 5% of runs, control near 0). A real gauge record is never exact. This script repeats
that check after adding observation error to the "observed" hydrograph, and asks whether the
tracking survives.

Pure post-hoc analysis. Runs no tRIBS. It re-reads the 500 stored hydrographs
(<run_id>_compare_obs_sim.csv), scores each against a noisy copy of every truth, picks the best
runs, and measures the tracking slope again. Writes only into its own new folder
(calibration_work/03_comparisons/summary_tables/recovery_noise_cc_110/).

Keep this file in the same folder as analyze_cc_recovery_110.py, analyze_cc_tests_110.py,
analyze_cc_pca_110.py and rescore_cc_truths_110.py. Run it from the lab/ directory.

WHAT IT DOES
------------
For each noise level (none, then 5%, 10%, 20% by default) and for each of 100 noisy replicates:
  1. For every cc-ON truth, build a noisy "observed" hydrograph from that truth (noise model below).
     Every truth gets its OWN independent noise draw, like a separate experiment with its own
     gauge error. (The earlier noise test in analyze_cc_tests_110.py used one draw shared by all
     truths; independent draws are the harder, more realistic case for a slope across truths.)
  2. Score all 250 cc-ON runs and all 250 control runs (cc inert) against that noisy hydrograph
     with KGE_2012, exactly as the main analysis does. The control sees the same noisy hydrograph.
  3. Take the best 20% / 10% / 5% by KGE_2012, find the median cc of each truth's best runs, and fit
     the tracking slope of (median best-run cc) against (true cc) on log scales.
  4. Shuffle cc among the runs (1000 shuffles, one shuffle shared by all truths) to get the p-value
     for "could a fit that ignores cc give a slope this large?".
Then it summarizes the slope over the replicates: its typical value, how much it varies with the
luck of the noise draw, and in what share of replicates it is still significant.

NOISE MODEL (the same as the earlier noise test)
------------------------------------------------
    observed = truth x (1 + b) x exp(e)
    b ~ Normal(0, sd)        one error for the whole hydrograph (a volume or rating-curve error)
    e ~ Normal(0, sd / 2)    independent at every 5-minute step (measurement scatter)
sd is 5%, 10% and 20% by default. --noise_parts volume uses only b; --noise_parts shape uses only e.
The whole-hydrograph error b is the part most likely to hurt: cc and Ks both change the runoff
volume, so a volume error can be traded against either of them.

VOCABULARY
-----------
  Replicate         One noisy "observed" data set. 100 replicates show how much the result changes
                    with the luck of the noise draw.
  Middle 80% of replicates   The range holding the central 80% of the replicate values (10th to 90th
                    percentile). It measures the effect of noise, not of the choice of runs.
  Tracking slope    As in analyze_cc_recovery_110.py: slope of (median best-run cc) against (true cc),
                    log scales. 1 = follows the truth one-for-one; 0 = ignores it. The edges of the
                    sampled range squash it below 1 even for a perfect method.
  Significant       The share of replicates in which shuffling cc almost never (p < 0.05) gives a slope
                    this large. With no noise it is a single yes or no.
  Typical error     The factor by which the best runs' median cc typically misses the true cc
                    (x1.00 = perfect), next to what a no-information fit gives.
  Control           Runs whose cc is inert. Its slope must stay near 0 at every noise level; if it
                    does not, noise alone is creating a pattern.

HOW TO READ THE RESULTS
------------------------
  * Slope stays near its no-noise value and the share significant stays high: the tracking survives
    an imperfect observed hydrograph at that noise level.
  * Slope shrinks toward 0 as noise grows: the tracking is fragile. The noise level where its middle
    80% first reaches 0 is the limit of what the data can tell you about cc.
  * Slope stays positive but the share significant falls: the average tracking survives but a
    single real data set could easily fail to show it.
  * Rates and slopes for 5/10/20% say what this one storm and these 250 runs can do. They do not
    make a real basin's gauge error equal to 5%, 10% or 20%.

OUTPUTS   (calibration_work/03_comparisons/summary_tables/recovery_noise_cc_110/)
-------------------------------------------------------------------------------
  recovery_noise_summary_110.csv       per cut x truth set x series x noise level: slope, spread over
                                       replicates, share significant, typical error
  recovery_noise_per_truth_110.csv     per cut x series x noise level x truth: median best-run cc over
                                       replicates and its middle 80%
  recovery_noise_realizations_110.csv  every replicate (slope, p-value, typical error)
  recovery_noise_input_checks_110.csv
  fig_recovery_noise_slope_110.png     slope against noise level, one panel per cut and truth set
  fig_recovery_noise_truth_110.png     best-run cc against true cc, one panel per noise level
  PROVENANCE_recovery_noise_cc_110.json
(With --noise_parts volume or shape the file names end in _volume or _shape.)

USAGE (run from the lab/ directory)
------------------------------------
    python analyze_cc_recovery_noise_110.py                      # defaults, a minute or two
    python analyze_cc_recovery_noise_110.py --noise_reps 200 --perm 2000
    python analyze_cc_recovery_noise_110.py --noise_levels 0.02 0.05 0.10
    python analyze_cc_recovery_noise_110.py --noise_parts volume   # only the whole-hydrograph error
    python analyze_cc_recovery_noise_110.py --no_plots

CAVEATS
--------
One storm, 250 runs, routing pinned at truth: a best case. Replicates reuse the same 250 runs, so
the control column is not the 5% calibration of the permutation test (that was checked on fresh
random designs). The noise is multiplicative and has no memory in time; real gauge errors can be
larger at high flows and correlated in time.
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
             "folder (it reuses their loader and input checks). Put the files in lab/ and run from there.")
try:
    import rescore_cc_truths_110 as RS
except ImportError:
    sys.exit("This script needs rescore_cc_truths_110.py in the same folder (it re-reads the stored hydrographs).")

_NEEDED_R = ["T", "P", "top_index", "slope_of", "tracking_sets", "truth_log", "fcc", "ffactor", "LOG_LO", "LOG_HI", "CENTER_LOG"]
_missing = [n for n in _NEEDED_R if not hasattr(R, n)]
if _missing:
    sys.exit("analyze_cc_recovery_110.py in this folder is an older version (missing: %s). Replace it with the current copy."
             % ", ".join(_missing))
T, P = R.T, R.P
_NEEDED_T = ["kge_matrix", "build_tensors", "make_perm_matrix", "p_low", "p_high", "fmt_p", "_style", "box_pos", "COL_CC", "COL_NULL",
             "INK", "INK2", "GRID", "SURFACE", "SER_ON", "SER_OFF", "OFF_TRUTH"]
_missing = [n for n in _NEEDED_T if not hasattr(T, n)]
if _missing:
    sys.exit("analyze_cc_tests_110.py in this folder is an older version (missing: %s). Replace it with the current copy."
             % ", ".join(_missing))
_NEEDED_RS = ["load_all_compares", "read_truth_5min", "md5_of", "CAND_DIRNAME"]
_missing = [n for n in _NEEDED_RS if not hasattr(RS, n)]
if _missing:
    sys.exit("rescore_cc_truths_110.py in this folder is an older version (missing: %s). Replace it with the current copy."
             % ", ".join(_missing))

SER_ON, SER_OFF, OFF_TRUTH = T.SER_ON, T.SER_OFF, T.OFF_TRUTH
OUT_DIRNAME = "recovery_noise_cc_110"
ALPHA = 0.05
NOISE_FREE_TOL = 1e-6
LEVEL_NAME = {0.0: "none"}


def level_label(lv):
    return LEVEL_NAME.get(lv, "%g%%" % round(100 * lv, 6))


# ------------------------------------------------------------------
# Stored hydrographs and truth series
# ------------------------------------------------------------------
def load_hydrographs(summary_dir, data, truths):
    """Simulated hydrographs of every run (runs x time) for both series, and the truth hydrograph of every truth,
    on the compare-CSV time index. The truths are verified against the checksums recorded at re-scoring."""
    calib = summary_dir.parent.parent
    csv_dir = calib / "03_comparisons" / "csv_exports"
    cand_dir = calib / RS.CAND_DIRNAME
    prov_path = summary_dir / P.RESCORE_DIR / P.RESCORE_PROV
    sims, idx0, obs_off = {}, None, None
    for s in (SER_ON, SER_OFF):
        res = pd.DataFrame({"run_id": data[s]["run_id"]})
        cache, prob = RS.load_all_compares(res, csv_dir)
        if prob:
            sys.exit("STOPPED: %d compare CSV(s) missing or unreadable, e.g. %s\n  Looked in: %s\n  This script re-reads every run's "
                     "stored hydrograph (<run_id>_compare_obs_sim.csv). Nothing was written." % (len(prob), list(prob.items())[0], csv_dir))
        mats = []
        for rid in data[s]["run_id"]:
            d = cache[rid]
            if idx0 is None:
                idx0 = d.index
            if not d.index.equals(idx0):
                sys.exit("STOPPED: compare CSV %s has a different time index from the others; refusing." % rid)
            mats.append(d["Simulated"].to_numpy(float))
        sims[s] = np.vstack(mats)
        if s == SER_ON:
            obs_off = cache[data[s]["run_id"][0]]["Observed"].to_numpy(float)

    prov = json.loads(prov_path.read_text()) if prov_path.exists() else {}
    cand = {}
    for name, info in prov.get("candidate_truths", {}).items():
        p = cand_dir / name
        if not p.exists():
            sys.exit("STOPPED: candidate truth file missing: %s" % p)
        if RS.md5_of(p) != info["md5"]:
            sys.exit("STOPPED: candidate truth %s no longer matches its recorded checksum; refusing." % name)
        s5 = RS.read_truth_5min(p).reindex(idx0)
        if s5.isna().any():
            sys.exit("STOPPED: truth %s has no value at some compare-CSV timestamps." % name)
        cand[float(info["cc_mmhr"])] = s5.to_numpy(float)

    truth_obs = {}
    for ti, t in enumerate(truths):
        if t == OFF_TRUTH:
            truth_obs[ti] = obs_off
            continue
        hit = [k for k in cand if np.isclose(t, k, rtol=1e-6)]
        if not hit:
            sys.exit("STOPPED: truth cc = %g has no verified candidate hydrograph in %s." % (t, cand_dir))
        truth_obs[ti] = cand[hit[0]]

    worst = 0.0
    for ti in truth_obs:
        for s in (SER_ON, SER_OFF):
            k = T.kge_matrix(sims[s], truth_obs[ti])
            stored = data[s]["K"][ti]
            ok = np.isfinite(stored)
            worst = max(worst, float(np.max(np.abs(k[ok] - stored[ok]))))
    if worst > NOISE_FREE_TOL:
        sys.exit("STOPPED: re-computing KGE_2012 without noise differs from the stored values by up to %.3g. The noisy results "
                 "would not be comparable to the main analysis." % worst)
    return sims, truth_obs, worst


def noisy_copy(truth, lv, parts, seed, li, rep, t):
    """One noisy observed hydrograph. The generator depends on (seed, level, replicate, truth value) only, so results do not
    change when --truths selects a subset."""
    rng = np.random.default_rng([seed, 88, li, rep, int(round(t * 1000))])
    b = rng.normal(0.0, lv)
    e = rng.normal(0.0, lv / 2.0, truth.shape[0])
    if parts == "shape":
        b = 0.0
    elif parts == "volume":
        e = np.zeros(truth.shape[0])
    return truth * max(1.0 + b, 0.05) * np.exp(e)


# ------------------------------------------------------------------
# One replicate
# ------------------------------------------------------------------
def one_replicate(Kmat, data, truths, fracs, sets, perm):
    """Tracking statistics for one noisy data set. Kmat[s] has shape (truths, runs). perm[s] is a (Bp, runs) shuffle matrix
    (the same shuffle is applied to every truth). Returns rows and the best-run median log10 cc per (series, cut, truth)."""
    rows, meds = [], {}
    for s in (SER_ON, SER_OFF):
        cc = data[s]["X"][:, 2]
        Pm = perm[s]
        for frac in fracs:
            pt = np.full(len(truths), np.nan)
            nl = np.full((Pm.shape[0], len(truths)), np.nan)
            for ti, t in enumerate(truths):
                if t <= 0:
                    continue
                idx = R.top_index(Kmat[s][ti], frac)
                pt[ti] = np.median(cc[idx])
                nl[:, ti] = np.median(cc[Pm[:, idx]], axis=1)
            meds[(s, frac)] = pt
            for label, sel in sets:
                x = np.log10(np.array([truths[i] for i in sel]))
                obs = pt[sel]
                sl = float(R.slope_of(x, obs))
                err = float(np.mean(np.abs(obs - x)))
                nsel = nl[:, sel]
                sl_null = R.slope_of(x, nsel)
                err_null = np.mean(np.abs(nsel - x), axis=1)
                rows.append({"series": s, "top_frac": frac, "truth_set": label, "n_truths": len(sel), "slope": sl,
                             "slope_null_mean": float(sl_null.mean()), "p_slope": T.p_high(sl_null, sl),
                             "typical_error_factor": float(10 ** err), "null_typical_error_factor": float(10 ** err_null.mean()),
                             "p_error_smaller": T.p_low(err_null, err)})
    return rows, meds


def run_replicates(args, data, truths, fracs, sets, sims, truth_obs):
    n = sims[SER_ON].shape[0]
    T_len = sims[SER_ON].shape[1]
    perm = {s: T.make_perm_matrix(n, args.perm, np.random.default_rng([args.seed, 99, si]))
            for si, s in enumerate((SER_ON, SER_OFF))}
    levels = [0.0] + list(args.noise_levels)
    rz, tz = [], []
    for li, lv in enumerate(levels):
        reps = 1 if lv == 0 else args.noise_reps
        for r in range(reps):
            if lv == 0:
                Kmat = {s: data[s]["K"] for s in (SER_ON, SER_OFF)}      # the stored values: identical to the main analysis
            else:
                Kmat = {s: np.full((len(truths), n), np.nan) for s in (SER_ON, SER_OFF)}
                for ti, t in enumerate(truths):
                    if t <= 0:
                        continue
                    o = noisy_copy(truth_obs[ti], lv, args.noise_parts, args.seed, li, r, t)
                    for s in (SER_ON, SER_OFF):
                        Kmat[s][ti] = T.kge_matrix(sims[s], o)
            rows, meds = one_replicate(Kmat, data, truths, fracs, sets, perm)
            for row in rows:
                row.update({"noise_level": lv, "rep": r})
                rz.append(row)
            for (s, frac), pt in meds.items():
                for ti, t in enumerate(truths):
                    if t > 0:
                        tz.append({"noise_level": lv, "rep": r, "series": s, "top_frac": frac, "truth_cc_mmhr": t,
                                   "median_log10_cc": float(pt[ti])})
        print("  noise %-5s done (%d replicate%s)" % (level_label(lv), reps, "" if reps == 1 else "s"))
    return pd.DataFrame(rz), pd.DataFrame(tz)


# ------------------------------------------------------------------
# Summaries
# ------------------------------------------------------------------
def summarize(rz):
    rows = []
    for (lv, frac, s, label), g in rz.groupby(["noise_level", "top_frac", "series", "truth_set"], sort=True):
        sl = g["slope"].to_numpy(float)
        rows.append({
            "noise_level": lv, "top_frac": frac, "series": s, "truth_set": label, "n_truths": int(g["n_truths"].iloc[0]), "n_reps": len(g),
            "slope_median": float(np.median(sl)), "slope_p10": float(np.percentile(sl, 10)), "slope_p90": float(np.percentile(sl, 90)),
            "share_significant": float(np.mean(g["p_slope"].to_numpy(float) < ALPHA)),
            "share_slope_positive": float(np.mean(sl > 0)),
            "typical_error_median": float(np.median(g["typical_error_factor"])),
            "null_typical_error_median": float(np.median(g["null_typical_error_factor"])),
            "share_error_smaller": float(np.mean(g["p_error_smaller"].to_numpy(float) < ALPHA)),
        })
    return pd.DataFrame(rows)


def per_truth_table(tz):
    rows = []
    for (lv, frac, s, t), g in tz.groupby(["noise_level", "top_frac", "series", "truth_cc_mmhr"], sort=True):
        v = g["median_log10_cc"].to_numpy(float)
        rows.append({"noise_level": lv, "top_frac": frac, "series": s, "truth_cc_mmhr": t, "n_reps": len(v),
                     "median_cc": float(10 ** np.median(v)), "cc_p10": float(10 ** np.percentile(v, 10)),
                     "cc_p90": float(10 ** np.percentile(v, 90)), "factor_vs_truth": float(10 ** (np.median(v) - np.log10(t)))})
    return pd.DataFrame(rows)


def compare_with_saved(summ, summary_dir):
    """The no-noise rows should reproduce the slopes saved by analyze_cc_recovery_110.py. Returns a message."""
    p = summary_dir / R.OUT_DIRNAME / "recovery_tracking_110.csv"
    if not p.exists():
        return "no saved recovery_tracking_110.csv found to compare with (run analyze_cc_recovery_110.py first if you want this check)."
    old = pd.read_csv(p)
    new = summ[summ.noise_level == 0.0][["series", "top_frac", "truth_set", "slope_median"]]
    m = new.merge(old[["series", "top_frac", "truth_set", "slope"]], on=["series", "top_frac", "truth_set"], how="inner")
    if m.empty:
        return "the saved recovery_tracking_110.csv has no matching rows (different cuts or truths), so no comparison was made."
    d = float(np.max(np.abs(m["slope_median"].to_numpy(float) - m["slope"].to_numpy(float))))
    if d < 1e-9:
        return "the no-noise rows reproduce the slopes saved by analyze_cc_recovery_110.py (largest difference %.1e, %d rows)." % (d, len(m))
    return ("NOTE: the no-noise slopes differ from the saved recovery_tracking_110.csv by up to %.3g (%d rows). That file may come from "
            "different settings (cuts, --truths) or different inputs." % (d, len(m)))


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------
def print_summary(summ, fracs, levels, n_by_frac):
    for frac in fracs:
        print("\n" + "=" * 132)
        print("RECOVERY WITH NOISY OBSERVATIONS -- best %d%% of runs by KGE_2012 (%d runs)" % (round(100 * frac), n_by_frac[frac]))
        print("=" * 132)
        print("%-18s %-6s | %-24s %-11s | %-20s | %-24s %-11s"
              % ("truths", "noise", "ON: tracking slope", "significant", "typ. error (no-info)", "control: slope", "significant"))
        print("%-18s %-6s | %-24s %-11s | %-20s | %-24s %-11s"
              % ("", "", "[middle 80% of repl.]", "", "", "[middle 80% of repl.]", ""))
        d = summ[summ.top_frac == frac]
        for label in d.truth_set.unique():
            for lv in levels:
                a = d[(d.truth_set == label) & (d.series == SER_ON) & (d.noise_level == lv)]
                c = d[(d.truth_set == label) & (d.series == SER_OFF) & (d.noise_level == lv)]
                if a.empty or c.empty:
                    continue
                a, c = a.iloc[0], c.iloc[0]

                def cell(r):
                    if r.n_reps == 1:
                        return "%+5.2f" % r.slope_median, "yes" if r.share_significant > 0 else "no"
                    return ("%+5.2f [%+5.2f, %+5.2f]" % (r.slope_median, r.slope_p10, r.slope_p90), "%3.0f%%" % (100 * r.share_significant))
                ca, sa = cell(a)
                cc_, sc = cell(c)
                terr = "x%.2f (x%.2f)" % (a.typical_error_median, a.null_typical_error_median)
                print("%-18s %-6s | %-24s %-11s | %-20s | %-24s %-11s"
                      % (label if lv == levels[0] else "", level_label(lv), ca, sa, terr, cc_, sc))
            print()
    print("  'none' = the exact truth, as in analyze_cc_recovery_110.py. Other rows: the median over the noisy replicates, with the "
          "middle 80% of the replicates in brackets.")
    print("  significant: with no noise a yes or no; with noise, the share of replicates in which shuffling cc rarely (p < 0.05) gives a "
          "slope this large.")
    print("  typ. error: the typical factor by which the best runs' median misses the true cc (x1.00 = perfect), and in brackets what a "
          "no-information fit gives.")
    print("  The control should keep a slope near 0 at every noise level and a low 'significant' share; if not, noise alone creates a pattern.")
    print("  Replicates reuse the same 250 runs, so the control share is not a 5% calibration. Noise: observed = truth x (1 + b) x exp(e), "
          "b one error per hydrograph (sd), e per 5-minute step (sd / 2); each truth gets its own independent draw.")


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def fig_slope(summ, fracs, levels, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    labels = list(summ.truth_set.unique())
    if not labels:
        return False
    nr, nc = len(labels), len(fracs)
    fig, axes = plt.subplots(nr, nc, figsize=(4.6 * nc + 0.4, 3.9 * nr + 1.1), facecolor=T.SURFACE, sharey=True, squeeze=False)
    xs = np.arange(len(levels))
    ymin = min(-0.3, float(np.nanmin(summ.slope_p10)) - 0.1)
    ymax = max(1.1, float(np.nanmax(summ.slope_p90)) + 0.1)
    for ri, label in enumerate(labels):
        for ci, frac in enumerate(fracs):
            ax = axes[ri][ci]
            T._style(ax)
            ax.axhline(1.0, color=T.INK, linestyle="--", linewidth=1.0, zorder=2)
            ax.axhline(0.0, color=T.COL_NULL, linestyle=":", linewidth=1.4, zorder=2)
            for s, col, off, mk, filled in ((SER_OFF, T.COL_NULL, +0.07, "s", False), (SER_ON, T.COL_CC, -0.07, "o", True)):
                g = summ[(summ.truth_set == label) & (summ.top_frac == frac) & (summ.series == s)].set_index("noise_level")
                lv_ok = [i for i, lv in enumerate(levels) if lv in g.index]
                x = xs[lv_ok] + off
                med = np.array([g.loc[levels[i], "slope_median"] for i in lv_ok])
                lo = np.array([g.loc[levels[i], "slope_p10"] for i in lv_ok])
                hi = np.array([g.loc[levels[i], "slope_p90"] for i in lv_ok])
                ax.vlines(x, lo, hi, color=col, linewidth=2.2 if s == SER_ON else 1.5, zorder=3)
                ax.plot(x, med, linestyle="-", color=col, linewidth=1.0, alpha=0.6, zorder=3)
                ax.plot(x, med, linestyle="none", marker=mk, markersize=6, markerfacecolor=col if filled else T.SURFACE,
                        markeredgecolor=col, markeredgewidth=1.5, zorder=4)
            ax.set_xticks(xs)
            ax.set_xticklabels([level_label(lv) for lv in levels], fontsize=8.5)
            ax.set_xlim(-0.5, len(levels) - 0.5)
            ax.set_ylim(ymin, ymax)
            nt = int(summ[(summ.truth_set == label) & (summ.top_frac == frac)].n_truths.iloc[0])
            ax.set_title("Best %d%% of runs, %s (%d)" % (round(100 * frac), label, nt), color=T.INK, fontsize=9, loc="left")
            if ri == nr - 1:
                ax.set_xlabel("observation error (sd)", color=T.INK2, fontsize=9)
            if ci == 0:
                ax.set_ylabel("tracking slope", color=T.INK2, fontsize=9)
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=6, linewidth=2.2, label="channel loss ON: median over replicates, bar = middle 80% of replicates"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=6, markerfacecolor=T.SURFACE, markeredgewidth=1.5, linewidth=1.5,
               label="control (cc inert): the same"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=1.0, label="1 = best runs follow the true cc one-for-one"),
        Line2D([0], [0], color=T.COL_NULL, linestyle=":", linewidth=1.4, label="0 = best runs ignore cc"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Do the best runs still follow the true cc when the observed hydrograph has error?", color=T.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.09 if nr == 1 else 0.07, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


def fig_truth(pt, frac, levels, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    d = pt[pt.top_frac == frac]
    if d.empty:
        return False
    tv = sorted(d.truth_cc_mmhr.unique())
    nl = len(levels)
    fig, axes = plt.subplots(1, nl, figsize=(4.0 * nl + 0.5, 5.6), facecolor=T.SURFACE, sharey=True, squeeze=False)
    lo, hi = 10 ** R.LOG_LO, 10 ** R.LOG_HI
    lim = (lo * 0.75, hi * 1.3)
    for ci, lv in enumerate(levels):
        ax = axes[0][ci]
        T._style(ax)
        ax.axhspan(10 ** (R.LOG_LO + 0.1 * (R.LOG_HI - R.LOG_LO)), 10 ** (R.LOG_LO + 0.9 * (R.LOG_HI - R.LOG_LO)),
                   color=T.COL_NULL, alpha=0.16, linewidth=0, zorder=1)
        ax.axhline(10 ** R.CENTER_LOG, color=T.COL_NULL, linestyle=":", linewidth=1.3, zorder=2)
        ax.plot(lim, lim, color=T.INK, linestyle="--", linewidth=1.0, zorder=2)
        for s, col, off, mk, filled in ((SER_OFF, T.COL_NULL, +1, "s", False), (SER_ON, T.COL_CC, -1, "o", True)):
            g = d[(d.series == s) & (d.noise_level == lv)].sort_values("truth_cc_mmhr")
            if g.empty:
                continue
            x = g.truth_cc_mmhr.values * 10 ** (off * 0.012)
            ax.vlines(x, g.cc_p10, g.cc_p90, color=col, linewidth=2.0 if s == SER_ON else 1.4, zorder=3)
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
        nrep = int(d[(d.series == SER_ON) & (d.noise_level == lv)].n_reps.iloc[0])
        ax.set_title("noise: %s" % level_label(lv) if lv == 0 else "noise sd %s (%d replicates)" % (level_label(lv), nrep),
                     color=T.INK, fontsize=9.5, loc="left")
        ax.set_xlabel("true cc of the truth (mm/hr)", color=T.INK2, fontsize=9)
        if ci == 0:
            ax.set_ylabel("median cc of the best runs (mm/hr)", color=T.INK2, fontsize=9)
    handles = [
        Line2D([0], [0], color=T.COL_CC, marker="o", markersize=5.5, linewidth=2.0,
               label="channel loss ON: median over replicates, bar = middle 80% of replicates"),
        Line2D([0], [0], color=T.COL_NULL, marker="s", markersize=5.5, markerfacecolor=T.SURFACE, markeredgewidth=1.5, linewidth=1.4,
               label="control (cc inert): the same"),
        Line2D([0], [0], color=T.INK, linestyle="--", linewidth=1.0, label="perfect recovery (median sits on the true cc)"),
        Patch(facecolor=T.COL_NULL, alpha=0.16, label="no information: middle 80% of the sampled range (dotted = its middle)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8, labelcolor=T.INK2)
    fig.suptitle("Best %d%% of runs: where does the median cc land when the observed hydrograph has error?  (no-noise panel: one value, no bar)"
                 % round(100 * frac), color=T.INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.11, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=T.SURFACE)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110 -- does the tracking of the true cc survive an imperfect observed hydrograph?")
    ap.add_argument("--top_fracs", type=float, nargs="+", default=[0.20, 0.10, 0.05],
                    help="Top fractions by KGE_2012; the first is the main one (default 0.20 0.10 0.05)")
    ap.add_argument("--noise_levels", type=float, nargs="+", default=[0.05, 0.10, 0.20],
                    help="Noise sd levels as fractions (default 0.05 0.10 0.20)")
    ap.add_argument("--noise_reps", type=int, default=100, help="Noisy replicates per level (default 100)")
    ap.add_argument("--noise_parts", choices=["both", "volume", "shape"], default="both",
                    help="both = whole-hydrograph error b and per-step scatter e (default); volume = only b; shape = only e")
    ap.add_argument("--perm", type=int, default=1000, help="Permutation shuffles per replicate (default 1000)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--expect_n", type=int, default=250)
    ap.add_argument("--truths", type=float, nargs="+", default=None,
                    help="Only these cc-ON truths (the cc-OFF reference is always kept for the checks)")
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--continue_on_fail", action="store_true")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()
    if any(lv <= 0 for lv in args.noise_levels):
        sys.exit("--noise_levels must be positive fractions, e.g. 0.05 0.10 0.20 (the no-noise case is always included).")
    if args.noise_reps < 2:
        sys.exit("--noise_reps must be at least 2.")
    if args.perm < 200:
        sys.exit("--perm must be at least 200.")

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    out_dir = summary_dir / OUT_DIRNAME
    tag = "" if args.noise_parts == "both" else "_" + args.noise_parts
    print("\n" + "=" * 78 + "\nSeries 110 -- does the tracking of the true cc survive observation error?\n" + "=" * 78)
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
    data = T.build_tensors(df, truths)
    sets = R.tracking_sets(truths)
    if not sets:
        sys.exit("STOPPED: the tracking slope needs at least %d cc-ON truths." % R.MIN_TRACK_TRUTHS)
    n_runs = data[SER_ON]["X"].shape[0]
    n_by_frac = {f: T.n_top_for(n_runs, f) for f in fracs}
    levels = [0.0] + list(args.noise_levels)
    print("\n%d cc-ON truth(s), %d runs per series, cuts %s, noise sd %s, %d replicates per level, %d shuffles per replicate, noise = %s"
          % (sum(1 for t in truths if t > 0), n_runs, ["%d%%" % round(100 * f) for f in fracs],
             ["%g%%" % round(100 * lv, 6) for lv in args.noise_levels], args.noise_reps, args.perm, args.noise_parts))

    print("\nReading the stored hydrographs ...")
    sims, truth_obs, worst = load_hydrographs(summary_dir, data, truths)
    print("  hydrographs: %d + %d runs x %d five-minute steps" % (sims[SER_ON].shape[0], sims[SER_OFF].shape[0], sims[SER_ON].shape[1]))
    print("  noise-free check: KGE_2012 recomputed here equals the stored values for every (truth, run); worst difference %.1e" % worst)

    print("\nAdding noise and re-picking the best runs ...")
    rz, tz = run_replicates(args, data, truths, fracs, sets, sims, truth_obs)
    summ = summarize(rz)
    ptab = per_truth_table(tz)
    print_summary(summ, fracs, levels, n_by_frac)
    print("\n  " + compare_with_saved(summ, summary_dir))

    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ck.rows).to_csv(out_dir / "recovery_noise_input_checks_110.csv", index=False)
    rz.to_csv(out_dir / ("recovery_noise_realizations_110%s.csv" % tag), index=False)
    summ.to_csv(out_dir / ("recovery_noise_summary_110%s.csv" % tag), index=False)
    ptab.to_csv(out_dir / ("recovery_noise_per_truth_110%s.csv" % tag), index=False)
    if not args.no_plots:
        try:
            fig_slope(summ, fracs, levels, out_dir / ("fig_recovery_noise_slope_110%s.png" % tag))
            fig_truth(ptab, fracs[0], levels, out_dir / ("fig_recovery_noise_truth_110%s.png" % tag))
        except Exception as e:
            print("\n  (figures skipped: %s -- the CSVs are saved regardless)" % e)

    def md5_of(p):
        return hashlib.md5(Path(p).read_bytes()).hexdigest()
    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "pandas": pd.__version__, "numpy": np.__version__,
            "inputs": {n: md5_of(summary_dir / n) for n in (P.STAGE2_NAME, P.CONTROL_NAME)},
            "noise_free_check_worst_diff": worst,
            "settings": {k: (v if not isinstance(v, Path) else str(v)) for k, v in vars(args).items()},
            "truths": [float(t) for t in truths], "checks": ck.rows}
    (out_dir / ("PROVENANCE_recovery_noise_cc_110%s.json" % tag)).write_text(json.dumps(prov, indent=2, default=str))

    print("\nSaved to: %s" % out_dir)
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results." % ck.n("FAIL"))
    print("  Caveat: one storm, %d runs, routing pinned at truth. Multiplicative noise with no memory in time; real gauge error can be larger "
          "at high flows and correlated in time." % n_runs)


if __name__ == "__main__":
    main()
