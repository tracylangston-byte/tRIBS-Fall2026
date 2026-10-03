"""
analyze_cc_tests_110.py
=======================
Series 110 -- statistical tests that go with analyze_cc_pca_110.py.

The PCA script answered "how tightly do the best runs pin each parameter?".
This script answers the follow-up questions the PCA could not:

  1. Is the cc signal real, or could it be chance?          -> PERMUTATION TEST
  2. How much do the trends across truths wobble, given that
     all 12 truths re-score the SAME 250 runs?              -> JOINT BOOTSTRAP
  3. How does the best fit change as cc (or Ks, f) is moved
     away from its best value?                              -> PROFILE CURVES
  4. (optional, --noise) Do the conclusions survive an
     imperfect "observed" hydrograph?                       -> NOISE TEST

Pure post-hoc analysis. Runs no tRIBS. Writes only into its own new folder
(calibration_work/03_comparisons/summary_tables/tests_cc_110/).

Keep this file in the same folder as analyze_cc_pca_110.py (it reuses that
script's loader and input checks). --noise also needs rescore_cc_truths_110.py
in the same folder.

VOCABULARY
-----------
  Bootstrap   Re-draw the runs (with replacement) many times to see how much a
      number wobbles from sample to sample. The middle 95% of the wobble is the
      "95% bootstrap interval". "Joint" means each re-draw is applied to all 12
      truths at once, because the truths share the same runs.
  Null hypothesis (H0)   The boring explanation we try to rule out. Here: "cc has
      no effect on how well a run fits; any pattern with cc is chance".
  Permutation test   Shuffle the cc values among the runs (keeping each run's Ks,
      f and fit as they are). A shuffled data set is what the world would look like
      if H0 were true. Do it thousands of times to get the "null distribution".
  p-value   The share of shuffled data sets that look at least as cc-informative as
      the real one. Small p (say < 0.05) = hard to explain by chance. Large p means
      "cannot rule out chance", NOT "cc has no effect".
  Profile curve   For each slice of one parameter (a "bin"), the best fit any run
      in that slice achieves. A flat curve = the data cannot tell those values
      apart. A peaked curve = the data can.
  Detection rate   (noise test) In what share of noisy "observed" hydrographs does
      the cc signal still come out significant (p < 0.05). Reported for the
      trade-off test (|r(Ks,cc)|) and for the pinning test (shrink_cc).

WHY A PERMUTATION TEST FITS THIS DESIGN
-----------------------------------------
The 250 (Ks, f, cc) points come from a Latin hypercube, so cc was assigned to
runs independently of Ks and f. If cc did not matter, the fit of a run would
depend only on its Ks and f, and which cc it happened to get would be arbitrary.
Shuffling cc across runs therefore reproduces H0 exactly (up to tiny chance
correlations in the design). The cc-OFF control series is the built-in check: cc
has no effect there, so its p-values should look random (about 5% below 0.05).

THE TESTS (all use the best top_frac of runs by KGE_2012, as in the PCA script)
---------------------------------------------------------------------------------
  shrink_cc      spread of cc among the best runs / spread of the sampled box.
                 Small = cc pinned.   (test: p = share of shuffles this small)
  |r(Ks,cc)|     rank correlation between Ks and cc among the best runs.
                 Large = Ks and cc offset each other.
  cc weight      how much cc takes part in the tightest parameter combination
                 (the last principal component). Large = cc is part of what the
                 data pin down.
  profile range  best-in-bin KGE: highest bin minus lowest bin along cc.
  bins near best number of cc bins that reach within --profile_tol of the best KGE.
  combined       the mean of shrink_cc, of |r(Ks,cc)| and of cc weight over all
                 cc-ON truths, with ONE shuffle shared by every truth. Each is a
                 single test ("is cc ever informative?"), so none suffers from 12
                 separate p-values.
12 truths x several tests give many p-values that are NOT independent (same runs).
Read the combined p first; treat single-truth p-values as descriptive.

WHICH TEST TO TRUST (from checks on stand-in data with a known cc effect)
---------------------------------------------------------------------------
cc often shows up as a TRADE-OFF with Ks (a different Ks can make up for a different
cc), not as a narrow cc range among the best runs. So |r(Ks,cc)| and cc weight detect
a cc effect much more often than shrink_cc or the profile curves do. A significant
|r(Ks,cc)| means cc is INFORMATIVE (the fit is sensitive to it); it does not mean cc
is individually identifiable. Only a small shrink_cc (cc pinned on its own) speaks to
identifiability. With cc inert, all of these tests came out below p = 0.05 about 5% of
the time, as they should.

JOINT BOOTSTRAP CONTRASTS (cc-ON truths only; control series drawn independently)
----------------------------------------------------------------------------------
  cc gap         mean over truths of (control shrink_cc - cc-ON shrink_cc).
                 Positive = cc more pinned than when it has no effect.
  slope / decade trend of a statistic with log10(true cc): Ks shrinkage,
                 r(Ks,cc), cc weight. No cut-point to choose.
  hi minus lo    mean for truths >= --hi_min minus mean for truths <= --lo_max.
                 EXPLORATORY: those split points were picked after the pattern was
                 seen, so treat these intervals as descriptive.
  Each contrast is also given for "interior" truths only (not within 25% of the
  cc box edge, where the box clips the ridge).

OUTPUTS   (calibration_work/03_comparisons/summary_tables/tests_cc_110/)
--------------------------------------------------------------------------
  tests_permutation_110.csv         per truth x series x top_frac: observed values and p-values
  tests_combined_110.csv            combined permutation tests
  tests_bootstrap_truth_110.csv     per truth x series x top_frac: estimates and joint-bootstrap 95% CIs
  tests_bootstrap_contrasts_110.csv across-truth contrasts with CIs
  tests_profiles_110.csv            profile curves (Ks, f, cc), every truth and series
  tests_profile_tests_110.csv       profile summaries and cc permutation p-values
  fig_tests_profiles_110.png  fig_tests_bootstrap_110.png  fig_tests_permutation_110.png
  with --noise:  tests_noise_summary_110.csv  tests_noise_realizations_110.csv  fig_tests_noise_110.png
  PROVENANCE_tests_cc_110.json

USAGE (run from the lab/ directory)
------------------------------------
    python analyze_cc_tests_110.py                    # tests 1-3, about a minute
    python analyze_cc_tests_110.py --noise            # adds the noise test
    python analyze_cc_tests_110.py --top_fracs 0.2 0.1 0.05 --boot 2000 --perm 10000
    python analyze_cc_tests_110.py --figure_truths 157 425
    python analyze_cc_tests_110.py --truths 157 425 --no_plots

CAVEATS
--------
One storm, 250 points in 3-D, routing pinned at truth: a best case. A permutation
test asks only "is cc informative at all?"; it does not say how well cc can be
estimated. Hi-minus-lo contrasts are exploratory. A non-significant result with 250
runs is not evidence that cc is unidentifiable: it may mean too few runs near the
best fit.
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
    import analyze_cc_pca_110 as P
except ImportError:
    sys.exit("This script needs analyze_cc_pca_110.py in the same folder (it reuses its loader and "
             "input checks). Put both files in lab/ and run from there.")

_NEEDED = ["BOX", "CONTROL_NAME", "STAGE2_NAME", "C_CC", "GRID", "INK", "INK2", "MIN_TOP", "NEAR_EDGE", "NULLC", "OFF_TRUTH",
           "PRIOR_STD", "RESCORE_DIR", "RESCORE_PROV", "SER_OFF", "SER_ON", "SURFACE", "load_all", "pick_figure_truths",
           "print_checks", "run_checks", "to_X", "truth_label"]
_missing = [n for n in _NEEDED if not hasattr(P, n)]
if _missing:
    sys.exit("analyze_cc_pca_110.py in this folder is an older version (missing: %s). Replace it with the current copy." % ", ".join(_missing))

OUT_DIRNAME = "tests_cc_110"
SER_ON, SER_OFF = P.SER_ON, P.SER_OFF
OFF_TRUTH = P.OFF_TRUTH
STAT_NAMES = ["shrink_Ks", "shrink_f", "shrink_cc", "rho_Ks_f", "rho_Ks_cc", "rho_f_cc", "eig3", "w_cc"]
TRUTH_KS, TRUTH_F = 7.0, 0.012
BIN_PARAMS = [("Ks", 0, "lin"), ("f", 1, "log"), ("cc", 2, "log")]
PARAM_BOX = [(4.0, 10.5), (0.004, 0.030), (30.0, 1000.0)]
COL_CC = P.C_CC
COL_NULL = P.NULLC
INK, INK2, GRID, SURFACE = P.INK, P.INK2, P.GRID, P.SURFACE


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def n_top_for(n_finite, frac):
    return max(P.MIN_TOP, int(round(frac * n_finite)))


def box_pos(t):
    lo, hi = np.log10(P.BOX["channelconductivity_mmhr"][0]), np.log10(P.BOX["channelconductivity_mmhr"][1])
    return (np.log10(t) - lo) / (hi - lo)


def _ranks(Xb):
    try:
        return stats.rankdata(Xb, axis=1)
    except TypeError:                                   # very old scipy
        return np.apply_along_axis(stats.rankdata, 1, Xb)


def _rho_cols(R, i, j):
    with np.errstate(all="ignore"):
        a = R[:, :, i] - R[:, :, i].mean(axis=1, keepdims=True)
        b = R[:, :, j] - R[:, :, j].mean(axis=1, keepdims=True)
        return (a * b).sum(axis=1) / np.sqrt((a ** 2).sum(axis=1) * (b ** 2).sum(axis=1))


def top_stats_batch(Xb):
    """Xb: (B, m, 3) parameters (Ks, log10 f, log10 cc) of the top set, one row of m runs
    per resample or shuffle. Returns a dict of (B,) arrays (STAT_NAMES)."""
    m = Xb.shape[1]
    Z = (Xb - Xb.mean(axis=1, keepdims=True)) / P.PRIOR_STD
    C = np.einsum("bmi,bmj->bij", Z, Z) / (m - 1)
    shrink = np.sqrt(np.clip(np.einsum("bii->bi", C), 0, None))
    evals, evecs = np.linalg.eigh(C)                    # ascending: column 0 = tightest combination
    R = _ranks(Xb)
    return {
        "shrink_Ks": shrink[:, 0], "shrink_f": shrink[:, 1], "shrink_cc": shrink[:, 2],
        "rho_Ks_f": _rho_cols(R, 0, 1), "rho_Ks_cc": _rho_cols(R, 0, 2), "rho_f_cc": _rho_cols(R, 1, 2),
        "eig3": evals[:, 0], "w_cc": np.abs(evecs[:, 2, 0]),
    }


def point_top(X, k, frac):
    """Top set (run indices) for the real data and its statistics."""
    ok = np.flatnonzero(np.isfinite(k))
    m = n_top_for(len(ok), frac)
    order = ok[np.argsort(-k[ok], kind="stable")][:m]
    st = top_stats_batch(X[order][None])
    return order, {s: float(v[0]) for s, v in st.items()}


# ------------------------------------------------------------------
# Data tensors
# ------------------------------------------------------------------
def build_tensors(df, truths):
    """Align runs across truths by run_id (the input checks guarantee identical run sets)."""
    out = {}
    for s in (SER_ON, SER_OFF):
        base = df[(df.truth_cc_mmhr == truths[0]) & (df.series == s)].reset_index(drop=True)
        order = base["run_id"].tolist()
        K = np.full((len(truths), len(order)), np.nan)
        for ti, t in enumerate(truths):
            g = df[(df.truth_cc_mmhr == t) & (df.series == s)].set_index("run_id")["kge_2012"]
            K[ti] = g.reindex(order).to_numpy(float)
        out[s] = {"run_id": order, "X": P.to_X(base), "K": K}
    return out


# ------------------------------------------------------------------
# 1. Permutation test
# ------------------------------------------------------------------
def make_perm_matrix(n, Bp, rng):
    return rng.permuted(np.tile(np.arange(n), (Bp, 1)), axis=1)


def perm_null_stats(X, top_idx, Pm):
    """Null statistics when cc is shuffled across ALL runs and the real top set is kept.
    Run i receives the cc of run Pm[b, i]."""
    cc = X[:, 2]
    m = len(top_idx)
    Bp = Pm.shape[0]
    Xb = np.empty((Bp, m, 3))
    Xb[:, :, 0] = X[top_idx, 0][None, :]
    Xb[:, :, 1] = X[top_idx, 1][None, :]
    Xb[:, :, 2] = cc[Pm[:, top_idx]]
    return top_stats_batch(Xb)


def p_low(null, obs):
    return (1.0 + np.sum(null <= obs + 1e-12)) / (len(null) + 1.0)


def p_high(null, obs):
    return (1.0 + np.sum(null >= obs - 1e-12)) / (len(null) + 1.0)


def bin_index(vals, scale, lo, hi, nb):
    v = np.log10(vals) if scale == "log" else np.asarray(vals, float)
    a, b = (np.log10(lo), np.log10(hi)) if scale == "log" else (lo, hi)
    return np.clip(np.floor((v - a) / (b - a) * nb).astype(int), 0, nb - 1)


def best_per_bin(K, bins, nb):
    out = np.full(nb, np.nan)
    for k in range(nb):
        v = K[(bins == k) & np.isfinite(K)]
        if len(v):
            out[k] = v.max()
    return out


def run_permutation(data, truths, fracs, Bp, nb, tol, rng_seed):
    """Returns (rows, combined_rows, profile_test_rows, null_store)."""
    rows, comb_rows, prof_rows, null_store = [], [], [], {}
    T = len(truths)
    Ton = [i for i, t in enumerate(truths) if t > 0]
    interior = [i for i in Ton if P.NEAR_EDGE <= box_pos(truths[i]) <= 1 - P.NEAR_EDGE]
    for si, s in enumerate((SER_ON, SER_OFF)):
        X, K = data[s]["X"], data[s]["K"]
        n = X.shape[0]
        rng = np.random.default_rng([rng_seed, 11, si])
        Pm = make_perm_matrix(n, Bp, rng)
        bins_pool = bin_index(10 ** X[:, 2], "log", *PARAM_BOX[2], nb)
        for frac in fracs:
            obs_all = {st: np.full(T, np.nan) for st in STAT_NAMES}
            null_all = {st: np.full((Bp, T), np.nan) for st in STAT_NAMES}
            for ti, t in enumerate(truths):
                top, obs = point_top(X, K[ti], frac)
                null = perm_null_stats(X, top, Pm)
                for st in STAT_NAMES:
                    obs_all[st][ti] = obs[st]
                    null_all[st][:, ti] = null[st]
                rows.append({
                    "truth_cc_mmhr": t, "series": s, "top_frac": frac, "n_top": len(top),
                    "shrink_cc": obs["shrink_cc"], "null_mean_shrink_cc": float(null["shrink_cc"].mean()),
                    "null_p05_shrink_cc": float(np.percentile(null["shrink_cc"], 5)),
                    "p_shrink_cc": p_low(null["shrink_cc"], obs["shrink_cc"]),
                    "rho_Ks_cc": obs["rho_Ks_cc"],
                    "p_abs_rho_Ks_cc": p_high(np.abs(null["rho_Ks_cc"]), abs(obs["rho_Ks_cc"])),
                    "rho_f_cc": obs["rho_f_cc"],
                    "p_abs_rho_f_cc": p_high(np.abs(null["rho_f_cc"]), abs(obs["rho_f_cc"])),
                    "w_cc": obs["w_cc"], "null_mean_w_cc": float(null["w_cc"].mean()),
                    "p_w_cc": p_high(null["w_cc"], obs["w_cc"]),
                    "box_pos_cc": box_pos(t) if t > 0 else np.nan,
                })
            null_store[(s, frac)] = null_all
            # combined tests: one shuffle shared by every cc-ON truth
            for label, idx in (("all cc-ON truths", Ton), ("interior truths", interior)):
                if len(idx) == 0:
                    continue
                for st, side in (("shrink_cc", "low"), ("rho_Ks_cc", "abs"), ("w_cc", "high")):
                    if st == "rho_Ks_cc":
                        o = float(np.mean(np.abs(obs_all[st][idx])))
                        nl = np.mean(np.abs(null_all[st][:, idx]), axis=1)
                        p = p_high(nl, o)
                    else:
                        o = float(np.mean(obs_all[st][idx]))
                        nl = np.mean(null_all[st][:, idx], axis=1)
                        p = p_low(nl, o) if side == "low" else p_high(nl, o)
                    comb_rows.append({"series": s, "top_frac": frac, "truth_set": label, "n_truths": len(idx),
                                      "statistic": ("mean " + ("|" + st + "|" if st == "rho_Ks_cc" else st)),
                                      "observed": o, "null_mean": float(nl.mean()),
                                      "null_p05": float(np.percentile(nl, 5)),
                                      "null_p95": float(np.percentile(nl, 95)), "p_value": p})
        # profile tests along cc (independent of top_frac)
        empty_bins = np.array([not np.any(bins_pool == k) for k in range(nb)])
        for ti, t in enumerate(truths):
            Kt = np.where(np.isfinite(K[ti]), K[ti], -np.inf)
            obs_best = best_per_bin(K[ti], bins_pool, nb)
            gbest = np.nanmax(obs_best)
            obs_range = float(np.nanmax(obs_best) - np.nanmin(obs_best))
            obs_near = int(np.sum(obs_best >= gbest - tol))
            bins_perm = bins_pool[Pm]
            bb = np.empty((Bp, nb))
            for k in range(nb):
                bb[:, k] = np.where(bins_perm == k, Kt[None, :], -np.inf).max(axis=1)
            if empty_bins.any():          # a bin with no runs has no profile value, in the data or in the null
                bb = bb[:, ~empty_bins]
            null_range = bb.max(axis=1) - bb.min(axis=1)
            null_near = (bb >= bb.max(axis=1, keepdims=True) - tol).sum(axis=1)
            prof_rows.append({
                "truth_cc_mmhr": t, "series": s, "profile_tol": tol, "n_bins": nb,
                "global_best_kge": float(gbest),
                "profile_range_cc": obs_range, "null_mean_range": float(null_range.mean()),
                "p_profile_range": p_high(null_range, obs_range),
                "bins_near_best_cc": obs_near, "null_mean_bins_near": float(null_near.mean()),
                "p_bins_near_best": p_low(null_near.astype(float), float(obs_near)),
            })
    return rows, comb_rows, prof_rows, null_store


# ------------------------------------------------------------------
# 2. Joint bootstrap
# ------------------------------------------------------------------
def run_bootstrap(data, truths, fracs, B, seed):
    """res[(series, frac)][stat] = (B, T). The same resampled run indices are used for every truth
    (and every top_frac) within a series; ON and control are drawn independently."""
    T = len(truths)
    res = {}
    for si, s in enumerate((SER_ON, SER_OFF)):
        X, K = data[s]["X"], data[s]["K"]
        n = X.shape[0]
        rng = np.random.default_rng([seed, 22, si])
        idx = rng.integers(0, n, (B, n))
        for frac in fracs:
            store = {st: np.full((B, T), np.nan) for st in STAT_NAMES}
            for ti in range(T):
                m = n_top_for(int(np.isfinite(K[ti]).sum()), frac)
                Kfill = np.where(np.isfinite(K[ti]), K[ti], -np.inf)
                Kb = Kfill[idx]
                order = np.argsort(-Kb, axis=1, kind="stable")[:, :m]
                sel = np.take_along_axis(idx, order, axis=1)
                st = top_stats_batch(X[sel])
                for name in STAT_NAMES:
                    store[name][:, ti] = st[name]
            res[(s, frac)] = store
    return res


def point_arrays(data, truths, fracs):
    out = {}
    T = len(truths)
    for s in (SER_ON, SER_OFF):
        X, K = data[s]["X"], data[s]["K"]
        for frac in fracs:
            store = {st: np.full((1, T), np.nan) for st in STAT_NAMES}
            for ti in range(T):
                _, obs = point_top(X, K[ti], frac)
                for st in STAT_NAMES:
                    store[st][0, ti] = obs[st]
            out[(s, frac)] = store
    return out


def contrast_values(on, off, truths, lo_max, hi_min):
    """on / off: dict stat -> (B, T). Returns {name: (B,)} for 'all' and 'interior' truth sets."""
    t_arr = np.asarray(truths, float)
    Ton = np.flatnonzero(t_arr > 0)
    interior = np.array([i for i in Ton if P.NEAR_EDGE <= box_pos(t_arr[i]) <= 1 - P.NEAR_EDGE], dtype=int)
    out = {}
    for tag, idx in (("all", Ton), ("interior", interior)):
        if len(idx) < 3:
            continue
        x = np.log10(t_arr[idx])
        xc = x - x.mean()
        out["cc_gap_mean|" + tag] = np.mean(off["shrink_cc"][:, idx] - on["shrink_cc"][:, idx], axis=1)
        for st in ("shrink_Ks", "rho_Ks_cc", "w_cc"):
            y = on[st][:, idx]
            out[st + "_slope_per_decade|" + tag] = ((y - y.mean(axis=1, keepdims=True)) * xc).sum(axis=1) / (xc ** 2).sum()
            lo = [j for j, i in enumerate(idx) if t_arr[i] <= lo_max]
            hi = [j for j, i in enumerate(idx) if t_arr[i] >= hi_min]
            if lo and hi:
                out[st + "_hi_minus_lo|" + tag] = y[:, hi].mean(axis=1) - y[:, lo].mean(axis=1)
    return out


# ------------------------------------------------------------------
# 3. Profile curves
# ------------------------------------------------------------------
def run_profiles(data, truths, nb):
    rows = []
    for s in (SER_ON, SER_OFF):
        X, K = data[s]["X"], data[s]["K"]
        for pname, col, scale in BIN_PARAMS:
            lo, hi = PARAM_BOX[col]
            vals = 10 ** X[:, col] if scale == "log" else X[:, col]
            bins = bin_index(vals, scale, lo, hi, nb)
            edges = (np.logspace(np.log10(lo), np.log10(hi), nb + 1) if scale == "log" else np.linspace(lo, hi, nb + 1))
            mid = np.sqrt(edges[:-1] * edges[1:]) if scale == "log" else 0.5 * (edges[:-1] + edges[1:])
            for ti, t in enumerate(truths):
                k = K[ti]
                for b in range(nb):
                    v = np.sort(k[(bins == b) & np.isfinite(k)])[::-1]
                    rows.append({"truth_cc_mmhr": t, "series": s, "parameter": pname, "bin": b,
                                 "bin_lo": edges[b], "bin_hi": edges[b + 1], "bin_mid": mid[b], "n_runs": len(v),
                                 "best_kge": float(v[0]) if len(v) else np.nan,
                                 "mean_top3_kge": float(v[:3].mean()) if len(v) else np.nan,
                                 "median_kge": float(np.median(v)) if len(v) else np.nan})
    return rows


# ------------------------------------------------------------------
# 4. Noise test (optional)
# ------------------------------------------------------------------
def kge_matrix(S, obs):
    """KGE_2012 of every row of S (runs x time) against one observed series (time,).
    Same formula as the scorer / rescore script: gamma = alpha / beta."""
    with np.errstate(all="ignore"):
        so = S - S.mean(axis=1, keepdims=True)
        oo = obs - obs.mean()
        r = (so @ oo) / np.sqrt((so ** 2).sum(axis=1) * (oo ** 2).sum())
        alpha = S.std(axis=1) / obs.std()
        beta = S.mean(axis=1) / obs.mean()
        gamma = alpha / beta
        return 1 - np.sqrt((r - 1) ** 2 + (gamma - 1) ** 2 + (beta - 1) ** 2)


def spearman_abs_null(m, Bn, rng):
    """Null for |Spearman r| between Ks and cc inside a top set of m runs when cc has no effect: the cc ranks are then a
    uniformly random ordering, whatever Ks values the top set has. r = 1 - 6 sum(d^2) / (m (m^2 - 1))."""
    perm = np.argsort(rng.random((Bn, m)), axis=1)
    d = perm - np.arange(m)[None, :]
    return np.abs(1.0 - 6.0 * (d ** 2).sum(axis=1) / (m * (m ** 2 - 1.0)))


def subset_std_null(cc_log, m, Bn, rng):
    """Spread of cc (scaled to the box) in a random subset of m runs: the permutation null for shrink_cc."""
    n = len(cc_log)
    sub = np.argsort(rng.random((Bn, n)), axis=1)[:, :m]
    v = cc_log[sub]
    return v.std(axis=1, ddof=1) / P.PRIOR_STD[2]


def run_noise(args, summary_dir, df, data, truths, fracs):
    try:
        import rescore_cc_truths_110 as RS
    except ImportError:
        sys.exit("--noise needs rescore_cc_truths_110.py in the same folder.")
    calib = summary_dir.parent.parent
    csv_dir = calib / "03_comparisons" / "csv_exports"
    cand_dir = calib / RS.CAND_DIRNAME
    prov_path = summary_dir / P.RESCORE_DIR / P.RESCORE_PROV
    print("\n" + "=" * 100 + "\nNOISE TEST -- re-scoring the stored hydrographs against noisy versions of each truth\n" + "=" * 100)

    # --- hydrographs ---------------------------------------------------------
    sims, idx0 = {}, None
    for s in (SER_ON, SER_OFF):
        res = pd.DataFrame({"run_id": data[s]["run_id"]})
        cache, prob = RS.load_all_compares(res, csv_dir)
        if prob:
            sys.exit("STOPPED: %d compare CSV(s) missing or unreadable, e.g. %s\n  Looked in: %s\n  --noise re-reads every run's stored "
                     "hydrograph (<run_id>_compare_obs_sim.csv). Everything else was already saved." % (len(prob), list(prob.items())[0], csv_dir))
        mats = []
        for rid in data[s]["run_id"]:
            d = cache[rid]
            if idx0 is None:
                idx0 = d.index
            if not d.index.equals(idx0):
                sys.exit("Compare CSV %s has a different time index from the others; refusing." % rid)
            mats.append(d["Simulated"].to_numpy(float))
        sims[s] = np.vstack(mats)
        if s == SER_ON:
            first = cache[data[s]["run_id"][0]]
            obs_off_truth = first["Observed"].to_numpy(float)
    T_len = sims[SER_ON].shape[1]
    print("  hydrographs: %d + %d runs x %d five-minute steps" % (sims[SER_ON].shape[0], sims[SER_OFF].shape[0], T_len))

    # --- truth series, verified against the re-scoring provenance ---------------
    prov = json.loads(prov_path.read_text()) if prov_path.exists() else {}
    cand_prov = prov.get("candidate_truths", {})
    truth_obs = {OFF_TRUTH: obs_off_truth}
    for name, info in cand_prov.items():
        p = cand_dir / name
        if not p.exists():
            sys.exit("Candidate truth file missing: %s" % p)
        if RS.md5_of(p) != info["md5"]:
            sys.exit("Candidate truth %s no longer matches its recorded checksum; refusing." % name)
        s5 = RS.read_truth_5min(p).reindex(idx0)
        if s5.isna().any():
            sys.exit("Truth %s has no value at some compare-CSV timestamps." % name)
        truth_obs[float(info["cc_mmhr"])] = s5.to_numpy(float)
    use = [t for t in truths if t in truth_obs or any(np.isclose(t, k, rtol=1e-6) for k in truth_obs)]
    key = {t: (t if t in truth_obs else [k for k in truth_obs if np.isclose(t, k, rtol=1e-6)][0]) for t in use}

    # --- noise-free equivalence check ---------------------------------------------
    worst = 0.0
    for ti, t in enumerate(truths):
        if t not in key:
            continue
        for s in (SER_ON, SER_OFF):
            k = kge_matrix(sims[s], truth_obs[key[t]])
            stored = data[s]["K"][ti]
            ok = np.isfinite(stored)
            worst = max(worst, float(np.max(np.abs(k[ok] - stored[ok]))))
    if worst > 1e-6:
        sys.exit("STOPPED: re-computing KGE_2012 without noise differs from the stored values by up to %.3g. "
                 "The noise test would not be comparable to the main analysis." % worst)
    print("  noise-free check: KGE_2012 recomputed here equals the stored values for every (truth, run); worst diff %.1e"
          % worst)

    # --- noise model -----------------------------------------------------------------
    rng = np.random.default_rng([args.seed, 33])
    levels = args.noise_levels
    R = args.noise_reps
    Bn = 20000
    cc_log = {s: data[s]["X"][:, 2] for s in (SER_ON, SER_OFF)}
    n = sims[SER_ON].shape[0]
    null_cache, rho_null = {}, {}
    for frac in fracs:
        m = n_top_for(n, frac)
        rho_null[m] = spearman_abs_null(m, Bn, np.random.default_rng([args.seed, 55, m]))
        for s in (SER_ON, SER_OFF):
            null_cache[(s, m)] = subset_std_null(cc_log[s], m, Bn, np.random.default_rng([args.seed, 44, m]))
    rows = []
    for lv in [0.0] + list(levels):
        reps = 1 if lv == 0 else R
        for r in range(reps):
            if lv == 0:
                mult = np.ones(T_len)
            else:
                b = rng.normal(0, lv)
                e = rng.normal(0, lv / 2.0, T_len)
                mult = np.clip(1 + b, 0.05, None) * np.exp(e)
            for ti, t in enumerate(truths):
                if t not in key:
                    continue
                o = truth_obs[key[t]] * mult
                K = {s: kge_matrix(sims[s], o) for s in (SER_ON, SER_OFF)}
                for frac in fracs:
                    rec = {"noise_level": lv, "rep": r, "truth_cc_mmhr": t, "top_frac": frac}
                    for s, tag in ((SER_ON, "on"), (SER_OFF, "off")):
                        top, st = point_top(data[s]["X"], K[s], frac)
                        rec["shrink_cc_" + tag] = st["shrink_cc"]
                        rec["p_shrink_cc_" + tag] = p_low(null_cache[(s, len(top))], st["shrink_cc"])
                        rec["rho_Ks_cc_" + tag] = st["rho_Ks_cc"]
                        rec["p_rho_Ks_cc_" + tag] = p_high(rho_null[len(top)], abs(st["rho_Ks_cc"]))
                        if tag == "on":
                            rec["shrink_Ks_on"] = st["shrink_Ks"]
                            rec["best_kge_on"] = float(np.nanmax(K[s]))
                    rows.append(rec)
    rz = pd.DataFrame(rows)
    rz["detect_on"] = rz["p_shrink_cc_on"] < 0.05
    rz["detect_off"] = rz["p_shrink_cc_off"] < 0.05
    rz["detect_rho_on"] = rz["p_rho_Ks_cc_on"] < 0.05
    rz["detect_rho_off"] = rz["p_rho_Ks_cc_off"] < 0.05
    summ = (rz.groupby(["truth_cc_mmhr", "top_frac", "noise_level"])
              .agg(n_reps=("rep", "size"), median_shrink_cc_on=("shrink_cc_on", "median"),
                   median_shrink_cc_off=("shrink_cc_off", "median"),
                   detection_rate_on=("detect_on", "mean"), false_positive_rate_off=("detect_off", "mean"),
                   detection_rate_rho_on=("detect_rho_on", "mean"), false_positive_rate_rho_off=("detect_rho_off", "mean"),
                   median_shrink_Ks_on=("shrink_Ks_on", "median"), median_rho_Ks_cc_on=("rho_Ks_cc_on", "median"),
                   median_best_kge_on=("best_kge_on", "median")).reset_index())
    return rz, summ


def print_noise(summ, main_frac, levels):
    d = summ[summ.top_frac == main_frac]
    print("\nNOISE TEST at the top %d%% cut. Detection rate = share of noisy 'observed' hydrographs in which the test is "
          "still significant (p < 0.05)." % round(100 * main_frac))
    for title, col, fpcol, note in (
            ("A. trade-off test: is |r(Ks, cc)| in the best runs larger than chance?  (the more powerful test)",
             "detection_rate_rho_on", "false_positive_rate_rho_off", "'no noise' = 'yes' means significant against the exact truth."),
            ("B. pinning test: is cc shrinkage in the best runs smaller than chance?  (the weaker test)",
             "detection_rate_on", "false_positive_rate_off", "")):
        print("\n  " + title)
        hdr = "  %12s | %8s | " % ("truth cc", "no noise") + " ".join("%10s" % ("sd %g%%" % (100 * lv)) for lv in levels)
        print(hdr + " | control (cc inert): share significant")
        for t in sorted(d.truth_cc_mmhr.unique()):
            g = d[d.truth_cc_mmhr == t].set_index("noise_level")
            base = "yes" if g.loc[0.0, col] > 0 else "no"
            cells = " ".join("%10s" % ("%3.0f%%" % (100 * g.loc[lv, col])) if lv in g.index else "%10s" % "n/a" for lv in levels)
            fp = np.mean([g.loc[lv, fpcol] for lv in levels if lv in g.index]) if len(levels) else np.nan
            print("  %12s | %8s | %s | %3.0f%%" % (P.truth_label(t), base, cells, 100 * fp))
        if note:
            print("  " + note)
    print("\n  Noise model: observed = truth x (1 + b) x exp(e); b ~ N(0, sd) is one error for the whole hydrograph "
          "(rating-curve or volume error), e ~ N(0, sd/2) is independent at every 5-minute step.")
    print("  How to read: only rows with 'no noise' = yes say anything about robustness. A rate that stays high as noise grows = "
          "the signal survives an imperfect observed hydrograph; one that collapses = it is fragile.")
    print("  If 'no noise' = no, a rate above 0% is just scatter that the noise adds to which runs make the top set, not a signal.")
    print("  The control column should stay low (noise alone must not create a cc signal). All noisy copies reuse the same runs, so "
          "it is not the 5% calibration of the permutation tests; that was checked separately on fresh random designs.")


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------
def fmt_p(p):
    """p-values from B shuffles can never be exactly 0 (smallest is 1/(B+1)), so tiny ones print as a bound."""
    return "<0.001" if p < 0.001 else "%.3f" % p


def print_permutation(perm_df, comb_df, main_frac, fracs):
    on = perm_df[(perm_df.series == SER_ON) & (perm_df.top_frac == main_frac)].sort_values("truth_cc_mmhr")
    print("\n" + "=" * 118)
    print("1. PERMUTATION TEST (cc-ON runs, best %d%% by KGE_2012).  p = share of cc-shuffled data sets that look at least "
          "this cc-informative" % round(100 * main_frac))
    print("=" * 118)
    print("%12s | %9s %10s %7s | %8s %7s | %8s %7s" % ("truth cc", "shrink_cc", "null mean", "p", "r(Ks,cc)", "p", "cc weight", "p"))
    for _, a in on.iterrows():
        print("%12s | %9.2f %10.2f %7s | %+8.2f %7s | %8.2f %7s"
              % (P.truth_label(a.truth_cc_mmhr), a.shrink_cc, a.null_mean_shrink_cc, fmt_p(a.p_shrink_cc),
                 a.rho_Ks_cc, fmt_p(a.p_abs_rho_Ks_cc), a.w_cc, fmt_p(a.p_w_cc)))
    print("  shrink_cc: small = cc pinned.  r(Ks,cc): p tests |r|.  cc weight: share of cc in the tightest combination.")
    if comb_df.empty:
        print("\n  (combined tests skipped: they need at least one cc-ON truth; only the cc-OFF truth is loaded)")
    else:
        print("\n  COMBINED TESTS (one shuffle shared by all truths; read these before the single-truth p-values)")
        print("  %-8s %-6s %-18s %-26s %9s %9s %7s" % ("series", "cut", "truths", "statistic", "observed", "null mean", "p"))
    for _, c in (comb_df if comb_df.empty else
                 comb_df.sort_values(["series", "top_frac", "truth_set", "statistic"], ascending=[False, False, True, True])).iterrows():
        print("  %-8s %-6s %-18s %-26s %9.3f %9.3f %7s"
              % (c.series.replace("cc ", ""), "%d%%" % round(100 * c.top_frac), c.truth_set, c.statistic, c.observed,
                 c.null_mean, fmt_p(c.p_value)))
    print("\n  HOW TO READ THESE (the three tests ask different questions)")
    print("    shrink_cc small, p small        -> cc is pinned ON ITS OWN: the best runs agree on cc (identifiable).")
    print("    r(Ks,cc) or cc weight p small,  -> cc matters, but only in COMBINATION with Ks (it can be offset by Ks):")
    print("      shrink_cc p not small            sensitive, not identifiable alone. This test is far more powerful than shrink_cc.")
    print("    none small                      -> cc cannot be detected at this sample size; NOT proof that cc is irrelevant.")
    off = perm_df[perm_df.series == SER_OFF]
    if len(off):
        tests = off[["p_shrink_cc", "p_abs_rho_Ks_cc", "p_abs_rho_f_cc", "p_w_cc"]].to_numpy().ravel()
        print("  Calibration: in the cc-OFF CONTROL series (cc has no effect) %d of %d single-truth p-values are below 0.05 "
              "(about 5%% expected; the tests are correlated, so this is a rough guide)." % (int((tests < 0.05).sum()), len(tests)))


def print_profile_tests(prof_df):
    on = prof_df[prof_df.series == SER_ON].sort_values("truth_cc_mmhr")
    off = prof_df[prof_df.series == SER_OFF].set_index("truth_cc_mmhr")
    print("\n" + "=" * 118)
    print("3. PROFILE CURVES along cc (best KGE_2012 in each of %d cc bins).  tolerance for 'near best' = %g"
          % (int(on.n_bins.iloc[0]), on.profile_tol.iloc[0]))
    print("=" * 118)
    print("%12s | %8s %9s %8s %7s | %9s %9s %7s | %s" % ("truth cc", "best KGE", "range", "null", "p", "bins near", "null", "p", "control p: range / bins"))
    for _, a in on.iterrows():
        b = off.loc[a.truth_cc_mmhr] if a.truth_cc_mmhr in off.index else None
        print("%12s | %8.4f %9.4f %8.4f %7s | %9d %9.1f %7s | %s"
              % (P.truth_label(a.truth_cc_mmhr), a.global_best_kge, a.profile_range_cc, a.null_mean_range, fmt_p(a.p_profile_range),
                 a.bins_near_best_cc, a.null_mean_bins_near, fmt_p(a.p_bins_near_best),
                 "%s / %s" % (fmt_p(b.p_profile_range), fmt_p(b.p_bins_near_best)) if b is not None else "n/a"))
    print("  range: best-bin KGE minus worst-bin KGE along cc (big = the fit depends on cc).  bins near: cc bins that still "
          "reach within the tolerance of the best KGE (few = cc is pinned).")
    print("  Profile curves are a weak test (see the figure for the shape); the control p-values should look random.")


def fmt_ci(v):
    return "%+.2f  [%+.2f, %+.2f]" % (v[0], v[1], v[2])


def print_contrasts(con_df, fracs):
    names = {
        "cc_gap_mean": "cc gap: control minus ON shrink_cc",
        "shrink_Ks_slope_per_decade": "Ks shrinkage, slope per decade of true cc",
        "rho_Ks_cc_slope_per_decade": "r(Ks,cc), slope per decade",
        "w_cc_slope_per_decade": "cc weight, slope per decade",
        "shrink_Ks_hi_minus_lo": "Ks shrinkage, hi minus lo (exploratory)",
        "rho_Ks_cc_hi_minus_lo": "r(Ks,cc), hi minus lo (exploratory)",
        "w_cc_hi_minus_lo": "cc weight, hi minus lo (exploratory)",
    }
    print("\n" + "=" * 118)
    print("2. JOINT BOOTSTRAP -- trends across cc-ON truths.  estimate [95% interval]; 'P>0' = share of re-draws above zero")
    print("=" * 118)
    for frac in fracs:
        print("\n  best %d%% of runs" % round(100 * frac))
        print("  %-44s | %-26s %6s | %-26s %6s" % ("contrast", "all truths", "P>0", "interior truths only", "P>0"))
        for key, label in names.items():
            cells = []
            for tag in ("all", "interior"):
                r = con_df[(con_df.top_frac == frac) & (con_df.contrast == key) & (con_df.truth_set == tag)]
                if r.empty:
                    cells.append(("n/a", np.nan))
                else:
                    r = r.iloc[0]
                    cells.append((fmt_ci((r.estimate, r.ci_lo, r.ci_hi)), r.share_positive))
            print("  %-44s | %-26s %6s | %-26s %6s"
                  % (label, cells[0][0], "%.2f" % cells[0][1] if np.isfinite(cells[0][1]) else "", cells[1][0],
                     "%.2f" % cells[1][1] if np.isfinite(cells[1][1]) else ""))
    print("\n  Reading: an interval that excludes 0 (and P>0 near 0 or 1) is a trend the re-draws rarely reverse. hi/lo split points "
          "were chosen after seeing the data; use the slopes as the headline.")


# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def fig_profiles(prof_curves, truths_show, nb, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(truths_show), 3, figsize=(11.5, 2.9 * len(truths_show)), facecolor=SURFACE, squeeze=False)
    ymin = 1.0
    for t in truths_show:
        ymin = min(ymin, np.nanmin(prof_curves[(prof_curves.truth_cc_mmhr == t)]["best_kge"]))
    ymin = max(0.0, np.floor(ymin * 20) / 20.0)
    for i, t in enumerate(truths_show):
        for j, (pname, col, scale) in enumerate(BIN_PARAMS):
            ax = axes[i][j]
            _style(ax)
            for s, color, ls, mk in ((SER_OFF, COL_NULL, "--", "s"), (SER_ON, COL_CC, "-", "o")):
                d = prof_curves[(prof_curves.truth_cc_mmhr == t) & (prof_curves.series == s) & (prof_curves.parameter == pname)].sort_values("bin")
                ax.plot(d["bin_mid"], d["best_kge"], color=color, linestyle=ls, linewidth=1.8, marker=mk, markersize=4.5,
                        markeredgecolor=SURFACE, markeredgewidth=0.8,
                        label=("channel loss ON" if s == SER_ON else "control: loss OFF (cc inert)"))
            if scale == "log":
                ax.set_xscale("log")
                ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
                ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
            tv = {"Ks": TRUTH_KS, "f": TRUTH_F, "cc": (t if t > 0 else None)}[pname]
            if tv is not None:
                ax.axvline(tv, color=INK2, linewidth=0.9, linestyle=":")
                ax.text(tv, ymin + 0.01, " truth", fontsize=7.5, color=INK2, va="bottom")
            ax.set_ylim(ymin, 1.005)
            ax.set_xlabel({"Ks": "Ks multiplier", "f": "f (decay)", "cc": "cc (mm/hr)"}[pname], color=INK2, fontsize=8.5)
            if j == 0:
                ax.set_ylabel("best KGE_2012 in bin", color=INK2, fontsize=8.5)
            ax.set_title("truth cc = %s" % P.truth_label(t) if j == 0 else "", color=INK, fontsize=9.5, loc="left")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper right", frameon=False, fontsize=8.5, ncol=2, labelcolor=INK2, bbox_to_anchor=(0.99, 0.995))
    fig.suptitle("Profile curves: the best fit reachable at each value of a parameter (%d bins; peaked = pinned, flat = not)" % nb,
                 color=INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_bootstrap(con_df, fracs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [("cc_gap_mean", "cc gap (control minus ON)\npositive = cc pinned"),
            ("shrink_Ks_slope_per_decade", "Ks shrinkage, slope\npositive = Ks less pinned as cc rises"),
            ("rho_Ks_cc_slope_per_decade", "r(Ks, cc), slope\nnegative = trade-off grows"),
            ("w_cc_slope_per_decade", "cc weight in tightest combo, slope\npositive = cc takes part more")]
    fig, axes = plt.subplots(1, len(fracs), figsize=(4.6 * len(fracs) + 1.8, 4.3), sharey=True, facecolor=SURFACE, squeeze=False)
    for ax, frac in zip(axes[0], fracs):
        _style(ax)
        ax.axvline(0, color=INK2, linewidth=1.0)
        for yi, (key, lab) in enumerate(rows):
            for k, (tag, filled, off) in enumerate((("all", True, 0.12), ("interior", False, -0.12))):
                r = con_df[(con_df.top_frac == frac) & (con_df.contrast == key) & (con_df.truth_set == tag)]
                if r.empty:
                    continue
                r = r.iloc[0]
                y = yi + off
                ax.plot([r.ci_lo, r.ci_hi], [y, y], color=COL_CC, linewidth=1.8, solid_capstyle="round")
                ax.plot([r.estimate], [y], marker="o", markersize=6.5, color=COL_CC if filled else SURFACE,
                        markeredgecolor=COL_CC, markeredgewidth=1.6, linestyle="none",
                        label=("all cc-ON truths" if tag == "all" else "interior truths only") if yi == 0 else None)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([lab for _, lab in rows], fontsize=8, color=INK2)
        ax.invert_yaxis()
        ax.set_title("Best %d%% of runs" % round(100 * frac), color=INK, fontsize=10, loc="left")
        ax.set_xlabel("estimate with 95% joint-bootstrap interval", color=INK2, fontsize=8.5)
    axes[0][0].legend(frameon=False, fontsize=8, loc="lower right", labelcolor=INK2)
    fig.suptitle("Do the trends across truths survive re-drawing the runs? (an interval that crosses 0 is not distinguishable from no trend)",
                 color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_permutation(null_store, perm_df, fracs, truths, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Ton = [i for i, t in enumerate(truths) if t > 0]
    if not Ton:
        return False
    fig, axes = plt.subplots(2, len(fracs), figsize=(4.4 * len(fracs) + 0.6, 6.2), facecolor=SURFACE, squeeze=False)
    for ci, frac in enumerate(fracs):
        for ri, s in enumerate((SER_ON, SER_OFF)):
            ax = axes[ri][ci]
            _style(ax)
            nl = np.mean(null_store[(s, frac)]["shrink_cc"][:, Ton], axis=1)
            sub = perm_df[(perm_df.series == s) & (perm_df.top_frac == frac) & (perm_df.truth_cc_mmhr > 0)]
            obs = float(sub["shrink_cc"].mean())
            p = (1.0 + np.sum(nl <= obs + 1e-12)) / (len(nl) + 1.0)
            ax.hist(nl, bins=40, color=COL_NULL, alpha=0.55, edgecolor=SURFACE, linewidth=0.4)
            ax.axvline(obs, color=INK, linewidth=2.0)
            # put the label on the side of the plot away from the black line so the two never overlap
            lo_, hi_ = ax.get_xlim()
            left_side = (obs - lo_) > (hi_ - obs)
            ax.text(0.03 if left_side else 0.97, 0.93,
                    "observed %.3f\nshuffled mean %.3f\np %s" % (obs, nl.mean(), ("< 0.001" if p < 0.001 else "= %.3f" % p)),
                    transform=ax.transAxes, fontsize=8, color=INK, va="top", ha="left" if left_side else "right")
            ax.set_title("Best %d%%  |  %s" % (round(100 * frac), "channel loss ON" if s == SER_ON else "control (cc inert)"),
                         color=INK, fontsize=9.5, loc="left")
            ax.set_xlabel("mean cc shrinkage over cc-ON truths", color=INK2, fontsize=8.5)
            if ci == 0:
                ax.set_ylabel("shuffled data sets", color=INK2, fontsize=8.5)
    fig.suptitle("What 'cc has no effect' looks like (grey: cc shuffled among runs) versus the real data (black line). "
                 "Left of the grey = more pinned than chance", color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True


def fig_noise(summ, main_frac, levels, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = summ[(summ.top_frac == main_frac) & (summ.truth_cc_mmhr > 0)]
    if d.empty:
        return False
    shades = ["#9bdcc2", "#4fc298", "#12805a", "#0a4a33"]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), facecolor=SURFACE, sharey=True)
    panels = (("detection_rate_rho_on", "false_positive_rate_rho_off", "Trade-off test: |r(Ks, cc)| larger than chance"),
              ("detection_rate_on", "false_positive_rate_off", "Pinning test: cc shrinkage smaller than chance"))
    for ax, (col, fpcol, ttl) in zip(axes, panels):
        _style(ax)
        for k, lv in enumerate([0.0] + list(levels)):
            g = d[d.noise_level == lv].sort_values("truth_cc_mmhr")
            if g.empty:
                continue
            c = shades[min(k, len(shades) - 1)]
            lab = "no noise" if lv == 0 else "noise sd %g%%" % (100 * lv)
            # thickest line first so curves that coincide (e.g. all at 100%) are still all visible
            ax.plot(g.truth_cc_mmhr, g[col], color=c, linewidth=[5.5, 4.0, 2.8, 1.8][min(k, 3)], marker="o" if lv > 0 else "D",
                    markersize=5, markeredgecolor=SURFACE, markeredgewidth=0.8, label=lab)
        fp = d[d.noise_level > 0].groupby("truth_cc_mmhr")[fpcol].mean()
        ax.plot(fp.index, fp.values, color=COL_NULL, linestyle="--", linewidth=1.6, label="control (cc inert): share significant")
        ax.set_xscale("log")
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
        ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
        tv = sorted(d.truth_cc_mmhr.unique())
        ax.set_xticks(tv)
        ax.set_xticklabels(["%g" % v for v in tv], rotation=60, fontsize=7.5)
        ax.set_ylim(-0.03, 1.05)
        ax.set_xlabel("true cc of the truth (mm/hr)", color=INK2, fontsize=9)
        ax.set_title(ttl, color=INK, fontsize=9.5, loc="left")
    axes[0].set_ylabel("share of noisy hydrographs with p < 0.05", color=INK2, fontsize=8.5)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="center left")
    fig.suptitle("Does the cc signal survive an imperfect observed hydrograph? (best %d%% of runs; read only truths that are detected with no noise)"
                 % round(100 * main_frac), color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110 -- permutation, joint-bootstrap, profile and noise tests for cc identifiability.")
    ap.add_argument("--top_fracs", type=float, nargs="+", default=[0.20, 0.10, 0.05],
                    help="Top fractions by KGE_2012; the first is the main one (default 0.20 0.10 0.05)")
    ap.add_argument("--boot", type=int, default=1000, help="Joint-bootstrap re-draws (default 1000)")
    ap.add_argument("--perm", type=int, default=5000, help="Permutation shuffles (default 5000)")
    ap.add_argument("--bins", type=int, default=10, help="Bins per parameter for the profile curves (default 10)")
    ap.add_argument("--profile_tol", type=float, default=0.02, help="'Near best' tolerance in KGE for the profile test (default 0.02)")
    ap.add_argument("--lo_max", type=float, default=95.1, help="Largest truth cc in the 'lo' group (default 95.1)")
    ap.add_argument("--hi_min", type=float, default=331.0, help="Smallest truth cc in the 'hi' group (default 331)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--expect_n", type=int, default=250)
    ap.add_argument("--truths", type=float, nargs="+", default=None, help="Only these cc-ON truths (the cc-OFF reference is always kept)")
    ap.add_argument("--figure_truths", type=float, nargs="+", default=[157.0, 425.0])
    ap.add_argument("--noise", action="store_true", help="Also run the noise test (reads the stored hydrographs)")
    ap.add_argument("--noise_levels", type=float, nargs="+", default=[0.05, 0.10, 0.20],
                    help="Noise sd levels as fractions (default 0.05 0.10 0.20)")
    ap.add_argument("--noise_reps", type=int, default=100, help="Noisy hydrographs per level (default 100)")
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--continue_on_fail", action="store_true")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    out_dir = summary_dir / OUT_DIRNAME
    print("\n" + "=" * 78 + "\nSeries 110 -- permutation, bootstrap, profile and noise tests\n" + "=" * 78)
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
    data = build_tensors(df, truths)
    Ton = [i for i, t in enumerate(truths) if t > 0]
    print("\n%d truth(s) (%d cc-ON + %s), %d runs per series, cuts %s, %d shuffles, %d bootstrap re-draws"
          % (len(truths), len(Ton), "cc-OFF reference" if OFF_TRUTH in truths else "no reference",
             data[SER_ON]["X"].shape[0], ["%d%%" % round(100 * f) for f in fracs], args.perm, args.boot))

    # ---- 1. permutation ---------------------------------------------------------
    print("\nRunning permutation tests ...")
    prows, crows, prof_rows, null_store = run_permutation(data, truths, fracs, args.perm, args.bins, args.profile_tol, args.seed)
    perm_df, comb_df, proft_df = pd.DataFrame(prows), pd.DataFrame(crows), pd.DataFrame(prof_rows)
    print_permutation(perm_df, comb_df, main_frac, fracs)

    # ---- 2. joint bootstrap -------------------------------------------------------
    con_df = pd.DataFrame()
    boot_rows = []
    if len(Ton) >= 3:
        print("\nRunning the joint bootstrap ...")
        res = run_bootstrap(data, truths, fracs, args.boot, args.seed)
        pts = point_arrays(data, truths, fracs)
        con_rows = []
        for frac in fracs:
            est = contrast_values(pts[(SER_ON, frac)], pts[(SER_OFF, frac)], truths, args.lo_max, args.hi_min)
            bst = contrast_values(res[(SER_ON, frac)], res[(SER_OFF, frac)], truths, args.lo_max, args.hi_min)
            for name, arr in bst.items():
                key, tag = name.split("|")
                lo, hi = np.nanpercentile(arr, [2.5, 97.5])
                con_rows.append({"top_frac": frac, "contrast": key, "truth_set": tag, "estimate": float(est[name][0]),
                                 "ci_lo": float(lo), "ci_hi": float(hi), "share_positive": float(np.mean(arr > 0)),
                                 "exploratory": key.endswith("hi_minus_lo")})
        con_df = pd.DataFrame(con_rows)
        for (s, frac), store in res.items():
            for ti, t in enumerate(truths):
                row = {"truth_cc_mmhr": t, "series": s, "top_frac": frac}
                for st in STAT_NAMES:
                    lo, hi = np.nanpercentile(store[st][:, ti], [2.5, 97.5])
                    row[st] = float(pts[(s, frac)][st][0, ti])
                    row[st + "_ci_lo"], row[st + "_ci_hi"] = float(lo), float(hi)
                boot_rows.append(row)
        print_contrasts(con_df, fracs)
    else:
        print("\n  (joint bootstrap skipped: needs at least 3 cc-ON truths)")

    # ---- 3. profiles ---------------------------------------------------------------
    prof_curves = pd.DataFrame(run_profiles(data, truths, args.bins))
    print_profile_tests(proft_df)

    # ---- outputs ----------------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ck.rows).to_csv(out_dir / "tests_input_checks_110.csv", index=False)
    perm_df.to_csv(out_dir / "tests_permutation_110.csv", index=False)
    comb_df.to_csv(out_dir / "tests_combined_110.csv", index=False)
    if boot_rows:
        pd.DataFrame(boot_rows).to_csv(out_dir / "tests_bootstrap_truth_110.csv", index=False)
        con_df.to_csv(out_dir / "tests_bootstrap_contrasts_110.csv", index=False)
    prof_curves.to_csv(out_dir / "tests_profiles_110.csv", index=False)
    proft_df.to_csv(out_dir / "tests_profile_tests_110.csv", index=False)

    # ---- 4. noise -----------------------------------------------------------------------
    noise_summ = None
    if args.noise:
        rz, noise_summ = run_noise(args, summary_dir, df, data, truths, fracs)
        rz.to_csv(out_dir / "tests_noise_realizations_110.csv", index=False)
        noise_summ.to_csv(out_dir / "tests_noise_summary_110.csv", index=False)
        print_noise(noise_summ, main_frac, args.noise_levels)

    if not args.no_plots:
        try:
            show = P.pick_figure_truths(truths, args.figure_truths)
            fig_profiles(prof_curves, show, args.bins, out_dir / "fig_tests_profiles_110.png")
            if len(con_df):
                fig_bootstrap(con_df, fracs, out_dir / "fig_tests_bootstrap_110.png")
            fig_permutation(null_store, perm_df, fracs, truths, out_dir / "fig_tests_permutation_110.png")
            if noise_summ is not None:
                fig_noise(noise_summ, main_frac, args.noise_levels, out_dir / "fig_tests_noise_110.png")
        except Exception as e:
            print("\n  (figures skipped: %s -- the CSVs are saved regardless)" % e)

    def md5_of(p):
        h = hashlib.md5()
        h.update(Path(p).read_bytes())
        return h.hexdigest()
    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "pandas": pd.__version__, "numpy": np.__version__,
            "inputs": {n: md5_of(summary_dir / n) for n in (P.STAGE2_NAME, P.CONTROL_NAME)},
            "settings": {k: (v if not isinstance(v, Path) else str(v)) for k, v in vars(args).items()},
            "truths": [float(t) for t in truths], "checks": ck.rows}
    (out_dir / "PROVENANCE_tests_cc_110.json").write_text(json.dumps(prov, indent=2, default=str))

    print("\nSaved to: %s" % out_dir)
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results." % ck.n("FAIL"))
    print("  Caveat: one storm, 250 points in 3-D, routing pinned at truth. A non-significant test with 250 runs does not show cc is "
          "unidentifiable; it may mean too few runs near the best fit.")


if __name__ == "__main__":
    main()
