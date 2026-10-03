"""
analyze_cc_pca_110.py
=====================
Series 110 -- PCA analysis of the recorded runs: how tightly do the best-fitting
runs pin down cc, Ks and f, which parameters trade off against each other, and
which part of the hydrograph each parameter controls.

Pure post-hoc analysis. Reads existing CSVs, runs nothing, writes only into its
own new folder (calibration_work/03_comparisons/summary_tables/pca_cc_110/).

VOCABULARY (short)
-------------------
  PCA (principal component analysis)  Finds the directions along which a cloud
      of points varies most. Each direction is a "principal component" (PC);
      PC1 is the direction of greatest spread, PC2 the next, and so on.
  Loading   The recipe for a PC: how much of each parameter it contains. A PC
      with big loadings on both f and cc is a direction in which f and cc move
      together (a trade-off).
  Eigenvalue   How much spread there is along a PC.
  Behavioral set   The runs that fit the truth best (here: the top 20% by
      KGE_2012). If the best runs are scattered over the whole box, the data did
      not constrain the parameters; if they cluster, it did.
  Equifinality   Different parameter combinations giving nearly equal fits.

WHY THIS DIFFERS FROM LAST WEEK'S PCA (analyze_cc_stage2_110.py, section 3)
-----------------------------------------------------------------------------
Read from that script's code (I have not seen its output):
  1. It scored cc-ON runs against the cc-OFF truth. The truth had no channel
     loss, so the "top performers" were runs where Ks, f and cc happened to
     cancel an unwanted cc. That measures compensation, not whether cc can be
     recovered. Here every cc-ON truth from rescore_cc_truths_110.py is used,
     and the cc-OFF truth is kept as a reference row (truth cc = 0).
  2. It standardized the parameters INSIDE the top-50 subset, so every axis was
     forced to the same spread. A parameter that the data pins down tightly and
     one it does not looked alike. Here each parameter is scaled by the spread of
     the WHOLE sampled box (the "prior"), so 1.0 means "as spread out as the box"
     and 0.2 means "squeezed to a fifth of it".
  3. It reported one PCA with no sense of its uncertainty. Here every number has
     a bootstrap interval, the top-10% cut is shown next to the top-20% cut as a
     sensitivity check, and the cc-OFF control series (where cc has no effect) is
     analysed the same way as the null case.

THE ANALYSES
-------------
  A. BEHAVIORAL PCA   top fraction by KGE_2012 in (Ks, log10 f, log10 cc), scaled
     by the box spread. Reports per parameter "shrinkage" (behavioral spread /
     box spread, with a 95% bootstrap interval), PC eigenvalues and loadings, the
     pairwise rank correlations among the best runs, and whether cc is
     constrained beyond the null (cc-ON interval entirely below cc-OFF interval).
     Note on scale: with a top-20% cut, a parameter that alone fully determined
     the fit would show about 0.2, so ~0.2 is the floor and 1.0 is "unconstrained".
     A second table (A2) gives the loosest and tightest parameter COMBINATIONS.
     Truths within 25% of the cc box edge are flagged "(edge)": the box itself
     clips the ridge there, so some of the shrinkage is not the data's doing.
  B. RESPONSE PCA   all runs of a series: PCA of the hydrograph error metrics
     (KGE r / beta / gamma, peak error, peak timing, first arrival, rising limb,
     time to peak, duration above threshold, recession; rank-transformed so
     outliers cannot dominate). Then each parameter's rank correlation with each
     PC says WHICH response mode it moves. PBIAS is left out because it is
     exactly 100*(beta-1), the same information as beta.

INPUT CHECKS (run first; printed as PASS / WARN / FAIL)
---------------------------------------------------------
  - the Stage 2 / control CSVs are unchanged since rescore_cc_truths_110.py
    gated them (checksums in its provenance file), and that run excluded nothing
  - every (truth, series) has the expected number of runs, unique run_ids, and
    the same run_ids across truths
  - the (Ks, f, cc) values equal the seed-42 Latin hypercube design the sweep
    was supposed to use, and the cc in each run_id label equals the cc column
  - optpercolation is 1 for every cc-ON run and 0 for every cc-OFF run
  - key metrics are finite, parameters lie inside the sampled box
  - null check: in the cc-OFF series, cc must not correlate with KGE_2012
FAIL stops the script (--continue_on_fail to override). WARN is printed and kept.

OUTPUTS   (calibration_work/03_comparisons/summary_tables/pca_cc_110/)
------------------------------------------------------------------------
  pca_input_checks_110.csv        each check, its status and detail
  pca_shrinkage_110.csv           one row per truth x series x top_frac
  pca_behavioral_110.csv          PC eigenvalues and loadings (long)
  pca_response_110.csv            response PCs: variance, metric loadings, parameter correlations
  pca_param_metric_corr_110.csv   each parameter's rank correlation with each hydrograph metric
  fig_pca_shrinkage_110.png       shrinkage vs truth cc, top 20% and top 10%
  fig_pca_loadings_110.png        behavioral PC loadings for the chosen truths
  fig_pca_response_circle_110.png response PCs: metrics and parameters on one plot
  PROVENANCE_pca_cc_110.json      input checksums and settings

USAGE (run from the lab/ directory)
------------------------------------
    python analyze_cc_pca_110.py
    python analyze_cc_pca_110.py --figure_truths 122 425
    python analyze_cc_pca_110.py --top_fracs 0.2 0.1 --boot 500
    python analyze_cc_pca_110.py --truths 157 425            # only these cc-ON truths
    python analyze_cc_pca_110.py --no_plots

CAVEATS (carried into every reading of the output)
----------------------------------------------------
One storm, 250 points in 3-D, routing pinned at truth, noise-free truth: a best
case. A PCA describes straight-line structure; a curved ridge looks wider than it
is. A top-20% set is 50 runs, so loadings are noisy to a few tenths; judge by the
intervals and by whether the 20% and 10% cuts agree.
"""

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ------------------------------------------------------------------
# Constants (match run_cc_stage2_joint_lhs_110.py)
# ------------------------------------------------------------------
SER_ON, SER_OFF = "cc ON", "cc OFF"
OFF_TRUTH = 0.0                       # label for the original cc-OFF truth
STAGE2_NAME   = "lhs_results_joint_Ks_f_cc_110.csv"
CONTROL_NAME  = "lhs_results_joint_Ks_f_cc_CONTROL_110.csv"
RESCORE_DIR   = "rescored_cc_truths_110"
LONG_NAME     = "rescored_long_110.csv"
RESCORE_PROV  = "PROVENANCE_rescored_cc_truths_110.json"
OUT_DIRNAME   = "pca_cc_110"

PARAMS  = ["Ks_mult", "f_RS_abs", "channelconductivity_mmhr"]
PLABEL  = ["Ks", "log f", "log cc"]
BOX = {"Ks_mult": (4.0, 10.5), "f_RS_abs": (0.004, 0.030),
       "channelconductivity_mmhr": (30.0, 1000.0)}
LHS_PARAMS = {                       # order matters: it sets how the random numbers are used
    "Ks_mult":                  {"lo": 4.0,   "hi": 10.5,   "scale": "linear"},
    "f_RS_abs":                 {"lo": 0.004, "hi": 0.030,  "scale": "log"},
    "channelconductivity_mmhr": {"lo": 30.0,  "hi": 1000.0, "scale": "log"},
}
DESIGN_SEED = 42

# Spread (standard deviation) of a uniform sample over the box, in the
# coordinates the PCA uses: Ks linear, f and cc on log10.
PRIOR_STD = np.array([
    (10.5 - 4.0),
    (np.log10(0.030) - np.log10(0.004)),
    (np.log10(1000.0) - np.log10(30.0)),
]) / np.sqrt(12.0)

KEY_METRICS = ["kge_2012", "pbias_pct", "kge_beta", "kge_gamma", "kge_r"]
RESPONSE_METRICS = [
    "kge_r", "kge_beta", "kge_gamma", "peak_error_pct", "peak_timing_error_hr",
    "first_arrival_error_min", "rising_limb_steepness_ratio",
    "time_to_peak_from_exc_min", "duration_above_thresh_error_min",
    "recession_rate_ratio",
]
METRIC_SHORT = {
    "kge_r": "r", "kge_beta": "beta", "kge_gamma": "gamma",
    "peak_error_pct": "peak err", "peak_timing_error_hr": "peak time",
    "first_arrival_error_min": "1st arrival", "rising_limb_steepness_ratio": "rise slope",
    "time_to_peak_from_exc_min": "time to peak", "duration_above_thresh_error_min": "duration",
    "recession_rate_ratio": "recession",
}
_EXTRA = [m for m in RESPONSE_METRICS if m not in KEY_METRICS]
ALL_COLS = ["truth_cc_mmhr", "series", "run_id"] + PARAMS + ["optpercolation"] + KEY_METRICS + _EXTRA
NUMERIC_COLS = [c for c in ALL_COLS if c not in ("series", "run_id")]

NEAR_EDGE   = 0.25          # truth cc within this fraction of the log-cc box edge -> "edge" flag
MIN_TOP     = 10
MIN_ROWS    = 30
NULL_RHO    = 0.20          # |rho(cc, KGE)| in the cc-OFF series above this -> WARN

# Okabe-Ito-style hues, validated for 3-series all-pairs use (CVD and normal vision).
C_KS, C_F, C_CC = "#2a78d6", "#eb6834", "#1baf7a"
PCOLORS = [C_KS, C_F, C_CC]
INK, INK2, GRID, SURFACE, NULLC = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb", "#8a8984"


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


SERIES_TITLE = {SER_ON: "channel loss ON (cc varies)", SER_OFF: "control: loss OFF (cc inert)"}


def truth_label(t):
    return "cc-OFF truth" if t == OFF_TRUTH else "%g" % t


def to_X(df):
    """Parameters in PCA coordinates: Ks, log10 f, log10 cc."""
    return np.column_stack([
        df["Ks_mult"].to_numpy(float),
        np.log10(df["f_RS_abs"].to_numpy(float)),
        np.log10(df["channelconductivity_mmhr"].to_numpy(float)),
    ])


def generate_lhs_samples(n, params, seed):
    """Copied from run_cc_stage2_joint_lhs_110.py so the intended design can be rebuilt."""
    rng = np.random.default_rng(seed)
    samples = {}
    for param, b in params.items():
        lo, hi, scale = b["lo"], b["hi"], b.get("scale", "linear")
        if scale == "log":
            edges = np.linspace(np.log10(lo), np.log10(hi), n + 1)
            pts = 10 ** rng.uniform(edges[:-1], edges[1:])
        else:
            edges = np.linspace(lo, hi, n + 1)
            pts = rng.uniform(edges[:-1], edges[1:])
        rng.shuffle(pts)
        samples[param] = pts
    return pd.DataFrame(samples)


# ------------------------------------------------------------------
# Loading
# ------------------------------------------------------------------
def load_all(summary_dir):
    notes = []
    frames = []
    for name, series in ((STAGE2_NAME, SER_ON), (CONTROL_NAME, SER_OFF)):
        p = summary_dir / name
        if not p.exists():
            sys.exit("Required file not found: %s" % p)
        d = pd.read_csv(p)
        d["truth_cc_mmhr"] = OFF_TRUTH
        d["series"] = series
        frames.append(d)
    long_path = summary_dir / RESCORE_DIR / LONG_NAME
    if long_path.exists():
        frames.append(pd.read_csv(long_path))
    else:
        notes.append("%s not found: only the cc-OFF truth can be analysed. Run "
                     "rescore_cc_truths_110.py first for the cc-ON truths." % long_path.name)
    df = pd.concat(frames, ignore_index=True, sort=False)
    for c in ALL_COLS:
        if c not in df.columns:
            df[c] = np.nan
    df = df[ALL_COLS].copy()
    for c in NUMERIC_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df, notes


# ------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------
class Checks:
    def __init__(self):
        self.rows = []

    def add(self, status, name, detail):
        self.rows.append({"status": status, "check": name, "detail": detail})

    def n(self, status):
        return sum(1 for r in self.rows if r["status"] == status)


def run_checks(df, summary_dir, expect_n):
    ck = Checks()
    truths = sorted(df["truth_cc_mmhr"].unique())

    # --- 1. provenance of the re-scoring ---------------------------------
    prov_path = summary_dir / RESCORE_DIR / RESCORE_PROV
    has_long = any(t != OFF_TRUTH for t in truths)
    if not has_long:
        ck.add("WARN", "re-scoring provenance", "no re-scored truths loaded; nothing to verify")
    elif not prov_path.exists():
        ck.add("WARN", "re-scoring provenance",
               "%s missing; cannot confirm the re-scoring gates passed" % RESCORE_PROV)
    else:
        prov = json.loads(prov_path.read_text())
        bad = []
        for name in (STAGE2_NAME, CONTROL_NAME):
            rec = prov.get("inputs", {}).get(name, {}).get("md5")
            now = md5_of(summary_dir / name)
            if rec != now:
                bad.append(name)
        if bad:
            ck.add("FAIL", "inputs unchanged since re-scoring",
                   "checksum differs from the one recorded by the re-scoring run: %s. "
                   "Re-run rescore_cc_truths_110.py." % ", ".join(bad))
        else:
            ck.add("PASS", "inputs unchanged since re-scoring",
                   "Stage 2 and control CSV checksums match the re-scoring provenance")
        gates = prov.get("gates", {})
        excl = gates.get("runs_excluded", {})
        if excl:
            ck.add("WARN", "re-scoring excluded runs",
                   "%d run(s) were excluded at re-scoring (--drop_bad_runs)" % len(excl))
        else:
            ck.add("PASS", "re-scoring gates", "all %s runs passed gates A and B at re-scoring"
                   % gates.get("runs_checked", "?"))
        if gates.get("phase_metrics_from_scorer") is False:
            ck.add("WARN", "phase metrics", "re-scoring left phase metrics blank; response PCA uses fewer metrics")
        prov_cc = [v["cc_mmhr"] for v in prov.get("candidate_truths", {}).values()]
        extra = [t for t in truths if t != OFF_TRUTH and not any(np.isclose(t, c, rtol=1e-6) for c in prov_cc)]
        if extra:
            ck.add("FAIL", "truths verified at re-scoring",
                   "truth(s) %s are not in the verified candidate list" % [float(x) for x in extra])
        else:
            ck.add("PASS", "truths verified at re-scoring",
                   "%d candidate truth(s), all checksum-verified at re-scoring" % len([t for t in truths if t != OFF_TRUTH]))

    # --- 2. counts and run identity ---------------------------------------
    problems = []
    ids_by_series = {}
    for (t, s), g in df.groupby(["truth_cc_mmhr", "series"]):
        if len(g) != expect_n:
            problems.append("truth %s %s: %d runs (expected %d)" % (truth_label(t), s, len(g), expect_n))
        if g["run_id"].duplicated().any():
            problems.append("truth %s %s: duplicate run_id" % (truth_label(t), s))
        ids_by_series.setdefault(s, []).append(frozenset(g["run_id"]))
    for s, sets in ids_by_series.items():
        if len(set(sets)) > 1:
            problems.append("%s: run_id sets differ between truths" % s)
    if problems:
        ck.add("FAIL", "run counts and identity", "; ".join(problems[:4]))
    else:
        ck.add("PASS", "run counts and identity",
               "%d truth(s) x 2 series, %d unique runs each, same runs for every truth" % (len(truths), expect_n))

    # --- 3. series flags --------------------------------------------------
    on_bad = int((df.loc[df.series == SER_ON, "optpercolation"] != 1).sum())
    off_bad = int((df.loc[df.series == SER_OFF, "optpercolation"] != 0).sum())
    if on_bad or off_bad:
        ck.add("FAIL", "optpercolation flags", "%d cc-ON rows not 1, %d cc-OFF rows not 0" % (on_bad, off_bad))
    else:
        ck.add("PASS", "optpercolation flags", "cc-ON rows all 1, cc-OFF rows all 0")

    # --- 4. design: parameters equal the intended LHS ---------------------
    design = generate_lhs_samples(expect_n, LHS_PARAMS, DESIGN_SEED)
    D = design[PARAMS].to_numpy(float)
    ref_truth = OFF_TRUTH if OFF_TRUTH in truths else truths[0]
    for s in (SER_ON, SER_OFF):
        g = df[(df.truth_cc_mmhr == ref_truth) & (df.series == s)]
        R = g[PARAMS].to_numpy(float)
        if len(R) == 0 or np.isnan(R).any():
            ck.add("FAIL", "design match (%s)" % s, "missing or non-numeric parameter values")
            continue
        close = np.all(np.isclose(R[:, None, :], D[None, :, :], rtol=1e-3, atol=0.0), axis=2)
        hits = close.sum(axis=1)
        unmatched = int((hits == 0).sum())
        used = close.any(axis=0)
        if unmatched == 0 and used.sum() == min(len(R), len(D)):
            ck.add("PASS", "design match (%s)" % s,
                   "all %d (Ks, f, cc) points equal the seed-%d, n=%d design (rtol 1e-3)" % (len(R), DESIGN_SEED, expect_n))
        else:
            ck.add("FAIL", "design match (%s)" % s,
                   "%d point(s) do not match the seed-%d design; %d of %d design points are used"
                   % (unmatched, DESIGN_SEED, int(used.sum()), len(D)))

    # --- 5. cc in the run_id label equals the cc column --------------------
    g = df[df.truth_cc_mmhr == ref_truth]
    parsed, mism, noparse = 0, 0, 0
    for rid, cc in zip(g["run_id"], g["channelconductivity_mmhr"]):
        m = re.search(r"_cc([0-9p]+)_(?:s2|ctl)$", str(rid))
        if not m:
            noparse += 1
            continue
        parsed += 1
        if not np.isclose(float(m.group(1).replace("p", ".")), cc, rtol=1e-6):
            mism += 1
    if parsed == 0:
        ck.add("WARN", "run_id cc label", "no run_id matched the expected pattern; label check skipped")
    elif mism:
        ck.add("FAIL", "run_id cc label", "%d of %d run_id labels disagree with the cc column" % (mism, parsed))
    else:
        ck.add("PASS", "run_id cc label",
               "cc in %d run_id labels equals the cc column%s" % (parsed, "" if not noparse else " (%d not parseable)" % noparse))

    # --- 6. finite metrics --------------------------------------------------
    worst = (1.0, "")
    for (t, s), gg in df.groupby(["truth_cc_mmhr", "series"]):
        for m in KEY_METRICS:
            frac = float(np.isfinite(gg[m]).mean())
            if frac < worst[0]:
                worst = (frac, "truth %s %s %s" % (truth_label(t), s, m))
    if worst[0] < 0.99:
        ck.add("FAIL", "finite metrics", "worst finite fraction %.3f (%s)" % worst)
    else:
        ck.add("PASS", "finite metrics", "all key metrics >= 99%% finite (worst %.3f)" % worst[0])

    # --- 7. parameters inside the box ----------------------------------------
    out = []
    for p in PARAMS:
        lo, hi = BOX[p]
        v = df[p].to_numpy(float)
        if np.nanmin(v) < lo * (1 - 1e-9) or np.nanmax(v) > hi * (1 + 1e-9):
            out.append(p)
    if out:
        ck.add("FAIL", "parameters inside sampled box", "outside the box: %s" % ", ".join(out))
    else:
        ck.add("PASS", "parameters inside sampled box", "Ks, f and cc all inside the Stage 2 box")

    # --- 8. null check: cc is inert in the control series ----------------------
    worst_rho, worst_t = 0.0, None
    for t in truths:
        g = df[(df.truth_cc_mmhr == t) & (df.series == SER_OFF)]
        ok = np.isfinite(g["kge_2012"]) & np.isfinite(g["channelconductivity_mmhr"])
        if ok.sum() > 10:
            rho = stats.spearmanr(np.log10(g.loc[ok, "channelconductivity_mmhr"]), g.loc[ok, "kge_2012"])[0]
            if abs(rho) > abs(worst_rho):
                worst_rho, worst_t = rho, t
    if worst_t is not None and abs(worst_rho) > NULL_RHO:
        ck.add("WARN", "null check (cc inert in control)",
               "rho(log cc, KGE_2012) = %+.2f at truth %s in the cc-OFF series (limit %.2f)"
               % (worst_rho, truth_label(worst_t), NULL_RHO))
    else:
        ck.add("PASS", "null check (cc inert in control)",
               "largest |rho(log cc, KGE_2012)| in the cc-OFF series is %.2f (limit %.2f)" % (abs(worst_rho), NULL_RHO))
    return ck


# ------------------------------------------------------------------
# A. Behavioral PCA
# ------------------------------------------------------------------
def _top_indices(kge, frac):
    ok = np.flatnonzero(np.isfinite(kge))
    n_top = max(MIN_TOP, int(round(frac * len(ok))))
    order = ok[np.argsort(-kge[ok], kind="stable")][:n_top]
    return order


def behavioral_pca(X, kge, frac):
    idx = _top_indices(kge, frac)
    Xt = X[idx]
    Z = (Xt - Xt.mean(axis=0)) / PRIOR_STD
    C = np.cov(Z, rowvar=False)
    evals, evecs = np.linalg.eigh(C)
    o = np.argsort(evals)[::-1]
    evals, evecs = evals[o], evecs[:, o]
    for j in range(3):                                   # sign: largest loading positive
        if evecs[np.argmax(np.abs(evecs[:, j])), j] < 0:
            evecs[:, j] *= -1
    rho = lambda a, b: float(stats.spearmanr(Xt[:, a], Xt[:, b])[0])
    return {
        "n_top": int(len(idx)), "kge_cutoff": float(kge[idx[-1]]),
        "shrink": np.sqrt(np.diag(C)), "evals": evals, "evecs": evecs,
        "vol_ratio": float(np.sqrt(max(np.linalg.det(C), 0.0))),
        "rho_Ks_f": rho(0, 1), "rho_Ks_cc": rho(0, 2), "rho_f_cc": rho(1, 2),
    }


def bootstrap_shrinkage(X, kge, frac, B, rng):
    ok = np.flatnonzero(np.isfinite(kge))
    X, kge = X[ok], kge[ok]
    n = len(kge)
    n_top = max(MIN_TOP, int(round(frac * n)))
    out = np.empty((B, 3))
    for b in range(B):
        idx = rng.integers(0, n, n)
        order = np.argsort(-kge[idx], kind="stable")[:n_top]
        out[b] = X[idx][order].std(axis=0, ddof=1) / PRIOR_STD
    return np.percentile(out, [2.5, 97.5], axis=0)       # (2, 3)


# ------------------------------------------------------------------
# B. Response PCA
# ------------------------------------------------------------------
def normal_scores(v):
    r = stats.rankdata(v)
    return stats.norm.ppf((r - 0.5) / len(v))


def response_pca(sub, X):
    use = []
    for m in RESPONSE_METRICS:
        v = sub[m].to_numpy(float)
        if np.isfinite(v).mean() >= 0.90 and np.nanstd(v) > 0:
            use.append(m)
    if len(use) < 3:
        return None
    M = sub[use].to_numpy(float)
    ok = np.isfinite(M).all(axis=1) & np.isfinite(X).all(axis=1)
    if ok.sum() < MIN_ROWS:
        return None
    M, Xp = M[ok], X[ok]
    Z = np.column_stack([normal_scores(M[:, j]) for j in range(M.shape[1])])
    sd = Z.std(axis=0, ddof=1)
    if (sd == 0).any():
        return None
    Z = (Z - Z.mean(axis=0)) / sd
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    n = Z.shape[0]
    lam = S ** 2 / (n - 1)
    scores = U * S
    for k in range(Vt.shape[0]):                          # sign: strongest metric loading positive
        if Vt[k, np.argmax(np.abs(Vt[k]))] < 0:
            Vt[k] *= -1
            scores[:, k] *= -1
    k_use = min(3, len(use))
    corr_metric = Vt[:k_use] * np.sqrt(lam[:k_use])[:, None]      # metric-PC correlations
    rho = np.array([[stats.spearmanr(Xp[:, p], scores[:, k])[0] for k in range(k_use)]
                    for p in range(3)])                            # (3 params, k_use)
    rho_metric = np.array([[stats.spearmanr(Xp[:, p], M[:, j])[0] for j in range(M.shape[1])] for p in range(3)])
    return {"metrics": use, "n": int(n), "evr": lam / lam.sum(), "loadings": Vt[:k_use],
            "corr_metric": corr_metric, "rho": rho, "rho_metric": rho_metric}


# ------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------
def print_checks(ck):
    print("\n" + "=" * 100)
    print("INPUT CHECKS")
    print("=" * 100)
    for r in ck.rows:
        print("  [%s] %-36s %s" % (r["status"], r["check"], r["detail"]))
    print("  -> %d PASS, %d WARN, %d FAIL" % (ck.n("PASS"), ck.n("WARN"), ck.n("FAIL")))


def print_behavioral(shr, main_frac):
    d = shr[shr["top_frac"] == main_frac]
    on = d[d.series == SER_ON].set_index("truth_cc_mmhr")
    off = d[d.series == SER_OFF].set_index("truth_cc_mmhr")
    print("\n" + "=" * 118)
    print("A. BEHAVIORAL PCA -- top %d%% of runs by KGE_2012.  Shrinkage = spread of the best runs / spread of the whole "
          "box (1.0 = unconstrained)" % round(100 * main_frac))
    print("=" * 118)
    print("%12s %4s %8s | %5s %5s %5s  %-15s | %-15s | %-12s | %s"
          % ("truth cc", "n", "KGE >=", "Ks", "f", "cc", "cc 95% interval", "cc-OFF control", "cc constrained",
             "PC1 eig; loadings (Ks, f, cc)"))
    for t in sorted(on.index):
        a = on.loc[t]
        b = off.loc[t] if t in off.index else None
        ci = "%.2f-%.2f" % (a["ci_lo_cc"], a["ci_hi_cc"])
        bci = ("%.2f (%.2f-%.2f)" % (b["shrink_cc"], b["ci_lo_cc"], b["ci_hi_cc"])) if b is not None else "n/a"
        cons = "n/a" if pd.isna(a.get("cc_constrained_vs_null")) else ("YES" if a["cc_constrained_vs_null"] else "no")
        if t == OFF_TRUTH and cons == "YES":
            cons = "YES (see *)"
        elif cons == "YES" and a.get("near_box_edge") is True:
            cons = "YES (edge)"
        print("%12s %4d %8.3f | %5.2f %5.2f %5.2f  %-15s | %-15s | %-12s | %.2f; %+.2f %+.2f %+.2f"
              % (truth_label(t), a["n_top"], a["kge_cutoff"], a["shrink_Ks"], a["shrink_f"], a["shrink_cc"], ci, bci,
                 cons, a["eig1"], a["pc1_Ks"], a["pc1_logf"], a["pc1_logcc"]))
    print("  KGE >= : lowest KGE_2012 inside the top set (how loose the cut is).  A single parameter that alone "
          "decided the fit would show about %.2f." % main_frac)
    print("  cc constrained: cc-ON interval lies entirely below the cc-OFF control interval (the control is the null: "
          "cc has no effect there).")
    print("  (edge): the truth sits within %d%% of the cc box edge, so part of the shrinkage is the box clipping the "
          "ridge, not the data. Treat as weaker evidence." % round(100 * NEAR_EDGE))
    if OFF_TRUTH in on.index:
        print("  * cc-OFF truth: the truth has no channel loss, so YES means the best runs are pushed toward LOW cc "
              "(compensation), not that a true cc was recovered.")


def print_directions(shr, main_frac):
    d = shr[(shr["top_frac"] == main_frac) & (shr.series == SER_ON)].sort_values("truth_cc_mmhr")
    print("\n" + "=" * 118)
    print("A2. TRADE-OFFS AMONG THE BEST %d%% -- loosest and tightest parameter combinations, and pairwise rank "
          "correlations" % round(100 * main_frac))
    print("=" * 118)
    print("%12s | %-34s | %-34s | %6s %6s %6s"
          % ("truth cc", "loosest combo PC1: eig; Ks f cc", "tightest combo PC3: eig; Ks f cc", "r Ks-f", "r Ks-cc", "r f-cc"))
    for _, a in d.iterrows():
        print("%12s | %-34s | %-34s | %+6.2f %+6.2f %+6.2f"
              % (truth_label(a["truth_cc_mmhr"]),
                 "%.2f; %+.2f %+.2f %+.2f" % (a["eig1"], a["pc1_Ks"], a["pc1_logf"], a["pc1_logcc"]),
                 "%.2f; %+.2f %+.2f %+.2f" % (a["eig3"], a["pc3_Ks"], a["pc3_logf"], a["pc3_logcc"]),
                 a["rho_Ks_f"], a["rho_Ks_cc"], a["rho_f_cc"]))
    print("  tightest combo: the mix of Ks, f and cc that the best runs agree on most (small eig = well pinned down).")
    print("  loosest combo: the mix they disagree on most (the trade-off direction). Signs only matter relative to each other.")
    print("  r = rank correlation among the best runs; a large |r| between two parameters means one can offset the other.")


def print_param_metric(resp_rows, truths):
    rows = {(r["truth_cc_mmhr"], r["series"]): r for r in resp_rows}
    for t in truths:
        r = rows.get((t, SER_ON))
        if r is None or "rho_metric" not in r:
            continue
        ms = r["metrics"]
        print("\n  Parameter -> metric rank correlation, truth cc = %s, channel loss ON (all %d runs)" % (truth_label(t), r["n"]))
        print("  %-8s " % "" + " ".join("%11s" % METRIC_SHORT[m] for m in ms))
        for p, lab in enumerate(PLABEL):
            print("  %-8s " % lab + " ".join("%+11.2f" % v for v in r["rho_metric"][p]))


def print_response(resp_rows):
    if not resp_rows:
        print("\nB. RESPONSE PCA: not enough usable metrics or runs to run.")
        return
    r = pd.DataFrame(resp_rows)
    on = r[r.series == SER_ON]
    off = r[r.series == SER_OFF].set_index("truth_cc_mmhr")
    print("\n" + "=" * 118)
    print("B. RESPONSE PCA -- all runs of the series; which response mode does each parameter move?   "
          "(rank correlation with PC scores)")
    print("=" * 118)
    print("%12s %3s | %6s %6s | %-16s %-16s %-16s | %s"
          % ("truth cc", "m", "PC1 %", "PC2 %", "Ks -> PC (rho)", "f -> PC (rho)", "cc -> PC (rho)", "cc-OFF control: max |rho| cc"))
    for _, a in on.sort_values("truth_cc_mmhr").iterrows():
        rho = np.array(a["rho"])
        cells = []
        for p in range(3):
            k = int(np.argmax(np.abs(rho[p])))
            cells.append("PC%d (%+.2f)" % (k + 1, rho[p, k]))
        t = a["truth_cc_mmhr"]
        offtxt = "n/a"
        if t in off.index:
            offtxt = "%.2f" % np.max(np.abs(np.array(off.loc[t, "rho"])[2]))
        print("%12s %3d | %6.1f %6.1f | %-16s %-16s %-16s | %s"
              % (truth_label(t), len(a["metrics"]), 100 * a["evr"][0], 100 * a["evr"][1], cells[0], cells[1], cells[2], offtxt))
    print("  m = number of response metrics used.  Compare cc's row with Ks and f: the same PC means they move the same "
          "part of the hydrograph.")


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


def fig_shrinkage(shr, fracs, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fr = fracs[:2]
    fig, axes = plt.subplots(1, len(fr), figsize=(6.2 * len(fr), 4.6), sharey=True, facecolor=SURFACE)
    axes = np.atleast_1d(axes)
    for ax, frac in zip(axes, fr):
        _style(ax)
        d = shr[(shr.top_frac == frac) & (shr.truth_cc_mmhr > 0)]
        on = d[d.series == SER_ON].sort_values("truth_cc_mmhr")
        off = d[d.series == SER_OFF].sort_values("truth_cc_mmhr")
        x = on["truth_cc_mmhr"].to_numpy()
        ax.axhline(1.0, color=INK2, linewidth=0.8, linestyle=":")
        ax.text(x.max(), 1.29, "1.0 = as spread out as the whole sampled box", fontsize=7.5, color=INK2,
                va="center", ha="right")
        ax.fill_between(off["truth_cc_mmhr"], off["ci_lo_cc"], off["ci_hi_cc"], color=NULLC, alpha=0.22, linewidth=0)
        ax.plot(off["truth_cc_mmhr"], off["shrink_cc"], color=NULLC, linewidth=1.6, linestyle="--")
        ax.fill_between(x, on["ci_lo_cc"], on["ci_hi_cc"], color=C_CC, alpha=0.18, linewidth=0)
        for col, key, mk, lab in ((C_KS, "shrink_Ks", "o", "Ks"), (C_F, "shrink_f", "s", "f"), (C_CC, "shrink_cc", "D", "cc")):
            ax.plot(x, on[key], color=col, linewidth=1.8, marker=mk, markersize=5, markeredgecolor=SURFACE,
                    markeredgewidth=1.0, label=lab + (" (channel loss on)" if lab == "cc" else ""))
        ax.plot([], [], color=NULLC, linestyle="--", linewidth=1.6, label="cc, control with loss off (null)")
        ax.set_xscale("log")
        ax.set_xticks(x[::2] if len(x) > 6 else x)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
        ax.get_xaxis().set_minor_locator(matplotlib.ticker.NullLocator())
        ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_ylim(0, 1.36)
        ax.set_xlabel("true cc of the truth (mm/hr)", color=INK2, fontsize=9)
        ax.set_title("Best %d%% of runs (n = %d)" % (round(100 * frac), int(on["n_top"].iloc[0])),
                     color=INK, fontsize=10, loc="left")
    axes[0].set_ylabel("shrinkage (spread of best runs / spread of box)", color=INK2, fontsize=9)
    axes[0].legend(fontsize=7.5, frameon=False, loc="lower right", labelcolor=INK2)
    fig.suptitle("How tightly the best-fitting runs pin each parameter (lower = better constrained; shaded = 95% bootstrap interval)",
                 color=INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_loadings(beh, shr, truths, frac, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(truths), 2, figsize=(9.0, 2.7 * len(truths)), sharey=True, facecolor=SURFACE, squeeze=False)
    for i, t in enumerate(truths):
        for j, s in enumerate((SER_ON, SER_OFF)):
            ax = axes[i][j]
            _style(ax)
            d = beh[(beh.truth_cc_mmhr == t) & (beh.series == s) & (beh.top_frac == frac)].sort_values("component")
            if d.empty:
                ax.set_visible(False)
                continue
            xs = np.arange(len(d))
            w = 0.26
            for k, (col, key, lab) in enumerate(zip(PCOLORS, ("loading_Ks", "loading_logf", "loading_logcc"), PLABEL)):
                ax.bar(xs + (k - 1) * w, d[key], w * 0.92, color=col, label=lab)
            for x0, lam in zip(xs, d["eigenvalue_ratio"]):
                ax.text(x0, 1.06, "eig %.2f" % lam, ha="center", fontsize=7.5, color=INK2)
            ax.axhline(0, color=INK2, linewidth=0.8)
            ax.set_xticks(xs)
            ax.set_xticklabels(d["component"], fontsize=8)
            ax.set_ylim(-1.0, 1.2)
            ax.set_title("truth cc = %s  |  %s" % (truth_label(t), SERIES_TITLE[s]), color=INK, fontsize=9, loc="left")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, frameon=False, loc="upper right", labelcolor=INK2, ncol=3,
               bbox_to_anchor=(0.99, 0.995))
    for i in range(len(truths)):
        axes[i][0].set_ylabel("loading", color=INK2, fontsize=8)
    fig.suptitle("What each principal component is made of, best %d%% of runs\n"
                 "(eig = spread along the PC relative to the whole box)" % round(100 * frac),
                 color=INK, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_circle(resp_rows, truths, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = {(r["truth_cc_mmhr"], r["series"]): r for r in resp_rows}
    use = [t for t in truths if (t, SER_ON) in rows]
    if not use:
        return False
    fig, axes = plt.subplots(1, len(use), figsize=(max(5.2 * len(use), 7.5), 6.6), facecolor=SURFACE, squeeze=False)
    for ax, t in zip(axes[0], use):
        r = rows[(t, SER_ON)]
        ax.set_facecolor(SURFACE)
        ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, color=GRID, linewidth=1.0))
        ax.axhline(0, color=GRID, linewidth=0.8)
        ax.axvline(0, color=GRID, linewidth=0.8)
        cm = r["corr_metric"]
        groups = []                                   # metrics whose arrow tips nearly coincide share one label
        for j, m in enumerate(r["metrics"]):
            ax.annotate("", xy=(cm[0, j], cm[1, j]), xytext=(0, 0),
                        arrowprops=dict(arrowstyle="-|>", color=NULLC, lw=0.9, shrinkA=0, shrinkB=0))
            for g in groups:
                if np.hypot(cm[0, g[0]] - cm[0, j], cm[1, g[0]] - cm[1, j]) < 0.22:
                    g.append(j)
                    break
            else:
                groups.append([j])
        for g in groups:
            gx, gy = float(np.mean(cm[0, g])), float(np.mean(cm[1, g]))
            names = "\n".join(METRIC_SHORT[r["metrics"][j]] for j in g)
            norm = max(np.hypot(gx, gy), 1e-9)
            ax.text(gx * 1.10 + 0.0, gy * 1.10 + 0.02 * len(g) * np.sign(gy if gy else 1), names, fontsize=7.2,
                    color=INK2, ha="center", va="center", linespacing=0.95)
        for p, (col, lab) in enumerate(zip(PCOLORS, PLABEL)):
            x, y = r["rho"][p, 0], r["rho"][p, 1]
            ax.plot([0, x], [0, y], color=col, linewidth=2.2, solid_capstyle="round")
            ax.plot([x], [y], marker=["o", "s", "D"][p], color=col, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.0)
            ax.text(x + (0.06 if x >= 0 else -0.06), y - 0.07, lab, fontsize=9, color=INK, fontweight="bold",
                    ha="left" if x >= 0 else "right", va="center")
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.set_aspect("equal")
        ax.set_xlabel("PC1 (%.0f%% of response variance)" % (100 * r["evr"][0]), color=INK2, fontsize=8.5)
        ax.set_ylabel("PC2 (%.0f%%)" % (100 * r["evr"][1]), color=INK2, fontsize=8.5)
        ax.tick_params(colors=INK2, labelsize=7.5)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title("truth cc = %s  |  %s" % (truth_label(t), SERIES_TITLE[SER_ON]), color=INK, fontsize=9.5, loc="left")
    fig.suptitle("Which response each parameter moves\n"
                 "gray arrows = hydrograph error metrics, coloured lines = parameters; lines pointing the same way "
                 "move the same part of the hydrograph.\n"
                 "A short coloured line only means the parameter barely moves PC1 and PC2; "
                 "check pca_param_metric_corr_110.csv for the rest.",
                 color=INK, fontsize=9.5, x=0.01, ha="left")
    fig.subplots_adjust(left=0.05, right=0.99, top=0.84, bottom=0.10, wspace=0.12)
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
def pick_figure_truths(available, wanted):
    out = []
    cc_av = [t for t in available if t != OFF_TRUTH]
    for w in wanted:
        if not cc_av:
            break
        t = min(cc_av, key=lambda c: abs(np.log(c / w)))
        if t not in out:
            out.append(t)
    if OFF_TRUTH in available:
        out.append(OFF_TRUTH)
    return out


def main():
    ap = argparse.ArgumentParser(description="Series 110 -- PCA of the recorded runs (cc, Ks, f and the hydrograph response).")
    ap.add_argument("--top_fracs", type=float, nargs="+", default=[0.20, 0.10],
                    help="Top fractions by KGE_2012 for the behavioral PCA; the first is the main one (default 0.20 0.10)")
    ap.add_argument("--boot", type=int, default=300, help="Bootstrap repeats for the intervals (default 300)")
    ap.add_argument("--seed", type=int, default=0, help="Random seed for the bootstrap (default 0)")
    ap.add_argument("--expect_n", type=int, default=250, help="Expected runs per series (default 250)")
    ap.add_argument("--truths", type=float, nargs="+", default=None,
                    help="Only these cc-ON truths (mm/hr); the cc-OFF reference is always kept")
    ap.add_argument("--figure_truths", type=float, nargs="+", default=[157.0, 425.0],
                    help="cc-ON truths to show in the loadings / circle figures (default 157 425; nearest available)")
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--continue_on_fail", action="store_true", help="Carry on even if an input check FAILs")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()

    summary_dir = args.summary_dir or (Path.cwd().parent / "calibration_work" / "03_comparisons" / "summary_tables")
    out_dir = summary_dir / OUT_DIRNAME
    print("\n" + "=" * 78 + "\nSeries 110 -- PCA of the recorded runs\n" + "=" * 78)
    print("Reading from: %s" % summary_dir)

    df, notes = load_all(summary_dir)
    for n in notes:
        print("  NOTE: %s" % n)

    ck = run_checks(df, summary_dir, args.expect_n)
    print_checks(ck)
    if ck.n("FAIL") and not args.continue_on_fail:
        sys.exit("\nSTOPPED: an input check failed. Nothing was analysed or written. "
                 "(--continue_on_fail to override.)")

    if args.truths:
        keep = [OFF_TRUTH] + [t for t in sorted(df.truth_cc_mmhr.unique())
                               if t != OFF_TRUTH and any(np.isclose(t, v, rtol=1e-3) for v in args.truths)]
        df = df[df.truth_cc_mmhr.isin(keep)]
        if len(keep) == 1:
            print("  NOTE: --truths matched no cc-ON truth; analysing the cc-OFF reference only.")
    truths = sorted(df["truth_cc_mmhr"].unique())

    shr_rows, beh_rows, resp_rows = [], [], []
    for t in truths:
        for s in (SER_ON, SER_OFF):
            sub = df[(df.truth_cc_mmhr == t) & (df.series == s)].reset_index(drop=True)
            X = to_X(sub)
            kge = sub["kge_2012"].to_numpy(float)
            if np.isfinite(kge).sum() < 2 * MIN_TOP:
                continue
            for frac in args.top_fracs:
                res = behavioral_pca(X, kge, frac)
                rng = np.random.default_rng([args.seed, int(round(t * 100)), 1 if s == SER_ON else 2, int(round(frac * 1000))])
                ci = bootstrap_shrinkage(X, kge, frac, args.boot, rng)
                ev, V = res["evals"], res["evecs"]
                row = {"truth_cc_mmhr": t, "series": s, "top_frac": frac, "n_top": res["n_top"],
                       "kge_cutoff": res["kge_cutoff"],
                       "shrink_Ks": res["shrink"][0], "shrink_f": res["shrink"][1], "shrink_cc": res["shrink"][2],
                       "ci_lo_Ks": ci[0, 0], "ci_hi_Ks": ci[1, 0], "ci_lo_f": ci[0, 1], "ci_hi_f": ci[1, 1],
                       "ci_lo_cc": ci[0, 2], "ci_hi_cc": ci[1, 2],
                       "ellipsoid_volume_ratio": res["vol_ratio"],
                       "rho_Ks_f": res["rho_Ks_f"], "rho_Ks_cc": res["rho_Ks_cc"], "rho_f_cc": res["rho_f_cc"],
                       "eig1": ev[0], "eig2": ev[1], "eig3": ev[2],
                       "pc1_Ks": V[0, 0], "pc1_logf": V[1, 0], "pc1_logcc": V[2, 0],
                       "pc3_Ks": V[0, 2], "pc3_logf": V[1, 2], "pc3_logcc": V[2, 2]}
                shr_rows.append(row)
                for k in range(3):
                    beh_rows.append({"truth_cc_mmhr": t, "series": s, "top_frac": frac, "n_top": res["n_top"],
                                     "kge_cutoff": res["kge_cutoff"], "component": "PC%d" % (k + 1),
                                     "eigenvalue_ratio": ev[k], "explained_variance_ratio": ev[k] / ev.sum(),
                                     "loading_Ks": V[0, k], "loading_logf": V[1, k], "loading_logcc": V[2, k]})
            rp = response_pca(sub, X)
            if rp is not None:
                rp.update({"truth_cc_mmhr": t, "series": s})
                resp_rows.append(rp)

    shr = pd.DataFrame(shr_rows)
    beh = pd.DataFrame(beh_rows)
    if shr.empty:
        sys.exit("No behavioral PCA could be computed (too few finite KGE_2012 values).")
    # is cc constrained beyond the null?
    flags = []
    for _, r in shr.iterrows():
        if r["series"] != SER_ON:
            flags.append(np.nan)
            continue
        o = shr[(shr.truth_cc_mmhr == r["truth_cc_mmhr"]) & (shr.series == SER_OFF) & (shr.top_frac == r["top_frac"])]
        flags.append(np.nan if o.empty else bool(r["ci_hi_cc"] < o["ci_lo_cc"].iloc[0]))
    shr["cc_constrained_vs_null"] = flags
    lo, hi = np.log10(BOX["channelconductivity_mmhr"][0]), np.log10(BOX["channelconductivity_mmhr"][1])
    pos = [(np.log10(t) - lo) / (hi - lo) if t > 0 else np.nan for t in shr["truth_cc_mmhr"]]
    shr["truth_box_pos"] = pos
    shr["near_box_edge"] = [bool(p < NEAR_EDGE or p > 1 - NEAR_EDGE) if np.isfinite(p) else np.nan for p in pos]

    main_frac = args.top_fracs[0]
    print_behavioral(shr, main_frac)
    print_directions(shr, main_frac)
    print_response(resp_rows)
    ft = pick_figure_truths(truths, args.figure_truths)
    print_param_metric(resp_rows, [t for t in ft if t != OFF_TRUTH] or ft)

    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ck.rows).to_csv(out_dir / "pca_input_checks_110.csv", index=False)
    shr.to_csv(out_dir / "pca_shrinkage_110.csv", index=False)
    beh.to_csv(out_dir / "pca_behavioral_110.csv", index=False)
    resp_out = []
    for r in resp_rows:
        for k in range(r["loadings"].shape[0]):
            row = {"truth_cc_mmhr": r["truth_cc_mmhr"], "series": r["series"], "component": "PC%d" % (k + 1),
                   "n_runs_used": r["n"], "explained_variance_ratio": r["evr"][k],
                   "rho_Ks": r["rho"][0, k], "rho_logf": r["rho"][1, k], "rho_logcc": r["rho"][2, k]}
            for j, m in enumerate(r["metrics"]):
                row["load_" + m] = r["loadings"][k, j]
            resp_out.append(row)
    pd.DataFrame(resp_out).to_csv(out_dir / "pca_response_110.csv", index=False)
    pm = []
    for r in resp_rows:
        for p, lab in enumerate(("Ks", "log_f", "log_cc")):
            row = {"truth_cc_mmhr": r["truth_cc_mmhr"], "series": r["series"], "parameter": lab, "n_runs_used": r["n"]}
            for j, m in enumerate(r["metrics"]):
                row[m] = r["rho_metric"][p, j]
            pm.append(row)
    pd.DataFrame(pm).to_csv(out_dir / "pca_param_metric_corr_110.csv", index=False)

    if not args.no_plots:
        try:
            fig_shrinkage(shr, args.top_fracs, out_dir / "fig_pca_shrinkage_110.png") if (shr.truth_cc_mmhr > 0).any() else None
            fig_loadings(beh, shr, ft, main_frac, out_dir / "fig_pca_loadings_110.png")
            fig_circle(resp_rows, ft, out_dir / "fig_pca_response_circle_110.png")
        except Exception as e:
            print("\n  (figures skipped: %s -- the CSVs are saved regardless)" % e)

    prov = {"script": Path(__file__).name, "created_local": datetime.now().isoformat(timespec="seconds"),
            "pandas": pd.__version__, "numpy": np.__version__,
            "inputs": {n: {"md5": md5_of(summary_dir / n)} for n in (STAGE2_NAME, CONTROL_NAME)},
            "rescored_long": ({"md5": md5_of(summary_dir / RESCORE_DIR / LONG_NAME)}
                              if (summary_dir / RESCORE_DIR / LONG_NAME).exists() else None),
            "top_fracs": args.top_fracs, "boot": args.boot, "seed": args.seed, "expect_n": args.expect_n,
            "truths": [float(t) for t in truths], "checks": ck.rows}
    (out_dir / "PROVENANCE_pca_cc_110.json").write_text(json.dumps(prov, indent=2, default=str))

    print("\nSaved to: %s" % out_dir)
    print("  pca_shrinkage_110.csv  pca_behavioral_110.csv  pca_response_110.csv  pca_param_metric_corr_110.csv  pca_input_checks_110.csv"
          "%s" % ("" if args.no_plots else "  + 3 figures"))
    if ck.n("FAIL"):
        print("\n  !! WARNING: this ran with %d FAILED input check(s) (--continue_on_fail). Do not trust or show these results."
              % ck.n("FAIL"))
    print("  Caveat: one storm, 250 points in 3-D, routing pinned at truth, noise-free truth -- a best case.")


if __name__ == "__main__":
    main()
