"""
analyze_series121_partials_and_dial.py
======================================
Series 121 follow-up analysis. It answers two questions that the first Series 121 script
(analyze_series121_cvrn_identifiability.py) could not:

  1. PARTIAL CORRELATIONS. cv looked like a weak parameter and n looked even weaker when each was
     correlated with the score on its own, because flowexp swings the score so much that it drowns
     everything else. A partial correlation asks: "with the other parameters held fixed, how strongly
     does THIS parameter move the score?" This script gives both views side by side, plus how much of the
     variance each parameter accounts for alone and uniquely.

  2. THE "DIAL" (single-index) TEST. Is the effect of the three routing parameters on the score really
     three separate effects, or one combined effect (one dial that is turned by a mix of cv, flowexp and n)?
     The test fits  score = g(w1*cv + w2*flowexp + w3*n)  (all three standardized, g a smooth curve) and
     compares it, on runs the model never saw (cross-validation), with a model that is free to treat the three
     parameters separately. If the one-dial model does as well, the three parameters act as one.
     It then checks whether the dial points the same way at every anchor, whether it agrees with the
     direction the best runs are pinned along (the PCA in the first script), whether different metrics
     use the same dial, and where the true parameter set sits on it.

It runs no tRIBS and changes none of your results files.

WHAT IT PRODUCES  (same folder as the first Series 121 script, with different file names)
-----------------------------------------------------------------------------------------
Tables (CSV, tidy):
    series121_partial_correlations.csv   scope x anchor x model x parameter: marginal r, partial r (with 95% interval),
                                         p-value, share of variance alone (marginal r squared) and unique share
    series121_model_r2.csv               scope x anchor x model: R squared (in-sample) of the linear models
    series121_dial_summary.csv           scope x anchor: the one-dial fit for KGE_2012: weights, held-out R squared of
                                         six models, dial values of the truth, natural-unit equivalences, the PCA check
    series121_dial_by_metric.csv         the same test for each metric (pooled independent runs)
    series121_dial_bins.csv              score against position on the dial (pooled independent runs)
    series121_timing_check.csv           per anchor: runs that score well, and how many of those also arrive on time
Figures (PNG):
    fig_series121_5_partials.png         alone vs with the others held fixed; per-anchor partial correlations
    fig_series121_6_dial.png             score vs the dial; held-out R squared; weights per anchor; weights per metric

HOW TO RUN (from the lab/ folder, like the other scripts)
---------------------------------------------------------
    python Series120_plots/analyze_series121_partials_and_dial.py
or, if you keep it directly in lab/:
    python analyze_series121_partials_and_dial.py

It takes about half a minute (the cross-validation refits the dial many times).

Options:
    --csv PATH           the results file (default calibration_work/03_comparisons/summary_tables/
                         lhs_results_anchor_ALL_121.csv, looked for next to lab/)
    --out_dir PATH       where to write (default calibration_work/03_comparisons/sensitivity_plots/series121_cvrn/)
    --top_frac X         share of each anchor's runs (best KGE_2012) used for the pinned-direction check (default 0.20)
    --good_kge X         a run "scores well" at or above this KGE_2012 (default 0.95; a choice, not a standard)
    --arrival_tol_min X  "arrives on time" = first-arrival error within +/- X minutes (default 5, the metric's resolution)
    --dial_bin X         width of the dial bins, in standardized units (default 0.5)
    --truth_window X     "near the truth's dial value" = within +/- X standardized units (default 0.25)
    --cv_reps N          how many times the cross-validation is repeated with different random splits (default 3)
    --cv_folds N         folds in each cross-validation (default 5)
    --seed N             random seed for the splits (default 0; the same seed gives the same numbers)
    --min_runs N         fewest runs a dial fit may use (default 40; the full cubic model has 20 coefficients)
    --dpi N              figure resolution (default 150)
    --no_plots           write the tables only

READ THIS BEFORE QUOTING ANY NUMBER
-----------------------------------
  * One synthetic truth, one storm (12 Aug 2014, 1.0x), channel loss off, OPTSNOW 0, five anchors.
  * "Held-out R squared" = the share of the variance in runs that were left out of the fit. It cannot be
    inflated by a model that has more terms, so a flexible model only wins if it truly captures more.
  * The dial is a description of these runs, within this sampling box (cv 2.5 to 6.5, flowexp 0.18 to 0.35,
    n 0.02 to 0.10). The truth n sits very close to the lower edge of that box, and the box is not centred on the
    truth for flowexp either, so the weights are for this box. A different box could give different weights.
  * Standardized weights are in units of each parameter's own spread in the sample, so they depend on the box.
    The natural-unit equivalences (printed in section 3) translate them back.
  * A dial that explains the score tells you the parameters trade off with each other. It does NOT tell you why
    inside tRIBS; that is a hypothesis to check in the model code.
  * Partial correlations hold the OTHER ROUTING parameters fixed. Ks and f are held fixed by the design.
  * p-values are not corrected for the number of tests. Use them to rank and to flag, not to claim.
"""

import argparse
import itertools
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize

# ======================================================================
# SETTINGS (same values as the Series 121 run script and the first analysis script)
# ======================================================================
ANCHOR_ORDER = ["f0p004", "f0p008", "f0p012_true", "f0p015", "f0p02"]     # low f -> high f
TRUE_ANCHOR = "f0p012_true"

PARAMS = ["kinemvelcoef", "flowexp", "channelroughness"]
PARAM_SHORT = {"kinemvelcoef": "cv", "flowexp": "flowexp", "channelroughness": "n"}
PARAM_LONG = {"kinemvelcoef": "cv (kinemvelcoef)", "flowexp": "flowexp (r)", "channelroughness": "n (channelroughness)"}
SHORT = [PARAM_SHORT[p] for p in PARAMS]
BOX = {
    "kinemvelcoef":     dict(lo=2.5,  hi=6.5,  truth=4.5),
    "flowexp":          dict(lo=0.18, hi=0.35, truth=0.24),
    "channelroughness": dict(lo=0.02, hi=0.10, truth=0.026),
}
BOX_TOL = 1e-6
TRUTH = np.array([BOX[p]["truth"] for p in PARAMS])
I_CV, I_R, I_N = 0, 1, 2                      # positions of the parameters in PARAMS

# metrics tested one by one in the dial-by-metric section (column, label)
DIAL_METRICS = [
    ("kge_2012", "KGE_2012"),
    ("peak_error_pct", "peak error %"),
    ("peak_timing_error_hr", "peak timing err (h)"),
    ("first_arrival_error_min", "first arrival err (min)"),
    ("rising_limb_steepness_ratio", "rising-limb steepness"),
    ("time_to_peak_from_exc_min", "time to peak (min)"),
    ("duration_above_thresh_error_min", "duration above thr. err (min)"),
    ("recession_rate_ratio", "recession rate ratio"),
    ("pbias_pct", "PBIAS %"),
]
METRIC_LABEL = dict(DIAL_METRICS)

# the six models compared by cross-validation
MODEL_LABELS = [
    ("single", "ONE dial of all three (cubic g)"),
    ("full", "full cubic in all three (20 terms)"),
    ("r_only", "flowexp alone (cubic)"),
    ("r_n", "flowexp + n (cubic)"),
    ("r_cv", "flowexp + cv (cubic)"),
    ("cv_n", "cv + n, no flowexp (cubic)"),
]
MODEL_LABEL = dict(MODEL_LABELS)
DEG = 3

DEFAULT_CSV_NAME = "lhs_results_anchor_ALL_121.csv"
OUT_SUBFOLDER = "series121_cvrn"
FILES = dict(
    partial="series121_partial_correlations.csv",
    r2="series121_model_r2.csv",
    dial="series121_dial_summary.csv",
    bymetric="series121_dial_by_metric.csv",
    bins="series121_dial_bins.csv",
    timing="series121_timing_check.csv",
    fig5="fig_series121_5_partials.png",
    fig6="fig_series121_6_dial.png",
)

# colours (the light-mode surface and ink of the other Series 120/121 figures)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
GREY_DOT = "#b4b3ab"
GOLD = "#f0a30a"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
PARAM_COLOR = {"kinemvelcoef": BLUE, "flowexp": ORANGE, "channelroughness": AQUA}   # palette slots 1, 2, 3


class DataProblem(Exception):
    """Something is wrong with the input file. Printed as a plain message, without a traceback."""


# ======================================================================
# Loading and checking the results file (same checks as the first Series 121 script)
# ======================================================================
def load_results(path):
    path = Path(path)
    if not path.exists():
        raise DataProblem(f"Could not find {path}")
    try:
        df = pd.read_csv(path)
    except Exception as exc:                                  # noqa: BLE001  (any read problem -> plain message)
        raise DataProblem(f"Could not read {path} as a CSV: {exc}")
    print(f"Read {len(df)} rows and {df.shape[1]} columns from {path.name}")

    need = ["run_id", "anchor_label", "batch", "draw_idx", "seed"] + PARAMS + ["kge_2012"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise DataProblem("The file is missing these column(s): " + ", ".join(missing)
                          + ". Is this the Series 121 ALL file (lhs_results_anchor_ALL_121.csv)?")
    return df


def check_and_clean(df):
    """Run the checks. Returns (clean frame, dict of flags). Stops (DataProblem) only if the data cannot be used."""
    print("\nChecks on the input")
    flags = {"pool_ok": True}

    for c in PARAMS + ["kge_2012"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    bad = df[PARAMS + ["kge_2012"]].isna().any(axis=1)
    if bad.any():
        print(f"  NOTE  {int(bad.sum())} row(s) have a missing or non-numeric parameter or KGE_2012; they are left out.")
        df = df[~bad].copy()
    if df.empty:
        raise DataProblem("No usable rows are left.")
    print(f"  ok    {len(df)} usable rows")

    if df["run_id"].duplicated().any():
        dup = df.loc[df["run_id"].duplicated(keep=False), "run_id"].unique()
        raise DataProblem(f"{len(dup)} run_id value(s) appear more than once (for example {dup[0]}). "
                          "Counting a run twice would distort every number; check the ALL file.")
    print("  ok    every run_id is unique")

    if "series" in df.columns and not (df["series"] == 121).all():
        print(f"  WARN  the 'series' column is not 121 everywhere: {sorted(df['series'].unique().tolist())}")
    if "storm_leg" in df.columns and df["storm_leg"].nunique() > 1:
        print(f"  WARN  more than one storm in the file: {sorted(df['storm_leg'].unique().tolist())}")
    if "optpercolation" in df.columns and not (df["optpercolation"] == 0).all():
        print("  WARN  optpercolation is not 0 everywhere (channel loss was on for some runs)")

    present = [a for a in ANCHOR_ORDER if a in set(df["anchor_label"])]
    extra = sorted(set(df["anchor_label"]) - set(ANCHOR_ORDER))
    if extra:
        print(f"  WARN  anchor label(s) not recognised and left out: {extra}")
        df = df[df["anchor_label"].isin(ANCHOR_ORDER)].copy()
    lost = [a for a in ANCHOR_ORDER if a not in present]
    if lost:
        print(f"  WARN  no runs for anchor(s): {lost}")
    if len(present) < 1:
        raise DataProblem("No recognised anchors are present.")
    flags["anchors"] = present
    flags["has_truth_anchor"] = TRUE_ANCHOR in present

    if not set(df["batch"]) <= {"shared", "indep"}:
        raise DataProblem(f"The 'batch' column should hold 'shared' or 'indep', found {sorted(df['batch'].unique().tolist())}")
    counts = df.groupby(["anchor_label", "batch"]).size().unstack(fill_value=0).reindex(present)
    for col in ("shared", "indep"):
        if col not in counts:
            counts[col] = 0
    print("  rows per anchor:   " + "   ".join(f"{a} {int(counts.loc[a, 'shared'])}+{int(counts.loc[a, 'indep'])}" for a in present)
          + "   (shared+independent)")
    odd = [a for a in present if counts.loc[a, "shared"] != 50 or counts.loc[a, "indep"] != 50]
    if odd:
        print(f"  WARN  not 50+50 at: {odd}; the numbers are still computed from what is there.")

    for p in PARAMS:
        b = BOX[p]
        outside = ((df[p] < b["lo"] - BOX_TOL) | (df[p] > b["hi"] + BOX_TOL)).sum()
        if outside:
            print(f"  WARN  {int(outside)} run(s) have {PARAM_SHORT[p]} outside the sampling box {b['lo']} to {b['hi']}")
    print("  ok    routing parameters checked against the sampling box")

    # the independent draws may be pooled only if they really are independent
    ind = df[df["batch"] == "indep"]
    seeds = ind.groupby("anchor_label")["seed"].agg(lambda s: tuple(sorted(set(s))))
    all_seeds = [s for ss in seeds.values for s in ss]
    seeds_distinct = len(all_seeds) == len(set(all_seeds))
    dup_sets = int(ind.duplicated(PARAMS).sum())
    if ind["anchor_label"].nunique() < 2:
        flags["pool_ok"] = False
        print("  NOTE  fewer than two anchors have independent draws; no pooled numbers are computed.")
    elif not seeds_distinct or dup_sets:
        flags["pool_ok"] = False
        print(f"  WARN  independent draws are not independent (shared seeds: {not seeds_distinct}; "
              f"repeated parameter sets: {dup_sets}); pooled numbers are skipped.")
    else:
        print("  ok    independent draws use a different seed at every anchor and no parameter set repeats;"
              " pooling them is fair")
    return df, flags


# ======================================================================
# Part 1: partial correlations and variance shares
# ======================================================================
def design(cols):
    """Design matrix with an intercept."""
    return np.column_stack([np.ones(len(cols[0]))] + list(cols))


def sse_of(y, cols):
    X = design(cols) if len(cols) else np.ones((len(y), 1))
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    return float(r @ r), b


def r2_of(y, cols):
    sst = float(((y - y.mean()) ** 2).sum())
    return 1.0 - sse_of(y, cols)[0] / sst


def block_partial(y, block, controls):
    """Partial R squared of a block of terms given the controls, with the F-test p-value and the sign of the block's
    first term in the full model. Partial R squared = share of the variance the controls leave unexplained
    that the block then explains."""
    n = len(y)
    sse_r, _ = sse_of(y, controls)
    sse_f, b = sse_of(y, list(controls) + list(block))
    k_full = len(controls) + len(block)
    pr2 = 0.0 if sse_r <= 0 else max(0.0, (sse_r - sse_f) / sse_r)
    dof = n - k_full - 1
    if dof <= 0 or sse_f <= 0:
        return pr2, np.nan, np.nan, 1.0
    F = ((sse_r - sse_f) / len(block)) / (sse_f / dof)
    p = float(stats.f.sf(F, len(block), dof))
    sign = 1.0 if b[1 + len(controls)] >= 0 else -1.0
    return pr2, float(F), p, sign


def corr_ci(r, n_eff):
    if not np.isfinite(r) or abs(r) >= 1 or n_eff <= 3:
        return np.nan, np.nan
    z, se = np.arctanh(r), 1.0 / np.sqrt(n_eff - 3)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def partial_scopes(df, flags):
    """List of (scope, anchor_label, frame, centre). Raw pooled, pooled with anchor means removed, and each anchor."""
    out = []
    if flags["pool_ok"]:
        ind = df[df["batch"] == "indep"]
        out.append(("pooled_indep", "ALL_INDEP", ind, False))
        out.append(("pooled_indep_centered", "ALL_INDEP", ind, True))
    for a in flags["anchors"]:
        out.append(("anchor_all", a, df[df["anchor_label"] == a], False))
    return out


def centred_xy(sub, centre, ycol="kge_2012"):
    X = sub[PARAMS].to_numpy(float)
    y = sub[ycol].to_numpy(float)
    if centre:
        X = X - sub.groupby("anchor_label")[PARAMS].transform("mean").to_numpy(float)
        y = y - sub.groupby("anchor_label")[ycol].transform("mean").to_numpy(float)
    return X, y


def partials_and_variance(df, flags):
    """Returns (partial-correlation table, model-R2 table). Models: 'linear' = cv + flowexp + n;
    'linear_plus_flowexp_sq' = the same plus flowexp squared (flowexp acts non-linearly)."""
    prow, mrow = [], []
    for scope, anchor, sub, centre in partial_scopes(df, flags):
        X, y = centred_xy(sub, centre)
        n = len(y)
        if n < 12 or np.ptp(y) == 0:
            print(f"  NOTE  {scope} {anchor}: too few runs or no spread in KGE_2012; partial correlations skipped.")
            continue
        cv, r, nn = X[:, I_CV], X[:, I_R], X[:, I_N]
        r_sq = r ** 2
        cols = {"cv": cv, "flowexp": r, "n": nn}
        models = {
            "linear": {"terms": {"cv": [cv], "flowexp": [r], "n": [nn]}},
            "linear_plus_flowexp_sq": {"terms": {"cv": [cv], "flowexp": [r, r_sq], "n": [nn]}},
        }
        # R squared of the simple models
        for name, cs in (("flowexp alone", [r]), ("cv alone", [cv]), ("n alone", [nn]),
                         ("flowexp + cv", [r, cv]), ("flowexp + n", [r, nn]), ("cv + n", [cv, nn]),
                         ("cv + flowexp + n (linear)", [cv, r, nn]),
                         ("cv + flowexp + flowexp squared + n", [cv, r, r_sq, nn])):
            mrow.append(dict(scope=scope, anchor_label=anchor, n=n, model=name, r2_in_sample=r2_of(y, cs)))
        # linear-model sign of each parameter (used to sign a block)
        for mname, spec in models.items():
            terms = spec["terms"]
            all_cols = [c for t in terms.values() for c in t]
            r2_full = r2_of(y, all_cols)
            for pname, p in zip(("cv", "flowexp", "n"), PARAMS):
                block = terms[pname]
                controls = [c for k, t in terms.items() if k != pname for c in t]
                pr2, F, pv, sign = block_partial(y, block, controls)
                if len(block) > 1:                                  # flowexp with its square: take the sign from the linear model
                    others = [c for k, t in models["linear"]["terms"].items() if k != pname for c in t]
                    sign = block_partial(y, models["linear"]["terms"][pname], others)[3]
                r_partial = sign * np.sqrt(pr2)
                lo, hi = (corr_ci(r_partial, n - len(controls)) if len(block) == 1 else (np.nan, np.nan))
                marg = corr_pair_simple(cols[pname], y)
                r2_without = r2_of(y, controls)
                prow.append(dict(scope=scope, anchor_label=anchor, n=n, model=mname, parameter=p,
                                 parameter_short=pname, term=("flowexp and flowexp squared (block)" if len(block) > 1 else "single term"),
                                 marginal_r=marg, marginal_r2=marg ** 2 if np.isfinite(marg) else np.nan,
                                 partial_r=float(r_partial), partial_ci_lo=lo, partial_ci_hi=hi, partial_p=pv,
                                 partial_r2=pr2, r2_full=r2_full, r2_without=r2_without, unique_r2=r2_full - r2_without))
    return pd.DataFrame(prow), pd.DataFrame(mrow)


def corr_pair_simple(x, y):
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


# ======================================================================
# Part 2: the single-index ("dial") model
# ======================================================================
def poly1(t, deg=DEG):
    return np.column_stack([t ** k for k in range(deg + 1)])


def polyfull(Z, deg=DEG):
    """Polynomial of every degree up to deg in all columns of Z, with every cross term (20 columns for 3 variables, deg 3)."""
    cols = [np.ones(len(Z))]
    for d in range(1, deg + 1):
        for idx in itertools.combinations_with_replacement(range(Z.shape[1]), d):
            cols.append(np.prod(Z[:, idx], axis=1))
    return np.column_stack(cols)


def unit_vector(theta, phi):
    return np.array([np.cos(phi) * np.cos(theta), np.cos(phi) * np.sin(theta), np.sin(phi)])


def sse_dir(w, Z, y):
    X = poly1(Z @ w)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    return float(r @ r), b


def fit_single_index(Z, y, ngrid=(30, 15)):
    """Best direction w (length 1) for  y = g(Z @ w),  g a cubic. Coarse search over directions, then a local polish
    from the three best. w and -w are the same dial turned the other way, so the sign is fixed by convention:
    the flowexp weight is positive."""
    thetas = np.linspace(0, np.pi, ngrid[0], endpoint=False)
    phis = np.linspace(-np.pi / 2, np.pi / 2, ngrid[1])
    cand = sorted(((sse_dir(unit_vector(t, p), Z, y)[0], t, p) for t in thetas for p in phis))[:3]
    best = None
    for _, t0, p0 in cand:
        res = minimize(lambda a: sse_dir(unit_vector(a[0], a[1]), Z, y)[0], [t0, p0], method="Nelder-Mead",
                       options=dict(xatol=1e-4, fatol=1e-10, maxiter=400))
        if best is None or res.fun < best.fun:
            best = res
    w = unit_vector(*best.x)
    if w[I_R] < 0:
        w = -w
    sse, b = sse_dir(w, Z, y)
    return w, b, sse


def predict_single(w, b, Z):
    return poly1(Z @ w) @ b


def cv_scores(Z, y, reps, k, seed):
    """Repeated k-fold cross-validated R squared for the six models (held-out runs only)."""
    rng = np.random.default_rng(seed)
    n = len(y)
    sst = float(((y - y.mean()) ** 2).sum())
    out = {key: [] for key, _ in MODEL_LABELS}
    subsets = {"r_only": [I_R], "r_n": [I_R, I_N], "r_cv": [I_CV, I_R], "cv_n": [I_CV, I_N]}
    for _ in range(reps):
        perm = rng.permutation(n)
        folds = np.array_split(perm, k)
        sse = {key: 0.0 for key in out}
        for f in folds:
            tr = np.setdiff1d(perm, f)
            Ztr, ytr, Zte, yte = Z[tr], y[tr], Z[f], y[f]
            w, b, _ = fit_single_index(Ztr, ytr)
            sse["single"] += float(((yte - predict_single(w, b, Zte)) ** 2).sum())
            bb, *_ = np.linalg.lstsq(polyfull(Ztr), ytr, rcond=None)
            sse["full"] += float(((yte - polyfull(Zte) @ bb) ** 2).sum())
            for key, cols in subsets.items():
                b2, *_ = np.linalg.lstsq(polyfull(Ztr[:, cols]), ytr, rcond=None)
                sse[key] += float(((yte - polyfull(Zte[:, cols]) @ b2) ** 2).sum())
        for key in out:
            out[key].append(1.0 - sse[key] / sst)
    return {key: float(np.mean(v)) for key, v in out.items()}


def scope_matrix(sub, centre):
    """Standardized parameters of a scope. With centre=True each anchor's own mean is removed first (and the overall mean put
    back), so differences in the typical parameters between anchors cannot leak in. Returns Z, mu, sd (mu and sd in natural units)."""
    X = sub[PARAMS].to_numpy(float)
    if centre:
        X = X - sub.groupby("anchor_label")[PARAMS].transform("mean").to_numpy(float) + X.mean(axis=0)
    mu, sd = X.mean(axis=0), X.std(axis=0, ddof=1)
    return (X - mu) / sd, mu, sd


def scope_metric(sub, m, centre):
    s = pd.to_numeric(sub[m], errors="coerce")
    if centre:
        s = s - s.groupby(sub["anchor_label"]).transform("mean")
    return s.to_numpy(float)


def fmt_w(w):
    return "  ".join(f"{s} {v:+.2f}" for s, v in zip(SHORT, w))


def natural_units(w, sd):
    """Translate the dial weights back into natural units. 'per_unit' = how far one natural unit of each parameter moves the dial."""
    per_unit = w / sd
    out = {}
    if abs(per_unit[I_R]) > 1e-9:
        out["flowexp_offset_per_plus1_cv"] = float(-per_unit[I_CV] / per_unit[I_R])
        out["flowexp_offset_per_plus0p01_n"] = float(-per_unit[I_N] / per_unit[I_R] * 0.01)
        for i, p in enumerate(PARAMS):
            width = BOX[p]["hi"] - BOX[p]["lo"]
            out[f"span_{PARAM_SHORT[p]}_in_flowexp_units"] = float(width * abs(per_unit[i]) / abs(per_unit[I_R]))
    else:
        out["flowexp_offset_per_plus1_cv"] = out["flowexp_offset_per_plus0p01_n"] = np.nan
        for p in PARAMS:
            out[f"span_{PARAM_SHORT[p]}_in_flowexp_units"] = np.nan
    return out


def fit_scope(Z, y, args):
    """Fit the dial and score the six models. Returns a dict, or None if the scope is too small or flat."""
    ok = np.isfinite(y)
    Z, y = Z[ok], y[ok]
    if len(y) < args.min_runs or np.ptp(y) == 0:
        return None
    w, b, sse = fit_single_index(Z, y)
    sst = float(((y - y.mean()) ** 2).sum())
    bb, *_ = np.linalg.lstsq(polyfull(Z), y, rcond=None)
    sse_full = float(((y - polyfull(Z) @ bb) ** 2).sum())
    cv = cv_scores(Z, y, args.cv_reps, args.cv_folds, args.seed)
    s = Z @ w
    rho = float(stats.spearmanr(s, y)[0])
    return dict(n=len(y), w=w, b=b, r2_in_single=1 - sse / sst, r2_in_full=1 - sse_full / sst, cv=cv,
                spearman_with_dial=rho, s_min=float(s.min()), s_max=float(s.max()))


def dial_analysis(df, flags, args):
    """Fits for KGE_2012 in the pooled centered scope and in each anchor; for each metric in the pooled centered scope."""
    t0 = time.time()
    summary_rows, fits = [], {}
    scopes = []
    if flags["pool_ok"]:
        scopes.append(("pooled_indep_centered", "ALL_INDEP", df[df["batch"] == "indep"], True))
    for a in flags["anchors"]:
        scopes.append(("anchor_all", a, df[df["anchor_label"] == a], False))
    for scope, anchor, sub, centre in scopes:
        Z, mu, sd = scope_matrix(sub, centre)
        y = scope_metric(sub, "kge_2012", centre)
        fit = fit_scope(Z, y, args)
        if fit is None:
            print(f"  NOTE  {scope} {anchor}: fewer than --min_runs ({args.min_runs}) runs, or no spread; no dial fit.")
            continue
        w = fit["w"]
        nat = natural_units(w, sd)
        s_truth = float(((TRUTH - mu) / sd) @ w)
        ts = np.linspace(fit["s_min"], fit["s_max"], 400)
        g_peak = float(ts[np.argmax(poly1(ts) @ fit["b"])])
        row = dict(scope=scope, anchor_label=anchor, n=fit["n"], w_cv=w[I_CV], w_flowexp=w[I_R], w_n=w[I_N],
                   r2_in_sample_single=fit["r2_in_single"], r2_in_sample_full=fit["r2_in_full"],
                   **{f"cv_r2_{k}": v for k, v in fit["cv"].items()},
                   truth_dial_value=s_truth, dial_min=fit["s_min"], dial_max=fit["s_max"], g_highest_at=g_peak,
                   spearman_dial_vs_kge=fit["spearman_with_dial"], **nat)
        summary_rows.append(row)
        fits[(scope, anchor)] = dict(fit=fit, mu=mu, sd=sd, row=row)
        print(f"  fitted the dial for {scope} {anchor if anchor != 'ALL_INDEP' else ''}".rstrip() + f"   ({time.time() - t0:.0f}s)")
    summary = pd.DataFrame(summary_rows)

    # the pinned direction (smallest-spread direction among the top runs) vs the fitted dial, per anchor
    pc_rows = []
    for a in flags["anchors"]:
        key = ("anchor_all", a)
        if key not in fits:
            continue
        sub = df[df["anchor_label"] == a]
        Za, _, _ = scope_matrix(sub, False)
        n_top = max(8, int(np.ceil(len(sub) * args.top_frac)))
        order = np.argsort(-sub["kge_2012"].to_numpy(float))[:n_top]
        if len(order) < 8:
            continue
        Xs = Za[order] - Za[order].mean(axis=0)
        _, S, Vt = np.linalg.svd(Xs, full_matrices=False)
        v3 = Vt[2]
        if v3[I_R] < 0:
            v3 = -v3
        pc_rows.append(dict(anchor_label=a, n_top_runs=len(order), pc3_variance_share=float(S[2] ** 2 / (S ** 2).sum()),
                            pc3_cv=v3[I_CV], pc3_flowexp=v3[I_R], pc3_n=v3[I_N],
                            abs_cosine_with_dial=abs(float(v3 @ fits[key]["fit"]["w"]))))
    pcs = pd.DataFrame(pc_rows)

    # the same test per metric (pooled independent runs, anchor means removed)
    metric_rows, metric_w = [], {}
    if flags["pool_ok"]:
        sub = df[df["batch"] == "indep"]
        Z, mu, sd = scope_matrix(sub, True)
        w_kge = fits[("pooled_indep_centered", "ALL_INDEP")]["fit"]["w"] if ("pooled_indep_centered", "ALL_INDEP") in fits else None
        for m, label in DIAL_METRICS:
            if m not in sub.columns:
                print(f"  NOTE  column {m} is not in the file; skipped in the by-metric test.")
                continue
            y = scope_metric(sub, m, True)
            ok = np.isfinite(y)
            fit = fit_scope(Z[ok], y[ok], args)
            if fit is None:
                print(f"  NOTE  {m}: too few runs or no spread; skipped in the by-metric test.")
                continue
            w = fit["w"]
            metric_w[m] = w
            metric_rows.append(dict(metric=m, metric_label=label, n=fit["n"], w_cv=w[I_CV], w_flowexp=w[I_R], w_n=w[I_N],
                                    **{f"cv_r2_{k}": v for k, v in fit["cv"].items()},
                                    spearman_metric_vs_dial=fit["spearman_with_dial"],
                                    abs_cosine_with_kge2012_dial=(abs(float(w @ w_kge)) if w_kge is not None else np.nan)))
        print(f"  finished the by-metric test   ({time.time() - t0:.0f}s)")
    bymetric = pd.DataFrame(metric_rows)

    # score against position on the dial (raw runs, pooled independent)
    bins, near_truth = pd.DataFrame(), {}
    key = ("pooled_indep_centered", "ALL_INDEP")
    if key in fits:
        f = fits[key]
        ind = df[df["batch"] == "indep"]
        s_run = ((ind[PARAMS].to_numpy(float) - f["mu"]) / f["sd"]) @ f["fit"]["w"]
        k = ind["kge_2012"].to_numpy(float)
        edges = np.arange(np.floor(s_run.min() / args.dial_bin) * args.dial_bin,
                          np.ceil(s_run.max() / args.dial_bin) * args.dial_bin + args.dial_bin / 2, args.dial_bin)
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (s_run >= lo) & (s_run < hi) if hi < edges[-1] else (s_run >= lo) & (s_run <= hi)
            if not m.any():
                continue
            kk = k[m]
            rows.append(dict(dial_from=float(lo), dial_to=float(hi), n_runs=int(m.sum()), kge2012_median=float(np.median(kk)),
                             kge2012_q25=float(np.quantile(kk, 0.25)), kge2012_q75=float(np.quantile(kk, 0.75)),
                             n_good=int((kk >= args.good_kge).sum()), share_good=float((kk >= args.good_kge).mean())))
        bins = pd.DataFrame(rows)
        s_truth = f["row"]["truth_dial_value"]
        near = np.abs(s_run - s_truth) <= args.truth_window
        near_truth = dict(s_truth=s_truth, n_near=int(near.sum()),
                          n_near_good=int((k[near] >= args.good_kge).sum()) if near.any() else 0,
                          share_good=float((k[near] >= args.good_kge).mean()) if near.any() else np.nan,
                          s_run=s_run, kge=k)
    return summary, pcs, bymetric, bins, near_truth, metric_w, fits


# ======================================================================
# Part 3: do the good runs also arrive on time?
# ======================================================================
def timing_check(df, flags, args):
    col = "first_arrival_error_min"
    if col not in df.columns:
        print(f"  NOTE  column {col} is not in the file; the timing check is skipped.")
        return pd.DataFrame()
    rows = []
    for a in flags["anchors"]:
        sub = df[df["anchor_label"] == a]
        err = pd.to_numeric(sub[col], errors="coerce")
        good = sub["kge_2012"] >= args.good_kge
        on_time = err.abs() <= args.arrival_tol_min
        both = good & on_time
        row = dict(anchor_label=a, n_runs=len(sub), n_arrival_missing=int(err.isna().sum()),
                   n_good=int(good.sum()), n_good_and_on_time=int(both.sum()),
                   share_of_good_on_time=(float(both.sum() / good.sum()) if good.sum() else np.nan),
                   arrival_tolerance_min=args.arrival_tol_min, good_kge_threshold=args.good_kge)
        for label, mask in (("good", good), ("good_on_time", both)):
            for p in PARAMS:
                vals = sub.loc[mask, p]
                row[f"{label}_{PARAM_SHORT[p]}_min"] = float(vals.min()) if len(vals) else np.nan
                row[f"{label}_{PARAM_SHORT[p]}_max"] = float(vals.max()) if len(vals) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ======================================================================
# Printing
# ======================================================================
def fnum(x, w=7, d=3, sign=True):
    if x is None or not np.isfinite(x):
        return "n/a".rjust(w)
    return (f"{x:+.{d}f}" if sign else f"{x:.{d}f}").rjust(w)


def fpv(p, w=9):
    return "n/a".rjust(w) if (p is None or not np.isfinite(p)) else f"{p:.1e}".rjust(w)


def sec(title):
    print("\n" + "=" * 86)
    print(title)
    print("=" * 86)


def scope_name(scope, anchor):
    return {"pooled_indep": "pooled independent", "pooled_indep_centered": "  anchor means removed"}.get(scope, anchor)


def print_report(df, flags, partial, r2tab, dial, pcs, bymetric, bins, near_truth, timing, args):
    present = flags["anchors"]

    # ---- 1 ----
    sec("1. One parameter at a time vs with the others held fixed (KGE_2012; linear model: cv + flowexp + n)")
    print("  'alone' = ordinary Pearson r between the parameter and the score. 'held fixed' = partial correlation: the same")
    print("  thing after the other two routing parameters have been accounted for. A big jump from 'alone' to 'held fixed'")
    print("  means the parameter matters, but another parameter's swings were hiding it.")
    if partial.empty:
        print("  (nothing computed: see the checks above.)")
    else:
        lin = partial[partial.model == "linear"]
        print(f"\n  {'':<24}{'n':>5}" + "".join(f"{PARAM_SHORT[p]:>26}" for p in PARAMS))
        print(f"  {'':<24}{'':>5}" + "".join(f"{'alone':>12}{'held fixed':>12}  " for _ in PARAMS))
        for (scope, anchor), g in lin.groupby(["scope", "anchor_label"], sort=False):
            line = f"  {scope_name(scope, anchor):<24}{int(g.n.iloc[0]):>5d}"
            for p in PARAMS:
                r = g[g.parameter == p].iloc[0]
                line += f"{fnum(r.marginal_r, 12)}{fnum(r.partial_r, 12)}  "
            print(line)
        print("\n  95% intervals (Fisher z) for the held-fixed correlation, pooled independent runs:")
        g = lin[lin.scope == "pooled_indep"]
        for p in PARAMS:
            r = g[g.parameter == p]
            if len(r):
                r = r.iloc[0]
                print(f"      {PARAM_SHORT[p]:<8} {r.partial_r:+.3f}   [{r.partial_ci_lo:+.3f}, {r.partial_ci_hi:+.3f}]   p = {fpv(r.partial_p).strip()}")
        print("\n  If flowexp acts non-linearly, the straight-line partial for the others is slightly too small. Adding flowexp squared:")
        q = partial[partial.model == "linear_plus_flowexp_sq"]
        print(f"  {'':<24}" + "".join(f"{PARAM_SHORT[p]:>12}" for p in PARAMS) + "     (flowexp = block of flowexp and flowexp squared)")
        for (scope, anchor), g in q.groupby(["scope", "anchor_label"], sort=False):
            line = f"  {scope_name(scope, anchor):<24}"
            for p in PARAMS:
                r = g[g.parameter == p].iloc[0]
                line += f"{fnum(r.partial_r, 12)}"
            print(line)

    # ---- 2 ----
    sec("2. How much of the score's variance each parameter accounts for: alone vs uniquely")
    print("  'alone' = r squared of the parameter by itself. 'unique' = how much R squared is lost if that parameter is removed from the")
    print("  model that has all three (and flowexp squared). The unique shares need not add up to the total, because parameters that")
    print("  trade off share some of the explained variance.")
    if not partial.empty:
        q = partial[partial.model == "linear_plus_flowexp_sq"]
        print(f"\n  {'':<24}{'n':>5}{'model R2':>10}" + "".join(f"{PARAM_SHORT[p] + ' alone':>14}{'unique':>8}" for p in PARAMS))
        for (scope, anchor), g in q.groupby(["scope", "anchor_label"], sort=False):
            line = f"  {scope_name(scope, anchor):<24}{int(g.n.iloc[0]):>5d}{g.r2_full.iloc[0]:>10.3f}"
            for p in PARAMS:
                r = g[g.parameter == p].iloc[0]
                line += f"{r.marginal_r2:>14.3f}{r.unique_r2:>8.3f}"
            print(line)
        print("  (For flowexp, 'unique' removes flowexp and flowexp squared together.)")
        print("\n  In-sample R squared of simpler models:")
        if "pooled_indep" in set(r2tab.scope):
            cols_ = [(("pooled_indep", "ALL_INDEP"), "pooled"), (("pooled_indep_centered", "ALL_INDEP"), "means removed")]
        else:
            cols_ = [(("anchor_all", a), a) for a in present]
        cols_ = [(k_, nme) for k_, nme in cols_ if ((r2tab.scope == k_[0]) & (r2tab.anchor_label == k_[1])).any()]
        print(f"  {'model':<40}" + "".join(f"{nme:>15}" for _, nme in cols_))
        for mname in r2tab.model.unique():
            line = f"  {mname:<40}"
            for (sc_, an_), _ in cols_:
                v = r2tab[(r2tab.scope == sc_) & (r2tab.anchor_label == an_) & (r2tab.model == mname)].r2_in_sample
                line += f"{v.iloc[0]:>15.3f}" if len(v) else f"{'n/a':>15}"
            print(line)

    # ---- 3 ----
    sec("3. The dial: can ONE combination of cv, flowexp and n stand in for all three? (KGE_2012)")
    print("  Held-out R squared = the share of the variance in runs left out of the fit (repeated cross-validation, "
          f"{args.cv_reps} x {args.cv_folds}-fold).")
    if dial.empty:
        print("  (no dial fit: see the notes above.)")
    else:
        short_col = {"single": "ONE dial", "full": "full cubic", "r_only": "flowexp", "r_n": "fl+n", "r_cv": "fl+cv", "cv_n": "cv+n"}
        print(f"\n  {'':<24}{'n':>5}" + "".join(f"{v:>11}" for v in short_col.values()))
        for r in dial.itertuples():
            line = f"  {scope_name(r.scope, r.anchor_label):<24}{r.n:>5d}"
            for k in short_col:
                line += f"{getattr(r, 'cv_r2_' + k):>11.3f}"
            print(line)
        print("  ONE dial = score depends on a single mix of the three parameters. Full cubic = the three parameters free to act separately")
        print("  (20 coefficients). fl = flowexp, n, cv = the other two. If ONE dial matches full cubic, the three act as one.")
        print("\n  Fitted dial weights on the standardized parameters (flowexp positive by convention; cv and n are relative to it):")
        for r in dial.itertuples():
            print(f"    {scope_name(r.scope, r.anchor_label):<24}cv {r.w_cv:+.2f}   flowexp {r.w_flowexp:+.2f}   n {r.w_n:+.2f}")
        pool = dial[dial.scope == "pooled_indep_centered"]
        ref = pool.iloc[0] if len(pool) else dial.iloc[0]
        print(f"\n  In natural units ({scope_name(ref.scope, ref.anchor_label).strip()} fit; a change in one parameter that leaves the dial where it was):")
        print(f"    raising cv by 1 is cancelled by changing flowexp by {ref.flowexp_offset_per_plus1_cv:+.4f}")
        print(f"    raising n by 0.01 is cancelled by changing flowexp by {ref.flowexp_offset_per_plus0p01_n:+.4f}")
        print("    Distance the dial moves across the whole sampled box, written as the equivalent change in flowexp:")
        print(f"      flowexp {ref.span_flowexp_in_flowexp_units:.3f}    cv {ref.span_cv_in_flowexp_units:.3f}    n {ref.span_n_in_flowexp_units:.3f}")
        print("    (a parameter whose span is small next to flowexp's has a smaller reach across this box, but it is still not zero)")
        print(f"\n  Where the truth sits and where the score peaks (standardized dial units):")
        print(f"    {'':<24}{'truth':>8}{'runs span':>20}{'score peaks at':>16}")
        for r in dial.itertuples():
            print(f"    {scope_name(r.scope, r.anchor_label):<24}{r.truth_dial_value:>+8.2f}{r.dial_min:>+10.2f} to {r.dial_max:<+6.2f}{r.g_highest_at:>+16.2f}")
        print("    'score peaks at' = where the AVERAGE score along the dial is highest. It need not equal the truth's position: runs at the same")
        print("    dial value still differ in the directions the dial does not capture, and the sampling box is not centred on the truth.")

    # ---- 4 ----
    sec("4. Score against position on the dial (pooled independent runs, raw KGE_2012)")
    if bins.empty:
        print("  (skipped: pooled numbers are not available.)")
    else:
        print(f"  Dial value of each run, from the pooled fit. 'good' = KGE_2012 >= {args.good_kge:g}.")
        print(f"    {'dial from':>10}{'to':>7}{'runs':>6}{'median':>8}{'q25-q75':>14}{'good':>7}{'share good':>12}")
        for r in bins.itertuples():
            print(f"    {r.dial_from:>+10.2f}{r.dial_to:>+7.2f}{r.n_runs:>6d}{r.kge2012_median:>8.3f}"
                  f"{r.kge2012_q25:>7.3f}-{r.kge2012_q75:<6.3f}{r.n_good:>7d}{r.share_good:>12.0%}")
        if near_truth:
            print(f"\n  Runs within +/-{args.truth_window:g} of the truth's dial value ({near_truth['s_truth']:+.2f}): {near_truth['n_near']} runs, "
                  f"{near_truth['n_near_good']} of them good ({near_truth['share_good']:.0%}).")
        print("  A run can sit at the truth's dial value with parameters far from the truth, because the parameters trade off.")

    # ---- 5 ----
    sec("5. The same test for each metric (pooled independent, anchor means removed): one dial or more?")
    if bymetric.empty:
        print("  (skipped: pooled numbers are not available.)")
    else:
        print(f"  {'metric':<32}{'ONE dial':>9}{'full cubic':>11}{'flowexp':>9}   {'weights: cv':>12}{'flowexp':>9}{'n':>7}   {'|cos| KGE dial':>15}")
        for r in bymetric.itertuples():
            print(f"  {r.metric_label:<32}{r.cv_r2_single:>9.3f}{r.cv_r2_full:>11.3f}{r.cv_r2_r_only:>9.3f}   "
                  f"{r.w_cv:>+12.2f}{r.w_flowexp:>+9.2f}{r.w_n:>+7.2f}   {fnum(r.abs_cosine_with_kge2012_dial, 15, 2, sign=False)}")
        print("  Columns ONE dial, full cubic, flowexp = held-out R squared. |cos| = how closely that metric's own dial points the same way as")
        print("  the KGE_2012 dial (1 = identical direction, 0 = at right angles). Metrics with a low R squared are not well described by")
        print("  any of these models; treat their weights as rough.")
        pair = ("peak_error_pct", "first_arrival_error_min")
        if set(pair) <= set(bymetric.metric):
            w1 = bymetric[bymetric.metric == pair[0]][["w_cv", "w_flowexp", "w_n"]].to_numpy(float)[0]
            w2 = bymetric[bymetric.metric == pair[1]][["w_cv", "w_flowexp", "w_n"]].to_numpy(float)[0]
            print(f"\n  Peak error and first arrival error: cosine between their dials = {abs(float(w1 @ w2)):.2f}. "
                  "Well below 1 means timing and shape are not turned by the same combination.")

    # ---- 6 ----
    sec(f"6. Does the dial point the same way as the direction the best runs are pinned along? (top {args.top_frac:.0%} of each anchor)")
    if pcs.empty:
        print("  (skipped: too few runs or no dial fit.)")
    else:
        print("  PC3 = the direction in which the top runs are LEAST spread out (the third principal component of the top runs).")
        print("  If the score depends on one combination, the best runs line up in a thin slab across it, so PC3 should match the dial.")
        print(f"    {'anchor':<14}{'top runs':>9}{'PC3 share':>11}   {'PC3 loadings (flowexp positive)':<36}{'|cosine| with dial':>20}")
        for r in pcs.itertuples():
            print(f"    {r.anchor_label:<14}{r.n_top_runs:>9d}{r.pc3_variance_share:>11.1%}   "
                  f"cv {r.pc3_cv:+.2f}  flowexp {r.pc3_flowexp:+.2f}  n {r.pc3_n:+.2f}      {r.abs_cosine_with_dial:>14.2f}")
        print("  With about 20 runs per anchor the PCA is a rough guide. A cosine near 1 at every anchor is the point.")

    # ---- 7 ----
    sec(f"7. Do the good runs also arrive on time? (KGE_2012 >= {args.good_kge:g} and first arrival within +/-{args.arrival_tol_min:g} min)")
    if timing.empty:
        print("  (skipped.)")
    else:
        print(f"  {'anchor':<14}{'runs':>6}{'good':>6}{'good + on time':>16}{'share on time':>15}")
        for r in timing.itertuples():
            print(f"  {r.anchor_label:<14}{r.n_runs:>6d}{r.n_good:>6d}{r.n_good_and_on_time:>16d}{fnum(r.share_of_good_on_time, 15, 2, sign=False)}")
        print("\n  Range of each parameter among the runs that score well, and among those that also arrive on time:")
        print(f"    {'anchor':<14}{'':<12}" + "".join(f"{PARAM_SHORT[p]:>20}" for p in PARAMS))
        for r in timing.itertuples():
            for label, nm in (("good", "good"), ("good_on_time", "+ on time")):
                line = f"    {r.anchor_label if label == 'good' else '':<14}{nm:<12}"
                for p in PARAMS:
                    lo, hi = getattr(r, f"{label}_{PARAM_SHORT[p]}_min"), getattr(r, f"{label}_{PARAM_SHORT[p]}_max")
                    line += (f"{lo:>9.3f} to {hi:<7.3f}" if np.isfinite(lo) else f"{'(no runs)':>20}")
                print(line)
        miss = int(timing["n_arrival_missing"].sum())
        if miss:
            print(f"  NOTE: {miss} run(s) have no first-arrival value; they count as 'not on time'.")
        print(f"  The first-arrival metric has {args.arrival_tol_min:g}-minute resolution at best, so a tolerance tighter than that is not meaningful.")
        print("  If requiring on-time arrival barely narrows the ranges, timing does not separate the good runs any further.")

    print("\nCAUTIONS: one synthetic truth, one storm; the dial is for THIS sampling box; p-values are not corrected for the number of tests;")
    print("          the reason the parameters trade off inside tRIBS is a hypothesis, not something these numbers show.")


# ======================================================================
# Figures
# ======================================================================
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "DejaVu Sans"
    return plt


def style_axes(ax, grid=True, axis="both"):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    if grid:
        ax.grid(True, color=GRID, linewidth=0.8, axis=axis)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_SOFT, labelsize=8.5)


def not_available(ax, msg="pooled numbers not available"):
    """Blank panel with a message (used when the pooled runs cannot be used)."""
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", color=INK_SOFT, fontsize=10)


def save_fig(fig, path, dpi):
    fig.savefig(path, dpi=dpi, facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"  saved {path.name}")


def figure_5(df, flags, partial, out_path, dpi):
    plt = _plt()
    from matplotlib.patches import Patch
    lin = partial[partial.model == "linear"]
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.3), gridspec_kw=dict(width_ratios=[1, 1, 1.7]))
    fig.patch.set_facecolor("white")
    have_pool = (lin.scope == "pooled_indep").any()
    pooled = lin[lin.scope == "pooled_indep"].set_index("parameter") if have_pool else None

    # A: alone vs held fixed (signed r)
    ax = axes[0]
    style_axes(ax, axis="y")
    if have_pool:
        x = np.arange(len(PARAMS))
        w = 0.36
        for i, p in enumerate(PARAMS):
            r = pooled.loc[p]
            ax.bar(i - w / 2, r.marginal_r, w, color=GREY_DOT, zorder=3)
            ax.bar(i + w / 2, r.partial_r, w, color=PARAM_COLOR[p], zorder=3)
            for xv, val in ((i - w / 2, r.marginal_r), (i + w / 2, r.partial_r)):
                ax.text(xv, val + (0.04 if val >= 0 else -0.04), f"{val:+.2f}", ha="center", va="bottom" if val >= 0 else "top",
                        fontsize=8.5, color=INK)
        ax.axhline(0, color=INK, linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([PARAM_SHORT[p] for p in PARAMS], fontsize=10)
        ax.set_ylim(-1.15, 1.15)
        ax.set_ylabel("correlation with KGE_2012", fontsize=9.5, color=INK_SOFT)
    else:
        not_available(ax)
    ax.set_title("Alone (grey) vs with the others\nheld fixed (colour); pooled independent runs", fontsize=10.5, color=INK, loc="left")

    # B: share of variance alone vs unique
    ax = axes[1]
    style_axes(ax, axis="y")
    q = partial[(partial.model == "linear_plus_flowexp_sq") & (partial.scope == "pooled_indep")].set_index("parameter")
    if len(q):
        w = 0.36
        for i, p in enumerate(PARAMS):
            r = q.loc[p]
            ax.bar(i - w / 2, r.marginal_r2, w, color=GREY_DOT, zorder=3)
            ax.bar(i + w / 2, r.unique_r2, w, color=PARAM_COLOR[p], zorder=3)
            for xv, val in ((i - w / 2, r.marginal_r2), (i + w / 2, r.unique_r2)):
                ax.text(xv, val + 0.015, f"{val:.2f}", ha="center", va="bottom", fontsize=8.5, color=INK)
        ax.set_xticks(range(len(PARAMS)))
        ax.set_xticklabels([PARAM_SHORT[p] for p in PARAMS], fontsize=10)
        ax.set_ylim(0, 1.0)
        ax.set_ylabel("share of the variance in KGE_2012", fontsize=9.5, color=INK_SOFT)
        ax.text(0.98, 0.97, f"all three together: {q.r2_full.iloc[0]:.2f}", transform=ax.transAxes, ha="right", va="top",
                fontsize=9, color=INK_SOFT)
    else:
        not_available(ax)
    ax.set_title("Variance accounted for:\nalone (grey) vs uniquely (colour)", fontsize=10.5, color=INK, loc="left")

    # C: per-anchor partials
    ax = axes[2]
    style_axes(ax, axis="y")
    groups = [(a, a) for a in flags["anchors"]]
    if have_pool:
        groups.append(("pooled_indep", "pooled"))
    n_g = len(groups)
    w = 0.26
    for gi, (key, lab) in enumerate(groups):
        sc = "pooled_indep" if key == "pooled_indep" else "anchor_all"
        an = "ALL_INDEP" if key == "pooled_indep" else key
        g = lin[(lin.scope == sc) & (lin.anchor_label == an)].set_index("parameter")
        if g.empty:
            continue
        for k, p in enumerate(PARAMS):
            val = g.loc[p, "partial_r"]
            ax.bar(gi + (k - 1) * w, val, w, color=PARAM_COLOR[p], zorder=3)
            ax.text(gi + (k - 1) * w, val + (0.03 if val >= 0 else -0.03), f"{val:+.2f}", ha="center", fontsize=7.3, color=INK,
                    va="bottom" if val >= 0 else "top", rotation=90)
    ax.axhline(0, color=INK, linewidth=0.8)
    if have_pool:
        ax.axvline(n_g - 1.5, color=AXIS, linewidth=1.0)
    ax.set_xticks(range(n_g))
    ax.set_xticklabels([("%s\n(truth)" % lab) if lab == TRUE_ANCHOR else lab for _, lab in groups], fontsize=8.5)
    ax.set_ylim(-1.3, 1.3)
    ax.set_ylabel("partial correlation with KGE_2012", fontsize=9.5, color=INK_SOFT)
    ax.legend(handles=[Patch(color=PARAM_COLOR[p], label=PARAM_LONG[p]) for p in PARAMS], loc="lower center", ncol=3,
              frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Partial correlation at each anchor (all 100 runs) and pooled", fontsize=10.5, color=INK, loc="left", pad=30)
    fig.suptitle("Series 121: do cv and n matter once flowexp is accounted for?", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "Linear model, partial correlation of each routing parameter with KGE_2012 holding the other two fixed. "
                         "Variance shares in the middle panel come from the model that also has flowexp squared.",
             ha="center", va="bottom", fontsize=8.3, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    save_fig(fig, out_path, dpi)


def figure_6(df, flags, dial, bymetric, bins, near_truth, args, out_path, dpi):
    plt = _plt()
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    fig, axes = plt.subplots(2, 2, figsize=(14.0, 10.2), gridspec_kw=dict(height_ratios=[1, 1.05]))
    fig.patch.set_facecolor("white")
    pool = dial[dial.scope == "pooled_indep_centered"]

    # A: score vs dial
    ax = axes[0][0]
    style_axes(ax)
    if near_truth and len(bins):
        s, k = near_truth["s_run"], near_truth["kge"]
        ax.scatter(s, k, s=16, color=GREY_DOT, alpha=0.9, linewidths=0, zorder=3)
        mid = (bins.dial_from + bins.dial_to) / 2
        ax.plot(mid, bins.kge2012_median, color=INK, linewidth=2.0, marker="o", markersize=4.5, zorder=4)
        ax.axvline(near_truth["s_truth"], color=GOLD, linewidth=2.0, zorder=2)
        ax.axhline(args.good_kge, color=INK_SOFT, linestyle=(0, (2, 3)), linewidth=1.0, zorder=2)
        ax.text(0.99, args.good_kge + 0.008, f"KGE_2012 = {args.good_kge:g}", fontsize=8.3, color=INK_SOFT, va="bottom", ha="right",
                transform=ax.get_yaxis_transform())
        ymin = max(0.0, np.floor((k.min() - 0.02) * 10) / 10)
        ax.set_ylim(ymin, 1.02)
        ax.text(near_truth["s_truth"], ymin + 0.02, " truth", color=INK, fontsize=9, va="bottom")
        ax.set_xlabel("position along the dial (standardized units; the flowexp end is positive)", fontsize=9, color=INK_SOFT)
        ax.set_ylabel("KGE_2012", fontsize=9.5, color=INK_SOFT)
        ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", color=GREY_DOT, markersize=6, label="runs (independent draws)"),
                           Line2D([0], [0], color=INK, linewidth=2, marker="o", markersize=4.5, label="median of each bin"),
                           Line2D([0], [0], color=GOLD, linewidth=2, label="the truth's position")],
                  loc="lower left", frameon=False, fontsize=8.3, bbox_to_anchor=(0.0, 0.07))
    else:
        not_available(ax)
    ax.set_title("Score against position on the dial", fontsize=10.5, color=INK, loc="left")

    # B: held-out R2 of the six models, pooled
    ax = axes[0][1]
    style_axes(ax, axis="x")
    if len(pool):
        r = pool.iloc[0]
        keys = [k for k, _ in MODEL_LABELS]
        vals = [r[f"cv_r2_{k}"] for k in keys]
        colors = [INK, "#8d8c85"] + [GREY_DOT] * 4
        ypos = np.arange(len(keys))[::-1]
        ax.barh(ypos, vals, 0.62, color=colors, zorder=3)
        for yv, v in zip(ypos, vals):
            ax.text(v + 0.01, yv, f"{v:.3f}", va="center", fontsize=9, color=INK)
        ax.set_yticks(ypos)
        ax.set_yticklabels([MODEL_LABEL[k] for k in keys], fontsize=9)
        ax.set_xlim(0, 1.12)
        ax.set_xlabel("held-out R squared (pooled independent runs, anchor means removed)", fontsize=9, color=INK_SOFT)
    else:
        not_available(ax)
    ax.set_title("Does one dial do as well as three separate parameters?", fontsize=10.5, color=INK, loc="left")

    # C: weights per scope
    ax = axes[1][0]
    style_axes(ax, axis="y")
    rows = list(dial.itertuples())
    w = 0.26
    for gi, r in enumerate(rows):
        for k, (p, val) in enumerate(zip(PARAMS, (r.w_cv, r.w_flowexp, r.w_n))):
            ax.bar(gi + (k - 1) * w, val, w, color=PARAM_COLOR[p], zorder=3)
            ax.text(gi + (k - 1) * w, val + (0.03 if val >= 0 else -0.03), f"{val:+.2f}", ha="center", fontsize=7.5, color=INK,
                    va="bottom" if val >= 0 else "top")
    ax.axhline(0, color=INK, linewidth=0.8)
    labels = []
    for r in rows:
        labels.append("pooled\n(anchor means\nremoved)" if r.scope == "pooled_indep_centered" else (r.anchor_label + ("\n(truth)" if r.anchor_label == TRUE_ANCHOR else "")))
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, fontsize=8.5)
    if len(rows) > 1 and rows[0].scope == "pooled_indep_centered":
        ax.axvline(0.5, color=AXIS, linewidth=1.0)
    ax.set_ylim(-1.2, 1.2)
    ax.set_ylabel("weight on the standardized parameter", fontsize=9.5, color=INK_SOFT)
    ax.legend(handles=[Patch(color=PARAM_COLOR[p], label=PARAM_LONG[p]) for p in PARAMS], loc="lower right", ncol=1, frameon=False, fontsize=8.5)
    ax.set_title("Is the dial turned by the same mix at every anchor?", fontsize=10.5, color=INK, loc="left")

    # D: weights per metric
    ax = axes[1][1]
    style_axes(ax, axis="x")
    if len(bymetric):
        n_m = len(bymetric)
        w = 0.26
        ylabels = []
        for mi, r in enumerate(bymetric.itertuples()):
            yc = n_m - 1 - mi
            for k, (p, val) in enumerate(zip(PARAMS, (r.w_cv, r.w_flowexp, r.w_n))):
                ax.barh(yc + (1 - k) * w, val, w, color=PARAM_COLOR[p], zorder=3)
            ylabels.append(f"{r.metric_label}  (R² {r.cv_r2_single:.2f})")
        ax.axvline(0, color=INK, linewidth=0.8)
        ax.set_yticks(range(n_m))
        ax.set_yticklabels(ylabels[::-1], fontsize=8.5)
        ax.set_xlim(-1.1, 1.1)
        ax.set_xlabel("dial weight (flowexp positive by convention); R² = held-out, one dial", fontsize=9, color=INK_SOFT)
        ax.legend(handles=[Patch(color=PARAM_COLOR[p], label=PARAM_SHORT[p]) for p in PARAMS], loc="lower left", ncol=3, frameon=False, fontsize=8.5)
    else:
        not_available(ax)
    ax.set_title("Does every metric use the same dial?", fontsize=10.5, color=INK, loc="left")
    fig.suptitle("Series 121: one dial, or three separate parameters?", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "One storm (1.0x), channel loss off, one sampling box. Held-out R squared comes from repeated cross-validation.",
             ha="center", va="bottom", fontsize=8.3, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.02, 1, 0.965))
    save_fig(fig, out_path, dpi)


# ======================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--csv", type=Path, default=None, help=f"the results file (default .../summary_tables/{DEFAULT_CSV_NAME})")
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--top_frac", type=float, default=0.20)
    ap.add_argument("--good_kge", type=float, default=0.95)
    ap.add_argument("--arrival_tol_min", type=float, default=5.0)
    ap.add_argument("--dial_bin", type=float, default=0.5)
    ap.add_argument("--truth_window", type=float, default=0.25)
    ap.add_argument("--cv_reps", type=int, default=3)
    ap.add_argument("--cv_folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min_runs", type=int, default=40)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()
    if not (0.05 <= args.top_frac <= 0.5):
        ap.error("--top_frac must be between 0.05 and 0.5.")
    if not (0.0 < args.good_kge < 1.0):
        ap.error("--good_kge must be between 0 and 1.")
    if args.arrival_tol_min <= 0:
        ap.error("--arrival_tol_min must be above 0.")
    if not (0.1 <= args.dial_bin <= 2.0):
        ap.error("--dial_bin must be between 0.1 and 2.")
    if not (0.05 <= args.truth_window <= 2.0):
        ap.error("--truth_window must be between 0.05 and 2.")
    if args.cv_reps < 1 or args.cv_reps > 50:
        ap.error("--cv_reps must be between 1 and 50.")
    if not (3 <= args.cv_folds <= 20):
        ap.error("--cv_folds must be between 3 and 20.")
    if args.min_runs < 25:
        ap.error("--min_runs must be at least 25 (the full cubic model has 20 coefficients).")
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")

    root = Path.cwd().parent
    calib_dir = root / "calibration_work"
    csv_path = args.csv or (calib_dir / "03_comparisons" / "summary_tables" / DEFAULT_CSV_NAME)
    out_dir = args.out_dir or (calib_dir / "03_comparisons" / "sensitivity_plots" / OUT_SUBFOLDER)
    if out_dir.resolve() == Path(csv_path).resolve().parent:
        print(f"The output folder {out_dir} is the folder holding the results file; refusing to write there.")
        return 1

    try:
        df = load_results(csv_path)
        df, flags = check_and_clean(df)
    except DataProblem as exc:
        print(f"\nSTOPPED: {exc}")
        if not Path(csv_path).exists():
            print("Run this from the lab/ folder, or point to the file with --csv.")
        return 1

    print("\nComputing (the cross-validation takes a little while)...")
    partial, r2tab = partials_and_variance(df, flags)
    dial, pcs, bymetric, bins, near_truth, metric_w, fits = dial_analysis(df, flags, args)
    timing = timing_check(df, flags, args)

    print_report(df, flags, partial, r2tab, dial, pcs, bymetric, bins, near_truth, timing, args)

    # attach the pinned-direction check to the dial summary so it is one table
    if len(dial) and len(pcs):
        extra = pcs.rename(columns={"abs_cosine_with_dial": "abs_cosine_pc3_with_dial"})
        extra = extra.assign(scope="anchor_all")
        dial = dial.merge(extra, on=["scope", "anchor_label"], how="left")

    out_dir.mkdir(parents=True, exist_ok=True)
    sec("Files written to " + str(out_dir))
    tables = [("partial", partial), ("r2", r2tab), ("dial", dial), ("bymetric", bymetric), ("bins", bins), ("timing", timing)]
    for key, table in tables:
        if table.empty:
            print(f"  (skipped {FILES[key]}: nothing to write)")
            continue
        table.to_csv(out_dir / FILES[key], index=False, float_format="%.6g")
        print(f"  saved {FILES[key]}  ({len(table)} rows)")

    if not args.no_plots:
        print("\nDrawing the figures...")
        if partial.empty:
            print(f"  (skipped {FILES['fig5']}: no partial correlations)")
        else:
            figure_5(df, flags, partial, out_dir / FILES["fig5"], args.dpi)
        if dial.empty:
            print(f"  (skipped {FILES['fig6']}: no dial fit)")
        else:
            figure_6(df, flags, dial, bymetric, bins, near_truth, args, out_dir / FILES["fig6"], args.dpi)

    print("\nDone. This script only reads the results file; it changes nothing in it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
