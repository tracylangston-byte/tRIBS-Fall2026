"""
analyze_series121_cvrn_identifiability.py
=========================================
Series 121 analysis: do the three routing parameters (cv, flowexp, n) change the score, and does
that change from one Ks/f anchor to the next? Reads the Series 121 results, writes tables and four
figures. It runs no tRIBS and changes none of your results files.

Adapted from the Series 101 script (analyze_series101_cvrn_identifiability.py, which stays as it is) for
the tRIBS 6.0.0 re-run. Differences from that script, all because of how Series 121 was designed:

  * Every anchor has 100 runs = 50 SHARED draws (the same 50 routing sets at every anchor, seed 42)
    + 50 INDEPENDENT draws (a different seed at each anchor). The results carry a 'batch' column.
  * The shared draws are NOT pooled across anchors. Pooling would count the same 50 numbers five
    times (pseudoreplication). They are used for what they are good for: comparing anchors, because
    the same routing set is run at every anchor.
  * Only the independent draws are pooled (250 runs). The pooled numbers are also given with the
    anchor means removed ("centered"), so a difference in typical score between anchors cannot
    leak into the pooled correlation.
  * The Series 101 PCA was underpowered and used shared draws. Here the PCA is run per anchor on the
    top runs and is labelled as the rough guide it is (20 runs by default).

WHAT IT PRODUCES
----------------
Tables (CSV, tidy, one row per thing counted):
    series121_correlation_summary.csv    scope x anchor x parameter x metric: Pearson r (with 95% interval),
                                         Spearman rho, p-values, n. Scopes: anchor_all (100 runs),
                                         anchor_shared (50), anchor_indep (50), pooled_indep (250),
                                         pooled_indep_centered (250, anchor means removed)
    series121_batch_consistency.csv      per anchor and parameter: Pearson r with KGE_2012 in the shared
                                         batch vs the independent batch, and a Fisher z test of whether they differ
    series121_anchor_summary.csv         per anchor: how the scores are spread, how many runs score well,
                                         where the best runs sit in routing-parameter space
    series121_shared_rank_stability.csv  the 50 shared draws: do they rank the same way at every anchor?
    series121_shared_vs_truth_anchor.csv the same routing set run at another anchor vs at the truth anchor
    series121_pca_summary.csv            PCA on the standardized (cv, flowexp, n) of each anchor's top runs
Figures (PNG):
    fig_series121_1_param_vs_kge.png         KGE_2012 against each routing parameter, one panel per anchor
    fig_series121_2_correlation_heatmaps.png Spearman rho: parameter x metric x anchor, plus the pooled row
    fig_series121_3_shared_draws.png         the same 50 routing sets at the five anchors; rank agreement
    fig_series121_4_top_runs_pca.png         where the best runs sit in (cv, flowexp, n); PCA loadings

HOW TO RUN (from the lab/ folder, like the other scripts)
---------------------------------------------------------
    python Series120_plots/analyze_series121_cvrn_identifiability.py
or, if you keep it directly in lab/:
    python analyze_series121_cvrn_identifiability.py

Options:
    --csv PATH         the results file (default calibration_work/03_comparisons/summary_tables/
                       lhs_results_anchor_ALL_121.csv, looked for next to lab/)
    --out_dir PATH     where to write (default calibration_work/03_comparisons/sensitivity_plots/series121_cvrn/)
    --top_frac X       share of each anchor's runs (best KGE_2012) used for "top runs" and the PCA (default 0.20)
    --min_n N          fewest top runs allowed for a PCA (default 8)
    --good_kge X       a run "scores well" at or above this KGE_2012 (default 0.95; a choice, not a standard)
    --p_flag X         p-values below this are marked with * and drawn bold (default 0.01)
    --dpi N            figure resolution (default 150)
    --no_plots         write the tables only

READ THIS BEFORE QUOTING ANY NUMBER
-----------------------------------
  * One synthetic truth, one storm (12 Aug 2014, 1.0x), channel loss off, OPTSNOW 0, five anchors.
  * Every correlation here is one parameter against one score. It does not show whether two parameters
    can trade off against each other; the top-run scatter and the PCA are a first look at that, with
    only about 20 runs per anchor.
  * The sampling box matters. The truth n (0.026) is close to the lower edge of the sampled n range
    (0.02 to 0.10), so "the best runs have low n" may partly be the box, not the model.
  * p-values are NOT corrected for the number of tests (several hundred are in the CSV). Use them to
    rank and to flag, not to claim.
  * 50 runs give a correlation an uncertainty of roughly +/-0.28 either way (95%). Differences between
    anchors, or between batches, smaller than that are not distinguishable from sampling noise.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ======================================================================
# SETTINGS (same values as the Series 121 run script)
# ======================================================================
ANCHOR_ORDER = ["f0p004", "f0p008", "f0p012_true", "f0p015", "f0p02"]     # low f -> high f
TRUE_ANCHOR = "f0p012_true"

PARAMS = ["kinemvelcoef", "flowexp", "channelroughness"]
PARAM_SHORT = {"kinemvelcoef": "cv", "flowexp": "flowexp", "channelroughness": "n"}
PARAM_LONG = {"kinemvelcoef": "cv (kinemvelcoef)", "flowexp": "flowexp (r)", "channelroughness": "n (channelroughness)"}
# sampling box and truth of each routing parameter
BOX = {
    "kinemvelcoef":     dict(lo=2.5,  hi=6.5,  truth=4.5),
    "flowexp":          dict(lo=0.18, hi=0.35, truth=0.24),
    "channelroughness": dict(lo=0.02, hi=0.10, truth=0.026),
}
BOX_TOL = 1e-6

# metric column -> label used in tables and figures
METRICS = [
    ("kge_2012", "KGE_2012"),
    ("kge", "KGE_2009"),
    ("pbias_pct", "PBIAS %"),
    ("peak_error_pct", "peak error %"),
    ("peak_timing_error_hr", "peak timing err (h)"),
    ("first_arrival_error_min", "first arrival err (min)"),
    ("rising_limb_steepness_ratio", "rising-limb steepness"),
    ("time_to_peak_from_exc_min", "time to peak (min)"),
    ("duration_above_thresh_error_min", "duration above thr. err (min)"),
    ("recession_rate_ratio", "recession rate ratio"),
]
METRIC_LABEL = dict(METRICS)
TOP_OVERLAP_K = 10          # "top 10 of the 50 shared draws"

DEFAULT_CSV_NAME = "lhs_results_anchor_ALL_121.csv"
OUT_SUBFOLDER = "series121_cvrn"
FILES = dict(
    corr="series121_correlation_summary.csv",
    batch="series121_batch_consistency.csv",
    anchors="series121_anchor_summary.csv",
    rank="series121_shared_rank_stability.csv",
    vs_truth="series121_shared_vs_truth_anchor.csv",
    pca="series121_pca_summary.csv",
    fig1="fig_series121_1_param_vs_kge.png",
    fig2="fig_series121_2_correlation_heatmaps.png",
    fig3="fig_series121_3_shared_draws.png",
    fig4="fig_series121_4_top_runs_pca.png",
)

# colours (the light-mode surface and ink of the other Series 120 figures)
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
RED = "#e34948"
NEUTRAL_MID = "#eeede8"
PARAM_COLOR = {"kinemvelcoef": BLUE, "flowexp": ORANGE, "channelroughness": AQUA}   # palette slots 1, 2, 3


class DataProblem(Exception):
    """Something is wrong with the input file. Printed as a plain message, without a traceback."""


# ======================================================================
# Loading and checking the results file
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
    flags = {"shared_ok": True, "pool_ok": True}
    notes = []

    # numbers must be numbers
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
                          "Counting a run twice would distort every correlation; check the ALL file.")
    print("  ok    every run_id is unique")

    if "series" in df.columns and not (df["series"] == 121).all():
        print(f"  WARN  the 'series' column is not 121 everywhere: {sorted(df['series'].unique().tolist())}")
    if "storm_leg" in df.columns and df["storm_leg"].nunique() > 1:
        print(f"  WARN  more than one storm in the file: {sorted(df['storm_leg'].unique().tolist())}")
    if "optpercolation" in df.columns and not (df["optpercolation"] == 0).all():
        print("  WARN  optpercolation is not 0 everywhere (channel loss was on for some runs)")

    # anchors present
    present = [a for a in ANCHOR_ORDER if a in set(df["anchor_label"])]
    extra = sorted(set(df["anchor_label"]) - set(ANCHOR_ORDER))
    if extra:
        print(f"  WARN  anchor label(s) not recognised and left out: {extra}")
        df = df[df["anchor_label"].isin(ANCHOR_ORDER)].copy()
    lost = [a for a in ANCHOR_ORDER if a not in present]
    if lost:
        print(f"  WARN  no runs for anchor(s): {lost}")
    if len(present) < 2:
        raise DataProblem("Fewer than two anchors are present, so there is nothing to compare.")
    flags["anchors"] = present
    flags["has_truth_anchor"] = TRUE_ANCHOR in present

    # batches
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
    flags["counts"] = counts

    # draws inside the box
    for p in PARAMS:
        b = BOX[p]
        outside = ((df[p] < b["lo"] - BOX_TOL) | (df[p] > b["hi"] + BOX_TOL)).sum()
        if outside:
            print(f"  WARN  {int(outside)} run(s) have {PARAM_SHORT[p]} outside the sampling box {b['lo']} to {b['hi']}")
    print("  ok    routing parameters checked against the sampling box")

    # shared draws really are shared
    sh = df[df["batch"] == "shared"]
    if sh.duplicated(["anchor_label", "draw_idx"]).any():
        raise DataProblem("Two shared runs have the same anchor and draw_idx; the paired comparison would be wrong.")
    n_sh_anchors = sh["anchor_label"].nunique()
    if n_sh_anchors < 2:
        flags["shared_ok"] = False
        print("  NOTE  fewer than two anchors have shared draws; the shared-draw comparison is skipped.")
    else:
        same = True
        for p in PARAMS:
            wide = sh.pivot(index="draw_idx", columns="anchor_label", values=p).dropna()
            if wide.empty or not (wide.nunique(axis=1) == 1).all():
                same = False
        if same:
            print("  ok    the shared draws carry identical (cv, flowexp, n) at every anchor")
        else:
            flags["shared_ok"] = False
            print("  WARN  the shared draws are NOT identical across anchors (or have no common draw_idx);"
                  " the shared-draw comparison is skipped.")

    # independent draws really are independent
    ind = df[df["batch"] == "indep"]
    seeds = ind.groupby("anchor_label")["seed"].agg(lambda s: tuple(sorted(set(s))))
    all_seeds = [s for ss in seeds.values for s in ss]
    seeds_distinct = len(all_seeds) == len(set(all_seeds))        # no seed used at two anchors
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
# Correlations
# ======================================================================
def corr_pair(x, y):
    """Pearson and Spearman for two arrays; NaN where it cannot be computed."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    out = dict(n=n, pearson_r=np.nan, pearson_p=np.nan, pearson_ci_lo=np.nan, pearson_ci_hi=np.nan,
               spearman_rho=np.nan, spearman_p=np.nan)
    if n < 4 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return out
    r, p = stats.pearsonr(x, y)
    rho, ps = stats.spearmanr(x, y)
    out.update(pearson_r=float(r), pearson_p=float(p), spearman_rho=float(rho), spearman_p=float(ps))
    if abs(r) < 1.0:
        z, se = np.arctanh(r), 1.0 / np.sqrt(n - 3)
        out["pearson_ci_lo"], out["pearson_ci_hi"] = float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))
    return out


def metric_columns(df):
    return [m for m, _ in METRICS if m in df.columns]


def build_scopes(df, flags):
    """List of (scope, anchor_label, frame). The shared draws are never pooled across anchors."""
    present = flags["anchors"]
    scopes = []
    for a in present:
        sub = df[df["anchor_label"] == a]
        scopes.append(("anchor_all", a, sub))
        scopes.append(("anchor_shared", a, sub[sub["batch"] == "shared"]))
        scopes.append(("anchor_indep", a, sub[sub["batch"] == "indep"]))
    if flags["pool_ok"]:
        ind = df[df["batch"] == "indep"].copy()
        scopes.append(("pooled_indep", "ALL_INDEP", ind))
        cols = PARAMS + metric_columns(df)
        cen = ind.copy()
        for c in cols:
            cen[c] = pd.to_numeric(cen[c], errors="coerce")
        cen[cols] = cen[cols] - cen.groupby("anchor_label")[cols].transform("mean")
        scopes.append(("pooled_indep_centered", "ALL_INDEP", cen))
    return scopes


def compute_correlations(df, flags):
    rows = []
    mcols = metric_columns(df)
    for scope, anchor, sub in build_scopes(df, flags):
        for p in PARAMS:
            for m in mcols:
                row = dict(scope=scope, anchor_label=anchor, parameter=p, parameter_short=PARAM_SHORT[p],
                           metric=m, metric_label=METRIC_LABEL[m])
                row.update(corr_pair(sub[p], sub[m]))
                rows.append(row)
    return pd.DataFrame(rows)


def fisher_two_sample(r1, n1, r2, n2):
    """Two-sided p that two correlations from independent samples differ (Fisher z)."""
    if not (np.isfinite(r1) and np.isfinite(r2)) or n1 < 4 or n2 < 4 or abs(r1) >= 1 or abs(r2) >= 1:
        return np.nan, np.nan
    z = (np.arctanh(r1) - np.arctanh(r2)) / np.sqrt(1.0 / (n1 - 3) + 1.0 / (n2 - 3))
    return float(z), float(2 * stats.norm.sf(abs(z)))


def batch_consistency(df, flags):
    rows = []
    for a in flags["anchors"]:
        sub = df[df["anchor_label"] == a]
        sh, ind = sub[sub["batch"] == "shared"], sub[sub["batch"] == "indep"]
        for p in PARAMS:
            c1, c2 = corr_pair(sh[p], sh["kge_2012"]), corr_pair(ind[p], ind["kge_2012"])
            z, pv = fisher_two_sample(c1["pearson_r"], c1["n"], c2["pearson_r"], c2["n"])
            rows.append(dict(anchor_label=a, parameter=p, parameter_short=PARAM_SHORT[p],
                             n_shared=c1["n"], r_shared=c1["pearson_r"], n_indep=c2["n"], r_indep=c2["pearson_r"],
                             difference=(c1["pearson_r"] - c2["pearson_r"]), fisher_z=z, fisher_p=pv))
    return pd.DataFrame(rows)


# ======================================================================
# Anchor summary, shared-draw comparison, PCA
# ======================================================================
def top_runs(sub, top_frac, min_n):
    n_keep = max(min_n, int(np.ceil(len(sub) * top_frac)))
    return sub.sort_values("kge_2012", ascending=False).head(n_keep)


def anchor_summary(df, flags, top_frac, min_n, good_kge):
    rows = []
    for a in flags["anchors"]:
        sub = df[df["anchor_label"] == a]
        k = sub["kge_2012"]
        top = top_runs(sub, top_frac, min_n)
        row = dict(anchor_label=a,
                   Ks_mult=float(sub["Ks_mult"].iloc[0]) if "Ks_mult" in sub else np.nan,
                   f_RS_abs=float(sub["f_RS_abs"].iloc[0]) if "f_RS_abs" in sub else np.nan,
                   n_runs=len(sub), n_shared=int((sub["batch"] == "shared").sum()), n_indep=int((sub["batch"] == "indep").sum()),
                   kge2012_median=float(k.median()), kge2012_q25=float(k.quantile(0.25)), kge2012_q75=float(k.quantile(0.75)),
                   kge2012_min=float(k.min()), kge2012_max=float(k.max()),
                   good_kge_threshold=good_kge, n_good=int((k >= good_kge).sum()), share_good=float((k >= good_kge).mean()))
        if "pbias_pct" in sub:
            pb = pd.to_numeric(sub["pbias_pct"], errors="coerce")
            row.update(pbias_median=float(pb.median()), pbias_min=float(pb.min()), pbias_max=float(pb.max()),
                       n_within_2pct=int((pb.abs() < 2.0).sum()))
        row.update(n_top=len(top), top_kge2012_min=float(top["kge_2012"].min()))
        for p in PARAMS:
            b = BOX[p]
            s = PARAM_SHORT[p]
            row[f"top_mean_{s}"] = float(top[p].mean())
            row[f"top_sd_{s}"] = float(top[p].std(ddof=1)) if len(top) > 1 else np.nan
            row[f"top_mean_{s}_minus_truth_in_box_widths"] = float((top[p].mean() - b["truth"]) / (b["hi"] - b["lo"]))
            row[f"all_sd_{s}"] = float(sub[p].std(ddof=1))
        for (p, q) in (("kinemvelcoef", "flowexp"), ("kinemvelcoef", "channelroughness"), ("flowexp", "channelroughness")):
            nm = f"{PARAM_SHORT[p]}_{PARAM_SHORT[q]}"
            row[f"all_corr_{nm}"] = corr_pair(sub[p], sub[q])["pearson_r"]
            row[f"top_corr_{nm}"] = corr_pair(top[p], top[q])["pearson_r"]
        rows.append(row)
    return pd.DataFrame(rows)


def shared_analysis(df, flags):
    """Returns (wide table of KGE_2012 by draw and anchor, rank-stability long table, vs-truth-anchor table)."""
    sh = df[df["batch"] == "shared"]
    wide = sh.pivot(index="draw_idx", columns="anchor_label", values="kge_2012")
    wide = wide[[a for a in flags["anchors"] if a in wide.columns]].dropna()
    anchors = list(wide.columns)
    rank_rows = []
    for i, a in enumerate(anchors):
        for b in anchors[i:]:
            rho, p = stats.spearmanr(wide[a], wide[b]) if a != b else (1.0, 0.0)
            rank_rows.append(dict(anchor_a=a, anchor_b=b, n_draws=len(wide), spearman_rho=float(rho), spearman_p=float(p)))
    vs_rows = []
    if TRUE_ANCHOR in anchors:
        top_true = set(wide[TRUE_ANCHOR].nlargest(TOP_OVERLAP_K).index)
        for a in anchors:
            if a == TRUE_ANCHOR:
                continue
            d = wide[a] - wide[TRUE_ANCHOR]
            try:
                pw = float(stats.wilcoxon(d).pvalue)
            except ValueError:                                  # all differences exactly zero
                pw = np.nan
            vs_rows.append(dict(anchor_label=a, n_draws=len(wide), median_diff=float(d.median()),
                                q25_diff=float(d.quantile(0.25)), q75_diff=float(d.quantile(0.75)),
                                share_better_than_truth_anchor=float((d > 0).mean()), wilcoxon_p=pw,
                                **{f"top{TOP_OVERLAP_K}_overlap_with_truth_anchor": len(top_true & set(wide[a].nlargest(TOP_OVERLAP_K).index))}))
    return wide, pd.DataFrame(rank_rows), pd.DataFrame(vs_rows)


def run_pca(df, flags, top_frac, min_n):
    """SVD-based PCA of the standardized (cv, flowexp, n) of each anchor's top runs.
    Sign of each component is fixed so its largest-|loading| parameter is positive (the raw sign is arbitrary)."""
    rows = []
    for a in flags["anchors"]:
        sub = df[df["anchor_label"] == a]
        top = top_runs(sub, top_frac, min_n)
        if len(top) < min_n:
            print(f"  Skipping PCA for {a}: only {len(top)} runs (less than --min_n {min_n}).")
            continue
        X = top[PARAMS].to_numpy(dtype=float)
        sd = X.std(axis=0, ddof=1)
        Xs = (X - X.mean(axis=0)) / np.where(sd == 0, 1.0, sd)
        _, S, Vt = np.linalg.svd(Xs, full_matrices=False)
        ratio = (S ** 2) / (S ** 2).sum()
        for i in range(Vt.shape[0]):
            v = Vt[i].copy()
            if v[np.argmax(np.abs(v))] < 0:
                v = -v
            rows.append(dict(anchor_label=a, n_top_runs=len(top), top_frac=top_frac, component=f"PC{i + 1}",
                             explained_variance_ratio=float(ratio[i]),
                             loading_cv=float(v[0]), loading_flowexp=float(v[1]), loading_n=float(v[2])))
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


def star(p, p_flag):
    return "*" if (np.isfinite(p) and p < p_flag) else " "


def sec(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def print_report(df, flags, corr, batch, asum, wide, rank, vs, pca, args):
    present = flags["anchors"]
    pf = args.p_flag

    sec("1. How the runs score, by anchor (KGE_2012; 1.0 is a perfect match to the truth)")
    print(f"  {'anchor':<12} {'Ks':>5} {'f':>6} {'runs':>5} {'median':>7} {'q25-q75':>13} {'best':>6} "
          f"{'>= ' + format(args.good_kge, 'g'):>7} {'PBIAS median':>13} {'|PBIAS|<2%':>11}")
    for _, r in asum.iterrows():
        pb = (f"{r['pbias_median']:>+13.2f} {int(r['n_within_2pct']):>8d}/{int(r['n_runs'])}"
              if "pbias_median" in asum.columns else "  (no pbias_pct column)")
        print(f"  {r.anchor_label:<12} {r.Ks_mult:>5.2f} {r.f_RS_abs:>6.3f} {r.n_runs:>5d} {r.kge2012_median:>7.3f} "
              f"{r.kge2012_q25:>6.3f}-{r.kge2012_q75:<6.3f} {r.kge2012_max:>6.3f} {int(r.n_good):>7d} " + pb)
    print("  (Ks and f are held fixed at each anchor; only cv, flowexp and n vary between the runs of one anchor.)")

    sec("2. KGE_2012 against each routing parameter, per anchor (all 100 runs at that anchor)")
    print(f"  {'anchor':<12}" + "".join(f"{PARAM_SHORT[p]:>26}" for p in PARAMS))
    print(f"  {'':<12}" + "".join(f"{'Pearson r':>13}{'Spearman':>13}" for _ in PARAMS))
    sub = corr[(corr.metric == "kge_2012") & (corr.scope == "anchor_all")]
    for a in present:
        line = f"  {a:<12}"
        for p in PARAMS:
            r = sub[(sub.anchor_label == a) & (sub.parameter == p)].iloc[0]
            line += f"{fnum(r.pearson_r, 9)}{star(r.pearson_p, pf)}   {fnum(r.spearman_rho, 9)}{star(r.spearman_p, pf)}  "
        print(line)
    if flags["pool_ok"]:
        print("  " + "-" * 76)
        for scope, name in (("pooled_indep", "pooled indep"), ("pooled_indep_centered", "  centered")):
            sub2 = corr[(corr.metric == "kge_2012") & (corr.scope == scope)]
            line = f"  {name:<12}"
            for p in PARAMS:
                r = sub2[sub2.parameter == p].iloc[0]
                line += f"{fnum(r.pearson_r, 9)}{star(r.pearson_p, pf)}   {fnum(r.spearman_rho, 9)}{star(r.spearman_p, pf)}  "
            print(line)
        r_un = corr[(corr.metric == "kge_2012") & (corr.scope == "pooled_indep")].set_index("parameter")
        print("  95% interval for the pooled Pearson r (250 independent runs):")
        for p in PARAMS:
            r = r_un.loc[p]
            print(f"      {PARAM_SHORT[p]:<8} {r.pearson_r:+.3f}   [{r.pearson_ci_lo:+.3f}, {r.pearson_ci_hi:+.3f}]")
        print("  'centered' = each anchor's own mean removed before pooling, so a different typical score at")
        print("  each anchor cannot leak into the pooled correlation. Pooled and centered agree when anchors do not matter.")
    else:
        print("  (Pooled numbers skipped: see the checks above.)")
    print(f"  * = p below {pf:g}, not corrected for the number of tests. With 100 runs a correlation is uncertain by about +/-0.2.")

    sec("3. Do the two batches agree? KGE_2012 correlation in the shared 50 vs the independent 50, per anchor")
    print("  Same anchor, two separate sets of 50 runs. A large difference would mean 50 runs is too few to trust.")
    print(f"  {'anchor':<12}{'param':<9}{'r shared':>10}{'r indep':>10}{'difference':>12}{'Fisher p':>10}")
    for _, r in batch.iterrows():
        print(f"  {r.anchor_label:<12}{r.parameter_short:<9}{fnum(r.r_shared, 10)}{fnum(r.r_indep, 10)}"
              f"{fnum(r.difference, 12)}{fnum(r.fisher_p, 10, 2, sign=False)}")
    pv = batch["fisher_p"].dropna()
    if len(pv):
        print(f"  {int((pv < 0.05).sum())} of {len(pv)} differences have p < 0.05; about {0.05 * len(pv):.1f} would by chance alone.")

    sec("4. All routing-parameter x metric correlations: pooled independent runs (Spearman rho)")
    if flags["pool_ok"]:
        print("  'pooled' = 250 independent runs. 'anchors' = lowest and highest rho among the five anchors (all 100 runs each).")
        for p in PARAMS:
            print(f"\n  {PARAM_LONG[p]}")
            print(f"    {'metric':<32}{'pooled rho':>11}{'p':>10}   {'anchors: lowest':>16}{'highest':>9}")
            for m, label in METRICS:
                if m not in df.columns:
                    continue
                r = corr[(corr.scope == "pooled_indep") & (corr.parameter == p) & (corr.metric == m)].iloc[0]
                per = corr[(corr.scope == "anchor_all") & (corr.parameter == p) & (corr.metric == m)]["spearman_rho"]
                print(f"    {label:<32}{fnum(r.spearman_rho, 11)}{star(r.spearman_p, pf)}{fpv(r.spearman_p)}   "
                      f"{fnum(per.min(), 16)}{fnum(per.max(), 9)}")
        print("\n  Sign reminder: for error metrics (peak error, first arrival error, PBIAS) the sign is the direction of the")
        print("  error, not 'better or worse'. A rho near zero for a ratio metric that is 1.0 at the truth is not the same as no effect.")
    else:
        print("  (skipped: pooling is not fair for this file, see the checks above.)")

    sec("5. The same 50 routing sets at every anchor (the shared draws)")
    if flags["shared_ok"] and len(wide):
        print(f"  {len(wide)} routing sets are in every anchor. Spearman rho of KGE_2012 between anchors (1.0 = the sets rank identically):")
        anchors = list(wide.columns)
        print("    " + " " * 12 + "".join(f"{a:>13}" for a in anchors))
        lookup = {(r.anchor_a, r.anchor_b): r.spearman_rho for r in rank.itertuples()}
        for a in anchors:
            line = f"    {a:<12}"
            for b in anchors:
                v = lookup.get((a, b), lookup.get((b, a), np.nan))
                line += f"{v:>13.3f}"
            print(line)
        print("  Best shared routing set at each anchor (draw_idx, KGE_2012):")
        for a in anchors:
            i = wide[a].idxmax()
            print(f"      {a:<12} draw {int(i):>2}   KGE_2012 {wide[a].max():.3f}")
        if len(vs):
            print(f"\n  The same routing set at another anchor, minus its score at the truth anchor ({TRUE_ANCHOR}):")
            print(f"    {'anchor':<12}{'median diff':>12}{'q25 to q75':>18}{'better':>9}{'Wilcoxon p':>12}{'top-' + str(TOP_OVERLAP_K) + ' overlap':>16}")
            for r in vs.itertuples():
                ov = getattr(r, f"top{TOP_OVERLAP_K}_overlap_with_truth_anchor")
                print(f"    {r.anchor_label:<12}{r.median_diff:>+12.3f}{r.q25_diff:>+9.3f} {r.q75_diff:>+8.3f}"
                      f"{r.share_better_than_truth_anchor:>9.0%}{fpv(r.wilcoxon_p, 12)}{ov:>12d} of {TOP_OVERLAP_K}")
            print("  'better' = the share of the 50 routing sets that score higher at that anchor than at the truth anchor.")
            print(f"  'top-{TOP_OVERLAP_K} overlap' = how many of the truth anchor's {TOP_OVERLAP_K} best routing sets are also among that anchor's {TOP_OVERLAP_K} best.")
    else:
        print("  (skipped: see the checks above.)")

    sec(f"6. Where the best runs sit (top {args.top_frac:.0%} of each anchor by KGE_2012) and the PCA")
    print("  Mean of the top runs, as a distance above (+) or below (-) the truth, in sampling-box widths (box width = upper minus lower limit):")
    print(f"    {'anchor':<12}{'top n':>6}{'KGE_2012 >=':>12}" + "".join(f"{PARAM_SHORT[p]:>10}" for p in PARAMS))
    for _, r in asum.iterrows():
        line = f"    {r.anchor_label:<12}{int(r.n_top):>6d}{r.top_kge2012_min:>12.3f}"
        for p in PARAMS:
            line += f"{r[f'top_mean_{PARAM_SHORT[p]}_minus_truth_in_box_widths']:>+10.2f}"
        print(line)
    mid = {p: ((BOX[p]["lo"] + BOX[p]["hi"]) / 2 - BOX[p]["truth"]) / (BOX[p]["hi"] - BOX[p]["lo"]) for p in PARAMS}
    print(f"    {'box middle':<12}{'':>6}{'':>12}" + "".join(f"{mid[p]:>+10.2f}" for p in PARAMS) + "   <- where a randomly drawn run sits on average")
    print("  A top-run mean close to the 'box middle' row means the best runs are not concentrated on that parameter;")
    print("  a top-run mean well away from it (towards 0, the truth) means the best runs are pulled towards the truth value.")
    bn = BOX["channelroughness"]
    print(f"  Note: the truth n ({bn['truth']}) is only {(bn['truth'] - bn['lo']) / (bn['hi'] - bn['lo']):.3f} of a box width above the lower edge of the sampled n range,")
    print("  so almost every sampled n is above the truth; the sampled range is not centred on it.")
    print("\n  Correlation between routing parameters among the top runs (about 20 runs, so +-0.4 uncertain; compare to the all-runs value):")
    print(f"    {'anchor':<12}{'cv-flowexp':>20}{'cv-n':>20}{'flowexp-n':>20}")
    print(f"    {'':<12}" + "".join(f"{'top':>10}{'all':>10}" for _ in range(3)))
    for _, r in asum.iterrows():
        print(f"    {r.anchor_label:<12}" + "".join(f"{fnum(r[f'top_corr_{nm}'], 10)}{fnum(r[f'all_corr_{nm}'], 10)}"
                                                     for nm in ("cv_flowexp", "cv_n", "flowexp_n")))
    print("  If two parameters trade off along a ridge, their correlation among the top runs is far from zero while the all-runs value is near zero.")
    print("\n  PCA of standardized (cv, flowexp, n) among the top runs (sign set so the biggest loading is positive):")
    if pca.empty:
        print("  (no PCA: too few top runs)")
    else:
        for a in present:
            s = pca[pca.anchor_label == a]
            if s.empty:
                continue
            print(f"    {a} ({int(s.n_top_runs.iloc[0])} runs)")
            for r in s.itertuples():
                print(f"      {r.component}  {r.explained_variance_ratio:>6.1%} of variance   cv {r.loading_cv:+.2f}   flowexp {r.loading_flowexp:+.2f}   n {r.loading_n:+.2f}")
        print("  A component that holds a large share of the variance is a direction in which the top runs are spread out (not pinned down).")
        print("  With about 20 runs per anchor the loadings are rough; treat them as a pointer to the scatter plots, not a result.")


# ======================================================================
# Figures
# ======================================================================
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "DejaVu Sans"
    return plt


def style_axes(ax, grid=True):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    if grid:
        ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_SOFT, labelsize=8.5)


def anchor_titles(df, anchors):
    out = {}
    for a in anchors:
        sub = df[df["anchor_label"] == a]
        ks = float(sub["Ks_mult"].iloc[0]) if "Ks_mult" in sub else np.nan
        f = float(sub["f_RS_abs"].iloc[0]) if "f_RS_abs" in sub else np.nan
        tag = "  (truth)" if a == TRUE_ANCHOR else ""
        out[a] = f"{a}{tag}\nKs {ks:g}, f {f:g}"
    return out


def luminance(rgba):
    r, g, b = rgba[:3]
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def save_fig(fig, path, dpi):
    fig.savefig(path, dpi=dpi, facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)
    print(f"  saved {path.name}")


def figure_1(df, flags, corr, titles, out_path, dpi):
    plt = _plt()
    from matplotlib.lines import Line2D
    anchors = flags["anchors"]
    nrow, ncol = len(PARAMS), len(anchors)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.9 * ncol, 2.75 * nrow + 1.0), sharey=True, squeeze=False)
    fig.patch.set_facecolor("white")
    ymin = max(0.0, np.floor((df["kge_2012"].min() - 0.02) * 10) / 10)
    sub_all = corr[(corr.scope == "anchor_all") & (corr.metric == "kge_2012")]
    for i, p in enumerate(PARAMS):
        b = BOX[p]
        pad = 0.04 * (b["hi"] - b["lo"])
        for j, a in enumerate(anchors):
            ax = axes[i][j]
            style_axes(ax)
            s = df[df["anchor_label"] == a]
            sh, ind = s[s["batch"] == "shared"], s[s["batch"] == "indep"]
            ax.scatter(sh[p], sh["kge_2012"], s=17, facecolors="none", edgecolors=INK_SOFT, linewidths=0.9, zorder=3)
            ax.scatter(ind[p], ind["kge_2012"], s=17, color=BLUE, alpha=0.78, linewidths=0, zorder=3.2)
            ax.axvline(b["truth"], color=INK, linestyle=(0, (4, 3)), linewidth=1.0, zorder=2)
            ax.set_xlim(b["lo"] - pad, b["hi"] + pad)
            ax.set_ylim(ymin, 1.01)
            r = sub_all[(sub_all.anchor_label == a) & (sub_all.parameter == p)].iloc[0]
            ax.text(0.04, 0.05, f"rho = {r.spearman_rho:+.2f}", transform=ax.transAxes, fontsize=9, color=INK,
                    ha="left", va="bottom", bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec="none", alpha=0.85))
            ax.set_xlabel(PARAM_LONG[p], fontsize=8.5, color=INK_SOFT)
            if j == 0:
                ax.set_ylabel("KGE_2012", fontsize=9.5, color=INK_SOFT)
            if i == 0:
                ax.set_title(titles[a], fontsize=9.5, color=INK, pad=6)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor="none", markeredgecolor=INK_SOFT, markersize=6,
               label="shared draws (the same 50 at every anchor)"),
        Line2D([0], [0], marker="o", linestyle="", color=BLUE, alpha=0.8, markersize=6, label="independent draws"),
        Line2D([0], [0], color=INK, linestyle=(0, (4, 3)), linewidth=1.0, label="truth value of the parameter"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle("Series 121: KGE_2012 against each routing parameter, one column per anchor", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.005, "rho = Spearman rank correlation over all 100 runs at that anchor. One storm (1.0x), channel loss off.",
             ha="center", va="bottom", fontsize=8.5, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.025, 1, 0.925))
    save_fig(fig, out_path, dpi)


def heat_panel(ax, mat, rows, cols, norm, cmap, bold, p_note=None):
    ax.imshow(mat, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=32, ha="right", fontsize=8.5, color=INK_SOFT)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=8.5, color=INK_SOFT)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if not np.isfinite(v):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=8, color=INK_SOFT)
                continue
            dark = luminance(cmap(norm(v))) < 0.42
            ax.text(j, i, f"{v:+.2f}" if norm.vmin < 0 else f"{v:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if dark else INK, fontweight="bold" if bold[i, j] else "normal")


def figure_2(df, flags, corr, titles, out_path, dpi, p_flag):
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    anchors = flags["anchors"]
    mcols = [m for m, _ in METRICS if m in df.columns]
    col_labels = [METRIC_LABEL[m] for m in mcols]
    row_keys = [("anchor_all", a) for a in anchors] + ([("pooled_indep", "ALL_INDEP")] if flags["pool_ok"] else [])
    counts = {a: int((df["anchor_label"] == a).sum()) for a in anchors}
    row_labels = [f"{a}  (n={counts[a]})" for a in anchors] + (
        [f"pooled independent  (n={int(corr[(corr.scope == 'pooled_indep')]['n'].max())})"] if flags["pool_ok"] else [])
    cmap = LinearSegmentedColormap.from_list("div", [BLUE, NEUTRAL_MID, RED])
    norm = Normalize(vmin=-1, vmax=1)
    fig, axes = plt.subplots(len(PARAMS), 1, figsize=(11.8, 3.15 * len(PARAMS) + 1.2), squeeze=False)
    fig.patch.set_facecolor("white")
    for k, p in enumerate(PARAMS):
        ax = axes[k][0]
        mat = np.full((len(row_keys), len(mcols)), np.nan)
        bold = np.zeros_like(mat, dtype=bool)
        for i, (scope, a) in enumerate(row_keys):
            for j, m in enumerate(mcols):
                r = corr[(corr.scope == scope) & (corr.anchor_label == a) & (corr.parameter == p) & (corr.metric == m)]
                if len(r):
                    mat[i, j] = r.iloc[0]["spearman_rho"]
                    bold[i, j] = bool(np.isfinite(r.iloc[0]["spearman_p"]) and r.iloc[0]["spearman_p"] < p_flag)
        heat_panel(ax, mat, row_labels, col_labels, norm, cmap, bold)
        if flags["pool_ok"]:
            ax.axhline(len(anchors) - 0.5, color=INK, linewidth=1.6)
        ax.set_title(PARAM_LONG[p], loc="left", fontsize=11, color=INK, pad=6)
    fig.subplots_adjust(left=0.2, right=0.9, top=0.93, bottom=0.12, hspace=0.95)
    cax = fig.add_axes([0.925, 0.25, 0.015, 0.5])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Spearman rho", fontsize=9, color=INK_SOFT)
    cb.ax.tick_params(labelsize=8, colors=INK_SOFT)
    cb.outline.set_visible(False)
    fig.suptitle("Series 121: Spearman rho between each routing parameter and each metric", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.5, 0.004, f"Bold = p below {p_flag:g} (not corrected for the number of tests). Rows above the line: all 100 runs at one anchor. "
                         "Bottom row: the 250 independent runs pooled.", ha="center", va="bottom", fontsize=8.3, color=INK_SOFT)
    save_fig(fig, out_path, dpi)


def figure_3(df, flags, wide, rank, titles, out_path, dpi):
    plt = _plt()
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    anchors = list(wide.columns)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.6), gridspec_kw=dict(width_ratios=[1.25, 1]))
    fig.patch.set_facecolor("white")
    style_axes(ax1)
    x = np.arange(len(anchors))
    for _, row in wide.iterrows():
        ax1.plot(x, row.to_numpy(), color=GREY_DOT, linewidth=0.8, alpha=0.7, zorder=2)
    med = wide.median().to_numpy()
    ax1.plot(x, med, color=INK, linewidth=2.4, marker="o", markersize=6, zorder=4)
    ax1.text(x[-1] + 0.08, med[-1], "median", fontsize=9, color=INK, va="center")
    ax1.set_xticks(x)
    ax1.set_xticklabels([titles[a] for a in anchors], fontsize=8.5)
    ax1.set_xlim(-0.3, len(anchors) - 1 + 0.75)
    ax1.set_ylabel("KGE_2012", fontsize=9.5, color=INK_SOFT)
    ax1.set_title(f"The same {len(wide)} routing sets run at each anchor (one grey line each)", fontsize=10.5, color=INK, loc="left")
    # rank-agreement heatmap
    style_axes(ax2, grid=False)
    mat = np.full((len(anchors), len(anchors)), np.nan)
    for r in rank.itertuples():
        i, j = anchors.index(r.anchor_a), anchors.index(r.anchor_b)
        mat[i, j] = mat[j, i] = r.spearman_rho
    lo = min(0.0, float(np.nanmin(mat)))
    cmap = LinearSegmentedColormap.from_list("seq", ["#f3f7fc", BLUE])
    norm = Normalize(vmin=lo, vmax=1.0)
    ax2.set_facecolor("white")
    heat_panel(ax2, mat, anchors, anchors, norm, cmap, np.zeros_like(mat, dtype=bool))
    ax2.set_title("Do the routing sets rank the same way at every anchor?", fontsize=10.5, color=INK, loc="left")
    ax2.set_xlabel("Spearman rho of KGE_2012 between anchors (colour scale 0 to 1)", fontsize=8.5, color=INK_SOFT)
    fig.suptitle("Series 121: the shared draws", fontsize=12.5, color=INK, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save_fig(fig, out_path, dpi)


def figure_4(df, flags, asum, pca, titles, out_path, dpi, top_frac, min_n):
    plt = _plt()
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    anchors = flags["anchors"]
    pairs = [("kinemvelcoef", "flowexp"), ("kinemvelcoef", "channelroughness"), ("flowexp", "channelroughness")]
    ncol = len(anchors)
    fig, axes = plt.subplots(4, ncol, figsize=(2.9 * ncol, 12.3), squeeze=False)
    fig.patch.set_facecolor("white")
    for j, a in enumerate(anchors):
        s = df[df["anchor_label"] == a]
        top = top_runs(s, top_frac, min_n)
        rest = s[~s["run_id"].isin(top["run_id"])]
        for i, (px, py) in enumerate(pairs):
            ax = axes[i][j]
            style_axes(ax)
            ax.scatter(rest[px], rest[py], s=14, color=GREY_DOT, linewidths=0, alpha=0.9, zorder=3)
            ax.scatter(top[px], top[py], s=24, color=INK, linewidths=0, zorder=3.4)
            ax.scatter([BOX[px]["truth"]], [BOX[py]["truth"]], s=150, marker="*", color=GOLD, edgecolors=INK, linewidths=0.9, zorder=5)
            for axis_name, key in (("x", px), ("y", py)):
                b = BOX[key]
                pad = 0.04 * (b["hi"] - b["lo"])
                (ax.set_xlim if axis_name == "x" else ax.set_ylim)(b["lo"] - pad, b["hi"] + pad)
            ax.set_xlabel(PARAM_SHORT[px], fontsize=9, color=INK_SOFT)
            ax.set_ylabel(PARAM_SHORT[py] if j == 0 else "", fontsize=9, color=INK_SOFT)
            if i == 0:
                ax.set_title(titles[a], fontsize=9.5, color=INK, pad=6)
        # PCA loadings row
        ax = axes[3][j]
        style_axes(ax)
        sp = pca[pca.anchor_label == a] if len(pca) else pca
        if len(sp):
            xs = np.arange(len(sp))
            w = 0.26
            for k, (col, p) in enumerate((("loading_cv", "kinemvelcoef"), ("loading_flowexp", "flowexp"), ("loading_n", "channelroughness"))):
                ax.bar(xs + (k - 1) * w, sp[col], w, color=PARAM_COLOR[p], zorder=3)
            ax.set_xticks(xs)
            ax.set_xticklabels([f"{r.component}\n{r.explained_variance_ratio:.0%}" for r in sp.itertuples()], fontsize=8.5)
            ax.axhline(0, color=INK, linewidth=0.8)
            ax.set_ylim(-1.05, 1.05)
        else:
            ax.text(0.5, 0.5, "no PCA\n(too few runs)", transform=ax.transAxes, ha="center", va="center", color=INK_SOFT)
        if j == 0:
            ax.set_ylabel("PCA loading", fontsize=9, color=INK_SOFT)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", color=INK, markersize=6, label=f"top {top_frac:.0%} of the anchor's runs (by KGE_2012)"),
        Line2D([0], [0], marker="o", linestyle="", color=GREY_DOT, markersize=6, label="the other runs"),
        Line2D([0], [0], marker="*", linestyle="", markerfacecolor=GOLD, markeredgecolor=INK, markersize=12, label="truth"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.962))
    bar_handles = [Patch(color=PARAM_COLOR[p], label=PARAM_SHORT[p]) for p in PARAMS]
    fig.legend(handles=bar_handles, loc="lower center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.012),
               title="PCA loadings (percent = share of variance in that component)", title_fontsize=8.5)
    fig.suptitle("Series 121: where the best runs sit in (cv, flowexp, n), and a PCA of those runs", fontsize=12.5, color=INK, y=0.997)
    fig.tight_layout(rect=(0, 0.045, 1, 0.945))
    save_fig(fig, out_path, dpi)


# ======================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--csv", type=Path, default=None, help=f"the results file (default .../summary_tables/{DEFAULT_CSV_NAME})")
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--top_frac", type=float, default=0.20)
    ap.add_argument("--min_n", type=int, default=8)
    ap.add_argument("--good_kge", type=float, default=0.95)
    ap.add_argument("--p_flag", type=float, default=0.01)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()
    if not (0.05 <= args.top_frac <= 0.5):
        ap.error("--top_frac must be between 0.05 and 0.5.")
    if args.min_n < 5:
        ap.error("--min_n must be at least 5 (a PCA on fewer runs means nothing).")
    if not (0.0 < args.good_kge < 1.0):
        ap.error("--good_kge must be between 0 and 1.")
    if not (0.0 < args.p_flag < 1.0):
        ap.error("--p_flag must be between 0 and 1.")
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

    print("\nComputing...")
    corr = compute_correlations(df, flags)
    batch = batch_consistency(df, flags)
    asum = anchor_summary(df, flags, args.top_frac, args.min_n, args.good_kge)
    if flags["shared_ok"]:
        wide, rank, vs = shared_analysis(df, flags)
    else:
        wide, rank, vs = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    pca = run_pca(df, flags, args.top_frac, args.min_n)

    print_report(df, flags, corr, batch, asum, wide, rank, vs, pca, args)

    out_dir.mkdir(parents=True, exist_ok=True)
    sec("Files written to " + str(out_dir))
    tables = [("corr", corr), ("batch", batch), ("anchors", asum), ("rank", rank), ("vs_truth", vs), ("pca", pca)]
    for key, table in tables:
        if table.empty:
            print(f"  (skipped {FILES[key]}: nothing to write)")
            continue
        table.to_csv(out_dir / FILES[key], index=False, float_format="%.6g")
        print(f"  saved {FILES[key]}  ({len(table)} rows)")

    if not args.no_plots:
        print("\nDrawing the figures...")
        titles = anchor_titles(df, flags["anchors"])
        figure_1(df, flags, corr, titles, out_dir / FILES["fig1"], args.dpi)
        figure_2(df, flags, corr, titles, out_dir / FILES["fig2"], args.dpi, args.p_flag)
        if flags["shared_ok"] and len(wide) >= 2 and wide.shape[1] >= 2:
            figure_3(df, flags, wide, rank, titles, out_dir / FILES["fig3"], args.dpi)
        else:
            print(f"  (skipped {FILES['fig3']}: no usable shared draws)")
        figure_4(df, flags, asum, pca, titles, out_dir / FILES["fig4"], args.dpi, args.top_frac, args.min_n)

    print("\nDone. This script only reads the results file; it changes nothing in it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
