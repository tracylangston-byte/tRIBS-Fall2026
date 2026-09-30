"""
analyze_cc_stage2_110.py
==========================
Series 110 -- Stage 2 analysis. Runs AFTER run_cc_stage2_joint_lhs_110.py
has produced lhs_results_joint_Ks_f_cc_110.csv. This script does not run
tRIBS or build anything -- it is pure post-hoc analysis, in the same spirit
as analyze_series101_kge_components.py sitting downstream of the Series 101
sweep script.

Answers the actual question Stage 2 exists to answer: does letting
channelconductivity_mmhr vary jointly with Ks_mult/f_RS_abs pull the
recovered "best" Ks/f away from the current truth (7.0x / 0.012), and by
how much? Four analyses, each producing one CSV (+ printed summary):

1. MARGINAL CORRELATIONS (Pearson r, Ks_mult/f_RS_abs/cc vs PBIAS/KGE_2012)
   -> lhs_results_joint_Ks_f_cc_correlations_110.csv
   Baseline context, same format as series101_correlation_summary.csv.

2. CC-BAND PEAK TRACKING (the sharp test)
   -> lhs_results_joint_Ks_f_cc_band_peaks_110.csv
   Splits the sweep into cc tertiles (low/mid/high) and finds the
   best-KGE_2012 (Ks_mult, f_RS_abs) WITHIN each band -- same move
   ridge_width_vs_f.csv already made for Ks vs. f slices in Series 100.
   If the best-fit Ks/f stays near (7.0, 0.012) across all three cc bands,
   Ks/f are robust to cc being present. If it visibly drifts as cc
   increases, that IS the compensation the handoff was worried about,
   quantified rather than suspected.

3. PCA ON TOP PERFORMERS (top 20% by KGE_2012, standardized Ks/f/cc)
   -> lhs_results_joint_Ks_f_cc_pca_110.csv
   Same method as series101_pca_summary.csv (which did this for cv/r/n).
   If cc loads heavily onto the same principal component as Ks or f,
   that names the equifinal trade-off axis explicitly.

4. FEASIBLE-VOLUME FRACTION
   Fraction of sampled points with |PBIAS| < threshold, compared to
   Series 100's 2D feasible-area fraction (storm_series_summary.csv:
   area_frac_PBIAS_feasible ~ 0.003 for the matched-bounds Ks/f-only
   sweep). If the 3D (Ks/f/cc) feasible fraction is proportionally
   larger, that's a direct numeric statement that equifinality got
   worse with cc added as a free dimension -- not just an impression
   from a scatter plot.

Also flags: any point in the CAUTION ZONE near the untrusted cross-version
outlier corner (Ks~9.83x/f~0.0485) that shows up among top performers,
per Handoff_ChannelLossCalibration_v1.md section 4.

Usage (run from the lab/ directory, after Stage 2 has produced results):
    python analyze_cc_stage2_110.py
    python analyze_cc_stage2_110.py --pbias_threshold 5.0
    python analyze_cc_stage2_110.py --top_frac 0.15
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

TRUTH_KS = 7.0
TRUTH_F  = 0.012

# Series 100's own 2D (Ks/f-only) feasible-area fraction, for comparison
# -- from storm_series_summary_080_100n_125.csv, "100_narrow (1.0x rain,
# matched bounds)" row (the LHS bounds closest in spirit to this sweep's
# Ks/f window; PBIAS-feasibility defined the same way there).
SERIES100_2D_FEASIBLE_FRAC = 0.036982363452723176

CAUTION_KS_LO, CAUTION_KS_HI = 8.85, 10.5
CAUTION_F_LO,  CAUTION_F_HI  = 0.0437, 0.030


def load_results(path):
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run run_cc_stage2_joint_lhs_110.py first "
            f"-- this script only analyzes its output, it doesn't generate "
            f"data itself."
        )
    df = pd.read_csv(path)
    required = {"Ks_mult", "f_RS_abs", "channelconductivity_mmhr",
                "pbias_pct", "kge_2012"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Results CSV is missing expected columns: {missing}")
    print(f"Loaded {len(df)} rows from {path.name}")
    return df


# ------------------------------------------------------------------
# 1. MARGINAL CORRELATIONS
# ------------------------------------------------------------------
def marginal_correlations(df):
    rows = []
    params  = ["Ks_mult", "f_RS_abs", "channelconductivity_mmhr"]
    metrics = ["pbias_pct", "kge_2012", "kge", "nse"]
    for p in params:
        for m in metrics:
            if p not in df.columns or m not in df.columns:
                continue
            sub = df[[p, m]].dropna()
            if len(sub) < 3:
                continue
            r, pval = stats.pearsonr(sub[p], sub[m])
            rho, spval = stats.spearmanr(sub[p], sub[m])
            rows.append({
                "parameter": p, "metric": m,
                "pearson_r": r, "pearson_p": pval,
                "spearman_rho": rho, "spearman_p": spval,
                "n": len(sub),
            })
    out = pd.DataFrame(rows)
    print("\n" + "=" * 70)
    print("1. MARGINAL CORRELATIONS")
    print("=" * 70)
    for p in params:
        sub = out[(out["parameter"] == p) & (out["metric"] == "pbias_pct")]
        if not sub.empty:
            r = sub.iloc[0]["pearson_r"]
            print(f"  {p} vs PBIAS: r={r:+.3f}  "
                  f"(compare: Series 100 Ks_mult vs PBIAS was r~-1.00, "
                  f"Stage 1 cc-only vs PBIAS was r=-0.994)")
    return out


# ------------------------------------------------------------------
# 2. CC-BAND PEAK TRACKING -- the sharp test
# ------------------------------------------------------------------
def cc_band_peak_tracking(df):
    # Tertile split on cc (log-spaced sweep, so quantile-based tertiles on
    # the raw values naturally respect that spacing).
    df = df.copy()
    df["cc_band"] = pd.qcut(df["channelconductivity_mmhr"], 3,
                             labels=["low", "mid", "high"])

    rows = []
    for band in ["low", "mid", "high"]:
        sub = df[df["cc_band"] == band]
        if sub.empty:
            continue
        best = sub.loc[sub["kge_2012"].idxmax()]
        cc_range = f"{sub['channelconductivity_mmhr'].min():.1f}-{sub['channelconductivity_mmhr'].max():.1f}"
        rows.append({
            "cc_band": band,
            "cc_range_mmhr": cc_range,
            "n_in_band": len(sub),
            "best_Ks_mult": best["Ks_mult"],
            "best_f_RS_abs": best["f_RS_abs"],
            "best_cc": best["channelconductivity_mmhr"],
            "best_kge_2012": best["kge_2012"],
            "best_pbias_pct": best["pbias_pct"],
            "Ks_drift_from_truth": best["Ks_mult"] - TRUTH_KS,
            "f_drift_from_truth": best["f_RS_abs"] - TRUTH_F,
            "in_caution_zone": bool(
                CAUTION_KS_LO <= best["Ks_mult"] <= CAUTION_KS_HI and
                CAUTION_F_LO  <= best["f_RS_abs"] <= CAUTION_F_HI
            ),
        })
    out = pd.DataFrame(rows)

    print("\n" + "=" * 70)
    print("2. CC-BAND PEAK TRACKING (the sharp equifinality test)")
    print("=" * 70)
    print(f"  Truth: Ks_mult={TRUTH_KS}, f_RS_abs={TRUTH_F}")
    for _, r in out.iterrows():
        caution = "  [CAUTION ZONE]" if r["in_caution_zone"] else ""
        print(f"  cc {r['cc_band']:>4} ({r['cc_range_mmhr']} mm/hr, n={r['n_in_band']}): "
              f"best Ks={r['best_Ks_mult']:.3f} (drift {r['Ks_drift_from_truth']:+.3f})  "
              f"f={r['best_f_RS_abs']:.4f} (drift {r['f_drift_from_truth']:+.4f})  "
              f"KGE_2012={r['best_kge_2012']:.3f}{caution}")

    if len(out) == 3:
        ks_spread = out["best_Ks_mult"].max() - out["best_Ks_mult"].min()
        f_spread  = out["best_f_RS_abs"].max() - out["best_f_RS_abs"].min()
        print(f"\n  Ks_mult spread across cc bands: {ks_spread:.3f}  "
              f"({'SUBSTANTIAL -- Ks is compensating for cc' if ks_spread > 1.0 else 'modest'})")
        print(f"  f_RS_abs spread across cc bands: {f_spread:.4f}  "
              f"({'SUBSTANTIAL -- f is compensating for cc' if f_spread > 0.006 else 'modest'})")

    return out


# ------------------------------------------------------------------
# 3. PCA ON TOP PERFORMERS
# ------------------------------------------------------------------
def pca_top_performers(df, top_frac):
    n_top = max(int(len(df) * top_frac), 10)
    top = df.nlargest(n_top, "kge_2012")

    X = top[["Ks_mult", "f_RS_abs", "channelconductivity_mmhr"]].copy()
    # log f and cc before standardizing -- both were sampled log-uniform,
    # so PCA on raw values would be dominated by scale, not structure.
    X["f_RS_abs"] = np.log10(X["f_RS_abs"])
    X["channelconductivity_mmhr"] = np.log10(X["channelconductivity_mmhr"])

    scaler = StandardScaler()
    X_std  = scaler.fit_transform(X)

    pca = PCA(n_components=3)
    pca.fit(X_std)

    rows = []
    for i, comp in enumerate(["PC1", "PC2", "PC3"]):
        rows.append({
            "n_runs_used": n_top,
            "top_frac_requested": top_frac,
            "component": comp,
            "explained_variance_ratio": pca.explained_variance_ratio_[i],
            "loading_Ks_mult": pca.components_[i][0],
            "loading_log_f_RS_abs": pca.components_[i][1],
            "loading_log_cc": pca.components_[i][2],
        })
    out = pd.DataFrame(rows)

    print("\n" + "=" * 70)
    print(f"3. PCA ON TOP {top_frac:.0%} PERFORMERS (n={n_top}, standardized Ks/log(f)/log(cc))")
    print("=" * 70)
    for _, r in out.iterrows():
        print(f"  {r['component']} (explains {r['explained_variance_ratio']:.1%}): "
              f"Ks={r['loading_Ks_mult']:+.3f}  log(f)={r['loading_log_f_RS_abs']:+.3f}  "
              f"log(cc)={r['loading_log_cc']:+.3f}")
    pc1 = out.iloc[0]
    dominant = max(
        [("Ks_mult", abs(pc1["loading_Ks_mult"])),
         ("f_RS_abs", abs(pc1["loading_log_f_RS_abs"])),
         ("cc", abs(pc1["loading_log_cc"]))],
        key=lambda t: t[1],
    )
    print(f"\n  PC1 (dominant trade-off axis, {pc1['explained_variance_ratio']:.1%} of variance) "
          f"is led by: {dominant[0]} (|loading|={dominant[1]:.3f})")
    if abs(pc1["loading_Ks_mult"]) > 0.4 and abs(pc1["loading_log_cc"]) > 0.4:
        print("  Ks_mult and cc both load heavily on PC1 -- consistent with "
          "the Ks<->cc compensation the handoff predicted.")
    if abs(pc1["loading_log_f_RS_abs"]) > 0.4 and abs(pc1["loading_log_cc"]) > 0.4:
        print("  f_RS_abs and cc both load heavily on PC1 -- consistent with "
          "an f<->cc compensation axis.")

    return out


# ------------------------------------------------------------------
# 4. FEASIBLE-VOLUME FRACTION
# ------------------------------------------------------------------
def feasible_volume_fraction(df, pbias_threshold):
    feasible = df[df["pbias_pct"].abs() < pbias_threshold]
    frac = len(feasible) / len(df)

    print("\n" + "=" * 70)
    print(f"4. FEASIBLE-VOLUME FRACTION (|PBIAS| < {pbias_threshold}%)")
    print("=" * 70)
    print(f"  This sweep (3D: Ks_mult x f_RS_abs x cc): "
          f"{len(feasible)}/{len(df)} = {frac:.4f} ({frac:.1%})")
    print(f"  Series 100 (2D: Ks_mult x f_RS_abs only, cc off): "
          f"{SERIES100_2D_FEASIBLE_FRAC:.4f} ({SERIES100_2D_FEASIBLE_FRAC:.1%})")

    if frac > 0:
        ratio = frac / SERIES100_2D_FEASIBLE_FRAC
        print(f"  Ratio: {ratio:.1f}x")
        if ratio > 1.5:
            print("  The feasible region got proportionally LARGER with cc "
                  "added as a free dimension -- direct evidence equifinality "
                  "got worse, not just a qualitative impression.")
        elif ratio < 0.67:
            print("  The feasible region got proportionally SMALLER -- cc "
                  "may be constraining the fit more than it's loosening it, "
                  "at least in this parameter window.")
        else:
            print("  Feasible fraction is roughly comparable to the 2D case "
                  "-- no strong evidence either way from this metric alone.")
    else:
        print("  No points met the feasibility threshold in this sweep -- "
              "consider a looser threshold or check sweep coverage.")

    return pd.DataFrame([{
        "sweep": "stage2_joint_Ks_f_cc_110", "pbias_threshold_pct": pbias_threshold,
        "n_total": len(df), "n_feasible": len(feasible), "area_frac_feasible": frac,
    }, {
        "sweep": "series100_2D_Ks_f_only", "pbias_threshold_pct": pbias_threshold,
        "n_total": None, "n_feasible": None,
        "area_frac_feasible": SERIES100_2D_FEASIBLE_FRAC,
    }])


def main():
    parser = argparse.ArgumentParser(
        description="Series 110 Stage 2 -- equifinality analysis against "
                    "the joint Ks/f/cc LHS results.")
    parser.add_argument("--pbias_threshold", type=float, default=2.0,
                        help="|PBIAS| threshold (%%) for feasibility "
                             "(default: 2.0, matching the storm-comparison "
                             "feasibility definition used elsewhere in this "
                             "project)")
    parser.add_argument("--top_frac", type=float, default=0.2,
                        help="Fraction of top KGE_2012 runs to use for PCA "
                             "(default: 0.2, matching series101_pca_summary.csv)")
    args = parser.parse_args()

    script_dir   = Path.cwd()
    project_root = script_dir.parent
    calib_dir    = project_root / "calibration_work"
    summary_dir  = calib_dir / "03_comparisons" / "summary_tables"

    in_path = summary_dir / "lhs_results_joint_Ks_f_cc_110.csv"
    df = load_results(in_path)

    caution_in_data = df.get("caution_zone")
    if caution_in_data is not None and caution_in_data.any():
        n_caution = int(caution_in_data.sum())
        top10 = df.nlargest(min(10, len(df)), "kge_2012")
        n_caution_top10 = int(top10.get("caution_zone", pd.Series(dtype=bool)).sum())
        print(f"\nNOTE: {n_caution} of {len(df)} sampled points fall in the "
              f"untrusted cross-version-outlier caution zone "
              f"(Ks {CAUTION_KS_LO}-{CAUTION_KS_HI}x, f {CAUTION_F_LO}-{CAUTION_F_HI}). "
              f"{n_caution_top10} of the top 10 KGE_2012 runs are in that zone -- "
              f"{'treat top-performer conclusions with extra skepticism until Josh confirms the anomaly is understood.' if n_caution_top10 > 0 else 'none of the current top performers are affected.'}")

    corr_out   = marginal_correlations(df)
    bands_out  = cc_band_peak_tracking(df)
    pca_out    = pca_top_performers(df, args.top_frac)
    feas_out   = feasible_volume_fraction(df, args.pbias_threshold)

    corr_out.to_csv(summary_dir / "lhs_results_joint_Ks_f_cc_correlations_110.csv", index=False)
    bands_out.to_csv(summary_dir / "lhs_results_joint_Ks_f_cc_band_peaks_110.csv", index=False)
    pca_out.to_csv(summary_dir / "lhs_results_joint_Ks_f_cc_pca_110.csv", index=False)
    feas_out.to_csv(summary_dir / "lhs_results_joint_Ks_f_cc_feasible_frac_110.csv", index=False)

    print("\n" + "=" * 70)
    print("Saved 4 output CSVs to", summary_dir)
    print("  lhs_results_joint_Ks_f_cc_correlations_110.csv")
    print("  lhs_results_joint_Ks_f_cc_band_peaks_110.csv")
    print("  lhs_results_joint_Ks_f_cc_pca_110.csv")
    print("  lhs_results_joint_Ks_f_cc_feasible_frac_110.csv")
    print("=" * 70)


if __name__ == "__main__":
    main()
