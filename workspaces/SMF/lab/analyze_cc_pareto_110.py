#!/usr/bin/env python
"""
analyze_cc_pareto_110.py -- PRELIMINARY component Pareto for the Stage 2 channel-loss (cc) study.

What it does
------------
Builds a two-objective Pareto analysis from the DECOMPOSED KGE_2012 components, with KGE_2012
itself deliberately left out (it is built from these components, so pairing it with them is
circular):

    vol = |kge_beta  - 1|    volume (bias) error       -> minimize
    var = |kge_gamma - 1|    variability (CV) error    -> minimize
    cor = 1 - kge_r          correlation error         -> reported only (tie-breaker, NOT an axis;
                                                          it is strongly coupled to vol, rho ~ 0.7-0.8)

It is run on two series that share the SAME 250 (Ks, f, cc) points and the same truth:
    cc ON   Stage 2   (optpercolation = 1)
    cc OFF  CONTROL   (optpercolation = 0)   <- the null: cc is inert here, so any structure in cc
                                                among its front members is chance.

Outputs (all written to --out_dir, default <summary_tables>/pareto_cc_110/; inputs are read-only)
    cc_pareto_110_runs.csv          per-run objectives, layer, near-front flags, box-edge flags
    cc_pareto_110_summary.csv       per series x epsilon: how cc, Ks, f are distributed among near-front runs
    cc_pareto_110_partial_corr.csv  shape metrics vs log10(cc), controlling for Ks and log10(f)
    cc_pareto_110_overview.png      objective space, Ks-f map of near-front runs, cc ECDFs

Definitions
-----------
  layer k        non-dominated sorting layer on (vol, var); layer 1 is the strict Pareto front.
  near-front(e)  run b is near-front at tolerance e if some layer-1 run a satisfies
                 b_vol <= a_vol + e AND b_var <= a_var + e, i.e. b is "indistinguishable from the
                 front at objective tolerance e". Nested in e; e = 0 reduces to the front itself.
  edge flag      run lies within --edge_frac of a sampled-box bound (Ks linear; f and cc log).
                 Front members at an edge mean the true ridge probably continues outside the box.

Caveats (printed again at the end of every run)
  * PRELIMINARY: one storm (Aug 12, 2014), n = 250 points, routing pinned at truth (cv, r, n).
  * The truth has channel loss OFF, so cc-ON runs are structurally mismatched to it by construction:
    the cc-ON front shows how well Ks and f can OFFSET an unwanted cc, not how identifiable cc is.
  * The sampled box clips the ridge (several best cc-ON runs sit at box edges).
  * Finite-sample fronts are noisy; prefer the near-front sets and layers to strict front membership.

Safety
  * Reads only the two CSVs below. Never globs per-run compare files (650 '_s2' files exist, only
    250 belong to Stage 2); run membership is defined solely by run_id rows in the CSVs.
  * Asserts that no output path equals an input path.

Usage (from anywhere under the SMF repo; paths are auto-discovered, or pass them explicitly):
    python analyze_cc_pareto_110.py
    python analyze_cc_pareto_110.py --stage2_csv PATH --control_csv PATH --out_dir PATH
    python analyze_cc_pareto_110.py --eps 0 0.005 0.01 0.02 --eps_plot 0.01
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy import stats
    HAVE_SCIPY = True
except Exception:  # pragma: no cover
    HAVE_SCIPY = False

STAGE2_NAME = "lhs_results_joint_Ks_f_cc_110.csv"
CONTROL_NAME = "lhs_results_joint_Ks_f_cc_CONTROL_110.csv"
REL_TABLES = Path("calibration_work/03_comparisons/summary_tables")

# Sampled box, as recorded in Handoff v3 (Ks linear; f and cc log). The Stage 2 script's LHS_PARAMS
# is the authority: the script warns if any sampled value falls outside these, which would mean
# the numbers here are stale. Override with --ks_range / --f_range / --cc_range.
DEFAULT_BOUNDS = {
    "Ks_mult": (4.0, 10.5, "lin"),
    "f_RS_abs": (0.004, 0.030, "log"),
    "channelconductivity_mmhr": (30.0, 1000.0, "log"),
}
SHORT = {"Ks_mult": "Ks", "f_RS_abs": "f", "channelconductivity_mmhr": "cc"}

REQUIRED = ["run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr",
            "kge_beta", "kge_gamma", "kge_r", "kge_2012", "pbias_pct", "optpercolation"]
SHAPE_COLS = ["kge_gamma", "recession_rate_ratio", "rising_limb_steepness_ratio"]  # optional ones skipped if absent

C_ON, C_OFF = "#0072B2", "#E69F00"  # colorblind-safe blue / orange


# ----------------------------------------------------------------------------------------------
# File discovery / loading
# ----------------------------------------------------------------------------------------------
def locate(name):
    here = Path.cwd().resolve()
    bases = [here] + list(here.parents)[:6]
    for b in bases:
        for pat in (REL_TABLES / name, Path("*") / REL_TABLES / name, Path("*/*") / REL_TABLES / name):
            hits = sorted(b.glob(str(pat)))
            if hits:
                return hits[0]
    return None


def load_series(path, label, expect_perc):
    df = pd.read_csv(path)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        sys.exit(f"[{label}] {path.name} is missing required columns: {missing}")
    n0 = len(df)
    if df["run_id"].duplicated().any():
        sys.exit(f"[{label}] duplicate run_id values in {path.name}: refusing to continue.")
    for c in REQUIRED[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["kge_beta", "kge_gamma", "kge_r", "Ks_mult", "f_RS_abs",
                           "channelconductivity_mmhr"]).copy()
    vals = sorted(df["optpercolation"].dropna().unique().tolist())
    if vals != [expect_perc]:
        print(f"  !! [{label}] optpercolation values = {vals}, expected [{expect_perc}]. "
              f"Check you passed the right file for this series.")
    df["series"] = label
    df["vol"] = (df["kge_beta"] - 1.0).abs()
    df["var"] = (df["kge_gamma"] - 1.0).abs()
    df["cor"] = 1.0 - df["kge_r"]
    df["log10_cc"] = np.log10(df["channelconductivity_mmhr"])
    df["log10_f"] = np.log10(df["f_RS_abs"])
    caution = int(df["caution_zone"].astype(bool).sum()) if "caution_zone" in df.columns else 0
    print(f"  [{label}] {path.name}: {len(df)}/{n0} rows usable, optpercolation={vals}, "
          f"caution-zone flagged={caution}")
    return df.reset_index(drop=True)


# ----------------------------------------------------------------------------------------------
# Pareto machinery
# ----------------------------------------------------------------------------------------------
def non_dominated_layers(P):
    """Non-dominated sorting for minimization. Returns int layer (1 = Pareto front)."""
    n = len(P)
    layer = np.zeros(n, dtype=int)
    remaining = np.arange(n)
    k = 0
    while len(remaining):
        k += 1
        sub = P[remaining]
        keep = []
        for i in remaining:
            dom = np.all(sub <= P[i], axis=1) & np.any(sub < P[i], axis=1)
            if not dom.any():
                keep.append(i)
        layer[keep] = k
        remaining = np.array([i for i in remaining if layer[i] == 0], dtype=int)
    return layer


def near_front_mask(P, layer, eps):
    """True where a run is within eps (every objective) of at least one layer-1 run."""
    F = P[layer == 1]
    return np.array([bool(np.any(np.all(P[i] <= F + eps, axis=1))) for i in range(len(P))])


def eps_tag(e):
    return ("%g" % e).replace(".", "p")


# ----------------------------------------------------------------------------------------------
# Box edges, partial correlation
# ----------------------------------------------------------------------------------------------
def add_edge_flags(df, bounds, frac):
    out_of_box = 0
    for col, (lo, hi, scale) in bounds.items():
        v = df[col].to_numpy(float)
        if scale == "log":
            u = (np.log10(v) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
        else:
            u = (v - lo) / (hi - lo)
        out_of_box += int(np.sum((u < -0.01) | (u > 1.01)))
        df["edge_" + SHORT[col]] = np.where(u < frac, "lo", np.where(u > 1 - frac, "hi", ""))
    df["any_edge"] = (df[["edge_Ks", "edge_f", "edge_cc"]] != "").any(axis=1)
    return out_of_box


def partial_spearman(df, y, x, controls):
    """Spearman partial correlation of y with x given controls (rank-transform, residualize)."""
    sub = df[[y, x] + controls].dropna()
    n, k = len(sub), len(controls)
    if n < k + 5:
        return np.nan, np.nan, n
    R = sub.rank()
    Z = np.column_stack([np.ones(n)] + [R[c].to_numpy() for c in controls])

    def resid(v):
        beta, *_ = np.linalg.lstsq(Z, v, rcond=None)
        return v - Z @ beta

    ry, rx = resid(R[y].to_numpy()), resid(R[x].to_numpy())
    if ry.std() == 0 or rx.std() == 0:
        return np.nan, np.nan, n
    r = float(np.corrcoef(ry, rx)[0, 1])
    p = np.nan
    if HAVE_SCIPY and abs(r) < 1:
        dof = n - 2 - k
        t = r * np.sqrt(dof / (1 - r * r))
        p = float(2 * stats.t.sf(abs(t), dof))
    return r, p, n


# ----------------------------------------------------------------------------------------------
# Summaries
# ----------------------------------------------------------------------------------------------
def summarize_near_front(df, label, eps_list, cc_bounds):
    lo, hi = np.log10(cc_bounds[0]), np.log10(cc_bounds[1])
    third_cut = lo + (hi - lo) / 3.0
    rows = []
    for e in eps_list:
        m = df["nf_" + eps_tag(e)].to_numpy()
        sub = df[m]
        k = len(sub)
        lcc = sub["log10_cc"].to_numpy()
        n_low = int(np.sum(lcc < third_cut)) if k else 0
        binom_p = ks_p = np.nan
        if HAVE_SCIPY and k:
            try:
                binom_p = float(stats.binomtest(n_low, k, 1.0 / 3.0).pvalue)
            except Exception:
                pass
            if k >= 3 and k < len(df):
                ks_p = float(stats.ks_2samp(lcc, df["log10_cc"].to_numpy()).pvalue)
        rows.append(dict(
            series=label, eps=e, n_members=k, n_total=len(df),
            cc_min=sub["channelconductivity_mmhr"].min() if k else np.nan,
            cc_median=sub["channelconductivity_mmhr"].median() if k else np.nan,
            cc_max=sub["channelconductivity_mmhr"].max() if k else np.nan,
            cc_logspan_frac=((lcc.max() - lcc.min()) / (hi - lo)) if k else np.nan,
            frac_in_lowest_third_cc=(n_low / k) if k else np.nan,
            binom_p_vs_one_third=binom_p, ks_p_vs_all_runs=ks_p,
            Ks_min=sub["Ks_mult"].min() if k else np.nan, Ks_max=sub["Ks_mult"].max() if k else np.nan,
            f_min=sub["f_RS_abs"].min() if k else np.nan, f_max=sub["f_RS_abs"].max() if k else np.nan,
            n_at_box_edge=int(sub["any_edge"].sum()) if k else 0,
        ))
    return pd.DataFrame(rows)


def partial_corr_table(df, label, controls=("Ks_mult", "log10_f")):
    subsets = {
        "all runs": df,
        "volume-matched (|PBIAS|<2%)": df[df["pbias_pct"].abs() < 2.0],
    }
    rows = []
    for sname, sub in subsets.items():
        for y in [c for c in SHAPE_COLS if c in df.columns]:
            r, p, n = partial_spearman(sub, y, "log10_cc", list(controls))
            r0 = sub[[y, "log10_cc"]].dropna().corr(method="spearman").iloc[0, 1] if len(sub) > 2 else np.nan
            rows.append(dict(series=label, subset=sname, metric=y, n=n,
                             spearman_raw=r0, partial_spearman_given_Ks_logf=r, partial_p=p))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------------------------
def make_figure(on, off, bounds, eps_plot, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    floor = 1e-4
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))

    # A: objective space
    ax = axs[0]
    for df, col, name, mk in ((off, C_OFF, "cc OFF (control, null)", "^"), (on, C_ON, "cc ON (Stage 2)", "o")):
        x = np.maximum(df["vol"], floor); y = np.maximum(df["var"], floor)
        ax.scatter(x, y, s=14, c=col, alpha=0.30, marker=mk, linewidths=0)
        f1 = df["layer"] == 1
        ax.scatter(x[f1], y[f1], s=60, c=col, marker=mk, edgecolors="k", linewidths=0.9,
                   label=f"{name}: layer 1 (n={int(f1.sum())})")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"volume error  $|\beta-1|$   (floored at 1e-4)")
    ax.set_ylabel(r"variability error  $|\gamma-1|$   (floored at 1e-4)")
    ax.set_title("A. Objective space (minimize both)")
    ax.legend(fontsize=8, loc="lower right")

    # B: Ks-f map of cc-ON near-front runs, colored by cc
    ax = axs[1]
    ax.scatter(on["Ks_mult"], on["f_RS_abs"], s=10, c="0.80", linewidths=0, label="all cc-ON runs")
    m = on["nf_" + eps_tag(eps_plot)]
    sc = ax.scatter(on.loc[m, "Ks_mult"], on.loc[m, "f_RS_abs"], s=48, c=on.loc[m, "log10_cc"],
                    cmap="viridis", edgecolors="k", linewidths=0.4,
                    label=f"near-front, eps={eps_plot:g} (n={int(m.sum())})")
    f1 = on["layer"] == 1
    ax.scatter(on.loc[f1, "Ks_mult"], on.loc[f1, "f_RS_abs"], s=130, facecolors="none",
               edgecolors="k", linewidths=1.6, label="layer 1")
    (klo, khi, _), (flo, fhi, _) = bounds["Ks_mult"], bounds["f_RS_abs"]
    ax.add_patch(plt.Rectangle((klo, flo), khi - klo, fhi - flo, fill=False, ls="--", ec="r", lw=1.2))
    ax.set_yscale("log")
    ax.set_xlabel("Ks multiplier"); ax.set_ylabel("f (RS, abs)")
    ax.set_title("B. Where cc-ON near-front runs sit (red = sampled box)")
    cb = plt.colorbar(sc, ax=ax); cb.set_label("log10 cc (mm/hr)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=1, frameon=False)

    # C: ECDFs of log10 cc
    ax = axs[2]
    for df, col, name in ((on, C_ON, "cc ON"), (off, C_OFF, "cc OFF")):
        for sel, ls, lab in ((np.ones(len(df), bool), ":", "all runs"),
                             (df["nf_" + eps_tag(eps_plot)].to_numpy(), "-", f"near-front (eps={eps_plot:g})")):
            v = np.sort(df.loc[sel, "log10_cc"].to_numpy())
            if len(v):
                ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=col, ls=ls, lw=2,
                        label=f"{name}, {lab}")
    ax.set_xlabel("log10 cc (mm/hr)"); ax.set_ylabel("ECDF")
    ax.set_title("C. cc among near-front runs vs all (OFF = null)")
    ax.legend(fontsize=8, loc="lower right")

    fig.suptitle("PRELIMINARY: single storm (Aug 12 2014), n=250, routing at truth; truth has cc OFF; "
                 "sampled box clips the ridge", fontsize=10, color="firebrick")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage2_csv", type=Path, default=None)
    ap.add_argument("--control_csv", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--eps", type=float, nargs="+", default=[0.0, 0.005, 0.01, 0.02],
                    help="near-front tolerances on the (vol, var) objectives")
    ap.add_argument("--eps_plot", type=float, default=0.01, help="tolerance used in the figure (must be in --eps)")
    ap.add_argument("--edge_frac", type=float, default=0.05, help="box-edge band as a fraction of range")
    ap.add_argument("--ks_range", type=float, nargs=2, default=None)
    ap.add_argument("--f_range", type=float, nargs=2, default=None)
    ap.add_argument("--cc_range", type=float, nargs=2, default=None)
    args = ap.parse_args()

    if args.eps_plot not in args.eps:
        args.eps.append(args.eps_plot)
    args.eps = sorted(set(args.eps))

    p_on = args.stage2_csv or locate(STAGE2_NAME)
    p_off = args.control_csv or locate(CONTROL_NAME)
    if not p_on or not p_off or not Path(p_on).exists() or not Path(p_off).exists():
        sys.exit("Could not find the input CSVs. Run from inside the SMF repo or pass "
                 "--stage2_csv and --control_csv explicitly.")
    p_on, p_off = Path(p_on).resolve(), Path(p_off).resolve()
    out_dir = (args.out_dir or (p_on.parent / "pareto_cc_110")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    outs = {k: out_dir / f"cc_pareto_110_{k}" for k in
            ("runs.csv", "summary.csv", "partial_corr.csv", "overview.png")}
    for o in outs.values():
        assert o not in (p_on, p_off), f"refusing to overwrite an input: {o}"

    bounds = dict(DEFAULT_BOUNDS)
    if args.ks_range: bounds["Ks_mult"] = (*args.ks_range, "lin")
    if args.f_range: bounds["f_RS_abs"] = (*args.f_range, "log")
    if args.cc_range: bounds["channelconductivity_mmhr"] = (*args.cc_range, "log")

    print("=" * 96)
    print("PRELIMINARY component Pareto: cc ON (Stage 2) vs cc OFF (control null)")
    print("=" * 96)
    print("Inputs (read-only):")
    on = load_series(p_on, "cc ON", expect_perc=1.0)
    off = load_series(p_off, "cc OFF", expect_perc=0.0)

    # Pairing check (same sample points => exact twins)
    key = ["Ks_mult", "f_RS_abs", "channelconductivity_mmhr"]
    A, B = off[key].to_numpy(float), on[key].to_numpy(float)
    matched = sum(bool(np.any(np.all(np.abs(B - a) <= 1e-3 * np.abs(a), axis=1))) for a in A)
    print(f"  pairing: {matched}/{len(off)} control points have a Stage 2 twin (rtol 1e-3)")
    if matched < 0.9 * len(off):
        print("  !! Under 90% matched: the two series are NOT paired; the cc-off null is weaker.")

    for df in (on, off):
        oob = add_edge_flags(df, bounds, args.edge_frac)
        if oob:
            print(f"  !! [{df['series'].iloc[0]}] {oob} sampled values lie outside the assumed box. "
                  f"The bounds in this script are stale; pass --ks_range/--f_range/--cc_range.")

    # Layers and near-front sets
    for df in (on, off):
        P = df[["vol", "var"]].to_numpy(float)
        df["layer"] = non_dominated_layers(P)
        for e in args.eps:
            df["nf_" + eps_tag(e)] = near_front_mask(P, df["layer"].to_numpy(), e)

    # Report: layers
    print("\n--- Non-dominated layers on (|beta-1|, |gamma-1|) ---")
    for df in (on, off):
        counts = df["layer"].value_counts().sort_index()
        print(f"  {df['series'].iloc[0]}: {len(counts)} layers; sizes of first 5 = "
              f"{counts.head(5).tolist()}")

    cols = ["Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "vol", "var", "cor", "kge_2012",
            "pbias_pct", "edge_Ks", "edge_f", "edge_cc"]
    for df in (on, off):
        print(f"\n--- {df['series'].iloc[0]}: layer-1 (Pareto front) members, sorted by vol + var ---")
        f1 = df[df["layer"] == 1].copy()
        f1["_s"] = f1["vol"] + f1["var"]
        print(f1.sort_values("_s")[cols].round(4).to_string(index=False))

    # Summaries
    summ = pd.concat([summarize_near_front(df, df["series"].iloc[0], args.eps,
                                           bounds["channelconductivity_mmhr"][:2]) for df in (on, off)],
                     ignore_index=True)
    print("\n--- Near-front sets: how cc is distributed (OFF is the null; expect ~1/3 in lowest third) ---")
    show = ["series", "eps", "n_members", "cc_min", "cc_median", "cc_max", "cc_logspan_frac",
            "frac_in_lowest_third_cc", "binom_p_vs_one_third", "ks_p_vs_all_runs", "n_at_box_edge"]
    print(summ[show].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    pc = pd.concat([partial_corr_table(df, df["series"].iloc[0]) for df in (on, off)], ignore_index=True)
    print("\n--- Shape metrics vs log10(cc), controlling for Ks and log10(f) (Spearman partial) ---")
    print("    (the volume-matched subsets are small, n ~ 20-30: treat p-values as indicative only)")
    print(pc.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    # Save
    keep = ["series", "run_id", "Ks_mult", "f_RS_abs", "channelconductivity_mmhr", "kge_2012", "pbias_pct",
            "kge_beta", "kge_gamma", "kge_r", "vol", "var", "cor", "layer"] + \
           ["nf_" + eps_tag(e) for e in args.eps] + ["edge_Ks", "edge_f", "edge_cc", "any_edge"]
    pd.concat([on[keep], off[keep]], ignore_index=True).to_csv(outs["runs.csv"], index=False)
    summ.to_csv(outs["summary.csv"], index=False)
    pc.to_csv(outs["partial_corr.csv"], index=False)
    make_figure(on, off, bounds, args.eps_plot, outs["overview.png"])

    print("\nWrote:")
    for o in outs.values():
        print("  ", o)

    print("\nHow to read this (and what NOT to conclude):")
    print("  * cc OFF is the null. If cc ON near-front runs bunch at one end of cc while cc OFF ones")
    print("    spread (binom_p / ks_p small for ON, large for OFF), cc is constrained on this storm.")
    print("  * Front members at box edges (n_at_box_edge, edge_* columns): the ridge likely continues")
    print("    outside the sampled box, so fronts and near-front spans are clipped.")
    print("  * The truth has cc OFF: the cc-ON front measures how well Ks/f OFFSET an unwanted cc,")
    print("    not how identifiable cc is in a realistic setting. Single storm, n=250: PRELIMINARY.")
    print("  * Routing (cv, r, n) was pinned at truth; shape signatures may be absorbed once they are free.")


if __name__ == "__main__":
    main()
