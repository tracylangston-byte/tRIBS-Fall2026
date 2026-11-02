#!/usr/bin/env python3
"""
plot_pbias_surface_sycamore.py  --  Stage 3: how percent bias (pbias) depends on Ks, f and PsiB.

Reads lhs_volume_summary.csv (the 150-run Latin hypercube on W1, Jan 18-25 2017) and draws the "error surface" of the
flood-week percent bias:

    pbias = 100 * (model volume - gauge volume) / gauge volume        (0 = perfect volume; +20 = model 20% too wet)

Why this is not one single 3-D surface. Each run has THREE inputs (Ks, f, PsiB) and ONE error (pbias). That is four
numbers, and a 3-D plot has only three axes. So the script shows it three ways:

  1. fig_pbias_3d_cube.png      Ks, f and PsiB are the three axes; each run is a dot colored by its pbias (blue = too dry,
                                red = too wet, dots outlined in black are within the band). The gray sheet is where the fitted
                                pbias is exactly zero: the "ridge" of equally good volume. Two viewing angles.
  2. fig_pbias_3d_surface.png   The real 3-D error surface: height = pbias, floor = Ks and f. Three sheets are drawn, one for
                                each of three PsiB values, so you can see PsiB slide the surface up or down.
  3. fig_pbias_contours.png     The same thing as a map (easiest to read numbers from): Ks across, f up, one panel per PsiB.
                                The black line is pbias = 0, the dashed lines are the +/- band. Runs with a PsiB close to the
                                panel's PsiB are shown as dots.
  4. fig_pbias_fit_check.png    The fitted surface against the real run results, so you can see how well it describes them.

How the surface is made. Percent bias is nearly a straight line in the LOGARITHMS of the three knobs, so the script fits
    ln(1 + pbias/100) = a + b ln(Ks) + c ln(f) + d ln(PsiB)
by least squares and draws that. The fit quality (R-squared, and the typical miss in pbias points) is printed and written
on the figures. The dots are always the real model results; the surfaces and contours are the fit.

It also writes pbias_zero_ks_table.csv: the Ks that gives pbias = 0 (and the +/- band edges) for a grid of f and PsiB values.

Run it from workspaces/Sycamore/lab/ (it only reads the CSV and writes picture and table files; it never touches results/):
    python -u plot_pbias_surface_sycamore.py --dry-run     # checks the CSV and the fit, writes nothing
    python -u plot_pbias_surface_sycamore.py               # makes the four figures and the table in lhs_volume_plots/

Options: --csv FILE (default lhs_volume_summary.csv), --metric COLUMN (default flood_week__pbias_pct), --band N (default 10),
--out-dir DIR (default lhs_volume_plots), --dpi N (default 150).
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                      # no screen needed (Codespace); files only
import matplotlib.pyplot as plt            # noqa: E402
import numpy as np                         # noqa: E402
import pandas as pd                        # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

KNOBS = ["KS_MULT", "F_MULT", "PSIB_MULT"]
KS_RANGE, F_RANGE, PSIB_RANGE = (0.15, 0.30), (0.5, 2.0), (0.25, 1.5)
BASELINE = {"KS_MULT": 0.2, "F_MULT": 1.0, "PSIB_MULT": 0.5}
KS_TICKS, F_TICKS, PSIB_TICKS = [0.15, 0.2, 0.25, 0.3], [0.5, 1, 2], [0.25, 0.5, 1, 1.5]
PSIB_PANELS = [0.25, 0.5, 1.0, 1.5]
CMAP = "RdBu_r"                            # blue = model too dry, red = model too wet


def log(msg: str = "") -> None:
    print(msg, flush=True)


# --------------------------------------------------------------------------------------------
# Data and fit
# --------------------------------------------------------------------------------------------
def load(csv_path: Path, metric: str) -> pd.DataFrame:
    if not csv_path.is_file():
        raise SystemExit(f"Cannot find {csv_path}. Run this from workspaces/Sycamore/lab/, or use --csv.")
    df = pd.read_csv(csv_path)
    missing = [c for c in KNOBS + [metric] if c not in df.columns]
    if missing:
        raise SystemExit(f"{csv_path.name} has no column(s) {missing}. Columns that contain 'pbias': "
                         f"{[c for c in df.columns if 'pbias' in c]}")
    if "status" in df.columns:
        df = df[df["status"].isin(["done", "unscored"])]
    df = df.dropna(subset=KNOBS + [metric]).copy()
    if len(df) < 20:
        raise SystemExit(f"Only {len(df)} scored runs in {csv_path.name}; need at least 20 for a fit.")
    if (df[metric] <= -99).any():
        raise SystemExit(f"{metric} has values at or below -99%, which the log fit cannot use.")
    return df.reset_index(drop=True)


def design(ks, f, psib):
    ks, f, psib = np.broadcast_arrays(np.asarray(ks, float), np.asarray(f, float), np.asarray(psib, float))
    return np.column_stack([np.ones(ks.size), np.log(ks).ravel(), np.log(f).ravel(), np.log(psib).ravel()]), ks.shape


def fit_surface(df: pd.DataFrame, metric: str) -> dict:
    X, _ = design(df["KS_MULT"], df["F_MULT"], df["PSIB_MULT"])
    y = np.log1p(df[metric].to_numpy() / 100.0)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = 100.0 * np.expm1(X @ coef)
    resid = df[metric].to_numpy() - pred
    r2 = 1 - float(((df[metric] - pred) ** 2).sum() / ((df[metric] - df[metric].mean()) ** 2).sum())
    return dict(coef=coef, r2=r2, rmse=float(np.sqrt((resid ** 2).mean())), pred=pred, resid=resid)


def predict(coef, ks, f, psib):
    X, shape = design(ks, f, psib)
    return (100.0 * np.expm1(X @ coef)).reshape(shape)


def ks_for(coef, target_pbias: float, f, psib):
    """The Ks at which the fitted pbias equals target_pbias, for the given f and PsiB."""
    a, b, c, d = coef
    return np.exp((math.log1p(target_pbias / 100.0) - a - c * np.log(f) - d * np.log(psib)) / b)


# --------------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------------
def norm_for(values, band):
    lim = max(float(np.nanmax(np.abs(values))), 2 * band)
    return TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)


def set_log_ticks(ax, axis, ticks):
    vals = [math.log10(t) for t in ticks]
    labels = [f"{t:g}" for t in ticks]
    getattr(ax, f"set_{axis}ticks")(vals)
    getattr(ax, f"set_{axis}ticklabels")(labels)


def fig_cube(df, fit, metric, band, out, dpi):
    p = df[metric].to_numpy()
    norm = norm_for(p, band)
    fig = plt.figure(figsize=(15, 7))
    fg, pg = np.meshgrid(np.geomspace(*F_RANGE, 25), np.geomspace(*PSIB_RANGE, 25))
    ks0 = np.clip(ks_for(fit["coef"], 0.0, fg, pg), *KS_RANGE)
    inband = np.abs(p) <= band
    for i, (elev, azim) in enumerate([(22, -62), (24, 38)]):
        ax = fig.add_subplot(1, 2, i + 1, projection="3d")
        ax.plot_surface(np.log10(ks0), np.log10(fg), np.log10(pg), color="gray", alpha=0.25, linewidth=0, shade=False)
        sc = ax.scatter(np.log10(df["KS_MULT"][~inband]), np.log10(df["F_MULT"][~inband]), np.log10(df["PSIB_MULT"][~inband]),
                        c=p[~inband], cmap=CMAP, norm=norm, s=34, edgecolor="none", depthshade=False)
        ax.scatter(np.log10(df["KS_MULT"][inband]), np.log10(df["F_MULT"][inband]), np.log10(df["PSIB_MULT"][inband]),
                   c=p[inband], cmap=CMAP, norm=norm, s=46, edgecolor="black", linewidth=1.1, depthshade=False)
        ax.scatter([math.log10(BASELINE["KS_MULT"])], [math.log10(BASELINE["F_MULT"])], [math.log10(BASELINE["PSIB_MULT"])],
                   marker="*", s=260, c="gold", edgecolor="black", linewidth=1.0, depthshade=False)
        set_log_ticks(ax, "x", KS_TICKS)
        set_log_ticks(ax, "y", F_TICKS)
        set_log_ticks(ax, "z", PSIB_TICKS)
        ax.set_xlabel("Ks multiplier", labelpad=8)
        ax.set_ylabel("f multiplier", labelpad=8)
        ax.set_zlabel("PsiB multiplier", labelpad=16)
        ax.set_xlim(math.log10(KS_RANGE[0]), math.log10(KS_RANGE[1]))
        ax.set_ylim(math.log10(F_RANGE[0]), math.log10(F_RANGE[1]))
        ax.set_zlim(math.log10(PSIB_RANGE[0]), math.log10(PSIB_RANGE[1]))
        ax.view_init(elev=elev, azim=azim)
    cax = fig.add_axes([0.92, 0.2, 0.012, 0.6])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("flood-week percent bias (%)\nblue = model too dry, red = too wet")
    fig.suptitle(f"The {len(df)} runs in Ks-f-PsiB space, colored by percent bias (two views of the same cloud).\n"
                 f"Black outline: within +/-{band:g}%.  Gray sheet: fitted pbias = 0.  Gold star: baseline.", fontsize=12)
    fig.subplots_adjust(left=0.03, right=0.88, top=0.92, bottom=0.04, wspace=0.12)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


def fig_surface(df, fit, metric, band, out, dpi):
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="3d")
    kg, fg = np.meshgrid(np.geomspace(*KS_RANGE, 40), np.geomspace(*F_RANGE, 40))
    colors = {0.25: "#2166ac", 0.5: "#555555", 1.5: "#b2182b"}
    lo, hi = -60.0, 60.0
    ax.plot_surface(np.log10(kg), np.log10(fg), np.zeros_like(kg), color="black", alpha=0.10, linewidth=0, shade=False)
    for psib, col in colors.items():
        z = np.clip(predict(fit["coef"], kg, fg, psib), lo, hi)
        ax.plot_surface(np.log10(kg), np.log10(fg), z, color=col, alpha=0.38, linewidth=0.3, edgecolor=col, rstride=4, cstride=4, shade=False)
        ax.contour(np.log10(kg), np.log10(fg), z, levels=[0], zdir="z", offset=lo, colors=[col], linewidths=2.2)
    ax.set_zlim(lo, hi)
    set_log_ticks(ax, "x", KS_TICKS)
    set_log_ticks(ax, "y", F_TICKS)
    ax.set_xlabel("Ks multiplier", labelpad=10)
    ax.set_ylabel("f multiplier", labelpad=10)
    ax.set_zlabel("flood-week percent bias (%)", labelpad=8)
    ax.view_init(elev=20, azim=-48)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], color=c, lw=6, alpha=0.6, label=f"PsiB x{p:g}") for p, c in colors.items()]
              + [Line2D([0], [0], color="black", lw=1, alpha=0.5, label="pbias = 0 (floor-level sheet)")],
              loc="upper left", fontsize=10)
    ax.set_title("Fitted error surface: percent bias over Ks and f, for three PsiB values\n"
                 "(colored lines on the floor: where each surface crosses pbias = 0)", fontsize=12)
    fig.text(0.5, 0.01, f"Fit: R-squared {fit['r2']:.3f}; typical miss {fit['rmse']:.1f} percentage points ({len(df)} runs). "
                        f"Surface heights are clipped at +/-60%.", ha="center", fontsize=9)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def fig_contours(df, fit, metric, band, out, dpi):
    p = df[metric].to_numpy()
    norm = norm_for(p, band)
    kg, fg = np.meshgrid(np.geomspace(*KS_RANGE, 200), np.geomspace(*F_RANGE, 200))
    fig, axes = plt.subplots(2, 2, figsize=(12, 9.5), sharex=True, sharey=True)
    near_tol = math.log(1.35)                       # a run is "near" a panel if its PsiB is within a factor of 1.35
    levels = np.linspace(norm.vmin, norm.vmax, 41)
    for ax, psib in zip(axes.ravel(), PSIB_PANELS):
        z = np.clip(predict(fit["coef"], kg, fg, psib), norm.vmin, norm.vmax)
        cf = ax.contourf(kg, fg, z, levels=levels, cmap=CMAP, norm=norm, alpha=0.85)
        ax.contour(kg, fg, z, levels=[0], colors="black", linewidths=2.4)
        ax.contour(kg, fg, z, levels=[-band, band], colors="black", linewidths=1.2, linestyles="--")
        near = np.abs(np.log(df["PSIB_MULT"].to_numpy() / psib)) <= near_tol
        ax.scatter(df["KS_MULT"][near], df["F_MULT"][near], c=p[near], cmap=CMAP, norm=norm, s=60, edgecolor="black", linewidth=0.8)
        ax.scatter([BASELINE["KS_MULT"]], [BASELINE["F_MULT"]], marker="*", s=240, c="gold", edgecolor="black", linewidth=1.0,
                   label="baseline (PsiB 0.5 panel)" if psib == 0.5 else None)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(KS_TICKS)
        ax.set_xticklabels([f"{t:g}" for t in KS_TICKS])
        ax.set_yticks(F_TICKS)
        ax.set_yticklabels([f"{t:g}" for t in F_TICKS])
        ax.minorticks_off()
        ax.set_title(f"PsiB x{psib:g}   ({int(near.sum())} nearby runs shown)")
        ax.grid(alpha=0.25)
        if psib == 0.5:
            ax.legend(loc="lower left", fontsize=9)
    for ax in axes[1]:
        ax.set_xlabel("Ks multiplier")
    for ax in axes[:, 0]:
        ax.set_ylabel("f multiplier")
    cb = fig.colorbar(cf, ax=axes.ravel().tolist(), shrink=0.85, pad=0.02)
    cb.set_label("flood-week percent bias (%): blue = model too dry, red = too wet")
    fig.suptitle(f"Percent bias map. Solid line: pbias = 0.  Dashed: +/-{band:g}%.  Dots: real runs whose PsiB is within a factor of 1.35 "
                 f"of the panel's PsiB (same colors).  Fit R-squared {fit['r2']:.2f}.", fontsize=11)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def fig_fitcheck(df, fit, metric, band, out, dpi):
    p = df[metric].to_numpy()
    fig, ax = plt.subplots(figsize=(6.4, 6))
    sc = ax.scatter(fit["pred"], p, c=np.log10(df["KS_MULT"]), cmap="viridis", s=40, edgecolor="none")
    lim = [min(p.min(), fit["pred"].min()) - 3, max(p.max(), fit["pred"].max()) + 3]
    ax.plot(lim, lim, color="black", lw=1)
    ax.fill_between(lim, [l - 3 for l in lim], [l + 3 for l in lim], color="gray", alpha=0.15, label="+/- 3 points")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("fitted percent bias (%)")
    ax.set_ylabel("percent bias of the real run (%)")
    ax.set_title(f"How well the surface describes the runs\nR-squared {fit['r2']:.3f}; typical miss {fit['rmse']:.1f} points; near 0 the fit is tighter")
    cb = fig.colorbar(sc)
    cb.set_ticks([math.log10(t) for t in KS_TICKS])
    cb.set_ticklabels([f"{t:g}" for t in KS_TICKS])
    cb.set_label("Ks multiplier")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.25)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_zero_table(fit, band, out) -> pd.DataFrame:
    rows = []
    for f in (0.5, 0.75, 1.0, 1.5, 2.0):
        for psib in (0.25, 0.5, 1.0, 1.5):
            rows.append(dict(F_MULT=f, PSIB_MULT=psib,
                             KS_pbias_minus_band=round(float(ks_for(fit["coef"], -band, f, psib)), 4),
                             KS_pbias_zero=round(float(ks_for(fit["coef"], 0.0, f, psib)), 4),
                             KS_pbias_plus_band=round(float(ks_for(fit["coef"], band, f, psib)), 4)))
    # Ks goes DOWN as pbias goes UP, so the "minus" column is the larger Ks. Name the columns by what they mean.
    tab = pd.DataFrame(rows).rename(columns={"KS_pbias_minus_band": f"KS_at_pbias_-{band:g}",
                                             "KS_pbias_zero": "KS_at_pbias_0",
                                             "KS_pbias_plus_band": f"KS_at_pbias_+{band:g}"})
    tab["inside_sampled_Ks_range"] = tab["KS_at_pbias_0"].between(*KS_RANGE)
    tab.to_csv(out, index=False)
    return tab


# --------------------------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="3-D / contour views of flood-week percent bias against Ks, f and PsiB.")
    ap.add_argument("--csv", default="lhs_volume_summary.csv")
    ap.add_argument("--metric", default="flood_week__pbias_pct")
    ap.add_argument("--band", type=float, default=10.0, help="percent bias band to outline (default 10)")
    ap.add_argument("--out-dir", default="lhs_volume_plots")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--dry-run", action="store_true", help="read the CSV and fit the surface, write nothing")
    args = ap.parse_args(argv)

    df = load(Path(args.csv), args.metric)
    fit = fit_surface(df, args.metric)
    b = fit["coef"]
    log(f"Read {len(df)} scored runs from {args.csv}; metric {args.metric} "
        f"(range {df[args.metric].min():.1f} to {df[args.metric].max():.1f} %).")
    log(f"Fit: ln(1 + pbias/100) = {b[0]:+.3f} {b[1]:+.3f} ln(Ks) {b[2]:+.3f} ln(f) {b[3]:+.3f} ln(PsiB)")
    log(f"     R-squared {fit['r2']:.3f}; typical miss {fit['rmse']:.1f} percentage points; worst miss {np.abs(fit['resid']).max():.1f}.")
    log("     Reading the coefficients: doubling Ks multiplies (1 + pbias/100) by "
        f"{2 ** b[1]:.2f}; doubling f by {2 ** b[2]:.2f}; doubling PsiB by {2 ** b[3]:.2f}.")
    inband = int((df[args.metric].abs() <= args.band).sum())
    log(f"     {inband} of {len(df)} runs are within +/-{args.band:g}%.")
    if args.dry_run:
        log("DRY RUN: nothing written.")
        return 0

    out = Path(args.out_dir)
    out.mkdir(exist_ok=True)
    jobs = [("fig_pbias_3d_cube.png", fig_cube), ("fig_pbias_3d_surface.png", fig_surface),
            ("fig_pbias_contours.png", fig_contours), ("fig_pbias_fit_check.png", fig_fitcheck)]
    for name, fn in jobs:
        fn(df, fit, args.metric, args.band, out / name, args.dpi)
        log(f"Wrote {out / name}")
    tab = write_zero_table(fit, args.band, out / "pbias_zero_ks_table.csv")
    log(f"Wrote {out / 'pbias_zero_ks_table.csv'} ({len(tab)} rows). Ks that gives pbias = 0, by f and PsiB:")
    piv = tab.pivot(index="F_MULT", columns="PSIB_MULT", values="KS_at_pbias_0")
    log(piv.to_string())
    outside = int((~tab["inside_sampled_Ks_range"]).sum())
    if outside:
        log(f"({outside} of {len(tab)} table rows fall outside the sampled Ks range {KS_RANGE}; treat those as extrapolation.)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
