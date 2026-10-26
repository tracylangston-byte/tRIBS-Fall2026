"""
plot_series120_ridge_old_vs_new.py
==================================
Makes the two-panel figure "Same 400 draws, two builds: where the model volume
matches the truth".

    left panel   old build  (tRIBS 5.3.1, Series 100)
    right panel  new build  (tRIBS 6.0.0, Series 120)

Each dot is one run, placed at its Ks multiplier (x) and f (y, log scale) and
coloured by PBIAS (percent volume error versus the synthetic truth):
    blue = the run made too little water, red = too much, near-white = about right.
A black line is the "ridge": a smooth curve fitted through all 400 runs, drawn
where it predicts zero volume error. Circled dots are runs within 2% of the
truth volume. The gold star is the truth (Ks 7.0, f 0.012).

What this script does NOT do: it does not run tRIBS and it does not change any
file you already have. It only READS two results CSVs and WRITES one PNG.

HOW TO RUN (from the lab/ folder, like your other scripts)
    python plot_series120_ridge_old_vs_new.py

It looks for these two files in calibration_work/03_comparisons/summary_tables/
    lhs_results_synth_Ks_f_120.csv            (new build, 400 rows)
    lhs_results_synth_Ks_f_100_RESCORED.csv   (old build, 400 rows)
and saves the picture to
    calibration_work/03_comparisons/sensitivity_plots/series120/
        fig_series120_ridge_old_vs_new.png

OPTIONS (all optional)
    --new_csv PATH   use a different new-build CSV
    --old_csv PATH   use a different old-build CSV
    --out_dir PATH   save the PNG somewhere else
    --dpi 150        picture sharpness (72 to 600)

HOW THE RIDGE LINE IS MADE (so you can explain it)
    1. For every run, take the log of the volume ratio: log(simulated / true volume).
       That is log(1 + PBIAS/100). It is 0 when the volume is exactly right.
    2. Fit a smooth degree-4 polynomial surface to those 400 numbers, as a
       function of Ks and log10(f). (One fit per build.)
    3. For each f going up the page, find the Ks where the fitted surface crosses 0.
       Joining those points gives the black line.
    It is a smooth summary, not a measurement. Near the edges of the box the fit
    is less trustworthy, and the old build's fit is rougher than the new one
    (see Series120_Storm100_Result_v1.md).

THE TWO HAND-WRITTEN NOTES under the panel titles (NOTE_OLD, NOTE_NEW below)
    are typed in by me from looking at the picture, not computed. If you re-run
    this on different data, check they are still true or edit them.
"""

import argparse
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # draw to a file; no window needed (Codespaces has no screen)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

# DejaVu Sans ships with matplotlib, so the picture looks the same on any machine.
plt.rcParams["font.family"] = "DejaVu Sans"

# ======================================================================
# SETTINGS you might want to change
# ======================================================================
TRUTH_KS = 7.0                 # synthetic truth
TRUTH_F = 0.012
KS_BOX = (3.0, 11.0)           # the sampled Ks range (the "box")
F_BOX = (0.003, 0.05)          # the sampled f range
FEASIBLE_PCT = 2.0             # "within 2% of the truth volume"
COLOR_CLIP_PCT = 60.0          # PBIAS beyond +/- this is drawn in the darkest colour
POLY_DEGREE = 4                # smoothness of the fitted surface

# The new build's fitted ridge is drawn only up to this f. Above it the
# polynomial starts to curl back on itself, which is the fit's fault, not the
# model's. An arrow marks that the ridge carries on past the Ks = 11 edge.
NEW_RIDGE_MAX_F = 0.042

NOTE_OLD = "Ridge curves back down above f 0.02"
NOTE_NEW = "Ridge keeps rising; leaves the Ks box near f 0.04 to 0.045"

OLD_NAME = "lhs_results_synth_Ks_f_100_RESCORED.csv"
NEW_NAME = "lhs_results_synth_Ks_f_120.csv"
OUT_SUBFOLDER = "series120"
OUT_NAME = "fig_series120_ridge_old_vs_new.png"

# Colours (light-mode chart surface, grey text, one diverging blue-to-red ramp)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e1e0d9"
GOLD = "#f0a30a"
CMAP = LinearSegmentedColormap.from_list(
    "pbias_div", ["#184f95", "#5598e7", "#f0efec", "#ee8c7f", "#e34948", "#8f1f1f"]
)


# ======================================================================
# Loading
# ======================================================================
def load_runs(path):
    """Read one results CSV and keep the columns this figure needs."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Could not find {path}\n"
            "Run this from the lab/ folder, or point to the file with --new_csv / --old_csv."
        )
    df = pd.read_csv(path)
    need = ["Ks_mult", "f_RS_abs", "pbias_pct"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name} is missing column(s) {missing}.")
    df = df.dropna(subset=need).copy()
    # log(1 + PBIAS/100) needs PBIAS above -100%; a run that made no water at all
    # would break the fit, so it is left out of the FIT (still drawn as a dot).
    df["log_ratio"] = np.log1p(np.clip(df["pbias_pct"].values, -99.0, None) / 100.0)
    df["log10_f"] = np.log10(df["f_RS_abs"].values)
    return df


# ======================================================================
# The smooth fit and the ridge
# ======================================================================
def design_matrix(ks, log10_f, degree):
    """Columns 1, x, z, x^2, x*z, z^2 ... up to the given total degree,
    where x and z are Ks and log10(f) rescaled to be of order 1."""
    x = (np.asarray(ks) - TRUTH_KS) / 4.0
    z = (np.asarray(log10_f) - np.log10(TRUTH_F)) / 0.6
    terms = [x**p * z**q for p, q in product(range(degree + 1), repeat=2) if 0 < p + q <= degree]
    return np.column_stack([np.ones_like(x)] + terms)


def fit_surface(df, degree=POLY_DEGREE):
    """Least-squares fit of log volume ratio versus (Ks, log10 f)."""
    A = design_matrix(df["Ks_mult"].values, df["log10_f"].values, degree)
    coef, *_ = np.linalg.lstsq(A, df["log_ratio"].values, rcond=None)
    return coef


def ridge_line(coef, f_max=None, degree=POLY_DEGREE, n_points=300):
    """For each f (log-spaced), find the Ks where the fitted surface equals 0.
    Keeps only f values with exactly one crossing inside the Ks box.
    Returns arrays (ks, f)."""
    ks_grid = np.linspace(KS_BOX[0], KS_BOX[1], 801)
    f_values = np.geomspace(F_BOX[0], F_BOX[1] if f_max is None else f_max, n_points)
    ks_out, f_out = [], []
    for f in f_values:
        y = design_matrix(ks_grid, np.full_like(ks_grid, np.log10(f)), degree) @ coef
        flips = np.where(np.sign(y[:-1]) * np.sign(y[1:]) < 0)[0]
        if len(flips) != 1:
            continue
        i = flips[0]
        t = y[i] / (y[i] - y[i + 1])               # straight-line interpolation to the zero
        ks_out.append(ks_grid[i] + t * (ks_grid[i + 1] - ks_grid[i]))
        f_out.append(f)
    return np.array(ks_out), np.array(f_out)


def ridge_ks_at(ks_line, f_line, f_target):
    """Ks on the ridge at one f (for the printed check table), or nan."""
    if len(f_line) == 0 or not (f_line.min() <= f_target <= f_line.max()):
        return float("nan")
    return float(np.interp(np.log10(f_target), np.log10(f_line), ks_line))


# ======================================================================
# Drawing
# ======================================================================
def style_axes(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c3c2b7")
    ax.grid(True, color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_SOFT, labelsize=11)
    ax.set_xlim(KS_BOX[0] - 0.25, KS_BOX[1] + 0.55)
    ax.set_ylim(F_BOX[0] / 1.07, F_BOX[1] * 1.10)
    ax.set_yscale("log")
    ticks = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:g}" for t in ticks])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xticks(range(3, 12))
    ax.set_xlabel("Ks multiplier", color=INK_SOFT, fontsize=12)


def draw_panel(ax, df, title, note, ridge_xy, arrow=False):
    style_axes(ax)
    norm = Normalize(-COLOR_CLIP_PCT, COLOR_CLIP_PCT)
    ax.scatter(df["Ks_mult"], df["f_RS_abs"], c=df["pbias_pct"].clip(-COLOR_CLIP_PCT, COLOR_CLIP_PCT),
               cmap=CMAP, norm=norm, s=36, linewidths=0, zorder=3)

    feas = df[df["pbias_pct"].abs() < FEASIBLE_PCT]
    ax.scatter(feas["Ks_mult"], feas["f_RS_abs"], s=150, facecolors="none",
               edgecolors=INK, linewidths=1.6, zorder=5)

    ks_line, f_line = ridge_xy
    ax.plot(ks_line, f_line, color=INK, linewidth=2.4, zorder=4, solid_capstyle="round")
    ax.scatter([TRUTH_KS], [TRUTH_F], marker="*", s=420, color=GOLD, edgecolors=INK,
               linewidths=1.2, zorder=6)

    if arrow:
        ax.annotate("", xy=(11.4, 0.0475), xytext=(10.7, 0.0445),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6), zorder=6)

    ax.text(0.0, 1.14, title, transform=ax.transAxes, fontsize=14, fontweight="bold",
            color=INK, va="bottom")
    ax.text(0.0, 1.015, f"{len(feas)} of {len(df)} runs within {FEASIBLE_PCT:g}% of the truth volume\n{note}",
            transform=ax.transAxes, fontsize=11, color=INK_SOFT, va="bottom", linespacing=1.3)
    return len(feas)


def make_figure(new_df, old_df, out_png, dpi):
    fig, (ax_old, ax_new) = plt.subplots(1, 2, figsize=(13.4, 7.4), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    fig.subplots_adjust(left=0.075, right=0.865, top=0.78, bottom=0.17, wspace=0.07)

    old_ridge = ridge_line(fit_surface(old_df))
    new_ridge = ridge_line(fit_surface(new_df), f_max=NEW_RIDGE_MAX_F)

    draw_panel(ax_old, old_df, "Old build: tRIBS 5.3.1 (Series 100)", NOTE_OLD, old_ridge)
    draw_panel(ax_new, new_df, "New build: tRIBS 6.0.0 (Series 120)", NOTE_NEW, new_ridge, arrow=True)
    ax_old.set_ylabel("f (1/mm, log scale)", color=INK_SOFT, fontsize=12)
    ax_new.tick_params(labelleft=False)

    fig.text(0.075, 0.955, "Same 400 draws, two builds: where the model volume matches the truth",
             fontsize=17, fontweight="bold", color=INK, ha="left", va="center")

    # colour bar
    cax = fig.add_axes([0.895, 0.30, 0.014, 0.40])
    sm = plt.cm.ScalarMappable(norm=Normalize(-COLOR_CLIP_PCT, COLOR_CLIP_PCT), cmap=CMAP)
    cb = fig.colorbar(sm, cax=cax, extend="both", ticks=[-60, -30, 0, 30, 60])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=INK_SOFT, labelsize=11)
    cb.set_label("PBIAS (%), volume error vs truth\nblue: too little water, red: too much (clipped at 60)",
                 color=INK_SOFT, fontsize=10)

    handles = [
        Line2D([], [], marker="*", color="none", markerfacecolor=GOLD, markeredgecolor=INK,
               markersize=17, label=f"truth (Ks {TRUTH_KS:g}, f {TRUTH_F:g})"),
        Line2D([], [], color=INK, linewidth=2.4, label="zero-bias ridge (smooth fit through the runs)"),
        Line2D([], [], marker="o", color="none", markerfacecolor="none", markeredgecolor=INK,
               markeredgewidth=1.6, markersize=12, label=f"run within {FEASIBLE_PCT:g}% of the truth volume"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=11,
               labelcolor=INK_SOFT, bbox_to_anchor=(0.47, 0.02))

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, facecolor=fig.get_facecolor())
    plt.close(fig)
    return old_ridge, new_ridge


# ======================================================================
# Main
# ======================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--new_csv", type=Path, default=None)
    ap.add_argument("--old_csv", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")

    # Same folder habit as your other scripts: run from lab/, data sits in ../calibration_work
    script_dir = Path.cwd()
    calib_dir = script_dir.parent / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    out_dir = args.out_dir or (calib_dir / "03_comparisons" / "sensitivity_plots" / OUT_SUBFOLDER)

    new_df = load_runs(args.new_csv or summary_dir / NEW_NAME)
    old_df = load_runs(args.old_csv or summary_dir / OLD_NAME)
    print(f"Loaded {len(new_df)} new-build runs and {len(old_df)} old-build runs.")
    if len(new_df) != len(old_df):
        print("  NOTE: the two builds do not have the same number of runs. The figure still draws, "
              "but it is no longer 'the same draws on both builds'.")

    out_png = out_dir / OUT_NAME
    old_ridge, new_ridge = make_figure(new_df, old_df, out_png, args.dpi)

    # A small table so you can check the picture against the numbers
    print("\nKs on the fitted ridge (blank = the fit has no single crossing there, or it is off the plot):")
    print(f"  {'f':>8}  {'old build':>10}  {'new build':>10}")
    for f in (0.003, 0.005, 0.0075, 0.012, 0.02, 0.03, 0.04):
        o = ridge_ks_at(*old_ridge, f)
        n = ridge_ks_at(*new_ridge, f)
        fmt = lambda v: f"{v:10.2f}" if np.isfinite(v) else f"{'':>10}"
        print(f"  {f:8.4f}  {fmt(o)}  {fmt(n)}")
    n_old = int((old_df["pbias_pct"].abs() < FEASIBLE_PCT).sum())
    n_new = int((new_df["pbias_pct"].abs() < FEASIBLE_PCT).sum())
    print(f"\nRuns within {FEASIBLE_PCT:g}% of the truth volume: old {n_old}, new {n_new}")
    print(f"\nSaved: {out_png}")


if __name__ == "__main__":
    main()
