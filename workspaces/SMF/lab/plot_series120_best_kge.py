"""
plot_series120_best_kge.py
==========================
Shows WHICH runs scored best on KGE, and where they sit in the Ks-f plane.
Uses only the new-build results (tRIBS 6.0.0, Series 120, storm 1.0x).

KGE (Kling-Gupta efficiency) is one number that scores how well a simulated
hydrograph matches the observed one: 1 is perfect. It combines three things:
    r      correlation (does the shape/timing line up?)
    beta   ratio of the means (is the total volume right?)
    gamma  ratio of the variability (are peaks and lows the right size?)
    KGE = 1 - sqrt( (r-1)^2 + (beta-1)^2 + (gamma-1)^2 )
"kge_2012" is the version that measures variability with the coefficient of
variation (this is your main metric). Here the "observed" hydrograph is the
synthetic truth (Ks 7.0, f 0.012), so a perfect score is possible.

THE FIGURE (three panels and a small table)
    A  every run at its Ks (x) and f (y, log scale), coloured by KGE.
       Grey = scored below 0.90. Blue = 0.90 or better (darker = better).
       Black rings = the best few runs. Gold star = the truth.
       The dashed box is the area shown in panel B.
    B  a zoom on the best few runs, numbered 1, 2, 3 ... in order of score.
    C  KGE against f for every run (the best score at each f is the top edge
       of the cloud). The dotted line marks KGE 0.95.
    table  the same best runs as numbers.

HOW TO RUN (from the lab/ folder, like your other scripts)
    python plot_series120_best_kge.py

Reads   calibration_work/03_comparisons/summary_tables/lhs_results_synth_Ks_f_120.csv
Writes  calibration_work/03_comparisons/sensitivity_plots/series120/
            fig_series120_best_kge.png
            series120_top_kge_runs.csv     (the best 20 runs, all their scores)
It only READS your results file; nothing you already have is changed.

OPTIONS (all optional)
    --csv PATH       a different results CSV
    --out_dir PATH   save somewhere else
    --metric NAME    kge_2012 (default) or kge (the original 2009 version)
    --top N          how many best runs to ring/number in the figure (default 5)
    --floor X        KGE below this is drawn grey (default 0.90)
    --dpi N          picture sharpness (default 150)

THINGS TO KEEP IN MIND WHEN READING IT
    * The best run is the draw that happened to land closest to the truth.
      With perfect (noise-free) synthetic data that is expected. It shows the
      model CAN recover the truth, not that a real calibration would.
    * KGE here is scored on one storm with routing parameters held at truth.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # draw to a file; no window needed (Codespaces has no screen)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

plt.rcParams["font.family"] = "DejaVu Sans"   # ships with matplotlib: same look everywhere

# ======================================================================
# SETTINGS
# ======================================================================
TRUTH_KS = 7.0
TRUTH_F = 0.012
KS_BOX = (3.0, 11.0)
F_BOX = (0.003, 0.05)
LINE_KGE = 0.95                # the dotted reference line in panel C
C_YMIN = 0.80                  # panel C shows KGE from here up to 1.0

CSV_NAME = "lhs_results_synth_Ks_f_120.csv"
OUT_SUBFOLDER = "series120"
PNG_NAME = "fig_series120_best_kge.png"
TOP_CSV_NAME = "series120_top_kge_runs.csv"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e1e0d9"
GREY_DOT = "#cfcec6"
GOLD = "#f0a30a"
# one-hue ramp (light blue -> dark blue) for "better and better"
CMAP = LinearSegmentedColormap.from_list("kge_seq", ["#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
CMAP.set_under(GREY_DOT)       # anything below the floor is grey


# ======================================================================
# Loading and ranking
# ======================================================================
def load_runs(path, metric):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Could not find {path}\nRun this from the lab/ folder, or give the file with --csv.")
    df = pd.read_csv(path)
    need = ["Ks_mult", "f_RS_abs", "pbias_pct", metric]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name} is missing column(s) {missing}.")
    df = df.dropna(subset=need).reset_index(drop=True)
    # rank 1 = best score (ties broken by the order in the file)
    df = df.sort_values(metric, ascending=False, kind="stable").reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    return df


def save_top_table(df, metric, n, out_csv):
    cols = ["rank", "draw_idx", "run_id", "Ks_mult", "f_RS_abs", metric, "kge", "kge_r", "kge_beta",
            "kge_gamma", "nse", "pbias_pct"]
    cols = [c for c in dict.fromkeys(cols) if c in df.columns]      # keep order, drop duplicates/missing
    top = df.loc[: n - 1, cols].copy()
    top["Ks_minus_truth"] = top["Ks_mult"] - TRUTH_KS
    top["f_over_truth"] = top["f_RS_abs"] / TRUTH_F
    top.to_csv(out_csv, index=False)
    return top


# ======================================================================
# Drawing helpers
# ======================================================================
F_TICKS = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05]


def style_axes(ax, logy=True):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c3c2b7")
    ax.grid(True, color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_SOFT, labelsize=10.5)
    if logy:
        ax.set_yscale("log")
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())


def set_f_ticks(ax, axis="y"):
    labels = [f"{t:g}" for t in F_TICKS]
    if axis == "y":
        ax.set_yticks(F_TICKS)
        ax.set_yticklabels(labels)
    else:
        ax.set_xticks(F_TICKS)
        ax.set_xticklabels(labels)
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())


def panel_title(ax, letter, title, sub=None):
    """Title and one grey line under it, placed a fixed number of points above the panel."""
    ax.annotate(f"{letter}   {title}", xy=(0, 1), xycoords="axes fraction", xytext=(0, 27 if sub else 8),
                textcoords="offset points", fontsize=12.5, fontweight="bold", color=INK, va="bottom")
    if sub:
        ax.annotate(sub, xy=(0, 1), xycoords="axes fraction", xytext=(0, 8), textcoords="offset points",
                    fontsize=10, color=INK_SOFT, va="bottom")


def draw_runs(ax, df, norm, size):
    """All runs, coloured by score. Draw the worst first so the best sit on top."""
    d = df.sort_values(METRIC_COL[0], ascending=True)
    ax.scatter(d["Ks_mult"], d["f_RS_abs"], c=d[METRIC_COL[0]], cmap=CMAP, norm=norm, s=size,
               linewidths=0, zorder=3)


METRIC_COL = ["kge_2012"]      # set in main(); lets the helpers know which column to colour by


# ======================================================================
# The figure
# ======================================================================
def make_figure(df, metric, top_n, floor, out_png, dpi):
    METRIC_COL[0] = metric
    label = metric.upper() if metric != "kge" else "KGE (2009)"
    norm = Normalize(vmin=floor, vmax=1.0)
    best = df.iloc[:top_n]

    fig = plt.figure(figsize=(14.6, 12.0))
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1.0], height_ratios=[1, 1],
                          left=0.065, right=0.975, top=0.865, bottom=0.33, wspace=0.17, hspace=0.62)
    axA = fig.add_subplot(gs[:, 0])
    axC = fig.add_subplot(gs[0, 1])
    axB = fig.add_subplot(gs[1, 1])

    n_floor = int((df[metric] >= floor).sum())
    n95 = int((df[metric] >= LINE_KGE).sum())
    n99 = int((df[metric] >= 0.99).sum())

    # ---------------- A: the map -----------------------------------------
    style_axes(axA)
    draw_runs(axA, df, norm, 40)
    axA.scatter(best["Ks_mult"], best["f_RS_abs"], s=190, facecolors="none", edgecolors=INK,
                linewidths=1.8, zorder=5)
    axA.scatter([TRUTH_KS], [TRUTH_F], marker="*", s=430, color=GOLD, edgecolors=INK, linewidths=1.2, zorder=6)
    axA.set_xlim(KS_BOX[0] - 0.25, KS_BOX[1] + 0.25)
    axA.set_ylim(F_BOX[0] / 1.07, F_BOX[1] * 1.10)
    axA.set_xticks(range(3, 12))
    set_f_ticks(axA, "y")
    axA.set_xlabel("Ks multiplier", color=INK_SOFT, fontsize=11.5)
    axA.set_ylabel("f (1/mm, log scale)", color=INK_SOFT, fontsize=11.5)

    # dashed box = the zoom in panel B
    bx0, bx1, by0, by1 = zoom_limits(best)
    axA.add_patch(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fill=False, edgecolor=INK, linestyle=(0, (4, 3)),
                            linewidth=1.3, zorder=7))
    axA.annotate("B", xy=(bx1, by1), xytext=(bx1 + 0.35, by1 * 1.16), fontsize=12, fontweight="bold",
                 color=INK, zorder=8, arrowprops=dict(arrowstyle="-", color=INK, lw=1.0))
    panel_title(axA, "A", "Where the high-scoring runs sit",
                f"{len(df)} runs. {n_floor} score {floor:g} or better (blue); the other {len(df) - n_floor} are grey")

    # colour bar: horizontal, under panel A, to the right of the legend
    cax = axA.inset_axes([0.60, -0.125, 0.38, 0.02])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    ticks = [floor, (floor + 1.0) / 2, 1.0]
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", extend="min", ticks=ticks)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=INK_SOFT, labelsize=10)
    cb.ax.set_xticklabels([f"{t:.2f}" for t in ticks])
    cb.set_label(f"{label}   (grey triangle: below {floor:g})", color=INK_SOFT, fontsize=10)

    # ---------------- C: score versus f -----------------------------------
    style_axes(axC, logy=False)
    axC.set_xscale("log")
    shown = df[df[metric] >= C_YMIN]
    draw_c = shown.sort_values(metric, ascending=True)
    axC.scatter(draw_c["f_RS_abs"], draw_c[metric], c=draw_c[metric], cmap=CMAP, norm=norm, s=34,
                linewidths=0, zorder=3)
    axC.scatter(best["f_RS_abs"], best[metric], s=150, facecolors="none", edgecolors=INK, linewidths=1.6, zorder=5)
    axC.axhline(LINE_KGE, color=INK_SOFT, linestyle=":", linewidth=1.3, zorder=2)
    axC.axvline(TRUTH_F, color=GOLD, linestyle="--", linewidth=1.6, zorder=2)
    axC.set_xlim(F_BOX[0] / 1.07, F_BOX[1] * 1.10)
    axC.set_ylim(C_YMIN, 1.012)
    set_f_ticks(axC, "x")
    axC.set_xlabel("f (1/mm, log scale)", color=INK_SOFT, fontsize=11)
    axC.set_ylabel(label, color=INK_SOFT, fontsize=11)
    axC.text(TRUTH_F * 1.05, C_YMIN + 0.006, "truth f", color="#7a5200", fontsize=10, va="bottom",
             bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.5))
    hi = df[df[metric] >= LINE_KGE]
    n_hidden = int((df[metric] < C_YMIN).sum())
    if len(hi):
        last = hi.loc[hi["f_RS_abs"].idxmax()]
        axC.annotate(f"no run above f {last['f_RS_abs']:.4f}\nreaches {LINE_KGE:g}",
                     xy=(last["f_RS_abs"], last[metric]), xytext=(F_BOX[1] * 1.06, 0.985),
                     fontsize=10, color=INK, va="center", ha="right",
                     arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0))
    panel_title(axC, "C", f"{label} against f",
                f"{n95} runs reach {LINE_KGE:g}, {n99} reach 0.99. ({n_hidden} runs below {C_YMIN:g} not shown)")

    # ---------------- B: zoom on the best few -----------------------------
    style_axes(axB)
    inside = df[(df["Ks_mult"].between(bx0, bx1)) & (df["f_RS_abs"].between(by0, by1))]
    draw_runs(axB, inside, norm, 60)
    axB.scatter([TRUTH_KS], [TRUTH_F], marker="*", s=900, color=GOLD, edgecolors=INK, linewidths=1.3, zorder=5)
    for _, r in best.iterrows():
        axB.scatter([r["Ks_mult"]], [r["f_RS_abs"]], s=380, color=CMAP(norm(r[metric])), edgecolors=INK,
                    linewidths=1.5, zorder=8)
        axB.text(r["Ks_mult"], r["f_RS_abs"], str(int(r["rank"])), ha="center", va="center",
                 fontsize=10.5, fontweight="bold", color="white", zorder=9)
    axB.set_xlim(bx0, bx1)
    axB.set_ylim(by0, by1)
    yt = [t for t in (0.011, 0.012, 0.013, 0.014, 0.015) if by0 <= t <= by1]
    axB.set_yticks(yt)
    axB.set_yticklabels([f"{t:g}" for t in yt])
    axB.set_xlabel("Ks multiplier", color=INK_SOFT, fontsize=11)
    axB.set_ylabel("f (1/mm, log scale)", color=INK_SOFT, fontsize=11)
    panel_title(axB, "B", f"Zoom: the {len(best)} best runs, numbered by rank",
                "star = truth; small dots = the other runs in this window")

    # ---------------- table of the best few (figure coordinates, so rows never overlap) ----
    mono = dict(family="DejaVu Sans Mono", fontsize=11.5, va="center", color=INK, ha="left")
    fig.text(0.065, 0.205, f"The {len(best)} best runs by {label}", fontsize=12.5, fontweight="bold",
             color=INK, va="center")
    head = f"{'rank':>4}  {'Ks':>6}  {'f':>8}  {label:>9}  {'PBIAS %':>8}  {'Ks - 7.0':>9}  {'f / 0.012':>9}"
    fig.text(0.065, 0.177, head, **{**mono, "color": INK_SOFT})
    for i, (_, r) in enumerate(best.iterrows()):
        row = (f"{int(r['rank']):>4}  {r['Ks_mult']:6.2f}  {r['f_RS_abs']:8.4f}  {r[metric]:9.4f}  "
               f"{r['pbias_pct']:+8.2f}  {r['Ks_mult'] - TRUTH_KS:+9.2f}  {r['f_RS_abs'] / TRUTH_F:9.2f}")
        fig.text(0.065, 0.177 - 0.0235 * (i + 1), row, **mono)

    # ---------------- a plain-words reminder of what KGE is (static text, not computed) ----
    kx = 0.54
    fig.text(kx, 0.205, "What goes into KGE", fontsize=12.5, fontweight="bold", color=INK, va="center")
    lines = [
        "KGE = 1 \u2212 \u221a[ (r \u2212 1)\u00b2 + (\u03b2 \u2212 1)\u00b2 + (\u03b3 \u2212 1)\u00b2 ]",
        "r   correlation of simulated and true flow (shape and timing)",
        "\u03b2   simulated mean \u00f7 true mean (volume; PBIAS is the same idea)",
        "\u03b3   variability ratio, simulated \u00f7 true (size of peaks and lows)",
        "1.0 is perfect. Each kind of miss pulls the score down.",
    ]
    for i, line in enumerate(lines):
        fig.text(kx, 0.177 - 0.0235 * i, line, fontsize=11, color=INK if i == 0 else INK_SOFT, va="center")

    # legend + titles
    handles = [
        Line2D([], [], marker="*", color="none", markerfacecolor=GOLD, markeredgecolor=INK, markersize=16,
               label=f"truth (Ks {TRUTH_KS:g}, f {TRUTH_F:g})"),
        Line2D([], [], marker="o", color="none", markerfacecolor="none", markeredgecolor=INK,
               markeredgewidth=1.7, markersize=12, label=f"the {len(best)} best runs"),
    ]
    axA.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.095), ncol=2, frameon=False,
               fontsize=10.5, labelcolor=INK_SOFT)

    fig.text(0.065, 0.958, f"Which runs score best? tRIBS 6.0.0, Series 120, {len(df)} runs, storm 1.0x",
             fontsize=17, fontweight="bold", color=INK, ha="left", va="center")
    fig.text(0.065, 0.926, "Scored against the synthetic truth (so a perfect 1.0 is possible). "
             "Routing parameters are held at the truth values.", fontsize=11, color=INK_SOFT, ha="left", va="center")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, facecolor=fig.get_facecolor())
    plt.close(fig)


def zoom_limits(best):
    """The window for panel B: wide enough for the best runs and the truth, with a margin."""
    xs = np.append(best["Ks_mult"].values, TRUTH_KS)
    fs = np.append(best["f_RS_abs"].values, TRUTH_F)
    pad_x = max(0.3, 0.12 * (xs.max() - xs.min()))
    return xs.min() - pad_x, xs.max() + pad_x, fs.min() / 1.12, fs.max() * 1.12


# ======================================================================
# Main
# ======================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--csv", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--metric", default="kge_2012", choices=["kge_2012", "kge"])
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--floor", type=float, default=0.90)
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    if not (1 <= args.top <= 6):
        ap.error("--top must be between 1 and 6 (the zoom and table get crowded beyond that).")
    if not (0.0 <= args.floor < 1.0):
        ap.error("--floor must be at least 0 and below 1.")
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")

    script_dir = Path.cwd()
    calib_dir = script_dir.parent / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    out_dir = args.out_dir or (calib_dir / "03_comparisons" / "sensitivity_plots" / OUT_SUBFOLDER)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_runs(args.csv or summary_dir / CSV_NAME, args.metric)
    print(f"Loaded {len(df)} runs. Ranking by {args.metric}.")

    top = save_top_table(df, args.metric, 20, out_dir / TOP_CSV_NAME)
    print(f"\nBest 10 runs by {args.metric}:")
    show = ["rank", "Ks_mult", "f_RS_abs", args.metric, "pbias_pct"]
    print(top.loc[:9, show].to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    out_png = out_dir / PNG_NAME
    make_figure(df, args.metric, args.top, args.floor, out_png, args.dpi)

    m = df[args.metric]
    print(f"\nRuns scoring >= {args.floor:g}: {int((m >= args.floor).sum())}   "
          f">= 0.95: {int((m >= 0.95).sum())}   >= 0.99: {int((m >= 0.99).sum())}")
    print(f"\nSaved: {out_png}")
    print(f"Saved: {out_dir / TOP_CSV_NAME}")


if __name__ == "__main__":
    main()
