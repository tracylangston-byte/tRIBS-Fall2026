"""
plot_series120_storms.py
========================
Makes three figures that compare the storm legs of Series 120 (tRIBS 6.0.0,
South Mountain Fan, synthetic truth Ks 7.0 / f 0.012):

    0.8x storm   (rain forcing multiplied by 0.8)
    1.0x storm   (rain as recorded)
    1.25x storm  (rain forcing multiplied by 1.25)

Every leg uses the SAME 400 (Ks, f) draws, and every run is scored against
THAT leg's own synthetic truth. So the leg-to-leg differences you see are caused
by the storm size, not by different draws.

THE THREE FIGURES (all saved as PNG)
    1  fig_series120_storms_1_pbias_maps.png
       One map per storm. Each dot is one run at its Ks multiplier (x) and f
       (y, log scale), coloured by PBIAS (volume error versus that storm's truth):
       blue = too little water, red = too much, near-white = about right.
       Black line = the zero-bias "ridge". Circled dots = runs within 2% of the
       truth volume. Gold star = the truth. A thick black bar on the edge of a
       map means the ridge leaves the sampled Ks box there (it is above Ks 11).

    2  fig_series120_storms_2_feasible_bands.png
       Left: the band of Ks/f pairs that match each storm's truth volume to
       within 2%, one coloured band per storm, and in black the part where ALL
       the storms agree (the "joint" region). Right: how wide each band is in Ks
       at each f. A narrow joint region means that using all the storms together
       pins the parameters down better than any one storm does.

    3  fig_series120_storms_3_pbias_pairs.png
       Top row: for each pair of storms, the run-by-run difference in PBIAS
       (bigger storm minus smaller storm), drawn on the Ks-f map.
       Bottom row: the same pairs as a scatter, PBIAS in the smaller storm (x)
       against PBIAS in the bigger storm (y). A slope below 1 means the same
       Ks and f give smaller percent errors in the bigger storm.

It also writes two small tables you can open in Excel:
    series120_storms_summary.csv   one row per storm
    series120_storms_pairs.csv     one row per pair of storms

WHAT THIS SCRIPT DOES NOT DO
    It does not run tRIBS, and it does not change any file you already have.
    It only READS the results CSVs and WRITES PNGs and two small CSVs.
    It never touches current_run_config.json or any run folder, so it cannot
    interfere with a sweep that is still going. It takes about 15 seconds of
    computer time, so a preview with --legs 080,100 is fine while a sweep runs.

HOW TO RUN (from the lab/ folder, like your other scripts)
    python plot_series120_storms.py

It looks in calibration_work/03_comparisons/summary_tables/ for
    lhs_results_synth_Ks_f_120_storm080.csv     (0.8x)
    lhs_results_synth_Ks_f_120.csv              (1.0x)
    lhs_results_synth_Ks_f_120_storm125.csv     (1.25x)
and saves into
    calibration_work/03_comparisons/sensitivity_plots/series120_storms/

OPTIONS (all optional)
    --legs 080,100,125   which storms to use (default all three). Use
                         --legs 080,100 for a two-storm preview while the 1.25x
                         sweep is still running. At least two are needed.
    --csv_080 PATH  --csv_100 PATH  --csv_125 PATH   use a different file
    --summary_dir PATH   look for the CSVs in a different folder
    --out_dir PATH       save the figures somewhere else
    --feasible_pct 2     the "within this many percent of the truth volume" band
    --dpi 150            picture sharpness (72 to 600)

HOW THE SMOOTH SURFACE IS MADE (so you can explain it)
    The runs are scattered dots; the bands and ridge are drawn from a smooth
    surface fitted through them. For each storm:
      1. For every run take z = log(simulated volume / true volume)
         = log(1 + PBIAS/100). It is 0 when the volume is exactly right.
      2. Fit a "thin-plate spline" to z as a function of Ks and log10(f).
         (A thin-plate spline is like a thin metal sheet bent to pass very
         near all 400 points. A tiny smoothing term, 0.001, stops it from
         chasing rounding noise.)
      3. Read the ridge off as the z = 0 line, and the band as where the
         surface says |PBIAS| is below 2%.
    WHY NOT THE DEGREE-4 POLYNOMIAL used in plot_series120_ridge_old_vs_new.py?
    I tested both by leaving runs out and predicting them. Near the ridge the
    degree-4 polynomial misses by about 2 to 3% in volume, which is as big as
    the 2% band we are trying to draw. The spline misses by a few tenths of a percent.
    The script repeats that leave-out test on YOUR data every time and prints
    the result (look for "fit check" below). If it ever prints a WARNING, the
    bands should not be trusted.
    A fitted surface is still a smooth summary, not a measurement. The runs
    are the evidence. Near the very edge of the box, where there are fewer runs,
    the surface is least certain. Nothing is drawn outside the sampled box.

THINGS TO KEEP IN MIND WHEN READING THEM
    * A band or ridge that reaches the Ks = 11 edge is cut off by the box, not
      by the model. The truth may need a bigger Ks box there.
    * "Joint region" here means: every storm's fitted PBIAS is inside the band
      at the same time. It uses PBIAS (volume) only, with routing parameters
      held at the truth values, so it says how well volume across storms pins
      down Ks and f, not how well the hydrograph shape does.
    * The synthetic truth has no noise, so these regions are the best case.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # draw to a file; no window needed (Codespaces has no screen)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

# DejaVu Sans ships with matplotlib, so the picture looks the same on any machine.
plt.rcParams["font.family"] = "DejaVu Sans"

# ======================================================================
# SETTINGS you might want to change
# ======================================================================
TRUTH_KS = 7.0
TRUTH_F = 0.012
KS_BOX = (3.0, 11.0)           # the sampled Ks range (the "box")
F_BOX = (0.003, 0.05)          # the sampled f range
FEASIBLE_PCT = 2.0             # "within 2% of the truth volume"
COLOR_CLIP_PCT = 60.0          # PBIAS beyond +/- this is drawn in the darkest colour
SMOOTHING = 1e-3               # spline smoothing (see the docstring)
GRID_N = 401                   # the surface is read off a GRID_N x GRID_N grid
CV_FOLDS = 10                  # leave-out test for the fit check
CV_WARN = 0.010                # warn if the leave-out miss near the ridge is above 1% in volume
EXPECTED_DRAWS = 400

# One entry per storm leg, smallest rain first. The colours are slots 1 to 3 of
# the reference categorical palette (blue, orange, aqua); I ran the palette checker
# on them (all-pairs, colour-vision-safe).
LEGS = {
    "080": dict(scale=0.8, name="0.8x storm", note="rain x 0.8", csv="lhs_results_synth_Ks_f_120_storm080.csv",
                color="#2a78d6"),
    "100": dict(scale=1.0, name="1.0x storm", note="rain as recorded", csv="lhs_results_synth_Ks_f_120.csv",
                color="#eb6834"),
    "125": dict(scale=1.25, name="1.25x storm", note="rain x 1.25", csv="lhs_results_synth_Ks_f_120_storm125.csv",
                color="#1baf7a"),
}
OUT_SUBFOLDER = "series120_storms"
FIG1 = "fig_series120_storms_1_pbias_maps.png"
FIG2 = "fig_series120_storms_2_feasible_bands.png"
FIG3 = "fig_series120_storms_3_pbias_pairs.png"
SUMMARY_CSV = "series120_storms_summary.csv"
PAIRS_CSV = "series120_storms_pairs.csv"

# Colours (light-mode chart surface, grey text, one diverging blue-to-red ramp)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SOFT = "#52514e"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
GOLD = "#f0a30a"
CMAP = LinearSegmentedColormap.from_list(
    "pbias_div", ["#184f95", "#5598e7", "#f0efec", "#ee8c7f", "#e34948", "#8f1f1f"]
)
F_TICKS = [0.003, 0.005, 0.01, 0.02, 0.03, 0.05]


class DataProblem(Exception):
    """Something is wrong with an input file. Printed as a plain message, without a traceback."""


# ======================================================================
# Loading and checking one storm's results
# ======================================================================
def load_leg(key, path):
    """Read one results CSV, check it really is this storm, and add the columns the plots need."""
    spec = LEGS[key]
    path = Path(path)
    if not path.exists():
        raise DataProblem(f"Could not find {path}")
    df = pd.read_csv(path)
    need = ["Ks_mult", "f_RS_abs", "pbias_pct", "draw_idx"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise DataProblem(f"{path.name} is missing column(s) {missing}.")
    n_raw = len(df)
    df = df.dropna(subset=need).copy()
    if len(df) < n_raw:
        print(f"  NOTE: {path.name}: dropped {n_raw - len(df)} row(s) with missing values.")
    if df["draw_idx"].duplicated().any():
        raise DataProblem(f"{path.name} has repeated draw_idx values: the same draw appears twice.")
    if "obs_mode" in df.columns and not (df["obs_mode"] == "synth").all():
        raise DataProblem(f"{path.name}: not every row is obs_mode = synth, so this is not a synthetic-truth file.")

    # Is this file really the storm we think it is?
    if "rain_scale" in df.columns:
        vals = df["rain_scale"].astype(float).unique()
        if len(vals) != 1 or abs(vals[0] - spec["scale"]) > 1e-6:
            raise DataProblem(f"{path.name} has rain_scale {[float(v) for v in vals]} but this slot is the {spec['name']} "
                             f"(rain_scale {spec['scale']}). Wrong file for this slot?")
        how = f"rain_scale = {vals[0]:g}"
    elif "storm_leg" in df.columns:
        vals = df["storm_leg"].astype(str).unique()
        if len(vals) != 1 or vals[0] != f"{spec['scale']:.2f}x":
            raise DataProblem(f"{path.name} has storm_leg {list(vals)} but this slot is the {spec['name']}.")
        how = f"storm_leg = {vals[0]}"
    else:
        how = "no rain_scale or storm_leg column, so the storm size could not be confirmed"
    print(f"  {spec['name']:<12} {len(df):3d} runs  ({path.name};  {how})")
    if len(df) < EXPECTED_DRAWS:
        print(f"  WARNING: {spec['name']} has {len(df)} of {EXPECTED_DRAWS} runs. This looks like a sweep that "
              f"is not finished, so every number and figure below is PARTIAL.")
    if len(df) < 100:
        raise DataProblem(f"{path.name} has only {len(df)} runs: too few to fit a surface. Wait for more runs.")

    # log(1 + PBIAS/100) needs PBIAS above -100%; a run that made no water at all would break
    # the fit, so its value is clipped at -99% for the FIT (it is still drawn as a dot).
    df["z"] = np.log1p(np.clip(df["pbias_pct"].values, -99.0, None) / 100.0)
    df["log10_f"] = np.log10(df["f_RS_abs"].values)
    df["draw_idx"] = df["draw_idx"].astype(int)
    return df.sort_values("draw_idx").reset_index(drop=True)


# ======================================================================
# The smooth surface (thin-plate spline, plain numpy)
# ======================================================================
def scale_xy(ks, log10_f):
    """Ks and log10(f) rescaled to be of order 1 (so both directions count equally in the fit)."""
    return np.column_stack([(np.asarray(ks, float) - TRUTH_KS) / 4.0,
                            (np.asarray(log10_f, float) - np.log10(TRUTH_F)) / 0.6])


class ThinPlateSurface:
    """z(x, y) bent to pass very near the given points. Same maths as scipy's
    RBFInterpolator(kernel='thin_plate_spline'); I checked the two agree to 1e-13."""

    def __init__(self, xy, z, smoothing=SMOOTHING):
        self.xy = np.asarray(xy, float)
        n = len(self.xy)
        K = self._phi(self._dist(self.xy, self.xy))
        P = np.column_stack([np.ones(n), self.xy])
        A = np.zeros((n + 3, n + 3))
        A[:n, :n] = K + smoothing * np.eye(n)
        A[:n, n:] = P
        A[n:, :n] = P.T
        sol = np.linalg.solve(A, np.concatenate([np.asarray(z, float), np.zeros(3)]))
        self.w, self.a = sol[:n], sol[n:]

    @staticmethod
    def _dist(a, b):
        return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))

    @staticmethod
    def _phi(r):
        out = np.zeros_like(r)
        m = r > 0
        out[m] = r[m] ** 2 * np.log(r[m])
        return out

    def predict(self, pts, chunk=20000):
        pts = np.asarray(pts, float)
        out = np.empty(len(pts))
        for i in range(0, len(pts), chunk):
            p = pts[i:i + chunk]
            out[i:i + chunk] = self._phi(self._dist(p, self.xy)) @ self.w + self.a[0] + p @ self.a[1:]
        return out


def fit_surface(df):
    return ThinPlateSurface(scale_xy(df["Ks_mult"], df["log10_f"]), df["z"].values)


def fit_check(df):
    """Leave-out test: hide a tenth of the runs, predict them from the rest, repeat.
    Returns the typical miss (RMS, in log volume ratio, about the same as a fraction of the volume)
    for all runs and for the runs near the ridge (|PBIAS| < 10%)."""
    xy = scale_xy(df["Ks_mult"], df["log10_f"])
    z = df["z"].values
    order = np.random.default_rng(0).permutation(len(df))
    pred = np.empty(len(df))
    for fold in np.array_split(order, CV_FOLDS):
        keep = np.setdiff1d(order, fold)
        pred[fold] = ThinPlateSurface(xy[keep], z[keep]).predict(xy[fold])
    err = pred - z
    near = np.abs(df["pbias_pct"].values) < 10.0
    return dict(rms_all=float(np.sqrt(np.mean(err ** 2))),
                rms_near=float(np.sqrt(np.mean(err[near] ** 2))) if near.any() else float("nan"),
                n_near=int(near.sum()))


def make_grid():
    ks = np.linspace(KS_BOX[0], KS_BOX[1], GRID_N)
    lf = np.linspace(np.log10(F_BOX[0]), np.log10(F_BOX[1]), GRID_N)
    KS, LF = np.meshgrid(ks, lf)                      # rows = f, columns = Ks
    return ks, lf, KS, LF


def evaluate_on_grid(surface, KS, LF):
    xy = scale_xy(KS.ravel(), LF.ravel())
    return surface.predict(xy).reshape(KS.shape)


def pbias_of(z):
    return 100.0 * np.expm1(z)


# ======================================================================
# What the fitted surface says: ridge, edge exits, bands
# ======================================================================
def ridge_exits(Z, lf):
    """For each f row: does the fitted surface cross zero somewhere inside the box?
    If not, is the whole row over-predicting (ridge is ABOVE the Ks box) or under-predicting (BELOW)?
    Returns two lists of (f_low, f_high) ranges: above, below."""
    s = np.sign(Z)
    has_cross = np.any(s[:, :-1] * s[:, 1:] <= 0, axis=1)
    above = (~has_cross) & (Z[:, 0] > 0)
    below = (~has_cross) & (Z[:, 0] < 0)

    def ranges(flags):
        out, start = [], None
        for i, flag in enumerate(flags):
            if flag and start is None:
                start = i
            if (not flag or i == len(flags) - 1) and start is not None:
                end = i if flag else i - 1
                out.append((10 ** lf[start], 10 ** lf[end]))
                start = None
        return out
    return ranges(above), ranges(below)


def ridge_note(above, below):
    parts = []
    for lo, hi in above:
        parts.append(f"above the Ks box for f {lo:.3f} to {hi:.3f}")
    for lo, hi in below:
        parts.append(f"below the Ks box for f {lo:.3f} to {hi:.3f}")
    if not parts:
        return "Zero-bias ridge stays inside the Ks box at every f"
    return "Ridge is " + "; ".join(parts)


def band_stats(mask, ks, lf):
    """Width of the feasible band in Ks at each f, and whether it is cut off by the box edge."""
    dks = ks[1] - ks[0]
    width = mask.sum(axis=1) * dks
    cut = mask[:, 0] | mask[:, -1]
    return width, cut


# ======================================================================
# Small drawing helpers
# ======================================================================
def style_axes(ax, xlim=True):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.grid(True, color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_SOFT, labelsize=11)


def style_map_axes(ax, right_pad=0.55):
    style_axes(ax)
    ax.set_xlim(KS_BOX[0] - 0.25, KS_BOX[1] + right_pad)
    ax.set_ylim(F_BOX[0] / 1.07, F_BOX[1] * 1.10)
    ax.set_yscale("log")
    ax.set_yticks(F_TICKS)
    ax.set_yticklabels([f"{t:g}" for t in F_TICKS])
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xticks(range(int(KS_BOX[0]), int(KS_BOX[1]) + 1))
    ax.set_xlabel("Ks multiplier", color=INK_SOFT, fontsize=12)


def draw_box_outline(ax):
    ax.add_patch(Rectangle((KS_BOX[0], F_BOX[0]), KS_BOX[1] - KS_BOX[0], F_BOX[1] - F_BOX[0],
                           fill=False, edgecolor=AXIS, linewidth=1.0, zorder=1))


def draw_truth(ax, size=420):
    ax.scatter([TRUTH_KS], [TRUTH_F], marker="*", s=size, color=GOLD, edgecolors=INK, linewidths=1.2, zorder=8)


def panel_heading(ax, title, sub, sub_lines=2):
    """A bold title and a grey line (or two) under it, a fixed number of points above the panel."""
    ax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, 12 + 14 * sub_lines + 4),
                textcoords="offset points", fontsize=13.5, fontweight="bold", color=INK, va="bottom")
    ax.annotate(sub, xy=(0, 1), xycoords="axes fraction", xytext=(0, 8), textcoords="offset points",
                fontsize=10.5, color=INK_SOFT, va="bottom", linespacing=1.3)


def draw_edge_bars(ax, above, below):
    """A thick black bar along the right (above) or left (below) edge where the ridge leaves the Ks box."""
    for lo, hi in above:
        ax.plot([KS_BOX[1] + 0.2] * 2, [lo, hi], color=INK, linewidth=5, solid_capstyle="butt", zorder=7,
                clip_on=False)
    for lo, hi in below:
        ax.plot([KS_BOX[0] - 0.2] * 2, [lo, hi], color=INK, linewidth=5, solid_capstyle="butt", zorder=7,
                clip_on=False)


def save(fig, out_png, dpi):
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, facecolor=fig.get_facecolor())
    plt.close(fig)


def partial_text(keys, data):
    """A warning line for the top of every figure when a sweep is not finished (empty text if all is complete)."""
    short = [f"{LEGS[k]['name']} has {len(data[k])} of {EXPECTED_DRAWS} runs" for k in keys
             if len(data[k]) < EXPECTED_DRAWS]
    return ("PARTIAL RESULTS: " + "; ".join(short) + ". Treat as preliminary.") if short else ""


def fig_title(fig, title, sub, left, warn=""):
    fig.text(left, 0.962, title, fontsize=17, fontweight="bold", color=INK, ha="left", va="center")
    fig.text(left, 0.925, sub, fontsize=11, color=INK_SOFT, ha="left", va="center")
    if warn:
        fig.text(left, 0.895, warn, fontsize=11.5, fontweight="bold", color=INK, ha="left", va="center")


# ======================================================================
# Figure 1: one PBIAS map per storm
# ======================================================================
def figure_maps(keys, data, info, out_png, dpi):
    n = len(keys)
    fig, axes = plt.subplots(1, n, figsize=(5.7 * n + 1.5, 7.9), sharey=True, squeeze=False)
    axes = axes[0]
    fig.patch.set_facecolor(SURFACE)
    right = 0.89 if n == 3 else 0.865
    fig.subplots_adjust(left=0.07 if n == 3 else 0.075, right=right, top=0.77, bottom=0.19, wspace=0.07)
    norm = Normalize(-COLOR_CLIP_PCT, COLOR_CLIP_PCT)

    for ax, key in zip(axes, keys):
        spec, df, ii = LEGS[key], data[key], info[key]
        style_map_axes(ax)
        draw_box_outline(ax)
        ax.scatter(df["Ks_mult"], df["f_RS_abs"], c=df["pbias_pct"].clip(-COLOR_CLIP_PCT, COLOR_CLIP_PCT),
                   cmap=CMAP, norm=norm, s=36, linewidths=0.5, edgecolors=SURFACE, zorder=3)
        feas = df[df["pbias_pct"].abs() < FEASIBLE_PCT]
        ax.scatter(feas["Ks_mult"], feas["f_RS_abs"], s=150, facecolors="none", edgecolors=INK,
                   linewidths=1.6, zorder=5)
        ax.contour(ii["KS"], 10 ** ii["LF"], ii["Z"], levels=[0.0], colors=INK, linewidths=2.2, zorder=4)
        draw_edge_bars(ax, ii["above"], ii["below"])
        draw_truth(ax)
        panel_heading(ax, f"{spec['name']}  ({spec['note']})",
                      f"{len(feas)} of {len(df)} runs within {FEASIBLE_PCT:g}% of the truth volume\n"
                      f"{ridge_note(ii['above'], ii['below'])}")
    axes[0].set_ylabel("f (1/mm, log scale)", color=INK_SOFT, fontsize=12)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)

    fig_title(fig, "Where the model volume matches each storm's truth: tRIBS 6.0.0, Series 120",
              f"Same draws in every panel. Each dot is scored against its own storm's synthetic truth "
              f"(Ks {TRUTH_KS:g}, f {TRUTH_F:g}).", 0.07 if n == 3 else 0.075, warn=partial_text(keys, data))

    cax = fig.add_axes([right + 0.025, 0.28, 0.013, 0.34])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    cb = fig.colorbar(sm, cax=cax, extend="both", ticks=[-60, -30, 0, 30, 60])
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=INK_SOFT, labelsize=11)
    cb.set_label("PBIAS (%), volume error vs truth\nblue: too little water, red: too much\n(clipped at 60)",
                 color=INK_SOFT, fontsize=10)

    handles = [
        Line2D([], [], marker="*", color="none", markerfacecolor=GOLD, markeredgecolor=INK, markersize=17,
               label=f"truth (Ks {TRUTH_KS:g}, f {TRUTH_F:g})"),
        Line2D([], [], color=INK, linewidth=2.2, label="zero-bias ridge (smooth fit through the runs)"),
        Line2D([], [], marker="o", color="none", markerfacecolor="none", markeredgecolor=INK,
               markeredgewidth=1.6, markersize=12, label=f"run within {FEASIBLE_PCT:g}% of the truth volume"),
        Line2D([], [], color=INK, linewidth=5, solid_capstyle="butt", label="ridge leaves the Ks box here"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=10.5, labelcolor=INK_SOFT,
               bbox_to_anchor=(0.46, 0.01))
    save(fig, out_png, dpi)


# ======================================================================
# Figure 2: feasible bands, joint region, band width
# ======================================================================
def figure_bands(keys, data, info, joint, out_png, dpi):
    ks, lf = joint["ks"], joint["lf"]
    KS, F = joint["KS"], 10 ** joint["LF"]
    n = len(keys)
    fig = plt.figure(figsize=(14.0, 8.9))
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1.0], left=0.07, right=0.975, top=0.77, bottom=0.2, wspace=0.07)
    axM = fig.add_subplot(gs[0, 0])
    axW = fig.add_subplot(gs[0, 1], sharey=axM)

    # ---- map of bands ----
    style_map_axes(axM, right_pad=0.25)
    axM.set_xlim(KS_BOX[0] - 0.25, KS_BOX[1] + 0.25)
    draw_box_outline(axM)
    for key in keys:
        c = LEGS[key]["color"]
        m = info[key]["mask"].astype(float)
        axM.contourf(KS, F, m, levels=[0.5, 1.5], colors=[c], alpha=0.22, zorder=2)
        axM.contour(KS, F, m, levels=[0.5], colors=[c], linewidths=1.8, zorder=3)
    jm = joint["mask"].astype(float)
    if jm.max() > 0:
        axM.contourf(KS, F, jm, levels=[0.5, 1.5], colors=[INK], alpha=0.80, zorder=4)
        axM.contour(KS, F, jm, levels=[0.5], colors=[INK], linewidths=1.8, zorder=5)
    draw_truth(axM, 520)
    axM.set_ylabel("f (1/mm, log scale)", color=INK_SOFT, fontsize=12)

    # ---- band width ----
    style_axes(axW)
    axW.set_yscale("log")
    axW.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axW.tick_params(labelleft=False)
    F1 = 10 ** lf
    for key in keys:
        c = LEGS[key]["color"]
        w, cut = info[key]["width"], info[key]["cut"]
        solid = np.where(cut, np.nan, w)
        near_cut = cut | np.r_[cut[1:], False] | np.r_[False, cut[:-1]]
        dotted = np.where(near_cut, w, np.nan)
        axW.plot(solid, F1, color=c, linewidth=2.0, solid_capstyle="round", zorder=3)
        axW.plot(dotted, F1, color=c, linewidth=2.0, linestyle=(0, (1, 2)), zorder=3)
    jw, jcut = joint["width"], joint["cut"]
    axW.plot(np.where(jcut, np.nan, jw), F1, color=INK, linewidth=2.6, solid_capstyle="round", zorder=4)
    axW.plot(np.where(jcut | np.r_[jcut[1:], False] | np.r_[False, jcut[:-1]], jw, np.nan), F1, color=INK,
             linewidth=2.6, linestyle=(0, (1, 2)), zorder=4)
    axW.axhline(TRUTH_F, color=GOLD, linewidth=1.6, linestyle="--", zorder=2)
    axW.text(0.98, TRUTH_F * 1.04, "truth f", color=INK_SOFT, fontsize=10, ha="right", va="bottom",
             transform=axW.get_yaxis_transform())
    wmax = max(np.nanmax(info[k]["width"]) for k in keys)
    axW.set_xlim(0, max(1.0, wmax * 1.08))
    axW.set_ylim(F_BOX[0] / 1.07, F_BOX[1] * 1.10)
    axW.set_xlabel(f"Width of the {FEASIBLE_PCT:g}% band in Ks units, at each f", color=INK_SOFT, fontsize=12)
    # the two panels share the f axis; set its ticks last so the second panel's setup cannot reset them
    axM.set_yticks(F_TICKS)
    axM.set_yticklabels([f"{t:g}" for t in F_TICKS])
    axM.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axW.tick_params(labelleft=False)

    # ---- headings (all computed) ----
    mean_single = np.mean([info[k]["area"] for k in keys])
    j_area = joint["area"]
    if j_area > 0:
        ratio = mean_single / j_area
        jtxt = (f"Joint region (all {n} storms): {100 * j_area:.2f}% of the box, {ratio:.1f} times smaller than the "
                f"average single storm ({100 * mean_single:.2f}%)\n"
                f"It spans Ks {joint['ks_min']:.2f} to {joint['ks_max']:.2f} and f {joint['f_min']:.4f} to "
                f"{joint['f_max']:.4f}")
    else:
        jtxt = (f"The {n} storms' bands do not overlap anywhere on the fitted surfaces\n"
                f"(average single-storm band: {100 * mean_single:.2f}% of the box)")
    panel_heading(axM, f"Where each storm matches its truth volume to within {FEASIBLE_PCT:g}%", jtxt)
    panel_heading(axW, "How narrow the bands are",
                  "Dotted: the band is cut off by the edge of the box\nBlack: all storms together")

    counts = ", ".join(str(len(data[k])) for k in keys)
    fig_title(fig, "Using several storms together narrows down Ks and f: tRIBS 6.0.0, Series 120",
              f"Bands are read off a smooth surface fitted through each storm's runs ({counts} runs). Black = every "
              f"storm's fitted PBIAS is inside the band.", 0.07, warn=partial_text(keys, data))

    handles = [Patch(facecolor=LEGS[k]["color"], alpha=0.5, edgecolor=LEGS[k]["color"], linewidth=1.8,
                     label=f"{LEGS[k]['name']}: {100 * info[k]['area']:.2f}% of the box") for k in keys]
    handles.append(Patch(facecolor=INK, alpha=0.8, edgecolor=INK, linewidth=1.8,
                         label=f"all {n} storms together: {100 * j_area:.2f}%"))
    handles.append(Line2D([], [], marker="*", color="none", markerfacecolor=GOLD, markeredgecolor=INK,
                          markersize=17, label=f"truth (Ks {TRUTH_KS:g}, f {TRUTH_F:g})"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=11, labelcolor=INK_SOFT,
               bbox_to_anchor=(0.5, 0.025))
    save(fig, out_png, dpi)


# ======================================================================
# Figure 3: pairs of storms
# ======================================================================
def pair_list(keys):
    """(smaller storm, bigger storm) for every pair, ordered 1.0-vs-0.8, 1.25-vs-1.0, 1.25-vs-0.8 when all three exist."""
    pairs = [(keys[i], keys[j]) for i in range(len(keys)) for j in range(i + 1, len(keys))]
    pairs.sort(key=lambda p: (keys.index(p[1]) - keys.index(p[0]), keys.index(p[0])))
    return pairs


def pair_stats(data, small, big):
    """Pair the runs by draw_idx (same Ks and f in both storms) and compare their PBIAS."""
    a = data[small][["draw_idx", "Ks_mult", "f_RS_abs", "pbias_pct"]]
    b = data[big][["draw_idx", "Ks_mult", "f_RS_abs", "pbias_pct"]]
    m = a.merge(b, on="draw_idx", suffixes=("_s", "_b"))
    bad = (np.abs(m["Ks_mult_s"] - m["Ks_mult_b"]) > 1e-9) | (np.abs(m["f_RS_abs_s"] - m["f_RS_abs_b"]) > 1e-12)
    if bad.any():
        raise DataProblem(f"{LEGS[small]['name']} and {LEGS[big]['name']}: {int(bad.sum())} draw(s) have the same "
                         f"draw_idx but different Ks or f. These files do not come from the same draws.")
    xs, yb = m["pbias_pct_s"].values, m["pbias_pct_b"].values
    delta = yb - xs
    slope, icpt = np.polyfit(xs, yb, 1)
    rank = float(np.corrcoef(pd.Series(xs).rank().values, pd.Series(yb).rank().values)[0, 1])
    fs, fb = np.abs(xs) < FEASIBLE_PCT, np.abs(yb) < FEASIBLE_PCT
    near = np.abs(xs) < 10.0
    return dict(merged=m, delta=delta, n_paired=len(m), mean_delta=float(delta.mean()),
                median_delta=float(np.median(delta)), std_delta=float(delta.std(ddof=1)), slope=float(slope),
                intercept=float(icpt), pearson=float(np.corrcoef(xs, yb)[0, 1]), spearman=rank,
                feas_small=int(fs.sum()), feas_big=int(fb.sum()), feas_both=int((fs & fb).sum()),
                near_n=int(near.sum()), near_median_delta=float(np.median(delta[near])) if near.any() else np.nan)


def figure_pairs(keys, data, stats, out_png, dpi):
    pairs = pair_list(keys)
    n = len(pairs)
    fig, axes = plt.subplots(2, n, figsize=(max(12.5, 5.4 * n + 1.4), 11.3), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    right = 0.89 if n > 1 else 0.80
    fig.subplots_adjust(left=0.07, right=right, top=0.81, bottom=0.11, wspace=0.14 if n > 1 else 0.2, hspace=0.46)

    lim = float(np.nanpercentile(np.concatenate([np.abs(stats[p]["delta"]) for p in pairs]), 98))
    lim = max(10.0, float(np.ceil(lim / 10.0) * 10.0))
    dnorm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)

    allp = np.concatenate([np.concatenate([stats[p]["merged"]["pbias_pct_s"], stats[p]["merged"]["pbias_pct_b"]])
                           for p in pairs])
    lo, hi = np.floor(allp.min() / 10.0) * 10.0 - 5.0, np.ceil(allp.max() / 10.0) * 10.0 + 5.0

    for col, (small, big) in enumerate(pairs):
        st = stats[(small, big)]
        m = st["merged"]
        sname, bname = LEGS[small]["name"], LEGS[big]["name"]

        # ---- top: map of the difference ----
        ax = axes[0, col]
        style_map_axes(ax)
        draw_box_outline(ax)
        order = np.argsort(np.abs(st["delta"]))                         # biggest differences drawn last, on top
        sc = ax.scatter(m["Ks_mult_s"].values[order], m["f_RS_abs_s"].values[order],
                        c=np.clip(st["delta"][order], -lim, lim), cmap=CMAP, norm=dnorm, s=36, linewidths=0.5,
                        edgecolors=SURFACE, zorder=3)
        draw_truth(ax, 300)
        panel_heading(ax, f"{bname} minus {sname}",
                      f"PBIAS difference, run by run ({st['n_paired']} paired runs)", sub_lines=1)
        if col == 0:
            ax.set_ylabel("f (1/mm, log scale)", color=INK_SOFT, fontsize=12)
        else:
            ax.tick_params(labelleft=False)

        # ---- bottom: scatter ----
        ax2 = axes[1, col]
        style_axes(ax2)
        ax2.set_xlim(lo, hi)
        ax2.set_ylim(lo, hi)
        ax2.set_aspect("equal", adjustable="box")
        ax2.axhline(0, color=AXIS, linewidth=1.0, zorder=1)
        ax2.axvline(0, color=AXIS, linewidth=1.0, zorder=1)
        ax2.plot([lo, hi], [lo, hi], color=INK_SOFT, linewidth=1.3, linestyle=(0, (4, 3)), zorder=2)
        ax2.scatter(m["pbias_pct_s"], m["pbias_pct_b"], s=26, color=INK, alpha=0.45, linewidths=0.5,
                    edgecolors=SURFACE, zorder=3)
        xx = np.array([lo, hi])
        ax2.plot(xx, st["slope"] * xx + st["intercept"], color=GOLD, linewidth=2.4, zorder=4)
        ax2.set_xlabel(f"PBIAS in the {sname} (%)", color=INK_SOFT, fontsize=11.5)
        ax2.set_ylabel(f"PBIAS in the {bname} (%)", color=INK_SOFT, fontsize=11.5)
        if st["slope"] < 0.95:
            read = "bigger storm: smaller percent errors"
        elif st["slope"] > 1.05:
            read = "bigger storm: larger percent errors"
        else:
            read = "about the same percent errors"
        panel_heading(ax2, f"{bname} against {sname}",
                      f"fitted line: y = {st['slope']:.2f} x {'+' if st['intercept'] >= 0 else '-'} {abs(st['intercept']):.1f}   ({read})\n"
                      f"rank correlation {st['spearman']:.3f}", sub_lines=2)

    cax = fig.add_axes([right + 0.025, 0.58, 0.013, 0.2])
    sm = plt.cm.ScalarMappable(norm=dnorm, cmap=CMAP)
    cb = fig.colorbar(sm, cax=cax, extend="both")
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=INK_SOFT, labelsize=10.5)
    cb.set_label(f"PBIAS difference (percentage points)\nbigger storm minus smaller storm\n(clipped at {lim:g})",
                 color=INK_SOFT, fontsize=10)

    fig_title(fig, "How the volume error changes with storm size: tRIBS 6.0.0, Series 120",
              "The same Ks and f are run against each storm. Top: where the difference sits. "
              "Bottom: the same difference as a scatter.", 0.07, warn=partial_text(keys, data))
    handles = [
        Line2D([], [], color=INK_SOFT, linewidth=1.3, linestyle=(0, (4, 3)), label="same PBIAS in both storms"),
        Line2D([], [], color=GOLD, linewidth=2.4, label="straight-line fit through the runs"),
        Line2D([], [], marker="*", color="none", markerfacecolor=GOLD, markeredgecolor=INK, markersize=15,
               label=f"truth (Ks {TRUTH_KS:g}, f {TRUTH_F:g}), top row"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=10.5, labelcolor=INK_SOFT,
               bbox_to_anchor=(0.46, 0.005))
    save(fig, out_png, dpi)


# ======================================================================
# Main
# ======================================================================
def main():
    global FEASIBLE_PCT                       # --feasible_pct changes the band used everywhere below
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--legs", default="080,100,125", help="comma-separated storms to use (default 080,100,125)")
    for k in LEGS:
        ap.add_argument(f"--csv_{k}", type=Path, default=None)
    ap.add_argument("--summary_dir", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    ap.add_argument("--feasible_pct", type=float, default=FEASIBLE_PCT)
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600.")
    if not (0.1 <= args.feasible_pct <= 20.0):
        ap.error("--feasible_pct must be between 0.1 and 20.")
    FEASIBLE_PCT = args.feasible_pct

    keys = [k.strip() for k in args.legs.split(",") if k.strip()]
    bad = [k for k in keys if k not in LEGS]
    if bad or len(set(keys)) != len(keys):
        ap.error(f"--legs must be a comma-separated list drawn from {list(LEGS)} with no repeats (got {args.legs!r}).")
    if len(keys) < 2:
        ap.error("--legs needs at least two storms.")
    keys = [k for k in LEGS if k in keys]                       # smallest rain first, whatever order was typed

    calib_dir = Path.cwd().parent / "calibration_work"
    summary_dir = args.summary_dir or (calib_dir / "03_comparisons" / "summary_tables")
    out_dir = args.out_dir or (calib_dir / "03_comparisons" / "sensitivity_plots" / OUT_SUBFOLDER)

    # ---- find and load the files ----
    paths = {k: (getattr(args, f"csv_{k}") or summary_dir / LEGS[k]["csv"]) for k in keys}
    missing = [p for p in paths.values() if not Path(p).exists()]
    if missing:
        print("Could not find these results file(s):")
        for p in missing:
            print(f"    {p}")
        print("Run this from the lab/ folder, or point to the file(s) with --csv_080 / --csv_100 / --csv_125 "
              "or --summary_dir.\nIf the 1.25x sweep is still running you can preview with:  "
              "python plot_series120_storms.py --legs 080,100")
        raise SystemExit(1)
    print(f"Loading {len(keys)} storm legs from {summary_dir}:")
    data = {k: load_leg(k, paths[k]) for k in keys}
    draws = [set(data[k]["draw_idx"]) for k in keys]
    common = set.intersection(*draws)
    print(f"  {len(common)} draws are present in every leg.")

    # ---- fit one surface per storm, check it, read it off a grid ----
    ks, lf, KS, LF = make_grid()
    info = {}
    print("\nFit check (leave 10% of the runs out, predict them; typical miss in volume, %):")
    print(f"  {'storm':<12} {'all runs':>9} {'near ridge':>11}   (near ridge = within 10% of zero bias)")
    for k in keys:
        df = data[k]
        surf = fit_surface(df)
        Z = evaluate_on_grid(surf, KS, LF)
        chk = fit_check(df)
        at_truth = float(pbias_of(surf.predict(scale_xy([TRUTH_KS], [np.log10(TRUTH_F)]))[0]))
        mask = np.abs(pbias_of(Z)) < FEASIBLE_PCT
        above, below = ridge_exits(Z, lf)
        width, cut = band_stats(mask, ks, lf)
        info[k] = dict(KS=KS, LF=LF, Z=Z, mask=mask, above=above, below=below, width=width, cut=cut,
                       area=float(mask.mean()), check=chk, at_truth=at_truth)
        flag = ""
        if not np.isfinite(chk["rms_near"]) or chk["rms_near"] > CV_WARN:
            flag = "   WARNING: the fit misses by more than 1% near the ridge, so the bands are not trustworthy"
        print(f"  {LEGS[k]['name']:<12} {100 * chk['rms_all']:8.2f}% {100 * chk['rms_near']:10.2f}%   "
              f"(n near = {chk['n_near']}){flag}")
        print(f"  {'':<12} fitted PBIAS at the truth point (7.0, 0.012), which was not in the fit: {at_truth:+.2f}%"
              + ("   <- should be near 0" if abs(at_truth) > FEASIBLE_PCT else ""))

    # ---- joint region ----
    jmask = np.logical_and.reduce([info[k]["mask"] for k in keys])
    jwidth, jcut = band_stats(jmask, ks, lf)
    joint = dict(ks=ks, lf=lf, KS=KS, LF=LF, mask=jmask, width=jwidth, cut=jcut, area=float(jmask.mean()))
    if jmask.any():
        rows, cols = np.where(jmask)
        joint.update(ks_min=float(ks[cols].min()), ks_max=float(ks[cols].max()),
                     f_min=float(10 ** lf[rows].min()), f_max=float(10 ** lf[rows].max()))

    # ---- pairs ----
    pairs = pair_list(keys)
    pstats = {p: pair_stats(data, *p) for p in pairs}

    # ---- figures ----
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_maps(keys, data, info, out_dir / FIG1, args.dpi)
    figure_bands(keys, data, info, joint, out_dir / FIG2, args.dpi)
    figure_pairs(keys, data, pstats, out_dir / FIG3, args.dpi)

    # ---- tables ----
    rows = []
    for k in keys:
        df, ii = data[k], info[k]
        best = df.loc[df["kge_2012"].idxmax()] if "kge_2012" in df.columns else None
        rows.append({
            "storm": LEGS[k]["name"], "rain_scale": LEGS[k]["scale"], "n_runs": len(df),
            f"n_within_{args.feasible_pct:g}pct": int((df["pbias_pct"].abs() < FEASIBLE_PCT).sum()),
            "band_share_of_box_pct": 100 * ii["area"],
            "ridge_above_box_f_ranges": "; ".join(f"{a:.4f}-{b:.4f}" for a, b in ii["above"]),
            "ridge_below_box_f_ranges": "; ".join(f"{a:.4f}-{b:.4f}" for a, b in ii["below"]),
            "best_kge_2012": float(best["kge_2012"]) if best is not None else np.nan,
            "best_run_Ks": float(best["Ks_mult"]) if best is not None else np.nan,
            "best_run_f": float(best["f_RS_abs"]) if best is not None else np.nan,
            "best_run_pbias_pct": float(best["pbias_pct"]) if best is not None else np.nan,
            "fit_check_rms_all_pct": 100 * ii["check"]["rms_all"],
            "fit_check_rms_near_ridge_pct": 100 * ii["check"]["rms_near"],
            "fitted_pbias_at_truth_pct": ii["at_truth"],
        })
    summary = pd.DataFrame(rows)
    jrow = {"storm": f"all {len(keys)} storms together (joint)", "n_runs": len(common),
            "band_share_of_box_pct": 100 * joint["area"]}
    jrow[f"n_within_{args.feasible_pct:g}pct"] = int(np.all(
        [data[k].set_index("draw_idx").loc[sorted(common), "pbias_pct"].abs().values < FEASIBLE_PCT for k in keys],
        axis=0).sum())
    summary = pd.concat([summary, pd.DataFrame([jrow])], ignore_index=True)
    summary.to_csv(out_dir / SUMMARY_CSV, index=False)

    prow = []
    for (s, b) in pairs:
        st = pstats[(s, b)]
        prow.append({"smaller_storm": LEGS[s]["name"], "bigger_storm": LEGS[b]["name"], "n_paired": st["n_paired"],
                     "mean_diff_pct_pts": st["mean_delta"], "median_diff_pct_pts": st["median_delta"],
                     "std_diff_pct_pts": st["std_delta"], "slope_big_vs_small": st["slope"],
                     "intercept_pct": st["intercept"], "pearson": st["pearson"], "rank_correlation": st["spearman"],
                     "n_within_band_smaller": st["feas_small"], "n_within_band_bigger": st["feas_big"],
                     "n_within_band_both": st["feas_both"],
                     "median_diff_for_runs_near_zero_in_smaller_storm": st["near_median_delta"]})
    pd.DataFrame(prow).to_csv(out_dir / PAIRS_CSV, index=False)

    # ---- printed report ----
    print(f"\nRuns within {FEASIBLE_PCT:g}% of each storm's truth volume (the 'feasible' runs):")
    for k in keys:
        print(f"  {LEGS[k]['name']:<12} {int((data[k]['pbias_pct'].abs() < FEASIBLE_PCT).sum()):3d} of {len(data[k])}")
    print(f"  in every storm at once: {jrow[f'n_within_{args.feasible_pct:g}pct']} of {len(common)} paired draws")
    print("\nShare of the sampled box (Ks x log f) where the fitted PBIAS is inside the band:")
    for k in keys:
        print(f"  {LEGS[k]['name']:<12} {100 * info[k]['area']:6.2f}%    {ridge_note(info[k]['above'], info[k]['below'])}")
    if joint["area"] > 0:
        print(f"  {'all together':<12} {100 * joint['area']:6.2f}%    spans Ks {joint['ks_min']:.2f} to {joint['ks_max']:.2f}, "
              f"f {joint['f_min']:.4f} to {joint['f_max']:.4f}")
    else:
        print(f"  {'all together':<12}   0.00%    the fitted bands do not overlap")
    print("\nPairs of storms (bigger minus smaller storm, same draws):")
    print(f"  {'pair':<28} {'paired':>6} {'mean':>7} {'median':>7} {'std':>6} {'slope':>6} {'rank corr':>9} {'both feasible':>14}")
    for (s, b) in pairs:
        st = pstats[(s, b)]
        print(f"  {LEGS[b]['name'] + ' - ' + LEGS[s]['name']:<28} {st['n_paired']:6d} {st['mean_delta']:+7.2f} "
              f"{st['median_delta']:+7.2f} {st['std_delta']:6.2f} {st['slope']:6.2f} {st['spearman']:9.3f} "
              f"{st['feas_both']:14d}")
    print("\nSaved:")
    for name in (FIG1, FIG2, FIG3, SUMMARY_CSV, PAIRS_CSV):
        print(f"  {out_dir / name}")


if __name__ == "__main__":
    try:
        main()
    except DataProblem as err:
        print(f"\nERROR: {err}")
        raise SystemExit(1)
