#!/usr/bin/env python3
"""
regression_effects_lhs_volume_sycamore.py  --  How do Ks, f and PsiB move storm volume? (empirical, from the Latin hypercube runs)

WHAT THIS DOES
  Fits straight lines to the LOGARITHMS of the three soil settings (Ks, f, PsiB multipliers) and the logarithms of the storm
  volume ratios (model volume / gauge volume). In logs, each coefficient is an "elasticity": the percent change in volume for a
  1 percent change in that setting. Example: a Ks coefficient of -0.93 means a 10 percent higher Ks gives about 9 percent less
  volume. Doubling a setting multiplies the volume by 2 ** coefficient.

  Because the Latin hypercube varied the three settings independently of each other, each coefficient is the effect of that one
  setting with the other two held fixed.

  Models fitted (each one is: ln(volume) = constant + a*ln(Ks) + b*ln(f) + c*ln(PsiB)):
    1. Flood-week volume, W1, all 150 runs
    2. Flood-week volume, W2 runs (batches 1 and 2 together)
    3. Feb 18-19 volume, W2 runs
    4. Feb 27-28 volume, W2 runs
    5. Feb 27-28 volume given the flood-week volume (ln flood, ln PsiB, ln f, ln Ks as the inputs): shows which setting changes
       HOW the volume splits between the big storm and the smaller one, rather than the total.

WHAT IT READS (nothing else, and it changes none of them)
  W1 table (150 runs):  lhs_volume_summary.csv   (or lhs_volume_plots/lhs_volume_summary.csv; or give --w1)
  W2 tables:            w2_summary.csv and w2b_summary.csv   (either one alone also works; or give --w2 FILE [FILE ...])
  Needed columns: KS_MULT, F_MULT, PSIB_MULT, flood_week__vol_ratio, and for the W2 files also feb18_19__vol_ratio and
  feb27_28__vol_ratio. Rows with no score (a run that is not finished) are skipped.

WHAT IT WRITES (all in regression_out/, nothing else is touched)
  regression_partial_effects.png   one panel per setting: the volume vs that setting, with the other two held at typical values
  regression_doubling_effects.png  what doubling each setting does to the volume, with 95 percent uncertainty bars
  regression_feb27_vs_flood.png    Feb 27-28 volume vs flood-week volume, colored by PsiB, with the new anchor-rule box
  regression_coefficients.csv      every coefficient, its uncertainty, R2, and a leave-one-out error
  The terminal also prints a plain-language summary.

HOW TO RUN (from workspaces/Sycamore/lab/)
    python -u regression_effects_lhs_volume_sycamore.py
  Options:  --w1 FILE   --w2 FILE [FILE ...]   --out-dir DIR (regression_out)   --dpi 150   --feb18 (add a Feb 18-19 row to
            the first figure; Feb 18-19 is not part of the new anchor rule, so it is off by default)

CAUTIONS (also printed at the end)
  * These are local fits over narrow ranges (Ks varies by a factor of about 2, f about 4, PsiB about 6). Do not trust them
    outside the sampled ranges.
  * Everything else was held fixed: routing at baseline, channel loss off. Volume only, not hydrograph shape or timing.
  * The W2 set is a chosen set of runs (19 or fewer), not a random Latin hypercube, so its numbers are rougher.
  * A regression shows how the MODEL responds. It does not prove why the physics does.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")                      # no screen needed (Codespace); files only
import matplotlib.pyplot as plt            # noqa: E402
import numpy as np                         # noqa: E402
import pandas as pd                        # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm   # noqa: E402
from matplotlib.ticker import FuncFormatter, NullFormatter       # noqa: E402

try:
    from scipy import stats as _stats
except ImportError:                        # scipy is nice to have, not required
    _stats = None

KS, F, PSIB = "KS_MULT", "F_MULT", "PSIB_MULT"
FLOOD, FEB18, FEB27 = "flood_week__vol_ratio", "feb18_19__vol_ratio", "feb27_28__vol_ratio"
SETTINGS = [(KS, "Ks"), (F, "f"), (PSIB, "PsiB")]
FLOOD_BAND = (0.90, 1.12)     # new anchor rule
FEB27_BAND = (0.90, 1.10)     # new anchor rule

# colors (validated reference palette; text in ink colors, never in series colors)
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID, BAND = "#0b0b0b", "#52514e", "#e4e3df", "#eceae4"
BLUE_RAMP = LinearSegmentedColormap.from_list("blue", ["#86b6ef", "#3987e5", "#1c5cab", "#104281"])

W1_DEFAULTS = ["lhs_volume_summary.csv", "lhs_volume_plots/lhs_volume_summary.csv", "lhs_volume_plots/lhs_volume_mainstorm.csv"]
W2_DEFAULTS = ["w2_summary.csv", "w2b_summary.csv"]


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    print(f"\nERROR: {msg}\n", file=sys.stderr, flush=True)
    sys.exit(2)


# --------------------------------------------------------------------------------------------
# Reading the tables
# --------------------------------------------------------------------------------------------
def read_scored(path: Path, need: list) -> pd.DataFrame:
    """Read one summary CSV; keep rows where every needed column is a positive number."""
    t = pd.read_csv(path)
    miss = [c for c in need if c not in t.columns]
    if miss:
        die(f"{path} is missing the columns {miss}")
    if "status" in t.columns and t["status"].notna().any():
        t = t[t["status"].astype(str).str.strip().eq("done")]
    for c in need:
        t[c] = pd.to_numeric(t[c], errors="coerce")
    t = t[np.isfinite(t[need]).all(axis=1) & (t[need] > 0).all(axis=1)].copy()
    label = "idx" if "idx" in t.columns else ("draw" if "draw" in t.columns else None)
    t["draw_no"] = t[label].astype(int) if label else -1
    return t.reset_index(drop=True)


def find_w1(given: str | None) -> Path:
    if given:
        p = Path(given)
        if not p.is_file():
            die(f"--w1 file not found: {p}")
        return p
    for c in W1_DEFAULTS:
        if Path(c).is_file():
            return Path(c)
    die("could not find the W1 table. Looked for: " + ", ".join(W1_DEFAULTS) + ". Run from lab/, or give --w1 FILE.")
    raise SystemExit(2)


# --------------------------------------------------------------------------------------------
# Regression (ordinary least squares, written out so every step is visible)
# --------------------------------------------------------------------------------------------
class Fit:
    """ln(y) = b0 + sum(b_i * x_i). x columns are already logs."""

    def __init__(self, name: str, terms: list, X: np.ndarray, y: np.ndarray):
        self.name, self.terms, self.X, self.y = name, terms, X, y
        n, k = X.shape
        A = np.column_stack([np.ones(n), X])
        self.A = A
        self.b = np.linalg.lstsq(A, y, rcond=None)[0]
        resid = y - A @ self.b
        self.n, self.dof = n, n - (k + 1)
        s2 = float(resid @ resid) / self.dof
        inv = np.linalg.inv(A.T @ A)
        self.se = np.sqrt(np.diag(s2 * inv))
        self.s = float(np.sqrt(s2))
        self.r2 = 1.0 - float(resid @ resid) / float(((y - y.mean()) ** 2).sum())
        hat = np.diag(A @ inv @ A.T)
        self.loo = float(np.sqrt(np.mean((resid / (1.0 - hat)) ** 2)))   # leave-one-out error, in ln units
        self.tcrit = float(_stats.t.ppf(0.975, self.dof)) if _stats else 1.96

    def coef(self, term: str) -> float:
        return float(self.b[1 + self.terms.index(term)])

    def err(self, term: str) -> float:
        return float(self.se[1 + self.terms.index(term)])

    def ci(self, term: str) -> tuple:
        c, e = self.coef(term), self.err(term)
        return c - self.tcrit * e, c + self.tcrit * e

    def r2_without(self, term: str) -> float:
        keep = [i for i, t in enumerate(self.terms) if t != term]
        return Fit("tmp", [self.terms[i] for i in keep], self.X[:, keep], self.y).r2


def fit_volume(name: str, df: pd.DataFrame, ycol: str) -> Fit:
    X = np.log(df[[KS, F, PSIB]].to_numpy(float))
    return Fit(name, ["lnKs", "lnf", "lnPsiB"], X, np.log(df[ycol].to_numpy(float)))


def fit_given_flood(df: pd.DataFrame, with_ks: bool = True) -> Fit:
    cols = [("lnFlood", np.log(df[FLOOD])), ("lnPsiB", np.log(df[PSIB])), ("lnf", np.log(df[F]))]
    if with_ks:
        cols.append(("lnKs", np.log(df[KS])))
    name = "Feb 27-28 volume given flood-week volume" + ("" if with_ks else " (no Ks term)")
    return Fit(name, [c[0] for c in cols], np.column_stack([c[1].to_numpy(float) for c in cols]), np.log(df[FEB27].to_numpy(float)))


# --------------------------------------------------------------------------------------------
# Plain-language printout and the coefficient CSV
# --------------------------------------------------------------------------------------------
def print_fit(m: Fit, df: pd.DataFrame | None = None) -> None:
    log(f"\n{m.name}   (n = {m.n} runs, R2 = {m.r2:.3f}, typical miss = {100 * (np.exp(m.s) - 1):.1f} percent, leave-one-out miss = {100 * (np.exp(m.loo) - 1):.1f} percent)")
    for t in m.terms:
        lo, hi = m.ci(t)
        extra = ""
        if t != "lnFlood":
            extra = f"  doubling it multiplies volume by {2 ** m.coef(t):.2f}"
        else:
            extra = f"  1 percent more flood-week volume goes with {m.coef(t):.2f} percent more Feb 27-28 volume"
        log(f"   {t:7s} {m.coef(t):+.3f}   (95 percent range {lo:+.3f} to {hi:+.3f}){extra}")


def coefficient_rows(m: Fit, ranges: dict | None) -> list:
    rows = []
    for t in m.terms:
        lo, hi = m.ci(t)
        row = dict(model=m.name, n=m.n, term=t, coef=m.coef(t), se=m.err(t), ci95_lo=lo, ci95_hi=hi,
                   doubling_factor=2 ** m.coef(t), r2=m.r2, resid_sd_ln=m.s, leave_one_out_rmse_ln=m.loo)
        if ranges and t in ranges:
            row["sampled_min"], row["sampled_max"] = ranges[t]
            row["pct_change_lowest_to_highest"] = 100 * (np.exp(m.coef(t) * np.log(ranges[t][1] / ranges[t][0])) - 1)
        rows.append(row)
    return rows


# --------------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------------
def style_axes(ax) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, lw=0.7, which="major")
    ax.set_axisbelow(True)


def pretty_ticks(ax, axis: str, lo: float, hi: float, candidates: list) -> None:
    ticks = [c for c in candidates if lo / 1.12 <= c <= hi * 1.12]
    if len(ticks) < 2:
        ticks = list(np.geomspace(lo, hi, 4))
    ax_ = ax.xaxis if axis == "x" else ax.yaxis
    (ax.set_xticks if axis == "x" else ax.set_yticks)(ticks)
    ax_.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax_.set_minor_formatter(NullFormatter())


X_TICKS = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1, 1.5, 2, 3]
Y_TICKS = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 2, 2.5, 3]


def fig_partial(rows: list, out: Path, dpi: int) -> None:
    """rows = [(title, df, ycol, fit, color, band)]. One row of three panels (Ks, f, PsiB) per entry."""
    nr = len(rows)
    fig, axes = plt.subplots(nr, 3, figsize=(11.5, 3.6 * nr + 0.9), sharey="row", squeeze=False)
    for r, (title, df, ycol, m, color, band) in enumerate(rows):
        X = np.log(df[[KS, F, PSIB]].to_numpy(float))
        y = np.log(df[ycol].to_numpy(float))
        mean = X.mean(axis=0)
        lo_y = hi_y = None
        for j, (col, nm) in enumerate(SETTINGS):
            ax = axes[r, j]
            others = [k for k in range(3) if k != j]
            part = y - sum(m.b[1 + k] * (X[:, k] - mean[k]) for k in others)       # volume with the other two settings held typical
            a = m.b[0] + sum(m.b[1 + k] * mean[k] for k in others)
            xs = np.exp(np.linspace(X[:, j].min(), X[:, j].max(), 50))
            ys = np.exp(a + m.b[1 + j] * np.log(xs))
            if band:
                ax.axhspan(band[0], band[1], color=BAND, zorder=0)
            ax.axhline(1.0, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=1)
            ax.scatter(df[col], np.exp(part), s=20 if len(df) > 60 else 38, color=color, edgecolor="white", lw=0.8, alpha=0.9, zorder=3)
            ax.plot(xs, ys, color=INK, lw=2, zorder=4)
            ax.set_xscale("log")
            ax.set_yscale("log")
            style_axes(ax)
            pretty_ticks(ax, "x", float(df[col].min()), float(df[col].max()), X_TICKS)
            c = m.coef(m.terms[j])
            ax.set_title(f"{nm}: doubling it gives volume x{2 ** c:.2f}", fontsize=9.5, loc="left", color=INK)
            ax.set_xlabel(f"{nm} multiplier (log scale)", fontsize=8.5, color=INK2)
            if j == 0:
                ax.set_ylabel(f"{title}\nmodel / gauge volume", fontsize=8.5, color=INK2)
            lo_y = np.exp(part).min() if lo_y is None else min(lo_y, np.exp(part).min())
            hi_y = np.exp(part).max() if hi_y is None else max(hi_y, np.exp(part).max())
        pretty_ticks(axes[r, 0], "y", float(lo_y), float(hi_y), Y_TICKS)
        axes[r, 0].text(0.02, 0.04, f"{m.n} runs, R2 = {m.r2:.2f}", transform=axes[r, 0].transAxes, fontsize=8, color=INK2)
    fig.suptitle("How each setting moves storm volume (the other two held at typical values)", fontsize=12, color=INK, x=0.01, y=0.99, ha="left", va="top")
    fig.text(0.01, 0.005, "Each dot is one run, shifted so the other two settings sit at their typical values. Line = fitted effect. "
                          "Dashed line = gauge volume (ratio 1). Gray band = anchor-rule band. "
                          "Dots that curve away from the line at the edges show where a straight line is too simple.", fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    fig.savefig(out, dpi=dpi, facecolor="white")
    plt.close(fig)


def fig_doubling(models: list, out: Path, dpi: int) -> None:
    """models = [(label, fit, color)]. Grouped bars: percent change in volume when each setting is doubled."""
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    nm = len(models)
    width = 0.34
    for i, (label, m, color) in enumerate(models):
        xs = np.arange(3) + (i - (nm - 1) / 2) * (width + 0.04)
        vals, lo, hi = [], [], []
        for t in m.terms[:3]:
            vals.append(100 * (2 ** m.coef(t) - 1))
            l, h = m.ci(t)
            lo.append(100 * (2 ** l - 1))
            hi.append(100 * (2 ** h - 1))
        vals = np.array(vals)
        ax.bar(xs, vals, width=width, color=color, label=label, zorder=3)
        ax.errorbar(xs, vals, yerr=[vals - np.array(lo), np.array(hi) - vals], fmt="none", ecolor=INK, elinewidth=1.2, capsize=3, zorder=4)
        for x, v, h_, l_ in zip(xs, vals, hi, lo):
            ax.text(x, (h_ if v >= 0 else l_) + (3 if v >= 0 else -3), f"{v:+.0f}%", ha="center", va="bottom" if v >= 0 else "top", fontsize=9, color=INK)
    ax.axhline(0, color=INK2, lw=1, zorder=2)
    all_lo = [100 * (2 ** m.ci(t)[0] - 1) for _, m, _ in models for t in m.terms[:3]]
    all_hi = [100 * (2 ** m.ci(t)[1] - 1) for _, m, _ in models for t in m.terms[:3]]
    ax.set_ylim(min(all_lo) - 14, max(all_hi) + 14)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Ks", "f", "PsiB"], fontsize=10, color=INK)
    ax.set_ylabel("change in storm volume when the setting is doubled (%)", fontsize=9, color=INK2)
    style_axes(ax)
    ax.grid(False, axis="x")
    ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK)
    ax.set_title("What doubling each setting does to storm volume", fontsize=12, loc="left", color=INK)
    fig.text(0.01, 0.01, "Bars = fitted effect. Black whiskers = 95 percent range. The smaller storm is more sensitive to every setting.",
             fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out, dpi=dpi, facecolor="white")
    plt.close(fig)


def fig_feb27_vs_flood(df: pd.DataFrame, m3: Fit, out: Path, dpi: int) -> None:
    """Feb 27-28 vs flood week, colored by PsiB, with the new anchor-rule box and fitted lines at three PsiB values."""
    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    norm = LogNorm(vmin=float(df[PSIB].min()), vmax=float(df[PSIB].max()))
    ax.add_patch(plt.Rectangle((FLOOD_BAND[0], FEB27_BAND[0]), FLOOD_BAND[1] - FLOOD_BAND[0], FEB27_BAND[1] - FEB27_BAND[0],
                               facecolor=BAND, edgecolor=INK2, lw=1, ls=(0, (4, 3)), zorder=1))
    ax.text(FLOOD_BAND[0] + 0.004, FEB27_BAND[1] - 0.012, "anchor-rule box", fontsize=8.5, color=INK2, va="top")
    f_med = float(df[F].median())
    xs = np.linspace(float(df[FLOOD].min()), float(df[FLOOD].max()), 40)
    for p in sorted({round(float(v), 2) for v in (df[PSIB].min(), np.sqrt(df[PSIB].min() * df[PSIB].max()), df[PSIB].max())}):
        ys = np.exp(m3.b[0] + m3.coef("lnFlood") * np.log(xs) + m3.coef("lnPsiB") * np.log(p) + m3.coef("lnf") * np.log(f_med))
        ax.plot(xs, ys, color=BLUE_RAMP(norm(p)), lw=1.8, zorder=2)
        ax.text(xs[-1] + 0.04, ys[-1], f"PsiB {p:g}", fontsize=8.5, color=INK2, va="center")
        ax.plot([xs[-1], xs[-1] + 0.035], [ys[-1], ys[-1]], color=BLUE_RAMP(norm(p)), lw=0.8, ls=":", zorder=2)
    sc = ax.scatter(df[FLOOD], df[FEB27], c=df[PSIB], cmap=BLUE_RAMP, norm=norm, s=62, edgecolor=INK, lw=0.9, zorder=3)
    for _, r in df.iterrows():
        if r["draw_no"] >= 0:
            ax.annotate(str(int(r["draw_no"])), (r[FLOOD], r[FEB27]), textcoords="offset points", xytext=(-4, 7), ha="right", fontsize=7.5, color=INK2, zorder=4)
    cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.045)
    cb.set_label("PsiB multiplier (log scale)", fontsize=9, color=INK2)
    cb.ax.tick_params(labelsize=8, colors=INK2)
    ticks = [t for t in (0.3, 0.5, 0.7, 1.0, 1.4) if norm.vmin <= t <= norm.vmax]
    cb.set_ticks(ticks)
    cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    cb.ax.yaxis.set_minor_formatter(NullFormatter())
    style_axes(ax)
    ax.set_xlim(0.70, 1.30)
    ax.set_ylim(0.20, 1.18)
    ax.set_xlabel("flood-week volume  (model / gauge)", fontsize=9.5, color=INK2)
    ax.set_ylabel("Feb 27-28 volume  (model / gauge)", fontsize=9.5, color=INK2)
    ax.set_title("Feb 27-28 rises faster than the flood week, and lower PsiB lifts it", fontsize=11.5, loc="left", color=INK)
    fig.text(0.01, 0.01, "Dots = W2 runs (numbers are draw numbers). Lines = fitted Feb 27-28 volume at three PsiB values, f at its median. "
                         "Lines are drawn only over the range of the data.", fontsize=7.5, color=INK2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out, dpi=dpi, facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Regression of storm volume on Ks, f and PsiB (see the top of this file).")
    ap.add_argument("--w1", default=None, help="the 150-run W1 table")
    ap.add_argument("--w2", nargs="+", default=None, help="one or more W2 summary tables (default w2_summary.csv w2b_summary.csv)")
    ap.add_argument("--out-dir", default="regression_out")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--feb18", action="store_true", help="add a Feb 18-19 row to the partial-effects figure")
    a = ap.parse_args()

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    log("regression_effects_lhs_volume_sycamore.py")
    if not _stats:
        log("  (scipy not found: 95 percent ranges use 1.96, which is slightly too narrow for the small W2 sets)")

    # ---- W1
    w1p = find_w1(a.w1)
    w1 = read_scored(w1p, [KS, F, PSIB, FLOOD])
    log(f"  W1: {len(w1)} runs from {w1p}")
    rng = {"lnKs": (w1[KS].min(), w1[KS].max()), "lnf": (w1[F].min(), w1[F].max()), "lnPsiB": (w1[PSIB].min(), w1[PSIB].max())}
    log(f"  sampled ranges (multipliers): Ks {rng['lnKs'][0]:.3f} to {rng['lnKs'][1]:.3f}, f {rng['lnf'][0]:.3f} to {rng['lnf'][1]:.3f}, "
        f"PsiB {rng['lnPsiB'][0]:.3f} to {rng['lnPsiB'][1]:.3f}")
    xin = np.log(w1[[KS, F, PSIB]].to_numpy(float))
    log("  correlations between the log settings (near 0 = the design varied them independently): "
        f"Ks-f {np.corrcoef(xin[:, 0], xin[:, 1])[0, 1]:+.2f}, Ks-PsiB {np.corrcoef(xin[:, 0], xin[:, 2])[0, 1]:+.2f}, "
        f"f-PsiB {np.corrcoef(xin[:, 1], xin[:, 2])[0, 1]:+.2f}")
    m1 = fit_volume("Flood-week volume, W1", w1, FLOOD)
    print_fit(m1)
    log("   How much of the run-to-run variation each setting explains (R2 lost if it is removed): "
        + ", ".join(f"{t} {m1.r2 - m1.r2_without(t):.2f}" for t in m1.terms))
    d = m1.coef("lnKs") * np.log(1.10)
    if m1.coef("lnf") > 0 and m1.coef("lnPsiB") < 0:
        log(f"   Trade rates at constant flood-week volume: a 10 percent higher Ks is undone by f about {100 * (np.exp(-d / m1.coef('lnf')) - 1):.0f} "
            f"percent higher, or PsiB about {100 * (1 - np.exp(-d / m1.coef('lnPsiB'))):.0f} percent lower.")
    all_rows = coefficient_rows(m1, rng)

    # ---- W2
    w2_paths = a.w2 if a.w2 else W2_DEFAULTS
    parts = []
    for p in w2_paths:
        if Path(p).is_file():
            t = read_scored(Path(p), [KS, F, PSIB, FLOOD, FEB18, FEB27])
            log(f"  W2: {len(t)} scored runs from {p}")
            parts.append(t)
        else:
            log(f"  W2: {p} not found (skipped)")
    figures = []
    models = []
    if parts:
        w2 = pd.concat(parts, ignore_index=True)
        key = "run_id" if "run_id" in w2.columns else "draw_no"
        w2 = w2.drop_duplicates(subset=key).reset_index(drop=True)
        log(f"  W2: {len(w2)} runs in total after combining")
        if len(w2) < 8:
            log("  Fewer than 8 W2 runs: skipping the W2 models (too few for a 3-setting fit).")
        else:
            rng2 = {"lnKs": (w2[KS].min(), w2[KS].max()), "lnf": (w2[F].min(), w2[F].max()), "lnPsiB": (w2[PSIB].min(), w2[PSIB].max())}
            m2 = fit_volume("Flood-week volume, W2", w2, FLOOD)
            m_f18 = fit_volume("Feb 18-19 volume, W2", w2, FEB18)
            m_f27 = fit_volume("Feb 27-28 volume, W2", w2, FEB27)
            m3 = fit_given_flood(w2, with_ks=True)
            m3b = fit_given_flood(w2, with_ks=False)
            for m in (m2, m_f18, m_f27, m3, m3b):
                print_fit(m)
                all_rows += coefficient_rows(m, rng2 if m.terms[0] == "lnKs" else None)
            klo, khi = m3.ci("lnKs")
            log("\n   Reading the last two models: once the flood-week volume is known, Ks shows no clear extra effect on Feb 27-28 "
                f"(coefficient {m3.coef('lnKs'):+.2f}, 95 percent range {klo:+.2f} to {khi:+.2f}, which includes zero), "
                "so Ks looks like the TOTAL-volume dial: it reaches Feb 27-28 through the flood-week volume. "
                f"PsiB ({m3b.coef('lnPsiB'):+.2f}) and f ({m3b.coef('lnf'):+.2f}) change HOW the water splits between the storms. "
                "With only 19 runs, 'no clear effect' is not the same as 'proven zero'.")
            if min(m3b.coef("lnPsiB"), 0) < 0:
                log(f"   At a fixed flood-week volume, halving PsiB changes Feb 27-28 volume by {100 * (0.5 ** m3b.coef('lnPsiB') - 1):+.0f} percent "
                    f"and halving f by {100 * (0.5 ** m3b.coef('lnf') - 1):+.0f} percent. PsiB is sampled down to {w2[PSIB].min():.2f} only, "
                    "so anything lower is a guess that needs real runs.")

            prow = [("Flood week (W1)", w1, FLOOD, m1, BLUE, FLOOD_BAND),
                    (f"Feb 27-28 (W2)", w2, FEB27, m_f27, ORANGE, FEB27_BAND)]
            if a.feb18:
                prow.append(("Feb 18-19 (W2)", w2, FEB18, m_f18, ORANGE, None))
            fig_partial(prow, out / "regression_partial_effects.png", a.dpi)
            fig_doubling([("Flood week (W1, %d runs)" % m1.n, m1, BLUE), ("Feb 27-28 (W2, %d runs)" % m_f27.n, m_f27, ORANGE)],
                         out / "regression_doubling_effects.png", a.dpi)
            fig_feb27_vs_flood(w2, m3b, out / "regression_feb27_vs_flood.png", a.dpi)
            figures = ["regression_partial_effects.png", "regression_doubling_effects.png", "regression_feb27_vs_flood.png"]
    else:
        log("  No W2 table found: only the W1 model was fitted, and only the first figure's flood-week row is drawn.")
        fig_partial([("Flood week (W1)", w1, FLOOD, m1, BLUE, FLOOD_BAND)], out / "regression_partial_effects.png", a.dpi)
        figures = ["regression_partial_effects.png"]

    pd.DataFrame(all_rows).to_csv(out / "regression_coefficients.csv", index=False)
    log("\nWROTE in " + str(out) + "/ : " + ", ".join(figures + ["regression_coefficients.csv"]))
    log("\nCAUTIONS: local fits over narrow ranges (do not trust them outside the sampled ranges); everything else held fixed "
        "(routing at baseline, channel loss off); volume only, not shape or timing; the W2 set is chosen, not random, so its "
        "numbers are rougher; a regression shows how the model responds, not why the physics does.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
