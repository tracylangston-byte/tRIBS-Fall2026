#!/usr/bin/env python3
"""
plot_hydrographs_w2_sycamore.py  --  Look at the hydrographs from the 13 Stage 3 runs promoted to W2 (Jan 18 - Mar 12 2017).

A HYDROGRAPH is a plot of streamflow (m3/s) against time at one point, here the watershed outlet. In every panel the BLACK line is
the USGS gauge (09510200), the DASHED GRAY line is the Stage 2 baseline (Ks 0.2, f 1, PsiB 0.5), and the other lines are the 13 runs.

Three figures, each answering ONE question. Read them in this order:

  w2_hydro_1_overview.png   QUESTION: over the whole 53 days, where do the 13 runs miss the gauge?
                            Three panels stacked, same dates:
                              A. streamflow, ordinary (linear) scale: shows the big storms.
                              B. the same flow on a LOG scale: each gridline is 10 times the one below it, so small flows
                                 (0.01 to 10 m3/s) become visible. This is where the recessions ("tails") show up: the gauge
                                 keeps carrying a few m3/s between storms, the model drops toward zero.
                              C. the running total of water past the outlet (cumulative volume, million m3). A steady climb
                                 is water leaving all the time; a staircase is water leaving only during storms.
                            Color = whole-window volume ratio (model / gauge), light = low, dark = high.

  w2_hydro_2_storms.png     QUESTION: how does each storm look up close?
                            Four panels, each with its own vertical scale: the January flood week, the second January storm,
                            the Feb 18-19 storm, the Feb 27-28 storm. Same colors as figure 1.

  w2_hydro_3_anchors.png    QUESTION: what do the proposed W4 anchors look like?
                            The same four storm panels, but only the anchors (default draws 130, 99, 85 and 5), each in its
                            own color. The legend gives each one's Ks, f, PsiB and its W2 volume ratios.
                            Draws 85 and 5 are the pair with the same W2 behavior but very different parameters.

Nothing here changes any model file: it only READS results/<run>_Outlet.qout, the gauge file, w2_summary.csv (and the baseline's
row of score_out_w2/scores_w2.csv, for the legend) and WRITES picture files into a new folder, w2_plots/.

The running total skips any gap longer than 3 hours (in the model or in the gauge file), the same rule the scorer uses, so a gauge
outage does not add made-up water.

Run it from workspaces/Sycamore/lab/ (the folder that holds score_run_sycamore.py and results/):
    python -u plot_hydrographs_w2_sycamore.py --dry-run     # checks every file it needs; writes nothing
    python -u plot_hydrographs_w2_sycamore.py               # makes the three figures

Options:
    --csv FILE         default w2_summary.csv
    --results-dir DIR  default results
    --obs FILE         default ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv
    --out-dir DIR      default w2_plots
    --anchors A,B,C,D  draw numbers (or full run ids) for figure 3; default 130,99,85,5
    --baseline ID      default W2_ks0p2_cv6_psib0p5_nocl (skipped if its file is not found)
    --dpi N            default 150
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                      # no screen needed (Codespace); files only
import matplotlib.dates as mdates          # noqa: E402
import matplotlib.pyplot as plt            # noqa: E402
import numpy as np                         # noqa: E402
import pandas as pd                        # noqa: E402
from matplotlib import cm                  # noqa: E402
from matplotlib.colors import Normalize    # noqa: E402

DEFAULT_OBS = "../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv"
DEFAULT_BASELINE = "W2_ks0p2_cv6_psib0p5_nocl"
DEFAULT_ANCHORS = "130,99,85,5"
W2_START = "2017-01-18 00:00"              # used only if the scorer has no W2 entry in RUN_WINDOWS
W2_END = "2017-03-12 00:00"
MAX_GAP_H = 3.0                            # same rule as the scorer: a gap longer than this adds no volume
LOG_FLOOR = 0.01                           # m3/s; flows below this are drawn at the bottom of the log panel

SPAN = "span__vol_ratio"
F18, F27 = "feb18_19__vol_ratio", "feb27_28__vol_ratio"

# time ranges (Arizona local time, same as the gauge file and the model output)
WHOLE = (W2_START, W2_END)
STORMS = [("A. January flood week (Jan 18-25)", ("2017-01-18 00:00", "2017-01-25 00:00")),
          ("B. Second January storm (Jan 23 12:00 to Jan 25)", ("2017-01-23 12:00", "2017-01-25 00:00")),
          ("C. Feb 18-19 storm (Feb 18 to Feb 21)", ("2017-02-18 00:00", "2017-02-21 00:00")),
          ("D. Feb 27-28 storm (Feb 27 to Mar 3)", ("2017-02-27 00:00", "2017-03-03 00:00"))]

C_GAUGE = "black"
C_BASE = "#666666"
ANCHOR_COLORS = ["#0077BB", "#EE7733", "#009988", "#AA3377"]   # Paul Tol "vibrant", colorblind-safe
ANCHOR_STYLES = ["-", "-", "--", "-."]                         # extra cue so overlapping lines can be told apart
CMAP = plt.get_cmap("Blues")                                    # one hue, light to dark = low to high volume
CMAP_LO, CMAP_HI = 0.30, 0.95                                   # skip the near-white end so the thin lines stay visible


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    print(f"\nERROR: {msg}\n", file=sys.stderr, flush=True)
    sys.exit(2)


# ------------------------------------------------------------------------------------------------------------------
# reading
# ------------------------------------------------------------------------------------------------------------------
def import_scorer():
    """Import Tracy's scorer for its file readers (read_obs, read_qout) so the plots read files exactly as the scorer does."""
    for p in (Path.cwd(), Path(__file__).resolve().parent):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    try:
        import score_run_sycamore as sc
    except ImportError:
        die("cannot import score_run_sycamore.py. Run this script from the lab/ folder (the one that holds score_run_sycamore.py).")
    for name in ("read_obs", "read_qout"):
        if not hasattr(sc, name):
            die(f"score_run_sycamore.py has no '{name}'. Is it the current copy from lab/?")
    return sc


def window_start(sc) -> pd.Timestamp:
    rw = getattr(sc, "RUN_WINDOWS", {}) or {}
    if "W2" in rw:
        return pd.Timestamp(rw["W2"][0])
    return pd.Timestamp(W2_START)


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        die(f"summary CSV not found: {path}. Run this from lab/ (after run_w2_candidates_sycamore.py --summary-only), or give --csv.")
    df = pd.read_csv(path)
    need = ["run_id", "idx", "group", "KS_MULT", "F_MULT", "PSIB_MULT", SPAN, F18, F27]
    miss = [c for c in need if c not in df.columns]
    if miss:
        die(f"{path} is missing columns {miss}")
    if "status" in df.columns:
        df = df[df["status"] == "done"]
    df = df[np.isfinite(df[SPAN])].copy()
    df["draw"] = df["idx"].astype(int)
    return df.reset_index(drop=True)


def find_qout(results_dir: Path, rid: str):
    f = results_dir / f"{rid}_Outlet.qout"
    if f.exists():
        return f
    hits = sorted(results_dir.glob(f"**/{rid}_Outlet.qout"))
    return hits[0] if len(hits) == 1 else None


def load_runs(sc, ids, results_dir: Path):
    """Read each run's outlet flow once. Returns ({run_id: Series}, [run_ids that could not be read])."""
    start = window_start(sc)
    flows, bad = {}, []
    for rid in ids:
        if rid in flows:
            continue
        f = find_qout(results_dir, rid)
        if f is None:
            bad.append(rid)
            continue
        try:
            q = sc.read_qout(f, start)
        except Exception as exc:  # noqa: BLE001
            log(f"  could not read {f}: {exc}")
            bad.append(rid)
            continue
        if q.empty:
            bad.append(rid)
            continue
        flows[rid] = q
    return flows, bad


# ------------------------------------------------------------------------------------------------------------------
# small drawing helpers
# ------------------------------------------------------------------------------------------------------------------
def clip(s: pd.Series, rng) -> pd.Series:
    return s[(s.index >= pd.Timestamp(rng[0])) & (s.index <= pd.Timestamp(rng[1]))]


def cum_volume(s: pd.Series) -> pd.Series:
    """Running total of water past the outlet (million m3), trapezoid rule, skipping gaps longer than MAX_GAP_H and missing values."""
    s = s.dropna()
    if len(s) < 2:
        return pd.Series(dtype=float)
    t = s.index.values.astype("datetime64[s]").astype("int64")
    v = s.to_numpy(float)
    dt = np.diff(t).astype(float)
    step = 0.5 * (v[1:] + v[:-1]) * dt
    step[dt > MAX_GAP_H * 3600] = 0.0
    return pd.Series(np.concatenate([[0.0], np.cumsum(step)]) / 1e6, index=s.index)


def style_time_axis(ax, rng):
    span_h = (pd.Timestamp(rng[1]) - pd.Timestamp(rng[0])) / pd.Timedelta(hours=1)
    if span_h <= 48:
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=4 if span_h <= 30 else 6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:%M"))
    elif span_h <= 200:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    else:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=7))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.set_xlim(pd.Timestamp(rng[0]), pd.Timestamp(rng[1]))
    ax.grid(alpha=0.25)
    ax.tick_params(labelsize=8)


def ylabel_flow(ax):
    ax.set_ylabel("Outlet streamflow (m³/s)", fontsize=9)


def log_clip(s: pd.Series) -> pd.Series:
    return s.clip(lower=LOG_FLOOR)


def color_for(v, norm):
    return CMAP(CMAP_LO + (CMAP_HI - CMAP_LO) * norm(v))


def add_colorbar(fig, norm, rect):
    cax = fig.add_axes(rect)
    sm = cm.ScalarMappable(norm=Normalize(0, 1), cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
        "blues_cut", [CMAP(CMAP_LO), CMAP(CMAP_HI)]))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    ticks = np.linspace(0, 1, 5)
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{norm.vmin + t * (norm.vmax - norm.vmin):.2f}" for t in ticks])
    cb.set_label("whole-window volume ratio\n(model / gauge)", fontsize=9)
    cb.ax.tick_params(labelsize=8)


def legend_handles(have_base: bool):
    h = [plt.Line2D([], [], color=C_GAUGE, lw=2.2, label="USGS gauge")]
    if have_base:
        h.append(plt.Line2D([], [], color=C_BASE, lw=1.6, ls="--", label="Stage 2 baseline (Ks 0.2, f 1, PsiB 0.5)"))
    h.append(plt.Line2D([], [], color=color_for(0.5, Normalize(0, 1)), lw=1.4, label="the 13 runs (color = volume ratio)"))
    return h


# ------------------------------------------------------------------------------------------------------------------
# figure 1: whole window
# ------------------------------------------------------------------------------------------------------------------
def fig_overview(df, flows, gauge, base_q, norm, out, dpi):
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(12.5, 10.4), sharex=True, gridspec_kw=dict(height_ratios=[2.0, 2.0, 1.6]))
    for rid, v in zip(df["run_id"], df[SPAN]):
        col = color_for(v, norm)
        q = clip(flows[rid], WHOLE)
        a1.plot(q.index, q.values, color=col, lw=0.8, alpha=0.85)
        a2.plot(q.index, log_clip(q).values, color=col, lw=0.8, alpha=0.85)
        c = cum_volume(q)
        if len(c):
            a3.plot(c.index, c.values, color=col, lw=0.9, alpha=0.9)
    g = clip(gauge.dropna(), WHOLE)
    cg = cum_volume(g)
    if base_q is not None:
        q = clip(base_q, WHOLE)
        a1.plot(q.index, q.values, color=C_BASE, lw=1.4, ls="--", zorder=8)
        a2.plot(q.index, log_clip(q).values, color=C_BASE, lw=1.4, ls="--", zorder=8)
        cb_ = cum_volume(q)
        a3.plot(cb_.index, cb_.values, color=C_BASE, lw=1.4, ls="--", zorder=8)
    a1.plot(g.index, g.values, color=C_GAUGE, lw=1.6, zorder=10)
    a2.plot(g.index, log_clip(g).values, color=C_GAUGE, lw=1.6, zorder=10)
    a3.plot(cg.index, cg.values, color=C_GAUGE, lw=1.8, zorder=10)
    a2.set_yscale("log")
    a2.set_ylim(LOG_FLOOR, 400)
    titles = [(a1, "A. Streamflow, ordinary scale"),
              (a2, f"B. Same flow on a LOG scale (each gridline is 10 times the one below; flows under {LOG_FLOOR:g} m³/s sit on the bottom edge)"),
              (a3, "C. Running total of water past the outlet")]
    for ax, t in titles:
        style_time_axis(ax, WHOLE)
        ax.set_title(t, fontsize=10, loc="left")
    ylabel_flow(a1)
    ylabel_flow(a2)
    a3.set_ylabel("Cumulative volume (million m³)", fontsize=9)
    a3.set_xlabel("Date (Arizona time)", fontsize=9)
    for ax in (a1, a2, a3):
        ax.tick_params(labelbottom=True)
    a1.legend(handles=legend_handles(base_q is not None), loc="upper right", fontsize=9, frameon=False)
    fig.suptitle("W2 (Jan 18 - Mar 12 2017): the 13 promoted runs against the gauge", fontsize=12)
    fig.tight_layout(rect=(0, 0, 0.9, 0.97))
    add_colorbar(fig, norm, [0.925, 0.30, 0.015, 0.40])
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
# figure 2: four storms, all runs
# ------------------------------------------------------------------------------------------------------------------
def fig_storms(df, flows, gauge, base_q, norm, out, dpi):
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.8))
    for ax, (title, rng) in zip(axes.ravel(), STORMS):
        for rid, v in zip(df["run_id"], df[SPAN]):
            q = clip(flows[rid], rng)
            ax.plot(q.index, q.values, color=color_for(v, norm), lw=0.9, alpha=0.85)
        if base_q is not None:
            q = clip(base_q, rng)
            ax.plot(q.index, q.values, color=C_BASE, lw=1.5, ls="--", zorder=8)
        g = clip(gauge.dropna(), rng)
        ax.plot(g.index, g.values, color=C_GAUGE, lw=2.0, zorder=10)
        style_time_axis(ax, rng)
        ax.set_title(title, fontsize=10, loc="left")
        ylabel_flow(ax)
    axes[0, 0].legend(handles=legend_handles(base_q is not None), loc="upper right", fontsize=8, frameon=False)
    fig.suptitle("Each storm up close (each panel has its own vertical scale)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 0.9, 0.96))
    add_colorbar(fig, norm, [0.925, 0.30, 0.015, 0.40])
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
# figure 3: anchors
# ------------------------------------------------------------------------------------------------------------------
def choose_anchors(df, arg):
    rows = []
    for tok in [t.strip() for t in arg.split(",") if t.strip()]:
        hit = df[df["draw"] == int(tok)] if tok.isdigit() else df[df["run_id"] == tok]
        if hit.empty:
            die(f"--anchors: '{tok}' is not a finished draw number or run id in the CSV")
        rows.append(hit.iloc[0])
    if not 1 <= len(rows) <= len(ANCHOR_COLORS):
        die(f"--anchors: give between 1 and {len(ANCHOR_COLORS)} runs")
    return pd.DataFrame(rows)


def anchor_label(r):
    return (f"draw {int(r['draw'])}: Ks {r['KS_MULT']:g}, f {r['F_MULT']:g}, PsiB {r['PSIB_MULT']:g}   |   "
            f"Feb 18-19 {r[F18]:.2f}, Feb 27-28 {r[F27]:.2f}, whole window {r[SPAN]:.2f}")


def fig_anchors(anch, flows, gauge, base_q, out, dpi):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.6))
    for ax, (title, rng) in zip(axes.ravel(), STORMS):
        g = clip(gauge.dropna(), rng)
        ax.plot(g.index, g.values, color=C_GAUGE, lw=2.2, zorder=10)
        if base_q is not None:
            q = clip(base_q, rng)
            ax.plot(q.index, q.values, color=C_BASE, lw=1.4, ls="--", zorder=5)
        for i, (_, r) in enumerate(anch.iterrows()):
            q = clip(flows[r["run_id"]], rng)
            ax.plot(q.index, q.values, color=ANCHOR_COLORS[i], ls=ANCHOR_STYLES[i], lw=1.7, zorder=6)
        style_time_axis(ax, rng)
        ax.set_title(title, fontsize=10, loc="left")
        ylabel_flow(ax)
    handles = [plt.Line2D([], [], color=C_GAUGE, lw=2.2, label="USGS gauge")]
    if base_q is not None:
        handles.append(plt.Line2D([], [], color=C_BASE, lw=1.4, ls="--", label="Stage 2 baseline (Ks 0.2, f 1, PsiB 0.5)"))
    for i, (_, r) in enumerate(anch.iterrows()):
        handles.append(plt.Line2D([], [], color=ANCHOR_COLORS[i], ls=ANCHOR_STYLES[i], lw=2, label=anchor_label(r)))
    fig.legend(handles=handles, loc="lower center", ncol=1, fontsize=9, frameon=False)
    fig.suptitle("The proposed W4 anchors on the four W2 storms", fontsize=12)
    fig.tight_layout(rect=(0, 0.02 + 0.024 * len(handles), 1, 0.96))
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Hydrograph figures for the 13 runs promoted to W2 (see the top of this file).")
    ap.add_argument("--csv", default="w2_summary.csv")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--obs", default=DEFAULT_OBS)
    ap.add_argument("--out-dir", default="w2_plots")
    ap.add_argument("--anchors", default=DEFAULT_ANCHORS)
    ap.add_argument("--baseline", default=DEFAULT_BASELINE)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    sc = import_scorer()
    df = load_table(Path(a.csv))
    log(f"plot_hydrographs_w2_sycamore.py   {len(df)} finished runs in {a.csv}")
    results = Path(a.results_dir)
    flows, bad = load_runs(sc, list(df["run_id"]) + [a.baseline], results)
    base_q = flows.pop(a.baseline, None)
    if a.baseline in bad:
        bad.remove(a.baseline)
        log(f"  baseline {a.baseline}: outlet file not found, so the baseline line is left out")
    if bad:
        die("cannot read these outlet files (nothing was drawn):\n  " + "\n  ".join(bad))
    if not Path(a.obs).exists():
        die(f"gauge file not found: {a.obs}. Run from lab/, or give --obs.")
    gauge = sc.read_obs(a.obs)
    anch = choose_anchors(df, a.anchors)
    log(f"  read {len(flows)} run files; gauge has {int(gauge.notna().sum())} values; anchors: " + ", ".join(str(int(d)) for d in anch["draw"]))
    log(f"  whole-window volume ratio range: {df[SPAN].min():.2f} to {df[SPAN].max():.2f}")
    if a.dry_run:
        log("\nDRY RUN: every file was found and read; nothing was written.")
        return 0

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    norm = Normalize(vmin=float(df[SPAN].min()), vmax=float(df[SPAN].max()))
    fig_overview(df, flows, gauge, base_q, norm, out / "w2_hydro_1_overview.png", a.dpi)
    fig_storms(df, flows, gauge, base_q, norm, out / "w2_hydro_2_storms.png", a.dpi)
    fig_anchors(anch, flows, gauge, base_q, out / "w2_hydro_3_anchors.png", a.dpi)
    log(f"\nWROTE {out / 'w2_hydro_1_overview.png'}\n      {out / 'w2_hydro_2_storms.png'}\n      {out / 'w2_hydro_3_anchors.png'}")
    log("Nothing in results/ was touched.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
