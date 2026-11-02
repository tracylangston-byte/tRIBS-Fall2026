#!/usr/bin/env python3
"""
plot_hydrographs_lhs_volume_sycamore.py  --  Stage 3: look at the hydrographs from the 150-run Latin hypercube (W1).

A HYDROGRAPH is a plot of streamflow (m3/s) against time at one point, here the watershed outlet. The black line in every
panel is the USGS gauge (09510200). The colored lines are model runs. Everything is the flood week, Jan 18-25 2017.

Three figures, each answering ONE question. Read them in this order:

  hydro_1_volume_groups.png   QUESTION: what does "volume ratio" look like on a hydrograph?
                              Three panels, same axes. Left: runs whose flood-week volume is TOO LOW (ratio under 0.9).
                              Middle: runs INSIDE the band (0.9 to 1.1). Right: runs TOO HIGH (over 1.1).
                              (Volume ratio = model volume / gauge volume. 1.0 = the model carried the same amount of water.)

  hydro_2_inband_detail.png   QUESTION: among the runs that pass the volume test, how different are the hydrographs?
                              Only the in-band runs, four panels: the whole week; a close-up of the main peak; the running
                              total of water (cumulative volume); the second, smaller storm. Color = volume ratio.

  hydro_3_three_picks.png     QUESTION: what do individual runs look like, with their numbers?
                              Three in-band runs (one low, one middle, one high volume) plus the gauge, whole week and
                              close-up. The numbers (Ks, f, PsiB, volume, peak error, arrival error) are in the legend.
                              Also writes hydrograph_picks.csv. Add the Stage 2 baseline run if its files are in results/.

Nothing here changes any model file: it only READS results/<run>_Outlet.qout, the gauge file and lhs_volume_summary.csv, and
WRITES picture files (and one small CSV) into a new folder, lhs_volume_plots/ (the same one the pbias figures used).

Run it from workspaces/Sycamore/lab/ (the folder that holds score_run_sycamore.py and results/):
    python -u plot_hydrographs_lhs_volume_sycamore.py --dry-run     # checks every file it needs; writes nothing
    python -u plot_hydrographs_lhs_volume_sycamore.py               # makes the three figures

Options:
    --csv FILE        default lhs_volume_summary.csv
    --results-dir DIR default results
    --obs FILE        default ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv
    --out-dir DIR     default lhs_volume_plots
    --band N          volume band half-width, default 0.10 (so 0.9 to 1.1)
    --picks A,B,C     choose the three runs for figure 3 yourself, by draw number (for example 105,12,77) or full run id.
                      Default: the in-band runs whose volume is closest to 0.92, 1.00 and 1.08.
    --baseline ID     the Stage 2 baseline run id to add to figure 3 (default W1_ks0p2_cv6_psib0p5_nocl; skipped if not found)
    --dpi N           default 150
"""
from __future__ import annotations

import argparse
import re
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
DEFAULT_BASELINE = "W1_ks0p2_cv6_psib0p5_nocl"
WINDOW = "W1"
VOL, PEAK, ARR = "flood_week__vol_ratio", "flood_week__peak_err_pct", "flood_week__arrival_err_h"

# time ranges for the panels (Arizona local time, same as the gauge file and the model output)
WEEK = ("2017-01-18 00:00", "2017-01-25 00:00")
PEAK_ZOOM = ("2017-01-20 18:00", "2017-01-21 12:00")
SECOND_STORM = ("2017-01-23 12:00", "2017-01-25 00:00")

# colors (colorblind-safe, from Paul Tol's palettes)
C_LOW, C_MID, C_HIGH = "#4477AA", "#228833", "#CC6677"
PICK_COLORS = ["#0077BB", "#EE7733", "#009988"]
C_GAUGE = "black"


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
    for name in ("read_obs", "read_qout", "RUN_WINDOWS", "cumulative"):
        if not hasattr(sc, name):
            die(f"score_run_sycamore.py has no '{name}'. Is it the current copy from lab/?")
    return sc


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        die(f"summary CSV not found: {path}. Run this from lab/, or give --csv.")
    df = pd.read_csv(path)
    need = ["run_id", "KS_MULT", "F_MULT", "PSIB_MULT", VOL, PEAK, ARR]
    miss = [c for c in need if c not in df.columns]
    if miss:
        die(f"{path} is missing columns {miss}")
    if "status" in df.columns:
        df = df[df["status"] == "done"]
    df = df[np.isfinite(df[VOL])].copy()
    df["draw"] = df["run_id"].str.extract(r"_d(\d+)_", expand=False).astype(float)
    return df.reset_index(drop=True)


def load_runs(sc, df: pd.DataFrame, results_dir: Path, extra_ids=()):
    """Read each run's outlet flow once. Returns ({run_id: Series}, [run_ids that could not be read])."""
    start = pd.Timestamp(sc.RUN_WINDOWS[WINDOW][0])
    flows, bad = {}, []
    for rid in list(df["run_id"]) + list(extra_ids):
        if rid in flows:
            continue
        f = results_dir / f"{rid}_Outlet.qout"
        if not f.exists():
            hits = sorted(results_dir.glob(f"**/{rid}_Outlet.qout"))
            f = hits[0] if len(hits) == 1 else None
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


def cum_volume(sc, s: pd.Series, rng) -> pd.Series:
    """Running total of water past the outlet (million m3) from the start of rng: trapezoid rule, as the scorer does."""
    s = clip(s.dropna(), rng)
    if len(s) < 2:
        return pd.Series(dtype=float)
    t = (s.index - s.index[0]).total_seconds().to_numpy(float)
    return pd.Series(sc.cumulative(t, s.to_numpy(float)) / 1e6, index=s.index)


def style_time_axis(ax, rng):
    span_h = (pd.Timestamp(rng[1]) - pd.Timestamp(rng[0])) / pd.Timedelta(hours=1)
    if span_h <= 48:
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=4 if span_h <= 30 else 6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:%M"))
    else:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.set_xlim(pd.Timestamp(rng[0]), pd.Timestamp(rng[1]))
    ax.grid(alpha=0.25)
    ax.tick_params(labelsize=8)


def ylabel_flow(ax):
    ax.set_ylabel("Outlet streamflow (m³/s)", fontsize=9)


def draw_gauge(ax, gauge, rng, lw=2.0, label=None):
    g = clip(gauge.dropna(), rng)
    ax.plot(g.index, g.values, color=C_GAUGE, lw=lw, zorder=10, label=label)


def band_group(v, lo, hi):
    return "low" if v < lo else ("high" if v > hi else "mid")


# ------------------------------------------------------------------------------------------------------------------
# figure 1: three groups by volume
# ------------------------------------------------------------------------------------------------------------------
def fig_groups(df, flows, gauge, lo, hi, out, dpi):
    groups = [("low", C_LOW, f"Volume TOO LOW\n(ratio under {lo:.1f})"),
              ("mid", C_MID, f"Volume IN THE BAND\n({lo:.1f} to {hi:.1f})"),
              ("high", C_HIGH, f"Volume TOO HIGH\n(ratio over {hi:.1f})")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.9), sharey=True)
    for ax, (g, col, title) in zip(axes, groups):
        sub = df[df[VOL].apply(lambda v: band_group(v, lo, hi)) == g]
        for rid in sub["run_id"]:
            q = clip(flows[rid], WEEK)
            ax.plot(q.index, q.values, color=col, lw=0.8, alpha=0.35)
        draw_gauge(ax, gauge, WEEK, label="USGS gauge")
        ax.set_title(f"{title}\n{len(sub)} runs", fontsize=11)
        style_time_axis(ax, WEEK)
        ax.set_xlabel("Date (Arizona time)", fontsize=9)
    ylabel_flow(axes[0])
    axes[0].legend(handles=[plt.Line2D([], [], color=C_GAUGE, lw=2, label="USGS gauge"),
                            plt.Line2D([], [], color=C_MID, lw=1.5, label="model runs")], loc="upper right", fontsize=9, frameon=False)
    fig.suptitle("Flood week (Jan 18-25 2017): the model hydrographs sorted by flood-week volume ratio (model / gauge)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
# figure 2: in-band detail
# ------------------------------------------------------------------------------------------------------------------
def fig_inband(sc, df, flows, gauge, lo, hi, out, dpi):
    sub = df[(df[VOL] >= lo) & (df[VOL] <= hi)].sort_values(VOL)
    norm = Normalize(vmin=lo, vmax=hi)
    cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.6))
    (a1, a2), (a3, a4) = axes
    for rid, v in zip(sub["run_id"], sub[VOL]):
        col = cmap(norm(v))
        for ax, rng in ((a1, WEEK), (a2, PEAK_ZOOM), (a4, SECOND_STORM)):
            q = clip(flows[rid], rng)
            ax.plot(q.index, q.values, color=col, lw=0.9, alpha=0.7)
        c = cum_volume(sc, flows[rid], WEEK)
        if len(c):
            a3.plot(c.index, c.values, color=col, lw=0.9, alpha=0.7)
    for ax, rng in ((a1, WEEK), (a2, PEAK_ZOOM), (a4, SECOND_STORM)):
        draw_gauge(ax, gauge, rng)
    cg = cum_volume(sc, gauge, WEEK)
    a3.plot(cg.index, cg.values, color=C_GAUGE, lw=2.0, zorder=10)
    titles = {a1: "A. The whole week", a2: "B. Close-up of the main peak (Jan 20 18:00 to Jan 21 12:00)",
              a3: "C. Running total of water past the outlet", a4: "D. The second, smaller storm (Jan 23 12:00 to Jan 25)"}
    for ax, rng in ((a1, WEEK), (a2, PEAK_ZOOM), (a3, WEEK), (a4, SECOND_STORM)):
        style_time_axis(ax, rng)
        ax.set_title(titles[ax], fontsize=10, loc="left")
    for ax in (a1, a2, a4):
        ylabel_flow(ax)
    a3.set_ylabel("Cumulative volume (million m³)", fontsize=9)
    a1.legend(handles=[plt.Line2D([], [], color=C_GAUGE, lw=2, label="USGS gauge")], loc="upper right", fontsize=9, frameon=False)
    fig.suptitle(f"The {len(sub)} runs inside the volume band: black = gauge, color = volume ratio (model / gauge)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 0.92, 0.96))
    cax = fig.add_axes([0.935, 0.18, 0.015, 0.62])
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label("volume ratio", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
# figure 3: three picks with numbers
# ------------------------------------------------------------------------------------------------------------------
def choose_picks(df, lo, hi, picks_arg):
    if picks_arg:
        rows = []
        for tok in [t.strip() for t in picks_arg.split(",") if t.strip()]:
            if tok.isdigit():
                hit = df[df["draw"] == int(tok)]
            else:
                hit = df[df["run_id"] == tok]
            if hit.empty:
                die(f"--picks: '{tok}' is not a finished draw number or run id in the CSV")
            rows.append(hit.iloc[0])
        return pd.DataFrame(rows)
    sub = df[(df[VOL] >= lo) & (df[VOL] <= hi)]
    if len(sub) < 3:
        die(f"only {len(sub)} runs are inside the volume band; cannot pick three. Use --picks.")
    chosen = []
    for target in (lo + 0.02, 1.0, hi - 0.02):
        pool = sub[~sub["run_id"].isin([c["run_id"] for c in chosen])]
        chosen.append(pool.iloc[(pool[VOL] - target).abs().argsort().iloc[0]])
    return pd.DataFrame(chosen)


def pick_label(letter, r):
    arr = r[ARR]
    arr_txt = f"{arr:+.1f} h" if np.isfinite(arr) else "n/a"
    return (f"{letter}: Ks {r['KS_MULT']:.3g}, f {r['F_MULT']:.3g}, PsiB {r['PSIB_MULT']:.3g}   |   "
            f"volume {r[VOL]:.2f}, peak {r[PEAK]:+.0f}%, arrival {arr_txt}")


def fig_picks(df_picks, flows, gauge, baseline_q, out, dpi):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.4))
    handles = [plt.Line2D([], [], color=C_GAUGE, lw=2.2, label="USGS gauge")]
    for (ax, rng) in zip(axes, (WEEK, PEAK_ZOOM)):
        draw_gauge(ax, gauge, rng, lw=2.4)
        if baseline_q is not None:
            q = clip(baseline_q, rng)
            ax.plot(q.index, q.values, color="gray", lw=1.6, ls="--", zorder=5)
        for i, (_, r) in enumerate(df_picks.iterrows()):
            q = clip(flows[r["run_id"]], rng)
            ax.plot(q.index, q.values, color=PICK_COLORS[i], lw=1.6, zorder=6)
        style_time_axis(ax, rng)
        ylabel_flow(ax)
    axes[0].set_title("Whole week", fontsize=10, loc="left")
    axes[1].set_title("Close-up of the main peak (Jan 20 18:00 to Jan 21 12:00)", fontsize=10, loc="left")
    for i, (_, r) in enumerate(df_picks.iterrows()):
        handles.append(plt.Line2D([], [], color=PICK_COLORS[i], lw=2, label=pick_label("ABC"[i], r)))
    if baseline_q is not None:
        handles.append(plt.Line2D([], [], color="gray", lw=1.6, ls="--", label="Stage 2 baseline (Ks 0.2, f 1, PsiB 0.5)"))
    fig.legend(handles=handles, loc="lower center", ncol=1, fontsize=9, frameon=False)
    fig.suptitle("Three runs inside the volume band, compared with the gauge", fontsize=12)
    fig.tight_layout(rect=(0, 0.17 if baseline_q is not None else 0.15, 1, 0.96))
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Hydrograph figures for the Stage 3 W1 Latin hypercube (see the top of this file).")
    ap.add_argument("--csv", default="lhs_volume_summary.csv")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--obs", default=DEFAULT_OBS)
    ap.add_argument("--out-dir", default="lhs_volume_plots")
    ap.add_argument("--band", type=float, default=0.10)
    ap.add_argument("--picks", default=None)
    ap.add_argument("--baseline", default=DEFAULT_BASELINE)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    lo, hi = 1 - a.band, 1 + a.band
    sc = import_scorer()
    df = load_table(Path(a.csv))
    log(f"plot_hydrographs_lhs_volume_sycamore.py   {len(df)} finished runs in {a.csv}")
    log(f"  volume band {lo:.2f} to {hi:.2f}:  too low {int((df[VOL] < lo).sum())},  in band {int(((df[VOL] >= lo) & (df[VOL] <= hi)).sum())},"
        f"  too high {int((df[VOL] > hi).sum())}")

    gauge = sc.read_obs(a.obs)
    g_week = clip(gauge.dropna(), WEEK)
    if g_week.empty:
        die(f"the gauge file {a.obs} has no data for {WEEK[0]} to {WEEK[1]}")
    log(f"  gauge: {len(g_week)} values in the flood week, peak {g_week.max():.1f} m3/s at {g_week.idxmax():%b %d %H:%M}")

    results_dir = Path(a.results_dir)
    flows, bad = load_runs(sc, df, results_dir, extra_ids=[a.baseline])
    baseline_q = flows.get(a.baseline)
    if a.baseline in bad:
        bad.remove(a.baseline)
    if baseline_q is None:
        log(f"  (baseline run {a.baseline} not found in {results_dir}; figure 3 will leave it out)")
    else:
        log(f"  baseline run {a.baseline} found")
    log(f"  read {len(df) - len(bad)} of {len(df)} outlet files from {results_dir}")
    if bad:
        log(f"  WARNING: {len(bad)} run(s) have no readable outlet file and will be left out, for example {bad[:3]}")
        df = df[~df["run_id"].isin(bad)].reset_index(drop=True)
    if df.empty:
        die(f"no outlet files could be read from {results_dir}. Run this from lab/, or give --results-dir.")
    if int(((df[VOL] >= lo) & (df[VOL] <= hi)).sum()) == 0:
        die("no run is inside the volume band, so figures 2 and 3 would be empty")

    picks = choose_picks(df, lo, hi, a.picks)
    log("  figure 3 will show:")
    for i, (_, r) in enumerate(picks.iterrows()):
        log(f"    {pick_label('ABC'[i], r)}   [{r['run_id']}]")

    out_dir = Path(a.out_dir)
    files = ["hydro_1_volume_groups.png", "hydro_2_inband_detail.png", "hydro_3_three_picks.png", "hydrograph_picks.csv"]
    if a.dry_run:
        log(f"\nDRY RUN: everything it needs was found. It would write to {out_dir}/: " + ", ".join(files) + ". Nothing was written.")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    fig_groups(df, flows, gauge, lo, hi, out_dir / files[0], a.dpi)
    log(f"  wrote {out_dir / files[0]}")
    fig_inband(sc, df, flows, gauge, lo, hi, out_dir / files[1], a.dpi)
    log(f"  wrote {out_dir / files[1]}")
    fig_picks(picks, flows, gauge, baseline_q, out_dir / files[2], a.dpi)
    log(f"  wrote {out_dir / files[2]}")
    keep = ["run_id", "draw", "KS_MULT", "F_MULT", "PSIB_MULT", VOL, PEAK, ARR]
    out_picks = picks[keep].assign(letter=list("ABC")[:len(picks)])
    out_picks["draw"] = out_picks["draw"].astype("Int64")
    out_picks.to_csv(out_dir / files[3], index=False)
    log(f"  wrote {out_dir / files[3]}")
    log("\nDone. Open the three PNG files in order (1, 2, 3). Nothing in results/ was touched.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
