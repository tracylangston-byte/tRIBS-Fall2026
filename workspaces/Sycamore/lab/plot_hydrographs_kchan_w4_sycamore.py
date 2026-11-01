#!/usr/bin/env python3
"""
plot_hydrographs_kchan_w4_sycamore.py  --  Look at the channel-loss (Kchan) check on W4 (Oct 1 2016 - Apr 1 2017).

A HYDROGRAPH is a plot of streamflow (m3/s) against time at one point, here the watershed outlet.
In every figure the BLACK line is the USGS gauge (09510200), the DASHED GRAY line is draw 99 with channel loss OFF, and the
BLUE lines are the same soil set with channel loss ON at Kchan 0.5, 1 and 2 mm/hr (light blue = lowest Kchan, dark blue = highest).

Four figures, each answering ONE question. Read them in this order:

  kchan_hydro_1_fall.png     QUESTION: which fall (Oct-Dec) flows does channel loss remove, and which stay?
                             Three panels stacked, Oct 1 to Dec 31 2016:
                               A. streamflow, ordinary (linear) scale. The dotted line is 0.13 m3/s, the false-flow threshold
                                  (it equals the gauge's own largest Oct-Dec flow). Anything the model draws above it, where the
                                  gauge is flat, is false flow.
                               B. the same flow on a LOG scale: each gridline is 10 times the one below it, so small flows
                                  (0.01 to 10 m3/s) become visible. Flow under 0.01 m3/s, including exactly zero, sits on the
                                  0.01 line at the bottom. A line riding that bottom line is a dry channel.
                               C. the running total of water past the outlet (cumulative volume, million m3).

  kchan_hydro_2_tails.png    QUESTION: what happens to the recession "tails"?
                             A TAIL is the slow decline in flow after a storm. The gauge keeps carrying flow for days after the
                             January flood and again in early March. Two columns: the tail after the January flood (Jan 26 - Feb 2)
                             and the March tail (Mar 3 - 11). Top row: log-scale flow. Bottom row: the running total of water
                             during that window only (it starts at zero at the left edge of each panel).

  kchan_hydro_3_storms.png   QUESTION: what does channel loss cost the storms themselves?
                             The four storm panels used for W2 (January flood week, second January storm, Feb 18-19, Feb 27-28),
                             ordinary scale, each panel with its own vertical scale.

  kchan_hydro_4_overview.png QUESTION: over the whole six months, where does the water go?
                             A. flow on a log scale for the whole window; B. the running total for the whole window.

The legends give the matching numbers from kchan_w4_summary.csv (volume ratios are model / gauge, so 1.00 is a perfect match).

Nothing here changes any model file: it only READS results/<run>_Outlet.qout, the gauge file and kchan_w4_summary.csv, and WRITES
picture files into a new folder, kchan_w4_plots/.

The running total skips any gap longer than 3 hours (in the model or in the gauge file), the same rule the scorer uses, so a gauge
outage does not add made-up water.

Run it from workspaces/Sycamore/lab/ (the folder that holds score_run_sycamore.py and results/):
    python -u plot_hydrographs_kchan_w4_sycamore.py --dry-run     # checks every file it needs; writes nothing
    python -u plot_hydrographs_kchan_w4_sycamore.py               # makes the four figures

If you added Kchan values later (5, 10, or lower ones), run the batch's --summary-only first so kchan_w4_summary.csv is current,
then run this again: the new runs are drawn as extra blue lines (up to 8 runs).

Options:
    --csv FILE         default kchan_w4_summary.csv
    --results-dir DIR  default results
    --obs FILE         default ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv
    --out-dir DIR      default kchan_w4_plots
    --dpi N            default 150
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                      # no screen needed (Codespace); files only
import matplotlib.dates as mdates          # noqa: E402
import matplotlib.pyplot as plt            # noqa: E402
import numpy as np                         # noqa: E402
import pandas as pd                        # noqa: E402

DEFAULT_OBS = "../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv"
W4_START = "2016-10-01 00:00"              # used only if the scorer has no W4 entry in RUN_WINDOWS
W4_END = "2017-04-01 00:00"
MAX_GAP_H = 3.0                            # same rule as the scorer: a gap longer than this adds no volume
LOG_FLOOR = 0.01                           # m3/s; flows below this are drawn ON the 0.01 line of the log panels
LOG_YMIN = LOG_FLOOR / 3                   # the log panels reach a little below that line so a dry channel shows as a line
FF_THRESHOLD = 0.13                        # m3/s; the false-flow threshold (the gauge's own Oct-Dec maximum)
MAX_RUNS = 8

# time ranges (Arizona local time, same as the gauge file and the model output)
WHOLE = (W4_START, W4_END)
FALL = ("2016-10-01 00:00", "2017-01-01 00:00")
TAILS = [("A. Tail after the January flood (Jan 26 - Feb 2)", ("2017-01-26 00:00", "2017-02-03 00:00")),
         ("B. March tail (Mar 3 - 11)", ("2017-03-03 00:00", "2017-03-12 00:00"))]
STORMS = [("A. January flood week (Jan 18-25)", ("2017-01-18 00:00", "2017-01-25 00:00")),
          ("B. Second January storm (Jan 23 12:00 to Jan 25)", ("2017-01-23 12:00", "2017-01-25 00:00")),
          ("C. Feb 18-19 storm (Feb 18 to Feb 21)", ("2017-02-18 00:00", "2017-02-21 00:00")),
          ("D. Feb 27-28 storm (Feb 27 to Mar 3)", ("2017-02-27 00:00", "2017-03-03 00:00"))]

C_GAUGE = "black"
C_OFF = "#666666"
CMAP = plt.get_cmap("Blues")               # one hue, light to dark = low to high Kchan

# score columns used in the legends (all are columns of kchan_w4_summary.csv; missing ones just show as "-")
C_FFV, C_FFDAYS = "octdec__vol_model_Mm3", "octdec__days_above_model"
C_FFGV, C_FFGMAX = "octdec__vol_gauge_Mm3", "octdec__max_gauge_m3s"
C_VOL, C_J2 = "flood_week__vol_ratio", "jan_second_storm__vol_ratio"
C_F18, C_F27 = "feb18_19__vol_ratio", "feb27_28__vol_ratio"
C_T1, C_T2 = "tail_jan26_feb2__tail_vol_ratio", "tail_mar3_mar11__tail_vol_ratio"
C_OM, C_DRY = "oct_mar__vol_ratio", "dry_pct"


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
    if "W4" in rw:
        return pd.Timestamp(rw["W4"][0])
    return pd.Timestamp(W4_START)


def num(row, col) -> float:
    try:
        v = float(row.get(col, float("nan")))
    except (TypeError, ValueError):
        return float("nan")
    return v if math.isfinite(v) else float("nan")


def fmt(x, nd=2) -> str:
    return "-" if not math.isfinite(x) else f"{x:.{nd}f}"


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        die(f"summary CSV not found: {path}. Run this from lab/ (after run_kchan_w4_sycamore.py --summary-only), or give --csv.")
    df = pd.read_csv(path)
    miss = [c for c in ("run_id", "KCHAN_MMHR") if c not in df.columns]
    if miss:
        die(f"{path} is missing columns {miss}")
    df = df[df["run_id"].notna()].copy()
    df["is_ref"] = df["KCHAN_MMHR"].astype(str).str.strip().str.lower().eq("off")
    if "status" in df.columns:
        df = df[df["status"].isin(["done", "reference"])]
    df["K"] = pd.to_numeric(df["KCHAN_MMHR"], errors="coerce")
    df.loc[df["is_ref"], "K"] = 0.0
    df = df[np.isfinite(df["K"])].sort_values("K").reset_index(drop=True)
    if int((~df["is_ref"]).sum()) < 1:
        die(f"{path} has no finished Kchan run yet.")
    if int((~df["is_ref"]).sum()) > MAX_RUNS:
        die(f"{path} has more than {MAX_RUNS} Kchan runs; this script draws at most {MAX_RUNS}.")
    return df


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
    elif span_h <= 250:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    elif span_h <= 1000:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=7))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    elif span_h <= 2500:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    else:
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.set_xlim(pd.Timestamp(rng[0]), pd.Timestamp(rng[1]))
    ax.grid(alpha=0.25)
    ax.tick_params(labelsize=8)


def ylabel_flow(ax):
    ax.set_ylabel("Outlet streamflow (m³/s)", fontsize=9)


def log_clip(s: pd.Series) -> pd.Series:
    return s.clip(lower=LOG_FLOOR)


def log_axis(ax, ymax):
    ax.set_yscale("log")
    ax.set_ylim(LOG_YMIN, ymax)


def styles(df):
    """One (color, linestyle, linewidth) per row. Loss off is gray dashed; Kchan runs run light to dark blue in order of Kchan."""
    n = int((~df["is_ref"]).sum())
    levels = np.linspace(0.40, 0.97, n) if n > 1 else [0.75]
    out, j = {}, 0
    for rid, is_ref in zip(df["run_id"], df["is_ref"]):
        if is_ref:
            out[rid] = (C_OFF, "--", 1.7)
        else:
            out[rid] = (CMAP(levels[j]), "-", 1.6)
            j += 1
    return out


def name_of(r) -> str:
    return "loss OFF (same soil set)" if r["is_ref"] else f"Kchan {r['K']:g} mm/hr"


def legend_handles(df, sty, labeler, gauge_label):
    h = [plt.Line2D([], [], color=C_GAUGE, lw=2.2, label=gauge_label)]
    for _, r in df.iterrows():
        c, ls, lw = sty[r["run_id"]]
        h.append(plt.Line2D([], [], color=c, ls=ls, lw=lw + 0.4, label=f"{name_of(r)}:  {labeler(r)}"))
    return h


def bottom_legend(fig, handles, top_rect):
    """Legend under the panels; the panels are packed into the space above it."""
    fig.legend(handles=handles, loc="lower center", ncol=1, fontsize=9, frameon=False)
    fig.tight_layout(rect=(0, 0.02 + 0.026 * len(handles), 1, top_rect))


def draw_lines(ax, df, flows, sty, rng, transform=None, zorder0=3):
    """Draw each run in order (off first, highest Kchan last); the gauge is added by the caller on top."""
    for i, (_, r) in enumerate(df.iterrows()):
        q = clip(flows[r["run_id"]], rng)
        if transform is not None:
            q = transform(q)
        if len(q) == 0:
            continue
        c, ls, lw = sty[r["run_id"]]
        ax.plot(q.index, q.values, color=c, ls=ls, lw=lw, zorder=zorder0 + i)


# ------------------------------------------------------------------------------------------------------------------
# figure 1: fall
# ------------------------------------------------------------------------------------------------------------------
def fig_fall(df, flows, gauge, sty, gauge_info, out, dpi):
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(12.5, 10.8), sharex=True, gridspec_kw=dict(height_ratios=[2.0, 2.0, 1.5]))
    draw_lines(a1, df, flows, sty, FALL)
    draw_lines(a2, df, flows, sty, FALL, transform=log_clip)
    for i, (_, r) in enumerate(df.iterrows()):
        c = cum_volume(clip(flows[r["run_id"]], FALL))
        if len(c):
            col, ls, lw = sty[r["run_id"]]
            a3.plot(c.index, c.values, color=col, ls=ls, lw=lw, zorder=3 + i)
    g = clip(gauge.dropna(), FALL)
    a1.plot(g.index, g.values, color=C_GAUGE, lw=2.0, zorder=20)
    a2.plot(g.index, log_clip(g).values, color=C_GAUGE, lw=2.0, zorder=20)
    cg = cum_volume(g)
    if len(cg):
        a3.plot(cg.index, cg.values, color=C_GAUGE, lw=2.0, zorder=20)
    for ax in (a1, a2):
        ax.axhline(FF_THRESHOLD, color="#CC3311", lw=1.0, ls=":", zorder=2)
    a1.text(pd.Timestamp(FALL[0]) + pd.Timedelta(days=1), FF_THRESHOLD, "0.13 m³/s false-flow threshold", color="#CC3311", fontsize=8,
            va="bottom", ha="left")
    log_axis(a2, 30)
    for ax, t in [(a1, "A. Streamflow, ordinary scale"),
                  (a2, f"B. Same flow on a LOG scale (each gridline is 10 times the one below; flow under {LOG_FLOOR:g} m³/s, including zero, sits on the {LOG_FLOOR:g} line at the bottom)"),
                  (a3, "C. Running total of water past the outlet")]:
        style_time_axis(ax, FALL)
        ax.set_title(t, fontsize=10, loc="left")
        ax.tick_params(labelbottom=True)
    ylabel_flow(a1)
    ylabel_flow(a2)
    a3.set_ylabel("Cumulative volume (million m³)", fontsize=9)
    a3.set_xlabel("Date (Arizona time)", fontsize=9)
    fig.suptitle("Oct 1 - Dec 31 2016: the gauge sees almost no flow; what does channel loss do to the model's false flow?", fontsize=12)

    def lab(r):
        return f"Oct-Dec {fmt(num(r, C_FFV))} M m³, {fmt(num(r, C_FFDAYS), 0)} days above 0.13 m³/s"
    gl = f"USGS gauge:  Oct-Dec {fmt(gauge_info['vol'], 3)} M m³, largest flow {fmt(gauge_info['max'])} m³/s"
    bottom_legend(fig, legend_handles(df, sty, lab, gl), 0.97)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
# figure 2: tails
# ------------------------------------------------------------------------------------------------------------------
def fig_tails(df, flows, gauge, sty, out, dpi):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.6))
    for col, (title, rng) in enumerate(TAILS):
        top, bot = axes[0, col], axes[1, col]
        draw_lines(top, df, flows, sty, rng, transform=log_clip)
        for i, (_, r) in enumerate(df.iterrows()):
            c = cum_volume(clip(flows[r["run_id"]], rng))
            if len(c):
                cc, ls, lw = sty[r["run_id"]]
                bot.plot(c.index, c.values, color=cc, ls=ls, lw=lw, zorder=3 + i)
        g = clip(gauge.dropna(), rng)
        top.plot(g.index, log_clip(g).values, color=C_GAUGE, lw=2.2, zorder=20)
        cg = cum_volume(g)
        if len(cg):
            bot.plot(cg.index, cg.values, color=C_GAUGE, lw=2.2, zorder=20)
        log_axis(top, 100)
        top.set_title(title + " - flow, LOG scale", fontsize=10, loc="left")
        bot.set_title("Water past the outlet during this window only", fontsize=10, loc="left")
        style_time_axis(top, rng)
        style_time_axis(bot, rng)
        ylabel_flow(top)
        bot.set_ylabel("Cumulative volume (million m³)", fontsize=9)
    fig.suptitle("The recession tails: the gauge keeps carrying flow after the storms; what does channel loss do to the model's?", fontsize=12)

    def lab(r):
        return f"tail ratio {fmt(num(r, C_T1))} (Jan 26 - Feb 2),  {fmt(num(r, C_T2), 3)} (Mar 3 - 11)"
    bottom_legend(fig, legend_handles(df, sty, lab, "USGS gauge (tail ratio = model tail volume / gauge tail volume; 1.00 = match)"), 0.95)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
# figure 3: storms
# ------------------------------------------------------------------------------------------------------------------
def fig_storms(df, flows, gauge, sty, out, dpi):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.8))
    for ax, (title, rng) in zip(axes.ravel(), STORMS):
        draw_lines(ax, df, flows, sty, rng)
        g = clip(gauge.dropna(), rng)
        ax.plot(g.index, g.values, color=C_GAUGE, lw=2.2, zorder=20)
        style_time_axis(ax, rng)
        ax.set_title(title, fontsize=10, loc="left")
        ylabel_flow(ax)
    fig.suptitle("The four storms up close (each panel has its own vertical scale)", fontsize=12)

    def lab(r):
        return (f"volume ratio: flood week {fmt(num(r, C_VOL))}, 2nd Jan storm {fmt(num(r, C_J2), 3)}, "
                f"Feb 18-19 {fmt(num(r, C_F18))}, Feb 27-28 {fmt(num(r, C_F27))}")
    bottom_legend(fig, legend_handles(df, sty, lab, "USGS gauge (volume ratio = model / gauge; 1.00 = match)"), 0.95)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
# figure 4: whole window
# ------------------------------------------------------------------------------------------------------------------
def fig_overview(df, flows, gauge, sty, out, dpi):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(12.5, 9.0), sharex=True, gridspec_kw=dict(height_ratios=[2.0, 1.7]))
    draw_lines(a1, df, flows, sty, WHOLE, transform=log_clip)
    for i, (_, r) in enumerate(df.iterrows()):
        c = cum_volume(clip(flows[r["run_id"]], WHOLE))
        if len(c):
            cc, ls, lw = sty[r["run_id"]]
            a2.plot(c.index, c.values, color=cc, ls=ls, lw=lw, zorder=3 + i)
    g = clip(gauge.dropna(), WHOLE)
    a1.plot(g.index, log_clip(g).values, color=C_GAUGE, lw=1.6, zorder=20)
    cg = cum_volume(g)
    if len(cg):
        a2.plot(cg.index, cg.values, color=C_GAUGE, lw=2.0, zorder=20)
    log_axis(a1, 400)
    for ax, t in [(a1, f"A. Streamflow on a LOG scale (flow under {LOG_FLOOR:g} m³/s, including zero, sits on the {LOG_FLOOR:g} line at the bottom)"),
                  (a2, "B. Running total of water past the outlet")]:
        style_time_axis(ax, WHOLE)
        ax.set_title(t, fontsize=10, loc="left")
        ax.tick_params(labelbottom=True)
    ylabel_flow(a1)
    a2.set_ylabel("Cumulative volume (million m³)", fontsize=9)
    a2.set_xlabel("Date (Arizona time)", fontsize=9)
    fig.suptitle("The whole window (Oct 1 2016 - Apr 1 2017): where does the water go?", fontsize=12)

    def lab(r):
        return f"whole-window volume ratio {fmt(num(r, C_OM))}, hours with zero flow {fmt(num(r, C_DRY), 1)}%"
    bottom_legend(fig, legend_handles(df, sty, lab, "USGS gauge (volume ratio = model / gauge; 1.00 = match)"), 0.96)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Hydrograph figures for the channel-loss (Kchan) check on W4 (see the top of this file).")
    ap.add_argument("--csv", default="kchan_w4_summary.csv")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--obs", default=DEFAULT_OBS)
    ap.add_argument("--out-dir", default="kchan_w4_plots")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    sc = import_scorer()
    df = load_table(Path(a.csv))
    n_ref = int(df["is_ref"].sum())
    log(f"plot_hydrographs_kchan_w4_sycamore.py   {len(df)} runs in {a.csv}: "
        + ", ".join("loss off" if r else f"Kchan {k:g}" for r, k in zip(df["is_ref"], df["K"])))
    if n_ref == 0:
        log("  (no loss-off reference row in the CSV, so the gray dashed line is left out)")
    flows, bad = load_runs(sc, list(df["run_id"]), Path(a.results_dir))
    if bad:
        die("cannot read these outlet files (nothing was drawn):\n  " + "\n  ".join(bad))
    if not Path(a.obs).exists():
        die(f"gauge file not found: {a.obs}. Run from lab/, or give --obs.")
    gauge = sc.read_obs(a.obs)
    ref = df[df["is_ref"]]
    src = ref.iloc[0] if len(ref) else df.iloc[0]
    gauge_info = dict(vol=num(src, C_FFGV), max=num(src, C_FFGMAX))
    log(f"  read {len(flows)} run files; gauge has {int(gauge.notna().sum())} values")
    if a.dry_run:
        log("\nDRY RUN: every file was found and read; nothing was written.")
        return 0

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sty = styles(df)
    names = ["kchan_hydro_1_fall.png", "kchan_hydro_2_tails.png", "kchan_hydro_3_storms.png", "kchan_hydro_4_overview.png"]
    fig_fall(df, flows, gauge, sty, gauge_info, out / names[0], a.dpi)
    fig_tails(df, flows, gauge, sty, out / names[1], a.dpi)
    fig_storms(df, flows, gauge, sty, out / names[2], a.dpi)
    fig_overview(df, flows, gauge, sty, out / names[3], a.dpi)
    log("\nWROTE " + "\n      ".join(str(out / n) for n in names))
    log("Nothing in results/ was touched.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
