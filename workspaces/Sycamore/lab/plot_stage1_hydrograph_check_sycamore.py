#!/usr/bin/env python3
"""
plot_stage1_hydrograph_check_sycamore.py  --  Stage 1 hydrograph check (the design's "hydrograph gate").

READ-ONLY. It runs no tRIBS and changes none of your files. It reads the outlet files already in results/,
the gauge file, and the score CSVs from Stage 1, and writes two figures plus two small text/CSV files into a
NEW folder, hydrograph_check_stage1/.

TWO FIGURES, ONE PER STAGE 1 DECISION
-------------------------------------
  Figure 1  fig_stage1_start_state.png    Decision: where the starting water table sits does not matter.
      Two panels, gauge in black in both.
        top     W3 window (Oct 1 2016 - Jan 25 2017), LOG flow axis: the baseline (thick grey), water table 0.70,
                water table 1.00, and bedrock +0 m (dashed; the one run that IS a different model, so it shows
                that the plot would reveal a difference if there were one). The log axis is the only way to see
                the low flows in Oct-Dec; the dotted line is 0.13 m3/s, the gauge's Oct-Dec maximum.
        bottom  the flood week (Jan 18-25), LINEAR flow axis, same lines.
      Under the figure: each run's flood peak, flood-week volume and Oct-Dec volume next to the baseline's.

  Figure 2  fig_stage1_proxy_flood_week.png    Decision: W1 (and W2) can stand in for W4 as a screening window.
      Six small panels, one per Ks / PsiB setting, ordered from the setting with the smallest W4 flood-week volume
      error to the largest. Each panel: the gauge (black), the W4 run (solid blue) and the W1 run (dashed orange),
      all on the SAME axes so the panels can be compared. W2 is not drawn: its flood week is identical to W1's
      (same start, same water table, and tRIBS is deterministic). Each panel prints the volume ratio and peak
      error for W4 and W1, and where each ranks among the six settings (1 = smallest error).

BEFORE IT DRAWS ANYTHING it checks itself, and stops if a check fails: it recomputes the flood-week volume ratio and
peak error of every plotted run from the outlet file, exactly the way score_run_sycamore.py does (it imports the
scorer's own readers), and compares them with the numbers already in your score CSVs. If they do not match to 1e-6,
the figures could not be trusted, so nothing is drawn.

HOW TO RUN  (from workspaces/Sycamore/lab/, where the other Stage 1 scripts are)
    python plot_stage1_hydrograph_check_sycamore.py --dry-run      # checks the inputs and the numbers, draws nothing
    python plot_stage1_hydrograph_check_sycamore.py                # draws and writes (takes well under a minute)

It needs score_run_sycamore.py in the same folder (it uses its file readers).

WHAT IT READS
    results/<run>_Outlet.qout                     the outlet files
    ../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv   the gauge (through the scorer)
    proxy_manifest.json                           which run is which setting and window
    score_out_proxy/scores_proxy.csv              for the self-check
    score_out_start_state/scores_start_state.csv  for the self-check (any run not in a CSV is simply not cross-checked)

WHAT IT WRITES  (new folder hydrograph_check_stage1/)
    fig_stage1_start_state.png   fig_stage1_proxy_flood_week.png
    hydrograph_check_numbers.csv   one row per plotted run (flood-week volume ratio, peak error, volumes)
    hydrograph_check_report.txt    a copy of everything printed

THE RUNS DRAWN IN FIGURE 1 are listed in START_RUNS below. To add the 0.82 water table or the other bedrock runs,
add a line there. The baseline is the Oct-Mar run used in the Stage 1 test (its first 2,784 hours are the W3 window).

READ THIS BEFORE QUOTING ANYTHING FROM THE FIGURES
    * These are model outputs against a gauge with known gaps (1 of 5 WY2017 peaks is USGS-estimated; resolution below
      0.5 m3/s is unchecked). The gauge in the low-flow part of figure 1 is the weakest evidence in the record.
    * Figure 2 only covers the January flood week, the same slice the rank check used. It cannot show whether W2
      can stand in for W4 on whole-season volume or Oct-Dec false flow (a caveat already recorded in the design doc).
    * W1 starts cold on Jan 18 (water table assumed 0.82 of bedrock, unconfirmed); W4 starts Oct 1 at 0.95 and
      reaches Jan 18 after about 3.5 months of simulated rain. Any difference between W1 and W4 in figure 2 is the
      effect of that cold start.
"""
from __future__ import annotations

import argparse
import io
import json
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

LAB = Path(__file__).resolve().parent
os.chdir(LAB)
sys.path.insert(0, str(LAB))
try:
    import score_run_sycamore as S
except ImportError as exc:                                                       # noqa: BLE001
    print(f"STOPPED: this script uses the file readers in score_run_sycamore.py and could not import it ({exc}).\n"
          f"Put both scripts in the same lab/ folder and run from there:\n    cd {LAB}")
    sys.exit(1)

# ----------------------------------------------------------------------------------------------------------
# What is drawn
# ----------------------------------------------------------------------------------------------------------
BASELINE_ID = "Sycamore_octmar_ks0p2_cv6_psib0p5_nocl"      # Ks x0.2, PsiB x0.5, water table 0.95, bedrock +1 m; Oct 1 start
BASELINE_LABEL = "baseline (water table 0.95, bedrock +1 m)"
# (short tag, label, run id, color, line style, line width)
START_RUNS = [
    ("gw0p70", "water table 0.70 of bedrock", "W3_ks0p2_cv6_psib0p5_gw0p70_nocl", "#d95f02", "-", 1.5),
    ("gw1p00", "water table 1.00 of bedrock", "W3_ks0p2_cv6_psib0p5_gw1p00_nocl", "#1b9e77", "-", 1.5),
    ("bed0p0", "bedrock +0 m (a different model)", "W3_ks0p2_cv6_psib0p5_bed0p0_nocl", "#c0392b", (0, (4, 2)), 1.7),
]
SETTING_LABELS = {          # same numbering as run_proxy_sycamore.py
    1: ("Ks x0.2, PsiB x0.5", "baseline"),
    2: ("Ks x0.2, PsiB x1.0", ""),
    3: ("Ks x0.1, PsiB x0.5", ""),
    4: ("Ks x0.3, PsiB x0.5", ""),
    5: ("Ks x0.3, PsiB x1.0", ""),
    6: ("Ks x0.1, PsiB x1.0", ""),
}
RESULTS = "results"
MANIFEST = "proxy_manifest.json"
PROXY_CSV = "score_out_proxy/scores_proxy.csv"
START_CSV = "score_out_start_state/scores_start_state.csv"
OUT_DIR = "hydrograph_check_stage1"
FIG1, FIG2 = "fig_stage1_start_state.png", "fig_stage1_proxy_flood_week.png"
NUMBERS_CSV, REPORT_TXT = "hydrograph_check_numbers.csv", "hydrograph_check_report.txt"
FLOOR = 1e-3                 # flows below this are drawn at this value on the log axis (a log axis cannot show zero)
CHECK_TOL = 1e-6
INK, INK_SOFT = "#222222", "#666666"
BLUE, ORANGE, GREY = "#1f77b4", "#d95f02", "#8c8c8c"

_LINES: list[str] = []


def say(text: str = "") -> None:
    print(text, flush=True)
    _LINES.append(text)


def sec(title: str) -> None:
    say("")
    say("=" * 100)
    say(title)
    say("=" * 100)


class Stop(Exception):
    """A problem the user can fix; carries the message to print."""


# ----------------------------------------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------------------------------------
def read_run(run_id: str, window: str) -> pd.Series:
    path, _, _, err = S.find_run_files(run_id, Path(RESULTS))
    if path is None:
        raise Stop(f"Outlet file for {run_id} not found: {err}")
    q = S.read_qout(path, S.RUN_WINDOWS[window][0])
    if q.empty:
        raise Stop(f"{path} has no readable rows.")
    return q


def load_proxy_runs() -> dict:
    """{(setting, window): run_id}, the same way analyze_proxy_sycamore.py picks them (a run made again beats a reused one)."""
    mp = Path(MANIFEST)
    if not mp.is_file():
        raise Stop(f"{MANIFEST} not found in {LAB}. Run this script from the lab folder where the proxy check was run.")
    man = json.loads(mp.read_text())
    pick: dict = {}
    for rid, rec in man["runs"].items():
        if rec.get("status") not in ("done", "reused") or "setting" not in rec:
            continue
        key = (int(rec["setting"]), rec["window"])
        rank = 1 if rec.get("kind") == "run" else 0
        if key not in pick or rank > pick[key][0]:
            pick[key] = (rank, rid)
    runs = {k: v[1] for k, v in pick.items()}
    missing = [(s, w) for s in SETTING_LABELS for w in ("W1", "W4") if (s, w) not in runs]
    if missing:
        raise Stop("The proxy manifest has no finished run for: " + ", ".join(f"setting {s} {w}" for s, w in missing)
                   + ".\nFigure 2 needs all six settings in W1 and W4.")
    return runs


def load_score_tables() -> pd.DataFrame:
    parts = []
    for p in (PROXY_CSV, START_CSV):
        if Path(p).is_file():
            df = pd.read_csv(p)
            if "run_id" in df.columns:
                parts.append(df.set_index("run_id"))
    if not parts:
        return pd.DataFrame()
    allr = pd.concat(parts)
    return allr[~allr.index.duplicated(keep="last")]


# ----------------------------------------------------------------------------------------------------------
# Numbers (the scorer's own method for the flood week)
# ----------------------------------------------------------------------------------------------------------
def flood_week_numbers(q: pd.Series, G, prof: dict) -> dict:
    a, b = prof["events"]["flood_week"]
    excl = [tuple(e) for e in prof["exclude"]]
    M = S.Flow(q)
    if not S.covers(M, a, b):
        raise Stop(f"A run does not cover the flood week ({a} to {b}); first sample {q.index[0]}, last {q.index[-1]}.")
    mm = M.volume(a, b, excl, match_to=G) / 1e6
    mr = M.volume(a, b) / 1e6
    g = G.volume(a, b, excl) / 1e6
    ev_m, ev_g = M.window(a, b, excl), G.window(a, b, excl)
    pm, pg = float(ev_m.max()), float(ev_g.max())
    return dict(vol_model_matched_Mm3=mm, vol_model_raw_Mm3=mr, vol_gauge_Mm3=g, vol_ratio=mm / g,
                peak_model_m3s=pm, peak_gauge_m3s=pg, peak_err_pct=(pm / pg - 1) * 100,
                peak_time=ev_m.idxmax())


def octdec_volume(q: pd.Series, prof: dict) -> float:
    a, b = prof["false_flow"]["octdec"]
    return S.Flow(q).volume(a, b) / 1e6


# ----------------------------------------------------------------------------------------------------------
# Self-check
# ----------------------------------------------------------------------------------------------------------
def self_check(table: pd.DataFrame, rows: list[dict]) -> None:
    sec("0. Self-check: do the outlet files give back the numbers already in the score CSVs?")
    if table.empty:
        say("  WARN  neither score CSV was found, so nothing could be cross-checked.")
        return
    say(f"  {'run':<44}{'vol ratio (here)':>17}{'(CSV)':>9}{'peak % (here)':>15}{'(CSV)':>9}   result")
    fails, n_ok, n_skip = [], 0, 0
    seen = set()
    for r in rows:
        rid = r["run_id"]
        if rid in seen:                       # the reused W4 baseline appears twice; check it once
            continue
        seen.add(rid)
        if rid not in table.index or "flood_week__vol_ratio" not in table.columns:
            say(f"  {rid:<44}{'-':>17}{'-':>9}{'-':>15}{'-':>9}   not in a score CSV, not cross-checked")
            n_skip += 1
            continue
        cv, cp = float(table.at[rid, "flood_week__vol_ratio"]), float(table.at[rid, "flood_week__peak_err_pct"])
        bad = not (abs(r["vol_ratio"] - cv) <= CHECK_TOL * (1 + abs(cv)) and abs(r["peak_err_pct"] - cp) <= CHECK_TOL * (1 + abs(cp)))
        if bad:
            fails.append(f"{rid}: volume ratio here {r['vol_ratio']:.8f}, CSV {cv:.8f}; peak error here {r['peak_err_pct']:.6f}, CSV {cp:.6f}")
        say(f"  {rid:<44}{r['vol_ratio']:>17.6f}{cv:>9.4f}{r['peak_err_pct']:>15.4f}{cp:>9.3f}   {'FAIL' if bad else 'ok'}")
        n_ok += 0 if bad else 1
    if fails:
        raise Stop("The outlet files do NOT reproduce the scores already on record, so nothing drawn from them can be trusted:\n    - "
                   + "\n    - ".join(fails)
                   + "\n  Possible causes: a run was rerun after it was scored, or your score_run_sycamore.py differs from the one that scored them."
                   + "\n  Nothing was drawn or written. Paste this section to Claude.")
    say(f"  PASS  {n_ok} run(s) reproduce their scored flood-week numbers" + (f"; {n_skip} not in a CSV." if n_skip else "."))
    if n_skip and n_ok == 0:
        say("  WARN  no plotted run could be cross-checked.")


# ----------------------------------------------------------------------------------------------------------
# Drawing
# ----------------------------------------------------------------------------------------------------------
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9.5, "axes.spines.top": False, "axes.spines.right": False})
    return plt


def _style(ax):
    ax.grid(True, color="#e5e5e5", lw=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_SOFT, labelsize=8.5)


def _floor(s: pd.Series) -> pd.Series:
    return s.clip(lower=FLOOR)


def figure_start_state(path: Path, dpi: int, G_q: pd.Series, base: pd.Series, runs: list, prof: dict, table_lines: list[str]) -> None:
    plt = _plt()
    import matplotlib.dates as mdates
    w3a, w3b = pd.Timestamp(S.RUN_WINDOWS["W3"][0]), pd.Timestamp(S.RUN_WINDOWS["W3"][1])
    fa, fb = (pd.Timestamp(x) for x in prof["events"]["flood_week"])
    thr = prof.get("false_flow_threshold_m3s")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11.5, 10.2), gridspec_kw=dict(height_ratios=[1.15, 1.0]))
    fig.patch.set_facecolor("white")
    for ax in (ax1, ax2):
        _style(ax)

    rd = RUN_SERIES          # {tag: outlet series}, filled in main()

    def draw(ax, lo, hi, log):
        g = G_q.loc[lo:hi]
        f = _floor if log else (lambda s: s)
        ax.plot(base.loc[lo:hi].index, f(base.loc[lo:hi]), color=GREY, lw=5.0, alpha=0.55, solid_capstyle="round", label=BASELINE_LABEL, zorder=2)
        for tag, label, rid, color, style, lw in runs:
            q = rd[tag].loc[lo:hi]
            ax.plot(q.index, f(q), color=color, lw=lw, linestyle=style, label=label, zorder=3)
        ax.plot(g.index, f(g), color=INK, lw=1.2, label="USGS gauge 09510200", zorder=4)

    draw(ax1, w3a, w3b, True)
    ax1.set_yscale("log")
    top = max(float(base.loc[w3a:w3b].max()), float(G_q.loc[w3a:w3b].max()), *(float(rd[t].loc[w3a:w3b].max()) for t, *_ in runs))
    ax1.set_ylim(FLOOR * 0.8, top * 3)
    if thr:
        ax1.axhline(thr, color=INK_SOFT, linestyle=(0, (1, 2)), lw=1.1)
        ax1.text(w3a, thr, f"  {thr:g} m3/s: the gauge's Oct-Dec maximum", fontsize=8.2, color=INK_SOFT, va="bottom", ha="left")
    ax1.axvspan(fa, fb, color="#f1e6b8", alpha=0.55, lw=0, zorder=1)
    ax1.set_xlim(w3a, w3b)
    ax1.xaxis.set_major_locator(mdates.MonthLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax1.set_ylabel("discharge (m3/s), log scale", color=INK_SOFT)
    ax1.set_title("W3, Oct 1 2016 to Jan 25 2017, log axis (shaded: the flood week).  Flows under 0.001 m3/s are drawn at 0.001.",
                  loc="left", fontsize=10.5, color=INK)
    ax1.legend(loc="upper left", frameon=True, framealpha=0.92, edgecolor="#dddddd", fontsize=8.8)

    draw(ax2, fa, fb, False)
    ax2.set_xlim(fa, fb)
    ax2.xaxis.set_major_locator(mdates.DayLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax2.set_ylim(0, None)
    ax2.set_ylabel("discharge (m3/s), linear scale", color=INK_SOFT)
    ax2.set_title(f"Flood week, {fa:%b %d} to {fb:%b %d} 2017, linear axis", loc="left", fontsize=10.5, color=INK)

    fig.suptitle("Stage 1 start state: does where the water table starts change the hydrograph?", fontsize=12.5, color=INK, y=0.995)
    fig.text(0.015, 0.012, "\n".join(table_lines), ha="left", va="bottom", fontsize=8.2, color=INK, family="monospace", linespacing=1.45)
    fig.tight_layout(rect=(0, 0.13, 1, 0.985))
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def figure_proxy(path: Path, dpi: int, G_q: pd.Series, order: list, panels: dict, prof: dict) -> None:
    plt = _plt()
    import matplotlib.dates as mdates
    fa, fb = (pd.Timestamp(x) for x in prof["events"]["flood_week"])
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.6), sharex=True, sharey=True)
    fig.patch.set_facecolor("white")
    g = G_q.loc[fa:fb]
    ymax = max(float(g.max()), *(float(panels[s][w]["q"].loc[fa:fb].max()) for s in order for w in ("W1", "W4"))) * 1.10
    for ax, s in zip(axes.ravel(), order):
        _style(ax)
        p = panels[s]
        ax.plot(g.index, g, color=INK, lw=1.8, label="USGS gauge", zorder=4)
        ax.plot(p["W4"]["q"].loc[fa:fb].index, p["W4"]["q"].loc[fa:fb], color=BLUE, lw=2.0, label="W4 (Oct 1 start, water table 0.95)", zorder=3)
        ax.plot(p["W1"]["q"].loc[fa:fb].index, p["W1"]["q"].loc[fa:fb], color=ORANGE, lw=1.7, linestyle=(0, (4, 2)),
                label="W1 (Jan 18 start, water table 0.82)", zorder=3)
        name, tag = SETTING_LABELS[s]
        ax.set_title(f"{name}" + (f"  ({tag})" if tag else "")
                     + f"\nvolume rank: W4 {p['W4']['rank_vol']}, W1 {p['W1']['rank_vol']}      "
                       f"peak rank: W4 {p['W4']['rank_peak']}, W1 {p['W1']['rank_peak']}",
                     loc="left", fontsize=9.2, color=INK)
        ax.text(0.98, 0.97,
                f"W4  volume {p['W4']['vol_ratio']:.3f}  peak {p['W4']['peak_err_pct']:+.1f}%\n"
                f"W1  volume {p['W1']['vol_ratio']:.3f}  peak {p['W1']['peak_err_pct']:+.1f}%",
                transform=ax.transAxes, ha="right", va="top", fontsize=8.6, color=INK, family="monospace",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#dddddd", alpha=0.92))
        ax.set_xlim(fa, fb)
        ax.set_ylim(0, ymax)
        ax.set_xticks(list(pd.date_range(fa, fb, freq="2D")))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    for ax in axes[:, 0]:
        ax.set_ylabel("discharge (m3/s)", color=INK_SOFT)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=9.5)
    fig.suptitle("Stage 1 proxy check: does the cheap window (W1) give the same flood-week picture as W4?  "
                 "(panels ordered by W4 volume error, smallest first)", fontsize=12, color=INK, y=0.995)
    fig.text(0.5, 0.045, "Volume = model / gauge over Jan 18-25 (1.000 is perfect); peak = model peak / gauge peak - 1. "
             "Ranks: 1 = smallest error among the six settings. W2's flood week is identical to W1's, so it is not drawn.\n"
             "Any W1-W4 difference is the effect of the cold start (W1 begins Jan 18 with the water table assumed at 0.82 of bedrock; "
             "W4 reaches Jan 18 after about 3.5 months of simulated rain).",
             ha="center", va="bottom", fontsize=8.3, color=INK_SOFT, linespacing=1.5)
    fig.tight_layout(rect=(0, 0.095, 1, 0.975))
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------------------
RUN_SERIES: dict = {}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage 1 hydrograph check: two figures, read-only.")
    ap.add_argument("--dry-run", action="store_true", help="check the inputs and the numbers, then stop (writes nothing)")
    ap.add_argument("--out-dir", default=OUT_DIR)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--obs", default=S.DEFAULT_OBS, help=f"gauge file (default {S.DEFAULT_OBS})")
    args = ap.parse_args(argv)
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be between 72 and 600")
    t0 = time.time()
    warnings.filterwarnings("ignore", category=FutureWarning)
    prof = S.load_profile("wy2017", None)
    say(f"Stage 1 hydrograph check   python {sys.version.split()[0]}, pandas {pd.__version__}   {pd.Timestamp.now():%Y-%m-%d %H:%M}")
    say(f"  scorer in use: score_run_sycamore.py v{S.VERSION}; excluded gauge periods: {prof['exclude'] or 'none'}")

    gauge_q = S.read_obs(args.obs)
    G = S.Flow(gauge_q)
    say(f"  gauge {args.obs}: {gauge_q.index[0]:%b %d %Y} to {gauge_q.index[-1]:%b %d %Y}, {int(gauge_q.notna().sum())} values")

    # ---- figure 1 inputs
    base = read_run(BASELINE_ID, "W4")
    for tag, _, rid, *_ in START_RUNS:
        RUN_SERIES[tag] = read_run(rid, "W3")
    # ---- figure 2 inputs
    proxy_ids = load_proxy_runs()
    proxy_q = {(s, w): read_run(proxy_ids[(s, w)], w) for s in SETTING_LABELS for w in ("W1", "W4")}

    # ---- numbers and self-check
    rows = []
    nums_base = flood_week_numbers(base, G, prof)
    rows.append(dict(run_id=BASELINE_ID, role="start-state baseline", window="W4", **nums_base))
    for tag, label, rid, *_ in START_RUNS:
        rows.append(dict(run_id=rid, role=f"start-state {tag}", window="W3", **flood_week_numbers(RUN_SERIES[tag], G, prof)))
    panels: dict = {}
    for s in SETTING_LABELS:
        panels[s] = {}
        for w in ("W1", "W4"):
            n = flood_week_numbers(proxy_q[(s, w)], G, prof)
            panels[s][w] = dict(q=proxy_q[(s, w)], run_id=proxy_ids[(s, w)], **n)
            rows.append(dict(run_id=proxy_ids[(s, w)], role=f"proxy setting {s}", window=w, **n))
    table = load_score_tables()
    try:
        self_check(table, rows)
    except Stop as exc:
        say("\n" + "!" * 94)
        say(f"STOPPED: {exc}")
        say("!" * 94)
        return 3

    # ---- ranks for figure 2 (1 = smallest error; ties share the smaller rank)
    for w in ("W1", "W4"):
        vol_err = pd.Series({s: abs(panels[s][w]["vol_ratio"] - 1) for s in panels})
        peak_err = pd.Series({s: abs(panels[s][w]["peak_err_pct"]) for s in panels})
        rv, rp = vol_err.rank(method="min").astype(int), peak_err.rank(method="min").astype(int)
        for s in panels:
            panels[s][w]["rank_vol"], panels[s][w]["rank_peak"] = int(rv[s]), int(rp[s])
    order = sorted(panels, key=lambda s: (abs(panels[s]["W4"]["vol_ratio"] - 1), s))

    sec("1. Start state (figure 1): each run against the baseline, flood week and Oct-Dec")
    ob = octdec_volume(base, prof)
    fb_vol, fb_peak = nums_base["vol_model_raw_Mm3"], nums_base["peak_model_m3s"]
    say(f"  gauge flood-week peak {nums_base['peak_gauge_m3s']:.1f} m3/s, flood-week volume {nums_base['vol_gauge_Mm3']:.3f} M m3")
    head = f"{'run':<34}{'flood peak':>12}{'vs base':>9}{'flood vol':>11}{'vs base':>9}{'Oct-Dec vol':>13}{'vs base':>9}"
    say("  " + head)
    table_lines = ["Model flood-week peak (m3/s), flood-week volume and Oct-Dec volume (M m3), each against the baseline's.",
                   head]

    def pct(x, ref):
        return f"{100 * (x / ref - 1):>+8.1f}%" if ref else f"{'n/a':>9}"

    def line(label, n, od):
        return (f"{label:<34}{n['peak_model_m3s']:>12.1f}{pct(n['peak_model_m3s'], fb_peak)}"
                f"{n['vol_model_raw_Mm3']:>11.3f}{pct(n['vol_model_raw_Mm3'], fb_vol)}"
                f"{od:>13.3f}{pct(od, ob)}")

    lb = line("baseline", nums_base, ob)
    say("  " + lb)
    table_lines.append(lb)
    for tag, label, rid, *_ in START_RUNS:
        n = next(r for r in rows if r["run_id"] == rid)
        od = octdec_volume(RUN_SERIES[tag], prof)
        ln = line(label[:33], n, od)
        say("  " + ln)
        table_lines.append(ln)

    sec("2. Proxy check (figure 2): flood week, W4 against W1, ordered by W4 volume error")
    say(f"  {'setting':<24}{'W4 vol':>9}{'W1 vol':>9}{'W4 peak%':>10}{'W1 peak%':>10}   ranks vol (W4/W1)   ranks peak (W4/W1)")
    for s in order:
        p = panels[s]
        say(f"  {SETTING_LABELS[s][0]:<24}{p['W4']['vol_ratio']:>9.3f}{p['W1']['vol_ratio']:>9.3f}{p['W4']['peak_err_pct']:>+10.1f}{p['W1']['peak_err_pct']:>+10.1f}"
            f"   {p['W4']['rank_vol']} / {p['W1']['rank_vol']}                 {p['W4']['rank_peak']} / {p['W1']['rank_peak']}")

    if args.dry_run:
        sec("DRY RUN: nothing was drawn or written")
        say("Ready:  python plot_stage1_hydrograph_check_sycamore.py")
        return 0

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sec(f"Files written to {out}/")
    figure_start_state(out / FIG1, args.dpi, gauge_q, base, START_RUNS, prof, table_lines)
    say(f"  saved {FIG1}")
    figure_proxy(out / FIG2, args.dpi, gauge_q, order, panels, prof)
    say(f"  saved {FIG2}")
    df = pd.DataFrame(rows)
    df["peak_time"] = df["peak_time"].map(lambda v: v.strftime("%Y-%m-%d %H:%M") if isinstance(v, pd.Timestamp) else v)
    df.to_csv(out / NUMBERS_CSV, index=False, float_format="%.8g")
    say(f"  saved {NUMBERS_CSV}  ({len(df)} rows)")
    say(f"\nDone in {time.time() - t0:.0f} s. This script only reads your results; it changed none of them.")
    (out / REPORT_TXT).write_text("\n".join(_LINES) + "\n")
    print(f"  saved {REPORT_TXT}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Stop as exc:
        print("\n" + "!" * 94)
        print(f"STOPPED: {exc}")
        print("!" * 94)
        print("Nothing was changed by this refusal. Fix the cause above, then run the command again.")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nStopped by you. Nothing in your results was changed.")
        sys.exit(130)
