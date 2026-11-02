#!/usr/bin/env python3
"""
select_w2_candidates_lhs_volume_sycamore.py  --  Stage 3 -> Stage 4: choose which finished W1 runs go on to W2.

THE RULE (so the choice is repeatable and nobody picked favorites):

  1. Candidates are ONLY runs from the 150-run Latin hypercube that already finished. No new Ks, f or PsiB values are made up.
  2. MAIN GROUP: runs whose MAIN-STORM volume ratio (Jan 18 to Jan 23 12:00) is inside the band (default 0.9 to 1.1).
     From those, pick N (default 10) that are as SPREAD OUT as possible in f and PsiB. Why spread, not "best score"? Inside the
     band, W1 cannot tell f and PsiB apart (the three knobs act as one volume dial), so ranking by W1 KGE would be picking noise.
     A longer window (W2) is what can separate them, so we give W2 a wide range of f and PsiB to work with.
  3. CONTROL GROUP: runs that pass the WHOLE-WEEK volume band but NOT the main-storm band (default 3), also spread out.
     These are the "right week total, main storm too wet" runs. If W2 (which includes the February storms) prefers the main group,
     that supports using the main-storm volume as the target; if it prefers the controls, we learn that early.

HOW "SPREAD OUT" IS DONE (farthest-point selection, no randomness):
  Put every run in the group on a plane with ln(f) across and ln(PsiB) up, each axis scaled by the group's own range so both count
  equally. Start with the run whose volume ratio is closest to 1.0. Then repeatedly add the run that is FARTHEST from all runs
  already chosen. Repeat until N are chosen.

Input is lhs_volume_mainstorm.csv, written by score_mainstorm_lhs_volume_sycamore.py. This script reads only that file and
writes only into lhs_volume_plots/:
    w2_candidates.csv      the chosen runs (group, draw, Ks, f, PsiB, volume ratios, peak, arrival, KGE)
    w2_candidates.png      where the chosen runs sit among all the runs in the group (f vs PsiB, and Ks vs f)

Run from workspaces/Sycamore/lab/:
    python -u select_w2_candidates_lhs_volume_sycamore.py --dry-run     # prints the choice, writes nothing
    python -u select_w2_candidates_lhs_volume_sycamore.py               # prints it and writes the CSV and picture

Options: --csv FILE (lhs_volume_plots/lhs_volume_mainstorm.csv)  --n-main 10  --n-control 3  --band 0.10
         --out-dir DIR (lhs_volume_plots)  --dpi 150
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

WEEK_VOL = "flood_week__vol_ratio"
NEED = ["run_id", "draw", "KS_MULT", "F_MULT", "PSIB_MULT", WEEK_VOL, "main_vol_ratio", "main_peak_err_pct",
        "main_arrival_err_h", "main_kge2012"]


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    print(f"\nERROR: {msg}\n", file=sys.stderr, flush=True)
    sys.exit(2)


def farthest_points(g: pd.DataFrame, n: int, seed_col: str) -> list:
    """Row labels of n runs spread out in (ln f, ln PsiB). Seed = the run whose seed_col is closest to 1.0. Deterministic."""
    if len(g) <= n:
        return list(g.index)
    x = np.log(g["F_MULT"].to_numpy(float))
    y = np.log(g["PSIB_MULT"].to_numpy(float))
    x = (x - x.min()) / (x.max() - x.min()) if x.max() > x.min() else np.zeros_like(x)
    y = (y - y.min()) / (y.max() - y.min()) if y.max() > y.min() else np.zeros_like(y)
    pts = np.column_stack([x, y])
    chosen = [int(np.argmin(np.abs(g[seed_col].to_numpy(float) - 1.0)))]
    mind = np.linalg.norm(pts - pts[chosen[0]], axis=1)
    while len(chosen) < n:
        nxt = int(np.argmax(mind))
        chosen.append(nxt)
        mind = np.minimum(mind, np.linalg.norm(pts - pts[nxt], axis=1))
    return [g.index[i] for i in chosen]


def draw_picture(t, main, ctrl, sel_main, sel_ctrl, out, dpi):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
    panels = [("F_MULT", "PSIB_MULT", "f multiplier (log scale)", "PsiB multiplier (log scale)", "f and PsiB"),
              ("F_MULT", "KS_MULT", "f multiplier (log scale)", "Ks multiplier (log scale)", "f and Ks")]
    for ax, (xc, yc, xl, yl, ttl) in zip(axes, panels):
        ax.scatter(t[xc], t[yc], s=14, color="#DDDDDD", label=f"all {len(t)} runs", zorder=1)
        ax.scatter(main[xc], main[yc], s=34, color="#9DB8D6", label=f"main-storm band ({len(main)})", zorder=2)
        ax.scatter(ctrl[xc], ctrl[yc], s=34, color="#E5B7BF", label=f"week band only ({len(ctrl)})", zorder=2)
        ax.scatter(sel_main[xc], sel_main[yc], s=95, color="#0077BB", edgecolor="black", lw=0.8, label=f"chosen: main group ({len(sel_main)})", zorder=4)
        ax.scatter(sel_ctrl[xc], sel_ctrl[yc], s=95, marker="s", color="#CC3311", edgecolor="black", lw=0.8, label=f"chosen: control group ({len(sel_ctrl)})", zorder=4)
        for _, r in pd.concat([sel_main, sel_ctrl]).iterrows():
            ax.annotate(str(int(r["draw"])), (r[xc], r[yc]), textcoords="offset points", xytext=(6, 5), fontsize=8, zorder=5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xl, fontsize=9)
        ax.set_ylabel(yl, fontsize=9)
        ax.set_title(ttl, fontsize=10, loc="left")
        ax.grid(alpha=0.25, which="both")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, fontsize=8, frameon=False)
    fig.suptitle("Runs chosen for W2 (numbers are draw numbers)", fontsize=12)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description="Choose the W1 runs that go on to W2 (see the top of this file).")
    ap.add_argument("--csv", default="lhs_volume_plots/lhs_volume_mainstorm.csv")
    ap.add_argument("--n-main", type=int, default=10)
    ap.add_argument("--n-control", type=int, default=3)
    ap.add_argument("--band", type=float, default=0.10)
    ap.add_argument("--out-dir", default="lhs_volume_plots")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    lo, hi = 1 - a.band, 1 + a.band

    path = Path(a.csv)
    if not path.exists():
        die(f"{path} not found. Run score_mainstorm_lhs_volume_sycamore.py first (from lab/), or give --csv.")
    t = pd.read_csv(path)
    miss = [c for c in NEED if c not in t.columns]
    if miss:
        die(f"{path} is missing columns {miss}")
    t = t[np.isfinite(t["main_vol_ratio"])].reset_index(drop=True)

    in_week = t[WEEK_VOL].between(lo, hi)
    in_main = t["main_vol_ratio"].between(lo, hi)
    main = t[in_main].copy()
    ctrl = t[in_week & ~in_main].copy()
    log(f"select_w2_candidates_lhs_volume_sycamore.py   {len(t)} runs in {path}")
    log(f"  volume band {lo:.2f} to {hi:.2f}:  main-storm band {len(main)} runs;  week band only {len(ctrl)} runs")
    if main.empty:
        die("no run is inside the main-storm band; nothing to choose from")

    sel_main = main.loc[farthest_points(main, a.n_main, "main_vol_ratio")].copy()
    sel_ctrl = ctrl.loc[farthest_points(ctrl, a.n_control, WEEK_VOL)].copy() if len(ctrl) else ctrl.iloc[0:0].copy()
    sel_main["group"], sel_ctrl["group"] = "main", "control"

    def show(title, g, full):
        log(f"\n{title}")
        log(f"  {'draw':>4}  {'Ks':>6}  {'f':>5}  {'PsiB':>5}  {'week vol':>8}  {'main vol':>8}  {'peak %':>7}  {'arrival h':>9}  {'KGE':>6}")
        for _, r in g.sort_values("F_MULT").iterrows():
            log(f"  {int(r['draw']):>4}  {r['KS_MULT']:>6.3f}  {r['F_MULT']:>5.2f}  {r['PSIB_MULT']:>5.2f}  {r[WEEK_VOL]:>8.2f}  {r['main_vol_ratio']:>8.2f}"
                f"  {r['main_peak_err_pct']:>+7.0f}  {r['main_arrival_err_h']:>+9.2f}  {r['main_kge2012']:>6.3f}")
        log(f"  range chosen vs range available:  f {g['F_MULT'].min():.2f} to {g['F_MULT'].max():.2f} (of {full['F_MULT'].min():.2f} to {full['F_MULT'].max():.2f});  "
            f"PsiB {g['PSIB_MULT'].min():.2f} to {g['PSIB_MULT'].max():.2f} (of {full['PSIB_MULT'].min():.2f} to {full['PSIB_MULT'].max():.2f});  "
            f"Ks {g['KS_MULT'].min():.3f} to {g['KS_MULT'].max():.3f} (of {full['KS_MULT'].min():.3f} to {full['KS_MULT'].max():.3f})")

    show(f"MAIN GROUP: {len(sel_main)} of {len(main)} runs inside the main-storm band, spread across f and PsiB", sel_main, main)
    if len(sel_ctrl):
        show(f"CONTROL GROUP: {len(sel_ctrl)} of {len(ctrl)} runs inside the week band only (main storm too wet)", sel_ctrl, ctrl)
    else:
        log("\nCONTROL GROUP: none (no run passes the week band without also passing the main-storm band)")
    draws = ",".join(str(int(d)) for d in list(sel_main["draw"]) + list(sel_ctrl["draw"]))
    log(f"\nDRAW NUMBERS, main group first then controls:  {draws}")
    log("  Nothing here ranks runs by W1 score: the choice is by spread in f and PsiB, starting from the run closest to volume 1.0.")

    if a.dry_run:
        log("\nDRY RUN: nothing was written.")
        return 0
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    keep = ["group", "draw", "run_id", "KS_MULT", "F_MULT", "PSIB_MULT", WEEK_VOL, "main_vol_ratio", "main_peak_err_pct", "main_arrival_err_h", "main_kge2012"]
    res = pd.concat([sel_main, sel_ctrl])[keep]
    res["draw"] = res["draw"].astype(int)
    res.to_csv(out / "w2_candidates.csv", index=False)
    draw_picture(t, main, ctrl, sel_main, sel_ctrl, out / "w2_candidates.png", a.dpi)
    log(f"\nWROTE {out / 'w2_candidates.csv'} and {out / 'w2_candidates.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
