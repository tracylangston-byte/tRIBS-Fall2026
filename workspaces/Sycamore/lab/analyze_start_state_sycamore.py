#!/usr/bin/env python3
"""
analyze_start_state_sycamore.py  --  Stage 1, step 1: how many days does each start-state run
take to "forget" its different start and match the baseline hydrograph?

READ-ONLY: it reads results/*_Outlet.qout and writes one CSV and one PNG. It never runs tRIBS.

Agreed rule (Tracy + Claude, 2026-10-28)
    Take the daily-mean outlet flow of a run and of the baseline.
    A run has CONVERGED on the first day after which its daily mean stays within
        max(5 % of the baseline's daily mean, 0.13 m3/s)
    of the baseline's daily mean, on every remaining day through Jan 24.
    "Days to converge" = that first day counted from Oct 1 (day 0 = Oct 1).
    A run that never gets there is reported as "never": that is itself a finding.

Because the Jan 21 flood sits at the end of the window, the table also shows the same rule applied only
to the dry days BEFORE the flood week (days 0-107, Oct 1 to Jan 16). That separates "the start state is
still remembered" from "the flood itself came out a little different".

Run from workspaces/Sycamore/lab/:
    python analyze_start_state_sycamore.py
    python analyze_start_state_sycamore.py --baseline results/<other baseline>_Outlet.qout
"""
from __future__ import annotations

import argparse
import csv
import sys
from datetime import date, timedelta
from pathlib import Path

import numpy as np

HOURS = 2784                       # W3: Oct 1 2016 00:00 to Jan 25 2017 00:00
START = date(2016, 10, 1)
OCTDEC_HOURS = 92 * 24             # Oct 1 to Dec 31 (92 days)
FLOOD_START_DAY = 108              # Jan 17: the flood week is days 108-115 (Jan 17-24)
REL_TOL, ABS_TOL = 0.05, 0.13      # the agreed rule: 5 % of baseline or 0.13 m3/s, whichever is larger

BASELINE_DEFAULT = "results/Sycamore_octmar_ks0p2_cv6_psib0p5_nocl_Outlet.qout"
RUNS = [                           # (tag, what changed, run id)
    ("gw0p70", "water table 0.70 of bedrock", "W3_ks0p2_cv6_psib0p5_gw0p70_nocl"),
    ("gw0p82", "water table 0.82 of bedrock", "W3_ks0p2_cv6_psib0p5_gw0p82_nocl"),
    ("gw1p00", "water table 1.00 of bedrock", "W3_ks0p2_cv6_psib0p5_gw1p00_nocl"),
    ("bed0p0", "bedrock +0 m (and its water table)", "W3_ks0p2_cv6_psib0p5_bed0p0_nocl"),
    ("bed0p5", "bedrock +0.5 m (and its water table)", "W3_ks0p2_cv6_psib0p5_bed0p5_nocl"),
    ("bed2p0", "bedrock +2 m (and its water table)", "W3_ks0p2_cv6_psib0p5_bed2p0_nocl"),
]


def read_outlet(path: Path):
    """Return (time_hr, q) arrays from a tRIBS outlet file (comma-delimited, header, Time_hr,Q)."""
    if not path.is_file():
        sys.exit(f"File not found: {path}")
    t, q = [], []
    with open(path) as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        for row in reader:
            if len(row) >= 2:
                try:
                    t.append(float(row[0]))
                    q.append(float(row[1]))
                except ValueError:
                    continue
    t, q = np.array(t), np.array(q)
    if len(t) < 100:
        sys.exit(f"{path.name}: only {len(t)} readable rows (header was {header}).")
    if np.any(np.diff(t) < 0):
        sys.exit(f"{path.name}: Time_hr is not increasing, so I will not guess.")
    return t, q


def daily_means(t, q, hours=HOURS):
    """Mean flow for each whole day from Oct 1, using only times in [0, hours)."""
    n_days = hours // 24
    day = np.floor(t / 24.0).astype(int)
    out = np.full(n_days, np.nan)
    for d in range(n_days):
        sel = day == d
        if sel.any():
            out[d] = q[sel].mean()
    return out


def volume_mm3(t, q, h0, h1):
    """Flow volume (million m3) between hour h0 and hour h1 (trapezoid rule)."""
    sel = (t >= h0) & (t <= h1)
    if sel.sum() < 2:
        return float("nan")
    tt, qq = t[sel] * 3600.0, q[sel]
    return float(np.sum((qq[1:] + qq[:-1]) / 2.0 * np.diff(tt)) / 1e6)


def convergence_day(run_d, base_d, last_day):
    """First day after which the run stays within tolerance through last_day. None = never."""
    n = last_day + 1
    diff = np.abs(run_d[:n] - base_d[:n])
    tol = np.maximum(REL_TOL * base_d[:n], ABS_TOL)
    bad = np.where(~(diff < tol))[0]
    if bad.size == 0:
        return 0
    if bad[-1] == n - 1:
        return None
    return int(bad[-1] + 1)


def label_day(d):
    return "never" if d is None else f"{d}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", default=BASELINE_DEFAULT, help="outlet file of the baseline run")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out-csv", default="start_state_convergence.csv")
    ap.add_argument("--out-png", default="start_state_convergence.png")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    bt, bq = read_outlet(Path(args.baseline))
    print(f"Baseline file: {args.baseline}")
    print(f"  {len(bt)} rows, Time_hr from {bt[0]:g} to {bt[-1]:g}; about {len(bt) / max(bt[-1] - bt[0], 1):.1f} rows per hour")
    if bt[-1] < HOURS - 2:
        sys.exit(f"The baseline stops at hour {bt[-1]:g}; it needs to reach hour {HOURS}.")
    if bt[-1] > HOURS + 2:
        print(f"  Baseline runs to hour {bt[-1]:g}; only the first {HOURS} h (Oct 1 to Jan 25) are used.")
    base_d = daily_means(bt, bq)
    if np.isnan(base_d).any():
        sys.exit("The baseline has missing days.")
    n_days = len(base_d)
    flood_end = n_days - 1

    results = []
    daily = {}
    for tag, what, rid in RUNS:
        t, q = read_outlet(Path(args.results_dir) / f"{rid}_Outlet.qout")
        rd = daily_means(t, q)
        if np.isnan(rd).any():
            sys.exit(f"{rid}: missing days in the outlet file.")
        daily[tag] = rd
        c_all = convergence_day(rd, base_d, flood_end)
        c_pre = convergence_day(rd, base_d, FLOOD_START_DAY - 1)
        diff = np.abs(rd - base_d)
        tol = np.maximum(REL_TOL * base_d, ABS_TOL)
        outside = int((~(diff < tol)).sum())
        res = {
            "tag": tag, "changed": what, "run_id": rid,
            "converge_day_all": c_all, "converge_day_before_flood": c_pre,
            "days_outside_tolerance": outside,
            "max_abs_diff_m3s": round(float(diff.max()), 3),
            "octdec_vol_Mm3": round(volume_mm3(t, q, 0, OCTDEC_HOURS), 3),
            "flood_peak_m3s": round(float(q[(t >= FLOOD_START_DAY * 24) & (t <= HOURS)].max()), 1),
            "flood_vol_Mm3": round(volume_mm3(t, q, FLOOD_START_DAY * 24, HOURS), 3),
        }
        results.append(res)

    b_od = volume_mm3(bt, bq, 0, OCTDEC_HOURS)
    b_pk = float(bq[(bt >= FLOOD_START_DAY * 24) & (bt <= HOURS)].max())
    b_fv = volume_mm3(bt, bq, FLOOD_START_DAY * 24, HOURS)
    for r in results:
        r["octdec_vol_vs_base_pct"] = round(100 * (r["octdec_vol_Mm3"] / b_od - 1), 1) if b_od else float("nan")
        r["flood_peak_vs_base_pct"] = round(100 * (r["flood_peak_m3s"] / b_pk - 1), 1)
        r["flood_vol_vs_base_pct"] = round(100 * (r["flood_vol_Mm3"] / b_fv - 1), 1)

    def when(d):
        return "-" if d is None else (START + timedelta(days=d)).strftime("%b %d")

    print()
    print(f"Rule: converged = first day after which |daily mean - baseline daily mean| stays below "
          f"max({REL_TOL:.0%} of baseline, {ABS_TOL} m3/s) through Jan 24 (day {flood_end}).")
    print(f"Day 0 = Oct 1. 'before flood' = the same rule through day {FLOOD_START_DAY - 1} (Jan 16) only.\n")
    head = f"{'run':8s} {'days (all)':>10s} {'date':>7s} {'days (pre-flood)':>17s} {'date':>7s} {'days off':>9s} {'max diff':>9s}"
    print(head)
    print("-" * len(head))
    for r in results:
        print(f"{r['tag']:8s} {label_day(r['converge_day_all']):>10s} {when(r['converge_day_all']):>7s} "
              f"{label_day(r['converge_day_before_flood']):>17s} {when(r['converge_day_before_flood']):>7s} "
              f"{r['days_outside_tolerance']:>9d} {r['max_abs_diff_m3s']:>9.3f}")
    print(f"\nBaseline: Oct-Dec volume {b_od:.3f} Mm3, flood-week peak {b_pk:.1f} m3/s, flood-week volume {b_fv:.3f} Mm3")
    print(f"\n{'run':8s} {'Oct-Dec Mm3':>12s} {'vs base':>8s} {'flood peak':>11s} {'vs base':>8s} {'flood Mm3':>10s} {'vs base':>8s}")
    for r in results:
        print(f"{r['tag']:8s} {r['octdec_vol_Mm3']:>12.3f} {r['octdec_vol_vs_base_pct']:>7.1f}% {r['flood_peak_m3s']:>11.1f} "
              f"{r['flood_peak_vs_base_pct']:>7.1f}% {r['flood_vol_Mm3']:>10.3f} {r['flood_vol_vs_base_pct']:>7.1f}%")

    with open(args.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        for r in results:
            row = dict(r)
            row["converge_day_all"] = label_day(r["converge_day_all"])
            row["converge_day_before_flood"] = label_day(r["converge_day_before_flood"])
            w.writerow(row)
    print(f"\nWrote {args.out_csv}")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        days = np.arange(n_days)
        fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
        cols = plt.get_cmap("tab10")
        tol = np.maximum(REL_TOL * base_d, ABS_TOL)
        ax = axes[0]
        ax.semilogy(days, base_d, "k-", lw=2.5, label="baseline", zorder=5)
        for i, r in enumerate(results):
            ax.semilogy(days, np.maximum(daily[r["tag"]], 1e-4), color=cols(i), lw=1.2, label=r["tag"])
        ax.set_ylabel("daily mean outlet flow (m3/s)")
        ax.set_title("Start-state runs vs baseline, daily mean outlet flow")
        ax.legend(ncol=4, fontsize=9)
        ax = axes[1]
        for i, r in enumerate(results):
            ax.semilogy(days, np.maximum(np.abs(daily[r["tag"]] - base_d), 1e-5), color=cols(i), lw=1.2, label=r["tag"])
        ax.semilogy(days, tol, "k--", lw=2, label="tolerance")
        ax.axvline(FLOOD_START_DAY, color="gray", ls=":")
        ax.text(FLOOD_START_DAY + 0.5, ax.get_ylim()[0] * 1.5, "flood week", color="gray", fontsize=9)
        ax.set_ylabel("|run - baseline| (m3/s)")
        ax.set_xlabel("days since Oct 1 2016")
        ax.set_title("Difference from baseline (a run has converged once it stays under the dashed line)")
        fig.tight_layout()
        fig.savefig(args.out_png, dpi=130)
        print(f"Wrote {args.out_png}")


if __name__ == "__main__":
    main()
