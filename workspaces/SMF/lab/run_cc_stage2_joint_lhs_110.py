"""
run_cc_stage2_joint_lhs_110.py
================================
Series 110 -- Stage 2: joint Ks_mult / f_RS_abs / channelconductivity_mmhr
LHS, per Handoff_ChannelLossCalibration_v1.md section 2 and Stage 1's own
findings.

WHY THIS STEP EXISTS
---------------------
Every Ks_mult/f_RS_abs value calibrated so far in this project (the
current truth, the full identifiability hierarchy) was calibrated with
OPTPERCOLATION=0 -- a model with NO channel-loss process represented at
all. Stage 1 (cc-only sweep, fixed Ks/f at truth) showed two things that
make trusting that old calibration risky once cc is turned on:
  1. cc is cleanly identifiable via PBIAS (r=-0.994 across 18 points,
     3-3000 mm/hr) -- as clean a signature as Ks_mult's r~-1.00 in
     Series 100.
  2. At the point where cc alone brings PBIAS back to ~0% (cc~384 mm/hr),
     the hydrograph SHAPE is still wrong (rise_ratio=1.64 vs. the ideal
     1.0, peak_err=+8.8%) -- i.e. cc can fake a volume match while Ks/f
     are still off. That is direct, observed evidence of the
     Ks/f<->cc compensation risk the handoff flagged, not just a
     theoretical concern.
A cc-only sweep at FIXED Ks/f cannot resolve that -- it can only show cc's
own marginal effect at a possibly-already-biased Ks/f operating point.
This script lets Ks_mult, f_RS_abs, and channelconductivity_mmhr vary
TOGETHER so the analysis script (analyze_cc_stage2_110.py) can check
whether the recovered best-fit Ks/f drifts away from the current truth
(7.0x / 0.012) once cc is free to absorb part of the volume signal.

RANGES (widened from the original ~2-hour-budget proposal, now that this
is an overnight run and Stage 1 showed cc's effect has real teeth --
see conversation, 2026-09-24):
    Ks_mult:                   4.0   - 10.5   (linear)
    f_RS_abs:                  0.004 - 0.030  (log)
    channelconductivity_mmhr:  30    - 1000   (log)  <- Stage 1's own
                                                          plateau/peak/
                                                          early-decline
                                                          region (KGE_2012
                                                          peaked at
                                                          cc=291.69;
                                                          PBIAS crossed
                                                          zero near
                                                          cc=384)
Routing (kinemvelcoef/flowexp/channelroughness) held fixed at truth
throughout -- a DELIBERATE DEFERRAL per the handoff, not a resolution.
The pre-existing cv/r/n equifinality question (methods-and-findings.md)
stays open, independent of this step. optpercolation=1 for every run.

CAUTION ZONE: per the migration handoff, (Ks~9.83x, f~0.0485) is
untrustworthy ground (a genuine cross-version outlier, PBIAS=+311%,
KGE=-4.68, not yet explained by Josh). This script's f_RS_abs upper bound
(0.030) stays under that f, and results landing within ~10% of that Ks/f
corner are flagged in the output for manual skepticism rather than
silently trusted.

REQUIRES SYNTHETIC TRUTH MODE -- same hard check as Stage 1: refuses to
run unless exactly one *.qout sits in calibration_work/synth_truth/.

DISK MANAGEMENT: at n in the hundreds, this generates far more raw
per-run tRIBS output (pixel files, spatial snapshots) than anything
Series 110 has produced so far. Raw outputs for Series 96/97/99 were
PERMANENTLY LOST in the past to exactly this kind of accumulation
(methods-and-findings.md, "Data loss caution"). After each successful
run's metrics are safely read back and appended to the results CSV, this
script deletes that run's raw results directory
(calibration_work/02_results/110_channel_conductivity/<run_id>/) -- the
.in file and the two CSVs (compare_obs_sim, metrics_summary) are
untouched; only the bulky per-run spatial/pixel output goes. Failed/hung
runs are left alone (their raw output may be useful for debugging).

Usage (run from the lab/ directory, same convention as every other
tribs6 script):
    python run_cc_stage2_joint_lhs_110.py                  # 400 samples
    python run_cc_stage2_joint_lhs_110.py --n 200           # smaller run
    python run_cc_stage2_joint_lhs_110.py --skip_existing   # resume
    python run_cc_stage2_joint_lhs_110.py --timeout 300
    python run_cc_stage2_joint_lhs_110.py --no_cleanup      # keep raw dirs

Output:
    calibration_work/03_comparisons/summary_tables/lhs_results_joint_Ks_f_cc_110.csv
    calibration_work/03_comparisons/summary_tables/lhs_results_joint_Ks_f_cc_FAILED_110.csv

Downstream: analyze_cc_stage2_110.py reads the first CSV above and runs
the actual equifinality analysis (marginal correlations, cc-band peak
tracking, PCA on top performers, feasible-volume-fraction comparison
against Series 100's 2D result) -- this script only generates the data.
"""

import argparse
import os
import shutil
import sys
import signal
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

import build_sensitivity_run_tribs6 as builder

# ------------------------------------------------------------------
# Routing pinned at truth throughout -- deliberate deferral, not a
# resolution (see module docstring / Handoff_ChannelLossCalibration_v1.md
# section 2).
# ------------------------------------------------------------------
ROUTING_TRUTH = {
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}

TRUTH_KS = 7.0
TRUTH_F  = 0.012

LHS_PARAMS = {
    "Ks_mult":                  {"lo": 4.0,   "hi": 10.5,   "scale": "linear"},
    "f_RS_abs":                 {"lo": 0.004, "hi": 0.030,  "scale": "log"},
    "channelconductivity_mmhr": {"lo": 30.0,  "hi": 1000.0, "scale": "log"},
}

LHS_SERIES = "110"

# Untrustworthy cross-version-outlier corner (migration handoff section 4)
# -- results landing within this box get flagged, not excluded.
CAUTION_KS_LO, CAUTION_KS_HI = 8.85, 10.5   # within ~10% of Ks=9.83
CAUTION_F_LO,  CAUTION_F_HI  = 0.0437, 0.030  # within ~10% of f=0.0485 (capped at our own upper bound)


def generate_lhs_samples(n, params, seed=None):
    """Log- or linear-stratified LHS. Same mechanics as Stage 1 /
    Series 100 -- ported unchanged, now over 3 parameters at once
    (each parameter independently stratified into n bins and shuffled;
    this is stratified-independent sampling, not a true orthogonal LHS
    design, matching the precedent set by every prior sweep script in
    this project)."""
    rng     = np.random.default_rng(seed)
    samples = {}
    for param, bounds in params.items():
        lo, hi = bounds["lo"], bounds["hi"]
        scale  = bounds.get("scale", "linear")

        if scale == "log":
            log_lo, log_hi = np.log10(lo), np.log10(hi)
            intervals  = np.linspace(log_lo, log_hi, n + 1)
            log_points = rng.uniform(intervals[:-1], intervals[1:])
            points     = 10 ** log_points
        else:
            intervals = np.linspace(lo, hi, n + 1)
            points    = rng.uniform(intervals[:-1], intervals[1:])

        rng.shuffle(points)
        samples[param] = points
    return pd.DataFrame(samples)


def load_existing_results(out_path):
    if out_path.exists():
        try:
            df = pd.read_csv(out_path)
            print(f"  Loaded existing results: {len(df)} rows from {out_path.name}")
            return df
        except Exception as e:
            print(f"  Warning: could not load existing results ({e}). Starting fresh.")
    return pd.DataFrame()


def csv_already_exists(run_id, calib_dir):
    csv_path = (calib_dir / "03_comparisons" / "csv_exports"
                / f"{run_id}_compare_obs_sim.csv")
    return csv_path.exists()


def in_caution_zone(ks_val, f_val):
    return (CAUTION_KS_LO <= ks_val <= CAUTION_KS_HI and
            CAUTION_F_LO  <= f_val  <= CAUTION_F_HI)


def cleanup_raw_results(run_id, calib_dir):
    """Delete a successful run's raw results directory (pixel files,
    spatial snapshots) -- the .in file and the two CSVs (compare,
    metrics) live elsewhere and are untouched. Never raises: a cleanup
    failure shouldn't abort an overnight sweep over a non-critical
    housekeeping issue."""
    raw_dir = calib_dir / "02_results" / "110_channel_conductivity" / run_id
    try:
        if raw_dir.exists():
            shutil.rmtree(raw_dir)
            return True
    except Exception as e:
        print(f"  (cleanup warning: could not remove {raw_dir}: {e})")
    return False


def run_with_timeout(timeout_sec):
    """Timeout-safe subprocess execution -- ported unchanged from
    Stage 1 / Series 100. Own process group so a hang can be killed as a
    unit; output scanned for tRIBS's own failure warning since a clean
    Python exit code alone isn't sufficient evidence of a good run."""
    proc = subprocess.Popen(
        [sys.executable, "run_sensitivity_single_interp_tribs6.py"],
        cwd=Path.cwd(),
        start_new_session=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    t0 = time.time()
    try:
        stdout, _ = proc.communicate(timeout=timeout_sec)
        timed_out = False
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        stdout, _ = proc.communicate()
        timed_out = True
    elapsed = time.time() - t0

    if stdout:
        print(stdout, end="" if stdout.endswith("\n") else "\n")

    tribs_warning = bool(stdout) and "WARNING: tRIBS may have failed" in stdout
    returncode    = None if timed_out else proc.returncode
    return returncode, elapsed, timed_out, tribs_warning


def main():
    parser = argparse.ArgumentParser(
        description="Series 110 Stage 2 -- joint Ks_mult/f_RS_abs/"
                    "channelconductivity_mmhr LHS, vs synthetic truth.")
    parser.add_argument("--n", type=int, default=400,
                        help="Number of joint LHS samples (default: 400, "
                             "sized for an overnight run)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip samples whose compare CSV already exists "
                             "(resume an interrupted overnight run)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Per-run hard timeout in seconds (default: 300)")
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep raw per-run results directories instead "
                             "of deleting them after a successful run -- "
                             "uses much more disk over a 400-run sweep")
    args = parser.parse_args()

    script_dir   = Path.cwd()
    project_root = script_dir.parent
    calib_dir    = project_root / "calibration_work"
    summary_dir  = calib_dir / "03_comparisons" / "summary_tables"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # HARD SAFETY CHECK -- same as Stage 1 / Series 100. Refuses to
    # start rather than silently repeating the gauge-mode mistake that
    # cost a wasted Step 0 run earlier in this series.
    # ------------------------------------------------------------------
    synth_dir  = calib_dir / "synth_truth"
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            f"Expected exactly one *.qout file in {synth_dir} to activate "
            f"synthetic-truth mode, found {len(qout_files)}: "
            f"{[f.name for f in qout_files]}. Do not run Stage 2 in gauge "
            f"mode -- it would be scored against a different target than "
            f"the truth this whole series is built around."
        )
    print(f"Synthetic truth mode confirmed active: {qout_files[0].name}")

    samples = generate_lhs_samples(args.n, LHS_PARAMS, seed=args.seed)

    out_path        = summary_dir / f"lhs_results_joint_Ks_f_cc_{LHS_SERIES}.csv"
    failed_log_path = summary_dir / f"lhs_results_joint_Ks_f_cc_FAILED_{LHS_SERIES}.csv"

    existing_df = load_existing_results(out_path)
    existing_ids = set(existing_df["run_id"].values) if not existing_df.empty else set()
    results     = existing_df.to_dict("records") if not existing_df.empty else []

    failed_log_rows = []

    print(f"\n{'='*70}")
    print(f"Stage 2 -- Series {LHS_SERIES} -- joint Ks/f/cc LHS vs SYNTHETIC TRUTH")
    print(f"  ({args.n} samples, seed={args.seed}, timeout={args.timeout}s, "
          f"cleanup={'off' if args.no_cleanup else 'on'})")
    for p, b in LHS_PARAMS.items():
        print(f"  {p}: {b['lo']} - {b['hi']}  ({b['scale']})")
    print(f"  PINNED (truth routing): kinemvelcoef={ROUTING_TRUTH['kinemvelcoef']}  "
          f"flowexp={ROUTING_TRUTH['flowexp']}  "
          f"channelroughness={ROUTING_TRUTH['channelroughness']}")
    print(f"  optpercolation=1 for every run in this sweep")
    print(f"  Caution zone (flagged, not excluded): Ks {CAUTION_KS_LO}-{CAUTION_KS_HI}x, "
          f"f {CAUTION_F_LO}-{CAUTION_F_HI}")
    print(f"{'='*70}\n")

    completed, skipped, hung, failed, cleaned, caution_hits = 0, 0, 0, 0, 0, 0
    sweep_start = time.time()

    for i, row in samples.iterrows():
        ks_val = row["Ks_mult"]
        f_val  = row["f_RS_abs"]
        cc_val = row["channelconductivity_mmhr"]
        flagged = in_caution_zone(ks_val, f_val)

        try:
            run_id, _, log_file = builder.build_input_file(
                "channelconductivity_mmhr", cc_val,
                overrides={
                    "Ks_mult":          ks_val,
                    "f_RS_abs":         f_val,
                    "kinemvelcoef":     ROUTING_TRUTH["kinemvelcoef"],
                    "flowexp":          ROUTING_TRUTH["flowexp"],
                    "channelroughness": ROUTING_TRUTH["channelroughness"],
                    "optpercolation":   1,
                },
                tag="s2",
            )
        except Exception as e:
            print(f"[{i+1:>4}/{args.n}]  BUILD FAILED "
                  f"(Ks={ks_val:.2f}, f={f_val:.4f}, cc={cc_val:.1f}): {e}")
            failed += 1
            continue

        flag_note = "  [CAUTION ZONE]" if flagged else ""
        print(f"\n[{i+1:>4}/{args.n}]  Ks={ks_val:.3f}x  f={f_val:.4f}  "
              f"cc={cc_val:.1f} mm/hr  -> {run_id}{flag_note}")
        if flagged:
            caution_hits += 1

        if args.skip_existing and csv_already_exists(run_id, calib_dir):
            print(f"  SKIP (CSV exists): {run_id}")
            skipped += 1
            metrics_file = summary_dir / f"{run_id}_metrics_summary.csv"
            if metrics_file.exists() and run_id not in existing_ids:
                try:
                    m = pd.read_csv(metrics_file).iloc[0].to_dict()
                    m["caution_zone"] = flagged
                    results.append(m)
                    existing_ids.add(run_id)
                except Exception:
                    pass
            continue

        returncode, elapsed, timed_out, tribs_warning = run_with_timeout(args.timeout)

        if timed_out:
            status = "HANG"
        elif returncode != 0 or tribs_warning:
            status = "FAILED"
        else:
            status = "SUCCESS"

        if status == "SUCCESS":
            metrics_file = summary_dir / f"{run_id}_metrics_summary.csv"
            metrics = pd.read_csv(metrics_file).iloc[0].to_dict()
            metrics["caution_zone"] = flagged
            results = [r for r in results if r.get("run_id") != run_id]
            results.append(metrics)
            completed += 1
            print(f"  KGE_2012={metrics.get('kge_2012', float('nan')):.3f}  "
                  f"PBIAS={metrics.get('pbias_pct', float('nan')):+.1f}%")

            if not args.no_cleanup:
                if cleanup_raw_results(run_id, calib_dir):
                    cleaned += 1
        else:
            reason = ("wall-clock timeout" if timed_out else
                       "tRIBS reported non-zero exit" if tribs_warning else
                       f"scorer exited {returncode}")
            print(f"  {status}: {run_id}  ({reason})")
            if status == "HANG":
                hung += 1
            else:
                failed += 1
            failed_log_rows.append({
                "run_id": run_id, "status": status, "reason": reason,
                "elapsed_min": elapsed / 60,
                "Ks_mult": ks_val, "f_RS_abs": f_val,
                "channelconductivity_mmhr": cc_val,
                "caution_zone": flagged,
            })
            pd.DataFrame(failed_log_rows).to_csv(failed_log_path, index=False)

        elapsed_total = time.time() - sweep_start
        remaining     = args.n - completed - skipped - hung - failed
        if completed > 0:
            avg_time = elapsed_total / completed
            eta_min  = (avg_time * remaining) / 60
            print(f"  Run time: {elapsed/60:.1f} min  |  "
                  f"ETA: {eta_min:.0f} min ({eta_min/60:.1f} h) remaining")

        # Save incrementally -- an overnight run must survive being
        # interrupted at any point without losing completed work.
        if results:
            pd.DataFrame(results).to_csv(out_path, index=False)

    print(f"\nStage 2 sweep complete: {completed} ran, {skipped} skipped, "
          f"{hung} hung, {failed} failed, {cleaned} raw dirs cleaned, "
          f"{caution_hits} landed in the caution zone")

    if results:
        final_df = pd.DataFrame(results)
        if "kge_2012" in final_df.columns:
            final_df = final_df.sort_values("kge_2012", ascending=False)
        final_df.to_csv(out_path, index=False)
        print(f"Saved: {out_path.name}  ({len(final_df)} rows)")

        for p in ["Ks_mult", "f_RS_abs", "channelconductivity_mmhr"]:
            if p in final_df.columns:
                print(f"  {p} coverage: {final_df[p].min():.4f} - {final_df[p].max():.4f}")

        if "kge_2012" in final_df.columns and len(final_df) > 0:
            best = final_df.iloc[0]
            best_ks = best.get("Ks_mult")
            best_f  = best.get("f_RS_abs")
            print(f"\nBest run: Ks_mult={best_ks:.3f}  f_RS_abs={best_f:.4f}  "
                  f"cc={best.get('channelconductivity_mmhr'):.1f}  "
                  f"KGE_2012={best.get('kge_2012'):.3f}  "
                  f"PBIAS={best.get('pbias_pct'):+.1f}%"
                  f"{'  [CAUTION ZONE]' if best.get('caution_zone') else ''}")
            ks_drift = best_ks - TRUTH_KS
            f_drift  = best_f - TRUTH_F
            print(f"  Drift from truth ({TRUTH_KS}x / {TRUTH_F}): "
                  f"Ks_mult {ks_drift:+.3f}  f_RS_abs {f_drift:+.4f}")
            print(f"  Full equifinality analysis (marginal correlations, "
                  f"cc-band peak tracking, PCA, feasible-volume fraction): "
                  f"run analyze_cc_stage2_110.py against this CSV.")
    else:
        print("No results to save.")


if __name__ == "__main__":
    main()
