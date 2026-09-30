"""
run_cc_stage1_sweep_110.py
============================
Series 110 -- Stage 1: cc-ONLY log-spaced sweep at fixed Ks/f/cv/r/n (truth
values), per Handoff_ChannelLossCalibration_v1.md section 2, stage 1.

Purpose: marginal sensitivity and direction of channelconductivity_mmhr,
and an empirically-grounded range to carry into Stage 2 (joint Ks/f/cc LHS).
Explicitly a first look, not a calibration result -- see the handoff for why
a cc-only sweep can't by itself answer the OPTPERCOLATION=0 confound
question (that's Stage 2's job).

REQUIRES SYNTHETIC TRUTH MODE. This sweep is only meaningful scored against
the same synthetic truth used throughout Series 100/101 (Ks_mult=7.0x,
f_RS_abs=0.012, generated with optpercolation=0) -- NOT against the real
gauge. The Step 0 smoke test (2026-09-24) found calibration_work/synth_truth/
does not currently contain the expected file on this (tribs6) build, so the
scorer silently fell back to real-gauge mode. Per Handoff section 0, a
new-build synthetic truth for this exact point was very likely already
generated during the Step 1 migration's smoke test 3a (self-consistency:
PBIAS=0.0000%, KGE=1.0000 "scored against its own truth") -- that file
needs to be the ONLY *.qout in calibration_work/synth_truth/ before this
script is run for real. The safety check below refuses to proceed
otherwise, exactly like run_lhs_synth_Ks_f_100.py did for Series 100.

Usage (run from the lab/ directory):
    python run_cc_stage1_sweep_110.py                    # 18 samples, seed=42
    python run_cc_stage1_sweep_110.py --n 12              # fewer samples
    python run_cc_stage1_sweep_110.py --skip_existing     # resume interrupted run
    python run_cc_stage1_sweep_110.py --timeout 300       # per-run hard timeout

Output:
    calibration_work/03_comparisons/summary_tables/lhs_results_cc_only_110.csv
    calibration_work/03_comparisons/summary_tables/lhs_results_cc_only_FAILED_110.csv
"""

import argparse
import os
import sys
import signal
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

import build_sensitivity_run_tribs6 as builder

# ------------------------------------------------------------------
# Truth values -- Ks/f/cv/r/n pinned throughout this sweep. Passed
# explicitly (not read from BASELINE, which is stale) on every build call.
# ------------------------------------------------------------------
TRUTH_VALUES = {
    "Ks_mult":          7.0,
    "f_RS_abs":         0.012,
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}

# ------------------------------------------------------------------
# cc sweep range -- provisional, log-spaced. See Handoff section 4 and
# Blasch et al. (2006): transient onset infiltration (~1000 m/d =
# ~41,700 mm/hr) represents only ~10-26% (avg ~18%) of total cumulative
# infiltration over an event, so this range deliberately stays well below
# that transient figure and is centered near the existing BASELINE
# placeholder (70 mm/hr), not an assumed-good value.
# ------------------------------------------------------------------
CC_LO, CC_HI = 3.0, 3000.0

LHS_SERIES   = "110"
LHS_CATEGORY = "110_cc_only_stage1"


# ------------------------------------------------------------------
# LOG-STRATIFIED SAMPLE GENERATION -- ported unchanged from
# run_lhs_synth_Ks_f_100.py's generate_lhs_samples(), specialized to one
# parameter here. Kept as a general function (not hardcoded to cc) since
# Stage 2's joint Ks/f/cc LHS will need the same mechanics.
# ------------------------------------------------------------------
def generate_lhs_samples(n, params, seed=None):
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


# ------------------------------------------------------------------
# Timeout-safe execution -- same pattern as run_cc_smoketest_110.py /
# run_lhs_synth_Ks_f_100.py. Matters more here than it did for Series 100:
# optpercolation=1 has never been exercised across a wide cc range before,
# so an unknown instability at some extreme is a real (if unlikely)
# possibility, not just boilerplate caution.
# ------------------------------------------------------------------
def run_with_timeout(timeout_sec):
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
        description="Series 110 Stage 1 -- cc-only log-spaced sweep at "
                    "truth Ks/f/cv/r/n, scored against synthetic truth.")
    parser.add_argument("--n", type=int, default=18,
                        help="Number of log-spaced cc samples (default: 18)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip samples whose compare CSV already exists")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Per-run hard timeout in seconds (default: 300)")
    args = parser.parse_args()

    script_dir   = Path.cwd()
    project_root = script_dir.parent
    calib_dir    = project_root / "calibration_work"
    summary_dir  = calib_dir / "03_comparisons" / "summary_tables"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # HARD SAFETY CHECK: refuse to run unless synthetic-truth mode will
    # actually activate. Same check run_lhs_synth_Ks_f_100.py uses for
    # Series 100 -- exactly one *.qout in calibration_work/synth_truth/.
    # This is the check the Step 0 smoke test's own run silently failed
    # (it has no --param sweep loop to protect, so it just ran in gauge
    # mode without complaint). This script does NOT get that pass.
    # ------------------------------------------------------------------
    synth_dir  = calib_dir / "synth_truth"
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            f"Expected exactly one *.qout file in {synth_dir} to activate "
            f"synthetic-truth mode, found {len(qout_files)}: "
            f"{[f.name for f in qout_files]}. The Step 0 smoke test ran in "
            f"GAUGE mode for this exact reason -- locate the new-build "
            f"synthetic truth .qout for Ks_mult=7.0x/f_RS_abs=0.012 "
            f"(likely already generated during the Step 1 migration's "
            f"smoke test 3a self-consistency check) and place it alone in "
            f"this directory before running Stage 1 for real."
        )
    print(f"Synthetic truth mode confirmed active: {qout_files[0].name}")
    if "Ks7p0x" not in qout_files[0].name:
        print(f"  NOTE: filename doesn't contain 'Ks7p0x' -- double-check "
              f"this is really the current truth (Ks_mult=7.0x/"
              f"f_RS_abs=0.012) before trusting results from this sweep.")

    samples = generate_lhs_samples(
        args.n, {"channelconductivity_mmhr": {"lo": CC_LO, "hi": CC_HI, "scale": "log"}},
        seed=args.seed,
    )

    out_path        = summary_dir / f"lhs_results_cc_only_{LHS_SERIES}.csv"
    failed_log_path = summary_dir / f"lhs_results_cc_only_FAILED_{LHS_SERIES}.csv"

    existing_df      = load_existing_results(out_path)
    existing_run_ids = (set(existing_df["run_id"].values)
                        if not existing_df.empty else set())

    results         = []
    failed_log_rows = []
    if not existing_df.empty:
        results.extend(existing_df.to_dict("records"))

    print(f"\n{'='*70}")
    print(f"Stage 1 -- Series {LHS_SERIES} -- cc-only sweep vs SYNTHETIC TRUTH")
    print(f"  ({args.n} samples, seed={args.seed}, timeout={args.timeout}s)")
    print(f"  channelconductivity_mmhr: {CC_LO} - {CC_HI} mm/hr  (log-stratified)")
    print(f"  PINNED (truth): Ks_mult={TRUTH_VALUES['Ks_mult']}  "
          f"f_RS_abs={TRUTH_VALUES['f_RS_abs']}  cv={TRUTH_VALUES['kinemvelcoef']}  "
          f"r={TRUTH_VALUES['flowexp']}  n={TRUTH_VALUES['channelroughness']}")
    print(f"  optpercolation=1 for every run in this sweep")
    print(f"{'='*70}\n")

    completed, skipped, hung, failed = 0, 0, 0, 0
    sweep_start = time.time()

    for i, row in samples.iterrows():
        cc_value = row["channelconductivity_mmhr"]

        try:
            run_id, _, log_file = builder.build_input_file(
                "channelconductivity_mmhr", cc_value,
                overrides={
                    "Ks_mult":          TRUTH_VALUES["Ks_mult"],
                    "f_RS_abs":         TRUTH_VALUES["f_RS_abs"],
                    "kinemvelcoef":     TRUTH_VALUES["kinemvelcoef"],
                    "flowexp":          TRUTH_VALUES["flowexp"],
                    "channelroughness": TRUTH_VALUES["channelroughness"],
                    "optpercolation":   1,
                },
                tag="s1",
            )
        except Exception as e:
            print(f"[{i+1:>3}/{args.n}]  BUILD FAILED (cc={cc_value:.2f}): {e}")
            failed += 1
            continue

        print(f"\n[{i+1:>3}/{args.n}]  cc={cc_value:.2f} mm/hr  -> {run_id}")

        if args.skip_existing and csv_already_exists(run_id, calib_dir):
            print(f"  SKIP (CSV exists): {run_id}")
            skipped += 1
            metrics_file = summary_dir / f"{run_id}_metrics_summary.csv"
            if metrics_file.exists() and run_id not in existing_run_ids:
                try:
                    m = pd.read_csv(metrics_file).iloc[0].to_dict()
                    m["channelconductivity_mmhr_sample"] = cc_value
                    results.append(m)
                except Exception:
                    pass
            continue

        t0 = time.time()
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
            metrics["channelconductivity_mmhr_sample"] = cc_value
            results = [r for r in results if r.get("run_id") != run_id]
            results.append(metrics)
            completed += 1
            print(f"  KGE_2012={metrics.get('kge_2012', float('nan')):.3f}  "
                  f"PBIAS={metrics.get('pbias_pct', float('nan')):+.1f}%")
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
                "channelconductivity_mmhr": cc_value,
            })
            pd.DataFrame(failed_log_rows).to_csv(failed_log_path, index=False)

        elapsed_total = time.time() - sweep_start
        remaining     = args.n - completed - skipped - hung - failed
        if completed > 0:
            avg_time = elapsed_total / completed
            eta_min  = (avg_time * remaining) / 60
            print(f"  Run time: {elapsed/60:.1f} min  |  ETA: {eta_min:.0f} min remaining")

        if results:
            pd.DataFrame(results).to_csv(out_path, index=False)

    print(f"\nStage 1 sweep complete: {completed} ran, {skipped} skipped, "
          f"{hung} hung, {failed} failed")

    if results:
        final_df = pd.DataFrame(results)
        if "kge_2012" in final_df.columns:
            final_df = final_df.sort_values("kge_2012", ascending=False)
        final_df.to_csv(out_path, index=False)
        print(f"Saved: {out_path.name}  ({len(final_df)} rows)")

        if "channelconductivity_mmhr_sample" in final_df.columns and "pbias_pct" in final_df.columns:
            corr = final_df["channelconductivity_mmhr_sample"].corr(final_df["pbias_pct"])
            print(f"\nPearson r (cc vs PBIAS): {corr:.3f}  "
                  f"-- compare to Ks_mult's r~-1.00 with PBIAS in Series 100 "
                  f"to gauge whether cc shows a similarly clean signature.")
    else:
        print("No results to save.")


if __name__ == "__main__":
    main()
