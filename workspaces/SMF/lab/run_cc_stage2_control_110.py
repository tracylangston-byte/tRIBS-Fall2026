"""
run_cc_stage2_control_110.py
============================
Series 110 -- Stage 2 CONTROL: the same joint Ks_mult / f_RS_abs /
channelconductivity_mmhr LHS as run_cc_stage2_joint_lhs_110.py, but with
channel loss switched OFF (optpercolation=0) on the current tRIBS 6.0.0 build.

WHY THIS SCRIPT EXISTS
-----------------------
The Stage 2 headline result ("feasible parameter space is 2.5x larger once cc
is free: 9.2% of samples within +/-2% PBIAS vs. 3.7% in the 2D world") compares
a tRIBS 6.0.0 result (Stage 2) against a baseline from the OLD build (Series
100). Smoke test 3b already showed the two builds deviate at the same
coordinates (Spearman rho ~0.8), and the Series 100 ranges are not necessarily
the Stage 2 ranges. This control removes both problems: same build, same
samples, same truth file -- the ONLY thing that differs from Stage 2 is
optpercolation (1 -> 0).

PAIRED DESIGN
--------------
Everything that determines the sampled points is identical to Stage 2:
  - same LHS_PARAMS (same ranges, same dict order -- order matters, it sets the
    order in which the random generator is consumed),
  - same generate_lhs_samples() (copied unchanged),
  - same seed (default 42) and same n (default 250 -- NOT 400).
So every point here is the exact (Ks_mult, f_RS_abs, cc) twin of a Stage 2
point, run with channel loss OFF. cc is still sampled and still written into
the .in file, but with OPTPERCOLATION=0 it should have no effect. That is
deliberate: it keeps the draws identical, and it gives a free check (see
"WHAT TO CHECK" below). The per-point difference in PBIAS between a Stage 2
run and its twin here is exactly the effect of turning cc on.

COLLISION AVOIDANCE (read before editing)
------------------------------------------
Run IDs, the results CSV, and the FAILED log are all renamed so nothing in
Stage 2 can be overwritten:
  run_id tag:   "_ctl"  (Stage 2 used "_s2")
  results CSV:  lhs_results_joint_Ks_f_cc_CONTROL_110.csv
  failed log:   lhs_results_joint_Ks_f_cc_CONTROL_FAILED_110.csv
The Stage 2 CSV is only ever READ (for the pairing check and end-of-run
comparison), never written.

WHAT TO CHECK
--------------
  1. Start of run: the PAIRING CHECK line should say all (or nearly all) of the
     n samples matched a Stage 2 point. If it says few or none matched, this
     run used a different --n or --seed than Stage 2 -- stop and fix before
     spending hours on it. (A few unmatched points are fine if some Stage 2 runs
     failed.)
  2. First few runs: each run prints its Stage 2 twin's KGE/PBIAS right under
     its own. With cc OFF, PBIAS should usually be HIGHER (less negative) than
     the twin's, because the twin loses water to channel infiltration.
  3. End of run: correlation of control PBIAS with log10(cc) should be ~0. If it
     is not, cc is NOT inert at OPTPERCOLATION=0 and the control is not a clean
     "cc off" baseline -- tell Josh.

NOT VALID FOR THIS CSV: analyze_cc_stage2_110.py. Its cc correlations, cc-band
tracking, and PCA assume cc has an effect; here it should not. The end-of-run
diagnostics below cover the control-vs-Stage-2 comparison.

Routing (kinemvelcoef/flowexp/channelroughness) pinned at truth, as in Stage 2.
optpercolation=0 for every run in this sweep.

REQUIRES SYNTHETIC TRUTH MODE -- same hard check as Stage 2: refuses to run
unless exactly one *.qout sits in calibration_work/synth_truth/.

DISK MANAGEMENT: same as Stage 2 -- each successful run's raw results
directory is deleted after its metrics are safely read back; failed/hung runs
are left alone. Use --no_cleanup to keep them.

DO NOT run any other tRIBS build/run script in lab/ while this is running --
all sweep scripts share calibration_work/current_run_config.json with no
locking.

Usage (run from the lab/ directory, same convention as every other script):
    python run_cc_stage2_control_110.py                   # n=250, seed=42
    python run_cc_stage2_control_110.py --skip_existing   # resume after interruption
    python run_cc_stage2_control_110.py --timeout 300
    python run_cc_stage2_control_110.py --no_cleanup

    # survives a browser/Codespace disconnect:
    nohup python run_cc_stage2_control_110.py > control_110.log 2>&1 &
    tail -f control_110.log

Output:
    calibration_work/03_comparisons/summary_tables/lhs_results_joint_Ks_f_cc_CONTROL_110.csv
    calibration_work/03_comparisons/summary_tables/lhs_results_joint_Ks_f_cc_CONTROL_FAILED_110.csv
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
# Routing pinned at truth throughout -- identical to Stage 2.
# ------------------------------------------------------------------
ROUTING_TRUTH = {
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}

TRUTH_KS = 7.0
TRUTH_F  = 0.012

# MUST match run_cc_stage2_joint_lhs_110.py exactly (ranges, scales, AND dict
# order) or the samples will not pair with Stage 2.
LHS_PARAMS = {
    "Ks_mult":                  {"lo": 4.0,   "hi": 10.5,   "scale": "linear"},
    "f_RS_abs":                 {"lo": 0.004, "hi": 0.030,  "scale": "log"},
    "channelconductivity_mmhr": {"lo": 30.0,  "hi": 1000.0, "scale": "log"},
}
PARAM_COLS = list(LHS_PARAMS.keys())

LHS_SERIES = "110"
OUT_TAG    = "CONTROL"   # in output filenames
RUN_TAG    = "ctl"       # in run_ids (Stage 2 used "s2")

STAGE2_CSV_NAME    = f"lhs_results_joint_Ks_f_cc_{LHS_SERIES}.csv"
FEASIBLE_PBIAS_PCT = 2.0   # |PBIAS| < 2% = "feasible", as in the Stage 2 analysis
MATCH_RTOL         = 1e-3  # tolerance when matching samples to Stage 2 rows

# Untrustworthy cross-version-outlier corner -- same flagging as Stage 2.
CAUTION_KS_LO, CAUTION_KS_HI = 8.85, 10.5
CAUTION_F_LO,  CAUTION_F_HI  = 0.0437, 0.030


def generate_lhs_samples(n, params, seed=None):
    """Log- or linear-stratified LHS. COPIED UNCHANGED from
    run_cc_stage2_joint_lhs_110.py so the draws are identical (same seed and
    n give the same points). Each parameter independently stratified into n
    bins and shuffled -- stratified-independent sampling, not a true
    orthogonal LHS design, matching every prior sweep script in this project.
    """
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
    """Delete a successful run's raw results directory (pixel files, spatial
    snapshots). The .in file and the two CSVs live elsewhere and are
    untouched. Never raises."""
    raw_dir = calib_dir / "02_results" / "110_channel_conductivity" / run_id
    try:
        if raw_dir.exists():
            shutil.rmtree(raw_dir)
            return True
    except Exception as e:
        print(f"  (cleanup warning: could not remove {raw_dir}: {e})")
    return False


def run_with_timeout(timeout_sec):
    """Timeout-safe subprocess execution -- unchanged from Stage 2. Own
    process group so a hang can be killed as a unit; output scanned for
    tRIBS's own failure warning."""
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


# ----------------------------------------------------------------------
# Pairing with Stage 2 (read-only)
# ----------------------------------------------------------------------
def load_stage2(stage2_path):
    """Read the Stage 2 results CSV (never written). Returns None if it is
    missing or lacks the parameter columns."""
    if not stage2_path.exists():
        return None
    try:
        df = pd.read_csv(stage2_path)
    except Exception as e:
        print(f"  Warning: could not read Stage 2 CSV ({e}).")
        return None
    if any(c not in df.columns for c in PARAM_COLS):
        print(f"  Warning: Stage 2 CSV lacks one of {PARAM_COLS}; cannot pair.")
        return None
    return df


def match_to_stage2(points_df, stage2_df, rtol=MATCH_RTOL):
    """For each row of points_df, return the positional index of the Stage 2
    row with the same (Ks_mult, f_RS_abs, cc) within rtol, or -1 if none."""
    out = np.full(len(points_df), -1, dtype=int)
    if stage2_df is None or any(c not in points_df.columns for c in PARAM_COLS):
        return out
    s2   = stage2_df[PARAM_COLS].to_numpy(dtype=float)
    pts  = points_df[PARAM_COLS].to_numpy(dtype=float)
    for i, row in enumerate(pts):
        ok   = np.all(np.isclose(s2, row, rtol=rtol, atol=0.0), axis=1)
        hits = np.flatnonzero(ok)
        if hits.size:
            out[i] = hits[0]
    return out


def report_control_diagnostics(final_df, stage2_df):
    """End-of-run checks. Wrapped by the caller in try/except -- a diagnostic
    failure must never lose the saved results."""
    print(f"\n{'='*70}")
    print("CONTROL DIAGNOSTICS (channel loss OFF)")
    print(f"{'='*70}")
    if "pbias_pct" not in final_df.columns:
        print("  No pbias_pct column in results; skipping diagnostics.")
        return
    pb = pd.to_numeric(final_df["pbias_pct"], errors="coerce").to_numpy(float)
    good = np.isfinite(pb)
    n_good = int(good.sum())
    print(f"  Scored runs: {n_good}")
    if n_good < 3:
        return

    # 1. cc should be inert
    if "channelconductivity_mmhr" in final_df.columns:
        logcc = np.log10(pd.to_numeric(
            final_df["channelconductivity_mmhr"], errors="coerce").to_numpy(float))
        both = good & np.isfinite(logcc)
        r = np.corrcoef(pb[both], logcc[both])[0, 1]
        verdict = ("OK -- cc looks inert, control is a clean cc-off baseline"
                   if abs(r) < 0.1 else
                   "NOT NEAR ZERO -- cc may NOT be inert at OPTPERCOLATION=0; tell Josh")
        print(f"  corr(PBIAS, log10 cc) = {r:+.3f}   -> {verdict}")

    # 2. feasible fraction
    n_feas = int((np.abs(pb[good]) < FEASIBLE_PBIAS_PCT).sum())
    print(f"  Feasible (|PBIAS| < {FEASIBLE_PBIAS_PCT:g}%): {n_feas}/{n_good} "
          f"= {100.0 * n_feas / n_good:.1f}%   "
          f"(Stage 2, cc ON: 9.2%;  Series 100, old build 2D: 3.7%)")

    # 3. paired comparison with Stage 2
    if stage2_df is None or "pbias_pct" not in stage2_df.columns:
        print("  (Stage 2 CSV not available -- no paired comparison.)")
        return
    tw = match_to_stage2(final_df, stage2_df)
    ok = (tw >= 0) & good
    n_ok = int(ok.sum())
    print(f"  Paired with a Stage 2 twin: {n_ok}/{n_good}")
    if n_ok == 0:
        return
    pb_ctl = pb[ok]
    pb_s2  = pd.to_numeric(stage2_df["pbias_pct"], errors="coerce"
                           ).to_numpy(float)[tw[ok]]
    delta  = pb_s2 - pb_ctl
    print(f"  Effect of turning cc ON (Stage 2 PBIAS minus control PBIAS): "
          f"mean {np.nanmean(delta):+.2f}, median {np.nanmedian(delta):+.2f} pct-pts")
    f_ctl = int((np.abs(pb_ctl) < FEASIBLE_PBIAS_PCT).sum())
    f_s2  = int((np.abs(pb_s2)  < FEASIBLE_PBIAS_PCT).sum())
    print(f"  Feasible on the SAME {n_ok} points:  cc OFF {f_ctl} "
          f"({100.0 * f_ctl / n_ok:.1f}%)   cc ON {f_s2} "
          f"({100.0 * f_s2 / n_ok:.1f}%)   ratio ON/OFF = "
          f"{(f_s2 / f_ctl) if f_ctl else float('nan'):.2f}")


def main():
    parser = argparse.ArgumentParser(
        description="Series 110 Stage 2 CONTROL -- joint Ks/f/cc LHS with "
                    "optpercolation=0, paired point-for-point with Stage 2.")
    parser.add_argument("--n", type=int, default=250,
                        help="Number of LHS samples (default: 250 -- must equal "
                             "the n used for Stage 2 or the points will not pair)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42 -- must equal Stage 2's seed)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip samples whose compare CSV already exists "
                             "(resume an interrupted run)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Per-run hard timeout in seconds (default: 300)")
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep raw per-run results directories instead of "
                             "deleting them after a successful run")
    args = parser.parse_args()

    script_dir   = Path.cwd()
    project_root = script_dir.parent
    calib_dir    = project_root / "calibration_work"
    summary_dir  = calib_dir / "03_comparisons" / "summary_tables"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # HARD SAFETY CHECK -- synthetic truth mode must be active.
    # ------------------------------------------------------------------
    synth_dir  = calib_dir / "synth_truth"
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            f"Expected exactly one *.qout file in {synth_dir} to activate "
            f"synthetic-truth mode, found {len(qout_files)}: "
            f"{[f.name for f in qout_files]}. Do not run the control in "
            f"gauge mode -- it would be scored against a different target "
            f"than Stage 2."
        )
    print(f"Synthetic truth mode confirmed active: {qout_files[0].name}")

    samples = generate_lhs_samples(args.n, LHS_PARAMS, seed=args.seed)

    out_path        = summary_dir / f"lhs_results_joint_Ks_f_cc_{OUT_TAG}_{LHS_SERIES}.csv"
    failed_log_path = summary_dir / f"lhs_results_joint_Ks_f_cc_{OUT_TAG}_FAILED_{LHS_SERIES}.csv"
    stage2_path     = summary_dir / STAGE2_CSV_NAME

    # Guard: the control must never write to Stage 2's files.
    assert out_path.name != STAGE2_CSV_NAME, "control output would overwrite Stage 2 CSV"

    # ------------------------------------------------------------------
    # PAIRING CHECK (read-only) -- do these samples match Stage 2's points?
    # ------------------------------------------------------------------
    stage2_df = load_stage2(stage2_path)
    twin_idx  = match_to_stage2(samples, stage2_df)
    if stage2_df is None:
        print("\nPAIRING CHECK: Stage 2 CSV not found/readable -- running UNPAIRED.")
    else:
        n_matched = int((twin_idx >= 0).sum())
        print(f"\nPAIRING CHECK: {n_matched}/{args.n} samples matched a Stage 2 "
              f"point (Stage 2 CSV has {len(stage2_df)} rows).")
        if n_matched < 0.9 * min(args.n, len(stage2_df)):
            print("  *** WARNING: poor pairing. This run probably used a different "
                  "--n or --seed than Stage 2. Press Ctrl+C now and fix it before "
                  "spending hours on an unpaired control. ***")

    existing_df  = load_existing_results(out_path)
    existing_ids = set(existing_df["run_id"].values) if not existing_df.empty else set()
    results      = existing_df.to_dict("records") if not existing_df.empty else []

    failed_log_rows = []

    print(f"\n{'='*70}")
    print(f"Stage 2 CONTROL -- Series {LHS_SERIES} -- joint Ks/f/cc LHS, "
          f"channel loss OFF, vs SYNTHETIC TRUTH")
    print(f"  ({args.n} samples, seed={args.seed}, timeout={args.timeout}s, "
          f"cleanup={'off' if args.no_cleanup else 'on'})")
    for p, b in LHS_PARAMS.items():
        print(f"  {p}: {b['lo']} - {b['hi']}  ({b['scale']})")
    print(f"  PINNED (truth routing): kinemvelcoef={ROUTING_TRUTH['kinemvelcoef']}  "
          f"flowexp={ROUTING_TRUTH['flowexp']}  "
          f"channelroughness={ROUTING_TRUTH['channelroughness']}")
    print(f"  optpercolation=0 for every run in this sweep "
          f"(cc is sampled but should be inert)")
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
                    "optpercolation":   0,          # <-- the ONLY substantive change vs Stage 2
                },
                tag=RUN_TAG,
            )
        except Exception as e:
            print(f"[{i+1:>4}/{args.n}]  BUILD FAILED "
                  f"(Ks={ks_val:.2f}, f={f_val:.4f}, cc={cc_val:.1f}): {e}")
            failed += 1
            continue

        flag_note = "  [CAUTION ZONE]" if flagged else ""
        print(f"\n[{i+1:>4}/{args.n}]  Ks={ks_val:.3f}x  f={f_val:.4f}  "
              f"cc={cc_val:.1f} mm/hr (inert)  -> {run_id}{flag_note}")
        if flagged:
            caution_hits += 1

        if args.skip_existing and csv_already_exists(run_id, calib_dir):
            print(f"  SKIP (CSV exists): {run_id}")
            skipped += 1
            metrics_file = summary_dir / f"{run_id}_metrics_summary.csv"
            if metrics_file.exists() and run_id not in existing_ids:
                try:
                    m = pd.read_csv(metrics_file).iloc[0].to_dict()
                    m.setdefault("Ks_mult", ks_val)
                    m.setdefault("f_RS_abs", f_val)
                    m.setdefault("channelconductivity_mmhr", cc_val)
                    m["optpercolation"] = 0
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
            # Make sure the three sampled parameters are in the row (needed for
            # pairing/diagnostics); keep the scorer's own value if it wrote one.
            metrics.setdefault("Ks_mult", ks_val)
            metrics.setdefault("f_RS_abs", f_val)
            metrics.setdefault("channelconductivity_mmhr", cc_val)
            metrics["optpercolation"] = 0
            metrics["caution_zone"] = flagged
            results = [r for r in results if r.get("run_id") != run_id]
            results.append(metrics)
            completed += 1
            print(f"  CONTROL (cc OFF): KGE_2012={metrics.get('kge_2012', float('nan')):.3f}  "
                  f"PBIAS={metrics.get('pbias_pct', float('nan')):+.1f}%")

            # Show the Stage 2 twin so the pairing is visible run by run.
            tw = twin_idx[i]
            if stage2_df is not None and tw >= 0:
                s2row = stage2_df.iloc[tw]
                print(f"  Stage 2 twin (cc ON):  KGE_2012={s2row.get('kge_2012', float('nan')):.3f}  "
                      f"PBIAS={s2row.get('pbias_pct', float('nan')):+.1f}%")
            elif stage2_df is not None:
                print("  (no Stage 2 twin found for this point)")

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

        # Save incrementally -- the run must survive being interrupted at any point.
        if results:
            pd.DataFrame(results).to_csv(out_path, index=False)

    print(f"\nControl sweep complete: {completed} ran, {skipped} skipped, "
          f"{hung} hung, {failed} failed, {cleaned} raw dirs cleaned, "
          f"{caution_hits} landed in the caution zone")

    if results:
        final_df = pd.DataFrame(results)
        if "kge_2012" in final_df.columns:
            final_df = final_df.sort_values("kge_2012", ascending=False)
        final_df.to_csv(out_path, index=False)
        print(f"Saved: {out_path.name}  ({len(final_df)} rows)")

        for p in PARAM_COLS:
            if p in final_df.columns:
                print(f"  {p} coverage: {final_df[p].min():.4f} - {final_df[p].max():.4f}")

        if "kge_2012" in final_df.columns and len(final_df) > 0:
            best = final_df.iloc[0]
            best_ks = best.get("Ks_mult")
            best_f  = best.get("f_RS_abs")
            print(f"\nBest control run: Ks_mult={best_ks:.3f}  f_RS_abs={best_f:.4f}  "
                  f"KGE_2012={best.get('kge_2012'):.3f}  "
                  f"PBIAS={best.get('pbias_pct'):+.1f}%"
                  f"{'  [CAUTION ZONE]' if best.get('caution_zone') else ''}")
            print(f"  Drift from truth ({TRUTH_KS}x / {TRUTH_F}): "
                  f"Ks_mult {best_ks - TRUTH_KS:+.3f}  f_RS_abs {best_f - TRUTH_F:+.4f}")

        try:
            report_control_diagnostics(final_df, stage2_df)
        except Exception as e:
            print(f"\n(diagnostics skipped: {e} -- results are saved regardless)")
    else:
        print("No results to save.")


if __name__ == "__main__":
    main()
