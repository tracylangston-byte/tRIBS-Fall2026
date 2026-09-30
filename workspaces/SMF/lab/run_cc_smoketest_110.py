"""
run_cc_smoketest_110.py
========================
Series 110 -- Step 0 (smoke test) for the channel-loss (cc) identifiability
investigation described in Handoff_ChannelLossCalibration_v1.md.

WHAT THIS DOES
--------------
Builds and runs exactly ONE tRIBS run on the new (tribs6) build:

    Ks_mult          = 7.0    (truth)
    f_RS_abs         = 0.012  (truth)
    kinemvelcoef     = 4.5    (truth)
    flowexp          = 0.24   (truth)
    channelroughness = 0.026  (truth)
    optpercolation   = 1      (ON for the first time this project)
    channelconductivity_mmhr = 70   (BASELINE placeholder value -- NOT a
                                      calibrated choice, just what happened
                                      to already be sitting in BASELINE)

WHY
---
Per ways-of-working.md, any new anchor/configuration must be individually
smoke-tested via build+run+score before being trusted in a batch sweep.
`optpercolation=1` has never been exercised end-to-end on the tribs6 build
in this project -- Step 1 (the pytRIBS 1.0.0 migration) confirmed the build
itself is trustworthy, but only ever ran with optpercolation=0 pinned. This
script is that first exercise: it validates (a) the run doesn't crash or
hang with cc active, (b) the override plumbing in build_sensitivity_run_
tribs6.py's `overrides` dict genuinely reaches the model (checked explicitly
below by reading back optpercolation/channelconductivity_mmhr from the
scored metrics row), and (c) gives a real per-run wall-clock time to size
the Stage 1 (cc-only, ~18 runs) and Stage 2 (joint Ks/f/cc LHS, ~80 runs)
budgets against the available ~2-hour window.

Ks_mult/f_RS_abs/kinemvelcoef/flowexp/channelroughness are passed as
EXPLICIT overrides rather than relying on BASELINE in build_sensitivity_run_
tribs6.py, because that file's BASELINE dict currently holds stale
pre-Series-100 values (Ks_mult=6.1, f_RS_abs=0.020) that were never updated
after the truth reset. Explicit overrides here sidestep that regardless of
whether/when BASELINE gets cleaned up.

Usage (run from the lab/ directory, same as every other tribs6 script):
    python run_cc_smoketest_110.py
    python run_cc_smoketest_110.py --cc_value 200      # test a different cc
    python run_cc_smoketest_110.py --timeout 600        # more generous timeout

Output:
    Normal run artifacts land wherever build_sensitivity_run_tribs6.py's
    build_input_file() puts them (run_category "40_multivariable", since
    series 110 falls outside the 59-69 single-param block -- see that
    file's PARAM_CONFIG comment). The scored metrics row is
    calibration_work/03_comparisons/summary_tables/<run_id>_metrics_summary.csv,
    written by run_sensitivity_single_interp_tribs6.py exactly as for any
    other run -- this script does not duplicate that file, only reads it
    back afterward to report and verify.
"""

import argparse
import os
import sys
import signal
import subprocess
import time
from pathlib import Path

import pandas as pd

import build_sensitivity_run_tribs6 as builder

# ------------------------------------------------------------------
# Truth values (confirmed synthetic truth, Series 100/101) -- pinned
# explicitly rather than trusting BASELINE (see module docstring).
# ------------------------------------------------------------------
TRUTH_VALUES = {
    "Ks_mult":          7.0,
    "f_RS_abs":         0.012,
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}

SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"


# ------------------------------------------------------------------
# Timeout-safe execution -- ported from run_lhs_synth_Ks_f_100.py's
# run_with_timeout(), only the target script name changed. Runs the
# scorer as a separate subprocess in its own process group so a hang
# can be killed as a unit instead of blocking this script forever.
# ------------------------------------------------------------------
def run_with_timeout(timeout_sec):
    proc = subprocess.Popen(
        [sys.executable, SCORER_SCRIPT],
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
        description="Series 110 Step 0 -- single smoke-test run with "
                    "optpercolation=1/channelconductivity_mmhr active, "
                    "Ks/f/cv/r/n pinned at truth, on the tribs6 build.")
    parser.add_argument("--cc_value", type=float, default=70.0,
                        help="channelconductivity_mmhr to test (default: "
                             "70, the existing BASELINE placeholder -- not "
                             "a calibrated value)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Hard wall-clock timeout in seconds (default: "
                             "300, matching Series 100/101 convention)")
    args = parser.parse_args()

    print("=" * 70)
    print("Series 110 -- Step 0 smoke test (channel-loss / cc)")
    print(f"  Ks_mult={TRUTH_VALUES['Ks_mult']}  f_RS_abs={TRUTH_VALUES['f_RS_abs']}  "
          f"cv={TRUTH_VALUES['kinemvelcoef']}  r={TRUTH_VALUES['flowexp']}  "
          f"n={TRUTH_VALUES['channelroughness']}  (all pinned at truth)")
    print(f"  optpercolation=1 (OFF everywhere else in this project until now)")
    print(f"  channelconductivity_mmhr={args.cc_value}  (placeholder value, "
          f"not calibrated)")
    print(f"  timeout={args.timeout}s")
    print("=" * 70)

    # --- Build ---
    t_build0 = time.time()
    try:
        run_id, input_file, log_file = builder.build_input_file(
            "channelconductivity_mmhr",
            args.cc_value,
            overrides={
                "Ks_mult":          TRUTH_VALUES["Ks_mult"],
                "f_RS_abs":         TRUTH_VALUES["f_RS_abs"],
                "kinemvelcoef":     TRUTH_VALUES["kinemvelcoef"],
                "flowexp":          TRUTH_VALUES["flowexp"],
                "channelroughness": TRUTH_VALUES["channelroughness"],
                "optpercolation":   1,
            },
            tag="smoketest",
        )
    except Exception as e:
        print(f"\nBUILD FAILED before tRIBS was even launched: {e}")
        print("This is itself useful information -- it means the overrides "
              "mechanism or the build step has a problem independent of "
              "whether cc/optpercolation actually work at runtime.")
        sys.exit(1)
    build_elapsed = time.time() - t_build0
    print(f"\nBuilt {run_id} in {build_elapsed:.1f}s. Launching tRIBS run + scoring...\n")

    # --- Run + score (timeout-safe) ---
    returncode, run_elapsed, timed_out, tribs_warning = run_with_timeout(args.timeout)

    if timed_out:
        status = "HANG"
    elif returncode != 0 or tribs_warning:
        status = "FAILED"
    else:
        status = "SUCCESS"

    total_elapsed = build_elapsed + run_elapsed

    print("\n" + "=" * 70)
    print(f"RESULT: {status}")
    print(f"  Build time:      {build_elapsed:.1f}s")
    print(f"  Run+score time:  {run_elapsed:.1f}s")
    print(f"  TOTAL:           {total_elapsed:.1f}s "
          f"({total_elapsed/60:.2f} min)")
    print("=" * 70)

    if status != "SUCCESS":
        reason = ("wall-clock timeout" if timed_out else
                   "tRIBS reported non-zero exit" if tribs_warning else
                   f"{SCORER_SCRIPT} exited {returncode}")
        print(f"\n  Reason: {reason}")
        print(f"  Check the log for details: {log_file}")
        print("\n  This is a real finding, not just a script bug to fix -- it "
              "means cc/optpercolation may not be safe to sweep as freely as "
              "planned. Worth a look before writing the Stage 1/2 LHS script.")
        sys.exit(1)

    # --- Read back the scored metrics and verify the overrides really landed ---
    script_dir   = Path.cwd()
    project_root = script_dir.parent
    calib_dir    = project_root / "calibration_work"
    summary_file = (calib_dir / "03_comparisons" / "summary_tables"
                     / f"{run_id}_metrics_summary.csv")

    if not summary_file.exists():
        print(f"\n  WARNING: run reported SUCCESS but no metrics file found "
              f"at {summary_file}. Investigate before trusting this result.")
        sys.exit(1)

    m = pd.read_csv(summary_file).iloc[0]

    print(f"\n  KGE={m['kge']:.3f}  KGE_2012={m['kge_2012']:.3f}  "
          f"NSE={m['nse']:.3f}  PBIAS={m['pbias_pct']:+.1f}%")
    print(f"  Peak error: {m['peak_error_pct']:+.1f}%   "
          f"Volume error: {m['volume_error_pct']:+.1f}%   "
          f"Peak timing error: {m['peak_timing_error_hr']:+.2f} hr")

    # Defensive check: confirm optpercolation/channelconductivity_mmhr
    # echoed back from the actual run config match what was requested --
    # this is the core thing being validated by this smoke test.
    optperc_ok = int(m["optpercolation"]) == 1
    cc_ok      = abs(float(m["channelconductivity_mmhr"]) - args.cc_value) < 1e-6

    print(f"\n  Override verification:")
    print(f"    optpercolation echoed back as {int(m['optpercolation'])} "
          f"(expected 1) -- {'PASS' if optperc_ok else 'FAIL'}")
    print(f"    channelconductivity_mmhr echoed back as "
          f"{m['channelconductivity_mmhr']} (expected {args.cc_value}) -- "
          f"{'PASS' if cc_ok else 'FAIL'}")

    if not (optperc_ok and cc_ok):
        print("\n  At least one override did NOT reach the scored run. Do "
              "not proceed to Stage 1/2 until this is understood -- a "
              "silently-ignored override would mean every downstream cc "
              "run this session is actually scoring optpercolation=0.")
        sys.exit(1)

    print("\n  Smoke test PASSED: cc ran cleanly end-to-end and the "
          "override plumbing is confirmed working.")
    print(f"\n  Use total_elapsed={total_elapsed:.1f}s to size the Stage 1 "
          f"(~18 runs) and Stage 2 (~80 runs) budgets against the time "
          f"actually available this session.")


if __name__ == "__main__":
    main()
