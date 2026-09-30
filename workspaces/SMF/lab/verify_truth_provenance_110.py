"""
verify_truth_provenance_110.py
================================
Answers one question directly: is the truth file currently sitting in
calibration_work/synth_truth/ (as of this session, ks7_f012_Outlet.qout)
actually consistent with the CURRENT (tribs6) build -- or could it be a
leftover from the old, pre-migration tRIBS build?

This matters because the migration handoff (Handoff_ChannelLossCalibration
_v1.md) records that smoke test 3a -- a self-consistency check, run against
its own truth on the new build -- already passed (PBIAS=+0.0000%,
KGE(2012)=1.0000). But that was never directly confirmed against the
specific file this session has been using: when it was first uploaded, the
question of its provenance was asked and never actually answered, and Stage
1 (18 runs) and Stage 2 (400 runs, overnight) have both been scored against
it since. One of the four bugs found during the tRIBS 6.0.0/pytRIBS 1.0.0
port was a *.qout FORMAT CHANGE, and smoke test 3b separately established
a real, statistically robust deviation between old- and new-build output at
the same parameters (Spearman rho~0.8 for Ks_mult vs PBIAS/KGE). So "is
this truth file new-build or old-build" is not a paranoid question -- it's
the exact thing 3a was designed to rule out, and it's worth re-confirming
directly rather than trusting that it happened.

WHAT THIS SCRIPT DOES
-----------------------
Builds and runs ONE simulation at exactly the parameters the current truth
file is supposed to represent -- Ks_mult=7.0x, f_RS_abs=0.012,
kinemvelcoef=4.5, flowexp=0.24, channelroughness=0.026, optpercolation=0
(no channel loss -- the same "no cc" configuration every truth in this
project has always been defined under) -- on the CURRENTLY ACTIVE build.
It then scores that run using the exact same synthetic-truth-mode scorer
(run_sensitivity_single_interp_tribs6.py) already used everywhere else in
Series 110, which will automatically compare it against whatever *.qout
currently sits in calibration_work/synth_truth/.

INTERPRETING THE RESULT
-------------------------
If the current build reproduces the current truth file almost exactly
(PBIAS ~ 0.0000%, KGE_2012 ~ 1.0000, matching smoke test 3a's own
precedent), that's strong direct evidence the file IS new-build truth --
whatever its exact origin, the current (trusted, already-validated) build
regenerates it essentially bit-for-bit at the parameters it's supposed to
represent.

If PBIAS/KGE come back meaningfully off from a perfect match, that is
strong evidence the current truth file does NOT match what this build
actually produces at those parameters -- i.e. it is NOT safe to keep using,
and both Stage 1 and Stage 2 (the full overnight run) would need to be
rerun against a corrected truth file.

Either way, this run's own raw output IS a freshly-generated, indisputably
new-build truth candidate. This script prints exactly where it landed
(<run_id>_Outlet.qout under calibration_work/02_results/60_sensitivity/)
so it's available to promote into synth_truth/ if the check fails.

Usage (run from the lab/ directory):
    python verify_truth_provenance_110.py
    python verify_truth_provenance_110.py --timeout 300
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

TRUTH_VALUES = {
    "Ks_mult":          7.0,
    "f_RS_abs":         0.012,
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}

SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"

# Match smoke test 3a's own precedent for what counts as "reproduces truth"
PBIAS_TOLERANCE = 0.01   # percent
KGE_TOLERANCE   = 0.9999


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
        description="Verify the current synth_truth file is consistent "
                    "with the current (tribs6) build, by regenerating it "
                    "and checking self-consistency directly.")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Hard wall-clock timeout in seconds (default: 300)")
    args = parser.parse_args()

    script_dir   = Path.cwd()
    project_root = script_dir.parent
    calib_dir    = project_root / "calibration_work"

    # ------------------------------------------------------------------
    # Confirm there IS a truth file to check against -- this script is
    # meaningless in gauge mode.
    # ------------------------------------------------------------------
    synth_dir  = calib_dir / "synth_truth"
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            f"Expected exactly one *.qout file in {synth_dir}, found "
            f"{len(qout_files)}: {[f.name for f in qout_files]}. Nothing "
            f"to verify against."
        )
    truth_file = qout_files[0]
    print("=" * 70)
    print("Truth provenance check -- Series 110")
    print(f"  Regenerating: Ks_mult={TRUTH_VALUES['Ks_mult']}  "
          f"f_RS_abs={TRUTH_VALUES['f_RS_abs']}  "
          f"cv={TRUTH_VALUES['kinemvelcoef']}  r={TRUTH_VALUES['flowexp']}  "
          f"n={TRUTH_VALUES['channelroughness']}  optpercolation=0")
    print(f"  Scoring against CURRENT truth file: {truth_file.name}")
    print("=" * 70)

    # --- Build (Ks_mult as the nominal swept param, at its own truth
    #     value -- overrides pin everything else explicitly rather than
    #     trusting BASELINE, which is stale) ---
    try:
        run_id, input_file, log_file = builder.build_input_file(
            "Ks_mult",
            TRUTH_VALUES["Ks_mult"],
            overrides={
                "f_RS_abs":         TRUTH_VALUES["f_RS_abs"],
                "kinemvelcoef":     TRUTH_VALUES["kinemvelcoef"],
                "flowexp":          TRUTH_VALUES["flowexp"],
                "channelroughness": TRUTH_VALUES["channelroughness"],
                "optpercolation":   0,
            },
            tag="truthcheck",
        )
    except Exception as e:
        print(f"\nBUILD FAILED: {e}")
        sys.exit(1)

    print(f"\nBuilt {run_id}. Launching tRIBS run + scoring against "
          f"{truth_file.name}...\n")

    # --- Run + score ---
    returncode, elapsed, timed_out, tribs_warning = run_with_timeout(args.timeout)

    if timed_out:
        status = "HANG"
    elif returncode != 0 or tribs_warning:
        status = "FAILED"
    else:
        status = "SUCCESS"

    # Where this run's own raw output landed -- a fresh, indisputably
    # new-build truth candidate, regardless of how the check below turns out.
    raw_qout_candidate = (calib_dir / "02_results" / "60_sensitivity"
                           / run_id / f"{run_id}_Outlet.qout")

    print("\n" + "=" * 70)
    print(f"RESULT: {status}")
    print(f"  Elapsed: {elapsed/60:.2f} min")
    print("=" * 70)

    if status != "SUCCESS":
        reason = ("wall-clock timeout" if timed_out else
                   "tRIBS reported non-zero exit" if tribs_warning else
                   f"{SCORER_SCRIPT} exited {returncode}")
        print(f"\n  Could not complete the check: {reason}. See {log_file}.")
        sys.exit(1)

    summary_file = (calib_dir / "03_comparisons" / "summary_tables"
                    / f"{run_id}_metrics_summary.csv")
    if not summary_file.exists():
        print(f"\n  WARNING: run reported SUCCESS but no metrics file found "
              f"at {summary_file}.")
        sys.exit(1)

    m = pd.read_csv(summary_file).iloc[0]
    pbias    = float(m["pbias_pct"])
    kge      = float(m["kge"])
    kge_2012 = float(m["kge_2012"])
    obs_mode = m.get("obs_mode", "unknown")

    print(f"\n  obs_mode={obs_mode}  (must be 'synth' for this check to mean anything)")
    print(f"  PBIAS={pbias:+.4f}%   KGE={kge:.4f}   KGE_2012={kge_2012:.4f}")
    print(f"  (smoke test 3a's own precedent: PBIAS=+0.0000%, KGE_2012=1.0000)")

    consistent = (obs_mode == "synth" and abs(pbias) < PBIAS_TOLERANCE
                  and kge_2012 > KGE_TOLERANCE)

    print("\n" + "=" * 70)
    if consistent:
        print("VERDICT: CONSISTENT -- the current truth file IS reproduced")
        print("  by this build at the parameters it's supposed to represent.")
        print("  Whatever its exact origin, it behaves as new-build truth.")
        print("  Stage 1 and Stage 2 results scored against it stand.")
    else:
        print("VERDICT: NOT CONSISTENT -- the current truth file does NOT")
        print("  match what this build produces at Ks=7.0x/f=0.012/truth")
        print("  routing/optpercolation=0. Do not keep using it.")
        print(f"\n  This run's own output is a valid new-build truth candidate, at:")
        print(f"    {raw_qout_candidate}")
        print(f"  (confirm it exists, then replace {truth_file.name} in")
        print(f"   {synth_dir}/ with it -- remove the old file first so")
        print(f"   exactly one *.qout remains, per the hard safety checks")
        print(f"   in every Stage 1/Stage 2 script.)")
        print(f"\n  Stage 1 (18 runs) and Stage 2 (400 runs, last night's full")
        print(f"  overnight sweep) were both scored against the file now shown")
        print(f"  inconsistent, and would need to be rerun against the")
        print(f"  corrected truth before their results can be trusted.")
    print("=" * 70)


if __name__ == "__main__":
    main()
