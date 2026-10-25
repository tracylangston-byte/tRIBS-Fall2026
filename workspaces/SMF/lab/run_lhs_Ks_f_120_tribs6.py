"""
run_lhs_Ks_f_120_tribs6.py
==========================
Series 120 -- Scenario 1 (the 1.0x storm leg) re-run on tRIBS 6.0.0 /
pytRIBS 1.0.0: a joint Ks_mult x f_RS_abs Latin-hypercube-style sweep scored
against the SYNTHETIC TRUTH (Ks_mult = 7.0x, f_RS_abs = 0.012).

WHAT THIS IS
------------
The AHS-deck "Scenario 1" experiment (old Series 100) was run on the old
tRIBS 5.3.1 build. This script repeats it on the current build, with the
same box, the same sampler, the same seed and the same n, so that every point
here has an exact twin in the old Series 100 results. The old/new difference
at each twin is then the effect of the build change and nothing else.

Settings (decided with Tracy, 2026-10-24):
  - Ks_mult   3.0 - 11.0, linear stratification        (truth 7.0)
  - f_RS_abs  0.003 - 0.05, LOG stratification         (truth 0.012)
  - n = 400, seed = 42   (identical to old Series 100, so draws pair)
  - Routing pinned at truth: kinemvelcoef 4.5, flowexp r 0.24,
    channelroughness 0.026
  - Channel loss OFF (optpercolation = 0). OPTSNOW = 0 is set inside the
    builder (build_sensitivity_run_tribs6.py), not by this script.
  - Storm: the original 12 Aug 2014 event (the 1.0x leg). The 0.8x and 1.25x
    legs come later, with their own truth files, using the SAME draws.

HOW THE RUNS ARE NAMED AND WHERE THEY GO
-----------------------------------------
build_sensitivity_run_tribs6.py has no entry for series 120, and this script
does not edit that file. Instead it teaches the builder about series 120 in
memory, for this one process only (see install_series_120_hooks()): the
builder is told that Ks_mult runs belong to series 120 and that series 120
goes in its own folder. Nothing on disk changes. A self-check runs at start-up
and stops the script if the naming does not come out as expected.

    run_id          SMF_20140812_120_Ks<Ks>x_f<f>        e.g. ..._120_Ks7p123456x_f0p012346
    raw results     calibration_work/02_results/120_lhs_synth_Ks_f/<run_id>/   (deleted after scoring)
    run inputs      calibration_work/01_run_inputs/120_lhs_synth_Ks_f/
    results CSV     calibration_work/03_comparisons/summary_tables/lhs_results_synth_Ks_f_120.csv
    failed log      calibration_work/03_comparisons/summary_tables/lhs_results_synth_Ks_f_FAILED_120.csv

If you change --n or --seed (e.g. for an extension batch), the output file
names gain a _n<N> and/or _seed<S> suffix so batches never overwrite each
other. The old Series 100 files are only ever READ.

THREE MODES
-----------
1. --dry_run      Generates the draws, runs the pairing check against the old
                  Series 100 CSV, prints coverage. Builds nothing, runs no
                  tRIBS, needs no truth file. Do this first; it takes seconds.
2. --smoke_test   Builds, runs and scores ONE run at the truth point
                  (Ks 7.0, f 0.012). Expect PBIAS +0.0000% and KGE_2012
                  1.0000. This proves the series-120 plumbing, the truth file
                  and the scorer all agree before you spend ~7 hours.
3. (default)      The full 400-run sweep.

SAFETY CHECKS BUILT IN
-----------------------
  * Refuses to run unless exactly ONE *.qout is in calibration_work/synth_truth/
    (otherwise the scorer silently falls back to the real gauge) AND its name
    contains "ks7p0x" (the run-id style used by the current truth file,
    SMF_20140812_60_Ks7p0x_truthcheck_Outlet.qout) or "ks7_f012" (the retired
    hand-made name). This script is the 1.0x leg only; override with
    --allow_other_truth if you know why. A name is only a hint -- the smoke
    test is what proves the file's contents.
  * Every scored run is READ BACK: the metrics row must say obs_mode=synth and
    must report the Ks, f, cv, flowexp r, n and optpercolation this script
    asked for. Any mismatch STOPS the sweep (it means the override plumbing is
    broken and the remaining runs would be wasted).
  * Each run is a separate subprocess with a hard timeout (default 300 s); the
    whole process group is killed on a hang. HANG/FAILED runs are logged to the
    FAILED csv and the sweep moves on. Ctrl+C also kills the in-flight run so
    no orphan tRIBS keeps writing current_run_config.json.
  * The raw results folder is deleted only if it sits exactly at
    02_results/120_lhs_synth_Ks_f/<run_id>.

DO NOT run any other tRIBS build/run script in lab/ while this is running --
all sweep scripts share calibration_work/current_run_config.json with no
locking.

USAGE (run from the lab/ directory, like every other script)
-------------------------------------------------------------
    python run_lhs_Ks_f_120_tribs6.py --dry_run
    python run_lhs_Ks_f_120_tribs6.py --smoke_test
    python run_lhs_Ks_f_120_tribs6.py --limit 5              # first 5 draws only
    python run_lhs_Ks_f_120_tribs6.py                        # full sweep
    python run_lhs_Ks_f_120_tribs6.py --skip_existing        # resume after interruption

    # survives a browser/Codespace disconnect:
    nohup python run_lhs_Ks_f_120_tribs6.py > series120_storm100.log 2>&1 &
    tail -f series120_storm100.log

WHAT TO LOOK AT AFTER THE SWEEP
--------------------------------
The end-of-run report prints: feasible count (|PBIAS| < 2%) against the
extension threshold of 20; a per-f-band table of how many runs fell on each
side of zero bias; any runs in the caution zone; and a paired comparison with
old Series 100. The best run is printed for reference only: on an equifinal
ridge the best run wanders along the ridge, so it is not the estimate.

STATUS (2026-10-24)
-------------------
Written without access to the Codespace. Tested only against stand-ins for the
builder and scorer: naming/hook logic, a full 400-draw pass, resume, hang and
failure handling, failed-log carry-forward, read-back abort, truth-file guards,
cleanup guard, Ctrl+C, and pairing against the real old Series 100 CSV (400/400).
NOT yet run against the real builder or real tRIBS. The --smoke_test is the first
real test. Update this paragraph once it has passed.

Revision (2026-10-24, same day): the first version of the truth-file name check
required "ks7_f012" in the file name. That was wrong -- the truth file in use is
SMF_20140812_60_Ks7p0x_truthcheck_Outlet.qout (made by verify_truth_provenance_110.py,
which replaced the old ks7_f012_Outlet.qout). Tracy's first real --smoke_test
stopped on that check. The check now accepts "ks7p0x" or "ks7_f012". The dry run
(400/400 draws paired with old Series 100, on the real Codespace) had already
passed.
"""

import argparse
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import build_sensitivity_run_tribs6 as builder

# ----------------------------------------------------------------------
# Series identity
# ----------------------------------------------------------------------
SERIES    = "120"
CATEGORY  = "120_lhs_synth_Ks_f"
STORM_LEG = "1.00x"

SCORER_SCRIPT = "run_sensitivity_single_interp_tribs6.py"

# ----------------------------------------------------------------------
# Truth and pinned values
# ----------------------------------------------------------------------
TRUTH_KS = 7.0
TRUTH_F  = 0.012
TRUTH_TOL = 1e-6
# A truth file's NAME can only ever be a hint (it cannot prove what is inside).
# These are the name fragments that mean "the Ks 7.0 / f 0.012 truth":
#   "ks7p0x"   -> the run-id style, e.g. SMF_20140812_60_Ks7p0x_truthcheck_Outlet.qout
#                 (this is the file verify_truth_provenance_110.py made and that is
#                 in synth_truth/ now)
#   "ks7_f012" -> the older hand-made name, ks7_f012_Outlet.qout (retired)
# The smoke test, not the name, is what proves the contents.
TRUTH_NAME_HINTS = ("ks7p0x", "ks7_f012")

ROUTING_TRUTH = {
    "kinemvelcoef":     4.5,
    "flowexp":          0.24,
    "channelroughness": 0.026,
}

# ----------------------------------------------------------------------
# Sampling box. MUST match old Series 100 exactly (ranges, scales, AND dict
# order -- order sets the order the random generator is consumed) or the draws
# will not pair with it.
# ----------------------------------------------------------------------
LHS_PARAMS = {
    "Ks_mult":  {"lo": 3.0,   "hi": 11.0, "scale": "linear"},
    "f_RS_abs": {"lo": 0.003, "hi": 0.05, "scale": "log"},
}
PARAM_COLS   = list(LHS_PARAMS.keys())
DEFAULT_N    = 400
DEFAULT_SEED = 42

# ----------------------------------------------------------------------
# Analysis thresholds
# ----------------------------------------------------------------------
FEASIBLE_PBIAS_PCT     = 2.0    # |PBIAS| < 2% = "feasible"
MIN_FEASIBLE_PER_STORM = 20     # fewer than this -> add a new-seed batch
N_F_BANDS              = 6      # log-spaced f bands for the bracket check
MATCH_RTOL             = 1e-3   # tolerance when matching draws to old Series 100

# Corner where the old-build vs new-build cross-version outlier sat
# (Ks 9.83 / f 0.0485: PBIAS +311%, KGE -4.68 on the old build, unexplained).
# Flagged, never excluded. Upper edges are the box edges.
CAUTION_KS_LO = 8.85
CAUTION_F_LO  = 0.0437

OLD_CSV_CANDIDATES = [
    "lhs_results_synth_Ks_f_100_RESCORED.csv",
    "lhs_results_synth_Ks_f_100.csv",
]

FAILED_COLS = ["run_id", "status", "reason", "elapsed_min", "draw_idx", "Ks_mult", "f_RS_abs"]

# Smoke test pass criteria (same precedent as verify_truth_provenance_110.py)
SMOKE_PBIAS_TOL = 0.01     # percent
SMOKE_KGE_MIN   = 0.9999


# ======================================================================
# Sampler -- COPIED UNCHANGED from run_cc_stage2_control_110.py /
# run_lhs_synth_Ks_f_100.py so the draws are identical.
# ======================================================================
def generate_lhs_samples(n, params, seed=None):
    """Log- or linear-stratified sampling. Each parameter is split into n equal
    slices (log slices for log-scale parameters), one uniform draw per slice,
    then shuffled independently. Stratified-independent sampling, not a true
    orthogonal LHS design -- matches every prior sweep script in this project.
    Same seed + same n + same bounds + same dict order = identical draws."""
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


# ======================================================================
# Teaching the (unedited) builder about series 120
# ======================================================================
_HOOKS_INSTALLED = False


def install_series_120_hooks():
    """In-memory only: tell the builder that Ks_mult runs are series 120 and
    that series 120 lives in its own category folder. The builder file on disk
    is not touched. Idempotent. Stops the script if the result is not exactly
    what the rest of this script assumes."""
    global _HOOKS_INSTALLED
    if _HOOKS_INSTALLED:
        return

    cfg = dict(builder.PARAM_CONFIG["Ks_mult"])
    cfg["series"] = SERIES
    builder.PARAM_CONFIG["Ks_mult"] = cfg

    original_get_run_category = builder.get_run_category

    def get_run_category_120(series_str):
        if str(series_str) == SERIES:
            return CATEGORY
        return original_get_run_category(series_str)

    builder.get_run_category = get_run_category_120

    # --- self-check ---
    run_id, _ = builder.build_run_id("Ks_mult", 7.0)
    expected  = f"{builder.LOCATION}_{builder.EVENT_DATE}_{SERIES}_Ks7p0x"
    if run_id != expected:
        raise RuntimeError(
            f"Series-120 hook self-check failed: build_run_id gave {run_id!r}, "
            f"expected {expected!r}. The builder's naming has changed; stop and "
            f"look at build_sensitivity_run_tribs6.py before running anything.")
    cat = builder.get_run_category(builder.PARAM_CONFIG["Ks_mult"]["series"])
    if cat != CATEGORY:
        raise RuntimeError(
            f"Series-120 hook self-check failed: category {cat!r}, "
            f"expected {CATEGORY!r}.")
    _HOOKS_INSTALLED = True


# ======================================================================
# Small helpers
# ======================================================================
def is_truth_point(ks, f):
    return abs(ks - TRUTH_KS) < TRUTH_TOL and abs(f - TRUTH_F) < TRUTH_TOL


def in_caution_zone(ks, f):
    return ks >= CAUTION_KS_LO and f >= CAUTION_F_LO


def variant_suffix(n, seed):
    s = ""
    if n != DEFAULT_N:
        s += f"_n{n}"
    if seed != DEFAULT_SEED:
        s += f"_seed{seed}"
    return s


def sweep_tag(f, seed):
    """Suffix appended to the builder's run_id. Carries f (the builder's own
    run_id only encodes Ks) and, for non-default seeds, the seed."""
    tag = f"f{builder.value_to_label(f)}"
    if seed != DEFAULT_SEED:
        tag += f"_sd{seed}"
    return tag


def predicted_run_id(ks, tag):
    run_id, _ = builder.build_run_id("Ks_mult", ks)
    return f"{run_id}_{tag}" if tag else run_id


def pinned_overrides(ks, f):
    """Everything that must be pinned for a Series 120 run."""
    return {
        "Ks_mult":          ks,
        "f_RS_abs":         f,
        "kinemvelcoef":     ROUTING_TRUTH["kinemvelcoef"],
        "flowexp":          ROUTING_TRUTH["flowexp"],
        "channelroughness": ROUTING_TRUTH["channelroughness"],
        "optpercolation":   0,
    }


def build_one(ks, f, tag):
    return builder.build_input_file(
        "Ks_mult", ks, overrides=pinned_overrides(ks, f), tag=tag)


def csv_already_exists(run_id, calib_dir):
    return (calib_dir / "03_comparisons" / "csv_exports"
            / f"{run_id}_compare_obs_sim.csv").exists()


def load_existing_results(out_path):
    if out_path.exists():
        try:
            df = pd.read_csv(out_path)
            print(f"  Loaded existing results: {len(df)} rows from {out_path.name}")
            return df
        except Exception as e:
            print(f"  Warning: could not load existing results ({e}). Starting fresh.")
    return pd.DataFrame()


def cleanup_raw_results(run_id, calib_dir):
    """Delete a successful run's raw results folder (pixel files, spatial
    snapshots). Only deletes a folder located EXACTLY at
    02_results/<CATEGORY>/<run_id>. Never raises."""
    try:
        parent = (calib_dir / "02_results" / CATEGORY).resolve()
        raw    = (parent / run_id).resolve()
        if raw.parent != parent or raw.name != run_id:
            print(f"  (cleanup skipped: {raw} is not directly inside {parent})")
            return False
        if raw.exists():
            shutil.rmtree(raw)
            return True
    except Exception as e:
        print(f"  (cleanup warning for {run_id}: {e})")
    return False


def upsert_failed(rows, new_row):
    """Add a FAILED-log record, replacing any earlier record for the same run_id
    (a draw that hangs again on resume must not appear twice)."""
    return [r for r in rows if r.get("run_id") != new_row["run_id"]] + [new_row]


def _kill_group(proc):
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


def run_with_timeout(timeout_sec):
    """Run the scorer in its own process group so a hang can be killed as a
    unit. Returns (returncode, elapsed_sec, timed_out, tribs_warning)."""
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
        _kill_group(proc)
        stdout, _ = proc.communicate()
        timed_out = True
    except KeyboardInterrupt:
        _kill_group(proc)
        proc.wait()
        raise
    elapsed = time.time() - t0

    if stdout:
        print(stdout, end="" if stdout.endswith("\n") else "\n")

    tribs_warning = bool(stdout) and "WARNING: tRIBS may have failed" in stdout
    returncode    = None if timed_out else proc.returncode
    return returncode, elapsed, timed_out, tribs_warning


def _num(m, key):
    try:
        return float(m.get(key, np.nan))
    except (TypeError, ValueError):
        return np.nan


def readback_problems(m, ks, f):
    """Compare the scored metrics row against what this script asked for.
    Returns a list of human-readable problems (empty list = all good)."""
    probs = []
    if str(m.get("obs_mode")) != "synth":
        probs.append(f"obs_mode={m.get('obs_mode')!r}, expected 'synth' "
                     f"(scorer fell back to the real gauge)")
    expected = {
        "Ks_mult":          ks,
        "f_RS_abs":         f,
        "kinemvelcoef":     ROUTING_TRUTH["kinemvelcoef"],
        "flowexp":          ROUTING_TRUTH["flowexp"],
        "channelroughness": ROUTING_TRUTH["channelroughness"],
        "optpercolation":   0,
    }
    for key, want in expected.items():
        got = _num(m, key)
        if not np.isfinite(got) or not np.isclose(got, want, rtol=1e-6, atol=1e-12):
            probs.append(f"{key}: scored row says {m.get(key)!r}, asked for {want!r}")
    return probs


# ======================================================================
# Pairing with old Series 100 (read-only)
# ======================================================================
def find_old_series100(summary_dir, explicit=None):
    candidates = [Path(explicit)] if explicit else [summary_dir / c for c in OLD_CSV_CANDIDATES]
    for p in candidates:
        if p.exists():
            try:
                df = pd.read_csv(p)
            except Exception as e:
                print(f"  Warning: could not read {p.name} ({e}).")
                continue
            if any(c not in df.columns for c in PARAM_COLS):
                print(f"  Warning: {p.name} lacks {PARAM_COLS}; cannot pair.")
                continue
            return df, p
    return None, None


def match_to_old(points_df, old_df, rtol=MATCH_RTOL):
    """For each draw, positional index of the old-Series-100 row with the same
    (Ks_mult, f_RS_abs) within rtol, or -1."""
    out = np.full(len(points_df), -1, dtype=int)
    if old_df is None:
        return out
    old = old_df[PARAM_COLS].to_numpy(dtype=float)
    pts = points_df[PARAM_COLS].to_numpy(dtype=float)
    for i, row in enumerate(pts):
        ok   = np.all(np.isclose(old, row, rtol=rtol, atol=0.0), axis=1)
        hits = np.flatnonzero(ok)
        if hits.size:
            out[i] = hits[0]
    return out


def print_pairing_check(samples, old_df, old_path, n, seed):
    if old_df is None:
        print("\nPAIRING CHECK vs old Series 100 (tRIBS 5.3.1): skipped -- old CSV not "
              f"found in summary_tables/ (looked for {OLD_CSV_CANDIDATES}). "
              "Pass --old_csv PATH to point at it. The sweep itself is unaffected.")
        return np.full(len(samples), -1, dtype=int)
    twin = match_to_old(samples, old_df)
    n_matched = int((twin >= 0).sum())
    print(f"\nPAIRING CHECK vs old Series 100 ({old_path.name}, {len(old_df)} rows): "
          f"{n_matched}/{len(samples)} draws matched a twin.")
    if n_matched < 0.9 * min(len(samples), len(old_df)):
        if n == DEFAULT_N and seed == DEFAULT_SEED:
            print("  *** WARNING: poor pairing with the default n and seed. The box, the "
                  "sampler, or the dict order in this script differs from old Series 100. "
                  "Stop and find out why before spending hours. ***")
        else:
            print("  (Expected: n or seed differs from old Series 100, so the draws are "
                  "new points, not twins.)")
    return twin


# ======================================================================
# End-of-run report
# ======================================================================
def _rank(a):
    return pd.Series(a).rank().to_numpy(dtype=float)


def report_final(final_df, old_df, failed_rows, n_planned, n_full):
    print(f"\n{'=' * 70}")
    print(f"SERIES {SERIES} REPORT (storm leg {STORM_LEG}, channel loss OFF)")
    print(f"{'=' * 70}")
    if "pbias_pct" not in final_df.columns:
        print("  No pbias_pct column in results; skipping report.")
        return
    pb = pd.to_numeric(final_df["pbias_pct"], errors="coerce").to_numpy(float)
    good = np.isfinite(pb)
    n_good = int(good.sum())
    print(f"  Scored runs: {n_good} of {n_planned} planned; "
          f"{len(failed_rows)} HANG/FAILED/mismatch record(s) in the FAILED log (see FAILED csv).")
    if n_good < 3:
        return

    for p in PARAM_COLS:
        v = pd.to_numeric(final_df[p], errors="coerce")
        print(f"  {p} coverage: {v.min():.4f} - {v.max():.4f}")

    # --- feasible set vs the extension threshold ---
    n_feas = int((np.abs(pb[good]) < FEASIBLE_PBIAS_PCT).sum())
    print(f"\n  Feasible (|PBIAS| < {FEASIBLE_PBIAS_PCT:g}%): {n_feas}/{n_good} "
          f"= {100.0 * n_feas / n_good:.1f}%")
    if n_good < 0.9 * n_full:
        print(f"  -> partial sweep ({n_good} of {n_full} draws scored): the extension "
              f"threshold of {MIN_FEASIBLE_PER_STORM} only applies to a complete sweep.")
    elif n_feas < MIN_FEASIBLE_PER_STORM:
        print(f"  -> BELOW the extension threshold of {MIN_FEASIBLE_PER_STORM}. "
              f"Add a new-seed batch (never a larger n).")
    else:
        print(f"  -> meets the extension threshold of {MIN_FEASIBLE_PER_STORM}.")

    # --- per-f-band bracket check (informational) ---
    f_all = pd.to_numeric(final_df["f_RS_abs"], errors="coerce").to_numpy(float)
    lo, hi = LHS_PARAMS["f_RS_abs"]["lo"], LHS_PARAMS["f_RS_abs"]["hi"]
    edges = np.logspace(np.log10(lo), np.log10(hi), N_F_BANDS + 1)
    print(f"\n  f-band check (PBIAS = 100*sum(sim-obs)/sum(obs); <0 under-predicts, >0 over-predicts)")
    print(f"  {'f band':>17}  {'runs':>5}  {'PBIAS<0':>8}  {'PBIAS>0':>8}  {'feasible':>8}  note")
    for b in range(N_F_BANDS):
        last = (b == N_F_BANDS - 1)
        m = good & (f_all >= edges[b]) & ((f_all <= edges[b + 1]) if last else (f_all < edges[b + 1]))
        n_b   = int(m.sum())
        n_neg = int((pb[m] < 0).sum())
        n_pos = int((pb[m] > 0).sum())
        n_fe  = int((np.abs(pb[m]) < FEASIBLE_PBIAS_PCT).sum())
        note = ""
        if n_b == 0:
            note = "EMPTY"
        elif n_neg == 0 or n_pos == 0:
            note = "no bracket (all one side of zero bias)"
        print(f"  {edges[b]:7.4f}-{edges[b + 1]:7.4f}  {n_b:5d}  {n_neg:8d}  {n_pos:8d}  {n_fe:8d}  {note}")
    print("  (A band with no bracket can mean the ridge leaves the Ks box there -- "
          "a finding to look at, not necessarily a sampling gap.)")

    # --- best run, for reference only ---
    if "kge_2012" in final_df.columns:
        best = final_df.sort_values("kge_2012", ascending=False).iloc[0]
        print(f"\n  Best run (reference only -- the best run wanders along an equifinal ridge, "
              f"it is not the estimate):")
        print(f"    Ks_mult={best['Ks_mult']:.3f}  f_RS_abs={best['f_RS_abs']:.4f}  "
              f"KGE_2012={best['kge_2012']:.3f}  PBIAS={best['pbias_pct']:+.2f}%")

    # --- caution zone ---
    if "caution_zone" in final_df.columns:
        cz = final_df[final_df["caution_zone"].astype(str) == "True"]
        print(f"\n  Caution-zone runs (Ks >= {CAUTION_KS_LO}, f >= {CAUTION_F_LO}): {len(cz)}")
        for _, r in cz.head(15).iterrows():
            print(f"    Ks={r['Ks_mult']:.3f}  f={r['f_RS_abs']:.4f}  "
                  f"KGE_2012={r.get('kge_2012', float('nan')):.3f}  "
                  f"PBIAS={r.get('pbias_pct', float('nan')):+.1f}%")

    # --- paired comparison with old Series 100 ---
    if old_df is None or "pbias_pct" not in old_df.columns:
        print("\n  (old Series 100 CSV not available -- no paired comparison.)")
        return
    twin = match_to_old(final_df, old_df)
    ok = (twin >= 0) & good
    n_ok = int(ok.sum())
    print(f"\n  PAIRED WITH OLD SERIES 100 (tRIBS 5.3.1): {n_ok}/{n_good} points have a twin")
    if n_ok < 3:
        return
    pb_new = pb[ok]
    pb_old = pd.to_numeric(old_df["pbias_pct"], errors="coerce").to_numpy(float)[twin[ok]]
    fin = np.isfinite(pb_old)
    pb_new, pb_old = pb_new[fin], pb_old[fin]
    delta = pb_new - pb_old
    rho = np.corrcoef(_rank(pb_new), _rank(pb_old))[0, 1]
    f_new = int((np.abs(pb_new) < FEASIBLE_PBIAS_PCT).sum())
    f_old = int((np.abs(pb_old) < FEASIBLE_PBIAS_PCT).sum())
    print(f"    PBIAS shift (new minus old): median {np.median(delta):+.2f}, "
          f"mean {np.mean(delta):+.2f} pct-pts")
    print(f"    Rank correlation of PBIAS, new vs old: {rho:.3f}  (cross-version check earlier gave ~0.8)")
    print(f"    Feasible on the same {len(pb_new)} points:  new build {f_new}   old build {f_old}")


# ======================================================================
# Smoke test
# ======================================================================
def run_smoke_test(args, calib_dir, summary_dir):
    ks, f = TRUTH_KS, TRUTH_F
    tag = "smoke"
    run_id = predicted_run_id(ks, tag)
    print(f"\n{'=' * 70}")
    print(f"SMOKE TEST -- Series {SERIES} plumbing at the truth point")
    print(f"  Ks_mult={ks}  f_RS_abs={f}  cv={ROUTING_TRUTH['kinemvelcoef']}  "
          f"flowexp r={ROUTING_TRUTH['flowexp']}  n={ROUTING_TRUTH['channelroughness']}  "
          f"optpercolation=0")
    print(f"  Expected: PBIAS +0.0000%, KGE_2012 1.0000, obs_mode=synth")
    print(f"  run_id: {run_id}")
    print(f"{'=' * 70}")

    try:
        built_id, _, log_file = build_one(ks, f, tag)
    except Exception as e:
        print(f"\nSMOKE TEST FAILED at BUILD: {e}")
        return 1
    if built_id != run_id:
        print(f"\nSMOKE TEST FAILED: builder made run_id {built_id!r} but this script "
              f"expected {run_id!r}.")
        return 1

    returncode, elapsed, timed_out, tribs_warning = run_with_timeout(args.timeout)
    if timed_out or returncode != 0 or tribs_warning:
        why = ("wall-clock timeout" if timed_out else
               "tRIBS reported non-zero exit" if tribs_warning else
               f"{SCORER_SCRIPT} exited {returncode}")
        print(f"\nSMOKE TEST FAILED: {why}. See {log_file}.")
        return 1

    metrics_file = summary_dir / f"{run_id}_metrics_summary.csv"
    if not metrics_file.exists():
        print(f"\nSMOKE TEST FAILED: run said SUCCESS but {metrics_file} does not exist.")
        return 1
    m = pd.read_csv(metrics_file).iloc[0].to_dict()

    probs    = readback_problems(m, ks, f)
    pbias    = _num(m, "pbias_pct")
    kge_2012 = _num(m, "kge_2012")
    print(f"\n  obs_mode={m.get('obs_mode')}   PBIAS={pbias:+.4f}%   KGE_2012={kge_2012:.4f}")
    print(f"  This one run took {elapsed / 60:.2f} min  ->  400 runs is roughly "
          f"{400 * elapsed / 3600:.1f} h")

    passed = (not probs and abs(pbias) < SMOKE_PBIAS_TOL and kge_2012 > SMOKE_KGE_MIN)
    print(f"\n{'=' * 70}")
    if passed:
        print("SMOKE TEST: PASS. Plumbing, truth file and scorer agree.")
        print("  Next: python run_lhs_Ks_f_120_tribs6.py --limit 5     (optional short trial)")
        print("   or:  nohup python run_lhs_Ks_f_120_tribs6.py > series120_storm100.log 2>&1 &")
        code = 0
    else:
        print("SMOKE TEST: FAIL. Do not start the sweep.")
        for p in probs:
            print(f"  - read-back: {p}")
        if abs(pbias) >= SMOKE_PBIAS_TOL or not (kge_2012 > SMOKE_KGE_MIN):
            print(f"  - the truth point did not reproduce the truth file "
                  f"(need |PBIAS| < {SMOKE_PBIAS_TOL}% and KGE_2012 > {SMOKE_KGE_MIN}).")
            print("    verify_truth_provenance_110.py asks the same question directly.")
        code = 1
    print(f"{'=' * 70}")
    return code


# ======================================================================
# Main
# ======================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Series 120 -- Scenario 1 (1.0x storm) joint Ks_mult x f_RS_abs "
                    "sweep on tRIBS 6.0.0, scored against the synthetic truth, "
                    "paired point-for-point with old Series 100.")
    parser.add_argument("--n", type=int, default=DEFAULT_N,
                        help=f"Number of draws (default {DEFAULT_N}; changing it changes "
                             f"EVERY draw, so the points will not pair with old Series 100)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help=f"Random seed (default {DEFAULT_SEED}). A new seed gives a new "
                             f"batch of points -- use it for extension batches")
    parser.add_argument("--dry_run", action="store_true",
                        help="Generate draws and print the pairing check, then stop. "
                             "Builds nothing, runs no tRIBS")
    parser.add_argument("--smoke_test", action="store_true",
                        help="Run ONE build+run+score at the truth point and stop")
    parser.add_argument("--limit", type=int, default=None,
                        help="Only process the first N draws (for a short trial; "
                             "resume later with --skip_existing)")
    parser.add_argument("--skip_existing", action="store_true",
                        help="Skip draws whose compare CSV already exists (resume)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Per-run hard timeout in seconds (default 300)")
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Keep raw per-run results folders instead of deleting them "
                             "after a successful run")
    parser.add_argument("--old_csv", type=str, default=None,
                        help="Path to the old Series 100 results CSV (for pairing); "
                             "default: look in summary_tables/")
    parser.add_argument("--allow_other_truth", action="store_true",
                        help="Allow a truth file whose name contains none of "
                             f"{', '.join(repr(h) for h in TRUTH_NAME_HINTS)} "
                             "(only if you are sure it is the Ks 7.0 / f 0.012 truth)")
    args = parser.parse_args()

    script_dir   = Path.cwd()
    project_root = script_dir.parent
    calib_dir    = project_root / "calibration_work"
    summary_dir  = calib_dir / "03_comparisons" / "summary_tables"
    summary_dir.mkdir(parents=True, exist_ok=True)

    samples = generate_lhs_samples(args.n, LHS_PARAMS, seed=args.seed)
    samples["caution"] = [in_caution_zone(k, f)
                          for k, f in zip(samples["Ks_mult"], samples["f_RS_abs"])]
    todo = samples.iloc[:args.limit] if args.limit else samples

    old_df, old_path = find_old_series100(summary_dir, args.old_csv)
    twin_idx = print_pairing_check(samples, old_df, old_path, args.n, args.seed)

    variant         = variant_suffix(args.n, args.seed)
    out_path        = summary_dir / f"lhs_results_synth_Ks_f_{SERIES}{variant}.csv"
    failed_log_path = summary_dir / f"lhs_results_synth_Ks_f_FAILED_{SERIES}{variant}.csv"
    assert "_100" not in out_path.name, "output would collide with an old Series 100 file"

    # ------------------------------------------------------------------
    # DRY RUN -- nothing built, no truth file needed
    # ------------------------------------------------------------------
    if args.dry_run:
        qouts = list((calib_dir / "synth_truth").glob("*.qout")) \
            if (calib_dir / "synth_truth").exists() else []
        print(f"\nDRY RUN -- {len(samples)} draws, seed {args.seed}. Nothing was built or run.")
        for p in PARAM_COLS:
            print(f"  {p}: draws span {samples[p].min():.4f} - {samples[p].max():.4f}  "
                  f"(box {LHS_PARAMS[p]['lo']} - {LHS_PARAMS[p]['hi']}, {LHS_PARAMS[p]['scale']})")
        print(f"  Draws in the caution zone (Ks >= {CAUTION_KS_LO}, f >= {CAUTION_F_LO}): "
              f"{int(samples['caution'].sum())}")
        print("  First 5 draws:")
        for i, r in samples.head(5).iterrows():
            print(f"    [{i + 1:>3}]  Ks={r['Ks_mult']:.4f}  f={r['f_RS_abs']:.5f}")
        print(f"  Output would go to: {out_path.name}")
        print(f"  Truth file(s) in synth_truth/ right now: {[q.name for q in qouts] or 'none'}")
        return 0

    # ------------------------------------------------------------------
    # HARD SAFETY CHECKS -- run from lab/, synthetic truth active and the right one
    # ------------------------------------------------------------------
    if not (script_dir / SCORER_SCRIPT).exists():
        raise RuntimeError(f"{SCORER_SCRIPT} not found in {script_dir}. Run this script "
                           f"from the lab/ directory.")
    synth_dir  = calib_dir / "synth_truth"
    qout_files = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
    if len(qout_files) != 1:
        raise RuntimeError(
            f"Expected exactly one *.qout file in {synth_dir} to activate synthetic-truth "
            f"mode, found {len(qout_files)}: {[f.name for f in qout_files]}. Without it the "
            f"scorer silently falls back to the real gauge.")
    truth_name = qout_files[0].name
    print(f"\nSynthetic truth mode confirmed active: {truth_name}")
    name_ok = any(h in truth_name.lower() for h in TRUTH_NAME_HINTS)
    if name_ok:
        print("  The name looks like the Ks 7.0 / f 0.012 truth. (A name is only a hint; "
              "the smoke test is what proves the file's contents.)")
    elif args.allow_other_truth:
        print(f"  NOTE: the name contains none of "
              f"{', '.join(repr(h) for h in TRUTH_NAME_HINTS)}, but --allow_other_truth was "
              f"given, so continuing. Only the smoke test can show this file is the "
              f"Ks 7.0 / f 0.012 truth.")
    else:
        raise RuntimeError(
            f"The truth file is {truth_name!r}, whose name contains none of "
            f"{', '.join(repr(h) for h in TRUTH_NAME_HINTS)}. This script is the 1.0x leg "
            f"(truth Ks 7.0, f 0.012). If you are sure this file is that truth, re-run with "
            f"--allow_other_truth.")

    install_series_120_hooks()

    if args.smoke_test:
        return run_smoke_test(args, calib_dir, summary_dir)

    # ------------------------------------------------------------------
    # SWEEP
    # ------------------------------------------------------------------
    existing_df  = load_existing_results(out_path)
    existing_ids = set(existing_df["run_id"].values) if not existing_df.empty else set()
    results      = existing_df.to_dict("records") if not existing_df.empty else []
    if not existing_df.empty and not args.skip_existing:
        print("  NOTE: existing results found and --skip_existing not given: every draw will "
              "be re-run (rows are replaced by run_id). Use --skip_existing to resume.")

    # The FAILED log is carried forward on a resume, so an interruption does not
    # erase the record of earlier HANG/FAILED draws. A draw that later succeeds
    # is dropped from it.
    failed_log_rows = []
    if args.skip_existing and failed_log_path.exists():
        try:
            failed_log_rows = pd.read_csv(failed_log_path).to_dict("records")
            print(f"  Loaded {len(failed_log_rows)} earlier HANG/FAILED record(s) from "
                  f"{failed_log_path.name}")
        except Exception as e:
            print(f"  Warning: could not read {failed_log_path.name} ({e}); starting a new log.")
    n_total = len(todo)

    print(f"\n{'=' * 70}")
    print(f"Series {SERIES} -- Ks_mult x f_RS_abs vs SYNTHETIC TRUTH ({truth_name}), "
          f"storm leg {STORM_LEG}, channel loss OFF")
    print(f"  {n_total} of {len(samples)} draws, seed={args.seed}, timeout={args.timeout}s, "
          f"cleanup={'off' if args.no_cleanup else 'on'}")
    for p, b in LHS_PARAMS.items():
        print(f"  {p}: {b['lo']} - {b['hi']}  ({b['scale']})")
    print(f"  PINNED (truth): kinemvelcoef={ROUTING_TRUTH['kinemvelcoef']}  "
          f"flowexp r={ROUTING_TRUTH['flowexp']}  channelroughness={ROUTING_TRUTH['channelroughness']}"
          f"  optpercolation=0")
    print(f"  Output: {out_path.name}")
    print(f"{'=' * 70}\n")

    completed = skipped = excluded = hung = failed = cleaned = caution_hits = 0
    aborted = False
    sweep_start = time.time()

    for i, row in todo.iterrows():
        ks      = float(row["Ks_mult"])
        f_val   = float(row["f_RS_abs"])
        flagged = bool(row["caution"])
        tag     = sweep_tag(f_val, args.seed)
        run_id  = predicted_run_id(ks, tag)

        if is_truth_point(ks, f_val):
            print(f"[{i + 1:>3}/{len(samples)}]  EXCLUDED: draw equals the truth point.")
            excluded += 1
            continue

        flag_note = "  [CAUTION ZONE]" if flagged else ""
        print(f"\n[{i + 1:>3}/{len(samples)}]  Ks={ks:.3f}x  f={f_val:.5f}  -> {run_id}{flag_note}")
        if flagged:
            caution_hits += 1

        # ---- resume ----
        if args.skip_existing and csv_already_exists(run_id, calib_dir):
            print(f"  SKIP (compare CSV exists): {run_id}")
            skipped += 1
            metrics_file = summary_dir / f"{run_id}_metrics_summary.csv"
            if metrics_file.exists() and run_id not in existing_ids:
                try:
                    m = pd.read_csv(metrics_file).iloc[0].to_dict()
                    m.update({"series": SERIES, "storm_leg": STORM_LEG, "seed": args.seed,
                              "draw_idx": int(i), "caution_zone": flagged,
                              "swept_param": f"lhs_synth_Ks_f_{SERIES}"})
                    results.append(m)
                    existing_ids.add(run_id)
                except Exception:
                    pass
            continue

        # ---- build ----
        t_run = time.time()
        try:
            built_id, _, log_file = build_one(ks, f_val, tag)
            if built_id != run_id:
                raise RuntimeError(f"builder made run_id {built_id!r}, expected {run_id!r}")
        except Exception as e:
            print(f"  BUILD FAILED: {e}")
            failed += 1
            failed_log_rows = upsert_failed(failed_log_rows, {
                "run_id": run_id, "status": "BUILD_ERROR", "reason": str(e),
                "elapsed_min": (time.time() - t_run) / 60,
                "draw_idx": int(i), "Ks_mult": ks, "f_RS_abs": f_val})
            pd.DataFrame(failed_log_rows).to_csv(failed_log_path, index=False)
            continue

        # ---- run + score ----
        returncode, elapsed, timed_out, tribs_warning = run_with_timeout(args.timeout)
        reason = ""
        if timed_out:
            status, reason = "HANG", "wall-clock timeout"
        elif tribs_warning:
            status, reason = "FAILED", "tRIBS reported non-zero exit"
        elif returncode != 0:
            status, reason = "FAILED", f"scorer exited {returncode}"
        else:
            status = "SUCCESS"

        if status == "SUCCESS":
            metrics_file = summary_dir / f"{run_id}_metrics_summary.csv"
            try:
                metrics = pd.read_csv(metrics_file).iloc[0].to_dict()
            except Exception as e:
                status, reason = "FAILED", f"metrics file unreadable: {e}"
            else:
                probs = readback_problems(metrics, ks, f_val)
                if probs:
                    print("\n" + "!" * 70)
                    print("READ-BACK MISMATCH -- STOPPING THE SWEEP")
                    for p in probs:
                        print(f"  - {p}")
                    print("The scored row does not match what this script asked for, so the "
                          "remaining runs would be wasted. Results so far are saved.")
                    print("!" * 70)
                    failed_log_rows = upsert_failed(failed_log_rows, {
                        "run_id": run_id, "status": "READBACK_MISMATCH",
                        "reason": "; ".join(probs),
                        "elapsed_min": elapsed / 60, "draw_idx": int(i),
                        "Ks_mult": ks, "f_RS_abs": f_val})
                    pd.DataFrame(failed_log_rows).to_csv(failed_log_path, index=False)
                    failed += 1
                    aborted = True
                    break

                metrics.update({"series": SERIES, "storm_leg": STORM_LEG, "seed": args.seed,
                                "draw_idx": int(i), "caution_zone": flagged,
                                "swept_param": f"lhs_synth_Ks_f_{SERIES}"})
                results = [r for r in results if r.get("run_id") != run_id]
                results.append(metrics)
                completed += 1
                if any(r.get("run_id") == run_id for r in failed_log_rows):
                    failed_log_rows = [r for r in failed_log_rows if r.get("run_id") != run_id]
                    pd.DataFrame(failed_log_rows, columns=FAILED_COLS).to_csv(
                        failed_log_path, index=False)
                print(f"  NEW (tRIBS 6.0):  KGE_2012={_num(metrics, 'kge_2012'):.3f}  "
                      f"PBIAS={_num(metrics, 'pbias_pct'):+.1f}%")
                tw = twin_idx[i] if i < len(twin_idx) else -1
                if old_df is not None and tw >= 0:
                    o = old_df.iloc[tw]
                    print(f"  OLD twin (5.3.1): KGE_2012={_num(o, 'kge_2012'):.3f}  "
                          f"PBIAS={_num(o, 'pbias_pct'):+.1f}%")
                if not args.no_cleanup and cleanup_raw_results(run_id, calib_dir):
                    cleaned += 1

        if status != "SUCCESS":
            print(f"  {status}: {run_id}  ({reason})")
            if status == "HANG":
                hung += 1
            else:
                failed += 1
            failed_log_rows = upsert_failed(failed_log_rows, {
                "run_id": run_id, "status": status, "reason": reason,
                "elapsed_min": elapsed / 60, "draw_idx": int(i),
                "Ks_mult": ks, "f_RS_abs": f_val})
            pd.DataFrame(failed_log_rows).to_csv(failed_log_path, index=False)

        attempted = completed + hung + failed
        remaining = n_total - attempted - skipped - excluded
        if attempted > 0:
            avg = (time.time() - sweep_start) / attempted
            eta_min = avg * max(remaining, 0) / 60
            print(f"  Run time: {(time.time() - t_run) / 60:.1f} min  |  "
                  f"ETA: {eta_min:.0f} min ({eta_min / 60:.1f} h) remaining")

        # save after every run
        if results:
            pd.DataFrame(results).to_csv(out_path, index=False)

    # ------------------------------------------------------------------
    # Wrap up
    # ------------------------------------------------------------------
    print(f"\nSweep {'STOPPED EARLY' if aborted else 'complete'}: {completed} ran, "
          f"{skipped} skipped, {excluded} excluded, {hung} hung, {failed} failed, "
          f"{cleaned} raw folders cleaned, {caution_hits} in the caution zone")

    if results:
        final_df = pd.DataFrame(results)
        if "kge_2012" in final_df.columns:
            final_df = final_df.sort_values("kge_2012", ascending=False)
        final_df.to_csv(out_path, index=False)
        print(f"Saved: {out_path.name}  ({len(final_df)} rows)")
        try:
            report_final(final_df, old_df, failed_log_rows, n_total, len(samples))
        except Exception as e:
            print(f"\n(report skipped: {e} -- results are saved regardless)")
    else:
        print("No results to save.")

    return 1 if aborted else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted (Ctrl+C). The run in flight was killed; every finished run is "
              "already saved. Resume with the same command plus --skip_existing.")
        sys.exit(130)
