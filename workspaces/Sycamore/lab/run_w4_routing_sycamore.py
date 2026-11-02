#!/usr/bin/env python3
"""
run_w4_routing_sycamore.py  --  Stage 5, Step 7: run draw 12 on the long window W4 (Oct 1 2016 - Apr 1 2017) with the NEW ROUTING
chosen by the W1B second-storm check.

Why. The W1B check (Stage5_Routing_Log.md, Section 8.5) picked one routing setting at both anchors by the approved rule: sweep point
k10 of rt1_list.csv,
        cv (KINEMVELCOEF) 6.39,  flowexp (FLOWEXP) 0.2973,  n (CHANNELROUGHNESS) 0.067       (baseline was 6, 0.3, 0.04).
It met arrival within 1 h and peak within 15% on the flood week (W1) and on Feb 27-28 (W1B). This run puts it on the whole season,
which is the first test of what it does to the Oct-Dec false flow, the Feb 18-19 overshoot, the tails and the by-storm table.

Which run. ONE run: anchor draw 12 (Ks 0.187, f 0.5874, PsiB 0.2939 x baseline), the same soil setting as the earlier W4 run
W4_ks0p187_cv6_psib0p2939_f0p5874_d012_nocl, but with the three routing values above. Everything else is identical to that run:
channel conductivity off (OPTPERCOLATION 0), snow on (OPTSNOW 1), width coefficient 2.33, thetaS, m, anisotropy and bedrock at baseline
(bedrock = SOLUS + 1 m), warm start (water table 0.95 x bedrock), Oct 1 2016, 4,368 h. The new run is named
W4_ks0p187_cv6p39_psib0p2939_f0p5874_d012_fe0p2973_n0p067_nocl.   (--anchors 25 also works: draw 25 was in the W1B check too.)

This script IMPORTS run_w4_draw12_sycamore.py (keep it in the same folder) and reuses its tested build, check, tRIBS and scoring
functions. It changes only: (1) the three routing values for the planned run; (2) every output file has its own name (w4rt_*), so the
earlier w4d12_* files are never touched; (3) the run name; (4) the summary, which prints the new run NEXT TO the earlier draw-12 run
at baseline routing (read only) and the W4 Stage 2 baseline; (5) two extra checks (below).

Extra checks added here:
  * The routing values are read from row k=10 of rt1_list.csv (the frozen sweep list, md5 c1e799f82ff16d30d91ff377b63356b7) and must
    equal the values above. If the list was changed or row 10 differs, nothing starts.
  * After the build, the new .in file is compared with the .in file of the earlier draw-12 W4 run at baseline routing (if it is in
    lab/): every keyword must be equal except cv, flowexp, n and the output name. (The usual check against the W4 Stage 2 baseline's
    .in file is also done, as before.)

NOTHING IS REPEATED OR OVERWRITTEN. The run has a unique name, and a planned run whose name already has files in results/ stops the
start-up check. The existing W4 baseline and the earlier draw-12 W4 run are only read (scored into the comparison table, never rerun).

Run it from workspaces/Sycamore/lab/, in this order:
    python -u run_w4_routing_sycamore.py --dry-run               # checks only, builds nothing
    python -u run_w4_routing_sycamore.py --build-only            # builds the run and checks it, no tRIBS (about 1 min)
    python -u run_w4_routing_sycamore.py                         # the full run (about 36 min), then paste the table
    python -u run_w4_routing_sycamore.py --summary-only          # just reprint the table from what is finished

Options (as in run_w4_draw12_sycamore.py): --anchors 12  --limit N  --redo ID  --timeout-min N  --heartbeat-min N  --notebook FILE

Files it makes (all in lab/): w4rt_list.csv (the frozen list), w4rt_manifest.json, w2.lock (only while running), w4rt_logs/,
score_out_w4rt/scores_w4rt.csv, w4rt_summary.csv, and <run>.in, simulation_<run>.log and results/<run>*.

Do NOT run the notebook (or any cell that rebuilds data/model/) while this runs: every build rewrites lab/data/model/. Keep the
terminal open: the progress lines also stop the Codespace from timing out for being idle. If the connection drops, run the same
command again: a finished run is skipped.
"""
from __future__ import annotations

import csv
import hashlib
import math
import sys
from pathlib import Path

import run_w4_draw12_sycamore as base      # the tested W4 runner: its functions are reused here

# --------------------------------------------------------------------------------------------
# What this batch does
# --------------------------------------------------------------------------------------------
K_POINT = 10                                                       # sweep point 10 of rt1_list.csv (the W1B winner at both anchors)
ROUTING = {"KINEMVELCOEF": 6.39, "FLOWEXP": 0.2973, "CHANNELROUGHNESS": 0.067}
ROUTING_KEYS = tuple(ROUTING)
RT1_LIST = "rt1_list.csv"
RT1_LIST_MD5 = "c1e799f82ff16d30d91ff377b63356b7"                  # the frozen sweep list (Output, Nov 1)
ALLOWED_DRAWS = (12, 25)                                           # the two anchors that went through the W1B check
EXPECTED_DRAW25 = (0.1852, 0.5055, 0.3358)                         # Ks, f, PsiB of draw 25 (handoff values)
GROUP = "k10_routing"
COMPARE_SCORES = "score_out_w4d12/scores_w4d12.csv"                # the earlier draw-12 run (baseline routing): read only
SCRIPT_NAME = "run_w4_routing_sycamore.py"

# Output names of THIS batch (the earlier batch's w4d12_* names are never used here).
base.LOG_DIR = "w4rt_logs"
base.SCORE_DIR = "score_out_w4rt"
base.SCORE_CSV = "scores_w4rt.csv"
base.SUMMARY_CSV = "w4rt_summary.csv"
base.DRAWS_CSV = "w4rt_list.csv"
base.MANIFEST = "w4rt_manifest.json"
base.EXPECTED_DRAWS[25] = EXPECTED_DRAW25

CRITERIA_ARRIVAL_H, CRITERIA_PEAK_PCT = 1.0, 15.0                  # the closing criteria (reported only; nothing is decided here)

_orig_read_anchors = base.read_anchors
_orig_settings_for = base.settings_for
_orig_describe = base.describe
_orig_check_in_file = base.check_in_file
_orig_preflight = base.preflight
_orig_save_manifest = base.save_manifest


def log(msg: str = "") -> None:
    base.log(msg)


def file_md5(path: Path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


# --------------------------------------------------------------------------------------------
# Names and the routing source
# --------------------------------------------------------------------------------------------
def baseline_routing_id(job: dict) -> str:
    """The earlier W4 run of the SAME soil at baseline routing, e.g. W4_ks0p187_cv6_psib0p2939_f0p5874_d012_nocl."""
    return base.run_id_for(job["ks"], job["f"], job["psib"], job["idx"])


def new_run_id(ks, f, psib, idx) -> str:
    v = base.vtag
    return (f"W4_ks{v(ks)}_cv{v(ROUTING['KINEMVELCOEF'])}_psib{v(psib)}_f{v(f)}_d{idx:03d}"
            f"_fe{v(ROUTING['FLOWEXP'])}_n{v(ROUTING['CHANNELROUGHNESS'])}_nocl")


def check_routing_source(lab: Path) -> str:
    """Row K_POINT of the frozen sweep list must hold exactly the routing values above. Raises base.Abort."""
    p = lab / RT1_LIST
    if not p.is_file():
        raise base.Abort(f"{RT1_LIST} not found in lab/. It is the frozen routing list the W1B check was made from.")
    md5 = file_md5(p)
    if md5 != RT1_LIST_MD5:
        raise base.Abort(f"{RT1_LIST} has md5 {md5}, not the frozen {RT1_LIST_MD5}. The list must not change; restore it.")
    with open(p, newline="") as fh:
        rows = {int(r["k"]): r for r in csv.DictReader(fh)}
    r = rows.get(K_POINT)
    if r is None:
        raise base.Abort(f"{RT1_LIST} has no row k = {K_POINT}.")
    got = {"KINEMVELCOEF": float(r["cv"]), "FLOWEXP": float(r["flowexp"]), "CHANNELROUGHNESS": float(r["n"])}
    bad = [f"{k}: list has {got[k]!r}, this script has {ROUTING[k]!r}" for k in ROUTING if abs(got[k] - ROUTING[k]) > 1e-9]
    if bad:
        raise base.Abort(f"row k = {K_POINT} of {RT1_LIST} does not match the routing values in this script:\n  " + "\n  ".join(bad))
    return f"routing = row k={K_POINT} of {RT1_LIST} (md5 {md5[:8]}...): cv {ROUTING['KINEMVELCOEF']:g}, flowexp {ROUTING['FLOWEXP']:g}, n {ROUTING['CHANNELROUGHNESS']:g}"


# --------------------------------------------------------------------------------------------
# Overrides of the base functions (the base module looks these up by name at call time)
# --------------------------------------------------------------------------------------------
def read_anchors(lab: Path, anchors_arg: str) -> list:
    """The base checks each draw against the Stage 3 design table and the handoff values (under the baseline-routing names).
    Here the run names are then replaced by the new ones, and the routing source is checked."""
    check_routing_source(lab)
    draws = _orig_read_anchors(lab, anchors_arg)
    bad = [d["idx"] for d in draws if d["idx"] not in ALLOWED_DRAWS]
    if bad:
        raise base.Abort(f"draw(s) {bad}: this routing went through the W1B check at draws {list(ALLOWED_DRAWS)} only.")
    for d in draws:
        d["run_id"] = new_run_id(d["ks"], d["f"], d["psib"], d["idx"])
        d["group"] = GROUP
    return draws


def settings_for(job: dict) -> dict:
    """The base settings, with the new routing for the planned run (the baseline reference job keeps the baseline routing)."""
    s = _orig_settings_for(job)
    if job.get("kind") == "run":
        s.update(ROUTING)
    return s


def describe(job: dict) -> str:
    text = _orig_describe(job)
    if job.get("kind") == "run":
        text += f"; routing cv {ROUTING['KINEMVELCOEF']:g}, flowexp {ROUTING['FLOWEXP']:g}, n {ROUTING['CHANNELROUGHNESS']:g}"
    return text


def save_manifest(man: dict) -> None:
    man["script"] = SCRIPT_NAME
    man["routing"] = dict(ROUTING, source=f"{RT1_LIST} row k={K_POINT}", rt1_list_md5=RT1_LIST_MD5)
    _orig_save_manifest(man)


def check_in_file(new_in: Path, ref_in: Path, job: dict) -> None:
    """The usual check against the W4 Stage 2 baseline's .in (which now also expects the new routing values), plus a comparison with
    the .in of the earlier draw-12 W4 run (same soil, baseline routing): equal in every keyword except cv, flowexp, n, output name."""
    _orig_check_in_file(new_in, ref_in, job)
    earlier = Path(f"{baseline_routing_id(job)}.in")
    if not earlier.is_file():
        log(f"{job['run_id']}: note: {earlier.name} is not in lab/, so the same-soil .in comparison was skipped")
        return
    new, old = base.read_in_file(new_in), base.read_in_file(earlier)
    allowed = set(ROUTING_KEYS) | {"OUTFILENAME"}
    problems = []
    for key in sorted(set(new) | set(old)):
        if key in allowed:
            continue
        a, b = old.get(key), new.get(key)
        if a is None or b is None or not base.same_value(a, b):
            problems.append(f"{key}: earlier run {a!r} but this run {b!r}")
    for key, want in ROUTING.items():
        if not base.same_value(new.get(key), want):
            problems.append(f"{key}: this run reads {new.get(key)!r}, wanted {want!r}")
    if problems:
        raise base.RunError(f"The .in file differs from {earlier.name} (same soil, baseline routing) in more than the routing and the "
                            f"output name:\n  " + "\n  ".join(problems))
    differs = sum(1 for k in ROUTING_KEYS if not base.same_value(old.get(k), new.get(k)))
    log(f"{job['run_id']}: same-soil comparison OK: the .in matches {earlier.name} in every keyword except cv, flowexp, n and the "
        f"output name ({differs} of the 3 routing values differ from baseline)")


def preflight(args, lab: Path, jobs: list, base_problems: list, nb_note: str, man: dict) -> list:
    out = _orig_preflight(args, lab, jobs, base_problems, nb_note, man)
    try:
        out.append(("OK", check_routing_source(lab)))
    except base.Abort as exc:
        out.append(("FAIL", str(exc)))
    cmp_scores = base_scores_for_compare(lab)
    for j in jobs:
        rid = baseline_routing_id(j)
        has_in = (lab / f"{rid}.in").is_file()
        has_q = (lab / "results" / f"{rid}_Outlet.qout").is_file()
        has_sc = rid in cmp_scores
        if has_in and has_q and has_sc:
            out.append(("OK", f"comparison run {rid} (same soil, baseline routing): .in, outlet file and scores are there (read only)"))
        else:
            missing = [n for n, ok in ((".in file", has_in), ("outlet file", has_q), (f"score row in {COMPARE_SCORES}", has_sc)) if not ok]
            out.append(("WARN", f"comparison run {rid} is missing: {', '.join(missing)}. The run itself is not affected; the table will show '-' for it "
                                f"and the same-soil .in comparison may be skipped."))
    return out


base.read_anchors = read_anchors
base.settings_for = settings_for
base.describe = describe
base.save_manifest = save_manifest
base.check_in_file = check_in_file
base.preflight = preflight


# --------------------------------------------------------------------------------------------
# The tables
# --------------------------------------------------------------------------------------------
def base_scores_for_compare(lab: Path) -> dict:
    out = {}
    p = lab / COMPARE_SCORES
    if p.is_file():
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                out[r["run_id"]] = r
    return out


def f_(x, nd=2, sign=False) -> str:
    return base.fm(x, nd, sign)


def n_(sc, col) -> float:
    return base.fnum(sc, col)


STORMS = [("flood_week", "flood week Jan 18-25"), ("jan_second_storm", "second Jan storm"), ("feb18_19", "Feb 18-19"),
          ("feb27_28", "Feb 27-28")]


def print_summary(man: dict, lab: Path, draws: list) -> None:
    mine = base.load_scores(lab)
    other = base_scores_for_compare(lab)
    get = lambda rid: mine.get(rid) or other.get(rid)
    log()
    log(f"W4 WITH THE W1B-CHOSEN ROUTING (Oct 1 2016 - Apr 1 2017). Routing: cv {ROUTING['KINEMVELCOEF']:g}, flowexp {ROUTING['FLOWEXP']:g}, "
        f"n {ROUTING['CHANNELROUGHNESS']:g} (baseline 6, 0.3, 0.04). Warm start (water table 0.95 x bedrock); channel loss off; snow on.")
    log("  Oct-Dec = model volume Oct 1 - Dec 31 in million m3 (gauge: 0). days = days the model flow is above 0.13 m3/s; max = largest Oct-Dec flow (m3/s).")
    log("  flood = flood week Jan 18-25; F27 = Feb 27-28. vol = volume ratio (model / gauge), pk = peak error %, arr = first-arrival error h (negative = early).")
    log("  Oct-Mar = whole-window volume ratio; KGE = KGE_2012 over the whole span (read only after the volume gate; not a pass/fail here).")
    log("  Nothing is promoted or decided here: this table shows what the routing does over the whole season.")
    head = (f"  {'run':<44} {'OctDec':>7} {'days':>5} {'max':>6} | {'flood vol':>9} {'pk %':>6} {'arr h':>6} | "
            f"{'F27 vol':>7} {'pk %':>6} {'arr h':>6} | {'Oct-Mar':>7} {'KGE':>6}")

    def line(label, sc):
        if sc is None:
            return f"  {label:<44} (no scores)"
        return (f"  {label:<44} {f_(n_(sc, base.C_FFV), 2):>7} {f_(n_(sc, base.C_FFDAYS), 0):>5} {f_(n_(sc, base.C_FFMAX), 1):>6} | "
                f"{f_(n_(sc, base.C_VOL), 3):>9} {f_(n_(sc, base.C_PK), 1, True):>6} {f_(n_(sc, base.C_ARR), 2, True):>6} | "
                f"{f_(n_(sc, base.C_F27), 3):>7} {f_(n_(sc, 'feb27_28__peak_err_pct'), 1, True):>6} {f_(n_(sc, 'feb27_28__arrival_err_h'), 2, True):>6} | "
                f"{f_(n_(sc, base.C_OM), 3):>7} {f_(n_(sc, 'span__kge2012'), 3):>6}")

    ref = get(base.BASELINE_ID)
    if ref is not None:
        log(f"  Gauge, Oct-Dec: volume {f_(n_(ref, base.C_FFG), 3)} M m3, largest flow {f_(n_(ref, base.C_FFGMAX), 2)} m3/s, "
            f"{f_(n_(ref, base.C_FFGDAYS), 0)} days above {base.FF_THRESHOLD:g}.")
    log()
    log(head)
    if ref is not None:
        log(line("W4 baseline (Ks 0.2, f 1, PsiB 0.5; cv 6)", ref))

    rows, finished = [], []
    for d in draws:
        rid = d["run_id"]
        rec = man["runs"].get(rid) or {}
        status = rec.get("status", "not run")
        cmp_id = base.run_id_for(d["ks"], d["f"], d["psib"], d["idx"])
        cmp_sc = get(cmp_id)
        new_sc = mine.get(rid) if status in ("done", "unscored") else None
        log(line(f"draw {d['idx']}, baseline routing (cv 6, 0.3, 0.04)", cmp_sc))
        if new_sc is None:
            log(f"  {'draw ' + str(d['idx']) + ', new routing':<44} ({status})")
        else:
            log(line(f"draw {d['idx']}, NEW routing ({ROUTING['KINEMVELCOEF']:g}, {ROUTING['FLOWEXP']:g}, {ROUTING['CHANNELROUGHNESS']:g})", new_sc))
            finished.append((d, new_sc, cmp_sc))
        row = dict(idx=d["idx"], run_id=rid, group=d["group"], KS_MULT=d["ks"], F_MULT=d["f"], PSIB_MULT=d["psib"],
                   cv=ROUTING["KINEMVELCOEF"], flowexp=ROUTING["FLOWEXP"], n=ROUTING["CHANNELROUGHNESS"], status=status,
                   compare_run_id=cmp_id)
        if new_sc is not None:
            for k, v in new_sc.items():
                if k not in ("run_id", "window"):
                    row[k] = v
            row["build_min"], row["tribs_min"], row["output_mb"] = rec.get("build_min", ""), rec.get("tribs_min", ""), rec.get("output_mb", "")
        rows.append(row)

    for d, sc, cmp_sc in finished:
        log()
        log(f"  BY STORM, draw {d['idx']}   (vol = volume ratio; pk = peak error %; arr = first-arrival error h; KGE = KGE_2012; gamma = KGE gamma)")
        log(f"  {'':<22} | {'baseline routing':^40} | {'NEW routing':^40}")
        log(f"  {'storm':<22} | {'vol':>6} {'pk %':>7} {'arr h':>6} {'KGE':>6} {'gamma':>6} | {'vol':>6} {'pk %':>7} {'arr h':>6} {'KGE':>6} {'gamma':>6}")

        def cells(s, key):
            if s is None:
                return f"{'-':>6} {'-':>7} {'-':>6} {'-':>6} {'-':>6}"
            return (f"{f_(n_(s, key + '__vol_ratio'), 3):>6} {f_(n_(s, key + '__peak_err_pct'), 1, True):>7} "
                    f"{f_(n_(s, key + '__arrival_err_h'), 2, True):>6} {f_(n_(s, key + '__kge2012'), 2):>6} {f_(n_(s, key + '__kge_gamma'), 2):>6}")
        for key, label in STORMS:
            log(f"  {label:<22} | {cells(cmp_sc, key)} | {cells(sc, key)}")
        for key, label in (("tail_jan26_feb2", "tail Jan 26-Feb 2"), ("tail_mar3_mar11", "tail Mar 3-11")):
            def tcell(s):
                if s is None:
                    return f"{'-':>6} (recession rate ratio {'-':>5})"
                return f"{f_(n_(s, key + '__tail_vol_ratio'), 3):>6} (recession rate ratio {f_(n_(s, key + '__recession_rate_ratio'), 2):>5})"
            log(f"  {label:<22} | tail vol {tcell(cmp_sc)} | tail vol {tcell(sc)}")
        mo = []
        for tag, label in (("2017_01", "Jan"), ("2017_02", "Feb"), ("2017_03", "Mar")):
            a = f_(n_(cmp_sc, f"month_{tag}__vol_ratio"), 3) if cmp_sc else "-"
            b = f_(n_(sc, f"month_{tag}__vol_ratio"), 3)
            mo.append(f"{label} {a} -> {b}")
        log("  monthly volume ratio (baseline routing -> NEW): " + "; ".join(mo))
        if cmp_sc is not None and math.isfinite(n_(cmp_sc, base.C_FFV)) and math.isfinite(n_(sc, base.C_FFV)):
            a, b = n_(cmp_sc, base.C_FFV), n_(sc, base.C_FFV)
            log(f"  Oct-Dec false flow: {f_(a, 2)} -> {f_(b, 2)} M m3 ({f_(b - a, 2, True)} M m3"
                + (f", {100 * (b / a - 1):+.0f}%" if a > 0 else "") + f"); days above 0.13: {f_(n_(cmp_sc, base.C_FFDAYS), 0)} -> {f_(n_(sc, base.C_FFDAYS), 0)}; "
                f"largest Oct-Dec flow {f_(n_(cmp_sc, base.C_FFMAX), 1)} -> {f_(n_(sc, base.C_FFMAX), 1)} m3/s")
        log()
        log(f"  CLOSING CRITERIA on this run (arrival within {CRITERIA_ARRIVAL_H:g} h, peak within {CRITERIA_PEAK_PCT:g}%): reported only")
        for key, label in (("flood_week", "flood week"), ("feb27_28", "Feb 27-28")):
            arr, pk = n_(sc, f"{key}__arrival_err_h"), n_(sc, f"{key}__peak_err_pct")
            ok_a = math.isfinite(arr) and abs(arr) <= CRITERIA_ARRIVAL_H
            ok_p = math.isfinite(pk) and abs(pk) <= CRITERIA_PEAK_PCT
            log(f"    {label:<11} arrival {f_(arr, 2, True):>6} h  {'met' if ok_a else 'NOT met'};  peak {f_(pk, 1, True):>6} %  {'met' if ok_p else 'NOT met'}")

    done = len(finished)
    log()
    log(f"  {done} of {len(draws)} run(s) scored; {sum(1 for r in rows if r['status'] == 'failed')} failed; "
        f"{sum(1 for r in rows if r['status'] == 'not run')} not run yet.")
    rejected = [rid for rid, rec in man["runs"].items() if rec.get("status") == "reuse_rejected"]
    if rejected:
        log("  NOT used (failed their checks) and NOT rerun: " + ", ".join(rejected))
    if rows:
        cols = []
        for r in rows:
            for c in r:
                if c not in cols:
                    cols.append(c)
        with open(lab / base.SUMMARY_CSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c, "") for c in cols})
        log()
        log(f"  Wrote {base.SUMMARY_CSV} (all score columns for the new run(s)).")


base.print_summary = print_summary


def main(argv=None) -> int:
    return base.main(argv)


if __name__ == "__main__":
    sys.exit(main())
