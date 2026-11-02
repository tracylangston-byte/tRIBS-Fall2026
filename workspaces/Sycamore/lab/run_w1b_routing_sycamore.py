#!/usr/bin/env python3
"""
run_w1b_routing_sycamore.py  --  Stage 5, Step 6: the SECOND-STORM CHECK of the routing sweep on W1B (Feb 21 - Mar 3 2017, 240 h).

WHAT IT DOES.  The routing sweep (run_sweep_routing_sycamore.py, batch rt1) scored 90 routing settings on W1, the flood week. This script takes
the best of them at each soil anchor (selected by the rule below, from the W1 scores ALONE), runs them on W1B, and scores the Feb 27-28 storm
(first-arrival error, peak error, volume, KGE) with score_run_sycamore.py. It uses the same .in-edit shortcut as the sweep: build the anchor's
soil ONCE (here on the W1B dates: start Feb 21, 240 h, cold start, water table 0.82 x bedrock, the same bedrock and water-table files as W1),
then for every routing point copy the .in, change only cv (KINEMVELCOEF), flowexp (FLOWEXP), n (CHANNELROUGHNESS) and the output name, and run tRIBS.
About 2.6 min per run.  It runs NOTHING on W4 and it does not promote anything: Tracy decides what goes to W4.

THE RULES (agreed by Tracy on Nov 1 2026, BEFORE any W1B result exists; they are written into the frozen list's header columns)
  Selection (--make-list; from rt1_summary.csv, the W1 scores only).  At each anchor:
      miss M = max( |first-arrival error| / 1 h ,  |peak error| / 15 % )         smaller is better; M <= 1 means both closing criteria are met on W1.
      eligible = the sweep points with M <= 1.
      chosen   = the 6 eligible points with the lowest M, plus the 2 eligible points with the highest KGE_2012 among the rest (so the best-shaped
                 runs are not lost).  Fewer than 8 eligible: take all of them.  Fewer than 5 chosen: add the points closest in arrival to reach 5.
      Ties are broken by the smaller sweep number k.  The baseline routing (cv 6, flowexp 0.3, n 0.04) is also run on W1B at each anchor, as the reference.
  Pass rule (W1B).  A run passes if it met both criteria on W1 (M <= 1) AND meets both on the W1B second storm (Feb 27-28: |arrival| <= 1 h, |peak| <= 15 %),
      with the same scorer definitions as W1.  Winner at each anchor = the passing run with the smallest worst-case miss max(M on W1, M on W1B); runs within
      0.05 of that best value are tied, and the tie goes to the smaller |gamma - 1| on W1.  NO LOOSENING after the fact: if nothing passes, the thresholds do
      not change; the result is reported and Tracy decides.  Shape (KGE gamma, KGE) is reported but is NOT a gate at this step.

THE LIST (frozen before any run, so no W1B result can influence it)
    --make-list reads rt1_summary.csv (your sweep summary) and rt1_list.csv, applies the selection rule, and writes {batch}_list.csv (default rt1b_list.csv):
    for each chosen point its routing values, the group it was chosen for, and its W1 numbers.  It is NEVER overwritten; its md5 is recorded in the manifest at
    the first run and checked on every start.  The values are cross-checked against rt1_list.csv (the frozen sweep list) before the list is written.

RUN 0.  At each anchor the baseline routing is run first through the same shortcut path.  If an earlier full-build W1B run of the same soil exists in results/
    (the Step A/B runs are named W1B_ks.._cv6_psib.._f.._d0NN_nocl), run 0 is compared with it byte for byte.  Identical = the W1B shortcut reproduces a full build.
    Not identical = a WARNING with the largest difference (the batch goes on).  No earlier run = a note (not an error).

SAFEGUARDS (inherited from the sweep script)
    heartbeat; shared w2.lock; free-disk guard; one tRIBS at a time; each edited .in checked line by line against the source .in (only cv, flowexp, n and the
    output name may differ) and against the read-only reference .in; the soil build checked against the SOLUS grids; outlet file checked; a fingerprint of the
    soil, bedrock and water-table files taken after each anchor's build and re-checked before EVERY run; list md5 checked on every start; a tRIBS failure on one
    point is recorded and the batch goes on, 3 failures in a row stop it; on a restart the anchor's soil is rebuilt (about 1 min) and finished runs are never rerun.
    NOTHING IS DELETED except a run's own half-finished files when that run is restarted.

HOW TO RUN (from workspaces/Sycamore/lab/; keep run_w4_draw12_sycamore.py, run_trial_inedit_sycamore.py, run_sweep_routing_sycamore.py and this file together,
and keep rt1_summary.csv and rt1_list.csv there):
    python -u run_w1b_routing_sycamore.py --make-list          # once: pick the points by the rule, freeze the list, print it
    python -u run_w1b_routing_sycamore.py --dry-run            # checks only, builds nothing
    python -u run_w1b_routing_sycamore.py --limit 2            # first check: 1 soil build + baseline run 0 + point 1 (about 5 min)
    python -u run_w1b_routing_sycamore.py                      # everything (about 55 min), then the pass table
    python -u run_w1b_routing_sycamore.py --summary-only
If the connection drops, run the same command again: finished runs are skipped.  Do NOT run the notebook (or anything that rebuilds data/model/) while this runs.

Files it makes (in lab/): {batch}_list.csv, {batch}_manifest.json, {batch}_logs/, score_out_{batch}/scores_{batch}.csv, {batch}_summary.csv, w2.lock (only while running),
and per run: <run>.in, simulation_<run>.log, results/<run>*.  Run names: W1B_{batch}_<soil>_r<k>_cv.._fe.._n..; source builds W1BSRC_{batch}_<soil>.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import run_w4_draw12_sycamore as base      # the tested W4 runner: build, check, tRIBS functions
import run_trial_inedit_sycamore as trial  # the tested shortcut helpers
import run_sweep_routing_sycamore as sw    # the routing sweep: edit, fingerprint and list-reading helpers

# --------------------------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------------------------
BATCH_DEFAULT = "rt1b"
SOURCE_BATCH = "rt1"                                  # the W1 sweep this batch is chosen from
SUMMARY_DEFAULT = "rt1_summary.csv"
ANCHORS_DEFAULT = "12,25"
W1B = dict(name="W1B", start="02/21/2017/00/00", hours=240, label="feb21_to_mar3", gw=0.82, scorer_start="2017-02-21 00:00")
N_LOW, N_KGE, MIN_TOTAL = 6, 2, 5                     # selection rule: lowest-miss points, best-KGE extras, minimum chosen
CLOSE_ARRIVAL_H, CLOSE_PEAK_PCT = 1.0, 15.0           # the closing criteria
TIE = 0.05                                            # winners within this of the best worst-case miss are tied
OFFSET = 0.0104                                       # Step A: Feb 27-28 volume ratio on W1B is lower than on W2 by this (information only)
MB_PER_RUN = 25.0                                     # W1B results files per run (W1 was about 17 MB for 168 h)
MAX_CONSECUTIVE_FAILS = 3
SELECTION_RULE = "6 lowest miss + 2 best KGE among eligible (miss<=1); min 5; baseline run 0 as reference"
GROUPS = ("low_miss", "best_kge", "padded_arrival", "baseline")
W1_ARR, W1_PK, W1_VOL = "flood_week__arrival_err_h", "flood_week__peak_err_pct", "flood_week__vol_ratio"
W1_GAM, W1_KGE, W1_R = "flood_week__kge_gamma", "flood_week__kge2012", "flood_week__kge_r"
B_ARR, B_PK, B_VOL = "feb27_28__arrival_err_h", "feb27_28__peak_err_pct", "feb27_28__vol_ratio"
B_GAM, B_KGE, B_R, B_FLAG = "feb27_28__kge_gamma", "feb27_28__kge2012", "feb27_28__kge_r", "feb27_28__arrival_flag"
OUTLET = sw.OUTLET
KEYS = sw.KEYS
LIST_COLS = ["anchor", "order", "k", "cv", "flowexp", "n", "group", "w1_arrival_h", "w1_peak_pct", "w1_miss", "w1_gamma", "w1_kge", "w1_vol",
             "summary_md5", "source_list_md5", "selection_rule"]


def log(msg: str = "") -> None:
    base.log(msg)


class RunError(base.RunError):
    pass


def nan() -> float:
    return float("nan")


def finite(x) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(x)


def miss(arr: float, pk: float) -> float:
    """Worst-case miss against the closing criteria: 1.0 = exactly on a limit; infinity if either number is missing."""
    if not (finite(arr) and finite(pk)):
        return float("inf")
    return max(abs(arr) / CLOSE_ARRIVAL_H, abs(pk) / CLOSE_PEAK_PCT)


def fmt_num(x: float, nd: int = 5) -> str:
    return "nan" if not finite(x) else f"{x:.{nd}g}"


# --------------------------------------------------------------------------------------------
# The selection rule (pure function: easy to test)
# --------------------------------------------------------------------------------------------
def select_points(rows: list) -> list:
    """rows: dicts with k, arr, pk, kge (W1 numbers) for the sweep points of ONE anchor. Returns [(row, group)] in run order."""
    for r in rows:
        r["miss"] = miss(r["arr"], r["pk"])
    eligible = sorted([r for r in rows if r["miss"] <= 1.0], key=lambda r: (r["miss"], r["k"]))
    low = eligible[:N_LOW]
    rest = eligible[N_LOW:]
    by_kge = sorted(rest, key=lambda r: (-(r["kge"] if finite(r["kge"]) else -1e9), r["k"]))
    extra = by_kge[:N_KGE]
    chosen = [(r, "low_miss") for r in low] + [(r, "best_kge") for r in extra]
    if len(chosen) < MIN_TOTAL:
        taken = {r["k"] for r, _ in chosen}
        pool = sorted([r for r in rows if r["k"] not in taken and finite(r["arr"])], key=lambda r: (abs(r["arr"]), r["k"]))
        for r in pool[: MIN_TOTAL - len(chosen)]:
            chosen.append((r, "padded_arrival"))
    return chosen


# --------------------------------------------------------------------------------------------
# The list: make and freeze, read and check
# --------------------------------------------------------------------------------------------
def list_path(lab: Path, batch: str) -> Path:
    return lab / f"{batch}_list.csv"


def routing_text(key: str, v: float) -> str:
    return f"{v:.{sw.DECIMALS[key]}f}".rstrip("0").rstrip(".")


def make_list(args, lab: Path) -> int:
    path = list_path(lab, args.batch)
    if path.exists():
        print(f"{path.name} already exists and is NEVER overwritten (it is the frozen list). For a new list, give it a new --batch name.")
        return 1
    spath, lpath = lab / args.summary, lab / f"{SOURCE_BATCH}_list.csv"
    for p in (spath, lpath):
        if not p.is_file():
            print(f"{p.name} not found in lab/. It is needed to pick the points (the sweep summary and the frozen sweep list).")
            return 1
    try:
        src_rows = sw.read_list(lpath)
    except base.Abort as exc:
        print(str(exc))
        return 1
    src = {r["k"]: r for r in src_rows}
    with open(spath, newline="") as fh:
        table = list(csv.DictReader(fh))
    anchors = [int(a) for a in args.anchors.replace(" ", "").split(",") if a]
    if not anchors or len(set(anchors)) != len(anchors):
        print("--anchors is empty or names a draw twice.")
        return 1
    smd5, lmd5 = sw.file_md5(spath), sw.file_md5(lpath)
    out_rows, problems = [], []
    print(f"Selection rule: {SELECTION_RULE}")
    print(f"  source: {spath.name} (md5 {smd5}) and {lpath.name} (md5 {lmd5})")
    for a in anchors:
        mine = [r for r in table if r.get("anchor", "").strip() == str(a)]
        sweep, base0 = [], None
        for r in mine:
            if r.get("status") != "done":
                continue
            rec = dict(k=int(r["k"]), cv=float(r["cv"]), fe=float(r["flowexp"]), n=float(r["n"]), arr=base.fnum(r, W1_ARR), pk=base.fnum(r, W1_PK),
                       kge=base.fnum(r, W1_KGE), gam=base.fnum(r, W1_GAM), vol=base.fnum(r, W1_VOL))
            if r.get("kind") == "sweep":
                sweep.append(rec)
            elif r.get("kind") == "baseline":
                base0 = rec
        if not sweep:
            problems.append(f"draw {a}: no scored sweep rows in {spath.name}")
            continue
        if base0 is None:
            problems.append(f"draw {a}: no scored baseline (run 0) row in {spath.name}")
            continue
        for rec in sweep:       # the summary's routing values must be the frozen list's values
            s = src.get(rec["k"])
            if s is None:
                problems.append(f"draw {a}, k={rec['k']}: not in {lpath.name}")
            elif any(abs(rec[x] - s[y]) > 1e-9 for x, y in (("cv", "KINEMVELCOEF"), ("fe", "FLOWEXP"), ("n", "CHANNELROUGHNESS"))):
                problems.append(f"draw {a}, k={rec['k']}: the summary's cv/flowexp/n differ from {lpath.name}")
        chosen = select_points(sweep)
        elig = sum(1 for r in sweep if r["miss"] <= 1.0)
        print(f"\n--- draw {a}: {len(sweep)} sweep points scored on W1; {elig} meet both closing criteria (miss <= 1); {len(chosen)} chosen ---")
        print(f"    {'order':>5} {'k':>3} {'group':<15} {'cv':>7} {'flowexp':>8} {'n':>7} {'W1 arr h':>9} {'W1 peak %':>10} {'miss':>6} {'gamma':>6} {'KGE':>6}")
        out_rows.append(dict(anchor=a, order=0, k=0, cv=base0["cv"], flowexp=base0["fe"], n=base0["n"], group="baseline", w1_arrival_h=base0["arr"],
                             w1_peak_pct=base0["pk"], w1_miss=miss(base0["arr"], base0["pk"]), w1_gamma=base0["gam"], w1_kge=base0["kge"], w1_vol=base0["vol"]))
        for i, (r, g) in enumerate(chosen, 1):
            out_rows.append(dict(anchor=a, order=i, k=r["k"], cv=r["cv"], flowexp=r["fe"], n=r["n"], group=g, w1_arrival_h=r["arr"], w1_peak_pct=r["pk"],
                                 w1_miss=r["miss"], w1_gamma=r["gam"], w1_kge=r["kge"], w1_vol=r["vol"]))
            print(f"    {i:>5} {r['k']:>3} {g:<15} {r['cv']:>7g} {r['fe']:>8g} {r['n']:>7g} {r['arr']:>+9.2f} {r['pk']:>+10.1f} {r['miss']:>6.2f} {r['gam']:>6.3f} {r['kge']:>6.3f}")
        if len(chosen) < N_LOW + N_KGE:
            print(f"    note: only {len(chosen)} chosen (fewer than {N_LOW + N_KGE}): " + ("padded with the closest in arrival" if any(g == 'padded_arrival' for _, g in chosen) else "all eligible points taken"))
    if problems:
        print("\nNothing was written. Problems:\n  " + "\n  ".join(problems))
        return 1
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(LIST_COLS)
        for r in out_rows:
            w.writerow([r["anchor"], r["order"], r["k"], routing_text("KINEMVELCOEF", r["cv"]), routing_text("FLOWEXP", r["flowexp"]), routing_text("CHANNELROUGHNESS", r["n"]),
                        r["group"], fmt_num(r["w1_arrival_h"]), fmt_num(r["w1_peak_pct"]), fmt_num(r["w1_miss"]), fmt_num(r["w1_gamma"]), fmt_num(r["w1_kge"]),
                        fmt_num(r["w1_vol"]), smd5, lmd5, SELECTION_RULE])
    read_list(path)         # the saved file must pass the same checks the batch applies
    n_runs = len(out_rows)
    print(f"\nWrote and froze {path.name}: {n_runs} rows ({n_runs - len(anchors)} chosen points + {len(anchors)} baseline runs), md5 {sw.file_md5(path)}.")
    print(f"Next: python -u run_w1b_routing_sycamore.py --dry-run   (about {n_runs} runs, roughly {base.fmt_hm(n_runs * 2.7 * 60 + len(anchors) * 1.1 * 60)} in all)")
    return 0


def read_list(path: Path) -> list:
    if not path.is_file():
        raise base.Abort(f"{path.name} not found. Make it first: python -u run_w1b_routing_sycamore.py --make-list")
    rows, problems = [], []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh)
        miss_cols = [c for c in LIST_COLS if c not in (rd.fieldnames or [])]
        if miss_cols:
            raise base.Abort(f"{path.name} lacks the columns {miss_cols}.")
        for r in rd:
            try:
                row = dict(anchor=int(r["anchor"]), order=int(r["order"]), k=int(r["k"]), group=r["group"].strip(),
                           KINEMVELCOEF=float(r["cv"]), FLOWEXP=float(r["flowexp"]), CHANNELROUGHNESS=float(r["n"]),
                           w1_arr=float(r["w1_arrival_h"]), w1_pk=float(r["w1_peak_pct"]), w1_miss=float(r["w1_miss"]), w1_gamma=float(r["w1_gamma"]),
                           w1_kge=float(r["w1_kge"]), w1_vol=float(r["w1_vol"]), summary_md5=r["summary_md5"], source_list_md5=r["source_list_md5"])
            except ValueError:
                problems.append(f"unreadable row {dict(r)}")
                continue
            rows.append(row)
    seen = set()
    for r in rows:
        tag = f"draw {r['anchor']}, k={r['k']}"
        if (r["anchor"], r["k"]) in seen:
            problems.append(f"{tag}: appears twice")
        seen.add((r["anchor"], r["k"]))
        if r["group"] not in GROUPS:
            problems.append(f"{tag}: group {r['group']!r} is not one of {GROUPS}")
        is_base = tuple(r[k] for k in KEYS) == tuple(float(sw.BASE_ROUTING[k]) for k in KEYS)
        if (r["group"] == "baseline") != is_base or (r["group"] == "baseline") != (r["k"] == 0):
            problems.append(f"{tag}: the baseline row and the baseline routing values (cv 6, flowexp 0.3, n 0.04, k=0) must go together")
        for key in KEYS:
            v, (lo, hi) = r[key], sw.SANITY[key]
            if not (lo <= v <= hi):
                problems.append(f"{tag}: {sw.SHORT[key]} = {v:g} is outside the sanity limits {lo:g}-{hi:g}")
            if round(v, sw.DECIMALS[key]) != v:
                problems.append(f"{tag}: {sw.SHORT[key]} = {v!r} has more than {sw.DECIMALS[key]} decimals")
    by_anchor = {}
    for r in rows:
        by_anchor.setdefault(r["anchor"], []).append(r)
    for a, rs in by_anchor.items():
        if sum(1 for r in rs if r["group"] == "baseline") != 1:
            problems.append(f"draw {a}: needs exactly one baseline row")
        trips = [tuple(r[k] for k in KEYS) for r in rs]
        if len(set(trips)) != len(trips):
            problems.append(f"draw {a}: the same cv, flowexp, n appear twice")
        if [r["order"] for r in sorted(rs, key=lambda x: x["order"])] != list(range(len(rs))):
            problems.append(f"draw {a}: the order column is not 0, 1, 2, ...")
    if not rows:
        problems.append("the list is empty")
    if problems:
        raise base.Abort(f"{path.name} has problems:\n  " + "\n  ".join(problems[:12]))
    return rows


# --------------------------------------------------------------------------------------------
# Names and small helpers
# --------------------------------------------------------------------------------------------
def run_name(batch: str, draw: dict, k: int, vals: dict) -> str:
    return (f"W1B_{batch}_{sw.soil_tag(draw)}_r{k:03d}_cv{base.vtag(vals['KINEMVELCOEF'])}_fe{base.vtag(vals['FLOWEXP'])}"
            f"_n{base.vtag(vals['CHANNELROUGHNESS'])}")


def src_name(batch: str, draw: dict) -> str:
    return f"W1BSRC_{batch}_{sw.soil_tag(draw)}"


def earlier_w1b_id(draw: dict) -> str:
    """The earlier full-build W1B run of this soil at baseline routing (Step A/B naming), if one exists."""
    return "W1B_" + draw["run_id"][len("W4_"):]


def run_scorer_w1b(rid: str, scorer: str, score_dir: Path, log_dir: Path, write_checksum: bool = True) -> dict:
    """Score one W1B run. The run id starts with W1B_, which the scorer does not know as a window, so the start time is given (as the Step B runner did)."""
    cmd = [sys.executable, scorer, "score", "--run", rid, "--start", W1B["scorer_start"], "--results-dir", "results",
           "--out-dir", str(score_dir), "--csv-name", base.SCORE_CSV]
    if write_checksum:
        cmd.append("--write-checksum")
    res = subprocess.run(cmd, capture_output=True, text=True)
    (log_dir / f"{rid}_score.txt").write_text((res.stdout or "") + (res.stderr or ""))
    if res.returncode != 0:
        raise base.RunError(f"The scorer exited with code {res.returncode}: {(res.stderr or res.stdout)[-400:]}")
    csv_path = score_dir / base.SCORE_CSV
    if csv_path.is_file():
        with open(csv_path, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("run_id") == rid:
                    return dict(row)
    raise base.RunError(f"The scorer ran but {rid} is not in {csv_path}.")


def one_line(row: dict) -> str:
    g = lambda c: base.fnum(row, c)
    return (f"Feb27-28: arrival {base.fm(g(B_ARR), 2, True)} h  peak {base.fm(g(B_PK), 1, True)} %  vol {base.fm(g(B_VOL), 3)}  "
            f"KGE r {base.fm(g(B_R), 3)}  gamma {base.fm(g(B_GAM), 3)}  KGE {base.fm(g(B_KGE), 3)}")


# --------------------------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------------------------
def preflight(args, lab: Path, plan: list, rows: list, man: dict, ref_in: Path) -> list:
    out = []

    def add(level, msg):
        out.append((level, msg))

    for f in (args.notebook, "run_w4_draw12_sycamore.py", "run_trial_inedit_sycamore.py", "run_sweep_routing_sycamore.py", base.SCORER):
        add("OK" if (lab / f).is_file() else "FAIL", f"file in lab/: {f}")
    add("OK" if (os.path.isfile(args.tribs) and os.access(args.tribs, os.X_OK)) else "FAIL", f"tRIBS program: {args.tribs}")
    n_new = sum(len(a["pending"]) for a in plan)
    need = max(args.min_free_gb, n_new * MB_PER_RUN / 1000.0 + 0.5)
    free = base.disk_free_gb(lab)
    add("OK" if free >= need else "FAIL", f"free disk space: {free:.1f} GB (need about {need:.1f}: {n_new} runs x {MB_PER_RUN:g} MB + margin)")
    comp = base.competing_tribs()
    add("FAIL" if comp else "OK", "no other tRIBS running" + (": " + comp[0] if comp else ""))
    for f in ("../init_data/bedrock_depth_m.asc", "../init_data/scgrid.gdf", "../init_data/soils.sdt", "../init_data/met/Master_Met.sdf",
              "../init_data/LandUse.asc", "../init_data/obs/sycamore_09510200_obs_flow_WY2017.csv"):
        add("OK" if (lab / f).is_file() else "FAIL", f"input: {f}")
    try:
        got = base.read_in_file(ref_in)
        need_kw = list(base.ROUTING_KWS) + ["OPTPERCOLATION", "OPTSNOW", "RUNTIME", "STARTDATE", "GWATERFILE", "BEDROCKFILE", "OUTFILENAME"]
        missing = [k for k in need_kw if k not in got]
        bad = [f"{k} is {got.get(k)!r}, wanted {v!r}" for k, v in sw.BASE_ROUTING.items() if k in got and not base.same_value(got[k], v)]
        add("FAIL" if (missing or bad) else "OK", f"reference .in ({ref_in.name}, the W2 baseline's; read only) has the needed keywords and baseline routing"
            + (f" (missing: {missing})" if missing else "") + (f" ({'; '.join(bad)})" if bad else ""))
    except Exception as exc:
        add("FAIL", f"reference .in {ref_in.name}: {exc}")
    lp = list_path(lab, args.batch)
    md5 = sw.file_md5(lp)
    if man.get("list_md5") and man["list_md5"] != md5:
        add("FAIL", f"{lp.name} has changed since the first start (md5 {md5} vs recorded {man['list_md5']}). The list is frozen: restore it, or use a new --batch name.")
    else:
        add("OK", f"W1B list {lp.name}: {len(rows)} valid rows, md5 {md5}" + ("" if man.get("list_md5") else " (will be recorded at the first start)"))
    src_list = lab / f"{SOURCE_BATCH}_list.csv"
    smd5s = {r["source_list_md5"] for r in rows}
    if src_list.is_file():
        now = sw.file_md5(src_list)
        add("OK" if smd5s == {now} else "WARN", f"{src_list.name} md5 {now} " + ("matches" if smd5s == {now} else f"differs from the one the list was made from ({', '.join(sorted(smd5s))})"))
    else:
        add("WARN", f"{src_list.name} not found; its md5 could not be compared with the one the list was made from")
    for a in plan:
        if a["draw"]["idx"] not in {r["anchor"] for r in rows}:
            add("FAIL", f"draw {a['draw']['idx']} is not in the list")
    try:
        import nbformat
        nb = nbformat.read(str(lab / args.notebook), as_version=4)
        idx = base.find_settings_cell(nb)
        for a in plan:
            sw.set_routing(sw.BASE_ROUTING)
            base.patched_source(nb.cells[idx].source, sw.make_job(a["draw"], a["src_id"]))
        add("OK", f"notebook Settings cell (cell {idx}) patches correctly for each anchor's soil build on the W1B dates ({W1B['start']}, {W1B['hours']} h)")
    except ImportError:
        add("FAIL", "nbformat is not installed: pip install nbformat nbclient ipykernel")
    except Exception as exc:
        add("FAIL", f"Settings cell check: {exc}")
    try:
        res = subprocess.run([sys.executable, str(lab / args.scorer), "score", "--help"], capture_output=True, text=True)
        flags = ["--run", "--start", "--results-dir", "--out-dir", "--csv-name", "--write-checksum"]
        miss_f = [f for f in flags if f not in res.stdout]
        add("FAIL" if (res.returncode != 0 or miss_f) else "OK", "scorer accepts " + " ".join(flags) + (f" (missing: {miss_f})" if miss_f else ""))
    except Exception as exc:
        add("FAIL", f"scorer check: {exc}")
    if not args.skip_kernel_check:
        try:
            from nbformat.v4 import new_code_cell, new_notebook
            nb = new_notebook(cells=[new_code_cell("import sys, os\nfrom pytRIBS.classes import *\nprint('KERNEL', sys.executable, os.getcwd())")])
            base.execute_nb(nb, lab, 300)
            add("OK", "notebook kernel imports pytRIBS")
        except Exception as exc:
            add("FAIL", f"notebook kernel could not import pytRIBS: {type(exc).__name__}: {str(exc)[:300]}")
    taken_bad = []
    for a in plan:
        for rid in [a["src_id"]] + [p["rid"] for p in a["pending"]]:
            files = [p.name for p in (lab / "results").glob(rid + "*")] if (lab / "results").is_dir() else []
            files += [p.name for p in lab.glob(rid + "*")]
            known = rid in man["runs"] or rid == a["src_id"]
            if files and not known:
                taken_bad.append(f"{rid} ({files[0]})")
    add("FAIL" if taken_bad else "OK", "run names unused" + (": already exist and are not in the manifest: " + "; ".join(taken_bad[:3]) if taken_bad else ""))
    for a in plan:
        d = a["draw"]
        e = lab / "results" / f"{earlier_w1b_id(d)}{OUTLET}"
        add("OK" if e.is_file() else "NOTE", f"draw {d['idx']}: " + (f"earlier full-build W1B result {e.name} found; run 0 will be compared with it byte for byte"
                                                                    if e.is_file() else "no earlier full-build W1B result of this soil in results/; run 0 will not be compared (not an error)"))
    return out


# --------------------------------------------------------------------------------------------
# One anchor, one run
# --------------------------------------------------------------------------------------------
def build_anchor(a: dict, args, lab: Path, man: dict, ref_in: Path) -> Path:
    """Build the anchor's soil at baseline routing on the W1B dates (no tRIBS), check it, take the fingerprint. Returns the source .in path."""
    d, sid = a["draw"], a["src_id"]
    log_dir = lab / base.LOG_DIR
    sw.set_routing(sw.BASE_ROUTING)
    job = sw.make_job(d, sid)
    base.remove_stale(lab, sid)
    log(f"=== anchor draw {d['idx']}: building its soil once for W1B (Ks x{d['ks']:g}, f x{d['f']:g}, PsiB x{d['psib']:g}) as '{sid}'; no tRIBS run ===")
    t = time.time()
    text, warns = base.build_run(job, args, lab, log_dir)
    src_in = lab / f"{sid}.in"
    base.check_in_file(src_in, ref_in, job)
    grid_note = base.check_grids(lab, job, t)
    fp = sw.soil_fingerprint(lab, src_in)
    man["anchors"][str(d["idx"])] = dict(src_id=sid, built=base.now_iso(), build_min=round((time.time() - t) / 60, 2), fingerprint=fp, grid_check=grid_note)
    base.save_manifest(man)
    log(f"build OK in {base.fmt_min(time.time() - t)}; .in matches the reference except the allowed lines (W1B start and length included). {grid_note}")
    for w in warns:
        log(f"notebook said: {w}")
    return src_in


def run_one(a: dict, p: dict, src_in: Path, args, lab: Path, man: dict, ref_in: Path) -> None:
    """One routing point on W1B: edit the .in, run tRIBS, check the outlet, score. RunError = this run failed (recorded, batch goes on)."""
    d, rid, vals, k = a["draw"], p["rid"], p["vals"], p["k"]
    score_dir, log_dir = lab / base.SCORE_DIR, lab / base.LOG_DIR
    prob = sw.fingerprint_problem(lab, src_in, man["anchors"][str(d["idx"])]["fingerprint"])
    if prob:
        raise base.Abort(f"The soil build was disturbed before {rid}: {prob.rstrip('.')}. Was the notebook run? Nothing was started; run the same command again "
                         f"(it rebuilds the anchor's soil first).")
    free = base.disk_free_gb(lab)
    if free < args.min_free_gb:
        raise base.Abort(f"Only {free:.1f} GB free; need {args.min_free_gb:g} GB. Free some space and run again to resume.")
    comp = base.competing_tribs()
    if comp:
        raise base.Abort("Another tRIBS is running, so the batch stops (two at once is not allowed): " + comp[0])
    rec = man["runs"].setdefault(rid, {})
    rec.clear()
    rec.update(status="running", started=base.now_iso(), draw=d["idx"], k=k, group=p["group"], cv=sw.num(vals["KINEMVELCOEF"]), flowexp=sw.num(vals["FLOWEXP"]),
               n=sw.num(vals["CHANNELROUGHNESS"]), kind=p["kind"])
    base.save_manifest(man)
    base.remove_stale(lab, rid)
    new_in = lab / f"{rid}.in"
    t0 = time.time()
    sw.write_edited_in(src_in, new_in, vals, rid)
    sw.set_routing(vals)
    job = sw.make_job(d, rid)
    try:
        base.check_in_file(new_in, ref_in, job)
    except base.RunError as exc:
        raise base.Abort(f"{rid}: the edited .in does not match the reference: {exc}")
    base.install_term_handler()
    rc, secs = base.run_tribs(args.tribs, new_in.name, f"simulation_{rid}.log", base.timeout_for(args),
                              progress_file=lab / "results" / f"{rid}{OUTLET}", hours=W1B["hours"],
                              heartbeat_min=args.heartbeat_min, label=rid)
    rec.update(tribs_min=round(secs / 60, 2), tribs_returncode=rc)
    if rc != 0:
        raise RunError(f"tRIBS exited with code {rc}; see simulation_{rid}.log")
    outlet = base.check_outlet(lab / "results" / f"{rid}{OUTLET}", W1B["hours"])
    rec.update(outlet=outlet, status="unscored", output_mb=base.output_mb(lab, rid))
    base.save_manifest(man)
    if p["kind"] == "baseline":
        ref = lab / "results" / f"{earlier_w1b_id(d)}{OUTLET}"
        if ref.is_file():
            same = sw.md5_or_none(ref) == outlet["md5"]
            rec["earlier_w1b_identical"] = same
            if same:
                log(f"{rid}: run 0 is BYTE-IDENTICAL to the earlier full-build W1B result ({ref.name}): the W1B shortcut path reproduces a full build.")
            else:
                gap = trial.outlet_gap(lab, rid, earlier_w1b_id(d))
                rec["earlier_w1b_gap"] = gap
                log(f"{rid}: WARNING: run 0 is NOT identical to the earlier W1B result ({ref.name}): {gap}. The batch goes on, but find out why before trusting "
                    f"comparisons with earlier W1B numbers (the earlier run may have used an older notebook or settings).")
        else:
            log(f"{rid}: no earlier full-build W1B result of this soil to compare run 0 with (not an error).")
    try:
        row = run_scorer_w1b(rid, str(lab / args.scorer), score_dir, log_dir)
        rec["scores"] = base.show_scores(row)
        rec["status"] = "done"
        log(f"{rid}: {base.fmt_min(time.time() - t0)}; {one_line(row)}")
    except base.RunError as exc:
        rec["score_error"] = str(exc)
        log(f"{rid}: scoring failed (the run itself is fine): {exc}")
    rec["finished"] = base.now_iso()
    base.save_manifest(man)


# --------------------------------------------------------------------------------------------
# The pass table
# --------------------------------------------------------------------------------------------
def load_scores(lab: Path) -> dict:
    p = lab / base.SCORE_DIR / base.SCORE_CSV
    if not p.is_file():
        return {}
    with open(p, newline="") as fh:
        return {r["run_id"]: r for r in csv.DictReader(fh)}


def evaluate(p: dict, row: dict) -> dict:
    """Everything the pass table needs for one run: its frozen W1 numbers and its W1B scores."""
    w1 = p["w1"]
    ba, bp = base.fnum(row, B_ARR), base.fnum(row, B_PK)
    mb = miss(ba, bp)
    m1 = w1["miss"]
    return dict(anchor=p["anchor"], k=p["k"], kind=p["kind"], group=p["group"], cv=sw.num(p["vals"]["KINEMVELCOEF"]), flowexp=sw.num(p["vals"]["FLOWEXP"]),
                n=sw.num(p["vals"]["CHANNELROUGHNESS"]), run_id=p["rid"], w1_arrival_h=w1["arr"], w1_peak_pct=w1["pk"], w1_miss=m1, w1_gamma=w1["gamma"], w1_kge=w1["kge"],
                b_arrival_h=ba, b_peak_pct=bp, b_vol=base.fnum(row, B_VOL), b_vol_w2eq=base.fnum(row, B_VOL) + OFFSET, b_gamma=base.fnum(row, B_GAM),
                b_kge=base.fnum(row, B_KGE), b_kge_r=base.fnum(row, B_R), b_arrival_flag=row.get(B_FLAG, ""), b_miss=mb,
                passes=bool(m1 <= 1.0 and mb <= 1.0), worst=max(m1, mb))


def pick_winner(evals: list):
    """Approved rule: among passing sweep points, smallest worst-case miss; within TIE of the best, the smaller |gamma - 1| on W1."""
    ok = [e for e in evals if e["kind"] == "sweep" and e["passes"]]
    if not ok:
        return None
    best = min(e["worst"] for e in ok)
    tied = [e for e in ok if e["worst"] <= best + TIE]
    return min(tied, key=lambda e: (abs(e["w1_gamma"] - 1.0) if finite(e["w1_gamma"]) else 9e9, e["worst"], e["k"]))


def print_summary(man: dict, lab: Path, plan_all: list, args) -> None:
    scores = load_scores(lab)
    all_rows, cols = [], []
    log()
    log("=" * 110)
    log(f"W1B SECOND-STORM CHECK, batch {args.batch} (Feb 27-28 storm on the Feb 21 - Mar 3 window, Kchan off).")
    log(f"Pass rule: met both criteria on W1 (miss <= 1) AND meets both on W1B (|arrival| <= {CLOSE_ARRIVAL_H:g} h, |peak| <= {CLOSE_PEAK_PCT:g}%). Miss = max(|arrival|/{CLOSE_ARRIVAL_H:g} h, |peak|/{CLOSE_PEAK_PCT:g}%).")
    log("Winner = smallest worst-case miss over W1 and W1B (ties within 0.05: smaller |gamma-1| on W1). Nothing is promoted here; Tracy decides.")
    for a in plan_all:
        d = a["draw"]
        evals = []
        for p in a["all"]:
            r = man["runs"].get(p["rid"])
            if not r or r.get("status") != "done":
                continue
            evals.append(evaluate(p, scores.get(p["rid"], {})))
        n_planned = len(a["all"])
        failed = [p for p in a["all"] if (man["runs"].get(p["rid"]) or {}).get("status") == "failed"]
        log("")
        log(f"--- anchor draw {d['idx']} (Ks x{d['ks']:g}, f x{d['f']:g}, PsiB x{d['psib']:g}): {len(evals)} scored of {n_planned} planned; {len(failed)} failed ---")
        if not evals:
            continue
        log(f"    {'k':>3} {'group':<15} {'cv':>7} {'flowexp':>8} {'n':>7} | {'W1 arr':>7} {'W1 pk%':>7} {'miss':>5} | {'W1B arr':>8} {'W1B pk%':>8} {'miss':>5} {'vol':>6} {'gamma':>6} {'KGE':>6} | {'worst':>5} pass")
        order = sorted(evals, key=lambda e: (e["kind"] != "baseline", e["worst"] if finite(e["worst"]) else 9e9, e["k"]))
        for e in order:
            tag = "baseline" if e["kind"] == "baseline" else ("PASS" if e["passes"] else "-")
            log(f"    {e['k']:>3} {e['group']:<15} {e['cv']:>7g} {e['flowexp']:>8g} {e['n']:>7g} | {base.fm(e['w1_arrival_h'], 2, True):>7} {base.fm(e['w1_peak_pct'], 1, True):>7} "
                f"{base.fm(e['w1_miss'], 2):>5} | {base.fm(e['b_arrival_h'], 2, True):>8} {base.fm(e['b_peak_pct'], 1, True):>8} {base.fm(e['b_miss'], 2):>5} "
                f"{base.fm(e['b_vol'], 3):>6} {base.fm(e['b_gamma'], 3):>6} {base.fm(e['b_kge'], 3):>6} | {base.fm(e['worst'], 2):>5} {tag}")
            if e["b_arrival_flag"] not in ("", "ok"):
                log(f"        note: no W1B arrival time for k={e['k']}: {e['b_arrival_flag']}")
        sweep = [e for e in evals if e["kind"] == "sweep"]
        npass = sum(1 for e in sweep if e["passes"])
        log(f"    {npass} of {len(sweep)} scored sweep points pass the rule on W1B.")
        b0 = next((e for e in evals if e["kind"] == "baseline"), None)
        if b0:
            log(f"    baseline routing on W1B: arrival {base.fm(b0['b_arrival_h'], 2, True)} h, peak {base.fm(b0['b_peak_pct'], 1, True)} % (miss {base.fm(b0['b_miss'], 2)}).")
        win = pick_winner(evals)
        if win:
            log(f"    WINNER by the approved rule: k={win['k']} (cv {win['cv']:g}, flowexp {win['flowexp']:g}, n {win['n']:g}): worst-case miss {win['worst']:.2f}, "
                f"W1 gamma {win['w1_gamma']:.3f}. Promotion to W4 is Tracy's decision.")
        elif sweep:
            log("    NO run passes. By the rule the thresholds are NOT loosened; report this and Tracy decides what comes next.")
        for p in failed[:5]:
            log(f"    FAILED k={p['k']} {p['rid']}: {(man['runs'][p['rid']].get('error') or '')[:160]}")
        for e in evals:
            e2 = dict(e)
            e2["winner"] = bool(win and win["run_id"] == e["run_id"])
            all_rows.append(e2)
            for c in e2:
                if c not in cols:
                    cols.append(c)
    if all_rows:
        path = lab / f"{args.batch}_summary.csv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(all_rows)
        log("")
        log(f"WROTE {path.name} ({len(all_rows)} rows: W1 numbers from the frozen list, W1B scores, miss, pass, winner flag).")
        log(f"The Feb 27-28 volume ratio is also given as W2-equivalent (+{OFFSET}, the Step A shift); routing is not expected to change volume.")


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="W1B second-storm check of the routing sweep (see the top of this file).")
    ap.add_argument("--batch", default=BATCH_DEFAULT, help=f"name of this batch (default {BATCH_DEFAULT}); names the list, manifest, scores and runs")
    ap.add_argument("--make-list", action="store_true", help="pick the points by the approved rule from the sweep summary, freeze the list, and stop")
    ap.add_argument("--summary", default=SUMMARY_DEFAULT, help=f"(--make-list) the W1 sweep summary (default {SUMMARY_DEFAULT})")
    ap.add_argument("--anchors", default=ANCHORS_DEFAULT, help=f"draw numbers in run order (default {ANCHORS_DEFAULT}); draws come from lhs_volume_list.csv")
    ap.add_argument("--dry-run", action="store_true", help="check everything, show the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the pass table from the runs already scored, then stop")
    ap.add_argument("--limit", type=int, default=None, help="start at most N tRIBS runs in this call (run 0 counts); finished runs do not count")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done (you almost never want this)")
    ap.add_argument("--timeout-min", type=float, default=None)
    ap.add_argument("--cell-timeout", type=int, default=1800)
    ap.add_argument("--min-free-gb", type=float, default=2.0)
    ap.add_argument("--heartbeat-min", type=float, default=base.HEARTBEAT_MIN)
    ap.add_argument("--notebook", default=None)
    ap.add_argument("--scorer", default=base.SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (base.TRIBS_DEFAULT if os.path.isfile(base.TRIBS_DEFAULT)
                                                                       else shutil.which("tRIBS") or base.TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    if args.make_list:
        return make_list(args, lab)

    b = args.batch
    base.LOG_DIR, base.SCORE_DIR, base.SCORE_CSV, base.MANIFEST = f"{b}_logs", f"score_out_{b}", f"scores_{b}.csv", f"{b}_manifest.json"
    (lab / base.LOG_DIR).mkdir(exist_ok=True)
    base._LOGFILE = str(lab / base.LOG_DIR / "batch_log.txt")
    args.redo_ids = {t.strip() for t in (args.redo or "").split(",") if t.strip()}
    args.notebook, nb_note = base.pick_notebook(lab, args.notebook)

    try:
        rows = read_list(list_path(lab, b))
        wanted = [int(x) for x in args.anchors.replace(" ", "").split(",") if x]
        in_list = []
        for r in rows:
            if r["anchor"] not in in_list:
                in_list.append(r["anchor"])
        anchors = [a for a in wanted if a in in_list]
        if not anchors or len(set(wanted)) != len(wanted) or len(anchors) != len(wanted):
            raise base.Abort(f"--anchors {args.anchors!r} must name draws that are in the list ({in_list}), each once.")
        # The draws are read FIRST, while the imported module still has its W4 settings (its checks name the design run with the W4 prefix).
        base.EXPECTED_DRAWS.update(sw.EXPECTED_DRAWS)
        draws = base.read_anchors(lab, ",".join(str(a) for a in anchors))
    except base.Abort as exc:
        print(str(exc))
        return 1
    base.WIN = dict(W1B)                                  # now switch to W1B: dates, length, cold-start water table
    base.BASE["GW_FRACTION"] = W1B["gw"]
    sw.set_routing(sw.BASE_ROUTING)
    man = base.load_manifest()
    man.setdefault("runs", {})
    man.setdefault("anchors", {})

    plan_all = []
    for a_idx, d in zip(anchors, draws):
        mine = sorted([r for r in rows if r["anchor"] == a_idx], key=lambda r: r["order"])
        a = dict(draw=d, src_id=src_name(b, d), all=[], pending=[])
        for r in mine:
            vals = {k: r[k] for k in KEYS}
            w1 = dict(arr=r["w1_arr"], pk=r["w1_pk"], miss=r["w1_miss"], gamma=r["w1_gamma"], kge=r["w1_kge"])
            p = dict(anchor=a_idx, k=r["k"], kind="baseline" if r["group"] == "baseline" else "sweep", group=r["group"], vals=vals, w1=w1)
            p["rid"] = run_name(b, d, r["k"], vals)
            a["all"].append(p)
        for p in a["all"]:
            st = (man["runs"].get(p["rid"]) or {}).get("status")
            if st != "done" or p["rid"] in args.redo_ids:
                a["pending"].append(p)
        plan_all.append(a)
    ref_in = lab / f"{sw.REF_ID}.in"

    if args.summary_only:
        print_summary(man, lab, plan_all, args)
        return 0

    log(f"Working folder: {lab}")
    log(f"Notebook: {args.notebook}" + (f"  (note: {nb_note})" if nb_note else ""))
    log(f"Batch {b}: W1B list {list_path(lab, b).name} ({len(rows)} rows); window W1B {W1B['start']}, {W1B['hours']} h, water table {W1B['gw']} x bedrock, Kchan off; "
        f"scored on the Feb 27-28 storm.")
    for a in plan_all:
        d = a["draw"]
        n_sweep = sum(1 for p in a["all"] if p["kind"] == "sweep")
        log(f"  anchor draw {d['idx']}: Ks x{d['ks']:g}, f x{d['f']:g}, PsiB x{d['psib']:g}: {len(a['all'])} planned runs (run 0 + {n_sweep} chosen points), {len(a['pending'])} still to do")
    n_todo = sum(len(a["pending"]) for a in plan_all)
    n_builds = sum(1 for a in plan_all if a["pending"])
    n_this = n_todo if args.limit is None else min(n_todo, args.limit)
    est = n_this * 2.7 + n_builds * 1.1
    log(f"  Time: about {base.fmt_hm(est * 60)} ({n_this} runs x about 2.7 min + {n_builds} soil build(s) x about 1.1 min).")

    results = preflight(args, lab, plan_all, rows, man, ref_in)
    log()
    for level, msg in results:
        log(f"  {level:5s} {msg}")
    fails = [m for lv, m in results if lv == "FAIL"]
    log()
    if args.dry_run:
        log("DRY RUN: nothing was built or run. " + ("All checks passed." if not fails else f"{len(fails)} check(s) FAILED."))
        return 1 if fails else 0
    if fails:
        log(f"{len(fails)} check(s) FAILED, so nothing was started. Fix them and run again.")
        return 1
    if n_todo == 0:
        log("Every planned run is already done. Printing the pass table.")
        print_summary(man, lab, plan_all, args)
        return 0

    try:
        base.take_lock()
    except base.Abort as exc:
        log(str(exc))
        return 1

    man.setdefault("list_md5", sw.file_md5(list_path(lab, b)))
    man["batch"], man["script"], man["window"] = b, "run_w1b_routing_sycamore.py", W1B
    man.setdefault("started", base.now_iso())
    base.save_manifest(man)
    code, fails_in_row, n_started = 0, 0, 0
    t_batch = time.time()
    try:
        for a in plan_all:
            if not a["pending"]:
                continue
            if args.limit is not None and n_started >= args.limit:
                break
            # runs that ran earlier but were never scored: score them (no tRIBS)
            for p in list(a["pending"]):
                r = man["runs"].get(p["rid"]) or {}
                if r.get("status") == "unscored" and p["rid"] not in args.redo_ids:
                    try:
                        row = run_scorer_w1b(p["rid"], str(lab / args.scorer), lab / base.SCORE_DIR, lab / base.LOG_DIR)
                        r["scores"], r["status"] = base.show_scores(row), "done"
                        log(f"{p['rid']}: ran earlier, scored now: {one_line(row)}")
                        a["pending"].remove(p)
                        base.save_manifest(man)
                    except base.RunError as exc:
                        log(f"{p['rid']}: scoring failed again: {exc}")
            if not a["pending"]:
                continue
            src_in = build_anchor(a, args, lab, man, ref_in)
            for p in a["pending"]:
                if args.limit is not None and n_started >= args.limit:
                    log(f"--limit {args.limit} reached; stopping. Run the same command again to continue.")
                    break
                n_started += 1
                try:
                    run_one(a, p, src_in, args, lab, man, ref_in)
                    fails_in_row = 0
                except base.RunError as exc:                       # this run failed (tRIBS error or timeout, bad outlet file); the others may be fine
                    fails_in_row += 1
                    man["runs"].setdefault(p["rid"], {}).update(status="failed", error=str(exc), finished=base.now_iso())
                    base.save_manifest(man)
                    log(f"{p['rid']}: FAILED ({exc}). Recorded; the batch goes on.")
                    if fails_in_row >= MAX_CONSECUTIVE_FAILS:
                        raise base.Abort(f"{MAX_CONSECUTIVE_FAILS} runs failed in a row, so the batch stops. Look at the last simulation_*.log files before restarting.")
                done_n = sum(1 for x in man["runs"].values() if x.get("status") == "done")
                left = sum(len([q for q in a2["pending"] if (man["runs"].get(q["rid"]) or {}).get("status") != "done"]) for a2 in plan_all)
                per = (time.time() - t_batch) / n_started
                log(f"    progress: {done_n} runs done in total; about {left} left; roughly {base.fmt_hm(per * left)} to go at this pace")
        log("Batch call finished.")
    except base.Abort as exc:
        log(f"BATCH STOPPED: {exc}")
        code = 1
    except KeyboardInterrupt:
        for rid, r in man["runs"].items():
            if r.get("status") == "running":
                r["status"] = "interrupted"
        base.save_manifest(man)
        log("BATCH INTERRUPTED. The lock is removed. Run the same command again to continue; the interrupted run is restarted from its own files.")
        code = 130
    finally:
        sw.set_routing(sw.BASE_ROUTING)
    print_summary(man, lab, plan_all, args)
    return code


if __name__ == "__main__":
    sys.exit(main())
