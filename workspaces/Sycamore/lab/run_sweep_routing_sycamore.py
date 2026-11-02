#!/usr/bin/env python3
"""
run_sweep_routing_sycamore.py  --  Stage 5, Steps 3 and 4: the ROUTING SWEEP (cv, flowexp, n) on W1 at the soil anchors, using the
.in-edit shortcut that Step 2 proved is exactly equal to a full build.

What it does.  For each soil anchor (draw 12 first, then draw 25) it builds the anchor's soil ONCE, then for every routing point on the
frozen list it copies the .in file, changes only KINEMVELCOEF (cv), FLOWEXP and CHANNELROUGHNESS (n) and the output name, runs tRIBS on
W1 (the flood week, Jan 18-25 2017, 168 h, cold start, water table 0.82 x bedrock, Kchan off) and scores the run with
score_run_sycamore.py. About 1.7 min per run (Step 2: 1.72 min with the shortcut vs 2.90 min with a full build).

Two modes besides the sweep:
    --make-list    makes and FREEZES the routing list (the 60 sampled points) and stops. Run it once, before any sweep run.
    --summary-only prints the tables from what is already scored.

THE ROUTING LIST (frozen before any run, so no result can influence it)
    * Sampling is the project's own sampler, copied unchanged from run_lhs_Ks_f_120_tribs6.py: "stratified-independent" sampling. Each
      parameter's range is split into n equal slices, one uniform draw is taken in each slice, and each parameter's draws are then
      shuffled independently. It is NOT a true orthogonal Latin hypercube (the project says so itself). Series 121 and 122 used the same
      method, on a LINEAR scale for cv, flowexp and n. Same seed + same n + same box = identical draws.
    * Default box: the Stage 2 box: cv 4-10, flowexp 0.26-0.32, n 0.04-0.08. Default seed 501 (Series 121 used 42 and 201-205, Series 122
      used 301). Values are rounded (cv 3 decimals, flowexp and n 4 decimals) BEFORE the list is saved, so the saved number is exactly
      the number written into the .in file.
    * 60 points. The second anchor (draw 25) runs the FIRST 30 rows of the same list, so the two anchors see the same routing points.
      The first 30 rows are a random half of the 60 slices, not a 30-slice design, so they cover the box a bit less evenly (make-list
      prints how evenly). Never pool the shared points across anchors as if they were independent.
    * The list file {batch}_list.csv is never overwritten; its md5 is recorded in the manifest at the first run and checked every time.

EXTRA RUN 0 (a check, not a sweep point).  Before each anchor's sweep points, the script runs the BASELINE routing (cv 6, flowexp 0.3,
n 0.04) through the same shortcut path. Stage 3 already ran that exact soil + routing on W1, so the result should be byte-identical to
the Stage 3 W1 outlet file (results/W1_<soil>_nocl_Outlet.qout). If it is not, the script prints a WARNING with the largest difference;
it does not stop. (--no-baseline-run skips it.)

SAFEGUARDS (all inherited from the W4 and trial scripts, plus the new ones marked *)
    heartbeat line while tRIBS runs; shared w2.lock (cannot run beside any other runner); free-disk guard; one tRIBS at a time;
    each edited .in file checked line by line against the source .in (only cv, flowexp, n and the output name may differ) and against
    the read-only reference .in; the soil build checked against the SOLUS grids; outlet file checked (rows, hours, no NaN);
    * a fingerprint of the soil, bedrock and water-table files taken right after each anchor's build and re-checked before EVERY run
      (if someone runs the notebook in the middle, the sweep stops instead of silently using other soil);
    * the list's md5 checked on every start; the manifest remembers every run;
    * a tRIBS failure on one routing point is recorded and the sweep goes on; three failures in a row stop it;
    * on a restart the anchor's soil is rebuilt (about 1 min) so the files are known to be right; finished runs are never rerun.
    NOTHING IS DELETED except a run's own half-finished files when that run is restarted. Finished runs keep all their files.

DISK.  A W1 run makes about 17 MB of results files. 91 + 31 runs is about 2.1 GB. The guard stops the sweep below --min-free-gb (2).

HOW TO RUN (from workspaces/Sycamore/lab/; keep run_w4_draw12_sycamore.py, run_trial_inedit_sycamore.py and this file together):
    python -u run_sweep_routing_sycamore.py --make-list                 # once: make and freeze the list, print its coverage
    python -u run_sweep_routing_sycamore.py --dry-run                   # checks only, builds nothing
    python -u run_sweep_routing_sycamore.py --limit 2                   # first check: 1 build + baseline run 0 + routing point 1 (about 5 min)
    python -u run_sweep_routing_sycamore.py                             # the whole sweep (about 3 h), then look at the summary
    python -u run_sweep_routing_sycamore.py --summary-only
If the connection drops, run the same command again: finished runs are skipped.
Do NOT run the notebook (or anything that rebuilds data/model/) while this runs.

Options: --batch NAME (default rt1; names every file and run, so a follow-up batch never collides)  --anchors 12:60,25:30  --seed N  --n N
         --cv-range LO HI  --flowexp-range LO HI  --n-range LO HI (only for --make-list)  --limit N  --redo ID  --timeout-min N

Files it makes (in lab/): {batch}_list.csv, {batch}_manifest.json, {batch}_logs/, score_out_{batch}/scores_{batch}.csv,
{batch}_summary.csv, w2.lock (only while running), and per run: <run>.in, simulation_<run>.log, results/<run>*.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np

import run_w4_draw12_sycamore as base      # the tested W4 runner: build, check, tRIBS and scorer functions
import run_trial_inedit_sycamore as trial  # the tested shortcut: edit and compare functions, and the W1 window

# --------------------------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------------------------
BATCH_DEFAULT = "rt1"
SEED_DEFAULT = 501
N_DEFAULT = 60
ANCHORS_DEFAULT = "12:60,25:30"               # draw number : how many rows of the list it runs (always rows 1..count)
KEYS = ("KINEMVELCOEF", "FLOWEXP", "CHANNELROUGHNESS")      # order matters: the sampler uses the random generator in this order
SHORT = {"KINEMVELCOEF": "cv", "FLOWEXP": "flowexp", "CHANNELROUGHNESS": "n"}
BOX_DEFAULT = {"KINEMVELCOEF": (4.0, 10.0), "FLOWEXP": (0.26, 0.32), "CHANNELROUGHNESS": (0.04, 0.08)}   # Stage 2 box
DECIMALS = {"KINEMVELCOEF": 3, "FLOWEXP": 4, "CHANNELROUGHNESS": 4}
SANITY = {"KINEMVELCOEF": (1.0, 20.0), "FLOWEXP": (0.15, 0.5), "CHANNELROUGHNESS": (0.01, 0.2)}   # a typo guard, not the box
BASE_ROUTING = {"KINEMVELCOEF": 6, "FLOWEXP": 0.3, "CHANNELROUGHNESS": 0.04}
EXPECTED_DRAWS = {12: (0.187, 0.5874, 0.2939), 25: (0.1852, 0.5055, 0.3358)}   # Ks, f, PsiB from the handoff / Stage 5 docs
REF_ID = trial.REF_ID                          # read-only reference .in (W2 baseline, same Jan 18 start)
OUTLET = "_Outlet.qout"
MB_PER_RUN = 17.0
MAX_CONSECUTIVE_FAILS = 3
CLOSE_ARRIVAL_H, CLOSE_PEAK_PCT = 1.0, 15.0    # Stage 5 closing criteria (looking only; the W1B selection rule is NOT applied here)
STAGE3_W1_REF = {12: 1.106, 25: 1.087}         # flood-week volume ratio at baseline routing, for the header of the summary


def log(msg: str = "") -> None:
    base.log(msg)


class RunError(base.RunError):
    pass


# --------------------------------------------------------------------------------------------
# The routing list: sampler (copied from run_lhs_Ks_f_120_tribs6.py, linear branch), rounding, freezing
# --------------------------------------------------------------------------------------------
def stratified_samples(n: int, box: dict, seed: int) -> dict:
    """Linear stratified-independent sampling: n equal slices per parameter, one uniform draw per slice, shuffled independently.
    Same seed + same n + same bounds + same key order = identical draws (same algorithm as the project's sampler)."""
    rng = np.random.default_rng(seed)
    out = {}
    for key in KEYS:
        lo, hi = box[key]
        edges = np.linspace(lo, hi, n + 1)
        pts = rng.uniform(edges[:-1], edges[1:])
        rng.shuffle(pts)
        out[key] = pts
    return out


def list_path(lab: Path, batch: str) -> Path:
    return lab / f"{batch}_list.csv"


def read_list(path: Path) -> list:
    if not path.is_file():
        raise base.Abort(f"{path.name} not found. Make it first: python -u run_sweep_routing_sycamore.py --make-list")
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append(dict(k=int(r["k"]), KINEMVELCOEF=float(r["cv"]), FLOWEXP=float(r["flowexp"]), CHANNELROUGHNESS=float(r["n"])))
    problems = []
    if [r["k"] for r in rows] != list(range(1, len(rows) + 1)):
        problems.append("the k column is not 1, 2, 3, ... in order")
    seen = set()
    for r in rows:
        trip = tuple(r[k] for k in KEYS)
        if trip in seen:
            problems.append(f"row {r['k']}: the same cv, flowexp, n appear twice")
        seen.add(trip)
        if trip == tuple(float(BASE_ROUTING[k]) for k in KEYS):
            problems.append(f"row {r['k']}: equals the baseline routing (run 0 already covers it)")
        for key in KEYS:
            v, (lo, hi) = r[key], SANITY[key]
            if not (lo <= v <= hi):
                problems.append(f"row {r['k']}: {SHORT[key]} = {v:g} is outside the sanity limits {lo:g}-{hi:g}")
            if round(v, DECIMALS[key]) != v:
                problems.append(f"row {r['k']}: {SHORT[key]} = {v!r} has more than {DECIMALS[key]} decimals")
    if problems:
        raise base.Abort(f"{path.name} has problems:\n  " + "\n  ".join(problems[:12]))
    return rows


def file_md5(path: Path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def coverage_report(cols: dict, n_sub: int) -> list:
    """How evenly the points cover each parameter's range (6 equal bins), for all rows and for the first n_sub rows; plus the largest
    pairwise correlation (the sampler shuffles each parameter on its own, so a chance correlation between parameters is possible)."""
    out = []
    for label, sl in (("all rows", slice(None)), (f"first {n_sub} rows", slice(0, n_sub))):
        out.append(f"  coverage, {label}:")
        for key in KEYS:
            v = np.asarray(cols[key][sl], float)
            lo, hi = BOX_USED[key]
            counts = np.histogram(v, bins=6, range=(lo, hi))[0]
            out.append(f"    {SHORT[key]:8s} range {v.min():.4g}-{v.max():.4g}; points per sixth of the box: {list(map(int, counts))}"
                       + ("   <-- an empty sixth" if (counts == 0).any() else ""))
        r = np.corrcoef([np.asarray(cols[k][sl], float) for k in KEYS])
        pairs = [("cv", "flowexp", r[0, 1]), ("cv", "n", r[0, 2]), ("flowexp", "n", r[1, 2])]
        worst = max(pairs, key=lambda p: abs(p[2]))
        out.append("    largest correlation between two parameters: " + f"{worst[0]}/{worst[1]} {worst[2]:+.2f}"
                   + ("   (high for this many points; consider another seed BEFORE freezing)" if abs(worst[2]) > 0.3 else ""))
    return out


BOX_USED = dict(BOX_DEFAULT)


def make_list(args, lab: Path) -> int:
    path = list_path(lab, args.batch)
    if path.exists():
        print(f"{path.name} already exists and is NEVER overwritten (it is the frozen list). If you really want a new list, give it a new --batch name.")
        return 1
    box = {"KINEMVELCOEF": tuple(args.cv_range), "FLOWEXP": tuple(args.flowexp_range), "CHANNELROUGHNESS": tuple(args.n_range)}
    for key, (lo, hi) in box.items():
        if not lo < hi:
            print(f"The {SHORT[key]} range must have lo < hi; got {lo:g}, {hi:g}.")
            return 1
    BOX_USED.update(box)
    raw = stratified_samples(args.n, box, args.seed)
    cols = {k: np.round(raw[k], DECIMALS[k]) for k in KEYS}
    trips = list(zip(*[cols[k] for k in KEYS]))
    if len(set(trips)) != len(trips):
        print("Two rows came out identical after rounding. Use another seed.")
        return 1
    if any(tuple(float(BASE_ROUTING[k]) for k in KEYS) == t for t in trips):
        print("A row equals the baseline routing exactly. Use another seed.")
        return 1
    n_sub = min(args.second_n, args.n)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["k", "cv", "flowexp", "n", "seed", "box_cv", "box_flowexp", "box_n"])
        for i in range(args.n):
            w.writerow([i + 1] + [f"{cols[k][i]:.{DECIMALS[k]}f}".rstrip("0").rstrip(".") for k in KEYS] + [args.seed]
                       + [f"{box[k][0]:g}-{box[k][1]:g}" for k in KEYS])
    read_list(path)         # the saved file must pass the same checks the sweep applies
    print(f"Wrote and froze {path.name}: {args.n} points, seed {args.seed}, linear stratified-independent sampling.")
    print(f"  box: cv {box['KINEMVELCOEF'][0]:g}-{box['KINEMVELCOEF'][1]:g}, flowexp {box['FLOWEXP'][0]:g}-{box['FLOWEXP'][1]:g}, "
          f"n {box['CHANNELROUGHNESS'][0]:g}-{box['CHANNELROUGHNESS'][1]:g}.  md5 {file_md5(path)}")
    print("\n".join(coverage_report(cols, n_sub)))
    print("  first 8 rows:")
    for i in range(min(8, args.n)):
        print(f"    k={i + 1:>3}  cv {cols['KINEMVELCOEF'][i]:g}  flowexp {cols['FLOWEXP'][i]:g}  n {cols['CHANNELROUGHNESS'][i]:g}")
    print("Next: python -u run_sweep_routing_sycamore.py --dry-run")
    return 0


# --------------------------------------------------------------------------------------------
# Names, jobs, small helpers
# --------------------------------------------------------------------------------------------
def num(x):
    return trial.num(x)


def soil_tag(draw: dict) -> str:
    return f"ks{base.vtag(draw['ks'])}_psib{base.vtag(draw['psib'])}_f{base.vtag(draw['f'])}_d{draw['idx']:03d}"


def run_name(batch: str, draw: dict, k: int, vals: dict) -> str:
    return (f"W1_{batch}_{soil_tag(draw)}_r{k:03d}_cv{base.vtag(vals['KINEMVELCOEF'])}_fe{base.vtag(vals['FLOWEXP'])}"
            f"_n{base.vtag(vals['CHANNELROUGHNESS'])}")


def src_name(batch: str, draw: dict) -> str:
    return f"W1SRC_{batch}_{soil_tag(draw)}"


def stage3_w1_id(draw: dict) -> str:
    """The Stage 3 W1 run of this soil: the design run id (W1_...), which the W4 name was made from."""
    return "W1_" + draw["run_id"][len("W4_"):]


def make_job(draw: dict, rid: str) -> dict:
    return trial.make_job(draw, rid)


def parse_anchors(text: str) -> list:
    out = []
    for part in text.replace(" ", "").split(","):
        if not part:
            continue
        d, _, c = part.partition(":")
        try:
            out.append((int(d), int(c) if c else None))
        except ValueError:
            raise base.Abort(f"--anchors must look like 12:60,25:30 (draw number, then how many list rows); got {text!r}.")
    if not out or len({d for d, _ in out}) != len(out):
        raise base.Abort("--anchors is empty or names a draw twice.")
    return out


def set_routing(vals: dict) -> None:
    base.BASE.update({k: num(vals[k]) for k in KEYS})


def soil_fingerprint(lab: Path, src_in: Path) -> dict:
    """Size and modified time of every file the shortcut runs depend on (soil grids, soil table, bedrock, water table); md5 of the two
    small tables. Taken right after a build; compared before every run."""
    soil = lab / base.MODEL_SOIL
    files = {soil / "scgrid.gdf", soil / "soils.sdt"}
    for v in base.read_gdf(soil / "scgrid.gdf").values():
        files.add(lab / (v + ".asc"))
    ini = base.read_in_file(src_in)
    for kw in ("GWATERFILE", "BEDROCKFILE"):
        if kw in ini:
            files.add(lab / ini[kw])
    fp = {}
    for p in sorted(files):
        if not p.is_file():
            raise RunError(f"{p} is missing after the build.")
        st = p.stat()
        fp[str(p.relative_to(lab))] = [st.st_size, st.st_mtime_ns]
    for name in ("scgrid.gdf", "soils.sdt"):
        fp["md5:" + name] = file_md5(soil / name)
    return fp


def fingerprint_problem(lab: Path, src_in: Path, saved: dict) -> str | None:
    try:
        now = soil_fingerprint(lab, src_in)
    except RunError as exc:
        return str(exc)
    changed = [k for k in sorted(set(saved) | set(now)) if saved.get(k) != now.get(k)]
    return ("these soil/bedrock/water-table files changed since the anchor's build: " + ", ".join(changed[:6])) if changed else None


def write_edited_in(src_in: Path, new_in: Path, vals: dict, rid: str) -> None:
    new_vals = {k: f"{num(vals[k]):g}" for k in KEYS}
    new_vals["OUTFILENAME"] = f"results/{rid}"
    try:
        new_in.write_text(trial.edit_in_text(src_in.read_text(), new_vals))
        trial.check_edit(src_in, new_in, new_vals)
    except base.RunError as exc:
        raise base.Abort(f"the .in edit is wrong (this is a script problem, not a model problem): {exc}")
    allowed = {trial.value_line_number(src_in, k) for k in new_vals}
    bad = [d for d in trial.diff_lines(src_in, new_in) if d[0] not in allowed]
    if bad:
        raise base.Abort(f"the edited .in differs from the source .in on lines other than the four that may change: {bad[:3]}")


def md5_or_none(p: Path):
    try:
        return file_md5(p)
    except OSError:
        return None


# --------------------------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------------------------
def preflight(args, lab: Path, plan: list, rows: list, man: dict, ref_in: Path) -> list:
    out = []

    def add(level, msg):
        out.append((level, msg))

    for f in (args.notebook, "run_w4_draw12_sycamore.py", "run_trial_inedit_sycamore.py", base.SCORER):
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
        miss = [k for k in list(base.ROUTING_KWS) + ["OPTPERCOLATION", "OPTSNOW", "RUNTIME", "STARTDATE", "GWATERFILE", "BEDROCKFILE", "OUTFILENAME"] if k not in got]
        bad = [f"{k} is {got.get(k)!r}, wanted {v!r}" for k, v in BASE_ROUTING.items() if k in got and not base.same_value(got[k], v)]
        add("FAIL" if (miss or bad) else "OK", f"reference .in ({ref_in.name}) is readable, has the needed keywords and baseline routing"
            + (f" (missing: {miss})" if miss else "") + (f" ({'; '.join(bad)})" if bad else ""))
        add("OK" if base.same_value(got.get("STARTDATE", ""), trial.W1["start"]) else "FAIL", f"reference .in starts {got.get('STARTDATE')} (W1 starts {trial.W1['start']})")
    except Exception as exc:
        add("FAIL", f"reference .in {ref_in.name}: {exc}")
    # the list: valid, and the same list as at the first start
    lp = list_path(lab, args.batch)
    md5 = file_md5(lp)
    if man.get("list_md5") and man["list_md5"] != md5:
        add("FAIL", f"{lp.name} has changed since the first start (md5 {md5} vs recorded {man['list_md5']}). The list is frozen: restore it, or use a new --batch name.")
    else:
        add("OK", f"routing list {lp.name}: {len(rows)} valid rows, md5 {md5}" + ("" if man.get("list_md5") else " (will be recorded at the first start)"))
    for a in plan:
        if a["count"] > len(rows):
            add("FAIL", f"draw {a['draw']['idx']}: --anchors asks for {a['count']} rows but the list has {len(rows)}")
    try:
        import nbformat
        nb = nbformat.read(str(lab / args.notebook), as_version=4)
        idx = base.find_settings_cell(nb)
        for a in plan:
            set_routing(BASE_ROUTING)
            base.patched_source(nb.cells[idx].source, make_job(a["draw"], a["src_id"]))
        add("OK", f"notebook Settings cell (cell {idx}) patches correctly for each anchor's soil build")
    except ImportError:
        add("FAIL", "nbformat is not installed: pip install nbformat nbclient ipykernel")
    except Exception as exc:
        add("FAIL", f"Settings cell check: {exc}")
    if not args.skip_kernel_check:
        try:
            from nbformat.v4 import new_code_cell, new_notebook
            nb = new_notebook(cells=[new_code_cell("import sys, os\nfrom pytRIBS.classes import *\nprint('KERNEL', sys.executable, os.getcwd())")])
            base.execute_nb(nb, lab, 300)
            add("OK", "notebook kernel imports pytRIBS")
        except Exception as exc:
            add("FAIL", f"notebook kernel could not import pytRIBS: {type(exc).__name__}: {str(exc)[:300]}")
    # names: a planned run may already have files only if the manifest knows it (half-finished runs are restarted)
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
        w1 = lab / "results" / f"{stage3_w1_id(d)}{OUTLET}"
        add("OK" if w1.is_file() else "WARN", f"draw {d['idx']}: Stage 3 W1 outlet file for the run-0 comparison " + ("found" if w1.is_file() else f"NOT found ({w1.name}); run 0 will not be compared"))
    return out


# --------------------------------------------------------------------------------------------
# One anchor
# --------------------------------------------------------------------------------------------
def build_anchor(a: dict, args, lab: Path, man: dict, ref_in: Path) -> Path:
    """Build the anchor's soil at baseline routing (no tRIBS), check it, take the fingerprint. Returns the source .in path."""
    d, sid = a["draw"], a["src_id"]
    log_dir = lab / base.LOG_DIR
    set_routing(BASE_ROUTING)
    job = make_job(d, sid)
    base.remove_stale(lab, sid)
    log(f"=== anchor draw {d['idx']}: building its soil once (Ks x{d['ks']:g}, f x{d['f']:g}, PsiB x{d['psib']:g}) as '{sid}'; no tRIBS run ===")
    t = time.time()
    text, warns = base.build_run(job, args, lab, log_dir)
    src_in = lab / f"{sid}.in"
    base.check_in_file(src_in, ref_in, job)
    grid_note = base.check_grids(lab, job, t)
    fp = soil_fingerprint(lab, src_in)
    man["anchors"][str(d["idx"])] = dict(src_id=sid, built=base.now_iso(), build_min=round((time.time() - t) / 60, 2), fingerprint=fp, grid_check=grid_note)
    base.save_manifest(man)
    log(f"build OK in {base.fmt_min(time.time() - t)}; .in matches the reference except the allowed lines. {grid_note}")
    for w in warns:
        log(f"notebook said: {w}")
    return src_in


def one_line(row: dict) -> str:
    g = lambda c: base.fnum(row, c)
    return (f"vol {base.fm(g('flood_week__vol_ratio'), 3)}  peak {base.fm(g('flood_week__peak_err_pct'), 1, True)} %  "
            f"arrival {base.fm(g('flood_week__arrival_err_h'), 2, True)} h  KGE r {base.fm(g('flood_week__kge_r'), 3)}  "
            f"gamma {base.fm(g('flood_week__kge_gamma'), 3)}  KGE {base.fm(g('flood_week__kge2012'), 3)}")


def run_one(a: dict, p: dict, src_in: Path, args, lab: Path, man: dict, ref_in: Path) -> None:
    """One routing point: edit the .in, run tRIBS, check the outlet, score. RunError = this run failed (recorded, sweep goes on)."""
    d, rid, vals, k = a["draw"], p["rid"], p["vals"], p["k"]
    score_dir, log_dir = lab / base.SCORE_DIR, lab / base.LOG_DIR
    prob = fingerprint_problem(lab, src_in, man["anchors"][str(d["idx"])]["fingerprint"])
    if prob:
        raise base.Abort(f"The soil build was disturbed before {rid}: {prob.rstrip('.')}. Was the notebook run? Nothing was started; run the same command again "
                         f"(it rebuilds the anchor's soil first).")
    free = base.disk_free_gb(lab)
    if free < args.min_free_gb:
        raise base.Abort(f"Only {free:.1f} GB free; need {args.min_free_gb:g} GB. Free some space and run again to resume.")
    comp = base.competing_tribs()
    if comp:
        raise base.Abort("Another tRIBS is running, so the sweep stops (two at once is not allowed): " + comp[0])
    rec = man["runs"].setdefault(rid, {})
    rec.clear()
    rec.update(status="running", started=base.now_iso(), draw=d["idx"], k=k, cv=num(vals["KINEMVELCOEF"]), flowexp=num(vals["FLOWEXP"]),
               n=num(vals["CHANNELROUGHNESS"]), kind=p["kind"])
    base.save_manifest(man)
    base.remove_stale(lab, rid)
    new_in = lab / f"{rid}.in"
    t0 = time.time()
    write_edited_in(src_in, new_in, vals, rid)
    set_routing(vals)
    job = make_job(d, rid)
    try:
        base.check_in_file(new_in, ref_in, job)
    except base.RunError as exc:
        raise base.Abort(f"{rid}: the edited .in does not match the reference: {exc}")
    base.install_term_handler()
    rc, secs = base.run_tribs(args.tribs, new_in.name, f"simulation_{rid}.log", base.timeout_for(args),
                              progress_file=lab / "results" / f"{rid}{OUTLET}", hours=trial.W1["hours"],
                              heartbeat_min=args.heartbeat_min, label=rid)
    rec.update(tribs_min=round(secs / 60, 2), tribs_returncode=rc)
    if rc != 0:
        raise RunError(f"tRIBS exited with code {rc}; see simulation_{rid}.log")
    outlet = base.check_outlet(lab / "results" / f"{rid}{OUTLET}", trial.W1["hours"])
    rec.update(outlet=outlet, status="unscored", output_mb=base.output_mb(lab, rid))
    base.save_manifest(man)
    if p["kind"] == "baseline":
        ref = lab / "results" / f"{stage3_w1_id(d)}{OUTLET}"
        if ref.is_file():
            same = md5_or_none(ref) == outlet["md5"]
            rec["stage3_w1_identical"] = same
            if same:
                log(f"{rid}: run 0 is BYTE-IDENTICAL to the Stage 3 W1 result ({ref.name}): the shortcut path reproduces the earlier work end to end.")
            else:
                gap = trial.outlet_gap(lab, rid, stage3_w1_id(d))
                rec["stage3_w1_gap"] = gap
                log(f"{rid}: WARNING: run 0 is NOT identical to the Stage 3 W1 result ({ref.name}): {gap}. The sweep goes on, but find out why before "
                    f"trusting comparisons with Stage 3 numbers.")
    try:
        row = base.run_scorer(rid, str(lab / args.scorer), score_dir, log_dir)
        rec["scores"] = base.show_scores(row)
        rec["status"] = "done"
        log(f"{rid}: {base.fmt_min(time.time() - t0)}; {one_line(row)}")
    except base.RunError as exc:
        rec["score_error"] = str(exc)
        log(f"{rid}: scoring failed (the run itself is fine): {exc}")
    rec["finished"] = base.now_iso()
    base.save_manifest(man)


# --------------------------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------------------------
def load_scores(lab: Path) -> dict:
    p = lab / base.SCORE_DIR / base.SCORE_CSV
    if not p.is_file():
        return {}
    with open(p, newline="") as fh:
        return {r["run_id"]: r for r in csv.DictReader(fh)}


def print_summary(man: dict, lab: Path, plan_all: list, args) -> None:
    scores = load_scores(lab)
    out_rows, cols = [], []
    log()
    log("=" * 100)
    log(f"ROUTING SWEEP SUMMARY, batch {args.batch} (W1, Kchan off). Closing criteria for LOOKING only: arrival within {CLOSE_ARRIVAL_H:g} h and peak within {CLOSE_PEAK_PCT:g}%.")
    log("The W1B selection rule is NOT applied here; it is fixed separately, before any W1B result is seen.")
    for a in plan_all:
        d = a["draw"]
        recs = []
        for p in a["all"]:
            r = man["runs"].get(p["rid"])
            if not r:
                continue
            row = scores.get(p["rid"], {})
            recs.append((p, r, row))
            flat = dict(anchor=d["idx"], k=p["k"], kind=p["kind"], cv=num(p["vals"]["KINEMVELCOEF"]), flowexp=num(p["vals"]["FLOWEXP"]),
                        n=num(p["vals"]["CHANNELROUGHNESS"]), run_id=p["rid"], status=r.get("status", ""), error=r.get("error", ""))
            flat.update({c: v for c, v in row.items() if c not in ("run_id", "window")})
            out_rows.append(flat)
            for c in flat:
                if c not in cols:
                    cols.append(c)
        done = [x for x in recs if x[1].get("status") == "done"]
        failed = [x for x in recs if x[1].get("status") == "failed"]
        log("")
        log(f"--- anchor draw {d['idx']} (Ks x{d['ks']:g}, f x{d['f']:g}, PsiB x{d['psib']:g}): {len(done)} scored of {len(a['all'])} planned; {len(failed)} failed ---")
        if not done:
            continue
        g = lambda row, c: base.fnum(row, c)
        sweep = [x for x in done if x[0]["kind"] == "sweep"]
        b0 = next((x for x in done if x[0]["kind"] == "baseline"), None)
        if b0:
            log(f"  run 0 (baseline routing 6 / 0.3 / 0.04): {one_line(b0[2])}   (Stage 3 W1 volume ratio at baseline routing: {STAGE3_W1_REF.get(d['idx'], '-')})"
                + ("   [byte-identical to Stage 3]" if b0[1].get("stage3_w1_identical") else ("   [NOT identical to Stage 3]" if b0[1].get("stage3_w1_identical") is False else "")))
        if sweep:
            arr = [g(x[2], "flood_week__arrival_err_h") for x in sweep]
            pk = [g(x[2], "flood_week__peak_err_pct") for x in sweep]
            vol = [g(x[2], "flood_week__vol_ratio") for x in sweep]
            fin = lambda v: [t for t in v if math.isfinite(t)]
            def rng(v, nd, sign=False):
                v = fin(v)
                return f"{base.fm(min(v), nd, sign)} to {base.fm(max(v), nd, sign)}" if v else "-"
            meets = [x for x in sweep if math.isfinite(g(x[2], "flood_week__arrival_err_h")) and math.isfinite(g(x[2], "flood_week__peak_err_pct"))
                     and abs(g(x[2], "flood_week__arrival_err_h")) <= CLOSE_ARRIVAL_H and abs(g(x[2], "flood_week__peak_err_pct")) <= CLOSE_PEAK_PCT]
            log(f"  sweep points scored: {len(sweep)};  arrival {rng(arr, 2, True)} h;  peak {rng(pk, 1, True)} %;  flood-week volume ratio {rng(vol, 3)}")
            log(f"  points meeting BOTH criteria (|arrival| <= {CLOSE_ARRIVAL_H:g} h and |peak| <= {CLOSE_PEAK_PCT:g}%): {len(meets)} of {len(sweep)}")
            noarr = [x for x in sweep if not math.isfinite(g(x[2], "flood_week__arrival_err_h"))]
            if noarr:
                log(f"  points with no arrival time (flag in the CSV): {len(noarr)}")
            ranked = sorted([x for x in sweep if math.isfinite(g(x[2], "flood_week__arrival_err_h"))], key=lambda x: abs(g(x[2], "flood_week__arrival_err_h")))
            log(f"  the {min(10, len(ranked))} points closest in arrival (for looking only):")
            log(f"    {'k':>3} {'cv':>7} {'flowexp':>8} {'n':>7} {'vol':>7} {'peak %':>8} {'arrival h':>10} {'KGE r':>7} {'gamma':>7} {'KGE':>7}")
            for p, r, row in ranked[:10]:
                v = p["vals"]
                log(f"    {p['k']:>3} {num(v['KINEMVELCOEF']):>7g} {num(v['FLOWEXP']):>8g} {num(v['CHANNELROUGHNESS']):>7g} "
                    f"{base.fm(g(row, 'flood_week__vol_ratio'), 3):>7} {base.fm(g(row, 'flood_week__peak_err_pct'), 1, True):>8} "
                    f"{base.fm(g(row, 'flood_week__arrival_err_h'), 2, True):>10} {base.fm(g(row, 'flood_week__kge_r'), 3):>7} "
                    f"{base.fm(g(row, 'flood_week__kge_gamma'), 3):>7} {base.fm(g(row, 'flood_week__kge2012'), 3):>7}")
        if failed:
            for p, r, row in failed[:5]:
                log(f"  FAILED k={p['k']} {p['rid']}: {r.get('error', '')[:160]}")
    if out_rows:
        path = lab / f"{args.batch}_summary.csv"
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(out_rows)
        log("")
        log(f"WROTE {path.name} ({len(out_rows)} rows: routing values, status and every score column). Rising-limb steepness is not in the scorer; it is computed later")
        log("from the kept outlet files (the analysis step), so no run needs repeating for it.")


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Routing sweep with the .in-edit shortcut (see the top of this file).")
    ap.add_argument("--batch", default=BATCH_DEFAULT, help=f"name of this batch (default {BATCH_DEFAULT}); names the list, manifest, scores and runs")
    ap.add_argument("--make-list", action="store_true", help="make and freeze the routing list, print its coverage, and stop")
    ap.add_argument("--seed", type=int, default=SEED_DEFAULT, help=f"(--make-list) random seed (default {SEED_DEFAULT}; used before: 42, 201-205, 301)")
    ap.add_argument("--n", type=int, default=N_DEFAULT, help=f"(--make-list) number of points (default {N_DEFAULT})")
    ap.add_argument("--second-n", type=int, default=30, help="(--make-list) rows the second anchor will use; only for the coverage printout")
    ap.add_argument("--cv-range", type=float, nargs=2, default=BOX_DEFAULT["KINEMVELCOEF"], metavar=("LO", "HI"))
    ap.add_argument("--flowexp-range", type=float, nargs=2, default=BOX_DEFAULT["FLOWEXP"], metavar=("LO", "HI"))
    ap.add_argument("--n-range", type=float, nargs=2, default=BOX_DEFAULT["CHANNELROUGHNESS"], metavar=("LO", "HI"))
    ap.add_argument("--anchors", default=ANCHORS_DEFAULT, help=f"draw:rows pairs in run order (default {ANCHORS_DEFAULT}); draws come from lhs_volume_list.csv")
    ap.add_argument("--no-baseline-run", action="store_true", help="skip run 0 (baseline routing) at each anchor")
    ap.add_argument("--dry-run", action="store_true", help="check everything, show the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the tables from the runs already scored, then stop")
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

    # ---- names for this batch (the imported module reads these globals) ----
    b = args.batch
    base.LOG_DIR, base.SCORE_DIR, base.SCORE_CSV, base.MANIFEST = f"{b}_logs", f"score_out_{b}", f"scores_{b}.csv", f"{b}_manifest.json"
    (lab / base.LOG_DIR).mkdir(exist_ok=True)
    base._LOGFILE = str(lab / base.LOG_DIR / "sweep_log.txt")
    args.redo_ids = {t.strip() for t in (args.redo or "").split(",") if t.strip()}
    args.notebook, nb_note = base.pick_notebook(lab, args.notebook)

    try:
        rows = read_list(list_path(lab, b))
        spec = parse_anchors(args.anchors)
        # The draws are read FIRST, while the imported module still has its W4 settings (its checks name the design run with the W4 prefix).
        base.EXPECTED_DRAWS.update(EXPECTED_DRAWS)
        draws = base.read_anchors(lab, ",".join(str(d) for d, _ in spec))
    except base.Abort as exc:
        print(str(exc))
        return 1
    base.WIN = dict(trial.W1)                       # now switch to W1: window, start, length, cold-start water table
    base.BASE["GW_FRACTION"] = trial.W1["gw"]
    set_routing(BASE_ROUTING)
    man = base.load_manifest()
    man.setdefault("runs", {})
    man.setdefault("anchors", {})

    plan_all = []
    for (didx, count), d in zip(spec, draws):
        count = len(rows) if count is None else count
        a = dict(draw=d, count=count, src_id=src_name(b, d), all=[], pending=[])
        pts = []
        if not args.no_baseline_run:
            pts.append(dict(k=0, kind="baseline", vals={k: BASE_ROUTING[k] for k in KEYS}))
        for r in rows[:count]:
            pts.append(dict(k=r["k"], kind="sweep", vals={k: r[k] for k in KEYS}))
        for p in pts:
            p["rid"] = run_name(b, d, p["k"], p["vals"])
        a["all"] = pts
        for p in pts:
            st = (man["runs"].get(p["rid"]) or {}).get("status")
            if st != "done" or p["rid"] in args.redo_ids:
                a["pending"].append(p)
        plan_all.append(a)
    ref_in = lab / f"{REF_ID}.in"

    if args.summary_only:
        print_summary(man, lab, plan_all, args)
        return 0

    log(f"Working folder: {lab}")
    log(f"Notebook: {args.notebook}" + (f"  (note: {nb_note})" if nb_note else ""))
    log(f"Batch {b}: routing list {list_path(lab, b).name} ({len(rows)} rows); window W1 {trial.W1['start']}, {trial.W1['hours']} h, water table {trial.W1['gw']} x bedrock, Kchan off.")
    for a in plan_all:
        d = a["draw"]
        log(f"  anchor draw {d['idx']}: Ks x{d['ks']:g}, f x{d['f']:g}, PsiB x{d['psib']:g}: {len(a['all'])} planned runs "
            f"({'run 0 + ' if not args.no_baseline_run else ''}rows 1-{a['count']}), {len(a['pending'])} still to do")
    n_todo = sum(len(a["pending"]) for a in plan_all)
    n_builds = sum(1 for a in plan_all if a["pending"])
    est = (n_todo * 1.8 + n_builds * 1.1) if args.limit is None else (min(n_todo, args.limit) * 1.8 + n_builds * 1.1)
    log(f"  Time: about {base.fmt_hm(est * 60)} ({min(n_todo, args.limit) if args.limit else n_todo} runs x about 1.8 min + {n_builds} soil build(s) x about 1.1 min).")

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
        log("Every planned run is already done. Printing the summary.")
        print_summary(man, lab, plan_all, args)
        return 0

    try:
        base.take_lock()
    except base.Abort as exc:
        log(str(exc))
        return 1

    man.setdefault("list_md5", file_md5(list_path(lab, b)))
    man["batch"], man["script"], man["window"] = b, "run_sweep_routing_sycamore.py", trial.W1
    man.setdefault("started", base.now_iso())
    base.save_manifest(man)
    code, started, fails_in_row, n_started = 0, 0, 0, 0
    t_sweep = time.time()
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
                        row = base.run_scorer(p["rid"], str(lab / args.scorer), lab / base.SCORE_DIR, lab / base.LOG_DIR)
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
                except base.RunError as exc:                     # this run failed (tRIBS error or timeout, bad outlet file); the others may be fine
                    fails_in_row += 1
                    man["runs"].setdefault(p["rid"], {}).update(status="failed", error=str(exc), finished=base.now_iso())
                    base.save_manifest(man)
                    log(f"{p['rid']}: FAILED ({exc}). Recorded; the sweep goes on.")
                    if fails_in_row >= MAX_CONSECUTIVE_FAILS:
                        raise base.Abort(f"{MAX_CONSECUTIVE_FAILS} runs failed in a row, so the sweep stops. Look at the last simulation_*.log files before restarting.")
                done_n = sum(1 for x in man["runs"].values() if x.get("status") == "done")
                left = sum(len([q for q in a2["pending"] if (man["runs"].get(q["rid"]) or {}).get("status") != "done"]) for a2 in plan_all)
                if n_started:
                    per = (time.time() - t_sweep) / n_started
                    log(f"    progress: {done_n} runs done in total; about {left} left; roughly {base.fmt_hm(per * left)} to go at this pace")
        log("Sweep call finished.")
    except base.Abort as exc:
        log(f"SWEEP STOPPED: {exc}")
        code = 1
    except KeyboardInterrupt:
        for rid, r in man["runs"].items():
            if r.get("status") == "running":
                r["status"] = "interrupted"
        base.save_manifest(man)
        log("SWEEP INTERRUPTED. The lock is removed. Run the same command again to continue; the interrupted run is restarted from its own files.")
        code = 130
    finally:
        set_routing(BASE_ROUTING)
    print_summary(man, lab, plan_all, args)
    return code


if __name__ == "__main__":
    sys.exit(main())
