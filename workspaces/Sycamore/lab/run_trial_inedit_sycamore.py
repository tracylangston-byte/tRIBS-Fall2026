#!/usr/bin/env python3
"""
run_trial_inedit_sycamore.py  --  Stage 5, Step 2: ONE-RUN TRIAL of the ".in-edit shortcut" for routing sweeps.

The question. A routing sweep changes only three numbers (cv = KINEMVELCOEF, flowexp = FLOWEXP, n = CHANNELROUGHNESS). Those three
numbers live in the model's input file (the .in file). The soil settings (Ks, f, PsiB) do NOT: they are baked into scaled soil grid
files that the notebook builds. So for a routing sweep at ONE soil setting, the soil grids never change, and the notebook build
(about 1 min) is repeated for nothing. The shortcut: build once, then for every routing run copy the .in file, change the three
numbers and the output name, and start tRIBS. This script tests, on ONE run, whether that gives exactly the same answer as a normal
full build.

What it does (all on W1, the flood week Jan 18-25 2017, 168 h, water table 0.82 x bedrock, Kchan off, draw 12 soil setting):
    1. BUILD "src":   normal notebook build of draw 12 at BASELINE routing (cv 6, flowexp 0.3, n 0.04). No tRIBS run. Its .in file is
                      the source for the shortcut.
    2. SHORTCUT run:  copies the src .in file, changes ONLY cv, flowexp, n and the output name (to the TEST values below), checks the
                      edit, and runs tRIBS on it. No notebook build.
    3. BUILD "full":  normal notebook build of draw 12 with the SAME test routing values patched into the notebook's Settings cell.
                      Runs tRIBS on it.
    4. COMPARES:      (a) the shortcut .in against the full-build .in, line by line (only the output name should differ);
                      (b) every results file, byte for byte (the main one is *_Outlet.qout).
The shortcut is VALID if the two outlet files are byte-identical. If they are not, the script says how far apart they are.

Test routing values (all inside the Stage 2 box, all different from baseline, so a no-change result cannot fool the test):
    cv 8, flowexp 0.29, n 0.055      (change with --cv, --flowexp, --n)

This script IMPORTS run_w4_draw12_sycamore.py (keep it in the same folder) and reuses its tested build, check and tRIBS functions,
changing only the window (to W1) and the water table (to 0.82).

NOTHING EXISTING IS OVERWRITTEN. The three runs have their own names (W1TRIAL_src_..., W1TRIAL_short_..., W1TRIAL_full_...); if any
of them already has files in lab/ or results/, the start-up check stops. The reference .in (W2_ks0p2_cv6_psib0p5_nocl.in) is only read.
Nothing is deleted at the end; the files stay so you can look at them.

Run it from workspaces/Sycamore/lab/, in this order:
    python -u run_trial_inedit_sycamore.py --dry-run     # checks only, builds nothing
    python -u run_trial_inedit_sycamore.py               # the whole trial: 2 builds + 2 tRIBS runs, about 6 min

Do NOT run the notebook (or anything that rebuilds data/model/) while this runs. It shares w2.lock with the other runners.

Files it makes (in lab/): the three .in files, simulation_<run>.log for the two tRIBS runs, results/<run>* for the two tRIBS runs,
inedit_trial_logs/ (build output, batch log) and inedit_trial_result.json (every number printed at the end).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import filecmp
import json
import os
import re
import shutil
import sys
import time
from pathlib import Path

import numpy as np

import run_w4_draw12_sycamore as base      # the tested W4 runner: its functions are reused here

# --------------------------------------------------------------------------------------------
# What the trial does
# --------------------------------------------------------------------------------------------
W1 = dict(name="W1", start="01/18/2017/00/00", hours=168, label="jan_storm_week", gw=0.82)   # from run_proxy_sycamore.py / Stage3_Result_v1.md
REF_ID = "W2_ks0p2_cv6_psib0p5_nocl"        # its .in file (same Jan 18 start) is the READ-ONLY reference for the .in comparison
DRAW = 12
LOG_DIR = "inedit_trial_logs"
RESULT_JSON = "inedit_trial_result.json"
OUTLET_SUFFIX = "_Outlet.qout"
ROUTING_KEYS = ("KINEMVELCOEF", "FLOWEXP", "CHANNELROUGHNESS")
CV_TEST, FLOWEXP_TEST, N_TEST = 8, 0.29, 0.055


def log(msg: str = "") -> None:
    base.log(msg)


class RunError(base.RunError):
    pass


# --------------------------------------------------------------------------------------------
# Names, jobs, settings
# --------------------------------------------------------------------------------------------
def num(x):
    """8.0 -> 8 (int), 0.29 -> 0.29, so the Settings cell shows 8 and not 8.0."""
    x = float(x)
    return int(x) if x.is_integer() else x


def trial_ids(draw: dict, cv, fe, n) -> dict:
    soil = f"ks{base.vtag(draw['ks'])}_psib{base.vtag(draw['psib'])}_f{base.vtag(draw['f'])}_d{draw['idx']:03d}"
    base_rt = f"cv{base.vtag(base.BASE['KINEMVELCOEF'])}_fe{base.vtag(base.BASE['FLOWEXP'])}_n{base.vtag(base.BASE['CHANNELROUGHNESS'])}"
    test_rt = f"cv{base.vtag(cv)}_fe{base.vtag(fe)}_n{base.vtag(n)}"
    return dict(src=f"W1TRIAL_src_{soil}_{base_rt}", short=f"W1TRIAL_short_{soil}_{test_rt}", full=f"W1TRIAL_full_{soil}_{test_rt}")


def make_job(draw: dict, run_id: str) -> dict:
    d = dict(draw)
    d["run_id"] = run_id
    return base.make_job(d)


def set_routing(values: dict) -> None:
    """The functions in run_w4_draw12_sycamore.py read the routing values from base.BASE, so a full build with other routing values
    needs BASE changed first (and changed back afterwards)."""
    base.BASE.update(values)


# --------------------------------------------------------------------------------------------
# Editing and comparing .in files
# --------------------------------------------------------------------------------------------
def edit_in_text(text: str, new_values: dict) -> str:
    """Change the VALUE line (the line after 'KEY: ...') of each keyword in new_values. Everything else stays byte for byte."""
    lines = text.splitlines(keepends=True)
    for key, value in new_values.items():
        hits = [i for i, ln in enumerate(lines) if re.match(rf"^{re.escape(key)}:\s", ln.rstrip("\r\n") + " ")]
        if len(hits) != 1:
            raise RunError(f"the .in file should have exactly one '{key}:' line; found {len(hits)}.")
        i = hits[0]
        if i + 1 >= len(lines) or not lines[i + 1].strip():
            raise RunError(f"the line after '{key}:' is missing or empty.")
        ending = lines[i + 1][len(lines[i + 1].rstrip("\r\n")):]
        lines[i + 1] = str(value) + ending
    return "".join(lines)


def check_edit(src_in: Path, new_in: Path, expected: dict) -> int:
    """The edited .in must equal the source .in in every keyword except the edited ones, and the edited ones must read as intended.
    Also counts the text lines that differ. Returns that count."""
    a, b = base.read_in_file(src_in), base.read_in_file(new_in)
    problems = []
    for k in sorted(set(a) | set(b)):
        if k in expected:
            if not base.same_value(b.get(k), expected[k]):
                problems.append(f"{k}: reads {b.get(k)!r}, wanted {expected[k]!r}")
        elif a.get(k) != b.get(k):
            problems.append(f"{k}: source {a.get(k)!r} but edited file {b.get(k)!r}  (this keyword may not change)")
    if problems:
        raise RunError("The edited .in file is not what was intended:\n  " + "\n  ".join(problems))
    return len(diff_lines(src_in, new_in))


def diff_lines(a_path: Path, b_path: Path) -> list:
    """Lines that differ between two text files, as (line number, text in a, text in b). A different line count is reported as one entry."""
    a, b = Path(a_path).read_text().splitlines(), Path(b_path).read_text().splitlines()
    if len(a) != len(b):
        return [(0, f"{len(a)} lines", f"{len(b)} lines")]
    return [(i + 1, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]


def value_line_number(path: Path, key: str) -> int | None:
    """1-based number of the value line (the line after 'KEY: ...')."""
    for i, ln in enumerate(Path(path).read_text().splitlines()):
        if re.match(rf"^{re.escape(key)}:\s", ln + " "):
            return i + 2
    return None


# --------------------------------------------------------------------------------------------
# Comparing results files
# --------------------------------------------------------------------------------------------
def result_files(lab: Path, rid: str) -> dict:
    return {p.name[len(rid):]: p for p in (lab / "results").glob(rid + "*") if p.is_file()}


def read_q(path: Path):
    rows = []
    with open(path) as fh:
        reader = csv.reader(fh)
        next(reader, None)
        for r in reader:
            if len(r) >= 2:
                try:
                    rows.append((float(r[0]), float(r[1])))
                except ValueError:
                    continue
    return np.array(rows)


def compare_results(lab: Path, rid_a: str, rid_b: str) -> list:
    fa, fb = result_files(lab, rid_a), result_files(lab, rid_b)
    out = []
    for suf in sorted(set(fa) | set(fb)):
        pa, pb = fa.get(suf), fb.get(suf)
        if pa is None or pb is None:
            out.append(dict(suffix=suf, status="only in " + ("full" if pa is None else "shortcut")))
        else:
            out.append(dict(suffix=suf, status="identical" if filecmp.cmp(pa, pb, shallow=False) else "DIFFERENT",
                            size_mb=round(pa.stat().st_size / 1e6, 2)))
    return out


def outlet_gap(lab: Path, rid_a: str, rid_b: str) -> dict:
    a, b = read_q(lab / "results" / (rid_a + OUTLET_SUFFIX)), read_q(lab / "results" / (rid_b + OUTLET_SUFFIX))
    if a.shape != b.shape:
        return dict(note=f"different shapes {a.shape} vs {b.shape}")
    d = np.abs(a[:, 1] - b[:, 1])
    i = int(d.argmax())
    return dict(max_abs_diff_m3s=float(d.max()), at_hour=float(a[i, 0]), n_rows_differing=int((d > 0).sum()),
                peak_a=float(a[:, 1].max()), peak_b=float(b[:, 1].max()))


# --------------------------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------------------------
def preflight(args, lab: Path, jobs: dict, ids: dict, ref_in: Path, test_routing: dict, base_routing: dict) -> list:
    out = []

    def add(level, msg):
        out.append((level, msg))

    add("OK" if (lab / args.notebook).is_file() else "FAIL", f"file in lab/: {args.notebook}")
    add("OK" if (lab / base.__file__.split(os.sep)[-1]).is_file() else "FAIL", "run_w4_draw12_sycamore.py is in lab/ (this script imports it)")
    tribs = args.tribs
    add("OK" if (os.path.isfile(tribs) and os.access(tribs, os.X_OK)) else "FAIL", f"tRIBS program: {tribs}")
    free = base.disk_free_gb(lab)
    add("OK" if free >= args.min_free_gb else "FAIL", f"free disk space: {free:.1f} GB (need {args.min_free_gb:g})")
    comp = base.competing_tribs()
    add("FAIL" if comp else "OK", "no other tRIBS running" + (": " + comp[0] if comp else ""))
    for f in ("../init_data/bedrock_depth_m.asc", "../init_data/scgrid.gdf", "../init_data/soils.sdt",
              "../init_data/met/Master_Met.sdf", "../init_data/LandUse.asc"):
        add("OK" if (lab / f).is_file() else "FAIL", f"input: {f}")
    try:
        got = base.read_in_file(ref_in)
        need = list(base.ROUTING_KWS) + ["OPTPERCOLATION", "OPTSNOW", "CHANNELCONDUCTIVITY", "RUNTIME", "STARTDATE", "GWATERFILE",
                                         "BEDROCKFILE", "OUTFILENAME"]
        miss = [k for k in need if k not in got]
        bad = [f"{k} is {got.get(k)!r}, wanted {v!r}" for k, v in base_routing.items() if k in got and not base.same_value(got[k], v)]
        add("FAIL" if (miss or bad) else "OK", f"reference .in ({ref_in.name}) is readable, has the needed keywords and baseline routing"
            + (f" (missing: {miss})" if miss else "") + (f" ({'; '.join(bad)})" if bad else ""))
        add("OK" if base.same_value(got.get("STARTDATE", ""), W1["start"]) else "WARN",
            f"reference .in starts {got.get('STARTDATE')} (W1 starts {W1['start']}; only RUNTIME and the output name are expected to differ)")
    except Exception as exc:
        add("FAIL", f"reference .in {ref_in.name}: {exc}")
    try:
        import nbformat
        nb = nbformat.read(str(lab / args.notebook), as_version=4)
        idx = base.find_settings_cell(nb)
        src = nb.cells[idx].source
        for key in ("src", "short", "full"):
            set_routing(test_routing if key == "full" else base_routing)
            base.patched_source(src, jobs[key])
        set_routing(base_routing)
        add("OK", f"notebook Settings cell (cell {idx}) patches correctly for the build runs (src and full)")
    except ImportError:
        add("FAIL", "nbformat is not installed: pip install nbformat nbclient ipykernel")
    except Exception as exc:
        set_routing(base_routing)
        add("FAIL", f"Settings cell check: {exc}")
    if not args.skip_kernel_check:
        try:
            from nbformat.v4 import new_code_cell, new_notebook
            nb = new_notebook(cells=[new_code_cell("import sys, os\nfrom pytRIBS.classes import *\nprint('KERNEL', sys.executable, os.getcwd())")])
            base.execute_nb(nb, lab, 300)
            add("OK", "notebook kernel imports pytRIBS")
        except Exception as exc:
            add("FAIL", f"notebook kernel could not import pytRIBS: {type(exc).__name__}: {str(exc)[:300]}")
    for key, rid in ids.items():
        taken = [p.name for p in (lab / "results").glob(rid + "*")] if (lab / "results").is_dir() else []
        taken += [p.name for p in lab.glob(rid + "*")]
        add("FAIL" if taken else "OK", f"name {rid} is unused" + (f" (already exists: {taken[:3]})" if taken else ""))
    return out


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="One-run trial of the .in-edit shortcut (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--cv", type=float, default=CV_TEST, help=f"test value of KINEMVELCOEF (default {CV_TEST})")
    ap.add_argument("--flowexp", type=float, default=FLOWEXP_TEST, help=f"test value of FLOWEXP (default {FLOWEXP_TEST})")
    ap.add_argument("--n", type=float, default=N_TEST, help=f"test value of CHANNELROUGHNESS (default {N_TEST})")
    ap.add_argument("--timeout-min", type=float, default=None)
    ap.add_argument("--cell-timeout", type=int, default=1800)
    ap.add_argument("--min-free-gb", type=float, default=2.0)
    ap.add_argument("--heartbeat-min", type=float, default=base.HEARTBEAT_MIN)
    ap.add_argument("--notebook", default=None)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (base.TRIBS_DEFAULT if os.path.isfile(base.TRIBS_DEFAULT)
                                                                       else shutil.which("tRIBS") or base.TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / LOG_DIR).mkdir(exist_ok=True)
    base._LOGFILE = str(lab / LOG_DIR / "trial_log.txt")
    args.notebook, nb_note = base.pick_notebook(lab, args.notebook)

    # The draw is read FIRST, while the imported module still has its W4 settings (its checks name the run with the W4 prefix).
    try:
        draw = base.read_anchors(lab, str(DRAW))[0]
    except base.Abort as exc:
        print(str(exc))
        return 1
    # Now switch the imported module to W1: window, start, length and the cold-start water table.
    base.WIN = dict(W1)
    base.BASE["GW_FRACTION"] = W1["gw"]
    base_routing = {k: base.BASE[k] for k in ROUTING_KEYS}
    test_routing = {"KINEMVELCOEF": num(args.cv), "FLOWEXP": num(args.flowexp), "CHANNELROUGHNESS": num(args.n)}
    if test_routing == base_routing:
        print("The test routing values equal the baseline values, so the trial could not tell a working shortcut from a broken one. Pick other values.")
        return 1
    ids = trial_ids(draw, test_routing["KINEMVELCOEF"], test_routing["FLOWEXP"], test_routing["CHANNELROUGHNESS"])
    jobs = {k: make_job(draw, rid) for k, rid in ids.items()}
    ref_in = lab / f"{REF_ID}.in"

    log(f"Working folder: {lab}")
    log(f"Notebook: {args.notebook}" + (f"  (note: {nb_note})" if nb_note else ""))
    log(f"Soil setting: draw {draw['idx']}: Ks x{draw['ks']:g}, f x{draw['f']:g}, PsiB x{draw['psib']:g}.  Window W1: {W1['start']}, {W1['hours']} h, water table {W1['gw']} x bedrock, Kchan off.")
    log(f"Baseline routing: {base_routing}.  TEST routing: {test_routing}")
    log("Plan: 1 build (src, baseline routing, no tRIBS) -> shortcut run (edit the .in, run tRIBS) -> 1 build (full, test routing) + tRIBS -> compare.")
    for k in ("src", "short", "full"):
        log(f"    {k:5s} {ids[k]}")
    log("  Time: 2 builds (about 1 min each) + 2 tRIBS runs (about 1.7 min each): roughly 6 min.")

    results = preflight(args, lab, jobs, ids, ref_in, test_routing, base_routing)
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

    try:
        base.take_lock()
    except base.Abort as exc:
        log(str(exc))
        return 1

    info, code = {"started": base.now_iso(), "draw": draw["idx"], "ids": ids, "baseline_routing": base_routing, "test_routing": test_routing}, 0
    log_dir = lab / LOG_DIR
    try:
        # ---- 1. build src (baseline routing), no tRIBS -------------------------------------------------
        set_routing(base_routing)
        log(f"=== 1 of 4: normal build of '{ids['src']}' at BASELINE routing (no tRIBS run) ===")
        t = time.time()
        text, warns = base.build_run(jobs["src"], args, lab, log_dir)
        base.check_in_file(lab / f"{ids['src']}.in", ref_in, jobs["src"])
        grid_note = base.check_grids(lab, jobs["src"], t)
        info["src_build_min"] = round((time.time() - t) / 60, 2)
        log(f"build OK in {base.fmt_min(time.time() - t)}; .in matches the reference except the allowed lines. {grid_note}")
        for w in warns:
            log(f"notebook said: {w}")

        # ---- 2. the shortcut: edit the .in, run tRIBS ------------------------------------------------------
        log(f"=== 2 of 4: SHORTCUT run '{ids['short']}': copy the src .in, change cv, flowexp, n and the output name, run tRIBS ===")
        t = time.time()
        src_in, short_in = lab / f"{ids['src']}.in", lab / f"{ids['short']}.in"
        new_vals = {"KINEMVELCOEF": f"{test_routing['KINEMVELCOEF']:g}", "FLOWEXP": f"{test_routing['FLOWEXP']:g}",
                    "CHANNELROUGHNESS": f"{test_routing['CHANNELROUGHNESS']:g}", "OUTFILENAME": f"results/{ids['short']}"}
        short_in.write_text(edit_in_text(src_in.read_text(), new_vals))
        n_lines = check_edit(src_in, short_in, new_vals)
        info["edit_seconds"] = round(time.time() - t, 3)
        log(f"edit OK in {info['edit_seconds']} s; {n_lines} text lines differ from the src .in (expected {len(new_vals)}: cv, flowexp, n, output name)")
        if n_lines != len(new_vals):
            log(f"WARNING: expected {len(new_vals)} differing lines, found {n_lines}.")
        set_routing(test_routing)
        base.check_in_file(short_in, ref_in, jobs["short"])        # the shortcut .in passes the same check as a normal build
        log("the shortcut .in also passes the normal-build check against the reference .in (routing values and output name as expected)")
        base.install_term_handler()
        rc, secs = base.run_tribs(args.tribs, f"{ids['short']}.in", f"simulation_{ids['short']}.log", base.timeout_for(args),
                                  progress_file=lab / "results" / f"{ids['short']}{OUTLET_SUFFIX}", hours=W1["hours"],
                                  heartbeat_min=args.heartbeat_min, label=ids["short"])
        if rc != 0:
            raise RunError(f"tRIBS exited with code {rc} on the shortcut run; see simulation_{ids['short']}.log")
        info["short_tribs_min"] = round(secs / 60, 2)
        info["short_outlet"] = base.check_outlet(lab / "results" / f"{ids['short']}{OUTLET_SUFFIX}", W1["hours"])
        log(f"shortcut tRIBS finished in {base.fmt_min(secs)}; outlet peak {info['short_outlet']['peak_m3s']:.2f} m3/s over {info['short_outlet']['rows']} rows")

        # ---- 3. full build with the test routing, run tRIBS -----------------------------------------------
        set_routing(test_routing)
        log(f"=== 3 of 4: normal FULL build of '{ids['full']}' with the test routing patched into the notebook, then tRIBS ===")
        t = time.time()
        text, warns = base.build_run(jobs["full"], args, lab, log_dir)
        base.check_in_file(lab / f"{ids['full']}.in", ref_in, jobs["full"])
        grid_note = base.check_grids(lab, jobs["full"], t)
        info["full_build_min"] = round((time.time() - t) / 60, 2)
        log(f"build OK in {base.fmt_min(time.time() - t)}; {grid_note}")
        base.install_term_handler()
        rc, secs = base.run_tribs(args.tribs, f"{ids['full']}.in", f"simulation_{ids['full']}.log", base.timeout_for(args),
                                  progress_file=lab / "results" / f"{ids['full']}{OUTLET_SUFFIX}", hours=W1["hours"],
                                  heartbeat_min=args.heartbeat_min, label=ids["full"])
        if rc != 0:
            raise RunError(f"tRIBS exited with code {rc} on the full-build run; see simulation_{ids['full']}.log")
        info["full_tribs_min"] = round(secs / 60, 2)
        info["full_outlet"] = base.check_outlet(lab / "results" / f"{ids['full']}{OUTLET_SUFFIX}", W1["hours"])
        log(f"full-build tRIBS finished in {base.fmt_min(secs)}; outlet peak {info['full_outlet']['peak_m3s']:.2f} m3/s over {info['full_outlet']['rows']} rows")

        # ---- 4. compare -----------------------------------------------------------------------------
        log("=== 4 of 4: comparison ===")
        full_in = lab / f"{ids['full']}.in"
        out_line = value_line_number(short_in, "OUTFILENAME")
        diffs = diff_lines(short_in, full_in)
        unexpected = [d for d in diffs if d[0] != out_line]
        info["in_diff_lines"] = [dict(line=d[0], shortcut=d[1], full=d[2]) for d in diffs]
        log(f".in files, shortcut vs full build: {len(diffs)} line(s) differ; the only expected one is the output-name line (line {out_line}).")
        for d in diffs:
            log(f"    line {d[0]}: shortcut {d[1]!r}   full {d[2]!r}" + ("" if d[0] == out_line else "   <-- UNEXPECTED"))
        files = compare_results(lab, ids["short"], ids["full"])
        info["results_files"] = files
        log("results files, shortcut vs full build (byte for byte):")
        for r in files:
            log(f"    {r['suffix']:34s} {r['status']:12s} {r.get('size_mb', ''):>8} MB")
        outlet = next((r for r in files if r["suffix"] == OUTLET_SUFFIX), None)
        outlet_same = bool(outlet and outlet["status"] == "identical")
        others_diff = [r["suffix"] for r in files if r["suffix"] != OUTLET_SUFFIX and r["status"] != "identical"]
        if not outlet_same:
            info["outlet_gap"] = outlet_gap(lab, ids["short"], ids["full"])
        info["outlet_identical"], info["in_identical_except_name"], info["other_files_differing"] = outlet_same, not unexpected, others_diff
        log()
        if outlet_same and not unexpected:
            verdict = "PASS: the shortcut gives a byte-identical outlet file and an .in file identical to the full build's except the output name."
            if others_diff:
                verdict += f"  Note: these other results files differ: {', '.join(others_diff)} (look at them before relying on the shortcut for anything but the outlet flow)."
        elif outlet_same:
            verdict = ("PASS WITH A NOTE: the outlet files are byte-identical, but the .in files differ in text beyond the output name "
                       "(see above), probably only in how a number is written. The shortcut result is valid; the number format in the edit may need matching.")
        else:
            verdict = f"FAIL: the outlet files differ ({info['outlet_gap']}). Do not use the shortcut until we know why."
        info["verdict"] = verdict
        log(verdict)
        t_short = info["short_tribs_min"] + info["edit_seconds"] / 60
        t_full = info["full_build_min"] + info["full_tribs_min"]
        log(f"Cost per routing run: shortcut {t_short:.2f} min (edit + tRIBS) vs full build {t_full:.2f} min (build {info['full_build_min']:.2f} + tRIBS {info['full_tribs_min']:.2f}).")
        info["cost_per_run_min"] = dict(shortcut=round(t_short, 2), full_build=round(t_full, 2))
        info["finished"] = base.now_iso()
        Path(RESULT_JSON).write_text(json.dumps(info, indent=2, default=str))
        log(f"Wrote {RESULT_JSON}. Nothing was deleted; the three .in files and the two results sets are still there.")
        code = 0 if outlet_same else 1
    except (RunError, base.RunError) as exc:
        log(f"TRIAL FAILED: {exc}")
        info["error"] = str(exc)
        Path(RESULT_JSON).write_text(json.dumps(info, indent=2, default=str))
        code = 1
    except base.Abort as exc:
        log(f"TRIAL STOPPED: {exc}")
        code = 1
    except KeyboardInterrupt:
        log("TRIAL INTERRUPTED. The lock is removed. Nothing is resumed: run the trial again after renaming or removing the W1TRIAL_ files it made.")
        code = 130
    finally:
        set_routing(base_routing)
    return code


if __name__ == "__main__":
    sys.exit(main())
