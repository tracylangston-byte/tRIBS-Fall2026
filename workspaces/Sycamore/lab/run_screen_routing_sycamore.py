#!/usr/bin/env python3
"""
run_screen_routing_sycamore.py  --  Stage 2, piece 1: the routing-knob screen on W1.

Question it answers: when ONE routing knob is moved and everything else stays at the baseline, how much does the
flood-week hydrograph change (volume, peak, arrival time)? The four routing knobs and the values screened (the design
doc's ranges table):

    cv        hillslope velocity coefficient (KINEMVELCOEF)    3   6   10   15      (baseline 6)
    rfe       flow exponent (FLOWEXP)                          0.18  0.24  0.30  0.35   (baseline 0.30)
    n         channel roughness (CHANNELROUGHNESS)             0.02  0.04  0.07  0.10   (baseline 0.04)
    wc        channel width coefficient (CHANNELWIDTHCOEFF)    1.5   2.33  3.5     (baseline 2.33)

Each knob is screened on W1 (Jan 18 - Jan 25 2017, 168 h, water table at 0.82 x bedrock), with the Oct 21 baseline for
everything else: Ks x0.2, PsiB x0.5, f x1, bedrock +1 m, channel loss off, snow on. The baseline value of each knob is
NOT run again: the existing W1 baseline run from the proxy check (W1_ks0p2_cv6_psib0p5_nocl) is the reference for all
four knobs. That leaves 11 new runs: 3 for cv, 3 for rfe, 3 for n, 2 for wc.

The script builds, runs and scores. It does not decide which knobs "matter": it prints, for every run, the flood-week
numbers and how far each moved from the baseline. Labeling a knob active, inert or trade-off against the design's
tolerances is the next step (the analysis), done on the CSV this script writes.

For each run to be made, the script:
    1. Copies the Settings cell of the Make_Sycamore_Model notebook, changes ONLY the one knob (plus the run names and
       the window), and runs the whole notebook headlessly (no browser) with nbclient. The notebook file is never changed.
    2. Checks the build: the new .in file must match the reference .in in every keyword except the ones this run may
       change, and the one screened keyword must have exactly the intended value. (The notebook's own read-back does
       not check flowexp, roughness or width, so this script does.) The bedrock / water-table / soil files are also checked.
    3. Runs tRIBS (one at a time, never two at once, in the foreground) and checks that it finished.
    4. Scores the run with score_run_sycamore.py (window W1).
    5. Writes everything to screen_routing_manifest.json, so a stopped batch can be resumed by running the same command.

The reference run is checked first (file present, complete, checksum and proxy manifest agree, its .in file has the
baseline routing values). If it fails any check the script makes a fresh baseline under a new name
(W1_ks0p2_cv6_psib0p5_nocl_screenbase); it never overwrites the old one.

Run it from workspaces/Sycamore/lab/. Typical use:
    python -u run_screen_routing_sycamore.py --dry-run                 # checks only, builds nothing
    python -u run_screen_routing_sycamore.py --limit 1 --build-only    # build the first run and check it, no tRIBS
    python -u run_screen_routing_sycamore.py --limit 1                 # first run, all steps
    python -u run_screen_routing_sycamore.py                           # everything; run it again to resume
    python -u run_screen_routing_sycamore.py --summary-only            # just print the table from what is finished

Files it makes (all in lab/): screen_routing_manifest.json, screen_routing.lock (only while running),
screen_routing_logs/, score_out_screen_routing/scores_screen_routing.csv, screen_routing_summary.csv, and for each run
<run>.in, simulation_<run>.log and results/<run>*.

Do NOT run the notebook (or any cell that rebuilds data/model/) while this runs: every build rewrites lab/data/model/,
and all the runs share it. Keep the terminal open: the progress line printed every 5 minutes also stops the Codespace
from timing out for being idle.
"""
from __future__ import annotations

import argparse
import ast
import atexit
import csv
import datetime as dt
import filecmp
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------------------------
# What the screen does
# --------------------------------------------------------------------------------------------
W1 = dict(start="01/18/2017/00/00", hours=168, label="jan_storm_week", gw=0.82)    # label = the notebook's RUNS key

# Everything that is not being screened (the Oct 21 baseline, with the W1 cold-start water table).
BASE = {
    "GW_FRACTION": W1["gw"], "BEDROCK_ADD_M": 1.0, "KS_MULT": 0.2, "PSIB_MULT": 0.5, "F_MULT": 1,
    "KINEMVELCOEF": 6, "FLOWEXP": 0.3, "CHANNELROUGHNESS": 0.04, "CHANNELWIDTHCOEFF": 2.33,
    "OPTSNOW": 1, "OPTPERCOLATION": 0, "CHANNELCONDUCTIVITY": 5,
    "REBUILD_SOIL": False,
}
# key = what you type after --knobs; tag = what goes in the run name; kw = the Settings name AND the .in keyword.
KNOBS = [
    dict(key="cv",  kw="KINEMVELCOEF",      tag="cv",  base=6,    values=[3, 6, 10, 15],
         label="hillslope velocity coefficient"),
    dict(key="rfe", kw="FLOWEXP",           tag="rfe", base=0.3,  values=[0.18, 0.24, 0.30, 0.35],
         label="flow exponent"),
    dict(key="n",   kw="CHANNELROUGHNESS",  tag="n",   base=0.04, values=[0.02, 0.04, 0.07, 0.10],
         label="channel roughness (Manning n)"),
    dict(key="wc",  kw="CHANNELWIDTHCOEFF", tag="wc",  base=2.33, values=[1.5, 2.33, 3.5],
         label="channel width coefficient"),
]
ROUTING_KWS = [k["kw"] for k in KNOBS]

BASELINE_ID = "W1_ks0p2_cv6_psib0p5_nocl"                  # the existing W1 baseline run (from the proxy check)
FRESH_BASELINE_ID = BASELINE_ID + "_screenbase"            # made only if the existing one fails its checks
W4_BASELINE_IN = "Sycamore_octmar_ks0p2_cv6_psib0p5_nocl.in"   # fallback reference .in if the W1 baseline is not usable
PROXY_MANIFEST = "proxy_manifest.json"
NOTEBOOK_CANDIDATES = ["Make_Sycamore_Model.ipynb", "Make_Sycamore_Model_clean.ipynb"]
SCORER = "score_run_sycamore.py"
TRIBS_DEFAULT = "/opt/tRIBS/build/tRIBS"
INIT_DATA = "../init_data"           # SOLUS_DIR in the notebook
MODEL_SOIL = "data/model/soil"
LOG_DIR = "screen_routing_logs"
SCORE_DIR = "score_out_screen_routing"
SCORE_CSV = "scores_screen_routing.csv"
SUMMARY_CSV = "screen_routing_summary.csv"
MANIFEST = "screen_routing_manifest.json"
LOCK = "screen_routing.lock"
MIN_PER_1000H = 7.8          # measured tRIBS speed on this machine: minutes per 1,000 model hours
HEARTBEAT_MIN = 5.0
ROWS_PER_MODEL_HOUR = 16     # the outlet file gets 16 rows per model hour

# Score columns shown in the table (all are columns of the scorer's CSV).
METRICS = [("vol", "flood_week__vol_ratio"), ("peak", "flood_week__peak_err_pct"),
           ("arr", "flood_week__arrival_err_h"), ("pkdt", "flood_week__peak_time_diff_h"),
           ("kge", "flood_week__kge2012"),
           ("s2vol", "jan_second_storm__vol_ratio"), ("s2peak", "jan_second_storm__peak_err_pct"),
           ("s2arr", "jan_second_storm__arrival_err_h")]
SHOW_COLS = [c for _, c in METRICS]

_LOGFILE = None


def log(msg: str = "") -> None:
    stamp = dt.datetime.now().strftime("%H:%M:%S")
    line = f"[{stamp}] {msg}" if msg else ""
    print(line, flush=True)
    if _LOGFILE is not None:
        try:
            with open(_LOGFILE, "a") as fh:
                fh.write(line + "\n")
        except OSError:
            pass


def now_iso() -> str:
    return dt.datetime.now().isoformat(timespec="seconds")


def fmt_min(seconds: float) -> str:
    return f"{seconds / 60:.1f} min"


class RunError(Exception):
    """A problem with one run. The batch records it and moves on to the next run."""


class Abort(Exception):
    """A problem with the whole batch (competing job, no disk). The batch stops."""


# --------------------------------------------------------------------------------------------
# Names and settings for each job
# --------------------------------------------------------------------------------------------
def vtag(x) -> str:
    """0.18 -> '0p18', 1.5 -> '1p5', 15 -> '15', 0.10 -> '0p1'."""
    return f"{x:g}".replace(".", "p")


def run_id_for(knob: dict | None, value=None) -> str:
    """cv goes in its usual slot (W1_ks0p2_cv3_psib0p5_nocl); the other knobs go after psib, before the loss setting."""
    if knob is None:
        return FRESH_BASELINE_ID
    if knob["key"] == "cv":
        return f"W1_ks0p2_cv{vtag(value)}_psib0p5_nocl"
    return f"W1_ks0p2_cv6_psib0p5_{knob['tag']}{vtag(value)}_nocl"


def make_job(knob: dict | None, value=None) -> dict:
    if knob is None:
        return dict(knob=None, kw=None, value=None, run_id=FRESH_BASELINE_ID, kind="baseline")
    return dict(knob=knob["key"], kw=knob["kw"], value=value, run_id=run_id_for(knob, value), kind="run")


def knob_by_key(key: str) -> dict:
    return next(k for k in KNOBS if k["key"] == key)


def screened_values(knob: dict) -> list:
    """The values that need a new run: every screened value except the baseline one."""
    return [v for v in knob["values"] if v != knob["base"]]


def settings_for(job: dict) -> dict:
    """Every NAME = value line of the Settings cell that this run sets."""
    s = dict(BASE)
    s.update(RUN_CHOICE=W1["label"], INPUT_FILE=f"{job['run_id']}.in", RESULTS_NAME=job["run_id"])
    if job["kw"]:
        s[job["kw"]] = job["value"]
    return s


def settings_record(job: dict) -> dict:
    """What the manifest remembers about the settings, to notice a changed setting on resume."""
    s = settings_for(job)
    rec = {k: s[k] for k in ("GW_FRACTION", "BEDROCK_ADD_M", "KS_MULT", "PSIB_MULT", "F_MULT", "KINEMVELCOEF", "FLOWEXP",
                             "CHANNELROUGHNESS", "CHANNELWIDTHCOEFF", "OPTSNOW", "OPTPERCOLATION", "CHANNELCONDUCTIVITY")}
    rec.update(START_DATE=W1["start"], RUN_HOURS=W1["hours"])
    return rec


def expected_files(s):
    """Names the notebook gives the bedrock and water-table grids for these settings."""
    bed = s["BEDROCK_ADD_M"]
    gw_name = f"InitGW_{int(round(s['GW_FRACTION'] * 100))}pct_plus{bed:g}m_mm.asc"
    bed_name = "bedrock_depth_m.asc" if bed == 0 else f"bedrock_depth_plus{bed:g}m.asc"
    return gw_name, bed_name


def timeout_for(args) -> float:
    if args.timeout_min:
        return args.timeout_min
    return max(10.0, 3.0 * W1["hours"] * MIN_PER_1000H / 1000.0)


# --------------------------------------------------------------------------------------------
# Patching the Settings cell
# --------------------------------------------------------------------------------------------
def find_settings_cell(nb):
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        src = cell.source
        if "# SETTINGS: everything specific" in src or ("GW_FRACTION = " in src and "INPUT_FILE = " in src):
            return i
    raise RunError("Could not find the Settings cell in the notebook.")


def patch_settings(src: str, values: dict) -> str:
    """Change only the listed NAME = value lines (keeps the trailing # comment)."""
    for name, value in values.items():
        pat = re.compile(rf"^({re.escape(name)}[ \t]*=[ \t]*)(.*?)([ \t]*(?:#.*)?)$", re.M)
        hits = pat.findall(src)
        if len(hits) != 1:
            raise RunError(f"Settings cell: expected exactly one line '{name} = ...' but found {len(hits)}.")
        src = pat.sub(lambda m, v=value: m.group(1) + repr(v) + m.group(3), src, count=1)
    return src


def patch_window(src: str, start: str, hours: int) -> str:
    """Replace the line 'START_DATE, RUN_HOURS = RUNS[RUN_CHOICE]' by the window's own start and length."""
    pat = re.compile(r"^(START_DATE[ \t]*,[ \t]*RUN_HOURS[ \t]*=[ \t]*)(.*?)([ \t]*(?:#.*)?)$", re.M)
    hits = pat.findall(src)
    if len(hits) != 1:
        raise RunError(f"Settings cell: expected exactly one line 'START_DATE, RUN_HOURS = ...' but found {len(hits)}.")
    return pat.sub(lambda m: m.group(1) + repr((start, hours)) + m.group(3), src, count=1)


def verify_patched(src: str, values: dict, window: tuple | None = None) -> None:
    """Parse the patched Settings cell and confirm each patched name is assigned once, to the right value."""
    tree = ast.parse(src)
    found, win = {}, []
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        tgt = node.targets[0]
        if isinstance(tgt, ast.Name) and tgt.id in values:
            found.setdefault(tgt.id, []).append(ast.literal_eval(node.value))
        elif isinstance(tgt, ast.Tuple) and [getattr(e, "id", None) for e in tgt.elts] == ["START_DATE", "RUN_HOURS"]:
            try:
                win.append(tuple(ast.literal_eval(node.value)))
            except ValueError:
                win.append(None)
    for name, want in values.items():
        got = found.get(name, [])
        if len(got) != 1 or got[0] != want or type(got[0]) is not type(want):
            raise RunError(f"Settings cell after patching: {name} = {got}, wanted exactly one value {want!r}.")
    if window is not None and win != [tuple(window)]:
        raise RunError(f"Settings cell after patching: START_DATE, RUN_HOURS = {win}, wanted exactly {tuple(window)!r}.")


def patched_source(src: str, job: dict) -> str:
    s = settings_for(job)
    out = patch_window(patch_settings(src, s), W1["start"], W1["hours"])
    verify_patched(out, s, (W1["start"], W1["hours"]))
    return out


# --------------------------------------------------------------------------------------------
# Running the notebook headlessly
# --------------------------------------------------------------------------------------------
def make_kernel_manager():
    """A kernel that uses THIS python (so it sees the same packages as this script)."""
    from jupyter_client.kernelspec import KernelSpec
    from jupyter_client.manager import KernelManager
    km = KernelManager(kernel_name="python3")
    km._kernel_spec = KernelSpec(argv=[sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
                                 display_name="screen", language="python")
    return km


def execute_nb(nb, lab: Path, cell_timeout: int):
    """Execute all cells of an nbformat notebook in a fresh kernel started in lab/. Returns the notebook."""
    import nbclient
    km = make_kernel_manager()
    client = nbclient.NotebookClient(nb, km=km, timeout=cell_timeout, startup_timeout=180,
                                     allow_errors=False, resources={"metadata": {"path": str(lab)}})
    try:
        client.execute()
    finally:
        try:
            km.shutdown_kernel(now=True)
        except Exception:
            pass
        try:
            km.cleanup_resources()
        except Exception:
            pass
    return nb


def outputs_text(nb) -> str:
    """All text a notebook printed (images dropped)."""
    chunks = []
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        for out in cell.get("outputs", []):
            t = out.get("output_type")
            if t == "stream":
                chunks.append(out.get("text", ""))
            elif t in ("execute_result", "display_data"):
                txt = out.get("data", {}).get("text/plain")
                if txt:
                    chunks.append(txt + "\n")
            elif t == "error":
                chunks.append("\n".join(out.get("traceback", [])) + "\n")
    return "".join(chunks)


def build_run(job: dict, args, lab: Path, log_dir: Path) -> tuple:
    """Patch + execute the notebook. Returns (output_text, warnings). Raises RunError on any problem."""
    import nbformat
    rid = job["run_id"]
    nb = nbformat.read(str(lab / args.notebook), as_version=4)
    for cell in nb.cells:                               # start from empty outputs
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    idx = find_settings_cell(nb)
    nb.cells[idx].source = patched_source(nb.cells[idx].source, job)
    try:
        execute_nb(nb, lab, args.cell_timeout)
    except Exception as exc:                            # a cell failed, or the kernel did not start
        failed = log_dir / f"{rid}_build_FAILED_executed.ipynb"
        try:
            nbformat.write(nb, str(failed))
        except Exception:
            pass
        partial = outputs_text(nb)
        (log_dir / f"{rid}_build_output.txt").write_text(partial)
        tail = re.sub(r"\x1b\[[0-9;]*m", "", "\n".join(partial.splitlines()[-12:]))
        raise RunError(f"The notebook stopped with an error: {type(exc).__name__}: {str(exc)[:600]}\n"
                       f"  (partial notebook saved as {failed.name}; last output lines follow)\n{tail}")
    text = outputs_text(nb)
    (log_dir / f"{rid}_build_output.txt").write_text(text)
    return text, check_build_text(text, job)


def check_build_text(text: str, job: dict) -> list:
    """Lines the notebook should have printed. Returns a list of warnings; raises RunError on a problem."""
    s = settings_for(job)
    gw_name, bed_name = expected_files(s)
    must = [
        f"Run choice: {s['RUN_CHOICE']} | start {W1['start']} | {W1['hours']} hours",
        "Every expected setting is in the input file.",
        "OK: every hourly file the run needs is there.",
        f"(KS_MULT = {s['KS_MULT']:g},", f"(PSIB_MULT = {s['PSIB_MULT']:g},", f"(F_MULT = {s['F_MULT']:g},",
        f"Bedrock file in use: {INIT_DATA}/{bed_name}",
        f"GWATERFILE: {MODEL_SOIL}/{gw_name}",
        f"BEDROCKFILE: {MODEL_SOIL}/{bed_name}",
    ]
    missing = [m for m in must if m not in text]
    if not re.search(r"All \d+ checks passed\.", text):
        missing.append("All <n> checks passed.")
    if missing:
        raise RunError("The notebook output is missing line(s) it should have printed:\n  " + "\n  ".join(missing))
    if re.search(r"^PROBLEM:", text, re.M):
        raise RunError("The notebook printed a PROBLEM line (see the build output file).")
    return [ln.strip() for ln in text.splitlines() if re.match(r"^CHECK\s", ln)]


# --------------------------------------------------------------------------------------------
# Checks on what the build produced
# --------------------------------------------------------------------------------------------
def read_in_file(path: Path) -> dict:
    lines = Path(path).read_text().splitlines()
    found = {}
    for i, line in enumerate(lines):
        m = re.match(r"^([A-Z0-9]+):\s", line + " ")
        if m and i + 1 < len(lines):
            found[m.group(1)] = lines[i + 1].strip()
    return found


def same_value(a, b) -> bool:
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return str(a).strip() == str(b).strip()


def check_in_file(new_in: Path, ref_in: Path, job: dict) -> None:
    """The new .in must equal the reference .in except for the run-specific lines, and the routing keywords must have
    exactly their intended values (the baseline value for three knobs, the screened value for the fourth)."""
    s = settings_for(job)
    gw_name, _ = expected_files(s)
    new, ref = read_in_file(new_in), read_in_file(ref_in)
    expect = {"RUNTIME": W1["hours"], "STARTDATE": W1["start"], "OUTFILENAME": f"results/{job['run_id']}",
              "GWATERFILE": f"{MODEL_SOIL}/{gw_name}"}
    for kw in ROUTING_KWS:
        expect[kw] = s[kw]
    problems = []
    for key in sorted(set(ref) | set(new)):
        b, n = ref.get(key), new.get(key)
        if key in expect:
            if n is None or not same_value(n, expect[key]):
                problems.append(f"{key}: written {n!r}, expected {expect[key]!r}")
        elif b is None or n is None or not same_value(b, n):
            problems.append(f"{key}: reference {b!r} but this run {n!r}  (this keyword may not change)")
    if problems:
        raise RunError("The .in file differs from the reference in more than it should:\n  " + "\n  ".join(problems))


def read_asc(path: Path):
    """Read an ESRI ASCII grid. Returns (header dict, 2-D array)."""
    header, nhead = {}, 0
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 2 and parts[0].lower() in ("ncols", "nrows", "xllcorner", "yllcorner", "xllcenter",
                                                        "yllcenter", "cellsize", "nodata_value"):
                header[parts[0].lower()] = float(parts[1])
                nhead += 1
            else:
                break
    data = np.loadtxt(path, skiprows=nhead, dtype=float)
    return header, data.reshape(int(header["nrows"]), int(header["ncols"]))


def check_grids(lab: Path, job: dict, t_build_start: float) -> str:
    """Bedrock, water-table and soil-multiplier files in data/model/soil must be fresh and correct."""
    s = settings_for(job)
    gw_name, bed_name = expected_files(s)
    soil = lab / MODEL_SOIL
    gw_path, bed_path = soil / gw_name, soil / bed_name
    for p in (gw_path, bed_path):
        if not p.is_file():
            raise RunError(f"Expected file not found: {p}")
    if gw_path.stat().st_mtime < t_build_start - 2:
        raise RunError(f"{gw_name} in data/model/soil is older than this build (stale file).")
    if not filecmp.cmp(bed_path, lab / INIT_DATA / bed_name, shallow=False):
        raise RunError(f"{bed_name} in data/model/soil differs from the copy in {INIT_DATA}.")
    hb, bed = read_asc(bed_path)
    hg, gw = read_asc(gw_path)
    nb_, ng_ = hb.get("nodata_value", -9999), hg.get("nodata_value", -9999)
    if bed.shape != gw.shape:
        raise RunError("Bedrock and water-table grids have different shapes.")
    ok = (bed != nb_) & (gw != ng_) & (bed > 0)
    if not ok.any():
        raise RunError("No valid cells to compare in the bedrock and water-table grids.")
    orig_h, orig = read_asc(lab / INIT_DATA / "bedrock_depth_m.asc")
    ov = (orig != orig_h.get("nodata_value", -9999)) & (bed != nb_)
    diff = (bed - orig)[ov]
    if diff.size == 0 or np.max(np.abs(diff - s["BEDROCK_ADD_M"])) > 0.011:
        raise RunError(f"Bedrock grid is not the SOLUS grid + {s['BEDROCK_ADD_M']:g} m "
                       f"(difference ranges {diff.min():.3f} to {diff.max():.3f} m).")
    med = float(np.median(gw[ok] / bed[ok]))
    if abs(med - s["GW_FRACTION"] * 1000) <= 5:
        unit = "mm"
    elif abs(med - s["GW_FRACTION"]) <= 0.005:
        unit = "m"
    else:
        raise RunError(f"Water-table grid / bedrock grid has median {med:.4f}; expected "
                       f"{s['GW_FRACTION'] * 1000:g} (mm per m) for GW_FRACTION = {s['GW_FRACTION']}.")
    with open(soil / "scgrid.gdf") as fh:
        rows = {r["Variable"].strip(): r["BasePath"].strip() for r in csv.DictReader(fh)}
    for var, mult in (("KS", s["KS_MULT"]), ("PB", s["PSIB_MULT"]), ("FD", s["F_MULT"])):
        base = rows.get(var, "")
        tg = f"_x{mult:g}"
        if (mult != 1 and not base.endswith(tg)) or (mult == 1 and "_x" in os.path.basename(base)):
            raise RunError(f"scgrid.gdf row {var} is '{base}', which does not match multiplier {mult:g}.")
    return (f"Ks x{s['KS_MULT']:g}, PsiB x{s['PSIB_MULT']:g}; water table = {s['GW_FRACTION']:g} x bedrock ({unit}), "
            f"median bedrock {np.median(bed[ok]):.2f} m, median water-table depth {np.median(gw[ok]):.0f} {unit}")


# --------------------------------------------------------------------------------------------
# Running tRIBS and scoring
# --------------------------------------------------------------------------------------------
def competing_tribs() -> list:
    out = subprocess.run(["ps", "-eo", "pid=,args="], capture_output=True, text=True).stdout
    hits = []
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        pid, _, rest = line.partition(" ")
        words = rest.split()
        if pid.isdigit() and int(pid) != os.getpid() and words and os.path.basename(words[0]) == "tRIBS":
            hits.append(line)
    return hits


def disk_free_gb(path=".") -> float:
    return shutil.disk_usage(path).free / 1e9


def count_rows(path: Path) -> int:
    try:
        with open(path, "rb") as fh:
            return sum(chunk.count(b"\n") for chunk in iter(lambda: fh.read(1 << 20), b""))
    except OSError:
        return 0


def kill_group(proc) -> None:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(timeout=10)
    except (ProcessLookupError, PermissionError):
        pass


def _on_term(signum, frame):
    raise KeyboardInterrupt


def install_term_handler() -> None:
    """Make 'kill' (SIGTERM) behave like Ctrl-C: stop tRIBS, mark the run interrupted, remove the lock. Running a notebook
    resets this handler (nbclient does it), so it is installed again before every tRIBS run."""
    signal.signal(signal.SIGTERM, _on_term)


def run_tribs(tribs: str, infile: str, logname: str, timeout_min: float,
              progress_file: Path | None = None, hours: int | None = None,
              heartbeat_min: float = HEARTBEAT_MIN, label: str = "") -> tuple:
    """Run tRIBS to completion. Returns (return code, seconds). Kills the whole process group on timeout.

    tRIBS prints nothing while it runs. Every `heartbeat_min` minutes this prints a progress line, which also counts as
    terminal activity, so a Codespace (30 minute idle timeout) does not shut down in the middle of a run.
    """
    t0 = time.time()
    deadline = t0 + timeout_min * 60
    beat = max(heartbeat_min, 0.01) * 60
    expected = hours * ROWS_PER_MODEL_HOUR if hours else None
    with open(logname, "w") as fh:
        proc = subprocess.Popen([tribs, infile], stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            while True:
                try:
                    rc = proc.wait(timeout=max(0.1, min(beat, deadline - time.time())))
                    break
                except subprocess.TimeoutExpired:
                    if time.time() >= deadline:
                        kill_group(proc)
                        raise RunError(f"tRIBS did not finish within {timeout_min:g} minutes and was stopped.")
                    msg = f"{label}: tRIBS still running, {(time.time() - t0) / 60:.1f} min so far"
                    if progress_file is not None:
                        n = count_rows(progress_file)
                        msg += f"; outlet file has {n} rows" + (f" of about {expected} ({100 * n / expected:.0f}%)" if expected else "")
                    log(msg)
        except BaseException:
            kill_group(proc)
            raise
    return rc, time.time() - t0


def check_outlet(path: Path, hours: int) -> dict:
    if not path.is_file():
        raise RunError(f"tRIBS did not write {path}.")
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
    if len(rows) < 0.99 * hours:
        raise RunError(f"{path.name} has only {len(rows)} rows; expected about {hours}.")
    arr = np.array(rows)
    if not np.isfinite(arr[:, 1]).all():
        raise RunError(f"{path.name} contains NaN or infinite flows.")
    if arr[-1, 0] < hours - 2:
        raise RunError(f"{path.name} stops at hour {arr[-1, 0]:g}; the run should reach hour {hours}.")
    return {"rows": len(rows), "last_hour": float(arr[-1, 0]), "peak_m3s": float(arr[:, 1].max()),
            "md5": hashlib.md5(path.read_bytes()).hexdigest()}


def run_scorer(rid: str, scorer: str, score_dir: Path, log_dir: Path, write_checksum: bool = True) -> dict:
    """Score one W1 run; returns its whole row of the scorer's CSV (all values as text)."""
    cmd = [sys.executable, scorer, "score", "--run", rid, "--window", "W1", "--results-dir", "results",
           "--out-dir", str(score_dir), "--csv-name", SCORE_CSV]
    if write_checksum:
        cmd.append("--write-checksum")
    res = subprocess.run(cmd, capture_output=True, text=True)
    (log_dir / f"{rid}_score.txt").write_text((res.stdout or "") + (res.stderr or ""))
    if res.returncode != 0:
        raise RunError(f"The scorer exited with code {res.returncode}: {(res.stderr or res.stdout)[-400:]}")
    csv_path = score_dir / SCORE_CSV
    if csv_path.is_file():
        with open(csv_path, newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("run_id") == rid:
                    return dict(row)
    raise RunError(f"The scorer ran but {rid} is not in {csv_path}.")


def fnum(row, col) -> float:
    try:
        v = float((row or {}).get(col, ""))
        return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def show_scores(row: dict) -> dict:
    out = {}
    for c in SHOW_COLS:
        v = fnum(row, c)
        if math.isfinite(v):
            out[c] = round(v, 3)
    return out


def one_line(scores: dict) -> str:
    names = {"flood_week__vol_ratio": "flood vol ratio", "flood_week__peak_err_pct": "peak err %",
             "flood_week__arrival_err_h": "arrival err h"}
    return ", ".join(f"{names[c]} {scores[c]}" for c in names if c in scores)


# --------------------------------------------------------------------------------------------
# Manifest, lock
# --------------------------------------------------------------------------------------------
def load_manifest() -> dict:
    if Path(MANIFEST).is_file():
        try:
            return json.loads(Path(MANIFEST).read_text())
        except json.JSONDecodeError:
            shutil.copy(MANIFEST, MANIFEST + ".corrupt")
    return {"runs": {}}


def save_manifest(man: dict) -> None:
    tmp = MANIFEST + ".tmp"
    Path(tmp).write_text(json.dumps(man, indent=2))
    os.replace(tmp, MANIFEST)


def take_lock() -> None:
    p = Path(LOCK)
    if p.is_file():
        try:
            pid = int(p.read_text().strip())
            os.kill(pid, 0)
            raise Abort(f"Another copy of this script is running (process {pid}). Lock file: {LOCK}")
        except (ValueError, ProcessLookupError, PermissionError):
            pass                                        # stale lock from a stopped run
    p.write_text(str(os.getpid()))
    atexit.register(lambda: p.unlink() if p.is_file() and p.read_text().strip() == str(os.getpid()) else None)


# --------------------------------------------------------------------------------------------
# The reference run (the existing W1 baseline)
# --------------------------------------------------------------------------------------------
def verify_baseline_in(path: Path) -> None:
    """The W1 baseline's .in file must carry the baseline routing values, loss off, snow on, and the W1 window."""
    got = read_in_file(path)
    gw_name, _ = expected_files(dict(BASE))
    want = {"STARTDATE": W1["start"], "RUNTIME": W1["hours"], "OPTPERCOLATION": BASE["OPTPERCOLATION"],
            "OPTSNOW": BASE["OPTSNOW"], "CHANNELCONDUCTIVITY": BASE["CHANNELCONDUCTIVITY"],
            "KINEMVELCOEF": BASE["KINEMVELCOEF"], "FLOWEXP": BASE["FLOWEXP"],
            "CHANNELROUGHNESS": BASE["CHANNELROUGHNESS"], "CHANNELWIDTHCOEFF": BASE["CHANNELWIDTHCOEFF"]}
    bad = [f"{k} is {got.get(k)!r}, wanted {v!r}" for k, v in want.items() if k not in got or not same_value(got[k], v)]
    gwf = got.get("GWATERFILE", "")
    if not gwf.endswith(gw_name):
        bad.append(f"GWATERFILE is {gwf!r}, wanted a file ending in {gw_name}")
    if bad:
        raise RunError(f"{path.name} is not the baseline: " + "; ".join(bad))


def baseline_check(lab: Path) -> tuple:
    """Read-only check of the existing W1 baseline. Returns (problems, notes). No problems = it can be the reference."""
    rid, problems, notes = BASELINE_ID, [], []
    q = lab / "results" / f"{rid}_Outlet.qout"
    if not q.is_file():
        return [f"results/{q.name} not found"], notes
    try:
        outlet = check_outlet(q, W1["hours"])
    except RunError as exc:
        return [str(exc)], notes
    side = q.with_suffix(q.suffix + ".sha256")
    if side.is_file():
        saved = side.read_text().split()[0]
        if saved != hashlib.sha256(q.read_bytes()).hexdigest():
            problems.append(f"{q.name} no longer matches its .sha256 checksum file (the kept file has changed)")
        else:
            notes.append("checksum file matches")
    else:
        notes.append("no .sha256 file yet (it is written when the baseline is scored)")
    pm = Path(PROXY_MANIFEST)
    if pm.is_file():
        try:
            prec = json.loads(pm.read_text())["runs"].get(rid)
        except Exception:                                                  # unreadable manifest is not a reason to reject
            prec = None
            notes.append(f"{PROXY_MANIFEST} could not be read")
        if prec:
            if prec.get("status") != "done":
                problems.append(f"{PROXY_MANIFEST} says this run's status is {prec.get('status')!r}, not 'done'")
            md5 = (prec.get("outlet") or {}).get("md5")
            if md5 and md5 != outlet["md5"]:
                problems.append(f"the outlet file's md5 ({outlet['md5']}) differs from the one recorded in {PROXY_MANIFEST} ({md5})")
            elif md5:
                notes.append(f"md5 matches {PROXY_MANIFEST}")
            want = {"GW_FRACTION": W1["gw"], "BEDROCK_ADD_M": BASE["BEDROCK_ADD_M"], "KS_MULT": BASE["KS_MULT"],
                    "PSIB_MULT": BASE["PSIB_MULT"], "START_DATE": W1["start"], "RUN_HOURS": W1["hours"]}
            if prec.get("settings") and prec["settings"] != want:
                problems.append(f"{PROXY_MANIFEST} records settings {prec['settings']}, expected {want}")
        elif prec is None and "could not be read" not in " ".join(notes):
            notes.append(f"{PROXY_MANIFEST} has no record of this run")
    else:
        notes.append(f"{PROXY_MANIFEST} not found, so the run's origin is checked only by its file, checksum and .in file")
    inp = lab / f"{rid}.in"
    if not inp.is_file():
        problems.append(f"{inp.name} not found in lab/")
    else:
        try:
            verify_baseline_in(inp)
            notes.append(f"{inp.name} has the baseline routing values")
        except RunError as exc:
            problems.append(str(exc))
    return problems, notes


def do_baseline(args, lab: Path, man: dict) -> bool:
    """Check and score the existing W1 baseline into the screen's own score CSV. Returns True if it can be the reference."""
    rid = BASELINE_ID
    log_dir, score_dir = lab / LOG_DIR, lab / SCORE_DIR
    log_dir.mkdir(exist_ok=True)
    score_dir.mkdir(exist_ok=True)
    rec = man["runs"].setdefault(rid, {})
    rec.clear()
    rec.update(kind="reuse", knob=None, status="checking", started=now_iso())
    problems, notes = baseline_check(lab)
    if not problems:
        try:
            q = lab / "results" / f"{rid}_Outlet.qout"
            has_sidecar = q.with_suffix(q.suffix + ".sha256").is_file()
            row = run_scorer(rid, str(lab / args.scorer), score_dir, log_dir, write_checksum=not has_sidecar)
            rec["scores"] = show_scores(row)
            rec["outlet"] = check_outlet(q, W1["hours"])
        except RunError as exc:
            problems.append(str(exc))
    if problems:
        rec.update(status="reuse_rejected", problems=problems, finished=now_iso())
        man["baseline_id"] = FRESH_BASELINE_ID
        save_manifest(man)
        log(f"{rid}: CANNOT be the reference run:")
        for p in problems:
            log(f"    {p}")
        log(f"    A fresh baseline will be made as {FRESH_BASELINE_ID}; {rid} itself is left untouched.")
        return False
    rec.update(status="reused", finished=now_iso(), settings=settings_record(dict(make_job(None), run_id=rid)))
    man["baseline_id"] = rid
    save_manifest(man)
    for n in notes:
        log(f"{rid}: {n}")
    log(f"{rid}: reused as the reference run. {one_line(rec['scores'])}")
    return True


# --------------------------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------------------------
def pick_notebook(lab: Path, given: str | None) -> tuple:
    """(notebook name, note). With no --notebook, Make_Sycamore_Model.ipynb if it exists (that is what built the baseline),
    otherwise Make_Sycamore_Model_clean.ipynb."""
    if given:
        return given, ""
    present = [n for n in NOTEBOOK_CANDIDATES if (lab / n).is_file()]
    if not present:
        return NOTEBOOK_CANDIDATES[0], ""
    note = ""
    if len(present) > 1:
        note = (f"both {present[0]} and {present[1]} are in lab/; using {present[0]} (the one the proxy runner uses). "
                f"They should have the same Settings cell. Use --notebook to choose the other.")
    return present[0], note


def preflight(args, lab: Path, jobs: list, base_ok: bool, nb_note: str) -> list:
    """Checks before anything is built. Returns a list of (level, message)."""
    out = []

    def add(level, msg):
        out.append((level, msg))

    for f in (args.notebook, args.scorer):
        add("OK" if (lab / f).is_file() else "FAIL", f"file in lab/: {f}")
    if nb_note:
        add("WARN", nb_note)
    w4 = lab / args.baseline_in
    if base_ok:
        add("OK" if w4.is_file() else "WARN", f"fallback reference .in (not needed, the W1 baseline is usable): {args.baseline_in}")
    else:
        add("OK" if w4.is_file() else "FAIL", f"reference .in for a fresh baseline: {args.baseline_in}")
    obs = lab / INIT_DATA / "obs" / "sycamore_09510200_obs_flow_WY2017.csv"
    add("OK" if obs.is_file() else "FAIL", f"gauge file: {INIT_DATA}/obs/{obs.name}")
    for f in ("../init_data/bedrock_depth_m.asc", "../init_data/met/Master_Met.sdf", "../init_data/LandUse.asc"):
        add("OK" if (lab / f).is_file() else "FAIL", f"input: {f}")
    tribs = args.tribs
    add("OK" if (os.path.isfile(tribs) and os.access(tribs, os.X_OK)) else "FAIL", f"tRIBS program: {tribs}")
    free = disk_free_gb(lab)
    add("OK" if free >= args.min_free_gb else "FAIL", f"free disk space: {free:.1f} GB (need {args.min_free_gb:g})")
    comp = competing_tribs()
    add("FAIL" if comp else "OK", "no other tRIBS running" + (": " + comp[0] if comp else ""))
    try:
        import nbformat
        nb = nbformat.read(str(lab / args.notebook), as_version=4)
        idx = find_settings_cell(nb)
        every = [make_job(knob_by_key(k["key"]), v) for k in KNOBS for v in screened_values(k)] + [make_job(None)]
        for j in every:
            patched_source(nb.cells[idx].source, j)
        add("OK", f"Settings cell found (cell {idx}) and patchable for all {len(every)} runs (11 screened values + a spare baseline)")
    except ImportError:
        add("FAIL", "nbformat is not installed: pip install nbformat nbclient ipykernel")
    except Exception as exc:
        add("FAIL", f"Settings cell check: {exc}")
    # the reference .in has every keyword the comparison needs
    ref = lab / (f"{BASELINE_ID}.in" if base_ok else args.baseline_in)
    try:
        got = read_in_file(ref)
        need = ROUTING_KWS + ["OPTPERCOLATION", "OPTSNOW", "CHANNELCONDUCTIVITY", "RUNTIME", "STARTDATE", "GWATERFILE", "OUTFILENAME"]
        miss = [k for k in need if k not in got]
        add("FAIL" if miss else "OK", f"reference .in ({ref.name}) has the keywords the comparison needs" + (f" (missing: {miss})" if miss else ""))
    except Exception as exc:
        add("FAIL", f"reference .in: {exc}")
    try:
        res = subprocess.run([sys.executable, str(lab / args.scorer), "score", "--help"], capture_output=True, text=True)
        flags = ["--run", "--window", "--results-dir", "--out-dir", "--csv-name", "--write-checksum"]
        miss = [f for f in flags if f not in res.stdout]
        add("FAIL" if (res.returncode != 0 or miss) else "OK",
            "scorer accepts " + " ".join(flags) + (f" (missing: {miss})" if miss else ""))
    except Exception as exc:
        add("FAIL", f"scorer check: {exc}")
    if not args.skip_kernel_check:
        try:
            from nbformat.v4 import new_code_cell, new_notebook
            nb = new_notebook(cells=[new_code_cell(
                "import sys, os\nfrom pytRIBS.classes import *\nprint('KERNEL', sys.executable, os.getcwd())")])
            execute_nb(nb, lab, 300)
            txt = outputs_text(nb).strip()
            add("OK", f"notebook kernel imports pytRIBS ({txt.replace('KERNEL ', '')})")
        except Exception as exc:
            add("FAIL", f"notebook kernel could not import pytRIBS: {type(exc).__name__}: {str(exc)[:300]}")
    man = load_manifest()
    for j in jobs:
        stale = list((lab / "results").glob(j["run_id"] + "*")) if (lab / "results").is_dir() else []
        if stale and j["run_id"] not in man["runs"]:
            add("WARN", f"results/ already has {len(stale)} file(s) named {j['run_id']}* (will be replaced)")
    kernels = subprocess.run(["ps", "-eo", "args="], capture_output=True, text=True).stdout.count("ipykernel_launcher")
    if kernels:
        add("WARN", f"{kernels} notebook kernel(s) are running (normal if VS Code is open). Do not run cells in "
                    f"{args.notebook} or anything that rebuilds data/model/ while the batch runs.")
    return out


# --------------------------------------------------------------------------------------------
# One run
# --------------------------------------------------------------------------------------------
def remove_stale(lab: Path, rid: str) -> None:
    for p in list((lab / "results").glob(rid + "*")):
        try:
            p.unlink() if p.is_file() else shutil.rmtree(p)
        except OSError:
            pass
    for p in (lab / f"{rid}.in", lab / f"simulation_{rid}.log"):
        if p.is_file():
            p.unlink()


def describe(job: dict) -> str:
    return "baseline (no change)" if job["kw"] is None else f"{job['knob']} = {job['value']:g}"


def do_run(job: dict, args, lab: Path, man: dict, ref_in: Path) -> None:
    rid, s = job["run_id"], settings_for(job)
    rec = man["runs"].setdefault(rid, {})
    log_dir, score_dir = lab / LOG_DIR, lab / SCORE_DIR
    log_dir.mkdir(exist_ok=True)
    score_dir.mkdir(exist_ok=True)
    cur = settings_record(job)

    if rec.get("status") in ("done", "unscored") and rec.get("settings") != cur and rid not in args.redo_ids:
        log(f"{rid}: the settings changed since it was run ({rec.get('settings')} -> {cur}); running it again")
        args.redo_ids.add(rid)
    if rec.get("status") == "done" and rid not in args.redo_ids:
        log(f"{rid}: already done, skipping")
        return
    if rec.get("status") == "unscored" and rid not in args.redo_ids:
        log(f"{rid}: ran earlier but was not scored; scoring again")
        rec["scores"] = show_scores(run_scorer(rid, str(lab / args.scorer), score_dir, log_dir))
        rec["status"] = "done"
        rec["finished"] = now_iso()
        save_manifest(man)
        return

    free = disk_free_gb(lab)
    if free < args.min_free_gb:
        raise Abort(f"Only {free:.1f} GB free; need {args.min_free_gb:g} GB. Free some space and run again to resume.")
    comp = competing_tribs()
    if comp:
        raise Abort("Another tRIBS is running, so the batch stops (two at once is not allowed): " + comp[0])

    rec.clear()
    rec.update(status="building", started=now_iso(), kind=job["kind"], knob=job["knob"], kw=job["kw"], value=job["value"],
               settings=cur)
    save_manifest(man)
    remove_stale(lab, rid)
    log(f"=== {rid}: {describe(job)} on W1 ({W1['hours']} h from {W1['start']}, water table {W1['gw']:g} x bedrock) ===")

    t_build = time.time()
    text, warns = build_run(job, args, lab, log_dir)
    check_in_file(lab / s["INPUT_FILE"], ref_in, job)
    grid_note = check_grids(lab, job, t_build)
    rec.update(build_min=round((time.time() - t_build) / 60, 2), warnings=warns, grid_check=grid_note, status="built")
    save_manifest(man)
    log(f"{rid}: build OK in {fmt_min(time.time() - t_build)}; .in file matches the reference except the allowed lines"
        + ("" if job["kw"] is None else f", and {job['kw']} = {job['value']:g}"))
    log(f"{rid}: {grid_note}")
    for wmsg in warns:
        log(f"{rid}: notebook said: {wmsg}")
    if args.build_only:
        log(f"{rid}: --build-only, not running tRIBS")
        return

    install_term_handler()
    log(f"{rid}: starting tRIBS ({W1['hours']} model hours, expect about {W1['hours'] * MIN_PER_1000H / 1000:.1f} min)")
    rec["status"] = "running"
    save_manifest(man)
    rc, secs = run_tribs(args.tribs, s["INPUT_FILE"], f"simulation_{rid}.log", timeout_for(args),
                         progress_file=lab / "results" / f"{rid}_Outlet.qout", hours=W1["hours"],
                         heartbeat_min=args.heartbeat_min, label=rid)
    rec["tribs_min"] = round(secs / 60, 2)
    rec["tribs_returncode"] = rc
    if rc != 0:
        raise RunError(f"tRIBS exited with code {rc}; see simulation_{rid}.log")
    outlet = check_outlet(lab / "results" / f"{rid}_Outlet.qout", W1["hours"])
    rec["outlet"] = outlet
    rec["status"] = "unscored"
    save_manifest(man)
    log(f"{rid}: tRIBS finished in {fmt_min(secs)}; outlet peak {outlet['peak_m3s']:.1f} m3/s over {outlet['rows']} rows")

    try:
        rec["scores"] = show_scores(run_scorer(rid, str(lab / args.scorer), score_dir, log_dir))
        rec["status"] = "done"
    except RunError as exc:
        rec["score_error"] = str(exc)
        log(f"{rid}: scoring failed (the run itself is fine): {exc}")
    rec["finished"] = now_iso()
    save_manifest(man)
    if rec.get("scores"):
        log(f"{rid}: {one_line(rec['scores'])}")


# --------------------------------------------------------------------------------------------
# The table and the summary CSV
# --------------------------------------------------------------------------------------------
def load_scores(lab: Path) -> dict:
    p = lab / SCORE_DIR / SCORE_CSV
    out = {}
    if p.is_file():
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                out[r["run_id"]] = r
    return out


def fm(x, nd=2, sign=False) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "-"
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def print_summary(man: dict, lab: Path, knob_keys: list) -> None:
    scores = load_scores(lab)
    base_id = man.get("baseline_id") or BASELINE_ID
    base = scores.get(base_id)
    out_rows = []
    log()
    log("ROUTING SCREEN ON W1 (flood week, Jan 18-25 2017). Everything not screened is at the baseline.")
    log("  vol = model volume / gauge volume.  peak = model peak vs gauge peak, %.  arr = first-arrival error, hours "
        "(model minus gauge; negative = early).")
    log("  'change' columns are the run minus the baseline row (volume ratio points, peak percentage points, hours).")
    log("  2nd storm = Jan 23 12:00 to Jan 25 00:00.  Proposed tolerances: volume within 10%, peak within 15%, arrival within 1 h.")
    if base is None:
        log(f"  (the baseline run {base_id} has no score yet, so the change columns are blank)")
    for k in KNOBS:
        if k["key"] not in knob_keys:
            continue
        log()
        log(f"  {k['key']:<4}{k['label']} ({k['kw']}), baseline {k['base']:g}")
        log(f"    {'value':>7} | {'vol':>6} {'peak %':>7} {'arr h':>6} | {'chg vol':>8} {'chg peak':>9} {'chg arr':>8} | "
            f"{'2nd vol':>8} {'2nd pk %':>9} {'2nd arr':>8} | note")
        for v in sorted(k["values"]):
            is_base = (v == k["base"])
            rid = base_id if is_base else run_id_for(k, v)
            sc = scores.get(rid)
            status = (man["runs"].get(rid) or {}).get("status", "not run")
            row = dict(knob=k["key"], kw=k["kw"], value=v, run_id=rid, is_baseline=is_base, status=status)
            if sc is None:
                log(f"    {v:>7g} | {'(' + ('baseline: ' if is_base else '') + status + ')':^60}")
                out_rows.append(row)
                continue
            vals = {n: fnum(sc, c) for n, c in METRICS}
            d = {n: (vals[n] - fnum(base, c) if (base is not None and not is_base) else float("nan")) for n, c in METRICS}
            flag = sc.get("flood_week__arrival_flag", "")
            note = "baseline" if is_base else ("" if flag in ("ok", "") else f"arrival: {flag}")
            log(f"    {v:>7g} | {fm(vals['vol'], 3):>6} {fm(vals['peak'], 1):>7} {fm(vals['arr']):>6} | "
                f"{fm(d['vol'], 3, True):>8} {fm(d['peak'], 1, True):>9} {fm(d['arr'], 2, True):>8} | "
                f"{fm(vals['s2vol'], 3):>8} {fm(vals['s2peak'], 1):>9} {fm(vals['s2arr']):>8} | {note}")
            for n, c in METRICS:
                row[c] = vals[n]
                row["change_" + c] = d[n]
            row["flood_week__arrival_flag"] = flag
            out_rows.append(row)
    if out_rows:
        cols = []
        for r in out_rows:
            for c in r:
                if c not in cols:
                    cols.append(c)
        with open(lab / SUMMARY_CSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for r in out_rows:
                w.writerow({c: ("" if (isinstance(r.get(c), float) and not math.isfinite(r[c])) else r.get(c, "")) for c in cols})
        log()
        log(f"  Wrote {SUMMARY_CSV}. This table does not label any knob active, inert or trade-off yet: that is the analysis step.")


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Stage 2 piece 1: the routing-knob screen on W1 (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--summary-only", action="store_true", help="print the table from the runs already scored, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished runs do not count)")
    ap.add_argument("--knobs", default=None, help="comma-separated knobs: cv,rfe,n,wc (default: all four)")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done (or 'all')")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop a tRIBS run after this many minutes (default: 10)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=1.5)
    ap.add_argument("--heartbeat-min", type=float, default=HEARTBEAT_MIN,
                    help="print a progress line every this many minutes while tRIBS runs (default 5; keeps a Codespace awake)")
    ap.add_argument("--notebook", default=None, help="notebook to run (default: Make_Sycamore_Model.ipynb, else Make_Sycamore_Model_clean.ipynb)")
    ap.add_argument("--baseline-in", default=W4_BASELINE_IN, help="reference .in used only if the W1 baseline cannot be")
    ap.add_argument("--scorer", default=SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (TRIBS_DEFAULT if os.path.isfile(TRIBS_DEFAULT)
                                                                       else shutil.which("tRIBS") or TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    return ap.parse_args(argv)


def pick(arg, allowed, label):
    if not arg:
        return list(allowed)
    want = [x.strip() for x in arg.split(",") if x.strip()]
    unknown = [x for x in want if x not in allowed]
    if unknown:
        raise SystemExit(f"Unknown {label} {unknown}. Choose from {list(allowed)}")
    return [a for a in allowed if a in want]


def main(argv=None) -> int:
    global _LOGFILE
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / LOG_DIR).mkdir(exist_ok=True)
    _LOGFILE = str(lab / LOG_DIR / "batch_log.txt")
    knob_keys = pick(args.knobs, [k["key"] for k in KNOBS], "knob")
    man = load_manifest()

    if args.summary_only:
        print_summary(man, lab, knob_keys)
        return 0

    args.notebook, nb_note = pick_notebook(lab, args.notebook)
    all_ids = [run_id_for(k, v) for k in KNOBS for v in screened_values(k)] + [FRESH_BASELINE_ID]
    if args.redo == "all":
        args.redo_ids = set(all_ids)
    else:
        args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())

    base_problems, base_notes = baseline_check(lab)
    base_ok = not base_problems

    def planned(base_usable: bool) -> list:
        jobs = [] if base_usable else [make_job(None)]
        for k in KNOBS:
            if k["key"] in knob_keys:
                jobs += [make_job(k, v) for v in screened_values(k)]
        return jobs

    preview = planned(base_ok)
    log(f"Working folder: {lab}")
    log(f"Notebook: {args.notebook}")
    log(f"Plan: knobs {knob_keys}, window W1 ({W1['hours']} h), {len(preview)} run(s) to make.")
    if base_ok:
        log(f"Reference run: existing {BASELINE_ID} ({'; '.join(base_notes)})")
    else:
        log(f"Reference run: existing {BASELINE_ID} CANNOT be used:")
        for p in base_problems:
            log(f"    {p}")
        log(f"    so a fresh baseline {FRESH_BASELINE_ID} will be made first.")
    for j in preview:
        gw_name, _ = expected_files(settings_for(j))
        log(f"  {j['run_id']:42s} {describe(j):22s} water-table file {gw_name}")
    log(f"  tRIBS time alone: about {len(preview) * W1['hours'] * MIN_PER_1000H / 1000:.0f} min. "
        f"Each run also needs a notebook build (the first build will show how long that takes).")

    results = preflight(args, lab, preview, base_ok, nb_note)
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
        take_lock()
    except Abort as exc:
        log(str(exc))
        return 1

    man["script"] = "run_screen_routing_sycamore.py"
    man["last_started"] = now_iso()
    save_manifest(man)

    install_term_handler()

    t_batch, consecutive_fail, code = time.time(), 0, 0
    try:
        base_ok = do_baseline(args, lab, man)
        ref_in = lab / (f"{BASELINE_ID}.in" if base_ok else args.baseline_in)
        todo = planned(base_ok)
        if args.limit is not None:
            def finished(j):
                rec = man["runs"].get(j["run_id"], {})
                return (rec.get("status") == "done" and rec.get("settings") == settings_record(j)
                        and j["run_id"] not in args.redo_ids)
            todo = [j for j in todo if not finished(j)][:args.limit]
            if not todo:
                log("Nothing left to run: every run in this selection is already done.")
        for j in todo:
            rid = j["run_id"]
            try:
                do_run(j, args, lab, man, ref_in)
                consecutive_fail = 0
            except RunError as exc:
                man["runs"].setdefault(rid, {}).update(status="failed", error=str(exc), finished=now_iso())
                save_manifest(man)
                log(f"{rid}: FAILED. {exc}")
                consecutive_fail += 1
                code = 1
                if consecutive_fail >= 2:
                    log("Two runs failed in a row, so the batch stops (the cause is probably the same).")
                    break
    except Abort as exc:
        log(f"BATCH STOPPED: {exc}")
        code = 1
    except KeyboardInterrupt:
        log("BATCH INTERRUPTED. Run the same command again to resume; finished runs are skipped.")
        for rid, rec in man["runs"].items():
            if rec.get("status") in ("building", "running", "checking"):
                rec["status"] = "interrupted"
        save_manifest(man)
        code = 130
    finally:
        save_manifest(man)

    log()
    log(f"Batch time: {fmt_min(time.time() - t_batch)}")
    log(f"{'run':44s} {'status':16s} {'build':>7s} {'tRIBS':>7s}")
    for rid, rec in man["runs"].items():
        log(f"{rid:44s} {rec.get('status', '?'):16s} {str(rec.get('build_min', '')):>7s} {str(rec.get('tribs_min', '')):>7s}")
    print_summary(man, lab, knob_keys)
    return code


if __name__ == "__main__":
    sys.exit(main())
