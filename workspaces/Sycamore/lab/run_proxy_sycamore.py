#!/usr/bin/env python3
"""
run_proxy_sycamore.py  --  Stage 1, step 2: the proxy check.

Question it answers: do the short test windows W1 (Jan 18-25) and W2 (Jan 18 - Mar 12) rank different soil
settings the same way the long window W4 (Oct 1 - Apr 1) does? If they do, the short windows are a valid, cheap
screen for calibration. (The ranking itself is done afterwards by analyze_proxy_sycamore.py; this script only
builds, runs and scores.)

Six settings (everything else as the Oct 21 baseline: channel loss off, cv 6, flowexp 0.3, roughness 0.04,
width 2.33, f x1, bedrock +1 m, snow on):

    1  Ks x0.2  PsiB x0.5   (the baseline)
    2  Ks x0.2  PsiB x1.0
    3  Ks x0.1  PsiB x0.5
    4  Ks x0.3  PsiB x0.5
    5  Ks x0.3  PsiB x1.0
    6  Ks x0.1  PsiB x1.0

Each setting is run on W1, W2 and W4. Two W4 runs already exist from earlier work (settings 1 and 2) and are
REUSED if they check out: they are not run again, only scored with the window W4. Setting 1 must be
byte-identical to the file whose checksum was recorded; setting 2 must reproduce the numbers in the handoff.
If a reused run fails its check, the script runs that W4 again under a new name.

For every run to be made, the script:
    1. Copies the Settings cell of Make_Sycamore_Model.ipynb, changes ONLY the values for this run (the soil
       multipliers, the starting water table, the window start and length, the two file names), and runs the
       whole notebook headlessly (no browser) with nbclient.
    2. Checks the build: the new .in file must match the baseline .in in every keyword except the ones this run
       may change, and the bedrock / water-table / soil files must be fresh and correct.
    3. Runs tRIBS (one at a time, never two at once) and checks that it finished.
    4. Scores the run with score_run_sycamore.py.
    5. Writes everything to proxy_manifest.json, so a stopped batch can be resumed by running the same command.

Order: reused W4 runs are checked first, then W1 (six short runs), then W2 (six medium runs), then the W4 runs
that still have to be made. So the short windows can be looked at long before the last W4 run finishes.

Run it from workspaces/Sycamore/lab/. Typical use:
    python -u run_proxy_sycamore.py --reuse-check               # only score and check the two old W4 runs
    python -u run_proxy_sycamore.py --dry-run                   # checks only, builds nothing
    python -u run_proxy_sycamore.py --limit 1 --build-only      # build the first run and check it, no tRIBS
    python -u run_proxy_sycamore.py --limit 1                   # first run, all steps
    python -u run_proxy_sycamore.py                             # everything; run it again to resume

Do NOT run Make_Sycamore_Model.ipynb (or any cell that rebuilds data/model/) while this runs: every build
rewrites lab/data/model/, and all the runs share it.
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
# What the batch does
# --------------------------------------------------------------------------------------------
# Settings shared by every run (the Oct 21 baseline). GW_FRACTION and the window differ by window, see WINDOWS.
BASE = {
    "BEDROCK_ADD_M": 1.0, "F_MULT": 1,
    "KINEMVELCOEF": 6, "FLOWEXP": 0.3, "CHANNELROUGHNESS": 0.04, "CHANNELWIDTHCOEFF": 2.33,
    "OPTSNOW": 1, "OPTPERCOLATION": 0, "CHANNELCONDUCTIVITY": 5,
    "REBUILD_SOIL": False,
}
# (number, Ks multiplier, PsiB multiplier)
SETTINGS = [
    (1, 0.2, 0.5),
    (2, 0.2, 1.0),
    (3, 0.1, 0.5),
    (4, 0.3, 0.5),
    (5, 0.3, 1.0),
    (6, 0.1, 1.0),
]
COLD_GW = 0.82      # W1 and W2 start cold on Jan 18: water table at 82% of bedrock depth (the design's value)
WARM_GW = 0.95      # W4 starts Oct 1: 95% of bedrock depth (the start-state rule from step 1)
WINDOWS = {
    # start text for tRIBS, run length in hours, label printed by the notebook, starting water table
    "W1": dict(start="01/18/2017/00/00", hours=168, label="jan_storm_week", gw=COLD_GW),
    "W2": dict(start="01/18/2017/00/00", hours=1272, label="jan_to_mar", gw=COLD_GW),     # Jan 18 to Mar 12, 53 days
    "W4": dict(start="10/01/2016/00/00", hours=4368, label="oct_to_mar", gw=WARM_GW),
}
WINDOW_ORDER = ["W1", "W2", "W4"]
MIN_PER_1000H = 7.8          # measured tRIBS speed on this machine: minutes per 1,000 model hours

# W4 runs that already exist. A reused run is only scored, never run again, and only if it passes its check.
REUSE = {
    1: dict(run_id="Sycamore_octmar_ks0p2_cv6_psib0p5_nocl", md5="aab95272a5473e66abf65d80ba61b1f9", numbers={}),
    2: dict(run_id="Sycamore_octmar_ks0p2_cv6_nocl", md5=None,
            # (value in the handoff, how far the new score may differ). All are columns of the scorer's CSV.
            numbers={"oct_mar__vol_model_Mm3_raw": (11.49, 0.03),     # Oct-Mar model volume, M m3
                     "flood_week__vol_ratio": (0.904, 0.004),         # flood-week volume, -9.6%
                     "flood_week__peak_err_pct": (-25.8, 0.4)}),      # flood-week peak, -25.8%
}

NOTEBOOK = "Make_Sycamore_Model.ipynb"
BASELINE_IN = "Sycamore_octmar_ks0p2_cv6_psib0p5_nocl.in"       # the W4 baseline; every new .in is compared with it
SCORER = "score_run_sycamore.py"
TRIBS_DEFAULT = "/opt/tRIBS/build/tRIBS"
INIT_DATA = "../init_data"           # SOLUS_DIR in the notebook
MODEL_SOIL = "data/model/soil"
LOG_DIR = "proxy_logs"
SCORE_DIR = "score_out_proxy"
SCORE_CSV = "scores_proxy.csv"
MANIFEST = "proxy_manifest.json"
LOCK = "proxy.lock"

# Keywords in the .in file that a run is allowed to change. Each is also checked against its expected value.
ALLOWED_DIFF_ALWAYS = {"RUNTIME", "STARTDATE", "OUTFILENAME", "GWATERFILE"}
# Score columns worth printing after each run (printed only if the column exists and has a value).
SHOW_COLS = ["flood_week__vol_ratio", "flood_week__peak_err_pct", "flood_week__arrival_err_h",
             "jan_flood_tail__vol_ratio", "feb18_19_tail__vol_ratio", "feb27_28_tail__vol_ratio",
             "oct_mar__vol_ratio", "span__vol_ratio"]

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
def tag(x: float) -> str:
    """0.1 -> '0p1', 1.0 -> '1p0'."""
    return f"{x:.1f}".replace(".", "p")


def run_id_for(window: str, ks: float, psib: float) -> str:
    return f"{window}_ks{tag(ks)}_cv6_psib{tag(psib)}_nocl"


def make_job(number: int, window: str, kind: str = "run", run_id: str | None = None) -> dict:
    _, ks, psib = next(s for s in SETTINGS if s[0] == number)
    return dict(setting=number, ks=ks, psib=psib, window=window, kind=kind,
                run_id=run_id or run_id_for(window, ks, psib))


def settings_for(job: dict) -> dict:
    """Every NAME = value line of the Settings cell that this run changes, plus the window."""
    w = WINDOWS[job["window"]]
    s = dict(BASE)
    s.update(GW_FRACTION=w["gw"], KS_MULT=job["ks"], PSIB_MULT=job["psib"], RUN_CHOICE=w["label"])
    s["INPUT_FILE"] = f"{job['run_id']}.in"
    s["RESULTS_NAME"] = job["run_id"]
    return s


def settings_record(job: dict) -> dict:
    """What the manifest remembers about the settings, to notice a changed setting on resume."""
    w = WINDOWS[job["window"]]
    return dict(GW_FRACTION=w["gw"], BEDROCK_ADD_M=BASE["BEDROCK_ADD_M"], KS_MULT=job["ks"], PSIB_MULT=job["psib"],
                START_DATE=w["start"], RUN_HOURS=w["hours"])


def expected_files(s):
    """Names the notebook gives the bedrock and water-table grids for these settings."""
    bed = s["BEDROCK_ADD_M"]
    gw_name = f"InitGW_{int(round(s['GW_FRACTION'] * 100))}pct_plus{bed:g}m_mm.asc"
    bed_name = "bedrock_depth_m.asc" if bed == 0 else f"bedrock_depth_plus{bed:g}m.asc"
    return gw_name, bed_name


def timeout_for(job: dict, args) -> float:
    if args.timeout_min:
        return args.timeout_min
    return max(10.0, 3.0 * WINDOWS[job["window"]]["hours"] * MIN_PER_1000H / 1000.0)


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
    """Replace the line 'START_DATE, RUN_HOURS = RUNS[RUN_CHOICE]' by the window's own start and length.
    (W2 is not in the notebook's RUNS list, so every window is written out explicitly.)"""
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
    w = WINDOWS[job["window"]]
    out = patch_window(patch_settings(src, s), w["start"], w["hours"])
    verify_patched(out, s, (w["start"], w["hours"]))
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
                                 display_name="proxy", language="python")
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
    w = WINDOWS[job["window"]]
    gw_name, bed_name = expected_files(s)
    must = [
        f"Run choice: {s['RUN_CHOICE']} | start {w['start']} | {w['hours']} hours",
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
    warns = [ln.strip() for ln in text.splitlines() if re.match(r"^CHECK\s", ln)]
    return warns


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


def check_in_file(new_in: Path, base_in: Path, job: dict) -> None:
    s = settings_for(job)
    w = WINDOWS[job["window"]]
    gw_name, _ = expected_files(s)
    new, base = read_in_file(new_in), read_in_file(base_in)
    expect = {"RUNTIME": w["hours"], "STARTDATE": w["start"], "OUTFILENAME": f"results/{job['run_id']}",
              "GWATERFILE": f"{MODEL_SOIL}/{gw_name}"}
    problems = []
    for key in sorted(set(base) | set(new)):
        b, n = base.get(key), new.get(key)
        if key in expect:
            if n is None or not same_value(n, expect[key]):
                problems.append(f"{key}: written {n!r}, expected {expect[key]!r}")
        elif b is None or n is None or not same_value(b, n):
            problems.append(f"{key}: baseline {b!r} but this run {n!r}  (this keyword may not change)")
    if problems:
        raise RunError("The .in file differs from the baseline in more than it should:\n  " + "\n  ".join(problems))


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
    # bedrock: the SOLUS original plus BEDROCK_ADD_M everywhere it has data
    orig_h, orig = read_asc(lab / INIT_DATA / "bedrock_depth_m.asc")
    ov = (orig != orig_h.get("nodata_value", -9999)) & (bed != nb_)
    diff = (bed - orig)[ov]
    if diff.size == 0 or np.max(np.abs(diff - s["BEDROCK_ADD_M"])) > 0.011:
        raise RunError(f"Bedrock grid is not the SOLUS grid + {s['BEDROCK_ADD_M']:g} m "
                       f"(difference ranges {diff.min():.3f} to {diff.max():.3f} m).")
    # water table: gw depth / bedrock depth should equal GW_FRACTION (grid is in mm, bedrock in m)
    med = float(np.median(gw[ok] / bed[ok]))
    if abs(med - s["GW_FRACTION"] * 1000) <= 5:
        unit = "mm"
    elif abs(med - s["GW_FRACTION"]) <= 0.005:
        unit = "m"
    else:
        raise RunError(f"Water-table grid / bedrock grid has median {med:.4f}; expected "
                       f"{s['GW_FRACTION'] * 1000:g} (mm per m) for GW_FRACTION = {s['GW_FRACTION']}.")
    # soil multipliers: scgrid.gdf in the model folder must point at the right multiplied grids
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


HEARTBEAT_MIN = 5.0          # how often to print a progress line while tRIBS runs (see --heartbeat-min)
ROWS_PER_MODEL_HOUR = 16     # the outlet file gets 16 rows per model hour (W2: 20352 rows for 1272 h)


def count_rows(path: Path) -> int:
    """Number of lines in a file, or 0 if it does not exist yet."""
    try:
        with open(path, "rb") as fh:
            return sum(chunk.count(b"\n") for chunk in iter(lambda: fh.read(1 << 20), b""))
    except OSError:
        return 0


def run_tribs(tribs: str, infile: str, logname: str, timeout_min: float,
              progress_file: Path | None = None, hours: int | None = None,
              heartbeat_min: float = HEARTBEAT_MIN, label: str = "") -> tuple:
    """Run tRIBS to completion. Returns (return code, seconds). Kills the whole process group on timeout.

    tRIBS prints nothing while it runs. Every `heartbeat_min` minutes this prints a progress line, which also counts as
    terminal activity, so a Codespace (30 minute idle timeout) does not shut down in the middle of a long run.
    """
    t0 = time.time()
    deadline = t0 + timeout_min * 60
    beat = max(heartbeat_min, 0.05) * 60
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
                    msg = f"{label}: tRIBS still running, {(time.time() - t0) / 60:.0f} min so far"
                    if progress_file is not None:
                        n = count_rows(progress_file)
                        msg += f"; outlet file has {n} rows" + (f" of about {expected} ({100 * n / expected:.0f}%)" if expected else "")
                    log(msg)
        except BaseException:
            kill_group(proc)
            raise
    return rc, time.time() - t0


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


def run_scorer(rid: str, window: str, scorer: str, score_dir: Path, log_dir: Path) -> dict:
    """Score one run; returns its whole row of the scorer's CSV (all values as text)."""
    cmd = [sys.executable, scorer, "score", "--run", rid, "--window", window, "--results-dir", "results",
           "--out-dir", str(score_dir), "--csv-name", SCORE_CSV, "--write-checksum"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    (log_dir / f"{rid}_score.txt").write_text((res.stdout or "") + (res.stderr or ""))
    if res.returncode != 0:
        raise RunError(f"The scorer exited with code {res.returncode}: {(res.stderr or res.stdout)[-400:]}")
    csv_path = score_dir / SCORE_CSV
    if csv_path.is_file():
        with open(csv_path) as fh:
            for row in csv.DictReader(fh):
                if row.get("run_id") == rid:
                    return dict(row)
    raise RunError(f"The scorer ran but {rid} is not in {csv_path}.")


def show_scores(row: dict) -> dict:
    out = {}
    for c in SHOW_COLS:
        v = row.get(c, "")
        if v not in ("", None):
            try:
                out[c] = round(float(v), 3)
            except ValueError:
                out[c] = v
    return out


# --------------------------------------------------------------------------------------------
# Manifest, lock, preflight
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


def preflight(args, lab: Path, jobs: list, reuse_jobs: list, scoring_only: bool = False) -> list:
    """Checks before anything is built. Returns a list of (level, message)."""
    out = []

    def add(level, msg):
        out.append((level, msg))

    files = [args.scorer] if scoring_only else [args.notebook, args.baseline_in, args.scorer]
    for f in files:
        add("OK" if (lab / f).is_file() else "FAIL", f"file in lab/: {f}")
    obs = lab / INIT_DATA / "obs" / "sycamore_09510200_obs_flow_WY2017.csv"
    add("OK" if obs.is_file() else "FAIL", f"gauge file: {INIT_DATA}/obs/{obs.name}")
    for j in reuse_jobs:
        p = lab / "results" / f"{j['run_id']}_Outlet.qout"
        add("OK" if p.is_file() else "WARN", f"old W4 run to reuse (setting {j['setting']}): {p.name}"
            + ("" if p.is_file() else "  -> missing, so this W4 run will be made again"))
    if scoring_only:
        return out
    for f in ("../init_data/bedrock_depth_m.asc", "../init_data/met/Master_Met.sdf", "../init_data/LandUse.asc"):
        add("OK" if (lab / f).is_file() else "FAIL", f"input: {f}")
    tribs = args.tribs
    add("OK" if (os.path.isfile(tribs) and os.access(tribs, os.X_OK)) else "FAIL", f"tRIBS program: {tribs}")
    free = disk_free_gb(lab)
    add("OK" if free >= args.min_free_gb else "FAIL", f"free disk space: {free:.1f} GB (need {args.min_free_gb:g})")
    comp = competing_tribs()
    add("FAIL" if comp else "OK", "no other tRIBS running" + (": " + comp[0] if comp else ""))
    # notebook: Settings cell present and patchable for every run, for every window
    try:
        import nbformat
        nb = nbformat.read(str(lab / args.notebook), as_version=4)
        idx = find_settings_cell(nb)
        every = [make_job(n, w) for n, _, _ in SETTINGS for w in WINDOW_ORDER]
        for j in every:
            patched_source(nb.cells[idx].source, j)
        add("OK", f"Settings cell found (cell {idx}) and patchable for all {len(every)} setting/window combinations")
    except ImportError:
        add("FAIL", "nbformat is not installed: pip install nbformat nbclient ipykernel")
    except Exception as exc:
        add("FAIL", f"Settings cell check: {exc}")
    # baseline .in has the keywords we compare
    try:
        base = read_in_file(lab / args.baseline_in)
        need = {"STARTDATE": WINDOWS["W4"]["start"], "RUNTIME": WINDOWS["W4"]["hours"], "OPTPERCOLATION": 0,
                "KINEMVELCOEF": 6, "FLOWEXP": 0.3, "CHANNELROUGHNESS": 0.04, "CHANNELWIDTHCOEFF": 2.33, "OPTSNOW": 1,
                "CHANNELCONDUCTIVITY": 5}
        bad = [k for k, v in need.items() if k not in base or not same_value(base[k], v)]
        add("FAIL" if bad else "OK", "baseline .in matches the W4 baseline settings" + (f" (differs: {bad})" if bad else ""))
    except Exception as exc:
        add("FAIL", f"baseline .in: {exc}")
    # scorer flags exist
    try:
        res = subprocess.run([sys.executable, str(lab / args.scorer), "score", "--help"], capture_output=True, text=True)
        flags = ["--run", "--window", "--results-dir", "--out-dir", "--csv-name", "--write-checksum"]
        miss = [f for f in flags if f not in res.stdout]
        add("FAIL" if (res.returncode != 0 or miss) else "OK",
            "scorer accepts " + " ".join(flags) + (f" (missing: {miss})" if miss else ""))
    except Exception as exc:
        add("FAIL", f"scorer check: {exc}")
    # kernel can import pytRIBS
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
    # stale output with these names
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
# Reusing an old W4 run
# --------------------------------------------------------------------------------------------
def do_reuse(job: dict, args, lab: Path, man: dict) -> bool:
    """Score an existing W4 run and check it is the run we think it is. Returns True if it can be reused."""
    n, rid = job["setting"], job["run_id"]
    rec = man["runs"].setdefault(rid, {})
    log_dir, score_dir = lab / LOG_DIR, lab / SCORE_DIR
    log_dir.mkdir(exist_ok=True)
    score_dir.mkdir(exist_ok=True)
    spec = REUSE[n]
    rec.clear()
    rec.update(setting=n, window="W4", kind="reuse", ks=job["ks"], psib=job["psib"], status="checking", started=now_iso())
    problems = []
    try:
        outlet = check_outlet(lab / "results" / f"{rid}_Outlet.qout", WINDOWS["W4"]["hours"])
        rec["outlet"] = outlet
        if spec["md5"] and outlet["md5"] != spec["md5"]:
            problems.append(f"the outlet file's md5 is {outlet['md5']} but the recorded one is {spec['md5']} "
                            f"(the file has changed since it was checked)")
        row = run_scorer(rid, "W4", str(lab / args.scorer), score_dir, log_dir)
        rec["scores"] = show_scores(row)
        for col, (want, tol) in spec["numbers"].items():
            try:
                got = float(row.get(col, ""))
            except ValueError:
                problems.append(f"{col} is blank in the score")
                continue
            if abs(got - want) > tol:
                problems.append(f"{col} = {got:.3f}, the handoff says {want} (allowed difference {tol})")
            else:
                log(f"{rid}: {col} = {got:.3f} matches the handoff value {want}")
    except RunError as exc:
        problems.append(str(exc))
    if problems:
        rec.update(status="reuse_rejected", problems=problems, finished=now_iso())
        save_manifest(man)
        log(f"{rid}: CANNOT be reused:")
        for p in problems:
            log(f"    {p}")
        return False
    rec.update(status="reused", finished=now_iso(), settings=settings_record(job))
    save_manifest(man)
    log(f"{rid}: reused as the W4 run of setting {n} (Ks x{job['ks']:g}, PsiB x{job['psib']:g}). "
        + ", ".join(f"{k}={v}" for k, v in rec["scores"].items() if k in ("flood_week__vol_ratio", "flood_week__peak_err_pct", "oct_mar__vol_ratio")))
    return True


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


def do_run(job: dict, args, lab: Path, man: dict) -> None:
    rid, s, w = job["run_id"], settings_for(job), WINDOWS[job["window"]]
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
        rec["scores"] = show_scores(run_scorer(rid, job["window"], str(lab / args.scorer), score_dir, log_dir))
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
    rec.update(status="building", started=now_iso(), setting=job["setting"], window=job["window"], kind="run",
               ks=job["ks"], psib=job["psib"], settings=cur)
    save_manifest(man)
    remove_stale(lab, rid)
    log(f"=== {rid}: setting {job['setting']} (Ks x{job['ks']:g}, PsiB x{job['psib']:g}) on {job['window']} "
        f"({w['hours']} h from {w['start']}, water table {w['gw']:g} x bedrock) ===")

    t_build = time.time()
    text, warns = build_run(job, args, lab, log_dir)
    check_in_file(lab / s["INPUT_FILE"], lab / args.baseline_in, job)
    grid_note = check_grids(lab, job, t_build)
    rec.update(build_min=round((time.time() - t_build) / 60, 2), warnings=warns, grid_check=grid_note, status="built")
    save_manifest(man)
    log(f"{rid}: build OK in {fmt_min(time.time() - t_build)}; .in file matches the baseline except the allowed lines")
    log(f"{rid}: {grid_note}")
    for wmsg in warns:
        log(f"{rid}: notebook said: {wmsg}")
    if args.build_only:
        log(f"{rid}: --build-only, not running tRIBS")
        return

    log(f"{rid}: starting tRIBS ({w['hours']} model hours, expect about {w['hours'] * MIN_PER_1000H / 1000:.0f} min)")
    rec["status"] = "running"
    save_manifest(man)
    rc, secs = run_tribs(args.tribs, s["INPUT_FILE"], f"simulation_{rid}.log", timeout_for(job, args),
                         progress_file=lab / "results" / f"{rid}_Outlet.qout", hours=w["hours"],
                         heartbeat_min=args.heartbeat_min, label=rid)
    rec["tribs_min"] = round(secs / 60, 2)
    rec["tribs_returncode"] = rc
    if rc != 0:
        raise RunError(f"tRIBS exited with code {rc}; see simulation_{rid}.log")
    outlet = check_outlet(lab / "results" / f"{rid}_Outlet.qout", w["hours"])
    rec["outlet"] = outlet
    rec["status"] = "unscored"
    save_manifest(man)
    log(f"{rid}: tRIBS finished in {fmt_min(secs)}; outlet peak {outlet['peak_m3s']:.1f} m3/s over {outlet['rows']} rows")

    try:
        rec["scores"] = show_scores(run_scorer(rid, job["window"], str(lab / args.scorer), score_dir, log_dir))
        rec["status"] = "done"
    except RunError as exc:
        rec["score_error"] = str(exc)
        log(f"{rid}: scoring failed (the run itself is fine): {exc}")
    rec["finished"] = now_iso()
    save_manifest(man)
    if rec.get("scores"):
        log(f"{rid}: " + ", ".join(f"{k}={v}" for k, v in rec["scores"].items()))


# --------------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Stage 1 step 2: the proxy check (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--reuse-check", action="store_true", help="only score and check the old W4 runs, then stop")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS")
    ap.add_argument("--limit", type=int, default=None, help="make only the first N runs that still need work (finished runs and reuse checks do not count)")
    ap.add_argument("--windows", default=None, help="comma-separated windows to make, e.g. W1 or W1,W2 (default: W1,W2,W4)")
    ap.add_argument("--settings", default=None, help="comma-separated setting numbers, e.g. 3,4 (default: all six)")
    ap.add_argument("--redo", default=None, help="comma-separated run ids to run again even if done (or 'all')")
    ap.add_argument("--timeout-min", type=float, default=None, help="stop a tRIBS run after this many minutes (default: 3 x the expected time, at least 10)")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=1.5)
    ap.add_argument("--heartbeat-min", type=float, default=HEARTBEAT_MIN,
                    help="print a progress line every this many minutes while tRIBS runs (default 5; keeps a Codespace awake)")
    ap.add_argument("--notebook", default=NOTEBOOK)
    ap.add_argument("--baseline-in", default=BASELINE_IN)
    ap.add_argument("--scorer", default=SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (TRIBS_DEFAULT if os.path.isfile(TRIBS_DEFAULT)
                                                                       else shutil.which("tRIBS") or TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    return ap.parse_args(argv)


def is_complete(man: dict) -> bool:
    """True when every setting has a scored result on every window (old W4 runs count if they were reused)."""
    for n, _, _ in SETTINGS:
        for w in WINDOW_ORDER:
            ok = man["runs"].get(make_job(n, w)["run_id"], {}).get("status") == "done"
            if w == "W4" and n in REUSE:
                ok = ok or man["runs"].get(REUSE[n]["run_id"], {}).get("status") == "reused"
            if not ok:
                return False
    return True


def pick(arg, allowed, label):
    if not arg:
        return list(allowed)
    want = [x.strip() for x in arg.split(",") if x.strip()]
    unknown = [x for x in want if x not in [str(a) for a in allowed]]
    if unknown:
        raise SystemExit(f"Unknown {label} {unknown}. Choose from {[str(a) for a in allowed]}")
    return [a for a in allowed if str(a) in want]


def main(argv=None) -> int:
    global _LOGFILE
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / LOG_DIR).mkdir(exist_ok=True)
    _LOGFILE = str(lab / LOG_DIR / "batch_log.txt")

    numbers = [int(x) for x in pick(args.settings, [s[0] for s in SETTINGS], "setting number")]
    windows = pick(args.windows, WINDOW_ORDER, "window")
    man = load_manifest()

    reuse_jobs = [make_job(n, "W4", "reuse", REUSE[n]["run_id"]) for n in numbers if n in REUSE] if "W4" in windows else []
    reuse_ok = {j["setting"]: man["runs"].get(j["run_id"], {}).get("status") == "reused" for j in reuse_jobs}

    def planned_runs():
        """Run jobs in order: W1, W2, then the W4 runs for settings that have no usable old run."""
        out = []
        for w in windows:
            for n in numbers:
                if w == "W4" and reuse_ok.get(n):
                    continue
                out.append(make_job(n, w))
        return out

    all_ids = [make_job(n, w)["run_id"] for n in [s[0] for s in SETTINGS] for w in WINDOW_ORDER]
    if args.redo == "all":
        args.redo_ids = set(all_ids)
    else:
        args.redo_ids = set(x.strip() for x in (args.redo or "").split(",") if x.strip())

    preview = planned_runs()
    log(f"Working folder: {lab}")
    log(f"Plan: settings {numbers}, windows {windows}. Old W4 runs to check for reuse: "
        + (", ".join(f"setting {j['setting']}" for j in reuse_jobs) or "none"))
    for j in ([] if args.reuse_check else preview):
        w = WINDOWS[j["window"]]
        gw_name, bed_name = expected_files(settings_for(j))
        log(f"  {j['run_id']:36s} {j['window']} {w['hours']:>5d} h  water-table file {gw_name}")
    # Same starting list both ways, so the two figures do not depend on whether the reuse check has already run.
    reuse_set = {j["setting"] for j in reuse_jobs}
    needed = [j for j in preview if not (j["window"] == "W4" and j["setting"] in reuse_set)]
    est_reuse = sum(WINDOWS[j["window"]]["hours"] for j in needed) * MIN_PER_1000H / 1000 / 60
    est = est_reuse + len(reuse_set) * WINDOWS["W4"]["hours"] * MIN_PER_1000H / 1000 / 60
    if not args.reuse_check:
        log(f"  tRIBS time alone: about {est_reuse:.1f} h if both old W4 runs are reused, {est:.1f} h if neither is. "
            f"Each run also needs a notebook build.")

    results = preflight(args, lab, preview, reuse_jobs, scoring_only=args.reuse_check)
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

    man["script"] = "run_proxy_sycamore.py"
    man["last_started"] = now_iso()
    save_manifest(man)

    def on_term(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, on_term)

    t_batch, consecutive_fail, code = time.time(), 0, 0
    try:
        for j in reuse_jobs:
            reuse_ok[j["setting"]] = do_reuse(j, args, lab, man)
        if args.reuse_check:
            log()
            log("REUSE CHECK finished. " + ("Both old W4 runs can be reused." if all(reuse_ok.values()) and reuse_ok
                                            else "At least one old W4 run cannot be reused; it will be made again in the full batch."))
            return 0 if all(reuse_ok.values()) else 1
        todo = planned_runs()
        if args.limit is not None:
            # --limit counts runs that still need work, so runs already finished must not use up the allowance
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
                do_run(j, args, lab, man)
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
    log(f"{'run':40s} {'status':16s} {'build':>7s} {'tRIBS':>7s}  flood-week volume ratio / peak error %")
    shown = [j["run_id"] for j in reuse_jobs] + [make_job(n, w)["run_id"] for w in windows for n in numbers]
    for rid in dict.fromkeys(shown):
        rec = man["runs"].get(rid)
        if rec:
            sc = rec.get("scores") or {}
            fw = f"{sc.get('flood_week__vol_ratio', '')} / {sc.get('flood_week__peak_err_pct', '')}" if sc else ""
            log(f"{rid:40s} {rec.get('status', '?'):16s} {str(rec.get('build_min', '')):>7s} {str(rec.get('tribs_min', '')):>7s}  {fw}")
    if code == 0 and is_complete(man):
        log("All 18 setting/window results are scored. Next: python analyze_proxy_sycamore.py")
    return code


if __name__ == "__main__":
    sys.exit(main())
