#!/usr/bin/env python3
"""
run_start_state_sycamore.py  --  Stage 1, step 1: the six start-state runs on window W3.

What it does, for each of six runs (one change at a time from the baseline):
    1. Copies the Settings cell of Make_Sycamore_Model.ipynb, changes ONLY the values for this run,
       and runs the whole notebook headlessly (no browser) with nbclient. The notebook does the
       building, exactly as when you press "Run All".
    2. Checks the build: the new .in file must match the baseline .in in every keyword except the
       ones this run is allowed to change, and the bedrock / water-table / soil files must be fresh
       and correct (this is the guard against the stale-water-table-file problem).
    3. Runs tRIBS (one at a time, never two at once) and checks that it finished.
    4. Scores the run with score_run_sycamore.py.
    5. Writes everything it did to start_state_manifest.json, so a stopped batch can be resumed.

Run it from workspaces/Sycamore/lab/ (it moves itself there). Typical use:
    python -u run_start_state_sycamore.py --dry-run          # checks only, builds nothing
    python -u run_start_state_sycamore.py --limit 1 --build-only
    python -u run_start_state_sycamore.py --limit 1          # first run, all steps
    nohup python -u run_start_state_sycamore.py > start_state_batch.out 2>&1 &   # the rest

Do NOT run Make_Sycamore_Model.ipynb (or any cell that rebuilds data/model/) while this runs:
every build rewrites lab/data/model/, and the six runs share it.
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
# Settings of the baseline run (Sycamore_octmar_ks0p2_cv6_psib0p5_nocl), with the W3 window.
BASE = {
    "GW_FRACTION": 0.95, "BEDROCK_ADD_M": 1.0,
    "KS_MULT": 0.2, "PSIB_MULT": 0.5, "F_MULT": 1,
    "KINEMVELCOEF": 6, "FLOWEXP": 0.3, "CHANNELROUGHNESS": 0.04, "CHANNELWIDTHCOEFF": 2.33,
    "OPTSNOW": 1, "OPTPERCOLATION": 0, "CHANNELCONDUCTIVITY": 5,
    "REBUILD_SOIL": False, "RUN_CHOICE": "oct_to_jan",
}
# (tag, setting changed, new value). One change at a time from the baseline.
TESTS = [
    ("gw0p70", "GW_FRACTION", 0.70),
    ("gw0p82", "GW_FRACTION", 0.82),
    ("gw1p00", "GW_FRACTION", 1.00),
    ("bed0p0", "BEDROCK_ADD_M", 0.0),
    ("bed0p5", "BEDROCK_ADD_M", 0.5),
    ("bed2p0", "BEDROCK_ADD_M", 2.0),
    # Control: the baseline itself (nothing changed) on the short W3 window. Added 2026-10-28 to check that the
    # long baseline cut at Jan 25 gives the same answer as a true W3 run. Last in the list, so the six runs
    # above are skipped on a resume and only this one runs.
    ("base", "GW_FRACTION", 0.95),
]
W3_HOURS = 2784                      # Oct 1 2016 00:00 to Jan 25 2017 00:00
RUN_PREFIX = "W3_ks0p2_cv6_psib0p5"

NOTEBOOK = "Make_Sycamore_Model.ipynb"
BASELINE_IN = "Sycamore_octmar_ks0p2_cv6_psib0p5_nocl.in"
SCORER = "score_run_sycamore.py"
TRIBS_DEFAULT = "/opt/tRIBS/build/tRIBS"
INIT_DATA = "../init_data"           # SOLUS_DIR in the notebook
MODEL_SOIL = "data/model/soil"
LOG_DIR = "start_state_logs"
SCORE_DIR = "score_out_start_state"
SCORE_CSV = "scores_start_state.csv"
MANIFEST = "start_state_manifest.json"
LOCK = "start_state.lock"

# Keywords in the .in file that a run is allowed to change, and nothing else.
ALLOWED_DIFF_ALWAYS = {"RUNTIME", "OUTFILENAME", "GWATERFILE"}
# Score columns worth printing after each run (printed only if the column exists).
SHOW_COLS = ["octdec__vol_model_Mm3", "flood_week__vol_ratio", "flood_week__peak_err_pct",
             "flood_week__arrival_err_h", "kge_volume_gate", "proc__mdgw_first_mm",
             "proc__mdgw_at_flood_week_start_mm", "proc__hsrf_share_flood_week"]

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
# Names and settings for each run
# --------------------------------------------------------------------------------------------
def settings_for(test):
    tag, key, val = test
    s = dict(BASE)
    s[key] = val
    rid = f"{RUN_PREFIX}_{tag}_nocl"
    s["INPUT_FILE"] = f"{rid}.in"
    s["RESULTS_NAME"] = rid
    return rid, s


def expected_files(s):
    """Names the notebook gives the bedrock and water-table grids for these settings."""
    bed = s["BEDROCK_ADD_M"]
    gw_name = f"InitGW_{int(round(s['GW_FRACTION'] * 100))}pct_plus{bed:g}m_mm.asc"
    bed_name = "bedrock_depth_m.asc" if bed == 0 else f"bedrock_depth_plus{bed:g}m.asc"
    return gw_name, bed_name


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


def verify_patched(src: str, values: dict) -> None:
    """Parse the patched Settings cell and confirm each patched name is assigned once, to the right value."""
    tree = ast.parse(src)
    found = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in values:
                found.setdefault(name, []).append(ast.literal_eval(node.value))
    for name, want in values.items():
        got = found.get(name, [])
        if len(got) != 1 or got[0] != want or type(got[0]) is not type(want):
            raise RunError(f"Settings cell after patching: {name} = {got}, wanted exactly one value {want!r}.")


# --------------------------------------------------------------------------------------------
# Running the notebook headlessly
# --------------------------------------------------------------------------------------------
def make_kernel_manager():
    """A kernel that uses THIS python (so it sees the same packages as this script)."""
    from jupyter_client.kernelspec import KernelSpec
    from jupyter_client.manager import KernelManager
    km = KernelManager(kernel_name="python3")
    km._kernel_spec = KernelSpec(argv=[sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
                                 display_name="start_state", language="python")
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


def build_run(rid: str, s: dict, args, lab: Path, log_dir: Path) -> tuple:
    """Patch + execute the notebook. Returns (output_text, warnings). Raises RunError on any problem."""
    import nbformat
    nb = nbformat.read(str(lab / args.notebook), as_version=4)
    for cell in nb.cells:                               # start from empty outputs
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    idx = find_settings_cell(nb)
    wanted = {k: v for k, v in s.items()}
    nb.cells[idx].source = patch_settings(nb.cells[idx].source, wanted)
    verify_patched(nb.cells[idx].source, wanted)
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
    return text, check_build_text(text, s)


def check_build_text(text: str, s: dict) -> list:
    """Lines the notebook should have printed. Returns a list of warnings; raises RunError on a problem."""
    gw_name, bed_name = expected_files(s)
    start_date, hours = "10/01/2016/00/00", W3_HOURS
    must = [
        f"Run choice: {s['RUN_CHOICE']} | start {start_date} | {hours} hours",
        "Every expected setting is in the input file.",
        f"(KS_MULT = {s['KS_MULT']:g},", f"(PSIB_MULT = {s['PSIB_MULT']:g},", f"(F_MULT = {s['F_MULT']:g},",
        f"Bedrock file in use: {INIT_DATA}/{bed_name}",
        f"GWATERFILE: {MODEL_SOIL}/{gw_name}",
        f"BEDROCKFILE: {MODEL_SOIL}/{bed_name}",
    ]
    missing = [m for m in must if m not in text]
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


def check_in_file(new_in: Path, base_in: Path, rid: str, s: dict) -> None:
    gw_name, bed_name = expected_files(s)
    new, base = read_in_file(new_in), read_in_file(base_in)
    expect = {"RUNTIME": W3_HOURS, "OUTFILENAME": f"results/{rid}", "GWATERFILE": f"{MODEL_SOIL}/{gw_name}"}
    allowed = set(ALLOWED_DIFF_ALWAYS)
    if s["BEDROCK_ADD_M"] != BASE["BEDROCK_ADD_M"]:
        allowed.add("BEDROCKFILE")
        expect["BEDROCKFILE"] = f"{MODEL_SOIL}/{bed_name}"
    problems = []
    for key in sorted(set(base) | set(new)):
        b, n = base.get(key), new.get(key)
        if key in expect:
            if n is None or not same_value(n, expect[key]):
                problems.append(f"{key}: written {n!r}, expected {expect[key]!r}")
        elif key in allowed:
            pass
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


def check_grids(lab: Path, s: dict, t_build_start: float) -> str:
    """Bedrock, water-table and soil-multiplier files in data/model/soil must be fresh and correct."""
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
    # soil multipliers: scgrid.gdf in the model folder must point at the x0.2 and x0.5 grids, f unscaled
    with open(soil / "scgrid.gdf") as fh:
        rows = {r["Variable"].strip(): r["BasePath"].strip() for r in csv.DictReader(fh)}
    for var, mult in (("KS", s["KS_MULT"]), ("PB", s["PSIB_MULT"]), ("FD", s["F_MULT"])):
        base = rows.get(var, "")
        tag = f"_x{mult:g}"
        if (mult != 1 and not base.endswith(tag)) or (mult == 1 and "_x" in os.path.basename(base)):
            raise RunError(f"scgrid.gdf row {var} is '{base}', which does not match multiplier {mult:g}.")
    return (f"water table = {s['GW_FRACTION']:g} x bedrock ({unit}), median bedrock {np.median(bed[ok]):.2f} m, "
            f"median water-table depth {np.median(gw[ok]):.0f} {unit}")


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


def run_tribs(tribs: str, infile: str, logname: str, timeout_min: float) -> tuple:
    """Run tRIBS to completion. Returns (return code, seconds). Kills the whole process group on timeout."""
    t0 = time.time()
    with open(logname, "w") as fh:
        proc = subprocess.Popen([tribs, infile], stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            rc = proc.wait(timeout=timeout_min * 60)
        except subprocess.TimeoutExpired:
            kill_group(proc)
            raise RunError(f"tRIBS did not finish within {timeout_min:g} minutes and was stopped.")
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


def run_scorer(rid: str, scorer: str, score_dir: Path, log_dir: Path) -> dict:
    cmd = [sys.executable, scorer, "score", "--run", rid, "--window", "W3", "--results-dir", "results",
           "--out-dir", str(score_dir), "--csv-name", SCORE_CSV, "--write-checksum"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    (log_dir / f"{rid}_score.txt").write_text((res.stdout or "") + (res.stderr or ""))
    if res.returncode != 0:
        raise RunError(f"The scorer exited with code {res.returncode}: {(res.stderr or res.stdout)[-400:]}")
    csv_path = score_dir / SCORE_CSV
    if csv_path.is_file():
        with open(csv_path) as fh:
            for row in csv.DictReader(fh):
                if rid in row.values():
                    return {c: row[c] for c in SHOW_COLS if c in row}
    return {}


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


def preflight(args, lab: Path, tests, fatal_only=False) -> list:
    """Checks before anything is built. Returns a list of (level, message)."""
    out = []

    def add(level, msg):
        out.append((level, msg))

    for f in (args.notebook, args.baseline_in, args.scorer):
        add("OK" if (lab / f).is_file() else "FAIL", f"file in lab/: {f}")
    for f in ("../init_data/bedrock_depth_m.asc", "../init_data/met/Master_Met.sdf", "../init_data/LandUse.asc"):
        add("OK" if (lab / f).is_file() else "FAIL", f"input: {f}")
    tribs = args.tribs
    add("OK" if (os.path.isfile(tribs) and os.access(tribs, os.X_OK)) else "FAIL", f"tRIBS program: {tribs}")
    free = disk_free_gb(lab)
    add("OK" if free >= args.min_free_gb else "FAIL", f"free disk space: {free:.1f} GB (need {args.min_free_gb:g})")
    comp = competing_tribs()
    add("FAIL" if comp else "OK", "no other tRIBS running" + (": " + comp[0] if comp else ""))
    # notebook: Settings cell present and patchable for every run
    try:
        import nbformat
        nb = nbformat.read(str(lab / args.notebook), as_version=4)
        idx = find_settings_cell(nb)
        for t in tests:
            rid, s = settings_for(t)
            verify_patched(patch_settings(nb.cells[idx].source, s), s)
        add("OK", f"Settings cell found (cell {idx}) and patchable for all {len(tests)} runs")
    except ImportError:
        add("FAIL", "nbformat is not installed: pip install nbformat nbclient ipykernel")
    except Exception as exc:
        add("FAIL", f"Settings cell check: {exc}")
    # baseline .in has the keywords we compare
    try:
        base = read_in_file(lab / args.baseline_in)
        need = {"STARTDATE": "10/01/2016/00/00", "OPTPERCOLATION": 0, "KINEMVELCOEF": 6, "FLOWEXP": 0.3,
                "CHANNELROUGHNESS": 0.04, "CHANNELWIDTHCOEFF": 2.33, "OPTSNOW": 1, "CHANNELCONDUCTIVITY": 5}
        bad = [k for k, v in need.items() if k not in base or not same_value(base[k], v)]
        add("FAIL" if bad else "OK", "baseline .in matches the W3 baseline settings" + (f" (differs: {bad})" if bad else ""))
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
            import nbformat
            from nbformat.v4 import new_code_cell, new_notebook
            nb = new_notebook(cells=[new_code_cell(
                "import sys, os\nfrom pytRIBS.classes import *\nprint('KERNEL', sys.executable, os.getcwd())")])
            execute_nb(nb, lab, 300)
            txt = outputs_text(nb).strip()
            add("OK", f"notebook kernel imports pytRIBS ({txt.replace('KERNEL ', '')})")
        except Exception as exc:
            add("FAIL", f"notebook kernel could not import pytRIBS: {type(exc).__name__}: {str(exc)[:300]}")
    # stale output with these names
    for t in tests:
        rid, _ = settings_for(t)
        stale = list((lab / "results").glob(rid + "*")) if (lab / "results").is_dir() else []
        if stale and rid not in load_manifest()["runs"]:
            add("WARN", f"results/ already has {len(stale)} file(s) named {rid}* (will be replaced)")
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


def do_run(test, args, lab: Path, man: dict) -> None:
    rid, s = settings_for(test)
    rec = man["runs"].setdefault(rid, {})
    log_dir, score_dir = lab / LOG_DIR, lab / SCORE_DIR
    log_dir.mkdir(exist_ok=True)
    score_dir.mkdir(exist_ok=True)

    if rec.get("status") == "done" and rid not in args.redo_ids:
        log(f"{rid}: already done, skipping")
        return
    if rec.get("status") == "unscored" and rid not in args.redo_ids:
        log(f"{rid}: ran earlier but was not scored; scoring again")
        rec["scores"] = run_scorer(rid, str(lab / args.scorer), score_dir, log_dir)
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
    rec.update(status="building", started=now_iso(), settings={k: s[k] for k in ("GW_FRACTION", "BEDROCK_ADD_M")},
               changed=f"{test[1]} = {test[2]}")
    save_manifest(man)
    remove_stale(lab, rid)
    log(f"=== {rid}: {test[1]} = {test[2]} ===")

    t_build = time.time()
    text, warns = build_run(rid, s, args, lab, log_dir)
    check_in_file(lab / s["INPUT_FILE"], lab / args.baseline_in, rid, s)
    grid_note = check_grids(lab, s, t_build)
    rec.update(build_min=round((time.time() - t_build) / 60, 2), warnings=warns, grid_check=grid_note,
               status="built")
    save_manifest(man)
    log(f"{rid}: build OK in {fmt_min(time.time() - t_build)}; .in file matches baseline except the allowed lines")
    log(f"{rid}: {grid_note}")
    for w in warns:
        log(f"{rid}: notebook said: {w}")
    if args.build_only:
        log(f"{rid}: --build-only, not running tRIBS")
        return

    log(f"{rid}: starting tRIBS ({W3_HOURS} model hours, expect about 22 min)")
    rec["status"] = "running"
    save_manifest(man)
    rc, secs = run_tribs(args.tribs, s["INPUT_FILE"], f"simulation_{rid}.log", args.timeout_min)
    rec["tribs_min"] = round(secs / 60, 2)
    rec["tribs_returncode"] = rc
    if rc != 0:
        raise RunError(f"tRIBS exited with code {rc}; see simulation_{rid}.log")
    outlet = check_outlet(lab / "results" / f"{rid}_Outlet.qout", W3_HOURS)
    rec["outlet"] = outlet
    rec["status"] = "unscored"
    save_manifest(man)
    log(f"{rid}: tRIBS finished in {fmt_min(secs)}; outlet peak {outlet['peak_m3s']:.1f} m3/s over {outlet['rows']} rows")

    try:
        rec["scores"] = run_scorer(rid, str(lab / args.scorer), score_dir, log_dir)
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
    ap = argparse.ArgumentParser(description="Stage 1 start-state runs on window W3 (see the top of this file).")
    ap.add_argument("--dry-run", action="store_true", help="check everything, list the plan, build nothing")
    ap.add_argument("--build-only", action="store_true", help="build and check the inputs but do not run tRIBS")
    ap.add_argument("--limit", type=int, default=None, help="do only the first N runs of the list")
    ap.add_argument("--only", default=None, help="comma-separated tags, e.g. gw0p70,bed2p0")
    ap.add_argument("--redo", default=None, help="comma-separated tags to run again even if done (or 'all')")
    ap.add_argument("--timeout-min", type=float, default=75.0, help="stop a tRIBS run after this many minutes")
    ap.add_argument("--cell-timeout", type=int, default=1800, help="seconds allowed for one notebook cell")
    ap.add_argument("--min-free-gb", type=float, default=1.5)
    ap.add_argument("--notebook", default=NOTEBOOK)
    ap.add_argument("--baseline-in", default=BASELINE_IN)
    ap.add_argument("--scorer", default=SCORER)
    ap.add_argument("--tribs", default=os.environ.get("TRIBS_BIN") or (TRIBS_DEFAULT if os.path.isfile(TRIBS_DEFAULT)
                                                                       else shutil.which("tRIBS") or TRIBS_DEFAULT))
    ap.add_argument("--skip-kernel-check", action="store_true", help="skip the pytRIBS import test in --dry-run")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    global _LOGFILE
    args = parse_args(argv)
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    (lab / LOG_DIR).mkdir(exist_ok=True)
    _LOGFILE = str(lab / LOG_DIR / "batch_log.txt")

    tests = list(TESTS)
    if args.only:
        want = [t.strip() for t in args.only.split(",")]
        unknown = [t for t in want if t not in [x[0] for x in TESTS]]
        if unknown:
            print(f"Unknown tag(s) {unknown}. Choose from {[x[0] for x in TESTS]}")
            return 2
        tests = [t for t in TESTS if t[0] in want]
    if args.limit is not None:
        tests = tests[:args.limit]
    all_ids = [settings_for(t)[0] for t in TESTS]
    if args.redo == "all":
        args.redo_ids = set(all_ids)
    else:
        args.redo_ids = {settings_for(t)[0] for t in TESTS if args.redo and t[0] in args.redo.split(",")}

    log(f"Working folder: {lab}")
    log(f"Plan: {len(tests)} run(s), window W3 (Oct 1 2016 to Jan 25 2017, {W3_HOURS} h), baseline = ks0p2 cv6 psib0p5 nocl")
    for t in tests:
        rid, s = settings_for(t)
        gw_name, bed_name = expected_files(s)
        log(f"  {rid:42s} {t[1]} = {t[2]:<5}  water-table file {gw_name}, bedrock file {bed_name}")

    results = preflight(args, lab, tests)
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

    man = load_manifest()
    man["script"] = "run_start_state_sycamore.py"
    man["last_started"] = now_iso()
    save_manifest(man)

    def on_term(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, on_term)

    t_batch, consecutive_fail, code = time.time(), 0, 0
    try:
        for t in tests:
            rid = settings_for(t)[0]
            try:
                do_run(t, args, lab, man)
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
            if rec.get("status") in ("building", "running"):
                rec["status"] = "interrupted"
        save_manifest(man)
        code = 130
    finally:
        save_manifest(man)

    log()
    log(f"Batch time: {fmt_min(time.time() - t_batch)}")
    log(f"{'run':44s} {'status':12s} {'build':>8s} {'tRIBS':>8s}")
    for t in TESTS:
        rid = settings_for(t)[0]
        rec = man["runs"].get(rid)
        if rec:
            log(f"{rid:44s} {rec.get('status', '?'):12s} {str(rec.get('build_min', '')):>8s} {str(rec.get('tribs_min', '')):>8s}")
    return code


if __name__ == "__main__":
    sys.exit(main())
