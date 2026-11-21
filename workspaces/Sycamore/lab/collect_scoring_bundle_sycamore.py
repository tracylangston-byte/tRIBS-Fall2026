#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collect_scoring_bundle_sycamore.py   (v1.0, 2026-11-21)
=======================================================

WHAT IT DOES
  Copies the files that Claude needs to score the Lower Sycamore runs by window into ONE NEW FOLDER, packed as one or more zip files
  you can download from Codespace and attach to the chat. It READS your files and WRITES ONLY inside that new folder. It runs no model,
  changes nothing in results/, and never overwrites anything (the new folder must not exist yet).

WHAT IT COLLECTS (for every run it finds)
  results/<run>_Outlet.qout            the outlet flow (needed)
  results/<run>*.mrf                   basin-mean rain and water table by hour (needed for the storm windows)
  results/<run>*.rft                   runoff made, by type (needed for the mechanisms breakdown)
  <run>.in                             the input file (small; lets Claude check the start date and run length)
  ../init_data/obs/sycamore_09510200_obs_flow*.csv     the USGS gauge files
  The runs are: the ten Sycamore_WY2017sep_* runs, the Sycamore_WY2020_* runs, and the two older WY2017 runs that started on Oct 1
  (W5_d012_k10_f0p15_as300_ksx1p03_nocl and W5_d012_k10_base_nocl). Runs with "kchan" in the name are left out unless you add
  --include-kchan. Anything in a results/_partial_* folder is skipped.
  It does NOT copy the .pixel files or any other big file.

HOW TO RUN (from workspaces/Sycamore/lab/, in the terminal)
  python collect_scoring_bundle_sycamore.py              # collect everything and write the zip file(s)
  python collect_scoring_bundle_sycamore.py --dry-run    # only list what it would collect; writes nothing

  Options (you will probably not need them):
    --results results             folder with the result files (default results, searched below too)
    --obs ../init_data/obs        folder with the gauge files
    --out-dir NAME                the new folder (default score_bundle_<date>_<time>, in the current folder)
    --max-mb 20                   split the zips so that none is larger than about this many megabytes
    --include-kchan               also collect runs with "kchan" in their name
    --extra RUN_NAME              also collect this run (name without _Outlet.qout); can be repeated

WHAT YOU DO AFTERWARD
  The script prints the zip file names. In the left-hand file list of Codespace, open the new folder, right-click each .zip, choose
  Download, then attach the downloaded files to the chat. The MANIFEST.txt in the folder (also inside the first zip) lists every file,
  its size and its checksum, so a damaged upload can be spotted.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import re
import sys
import zipfile
import zlib
from pathlib import Path

VERSION = "1.0"
PREFIX_2017SEP = "Sycamore_WY2017sep_"
PREFIX_2020 = "Sycamore_WY2020_"
EXPECT_2017SEP = ["Retuned", "Original", "S1_Ks1p0", "S5_Ks1p0_Bed2p0", "S6_Ks0p5_Bed2p0", "S7_Ks1p0_Bed2p0_As500",
                  "S8_Ks1p0_Bed2p0_f0p02", "S9_Ks1p0_Bed0p5", "S10_Ks1p0_GW0p5", "S11_Ks1p0_PsiB1p0"]
OLD_RUNS = ["W5_d012_k10_f0p15_as300_ksx1p03_nocl", "W5_d012_k10_base_nocl"]
GAUGE_GLOB = "sycamore_09510200_obs_flow*.csv"


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str) -> None:
    sys.exit("\nSTOPPED: " + msg + "\n")


def md5_of(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()


def describe(path: Path) -> dict:
    """Size, checksum, line count and first/last lines of a text file (read once, read-only)."""
    data = path.read_bytes()
    lines = data.splitlines()
    d = {"path": path, "size": len(data), "md5": md5_of(data), "lines": len(lines), "data": data}
    d["first"] = lines[0].decode("utf-8", "replace")[:160] if lines else ""
    d["second"] = lines[1].decode("utf-8", "replace")[:160] if len(lines) > 1 else ""
    d["last"] = lines[-1].decode("utf-8", "replace")[:160] if lines else ""
    return d


def companions(outlet: Path, stem: str, lab: Path):
    """The .mrf, .rft and .in files that go with an outlet file. The result files carry the run hours after the run name."""
    pat = re.compile("^" + re.escape(stem) + r"(\d.*)?\.(mrf|rft)$")
    found = {"mrf": [], "rft": []}
    for p in sorted(outlet.parent.iterdir()):
        m = pat.match(p.name)
        if m and p.is_file():
            found[m.group(2)].append(p)
    ins = [p for p in (lab / (stem + ".in"), outlet.parent / (stem + ".in")) if p.is_file()]
    return found["mrf"], found["rft"], ins[:1]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Collect the files Claude needs for window scoring into a new folder of zip files (see the top of this file).")
    ap.add_argument("--results", default="results")
    ap.add_argument("--obs", default="../init_data/obs")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--max-mb", type=float, default=20.0)
    ap.add_argument("--include-kchan", action="store_true")
    ap.add_argument("--extra", action="append", default=[])
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    here = Path.cwd()
    results = Path(a.results)
    if not results.is_dir():
        die("cannot find the results folder '{}' from {}.\n  Run this from workspaces/Sycamore/lab/ (the folder that contains results/), "
            "or point to it with --results.".format(a.results, here))
    obs_dir = Path(a.obs)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M")
    out_dir = Path(a.out_dir) if a.out_dir else here / "score_bundle_{}".format(stamp)
    if out_dir.exists():
        die("{} already exists. This script never writes into an existing folder; use --out-dir NAME for a new one.".format(out_dir))

    log("\ncollect_scoring_bundle_sycamore.py v{}   (read-only; writes only into a new folder)".format(VERSION))
    log("  working folder: {}".format(here))
    log("  now: {} local, {} UTC (Phoenix is UTC-7)".format(dt.datetime.now().strftime("%Y-%m-%d %H:%M"), dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M")))

    # ---- which runs
    outlets = sorted(p for p in results.rglob("*_Outlet.qout") if "_partial_" not in str(p) and p.is_file())
    by_stem = {}
    for p in outlets:
        by_stem.setdefault(p.name[:-len("_Outlet.qout")], []).append(p)
    dup = {s: ps for s, ps in by_stem.items() if len(ps) > 1}
    for s, ps in dup.items():
        log("  NOTE: more than one outlet file for {}; using {}".format(s, ps[0]))
    want, skipped = [], []
    extras = set(a.extra)
    for stem in sorted(by_stem):
        is_ours = stem.startswith(PREFIX_2017SEP) or stem.startswith(PREFIX_2020) or stem in OLD_RUNS or stem in extras
        if not is_ours:
            skipped.append((stem, "not one of the runs listed at the top of this script"))
        elif "kchan" in stem.lower() and not a.include_kchan and stem not in extras:
            skipped.append((stem, "Kchan run (use --include-kchan to add it)"))
        else:
            want.append(stem)
    expected = [PREFIX_2017SEP + t for t in EXPECT_2017SEP] + OLD_RUNS + sorted(extras)
    missing_expected = [s for s in expected if s not in by_stem]

    log("\nRuns found ({}):".format(len(want)))
    log("  {:<56} {:>6} {:>6} {:>6}".format("run", "mrf", "rft", "in"))
    plan = []          # (arcname, path)
    notes = []
    sel_info = []
    for stem in want:
        outlet = by_stem[stem][0]
        mrf, rft, ins = companions(outlet, stem, here)
        log("  {:<56} {:>6} {:>6} {:>6}".format(stem, len(mrf), len(rft), len(ins)))
        if len(mrf) != 1 or len(rft) != 1:
            notes.append("{}: found {} .mrf and {} .rft file(s) (one of each is expected); {}".format(
                stem, len(mrf), len(rft), "all of them are collected" if (mrf or rft) else "none to collect"))
        plan.append(("results/" + outlet.name, outlet))
        for p in mrf + rft:
            plan.append(("results/" + p.name, p))
        for p in ins:
            plan.append(("inputs/" + p.name, p))
        sel_info.append((stem, outlet, mrf, rft, ins))
    gauges = sorted(obs_dir.glob(GAUGE_GLOB)) if obs_dir.is_dir() else []
    log("\nGauge files found in {} ({}):".format(obs_dir, len(gauges)))
    for p in gauges:
        log("  " + p.name)
        plan.append(("gauge/" + p.name, p))
    if not gauges:
        notes.append("NO gauge file found in {} (looked for {}). Use --obs to point to the folder.".format(obs_dir, GAUGE_GLOB))
    if missing_expected:
        log("\nExpected but NOT found (so they will not be in the bundle):")
        for s in missing_expected:
            log("  " + s)
    if skipped:
        log("\nOutlet files seen but left out ({}):".format(len(skipped)))
        for s, why in skipped:
            log("  {}  [{}]".format(s, why))
    for n in notes:
        log("  NOTE: " + n)
    if not want:
        die("no runs to collect. Is this the right results folder ({})?".format(results.resolve()))
    if a.dry_run:
        log("\nDRY RUN: nothing was written.")
        return 0

    # ---- read every file once (checksums, line counts) and estimate its zipped size
    log("\nReading {} files (read-only)...".format(len(plan)))
    infos = []
    for arc, p in plan:
        d = describe(p)
        d["arc"] = arc
        d["zsize"] = len(zlib.compress(d["data"], 6))
        infos.append(d)
    total_raw = sum(d["size"] for d in infos)
    total_z = sum(d["zsize"] for d in infos)
    log("  {} files, {:.1f} MB as they are, about {:.1f} MB zipped".format(len(infos), total_raw / 1e6, total_z / 1e6))

    # ---- manifest
    M = []
    M.append("MANIFEST for the Lower Sycamore scoring bundle   (collect_scoring_bundle_sycamore.py v{})".format(VERSION))
    M.append("made {} local / {} UTC in {}".format(dt.datetime.now().strftime("%Y-%m-%d %H:%M"), dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M"), here))
    M.append("python {}   command: {}".format(sys.version.split()[0], " ".join(sys.argv)))
    M.append("")
    M.append("RUNS COLLECTED (outlet rows, first and last Time_hr, last outlet line)")
    for stem, outlet, mrf, rft, ins in sel_info:
        d = next(x for x in infos if x["path"] == outlet)
        rows = d["lines"] - 1
        try:
            tfirst = d["second"].split(",")[0].split()[0]
            tlast = d["last"].split(",")[0].split()[0]
        except Exception:
            tfirst = tlast = "?"
        M.append("  {}: outlet rows {}, Time_hr {} to {}; last line: {}; mrf {} rft {} in {}".format(
            stem, rows, tfirst, tlast, d["last"], len(mrf), len(rft), len(ins)))
    M.append("")
    if missing_expected:
        M.append("EXPECTED BUT NOT FOUND: " + ", ".join(missing_expected))
    if skipped:
        M.append("OUTLET FILES LEFT OUT: " + ", ".join("{} [{}]".format(s, w) for s, w in skipped))
    for n in notes:
        M.append("NOTE: " + n)
    M.append("")
    M.append("FILES (path in the zip | bytes | lines | md5 | first line | second line)")
    for d in infos:
        M.append("  {} | {} | {} | {} | {} | {}".format(d["arc"], d["size"], d["lines"], d["md5"], d["first"], d["second"]))
    manifest = "\n".join(M) + "\n"

    # ---- plan the zip parts (files are never split between zips)
    limit = max(1.0, a.max_mb - 2.5) * 1e6     # leave room for the manifest and for the estimate being a little off
    parts, cur, cur_z = [], [], 0
    for d in infos:
        if cur and cur_z + d["zsize"] > limit:
            parts.append(cur)
            cur, cur_z = [], 0
        cur.append(d)
        cur_z += d["zsize"]
    if cur:
        parts.append(cur)
    n = len(parts)

    out_dir.mkdir(parents=True)
    (out_dir / "MANIFEST.txt").write_text(manifest, encoding="utf-8")
    written = []
    for i, part in enumerate(parts, 1):
        zp = out_dir / "sycamore_scoring_bundle_part{}of{}.zip".format(i, n)
        with zipfile.ZipFile(zp, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
            if i == 1:
                zf.writestr("MANIFEST.txt", manifest)
            for d in part:
                zf.writestr(d["arc"], d["data"])
        written.append(zp)

    # ---- verify what was written
    log("\nChecking the zip file(s)...")
    names = set()
    bad = []
    for zp in written:
        with zipfile.ZipFile(zp) as zf:
            if zf.testzip() is not None:
                bad.append(zp.name)
            for zi in zf.infolist():
                names.add(zi.filename)
    want_names = {d["arc"] for d in infos} | {"MANIFEST.txt"}
    if bad or names != want_names:
        die("the zip check failed (bad: {}; missing: {}). Nothing in your results was changed; delete the folder {} and run again.".format(
            bad, sorted(want_names - names), out_dir))
    log("  OK: every file is in a zip and every zip passes its check.")

    log("\nDONE. New folder: {}".format(out_dir))
    for zp in written:
        log("  {}   {:.1f} MB".format(zp.name, zp.stat().st_size / 1e6))
    log("  MANIFEST.txt   (also inside the first zip)")
    log("\nNEXT: in the Codespace file list, open that folder, right-click each .zip file, choose Download, and attach the downloaded")
    log("file(s) to the chat. Nothing else is needed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
