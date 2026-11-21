"""
apply_wy2017sep_patch.py  (v1.0)  --  makes the WY2017 "September 2016 spin-up" notebook from your S11 notebook

WHAT IT MAKES
  Reads  Make_Sycamore_Model_WY2020_screens_S11.ipynb  (the notebook that holds all of your settings: Re-tuned, Original, the two Kchan
  runs and S1 to S11) and writes a NEW file,  Make_Sycamore_Model_WY2017sep.ipynb.  Your S11 notebook is never changed, and the script
  will not overwrite an existing output file unless you add --force.

  The new notebook is the S11 notebook with the year changed. It changes ONLY these things:
      run window     Sep 1 2016 00:00, 6,551 hours (to May 31 2017 23:00)        (WY2020: Sep 1 2019, 6,575 hours)
      radar rain     first file radar_0901201601.asc, last file radar_0531201723.asc, folder ~/mrms_wk_wy2017_sep
      weather        ../init_data/met_wy2017_sep/Master_Met.sdf
      run names      Sycamore_WY2017sep_<tag>.in  and results/Sycamore_WY2017sep_<tag>*   (WY2020: Sycamore_WY2020_<tag>)
      some wording   the title, the notes and a few comments, so they no longer say "WY2020"
  It does NOT change any setting. Every setting (Ks, f, PsiB, As, Au, bedrock, starting water table, channel loss, routing, snow)
  keeps exactly the value it has in your S11 notebook, so the ten runs are the same ten runs you did for WY2020.

  September 2016 is spin-up: the run starts with the same starting state as every other run (water table 0.95 x bedrock depth) one
  month earlier, so the model has a month to settle before the scored window opens on Oct 1 2016. Nothing here scores anything.

  Why 6,551 hours: Sep 1 2016 00:00 to May 31 2017 23:00 is 273 days = 6,552 hourly weather rows, and the run is one hour shorter
  (the WY2020 run was 6,575 of 6,576, and the WY2017 full-year runs were 8,759 of 8,760). Your weather files run on to Sep 30 2017, so
  6,552 hours would also work; 6,551 keeps the same end time (23:00) as the WY2020 runs.

WHY A PATCH AND NOT A WHOLE NEW NOTEBOOK
  Same reason as before: this edits your own copy and checks the result, instead of retyping a notebook this size by hand.

HOW IT PROTECTS YOU
  - Every edit that matters looks for one exact piece of text in the Settings cell and must find it exactly once. If your notebook
    differs from the copy I read, the script says which edit failed and writes nothing. (Edits marked "wording only" are skipped, not
    fatal: they are title and note text.)
  - It stops if the notebook is not the S11 notebook (it needs all 15 settings), is not a WY2020 notebook, or is already a WY2017sep one.
  - After editing, it runs ONLY the Settings cell (no files, no tRIBS) once per setting, before and after the patch, and checks:
        * all 15 settings give exactly the same values after the patch as before it (Ks, f, PsiB, As, Au, bedrock, water table,
          channel loss, channel conductivity); only the input and results names change,
        * the ten settings you will run have the values written in the WY2020 notes (Re-tuned, Original, S1, S5 to S11),
        * every setting writes its own input and results name, of the form Sycamore_WY2017sep_<tag>,
        * the run window, radar file names, radar folder and weather list are the new ones, and the first and last radar file
          names agree with the start date and run length (so a typo in one of them is caught),
        * no executable line of the Settings cell still says 2019 or 2020,
        * no cell other than the Settings cell and the five text cells it edits has changed.
    It prints the table so you can read it. If any check fails, nothing is written.

RUN IT from workspaces/Sycamore/lab/ (nothing else is touched):
    python apply_wy2017sep_patch.py
    python apply_wy2017sep_patch.py --in Make_Sycamore_Model_WY2020_screens_S11.ipynb --out Make_Sycamore_Model_WY2017sep.ipynb
    python apply_wy2017sep_patch.py --force        # allow overwriting the output file

THEN, TO RUN THE TEN RUNS: use run_wy2017sep_series.sh (it builds each input file from the new notebook, checks it, and runs tRIBS,
one run after the other). Start with  CHECK_ONLY=1 ./run_wy2017sep_series.sh  (builds and checks all ten, runs no tRIBS).

VERSION HISTORY
  v1.0 (2026-11-21)  first version, based on apply_wy2020_s11_patch.py v1.0. Tested on a stand-in notebook (the real Settings cell and
                     notes as the earlier patches leave them), not on your own file.
"""
import argparse
import ast
import contextlib
import io
import json
import os
import sys
from datetime import datetime, timedelta

CELL_INTRO = '69ebac8e'        # markdown: how to run it
CELL_SETTINGS_MD = '2fb78b48'  # markdown: notes on the Settings cell
CELL_SETTINGS = 'b5dfbf73'     # code: the Settings cell
CELL_WEATHER_MD = '20cb6f72'   # markdown: check the weather files against the run window
CELL_RADAR_MD = '32618ba7'     # markdown: check the radar files against the run window
CELL_DONE_MD = 'c853aa7e'      # markdown: the "Done" cell
EDITED_CELLS = {CELL_INTRO, CELL_SETTINGS_MD, CELL_SETTINGS, CELL_WEATHER_MD, CELL_RADAR_MD, CELL_DONE_MD}

# ---- the new run window and inputs -----------------------------------------------------------------------------------
NEW_START = '09/01/2016/00/00'
NEW_HOURS = 6551
NEW_FIRST = 'radar_0901201601.asc'
NEW_LAST = 'radar_0531201723.asc'
NEW_MET = '../init_data/met_wy2017_sep/Master_Met.sdf'
NEW_ROOT = '~/mrms_wk_wy2017_sep'
NEW_RUN_KEY = 'wy2017sep'
NEW_PREFIX = 'Sycamore_WY2017sep_'
OLD_START = '09/01/2019/00/00'
OLD_HOURS = 6575
OLD_PREFIX = 'Sycamore_WY2020_'

# (cell id, text to find, text to put in its place, required, what it does).  Each must be found exactly once in its cell.
EDITS = [
    # --- the Settings cell: the edits that matter (required) ---
    (CELL_SETTINGS, "'../init_data/met_wy2020/Master_Met.sdf'", f"'{NEW_MET}'", True, 'weather station list'),
    (CELL_SETTINGS, "'~/mrms_wk_wy2020'", f"'{NEW_ROOT}'", True, 'folder searched for the radar files'),
    (CELL_SETTINGS, "'radar_0901201901.asc'", f"'{NEW_FIRST}'", True, 'first radar file'),
    (CELL_SETTINGS, "'radar_0531202023.asc'", f"'{NEW_LAST}'", True, 'last radar file'),
    (CELL_SETTINGS, "'wy2020': ('09/01/2019/00/00', 6575),", f"'{NEW_RUN_KEY}': ('{NEW_START}', {NEW_HOURS}),", True, 'run window'),
    (CELL_SETTINGS, "RUN_CHOICE = 'wy2020'", f"RUN_CHOICE = '{NEW_RUN_KEY}'", True, 'which run window is used'),
    (CELL_SETTINGS, 'INPUT_FILE = f"Sycamore_WY2020_', f'INPUT_FILE = f"{NEW_PREFIX}', True, 'input file name'),
    (CELL_SETTINGS, 'RESULTS_NAME = f"Sycamore_WY2020_', f'RESULTS_NAME = f"{NEW_PREFIX}', True, 'results name'),
    # --- the Settings cell: comments (wording only) ---
    (CELL_SETTINGS, "# WY2020 weather, made by Generate_Met_Forcing_Sycamore_WY2020.ipynb",
     "# WY2017 weather with a Sep 2016 start (Sep 1 2016 to Sep 30 2017), made by a copy of Generate_Met_Forcing_Sycamore.ipynb", False,
     'comment on the weather list'),
    (CELL_SETTINGS, "# first end-labeled hour: 09/01/2019 00:00 to 01:00", "# first end-labeled hour: 09/01/2016 00:00 to 01:00", False,
     'comment on the first radar file'),
    (CELL_SETTINGS, "# last hour the run uses: the hour ending 05/31/2020 23:00 (the weather files end at 23:00)",
     "# last hour the run uses: the hour ending 05/31/2017 23:00 (the weather files run on to 09/30/2017 23:00)", False,
     'comment on the last radar file'),
    (CELL_SETTINGS, "# Starts Sep 1 2019 so the Sep 10 and Sep 23 storms are in the run (Sep is spin-up, not scored; scored window Oct 1 - May 31).",
     "# Starts Sep 1 2016: September is spin-up, not scored (the scored window is Oct 1 2016 - May 31 2017).", False, 'comment on the run window (1)'),
    (CELL_SETTINGS, "# 6575 hours, not 6576: the weather files end 05/31/2020 23:00, and a run to 06/01 00:00 would go one hour past them",
     "# 6551 hours (273 days less one hour), so the run ends 05/31/2017 23:00, the same clock time the WY2020 run ended. The weather files", False,
     'comment on the run window (2)'),
    (CELL_SETTINGS, "# (the same reason the WY2017 full-year runs used 8759).",
     "# run on to 09/30/2017 23:00, so 6552 hours would also work; 6551 keeps the same end time as WY2020.", False,
     'comment on the run window (3)'),
    (CELL_SETTINGS, "# ---- WHICH SETTING: Re-tuned or Original (WY2020 validation run) ---",
     "# ---- WHICH SETTING: Re-tuned, Original or a screen (WY2017 run with a Sep 2016 spin-up) ---", False, 'comment above SETTING (1)'),
    (CELL_SETTINGS, "# run_wy2020_both.sh sets SYCAMORE_SETTING for you.", "# run_wy2017sep_series.sh sets SYCAMORE_SETTING for you.", False,
     'comment above SETTING (2)'),
    # --- the text cells (wording only) ---
    (CELL_INTRO, "water year 2020 (validation year, Sep 1 2019 start)",
     "water year 2017 with a September 2016 spin-up (Sep 1 2016 to May 31 2017; September is spin-up, not scored)", False, 'title'),
    (CELL_INTRO, "`../init_data/met_wy2020/Master_Met.sdf`, plus the `.mdf` files it lists | `Generate_Met_Forcing_Sycamore_WY2020.ipynb` |",
     "`../init_data/met_wy2017_sep/Master_Met.sdf`, plus the `.mdf` files it lists | a copy of `Generate_Met_Forcing_Sycamore.ipynb` with the Sep 2016 start |",
     False, 'inputs table: weather'),
    (CELL_INTRO, "or let `run_wy2020_both.sh` do both settings. WY2020 flow and model results stay sealed until the scoring plan is written: "
                 "do not open the results or compare with the gauge before then.",
     "or let `run_wy2017sep_series.sh` run all ten settings one after the other.", False, 'how to run tRIBS'),
    (CELL_INTRO, "**WY2020 setting:** this copy has one run window (`wy2020`)",
     "**WY2017 (Sep spin-up) setting:** this copy has one run window (`wy2017sep`)", False, 'how to run: setting'),
    (CELL_INTRO, "which `run_wy2020_both.sh` sets", "which `run_wy2017sep_series.sh` sets", False, 'how to run: environment variable'),
    (CELL_SETTINGS_MD,
     "`wy2020` starts 09/01/2019 00:00 and runs 6,575 hours. The radar files are labeled by the **end** of each hour, so the run first asks for "
     "`radar_0901201901.asc` and last for `radar_0531202023.asc`. The weather files cover 09/01/2019 00:00 to 05/31/2020 23:00.",
     "`wy2017sep` starts 09/01/2016 00:00 and runs 6,551 hours (to 05/31/2017 23:00). September 2016 is spin-up and is not scored. The radar files "
     "are labeled by the **end** of each hour, so the run first asks for `radar_0901201601.asc` and last for `radar_0531201723.asc`. The weather files "
     "cover 09/01/2016 00:00 to 09/30/2017 23:00.", False, 'Settings notes: run window'),
    (CELL_SETTINGS_MD, "`Sycamore_WY2020_Retuned` and `Sycamore_WY2020_Original`",
     "`Sycamore_WY2017sep_Retuned` and `Sycamore_WY2017sep_Original`", False, 'Settings notes: run names'),
    (CELL_WEATHER_MD,
     "The WY2020 weather files end at 05/31/2020 23:00, so this run is set to 6,575 hours (ending 05/31/2020 23:00). "
     "A run of 6,576 hours would go one hour past the weather and this check would say so.",
     "The WY2017 weather files (made with a September 2016 start) cover 09/01/2016 00:00 to 09/30/2017 23:00. This run is set to 6,551 hours "
     "(ending 05/31/2017 23:00), so the weather covers it with room to spare.", False, 'weather check text'),
    (CELL_RADAR_MD,
     "(The 17 hours that were missing from the source archive were written as zero-rain grids, so no file should be missing.)",
     "(Hours missing from the source archive were written as zero-rain grids, so no file should be missing. One of them is 09/11/2016 13:00, "
     "in the September spin-up; the hours on both sides of it are dry.)", False, 'radar check text'),
    (CELL_DONE_MD, "this 6,575-hour run", "this 6,551-hour run", False, 'Done cell: run length'),
    (CELL_DONE_MD,
     "`run_wy2020_both.sh` does this for both settings. WY2020 results are sealed: do not open them or compare with the gauge until the scoring plan is written.",
     "`run_wy2017sep_series.sh` does this for all ten settings, one after the other. September 2016 is spin-up: do not score it.", False,
     'Done cell: script and sealed note'),
]

KEYS = ['KS_MULT', 'F_MULT', 'PSIB_MULT', 'ANISOTROPY', 'ANISOTROPY_U', 'BEDROCK_ADD_M', 'GW_FRACTION', 'OPTPERCOLATION', 'CHANNELCONDUCTIVITY']
EXTRA = ['INPUT_FILE', 'RESULTS_NAME', 'START_DATE', 'RUN_HOURS', 'RUN_CHOICE', 'MET_STATIONS', 'RADAR_ROOT', 'RADAR_FIRST_FILE', 'RADAR_LAST_FILE',
         'OPTSNOW', 'UTC_OFFSET']
ALL15 = ['Re-tuned', 'Original', 'Original-Kchan0.5', 'Original-Kchan0.2', 'S1-Ks1.0', 'S2-Ks0.5', 'S3-Bedrock0.5', 'S4-GW0.5',
         'S5-Ks1.0-Bedrock2.0', 'S6-Ks0.5-Bedrock2.0', 'S7-Ks1.0-Bedrock2.0-As500', 'S8-Ks1.0-Bedrock2.0-f0.02',
         'S9-Ks1.0-Bedrock0.5', 'S10-Ks1.0-GW0.5', 'S11-Ks1.0-PsiB1.0']

# The ten settings that will be run, with the values written in the WY2020 notes (Re-tuned and Original from the Settings cell notes;
# S1 and S5 to S11 from apply_wy2020_screens_patch.py, apply_wy2020_s5_patch.py, apply_wy2020_s9_s10_patch.py and apply_wy2020_s11_patch.py).
_RT = dict(KS_MULT=0.19261, F_MULT=0.15, PSIB_MULT=0.2939, ANISOTROPY=300, ANISOTROPY_U=1, BEDROCK_ADD_M=1.0, GW_FRACTION=0.95, OPTPERCOLATION=0)
_S1 = dict(_RT, KS_MULT=1.0)
_S5 = dict(_S1, BEDROCK_ADD_M=2.0)
EXPECT = {
    'Re-tuned': (_RT, 'Retuned'),
    'Original': (dict(_RT, KS_MULT=0.187, F_MULT=0.5874, ANISOTROPY=1), 'Original'),
    'S1-Ks1.0': (_S1, 'S1_Ks1p0'),
    'S5-Ks1.0-Bedrock2.0': (_S5, 'S5_Ks1p0_Bed2p0'),
    'S6-Ks0.5-Bedrock2.0': (dict(_S5, KS_MULT=0.5), 'S6_Ks0p5_Bed2p0'),
    'S7-Ks1.0-Bedrock2.0-As500': (dict(_S5, ANISOTROPY=500), 'S7_Ks1p0_Bed2p0_As500'),
    'S8-Ks1.0-Bedrock2.0-f0.02': (dict(_S5, F_MULT=0.02), 'S8_Ks1p0_Bed2p0_f0p02'),
    'S9-Ks1.0-Bedrock0.5': (dict(_S1, BEDROCK_ADD_M=0.5), 'S9_Ks1p0_Bed0p5'),
    'S10-Ks1.0-GW0.5': (dict(_S1, GW_FRACTION=0.5), 'S10_Ks1p0_GW0p5'),
    'S11-Ks1.0-PsiB1.0': (dict(_S1, PSIB_MULT=1.0), 'S11_Ks1p0_PsiB1p0'),
}


def cell_by_id(nb, cid):
    hits = [c for c in nb['cells'] if c.get('id') == cid]
    if len(hits) != 1:
        raise SystemExit(f'STOP: expected exactly one cell with id {cid}, found {len(hits)}. Nothing was written.')
    return hits[0]


def get_text(cell):
    return ''.join(cell['source'])


def set_text(cell, text):
    """Store text the way Jupyter does: a list of lines, each ending in a newline except the last."""
    parts = text.split('\n')
    lines = [p + '\n' for p in parts[:-1]]
    if parts[-1] != '':
        lines.append(parts[-1])
    cell['source'] = lines


def settings_table(code, names):
    """Run ONLY the Settings cell (it needs nothing but os) once per setting; return {name: values or None}. No files are read."""
    out = {}
    saved = {k: os.environ.pop(k, None) for k in ('SYCAMORE_RADAR_ROOT', 'SYCAMORE_RADAR_DIR')}   # so the cell's own defaults are read
    try:
        for name in names:
            os.environ['SYCAMORE_SETTING'] = name
            ns = {'os': os}
            buf = io.StringIO()
            try:
                with contextlib.redirect_stdout(buf):
                    exec(compile(code, '<settings cell>', 'exec'), ns)
            except ValueError:
                out[name] = None            # the setting does not exist in this version of the cell
                continue
            out[name] = {k: ns[k] for k in KEYS + EXTRA}
    finally:
        os.environ.pop('SYCAMORE_SETTING', None)
        for k, v in saved.items():
            if v is not None:
                os.environ[k] = v
    return out


def tag_of(row, prefix):
    name = row['INPUT_FILE']
    if not (name.startswith(prefix) and name.endswith('.in')):
        return None
    return name[len(prefix):-len('.in')]


def executable_part(line):
    """The part of a line of Python before any # comment (a # inside a quoted string is kept)."""
    quote = None
    i = 0
    while i < len(line):
        ch = line[i]
        if quote:
            if ch == '\\':
                i += 2
                continue
            if ch == quote:
                quote = None
        elif ch in ('"', "'"):
            quote = ch
        elif ch == '#':
            return line[:i]
        i += 1
    return line


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--in', dest='src', default='Make_Sycamore_Model_WY2020_screens_S11.ipynb', help='the S11 notebook to read (never changed)')
    ap.add_argument('--out', dest='dst', default='Make_Sycamore_Model_WY2017sep.ipynb', help='the new notebook to write')
    ap.add_argument('--force', action='store_true', help='allow overwriting the output file')
    args = ap.parse_args()

    if not os.path.isfile(args.src):
        sys.exit(f'STOP: {args.src} not found. Run this from workspaces/Sycamore/lab/, or give the path with --in.')
    if os.path.exists(args.dst) and not args.force:
        sys.exit(f'STOP: {args.dst} already exists. Use --out for another name, or --force to overwrite it.')
    if os.path.abspath(args.src) == os.path.abspath(args.dst):
        sys.exit('STOP: the input and output are the same file. This script never overwrites the notebook it reads.')

    with open(args.src, encoding='utf-8') as fh:
        nb = json.load(fh)
    original_text = {c.get('id'): get_text(c) for c in nb['cells']}

    # ---- before: what the settings give today -----------------------------------------------------------------
    before = settings_table(get_text(cell_by_id(nb, CELL_SETTINGS)), ALL15)
    for name in ALL15:
        if before[name] is None:
            sys.exit(f"STOP: the setting '{name}' is not in this notebook, so it is not the S11 notebook I read "
                     '(did the S9, S10 and S11 patches run?). Nothing was written.')
    rt = before['Re-tuned']
    if rt['RUN_CHOICE'] == NEW_RUN_KEY or rt['INPUT_FILE'].startswith(NEW_PREFIX):
        sys.exit('STOP: this notebook is already a WY2017sep notebook. Nothing was written.')
    if rt['START_DATE'] != OLD_START or rt['RUN_HOURS'] != OLD_HOURS or not rt['INPUT_FILE'].startswith(OLD_PREFIX):
        sys.exit(f"STOP: this is not the WY2020 notebook I read (start {rt['START_DATE']}, {rt['RUN_HOURS']} hours, input file "
                 f"{rt['INPUT_FILE']}; I expected {OLD_START}, {OLD_HOURS} hours, {OLD_PREFIX}...). Nothing was written.")
    # The settings must already be right in the notebook you give me, so the patch only changes the year.
    failures = []
    for name, (want, tag) in EXPECT.items():
        for k, v in want.items():
            if before[name][k] != v:
                failures.append(f'(before the patch) {name}: {k} is {before[name][k]} (the WY2020 notes say {v})')
        if tag_of(before[name], OLD_PREFIX) != tag:
            failures.append(f'(before the patch) {name}: tag is {tag_of(before[name], OLD_PREFIX)} (the WY2020 notes say {tag})')
    if failures:
        print('The settings in the notebook you gave me are not the ones written in the WY2020 notes:')
        for f in failures:
            print('  - ' + f)
        sys.exit('\nSTOP: nothing was written. Send me these lines.')

    # ---- edits ------------------------------------------------------------------------------------------------
    print('Edits:')
    problems = []
    n_done = n_skipped = 0
    for cid, old, new, required, what in EDITS:
        cell = cell_by_id(nb, cid)
        text = get_text(cell)
        n = text.count(old)
        if n != 1:
            msg = f'cell {cid}: looked for the text for "{what}" and found it {n} times (needed exactly 1)'
            if required:
                problems.append(msg)
                print('  FAILED  ' + msg)
            else:
                n_skipped += 1
                print('  skipped ' + msg + ' (wording only; the patch still works)')
            continue
        set_text(cell, text.replace(old, new))
        n_done += 1
        print(f'  done    cell {cid}: {what}')
    if problems:
        sys.exit('\nSTOP: the notebook is not the copy I read, so nothing was written. Send me the lines marked FAILED.')
    print(f'\n{n_done} edits made, {n_skipped} wording-only edits skipped.')

    # ---- checks -----------------------------------------------------------------------------------------------
    failures = []
    settings_code = get_text(cell_by_id(nb, CELL_SETTINGS))
    ast.parse(settings_code)                                     # raises SyntaxError if an edit broke the cell
    print('The edited Settings cell still parses.')

    after = settings_table(settings_code, ALL15)
    for name in ALL15:
        if after[name] is None:
            failures.append(f'{name} is missing after the patch')
            continue
        for k in KEYS:
            if after[name][k] != before[name][k]:
                failures.append(f'{name}: {k} changed from {before[name][k]} to {after[name][k]} (a setting must not change)')
        for k in ('OPTSNOW', 'UTC_OFFSET'):
            if after[name][k] != before[name][k]:
                failures.append(f'{name}: {k} changed from {before[name][k]} to {after[name][k]}')
        t_old, t_new = tag_of(before[name], OLD_PREFIX), tag_of(after[name], NEW_PREFIX)
        if t_new is None or t_new != t_old:
            failures.append(f"{name}: input file is {after[name]['INPUT_FILE']} (needed {NEW_PREFIX}{t_old}.in)")
        elif after[name]['RESULTS_NAME'] != f'{NEW_PREFIX}{t_old}':
            failures.append(f"{name}: results name is {after[name]['RESULTS_NAME']} (needed {NEW_PREFIX}{t_old})")
    inputs = [after[n]['INPUT_FILE'] for n in ALL15 if after[n]]
    if len(set(inputs)) != len(inputs):
        failures.append('two settings would write the same input file name')

    for name, (want, tag) in EXPECT.items():
        row = after[name]
        if row is None:
            continue
        for k, v in want.items():
            if row[k] != v:
                failures.append(f'{name}: {k} is {row[k]} (the WY2020 notes say {v})')
        if tag_of(row, NEW_PREFIX) != tag:
            failures.append(f'{name}: tag is {tag_of(row, NEW_PREFIX)} (expected {tag})')

    r = after['Re-tuned']
    for key, want in [('START_DATE', NEW_START), ('RUN_HOURS', NEW_HOURS), ('RUN_CHOICE', NEW_RUN_KEY), ('MET_STATIONS', NEW_MET),
                      ('RADAR_ROOT', NEW_ROOT), ('RADAR_FIRST_FILE', NEW_FIRST), ('RADAR_LAST_FILE', NEW_LAST)]:
        if r[key] != want:
            failures.append(f'{key} is {r[key]!r} (needed {want!r})')
    # the radar file names must agree with the start date and run length (files are labeled by the END of each hour)
    t0 = datetime.strptime(r['START_DATE'], '%m/%d/%Y/%H/%M')
    first = f"radar_{(t0 + timedelta(hours=1)):%m%d%Y%H}.asc"
    last = f"radar_{(t0 + timedelta(hours=r['RUN_HOURS'])):%m%d%Y%H}.asc"
    if first != r['RADAR_FIRST_FILE'] or last != r['RADAR_LAST_FILE']:
        failures.append(f"the start date and run length need radar files {first} to {last}, but Settings says {r['RADAR_FIRST_FILE']} to {r['RADAR_LAST_FILE']}")
    end = t0 + timedelta(hours=r['RUN_HOURS'])
    if end != datetime(2017, 5, 31, 23, 0):
        failures.append(f'the run ends {end:%Y-%m-%d %H:%M} (needed 2017-05-31 23:00)')
    for name in ALL15:
        if after[name] and any(after[name][k] != r[k] for k in ('START_DATE', 'RUN_HOURS', 'RUN_CHOICE', 'MET_STATIONS', 'RADAR_ROOT',
                                                                  'RADAR_FIRST_FILE', 'RADAR_LAST_FILE')):
            failures.append(f'{name} has a different run window or input from Re-tuned')

    # no executable line of the Settings cell may still say 2019 or 2020 (comments are not held to this)
    for i, line in enumerate(settings_code.splitlines(), start=1):
        code_part = executable_part(line)
        if '2019' in code_part or '2020' in code_part:
            failures.append(f'Settings cell, line {i}, still says 2019 or 2020: {line.strip()[:100]}')

    # only the cells I meant to edit changed
    changed = [c.get('id') for c in nb['cells'] if get_text(c) != original_text.get(c.get('id'))]
    unexpected = [cid for cid in changed if cid not in EDITED_CELLS]
    if unexpected:
        failures.append(f'cells changed that should not have: {unexpected}')
    if len(nb['cells']) != len(original_text):
        failures.append('the number of cells changed')

    # leftovers in the text and comments, for you to read (not a failure)
    left = []
    for c in nb['cells']:
        for i, line in enumerate(get_text(c).splitlines(), start=1):
            if '2019' in line or '2020' in line:
                left.append((c.get('id'), c['cell_type'], i, line.strip()))
    print(f'\nLines that still mention 2019 or 2020 anywhere in the notebook (text and comments; you may want to read them): {len(left)}')
    for cid, ctype, i, line in left[:12]:
        print(f'  {ctype:8s} cell {cid} line {i}: {line[:110]}')
    if len(left) > 12:
        print(f'  ... and {len(left) - 12} more')

    names = [n for n in EXPECT]
    print('\nThe ten settings you will run, after the patch (every value is the same as in your S11 notebook):')
    width = max(len(n) for n in names)
    cols = ['KS_MULT', 'F_MULT', 'PSIB_MULT', 'ANISOTROPY', 'ANISOTROPY_U', 'BEDROCK_ADD_M', 'GW_FRACTION', 'OPTPERCOLATION']
    print(f"{'setting':{width}s}  " + '  '.join(f'{k:>13s}' for k in cols))
    for name in names:
        row = after[name]
        print(f'{name:{width}s}  ' + '  '.join(f'{str(row[k]):>13s}' for k in cols))
    print('\nInput files and results names:')
    for name in names:
        print(f"  {name:{width}s} -> {after[name]['INPUT_FILE']}   results/{after[name]['RESULTS_NAME']}*")
    print(f"\nRun window: {r['START_DATE']}, {r['RUN_HOURS']} hours (ends {end:%Y-%m-%d %H:%M}); radar {r['RADAR_FIRST_FILE']} to {r['RADAR_LAST_FILE']} "
          f"in {r['RADAR_ROOT']}; weather {r['MET_STATIONS']}")

    if failures:
        print('\nCHECKS FAILED, nothing was written:')
        for f in failures:
            print('  - ' + f)
        sys.exit(1)
    print('\nChecks passed: all 15 settings give the same values as before (only the input and results names changed); the ten settings you will run '
          'match the WY2020 notes; the run window, radar files and weather list are the WY2017 September-spin-up ones and agree with each other; '
          'only the Settings cell and the five text cells changed.')

    # ---- write ------------------------------------------------------------------------------------------------
    try:
        import nbformat
        nbformat.validate(nbformat.from_dict(nb))
        print('nbformat says the notebook is valid.')
    except ImportError:
        print('(nbformat is not installed, so the notebook structure was not validated that way.)')
    with open(args.dst, 'w', encoding='utf-8') as fh:
        json.dump(nb, fh, indent=1, ensure_ascii=False, sort_keys=True)
        fh.write('\n')
    print(f'\nWrote {args.dst}. Your {args.src} was not changed.')
    print('To run the ten runs: chmod +x run_wy2017sep_series.sh, then  CHECK_ONLY=1 ./run_wy2017sep_series.sh  first (builds and checks all ten, no tRIBS).')


if __name__ == '__main__':
    main()
