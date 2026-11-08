"""
apply_wy2020_s5_patch.py  (v1.0)  --  adds ONE more screen, S5, to Make_Sycamore_Model_WY2020_screens.ipynb

WHAT S5 IS
  S5 is S1 (Ks = SOLUS as is) with the bedrock 1 m deeper:   bedrock = SOLUS + 2 m   (S1 and Re-tuned: SOLUS + 1 m).
  So S5 is ONE change from S1. It is NOT one change from Re-tuned (it differs from Re-tuned in two values: Ks and bedrock).
  The starting water table is still 0.95 x the bedrock depth, so it moves down with the bedrock; that is how the notebook is built.

  Setting name :  S5-Ks1.0-Bedrock2.0
  Input file   :  Sycamore_WY2020_S5_Ks1p0_Bed2p0.in           results name  Sycamore_WY2020_S5_Ks1p0_Bed2p0
  New soil files: bedrock_depth_plus2m.asc  and  InitGW_95pct_plus2m_mm.asc      (the SOLUS file itself is never changed)

WHAT IT DOES
  Reads your screens notebook (the one the first patch made) and writes a NEW file, Make_Sycamore_Model_WY2020_screens_S5.ipynb.
  Your input notebook is never changed, and the script will not overwrite an existing output file unless you add --force.
  It adds one line to the SETTINGS table in the Settings cell, plus a few words in the notes. Nothing else is edited: the bedrock
  and starting-water-table cells already build their file names from the setting (the first patch did that).

HOW IT PROTECTS YOU
  - Every edit looks for one exact piece of text in one named cell and must find it exactly once. If your notebook differs from
    the copy I read, the script says which edit failed and writes nothing.
  - It stops if the notebook is not the screens notebook (no S1 to S4) or already has S5.
  - After editing, it checks that the Settings cell still parses, then runs ONLY the Settings cell (no files, no tRIBS) once per
    setting, and checks that:
        * the eight settings you already have (Re-tuned, Original, the two Kchan runs, S1 to S4) give the same values as before,
        * S5 differs from S1 in exactly one value (BEDROCK_ADD_M, 1.0 -> 2.0),
        * S5 writes input and soil file names that no other setting uses.
    It prints the table so you can read it. If any check fails, nothing is written.

RUN IT from workspaces/Sycamore/lab/ (nothing else is touched):
    python apply_wy2020_s5_patch.py
    python apply_wy2020_s5_patch.py --in Make_Sycamore_Model_WY2020_screens.ipynb --out Make_Sycamore_Model_WY2020_screens_S5.ipynb
    python apply_wy2020_s5_patch.py --force        # allow overwriting the output file

THEN, TO RUN S5
  Open the new notebook. In the Settings cell either change the word after the comma to 'S5-Ks1.0-Bedrock2.0'
  (SETTING = os.environ.get('SYCAMORE_SETTING', 'S5-Ks1.0-Bedrock2.0')) or set SYCAMORE_SETTING in the shell. Because the soil folder
  data/model/ is shared by all settings, run the soil cells again (Execute Cell and Below from the Settings cell) before the model run.
  The Settings cell prints a line with 'bedrock SOLUS + 2 m'; check it says 2 before you start.

VERSION HISTORY
  v1.0 (2026-11-07)  first version. Tested on a copy of the five cells the first patch edits, not on your own file.
"""
import argparse
import ast
import contextlib
import io
import json
import os
import sys

# The cells that get edited, by their notebook cell ids (the same ids the first patch used).
CELL_INTRO = '69ebac8e'        # markdown: how to run it
CELL_SETTINGS_MD = '2fb78b48'  # markdown: notes on the Settings cell
CELL_SETTINGS = 'b5dfbf73'     # code: the Settings cell
CELL_BEDROCK = '7e90bc2a'      # code: deepen the bedrock grid (read only here, to check it builds its name from the setting)
CELL_GW = '2ee2f5dd'           # code: initial groundwater (read only here, same reason)

NEW_NAME = 'S5-Ks1.0-Bedrock2.0'

S5_BLOCK = """# One more screen (added 2026-11-07): S5 is S1 (Ks = SOLUS as is) with the bedrock 1 m deeper (SOLUS + 2 m instead of + 1 m).
# It is ONE change from S1, not from Re-tuned. The starting water table stays at GW_FRACTION x bedrock depth, so it moves down with the bedrock.
SETTINGS['S5-Ks1.0-Bedrock2.0'] = dict(SETTINGS['S1-Ks1.0'], BEDROCK_ADD_M=2.0, tag='S5_Ks1p0_Bed2p0')
"""

# (cell id, text to find, text to put in its place, required, what it does)
EDITS = [
    (CELL_SETTINGS, "if SETTING not in SETTINGS:",
     S5_BLOCK + "if SETTING not in SETTINGS:", True,
     "add S5 to the SETTINGS table"),
    (CELL_SETTINGS, "are one-change screens from Re-tuned (see below).",
     "are one-change screens from Re-tuned (see below).\n"
     "# One more, S5-Ks1.0-Bedrock2.0 (added 2026-11-07), is S1 with the bedrock 1 m deeper (see below).", False,
     "comment above SETTING"),
    (CELL_INTRO, "`S3-Bedrock0.5` and `S4-GW0.5`)",
     "`S3-Bedrock0.5` and `S4-GW0.5`, and `S5-Ks1.0-Bedrock2.0`, which is S1 with the bedrock 1 m deeper)", False,
     "intro text"),
    (CELL_SETTINGS_MD, "or a one-change screen S1 to S4)",
     "a one-change screen S1 to S4, or S5, which is S1 with the bedrock 1 m deeper)", False,
     "Settings notes"),
    (CELL_SETTINGS_MD, "off for Re-tuned, Original and the S1 to S4 screens",
     "off for Re-tuned, Original and the S1 to S5 screens", False,
     "Settings notes (channel loss line)"),
]

# The values compared between settings.
KEYS = ['KS_MULT', 'F_MULT', 'ANISOTROPY', 'ANISOTROPY_U', 'BEDROCK_ADD_M', 'GW_FRACTION', 'OPTPERCOLATION', 'CHANNELCONDUCTIVITY']
OLD = ['Re-tuned', 'Original', 'Original-Kchan0.5', 'Original-Kchan0.2', 'S1-Ks1.0', 'S2-Ks0.5', 'S3-Bedrock0.5', 'S4-GW0.5']


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


def file_names(bedrock_add_m, gw_fraction):
    """The soil file names the bedrock and water table cells build (same expressions as those cells)."""
    bed_tag = f'{bedrock_add_m:g}'.replace('.', 'p')
    return f'bedrock_depth_plus{bed_tag}m.asc', f'InitGW_{int(round(gw_fraction * 100))}pct_plus{bed_tag}m_mm.asc'


def settings_table(code, names):
    """Run ONLY the Settings cell (it needs nothing but os) once per setting; return {name: values or None}. No files are read."""
    out = {}
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
        out[name] = {k: ns[k] for k in KEYS}
        out[name]['INPUT_FILE'] = ns['INPUT_FILE']
        out[name]['RESULTS_NAME'] = ns['RESULTS_NAME']
        out[name]['SOIL_FILES'] = file_names(ns['BEDROCK_ADD_M'], ns['GW_FRACTION'])
    os.environ.pop('SYCAMORE_SETTING', None)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--in', dest='src', default='Make_Sycamore_Model_WY2020_screens.ipynb', help='the screens notebook to read (never changed)')
    ap.add_argument('--out', dest='dst', default='Make_Sycamore_Model_WY2020_screens_S5.ipynb', help='the new notebook to write')
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

    # ---- before: what the existing settings give today -------------------------------------------------------
    before = settings_table(get_text(cell_by_id(nb, CELL_SETTINGS)), OLD + [NEW_NAME])
    for name in OLD:
        if before[name] is None:
            sys.exit(f"STOP: the setting '{name}' is not in this notebook, so it is not the screens notebook I read "
                     "(did the first patch run?). Nothing was written.")
    if before[NEW_NAME] is not None:
        sys.exit(f'STOP: this notebook already has the setting {NEW_NAME}. Nothing was written.')

    # The bedrock and water table cells must build their file names from the setting (the first patch did that).
    if "bed_tag = f'{BEDROCK_ADD_M:g}'.replace('.', 'p')" not in get_text(cell_by_id(nb, CELL_BEDROCK)):
        sys.exit('STOP: the bedrock cell does not build its file name from the setting the way the first patch does. Nothing was written.')
    if "bed_tag = f'{BEDROCK_ADD_M:g}'.replace('.', 'p')" not in get_text(cell_by_id(nb, CELL_GW)):
        sys.exit('STOP: the starting-water-table cell does not build its file name from the setting the way the first patch does. Nothing was written.')

    # ---- edits ----------------------------------------------------------------------------------------------
    print('Edits:')
    problems = []
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
                print('  skipped ' + msg + ' (wording only; the patch still works)')
            continue
        set_text(cell, text.replace(old, new))
        print(f'  done    cell {cid}: {what}')
    if problems:
        sys.exit('\nSTOP: the notebook is not the copy I read, so nothing was written. Send me the lines marked FAILED.')

    # ---- checks ---------------------------------------------------------------------------------------------
    ast.parse(get_text(cell_by_id(nb, CELL_SETTINGS)))          # raises SyntaxError if an edit broke the cell
    print('\nThe edited Settings cell still parses.')

    after = settings_table(get_text(cell_by_id(nb, CELL_SETTINGS)), OLD + [NEW_NAME])
    failures = []
    for name in OLD:
        if after[name] != before[name]:
            failures.append(f'{name} changed: before {before[name]} | after {after[name]}')
    s5, s1, rt = after[NEW_NAME], after['S1-Ks1.0'], after['Re-tuned']
    if s5 is None:
        failures.append(f'{NEW_NAME} is missing after the patch')
    else:
        diff1 = [k for k in KEYS if s5[k] != s1[k]]
        if diff1 != ['BEDROCK_ADD_M']:
            failures.append(f'S5 differs from S1 in {diff1} (needed exactly BEDROCK_ADD_M)')
        elif s5['BEDROCK_ADD_M'] != 2.0 or s1['BEDROCK_ADD_M'] != 1.0:
            failures.append(f"BEDROCK_ADD_M is S1 {s1['BEDROCK_ADD_M']}, S5 {s5['BEDROCK_ADD_M']} (needed 1.0 and 2.0)")
        others_in = [after[n]['INPUT_FILE'] for n in OLD]
        others_soil = [after[n]['SOIL_FILES'] for n in OLD]
        if s5['INPUT_FILE'] in others_in or s5['RESULTS_NAME'] in [after[n]['RESULTS_NAME'] for n in OLD]:
            failures.append('S5 would write an input or results name another setting already uses')
        if s5['SOIL_FILES'][0] in [x[0] for x in others_soil] or s5['SOIL_FILES'][1] in [x[1] for x in others_soil]:
            failures.append('S5 would write a soil file name another setting already uses')

    names = OLD + [NEW_NAME]
    print('\nSettings after the patch (a * marks a value that differs from S1, the setting S5 is built from):')
    width = max(len(n) for n in names)
    print(f"{'setting':{width}s}  " + '  '.join(f'{k:>13s}' for k in KEYS))
    for name in names:
        row = after[name]
        if row is None:
            continue
        cells = []
        for k in KEYS:
            mark = '*' if row[k] != s1[k] else ' '
            cells.append(f'{str(row[k]):>12s}{mark}')
        print(f'{name:{width}s}  ' + '  '.join(cells))
    print('\nInput files and soil files the settings write:')
    for name in names:
        if after[name]:
            print(f'  {name:{width}s} -> {after[name]["INPUT_FILE"]}   {after[name]["SOIL_FILES"][0]}   {after[name]["SOIL_FILES"][1]}')

    if failures:
        print('\nCHECKS FAILED, nothing was written:')
        for f in failures:
            print('  - ' + f)
        sys.exit(1)
    print('\nChecks passed: the eight settings you already have are unchanged, S5 differs from S1 in exactly one value '
          '(bedrock 1 m -> 2 m), and S5 writes file names no other setting uses.')

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
    print("To run S5: in the Settings cell use SETTING = ... 'S5-Ks1.0-Bedrock2.0', rerun the soil cells (Execute Cell and Below), "
          "and check the printed line says 'bedrock SOLUS + 2 m'.")


if __name__ == '__main__':
    main()
