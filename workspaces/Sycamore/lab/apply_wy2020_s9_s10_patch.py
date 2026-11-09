"""
apply_wy2020_s9_s10_patch.py  (v1.0)  --  adds two more screens, S9 and S10, to the S8 notebook

WHAT THEY ARE (each is S1 with exactly ONE value changed; S1 is Ks x1.0, f x0.15, As 300, bedrock SOLUS + 1 m, starting water table 0.95 x bedrock)

  S9-Ks1.0-Bedrock0.5     bedrock SOLUS + 0.5 m   (S1: + 1 m)        tag S9_Ks1p0_Bed0p5     the "thinner soil" run
  S10-Ks1.0-GW0.5         starting water table 0.5 x bedrock depth   (S1: 0.95)     tag S10_Ks1p0_GW0p5    the "wetter start" run

  Everything else (Ks x1.0, f x0.15, As 300, Au 1, routing, channel loss off, snow on) is the same as S1.
  Each setting is built from the S1 entry with dict(..., one_value=...), so nothing else can drift.

  Do not confuse these with S3 and S4 (already in the notebook, never run). S3 and S4 are built from Re-tuned (Ks x0.19261), so they
  carry the old Ks too. S9 and S10 are built from S1 (Ks x1.0), because S1 is the baseline the later screens were compared against.

  Read S9 with one caution: the starting water table is 0.95 x the bedrock depth, so when the bedrock moves up 0.5 m the starting water
  table moves up with it. S9 therefore changes the soil depth AND the starting water table together. S10 changes only the starting water table.
  The pair lets you separate them, but only roughly (S9 is not a pure soil-depth test).

WHAT THE SCRIPT DOES
  Reads your S8 notebook (Make_Sycamore_Model_WY2020_screens_S8.ipynb) and writes a NEW file, Make_Sycamore_Model_WY2020_screens_S10.ipynb.
  Your input notebook is never changed, and the script will not overwrite an existing output file unless you add --force.
  It adds two lines to the SETTINGS table in the Settings cell, plus a few words in the notes. Nothing else is edited.

HOW IT PROTECTS YOU
  - Every edit looks for one exact piece of text in one named cell and must find it exactly once. If your notebook differs from the
    copy I read, the script says which edit failed and writes nothing. (Edits marked "wording only" are skipped, not fatal.)
  - It stops if the notebook is not the S8 notebook (it needs Re-tuned, Original, the two Kchan runs and S1 to S8) or already has S9 or S10.
  - After editing, it checks that the Settings cell still parses, then runs ONLY the Settings cell (no files, no tRIBS) once per setting, and checks:
        * the twelve settings you already have give the same values as before,
        * S9 differs from S1 in exactly one value (BEDROCK_ADD_M, 1.0 -> 0.5) and S10 differs from S1 in exactly one value (GW_FRACTION, 0.95 -> 0.5),
        * the two new settings write input and results names that no other setting uses,
        * S9 writes bedrock_depth_plus0p5m.asc and InitGW_95pct_plus0p5m_mm.asc; S10 writes bedrock_depth_plus1m.asc (the same as S1) and
          InitGW_50pct_plus1m_mm.asc.
    It prints the table so you can read it. If any check fails, nothing is written.

RUN IT from workspaces/Sycamore/lab/ (nothing else is touched):
    python apply_wy2020_s9_s10_patch.py
    python apply_wy2020_s9_s10_patch.py --in Make_Sycamore_Model_WY2020_screens_S8.ipynb --out Make_Sycamore_Model_WY2020_screens_S10.ipynb
    python apply_wy2020_s9_s10_patch.py --force        # allow overwriting the output file

THEN, TO RUN THEM: use run_wy2020_screens_s9_s10.sh (it builds each input file from the new notebook and runs tRIBS, one run after the other).

VERSION HISTORY
  v1.0 (2026-11-09)  first version. Tested on a stand-in notebook (the Settings cell and the cells around it), not on your own file.
"""
import argparse
import ast
import contextlib
import io
import json
import os
import sys

CELL_INTRO = '69ebac8e'        # markdown: how to run it
CELL_SETTINGS_MD = '2fb78b48'  # markdown: notes on the Settings cell
CELL_SETTINGS = 'b5dfbf73'     # code: the Settings cell

S1 = 'S1-Ks1.0'
# new setting name -> (the one key that changes, its new value, tag)
NEW = {
    'S9-Ks1.0-Bedrock0.5': ('BEDROCK_ADD_M', 0.5, 'S9_Ks1p0_Bed0p5'),
    'S10-Ks1.0-GW0.5':     ('GW_FRACTION',   0.5, 'S10_Ks1p0_GW0p5'),
}
# the soil files each new setting must write (bedrock file, starting water table file)
EXPECTED_SOIL_FILES = {
    'S9-Ks1.0-Bedrock0.5': ('bedrock_depth_plus0p5m.asc', 'InitGW_95pct_plus0p5m_mm.asc'),
    'S10-Ks1.0-GW0.5':     ('bedrock_depth_plus1m.asc',   'InitGW_50pct_plus1m_mm.asc'),
}

S910_BLOCK = """# Two more screens (added 2026-11-09). Each is S1 (Ks x1.0, f x0.15, As 300, bedrock SOLUS + 1 m, starting water table 0.95 x bedrock) with ONE value changed:
# S9 has the bedrock at SOLUS + 0.5 m (the starting water table moves up with it, because it is 0.95 x bedrock depth); S10 has the starting
# water table at 0.5 x bedrock depth. These are built from S1, not from Re-tuned like S3 and S4 (which stay in the table, never run).
SETTINGS['S9-Ks1.0-Bedrock0.5'] = dict(SETTINGS['S1-Ks1.0'], BEDROCK_ADD_M=0.5, tag='S9_Ks1p0_Bed0p5')
SETTINGS['S10-Ks1.0-GW0.5'] = dict(SETTINGS['S1-Ks1.0'], GW_FRACTION=0.5, tag='S10_Ks1p0_GW0p5')
"""

# (cell id, text to find, text to put in its place, required, what it does)
EDITS = [
    (CELL_SETTINGS, "if SETTING not in SETTINGS:",
     S910_BLOCK + "if SETTING not in SETTINGS:", True,
     "add S9 and S10 to the SETTINGS table"),
    (CELL_SETTINGS, "are S5 with one value changed each (see below).",
     "are S5 with one value changed each (see below).\n"
     "# Two more (added 2026-11-09), S9 and S10, are S1 with one value changed each (see below).", False,
     "comment above SETTING"),
    (CELL_INTRO, "which are S5 with Ks x0.5, As 500 and f x0.02)",
     "which are S5 with Ks x0.5, As 500 and f x0.02, and S9 and S10, which are S1 with the bedrock 0.5 m shallower and the starting water table at 0.5)", False,
     "intro text"),
    (CELL_SETTINGS_MD, "or S6 to S8, which are S5 with one value changed)",
     "S6 to S8, which are S5 with one value changed, or S9 and S10, which are S1 with one value changed)", False,
     "Settings notes"),
    (CELL_SETTINGS_MD, "off for Re-tuned, Original and the S1 to S8 screens",
     "off for Re-tuned, Original and the S1 to S10 screens", False,
     "Settings notes (channel loss line)"),
]

KEYS = ['KS_MULT', 'F_MULT', 'ANISOTROPY', 'ANISOTROPY_U', 'BEDROCK_ADD_M', 'GW_FRACTION', 'OPTPERCOLATION', 'CHANNELCONDUCTIVITY']
OLD = ['Re-tuned', 'Original', 'Original-Kchan0.5', 'Original-Kchan0.2', 'S1-Ks1.0', 'S2-Ks0.5', 'S3-Bedrock0.5', 'S4-GW0.5',
       'S5-Ks1.0-Bedrock2.0', 'S6-Ks0.5-Bedrock2.0', 'S7-Ks1.0-Bedrock2.0-As500', 'S8-Ks1.0-Bedrock2.0-f0.02']


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
    ap.add_argument('--in', dest='src', default='Make_Sycamore_Model_WY2020_screens_S8.ipynb', help='the S8 notebook to read (never changed)')
    ap.add_argument('--out', dest='dst', default='Make_Sycamore_Model_WY2020_screens_S10.ipynb', help='the new notebook to write')
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
    names = OLD + list(NEW)
    before = settings_table(get_text(cell_by_id(nb, CELL_SETTINGS)), names)
    for name in OLD:
        if before[name] is None:
            sys.exit(f"STOP: the setting '{name}' is not in this notebook, so it is not the S8 notebook I read "
                     "(did the S6 to S8 patch run?). Nothing was written.")
    for name in NEW:
        if before[name] is not None:
            sys.exit(f'STOP: this notebook already has the setting {name}. Nothing was written.')

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

    after = settings_table(get_text(cell_by_id(nb, CELL_SETTINGS)), names)
    failures = []
    for name in OLD:
        if after[name] != before[name]:
            failures.append(f'{name} changed: before {before[name]} | after {after[name]}')
    s1 = after[S1]
    old_inputs = [after[n]['INPUT_FILE'] for n in OLD]
    old_results = [after[n]['RESULTS_NAME'] for n in OLD]
    new_inputs = []
    for name, (key, value, tag) in NEW.items():
        row = after[name]
        if row is None:
            failures.append(f'{name} is missing after the patch')
            continue
        diff = [k for k in KEYS if row[k] != s1[k]]
        if diff != [key]:
            failures.append(f'{name} differs from S1 in {diff} (needed exactly [{key!r}])')
        elif row[key] != value:
            failures.append(f'{name}: {key} is {row[key]} (needed {value})')
        if row['INPUT_FILE'] != f'Sycamore_WY2020_{tag}.in' or row['RESULTS_NAME'] != f'Sycamore_WY2020_{tag}':
            failures.append(f'{name}: input/results names are {row["INPUT_FILE"]} / {row["RESULTS_NAME"]} (needed tag {tag})')
        if row['INPUT_FILE'] in old_inputs or row['RESULTS_NAME'] in old_results:
            failures.append(f'{name} would write an input or results name an existing setting already uses')
        if row['SOIL_FILES'] != EXPECTED_SOIL_FILES[name]:
            failures.append(f'{name} would write soil files {row["SOIL_FILES"]} (needed {EXPECTED_SOIL_FILES[name]})')
        new_inputs.append(row['INPUT_FILE'])
    if len(set(new_inputs)) != len(new_inputs):
        failures.append('the two new settings would write the same input file')

    print('\nSettings after the patch (a * marks a value that differs from S1, the setting S9 and S10 are built from):')
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
    print('\nChecks passed: the twelve settings you already have are unchanged; S9 and S10 each differ from S1 in exactly one value '
          '(bedrock + 0.5 m, starting water table 0.5); their file names are not used by any other setting.')

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
    print('To run them: ./run_wy2020_screens_s9_s10.sh  (CHECK_ONLY=1 first to build the input files without running tRIBS).')


if __name__ == '__main__':
    main()
