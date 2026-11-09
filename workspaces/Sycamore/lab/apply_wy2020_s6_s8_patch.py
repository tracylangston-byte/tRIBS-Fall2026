"""
apply_wy2020_s6_s8_patch.py  (v1.0)  --  adds three more screens, S6, S7 and S8, to the S5 notebook

WHAT THEY ARE (each is S5 with exactly ONE value changed)
  S5 = Re-tuned values (f x0.15, As 300, Au 1) but Ks = SOLUS as is (x1.0) and the bedrock SOLUS + 2 m.

  S6-Ks0.5-Bedrock2.0         Ks x0.5            (S5: x1.0)          tag S6_Ks0p5_Bed2p0
  S7-Ks1.0-Bedrock2.0-As500   As 500             (S5: 300)           tag S7_Ks1p0_Bed2p0_As500     (Au stays 1, as in S5)
  S8-Ks1.0-Bedrock2.0-f0.02   f x0.02            (S5: x0.15)         tag S8_Ks1p0_Bed2p0_f0p02

  Everything else (bedrock + 2 m, starting water table 0.95 x bedrock, channel loss off, routing, snow on) is the same as S5.
  Each setting is built from the S5 entry with dict(..., one_value=...), so nothing else can drift.

WHAT THE SCRIPT DOES
  Reads your S5 notebook (Make_Sycamore_Model_WY2020_screens_S5.ipynb) and writes a NEW file, Make_Sycamore_Model_WY2020_screens_S8.ipynb.
  Your input notebook is never changed, and the script will not overwrite an existing output file unless you add --force.
  It adds three lines to the SETTINGS table in the Settings cell, plus a few words in the notes. Nothing else is edited.

HOW IT PROTECTS YOU
  - Every edit looks for one exact piece of text in one named cell and must find it exactly once. If your notebook differs from the
    copy I read, the script says which edit failed and writes nothing.
  - It stops if the notebook is not the S5 notebook (it needs Re-tuned, Original, the two Kchan runs and S1 to S5) or already has S6, S7 or S8.
  - After editing, it checks that the Settings cell still parses, then runs ONLY the Settings cell (no files, no tRIBS) once per setting, and checks:
        * the nine settings you already have give the same values as before,
        * S6, S7 and S8 each differ from S5 in exactly one value, and that value is the intended one,
        * the three new settings write input and results names that no other setting uses,
        * they use the same bedrock and starting-water-table files as S5 (same bedrock depth and same water-table fraction, so the same files).
    It prints the table so you can read it. If any check fails, nothing is written.

RUN IT from workspaces/Sycamore/lab/ (nothing else is touched):
    python apply_wy2020_s6_s8_patch.py
    python apply_wy2020_s6_s8_patch.py --in Make_Sycamore_Model_WY2020_screens_S5.ipynb --out Make_Sycamore_Model_WY2020_screens_S8.ipynb
    python apply_wy2020_s6_s8_patch.py --force        # allow overwriting the output file

THEN, TO RUN THEM: use run_wy2020_screens.sh (it builds each input file from the new notebook and runs tRIBS, one run after the other).

VERSION HISTORY
  v1.0 (2026-11-08)  first version. Tested on a stand-in notebook (the Settings cell and the cells around it), not on your own file.
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

S5 = 'S5-Ks1.0-Bedrock2.0'
# new setting name -> (the one key that changes, its new value, tag)
NEW = {
    'S6-Ks0.5-Bedrock2.0':       ('KS_MULT',    0.5,  'S6_Ks0p5_Bed2p0'),
    'S7-Ks1.0-Bedrock2.0-As500': ('ANISOTROPY', 500,  'S7_Ks1p0_Bed2p0_As500'),
    'S8-Ks1.0-Bedrock2.0-f0.02': ('F_MULT',     0.02, 'S8_Ks1p0_Bed2p0_f0p02'),
}

S68_BLOCK = """# Three more screens (added 2026-11-08). Each is S5 (Ks x1.0, bedrock SOLUS + 2 m, f x0.15, As 300) with ONE value changed:
# S6 has Ks x0.5, S7 has As 500, S8 has f x0.02. Nothing else differs from S5, so nothing else can drift.
SETTINGS['S6-Ks0.5-Bedrock2.0'] = dict(SETTINGS['S5-Ks1.0-Bedrock2.0'], KS_MULT=0.5, tag='S6_Ks0p5_Bed2p0')
SETTINGS['S7-Ks1.0-Bedrock2.0-As500'] = dict(SETTINGS['S5-Ks1.0-Bedrock2.0'], ANISOTROPY=500, tag='S7_Ks1p0_Bed2p0_As500')   # Au stays 1
SETTINGS['S8-Ks1.0-Bedrock2.0-f0.02'] = dict(SETTINGS['S5-Ks1.0-Bedrock2.0'], F_MULT=0.02, tag='S8_Ks1p0_Bed2p0_f0p02')
"""

# (cell id, text to find, text to put in its place, required, what it does)
EDITS = [
    (CELL_SETTINGS, "if SETTING not in SETTINGS:",
     S68_BLOCK + "if SETTING not in SETTINGS:", True,
     "add S6, S7 and S8 to the SETTINGS table"),
    (CELL_SETTINGS, "is S1 with the bedrock 1 m deeper (see below).",
     "is S1 with the bedrock 1 m deeper (see below).\n"
     "# Three more (added 2026-11-08), S6, S7 and S8, are S5 with one value changed each (see below).", False,
     "comment above SETTING"),
    (CELL_INTRO, "which is S1 with the bedrock 1 m deeper)",
     "which is S1 with the bedrock 1 m deeper, and S6, S7 and S8, which are S5 with Ks x0.5, As 500 and f x0.02)", False,
     "intro text"),
    (CELL_SETTINGS_MD, "or S5, which is S1 with the bedrock 1 m deeper)",
     "S5, which is S1 with the bedrock 1 m deeper, or S6 to S8, which are S5 with one value changed)", False,
     "Settings notes"),
    (CELL_SETTINGS_MD, "off for Re-tuned, Original and the S1 to S5 screens",
     "off for Re-tuned, Original and the S1 to S8 screens", False,
     "Settings notes (channel loss line)"),
]

KEYS = ['KS_MULT', 'F_MULT', 'ANISOTROPY', 'ANISOTROPY_U', 'BEDROCK_ADD_M', 'GW_FRACTION', 'OPTPERCOLATION', 'CHANNELCONDUCTIVITY']
OLD = ['Re-tuned', 'Original', 'Original-Kchan0.5', 'Original-Kchan0.2', 'S1-Ks1.0', 'S2-Ks0.5', 'S3-Bedrock0.5', 'S4-GW0.5', S5]


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
    ap.add_argument('--in', dest='src', default='Make_Sycamore_Model_WY2020_screens_S5.ipynb', help='the S5 notebook to read (never changed)')
    ap.add_argument('--out', dest='dst', default='Make_Sycamore_Model_WY2020_screens_S8.ipynb', help='the new notebook to write')
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
            sys.exit(f"STOP: the setting '{name}' is not in this notebook, so it is not the S5 notebook I read "
                     "(did the S5 patch run?). Nothing was written.")
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
    s5 = after[S5]
    old_inputs = [after[n]['INPUT_FILE'] for n in OLD]
    old_results = [after[n]['RESULTS_NAME'] for n in OLD]
    new_inputs = []
    for name, (key, value, tag) in NEW.items():
        row = after[name]
        if row is None:
            failures.append(f'{name} is missing after the patch')
            continue
        diff = [k for k in KEYS if row[k] != s5[k]]
        if diff != [key]:
            failures.append(f'{name} differs from S5 in {diff} (needed exactly [{key!r}])')
        elif row[key] != value:
            failures.append(f'{name}: {key} is {row[key]} (needed {value})')
        if row['INPUT_FILE'] != f'Sycamore_WY2020_{tag}.in' or row['RESULTS_NAME'] != f'Sycamore_WY2020_{tag}':
            failures.append(f'{name}: input/results names are {row["INPUT_FILE"]} / {row["RESULTS_NAME"]} (needed tag {tag})')
        if row['INPUT_FILE'] in old_inputs or row['RESULTS_NAME'] in old_results:
            failures.append(f'{name} would write an input or results name an existing setting already uses')
        if row['SOIL_FILES'] != s5['SOIL_FILES']:
            failures.append(f'{name} uses different bedrock/water-table files than S5: {row["SOIL_FILES"]} vs {s5["SOIL_FILES"]}')
        new_inputs.append(row['INPUT_FILE'])
    if len(set(new_inputs)) != len(new_inputs):
        failures.append('two of the new settings would write the same input file')

    print('\nSettings after the patch (a * marks a value that differs from S5, the setting S6 to S8 are built from):')
    width = max(len(n) for n in names)
    print(f"{'setting':{width}s}  " + '  '.join(f'{k:>13s}' for k in KEYS))
    for name in names:
        row = after[name]
        if row is None:
            continue
        cells = []
        for k in KEYS:
            mark = '*' if row[k] != s5[k] else ' '
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
    print('\nChecks passed: the nine settings you already have are unchanged; S6, S7 and S8 each differ from S5 in exactly one value '
          '(Ks x0.5, As 500, f x0.02); their file names are not used by any other setting.')

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
    print('To run them: ./run_wy2020_screens.sh  (CHECK_ONLY=1 first to build the input files without running tRIBS).')


if __name__ == '__main__':
    main()
