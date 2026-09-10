#!/usr/bin/env bash
# Export the lab notebooks to self-contained HTML for offline / low-bandwidth use.
#
# Why this exists: when the venue network is slow, the Codespaces web client can
# fail to load the notebook editor and VS Code falls back to offering the .ipynb
# as a text file. These HTML exports let someone read and follow along.
#
# IMPORTANT: run the notebooks first. tools/clean_notebook.py (the git clean
# filter) strips outputs from every .ipynb on the way into git, so a fresh clone
# has no figures to export. Run Make_SMF_Model.ipynb then Run_Model.ipynb in the
# codespace, save them, and only then run this script. The filter does not touch
# HTML, so the exports keep their figures while the committed .ipynb stays clean.
#
# Usage:
#   bash tools/export_notebooks.sh            # export working-tree notebooks as-is
#   bash tools/export_notebooks.sh --execute  # run each notebook first, then export

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LAB_DIR="$REPO_ROOT/workspaces/SMF/lab"
OUT_DIR="$REPO_ROOT/docs/notebooks"

EXECUTE=0
if [[ "${1:-}" == "--execute" ]]; then
    EXECUTE=1
elif [[ $# -gt 0 ]]; then
    echo "usage: $0 [--execute]" >&2
    exit 2
fi

if ! python3 -c "import nbconvert" 2>/dev/null; then
    echo "ERROR: nbconvert is not installed. Install it with:" >&2
    echo "  pip install nbconvert" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

# Ordered: Make_SMF_Model builds the model that Run_Model then runs.
NOTEBOOKS=(
    "Make_SMF_Model.ipynb"
    "Run_Model.ipynb"
)

missing_outputs=()

for nb in "${NOTEBOOKS[@]}"; do
    src="$LAB_DIR/$nb"
    if [[ ! -f "$src" ]]; then
        echo "skip: $nb not found at $src"
        continue
    fi

    if [[ $EXECUTE -eq 0 ]]; then
        # Warn loudly rather than silently shipping an export with no figures.
        if ! python3 - "$src" <<'PY'
import json, sys
nb = json.load(open(sys.argv[1]))
sys.exit(0 if any(c.get("outputs") for c in nb["cells"]) else 1)
PY
        then
            missing_outputs+=("$nb")
        fi
    fi

    echo "exporting $nb ..."
    args=(--to html --embed-images --output-dir "$OUT_DIR")
    if [[ $EXECUTE -eq 1 ]]; then
        # Run in the notebook's own directory; the lab reads and writes relative paths.
        args+=(--execute --ExecutePreprocessor.timeout=1800)
    fi

    ( cd "$LAB_DIR" && python3 -m nbconvert "${args[@]}" "$nb" )

    # nbconvert's --embed-images does not reliably inline images that markdown
    # cells reference by relative path, which would leave broken images once the
    # export is served from docs/. Inline them ourselves.
    python3 "$REPO_ROOT/tools/inline_images.py" "$OUT_DIR/${nb%.ipynb}.html" "$LAB_DIR"
done

echo
echo "HTML written to docs/notebooks/"

if [[ ${#missing_outputs[@]} -gt 0 ]]; then
    echo
    echo "WARNING: these notebooks had no cell outputs, so their exports have no figures:"
    for nb in "${missing_outputs[@]}"; do
        echo "  - $nb"
    done
    echo "Run them in the codespace and save, then re-run this script (or pass --execute)."
fi
