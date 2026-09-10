# Offline notebook exports

Self-contained HTML copies of the lab notebooks, figures embedded, for anyone whose
network cannot load the Codespaces notebook editor.

These are generated, not hand-edited. To refresh them:

1. In a codespace, run `Make_SMF_Model.ipynb` and then `Run_Model.ipynb` all the way
   through, and save both.
2. Run `bash tools/export_notebooks.sh`.
3. Commit the resulting `.html` files.

Step 1 is not optional. The git clean filter (`tools/clean_notebook.py`) strips outputs
from every `.ipynb` on the way into git, so a fresh clone has no figures to export. The
filter matches `*.ipynb` only, so these HTML exports keep their figures while the
committed notebooks stay clean.
