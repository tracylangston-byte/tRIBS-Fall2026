#!/usr/bin/env python3
"""Inline local <img> sources in an exported notebook as base64 data URIs.

nbconvert's --embed-images does not reliably inline images that markdown cells
reference by relative path (e.g. ../assets/diagram.png), which leaves the export
with broken images once it is served from docs/. This pass makes the HTML truly
self-contained, so a single file can be served, downloaded, or emailed and still
render every figure.

Usage: inline_images.py <html-file> <base-dir>

<base-dir> is the directory that relative src paths resolve against, i.e. the
directory the notebook lives in. Already-inlined (data:) and remote (http:) srcs
are left alone. Exits non-zero if a local image could not be resolved.
"""

import base64
import mimetypes
import re
import sys
from pathlib import Path

IMG_SRC = re.compile(r'(<img\b[^>]*?\bsrc=")([^"]+)(")', re.IGNORECASE)


def main() -> int:
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} <html-file> <base-dir>", file=sys.stderr)
        return 2

    html_path = Path(sys.argv[1])
    base_dir = Path(sys.argv[2])
    html = html_path.read_text(encoding="utf-8")

    inlined, missing = [], []

    def replace(match: re.Match) -> str:
        prefix, src, suffix = match.groups()
        if src.startswith(("data:", "http://", "https://", "//")):
            return match.group(0)

        target = (base_dir / src).resolve()
        if not target.is_file():
            missing.append(src)
            return match.group(0)

        mime = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
        payload = base64.b64encode(target.read_bytes()).decode("ascii")
        inlined.append(src)
        return f"{prefix}data:{mime};base64,{payload}{suffix}"

    html = IMG_SRC.sub(replace, html)

    if inlined:
        html_path.write_text(html, encoding="utf-8")
        for src in inlined:
            print(f"    inlined {src}")

    for src in missing:
        print(f"    MISSING {src}", file=sys.stderr)

    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
