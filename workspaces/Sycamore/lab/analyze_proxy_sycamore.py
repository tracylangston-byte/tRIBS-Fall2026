#!/usr/bin/env python3
"""
analyze_proxy_sycamore.py  --  Stage 1, step 2: do the short windows rank the settings the way W4 does?

Reads what run_proxy_sycamore.py wrote (proxy_manifest.json and score_out_proxy/scores_proxy.csv) and, for the six
soil settings, compares the ORDER in which each short window (W1, W2) puts them with the order W4 puts them in.
It changes no run and writes only three small files in score_out_proxy/.

How the comparison works
------------------------
For each setting we take one error number from each window, rank the six settings from smallest error (best) to
largest, and measure how well two rankings agree with the Spearman rank correlation (rho):
    rho = +1   identical order        rho = 0   unrelated order        rho = -1   exactly reversed
The design's rule: a window is a valid screen if rho is 0.8 or more.
With 6 settings, rho = 0.943 is one swap of two neighbours, 0.886 is two such swaps and 0.829 is three (all pass).
A swap of two settings that are two places apart gives 0.771, and four neighbour swaps give 0.771 too (both fail).

Like for like (the approved method, "A"): the January flood week (Jan 18-25) is scored in all three windows, so W1,
W2 and the Jan 18-25 slice of each W4 run are compared on the same eight days. This isolates the effect of the cold
start. W2 vs W4 is also compared on the two February storms and on the January flood with its recession tail, which
W1 cannot see.
Whole window ("B", secondary): each window's own whole-window volume error, W1 and W2 against W4's Oct-Mar error.

Run it from workspaces/Sycamore/lab/ after run_proxy_sycamore.py has finished (or part way: it uses what is there):
    python analyze_proxy_sycamore.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MANIFEST = "proxy_manifest.json"
SCORE_DIR = "score_out_proxy"
SCORE_CSV = "scores_proxy.csv"
PASS_RHO = 0.8
MIN_N = 4                        # fewer than this many settings and a rank correlation means nothing

LABELS = {
    1: "Ks 0.2  PsiB 0.5 (baseline)",
    2: "Ks 0.2  PsiB 1.0",
    3: "Ks 0.1  PsiB 0.5",
    4: "Ks 0.3  PsiB 0.5",
    5: "Ks 0.3  PsiB 1.0",
    6: "Ks 0.1  PsiB 1.0",
}

# (key, title, scorer column, how to turn the column into an error, windows compared with W4)
#   "abs1"   |ratio - 1|      (a volume ratio, 1 = perfect)
#   "abs"    |value|          (a peak error in percent, 0 = perfect)
#   "signed" the value as is  (orders the settings by how much water the model makes, not by error)
PRIMARY = [
    ("flood_vol", "Flood week (Jan 18-25): volume error |model/gauge - 1|", "flood_week__vol_ratio", "abs1", ("W1", "W2")),
    ("flood_peak", "Flood week (Jan 18-25): peak error |%|", "flood_week__peak_err_pct", "abs", ("W1", "W2")),
]
EXTRA = [
    ("jan_tail_vol", "January flood with its tail (Jan 18 - Feb 3): volume error", "jan_flood_tail__vol_ratio", "abs1", ("W2",)),
    ("feb18_vol", "February 18-19 storm with tail: volume error", "feb18_19_tail__vol_ratio", "abs1", ("W2",)),
    ("feb27_vol", "February 27-28 storm with tail: volume error", "feb27_28_tail__vol_ratio", "abs1", ("W2",)),
]
WHOLE = [
    ("whole_vol", "Whole window volume error (W4 = Oct 1 - Apr 1)", "span__vol_ratio", "abs1", ("W1", "W2")),
]
SIGNED = [
    ("flood_vol_signed", "Flood week: volume ratio, signed (orders by water made)", "flood_week__vol_ratio", "signed", ("W1", "W2")),
    ("flood_peak_signed", "Flood week: peak error %, signed (orders by peak size)", "flood_week__peak_err_pct", "signed", ("W1", "W2")),
]


def spearman(a: pd.Series, b: pd.Series) -> tuple[float, int]:
    """Spearman rank correlation over the settings where both have a value. Ties get the average rank."""
    both = pd.concat([a, b], axis=1).dropna()
    n = len(both)
    if n < MIN_N:
        return np.nan, n
    ra, rb = both.iloc[:, 0].rank(), both.iloc[:, 1].rank()
    if ra.std() == 0 or rb.std() == 0:
        return np.nan, n
    return float(ra.corr(rb)), n


def transform(x: float, how: str) -> float:
    if x is None or not np.isfinite(x):
        return np.nan
    return {"abs1": abs(x - 1.0), "abs": abs(x), "signed": x}[how]


def order_text(values: pd.Series) -> str:
    """Settings from smallest to largest value, e.g. '2 < 1 < 3'; equal values (to 4 places) joined with '='."""
    v = values.dropna().sort_values(kind="mergesort")
    out, prev = [], None
    for k, x in v.items():
        if prev is not None:
            out.append("=" if round(x, 4) == round(prev, 4) else "<")
        out.append(str(k))
        prev = x
    return " ".join(out) if out else "-"


def load(lab: Path):
    mpath = lab / MANIFEST
    cpath = lab / SCORE_DIR / SCORE_CSV
    if not mpath.is_file() or not cpath.is_file():
        raise SystemExit(f"Need {mpath.name} and {SCORE_DIR}/{SCORE_CSV}: run run_proxy_sycamore.py first.")
    man = json.loads(mpath.read_text())
    scores = pd.read_csv(cpath).set_index("run_id")
    # one run per (setting, window). A W4 run made again by the batch beats a reused one if both exist.
    pick: dict = {}
    for rid, rec in man["runs"].items():
        if rec.get("status") not in ("done", "reused") or "setting" not in rec or rid not in scores.index:
            continue
        key = (int(rec["setting"]), rec["window"])
        rank = 1 if rec.get("kind") == "run" else 0
        if key not in pick or rank > pick[key][0]:
            pick[key] = (rank, rid)
    return {k: v[1] for k, v in pick.items()}, scores


def main() -> int:
    lab = Path(__file__).resolve().parent
    os.chdir(lab)
    runs, scores = load(lab)
    have = sorted({s for s, _ in runs})
    lines: list[str] = []

    def say(text: str = "") -> None:
        print(text)
        lines.append(text)

    say("PROXY CHECK: do W1 and W2 rank the soil settings the way W4 does?")
    say(f"A window passes if its ranking agrees with W4's at Spearman rho >= {PASS_RHO}. Six settings, so rho moves in coarse steps.")
    say()
    say("Runs found (setting: windows with a scored run):")
    for s in sorted(LABELS):
        w = [w for w in ("W1", "W2", "W4") if (s, w) in runs]
        say(f"  {s}  {LABELS[s]:<30} {', '.join(w) if w else 'none'}")
    missing = [(s, w) for s in LABELS for w in ("W1", "W2", "W4") if (s, w) not in runs]
    if missing:
        say(f"  ({len(missing)} of 18 setting/window results are missing; the comparisons below use the settings that have both windows.)")
    say()

    long_rows, rho_rows = [], []

    def table(spec):
        key, title, col, how, wins = spec
        vals = pd.DataFrame(index=sorted(LABELS), columns=["W1", "W2", "W4"], dtype=float)
        for (s, w), rid in runs.items():
            raw = scores.at[rid, col] if col in scores.columns else np.nan
            vals.loc[s, w] = transform(float(raw) if pd.notna(raw) else np.nan, how)
            long_rows.append(dict(metric=key, setting=s, window=w, run_id=rid, column=col, error_kind=how,
                                  value=vals.loc[s, w]))
        return vals

    def block(spec, heading=None):
        key, title, col, how, wins = spec
        vals = table(spec)
        say(heading or title)
        shown = [w for w in ("W1", "W2", "W4") if vals[w].notna().any()]
        say("  setting  " + "".join(f"{w:>9}" for w in shown) + "      " + "label")
        for s in vals.index:
            say(f"  {s:>5}    " + "".join(f"{vals.loc[s, w]:>9.3f}" if pd.notna(vals.loc[s, w]) else f"{'-':>9}" for w in shown)
                + f"      {LABELS[s]}")
        w4 = vals["W4"].dropna()
        out = {}
        for w in wins:
            rho, n = spearman(vals[w], vals["W4"])
            out[w] = rho
            verdict = ("not enough settings" if not np.isfinite(rho) else ("PASS" if rho >= PASS_RHO else "FAIL"))
            say(f"  {w} vs W4: rho = {'-' if not np.isfinite(rho) else f'{rho:.3f}'} over {n} settings  -> {verdict}"
                f"      order {w}: {order_text(vals[w])}   |   order W4: {order_text(vals['W4'])}")
            rho_rows.append(dict(metric=key, title=title, window=w, rho=rho, n=n, verdict=verdict,
                                 order_window=order_text(vals[w]), order_W4=order_text(vals["W4"])))
        if len(w4) >= 2:
            say(f"  (W4 values run from {w4.min():.3f} to {w4.max():.3f}; if they are nearly equal, their order is noise)")
        say()
        return out

    say("=" * 100)
    say("PRIMARY TEST (approved method A): the January flood week, scored the same way in every window")
    say("=" * 100)
    prim = {s[0]: block(s) for s in PRIMARY}
    for w in ("W1", "W2"):
        v, p = prim["flood_vol"].get(w, np.nan), prim["flood_peak"].get(w, np.nan)
        ok = [np.isfinite(x) and x >= PASS_RHO for x in (v, p)]
        word = "PASS on both" if all(ok) else ("PASS on volume only" if ok[0] else ("PASS on peak only" if ok[1] else "FAIL on both"))
        if not (np.isfinite(v) or np.isfinite(p)):
            word = "no result yet"
        say(f"  >> {w}: {word}   (volume rho {v:.3f}, peak rho {p:.3f})" if np.isfinite(v) and np.isfinite(p) else f"  >> {w}: {word}")
    say()

    say("=" * 100)
    say("W2 vs W4 on what W1 cannot see: the recession after the flood and the two February storms")
    say("=" * 100)
    for s in EXTRA:
        block(s)

    say("=" * 100)
    say("SECONDARY (method B): each window's own whole-window volume error against W4's Oct-Mar error")
    say("=" * 100)
    for s in WHOLE:
        block(s)

    say("=" * 100)
    say("DIAGNOSTIC: signed values. Do the windows at least agree on which setting makes more water / a bigger peak?")
    say("=" * 100)
    for s in SIGNED:
        block(s)

    say("How to read this: the first table is the test the design asks for. The extra blocks say whether W2 can also be trusted")
    say("for the February storms and the recession, which decide whether W2 can stand in for W4 in later stages. If W1 or W2")
    say("fails, look at the order lines to see WHICH settings swapped: a swap between two settings whose W4 errors are almost")
    say("equal is noise; a swap between settings far apart in W4 is a real disagreement.")

    out = lab / SCORE_DIR
    out.mkdir(exist_ok=True)
    pd.DataFrame(long_rows).drop_duplicates(["metric", "setting", "window"]).to_csv(out / "proxy_check_values.csv", index=False)
    pd.DataFrame(rho_rows).to_csv(out / "proxy_check_spearman.csv", index=False)
    (out / "proxy_check_report.txt").write_text("\n".join(lines) + "\n")
    print(f"\nWROTE {SCORE_DIR}/proxy_check_report.txt, proxy_check_values.csv, proxy_check_spearman.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
