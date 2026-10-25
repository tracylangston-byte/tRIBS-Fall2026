"""
make_poster_hydrographs_110.py   (version 4.2)
=============================================
Series 110 -- the POSTER hydrograph figure: the real gauge and the best run for each of the three big storms,
three panels side by side, with a legend and a caption underneath. No tRIBS runs; it only reads files and
draws one picture.

WHAT CHANGED IN VERSION 4.2
---------------------------
  * The bold lead-in of the caption now reads: "tRIBS 6.0: Preliminary results show that no single Ks-f-Kchan set fits all
    three real 2014 storms, which points to forcing or initial-state limits." (your wording; DEFAULT_LEAD_OK near the top of
    the script). As before, the script checks the sweep tables first: if some run DOES have the real peak and volume within a
    factor of 2 on all three storms, the lead-in switches to a count instead of printing a claim the data do not support.

WHAT CHANGED IN VERSION 4 (the type now matches the Series 100 figure)
-----------------------------------------------------------------------
  * ONE text size, 24 pt at printed size, for ticks, axis label, panel titles, subtitles, legend AND caption (the Series 100
    figure measures 24 pt throughout; the caption used to be 20 pt and the titles 28 pt).
  * Panel titles are bold navy (#0C234B). Everything else is dark ink (#0F1825), the same ink as the Series 100 figure.
  * Caption = a bold lead-in sentence followed by regular text, like the Series 100 captions. Version 4.1: the first paragraph ends
    with "Not comparable with tRIBS 5.3.1 results." and the "All panels:" note is gone (--all_panels brings it back: build label,
    "Preliminary", the KGE definition, the symbol definitions and the fixed routing parameters).
  * Symbols are italic with real subscripts and superscripts: Ks (italic K, subscript s), f, Kchan (italic K, subscript chan),
    cv, r, n, mm^-1. Channel conductivity is now called Kchan everywhere in the figure (it was "cc").
  * The lead-in is your sentence (see version 4.2 above for the current wording). The script counts, from
    the sweep tables, how many runs have the real peak AND volume within a factor of 2 on all three storms (your 0.5-2.0
    "ballpark" band) and prints that count in the caption. If the count is not 0, the lead-in changes to say so instead of
    printing a claim the data do not support. Use --lead "..." to write your own lead-in.

WHAT IT DRAWS, IN PLAIN WORDS
------------------------------
Three panels: 12 Aug, 8 Sep, 27 Sep. Three runs are picked, one per storm: the run with the highest KGE_2012 on that
storm alone (among runs that made at least 1 cfs). Each run has its own colour in every panel:
  * blue   the best run for 12 Aug          * green  the best run for 8 Sep          * orange  the best run for 27 Sep
In each panel:
  * black line           the real gauge (gaps filled the same way as in score_storms_110.py)
  * thick coloured line  the run that was picked for THAT panel's storm
  * dotted coloured line the same run from another storm, shown in this storm (to see how it does there)
The dotted lines are: green and orange on 12 Aug; blue and orange on 8 Sep; blue and green on 27 Sep (the ms1 runs make no
flow on 27 Sep, so they lie along zero; leave them out with --no_dotted_on_27sep). Each panel is zoomed in on its event and has its own flow axis. Under each panel title a second line
says which sweep that panel's thick run came from:
  * original Ks range  = sweep ms1 (Ks multiplier 4.5 to 8.0)
  * low Ks range       = sweep ms2 (Ks multiplier 1.5 to 4.5)
By default 12 Aug and 8 Sep come from ms1 and 27 Sep from ms2 (no ms1 run makes any flow on 27 Sep).
Under the panels: a legend (two short rows) and a caption (24 pt) that lists the Ks, f and Kchan of
each of the three runs. The "best run, storms combined" line of version 1 is NOT drawn unless you add --show_combined.

HOW IT USES YOUR SCORER
------------------------
It does not copy the scoring: it imports score_storms_110.py (it must sit in the same folder) and calls its own
functions: the same windows file, the same gauge reading and gap filling, the same .qout reading, the same KGE_2012,
the same rescaling and the same two self-checks (GATE 1 and GATE 2). So the numbers match the scorer's tables. It never
changes score_storms_110.py, and it never writes outside its own output folder.

IMPORTANT TO KNOW (also printed at the end of every run)
---------------------------------------------------------
  * Each "best run" is picked by hindsight on one storm. It shows what the model CAN do on that storm. The dotted lines
    show what that same parameter set does on the other storms.
  * The 27 Sep run comes from a different sweep (ms2) than the other two (ms1), so its Ks range differs by design.
  * 27 Sep is flagged WEAK by the scorer's own rule (best raw KGE_2012 under 0.5); the script prints this.
  * Same two rain gauges for every storm, one monsoon season (2014), routing pinned, gauge error not quantified.

THE FONT
--------
All text in the figure (titles, axis labels, legend and caption) is set in Liberation Sans, which has the same letter widths as
Arial (so it matches an Arial poster) and is usually installed on Linux. The script says which font it found when it starts.
If Liberation Sans is missing it tries Arial, Arimo, Helvetica and Nimbus Sans; if none exists it uses DejaVu Sans (wider) and
says so, with the command to install it. To use a particular font file: --font_file path/to/font.ttf (a single file has no
bold or italic, so bold and italic text would look regular; the script says so).
The caption is laid out by the script itself (matplotlib cannot mix bold, italic and subscripts in one wrapped paragraph), by
measuring each word in the real font. Markup used inside the script's caption text: **bold**, *italic*, _{subscript}, ^{superscript}.

USAGE (from the lab/ folder)
-----------------------------
    python make_poster_hydrographs_110.py --check_only        # gates only; writes nothing
    python make_poster_hydrographs_110.py                     # gates, then the figure
    python make_poster_hydrographs_110.py --font_file Arial.ttf     # use a particular font file instead
    python make_poster_hydrographs_110.py --no_dotted_on_27sep # leave the blue and green runs off the 27 Sep panel
    python make_poster_hydrographs_110.py --show_combined     # add the best run of the combined score back in
    python make_poster_hydrographs_110.py --height_in 7.0     # fix the total height (the panels shrink to fit)
    python make_poster_hydrographs_110.py --caption_kge       # put each run's KGE in the caption
    python make_poster_hydrographs_110.py --all_panels        # add the "All panels:" note back under the caption
    python make_poster_hydrographs_110.py --all_panels --brief   # ... without the definitions of Ks, f and Kchan
    python make_poster_hydrographs_110.py --lead "Your own bold lead-in sentence."   # replace the bold lead-in (markup above, without **)

WHAT YOU GET   (calibration_work/03_comparisons/summary_tables/poster_110/)
-----------------------------------------------------------------------------
    fig_poster_hydrographs_110.png   22.5 in wide, 300 dpi; the height is worked out from the panels plus the legend and
                                     caption (printed when it runs)
    poster_best_runs_110.csv         for each storm: which run is drawn, its Ks, f and cc, and how it does on all three storms
    PROVENANCE_poster_110.json       checksums, versions, gate results, settings, the caption text
"""

import argparse
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import score_storms_110 as ss
except ImportError:
    sys.exit("STOPPED: cannot import score_storms_110.py. Put this script in the same folder as score_storms_110.py "
             "(the lab/ folder) and run it from there.")

VERSION = "4.2"

# ------------------------------------------------------------------
# Look of the figure (change these before a run if you want a different look)
# ------------------------------------------------------------------
# One text size for everything, as in the Series 100 figure (measured there: 24 pt, line pitch 1.25 x the "lp" height).
POSTER_PT = 24                  # ticks, axis label, subtitles, legend: points at printed size
TITLE_PT = 24                   # panel titles (bold navy; they were 28 pt before version 4)
CAPTION_PT = 24                 # caption (it was 20 pt before version 4)
TEXT_INK = "#0F1825"            # dark ink of all text (colour measured in the Series 100 figure)
HEAD_NAVY = "#0C234B"           # Arizona Navy: panel titles
SCRIPT_SCALE, SUB_DROP, SUP_RISE = 0.70, 0.18, 0.38     # subscript/superscript size, and shift as a share of the text size
LINE_PITCH, PARA_GAP = 1.25, 0.30                       # caption line spacing (x the height of "lp"); extra gap between paragraphs (x one line)
BALLPARK_LO, BALLPARK_HI = 0.5, 2.0   # your "ballpark" band (Series110_Compiled_Record.md 4.16): model peak and volume each 0.5 to 2.0 x real
RUN_COLOR = {"12 Aug": "#36a2ff", "8 Sep": "#4cd05a", "27 Sep": "#ff9a3d"}     # light but bright blue, green, orange
RUN_COLOR_NAME = {"12 Aug": "blue", "8 Sep": "green", "27 Sep": "orange"}
C_COMBINED = "#6b3fa0"          # only used with --show_combined
LW_REAL, LW_OWN, LW_CROSS, LW_COMBINED = 3.0, 5.5, 4.5, 4.0
DOTTED = (0, (0.05, 1.8))       # round dots: dash length ~0, gap 1.8 line-widths (used with dash_capstyle="round")
DOTTED_PERIOD = 0.05 + 1.8       # one dot plus its gap, in line-widths
CROSS_DEFAULT = {"12 Aug": ["8 Sep", "27 Sep"], "8 Sep": ["12 Aug", "27 Sep"], "27 Sep": ["12 Aug", "8 Sep"]}   # whose run is dotted where
# Routing was held fixed in every ms1 and ms2 run (Series110_Compiled_Record.md, 4.14 and 4.16: cv 4.5, flowexp 0.24, n 0.026).
ROUTING_TEXT = "*cv* (velocity coefficient) 4.5, *r* (velocity exponent) 0.24, *n* (channel roughness) 0.026"
FONT_CANDIDATES = ["Liberation Sans", "Arial", "Arimo", "Helvetica", "Nimbus Sans"]   # in this order
DEFAULT_PANELS = "12 Aug=ms1,8 Sep=ms1,27 Sep=ms2"
SWEEP_NAME = {"ms1": "original *K*_{s} range", "ms2": "low *K*_{s} range"}      # caption/subtitle markup (see strip_markup)
KS_M, F_M, KCHAN_M = "*K*_{s}", "*f*", "*K*_{chan}"      # symbols as markup: Ks, f, Kchan (channel conductivity; "cc" in the other scripts)
NBSP = "\u00a0"                 # a space that never breaks across lines (keeps "tRIBS 5.3.1", "12 Aug" and "168 mm/hr" together)
DEFAULT_LEAD_OK = ("tRIBS" + NBSP + "6.0: Preliminary results show that no single %s-%s-%s set fits all three real 2014 storms, "
                   "which points to forcing or initial-state limits." % (KS_M, F_M, KCHAN_M))
OUT_SUBFOLDER = "poster_110"
SIDE_MARGIN_IN = 0.4            # caption starts this far from the left edge and stops this far from the right


# ------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------
def parse_panels(text):
    """'12 Aug=ms1,8 Sep=ms1,27 Sep=ms2' -> [('12 Aug', 'ms1'), ...]. Each storm must be one of the three big ones."""
    out = []
    for part in text.split(","):
        if "=" not in part:
            sys.exit("STOPPED: --panels entry '%s' needs the form 'storm=sweep', for example '27 Sep=ms2'." % part.strip())
        storm, sweep = [x.strip() for x in part.split("=", 1)]
        if storm not in ss.BIG:
            sys.exit("STOPPED: '%s' is not one of the big storms %s." % (storm, ss.BIG))
        if not sweep:
            sys.exit("STOPPED: --panels entry for '%s' has no sweep label." % storm)
        out.append((storm, sweep))
    if len({s for s, _ in out}) != len(out):
        sys.exit("STOPPED: a storm appears twice in --panels.")
    return out


def show_gates(ck, title):
    print("\nGATES: %s" % title)
    ck.show()


def stop_if_failed(ck, title):
    if ck.failures():
        show_gates(ck, title)
        sys.exit("\nSTOPPED: nothing was drawn. Fix the cause above.")


def clock_label(x, pos=None):
    """Tick label hh:mm; at midnight the date as well (the 12 Aug event crosses midnight)."""
    import matplotlib.dates as mdates
    d = mdates.num2date(x)
    if d.hour == 0 and d.minute == 0:
        return d.strftime("00:00\n%d %b")
    return d.strftime("%H:%M")


def hour_step(span_hours, max_ticks=5):
    """Smallest whole-hour step that divides 24 and gives at most max_ticks ticks (ticks land on 00:00, 02:00, ...)."""
    for k in (1, 2, 3, 4, 6, 8, 12, 24):
        if span_hours / float(k) <= max_ticks:
            return k
    return 24


# ------------------------------------------------------------------
# The real gauge (same recipe and same checks as score_storms_110.py main())
# ------------------------------------------------------------------
def build_real(storms, xlsx_path, summary_dir, ck):
    raw = ss.read_gauge(xlsx_path)
    real = {}
    for s in storms:
        asis, filled, ftype = ss.build_real_window(raw, s["ws"], s["we"], ss.DEFAULT_MAX_GAP_MIN, ss.DEFAULT_MAX_HOLD_HR)
        n_asis, n_unc = int(asis.notna().sum()), int((ftype == "uncovered").sum())
        real[s["tag"]] = {"asis": asis if n_asis >= ss.MIN_ASIS_POINTS else None, "asis_full": asis, "filled": filled,
                          "ftype": ftype,
                          "real_peak_cfs": float(np.nanmax(filled.to_numpy())) / ss.CFS_TO_CMS,
                          "real_vol_acft": float(np.nansum(filled.to_numpy())) * ss.GRID_SEC / ss.M3_PER_ACFT}
        if n_unc:
            ck.add("WARN", "gauge, " + s["label"], "%d of %d 5-minute bins are uncovered by any reading" % (n_unc, len(filled)))
    stage0 = summary_dir / "real_gauge_110" / "real_gauge_observed_110.csv"
    if stage0.exists() and (storms[0]["ws"], storms[0]["we"]) == (ss.STAGE_B_START, ss.STAGE_B_END):
        try:
            o = pd.read_csv(stage0, parse_dates=["datetime"]).set_index("datetime")
            r12 = real[storms[0]["tag"]]
            same_idx = o.index.equals(r12["filled"].index)
            a_ok = np.allclose(o["Observed_asis_m3s"].to_numpy(float), r12["asis_full"].to_numpy(float), rtol=1e-9, atol=1e-12, equal_nan=True)
            f_ok = np.allclose(o["Observed_filled_m3s"].to_numpy(float), r12["filled"].to_numpy(float), rtol=1e-9, atol=1e-12, equal_nan=True)
            ck.add("PASS" if (same_idx and a_ok and f_ok) else "FAIL", "GATE 1",
                   "the 12 Aug real series built here %s %s from Stage 0" % ("equals" if (same_idx and a_ok and f_ok) else "DIFFERS from", stage0.name))
        except Exception as exc:
            ck.add("WARN", "GATE 1", "could not compare with %s: %s" % (stage0.name, str(exc)[:80]))
    else:
        ck.add("NOTE", "GATE 1", "%s not found: the 12 Aug real series was not compared with Stage 0" % stage0.name)
    return real


# ------------------------------------------------------------------
# One sweep: gates, scoring, combined score, best runs
# ------------------------------------------------------------------
def load_sweep(label, ctx):
    ck = ss.Checks()
    storms, runtime, real = ctx["storms"], ctx["runtime"], ctx["real"]
    design = ss.read_design(ctx["summary_dir"] / ("lhs_design_lowf_%s_ON_110.json" % label), label, ck)
    stop_if_failed(ck, label)
    if design["runtime_hours"] is None or int(design["runtime_hours"]) != runtime:
        ck.add("FAIL", "run length", "the %s design says %s h, the windows file says %d h" % (label, design["runtime_hours"], runtime))
    res_path = ctx["summary_dir"] / ("lhs_results_lowf_%s_ON_110.csv" % label)
    res, partial = ss.load_results(res_path, design, ck, False)
    stop_if_failed(ck, label)
    if not ctx["keep_dir"].exists():
        ck.add("FAIL", "kept .qout folder", "not found: %s" % ctx["keep_dir"])
    stop_if_failed(ck, label)
    print("\nScoring the %d runs of %s (reading each kept .qout once) ..." % (len(res), label))
    recs, sim_cfs, bad = ss.score_runs(res, ctx["keep_dir"], storms, real, runtime, ctx["synth5"], ck, False)
    show_gates(ck, label)
    if bad:
        print("\n  %d run(s) failed a gate or could not be read. First few:" % len(bad))
        for rid, why in list(bad.items())[:8]:
            print("    %s: %s" % (rid, why))
        sys.exit("\nSTOPPED: nothing was drawn. Fix the cause above.")
    if ck.failures():
        sys.exit("\nSTOPPED: a gate failed in %s. Nothing was drawn." % label)
    if len(recs) < 20:
        sys.exit("\nSTOPPED: only %d usable runs in %s; too few." % (len(recs), label))
    df_all = pd.DataFrame(recs)
    t = ss.arm_table(df_all, "filled", storms)
    infos, included = ss.rescale_storms(t)
    versions = ss.make_versions(included)
    ss.add_combined(t, versions)
    comb_rid, comb_score = None, np.nan
    if "ALL" in versions:
        pos = int(ss.top_positions(t, "ALL", 1)[0])
        comb_rid, comb_score = t.iloc[pos]["run_id"], float(t.iloc[pos][ss.ckey("ALL")])
    best = {}
    for s in storms:
        tg = s["tag"]
        raw = t[tg + "_kge_2012"].to_numpy(dtype=float)
        mf = t[tg + "_made_flow"].to_numpy(dtype=bool)
        ok = np.flatnonzero(mf & np.isfinite(raw))
        if len(ok):
            sub = t.iloc[ok].assign(_neg=-raw[ok]).sort_values(["_neg", "run_id"])
            best[tg] = int(sub.index[0])
        else:
            best[tg] = None
    return {"label": label, "design": design, "t": t, "sim_cfs": sim_cfs, "infos": infos, "included": included,
            "comb_rid": comb_rid, "comb_score": comb_score, "best": best, "n": len(t), "results_path": res_path,
            "gates": [{"status": s_, "name": n_, "text": x_} for s_, n_, x_ in ck.rows]}




# ------------------------------------------------------------------
# Fonts, caption, legends
# ------------------------------------------------------------------
def resolve_font(font_file):
    """Returns (family name, note). Liberation Sans if it can be found, else a look-alike, else DejaVu Sans."""
    import matplotlib.font_manager as fm
    if font_file:
        p = Path(font_file)
        if not p.exists():
            sys.exit("STOPPED: --font_file %s not found." % font_file)
        fm.fontManager.addfont(str(p))
        name = fm.FontProperties(fname=str(p)).get_name()
        return name, "font: %s (from %s)" % (name, p.name)
    for name in FONT_CANDIDATES:
        try:
            path = fm.findfont(fm.FontProperties(family=name), fallback_to_default=False)
        except Exception:
            continue
        return name, "font: %s (%s)" % (name, path)
    return "DejaVu Sans", ("*** font: Liberation Sans (and every look-alike) is not installed, so DejaVu Sans is used; its letters are wider "
                           "than Arial's. To install Liberation Sans: sudo apt-get update && sudo apt-get install -y fonts-liberation "
                           "then rm -rf ~/.cache/matplotlib and run this again ***")


def font_variant_note(font_family):
    """None if bold and italic of this family are separate font files (they are for Liberation Sans); otherwise a warning line."""
    import matplotlib.font_manager as fm
    try:
        reg = fm.findfont(fm.FontProperties(family=font_family), fallback_to_default=True)
        bold = fm.findfont(fm.FontProperties(family=font_family, weight="bold"), fallback_to_default=True)
        ital = fm.findfont(fm.FontProperties(family=font_family, style="italic"), fallback_to_default=True)
    except Exception:
        return None
    missing = [name for name, p in (("bold", bold), ("italic", ital)) if p == reg]
    if missing:
        return ("*** font: no separate %s file found for %s, so the %s text in the figure will look regular. "
                "Liberation Sans has all four files (regular, bold, italic, bold italic). ***" % (" or ".join(missing), font_family, " and ".join(missing)))
    return None


# ------------------------------------------------------------------
# Rich text: bold, italic, subscripts and superscripts in one wrapped paragraph.
# matplotlib cannot do this in a single text object, so the words are measured in the real font and each piece is placed by hand.
# Markup:  **bold**   *italic*   _{subscript}   ^{superscript}   and a new line starts a new paragraph.
# ------------------------------------------------------------------
def parse_markup(markup):
    """-> paragraphs; each is a list of words; each word is a list of atoms (text, bold, italic, script), script None/'sub'/'sup'."""
    out = []
    for para in markup.split("\n"):
        words, cur, buf = [], [], ""
        bold = italic = False
        i, n = 0, len(para)
        while i < n:
            ch = para[i]
            if ch == "*":
                if buf:
                    cur.append((buf, bold, italic, None))
                    buf = ""
                if para.startswith("**", i):
                    bold = not bold
                    i += 2
                else:
                    italic = not italic
                    i += 1
            elif para.startswith("_{", i) or para.startswith("^{", i):
                j = para.find("}", i)
                if j < 0:
                    raise ValueError("a { is not closed in the caption text: %r" % para)
                if buf:
                    cur.append((buf, bold, italic, None))
                    buf = ""
                cur.append((para[i + 2:j], bold, italic, "sub" if ch == "_" else "sup"))
                i = j + 1
            elif ch == " ":
                if buf:
                    cur.append((buf, bold, italic, None))
                    buf = ""
                if cur:
                    words.append(cur)
                    cur = []
                i += 1
            else:
                buf += ch
                i += 1
        if buf:
            cur.append((buf, bold, italic, None))
        if cur:
            words.append(cur)
        if bold or italic:
            raise ValueError("a ** or * is not closed in the caption text: %r" % para)
        if words:
            out.append(words)
    return out


def strip_markup(markup):
    """The same text without the markup: Ks, f, Kchan, mm^-1 ... (for the console, the CSV and the provenance file)."""
    import re
    s = re.sub(r"\^\{([^}]*)\}", r"^\1", markup)
    s = re.sub(r"_\{([^}]*)\}", r"\1", s)
    return s.replace("*", "").replace("\u00a0", " ")


def _fp(font_family, size_pt, bold, italic):
    from matplotlib.font_manager import FontProperties
    return FontProperties(family=font_family, size=size_pt, weight="bold" if bold else "normal", style="italic" if italic else "normal")


def layout_rich(markup, font_family, size_pt, max_w_in, dpi=300):
    """Word-wraps the markup to max_w_in inches, measured in the real font at the printed dpi.
    Returns {"lines": [{"atoms": [(dx_in, text, bold, italic, script)], "para_first": bool, "word_starts": [atom indices]}], "down": [first-baseline-to-baseline
    distances from the first baseline, in inches], "height_in", "h_lp", "d_lp", "pitch"}."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6.0, 1.0), dpi=dpi)
    try:
        r = fig.canvas.get_renderer()

        def width(text, bold, italic, scale):
            return r.get_text_width_height_descent(text, _fp(font_family, size_pt * scale, bold, italic), False)[0] / dpi

        space = width("a a", False, False, 1.0) - width("aa", False, False, 1.0)
        _, h_px, d_px = r.get_text_width_height_descent("lp", _fp(font_family, size_pt, False, False), False)
        h_lp, d_lp = h_px / dpi, d_px / dpi
        pitch = LINE_PITCH * h_lp
        lines = []
        for words in parse_markup(markup):
            line, x, started, starts = [], 0.0, False, []
            for word in words:
                pieces = [(a, width(a[0], a[1], a[2], SCRIPT_SCALE if a[3] else 1.0)) for a in word]
                ww = sum(w for _, w in pieces)
                if line and x + space + ww > max_w_in:
                    lines.append({"atoms": line, "para_first": not started, "word_starts": starts})
                    started = True
                    line, x, starts = [], 0.0, []
                if line:
                    x += space
                starts.append(len(line))
                for a, w in pieces:
                    line.append((x, a[0], a[1], a[2], a[3]))
                    x += w
            if line:
                lines.append({"atoms": line, "para_first": not started, "word_starts": starts})
    finally:
        plt.close(fig)
    down, y = [], 0.0
    for k, ln in enumerate(lines):
        if k > 0:
            y += pitch * (1.0 + (PARA_GAP if ln["para_first"] else 0.0))
        down.append(y)
    height = (down[-1] if down else 0.0) + h_lp
    return {"lines": lines, "down": down, "height_in": height, "h_lp": h_lp, "d_lp": d_lp, "pitch": pitch}


def layout_plain(layout):
    """The wrapped lines as plain text (for the console)."""
    out = []
    for ln in layout["lines"]:
        s = ""
        for k, (dx, text, bold, italic, script) in enumerate(ln["atoms"]):
            s += (" " if (k in ln["word_starts"] and k > 0) else "") + ("^" if script == "sup" else "") + text
        out.append(s)
    return out


def draw_rich(fig, layout, font_family, size_pt, x_in, first_baseline_in, color):
    """Draws a laid-out block with its first baseline at first_baseline_in (inches from the bottom of the figure)."""
    for ln, down in zip(layout["lines"], layout["down"]):
        y = first_baseline_in - down
        for (dx, text, bold, italic, script) in ln["atoms"]:
            dy = -SUB_DROP * size_pt / 72.0 if script == "sub" else (SUP_RISE * size_pt / 72.0 if script == "sup" else 0.0)
            fig.text(x_in + dx, y + dy, text, fontproperties=_fp(font_family, size_pt * (SCRIPT_SCALE if script else 1.0), bold, italic),
                     color=color, ha="left", va="baseline", transform=fig.dpi_scale_trans)


def count_fits(D):
    """How many runs of this sweep have the model's peak AND volume within BALLPARK_LO..HI x the real ones on all three big storms
    (and make at least FLOW_CFS there)? The same 'ballpark' rule as Series110_Compiled_Record.md 4.16. Returns (count, runs)."""
    t = D["t"]
    ok = np.ones(len(t), dtype=bool)
    for b in ss.BIG:
        tb = ss.tag_of(b)
        mf = t[tb + "_made_flow"].to_numpy(dtype=bool)
        pr = t[tb + "_peak_ratio"].to_numpy(dtype=float)
        vr = t[tb + "_vol_ratio"].to_numpy(dtype=float)
        ok &= mf & (pr >= BALLPARK_LO) & (pr <= BALLPARK_HI) & (vr >= BALLPARK_LO) & (vr <= BALLPARK_HI)
    return int(ok.sum()), len(t)


def lead_and_fits_sentence(fits, lead_override):
    """fits: [(sweep label, runs that fit all three, runs)]. Returns (lead-in markup, regular sentence markup, total fitting)."""
    total_ok = sum(k for _, k, _ in fits)
    total_n = sum(n for _, _, n in fits)
    where = [(SWEEP_NAME.get(lab, lab), k, n) for lab, k, n in fits]
    if total_ok == 0:
        lead = DEFAULT_LEAD_OK
        sent = ("None of the %d runs (%s) is within a factor of 2 of the real peak and volume on all three storms "
                "(a band we chose)." % (total_n, ", ".join("%d in the %s" % (n, nm) for nm, _, n in where)))
    else:
        lead = "tRIBS" + NBSP + "6.0: %d of %d runs are within a factor of 2 of the real peak and volume on all three real 2014 storms." % (total_ok, total_n)
        sent = "By sweep: %s (a band we chose)." % ", ".join("%d of %d in the %s" % (k, n, nm) for nm, k, n in where)
    if lead_override:
        lead = lead_override
    return lead, sent, total_ok


def build_caption(panels, rows_by_storm, sweeps, show_kge, lead, fits_sentence, all_panels, brief):
    """The caption as markup: two paragraphs (lead-in, the factor-of-2 count and the not-comparable sentence; the best runs), plus an
    'All panels:' paragraph only when all_panels is True."""
    p1 = "**%s** %s Not comparable with tRIBS" % (lead, fits_sentence) + NBSP + "5.3.1 results."
    entries = []
    for storm, sweep in panels:
        r = rows_by_storm[storm]
        if r.get("best_run_id") is None:
            entries.append("%s (%s): no run reaches %g cfs" % (storm.replace(" ", NBSP), RUN_COLOR_NAME[storm], ss.FLOW_CFS))
            continue
        e = "%s (%s) %s %s×, %s %s mm^{−1}, %s %s mm/hr" % (
            storm.replace(" ", NBSP), RUN_COLOR_NAME[storm], KS_M, ss.fval("Ks_mult", r["best_Ks_mult"]), F_M, ss.fval("f_RS_abs", r["best_f_RS_abs"]),
            KCHAN_M, ss.fval("channelconductivity_mmhr", r["best_Kchan_mmhr"]))
        e = e.replace(" mm^{", NBSP + "mm^{").replace(" mm/hr", NBSP + "mm/hr")
        if show_kge:
            e += ", KGE" + NBSP + "%.2f" % r["best_kge_2012"]
        entries.append(e)
    p2 = ("**Best run for each storm** (highest KGE among the runs that made at least %g cfs): %s." % (ss.FLOW_CFS, "; ".join(entries)))
    if not all_panels:
        return "\n".join([p1, p2])
    p3 = ["**All panels:** tRIBS" + NBSP + "6.0.0, real SMF gauge, 2014 (rain: SMF, SMPHQ gauges). Preliminary. "
          "KGE = Kling–Gupta efficiency (2012)."]
    if not brief:
        p3.append("%s= saturated hydraulic conductivity, as a multiple of baseline soil values; %s= conductivity decay with depth "
                  "(mm^{−1}); %s= channel conductivity (mm/hr)." % (KS_M + NBSP, F_M + NBSP, KCHAN_M + NBSP))
    p3.append("Routing held fixed: %s." % ROUTING_TEXT)
    return "\n".join([p1, p2, " ".join(p3)])


def legend_handles(panels, any_cross, show_combined):
    from matplotlib.lines import Line2D
    h1 = [Line2D([0], [0], color=ss.INK, lw=LW_REAL, label="real gauge, 2014")]
    for storm, _ in panels:
        h1.append(Line2D([0], [0], color=RUN_COLOR[storm], lw=LW_OWN, solid_capstyle="round", label="best run, %s" % storm))
    if show_combined:
        h1.append(Line2D([0], [0], color=C_COMBINED, lw=LW_COMBINED, ls=(0, (3, 2)), label="best run, storms combined"))
    h2 = []
    if any_cross:
        h2 = [Line2D([0], [0], color=ss.INK2, lw=LW_OWN, solid_capstyle="round", label="solid: the storm the run was chosen for"),
              Line2D([0], [0], color=ss.INK2, lw=LW_CROSS, ls=DOTTED, dash_capstyle="round", label="dotted: the same run in another storm")]
    return h1, h2


def _legend(fig, handles, ncol, anchor_xy):
    return fig.legend(handles=handles, loc="lower center", bbox_to_anchor=anchor_xy, bbox_transform=fig.transFigure, ncol=ncol,
                      frameon=False, labelcolor=TEXT_INK, fontsize=POSTER_PT, handlelength=2.6, columnspacing=2.4, borderpad=0.15)


def measure_bottom(plt, rc, width_in, h1, h2, caption, font_family):
    """Sizes (in inches) of the two legend rows and the wrapped caption, so the figure height can be worked out."""
    with plt.rc_context(rc):
        fig = plt.figure(figsize=(width_in, 6.0), dpi=100)
        renderer = fig.canvas.get_renderer()
        out = {}
        for key, handles in (("1", h1), ("2", h2)):
            out["ncol" + key], out["h" + key] = 1, 0.0
            if not handles:
                continue
            for ncol in range(len(handles), 0, -1):
                leg = _legend(fig, handles, ncol, (0.5, 0.0))
                ext = leg.get_window_extent(renderer)
                w_in, h_in = ext.width / fig.dpi, ext.height / fig.dpi
                leg.remove()
                if w_in <= width_in - 2 * SIDE_MARGIN_IN or ncol == 1:
                    out["ncol" + key], out["h" + key] = ncol, h_in
                    break
        plt.close(fig)
    out["layout"] = layout_rich(caption, font_family, CAPTION_PT, width_in - 2 * SIDE_MARGIN_IN)
    out["h_cap"] = out["layout"]["height_in"]
    return out


# ------------------------------------------------------------------
# The figure
# ------------------------------------------------------------------
def zoom_limits(x, curves, floor_cfs, pad_min):
    """Time limits that cover every curve while it is at or above floor_cfs, plus a margin on each side."""
    on = np.zeros(len(x), dtype=bool)
    for y in curves:
        if y is not None:
            on |= np.nan_to_num(np.asarray(y, dtype=float), nan=0.0) >= floor_cfs
    if not on.any():
        return x[0], x[-1]
    idx = np.flatnonzero(on)
    pad = pd.Timedelta(minutes=pad_min)
    return max(x[0], x[idx[0]] - pad), min(x[-1], x[idx[-1]] + pad)


def draw_figure(path, panels, sweeps, storms, real, runs, cross, show_combined, caption, font_family,
                width_in, height_in, panel_block_in, zoom, pad_min, zoom_frac):
    """runs: {storm: (sweep, run_id) or None}; cross: {storm: [other storms whose run is drawn dotted here]}.
    Returns the total figure height in inches."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    byl = {s["label"]: s for s in storms}
    rc = {"font.size": POSTER_PT, "axes.titlesize": TITLE_PT, "axes.labelsize": POSTER_PT,
          "xtick.labelsize": POSTER_PT, "ytick.labelsize": POSTER_PT, "legend.fontsize": POSTER_PT,
          "font.family": [font_family], "text.color": TEXT_INK, "axes.labelcolor": TEXT_INK}
    any_cross = any(cross.get(st) for st, _ in panels)
    h1, h2 = legend_handles(panels, any_cross, show_combined)
    m = measure_bottom(plt, rc, width_in, h1, h2, caption, font_family)
    sub_layouts = []                              # the one-line italic/subscript subtitle under each panel title
    for storm, sweep in panels:
        b = sweeps[sweep]["design"]["box"]["Ks_mult"]
        sub_layouts.append(layout_rich("%s, %g–%g" % (SWEEP_NAME.get(sweep, sweep), b[0], b[1]), font_family, POSTER_PT, width_in))
    margin, gap = 0.12, 0.06                      # inches: below the caption; between caption and legends
    y_cap = margin
    y_l2 = y_cap + m["h_cap"] + gap
    y_l1 = y_l2 + (m["h2"] if h2 else 0.0)
    bottom_in = y_l1 + m["h1"] + 0.12             # a little air between the tick labels and the legend
    if height_in is None:
        total_in = panel_block_in + bottom_in
    else:
        total_in = height_in
        if total_in - bottom_in < 3.0:
            sys.exit("STOPPED: --height_in %.1f leaves only %.1f in for the panels (the legend and caption need %.1f in). "
                     "Use at least %.1f." % (height_in, total_in - bottom_in, bottom_in, bottom_in + 3.0))
    with plt.rc_context(rc):
        fig, axes = plt.subplots(1, len(panels), figsize=(width_in, total_in))
        axes = np.atleast_1d(axes)
        fig.patch.set_facecolor(ss.SURFACE)
        for ax, (storm, sweep) in zip(axes, panels):
            s = byl[storm]
            tg, x = s["tag"], s["grid"]
            D = sweeps[sweep]
            ss._style(ax)
            ax.tick_params(labelsize=POSTER_PT, length=6, width=1.2, labelcolor=TEXT_INK)
            ax.grid(True, color=ss.GRID, lw=1.0)
            real_cfs = real[tg]["filled"].to_numpy() / ss.CFS_TO_CMS
            drawn = [real_cfs]
            # dotted: other storms' runs, shown in this storm
            todo = [o for o in cross.get(storm, []) if runs.get(o) is not None]
            for k, other in enumerate(todo):
                sw2, rid2 = runs[other]
                y = sweeps[sw2]["sim_cfs"][rid2][tg]
                # shift each dotted line's dots by a share of the dot spacing, so lines that coincide (two runs that both
                # make no flow lie on top of each other along zero) show alternating dots instead of hiding one another
                ls = (k * DOTTED_PERIOD / max(len(todo), 1), DOTTED[1])
                ax.plot(x, y, color=RUN_COLOR[other], lw=LW_CROSS, ls=ls, dash_capstyle="round", zorder=3)
                drawn.append(y)
            if show_combined and D["comb_rid"] is not None:
                y = D["sim_cfs"][D["comb_rid"]][tg]
                ax.plot(x, y, color=C_COMBINED, lw=LW_COMBINED, ls=(0, (3, 2)), zorder=4)
                drawn.append(y)
            own = runs.get(storm)
            if own is not None:
                y = sweeps[own[0]]["sim_cfs"][own[1]][tg]
                ax.plot(x, y, color=RUN_COLOR[storm], lw=LW_OWN, solid_capstyle="round", zorder=4.5)
                drawn.append(y)
            ax.plot(x, real_cfs, color=ss.INK, lw=LW_REAL, zorder=5)
            if own is None:
                ax.text(0.98, 0.95, "no run reaches\n%g cfs" % ss.FLOW_CFS, transform=ax.transAxes, ha="right", va="top",
                        fontsize=POSTER_PT, fontweight="bold", color=TEXT_INK)
            if zoom:
                floor = max(ss.FLOW_CFS, zoom_frac * real[tg]["real_peak_cfs"])
                lo, hi = zoom_limits(x, drawn, floor, pad_min)
            else:
                lo, hi = x[0], x[-1]
            ax.set_xlim(lo, hi)
            win = (x >= lo) & (x <= hi)
            top = max(float(np.nanmax(np.asarray(y, dtype=float)[win])) for y in drawn)
            ax.set_ylim(0, top * 1.08)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: "{:,.0f}".format(v)))
            step = hour_step((hi - lo).total_seconds() / 3600.0)
            ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, step)))
            ax.xaxis.set_major_formatter(FuncFormatter(clock_label))
            ax.set_title("%s · real peak %.0f cfs" % (storm, real[tg]["real_peak_cfs"]), fontsize=TITLE_PT, color=HEAD_NAVY,
                         fontweight="bold", loc="left", pad=POSTER_PT * 1.25 + 14)
        axes[0].set_ylabel("flow (cfs)", fontsize=POSTER_PT, color=TEXT_INK)
        _legend(fig, h1, m["ncol1"], (0.5, y_l1 / total_in))
        if h2:
            _legend(fig, h2, m["ncol2"], (0.5, y_l2 / total_in))
        cap = m["layout"]                          # bold lead-in + regular text, laid out by hand (see layout_rich)
        draw_rich(fig, cap, font_family, CAPTION_PT, SIDE_MARGIN_IN, y_cap + cap["height_in"] - (cap["h_lp"] - cap["d_lp"]), TEXT_INK)
        fig.tight_layout(rect=(0, bottom_in / total_in, 1, 1), pad=0.4, w_pad=2.0)
        for ax, lay in zip(axes, sub_layouts):     # subtitles go on after the layout is fixed: 6 pt above each panel, left-aligned
            pos = ax.get_position()
            draw_rich(fig, lay, font_family, POSTER_PT, pos.x0 * width_in, pos.y1 * total_in + 6.0 / 72.0 + lay["d_lp"], TEXT_INK)
        fig.savefig(path, dpi=300, facecolor=ss.SURFACE)
        plt.close(fig)
    return total_in, bottom_in


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Series 110: the poster hydrograph figure (no tRIBS runs).")
    ap.add_argument("--panels", default=DEFAULT_PANELS,
                    help="which sweep each storm's panel comes from (default '%s')" % DEFAULT_PANELS)
    ap.add_argument("--width_in", type=float, default=22.5, help="figure width in inches (default 22.5)")
    ap.add_argument("--height_in", type=float, default=None,
                    help="total figure height in inches. Default: the panels' height (--panel_block_in) plus whatever the legend and caption need")
    ap.add_argument("--panel_block_in", type=float, default=5.7,
                    help="height in inches of the panels with their titles and tick labels (default 5.7)")
    ap.add_argument("--no_dotted_on_27sep", action="store_true", help="leave the blue and green runs (dotted, along zero) off the 27 Sep panel")
    ap.add_argument("--show_combined", action="store_true", help="also draw the best run of the combined score (purple dashed)")
    ap.add_argument("--caption_kge", action="store_true", help="put each run's KGE in the caption")
    ap.add_argument("--all_panels", action="store_true",
                    help="add the 'All panels:' note under the caption (build label, Preliminary, KGE definition, symbol definitions, fixed routing)")
    ap.add_argument("--brief", action="store_true", help="with --all_panels: leave out the definitions of Ks, f and Kchan")
    ap.add_argument("--lead", default=None,
                    help="your own bold lead-in sentence for the caption (markup: *italic*, _{sub}, ^{sup}); default: the sentence the "
                         "factor-of-2 count supports")
    ap.add_argument("--font_file", default=None, help="path to a .ttf font to use for all text instead of Liberation Sans")
    ap.add_argument("--no_zoom", action="store_true", help="show each full storm window instead of zooming in on the event")
    ap.add_argument("--zoom_pad_min", type=float, default=30.0, help="minutes of margin around the event when zoomed (default 30)")
    ap.add_argument("--zoom_frac", type=float, default=0.02,
                    help="zoom to where a curve is above this share of the real peak (and at least 1 cfs) (default 0.02)")
    ap.add_argument("--check_only", action="store_true", help="run the gates and stop; writes nothing")
    ap.add_argument("--skip_gate2", action="store_true", help="skip the stored-metric check (not advised)")
    ap.add_argument("--windows", type=Path, default=None)
    ap.add_argument("--gauge_xlsx", type=Path, default=None)
    ap.add_argument("--out_dir", type=Path, default=None)
    args = ap.parse_args()
    if not (6.0 <= args.width_in <= 60.0 and (args.height_in is None or 2.0 <= args.height_in <= 30.0) and 3.0 <= args.panel_block_in <= 20.0):
        ap.error("--width_in must be 6 to 60, --height_in 2 to 30 and --panel_block_in 3 to 20.")
    if args.lead:
        try:
            parse_markup("**%s**" % args.lead)
        except ValueError as exc:
            sys.exit("STOPPED: --lead is not valid: %s. Use *text* for italic, _{text} for a subscript, ^{text} for a superscript, and do not "
                     "put ** in it (the whole lead-in is bold already)." % exc)
    panels = parse_panels(args.panels)
    cross = {k: list(v) for k, v in CROSS_DEFAULT.items()}
    if args.no_dotted_on_27sep:
        cross["27 Sep"] = []

    script_dir = Path.cwd()
    if not (script_dir / "score_storms_110.py").exists():
        sys.exit("STOPPED: run this from the lab/ folder (score_storms_110.py is not in %s)." % script_dir)
    project_root = script_dir.parent
    calib_dir = project_root / "calibration_work"
    summary_dir = calib_dir / "03_comparisons" / "summary_tables"
    synth_dir = calib_dir / "synth_truth"
    keep_dir = calib_dir / "kept_qout_110"
    csv_dir = calib_dir / "03_comparisons" / "csv_exports"
    win_path = args.windows or (script_dir / ss.WINDOWS_NAME)
    xlsx_path = args.gauge_xlsx or project_root.joinpath(*ss.GAUGE_REL)
    out_dir = args.out_dir or (summary_dir / OUT_SUBFOLDER)

    print("\n" + "=" * 100)
    print("Series 110 -- poster hydrograph figure (version %s; scorer functions from score_storms_110.py version %s)" % (VERSION, ss.VERSION))
    print("panels: %s   python %s  pandas %s  numpy %s" % ("; ".join("%s from %s" % p for p in panels), sys.version.split()[0], pd.__version__, np.__version__))
    print("=" * 100)
    if ".".join(pd.__version__.split(".")[:2]) != ss.EXPECT_PANDAS_MINOR:
        print("  *** NOTE: pandas here is %s; the runs were scored with pandas %s.x. Gate 2 will tell. ***" % (pd.__version__, ss.EXPECT_PANDAS_MINOR))
    font_family, font_note = resolve_font(args.font_file)
    print("  " + font_note)
    variant_note = font_variant_note(font_family)
    if variant_note:
        print("  " + variant_note)
    out_res = out_dir.resolve()
    for guarded in (synth_dir, csv_dir, keep_dir):
        g = guarded.resolve()
        if out_res == g or g in out_res.parents or out_res in g.parents:
            sys.exit("Output folder %s overlaps an input folder (%s); refusing." % (out_dir, guarded))

    ck = ss.Checks()
    storms, runtime = ss.load_windows(win_path, ck)
    if not xlsx_path.exists():
        ck.add("FAIL", "gauge workbook", "not found: %s" % xlsx_path)
    stop_if_failed(ck, "windows and gauge")
    print("\nReading the real gauge ...")
    real = build_real(storms, xlsx_path, summary_dir, ck)
    synth5 = None
    if not args.skip_gate2:
        qs = list(synth_dir.glob("*.qout")) if synth_dir.exists() else []
        if len(qs) != 1:
            ck.add("FAIL", "GATE 2", "expected exactly one *.qout in %s (the synthetic truth), found %d" % (synth_dir, len(qs)))
        else:
            synth5 = ss.read_qout_5min(qs[0])
            ck.add("NOTE", "GATE 2", "synthetic truth (read-only): %s" % qs[0].name)
    else:
        ck.add("WARN", "GATE 2", "skipped by --skip_gate2")
    stop_if_failed(ck, "windows, gauge, synthetic truth")
    ctx = {"storms": storms, "runtime": runtime, "real": real, "summary_dir": summary_dir, "keep_dir": keep_dir, "synth5": synth5}

    sweeps = {}
    for label in sorted({sw for _, sw in panels}):
        sweeps[label] = load_sweep(label, ctx)
    show_gates(ck, "windows, gauge, synthetic truth")

    # ---- which run is drawn for each storm ---------------------------------------------------------
    byl = {s["label"]: s for s in storms}
    rows, rows_by_storm, notes, runs = [], {}, [], {}
    print("\n" + "=" * 100)
    print("WHAT IS DRAWN   (one run per storm; its colour in every panel: %s)" % ", ".join(
        "%s = %s" % (st, RUN_COLOR_NAME[st]) for st, _ in panels))
    print("=" * 100)
    for storm, sw in panels:
        tg = byl[storm]["tag"]
        D = sweeps[sw]
        t = D["t"]
        n_flow = int(t[tg + "_made_flow"].sum())
        bpos = D["best"][tg]
        info = [i for i in D["infos"] if i["storm"] == storm][0]
        row = {"storm": storm, "colour": RUN_COLOR_NAME[storm], "sweep": sw, "sweep_name": strip_markup(SWEEP_NAME.get(sw, sw)),
               "ks_box_lo": D["design"]["box"]["Ks_mult"][0], "ks_box_hi": D["design"]["box"]["Ks_mult"][1],
               "runs_in_sweep": D["n"], "runs_making_flow_on_storm": n_flow,
               "real_peak_cfs": real[tg]["real_peak_cfs"], "real_volume_acft": real[tg]["real_vol_acft"],
               "storm_in_combined_mean": storm in D["included"], "storm_status": info["status"], "storm_flag": info["flag"],
               "routing_pinned_in_design_file": str(D["design"]["routing"])}
        if bpos is not None:
            r = t.iloc[bpos]
            runs[storm] = (sw, r["run_id"])
            row.update({"best_run_id": r["run_id"], "best_Ks_mult": r["Ks_mult"], "best_f_RS_abs": r["f_RS_abs"],
                        "best_Kchan_mmhr": r["channelconductivity_mmhr"], "best_kge_2012": r[tg + "_kge_2012"],
                        "best_peak_cfs": r[tg + "_sim_peak_cfs"], "best_peak_ratio": r[tg + "_peak_ratio"],
                        "best_volume_ratio": r[tg + "_vol_ratio"], "best_peak_dt_min": r[tg + "_peak_dt_min"]})
            for b in ss.BIG:           # how this run does on every big storm
                tb = ss.tag_of(b)
                row["run_on_%s_made_flow" % tb] = bool(r[tb + "_made_flow"])
                row["run_on_%s_peak_cfs" % tb] = r[tb + "_sim_peak_cfs"]
                row["run_on_%s_kge_2012" % tb] = r[tb + "_kge_2012"]
            print("  %-7s %-6s %-4s run %s   Ks %s  f %s  Kchan %s   (%d of %d runs make flow)"
                  % (storm, "(%s)" % RUN_COLOR_NAME[storm], sw, r["run_id"], ss.fval("Ks_mult", r["Ks_mult"]),
                     ss.fval("f_RS_abs", r["f_RS_abs"]), ss.fval("channelconductivity_mmhr", r["channelconductivity_mmhr"]), n_flow, D["n"]))
            if r[tg + "_kge_2012"] < ss.WEAK_KGE:
                notes.append("%s (%s): the best run's KGE_2012 is %.2f, under the scorer's WEAK line of %g; the figure shows the best this sweep can do, and it is a poor fit."
                             % (storm, sw, r[tg + "_kge_2012"], ss.WEAK_KGE))
        else:
            runs[storm] = None
            print("  %-7s %-6s %-4s NO run makes flow (%d of %d): the panel shows only the real gauge." % (storm, "(%s)" % RUN_COLOR_NAME[storm], sw, n_flow, D["n"]))
            notes.append("%s (%s): no run makes flow, so the panel shows the real gauge and a note." % (storm, sw))
        if D["comb_rid"] is not None:      # kept in the CSV for reference; drawn only with --show_combined
            c = t[t["run_id"] == D["comb_rid"]].iloc[0]
            row.update({"combined_run_id": D["comb_rid"], "combined_storms": "+".join(D["included"]), "combined_score": D["comb_score"],
                        "combined_Ks_mult": c["Ks_mult"], "combined_f_RS_abs": c["f_RS_abs"], "combined_Kchan_mmhr": c["channelconductivity_mmhr"],
                        "combined_kge_2012_this_storm": c[tg + "_kge_2012"], "combined_peak_cfs_this_storm": c[tg + "_sim_peak_cfs"]})
        elif args.show_combined:
            notes.append("%s: %s has no combined score, so no combined line in this panel." % (storm, sw))
        rows.append(row)
        rows_by_storm[storm] = row
    ids = [v for v in runs.values() if v is not None]
    if len(ids) != len(set(ids)):
        notes.append("the same run is the best for two storms, so two colours lie exactly on top of each other.")

    print("\nHOW EACH DRAWN RUN DOES ON EACH STORM   (model peak in cfs, KGE_2012; real peaks: %s)" % ", ".join(
        "%s %.0f" % (b, real[ss.tag_of(b)]["real_peak_cfs"]) for b in ss.BIG))
    for storm, _ in panels:
        r = rows_by_storm[storm]
        if r.get("best_run_id") is None:
            continue
        cells = []
        for b in ss.BIG:
            tb = ss.tag_of(b)
            cells.append("%s: %s" % (b, ("%.0f cfs, KGE %s" % (r["run_on_%s_peak_cfs" % tb], ss.fnum(r["run_on_%s_kge_2012" % tb], ".2f")))
                                       if r["run_on_%s_made_flow" % tb] else "no flow"))
        print("  run chosen for %-7s %-9s  %s" % (storm, r["best_run_id"], "   ".join(cells)))

    # ---- how many runs fit all three storms (your factor-of-2 band) -----------------------------
    fits = []
    for lab in dict.fromkeys(sw for _, sw in panels):          # each sweep once, in panel order
        k, n = count_fits(sweeps[lab])
        fits.append((lab, k, n))
    print("\nRUNS WITHIN A FACTOR OF 2 (%g to %g x) OF THE REAL PEAK AND VOLUME ON ALL THREE BIG STORMS (and at least %g cfs on each):"
          % (BALLPARK_LO, BALLPARK_HI, ss.FLOW_CFS))
    for lab, k, n in fits:
        print("  %-4s %d of %d runs" % (lab, k, n))
        if lab == "ms2" and n == 60 and k != 0:
            notes.append("ms2 shows %d of 60 runs within a factor of 2 on all three storms, but Series110_Compiled_Record.md 4.16 says 0 of 60. "
                         "Check before using the caption." % k)
    lead, fits_sentence, total_fit = lead_and_fits_sentence(fits, args.lead)
    if args.lead and total_fit > 0:
        notes.append("--lead is used as given, but %d run(s) are within a factor of 2 on all three storms; check that your lead-in still says what the data say." % total_fit)
    if total_fit > 0 and not args.lead:
        notes.append("%d run(s) are within a factor of 2 on all three storms, so the lead-in says that instead of 'no single ... set fits'." % total_fit)

    if args.brief and not args.all_panels:
        notes.append("--brief only changes the 'All panels:' note; add --all_panels to see it.")
    caption = build_caption(panels, rows_by_storm, sweeps, args.caption_kge, lead, fits_sentence, args.all_panels, args.brief)
    print("\nROUTING (in the caption only with --all_panels): \"%s\". The design files say: %s" % (
        strip_markup(ROUTING_TEXT), "; ".join("%s: %s" % (k, v["design"]["routing"]) for k, v in sweeps.items())))

    if args.check_only:
        lay = layout_rich(caption, font_family, CAPTION_PT, args.width_in - 2 * SIDE_MARGIN_IN)
        print("\nCAPTION (as it will be wrapped at %d pt in %s; %d lines, %.2f in tall):" % (CAPTION_PT, font_family, len(lay["lines"]), lay["height_in"]))
        for line in layout_plain(lay):
            print("  " + line)
        print("\n--check_only: gates finished; nothing was written.")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    fig_path = out_dir / "fig_poster_hydrographs_110.png"
    tmp = fig_path.with_name("tmp_" + fig_path.name)
    total_in, bottom_in = draw_figure(tmp, panels, sweeps, storms, real, runs, cross, args.show_combined, caption, font_family,
                                      args.width_in, args.height_in, args.panel_block_in, not args.no_zoom,
                                      args.zoom_pad_min, args.zoom_frac)
    os.replace(str(tmp), str(fig_path))
    ss.write_csv_atomic(pd.DataFrame(rows), out_dir / "poster_best_runs_110.csv")
    prov = {"script": Path(__file__).name, "version": VERSION, "script_md5": ss.md5_of(Path(__file__)),
            "scorer": {"file": "score_storms_110.py", "version": ss.VERSION, "md5": ss.md5_of(Path(ss.__file__))},
            "created_local": datetime.now().isoformat(timespec="seconds"), "python": sys.version.split()[0],
            "pandas": pd.__version__, "numpy": np.__version__, "panels": [{"storm": a, "sweep": b} for a, b in panels],
            "runs_drawn": {k: (None if v is None else {"sweep": v[0], "run_id": v[1]}) for k, v in runs.items()},
            "dotted": cross, "colours": RUN_COLOR, "show_combined": bool(args.show_combined), "font": font_family,
            "text_style": {"text_pt": POSTER_PT, "title_pt": TITLE_PT, "caption_pt": CAPTION_PT, "ink": TEXT_INK, "title_colour": HEAD_NAVY,
                           "line_pitch": LINE_PITCH, "bold_italic_files_found": variant_note is None},
            "caption_markup": caption, "caption": strip_markup(caption), "caption_lead": strip_markup(lead),
            "runs_within_factor_2_on_all_three_storms": {lab: {"runs": n, "fit": k} for lab, k, n in fits},
            "ballpark_band": [BALLPARK_LO, BALLPARK_HI],
            "figure": {"width_in": args.width_in, "height_in": round(total_in, 2), "legend_and_caption_block_in": round(bottom_in, 2),
                       "dpi": 300, "zoom": not args.no_zoom, "zoom_pad_min": args.zoom_pad_min, "zoom_frac": args.zoom_frac},
            "inputs": {"windows": {"file": win_path.name, "md5": ss.md5_of(win_path)},
                       "gauge": {"file": xlsx_path.name, "md5": ss.md5_of(xlsx_path)},
                       **{"sweep_" + k: {"design": v["design"]["path"].name, "results": v["results_path"].name,
                                          "results_md5": ss.md5_of(v["results_path"]), "runs": v["n"],
                                          "storms_in_combined_mean": v["included"]} for k, v in sweeps.items()}},
            "gates": {k: v["gates"] for k, v in sweeps.items()}, "notes": notes}
    ss.write_json_atomic(out_dir / "PROVENANCE_poster_110.json", prov)

    print("\nFIGURE SIZE: %.1f in wide x %.1f in tall (panels %.1f in, legend and caption %.1f in). Change the height with --height_in."
          % (args.width_in, total_in, total_in - bottom_in, bottom_in))
    print("\nCAVEATS")
    print("  - Each best run is picked by hindsight on one storm: it shows what the model can do there. The dotted lines show what the same parameters do elsewhere.")
    print("  - The 27 Sep run comes from the low Ks sweep (ms2), the other two from the original Ks sweep (ms1).")
    for c in ss.CAVEATS[:3]:
        print("  - " + c)
    for n in notes:
        print("  * NOTE: " + n)
    print("\nSaved to: %s" % out_dir)
    print("  fig_poster_hydrographs_110.png   poster_best_runs_110.csv   PROVENANCE_poster_110.json")


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        main()
