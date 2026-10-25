"""
make_poster_scoring_figure_110.py   (version 1)
===============================================
Poster panel "How runs are scored": ONE picture of an observed and a simulated hydrograph, with three coloured
callouts (timing, volume, shape) that match the three terms of the KGE equation under it. No tRIBS runs, no data
files: it only draws. Size is 10.5 x 3.5 in at 300 dpi (3150 x 1050 pixels) so it drops straight into the panel.

IT IS A SCHEMATIC, NOT DATA
---------------------------
Both curves are made from a simple formula (a skewed hump), not from a tRIBS run or the gauge. The simulated curve is
drawn late, flatter and broader than the observed one, and its AREA is set equal to the observed area (see
SIM_VOLUME_RATIO), so the picture shows that a run can match volume and still lose points for timing and shape.
The picture says "Schematic, not real data" in its corner unless you add --no_note.

WHAT EACH CALLOUT MEANS
-----------------------
  timing  (blue)    a double arrow between the observed peak and the simulated peak: the simulated peak arrives late.
                    This is what the KGE correlation term (KGE r) responds to.
  volume  (green)   the area under the observed curve is shaded: total flow. This is the KGE beta term
                    (mean simulated flow / mean observed flow).
  shape   (orange)  the part of the observed curve above half its peak has an orange halo: how peaked or spread out the
                    hydrograph is. This is a visual cue for the KGE gamma term, which is the ratio of the coefficients
                    of variation (CV simulated / CV observed). The halo is NOT a literal measurement of gamma.

COLOURS
-------
The three callout colours are the same three hex values that make_poster_hydrographs_110.py uses for the storms
(RUN_COLOR: 12 Aug #36a2ff, 8 Sep #4cd05a, 27 Sep #ff9a3d). They are bright, so the words on them are dark ink
(#0F1825) and never the colour itself. To change a colour, edit CALLOUT_COLOR below (one line each); nothing else
needs to change. Ink, grey and background are the same as in the other poster figures (#0b0b0b, #52514e, #fcfcfb).

TEXT SIZE
---------
Every piece of text is 24 pt at printed size (the poster floor). The script checks this after drawing and stops if
anything is smaller. Font: Liberation Sans (same letter widths as Arial), found the same way as in
make_poster_hydrographs_110.py; --font_file lets you point at a particular .ttf.

USAGE (from any folder; matplotlib and numpy are enough)
--------------------------------------------------------
    python make_poster_scoring_figure_110.py                    # the figure, 10.5 x 3.5 in
    python make_poster_scoring_figure_110.py --svg              # also write a vector .svg (text turned into outlines)
    python make_poster_scoring_figure_110.py --transparent      # no background, so it takes the panel's colour
    python make_poster_scoring_figure_110.py --no_note          # leave out "Schematic, not real data"
    python make_poster_scoring_figure_110.py --height_in 3.3    # other height (3.2 to 4.5 works; width 8 to 14)
    python make_poster_scoring_figure_110.py --font_file Arial.ttf

WHAT YOU GET
------------
    fig_poster_scoring_110.png   300 dpi, exactly width_in x height_in
    fig_poster_scoring_110.svg   only with --svg
Saved in calibration_work/03_comparisons/summary_tables/poster_110/ when that folder's parent exists (run from lab/),
otherwise in ./poster_110/. --out_dir picks another folder.
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

VERSION = "1"

# ------------------------------------------------------------------
# Look of the figure (change these before a run if you want a different look)
# ------------------------------------------------------------------
POSTER_PT = 24                  # ALL text, points at printed size (the poster floor)
TEXT_INK = "#0F1825"            # words (same as make_poster_hydrographs_110.py)
OBS_INK = "#0b0b0b"             # observed curve (the real-gauge ink of the other poster figures)
SIM_INK = "#52514e"             # simulated curve, axes, the note
SURFACE = "#fcfcfb"             # background (same as the other poster figures)
# the same hex values as RUN_COLOR in make_poster_hydrographs_110.py ("12 Aug", "8 Sep", "27 Sep")
CALLOUT_COLOR = {"timing": "#36a2ff", "volume": "#4cd05a", "shape": "#ff9a3d"}
FONT_CANDIDATES = ["Liberation Sans", "Arial", "Arimo", "Helvetica", "Nimbus Sans"]   # tried in this order

X_MAX = 10.0                    # the time axis runs 0 to 10 (no units: it is a schematic)
OBS = dict(t_peak=3.0, q_peak=1.0, k=3.2)    # observed hump: peak at 3.0, height 1, skew k
SIM = dict(t_peak=4.0, k=2.8)                # simulated hump: peak 1.0 later, flatter and broader; its height is set by SIM_VOLUME_RATIO
SIM_VOLUME_RATIO = 1.00         # simulated AREA / observed AREA (1.00 = same total volume, so only timing and shape differ)
HALF = 0.5                      # the "shape" halo covers the observed curve where it is above this share of its peak
TOP_ROOM_IN = 0.68              # inches kept free above the observed peak for the timing arrow and its tag
LW_CURVE, LW_HALO, LW_ARROW, LW_DOT = 4.2, 13.0, 3.4, 2.8     # line widths in points
MARGINS_IN = dict(left=0.95, bottom=0.62, right=0.15, top=0.12)   # room for the y label, x label, and a little air
SHAPE_TAG_LEFT_IN = 0.34        # the "shape" tag's text starts this far (inches) right of the y axis, so its box clears the axis line
NOTE_TEXT = "Schematic, not real data"
OUT_SUBFOLDER = "poster_110"


# ------------------------------------------------------------------
# The curves
# ------------------------------------------------------------------
def hump(t, t_peak, q_peak, k):
    """A skewed hump: rises quickly, recedes slowly; its maximum is q_peak at t = t_peak. Larger k = narrower and more peaked."""
    x = np.clip(t / t_peak, 0.0, None)
    return q_peak * (x * np.exp(1.0 - x)) ** k


def area(t, y):
    """Area under a curve (trapezoid rule), kept by hand so it works in every numpy version."""
    return float(np.sum((y[1:] + y[:-1]) * np.diff(t)) / 2.0)


def half_peak_crossings(t, y, frac):
    """Times where the curve crosses frac x its peak: (rising side, falling side)."""
    i = int(np.argmax(y))
    level = frac * float(y[i])
    t_up = float(np.interp(level, y[:i + 1], t[:i + 1]))                 # y rises along t[:i+1]
    t_dn = float(np.interp(level, y[i:][::-1], t[i:][::-1]))             # y falls along t[i:], so reverse both to make it increase
    return t_up, t_dn


# ------------------------------------------------------------------
# Fonts and checks
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
    return "DejaVu Sans", ("*** font: Liberation Sans (and every look-alike) is not installed, so DejaVu Sans is used; its letters are "
                           "wider than Arial's. To install Liberation Sans: sudo apt-get update && sudo apt-get install -y fonts-liberation "
                           "then rm -rf ~/.cache/matplotlib and run this again ***")


def luminance(hex_color):
    r, g, b = [int(hex_color[i:i + 2], 16) / 255.0 for i in (1, 3, 5)]
    lin = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)


def contrast(a, b):
    """WCAG contrast ratio of two hex colours (4.5 or more is the usual bar for text; poster text this large needs 3)."""
    hi, lo = sorted((luminance(a), luminance(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


# ------------------------------------------------------------------
# The figure
# ------------------------------------------------------------------
def draw_figure(png_path, svg_path, width_in, height_in, dpi, font_family, transparent, show_note):
    """Draws and saves the figure. Returns a dict of numbers for the printout."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.text import Text

    t = np.linspace(0.0, X_MAX, 1001)
    obs = hump(t, **OBS)
    sim_unit = hump(t, SIM["t_peak"], 1.0, SIM["k"])
    sim = sim_unit * (SIM_VOLUME_RATIO * area(t, obs) / area(t, sim_unit))     # same area as the observed curve (times the ratio)
    i_o, i_s = int(np.argmax(obs)), int(np.argmax(sim))
    x_o, x_s, q_o, q_s = float(t[i_o]), float(t[i_s]), float(obs[i_o]), float(sim[i_s])
    x_up, x_dn = half_peak_crossings(t, obs, HALF)

    rc = {"font.family": [font_family], "font.size": POSTER_PT, "axes.labelsize": POSTER_PT, "legend.fontsize": POSTER_PT,
          "text.color": TEXT_INK, "axes.labelcolor": TEXT_INK}
    with plt.rc_context(rc):
        fig = plt.figure(figsize=(width_in, height_in))
        fig.patch.set_facecolor(SURFACE)
        m = MARGINS_IN
        ax_w, ax_h = width_in - m["left"] - m["right"], height_in - m["bottom"] - m["top"]
        ax = fig.add_axes([m["left"] / width_in, m["bottom"] / height_in, ax_w / width_in, ax_h / height_in])
        y_per_in = q_o / (ax_h - TOP_ROOM_IN)            # data units per inch up the axis, so the observed peak sits TOP_ROOM_IN below the top
        x_per_in = X_MAX / ax_w
        ax.set_xlim(0.0, X_MAX)
        ax.set_ylim(0.0, ax_h * y_per_in)
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(SIM_INK)
            ax.spines[side].set_linewidth(2.2)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("time →", loc="right", labelpad=6)
        ax.set_ylabel("streamflow", labelpad=10)

        # volume: the area under the observed curve
        ax.fill_between(t, 0.0, obs, color=CALLOUT_COLOR["volume"], alpha=0.38, lw=0, zorder=1)
        # shape: a halo on the part of the observed curve above half its peak
        up = obs >= HALF * q_o
        ax.plot(t[up], obs[up], color=CALLOUT_COLOR["shape"], lw=LW_HALO, solid_capstyle="butt", zorder=3.2)
        # the two curves
        ax.plot(t, sim, color=SIM_INK, lw=LW_CURVE, ls=(0, (4.5, 2.6)), zorder=3.6)
        ax.plot(t, obs, color=OBS_INK, lw=LW_CURVE, zorder=4)

        # timing: dotted drops from the two peaks to a double arrow
        y_arrow = q_o + 0.36 * y_per_in
        dot = dict(color=CALLOUT_COLOR["timing"], lw=LW_DOT, ls=(0, (1, 2.2)), dash_capstyle="round", zorder=5)
        ax.plot([x_o, x_o], [q_o + 0.07 * y_per_in, y_arrow], **dot)
        ax.plot([x_s, x_s], [q_s + 0.07 * y_per_in, y_arrow], **dot)
        ax.annotate("", xy=(x_s, y_arrow), xytext=(x_o, y_arrow), zorder=5,
                    arrowprops=dict(arrowstyle="<|-|>", color=CALLOUT_COLOR["timing"], lw=LW_ARROW, shrinkA=0, shrinkB=0, mutation_scale=26))

        # tags: dark words on a solid patch of the callout colour
        def box(color):
            return dict(boxstyle="round,pad=0.28,rounding_size=0.3", fc=color, ec="none")
        tag = dict(fontsize=POSTER_PT, fontweight="bold", color=TEXT_INK, zorder=6)
        # timing tag: just right of the arrow's end
        ax.text(x_s + 0.42 * x_per_in, y_arrow, "timing", ha="left", va="center", bbox=box(CALLOUT_COLOR["timing"]), **tag)
        # volume tag: inside the shaded area, near the bottom
        ax.text(x_o, 0.37 * y_per_in, "volume", ha="center", va="center", bbox=box(CALLOUT_COLOR["volume"]), **tag)
        # shape tag: up and to the left of the rising limb, joined to the halo by a short line
        x_h = float(np.interp(0.75 * q_o, obs[:i_o + 1], t[:i_o + 1]))
        y_h = 0.75 * q_o
        ax.annotate("shape", xy=(x_h, y_h), xytext=(SHAPE_TAG_LEFT_IN * x_per_in, y_h + 0.30 * y_per_in), ha="left", va="center",
                    bbox=box(CALLOUT_COLOR["shape"]), arrowprops=dict(arrowstyle="-", color=CALLOUT_COLOR["shape"], lw=LW_ARROW + 0.4,
                                                                      shrinkA=0, shrinkB=0), **tag)

        # legend (top right) and the note under it
        handles = [Line2D([0], [0], color=OBS_INK, lw=LW_CURVE, label="observed"),
                   Line2D([0], [0], color=SIM_INK, lw=LW_CURVE, ls=(0, (4.5, 2.6)), label="simulated")]
        ax.legend(handles=handles, loc="upper right", frameon=False, handlelength=2.2, labelcolor=TEXT_INK, borderaxespad=0.0)
        if show_note:
            ax.text(1.0, 0.60, NOTE_TEXT, transform=ax.transAxes, ha="right", va="top", fontsize=POSTER_PT, color=SIM_INK, zorder=6)

        # check: nothing below the poster text floor
        sizes = [(tx.get_text(), tx.get_fontsize()) for tx in fig.findobj(Text) if tx.get_text().strip()]
        small = [(s, z) for s, z in sizes if z < POSTER_PT - 1e-9]
        if small:
            sys.exit("STOPPED: text smaller than %d pt in the figure: %s" % (POSTER_PT, small))

        fig.savefig(png_path, dpi=dpi, facecolor="none" if transparent else SURFACE, transparent=transparent)
        if svg_path is not None:
            fig.savefig(svg_path, facecolor="none" if transparent else SURFACE, transparent=transparent)
        plt.close(fig)
    return {"min_pt": min(z for _, z in sizes), "n_text": len(sizes), "area_ratio": area(t, sim) / area(t, obs),
            "peak_ratio": q_s / q_o, "peak_shift": x_s - x_o, "axes_in": (ax_w, ax_h), "half_width": x_dn - x_up}


# ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Poster: observed vs simulated hydrograph with timing / volume / shape callouts (schematic).")
    ap.add_argument("--width_in", type=float, default=10.5, help="figure width in inches (default 10.5)")
    ap.add_argument("--height_in", type=float, default=3.5, help="figure height in inches (default 3.5)")
    ap.add_argument("--dpi", type=int, default=300, help="pixels per inch (default 300)")
    ap.add_argument("--svg", action="store_true", help="also save a vector .svg (text becomes outlines)")
    ap.add_argument("--transparent", action="store_true", help="no background colour")
    ap.add_argument("--no_note", action="store_true", help="leave out the 'Schematic, not real data' note")
    ap.add_argument("--font_file", default=None, help="path to a .ttf font to use for all text instead of Liberation Sans")
    ap.add_argument("--out_dir", type=Path, default=None, help="folder for the files (default: see the top of this file)")
    args = ap.parse_args()
    if not (8.0 <= args.width_in <= 14.0 and 3.2 <= args.height_in <= 4.5):
        ap.error("--width_in must be 8 to 14 and --height_in 3.2 to 4.5 (the callouts are a fixed physical size, so a smaller figure crowds them).")
    if not (72 <= args.dpi <= 600):
        ap.error("--dpi must be 72 to 600.")

    cwd = Path.cwd()
    summary_dir = cwd.parent / "calibration_work" / "03_comparisons" / "summary_tables"
    out_dir = args.out_dir or ((summary_dir / OUT_SUBFOLDER) if summary_dir.exists() else (cwd / OUT_SUBFOLDER))

    print("\n" + "=" * 100)
    print("Poster scoring figure (version %s): %.1f x %.1f in at %d dpi = %d x %d pixels"
          % (VERSION, args.width_in, args.height_in, args.dpi, round(args.width_in * args.dpi), round(args.height_in * args.dpi)))
    print("=" * 100)
    font_family, font_note = resolve_font(args.font_file)
    print("  " + font_note)

    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "fig_poster_scoring_110.png"
    svg = out_dir / "fig_poster_scoring_110.svg" if args.svg else None
    info = draw_figure(png, svg, args.width_in, args.height_in, args.dpi, font_family, args.transparent, not args.no_note)

    print("\nCHECKS")
    print("  smallest text in the figure: %g pt (floor %d pt), %d text items" % (info["min_pt"], POSTER_PT, info["n_text"]))
    for name, color in CALLOUT_COLOR.items():
        print("  %-7s tag %s: contrast of the dark words on it %.1f : 1" % (name, color, contrast(TEXT_INK, color)))
    print("  schematic: simulated peak %.2f x the observed peak, %.2f time units later; simulated area %.2f x the observed area"
          % (info["peak_ratio"], info["peak_shift"], info["area_ratio"]))
    print("  plot area %.2f x %.2f in" % info["axes_in"])
    print("\nSaved to: %s" % out_dir)
    print("  " + png.name + ("   " + svg.name if svg else ""))


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        main()
