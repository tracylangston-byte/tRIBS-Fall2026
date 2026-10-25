#!/usr/bin/env python
"""Slide: data inventory for the Lower Sycamore Wash model (what each input is, where it came from, whether it
transfers to a basin with no gauge, its status and known gaps).

Content is typed in from init_data/data_sources.md (checked 2026-10-14). Nothing is read from the model files, so this
script has no inputs. Writes out/fig6_data_inventory.png (16:9 slide) and out/data_inventory.csv.
"""
import argparse, textwrap
import pandas as pd
import matplotlib.pyplot as plt
import srp_common as c

NAT = ("National", "#cde2fb", "#184f95")          # works for any U.S. basin
LOC = ("Local", "#fbe0d3", "#8f3411")             # needs a gauge or a county network
NONE = ("", None, None)

# (input, source, resolution, used for, status key, transfers tag, transfers label, checks and gaps)
ROWS = [
    ("Terrain (DEM)", "USGS 3DEP 1/3 arc-second tiles, The National Map", "10 m",
     "Basin outline, stream network, mesh", "done", NAT, "National",
     "Reads 537.6 m at the gauge (USGS: about 537 m). The 54 MB file is not in the repo; re-download the two tiles."),
    ("Streamflow", "USGS gauge 09510200 (NWIS), 15-minute and daily", "15 min",
     "Calibration target: storm timing, peaks, volume", "done", LOC, "Needs a gauge",
     "Outage Jan 1-2, 2017 lost the first peak. Mar 12 interpolated; five dry days in December zero-filled (USGS daily means are 0)."),
    ("Basin boundary", "StreamStats basin; USGS Watershed Boundary Dataset (HUC-10, HUC-12)", "Vector, 1:24,000",
     "Cross-check on the delineated area", "done", NAT, "National",
     "Model uses the DEM outline, 424.6 km² (0.06% from StreamStats). The HUC-10 is 74 km² larger: it runs past the gauge."),
    ("Rainfall (radar)", "MRMS gauge-corrected hourly radar, Iowa State Mesonet archive", "1 km, hourly",
     "Rain forcing (WY2017: 567 mm basin mean)", "done", NAT, "National",
     "Reads 69-87% of FCDMC gauges (four storms). 17 hours missing (0.19%). Hour label assumed end-of-hour. Older than SRP's post-2020 product."),
    ("Rain gauges", "FCDMC ALERT: Mt. Ord, Sycamore Cr Upper, Iron Dike", "Point",
     "Plausibility check on the radar only", "partly", LOC, "County network",
     "Five hand-pulled storm reports so far. Mt. Ord read near zero on Jan 20 (outage or snow)."),
    ("Weather", "NLDAS-2 via NASA GES DISC and pytRIBS", "0.125° (about 12 km), hourly",
     "Temperature, humidity, wind, solar radiation, pressure", "partly", NAT, "National",
     "Files written, not yet read by tRIBS. Three stations (657-1,435 m) for ground at 537-2,330 m; no extra lapse-rate correction."),
    ("Soils and bedrock", "SOLUS100 via pytRIBS; Ks and f from the Rosetta v3 pedotransfer function", "100 m, 0-30 cm",
     "Ks, f, θs, θr, m, ψb; bedrock depth sets the starting water table", "partly", NAT, "National",
     "Built, not yet run. Ks and f are estimates, not measurements; about 9% of cells sit at the f floor."),
    ("Land cover", "LANDFIRE 2016 Remap: vegetation type, cover and height", "30 m, resampled to 10 m",
     "8 land classes; vegetation fraction and height", "partly", NAT, "National",
     "Map built; table script written Oct 14. Other parameters from South Mountain. 2.6% mapped burned or logged, unchecked."),
    ("Fire history", "USFS Bush Fire reports, Wikipedia, local news", "Acres as reported",
     "Choose a pre-fire year (2020 Bush Fire)", "partly", NONE, "–",
     "Not yet clipped to the basin (at most about 30%). 2012 Sunflower Fire and WY2025 fires unchecked."),
]
STATUS = {"done": "Done", "partly": "In progress", "open": "Open"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(c.OUT / "fig6_data_inventory.png"))
    args = ap.parse_args()
    c.style()

    pd.DataFrame([dict(input=r[0], source=r[1], resolution=r[2], used_for=r[3], status=STATUS[r[4]], transfers=r[6], checks_and_gaps=r[7])
                  for r in ROWS]).to_csv(c.OUT / "data_inventory.csv", index=False)

    n_nat = sum(1 for r in ROWS if r[5] is NAT)
    n_loc = sum(1 for r in ROWS if r[5] is LOC)
    words = {6: "Six", 7: "Seven", 5: "Five", 4: "Four"}
    nums = {8: "eight", 7: "seven", 9: "nine", 6: "six"}
    head = (f"{words.get(n_nat, n_nat)} of the {nums.get(n_nat + n_loc, n_nat + n_loc)} inputs so far come from national datasets; only the flow record "
            f"and the rain-gauge check need local data")
    W, H = 13.333, 7.5
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.05, 0.935, "\n".join(textwrap.wrap(head, 80)), fontsize=21, fontweight="bold", va="center", linespacing=1.25)
    ax.text(0.05, 0.848, "What goes into the Lower Sycamore Wash model, where it came from, and which parts carry over to a basin with no gauge",
            fontsize=12.5, color=c.INK2, va="center")

    X = dict(inp=0.05, src=0.148, res=0.32, use=0.395, st=0.53, tr=0.625, note=0.72)
    WRAP = dict(src=28, res=13, use=23, note=49)
    y_head = 0.792
    for k, label in [("inp", "Input"), ("src", "Source"), ("res", "Resolution"), ("use", "Used for"), ("st", "Status"),
                     ("tr", "Transfers?"), ("note", "Checks and gaps")]:
        ax.text(X[k], y_head, label, fontsize=11.5, fontweight="bold", color=c.INK2, ha="left", va="center")
    ax.plot([0.05, 0.97], [y_head - 0.027, y_head - 0.027], color=c.INK2, lw=1.2)

    dy = 0.0775
    y0 = y_head - 0.027 - dy / 2
    for i, (inp, src, res, use, st, tag, tlabel, note) in enumerate(ROWS):
        y = y0 - i * dy
        ax.text(X["inp"], y, "\n".join(textwrap.wrap(inp, 13)), fontsize=11.5, fontweight="bold", ha="left", va="center", linespacing=1.2)
        ax.text(X["src"], y, "\n".join(textwrap.wrap(src, WRAP["src"])), fontsize=9.5, ha="left", va="center", linespacing=1.22)
        ax.text(X["res"], y, "\n".join(textwrap.wrap(res, WRAP["res"])), fontsize=9.5, ha="left", va="center", linespacing=1.22)
        ax.text(X["use"], y, "\n".join(textwrap.wrap(use, WRAP["use"])), fontsize=9.5, ha="left", va="center", linespacing=1.22)
        # status marker: filled = done, half = in progress, hollow = open
        kw = dict(marker="o", ms=10, mec=c.BLUE, mew=1.6, ls="none", clip_on=False)
        if st == "done":
            ax.plot([X["st"] + 0.007], [y], mfc=c.BLUE, **kw)
        elif st == "partly":
            ax.plot([X["st"] + 0.007], [y], fillstyle="left", mfc=c.BLUE, mfcalt=c.SURFACE, **kw)
        else:
            ax.plot([X["st"] + 0.007], [y], mfc=c.SURFACE, **kw)
        ax.text(X["st"] + 0.02, y, STATUS[st], fontsize=10, ha="left", va="center")
        if tag[1]:
            ax.text(X["tr"], y, tlabel, fontsize=9.5, fontweight="bold", color=tag[2], ha="left", va="center",
                    bbox=dict(boxstyle="round,pad=0.32,rounding_size=0.5", fc=tag[1], ec="none"))
        else:
            ax.text(X["tr"], y, tlabel, fontsize=10, color=c.MUTED, ha="left", va="center")
        ax.text(X["note"], y, "\n".join(textwrap.wrap(note, WRAP["note"])), fontsize=9, color=c.INK2, ha="left", va="center", linespacing=1.22)
        ax.plot([0.05, 0.97], [y - dy / 2, y - dy / 2], color=c.GRID, lw=1.0)

    foot = ("From the project's data log (data_sources.md), checked Oct 14, 2026. Filled dot = done, half = in progress, hollow = open. "
            "National = the same dataset exists for any U.S. basin; the others need a gauge or a county network. Numbers are for water year 2017.")
    ax.text(0.05, 0.034, "\n".join(textwrap.wrap(foot, 170)), fontsize=9, color=c.INK2, va="center", linespacing=1.4)
    fig.savefig(args.out, dpi=200)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
