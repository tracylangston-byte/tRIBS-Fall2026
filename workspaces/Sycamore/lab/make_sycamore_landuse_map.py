#!/usr/bin/env python
"""Collapse LANDFIRE 2016 Remap EVT into 8 tRIBS land use classes for Lower Sycamore Wash.

Adapted from Josh Cederstrom's make_landuse_map.py (South Mountain).  What changed:
  * 8 Sycamore classes instead of 4 (woodland is ~47% of this basin, so it needs its own class);
  * target grid is Sycamore_DEM.tif (10 m, NAD83 UTM 12N), because there is no soil grid yet;
  * the whole DEM extent is classified, including the 5 km buffer, so mesh boundary nodes that
    sit just outside the basin line never land on NODATA;
  * the output is init_data/LandUse.asc, the file name the Sycamore notebook loads;
  * class areas are also reported inside the basin boundary, to compare with the ArcGIS Pro numbers.

What it does
  1. Reads the 30 m LANDFIRE Existing Vegetation Type (EVT) raster and its class table (.vat.dbf).
  2. Resamples the EVT codes onto the DEM grid with NEAREST NEIGHBOUR (codes are categories, so
     they must never be averaged), so land cover lines up cell for cell with the DEM.
  3. Sorts every EVT type into one of the 8 classes below.
  4. Writes LandUse.asc (ASCII grid, NODATA -9999) and a crosswalk CSV that says which EVT type
     went into which class.

How an EVT type gets its class, first match wins:
  BY_VALUE   the explicit list of every EVT code found in the Sycamore basin (hand-checked)
  RULES_GP   the EVT group name (EVT_GP_N), for codes that only occur in the 5 km buffer
  RULES_PHYS the EVT physiognomy (EVT_PHYS)
  RULES_LF   the EVT lifeform (EVT_LF)
  Codes that fall through to nothing stop the script and are listed, so none is dropped silently.
  Codes classed by a rule rather than by BY_VALUE are flagged in the printout and the CSV.

Run:
    python make_sycamore_landuse_map.py --dry-run     # print class areas, write nothing
    python make_sycamore_landuse_map.py               # write LandUse.asc and the crosswalk CSV

Default locations (all inside workspaces/Sycamore/, change with the options below):
    init_data/landcover_raw/LF2016_EVT_CONUS.tif       (+ LF2016_EVT_CONUS.tif.vat.dbf beside it)
    init_data/Sycamore_DEM.tif                         target grid
    lab/data/preprocessing/*boundary*.shp              basin outline, used only for the area report;
                                                       found automatically if exactly one file matches
                                                       (boundary / watershed / basin), else pass --boundary
    -> init_data/LandUse.asc, init_data/LandUse_crosswalk.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.warp import Resampling, reproject
from scipy.ndimage import distance_transform_edt

def _find_workspace() -> Path:
    """workspaces/Sycamore: the nearest folder above this script that contains init_data/."""
    here = Path(__file__).resolve().parent
    for p in (here, *here.parents):
        if (p / "init_data").is_dir():
            return p
    return here.parent


WS = _find_workspace()                           # .../workspaces/Sycamore
INIT = WS / "init_data"
PREP = WS / "lab" / "data" / "preprocessing"     # where the delineation outputs live
EVT_TIF = INIT / "landcover_raw" / "LF2016_EVT_CONUS.tif"
DEM_TIF = INIT / "Sycamore_DEM.tif"
OUT_ASC = INIT / "LandUse.asc"
OUT_CSV = INIT / "LandUse_crosswalk.csv"

NODATA = -9999
TARGET_EPSG = 26912          # NAD83 / UTM zone 12N
MAX_FILL_FRACTION = 0.01     # stop if more than 1% of the grid falls outside the EVT download

# tRIBS class ID -> name
CLASSES = {
    1: "Woodland",                       # pinyon-juniper, Madrean oak, pine
    2: "Desert scrub",                   # Sonoran / creosote scrub, ruderal and planted scrub
    3: "Chaparral and oak shrub",
    4: "Bare rock",                      # bedrock outcrop, desert pavement
    5: "Recently burned or logged",      # as mapped in the 2016 base map
    6: "Developed",                      # roads, buildings, urban vegetation
    7: "Riparian and water",
    8: "Grassland",                      # semi-desert grassland, ruderal and planted grass
}

# 1. EVT VALUE -> class, for every code found inside the Sycamore basin (45 codes).
_BY_CLASS = {
    1: [7023,   # Madrean Encinal
        7024,   # Madrean Lower Montane Pine-Oak Forest and Woodland
        7025,   # Madrean Pinyon-Juniper Woodland
        7026,   # Madrean Upper Montane Conifer-Oak Forest and Woodland
        7051,   # S. Rocky Mtn Dry-Mesic Montane Mixed Conifer Forest and Woodland
        7054,   # S. Rocky Mtn Ponderosa Pine Woodland
        7116],  # Madrean Juniper Savanna
    2: [7087,   # Sonora-Mojave Creosotebush-White Bursage Desert Scrub
        7091,   # Sonoran Mid-Elevation Desert Scrub
        7109,   # Sonoran Paloverde-Mixed Cacti Desert Scrub
        9310],  # N. American Warm Desert Ruderal & Planted Scrub
    3: [7104,   # Mogollon Chaparral
        7107,   # Rocky Mountain Gambel Oak-Mixed Montane Shrubland
        7108],  # Sonora-Mojave Semi-Desert Chaparral
    4: [9147,   # N. American Warm Desert Bedrock Cliff and Outcrop
        9150],  # N. American Warm Desert Pavement
    5: [7192,   # Recently Logged-Shrub Cover
        7193,   # Recently Logged-Tree Cover
        7195,   # Recently Burned-Herb and Grass Cover
        7196,   # Recently Burned-Shrub Cover
        7197],  # Recently Burned-Tree Cover
    6: [7296, 7297, 7298, 7299,           # Developed low / medium / high intensity, roads
        7910, 7911, 7912, 7913, 7914],    # Urban deciduous / evergreen / mixed forest, herb, shrub
    7: [7292,   # Open Water
        9011,   # Arid West Emergent Marsh
        9019,   # Rocky Mountain Lower Montane-Foothill Riparian Woodland
        9034,   # N. American Warm Desert Riparian Woodland
        9035,   # N. American Warm Desert Lower Montane Riparian Woodland
        9148,   # N. American Warm Desert Cienega
        9152,   # N. American Warm Desert Riparian Mesquite Bosque Woodland
        9329,   # Western N. American Ruderal Wet Shrubland
        9533,   # N. American Warm Desert Riparian Herbaceous
        9534,   # N. American Warm Desert Riparian Shrubland
        9652,   # N. American Warm Desert Riparian Mesquite Bosque Shrubland
        9829],  # Western N. American Ruderal Wet Meadow & Marsh
    8: [7121,   # Apacherian-Chihuahuan Semi-Desert Shrub-Steppe
        7256,   # Apacherian-Chihuahuan Semi-Desert Grassland
        9810],  # N. American Warm Desert Ruderal & Planted Grassland
}
BY_VALUE = {v: c for c, vals in _BY_CLASS.items() for v in vals}

# 2. EVT_GP_N (group name) -> class, for codes that are not in BY_VALUE.
RULES_GP = {
    "Chaparral": 3,
    "Deciduous Shrubland": 3,
    "Desert Scrub": 2,
    "Creosotebush Desert Scrub": 2,
    "Introduced Upland Vegetation-Shrub": 2,
    "Grassland and Steppe": 8,
    "Introduced Perennial Grassland and Forbland": 8,
    "Transitional Herbaceous Vegetation": 5,     # recently burned / logged
    "Transitional Shrub Vegetation": 5,
    "Transitional Forest Vegetation": 5,
}

# 3. EVT_PHYS (physiognomy) -> class.
RULES_PHYS = {
    "Conifer": 1, "Conifer-Hardwood": 1, "Hardwood": 1,
    "Shrubland": 2, "Exotic Tree-Shrub": 2,
    "Grassland": 8, "Exotic Herbaceous": 8, "Agricultural": 8,
    "Sparsely Vegetated": 4,
    "Quarries-Strip Mines-Gravel Pits-Well and Wind Pads": 4,
    "Riparian": 7, "Open Water": 7,
    "Developed": 6, "Developed-Roads": 6, "Developed-Low Intensity": 6,
    "Developed-Medium Intensity": 6, "Developed-High Intensity": 6,
}

# 4. EVT_LF (lifeform) -> class, the last resort.
RULES_LF = {"Tree": 1, "Shrub": 2, "Herb": 8, "Agriculture": 8,
            "Sparse": 4, "Barren": 4, "Water": 7, "Developed": 6}


def crosswalk(vat: pd.DataFrame) -> dict[int, tuple[int, str]]:
    """EVT VALUE -> (class ID, how it was assigned)."""
    cw: dict[int, tuple[int, str]] = {}
    for r in vat.itertuples():
        v = int(r.VALUE)
        if v in BY_VALUE:
            cw[v] = (BY_VALUE[v], "value")
        elif r.EVT_GP_N in RULES_GP:
            cw[v] = (RULES_GP[r.EVT_GP_N], "rule: group")
        elif r.EVT_PHYS in RULES_PHYS:
            cw[v] = (RULES_PHYS[r.EVT_PHYS], "rule: physiognomy")
        elif r.EVT_LF in RULES_LF:
            cw[v] = (RULES_LF[r.EVT_LF], "rule: lifeform")
    return cw


def _find_boundary() -> Path | None:
    hits = [p for p in sorted(PREP.glob("*.shp"))
            if any(k in p.stem.lower() for k in ("boundary", "watershed", "basin"))]
    if len(hits) == 1:
        return hits[0]
    if hits:
        print("note: several possible basin outlines in lab/data/preprocessing/, none used. "
              "Pick one with --boundary:\n  " + "\n  ".join(p.name for p in hits) + "\n")
    return None


def build(evt_tif: Path, dem_tif: Path, boundary: Path | None, out_asc: Path, dry_run: bool) -> None:
    evt_vat = evt_tif.with_suffix(".tif.vat.dbf")
    for p in (evt_tif, evt_vat, dem_tif):
        if not p.exists():
            raise SystemExit(f"missing input: {p}")

    vat = gpd.read_file(evt_vat)
    cw = crosswalk(vat)
    names = vat.assign(VALUE=vat.VALUE.astype(int)).set_index("VALUE").EVT_NAME

    stale = sorted(set(BY_VALUE) - set(names.index))
    if stale:
        print(f"note: BY_VALUE codes that are not in this tile's class table: {stale}\n")

    # Target grid = the DEM grid.
    with rasterio.open(dem_tif) as dem:
        dst_tf, dst_crs, ny, nx = dem.transform, dem.crs, dem.height, dem.width
    if dst_crs is None or dst_crs.to_epsg() != TARGET_EPSG:
        raise SystemExit(f"{dem_tif.name} is not in EPSG:{TARGET_EPSG} (found {dst_crs})")
    if abs(abs(dst_tf.a) - abs(dst_tf.e)) > 1e-6 or dst_tf.b != 0 or dst_tf.d != 0 or dst_tf.e > 0:
        raise SystemExit("the DEM grid must have square cells and be north-up")
    cs = abs(dst_tf.a)

    # Resample the EVT *codes* (not classes) onto the DEM grid, then look the classes up.
    with rasterio.open(evt_tif) as src:
        evt = src.read(1).astype(np.int32)
        src_nodata = src.nodata if src.nodata is not None else NODATA
        src_crs, src_tf = src.crs, src.transform
    evt10 = np.full((ny, nx), NODATA, dtype=np.int32)
    reproject(evt, evt10, src_transform=src_tf, src_crs=src_crs, src_nodata=src_nodata,
              dst_transform=dst_tf, dst_crs=dst_crs, dst_nodata=NODATA,
              resampling=Resampling.nearest)

    ok = evt10 >= 0
    present = sorted(int(v) for v in np.unique(evt10[ok]))
    unmapped = [v for v in present if v not in cw]
    if unmapped:
        lines = [f"  {v}  {names.get(v, '(not in class table)')}" for v in unmapped]
        raise SystemExit("EVT codes with no class (add them to BY_VALUE):\n" + "\n".join(lines))

    lut = np.full(max(present) + 1, NODATA, dtype=np.int16)
    for v in present:
        lut[v] = cw[v][0]
    cls = np.full((ny, nx), NODATA, dtype=np.int16)
    cls[ok] = lut[evt10[ok]]

    # Cells that fall outside the EVT download (should be none): fill from the nearest class.
    gaps = ~ok
    n_gap = int(gaps.sum())
    if n_gap:
        if n_gap / cls.size > MAX_FILL_FRACTION:
            raise SystemExit(f"{n_gap} DEM cells ({100 * n_gap / cls.size:.1f}%) fall outside the "
                             "EVT download - check the download box and the projection")
        _, (ir, ic) = distance_transform_edt(gaps, return_indices=True)
        cls[gaps] = cls[ir[gaps], ic[gaps]]
        print(f"filled {n_gap} DEM cells outside the EVT tile from the nearest class\n")

    # Basin outline, used only to report areas inside the basin.
    basin = None
    if boundary is None:
        boundary = _find_boundary()
    if boundary is not None and not boundary.exists():
        raise SystemExit(f"missing basin outline: {boundary}")
    if boundary is not None:
        print(f"basin outline for the area report: {boundary.name}\n")
        b = gpd.read_file(boundary)
        b = b.set_crs(epsg=TARGET_EPSG) if b.crs is None else b.to_crs(epsg=TARGET_EPSG)
        shapes = [(g, 1) for g in b.geometry if g is not None and not g.is_empty]
        basin = rasterize(shapes, out_shape=(ny, nx), transform=dst_tf, fill=0,
                          dtype="uint8").astype(bool)
    else:
        print("note: no basin outline found; reporting the whole DEM grid only\n")

    _report(names, cw, cls, evt10, basin, cs, evt_tif, dem_tif)

    if dry_run:
        print("\n-- dry run, nothing written --")
        return

    out_asc.parent.mkdir(parents=True, exist_ok=True)
    xll, yll = dst_tf.c, dst_tf.f + ny * dst_tf.e        # lower-left corner of the grid
    with open(out_asc, "w") as f:
        f.write(f"ncols        {nx}\nnrows        {ny}\n"
                f"xllcorner    {xll:.6f}\nyllcorner    {yll:.6f}\n"
                f"cellsize     {cs:.6f}\nNODATA_value {NODATA}\n")
        np.savetxt(f, cls, fmt="%d")

    out_csv = out_asc.with_name(OUT_CSV.name)
    vals, counts = np.unique(evt10[ok], return_counts=True)
    grid_cells = dict(zip(vals.tolist(), counts.tolist()))
    rows = pd.DataFrame({"VALUE": present})
    rows["EVT_NAME"] = rows.VALUE.map(names)
    for col in ("EVT_PHYS", "EVT_LF", "EVT_GP_N"):
        rows[col] = rows.VALUE.map(vat.assign(VALUE=vat.VALUE.astype(int)).set_index("VALUE")[col])
    rows["TRIBS_ID"] = rows.VALUE.map(lambda v: cw[v][0])
    rows["TRIBS_CLASS"] = rows.TRIBS_ID.map(CLASSES)
    rows["ASSIGNED_BY"] = rows.VALUE.map(lambda v: cw[v][1])
    rows["GRID_CELLS"] = rows.VALUE.map(grid_cells)
    if basin is not None:
        bv, bc = np.unique(evt10[ok & basin], return_counts=True)
        in_basin = dict(zip(bv.tolist(), bc.tolist()))
        rows["BASIN_CELLS"] = rows.VALUE.map(in_basin).fillna(0).astype(int)
    rows.sort_values(["TRIBS_ID", "VALUE"]).to_csv(out_csv, index=False)
    out_asc.with_suffix(".prj").write_text(dst_crs.to_wkt())     # so GIS programs know the CRS
    print(f"\nwrote {out_asc}\nwrote {out_asc.with_suffix('.prj').name}\nwrote {out_csv}")


def _report(names, cw, cls, evt10, basin, cs, evt_tif, dem_tif) -> None:
    cell_km2 = cs ** 2 / 1e6
    print(f"source: {evt_tif.name}  ->  grid of {dem_tif.name} "
          f"({cls.shape[1]} x {cls.shape[0]} cells of {cs:g} m)\n")

    in_grid = {int(v) for v in np.unique(evt10[evt10 >= 0])}
    for cid, cname in CLASSES.items():
        members = sorted(v for v, (c, _) in cw.items() if c == cid and v in in_grid)
        print(f"{cid}  {cname}  ({len(members)} EVT types in the grid)")
        for v in members:
            flag = "" if cw[v][1] == "value" else f"   <- by {cw[v][1]}"
            print(f"      {v}  {names[v]}{flag}")

    n_grid = int((cls != NODATA).sum())
    head = "\nID  class                      "
    if basin is not None:
        n_basin = int(((cls != NODATA) & basin).sum())
        print(head + "basin km2  basin %   whole grid km2")
    else:
        n_basin = 0
        print(head + "whole grid km2")
    for cid, cname in CLASSES.items():
        n_g = int((cls == cid).sum())
        line = f"{cid:>2}  {cname:<26}"
        if basin is not None:
            n_b = int(((cls == cid) & basin).sum())
            line += f"{n_b * cell_km2:>9.2f}  {100 * n_b / n_basin:>7.1f}   "
        print(line + f"{n_g * cell_km2:>14.2f}")
    if basin is not None:
        print(f"\nbasin: {n_basin:,} cells = {n_basin * cell_km2:.2f} km2   "
              f"(whole grid {n_grid * cell_km2:.1f} km2)")

    by_rule = sorted(v for v, (_, how) in cw.items() if how != "value" and v in in_grid)
    if by_rule:
        print("\nEVT types classed by a rule, not by an explicit BY_VALUE entry (please review):")
        for v in by_rule:
            print(f"  {v}  {names[v]}  -> {cw[v][0]} {CLASSES[cw[v][0]]}  ({cw[v][1]})")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="print the summary without writing")
    ap.add_argument("--evt-tif", type=Path, default=EVT_TIF, help="LANDFIRE EVT .tif (its .tif.vat.dbf must sit beside it)")
    ap.add_argument("--dem", type=Path, default=DEM_TIF, help="DEM that defines the target grid")
    ap.add_argument("--boundary", type=Path, default=None, help="basin outline shapefile (area report only)")
    ap.add_argument("--out", type=Path, default=OUT_ASC, help="output ASCII grid")
    args = ap.parse_args(argv)
    build(args.evt_tif, args.dem, args.boundary, args.out, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
