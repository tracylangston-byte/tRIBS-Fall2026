#!/usr/bin/env python
"""Summarize LANDFIRE vegetation cover and height for each Sycamore land use class.

READ-ONLY: this script prints tables and writes nothing.

What it does
  1. Reads LandUse.asc (the 8-class land use map, 10 m, EPSG 26912).
  2. Reads LANDFIRE 2016 Remap EVC (cover) and EVH (height), 30 m, EPSG 5070, and
     resamples each onto the LandUse grid (nearest neighbour, the same method the
     land use map script used for EVT).
  3. Keeps only the cells inside the basin outline.
  4. Decodes the LANDFIRE codes and prints, for each class, the typical plant cover
     (a data-based starting value for V) and plant height (for h).

LANDFIRE codes (from the .tif.vat.dbf tables of the 2016 Remap release)
  EVC  tree 110-199   cover % = code - 100      EVH  tree  101-199   height m = code - 100
       shrub 210-299  cover % = code - 200           shrub 201-299   height m = (code - 200) / 10
       herb 310-399   cover % = code - 300           herb  301-399   height m = (code - 300) / 10
       100 = sparse vegetation (under 10% cover; counted here as --sparse-cover percent)
       anything else (water, developed, barren, crops) = not vegetation, 0% cover
  Nodata is 32767.

Run it from anywhere (it finds its own folder):
    python3 workspaces/Sycamore/lab/summarize_landuse_cover_height.py
Options:  --sparse-cover 5   --boundary PATH   --landuse PATH   --evc PATH   --evh PATH
"""
from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.features import rasterize
from rasterio.warp import Resampling, reproject

LAB = Path(__file__).resolve().parent          # .../workspaces/Sycamore/lab
INIT = LAB.parent / "init_data"
LANDUSE = INIT / "LandUse.asc"
EVC_TIF = INIT / "landcover_raw" / "LF2016_EVC_CONUS.tif"
EVH_TIF = INIT / "landcover_raw" / "LF2016_EVH_CONUS.tif"
BOUNDARY_DIR = LAB / "data" / "preprocessing"

DST_CRS = CRS.from_epsg(26912)
SRC_NODATA = 32767
LU_NODATA = -9999

CLASS_NAMES = {1: "Woodland", 2: "Desert scrub", 3: "Chaparral", 4: "Bare rock",
               5: "Burned/logged", 6: "Developed", 7: "Riparian/water", 8: "Grassland"}

# Lifeform codes used inside this script
NODATA, TREE, SHRUB, HERB, SPARSE, NONVEG = 0, 1, 2, 3, 4, 5
LIFE_NAMES = {TREE: "tree", SHRUB: "shrub", HERB: "herb"}
# EVC codes that are legitimately "not vegetation" (water, developed, barren, crops)
KNOWN_NONVEG = {11, 13, 14, 15, 16, 17, 22, 23, 24, 25, 31, 32, 63, 64, 65, 68}


def find_boundary(arg: str | None) -> Path:
    if arg:
        return Path(arg)
    hits = sorted(BOUNDARY_DIR.glob("*boundary*.shp"))
    if len(hits) != 1:
        raise SystemExit(f"Expected exactly one *boundary*.shp in {BOUNDARY_DIR}, found {len(hits)}: "
                         f"{[h.name for h in hits]}. Pass the one you want with --boundary.")
    return hits[0]


def read_landuse(path: Path):
    if not path.exists():
        raise SystemExit(f"missing input: {path}")
    with rasterio.open(path) as src:
        arr = src.read(1)
        nd = src.nodata
        transform = src.transform
        crs = src.crs
    if crs is not None and crs.to_epsg() not in (None, 26912):
        raise SystemExit(f"{path.name} is in {crs}, expected EPSG 26912")
    lu = arr.copy()
    if nd is not None:
        lu = np.where(arr == nd, LU_NODATA, arr)
    return lu.astype(np.int32), transform, arr.shape


def resample_to_grid(tif: Path, shape, transform) -> np.ndarray:
    if not tif.exists():
        raise SystemExit(f"missing input: {tif}")
    dst = np.full(shape, SRC_NODATA, dtype=np.int16)
    with rasterio.open(tif) as src:
        src_nd = int(src.nodata) if src.nodata is not None else SRC_NODATA
        reproject(source=rasterio.band(src, 1), destination=dst,
                  dst_transform=transform, dst_crs=DST_CRS,
                  resampling=Resampling.nearest, src_nodata=src_nd, dst_nodata=SRC_NODATA)
    return dst


def basin_mask(boundary: Path, shape, transform) -> np.ndarray:
    if not boundary.exists():
        raise SystemExit(f"missing input: {boundary}")
    gdf = gpd.read_file(boundary)
    if gdf.crs is not None and gdf.crs.to_epsg() != 26912:
        gdf = gdf.to_crs(26912)
    geoms = [g for g in gdf.geometry if g is not None and not g.is_empty]
    burned = rasterize(((g, 1) for g in geoms), out_shape=shape, transform=transform,
                       fill=0, dtype="uint8")
    return burned.astype(bool)


def decode_cover(code: np.ndarray, sparse_cover_pct: float):
    """EVC code -> (lifeform, cover as a 0-1 fraction). Cover is NaN for nodata."""
    life = np.full(code.shape, NONVEG, dtype=np.int8)
    cover = np.zeros(code.shape, dtype=np.float32)
    tree = (code >= 110) & (code <= 199)
    shrub = (code >= 210) & (code <= 299)
    herb = (code >= 310) & (code <= 399)
    sparse = code == 100
    nodata = (code == SRC_NODATA) | (code < 0)
    life[tree], cover[tree] = TREE, code[tree] - 100
    life[shrub], cover[shrub] = SHRUB, code[shrub] - 200
    life[herb], cover[herb] = HERB, code[herb] - 300
    life[sparse], cover[sparse] = SPARSE, sparse_cover_pct
    life[nodata], cover[nodata] = NODATA, np.nan
    unknown = ~(tree | shrub | herb | sparse | nodata) & ~np.isin(code, list(KNOWN_NONVEG))
    return life, cover / 100.0, unknown


def decode_height(code: np.ndarray):
    """EVH code -> (lifeform, height in metres). Height is NaN where there is none."""
    life = np.full(code.shape, NONVEG, dtype=np.int8)
    h = np.full(code.shape, np.nan, dtype=np.float32)
    tree = (code >= 101) & (code <= 199)
    shrub = (code >= 201) & (code <= 299)
    herb = (code >= 301) & (code <= 399)
    life[tree], h[tree] = TREE, code[tree] - 100
    life[shrub], h[shrub] = SHRUB, (code[shrub] - 200) / 10.0
    life[herb], h[herb] = HERB, (code[herb] - 300) / 10.0
    return life, h


def fmt(x, width, nd=2):
    return f"{'-':>{width}}" if x is None or not np.isfinite(x) else f"{x:{width}.{nd}f}"


def quartiles(x):
    if x.size == 0:
        return (np.nan,) * 4
    q = np.percentile(x, [25, 50, 75])
    return float(np.mean(x)), float(q[0]), float(q[1]), float(q[2])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--landuse", default=str(LANDUSE))
    ap.add_argument("--evc", default=str(EVC_TIF))
    ap.add_argument("--evh", default=str(EVH_TIF))
    ap.add_argument("--boundary", default=None, help="basin outline .shp (default: the one *boundary*.shp in lab/data/preprocessing)")
    ap.add_argument("--sparse-cover", type=float, default=5.0,
                    help="percent cover to assign to LANDFIRE 'sparse vegetation' cells (default 5)")
    args = ap.parse_args(argv)

    lu, transform, shape = read_landuse(Path(args.landuse))
    boundary = find_boundary(args.boundary)
    print(f"land use map : {Path(args.landuse).name}  ({shape[1]} x {shape[0]} cells, "
          f"{abs(transform.a):g} x {abs(transform.e):g} m)")
    print(f"basin outline: {boundary.name}")
    print("resampling EVC and EVH onto the land use grid (nearest neighbour)...")
    evc = resample_to_grid(Path(args.evc), shape, transform)
    evh = resample_to_grid(Path(args.evh), shape, transform)
    mask = basin_mask(boundary, shape, transform)

    cell_km2 = abs(transform.a * transform.e) / 1e6
    n_basin = int(mask.sum())
    print(f"basin cells  : {n_basin:,}  =  {n_basin * cell_km2:.1f} km2 (expected about 424.6)\n")

    lu_b, evc_b, evh_b = lu[mask], evc[mask], evh[mask]
    life_c, cover, unknown = decode_cover(evc_b, args.sparse_cover)
    life_h, height = decode_height(evh_b)

    # ---- checks ----
    print("Checks")
    print(f"  basin cells with no land use class : {(lu_b <= 0).sum():,}   (expect 0)")
    print(f"  basin cells with LANDFIRE no data  : {(life_c == NODATA).sum():,}   (expect 0)")
    if unknown.any():
        codes, counts = np.unique(evc_b[unknown], return_counts=True)
        print(f"  EVC codes the script does not know : {dict(zip(codes.tolist(), counts.tolist()))}   (counted as not vegetation)")
    else:
        print("  EVC codes the script does not know : none")
    both_veg = np.isin(life_c, [TREE, SHRUB, HERB])
    mismatch = both_veg & (life_h != life_c)
    print(f"  cells where EVC and EVH disagree on tree/shrub/herb : {mismatch.sum():,}   (expect 0)")
    print()

    classes = sorted(int(c) for c in np.unique(lu_b) if c > 0)

    # ---- Table 1 ----
    print("TABLE 1. Lifeform mix, and average cover and height, for each class (inside the basin)")
    print("  V_veg = average cover of the vegetated cells; V_all = same, but non-vegetated cells count as 0 cover.")
    print("  h = LANDFIRE height of the cell's dominant lifeform (trees, shrubs or herbs), in metres.")
    print("  %w/h = share of the class's cells that have a height (sparse, barren and developed cells have none).")
    print("  If a class is mostly non-vegetated, V_veg and h describe only its few vegetated cells; read V_all instead.")
    head = (f"{'ID':>2} {'class':<15}{'km2':>7}{'%bas':>6} |{'tree%':>6}{'shrub%':>7}{'herb%':>6}{'spars%':>7}"
            f"{'nonveg%':>8}{'nodata%':>8} |{'V_veg':>6}{'V_all':>6} |{'h_mean':>7}{'h_med':>6}{'%w/h':>6}")
    print(head)
    print("-" * len(head))
    for c in classes:
        sel = lu_b == c
        n = int(sel.sum())
        counts = np.bincount(life_c[sel], minlength=6) / n * 100
        veg = sel & np.isin(life_c, [TREE, SHRUB, HERB, SPARSE])
        valid = sel & (life_c != NODATA)
        v_veg = float(np.mean(cover[veg])) if veg.any() else np.nan
        v_all = float(np.sum(cover[veg]) / valid.sum()) if valid.any() else np.nan
        hs = height[sel & np.isfinite(height)]
        h_mean = float(np.mean(hs)) if hs.size else np.nan
        h_med = float(np.median(hs)) if hs.size else np.nan
        h_share = hs.size / n * 100
        print(f"{c:>2} {CLASS_NAMES.get(c, 'class ' + str(c)):<15}{n * cell_km2:7.1f}{n / n_basin * 100:6.1f} |"
              f"{counts[TREE]:6.1f}{counts[SHRUB]:7.1f}{counts[HERB]:6.1f}{counts[SPARSE]:7.1f}"
              f"{counts[NONVEG]:8.1f}{counts[NODATA]:8.1f} |{fmt(v_veg, 6)}{fmt(v_all, 6)} |{fmt(h_mean, 7)}{fmt(h_med, 6)}{h_share:6.1f}")
    print()

    # ---- Table 2 ----
    print("TABLE 2. Cover and height by dominant lifeform (only lifeforms with at least 1% of the class, and 500+ cells)")
    print("  cover in % (a cell's own cover value), height in metres; p25 and p75 are the quartiles.")
    head2 = (f"{'ID':>2} {'class':<15}{'lifeform':<9}{'cells':>10}{'%class':>7} |{'cov_mean':>9}{'p25':>5}{'med':>5}{'p75':>5}"
             f" |{'h_mean':>7}{'p25':>6}{'med':>6}{'p75':>6}")
    print(head2)
    print("-" * len(head2))
    for c in classes:
        sel = lu_b == c
        n = int(sel.sum())
        for lf in (TREE, SHRUB, HERB):
            m = sel & (life_c == lf)
            k = int(m.sum())
            if k < 500 or k / n < 0.01:
                continue
            cm, c25, c50, c75 = quartiles(cover[m] * 100)
            hv = height[m & np.isfinite(height)]
            hm, h25, h50, h75 = quartiles(hv)
            print(f"{c:>2} {CLASS_NAMES.get(c, 'class ' + str(c)):<15}{LIFE_NAMES[lf]:<9}{k:10,d}{k / n * 100:7.1f} |"
                  f"{fmt(cm, 9, 1)}{fmt(c25, 5, 0)}{fmt(c50, 5, 0)}{fmt(c75, 5, 0)} |"
                  f"{fmt(hm, 7)}{fmt(h25, 6, 1)}{fmt(h50, 6, 1)}{fmt(h75, 6, 1)}")
    print("\nNothing was written. Paste both tables and the Checks block back into the chat.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
