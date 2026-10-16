#!/usr/bin/env python
"""Write the tRIBS land use reclassification table (.ldt) for the Lower Sycamore Wash model.

Adapted from Josh's make_landuse_table.py (South Mountain). Pairs with LandUse.asc from
make_sycamore_landuse_map.py: one map and one table, 8 classes.

Table format is tRIBS v6.0.0: 14 columns, no `a` or `b1`, with RZD_m added (9999.99 = tRIBS default).
    ID, P, S, K, b2, Al, h, Kt, Rs, V, LAI, theta*_s, theta*_t, RZD_m

Parameters (Table 3.2, tRIBS docs):
    P         free throughfall coefficient, Rutter (-)
    S         canopy field capacity, Rutter (mm)
    K         drainage coefficient, Rutter (mm/hr)
    b2        drainage exponent, Rutter (1/mm)
    Al        albedo (-)
    h         vegetation height (m)
    Kt        optical transmission coefficient (-)
    Rs        stomatal resistance (s/m)
    V         vegetation fraction (-)
    LAI       leaf area index (-)
    theta*_s  stress threshold for evaporation   (absolute soil moisture, between thetaR and thetaS)
    theta*_t  stress threshold for transpiration (same range)

Where the numbers come from is written next to each class below. V and h are the LANDFIRE 2016
Remap EVC (cover) and EVH (height) averages inside the basin, from
summarize_landuse_cover_height.py (run 2026-10-14).

    python3 make_sycamore_landuse_table.py --dry-run     # print the table and run the checks, write nothing
    python3 make_sycamore_landuse_table.py               # same, then write init_data/land_use_params.ldt
"""
from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import rasterize

LAB = Path(__file__).resolve().parent                 # .../workspaces/Sycamore/lab
INIT = LAB.parent / "init_data"
LAND_RAS = INIT / "LandUse.asc"
THETA_S_GRID = INIT / "solus" / "theta_s_0-5cm.asc"   # the grids scgrid.gdf points tRIBS at (TS and TR)
THETA_R_GRID = INIT / "solus" / "theta_r_0-5cm.asc"
BOUNDARY_DIR = LAB / "data" / "preprocessing"
OUT_LDT = INIT / "land_use_params.ldt"

# Soil moisture stress thresholds, the same for every class (as in Josh's script).
# They are absolute soil moisture values and must sit between the largest thetaR and the
# smallest thetaS anywhere in the basin (the basin's thetaS runs 0.34-0.43, thetaR 0.06-0.10).
# theta*_t = 0.30 is the SMF value. theta*_s = 0.33 is as close to SMF's 0.37 as the grids allow
# (SMF's 0.37 is above Sycamore's smallest thetaS of 0.34 and would stop tRIBS).
THETA_S_STAR = 0.33
THETA_T_STAR = 0.30

# Rutter canopy values (P, S, K, b2), albedo, Kt, Rs and LAI are Josh's South Mountain starting
# values; the sensitivity studies (Huner 2025) found them to matter little for outlet flow.
# V and h are Sycamore-specific, from LANDFIRE.
LAND_PARAMS = {
    # 1 Woodland (Madrean pinyon-juniper, encinal, pine-oak). 198.4 km2.
    #   V = 0.28 and h = 6.3 m: LANDFIRE tree cover 28.4% and tree height 6.28 m, basin average.
    #   LAI 1.8: pinyon-juniper stands in north-central Arizona measured 1.72-1.85 (Grier et al. 1992).
    1: {"P": 0.4, "S": 1.5, "K": 0.12, "b2": 4.7, "Al": 0.18, "h": 6.3,
        "Kt": 0.4, "Rs": 120, "V": 0.28, "LAI": 1.8},
    # 2 Desert scrub (Sonoran paloverde-cacti, creosotebush). 150.2 km2.
    #   Josh's desert shrubland row, with V = 0.23 and h = 0.7 m from LANDFIRE (shrub cover 22.8%,
    #   shrub height 0.67 m). Josh used V 0.25 and h 1.0 m for South Mountain.
    2: {"P": 0.4, "S": 1.5, "K": 0.12, "b2": 4.7, "Al": 0.18, "h": 0.7,
        "Kt": 0.4, "Rs": 120, "V": 0.23, "LAI": 1.5},
    # 3 Chaparral and oak shrub. 43.7 km2.
    #   Shrub row with LANDFIRE V = 0.24 (cover 23.6%) and h = 0.5 m (0.46 m). NOTE: the LANDFIRE
    #   description of Mogollon chaparral says 40-80% cover and 1-2.5 m shrubs, so the mapped values may
    #   be low. First candidate to revisit in calibration.
    3: {"P": 0.4, "S": 1.5, "K": 0.12, "b2": 4.7, "Al": 0.18, "h": 0.5,
        "Kt": 0.4, "Rs": 120, "V": 0.24, "LAI": 1.5},
    # 4 Bare rock (bedrock outcrop, desert pavement). 12.1 km2.
    #   Josh's barren row. LANDFIRE codes 99% of these cells "sparse vegetation" (under 10% cover).
    4: {"P": 0.99, "S": 0.01, "K": 0.001, "b2": 0.001, "Al": 0.25, "h": 0.01,
        "Kt": 0.99, "Rs": 9999, "V": 0.01, "LAI": 0.01},
    # 5 Recently burned or logged (as mapped in the 2016 base map). 11.0 km2.
    #   LANDFIRE codes 82% of these cells as trees: V = 0.29 (cover 28.9%), h = 4.7 m (mean height), so
    #   the woodland row with its own V and h. Pre-2017 fire perimeter check still open.
    5: {"P": 0.4, "S": 1.5, "K": 0.12, "b2": 4.7, "Al": 0.18, "h": 4.7,
        "Kt": 0.4, "Rs": 120, "V": 0.29, "LAI": 1.8},
    # 6 Developed (roads, developed, urban vegetation). 3.4 km2. Josh's developed row; LANDFIRE shows
    #   no vegetation cover in these cells.
    6: {"P": 0.99, "S": 0.01, "K": 0.001, "b2": 0.001, "Al": 0.15, "h": 0.01,
        "Kt": 0.99, "Rs": 9999, "V": 0.01, "LAI": 0.01},
    # 7 Riparian and water. 3.0 km2.
    #   Shrub row with LANDFIRE V = 0.31 (cover 31%) and h = 4.3 m: the class-average height of 42% trees
    #   (about 8 m) and 54% shrubs (about 1 m). Water cells are only 0.2% of the class.
    7: {"P": 0.4, "S": 1.5, "K": 0.12, "b2": 4.7, "Al": 0.18, "h": 4.3,
        "Kt": 0.4, "Rs": 120, "V": 0.31, "LAI": 1.5},
    # 8 Grassland. 2.8 km2.
    #   Josh's grassland row with LANDFIRE V = 0.27 (cover 27%) and h = 0.3 m (0.32 m). Josh used V 0.20.
    8: {"P": 0.6, "S": 0.5, "K": 0.12, "b2": 4.7, "Al": 0.22, "h": 0.3,
        "Kt": 0.6, "Rs": 100, "V": 0.27, "LAI": 1.0},
}
CLASS_NAMES = {1: "Woodland", 2: "Desert scrub", 3: "Chaparral", 4: "Bare rock",
               5: "Burned/logged", 6: "Developed", 7: "Riparian/water", 8: "Grassland"}

COLS = ["P", "S", "K", "b2", "Al", "h", "Kt", "Rs", "V", "LAI", "theta*_s", "theta*_t"]
# Allowed ranges, as the tRIBS 6.0.0 land table validator checks them (everything else just >= 0)
UNIT_INTERVAL = {"P", "Al", "Kt", "V", "theta*_s", "theta*_t"}


def find_boundary(arg: str | None) -> Path:
    if arg:
        return Path(arg)
    hits = sorted(BOUNDARY_DIR.glob("*boundary*.shp"))
    if len(hits) != 1:
        raise SystemExit(f"Expected exactly one *boundary*.shp in {BOUNDARY_DIR}, found {len(hits)}: "
                         f"{[h.name for h in hits]}. Pass the one you want with --boundary.")
    return hits[0]


def check_map_ids(land_ras: Path) -> None:
    """The table must cover exactly the classes present in the map."""
    if not land_ras.exists():
        raise SystemExit(f"missing input: {land_ras}")
    with rasterio.open(land_ras) as src:
        arr = src.read(1)
        nd = src.nodata
    ids = {int(v) for v in np.unique(arr) if v > 0 and (nd is None or v != nd)}
    if ids != set(LAND_PARAMS):
        raise SystemExit(f"land use map IDs {sorted(ids)} do not match table IDs {sorted(LAND_PARAMS)}")
    print(f"land use map : {land_ras.name} has classes {sorted(ids)}; the table covers exactly these")


def theta_range(grid: Path, geoms) -> tuple[float, float, int]:
    """Smallest and largest value of a soil moisture grid under the given shapes (valid cells only)."""
    if not grid.exists():
        raise SystemExit(f"missing input: {grid}")
    with rasterio.open(grid) as src:
        a = src.read(1).astype("float64")
        tr = src.transform
        crs = src.crs
    if crs is not None and crs.to_epsg() not in (None, 26912):
        raise SystemExit(f"{grid.name} is in {crs}, expected EPSG 26912")
    inside = rasterize(((g, 1) for g in geoms), out_shape=a.shape, transform=tr, fill=0, dtype="uint8").astype(bool)
    ok = inside & np.isfinite(a) & (a > 0) & (a < 1)       # drops no-data and zero padding
    if not ok.any():
        raise SystemExit(f"no valid cells of {grid.name} fall inside the basin; check the boundary and the grid")
    return float(a[ok].min()), float(a[ok].max()), int(ok.sum())


def check_thresholds(table: list[dict], boundary: Path, buffer_m: float) -> None:
    """Stress thresholds must sit between the largest thetaR and the smallest thetaS near the basin."""
    gdf = gpd.read_file(boundary)
    if gdf.crs is not None and gdf.crs.to_epsg() != 26912:
        gdf = gdf.to_crs(26912)
    geoms = [g.buffer(buffer_m) if buffer_m > 0 else g for g in gdf.geometry if g is not None and not g.is_empty]
    ts_min, ts_max, n_s = theta_range(THETA_S_GRID, geoms)
    tr_min, tr_max, n_r = theta_range(THETA_R_GRID, geoms)
    where = f"the basin plus a {buffer_m:g} m margin" if buffer_m > 0 else "the basin"
    print(f"soil grids   : over {where}: thetaS {ts_min:.3f} to {ts_max:.3f} ({n_s:,} cells), "
          f"thetaR {tr_min:.3f} to {tr_max:.3f} ({n_r:,} cells)")
    print(f"               thresholds must be above {tr_max:.3f} and below {ts_min:.3f}")
    for cls in table:
        for k in ("theta*_s", "theta*_t"):
            v = cls[k]
            if not (tr_max < v < ts_min):
                safe = np.floor((ts_min - 1e-6) * 100) / 100       # largest 2-decimal value below the smallest thetaS
                flag = "--theta-s" if k == "theta*_s" else "--theta-t"
                raise SystemExit(f"class {cls['ID']} {k} = {v} is outside thetaR..thetaS ({tr_max:.3f}..{ts_min:.3f}); "
                                 f"tRIBS would stop at the first cell where it does not fit.\n"
                                 f"Largest 2-decimal value that fits: {safe:.2f}. Re-run with, for example:  {flag} {safe:.2f}")
    print("               both thresholds fit")


def check_ranges(table: list[dict]) -> None:
    """Same ranges the tRIBS 6.0.0 land table validator enforces."""
    for cls in table:
        for k in COLS:
            v = cls[k]
            if v == 9999.99:
                continue
            if k in UNIT_INTERVAL and not (0 <= v <= 1):
                raise SystemExit(f"class {cls['ID']}: {k} = {v} must be between 0 and 1")
            if v < 0:
                raise SystemExit(f"class {cls['ID']}: {k} = {v} must not be negative")


def report(table: list[dict]) -> None:
    print("\nID class          " + "".join(f"{c:>9}" for c in COLS))
    for cls in table:
        print(f"{cls['ID']:>2} {CLASS_NAMES[cls['ID']]:<15}" + "".join(f"{cls[c]:>9g}" for c in COLS))
    print("(RZD_m is left out of the dictionaries, so 9999.99 is written: the tRIBS default root zone depth)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="print the table and run the checks, write nothing")
    ap.add_argument("--theta-s", type=float, default=THETA_S_STAR, help=f"evaporation stress threshold (default {THETA_S_STAR})")
    ap.add_argument("--theta-t", type=float, default=THETA_T_STAR, help=f"transpiration stress threshold (default {THETA_T_STAR})")
    ap.add_argument("--buffer-m", type=float, default=200.0,
                    help="check the soil grids this far beyond the basin outline too (default 200 m; 0 = basin only)")
    ap.add_argument("--boundary", default=None, help="basin outline .shp (default: the one *boundary*.shp in lab/data/preprocessing)")
    ap.add_argument("--landuse", default=str(LAND_RAS))
    ap.add_argument("--out", default=str(OUT_LDT))
    args = ap.parse_args(argv)

    table = [{"ID": cid, **p, "theta*_s": args.theta_s, "theta*_t": args.theta_t}
             for cid, p in sorted(LAND_PARAMS.items())]
    check_map_ids(Path(args.landuse))
    check_ranges(table)
    check_thresholds(table, find_boundary(args.boundary), args.buffer_m)
    report(table)

    if args.dry_run:
        print("\n-- dry run, nothing written --")
        return 0

    from pytRIBS.shared.inout import InOut             # imported here so --dry-run works anywhere
    out = Path(args.out)
    existed = out.exists()
    out.parent.mkdir(parents=True, exist_ok=True)
    InOut.write_landuse_table(table, str(out))
    print(f"\n{'replaced' if existed else 'wrote'} {out}\n")
    print(out.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
