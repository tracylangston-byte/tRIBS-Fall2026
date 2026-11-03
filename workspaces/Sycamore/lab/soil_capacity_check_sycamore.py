"""
soil_capacity_check_sycamore.py  --  Read-only look at the soil grids the Sycamore model uses (no model run, nothing written except to the screen).

QUESTION: on Jan 20 the gauges and radar show 0.5 to 4 mm/h of rain from noon to 19:00, and the real outlet rose from 3 to 84 m3/s while the
model made almost no flow. Can the model's soil take in 0.5 to 4 mm/h everywhere? If the soil's saturated conductivity (Ks) is far above the rain
rate in almost every cell, light rain cannot run off by infiltration excess and the early flow must come from somewhere else.

WHAT IT READS (only these; it never opens data/model, so it is safe while a run is building or running):
    ../init_data/scgrid.gdf          the list of SOLUS grids (KS, PB, TS, TR, PI, FD ...) -- the SOLUS files themselves, before any multiplier
    the grids that list names        ESRI ASCII files (.asc)
    optional --boundary <shp>        the basin outline (e.g. data/preprocessing/Sycamore_boundary.shp); without it every non-nodata cell counts

WHAT IT PRINTS: for each grid, the median and the 10th/90th percentile over the basin, as SOLUS gives it and as the working setting uses it
(Ks x0.19261, PsiB x0.2939, f x0.15 with the 1e-4 floor, as in the notebook). Then, for rain rates of 0.5 to 16 mm/h, the share of cells whose
in-use surface Ks is BELOW that rate.
LIMITS: Ks here is the surface (0-5 cm) value. Early in a storm the soil can take in MORE than Ks (the wetting-front suction helps), and the
capacity falls toward Ks as the soil wets, so "rain above Ks" is a floor for runoff, not the full answer.

Run (from lab/):   python soil_capacity_check_sycamore.py --boundary <path to Sycamore_boundary.shp>
"""
import argparse, csv, os
import numpy as np

NAMES = {"KS": "Ks 0-5 cm (mm/hr)", "PB": "PsiB (mm, negative)", "TS": "thetaS", "TR": "thetaR", "PI": "m (pore-size index)", "FD": "f (1/mm)", "PO": "porosity"}
F_FLOOR = 1e-4


def read_ascii(path):
    head = {}
    with open(path) as fh:
        for _ in range(6):
            k, v = fh.readline().split()
            head[k.lower()] = float(v)
        data = np.loadtxt(fh)
    return head, data.reshape(int(head["nrows"]), int(head["ncols"]))


def inside_mask(head, shape, boundary):
    """True for cells whose centre is inside the basin outline (needs geopandas + rasterio)."""
    import geopandas as gpd
    from rasterio.features import geometry_mask
    from affine import Affine
    cs = head["cellsize"]
    x0 = head["xllcorner"] if "xllcorner" in head else head["xllcenter"] - cs / 2
    y0 = head["yllcorner"] if "yllcorner" in head else head["yllcenter"] - cs / 2
    tr = Affine(cs, 0, x0, 0, -cs, y0 + cs * shape[0])
    shapes = list(gpd.read_file(boundary).geometry)
    return geometry_mask(shapes, out_shape=shape, transform=tr, invert=True)


def pct(a, q):
    return float(np.percentile(a, q))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--init-data", default="../init_data")
    ap.add_argument("--boundary", default=None)
    ap.add_argument("--ks-mult", type=float, default=0.19261)
    ap.add_argument("--psib-mult", type=float, default=0.2939)
    ap.add_argument("--f-mult", type=float, default=0.15)
    ap.add_argument("--rain", default="0.5,1,2,4,8,16")
    a = ap.parse_args()
    mult = {"KS": a.ks_mult, "PB": a.psib_mult, "FD": a.f_mult}
    with open(os.path.join(a.init_data, "scgrid.gdf")) as fh:
        rows = list(csv.DictReader(fh))
    print(f"grids listed in {a.init_data}/scgrid.gdf: {[r['Variable'].strip() for r in rows]}")
    vals, used, mask = {}, {}, None
    for r in rows:
        v = r["Variable"].strip()
        path = f"{r['BasePath'].strip()}.{r['FileExtension'].strip()}"
        head, grid = read_ascii(path)
        nd = head.get("nodata_value", -9999.0)
        good = grid != nd
        if a.boundary and mask is None:
            mask = inside_mask(head, grid.shape, a.boundary)
        keep = good & (mask if mask is not None else True)
        x = grid[keep].astype(float)
        vals[v] = x
        u = x * mult.get(v, 1.0)
        if v == "FD":
            u = np.maximum(u, F_FLOOR)
        used[v] = u
        print(f"  {v}: {path}  cells used {x.size}")
    print(f"\nbasin outline used: {'yes (' + a.boundary + ')' if mask is not None else 'NO - every non-nodata cell counts'}")
    print(f"\n{'grid':<24}{'SOLUS p10':>11}{'median':>10}{'p90':>10} | {'in use p10':>11}{'median':>10}{'p90':>10}")
    for v, x in vals.items():
        u = used[v]
        print(f"{NAMES.get(v, v):<24}{pct(x,10):>11.4g}{pct(x,50):>10.4g}{pct(x,90):>10.4g} | {pct(u,10):>11.4g}{pct(u,50):>10.4g}{pct(u,90):>10.4g}")
    if "KS" in used:
        ks, ks0 = used["KS"], vals["KS"]
        print("\nShare of basin cells whose surface Ks is BELOW the rain rate (rain above Ks can run off by infiltration excess):")
        print(f"{'rain mm/h':>10}{'in-use Ks':>12}{'SOLUS Ks':>12}")
        for rate in [float(s) for s in a.rain.split(",")]:
            print(f"{rate:>10g}{(ks < rate).mean():>12.1%}{(ks0 < rate).mean():>12.1%}")


if __name__ == "__main__":
    main()
