"""Shared helpers for the SRP meeting figures (Lower Sycamore Wash, WY2017)."""
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def _find_ws():
    """workspaces/Sycamore: the nearest folder above this script that contains init_data/ (e.g. when the
    scripts sit in workspaces/Sycamore/lab). Falls back to the sandbox copy of the repo."""
    for p in (HERE, *HERE.parents):
        if (p / "init_data").is_dir():
            return p
    return Path("/home/claude/tracylangston-byte/tribs-fall2026/workspaces/Sycamore")


WS = _find_ws()
INIT = WS / "init_data"
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
_b = [WS / "lab" / "data" / "preprocessing" / "Sycamore_boundary.shp", HERE / "boundary_in" / "Sycamore_boundary.shp"]
BOUNDARY = next((p for p in _b if p.exists()), _b[0])      # basin outline (NAD83 / UTM 12N shapefile)

GAUGE_XY = (449721.0, 3728390.0)      # USGS 09510200, EPSG:26912 (data_sources.md)

# --- look and feel (reference palette: light surface, ink, one blue) ---
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
BLUE = "#2a78d6"
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "axes.edgecolor": GRID, "axes.linewidth": 1.0, "axes.spines.top": False,
        "axes.spines.right": False, "grid.color": GRID, "grid.linewidth": 1.0,
        "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.titlecolor": INK, "legend.frameon": False, "savefig.dpi": 200,
    })


def read_asc(path):
    """ESRI ASCII grid -> (header dict, 2-D float array with NaN for nodata, extent L,R,B,T)."""
    hdr = {}
    n_hdr = 0
    with open(path) as f:
        while True:
            pos = f.tell()
            parts = f.readline().split()
            if len(parts) == 2 and parts[0].replace("_", "").isalpha():
                hdr[parts[0].lower()] = float(parts[1]); n_hdr += 1
            else:
                break
    a = pd.read_csv(path, sep=r"\s+", header=None, skiprows=n_hdr, dtype=np.float64,
                    engine="c").to_numpy()
    nd = hdr.get("nodata_value")
    if nd is not None:
        a = np.where(a == nd, np.nan, a)
    cs = hdr["cellsize"]
    L, B = hdr["xllcorner"], hdr["yllcorner"]
    R, T = L + hdr["ncols"] * cs, B + hdr["nrows"] * cs
    assert a.shape == (int(hdr["nrows"]), int(hdr["ncols"])), (a.shape, hdr)
    return hdr, a, (L, R, B, T)


def load_boundary(path):
    """Read polygon rings from a shapefile (pyshp). Returns list of (N,2) arrays in the file's CRS.
    Refuses files whose .prj is not NAD83 / UTM 12N, so a wrong-CRS outline cannot slip through."""
    import shapefile
    path = Path(path)
    prj = path.with_suffix(".prj")
    if prj.exists():
        w = prj.read_text()
        if "UTM zone 12N" not in w and "UTM_Zone_12N" not in w and "26912" not in w:
            raise SystemExit(f"{prj.name} is not NAD83 / UTM zone 12N; reproject the outline first.")
    else:
        print(f"warning: no .prj next to {path.name}; assuming EPSG:26912")
    sf = shapefile.Reader(str(path))
    rings = []
    for shp in sf.shapes():
        pts = np.array(shp.points)
        parts = list(shp.parts) + [len(pts)]
        for i in range(len(parts) - 1):
            rings.append(pts[parts[i]:parts[i + 1]])
    return rings


def inside_mask(rings, X, Y):
    """Boolean mask of points (X, Y arrays of cell centres) inside any ring (even-odd rule)."""
    from matplotlib.path import Path as MPath
    pts = np.column_stack([X.ravel(), Y.ravel()])
    m = np.zeros(len(pts), dtype=bool)
    for r in rings:
        m ^= MPath(r).contains_points(pts)
    return m.reshape(X.shape)


def grid_cover(rings, extent, shape, ss=4):
    """Fraction (0-1) of each cell of a north-up grid that lies inside the polygon.
    extent = (L, R, B, T); shape = (ny, nx). The polygon is drawn with PIL on a grid `ss` times finer
    (even-odd, so holes work) and block-averaged back down. PIL fills every pixel the outline touches,
    so the finer the sub-grid, the smaller that edge bias gets."""
    from PIL import Image, ImageDraw
    L, R, B, T = extent
    ny, nx = shape
    fine = np.zeros((ny * ss, nx * ss), dtype=bool)
    for r in rings:
        if len(r) < 3:
            continue
        px = (r[:, 0] - L) / (R - L) * nx * ss
        py = (T - r[:, 1]) / (T - B) * ny * ss
        img = Image.new("1", (nx * ss, ny * ss), 0)
        ImageDraw.Draw(img).polygon(list(zip(px.tolist(), py.tolist())), fill=1)
        fine ^= np.asarray(img, dtype=bool)
    return fine.reshape(ny, ss, nx, ss).mean(axis=(1, 3))


def grid_mask(rings, extent, shape, ss=4):
    """Boolean mask: cells at least half inside the polygon (close to the cell-centre rule)."""
    return grid_cover(rings, extent, shape, ss) >= 0.5


# The five WY2017 winter storms, shared by figure 1 and the event table. Each window is 72 hours and starts about
# 4 hours before the storm's first wet hour (hourly basin-mean radar rain above 0.1 mm), so the whole storm is inside it.
# (label, short name, window start, window end), Arizona time (MST).
EVENTS = [
    ("1", "Jan 1-2",  "2016-12-30 12:00", "2017-01-02 12:00"),
    ("2", "Jan 21",   "2017-01-19 06:00", "2017-01-22 06:00"),
    ("3", "Jan 24",   "2017-01-23 03:00", "2017-01-26 03:00"),
    ("4", "Feb 19",   "2017-02-17 21:00", "2017-02-20 21:00"),
    ("5", "Feb 28",   "2017-02-27 09:00", "2017-03-02 09:00"),
]
