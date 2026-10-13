#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_mrms_radar_rainfall_tribs6.py
===================================

Build hourly gridded radar rainfall input for tRIBS 6.0.0 from NOAA MRMS.

WHAT THIS SCRIPT DOES
---------------------
1. DOWNLOAD  - fetches hourly MRMS rainfall, one file per hour, from one of two
               public archives (choose with --source):
                 aws  NOAA bucket "noaa-mrms-pds": MultiSensor QPE Pass 2 (the
                      default). The archive starts 2020-10-14.
                 iem  Iowa State "MTArchive": GaugeCorr QPE, the older
                      gauge-corrected radar product. Use it for water years
                      before Oct 2020 (files confirmed to exist for 2017).
               Each file is a gzipped GRIB2 covering the whole US (a few MB), so the
               script immediately CLIPS it to a small window around your basin
               and (by default) deletes the big file. Only the small clips are
               kept, so a full water year fits comfortably in a Codespace.
2. PROCESS   - reprojects each clip onto a regular square-cell grid in your
               model's projected CRS (e.g. EPSG:26912), converts UTC to your
               local model clock, and writes one ESRI ASCII grid per hour named
               <basename><MMDDYYYYHH>.asc  (e.g. radar_0601200203.asc).
               It also writes a manifest CSV and a printable snippet of the
               tRIBS .in keywords you need.
3. AUDIT     - re-reads the files it wrote and checks them (file count, headers,
               NaN/negatives, spikes, coverage of the basin, zero-filled hours),
               writes a report + basin-mean hyetograph, and optionally compares
               against a rain gauge (including a lag test).

QUICK START
-----------
    # 0. one-time environment check (packages, GRIB support, network)
    python build_mrms_radar_rainfall_tribs6.py check-env

    # 1. see the plan without downloading anything (also prints the
    #    local-time -> MRMS-UTC mapping so you can eyeball it)
    python build_mrms_radar_rainfall_tribs6.py download --dry-run \\
        --basin path/to/sycamore_basin.shp \\
        --start "2023-10-01 00:00" --end "2023-10-07 23:00"

    # 2. try ONE WEEK first, end to end
    python build_mrms_radar_rainfall_tribs6.py all \\
        --basin path/to/sycamore_basin.shp \\
        --start "2023-10-01 00:00" --end "2023-10-07 23:00"

    # 3. then the full water year
    python build_mrms_radar_rainfall_tribs6.py all \\
        --basin path/to/sycamore_basin.shp \\
        --start "2023-10-01 00:00" --end "2024-09-30 23:00"

    # 3b. water years before Oct 2020: the same commands plus  --source iem
    #     (WY2017 = --start "2016-10-01 00:00" --end "2017-09-30 23:00")

    Every stage can also be run alone: download | process | audit.
    Use the SAME flags for each stage (or edit USER_SETTINGS below once).
    Downloads are resumable: re-running skips hours already clipped.

TIME CONVENTIONS  (READ THIS - it is the easiest place to be off by an hour)
---------------------------------------------------------------------------
* MRMS files are stamped in UTC with the time at the END of the accumulation
  hour. MRMS_..._20230715-180000 is the rain that fell 17:00-18:00 UTC.
  (This is the assumed convention for the AWS product. For the Iowa State
  GaugeCorr product it is also UNVERIFIED: the GRIB metadata does not say.
  Check it with the gauge lag test in the audit on a real storm.)
* Arizona has no daylight saving: local time = UTC - 7 h all year. The
  --start/--end you give are LOCAL model-clock hours, inclusive.
* The hour in a tRIBS filename (the HH in <base>MMDDYYYYHH.asc) is a LABEL.
  Whether tRIBS applies file "..HH" to the interval starting at HH or ending
  at HH is NOT documented anywhere I could find, so it is a setting:
      --stamp-convention start   label = START of the interval (DEFAULT)
                                  local 00:00 label <- MRMS 01:00 local stamp
      --stamp-convention end     label = END of the interval
                                  local 01:00 label <- MRMS 01:00 local stamp
  The default 'start' is my best inference from how a model reads "the file for
  the current time and applies it forward", but it is UNVERIFIED. Confirm with
  Josh Cederstrom, or with the one-week test in the audit output (compare the
  hyetograph to what tRIBS reports, and use the gauge lag test).
  The convention is recorded in domain.json and manifest.csv. Changing it only
  needs the PROCESS stage re-run (no re-download).

OUTPUT LAYOUT (under --workdir, default ./mrms_work)
----------------------------------------------------
    domain.json                 grid definition + all settings used
    clips/<product>/*.tif       small per-hour clips in lon/lat (UTC-stamped)
    raw/                        full-CONUS .grib2.gz kept here ONLY with --keep-raw
    download_log.csv            per-hour download result
    rain/<base>MMDDYYYYHH.asc   THE FILES tRIBS READS (folder = --outdir)
    manifest.csv                per-hour provenance + statistics
    tribs_input_snippet.txt     the .in keywords to copy into your model
    audit/                      audit_report.txt, hyetograph CSV, plots

DESIGN CHOICES WORTH KNOWING
----------------------------
* Units: MRMS QPE_01H is mm accumulated in one hour, which is numerically
  mm/hr. tRIBS 6.0.0 wants mm/hr. No conversion is applied.
* Gridded rain is HOURLY ONLY in tRIBS 6.0.0. Hourly averaging smooths
  convective peaks; expect lower peak intensities than a 5-minute gauge.
* Missing or "no coverage" cells (MRMS uses negative sentinel values) are
  filled with 0 by default and COUNTED (manifest + audit). Missing whole hours
  are written as all-zero grids by default and flagged. Use --on-missing-cells
  error / --on-missing-hours error|skip to change this. Nothing is silently
  invented: every filled hour is listed in the manifest and audit.
* The tRIBS NODATA value (-9999) appears only in the header; no data cell is
  ever -9999, because we don't know how tRIBS treats it as a rain rate.
* Resampling defaults to 'average' (area-weighted; conserves mean rain depth).
* The script does not use pytRIBS: pytRIBS 1.0.0 has no rainfall-grid writer
  that I could find, so the ASCII files are written directly.

REQUIRES: numpy, pandas, requests, rasterio (with GRIB driver; the pip wheel
has it). geopandas is needed only if you pass --basin (a polygon file);
--bbox W S E N needs no geopandas. matplotlib is optional (plots).

VERSION HISTORY
---------------
v1.0 (2026-09-29)  first version. Tested on synthetic MRMS-style GRIB2 and a
                   mock S3 server; NOT yet run against the real AWS bucket.
v1.1 (2026-10-09)  added --source iem (Iowa State archive, GaugeCorr product) for
                   water years before Oct 2020. Downloads are now recognised
                   whether or not the server already decompressed them (the
                   first bytes are checked, not the .gz name). --s3-base-url
                   is now --base-url (the old name still works). Tested on
                   synthetic GRIB2 and a mock server for both sources; NOT yet
                   run against the real Iowa State or AWS servers.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import shutil
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import requests

SCRIPT_VERSION = "1.1"

# =============================================================================
# USER SETTINGS  -  edit these once, or override any of them on the command
# line (e.g. --basin file.shp).  None means "must be given".
# =============================================================================
USER_SETTINGS = {
    # --- where things live -------------------------------------------------
    "workdir": "mrms_work",        # everything is created under here
    "outdir": None,                # None -> <workdir>/rain  (the .asc files)
    "basename": "radar_",          # file prefix; becomes RAINFILE base name
    "extension": "asc",            # RAINEXTENSION (no period)
    # --- where is the basin (give ONE of these) ------------------------------
    "basin": None,                 # polygon file (.shp/.geojson/.gpkg), any CRS
    "bbox": None,                  # or (west, south, east, north) in lon/lat
    # --- target grid --------------------------------------------------------
    "out_crs": "EPSG:26912",       # must match your tRIBS mesh CRS (UTM 12N)
    "cellsize": 1000.0,            # metres; MRMS is ~0.93 x 1.11 km here
    "buffer_km": 5.0,              # padding around the basin
    "resampling": "average",       # average | bilinear | nearest
    # --- time ---------------------------------------------------------------
    "start": None,                 # LOCAL clock, first hour label, inclusive
    "end": None,                   # LOCAL clock, last hour label, inclusive
    "utc_offset": -7.0,            # hours from UTC; Arizona = -7 all year
    "stamp_convention": "start",   # start | end   (see TIME CONVENTIONS)
    # --- MRMS source --------------------------------------------------------
    "source": "aws",               # aws (Oct 2020 on) | iem (older years)
    "product": None,               # None -> the source's default (see SOURCES)
    "fallback_product": None,      # e.g. "MultiSensor_QPE_01H_Pass1_00.00"
    "base_url": None,              # None -> the source's default (see SOURCES)
    "region": "CONUS",             # aws only
    # --- behaviour ----------------------------------------------------------
    "on_missing_cells": "zero",    # zero | error
    "on_missing_hours": "zero",    # zero | skip | error
    "workers": 4,                  # parallel downloads
    "retries": 4,
    "retry_wait": 2.0,             # seconds; doubles each retry
    "timeout": 60,                 # seconds per request
}

MRMS_MAX_VALID_MM = 1000.0     # anything above is treated as a bad value
NODATA_HEADER_VALUE = -9999    # only ever written in the header
CLIP_MARGIN_DEG = 0.03         # extra lon/lat margin around the target grid

_GRIB_LOCK = threading.Lock()  # GDAL's GRIB reader is not reliably thread-safe
_TLS = threading.local()

# The two places hourly MRMS files can come from. "probe_stamp" is an hour that is
# known to exist, used only by check-env.
SOURCES = {
    "aws": {"label": "NOAA MRMS archive on AWS (noaa-mrms-pds), starts 2020-10-14",
            "product": "MultiSensor_QPE_01H_Pass2_00.00",
            "base_url": "https://noaa-mrms-pds.s3.amazonaws.com",
            "probe_stamp": datetime(2023, 7, 15, 18)},
    "iem": {"label": "Iowa State IEM MTArchive (NCEP MRMS mirror), GaugeCorr product",
            "product": "GaugeCorr_QPE_01H_00.00",
            "base_url": "https://mtarchive.geol.iastate.edu",
            "probe_stamp": datetime(2017, 8, 15, 18)},
}


# =============================================================================
# small helpers
# =============================================================================
def die(msg: str, code: int = 2) -> None:
    """Print a plain-English error and stop the script."""
    print(f"\nERROR: {msg}\n", file=sys.stderr)
    sys.exit(code)


def parse_local_hour(text: str, name: str) -> datetime:
    """Parse '2023-10-01 00:00' (local clock) and insist on a whole hour."""
    ts = pd.Timestamp(text)
    if ts.tzinfo is not None:
        die(f"{name} must be a plain local time without a timezone, got {text!r}")
    if ts.minute or ts.second or ts.microsecond:
        die(f"{name} must be exactly on the hour, got {text!r}")
    return ts.to_pydatetime()


def hour_labels(start: datetime, end: datetime) -> list:
    """All local hour labels from start to end, inclusive."""
    if end < start:
        die(f"--end ({end}) is before --start ({start})")
    return [t.to_pydatetime() for t in pd.date_range(start, end, freq="h")]


def label_to_mrms_utc(label: datetime, utc_offset_hours: float, convention: str) -> datetime:
    """Convert a local tRIBS hour label to the MRMS file stamp (UTC).

    MRMS stamps mark the END of the accumulation hour.
      convention 'start': label is the interval START, so interval end = label + 1 h
      convention 'end'  : label is the interval END,   so interval end = label
    UTC = local - utc_offset  (Arizona: offset -7, so UTC = local + 7 h).
    """
    if convention == "start":
        interval_end_local = label + timedelta(hours=1)
    elif convention == "end":
        interval_end_local = label
    else:
        raise ValueError(f"unknown stamp convention {convention!r}")
    return interval_end_local - timedelta(hours=utc_offset_hours)


def asc_filename(basename: str, ext: str, label: datetime) -> str:
    """tRIBS gridded-rain file name: <base><MMDDYYYYHH>.<ext>."""
    return f"{basename}{label:%m%d%Y%H}.{ext}"


def s3_key(region: str, product: str, stamp_utc: datetime) -> str:
    day = f"{stamp_utc:%Y%m%d}"
    return f"{region}/{product}/{day}/MRMS_{product}_{day}-{stamp_utc:%H%M%S}.grib2.gz"


def source_url(source: str, base_url: str, region: str, product: str, stamp_utc: datetime) -> str:
    """Full download URL of one MRMS hour for the chosen source.

    aws: <base>/CONUS/<product>/<YYYYMMDD>/MRMS_<product>_<YYYYMMDD>-<HHMMSS>.grib2.gz
    iem: <base>/<YYYY>/<MM>/<DD>/mrms/ncep/<product folder>/<product>_<YYYYMMDD>-<HHMMSS>.grib2.gz
         (the folder name is the product without its trailing "_00.00")
    """
    base = base_url.rstrip("/")
    day = f"{stamp_utc:%Y%m%d}"
    if source == "aws":
        return f"{base}/{s3_key(region, product, stamp_utc)}"
    if source == "iem":
        folder = re.sub(r"_\d\d\.\d\d$", "", product)
        return (f"{base}/{stamp_utc:%Y/%m/%d}/mrms/ncep/{folder}/"
                f"{product}_{day}-{stamp_utc:%H%M%S}.grib2.gz")
    raise ValueError(f"unknown source {source!r}")


def resolve_source(args) -> None:
    """Fill in the product and base URL that depend on --source."""
    defaults = SOURCES[args.source]
    if args.product is None:
        args.product = defaults["product"]
    if args.base_url is None:
        args.base_url = defaults["base_url"]
    if args.source == "iem" and not re.search(r"_\d\d\.\d\d$", args.product):
        args.product += "_00.00"      # IEM file names carry the level suffix


def clip_path_for(clips_dir: Path, product: str, stamp_utc: datetime) -> Path:
    return clips_dir / product / f"mrms_{stamp_utc:%Y%m%d_%H}.tif"


def fmt_num(x: float) -> str:
    """Compact number for the ASCII header (400000 not 400000.000)."""
    s = f"{x:.3f}".rstrip("0").rstrip(".")
    return s if s else "0"


# =============================================================================
# grid / domain definition
# =============================================================================
@dataclass
class Grid:
    """A regular square-cell grid in a projected CRS (north-up)."""
    crs: str
    cellsize: float
    xmin: float
    ymin: float
    ncols: int
    nrows: int

    @property
    def xmax(self) -> float:
        return self.xmin + self.ncols * self.cellsize

    @property
    def ymax(self) -> float:
        return self.ymin + self.nrows * self.cellsize

    @property
    def transform(self):
        from rasterio.transform import from_origin
        return from_origin(self.xmin, self.ymax, self.cellsize, self.cellsize)


def load_basin_geometry(path: str, out_crs: str):
    """Read a basin polygon file and return ONE geometry in the model CRS."""
    try:
        import geopandas as gpd
    except ImportError:
        die("--basin needs geopandas (pip install geopandas). "
            "Alternatively give --bbox WEST SOUTH EAST NORTH in lon/lat.")
    if not Path(path).exists():
        die(f"basin file not found: {path}")
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        die(f"basin file {path} has no CRS defined; re-save it with one")
    gdf = gdf.to_crs(out_crs)
    union = getattr(gdf.geometry, "union_all", None)
    return union() if union is not None else gdf.geometry.unary_union


def build_grid(basin, bbox, out_crs: str, cellsize: float, buffer_km: float):
    """Make the target grid (snapped to multiples of cellsize) + basin geometry."""
    from rasterio.warp import transform_bounds
    geom = None
    if basin:
        geom = load_basin_geometry(basin, out_crs)
        xmin, ymin, xmax, ymax = geom.bounds
    elif bbox:
        west, south, east, north = bbox
        xmin, ymin, xmax, ymax = transform_bounds(
            "EPSG:4326", out_crs, west, south, east, north, densify_pts=21)
    else:
        die("give the basin location: --basin <polygon file>  or  "
            "--bbox WEST SOUTH EAST NORTH (lon/lat)")
    buf = buffer_km * 1000.0
    xmin, ymin, xmax, ymax = xmin - buf, ymin - buf, xmax + buf, ymax + buf
    gx0 = math.floor(xmin / cellsize) * cellsize
    gy0 = math.floor(ymin / cellsize) * cellsize
    gx1 = math.ceil(xmax / cellsize) * cellsize
    gy1 = math.ceil(ymax / cellsize) * cellsize
    grid = Grid(out_crs, float(cellsize), float(gx0), float(gy0),
                int(round((gx1 - gx0) / cellsize)), int(round((gy1 - gy0) / cellsize)))
    return grid, geom


def clip_bbox_lonlat(grid: Grid):
    """Lon/lat window (with margin) that fully covers the target grid."""
    from rasterio.warp import transform_bounds
    w, s, e, n = transform_bounds(grid.crs, "EPSG:4326", grid.xmin, grid.ymin,
                                  grid.xmax, grid.ymax, densify_pts=21)
    m = CLIP_MARGIN_DEG
    return (w - m, s - m, e + m, n + m)


def basin_mask(geom, grid: Grid):
    """Boolean array (nrows, ncols): True where a cell centre is in the basin."""
    if geom is None:
        return np.ones((grid.nrows, grid.ncols), dtype=bool)
    from rasterio import features
    m = features.rasterize([(geom, 1)], out_shape=(grid.nrows, grid.ncols),
                           transform=grid.transform, fill=0, all_touched=False,
                           dtype="uint8").astype(bool)
    if not m.any():  # basin smaller than a cell: fall back to touched cells
        m = features.rasterize([(geom, 1)], out_shape=(grid.nrows, grid.ncols),
                               transform=grid.transform, fill=0, all_touched=True,
                               dtype="uint8").astype(bool)
    return m


# =============================================================================
# run description shared by all stages
# =============================================================================
@dataclass
class Dirs:
    work: Path
    clips: Path
    raw: Path
    out: Path
    audit: Path


def make_dirs(args) -> Dirs:
    work = Path(args.workdir).resolve()
    out = Path(args.outdir).resolve() if args.outdir else work / "rain"
    return Dirs(work, work / "clips", work / "raw", out, work / "audit")


def print_time_mapping(labels, args, n: int = 3) -> None:
    """Show how local labels map to MRMS stamps - eyeball this!"""
    print("\nTime mapping (stamp convention = "
          f"'{args.stamp_convention}', UTC offset = {args.utc_offset:+g} h)")
    print(f"  {'tRIBS label (local)':<20} {'rain interval (local)':<27} "
          f"{'MRMS stamp (UTC)':<18} file name")
    if len(labels) > 2 * n:
        show = labels[:n] + [None] + labels[-n:]   # None prints a "..." row
    else:
        show = labels
    for lab in show:
        if lab is None:
            print("  ...")
            continue
        end_local = lab + timedelta(hours=1) if args.stamp_convention == "start" else lab
        start_local = end_local - timedelta(hours=1)
        stamp = label_to_mrms_utc(lab, args.utc_offset, args.stamp_convention)
        print(f"  {lab:%Y-%m-%d %H:%M}        "
              f"{start_local:%m-%d %H:%M} -> {end_local:%m-%d %H:%M}      "
              f"{stamp:%Y-%m-%d %H:%M}    "
              f"{asc_filename(args.basename, args.extension, lab)}")
    print("  MRMS stamps mark the END of each accumulation hour.\n")


def stamps_needed(labels, args) -> list:
    """Unique sorted MRMS UTC stamps needed for these labels."""
    return sorted({label_to_mrms_utc(l, args.utc_offset, args.stamp_convention)
                   for l in labels})


# =============================================================================
# STAGE 1: DOWNLOAD + CLIP
# =============================================================================
def _session() -> requests.Session:
    """One requests.Session per thread (Sessions aren't guaranteed thread-safe)."""
    if not hasattr(_TLS, "session"):
        s = requests.Session()
        s.headers["User-Agent"] = f"mrms-tribs6-script/{SCRIPT_VERSION}"
        _TLS.session = s
    return _TLS.session


def clip_grib_to_tif(grib_path: Path, bbox, out_tif: Path, tags: dict) -> None:
    """Read a full-CONUS MRMS GRIB2, cut out `bbox`, mask bad values, save a GeoTIFF.

    Bad values (MRMS uses negatives such as -3 for 'no coverage', or huge fill
    values) become NaN so later stages never see them as rain.
    """
    import rasterio
    from rasterio.transform import Affine
    from rasterio.windows import Window
    with _GRIB_LOCK:
        with rasterio.open(grib_path) as src:
            if src.crs is None or not src.crs.is_geographic:
                raise RuntimeError(f"unexpected MRMS grid CRS: {src.crs}")
            transform = src.transform
            if transform.c > 180.0:  # 0-360 longitudes -> -180..180
                transform = Affine.translation(-360.0, 0.0) * transform
            west, south, east, north = bbox
            col0 = max(0, math.floor((west - transform.c) / transform.a))
            col1 = min(src.width, math.ceil((east - transform.c) / transform.a))
            row0 = max(0, math.floor((north - transform.f) / transform.e))
            row1 = min(src.height, math.ceil((south - transform.f) / transform.e))
            if col1 <= col0 or row1 <= row0:
                raise RuntimeError("basin window lies outside the MRMS grid")
            win = Window(col0, row0, col1 - col0, row1 - row0)
            arr = src.read(1, window=win).astype("float32")
            nodata = src.nodata
            win_transform = transform * Affine.translation(col0, row0)
    bad = ~np.isfinite(arr) | (arr < 0) | (arr > MRMS_MAX_VALID_MM)
    if nodata is not None and np.isfinite(nodata):
        bad |= (arr == nodata)
    arr[bad] = np.nan
    out_tif.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_tif.with_suffix(".tmp.tif")
    with rasterio.open(tmp, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype="float32", crs="EPSG:4326", transform=win_transform,
                       nodata=float("nan"), compress="deflate") as dst:
        dst.write(arr, 1)
        dst.update_tags(**{k: str(v) for k, v in tags.items()})
    tmp.replace(out_tif)


def fetch_and_clip(args, dirs: Dirs, product: str, stamp: datetime, bbox) -> dict:
    """Download one MRMS hour, verify it, clip it, and clean up.

    Returns a dict describing what happened (status: ok | cached | missing | error).
    """
    result = {"mrms_stamp_utc": f"{stamp:%Y-%m-%d %H:%M}", "product": product,
              "status": "", "http": "", "note": "", "mb": 0.0}
    clip = clip_path_for(dirs.clips, product, stamp)
    if clip.exists() and not args.force:
        result["status"] = "cached"
        return result

    url = source_url(args.source, args.base_url, args.region, product, stamp)
    raw_dir = dirs.raw / product
    raw_dir.mkdir(parents=True, exist_ok=True)
    gz_path = raw_dir / Path(url).name
    part = gz_path.with_suffix(gz_path.suffix + ".part")

    last_note = ""
    for attempt in range(args.retries + 1):
        if attempt:
            time.sleep(args.retry_wait * (2 ** (attempt - 1)))
        try:
            with _session().get(url, stream=True, timeout=args.timeout) as r:
                result["http"] = r.status_code
                if r.status_code in (403, 404):
                    # Public buckets return 404 for a key that doesn't exist.
                    # 403 is treated the same but flagged, because it can also
                    # mean the network/proxy is blocking access.
                    result["status"] = "missing"
                    result["note"] = f"HTTP {r.status_code}"
                    return result
                if r.status_code != 200:
                    last_note = f"HTTP {r.status_code}"
                    continue
                with open(part, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
            part.replace(gz_path)
            result["mb"] = round(gz_path.stat().st_size / 1e6, 2)

            # Decompress fully: this also verifies the gzip checksum. Some
            # servers decompress on the fly, so look at the first bytes
            # instead of trusting the .gz in the file name.
            with tempfile.NamedTemporaryFile(dir=raw_dir, suffix=".grib2", delete=False) as tmp:
                tmp_path = Path(tmp.name)
            try:
                with open(gz_path, "rb") as probe:
                    magic = probe.read(4)
                if magic[:2] == b"\x1f\x8b":
                    with gzip.open(gz_path, "rb") as fin, open(tmp_path, "wb") as fout:
                        shutil.copyfileobj(fin, fout)
                elif magic == b"GRIB":
                    shutil.copyfile(gz_path, tmp_path)
                else:
                    raise OSError(f"download is neither gzip nor GRIB (first bytes {magic!r}); "
                                  "the server may have sent an error page")
                clip_grib_to_tif(tmp_path, bbox, clip,
                                 {"source_url": url, "mrms_stamp_utc": f"{stamp:%Y-%m-%dT%H:%M:%SZ}",
                                  "product": product})
            finally:
                tmp_path.unlink(missing_ok=True)
            if not args.keep_raw:
                gz_path.unlink(missing_ok=True)
            result["status"] = "ok"
            return result
        except (requests.RequestException, OSError, EOFError, gzip.BadGzipFile,
                RuntimeError) as exc:
            last_note = f"{type(exc).__name__}: {exc}"
            gz_path.unlink(missing_ok=True)
            part.unlink(missing_ok=True)
            # A clip-stage RuntimeError (bad grid/window) will not fix itself.
            if isinstance(exc, RuntimeError):
                break
    result["status"] = "error"
    result["note"] = last_note
    return result


def stage_download(args, dirs: Dirs, grid: Grid, labels) -> pd.DataFrame:
    if args.source == "iem" and args.workers > 4:
        print("NOTE: the Iowa State archive is a free university service; "
              "limiting to 4 parallel downloads.")
        args.workers = 4
    stamps = stamps_needed(labels, args)
    bbox = clip_bbox_lonlat(grid)
    print(f"DOWNLOAD: {len(stamps)} MRMS hours of '{args.product}' "
          f"({stamps[0]:%Y-%m-%d %H:%M} to {stamps[-1]:%Y-%m-%d %H:%M} UTC)")
    print(f"  clip window (lon/lat): W {bbox[0]:.3f}  S {bbox[1]:.3f}  E {bbox[2]:.3f}  N {bbox[3]:.3f}")
    print(f"  source: {args.source} - {SOURCES[args.source]['label']}")
    print(f"  base URL: {args.base_url}   workers: {args.workers}")
    if args.dry_run:
        print("  --dry-run: nothing downloaded.")
        for tag, s in (("first", stamps[0]), ("last ", stamps[-1])):
            print(f"  {tag} file: {source_url(args.source, args.base_url, args.region, args.product, s)}")
        return pd.DataFrame()

    dirs.clips.mkdir(parents=True, exist_ok=True)
    rows, t0, done = [], time.time(), 0
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futs = {pool.submit(fetch_and_clip, args, dirs, args.product, s, bbox): s for s in stamps}
        for fut in as_completed(futs):
            rows.append(fut.result())
            done += 1
            if done % 100 == 0 or done == len(stamps):
                rate = done / max(time.time() - t0, 1e-9)
                print(f"  {done}/{len(stamps)} done  ({rate:.1f} files/s, "
                      f"~{(len(stamps) - done) / max(rate, 1e-9) / 60:.1f} min left)")
    log = pd.DataFrame(rows)

    # Optional fallback product for hours the primary product doesn't have.
    if args.fallback_product:
        need = log[log["status"].isin(["missing", "error"])]["mrms_stamp_utc"].tolist()
        if need:
            print(f"  trying fallback '{args.fallback_product}' for {len(need)} hours ...")
            with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
                futs = [pool.submit(fetch_and_clip, args, dirs, args.fallback_product,
                                    datetime.strptime(s, "%Y-%m-%d %H:%M"), bbox) for s in need]
                fb_rows = [f.result() for f in futs]
            log = pd.concat([log, pd.DataFrame(fb_rows)], ignore_index=True)

    log = log.sort_values(["mrms_stamp_utc", "product"]).reset_index(drop=True)
    dirs.work.mkdir(parents=True, exist_ok=True)
    log.to_csv(dirs.work / "download_log.csv", index=False)

    # Summary + sanity checks on the results.
    counts = log[log["product"] == args.product]["status"].value_counts().to_dict()
    print(f"  result for primary product: {counts}")
    prim = log[log["product"] == args.product]
    got = prim["status"].isin(["ok", "cached"]).sum()
    if got == 0:
        http = prim["http"].astype(str).value_counts().to_dict()
        die("no MRMS files could be downloaded (HTTP codes seen: "
            f"{http}). If these are 403/407, your network is blocking "
            f"{args.base_url}; run 'check-env' to test. If they are all 404, "
            "the dates may be outside this source's archive (aws starts "
            "2020-10-14; iem GaugeCorr starts in 2015). "
            "Log: " + str(dirs.work / "download_log.csv"))
    n_bad = (prim["status"] == "error").sum()
    if n_bad:
        print(f"  WARNING: {n_bad} hours failed after retries (see download_log.csv). "
              "Re-running the stage will retry them.")
    n_miss = (prim["status"] == "missing").sum()
    if n_miss:
        print(f"  NOTE: {n_miss} hours are not in the archive (see download_log.csv); "
              "PROCESS will handle them per --on-missing-hours.")
    return log


# =============================================================================
# STAGE 2: PROCESS (reproject -> ASCII grids -> manifest -> .in snippet)
# =============================================================================
def write_ascii(path: Path, arr: np.ndarray, grid: Grid) -> None:
    """Write an ESRI ASCII grid (6-line header, rows top to bottom)."""
    header = (f"ncols {grid.ncols}\nnrows {grid.nrows}\n"
              f"xllcorner {fmt_num(grid.xmin)}\nyllcorner {fmt_num(grid.ymin)}\n"
              f"cellsize {fmt_num(grid.cellsize)}\nNODATA_value {NODATA_HEADER_VALUE}\n")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        f.write(header)
        np.savetxt(f, arr, fmt="%.3f")
    tmp.replace(path)


def reproject_clip(clip_path: Path, grid: Grid, resampling_name: str) -> np.ndarray:
    """Reproject one lon/lat clip onto the target grid; NaN where no valid data."""
    import rasterio
    from rasterio.warp import Resampling, reproject
    resampling = {"average": Resampling.average, "bilinear": Resampling.bilinear,
                  "nearest": Resampling.nearest}[resampling_name]
    with rasterio.open(clip_path) as src:
        arr = src.read(1).astype("float32")
        src_transform, src_crs = src.transform, src.crs
    dst = np.full((grid.nrows, grid.ncols), np.nan, dtype="float32")
    reproject(source=arr, destination=dst, src_transform=src_transform, src_crs=src_crs,
              dst_transform=grid.transform, dst_crs=grid.crs,
              src_nodata=np.nan, dst_nodata=np.nan, resampling=resampling)
    return dst


def find_clip(dirs: Dirs, args, stamp: datetime):
    """Return (source_label, path) for the best available clip, or (None, None)."""
    p = clip_path_for(dirs.clips, args.product, stamp)
    if p.exists():
        return "primary", p
    if args.fallback_product:
        p2 = clip_path_for(dirs.clips, args.fallback_product, stamp)
        if p2.exists():
            return "fallback", p2
    return None, None


def stage_process(args, dirs: Dirs, grid: Grid, geom, labels) -> pd.DataFrame:
    mask = basin_mask(geom, grid)
    n_basin = int(mask.sum())
    print(f"PROCESS: writing {len(labels)} hourly grids to {dirs.out}")
    print(f"  grid: {grid.ncols} x {grid.nrows} cells of {grid.cellsize:g} m, CRS {grid.crs}")
    print(f"  extent: x {grid.xmin:.0f}-{grid.xmax:.0f}, y {grid.ymin:.0f}-{grid.ymax:.0f}; "
          f"basin cells: {n_basin}" + ("" if geom is not None else "  (no basin polygon: stats use whole grid)"))
    if args.dry_run:
        print("  --dry-run: nothing written.")
        return pd.DataFrame()

    dirs.out.mkdir(parents=True, exist_ok=True)
    zero = np.zeros((grid.nrows, grid.ncols), dtype="float32")
    rows = []
    for i, lab in enumerate(labels):
        stamp = label_to_mrms_utc(lab, args.utc_offset, args.stamp_convention)
        src_kind, clip = find_clip(dirs, args, stamp)
        fname = asc_filename(args.basename, args.extension, lab)
        row = {"label_local": f"{lab:%Y-%m-%d %H:%M}",
               "interval_start_local": f"{(lab if args.stamp_convention == 'start' else lab - timedelta(hours=1)):%Y-%m-%d %H:%M}",
               "mrms_stamp_utc": f"{stamp:%Y-%m-%d %H:%M}", "file": fname}
        if clip is None:
            if args.on_missing_hours == "error":
                die(f"no MRMS data for label {lab} (UTC stamp {stamp}); "
                    "use --on-missing-hours zero|skip to continue past gaps")
            if args.on_missing_hours == "skip":
                row.update(source="none", status="skipped_missing_hour",
                           n_invalid_domain=grid.ncols * grid.nrows, n_invalid_basin=n_basin,
                           domain_mean_mm=np.nan, basin_mean_mm=np.nan, basin_max_mm=np.nan,
                           domain_max_mm=np.nan)
                rows.append(row)
                continue
            data = zero.copy()
            row.update(source="none", status="filled_zero_missing_hour",
                       n_invalid_domain=grid.ncols * grid.nrows, n_invalid_basin=n_basin)
        else:
            raw = reproject_clip(clip, grid, args.resampling)
            invalid = ~np.isfinite(raw)
            n_inv, n_inv_b = int(invalid.sum()), int((invalid & mask).sum())
            if n_inv and args.on_missing_cells == "error":
                die(f"{n_inv} cells have no valid MRMS data at label {lab} "
                    "(--on-missing-cells error)")
            data = np.where(invalid, 0.0, raw).astype("float32")
            row.update(source=src_kind + ("_" + args.fallback_product if src_kind == "fallback" else ""),
                       status="ok" if n_inv == 0 else "ok_cells_zero_filled",
                       n_invalid_domain=n_inv, n_invalid_basin=n_inv_b)
        write_ascii(dirs.out / fname, data, grid)
        bm = data[mask]
        row.update(domain_mean_mm=float(data.mean()), basin_mean_mm=float(bm.mean()),
                   basin_max_mm=float(bm.max()), domain_max_mm=float(data.max()))
        rows.append(row)
        if (i + 1) % 1000 == 0:
            print(f"  {i + 1}/{len(labels)} written")

    man = pd.DataFrame(rows)
    man.to_csv(dirs.work / "manifest.csv", index=False)
    write_domain_json(args, dirs, grid, labels)
    write_snippet(args, dirs, grid, labels)

    st = man["status"].value_counts().to_dict()
    print(f"  wrote {int((man['status'] != 'skipped_missing_hour').sum())} files; status counts: {st}")
    n_fill = int((man["status"] == "filled_zero_missing_hour").sum())
    if n_fill:
        print(f"  WARNING: {n_fill} hours had NO MRMS data and were written as all-zero grids "
              "(listed in manifest.csv).")
    n_cz = int((man["status"] == "ok_cells_zero_filled").sum())
    if n_cz:
        print(f"  NOTE: {n_cz} hours had some no-coverage cells set to 0 (n_invalid_* in manifest.csv).")
    print(f"  manifest: {dirs.work / 'manifest.csv'}")
    print(f"  .in snippet: {dirs.work / 'tribs_input_snippet.txt'}")
    return man


def write_domain_json(args, dirs: Dirs, grid: Grid, labels) -> None:
    info = {
        "script": "build_mrms_radar_rainfall_tribs6.py", "script_version": SCRIPT_VERSION,
        "created": datetime.now().isoformat(timespec="seconds"),
        "grid": asdict(grid),
        "clip_bbox_lonlat": clip_bbox_lonlat(grid),
        "run": {"start": f"{labels[0]:%Y-%m-%d %H:%M}", "end": f"{labels[-1]:%Y-%m-%d %H:%M}",
                "n_files": len(labels), "utc_offset": args.utc_offset,
                "stamp_convention": args.stamp_convention},
        "files": {"outdir": str(dirs.out), "basename": args.basename, "extension": args.extension},
        "source": {"product": args.product, "fallback_product": args.fallback_product,
                   "source_kind": args.source, "base_url": args.base_url,
                   "resampling": args.resampling,
                   "on_missing_cells": args.on_missing_cells, "on_missing_hours": args.on_missing_hours},
        "basin_file": str(Path(args.basin).resolve()) if args.basin else None,
        "bbox_lonlat": list(args.bbox) if args.bbox else None,
    }
    with open(dirs.work / "domain.json", "w") as f:
        json.dump(info, f, indent=2)


def write_snippet(args, dirs: Dirs, grid: Grid, labels) -> None:
    base_path = dirs.out / args.basename
    product_note = ""
    if args.source == "iem":
        product_note = f"""
Radar product note
  * This run used the Iowa State archive product {args.product}: radar rainfall
    adjusted with a local gauge bias correction. It is an OLDER algorithm than the
    AWS MultiSensor Pass 2 product used from Oct 2020 (MRMS v12). Do not splice
    the two products together inside one calibration without testing the difference.
"""
    text = f"""tRIBS 6.0.0 input-file (.in) settings for the gridded radar rainfall
generated by build_mrms_radar_rainfall_tribs6.py on {datetime.now():%Y-%m-%d %H:%M}

Copy these values into your .in file, matching the style of the neighbouring
lines (compare with big_spring.in from the tRIBS-benchmarks repo).

  RAINSOURCE     1
  RAINFILE       {base_path}
  RAINEXTENSION  {args.extension}
  RAININTRVL     1

Notes
  * RAINFILE = folder + base name only: NO date/hour and NO extension. tRIBS
    appends MMDDYYYYHH and the extension itself. The path above is absolute;
    convert it to the same relative style as your other paths (relative to
    where you run tRIBS).
  * RAINSOURCE: the tRIBS 6.0.0 Forcings page says 1 = gridded radar and
    2 = rain gauge. An older Model Input File page still says 1 = NEXRAD
    Stage III. VERIFY against big_spring.in before trusting this value.
  * RAININTRVL is in hours; gridded input is hourly only.
  * RAINSEARCH (hours to look ahead for a missing grid): use what your template
    uses. This script writes a file for every hour (missing hours are zero
    grids unless you chose --on-missing-hours skip), so it should not be needed.
  * Point met forcing (e.g. NLDAS) can stay as is; rainfall is configured
    independently. Confirm the NLDAS precip column is ignored when RAINSOURCE = 1.
  * Set your run start (STARTDATE) and RUNTIME so every hour is inside the
    files written:  {labels[0]:%m/%d/%Y %H:%M}  through  {labels[-1]:%m/%d/%Y %H:%M}  (local).
    If tRIBS complains about a missing final grid, rerun with --end one hour later.
  * UTCOFFSET in the .in file should be {args.utc_offset:+g} (local = UTC {args.utc_offset:+g} h).

Grid (must be in the same CRS as your mesh: {grid.crs})
  ncols {grid.ncols}   nrows {grid.nrows}   cellsize {grid.cellsize:g} m
  xllcorner {grid.xmin:.0f}   yllcorner {grid.ymin:.0f}

TIME CONVENTION USED: label = {args.stamp_convention.upper()} of the rain interval  (UNVERIFIED
in tRIBS - see the TIME CONVENTIONS section at the top of the script).
{product_note}"""
    (dirs.work / "tribs_input_snippet.txt").write_text(text)


# =============================================================================
# STAGE 3: AUDIT
# =============================================================================
def read_ascii(path: Path):
    """Return (header dict, 2-D data array). Raises ValueError on a malformed file."""
    with open(path) as f:
        hdr = {}
        for _ in range(6):
            parts = f.readline().split()
            if len(parts) != 2:
                raise ValueError("malformed header line")
            hdr[parts[0].lower()] = float(parts[1])
        body = f.read().split()
    need = {"ncols", "nrows", "xllcorner", "yllcorner", "cellsize", "nodata_value"}
    if not need.issubset(hdr):
        raise ValueError(f"header missing {sorted(need - set(hdr))}")
    n = int(hdr["ncols"]) * int(hdr["nrows"])
    if len(body) != n:
        raise ValueError(f"expected {n} values, found {len(body)}")
    data = np.array(body, dtype="float64").reshape(int(hdr["nrows"]), int(hdr["ncols"]))
    return hdr, data


class Report:
    """Collects PASS/WARN/FAIL/INFO lines and prints them together at the end."""
    def __init__(self):
        self.lines, self.n_fail, self.n_warn = [], 0, 0

    def add(self, level, msg):
        self.lines.append(f"[{level:<4}] {msg}")
        if level == "FAIL":
            self.n_fail += 1
        elif level == "WARN":
            self.n_warn += 1

    def text(self):
        return "\n".join(self.lines)


def stage_audit(args, dirs: Dirs, labels, grid_hint: Grid | None, geom) -> int:
    rep = Report()
    expected = {asc_filename(args.basename, args.extension, l): l for l in labels}
    man = None
    man_path = dirs.work / "manifest.csv"
    if man_path.exists():
        man = pd.read_csv(man_path)
    n_skip = 0 if man is None else int((man["status"] == "skipped_missing_hour").sum())

    print(f"AUDIT: checking {dirs.out}")
    rep.add("INFO", f"folder: {dirs.out}")
    rep.add("INFO", f"expected {len(expected)} files: {labels[0]:%Y-%m-%d %H:%M} to "
                    f"{labels[-1]:%Y-%m-%d %H:%M} local; stamp convention '{args.stamp_convention}'")

    # ---- 1. file inventory ---------------------------------------------------
    pat = re.compile(rf"^{re.escape(args.basename)}(\d{{10}})\.{re.escape(args.extension)}$")
    found = {}
    malformed = []
    for p in sorted(dirs.out.glob(f"{args.basename}*.{args.extension}")):
        m = pat.match(p.name)
        if not m:
            malformed.append(p.name)
            continue
        try:
            datetime.strptime(m.group(1), "%m%d%Y%H")
        except ValueError:
            malformed.append(p.name)
            continue
        found[p.name] = p
    if malformed:
        rep.add("FAIL", f"{len(malformed)} file names don't match <base>MMDDYYYYHH.{args.extension}: "
                        f"{malformed[:5]}")
    missing = [n for n in expected if n not in found]
    extra = [n for n in found if n not in expected]
    if missing:
        allowed = n_skip and len(missing) == n_skip
        rep.add("WARN" if allowed else "FAIL",
                f"{len(missing)} expected files are missing (first: {missing[:3]})"
                + (" - these were skipped on purpose (--on-missing-hours skip)" if allowed else ""))
    if extra:
        rep.add("WARN", f"{len(extra)} files outside the expected range are present "
                        f"(stale from an earlier run?) e.g. {extra[:3]}")
    if not missing and not malformed:
        rep.add("PASS", f"all {len(expected)} expected files are present")
    if not found:
        rep.add("FAIL", "no rainfall files found - nothing more to check")
        return _finish_audit(rep, dirs)

    # ---- 2. read every file: headers, values --------------------------------
    ref_hdr = None
    bad_hdr, bad_read, bad_nan, bad_neg, spikes_warn, spikes_fail = [], [], [], [], [], []
    depth_sum = None
    recs = []
    gauge_pix = None
    stats_geom_mask = None
    order = [n for n in expected if n in found]
    for name in order:
        try:
            hdr, data = read_ascii(found[name])
        except (ValueError, OSError) as exc:
            bad_read.append(f"{name}: {exc}")
            continue
        if ref_hdr is None:
            ref_hdr = hdr
            depth_sum = np.zeros_like(data)
            g = Grid(args.out_crs, hdr["cellsize"], hdr["xllcorner"], hdr["yllcorner"],
                     int(hdr["ncols"]), int(hdr["nrows"]))
            stats_geom_mask = basin_mask(geom, g)
            if args.gauge_csv:
                gauge_pix = _gauge_pixel(args, g)
        elif hdr != ref_hdr:
            bad_hdr.append(name)
            continue
        if not np.isfinite(data).all():
            bad_nan.append(name)
            continue
        if (data < 0).any():
            bad_neg.append(name)
            continue
        mx = float(data.max())
        if mx > args.rate_fail:
            spikes_fail.append((name, mx))
        elif mx > args.rate_warn:
            spikes_warn.append((name, mx))
        depth_sum += data
        rec = {"label_local": f"{expected[name]:%Y-%m-%d %H:%M}",
               "basin_mean_mm_hr": float(data[stats_geom_mask].mean()),
               "domain_max_mm_hr": mx}
        if gauge_pix is not None:
            rec["pixel_at_gauge_mm_hr"] = float(data[gauge_pix])
        recs.append(rec)

    if bad_read:
        rep.add("FAIL", f"{len(bad_read)} files could not be parsed, e.g. {bad_read[:3]}")
    if bad_hdr:
        rep.add("FAIL", f"{len(bad_hdr)} files have a header different from the first file, e.g. {bad_hdr[:3]}")
    if bad_nan:
        rep.add("FAIL", f"{len(bad_nan)} files contain NaN/inf, e.g. {bad_nan[:3]}")
    if bad_neg:
        rep.add("FAIL", f"{len(bad_neg)} files contain negative values (incl. -9999 in data), e.g. {bad_neg[:3]}")
    if spikes_fail:
        rep.add("FAIL", f"{len(spikes_fail)} files exceed {args.rate_fail:g} mm/hr, e.g. {spikes_fail[:3]}")
    if spikes_warn:
        rep.add("WARN", f"{len(spikes_warn)} files have a cell above {args.rate_warn:g} mm/hr "
                        f"(check they are real storms): {[(n, round(v, 1)) for n, v in spikes_warn[:5]]}")
    if not (bad_read or bad_hdr or bad_nan or bad_neg):
        rep.add("PASS", "all files parse, share one header, and contain no NaN/negative values")
    if not recs:
        rep.add("FAIL", "no usable files - nothing more to check")
        return _finish_audit(rep, dirs)

    # ---- 3. header vs plan, square cells, basin coverage -----------------------
    rep.add("INFO", f"grid: {int(ref_hdr['ncols'])} x {int(ref_hdr['nrows'])} cells, cellsize "
                    f"{ref_hdr['cellsize']:g} m, xll {ref_hdr['xllcorner']:.0f}, yll {ref_hdr['yllcorner']:.0f}")
    if grid_hint is not None:
        same = (int(ref_hdr["ncols"]) == grid_hint.ncols and int(ref_hdr["nrows"]) == grid_hint.nrows
                and abs(ref_hdr["xllcorner"] - grid_hint.xmin) < 1e-3
                and abs(ref_hdr["yllcorner"] - grid_hint.ymin) < 1e-3
                and abs(ref_hdr["cellsize"] - grid_hint.cellsize) < 1e-6)
        rep.add("PASS" if same else "FAIL",
                "file header matches the planned grid (domain.json / current settings)" if same
                else "file header does NOT match the planned grid - files are from a different run/settings")
    if geom is not None:
        gx0, gy0 = ref_hdr["xllcorner"], ref_hdr["yllcorner"]
        gx1 = gx0 + ref_hdr["ncols"] * ref_hdr["cellsize"]
        gy1 = gy0 + ref_hdr["nrows"] * ref_hdr["cellsize"]
        bx0, by0, bx1, by1 = geom.bounds
        inside = bx0 >= gx0 and by0 >= gy0 and bx1 <= gx1 and by1 <= gy1
        rep.add("PASS" if inside else "FAIL",
                "grid fully covers the basin polygon" if inside
                else "grid does NOT fully cover the basin polygon - enlarge --buffer-km and rerun")

    # ---- 4. manifest-derived flags -----------------------------------------
    if man is not None:
        n_fill = int((man["status"] == "filled_zero_missing_hour").sum())
        n_cz = int((man["status"] == "ok_cells_zero_filled").sum())
        n_fb = int(man["source"].astype(str).str.startswith("fallback").sum())
        if n_fill:
            rep.add("WARN", f"{n_fill} hours are zero-filled because MRMS had no data "
                            f"({100 * n_fill / len(man):.2f}% of hours); see manifest.csv")
        if n_cz:
            worst = int(man["n_invalid_basin"].max())
            rep.add("WARN", f"{n_cz} hours had some no-coverage cells set to 0 "
                            f"(max {worst} basin cells in one hour)")
        if n_fb:
            rep.add("WARN", f"{n_fb} hours came from the fallback product, not '{args.product}'")
        if not (n_fill or n_cz or n_fb):
            rep.add("PASS", "manifest: every hour has full primary-product data")
        conv = None
        dj = dirs.work / "domain.json"
        if dj.exists():
            conv = json.loads(dj.read_text()).get("run", {}).get("stamp_convention")
        if conv and conv != args.stamp_convention:
            rep.add("FAIL", f"files were built with stamp convention '{conv}' but audit was given "
                            f"'{args.stamp_convention}'")

    # ---- 5. hyetograph statistics ---------------------------------------------
    hy = pd.DataFrame(recs)
    hy["label_local"] = pd.to_datetime(hy["label_local"])
    total = float(hy["basin_mean_mm_hr"].sum())
    wet = int((hy["basin_mean_mm_hr"] > 0.1).sum())
    rep.add("INFO", f"basin-mean rainfall over the run: {total:.1f} mm  "
                    f"({wet} hours > 0.1 mm/hr, {100 * wet / len(hy):.1f}% of hours)")
    if total <= 0:
        rep.add("FAIL", "total rainfall is zero everywhere - wrong data or wrong time window?")
    top = hy[hy["basin_mean_mm_hr"] > 0].nlargest(5, "basin_mean_mm_hr")
    if len(top):
        rep.add("INFO", "wettest hours (basin mean mm/hr, local label): " + "; ".join(
            f"{r.label_local:%Y-%m-%d %H:%M} = {r.basin_mean_mm_hr:.2f}" for r in top.itertuples()))
    monthly = hy.set_index("label_local")["basin_mean_mm_hr"].resample("MS").sum()
    rep.add("INFO", "monthly basin-mean depth (mm): " + ", ".join(
        f"{i:%Y-%m}={v:.1f}" for i, v in monthly.items()))

    # ---- 6. optional gauge comparison ----------------------------------------
    if gauge_pix is not None:
        _gauge_compare(args, hy, rep)

    # ---- outputs --------------------------------------------------------------
    dirs.audit.mkdir(parents=True, exist_ok=True)
    hy.to_csv(dirs.audit / "basin_mean_hyetograph.csv", index=False)
    if not args.no_plots:
        _make_plots(args, dirs, hy, depth_sum, ref_hdr, geom)
    return _finish_audit(rep, dirs)


def _finish_audit(rep: Report, dirs: Dirs) -> int:
    verdict = "FAIL" if rep.n_fail else ("PASS with warnings" if rep.n_warn else "PASS")
    rep.lines.append(f"\nOVERALL: {verdict}  ({rep.n_fail} failures, {rep.n_warn} warnings)")
    txt = rep.text()
    dirs.audit.mkdir(parents=True, exist_ok=True)
    (dirs.audit / "audit_report.txt").write_text(txt + "\n")
    print(txt)
    print(f"\n  report: {dirs.audit / 'audit_report.txt'}")
    return 1 if rep.n_fail else 0


def _gauge_pixel(args, g: Grid):
    """(row, col) of the grid cell containing the gauge, or None if outside."""
    from pyproj import Transformer
    if args.gauge_lon is None or args.gauge_lat is None:
        die("--gauge-csv needs --gauge-lon and --gauge-lat (decimal degrees)")
    x, y = Transformer.from_crs("EPSG:4326", g.crs, always_xy=True).transform(
        args.gauge_lon, args.gauge_lat)
    col = math.floor((x - g.xmin) / g.cellsize)
    row = g.nrows - 1 - math.floor((y - g.ymin) / g.cellsize)
    if not (0 <= col < g.ncols and 0 <= row < g.nrows):
        die("the gauge location is outside the rainfall grid")
    return (row, col)


def _gauge_compare(args, hy: pd.DataFrame, rep: Report) -> None:
    """Compare the MRMS pixel over a gauge to the gauge's hourly totals."""
    df = pd.read_csv(args.gauge_csv)
    tcol = args.gauge_time_col or df.columns[0]
    vcol = args.gauge_value_col or df.columns[1]
    ts = pd.to_datetime(df[tcol], errors="coerce")
    vals = pd.to_numeric(df[vcol], errors="coerce") * args.gauge_scale
    ok = ts.notna() & vals.notna()
    s = pd.Series(vals[ok].values, index=pd.DatetimeIndex(ts[ok])).sort_index()
    if args.stamp_convention == "start":
        hourly = s.resample("1h", label="left", closed="left").sum(min_count=1)
    else:
        hourly = s.resample("1h", label="right", closed="right").sum(min_count=1)
    m = hy.set_index("label_local")["pixel_at_gauge_mm_hr"]
    joined = pd.concat([hourly.rename("gauge"), m.rename("mrms")], axis=1, join="inner").dropna()
    if len(joined) < 24:
        rep.add("WARN", f"gauge comparison skipped: only {len(joined)} overlapping hours "
                        "(check the gauge file's time column and period)")
        return
    g_tot, m_tot = joined["gauge"].sum(), joined["mrms"].sum()
    rep.add("INFO", f"gauge '{args.gauge_name}' vs MRMS pixel, {len(joined)} overlapping hours: "
                    f"gauge total {g_tot:.1f} mm, MRMS total {m_tot:.1f} mm, "
                    f"ratio MRMS/gauge = {m_tot / g_tot:.2f}" if g_tot > 0 else
                    f"gauge total is 0 mm over {len(joined)} hours")
    if g_tot <= 0:
        return
    lags = {}
    for k in range(-3, 4):
        # k > 0: compare gauge at time t with MRMS at time t + k hours
        shifted = m.shift(-k, freq="h")
        pair = pd.concat([hourly.rename("g"), shifted.rename("m")], axis=1, join="inner").dropna()
        lags[k] = pair["g"].corr(pair["m"]) if len(pair) > 24 else float("nan")
    best = max(lags, key=lambda k: -1 if np.isnan(lags[k]) else lags[k])
    rep.add("INFO", "hourly correlation gauge vs MRMS by lag (hours; + means the MRMS label is later): "
            + ", ".join(f"{k:+d}: {v:.2f}" for k, v in lags.items()))
    if best == 0:
        rep.add("PASS", f"gauge correlation peaks at lag 0 (r = {lags[0]:.2f}) - timing looks consistent")
    else:
        rep.add("WARN", f"gauge correlation peaks at lag {best:+d} h (r = {lags[best]:.2f}) rather than 0 "
                        f"(r = {lags[0]:.2f}). This may mean the stamp convention or the gauge's own "
                        "time convention is off by an hour. MRMS Pass 2 may also already include this "
                        "gauge, so this is not fully independent.")


def _make_plots(args, dirs: Dirs, hy, depth_sum, hdr, geom) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not installed - skipping plots)")
        return
    fig, ax = plt.subplots(figsize=(11, 3.5))
    # Draw each bar over the hour it represents: a 'start' label marks the
    # beginning of the interval (bar extends right); an 'end' label marks the
    # end (bar extends left).
    width = 1 / 24 if args.stamp_convention == "start" else -1 / 24
    ax.bar(hy["label_local"], hy["basin_mean_mm_hr"], width=width, align="edge", color="#378DBD")
    ax.set_ylabel("basin-mean rain (mm/hr)")
    ax.set_title("MRMS basin-mean hourly rainfall (tRIBS input, local time; "
                 f"file label = {args.stamp_convention.upper()} of interval)")
    fig.tight_layout()
    fig.savefig(dirs.audit / "basin_mean_hyetograph.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 5))
    x0, y0, cs = hdr["xllcorner"], hdr["yllcorner"], hdr["cellsize"]
    ext = (x0, x0 + hdr["ncols"] * cs, y0, y0 + hdr["nrows"] * cs)
    im = ax.imshow(depth_sum, extent=ext, cmap="YlGnBu", origin="upper")
    plt.colorbar(im, ax=ax, label="total depth over run (mm)")
    if geom is not None:
        for poly in getattr(geom, "geoms", [geom]):
            xs, ys = poly.exterior.xy
            ax.plot(xs, ys, color="#AB0520", lw=1.2)
    ax.set_title("Total MRMS rainfall depth in the tRIBS input grid")
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(axis="x", labelrotation=45)
    ax.set_xlabel("easting (m)")
    ax.set_ylabel("northing (m)")
    fig.tight_layout()
    fig.savefig(dirs.audit / "total_depth_map.png", dpi=130)
    plt.close(fig)


# =============================================================================
# check-env
# =============================================================================
def stage_check_env(args) -> int:
    ok = True
    print(f"python {sys.version.split()[0]}")
    for mod in ("numpy", "pandas", "requests", "rasterio", "pyproj", "geopandas", "matplotlib"):
        try:
            m = __import__(mod)
            print(f"  {mod:<11} {getattr(m, '__version__', '?')}")
        except ImportError:
            need = mod in ("numpy", "pandas", "requests", "rasterio", "pyproj")
            print(f"  {mod:<11} NOT INSTALLED" + ("  <-- required" if need else "  (optional)"))
            ok = ok and not need
    try:
        import rasterio
        with rasterio.Env() as env:
            drivers = env.drivers()
        print(f"  GDAL {rasterio.__gdal_version__}: GRIB driver "
              f"{'present' if 'GRIB' in drivers else 'MISSING  <-- required'}")
        ok = ok and ("GRIB" in drivers)
    except Exception as exc:  # noqa: BLE001
        print(f"  rasterio problem: {exc}")
        ok = False
    if not args.no_network:
        url = source_url(args.source, args.base_url, args.region, args.product,
                         SOURCES[args.source]["probe_stamp"])
        try:
            r = requests.head(url, timeout=args.timeout, allow_redirects=True)
            print(f"  network: HEAD {url}\n           -> HTTP {r.status_code}"
                  f"{'  (OK)' if r.status_code == 200 else '  <-- expected 200'}")
            ok = ok and r.status_code == 200
        except requests.RequestException as exc:
            print(f"  network: could not reach {args.base_url}: {exc}")
            ok = False
    print("\nENVIRONMENT OK" if ok else "\nENVIRONMENT HAS PROBLEMS (see above)")
    return 0 if ok else 1


# =============================================================================
# command line
# =============================================================================
def build_parser() -> argparse.ArgumentParser:
    U = USER_SETTINGS
    common = argparse.ArgumentParser(add_help=False)
    g = common.add_argument_group("where things live")
    g.add_argument("--workdir", default=U["workdir"])
    g.add_argument("--outdir", default=U["outdir"], help="folder for the .asc files (default <workdir>/rain)")
    g.add_argument("--basename", default=U["basename"], help="file prefix, e.g. radar_")
    g.add_argument("--extension", default=U["extension"], help="file extension without the period")
    g = common.add_argument_group("basin (give one)")
    g.add_argument("--basin", default=U["basin"], help="basin polygon file (.shp/.geojson/.gpkg)")
    g.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"), default=U["bbox"],
                   help="basin extent in lon/lat degrees")
    g = common.add_argument_group("target grid")
    g.add_argument("--out-crs", default=U["out_crs"], help="model CRS, e.g. EPSG:26912")
    g.add_argument("--cellsize", type=float, default=U["cellsize"], help="metres")
    g.add_argument("--buffer-km", type=float, default=U["buffer_km"])
    g.add_argument("--resampling", choices=["average", "bilinear", "nearest"], default=U["resampling"])
    g = common.add_argument_group("time (LOCAL model clock)")
    g.add_argument("--start", default=U["start"], help="first hour label, e.g. '2023-10-01 00:00'")
    g.add_argument("--end", default=U["end"], help="last hour label (inclusive)")
    g.add_argument("--utc-offset", type=float, default=U["utc_offset"], help="hours from UTC (Arizona -7)")
    g.add_argument("--stamp-convention", choices=["start", "end"], default=U["stamp_convention"],
                   help="does the hour in the file name mark the START or END of the rain interval? "
                        "(UNVERIFIED - see script header)")
    g = common.add_argument_group("MRMS source")
    g.add_argument("--source", choices=sorted(SOURCES), default=U["source"],
                   help="aws = NOAA bucket, MultiSensor Pass 2, starts 2020-10-14 (default); "
                        "iem = Iowa State archive, GaugeCorr product, for older years")
    g.add_argument("--product", default=U["product"],
                   help="default depends on --source (aws: MultiSensor_QPE_01H_Pass2_00.00, "
                        "iem: GaugeCorr_QPE_01H_00.00)")
    g.add_argument("--fallback-product", default=U["fallback_product"],
                   help="e.g. MultiSensor_QPE_01H_Pass1_00.00, used only for hours the main product lacks")
    g.add_argument("--base-url", "--s3-base-url", dest="base_url", default=U["base_url"],
                   help="default depends on --source (--s3-base-url is the old name)")
    g.add_argument("--region", default=U["region"], help="aws only")
    g = common.add_argument_group("behaviour")
    g.add_argument("--on-missing-cells", choices=["zero", "error"], default=U["on_missing_cells"])
    g.add_argument("--on-missing-hours", choices=["zero", "skip", "error"], default=U["on_missing_hours"])
    g.add_argument("--workers", type=int, default=U["workers"])
    g.add_argument("--retries", type=int, default=U["retries"])
    g.add_argument("--retry-wait", type=float, default=U["retry_wait"])
    g.add_argument("--timeout", type=float, default=U["timeout"])
    g.add_argument("--keep-raw", action="store_true", help="keep the full-CONUS .grib2.gz files (large!)")
    g.add_argument("--force", action="store_true", help="re-download even if a clip exists")
    g.add_argument("--dry-run", action="store_true", help="print the plan; download/write nothing")

    audit = argparse.ArgumentParser(add_help=False)
    g = audit.add_argument_group("audit options")
    g.add_argument("--rate-warn", type=float, default=100.0, help="mm/hr: warn above this")
    g.add_argument("--rate-fail", type=float, default=300.0, help="mm/hr: fail above this")
    g.add_argument("--no-plots", action="store_true")
    g.add_argument("--gauge-csv", help="hourly-or-finer gauge file (time in LOCAL clock, "
                                       "value = rain DEPTH per record, not a rate)")
    g.add_argument("--gauge-lon", type=float)
    g.add_argument("--gauge-lat", type=float)
    g.add_argument("--gauge-name", default="gauge")
    g.add_argument("--gauge-time-col", help="default: first column")
    g.add_argument("--gauge-value-col", help="default: second column")
    g.add_argument("--gauge-scale", type=float, default=1.0, help="multiply gauge values to get mm (25.4 for inches)")

    p = argparse.ArgumentParser(
        description="Build hourly MRMS radar rainfall grids for tRIBS 6.0.0 "
                    "(see the docstring at the top of this file).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="stage", required=True)
    sub.add_parser("download", parents=[common], help="download + clip MRMS hours")
    sub.add_parser("process", parents=[common], help="reproject and write tRIBS ASCII grids")
    sub.add_parser("audit", parents=[common, audit], help="check the written grids")
    sub.add_parser("all", parents=[common, audit], help="download, process, audit")
    ce = sub.add_parser("check-env", parents=[common], help="check packages, GRIB support, network")
    ce.add_argument("--no-network", action="store_true")
    return p


def run_period(args):
    """Resolve start/end (audit may fall back to domain.json) into hour labels."""
    if args.stage == "audit" and (args.start is None or args.end is None):
        dj = Path(args.workdir) / "domain.json"
        if not dj.exists():
            die("audit needs --start and --end (or a domain.json from an earlier process run)")
        run = json.loads(dj.read_text())["run"]
        args.start = args.start or run["start"]
        args.end = args.end or run["end"]
    if args.start is None or args.end is None:
        die("--start and --end are required (LOCAL clock, e.g. --start \"2023-10-01 00:00\" "
            "--end \"2024-09-30 23:00\")")
    return hour_labels(parse_local_hour(args.start, "--start"), parse_local_hour(args.end, "--end"))


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    resolve_source(args)
    if args.stage == "check-env":
        return stage_check_env(args)
    if args.cellsize < 900:
        print(f"NOTE: --cellsize {args.cellsize:g} m is finer than MRMS (~1 km); "
              "it adds no real detail.")
    labels = run_period(args)
    if args.source == "aws" and labels[0] < datetime(2020, 10, 14):
        print("NOTE: the AWS archive starts 2020-10-14; for earlier dates use --source iem.")
    if args.source == "iem" and labels[0] < datetime(2015, 7, 1):
        print("NOTE: the Iowa State GaugeCorr files were only confirmed from mid-2015 on; "
              "earlier hours may not exist.")
    dirs = make_dirs(args)
    grid, geom = build_grid(args.basin, args.bbox, args.out_crs, args.cellsize, args.buffer_km)
    print(f"build_mrms_radar_rainfall_tribs6.py v{SCRIPT_VERSION}   stage: {args.stage}")
    print(f"  period: {labels[0]:%Y-%m-%d %H:%M} to {labels[-1]:%Y-%m-%d %H:%M} local "
          f"({len(labels)} hours)   product: {args.product}   source: {args.source}")
    print_time_mapping(labels, args)
    if args.stage in ("download", "all"):
        stage_download(args, dirs, grid, labels)
    if args.stage in ("process", "all"):
        stage_process(args, dirs, grid, geom, labels)
    if args.stage in ("audit", "all") and not args.dry_run:
        return stage_audit(args, dirs, labels, grid, geom)
    return 0


if __name__ == "__main__":
    sys.exit(main())
