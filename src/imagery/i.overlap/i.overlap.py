#!/usr/bin/env python

##############################################################################
# MODULE:    i.overlap
#
# AUTHOR(S): Corey T. White <smortopahri@gmail.com>
#
# PURPOSE:   Reports image overlap from a directory of aerial photos from EXIF metadata.
#
# COPYRIGHT: (C) 2025 by Corey T. White and the GRASS Development Team
#
#            This program is free software under the GNU General Public
#            License (>=v2). Read the file COPYING that comes with GRASS
#            for details.
##############################################################################

"""Calculates overlap statistics from aerial photo EXIF data."""

# %module
# % description: Calculates overlap statistics from aerial photo EXIF data.
# % keyword: imagery
# % keyword: overlap
# % keyword: UAV
# %end
#
# %option G_OPT_F_INPUT
# % description: Directory of aerial photos
# % required: yes
# %end
#
# %option G_OPT_R_ELEV
# % key: elevation
# % required: yes
# %end
#
# %option
# % key: output
# % type: string
# % required: no
# % description: Output CSV file with overlap statistics
# %end
#
# %option
# % key: vector
# % type: string
# % required: no
# % description: Output vector map of image footprints
# %end

import sys
import os
import glob
import math
import csv
import numpy as np
from PIL import Image
from PIL.ExifTags import TAGS, GPSTAGS
from pyproj import CRS, Transformer
import grass.script as gs
from grass.pygrass.vector import VectorTopo
from grass.pygrass.vector.geometry import Area, Boundary, Centroid, Point

try:
    from pyproj import Geod

    geod = Geod(ellps="WGS84")
    USE_GEOD = True
except ImportError:
    USE_GEOD = False


def get_exif(image_path):
    """Extract relevant EXIF data efficiently with Pillow."""
    img = Image.open(image_path)
    exif = img._getexif()
    if not exif:
        return {}
    exif_data = {}
    for tag, val in exif.items():
        decoded = TAGS.get(tag, tag)
        if decoded == "GPSInfo":
            gps_data = {GPSTAGS.get(t, t): v for t, v in val.items()}
            exif_data["GPSInfo"] = gps_data
        else:
            exif_data[decoded] = val
    gs.debug(_("EXIF Data %s") % exif_data)
    return exif_data


def dms_to_dd(dms, ref):
    d, m, s = dms
    gs.debug(_("GPS data found in %s %s %s with ref %s") % (d, m, s, ref))
    dd = d + (m / 60.0 + (s / 3600.0))
    if ref in ["S", "W"]:
        dd = -dd
    return dd


def get_coords(exif):
    """Return lon, lat, alt if available."""
    gps = exif.get("GPSInfo")
    if not gps:
        return None
    lat = dms_to_dd(gps["GPSLatitude"], gps["GPSLatitudeRef"])
    lon = dms_to_dd(gps["GPSLongitude"], gps["GPSLongitudeRef"])
    alt = gps.get("GPSAltitude")
    if not alt:
        return None
    return lon, lat, alt


def get_focal_length(exif):
    """Return focal length in mm from EXIF data."""
    focal = exif.get("FocalLength")
    if not focal:
        return 0.1
    focal_mm = focal
    return focal_mm


def compute_sensor_size(exif):
    # image dimensions
    img_w = exif.get("ExifImageWidth")
    img_h = exif.get("ExifImageHeight")
    if not img_w or not img_h:
        return 0.1

    res_unit = exif.get("ResolutionUnit")  # 2=inches
    if "ResolutionUnit" not in exif:
        gs.warning("No ResolutionUnit found in EXIF data, assuming 1.0 (mm)")
        res_unit = 1.0

    conv = 25.4 if res_unit == 2 else 1.0
    # try to compute sensor size from FocalPlaneResolution
    if "FocalPlaneXResolution" in exif and "FocalPlaneYResolution" in exif:
        sensor_w_mm = (img_w / exif["FocalPlaneXResolution"]) * conv
        sensor_h_mm = (img_h / exif["FocalPlaneYResolution"]) * conv
    else:
        # fallback: common compact sensor (DJI, ~6.3mm or 13.2mm width)
        sensor_w_mm = 13.2
        sensor_h_mm = 8.8  # common height for compact sensors
    gs.debug(_("Sensor size: %smm x %smm") % (sensor_w_mm, sensor_h_mm))
    return (sensor_w_mm, sensor_h_mm)


def compute_gsd(exif, alt, focal_mm, sensor_size):
    """Estimate GSD (m/pixel) from EXIF and altitude."""
    # altitude in meters (if AGL, not AMSL!)
    if not alt:
        return 0.1  # fallback: 10 cm/px

    # image dimensions
    img_w = exif.get("ExifImageWidth")  # px
    img_h = exif.get("ExifImageHeight")  # px
    if not img_w or not img_h:
        return 0.1

    sensor_w_mm, sensor_h_mm = sensor_size
    gsd_w = (alt * sensor_w_mm) / (focal_mm * img_w)
    gsd_h = (alt * sensor_h_mm) / (focal_mm * img_h)
    gsd = (gsd_w + gsd_h) / 2  # average
    return gsd


def distance(p1, p2):
    """Meters between two lon/lat points."""
    if USE_GEOD:
        _, _, dist = geod.inv(p1[0], p1[1], p2[0], p2[1])
        return dist
    dx = (p2[0] - p1[0]) * 111320
    dy = (p2[1] - p1[1]) * 110540
    return math.sqrt(dx**2 + dy**2)


def calculate_overlaps(coords, gsd):
    """Forward overlaps based on consecutive image distances."""
    overlaps = []
    for i in range(1, len(coords)):
        dist = distance(coords[i - 1], coords[i])
        footprint = gsd * 1000  # assume ~1000 px footprint in track dir
        overlap = max(0, 1 - (dist / footprint))
        overlaps.append(overlap)
    return overlaps


def rotation_matrix(yaw, pitch, roll):
    """Return 3D rotation matrix from yaw/pitch/roll (degrees)."""
    yaw, pitch, roll = np.radians([yaw, pitch, roll])
    Ry = np.array(
        [[np.cos(yaw), -np.sin(yaw), 0], [np.sin(yaw), np.cos(yaw), 0], [0, 0, 1]]
    )
    Rp = np.array(
        [
            [np.cos(pitch), 0, np.sin(pitch)],
            [0, 1, 0],
            [-np.sin(pitch), 0, np.cos(pitch)],
        ]
    )
    Rr = np.array(
        [[1, 0, 0], [0, np.cos(roll), -np.sin(roll)], [0, np.sin(roll), np.cos(roll)]]
    )
    return Ry @ Rp @ Rr


def get_orientation(exif):
    """Extract yaw, pitch, roll; return defaults if not found."""
    yaw = exif.get("FlightYawDegree", 0.0)  # DJI tag
    pitch = exif.get("GimbalPitchDegree", -90.0)  # Default: nadir
    roll = exif.get("GimbalRollDegree", 0.0)

    # Some cameras use MakerNotes instead of EXIF
    # Ensure defaults are applied if None
    return yaw, pitch, roll


def intersect_ray_dem(x0, y0, z0, dir_vec, dem, step=1.0, max_dist=2000.0):
    """
    Trace a ray from camera through DEM until it hits ground.

    Args:
        x0, y0, z0 : camera position (in GRASS CRS)
        dir_vec    : ray direction vector (3D)
        dem        : DEM raster name
        step       : step size in meters along ray
        max_dist   : max distance to trace (m)

    Returns:
        (x, y) intersection coords in GRASS CRS, or None if no hit
    """
    # normalize direction
    dir_vec = dir_vec / np.linalg.norm(dir_vec)

    dist = 0.0
    while dist < max_dist:
        # advance along ray
        gx = x0 + dir_vec[0] * dist
        gy = y0 + dir_vec[1] * dist
        gz = z0 + dir_vec[2] * dist

        # sample DEM at this point
        gs.debug(_("Checking DEM at (%s, %s)") % (gx, gy))
        result = gs.raster_what(dem, coord=[[gx, gy]])
        val = result[0][dem]["value"]
        if val not in (None, "null", "No data"):
            ground_z = float(val)
            if gz <= ground_z:
                return gx, gy  # intersection point on ground

        dist += step

    gs.warning(f"No DEM intersection found for ray at ({x0},{y0})")
    return None


def make_footprint(lon, lat, alt, focal_length, sensor_size, yaw, pitch, roll, dem):
    """Return rectangle footprint coords around image center."""
    # footprint size in meters
    sensor_w, sensor_h = sensor_size
    # approximate offsets in degrees (small area assumption)
    # ground field of view (footprint) in meters
    # gfov_w = (sensor_w * alt) / focal_length
    # gfov_h = (sensor_h * alt) / focal_length

    corners = [
        (-sensor_w / 2, -sensor_h / 2),
        (sensor_w / 2, -sensor_h / 2),
        (sensor_w / 2, sensor_h / 2),
        (-sensor_w / 2, sensor_h / 2),
    ]

    # rotation
    R = rotation_matrix(yaw, pitch, roll)
    x0, y0 = lon, lat
    z0 = alt  # camera height above reference plane
    footprint = []
    for cx, cy in corners:
        # direction vector
        dir_vec = R @ np.array([cx, cy, focal_length])
        hit = intersect_ray_dem(x0, y0, z0, dir_vec, dem)
        if hit:
            footprint.append(hit)

    if not footprint:
        gs.warning("No DEM intersections found, footprint empty")
        return []

    footprint.append(footprint[0])  # close polygon

    return footprint


def write_vector(footprints, outmap):
    """Export footprint polygons to a GRASS vector map with pygrass."""
    COLS_TYPES = {"cat": "INTEGER PRIMARY KEY"}
    with VectorTopo(
        outmap, mode="w", with_z=False, tab_cols=COLS_TYPES, layer=1, overwrite=True
    ) as vect:
        for i, poly in enumerate(footprints):
            # Create boundary from points
            boundary = Boundary(points=[Point(x, y) for x, y in poly])
            # Centroid with category ID
            centroid = Centroid(x=poly[0][0], y=poly[0][1])
            # Add area
            vect.write(centroid)
            vect.write(geo_obj=boundary, cat=i + 1, attrs=())
        vect.table.conn.commit()
        vect.build()


def create_transformer():
    """Reproject list of (lon,lat) coords from WGS84 to GRASS CRS."""
    grass_proj = gs.read_command("g.proj", flags="jf")  # PROJ JSON string
    grass_crs = CRS.from_string(grass_proj.strip())
    wgs84 = CRS.from_epsg(4326)

    # Build transformer (lon/lat WGS84 → GRASS CRS)
    transformer = Transformer.from_crs(wgs84, grass_crs, always_xy=True)
    return transformer


def get_height_above_ground(lon, lat, alt, elevation):
    """Calculate height above ground level."""
    result = gs.raster_what(elevation, coord=[[lon, lat]])
    ground_elev = None
    if elevation in result[0]:
        val = result[0][elevation]["value"]
        if val not in (None, "null", "No data"):
            ground_elev = float(val)

    if ground_elev is None:
        gs.warning(
            _("Elevation raster %s not found at coordinates (%s, %s)")
            % (elevation, lon, lat)
        )
        return alt

    return alt - ground_elev


def main():
    options, flags = gs.parser()
    indir = options["input"]
    elevation = options["elevation"]
    outcsv = options["output"]
    outvec = options["vector"]

    photos = sorted(glob.glob(os.path.join(indir, "*.jpg")))
    coords, footprints, rows = [], [], []
    transformer = create_transformer()
    for img in photos:
        exif = get_exif(img)
        gps = get_coords(exif)
        if not gps:
            continue
        lon, lat, alt = gps
        focal_length_mm = get_focal_length(exif)
        sensor_size = compute_sensor_size(exif)
        gsd = compute_gsd(exif, alt, focal_length_mm, sensor_size)
        rlon, rlat = transformer.transform(lon, lat)  # reproject lon/lat
        coords.append((rlon, rlat))
        height_above_ground = get_height_above_ground(rlon, rlat, alt, elevation)
        yaw, pitch, roll = get_orientation(exif)
        footprint = make_footprint(
            lon=rlon,
            lat=rlat,
            alt=alt,
            focal_length=focal_length_mm,
            sensor_size=sensor_size,
            yaw=yaw,
            pitch=pitch,
            roll=roll,
            dem=elevation,
        )
        footprints.append(footprint)

    if not coords:
        gs.fatal("No GPS data found")

    overlaps = calculate_overlaps(coords, gsd)
    avg_overlap = sum(overlaps) / len(overlaps)
    gs.message(f"Average overlap: {avg_overlap:.2f}")

    if outcsv:
        with open(outcsv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["photo1", "photo2", "overlap"])
            for i, ov in enumerate(overlaps):
                writer.writerow([photos[i], photos[i + 1], f"{ov:.2f}"])

    if outvec:
        write_vector(footprints, outvec)


if __name__ == "__main__":
    sys.exit(main())
