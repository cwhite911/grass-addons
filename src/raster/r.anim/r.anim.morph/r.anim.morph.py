#!/usr/bin/env python3
"""
MODULE:    r.anim.morph
AUTHOR(S): Implemented from Baia framework by Lobo, Appert & Pietriga (2018)
           "Animation Plans for Before-and-After Satellite Images"
           IEEE Transactions on Visualization and Computer Graphics
PURPOSE:   Create pixel-based animated transitions between before-and-after
           raster maps using animation plans. Supports CONTRACTION, EXPANSION,
           DEFORMATION, APPEARANCE, DISAPPEARANCE, DIRECTIONAL, RADIAL,
           DEM-based, and monolithic BLEND primitives.
COPYRIGHT: (C) 2024 by the GRASS Development Team
           This program is free software under the GNU General Public License
           (>=v2). Read the file COPYING that comes with GRASS for details.
"""

# --------------------------------------------------------------------------
# GRASS module parameter block
# --------------------------------------------------------------------------

# %module
# % description: Morphs a before raster into an after raster as a spatially-aware animated transition using Baia animation plans
# % keyword: raster
# % keyword: animation
# % keyword: animation
# % keyword: morph
# % keyword: temporal
# % keyword: imagery
# % keyword: change detection
# % keyword: before-and-after
# % keyword: satellite
# %end

# --- Input/output -----------------------------------------------------------

# %option G_OPT_R_INPUT
# % key: before
# % label: Before raster map(s)
# % description: Single raster or comma-separated list of band rasters (R,G,B) for the before image
# % required: yes
# %end

# %option G_OPT_R_INPUT
# % key: after
# % label: After raster map(s)
# % description: Single raster or comma-separated list of band rasters (R,G,B) for the after image
# % required: yes
# %end

# %option G_OPT_R_OUTPUT
# % key: output
# % label: Output prefix for animation frames
# % description: Frames will be named <output>_001, <output>_002, etc. For multi-band, <output>_<band>_<frame>
# % required: yes
# %end

# --- Animation parameters ---------------------------------------------------

# %option
# % key: frames
# % type: integer
# % label: Number of animation frames
# % description: Total number of output raster frames to generate
# % answer: 30
# % options: 2-1000
# % required: yes
# %end

# %option
# % key: primitive
# % type: string
# % label: Animation primitive
# % description: Type of spatial transition to apply to the region of interest
# % options: blend,appearance,disappearance,contraction,expansion,deformation,radial,directional,dem,plan
# % descriptions: blend;Monolithic blend (uniform across all pixels);appearance;Region appears from nothing;disappearance;Region disappears;contraction;Shape shrinks inward (e.g. shrinking lake);expansion;Shape grows outward (e.g. flooding);deformation;Shape changes form (combines contraction and expansion);radial;Radial gradient from center point outward;directional;Linear progression in a compass direction;dem;Transition ordered by elevation (e.g. snow accumulation);plan;Use a pre-computed animation plan raster
# % answer: blend
# % required: yes
# %end

# --- Mask options -----------------------------------------------------------

# %option G_OPT_R_INPUT
# % key: mask_before
# % label: ROI mask in before image
# % description: Binary raster mask for the region of interest in the before image (1=ROI, 0=background). Required for contraction, expansion, deformation, appearance, disappearance.
# % required: no
# %end

# %option G_OPT_R_INPUT
# % key: mask_after
# % label: ROI mask in after image
# % description: Binary raster mask for the region of interest in the after image (1=ROI, 0=background). Required for contraction, expansion, deformation.
# % required: no
# %end

# --- Primitive-specific options ---------------------------------------------

# %option G_OPT_R_INPUT
# % key: dem
# % label: Digital elevation model raster
# % description: DEM used to derive transition timing (required for dem primitive). Higher elevations transition first by default.
# % required: no
# %end

# %option
# % key: direction
# % type: string
# % label: Direction for directional primitive
# % description: Compass direction of animation progression (N=animation moves northward, etc.)
# % options: N,S,E,W,NE,NW,SE,SW
# % answer: N
# % required: no
# %end

# %option G_OPT_R_INPUT
# % key: plan
# % label: Pre-computed animation plan (S matrix)
# % description: Raster holding normalized start times S[i,j] in [0,1]. A second raster named <plan>_E must also exist for end times.
# % required: no
# %end

# --- Staging / timing -------------------------------------------------------

# %option
# % key: roi_start
# % type: double
# % label: ROI animation start time
# % description: Normalized time [0,1] at which the ROI begins transitioning (enables staging)
# % answer: 0.0
# % required: no
# %end

# %option
# % key: roi_end
# % type: double
# % label: ROI animation end time
# % description: Normalized time [0,1] at which the ROI finishes transitioning
# % answer: 1.0
# % required: no
# %end

# %option
# % key: bg_start
# % type: double
# % label: Background animation start time
# % description: Normalized time [0,1] at which background begins transitioning
# % answer: 0.0
# % required: no
# %end

# %option
# % key: bg_end
# % type: double
# % label: Background animation end time
# % description: Normalized time [0,1] at which background finishes transitioning
# % answer: 1.0
# % required: no
# %end

# --- Additional options -----------------------------------------------------

# %option
# % key: stages
# % type: string
# % label: Multi-stage JSON definition file
# % description: Path to a JSON file defining a sequence of animation stages (see documentation for format)
# % required: no
# %end

# %option
# % key: output_plan
# % type: string
# % label: Output animation plan prefix
# % description: If set, saves the computed S and E matrices as rasters named <output_plan>_S and <output_plan>_E
# % required: no
# %end

# %option
# % key: radial_x
# % type: double
# % label: Radial center X coordinate
# % description: X (easting) coordinate for the center of radial animation. Defaults to map center.
# % required: no
# %end

# %option
# % key: radial_y
# % type: double
# % label: Radial center Y coordinate
# % description: Y (northing) coordinate for the center of radial animation. Defaults to map center.
# % required: no
# %end

# %option
# % key: blend_duration
# % type: double
# % label: Blend duration per pixel
# % description: Duration [0,1] of the blend window for each pixel (E[i,j] - S[i,j]). Use 0 for instantaneous swap.
# % answer: 0.3
# % required: no
# %end

# %flag
# % key: i
# % description: Invert DEM or distance gradient (e.g. snow accumulates from low to high elevation)
# %end

# %flag
# % key: c
# % description: Apply color-transfer preprocessing to match before image histogram to after image
# %end

# %flag
# % key: p
# % description: Print animation plan statistics and exit without generating frames
# %end

# --------------------------------------------------------------------------
# Implementation
# --------------------------------------------------------------------------

import os
import sys
import json

import grass.script as gs
from grass.script import array as garray


def read_raster_array(name):
    """Read a GRASS raster into a numpy float64 array."""
    arr = garray.array(name)
    return arr.astype("float64")


def write_raster_array(arr, name, overwrite=False):
    """Write a numpy array back to a GRASS raster map."""
    out = garray.array()
    out[...] = arr
    out.write(name, overwrite=overwrite)
    gs.run_command(
        "r.support", map=name, quiet=True, description="r.anim.morph animation frame"
    )


# --------------------------------------------------------------------------
# Animation plan primitives
# --------------------------------------------------------------------------


def plan_blend(shape, roi_window, bg_window):
    """
    Monolithic blend: every pixel transitions uniformly.
    roi_window / bg_window are (start, end) tuples in [0,1] but since
    there is no mask distinction here both map to the full image.
    """
    S = np.full(shape, float(roi_window[0]))
    E = np.full(shape, float(roi_window[1]))
    return S, E


def plan_appearance(mask_after, roi_window, bg_window, blend_duration):
    """
    APPEARANCE: A new object appears in the after image.
    Inside mask: progressive reveal from edge outward.
    Outside mask: background transitions during bg_window.
    """
    import numpy as np
    from scipy import ndimage

    shape = mask_after.shape
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))

    roi = mask_after > 0
    if not roi.any():
        gs.warning("mask_after is empty; returning uniform blend plan.")
        return S, E

    # Distance from edge of appearing region inward
    dist = ndimage.distance_transform_edt(roi)
    max_d = dist.max()
    if max_d > 0:
        norm_dist = dist / max_d  # 0=edge, 1=center
    else:
        norm_dist = np.zeros(shape)

    roi_dur = roi_window[1] - roi_window[0]
    # Edge pixels transition first, center last
    S[roi] = roi_window[0] + norm_dist[roi] * (roi_dur - blend_duration)
    E[roi] = S[roi] + blend_duration
    E[roi] = np.clip(E[roi], 0, 1)
    return S, E


def plan_disappearance(mask_before, roi_window, bg_window, blend_duration):
    """
    DISAPPEARANCE: An object that exists in before vanishes.
    Inside mask: progressive fade from edge inward.
    Outside mask: background transitions during bg_window.
    """
    import numpy as np
    from scipy import ndimage

    shape = mask_before.shape
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))

    roi = mask_before > 0
    if not roi.any():
        gs.warning("mask_before is empty; returning uniform blend plan.")
        return S, E

    # Distance from outside (background) into the disappearing region
    dist = ndimage.distance_transform_edt(roi)
    max_d = dist.max()
    if max_d > 0:
        norm_dist = dist / max_d
    else:
        norm_dist = np.zeros(shape)

    roi_dur = roi_window[1] - roi_window[0]
    # Outer ring disappears first, center last
    S[roi] = roi_window[0] + (1.0 - norm_dist[roi]) * (roi_dur - blend_duration)
    E[roi] = S[roi] + blend_duration
    E[roi] = np.clip(E[roi], 0, 1)
    return S, E


def plan_contraction(mask_before, mask_after, roi_window, bg_window, blend_duration):
    """
    CONTRACTION: Shape shrinks from before mask to after mask.
    Pixels in (before AND NOT after) transition ordered by distance to
    the new contour (after mask edge), far pixels first.
    The intersection transitions uniformly within the roi_window.
    Background transitions in bg_window.
    """
    import numpy as np
    from scipy import ndimage

    shape = mask_before.shape
    mb = mask_before > 0
    ma = mask_after > 0 if mask_after is not None else np.zeros(shape, bool)

    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))

    intersection = mb & ma
    disappearing = mb & ~ma

    if not disappearing.any():
        gs.warning(
            "No pixels to contract (masks may be identical or mask_after larger). "
            "Returning blend plan."
        )
        S[mb] = roi_window[0]
        E[mb] = roi_window[1]
        return S, E

    # Distance from each disappearing pixel to the after-mask boundary
    # We want distance measured *inside* the disappearing region to the
    # nearest after-mask pixel.
    if ma.any():
        # Pixels inside disappearing region: distance to nearest after-mask pixel
        dist_to_after = ndimage.distance_transform_edt(~ma)
    else:
        # No after mask -> full contraction to nothing
        dist_to_after = ndimage.distance_transform_edt(~np.zeros_like(mb))
        dist_to_after = ndimage.distance_transform_edt(mb)

    d = dist_to_after[disappearing]
    dmax = d.max()
    if dmax > 0:
        norm = (
            d / dmax
        )  # 0=close to new contour (transitions last), 1=far (transitions first)
    else:
        norm = np.zeros_like(d, dtype=float)

    roi_dur = roi_window[1] - roi_window[0]
    # Far pixels (norm≈1) transition early; close pixels (norm≈0) transition late
    S_vals = roi_window[0] + (1.0 - norm) * (roi_dur - blend_duration)
    E_vals = S_vals + blend_duration
    S[disappearing] = np.clip(S_vals, 0, 1)
    E[disappearing] = np.clip(E_vals, 0, 1)

    # Intersection: uniform transition across full roi_window
    S[intersection] = roi_window[0]
    E[intersection] = roi_window[1]

    return S, E


def plan_expansion(mask_before, mask_after, roi_window, bg_window, blend_duration):
    """
    EXPANSION: Shape grows from before mask to after mask.
    New pixels (in after AND NOT before) transition ordered by distance
    from the old contour, close pixels first, far pixels last.
    """
    import numpy as np
    from scipy import ndimage

    shape = mask_before.shape
    mb = mask_before > 0 if mask_before is not None else np.zeros(shape, bool)
    ma = mask_after > 0

    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))

    intersection = mb & ma
    appearing = ma & ~mb

    if not appearing.any():
        gs.warning(
            "No pixels to expand (masks may be identical or mask_before larger). "
            "Returning blend plan."
        )
        S[ma] = roi_window[0]
        E[ma] = roi_window[1]
        return S, E

    # Distance of each appearing pixel to the nearest before-mask pixel
    if mb.any():
        dist_to_before = ndimage.distance_transform_edt(~mb)
    else:
        dist_to_before = ndimage.distance_transform_edt(~np.zeros_like(ma))

    d = dist_to_before[appearing]
    dmax = d.max()
    if dmax > 0:
        norm = d / dmax  # 0=close to old contour (first), 1=far (last)
    else:
        norm = np.zeros_like(d, dtype=float)

    roi_dur = roi_window[1] - roi_window[0]
    S_vals = roi_window[0] + norm * (roi_dur - blend_duration)
    E_vals = S_vals + blend_duration
    S[appearing] = np.clip(S_vals, 0, 1)
    E[appearing] = np.clip(E_vals, 0, 1)

    # Intersection: uniform
    S[intersection] = roi_window[0]
    E[intersection] = roi_window[1]

    return S, E


def plan_deformation(mask_before, mask_after, roi_window, bg_window, blend_duration):
    """
    DEFORMATION: Shape changes form.
    Combines CONTRACTION (pixels leaving) and EXPANSION (pixels entering).
    Intersection transitions uniformly. Works by superimposing both plans.
    """
    import numpy as np
    from scipy import ndimage

    shape = mask_before.shape
    mb = mask_before > 0
    ma = mask_after > 0

    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))

    intersection = mb & ma
    contracting = mb & ~ma
    expanding = ma & ~mb

    # Intersection: uniform
    S[intersection] = roi_window[0]
    E[intersection] = roi_window[1]

    # Contracting pixels (disappearing from before)
    if contracting.any():
        if ma.any():
            dist_to_after = ndimage.distance_transform_edt(~ma)
        else:
            dist_to_after = ndimage.distance_transform_edt(mb)
        d = dist_to_after[contracting]
        dmax = d.max()
        norm = d / dmax if dmax > 0 else np.zeros_like(d, float)
        roi_dur = roi_window[1] - roi_window[0]
        S[contracting] = np.clip(
            roi_window[0] + (1.0 - norm) * (roi_dur - blend_duration), 0, 1
        )
        E[contracting] = np.clip(S[contracting] + blend_duration, 0, 1)

    # Expanding pixels (new in after)
    if expanding.any():
        if mb.any():
            dist_to_before = ndimage.distance_transform_edt(~mb)
        else:
            dist_to_before = ndimage.distance_transform_edt(ma)
        d = dist_to_before[expanding]
        dmax = d.max()
        norm = d / dmax if dmax > 0 else np.zeros_like(d, float)
        roi_dur = roi_window[1] - roi_window[0]
        S[expanding] = np.clip(roi_window[0] + norm * (roi_dur - blend_duration), 0, 1)
        E[expanding] = np.clip(S[expanding] + blend_duration, 0, 1)

    return S, E


def plan_radial(
    shape, cx, cy, roi_window, bg_window, blend_duration, mask=None, invert=False
):
    """
    RADIAL: Animation radiates outward from a center point.
    Pixels close to center transition first; distant pixels last.
    """
    import numpy as np

    rows, cols = shape
    Y, X = np.mgrid[0:rows, 0:cols]
    dist = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    max_d = dist.max()
    norm = dist / max_d if max_d > 0 else np.zeros(shape)

    if invert:
        norm = 1.0 - norm

    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))

    if mask is not None:
        roi = mask > 0
    else:
        roi = np.ones(shape, bool)

    roi_dur = roi_window[1] - roi_window[0]
    S[roi] = roi_window[0] + norm[roi] * (roi_dur - blend_duration)
    E[roi] = np.clip(S[roi] + blend_duration, 0, 1)

    return S, E


def plan_directional(
    shape, direction, roi_window, bg_window, blend_duration, mask=None, invert=False
):
    """
    DIRECTIONAL: Linear progression along a compass bearing.
    E.g., direction='N' means animation sweeps from south to north.
    """
    import numpy as np

    rows, cols = shape
    Y, X = np.mgrid[0:rows, 0:cols]  # Y: row 0=north, row N=south in raster convention

    dir_map = {
        "N": (-Y, None),
        "S": (Y, None),
        "E": (X, None),
        "W": (-X, None),
        "NE": (-Y + X, None),
        "NW": (-Y - X, None),
        "SE": (Y + X, None),
        "SW": (Y - X, None),
    }
    raw, _ = dir_map[direction.upper()]
    mn, mx = raw.min(), raw.max()
    norm = (raw - mn) / (mx - mn) if mx > mn else np.zeros(shape)

    if invert:
        norm = 1.0 - norm

    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))

    if mask is not None:
        roi = mask > 0
    else:
        roi = np.ones(shape, bool)

    roi_dur = roi_window[1] - roi_window[0]
    S[roi] = roi_window[0] + norm[roi] * (roi_dur - blend_duration)
    E[roi] = np.clip(S[roi] + blend_duration, 0, 1)

    return S, E


def plan_dem(dem_arr, roi_window, bg_window, blend_duration, mask=None, invert=False):
    """
    DEM-based: Transition timing derived from terrain elevation.
    Default: high elevations transition first (e.g. snow falls from peaks).
    invert=True: low elevations first (e.g. flood waters rise).
    """
    import numpy as np

    shape = dem_arr.shape
    d = dem_arr.copy()
    valid = np.isfinite(d)
    mn = d[valid].min()
    mx = d[valid].max()

    norm = np.zeros(shape)
    if mx > mn:
        norm[valid] = (d[valid] - mn) / (mx - mn)

    # High elevation -> norm=1. By default high goes first -> invert for S
    if not invert:
        norm = 1.0 - norm  # high elev gets low S (starts early)

    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))

    if mask is not None:
        roi = mask > 0
    else:
        roi = np.ones(shape, bool)

    roi_dur = roi_window[1] - roi_window[0]
    S[roi] = roi_window[0] + norm[roi] * (roi_dur - blend_duration)
    E[roi] = np.clip(S[roi] + blend_duration, 0, 1)

    return S, E


# --------------------------------------------------------------------------
# Frame rendering
# --------------------------------------------------------------------------


def render_frame(before_arr, after_arr, S, E, t):
    """
    Render a single animation frame at normalized time t ∈ [0,1].

    α(i,j) = clamp((t - S[i,j]) / max(E[i,j] - S[i,j], ε), 0, 1)
    frame   = (1-α) * before + α * after
    """
    import numpy as np

    duration = np.maximum(E - S, 1e-10)
    alpha = np.clip((t - S) / duration, 0.0, 1.0)
    return (1.0 - alpha) * before_arr + alpha * after_arr


# --------------------------------------------------------------------------
# Color transfer (histogram matching, optional)
# --------------------------------------------------------------------------


def color_transfer(source, target):
    """
    Match source raster statistics to target raster (mean/std transfer).
    Helps when before/after images have different illumination conditions.
    """
    import numpy as np

    valid_s = source[np.isfinite(source)]
    valid_t = target[np.isfinite(target)]

    if valid_s.size == 0 or valid_t.size == 0:
        return source

    mu_s, std_s = valid_s.mean(), valid_s.std()
    mu_t, std_t = valid_t.mean(), valid_t.std()

    if std_s < 1e-10:
        return source

    result = (source - mu_s) * (std_t / std_s) + mu_t
    return result


# --------------------------------------------------------------------------
# Multi-stage support
# --------------------------------------------------------------------------


def compose_plans(stages_list, shape):
    """
    Compose multiple (S, E) animation plan stages into a single (S, E) pair
    by taking the element-wise union (minimum S, maximum E for each pixel).

    stages_list: list of (S_arr, E_arr) tuples
    """
    import numpy as np

    S_final = np.ones(shape)
    E_final = np.zeros(shape)

    for S, E in stages_list:
        S_final = np.minimum(S_final, S)
        E_final = np.maximum(E_final, E)

    return S_final, E_final


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main():
    import numpy as np

    options, flags = gs.parser()

    # ---- Parse inputs -------------------------------------------------------
    before_names = [b.strip() for b in options["before"].split(",")]
    after_names = [a.strip() for a in options["after"].split(",")]

    if len(before_names) != len(after_names):
        gs.fatal("Number of 'before' and 'after' maps must match.")

    n_bands = len(before_names)
    output = options["output"]
    n_frames = int(options["frames"])
    primitive = options["primitive"].lower()
    overwrite = gs.overwrite()

    roi_window = (float(options["roi_start"]), float(options["roi_end"]))
    bg_window = (float(options["bg_start"]), float(options["bg_end"]))
    blend_dur = float(options["blend_duration"])
    invert = flags["i"]
    do_color = flags["c"]
    plan_only = flags["p"]

    # Validate timing
    for win, label in [(roi_window, "roi"), (bg_window, "bg")]:
        if not (0.0 <= win[0] <= win[1] <= 1.0):
            gs.fatal(f"{label}_start and {label}_end must satisfy 0 ≤ start ≤ end ≤ 1")

    # ---- Read first band to get raster shape --------------------------------
    gs.message("Reading raster metadata ...")
    before_arrays = [read_raster_array(n) for n in before_names]
    after_arrays = [read_raster_array(n) for n in after_names]
    shape = before_arrays[0].shape

    gs.message(
        f"Raster shape: {shape[0]} rows × {shape[1]} cols | "
        f"{n_bands} band(s) | {n_frames} frames"
    )

    # ---- Optional: color transfer -------------------------------------------
    if do_color:
        gs.message("Applying color transfer (histogram matching) ...")
        before_arrays = [
            color_transfer(b, a) for b, a in zip(before_arrays, after_arrays)
        ]

    # ---- Read optional masks -------------------------------------------------
    mask_before_arr = None
    mask_after_arr = None
    if options["mask_before"]:
        mask_before_arr = read_raster_array(options["mask_before"])
    if options["mask_after"]:
        mask_after_arr = read_raster_array(options["mask_after"])

    # ---- Handle multi-stage JSON file ---------------------------------------
    if options["stages"]:
        gs.message("Loading multi-stage animation plan from JSON ...")
        with open(options["stages"]) as f:
            stage_defs = json.load(f)
        stages_plans = []
        for sd in stage_defs:
            S_s, E_s = build_plan_from_def(
                sd, shape, before_arrays, after_arrays, invert, blend_dur
            )
            stages_plans.append((S_s, E_s))
        S, E = compose_plans(stages_plans, shape)

    else:
        # ---- Build animation plan from selected primitive -------------------
        gs.message(f"Computing animation plan: primitive='{primitive}' ...")

        if primitive == "blend":
            S, E = plan_blend(shape, roi_window, bg_window)

        elif primitive == "appearance":
            if mask_after_arr is None:
                gs.fatal("Primitive 'appearance' requires mask_after.")
            S, E = plan_appearance(mask_after_arr, roi_window, bg_window, blend_dur)

        elif primitive == "disappearance":
            if mask_before_arr is None:
                gs.fatal("Primitive 'disappearance' requires mask_before.")
            S, E = plan_disappearance(mask_before_arr, roi_window, bg_window, blend_dur)

        elif primitive == "contraction":
            if mask_before_arr is None:
                gs.fatal("Primitive 'contraction' requires mask_before.")
            S, E = plan_contraction(
                mask_before_arr, mask_after_arr, roi_window, bg_window, blend_dur
            )

        elif primitive == "expansion":
            if mask_after_arr is None:
                gs.fatal("Primitive 'expansion' requires mask_after.")
            S, E = plan_expansion(
                mask_before_arr, mask_after_arr, roi_window, bg_window, blend_dur
            )

        elif primitive == "deformation":
            if mask_before_arr is None or mask_after_arr is None:
                gs.fatal(
                    "Primitive 'deformation' requires both mask_before and mask_after."
                )
            S, E = plan_deformation(
                mask_before_arr, mask_after_arr, roi_window, bg_window, blend_dur
            )

        elif primitive == "radial":
            region = gs.region()
            # Convert geographic center to pixel coordinates
            if options["radial_x"] and options["radial_y"]:
                rx = float(options["radial_x"])
                ry = float(options["radial_y"])
                cx = (rx - region["w"]) / region["ewres"]
                cy = (region["n"] - ry) / region["nsres"]
            else:
                cx = shape[1] / 2.0
                cy = shape[0] / 2.0
            mask_arr = (
                mask_before_arr if mask_before_arr is not None else mask_after_arr
            )
            S, E = plan_radial(
                shape,
                cx,
                cy,
                roi_window,
                bg_window,
                blend_dur,
                mask=mask_arr,
                invert=invert,
            )

        elif primitive == "directional":
            mask_arr = (
                mask_before_arr if mask_before_arr is not None else mask_after_arr
            )
            S, E = plan_directional(
                shape,
                options["direction"],
                roi_window,
                bg_window,
                blend_dur,
                mask=mask_arr,
                invert=invert,
            )

        elif primitive == "dem":
            if not options["dem"]:
                gs.fatal("Primitive 'dem' requires a dem raster.")
            dem_arr = read_raster_array(options["dem"])
            mask_arr = (
                mask_before_arr if mask_before_arr is not None else mask_after_arr
            )
            S, E = plan_dem(
                dem_arr, roi_window, bg_window, blend_dur, mask=mask_arr, invert=invert
            )

        elif primitive == "plan":
            if not options["plan"]:
                gs.fatal(
                    "Primitive 'plan' requires the plan option (raster name for S matrix)."
                )
            plan_name = options["plan"]
            S = read_raster_array(plan_name)
            E_name = plan_name + "_E"
            E = read_raster_array(E_name)
            if S.shape != shape or E.shape != shape:
                gs.fatal(
                    "Animation plan rasters must match extent of before/after rasters."
                )

        else:
            gs.fatal(f"Unknown primitive: '{primitive}'")

    # ---- Print plan statistics ----------------------------------------------
    gs.message("Animation plan stats:")
    gs.message(f"  S:  min={S.min():.3f}  max={S.max():.3f}  mean={S.mean():.3f}")
    gs.message(f"  E:  min={E.min():.3f}  max={E.max():.3f}  mean={E.mean():.3f}")

    if plan_only:
        gs.message("Flag -p set: skipping frame generation.")
        _maybe_save_plan(S, E, options, shape, overwrite)
        return

    # ---- Optionally save animation plan rasters -----------------------------
    _maybe_save_plan(S, E, options, shape, overwrite)

    # ---- Generate animation frames ------------------------------------------
    gs.message(f"Generating {n_frames} animation frames ...")

    frame_times = np.linspace(0.0, 1.0, n_frames)
    pad = len(str(n_frames))

    for fi, t in enumerate(frame_times):
        frame_label = str(fi + 1).zfill(pad)
        gs.percent(fi, n_frames, 2)

        for bi, (before_arr, after_arr) in enumerate(zip(before_arrays, after_arrays)):
            frame = render_frame(before_arr, after_arr, S, E, t)

            if n_bands == 1:
                out_name = f"{output}_{frame_label}"
            else:
                out_name = f"{output}_b{bi + 1}_{frame_label}"

            write_raster_array(frame, out_name, overwrite=overwrite)
            gs.run_command(
                "r.support",
                map=out_name,
                quiet=True,
                title=f"Baia frame {frame_label} (t={t:.3f})",
                description=f"r.anim.morph {primitive} primitive, "
                f"band {bi + 1}, frame {fi + 1}/{n_frames}",
            )

    gs.percent(n_frames, n_frames, 2)
    gs.message(f"Done. Generated {n_frames} frame(s) with prefix '{output}'.")

    if n_bands == 1:
        gs.message(
            f"  Frames: {output}_{str(1).zfill(pad)} .. "
            f"{output}_{str(n_frames).zfill(pad)}"
        )
    else:
        gs.message(
            f"  Frames: {output}_b1_{str(1).zfill(pad)} .. "
            f"{output}_b{n_bands}_{str(n_frames).zfill(pad)}"
        )


def _maybe_save_plan(S, E, options, shape, overwrite):
    """Save S and E matrices as rasters if output_plan is specified."""
    if options.get("output_plan"):
        import numpy as np

        prefix = options["output_plan"]
        gs.message(f"Saving animation plan to '{prefix}_S' and '{prefix}_E' ...")
        write_raster_array(S, prefix + "_S", overwrite=overwrite)
        write_raster_array(E, prefix + "_E", overwrite=overwrite)
        gs.run_command(
            "r.support",
            map=prefix + "_S",
            quiet=True,
            title="Baia animation plan: S (start times)",
            units="normalized [0,1]",
        )
        gs.run_command(
            "r.support",
            map=prefix + "_E",
            quiet=True,
            title="Baia animation plan: E (end times)",
            units="normalized [0,1]",
        )


def build_plan_from_def(sd, shape, before_arrays, after_arrays, invert, blend_dur):
    """
    Build an (S, E) plan from a JSON stage definition dict.
    Used for multi-stage animation plans loaded from a JSON file.

    Expected JSON structure for each stage:
    {
        "primitive": "contraction",
        "mask_before": "lake_mask_2000",
        "mask_after":  "lake_mask_2010",
        "roi_start": 0.0,
        "roi_end":   0.6,
        "bg_start":  0.6,
        "bg_end":    1.0,
        "blend_duration": 0.2,
        "direction": "N",    // for directional
        "dem": "aster_dem",  // for dem
        "invert": false
    }
    """
    import numpy as np

    prim = sd.get("primitive", "blend").lower()
    roi_win = (float(sd.get("roi_start", 0.0)), float(sd.get("roi_end", 1.0)))
    bg_win = (float(sd.get("bg_start", 0.0)), float(sd.get("bg_end", 1.0)))
    bd = float(sd.get("blend_duration", blend_dur))
    inv = sd.get("invert", invert)

    mb = read_raster_array(sd["mask_before"]) if sd.get("mask_before") else None
    ma = read_raster_array(sd["mask_after"]) if sd.get("mask_after") else None

    if prim == "blend":
        return plan_blend(shape, roi_win, bg_win)
    elif prim == "appearance":
        return plan_appearance(ma, roi_win, bg_win, bd)
    elif prim == "disappearance":
        return plan_disappearance(mb, roi_win, bg_win, bd)
    elif prim == "contraction":
        return plan_contraction(mb, ma, roi_win, bg_win, bd)
    elif prim == "expansion":
        return plan_expansion(mb, ma, roi_win, bg_win, bd)
    elif prim == "deformation":
        return plan_deformation(mb, ma, roi_win, bg_win, bd)
    elif prim == "radial":
        cx = shape[1] / 2.0
        cy = shape[0] / 2.0
        mask = mb if mb is not None else ma
        return plan_radial(shape, cx, cy, roi_win, bg_win, bd, mask=mask, invert=inv)
    elif prim == "directional":
        mask = mb if mb is not None else ma
        return plan_directional(
            shape, sd.get("direction", "N"), roi_win, bg_win, bd, mask=mask, invert=inv
        )
    elif prim == "dem":
        dem_arr = read_raster_array(sd["dem"])
        mask = mb if mb is not None else ma
        return plan_dem(dem_arr, roi_win, bg_win, bd, mask=mask, invert=inv)
    else:
        gs.fatal(f"Unknown primitive in stage definition: '{prim}'")


if __name__ == "__main__":
    main()
