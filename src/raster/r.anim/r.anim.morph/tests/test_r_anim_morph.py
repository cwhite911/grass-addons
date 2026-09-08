#!/usr/bin/env python3
"""
r.anim.morph test & demo script
Validates all animation plan primitives and generates demo GIF animations.
Runs standalone (no GRASS required) using numpy/scipy/Pillow.

Usage:
    python test_r_anim_morph.py

Outputs PNG frame sequences and animated GIFs in ./demo_output/
"""

import os
import sys
import numpy as np
from scipy import ndimage

# ── Require Pillow for image output ─────────────────────────────────────────
try:
    from PIL import Image
except ImportError:
    sys.exit("Pillow required: pip install Pillow")

OUT_DIR = "./demo_output"
os.makedirs(OUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# Core Baia rendering engine (mirror of r.anim.morph.py, no GRASS deps)
# ─────────────────────────────────────────────────────────────────────────────


def render_frame(before, after, S, E, t):
    dur = np.maximum(E - S, 1e-10)
    alpha = np.clip((t - S) / dur, 0.0, 1.0)
    return (1.0 - alpha) * before + alpha * after


def plan_blend(shape, roi_window=(0, 1), bg_window=(0, 1)):
    return np.full(shape, float(roi_window[0])), np.full(shape, float(roi_window[1]))


def plan_appearance(mask_after, roi_window=(0, 1), bg_window=(0, 1), blend_dur=0.3):
    shape = mask_after.shape
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))
    roi = mask_after > 0
    dist = ndimage.distance_transform_edt(roi)
    mx = dist.max()
    norm = dist / mx if mx > 0 else np.zeros(shape)
    roi_dur = roi_window[1] - roi_window[0]
    S[roi] = roi_window[0] + norm[roi] * (roi_dur - blend_dur)
    E[roi] = np.clip(S[roi] + blend_dur, 0, 1)
    return S, E


def plan_disappearance(mask_before, roi_window=(0, 1), bg_window=(0, 1), blend_dur=0.3):
    shape = mask_before.shape
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))
    roi = mask_before > 0
    dist = ndimage.distance_transform_edt(roi)
    mx = dist.max()
    norm = dist / mx if mx > 0 else np.zeros(shape)
    roi_dur = roi_window[1] - roi_window[0]
    S[roi] = roi_window[0] + (1.0 - norm[roi]) * (roi_dur - blend_dur)
    E[roi] = np.clip(S[roi] + blend_dur, 0, 1)
    return S, E


def plan_contraction(
    mask_before, mask_after, roi_window=(0, 1), bg_window=(0, 1), blend_dur=0.3
):
    shape = mask_before.shape
    mb, ma = mask_before > 0, mask_after > 0
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))
    intersection = mb & ma
    disappearing = mb & ~ma
    S[intersection] = roi_window[0]
    E[intersection] = roi_window[1]
    if disappearing.any():
        d2a = (
            ndimage.distance_transform_edt(~ma)
            if ma.any()
            else ndimage.distance_transform_edt(mb)
        )
        d = d2a[disappearing]
        dmax = d.max()
        norm = d / dmax if dmax > 0 else np.zeros_like(d, float)
        roi_dur = roi_window[1] - roi_window[0]
        S[disappearing] = np.clip(
            roi_window[0] + (1 - norm) * (roi_dur - blend_dur), 0, 1
        )
        E[disappearing] = np.clip(S[disappearing] + blend_dur, 0, 1)
    return S, E


def plan_expansion(
    mask_before, mask_after, roi_window=(0, 1), bg_window=(0, 1), blend_dur=0.3
):
    shape = mask_after.shape
    mb = mask_before > 0 if mask_before is not None else np.zeros(shape, bool)
    ma = mask_after > 0
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))
    intersection = mb & ma
    appearing = ma & ~mb
    S[intersection] = roi_window[0]
    E[intersection] = roi_window[1]
    if appearing.any():
        d2b = (
            ndimage.distance_transform_edt(~mb)
            if mb.any()
            else ndimage.distance_transform_edt(ma)
        )
        d = d2b[appearing]
        dmax = d.max()
        norm = d / dmax if dmax > 0 else np.zeros_like(d, float)
        roi_dur = roi_window[1] - roi_window[0]
        S[appearing] = np.clip(roi_window[0] + norm * (roi_dur - blend_dur), 0, 1)
        E[appearing] = np.clip(S[appearing] + blend_dur, 0, 1)
    return S, E


def plan_deformation(
    mask_before, mask_after, roi_window=(0, 1), bg_window=(0, 1), blend_dur=0.3
):
    shape = mask_before.shape
    mb, ma = mask_before > 0, mask_after > 0
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))
    intersection = mb & ma
    S[intersection] = roi_window[0]
    E[intersection] = roi_window[1]
    # Contracting
    contracting = mb & ~ma
    if contracting.any():
        d2a = (
            ndimage.distance_transform_edt(~ma)
            if ma.any()
            else ndimage.distance_transform_edt(mb)
        )
        d = d2a[contracting]
        dmax = d.max()
        norm = d / dmax if dmax > 0 else np.zeros_like(d, float)
        roi_dur = roi_window[1] - roi_window[0]
        S[contracting] = np.clip(
            roi_window[0] + (1 - norm) * (roi_dur - blend_dur), 0, 1
        )
        E[contracting] = np.clip(S[contracting] + blend_dur, 0, 1)
    # Expanding
    expanding = ma & ~mb
    if expanding.any():
        d2b = (
            ndimage.distance_transform_edt(~mb)
            if mb.any()
            else ndimage.distance_transform_edt(ma)
        )
        d = d2b[expanding]
        dmax = d.max()
        norm = d / dmax if dmax > 0 else np.zeros_like(d, float)
        roi_dur = roi_window[1] - roi_window[0]
        S[expanding] = np.clip(roi_window[0] + norm * (roi_dur - blend_dur), 0, 1)
        E[expanding] = np.clip(S[expanding] + blend_dur, 0, 1)
    return S, E


def plan_radial(
    shape,
    cx=None,
    cy=None,
    roi_window=(0, 1),
    bg_window=(0, 1),
    blend_dur=0.3,
    mask=None,
    invert=False,
):
    rows, cols = shape
    if cx is None:
        cx = cols / 2.0
    if cy is None:
        cy = rows / 2.0
    Y, X = np.mgrid[0:rows, 0:cols]
    dist = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    mx = dist.max()
    norm = dist / mx if mx > 0 else np.zeros(shape)
    if invert:
        norm = 1 - norm
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))
    roi = mask > 0 if mask is not None else np.ones(shape, bool)
    roi_dur = roi_window[1] - roi_window[0]
    S[roi] = roi_window[0] + norm[roi] * (roi_dur - blend_dur)
    E[roi] = np.clip(S[roi] + blend_dur, 0, 1)
    return S, E


def plan_directional(
    shape,
    direction="N",
    roi_window=(0, 1),
    bg_window=(0, 1),
    blend_dur=0.3,
    mask=None,
    invert=False,
):
    rows, cols = shape
    Y, X = np.mgrid[0:rows, 0:cols]
    dir_map = {
        "N": -Y,
        "S": Y,
        "E": X,
        "W": -X,
        "NE": -Y + X,
        "NW": -Y - X,
        "SE": Y + X,
        "SW": Y - X,
    }
    raw = dir_map[direction.upper()]
    mn, mx = raw.min(), raw.max()
    norm = (raw - mn) / (mx - mn) if mx > mn else np.zeros(shape)
    if invert:
        norm = 1 - norm
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))
    roi = mask > 0 if mask is not None else np.ones(shape, bool)
    roi_dur = roi_window[1] - roi_window[0]
    S[roi] = roi_window[0] + norm[roi] * (roi_dur - blend_dur)
    E[roi] = np.clip(S[roi] + blend_dur, 0, 1)
    return S, E


def plan_dem(
    dem_arr, roi_window=(0, 1), bg_window=(0, 1), blend_dur=0.1, mask=None, invert=False
):
    shape = dem_arr.shape
    valid = np.isfinite(dem_arr)
    mn, mx = dem_arr[valid].min(), dem_arr[valid].max()
    norm = np.zeros(shape)
    if mx > mn:
        norm[valid] = (dem_arr[valid] - mn) / (mx - mn)
    if not invert:
        norm = 1 - norm  # high elevation -> low S (transitions early by default)
    S = np.full(shape, float(bg_window[0]))
    E = np.full(shape, float(bg_window[1]))
    roi = mask > 0 if mask is not None else np.ones(shape, bool)
    roi_dur = roi_window[1] - roi_window[0]
    S[roi] = roi_window[0] + norm[roi] * (roi_dur - blend_dur)
    E[roi] = np.clip(S[roi] + blend_dur, 0, 1)
    return S, E


# ─────────────────────────────────────────────────────────────────────────────
# Demo helpers
# ─────────────────────────────────────────────────────────────────────────────


def arr_to_rgb_img(arr, vmin=None, vmax=None):
    """Convert float array [0-255 or arbitrary] to a PIL RGB image."""
    if vmin is None:
        vmin = arr.min()
    if vmax is None:
        vmax = arr.max()
    if vmax > vmin:
        norm = np.clip((arr - vmin) / (vmax - vmin), 0, 1)
    else:
        norm = np.zeros_like(arr)
    rgb = (norm * 255).astype(np.uint8)
    return Image.fromarray(np.stack([rgb, rgb, rgb], axis=-1))


def save_gif(frames_pil, path, duration_ms=80):
    frames_pil[0].save(
        path, save_all=True, append_images=frames_pil[1:], loop=0, duration=duration_ms
    )
    print(f"  Saved: {path}")


def generate_and_save(name, before, after, S, E, n_frames=40, colormap=None):
    """Generate animation frames and save as animated GIF and plan PNGs."""
    print(f"\n[{name}]")
    times = np.linspace(0, 1, n_frames)
    pil_frames = []
    vmin = min(before.min(), after.min())
    vmax = max(before.max(), after.max())

    for t in times:
        frame = render_frame(before, after, S, E, t)
        pil_frames.append(arr_to_rgb_img(frame, vmin, vmax))

    save_gif(pil_frames, os.path.join(OUT_DIR, f"{name}.gif"))

    # Save plan visualization (S and E as grayscale images)
    Image.fromarray((S * 255).clip(0, 255).astype(np.uint8)).save(
        os.path.join(OUT_DIR, f"{name}_S.png")
    )
    Image.fromarray((E * 255).clip(0, 255).astype(np.uint8)).save(
        os.path.join(OUT_DIR, f"{name}_E.png")
    )
    print("  Plan (S, E) saved as PNG")

    # Quick stats
    print(f"  S: min={S.min():.3f}  max={S.max():.3f}  mean={S.mean():.3f}")
    print(f"  E: min={E.min():.3f}  max={E.max():.3f}  mean={E.mean():.3f}")


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic scenes for demonstration
# ─────────────────────────────────────────────────────────────────────────────


def make_lake_scene(rows=200, cols=300):
    """
    Simulate a shrinking lake (CONTRACTION scenario).
    Before: large lake.   After: smaller lake.
    Colors: deep blue = deep water, light blue = shallow, green = land.
    """
    # Background: green land
    before = np.full((rows, cols, 3), [60, 120, 40], dtype=float)
    after = np.full((rows, cols, 3), [65, 125, 45], dtype=float)

    # Large lake (before)
    Y, X = np.ogrid[:rows, :cols]
    cx_b, cy_b, rx_b, ry_b = cols // 2, rows // 2, cols // 3, rows // 3
    lake_b = ((X - cx_b) / rx_b) ** 2 + ((Y - cy_b) / ry_b) ** 2 <= 1.0

    # Small lake (after), shifted and smaller, simulating retreat
    cx_a, cy_a, rx_a, ry_a = cols // 2 + 10, rows // 2 - 5, cols // 5, rows // 5
    lake_a = ((X - cx_a) / rx_a) ** 2 + ((Y - cy_a) / ry_a) ** 2 <= 1.0

    # Paint water
    depth_b = np.sqrt(
        np.maximum(1.0 - ((X - cx_b) / rx_b) ** 2 - ((Y - cy_b) / ry_b) ** 2, 0)
    )
    before[lake_b] = np.column_stack(
        [
            20 + 30 * (1 - depth_b[lake_b]),
            100 + 50 * depth_b[lake_b],
            200 + 40 * depth_b[lake_b],
        ]
    )

    depth_a = np.sqrt(
        np.maximum(1.0 - ((X - cx_a) / rx_a) ** 2 - ((Y - cy_a) / ry_a) ** 2, 0)
    )
    after[lake_a] = np.column_stack(
        [
            20 + 30 * (1 - depth_a[lake_a]),
            100 + 50 * depth_a[lake_a],
            200 + 40 * depth_a[lake_a],
        ]
    )

    return before, after, (lake_b.astype(float), lake_a.astype(float))


def make_flood_scene(rows=200, cols=300):
    """Simulate a river flooding its floodplain (EXPANSION)."""
    before = np.full((rows, cols, 3), [100, 150, 80], dtype=float)
    after = np.full((rows, cols, 3), [90, 140, 75], dtype=float)
    Y, X = np.mgrid[:rows, :cols]
    # Narrow river before
    river = abs(X - cols // 2) < 15
    before[river] = [30, 80, 180]
    after[river] = [30, 80, 180]
    # Flood zone after (wider)
    flood = abs(X - cols // 2) < 60
    after[flood & ~river] = [50, 100, 170]
    return before, after, (river.astype(float), flood.astype(float))


def make_snow_scene(rows=200, cols=300):
    """Simulate DEM-controlled snow accumulation."""
    Y, X = np.mgrid[:rows, :cols]
    # Synthetic DEM: mountain in center
    dem = (
        np.exp(
            -(
                (X - cols // 2) ** 2 / (cols // 3) ** 2
                + (Y - rows // 2) ** 2 / (rows // 3) ** 2
            )
        )
        * 3000
    )

    # Before = summer (green valleys, grey peaks)
    before = np.zeros((rows, cols, 3), dtype=float)
    elev_norm = (dem - dem.min()) / (dem.max() - dem.min())
    before[..., 0] = 80 - 40 * elev_norm
    before[..., 1] = 120 - 60 * elev_norm
    before[..., 2] = 60 - 20 * elev_norm

    # After = winter (white snow everywhere)
    after = np.full((rows, cols, 3), 220.0)
    after[..., 0] -= 20 * (1 - elev_norm)
    after[..., 1] -= 10 * (1 - elev_norm)

    return before, after, dem


def main():
    print("=" * 60)
    print("r.anim.morph Animation Plan Demo")
    print("=" * 60)
    print(f"Output directory: {os.path.abspath(OUT_DIR)}\n")

    N_FRAMES = 50

    # ── 1. BLEND ─────────────────────────────────────────────────────────────
    before_l, after_l, (mb_l, ma_l) = make_lake_scene()
    shape = before_l.shape[:2]

    S, E = plan_blend(shape)
    before_g = before_l.mean(axis=2)
    after_g = after_l.mean(axis=2)
    generate_and_save("01_blend", before_g, after_g, S, E, N_FRAMES)

    # ── 2. CONTRACTION (shrinking lake) ───────────────────────────────────────
    S, E = plan_contraction(
        mb_l, ma_l, roi_window=(0.0, 0.75), bg_window=(0.75, 1.0), blend_dur=0.2
    )
    generate_and_save("02_contraction_lake", before_g, after_g, S, E, N_FRAMES)

    # ── 3. CONTRACTION with staging (ROI then BG) ─────────────────────────────
    S, E = plan_contraction(
        mb_l,
        ma_l,
        roi_window=(0.0, 0.6),
        bg_window=(0.0, 0.0),  # BG stays frozen during ROI
        blend_dur=0.15,
    )
    # Compose with BG stage
    S_bg, E_bg = plan_blend(shape, roi_window=(0.0, 0.0), bg_window=(0.6, 1.0))
    # BG pixels are where ma and mb are both zero
    bg_mask = (mb_l == 0) & (ma_l == 0)
    S_combined = S.copy()
    E_combined = E.copy()
    S_combined[bg_mask] = S_bg[bg_mask]
    E_combined[bg_mask] = E_bg[bg_mask]
    generate_and_save(
        "03_contraction_staged", before_g, after_g, S_combined, E_combined, N_FRAMES
    )

    # ── 4. EXPANSION (flooding) ───────────────────────────────────────────────
    before_f, after_f, (mb_f, ma_f) = make_flood_scene()
    shape_f = before_f.shape[:2]
    before_fg = before_f.mean(axis=2)
    after_fg = after_f.mean(axis=2)
    S, E = plan_expansion(
        mb_f, ma_f, roi_window=(0.0, 0.8), bg_window=(0.0, 0.0), blend_dur=0.2
    )
    generate_and_save("04_expansion_flood", before_fg, after_fg, S, E, N_FRAMES)

    # ── 5. DEFORMATION ─────────────────────────────────────────────────────────
    # Shift lake slightly and resize for deformation
    Y, X = np.ogrid[: shape[0], : shape[1]]
    rows, cols = shape
    mb_d = (
        ((X - cols // 2) / (cols // 4)) ** 2 + ((Y - rows // 2) / (rows // 4)) ** 2 <= 1
    ).astype(float)
    ma_d = (
        ((X - cols // 2 + 20) / (cols // 5)) ** 2
        + ((Y - rows // 2 - 10) / (rows // 3)) ** 2
        <= 1
    ).astype(float)
    S, E = plan_deformation(
        mb_d, ma_d, roi_window=(0.0, 0.85), bg_window=(0.85, 1.0), blend_dur=0.15
    )
    generate_and_save("05_deformation", before_g, after_g, S, E, N_FRAMES)

    # ── 6. APPEARANCE ─────────────────────────────────────────────────────────
    before_blank = np.full(shape, 80.0)
    S, E = plan_appearance(
        ma_l, roi_window=(0.0, 0.8), bg_window=(0.8, 1.0), blend_dur=0.2
    )
    generate_and_save("06_appearance", before_blank, before_g, S, E, N_FRAMES)

    # ── 7. DISAPPEARANCE ──────────────────────────────────────────────────────
    S, E = plan_disappearance(
        mb_l, roi_window=(0.0, 0.8), bg_window=(0.8, 1.0), blend_dur=0.2
    )
    after_blank = np.full(shape, 80.0)
    generate_and_save("07_disappearance", before_g, after_blank, S, E, N_FRAMES)

    # ── 8. RADIAL ─────────────────────────────────────────────────────────────
    S, E = plan_radial(
        shape, roi_window=(0.0, 0.9), bg_window=(0.9, 1.0), blend_dur=0.1
    )
    generate_and_save("08_radial", before_g, after_g, S, E, N_FRAMES)

    # ── 9. DIRECTIONAL (N) ───────────────────────────────────────────────────
    S, E = plan_directional(
        shape, "N", roi_window=(0.0, 0.9), bg_window=(0.9, 1.0), blend_dur=0.1
    )
    generate_and_save("09_directional_N", before_g, after_g, S, E, N_FRAMES)

    # ── 10. DIRECTIONAL (SE) ─────────────────────────────────────────────────
    S, E = plan_directional(
        shape, "SE", roi_window=(0.0, 0.9), bg_window=(0.9, 1.0), blend_dur=0.1
    )
    generate_and_save("10_directional_SE", before_g, after_g, S, E, N_FRAMES)

    # ── 11. DEM-based (snow from peaks) ──────────────────────────────────────
    before_s, after_s, dem = make_snow_scene()
    before_sg = before_s.mean(axis=2)
    after_sg = after_s.mean(axis=2)
    S, E = plan_dem(dem, roi_window=(0.0, 0.9), blend_dur=0.1, invert=False)
    generate_and_save("11_dem_snow_peaks_first", before_sg, after_sg, S, E, N_FRAMES)

    # ── 12. DEM-based inverted (snow from valleys) ────────────────────────────
    S, E = plan_dem(dem, roi_window=(0.0, 0.9), blend_dur=0.1, invert=True)
    generate_and_save("12_dem_snow_valleys_first", before_sg, after_sg, S, E, N_FRAMES)

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("All demos complete.")
    print(f"Animated GIFs and plan PNGs saved to: {os.path.abspath(OUT_DIR)}/")
    print("\nDemo inventory:")
    demos = [
        ("01_blend.gif", "Monolithic blend (baseline)"),
        ("02_contraction_lake.gif", "Contraction – shrinking lake"),
        ("03_contraction_staged.gif", "Contraction – staged (ROI then BG)"),
        ("04_expansion_flood.gif", "Expansion – river flooding"),
        ("05_deformation.gif", "Deformation – shape change"),
        ("06_appearance.gif", "Appearance – object materializes"),
        ("07_disappearance.gif", "Disappearance – object vanishes"),
        ("08_radial.gif", "Radial progression from center"),
        ("09_directional_N.gif", "Directional progression – N"),
        ("10_directional_SE.gif", "Directional progression – SE"),
        ("11_dem_snow_peaks_first.gif", "DEM – snow from peaks downward"),
        ("12_dem_snow_valleys_first.gif", "DEM – snow from valleys upward"),
    ]
    for fname, desc in demos:
        print(f"  {fname:<40}  {desc}")


if __name__ == "__main__":
    main()
