## DESCRIPTION

*r.anim.morph* implements the **Baia animation plan framework** ([Lobo, Appert &
Pietriga 2018](https://doi.org/10.1109/TVCG.2018.2796557)) for GRASS. Given a
*before* and *after* raster (or multi-band imagery group defined as a
comma-separated list), it computes a pixel-level **animation plan**: two rasters
S and E that specify, for each pixel, when its transition from before to after
begins and ends, and outputs a numbered sequence of raster frames representing
the animated transition.

The key insight of the Baia model is that meaningful before-and-after animations
require **different pixels to transition at different times**. For example, when
illustrating a shrinking lake, pixels close to the final shoreline should
disappear last, while pixels far from the final shoreline should disappear
first, giving the visual impression of the lake actually shrinking rather than
simply fading.

### Animation Plan Model

For each pixel (i,j) and animation time t ∈ [0,1]:

```text
α(i,j,t) = clamp( (t − S[i,j]) / (E[i,j] − S[i,j]) , 0, 1 )
frame(i,j,t) = (1 − α) · before(i,j) + α · after(i,j)
```

- S[i,j] = 0: pixel starts blending immediately
- E[i,j] = 1: pixel finishes blending at the very end
- S[i,j] = E[i,j]: instantaneous swap at that moment

### Animation Primitives

- **blend**  
  Monolithic blend: all pixels transition uniformly from start to end.
  Equivalent to classic cross-dissolve. Requires no masks.

- **appearance**  
  An entity appears in the after image that was not present before. Pixels
  inside *mask_after* fade in progressively from edge to center. Requires
  *mask_after*.

- **disappearance**  
  An entity present in the before image vanishes. Pixels inside *mask_before*
  fade out progressively from edge inward. Requires *mask_before*.

- **contraction**  
  A shape shrinks from *mask_before* to *mask_after* (e.g. a drying lake).
  Pixels farthest from the new contour transition first; pixels closest to the
  new contour transition last, giving the impression of the shape gradually
  pulling inward. Requires *mask_before*; *mask_after* optional (if absent,
  shape contracts to nothing).

- **expansion**  
  A shape grows from *mask_before* to *mask_after* (e.g. flooding). Pixels
  closest to the old contour appear first; pixels farthest appear last. Requires
  *mask_after*; *mask_before* optional.

- **deformation**  
  Shape changes form: combines contraction of pixels leaving and expansion of
  pixels entering. Both *mask_before* and *mask_after* are required.

- **radial**  
  Animation radiates outward from a center point. Pixels closest to the center
  transition first. Center defaults to image center but can be set with
  *radial_x* / *radial_y* (map coordinates). Use *-i* flag to invert (animate
  from outside in). An optional mask restricts the effect.

- **directional**  
  Linear progression along a compass bearing (N/S/E/W/NE/NW/SE/SW). Pixels at
  the leading edge transition first. Use *-i* to reverse direction.

- **dem**  
  Transition timing derived from a Digital Elevation Model. By default,
  high-elevation pixels transition first (e.g. snow accumulating from mountain
  peaks downward). Use *-i* to reverse (low elevations first, e.g. flooding).
  Requires *dem* raster.

- **plan**  
  Use a pre-computed animation plan. Supply the name of the S raster via *plan*;
  a companion raster named `<plan>_E` must also exist in the mapset. This allows
  animation plans derived from external data (simulations, expert masks, etc.)
  to drive the transition.

### Staging

Staging allows the Region of Interest (ROI) and the background to animate at
different times, preventing distracting background changes from competing with
the focal change. Use *roi_start*/*roi_end* and *bg_start*/*bg_end* to set
independent timing windows.

Example: animate ROI first, then background:

```text
roi_start=0.0  roi_end=0.5
bg_start=0.5   bg_end=1.0
```

Example: animate ROI and background concurrently (default):

```text
roi_start=0.0  roi_end=1.0
bg_start=0.0   bg_end=1.0
```

### Multi-stage Animations (JSON)

For complex animations involving multiple sequential or parallel stages, supply
a JSON file via *stages*. Each entry in the JSON array defines one stage with
its own primitive, masks, and timing windows. The stages are composed by taking
the element-wise minimum of start times and maximum of end times across all
stages, allowing both sequential and overlapping transitions.

Example JSON (Aral Sea shrinkage, two stages: ROI contracts, then BG fades):

```json
[
  {
    "primitive": "contraction",
    "mask_before": "aral_sea_2000",
    "mask_after": "aral_sea_2010",
    "roi_start": 0.0,
    "roi_end": 0.7,
    "bg_start": 0.0,
    "bg_end": 0.0,
    "blend_duration": 0.2
  },
  {
    "primitive": "blend",
    "roi_start": 0.7,
    "roi_end": 1.0,
    "bg_start": 0.7,
    "bg_end": 1.0,
    "blend_duration": 0.3
  }
]
```

## NOTES

### Multi-band / RGB Imagery

Supply comma-separated band raster names for *before* and *after*. The same
animation plan is applied to all bands. Output frames are named
`<output>_b1_001`, `<output>_b2_001`, etc. Use *i.group* to assemble output
frames into imagery groups for display.

### Color Transfer (*-c* flag)

When before and after images were acquired under different illumination
conditions, the *-c* flag applies a simple mean/standard-deviation color
transfer to match the before image statistics to the after image before
computing the blend. This reduces distraction from color histogram differences.

### Saving Animation Plans

Use *output_plan* to save the computed S and E matrices as named rasters. These
can be re-used with `primitive=plan`, visualized in GRASS, or exported to
GeoTIFF and used in external tools. The paper's original implementation stores S
in band R and E in band G of a TIFF file; *r.anim.morph* stores them as separate
GRASS rasters for consistency with GRASS data model conventions.

### Creating Animations

Use *g.gui.animation* to play back the generated frame sequence as an animation
within GRASS. Alternatively, export frames with *r.out.png* and compose them
into an animated GIF or video with external tools such as `ffmpeg` or `convert`
(ImageMagick).

```sh
# Export frames to PNG
for frame in $(g.list type=rast pattern="anim_*" mapset=.); do
    r.out.png input=$frame output=${frame}.png
done

# Assemble into animated GIF (ImageMagick)
convert -delay 5 -loop 0 anim_*.png aral_sea.gif

# Assemble into MP4 (ffmpeg)
ffmpeg -framerate 15 -pattern_type glob -i 'anim_*.png' -c:v libx264 aral_sea.mp4
```

## EXAMPLES

### Simple monolithic blend (baseline)

```sh
r.anim.morph before=lake_2000 after=lake_2010 \
    output=anim primitive=blend frames=30
```

### Contracting lake (Aral Sea / shrinking lake scenario)

```sh
# Create binary masks with r.mapcalc or r.threshold
r.mapcalc "lake_mask_2000 = if(lake_extent_2000 > 0, 1, 0)"
r.mapcalc "lake_mask_2010 = if(lake_extent_2010 > 0, 1, 0)"

r.anim.morph before=lake_2000 after=lake_2010 \
    output=lake_anim \
    primitive=contraction \
    mask_before=lake_mask_2000 \
    mask_after=lake_mask_2010 \
    roi_start=0.0 roi_end=0.8 \
    bg_start=0.8  bg_end=1.0 \
    blend_duration=0.2 \
    frames=60 \
    output_plan=lake_plan
```

### Flood expansion

```sh
r.anim.morph before=area_pre_flood after=area_post_flood \
    output=flood_anim \
    primitive=expansion \
    mask_before=water_mask_before \
    mask_after=water_mask_after \
    roi_start=0.0 roi_end=0.9 \
    bg_start=0.9  bg_end=1.0 \
    frames=45
```

### DEM-based snow accumulation

```sh
r.anim.morph before=ndsi_summer after=ndsi_winter \
    output=snow_anim \
    primitive=dem \
    dem=aster_dem \
    blend_duration=0.1 \
    frames=30
```

### DEM-based snow melt (low elevations first)

```sh
r.anim.morph before=ndsi_winter after=ndsi_summer \
    output=melt_anim \
    primitive=dem \
    dem=aster_dem \
    blend_duration=0.1 \
    frames=30 -i
```

### Northward directional progression

```sh
r.anim.morph before=before_img after=after_img \
    output=dir_anim \
    primitive=directional \
    direction=N \
    frames=30
```

### Multi-band RGB Sentinel-2 animation

```sh
# Sentinel-2 bands: B4=red, B3=green, B2=blue
r.anim.morph \
    before=S2_2020_B4,S2_2020_B3,S2_2020_B2 \
    after=S2_2024_B4,S2_2024_B3,S2_2024_B2 \
    output=s2_anim \
    primitive=contraction \
    mask_before=water_mask_2020 \
    mask_after=water_mask_2024 \
    frames=60 -c

# Assemble into RGB imagery groups for display
for f in $(seq -w 1 60); do
    i.group group=frame_${f} \
        input=s2_anim_b1_${f},s2_anim_b2_${f},s2_anim_b3_${f}
done
```

### Multi-stage JSON animation

```sh
r.anim.morph before=before after=after \
    output=staged_anim \
    stages=/path/to/animation_stages.json \
    frames=60
```

### Visualize animation in GRASS

```sh
g.gui.animation strds=my_animation_strds

# Or using the frame list directly:
g.list type=rast pattern="lake_anim_*" mapset=. output=/tmp/frames.txt
g.gui.animation rast=$(cat /tmp/frames.txt | tr '\n' ',')
```

## REFERENCES

Lobo MJ, Appert C, Pietriga E (2019).  
*Animation Plans for Before-and-After Satellite Images.*  
IEEE Transactions on Visualization and Computer Graphics, 25(2):1347–1360.  
[doi:10.1109/TVCG.2018.2796557](https://doi.org/10.1109/TVCG.2018.2796557)

Claramunt C, Thériault M (1995).  
Managing time in GIS an event-oriented approach.  
In: Recent Advances in Temporal Databases. Springer.

## SEE ALSO

*[g.gui.animation](g.gui.animation.html),  
[i.group](i.group.html),  
[r.buffer](r.buffer.html),  
[r.grow](r.grow.html),  
[r.mapcalc](r.mapcalc.html),  
[r.out.png](r.out.png.html),  
[r.series](r.series.html),  
[r.anim](r.anim.html),  
[t.rast.series](t.rast.series.html)*

## AUTHOR

Implemented for GRASS from the Baia framework (Lobo, Appert & Pietriga 2018).
Corey T. White, Center for Geospatial Analytics, NC State University
