#!/usr/bin/env python3

############################################################################
#
# MODULE:       r.out.stereogram
# AUTHOR:       Corey T. White, OpenPlains Inc.
# PURPOSE:      Get items from a STAC API server
# COPYRIGHT:    (C) 2023-2024 Corey White
#               This program is free software under the GNU General
#               Public License (>=v2). Read the file COPYING that
#               comes with GRASS for details.
#
#############################################################################

# %module
# % description: Create Magic Eye stereograms from SVG or text.
# % keyword: raster
# % keyword: stereogram
# %end

# %option
# % key: input
# % type: string
# % description: Text input for stereogram
# % required: no
# % guisection: Stereogram
# %end

# %option G_OPT_F_INPUT
# % key: svg
# % type: string
# % description: SVG file path input for stereogram
# % required: no
# % guisection: Stereogram
# %end

# %option G_OPT_R_INPUT
# % key: input_raster
# % type: string
# % description: Input depth raster for stereogram
# % required: no
# % guisection: Stereogram
# %end

# %option
# % key: pattern_size
# % type: integer
# % description: Size of the random pattern in pixels
# % required: no
# % answer: 60
# % guisection: Stereogram
# %end

# %option G_OPT_M_COLR
# % key: color_scheme
# % description: Color scheme for the stereogram
# % required: no
# % guisection: Stereogram
# %end

# %option G_OPT_F_OUTPUT
# % key: output
# % description: Name of the output stereogram png
# % required: yes
# %end

# %option G_OPT_R_INPUT
# % key: raster
# % description: Name of the base raster for dimensions
# % required: yes
# %end

# %flag
# % key: r
# % description: Save output as raster in addition to PNG
# % guisection: Stereogram
# %end

from __future__ import annotations
import grass.script as gs
import grass.script.array as garray
from grass.script import raster as grast
from grass.exceptions import CalledModuleError
import numpy as np
import numpy.typing as npt
import cv2
from io import BytesIO
from PIL import Image
from pathlib import Path
import sys
import uuid
import atexit

TMP_RASTER_NAME = None


def cleanup():
    """Cleanup temporary raster"""
    if TMP_RASTER_NAME:
        gs.run_command("g.remove", flags="f", type="raster", name=TMP_RASTER_NAME)
        gs.message(_("Temporary raster '%s' removed.") % TMP_RASTER_NAME)


def lazy_import_module(module_name: str):
    try:
        module = __import__(module_name)
        return module
    except ImportError as e:
        gs.fatal(f"Failed to import {module_name}: {e}")


def generate_temp_raster_name(raster_name: str) -> str:
    """Generate a temporary raster name"""
    uuid_str = str(uuid.uuid4())
    TMP_RASTER_NAME = f"tmp_{raster_name}_{uuid_str}"
    gs.debug(_("Temporary raster name: %s") % TMP_RASTER_NAME)
    return TMP_RASTER_NAME


def svg_to_image(svg_path: Path, width: int, height: int):
    cairosvg = lazy_import_module("cairosvg")
    if cairosvg is None:
        gs.fatal("CairoSVG is not available. Please install it to use SVG input.")
    if not svg_path.with_suffix(".svg"):
        gs.fatal("Input file must have a '.svg' extension.")
    if not svg_path.exists():
        gs.fatal(f"SVG file '{svg_path}' does not exist.")
    png_data = cairosvg.svg2png(
        url=str(svg_path), output_width=width, output_height=height
    )
    image = Image.open(BytesIO(png_data)).convert("L")
    return np.array(image)


def text_to_image(text, width, height):
    image = np.zeros((height, width), dtype=np.uint8)
    font = cv2.FONT_HERSHEY_SIMPLEX
    scale = min(width, height) / 200
    thickness = int(scale)
    textsize = cv2.getTextSize(text, font, scale, thickness)[0]
    textX = (image.shape[1] - textsize[0]) // 2
    textY = (image.shape[0] + textsize[1]) // 2
    cv2.putText(image, text, (textX, textY), font, scale, 255, thickness)
    cv2.imwrite(".depth_text.png", image)
    return image


def generate_stereogram(base_raster, depth_map, pattern_size=60):
    gs.message(_("Generating stereogram..."))
    stereogram = garray.array(base_raster).astype(np.uint8)
    h, w = depth_map.shape
    pattern = np.random.randint(0, 256, (h, pattern_size), dtype=np.uint8)
    offsets = (depth_map / 255.0 * pattern_size / 2).astype(int)
    x_indices = (np.arange(w) - offsets) % pattern_size
    stereogram[:] = pattern[np.arange(h)[:, None], x_indices]

    return stereogram


def write_stereogram_raster(stereogram: npt.NDArray[np.uint8], output: str) -> str:
    """Write the stereogram to a raster file"""
    gs.message(_("Writing stereogram raster..."))
    new_raster = garray.array()
    # calculate new map from input map and store as GRASS raster map
    new_raster[...] = stereogram.astype(np.uint8)

    # Strip the extension and path from the output name
    raster_output = Path(output).name

    output_name = (
        raster_output if raster_output else generate_temp_raster_name(raster_output)
    )
    new_raster.write(output_name)
    return output_name


def set_color_map(stereogram: str, base_raster: str, color_scheme: str | None = None):
    """Set the color map for the stereogram"""
    gs.message(_("Writing colormap..."))

    def __get_color_map_scale_offset(base_raster):
        try:
            base_info = grast.raster_info(base_raster)
            min_val = base_info["min"]
            max_val = base_info["max"]
            gs.message(_("Base raster min: %s, max: %s") % (min_val, max_val))
            scale = float(255 / (max_val - min_val))
            offset = float(-min_val * scale)
            gs.message(_("Color map scale: %s, offset: %s") % (scale, offset))
            return (scale, offset)
        except Exception as e:
            gs.fatal(_("Error getting color map scale and offset: %s") % e)

    try:
        if color_scheme:
            gs.run_command("r.colors", map=stereogram, color=color_scheme)
        else:
            scale, offset = __get_color_map_scale_offset(base_raster)
            gs.run_command(
                "r.colors",
                map=stereogram,
                raster=base_raster,
                offset=offset,
                scale=scale,
            )

    except CalledModuleError as e:
        gs.fatal(_("Error setting color map: %s") % e.stderr)


def export_as_png(stereogram_name: str, output: str):
    gs.message(_("Exporting stereogram as PNG..."))
    png_path = Path(output)
    if not png_path.with_suffix(".png"):
        gs.warning("Output file must have a '.png' extension.")
        png_path = Path(f"{output}.png")

    try:
        gs.run_command(
            "r.out.png",
            input=stereogram_name,
            output=png_path,
            flags="w",
            overwrite=True,
        )
    except CalledModuleError as e:
        gs.fatal(_("Error saving png: %s") % e.stderr)


def main():
    options, flags = gs.parser()

    # Check for input options
    input_text = options.get("input")
    svg_input = options.get("svg")
    input_raster = options.get("input_raster")
    base_raster = options["raster"]

    # Stereogram parameters
    color_scheme = options.get("color_scheme")
    pattern_size = int(options.get("pattern_size", 60))

    # Output parameters
    output_png = options["output"]
    output_raster = flags["r"]

    # Get the region dimensions
    region = gs.region()
    width = region["cols"]
    height = region["rows"]

    if input_raster:
        depth_image = garray.array(input_raster).astype(np.uint8)
    elif svg_input:
        svg_path = Path(svg_input)
        depth_image = svg_to_image(svg_path, width, height)
    elif input_text:
        depth_image = text_to_image(input_text, width, height)
    else:
        gs.fatal("Either input text or SVG image must be provided.")

    stereogram = generate_stereogram(base_raster, depth_image, pattern_size)

    output_raster_name = write_stereogram_raster(stereogram, output_png, output_raster)
    set_color_map(output_raster_name, base_raster, color_scheme)

    export_as_png(output_raster_name, output_png)

    gs.message(f"Magic eye stereogram '{output_png}' created successfully.")


if __name__ == "__main__":
    options, flags = gs.parser()
    atexit.register(cleanup)
    sys.exit(main())
