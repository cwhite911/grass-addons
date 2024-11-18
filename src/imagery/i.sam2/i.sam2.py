#!/usr/bin/env python3

############################################################################
#
# MODULE:       i.sam2
# AUTHOR:       Corey T. White, OpenPlains Inc.
# PURPOSE:      Uses the SAMGeo model for segmentation in GRASS GIS.
# COPYRIGHT:    (C) 2023-2024 Corey White
#               This program is free software under the GNU General
#               Public License (>=v2). Read the file COPYING that
#               comes with GRASS for details.
#
#############################################################################

# %module
# % description: Integrates SAMGeo model with text prompt for segmentation in GRASS GIS.
# % keyword: imagery
# % keyword: segmentation
# % keyword: object recognition
# % keyword: deep learning
# %end

# %option G_OPT_I_GROUP
# % key: group
# % description: Name of input imagery group
# % required: yes
# %end

# %option
# % key: output
# % type: string
# % description: Name of output segmented raster map
# % gisprompt: new,cell,raster
# % required: yes
# %end

# %option
# % key: model_path
# % type: string
# % description: Path to the SAMGeo model file (optional if using default model)
# % required: no
# %end

# %option
# % key: text_prompt
# % type: string
# % description: Optional text prompt to guide segmentation
# % required: yes
# %end

# %option
# % key: text_threshold
# % type: double
# % answer: 0.24
# % description: Text threshold for text segmentation
# % required: no
# % multiple: no
# %end

# %option
# % key: box_threshold
# % type: double
# % answer: 0.24
# % description: Box threshold for text segmentation
# % required: no
# % multiple: no
# %end


import os
import sys
import grass.script as gs
import torch
import numpy as np
from PIL import Image
from grass.script import array as garray


def main():
    from samgeo import SamGeo
    from samgeo.text_sam import LangSAM

    group = options["group"]

    output_raster = options["output"]
    model_path = options.get("model_path")
    text_prompt = options.get("text_prompt")
    text_threshold = float(options.get("text_threshold"))
    box_threshold = float(options.get("box_threshold"))

    # Set default model path if not provided
    default_model_path = os.path.join(
        gs.gisenv()["GISDBASE"], "samgeo_default_model.pth"
    )
    if not model_path:
        model_path = default_model_path

    # Set up paths to access the raster files
    tmp_dir = gs.tempdir()
    temp_input_path = os.path.join(tmp_dir, "input.tif")

    temp_output_path = tmp_dir
    guide_input_path = None

    rasters = gs.read_command("i.group", group=group, flags="lg").strip().split("\n")
    input_image_np = list([garray.array(raster, dtype=np.uint8) for raster in rasters])

    rgb_array = np.stack(input_image_np, axis=-1)

    if rgb_array.dtype != np.uint8:
        rgb_array = (
            (rgb_array - rgb_array.min()) / (rgb_array.max() - rgb_array.min()) * 255
        )
        rgb_array = rgb_array.astype(np.uint8)

    np_image = Image.fromarray(rgb_array[:, :, :3])

    # Get device
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        torch.cuda.empty_cache()

    try:
        if text_prompt:
            gs.message("Running LangSAM segmentation...")
            sam = LangSAM(
                model_type="sam2-hiera-large",
            )
            from torch.amp.autocast_mode import autocast

            with autocast(device_type=device):
                masks, boxes, phrases, logits = sam.predict(
                    image=np_image,
                    text_prompt=text_prompt,
                    box_threshold=box_threshold,
                    text_threshold=text_threshold,
                    return_results=True,
                )
        else:
            gs.message("Running SAMGeo segmentation...")
            sam = SamGeo(model_type="vit_h", model_path=model_path, device=device)
            sam.generate(
                input_image=temp_input_path,
                output=temp_output_path,
                guide=guide_input_path,
            )
    except Exception as e:
        gs.message(torch.cuda.memory_summary())
        gs.fatal(f"Error while running SAMGeo: {e}")
        return 1

    gs.message("Segmentation complete.")

    write_raster(masks, output_raster)
    return 0


def write_raster(input_np_array, output_raster):
    gs.message("Importing the segmented raster into GRASS GIS...")

    if len(input_np_array) == 0:
        gs.fatal("No masks found.")

    # Initialize the merged raster with zeros
    merged_raster = np.zeros_like(input_np_array[0], dtype=np.int32)

    for idx, band in enumerate(input_np_array):
        if band.shape != input_np_array[0].shape:
            gs.fatal("All masks must have the same shape.")

        unique_value = idx + 1  # Start unique values from 1

        # Use NumPy's vectorized operations to assign unique values
        mask = band != 0  # Create a mask where band is not zero
        merged_raster[mask] = unique_value

    # Convert the merged raster to a GRASS array
    mask_raster = garray.array()
    mask_raster[...] = merged_raster

    # Write the merged raster to the output file
    mask_raster.write(mapname=output_raster)


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main())
