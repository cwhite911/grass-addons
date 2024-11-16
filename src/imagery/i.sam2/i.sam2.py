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
# % description: Integrates SAMGeo model for segmentation in GRASS GIS.
# % keyword: raster
# % keyword: segmentation
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

# %flag
# % key: u
# % description: Update the default SAMGeo model
# %end

import os
import sys
import grass.script as gs
import torch
import requests
import numpy as np
from PIL import Image
from PIL.Image import Image as PILImage
from grass.pygrass import raster as r
from grass.script import array as garray


def update_model(default_model_path):
    gs.message("Updating the SAMGeo model...")
    model_url = "https://example.com/samgeo/default_model.pth"  # Replace with the actual model URL
    response = requests.get(model_url)
    if response.status_code == 200:
        with open(default_model_path, "wb") as model_file:
            model_file.write(response.content)
        gs.message("Model updated successfully.")
    else:
        gs.fatal(
            "Failed to update the model. Please check the URL or your internet connection."
        )


def main():
    from samgeo import SamGeo
    from samgeo.text_sam import LangSAM

    gisenv = gs.gisenv()
    group = options["group"]

    output_raster = options["output"]
    model_path = options.get("model_path")
    text_prompt = options.get("text_prompt")
    text_threshold = float(options.get("text_threshold"))
    box_threshold = float(options.get("box_threshold"))
    update_flag = flags["u"]

    # Set default model path if not provided
    default_model_path = os.path.join(
        gs.gisenv()["GISDBASE"], "samgeo_default_model.pth"
    )
    if not model_path:
        model_path = default_model_path

    # Update model if flag is set
    if update_flag:
        update_model(default_model_path)

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

    image_size = rgb_array.shape

    np_image = Image.fromarray(rgb_array[:, :, :3])
    image_np = np.array(np_image)

    gs.message(
        f"np_image type: {type(np_image)}, {isinstance(np_image, PILImage)}, image_np size: {image_np.size}"
    )

    # Get device
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        torch.cuda.empty_cache()

    try:
        if text_prompt:
            gs.message("Running LangSAM segmentation...")
            sam = LangSAM(
                # model_type="vit_h",
                model_type="sam2-hiera-large",
                # checkpoint=model_path
            )
            from torch.amp.autocast_mode import autocast

            with autocast(device_type=device):
                masks, boxes, phrases, logits = sam.predict(
                    image=np_image,
                    # output=temp_output_path,
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
        # Check GPU memory usage
        # allocated = torch.cuda.memory_allocated()
        # reserved = torch.cuda.memory_reserved()
        # free = torch.cuda.get_device_properties(0).total_memory - reserved

        # gs.message(f"Allocated Memory: {allocated / 1024**2:.2f} MB")
        # gs.message(f"Reserved Memory: {reserved / 1024**2:.2f} MB")
        # gs.message(f"Free Memory: {free / 1024**2:.2f} MB")
        return 1

    gs.message("Segmentation complete.")

    write_raster(masks, output_raster)
    return 0


def write_raster(input_np_array, output_raster):
    gs.message("Importing the segmented raster into GRASS GIS...")
    if input_np_array.shape[0] == 1:
        gs.message("Writing single-band raster...")
        mask_raster = garray.array()
        for y in range(mask_raster.shape[0]):
            for x in range(mask_raster.shape[1]):
                mask_raster[y][x] = input_np_array[0][y][x]
        mask_raster.write(mapname=output_raster)
    else:
        gs.message("Writing multi-band raster...")
        # TODO: Merge the bands into a single raster
        for idx, band in enumerate(input_np_array):
            if band.shape != input_np_array[0].shape:
                gs.fatal("All bands must have the same shape.")
            mask_raster = garray.array()
            for y in range(mask_raster.shape[0]):
                for x in range(mask_raster.shape[1]):
                    mask_raster[y][x] = band[y][x]

            mask_raster.write(mapname=f"{output_raster}.{idx}")


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main())
