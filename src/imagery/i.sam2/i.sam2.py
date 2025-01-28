#!/usr/bin/env python3

############################################################################
#
# MODULE:       i.sam2
# AUTHOR:       Corey T. White, OpenPlains Inc.
# PURPOSE:      Uses the SAMGeo model for segmentation in GRASS GIS.
# COPYRIGHT:    (C) 2023-2025 Corey White
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
# % description: Text prompt to guide segmentation
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

DEFAULT_MODEL_PATH = "samgeo_default_model.pth"


def get_device():
    """
    Determines the available device for computation (CUDA or CPU).

    This function checks if a CUDA-enabled GPU is available and returns "cuda" if it is,
    otherwise it returns "cpu". If CUDA is available, it also clears the CUDA cache.

    Returns:
        str: "cuda" if a CUDA-enabled GPU is available, otherwise "cpu".
    """
    device = "cuda" if torch.cuda.is_available() else "cpu"
    gs.message(_(f"Running computation on {device}..."))
    if device == "cuda":
        torch.cuda.empty_cache()
    return device


def read_raster_group(group):
    """
    Reads a group of raster maps and returns them as a list of numpy arrays.

    Parameters:
    group (str): The name of the raster group to read.

    Returns:
    list: A list of numpy arrays, each representing a raster map in the group.
    """
    gs.message(_("Reading imagery group..."))
    rasters = gs.read_command("i.group", group=group, flags="lg")
    # raster_list = [raster.strip() for raster in rasters.splitlines()]
    raster_list = map(str.split, rasters.splitlines())
    return [garray.array(raster, dtype=np.uint8) for raster in raster_list]


def normalize_rgb_array(rgb_array):
    """
    Normalize an RGB array to the range [0, 255].

    This function takes an RGB array and normalizes its values to the range [0, 255].
    If the input array is not of type np.uint8, it scales the values to fit within
    this range and converts the array to np.uint8.

    Parameters:
    rgb_array (numpy.ndarray): The input RGB array to be normalized.

    Returns:
    numpy.ndarray: The normalized RGB array with values in the range [0, 255] and type np.uint8.
    """
    if rgb_array.dtype != np.uint8:
        gs.message(_("Converting RGB array to uint8..."))
        min_val = rgb_array.min()
        max_val = rgb_array.max()
        scale = 255 / (max_val - min_val)
        rgb_array = ((rgb_array - min_val) * scale).astype(np.uint8)
        # rgb_array = ((rgb_array - rgb_array.min()) / (rgb_array.max() - rgb_array.min()) * 255).astype(np.uint8)
    return rgb_array


def run_langsam_segmentation(
    np_image, text_prompt, box_threshold, text_threshold, device
):
    """
    Runs LangSAM segmentation on the given image using the specified text prompt and thresholds.

    Parameters:
    np_image (numpy.ndarray): The input image as a NumPy array.
    text_prompt (str): The text prompt to guide the segmentation.
    box_threshold (float): The threshold for box predictions.
    text_threshold (float): The threshold for text predictions.
    device (str): The device to run the segmentation on (e.g., 'cpu' or 'cuda').

    Returns:
    list: A list of masks generated by the segmentation.
    """
    from samgeo.text_sam import LangSAM
    from torch.amp.autocast_mode import autocast

    gs.message(_("Running LangSAM segmentation..."))
    sam = LangSAM(model_type="sam2-hiera-large")
    with autocast(device_type=device):
        masks, boxes, phrases, logits = sam.predict(
            image=np_image,
            text_prompt=text_prompt,
            box_threshold=box_threshold,
            text_threshold=text_threshold,
            return_results=True,
        )
    return masks


def run_samgeo_segmentation(temp_input_path, temp_output_path, model_path, device):
    """
    Run SAMGeo segmentation on an input image and save the output.

    Parameters:
    temp_input_path (str): The file path to the input image.
    temp_output_path (str): The file path to save the segmented output image.
    model_path (str): The file path to the SAMGeo model.
    device (str): The device to run the model on (e.g., 'cpu', 'cuda').

    Returns:
    None
    """
    from samgeo import SamGeo

    gs.message(_("Running SAMGeo segmentation..."))
    sam = SamGeo(model_type="vit_h", model_path=model_path, device=device)
    sam.generate(input_image=temp_input_path, output=temp_output_path)


def write_raster(input_np_array, output_raster):
    """
    Writes a segmented raster into GRASS GIS.

    Parameters:
    input_np_array (list of numpy.ndarray): A list of numpy arrays representing the input raster bands.
    output_raster (str): The name of the output raster map to be created in GRASS GIS.

    Raises:
    ValueError: If the input array is empty or if the masks do not have the same shape.

    This function merges multiple raster bands into a single raster, where each band is assigned a unique value.
    The merged raster is then written to a GRASS GIS raster map.
    """

    gs.message(_("Importing the segmented raster into GRASS GIS..."))

    if len(input_np_array) == 0:
        gs.fatal("No masks found.")

    merged_raster = np.zeros_like(input_np_array[0], dtype=np.int32)
    for idx, band in enumerate(input_np_array):
        if band.shape != input_np_array[0].shape:
            gs.fatal(_("All masks must have the same shape."))
        unique_value = idx + 1
        mask = band != 0
        merged_raster[mask] = unique_value

    mask_raster = garray.array()
    mask_raster[...] = merged_raster
    mask_raster.write(mapname=output_raster)


def main():
    group = options["group"]
    output_raster = options["output"]
    model_path = options.get("model_path") or os.path.join(
        gs.gisenv()["GISDBASE"], DEFAULT_MODEL_PATH
    )
    text_prompt = options.get("text_prompt")
    text_threshold = float(options.get("text_threshold"))
    box_threshold = float(options.get("box_threshold"))

    tmp_dir = gs.tempdir()
    temp_input_path = os.path.join(tmp_dir, "input.tif")
    temp_output_path = tmp_dir

    input_image_np = read_raster_group(group)
    rgb_array = normalize_rgb_array(np.stack(input_image_np, axis=-1))
    np_image = Image.fromarray(rgb_array[:, :, :3])

    device = get_device()

    try:
        if text_prompt:
            masks = run_langsam_segmentation(
                np_image, text_prompt, box_threshold, text_threshold, device
            )
        else:
            run_samgeo_segmentation(
                temp_input_path, temp_output_path, model_path, device
            )
    except Exception as e:
        gs.fatal(_(f"Error while running SAMGeo: {e}"))
        return 1

    gs.message(_("Segmentation complete."))
    write_raster(masks, output_raster)
    return 0


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main())
