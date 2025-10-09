from collections import Counter
import logging
from time import time
from typing import Union, TextIO

from highdicom.volume import ChannelDescriptor
import numpy as np
import pydicom
from pydicom.uid import (
    ExplicitVRLittleEndian,
    JPEG2000Lossless,
)
import pandas as pd
import highdicom as hd

from idc_annotation_conversion.tcga_til_maps import metadata_config


def convert_segmentation(
    segmentation_array: np.ndarray,
    source_image: pydicom.Dataset,
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str],
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str],
    include_lut: bool = False,
) -> hd.seg.Segmentation:
    """Store segmentation masks as DICOM segmentations.

    Parameters
    ----------
    source_image: pydicom.Dataset
        Pydicom dataset containing the metadata of the image (already converted
        to DICOM format). This can be the full image datasets, but the
        PixelData attribute is not required.
    segmentation_array: np.ndarray
        Segmentation mask image as the original RGB array.
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str], optional
        Segmentation type (BINARY or FRACTIONAL) for the Segmentation Image
        (if any).
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str], optional
        Dimension organization type of the output segmentations.
    include_lut: bool, optional
        Whether to include a LUT transformation in the segmentation to store as
        a PALETTE COLOR image. Ignored if segmentation_type is not LABELMAP.

    Returns
    -------
    segmentations: hd.seg.Segmentation
        DICOM segmentation image encoding the original annotations

    """
    seg_start_time = time()

    if include_lut and segmentation_type == hd.seg.SegmentationTypeValues.LABELMAP:
        lut_transform = metadata_config.labelmap_lut
    else:
        lut_transform = None

    container_id = source_image.ContainerIdentifier
    segment_descriptions = []
    for (number, label) in enumerate(
        metadata_config.segmentation_channel_order_2018,
        start=1
    ):
        color = (
            metadata_config.display_colors[label]
            if lut_transform is None else None
        )

        desc = hd.seg.SegmentDescription(
            segment_number=number,
            segment_label=label,
            segmented_property_category=metadata_config.finding_codes_2018[label][0],
            segmented_property_type=metadata_config.finding_codes_2018[label][1],
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
            algorithm_identification=metadata_config.algorithm_identification_2018,
            tracking_id=f"{container_id}-{label}",
            tracking_uid=hd.UID(),
            display_color=color,
        )
        segment_descriptions.append(desc)

    segmentation_type = hd.seg.SegmentationTypeValues(segmentation_type)
    dimension_organization_type = hd.DimensionOrganizationTypeValues(
        dimension_organization_type
    )

    # Compression method depends on what is possible given the chosen
    # segmentation type
    transfer_syntax_uid = {
        hd.seg.SegmentationTypeValues.BINARY: ExplicitVRLittleEndian,
        hd.seg.SegmentationTypeValues.FRACTIONAL: JPEG2000Lossless,
        hd.seg.SegmentationTypeValues.LABELMAP: JPEG2000Lossless,
    }[segmentation_type]

    omit_empty_frames = dimension_organization_type.value != "TILED_FULL"

    mask = np.stack(
        [
            segmentation_array[:, :, 2] > 0,  # B channel encodes TIL negative
            segmentation_array[:, :, 0] > 0,  # R channel encodes TIL positive
        ],
        axis=-1
    )[None]

    source_image = hd.Image.from_dataset(source_image, copy=False)

    source_geometry = source_image.get_volume_geometry()

    # Account for the shift in position due to the resampling, assuming that
    # the corner of the two images remains the same
    # TODO factor this into highdicom as resampling
    new_pix_spacing = 0.05  # from description in paper
    # new_position = (
    #     np.array(source_geometry.position) +
    #     (new_pix_spacing / 2.0 - source_geometry.spacing[1] / 2.0) * np.array(source_geometry.unit_vectors()[1]) +
    #     (new_pix_spacing / 2.0 - source_geometry.spacing[2] / 2.0) * np.array(source_geometry.unit_vectors()[2])
    # )
    new_position = np.array(source_geometry.position)

    mask_volume = hd.Volume.from_components(
        direction=source_geometry.direction,
        position=new_position,
        spacing=np.array([1.0, new_pix_spacing, new_pix_spacing]),
        coordinate_system="SLIDE",
        array=mask,
        channels={ChannelDescriptor('SegmentNumber'): [1, 2]},
    )

    logging.info("Creating DICOM Segmentation")
    seg_start_time = time()
    segmentation = hd.seg.Segmentation(
        source_images=[source_image],
        pixel_array=mask_volume,
        segmentation_type=segmentation_type,
        segment_descriptions=segment_descriptions,
        series_instance_uid=hd.UID(),
        series_number=20,
        sop_instance_uid=hd.UID(),
        instance_number=1,
        manufacturer=metadata_config.seg_manufacturer_2018,
        manufacturer_model_name=metadata_config.seg_manufacturer_model_name_2018,
        software_versions=metadata_config.software_versions,
        device_serial_number=metadata_config.device_serial_number,
        transfer_syntax_uid=transfer_syntax_uid,
        dimension_organization_type=dimension_organization_type,
        omit_empty_frames=omit_empty_frames,
        series_description=metadata_config.segmentation_series_description_2018,
        palette_color_lut_transformation=lut_transform,
    )
    segmentation.add(metadata_config.other_trials_seq_element)
    seg_time = time() - seg_start_time
    logging.info(f"Created DICOM Segmentation in {seg_time:.1f}s.")

    return segmentation


def convert_txt_file(
    text_file: TextIO,
    source_image: hd.Image,
    container_id: str,
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str],
    is_probability: bool = False,
    dimension_organization_type: Union[
        hd.DimensionOrganizationTypeValues,
        str
    ] = "TILED_FULL",
) -> hd.seg.Segmentation:
    """Convert a TIL map in text file form to a DICOM segmentation.

    Parameters
    ----------
    text_file: TextIO
        Text file or file-like object containing the TIL map information.
    source_image: hd.Image
        Source image for the segmentation as a DICOM image.
    container_id: str
        Container identifier of the source image.
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str]
        Segmentation type (must be FRACTIONAL if `is_probability` is True).
    is_probability: bool = False
        Whether the text file contains a binary or probability TIL map.
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str]
        Dimension organization type to use in the Segmentation image.

    Returns
    -------
    highdicom.seg.Segmentation:
        Converted segmentation image.

    """
    segmentation_type = hd.seg.SegmentationTypeValues(segmentation_type)
    dimension_organization_type = hd.DimensionOrganizationTypeValues(
        dimension_organization_type
    )

    if is_probability:
        series_number = 22

        if segmentation_type != hd.seg.SegmentationTypeValues.FRACTIONAL:
            raise ValueError("Fractional type required for fractional segmentations")

        transfer_syntax_uid = JPEG2000Lossless
        array_dtype = np.float64
        series_description = metadata_config.segmentation_series_description_2022_fractional
    else:
        series_number = 21
        if segmentation_type == hd.seg.SegmentationTypeValues.LABELMAP:
            transfer_syntax_uid = JPEG2000Lossless
        else:
            transfer_syntax_uid = ExplicitVRLittleEndian
        array_dtype = bool
        series_description = metadata_config.segmentation_series_description_2022_binary

    omit_empty_frames = dimension_organization_type.value != "TILED_FULL"

    df = pd.read_csv(
        text_file,
        sep=" ",
        names=["x", "y", "value", "0"],
    )

    sorted_x_vals = df.sort_values("x").x.unique()
    diff_x = np.diff(sorted_x_vals)
    sorted_y_vals = df.sort_values("y").y.unique()
    diff_y = np.diff(sorted_y_vals)
    x_start = sorted_x_vals[0]
    y_start = sorted_y_vals[0]

    # Due to resampling, the coordinates of the segmentation are not quite
    # regularly spaced
    spacing_x = Counter(diff_x).most_common(1)[0][0]
    spacing_y = Counter(diff_y).most_common(1)[0][0]

    assert abs(spacing_x // 2 - x_start) < 2
    assert abs(spacing_y // 2 - y_start) < 2

    # Check the spacing is "near" regular: the gap between either neighbouring
    # pixels is either the most common spacing value, or within 1 pixel of it
    if np.any(np.abs(diff_x - spacing_x) > 1):
        raise RuntimeError("Missing pixels")
    if np.any(np.abs(diff_y - spacing_y) > 1):
        raise RuntimeError("Missing pixels")

    x_indices = ((df.x.values - x_start) / spacing_x).round().astype(np.uint32)
    y_indices = ((df.y.values - y_start) / spacing_y).round().astype(np.uint32)

    shape = (1, y_indices.max() + 1, x_indices.max() + 1)
    mask = np.zeros(shape, dtype=array_dtype)

    for x, y, v in zip(x_indices, y_indices, df.value):
        if not is_probability:
            v = bool(v)

        mask[0, y, x] = v

    # Create volume with correct spatial metadata
    source_geometry = source_image.get_volume_geometry()

    # Account for the shift in position due to the resampling, assuming that
    # the corner of the two images remains the same
    # TODO factor this into highdicom as resampling
    new_pix_spacing = 0.05  # from description in paper
    # new_position = (
    #     np.array(source_geometry.position) +
    #     (new_pix_spacing / 2.0 - source_geometry.spacing[1] / 2.0) * np.array(source_geometry.unit_vectors()[1]) +
    #     (new_pix_spacing / 2.0 - source_geometry.spacing[2] / 2.0) * np.array(source_geometry.unit_vectors()[2])
    # )
    new_position = np.array(source_geometry.position)

    mask_volume = hd.Volume.from_components(
        direction=source_geometry.direction,
        position=new_position,
        spacing=np.array([1.0, new_pix_spacing, new_pix_spacing]),
        coordinate_system="SLIDE",
        array=mask,
    )

    # If the image size is smaller than 32 in either dimension, there will be a
    # compression error, so switch to uncompressed
    if mask.shape[0] < 33 or mask.shape[2] < 33:
        transfer_syntax_uid = ExplicitVRLittleEndian

    segment_descriptions = []
    for (number, label) in enumerate(
        metadata_config.segmentation_channel_order_2022,
        start=1
    ):
        color = (
            metadata_config.display_colors[label]
            if segmentation_type != hd.seg.SegmentationTypeValues.LABELMAP
            else None
        )

        desc = hd.seg.SegmentDescription(
            segment_number=number,
            segment_label=label,
            segmented_property_category=metadata_config.finding_codes_2022[label][0],
            segmented_property_type=metadata_config.finding_codes_2022[label][1],
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
            algorithm_identification=metadata_config.algorithm_identification_2022,
            tracking_id=f"{container_id}-{label}",
            tracking_uid=hd.UID(),
            display_color=color,
        )
        segment_descriptions.append(desc)

    logging.info("Creating DICOM Segmentation")
    seg_start_time = time()
    segmentation = hd.seg.Segmentation(
        source_images=[source_image],
        pixel_array=mask_volume,
        segmentation_type=segmentation_type,
        segment_descriptions=segment_descriptions,
        series_instance_uid=hd.UID(),
        series_number=series_number,
        sop_instance_uid=hd.UID(),
        instance_number=1,
        manufacturer=metadata_config.seg_manufacturer_2022,
        manufacturer_model_name=metadata_config.seg_manufacturer_model_name_2022,
        software_versions=metadata_config.software_versions,
        device_serial_number=metadata_config.device_serial_number,
        transfer_syntax_uid=transfer_syntax_uid,
        series_description=series_description,
        dimension_organization_type=dimension_organization_type,
        omit_empty_frames=omit_empty_frames,
    )
    segmentation.add(metadata_config.other_trials_seq_element)
    seg_time = time() - seg_start_time
    logging.info(f"Created DICOM Segmentation in {seg_time:.1f}s.")

    return segmentation
