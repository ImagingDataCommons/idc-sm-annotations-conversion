from collections import Counter
import logging
from time import time
from typing import Union, BinaryIO

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

    container_id = source_image.ContainerIdentifier
    segment_descriptions = []
    for (number, label) in enumerate(
        metadata_config.segmentation_channel_order,
        start=1
    ):
        desc = hd.seg.SegmentDescription(
            segment_number=number,
            segment_label=label,
            segmented_property_category=metadata_config.finding_codes_2018[label][0],
            segmented_property_type=metadata_config.finding_codes_2018[label][1],
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
            algorithm_identification=metadata_config.algorithm_identification_2018,
            tracking_id=f"{container_id}-{label}",
            tracking_uid=hd.UID(),
        )
        segment_descriptions.append(desc)

    segmentation_type = hd.seg.SegmentationTypeValues(segmentation_type)
    dimension_organization_type = hd.DimensionOrganizationTypeValues(
        dimension_organization_type
    )

    if include_lut and segmentation_type == hd.seg.SegmentationTypeValues.LABELMAP:
        lut_transform = metadata_config.labelmap_lut
    else:
        lut_transform = None

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

    mask_volume = hd.Volume.from_components(
        direction=source_geometry.direction,
        position=source_geometry.position,
        spacing=np.array([1.0, 0.05, 0.05]),  # TODO this is an approximation. CHECK
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
    seg_time = time() - seg_start_time
    logging.info(f"Created DICOM Segmentation in {seg_time:.1f}s.")

    return segmentation


def convert_txt_file(
    text_file: BinaryIO,
    source_image: pydicom.Dataset,
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str],
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str],
    is_probability: bool = False,
) -> hd.seg.Segmentation:

    if is_probability:
        series_number = 22
        segmentation_type = hd.seg.SegmentationTypeValues.FRACTIONAL
        transfer_syntax_uid = JPEG2000Lossless
        array_dtype = np.float64
    else:
        series_number = 21
        segmentation_type = hd.seg.SegmentationTypeValues.BINARY
        transfer_syntax_uid = ExplicitVRLittleEndian
        array_dtype = bool


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

    assert spacing_x == x_start * 2
    assert spacing_y == y_start * 2

    # Check the spacing is "near" regular: the gap between either neighbouring
    # pixels is either the most common spacing value, or within 1 pixel of it
    if np.any(np.abs(diff_x - spacing_x) > 1):
        raise RuntimeError("Missing pixels")
    if np.any(np.abs(diff_y - spacing_y) > 1):
        raise RuntimeError("Missing pixels")

    x_indices = ((df.x.values - x_start) / spacing_x).round().astype(np.uint32)
    y_indices = ((df.y.values - y_start) / spacing_y).round().astype(np.uint32)

    shape = (y_indices.max() + 1, x_indices.max() + 1)
    mask = np.zeros(shape, dtype=array_dtype)

    for x, y, v in zip(x_indices, y_indices, df.value):
        if not is_probability:
            v = bool(v)

        mask[y, x] = bool(v)

    # TODO create volume with correct spatial metadata
    source_geometry = source_image.get_volume_geometry()

    mask_volume = hd.Volume.from_components(
        direction=source_geometry.direction,
        position=source_geometry.position,
        spacing=np.array([1.0, 0.05, 0.05]),  # TODO this is an approximation. CHECK
        coordinate_system="SLIDE",
        array=mask,
    )

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
        manufacturer=metadata_config.seg_manufacturer,
        manufacturer_model_name=metadata_config.seg_manufacturer_model_name,
        software_versions=metadata_config.software_versions,
        device_serial_number=metadata_config.device_serial_number,
        transfer_syntax_uid=transfer_syntax_uid,
        dimension_organization_type=dimension_organization_type,
        series_description=metadata_config.segmentation_series_description,
    )
    seg_time = time() - seg_start_time
    logging.info(f"Created DICOM Segmentation in {seg_time:.1f}s.")

    return segmentation
