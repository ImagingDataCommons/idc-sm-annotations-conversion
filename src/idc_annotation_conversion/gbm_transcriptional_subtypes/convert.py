import logging
from time import time
from typing import Union
from collections import Counter

import numpy as np
import pydicom
from pydicom.uid import (
    ExplicitVRLittleEndian,
    JPEG2000Lossless,
)
import pandas as pd
import highdicom as hd

from idc_annotation_conversion.gbm_transcriptional_subtypes import metadata_config


def convert_segmentation(
    segmentation_dataframe: pd.DataFrame,
    source_image: pydicom.Dataset,
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str],
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str],
    include_lut: bool = False,
) -> hd.seg.Segmentation:
    """Convert segmentation dataframe to DICOM segmentation.

    Parameters
    ----------
    source_image: pydicom.Dataset
        Pydicom dataset containing the metadata of the image (already converted
        to DICOM format). This can be the full image datasets, but the
        PixelData attribute is not required.
    segmentation_dataframe: pandas.DataFrame
        Patch-level subtype information contained in a dataframe mapping
        coordinates to subtype.
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
            segmented_property_category=metadata_config.finding_codes[label][0],
            segmented_property_type=metadata_config.finding_codes[label][1],
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
            algorithm_identification=metadata_config.algorithm_identification,
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

    source_image = hd.Image.from_dataset(source_image, copy=False)

    sorted_x_vals = segmentation_dataframe.sort_values("x").x.unique()
    diff_x = np.diff(sorted_x_vals)
    sorted_y_vals = segmentation_dataframe.sort_values("y").y.unique()
    diff_y = np.diff(sorted_y_vals)
    x_start = sorted_x_vals[0]
    y_start = sorted_y_vals[0]

    spacing_x = diff_x[0]
    spacing_y = diff_y[0]

    assert np.all(diff_x % spacing_x == 0)
    assert np.all(diff_y % spacing_y == 0)

    x_indices = ((segmentation_dataframe.x.values - x_start) / spacing_x).round().astype(np.uint32)
    y_indices = ((segmentation_dataframe.y.values - y_start) / spacing_y).round().astype(np.uint32)

    shape = (1, y_indices.max() + 1, x_indices.max() + 1)
    mask = np.zeros(shape, dtype=np.uint8)

    for x, y, t in zip(x_indices, y_indices, segmentation_dataframe.cell_type):
        mask[0, y, x] = metadata_config.segmentation_channel_order.index(t) + 1

    source_geometry = source_image.get_volume_geometry()

    # Account for the shift in position due to the resampling, assuming that
    # the corner of the two images remains the same
    # TODO factor this into highdicom as resampling
    new_pix_spacing = 0.056  # from description in paper (56 microns)
    new_position = (
        np.array(source_geometry.position) +
        (new_pix_spacing / 2.0 - source_geometry.spacing[1] / 2.0) * np.array(source_geometry.unit_vectors()[1]) +
        (new_pix_spacing / 2.0 - source_geometry.spacing[2] / 2.0) * np.array(source_geometry.unit_vectors()[2])
    )

    mask_volume = hd.Volume.from_components(
        direction=source_geometry.direction,
        position=new_position,
        spacing=np.array([1.0, new_pix_spacing, new_pix_spacing]),
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
        series_number=20,
        sop_instance_uid=hd.UID(),
        instance_number=1,
        manufacturer=metadata_config.seg_manufacturer,
        manufacturer_model_name=metadata_config.seg_manufacturer_model_name,
        software_versions=metadata_config.software_versions,
        device_serial_number=metadata_config.device_serial_number,
        transfer_syntax_uid=transfer_syntax_uid,
        dimension_organization_type=dimension_organization_type,
        omit_empty_frames=omit_empty_frames,
        series_description=metadata_config.segmentation_series_description,
        palette_color_lut_transformation=lut_transform,
    )
    seg_time = time() - seg_start_time
    logging.info(f"Created DICOM Segmentation in {seg_time:.1f}s.")

    return segmentation
