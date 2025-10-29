import logging
from time import time
from typing import cast, Union

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
    hghdicom.seg.Segmentation
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
        metadata_config.segmentation_channel_order,
        start=1
    ):
        color = (
            metadata_config.display_colors[label]
            if lut_transform is None else None
        )

        desc = hd.seg.SegmentDescription(
            segment_number=number,
            segment_label=label,
            segmented_property_category=metadata_config.finding_codes[label][0],
            segmented_property_type=metadata_config.finding_codes[label][1],
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
            algorithm_identification=metadata_config.algorithm_identification,
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
    assert spacing_x == spacing_y

    x_indices = ((segmentation_dataframe.x.values - x_start) / spacing_x).round().astype(np.uint32)
    y_indices = ((segmentation_dataframe.y.values - y_start) / spacing_y).round().astype(np.uint32)

    shape = (1, y_indices.max() + 1, x_indices.max() + 1)
    mask = np.zeros(shape, dtype=np.uint8)

    for x, y, t in zip(x_indices, y_indices, segmentation_dataframe.cell_type):
        mask[0, y, x] = metadata_config.segmentation_channel_order.index(t) + 1

    source_geometry = cast(hd.VolumeGeometry, source_image.get_volume_geometry())

    # Account for the shift in position due to the resampling, assuming that
    # the corner of the two images remains the same
    # TODO factor this into highdicom as resampling
    new_pix_spacing = spacing_x * source_geometry.pixel_spacing[0]
    # new_position = (
    #     np.array(source_geometry.position) +
    #     (new_pix_spacing / 2.0 - source_geometry.spacing[1] / 2.0) * np.array(source_geometry.unit_vectors()[1]) +
    #     (new_pix_spacing / 2.0 - source_geometry.spacing[2] / 2.0) * np.array(source_geometry.unit_vectors()[2])
    # )

    mask_volume = hd.Volume.from_components(
        direction=source_geometry.direction,
        position=source_geometry.position,
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

    segmentation.add(metadata_config.other_trials_seq_element)

    return segmentation


def convert_aggressiveness_map(
    scores: np.ndarray,
    coords: np.ndarray,
    source_image: pydicom.Dataset,
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str],
) -> hd.pm.ParametricMap:
    """Convert segmentation dataframe to DICOM segmentation.

    Parameters
    ----------
    scores: numpy.ndarray
        1-d numpy array giving aggressiveness scores.
    coords: numpy.ndarray
        2-d numpy array of shape (n, 2), where coords[i, :] gives the [x, y]
        total pixel matrix coordinates of scores[i].
    source_image: pydicom.Dataset
        Pydicom dataset containing the metadata of the image (already converted
        to DICOM format). This can be the full image datasets, but the
        PixelData attribute is not required.
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str], optional
        Dimension organization type of the output parametric map.

    Returns
    -------
    highdicom.pmap.ParametricMap
        DICOM parametric map image encoding the original annotations

    """
    dimension_organization_type = hd.DimensionOrganizationTypeValues(
        dimension_organization_type
    )
    source_image_hd = hd.Image.from_dataset(source_image, copy=False)

    sorted_x_vals = np.unique(coords[:, 0])
    sorted_y_vals = np.unique(coords[:, 1])
    diff_x = np.diff(sorted_x_vals)
    diff_y = np.diff(sorted_y_vals)

    spacing_x = diff_x.min()
    spacing_y = diff_y.min()

    assert np.all(diff_x % spacing_x == 0)
    assert np.all(diff_y % spacing_y == 0)
    assert spacing_x == spacing_y

    x_indices = ((coords[:, 0]) / spacing_x).round().astype(np.uint32)
    y_indices = ((coords[:, 1]) / spacing_y).round().astype(np.uint32)

    shape = (1, y_indices.max() + 1, x_indices.max() + 1)
    score_map = np.zeros(shape, dtype=scores.dtype)

    for x, y, score in zip(x_indices, y_indices, scores):
        score_map[0, y, x] = score

    source_geometry = cast(hd.VolumeGeometry, source_image_hd.get_volume_geometry())

    new_pix_spacing = spacing_x * source_geometry.pixel_spacing[0]

    position_sequence = hd.PlanePositionSequence(
        coordinate_system="SLIDE",
        image_position=source_geometry.position,
        pixel_matrix_position=[1, 1],
    )
    measures_sequence = hd.PixelMeasuresSequence(
        pixel_spacing=[new_pix_spacing] * 2,
        slice_thickness=source_geometry.spacing_between_slices,
    )

    pmap = hd.pm.ParametricMap(
        source_images=[source_image_hd],
        pixel_array=score_map,
        series_instance_uid=hd.UID(),
        sop_instance_uid=hd.UID(),
        series_number=21,
        instance_number=1,
        manufacturer=metadata_config.pmap_manufacturer,
        manufacturer_model_name=metadata_config.pmap_manufacturer_model_name,
        device_serial_number=metadata_config.pmap_device_serial_number,
        software_versions=metadata_config.pmap_software_versions,
        series_description=metadata_config.pmap_series_description,
        contains_recognizable_visual_features=False,
        real_world_value_mappings=metadata_config.pmap_real_world_value_mappings,
        window_center=0.5,
        window_width=1.0,
        content_label=metadata_config.pmap_content_label,
        content_description=metadata_config.pmap_content_description,
        pixel_measures=measures_sequence,
        plane_orientation=source_geometry.get_plane_orientation(),
        plane_positions=[position_sequence],
        contributing_equipment=metadata_config.pmap_contributing_equipment,
    )

    # TODO fold this into highdicom
    pmap.ImageType[2] = metadata_config.pmap_image_flavor
    pmap.ImageType[3] = metadata_config.derived_pixel_contrast

    # Current highdicom release (0.27.0) doesn't deal with dimension
    # organization type.  This should be fixed in the next release of highdicom
    pmap.DimensionOrganizationType = dimension_organization_type.value
    if (
        dimension_organization_type ==
        hd.DimensionOrganizationTypeValues.TILED_FULL
    ):
        if "DimensionIndexSequence" in pmap:
            del pmap.DimensionIndexSequence

        if "PerFrameFunctionalGroupsSequence" in pmap:
            del pmap.PerFrameFunctionalGroupsSequence

        pmap.TotalPixelMatrixFocalPlanes = 1

    pmap.add(metadata_config.other_trials_seq_element)

    return pmap
