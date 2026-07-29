import highdicom as hd
import numpy as np
from typing import cast
from pydicom.uid import ExplicitVRLittleEndian

import pandas as pd
from idc_annotation_conversion.panoptils import metadata_config


# Mapping from graphic type ('type' column) in the CSVs to DICOM graphic type
GRAPHIC_TYPE_MAPPING = {
    "polyline": hd.ann.GraphicTypeValues.POLYGON,
    "rectangle": hd.ann.GraphicTypeValues.RECTANGLE,
    "point": hd.ann.GraphicTypeValues.POINT,
}


# The pixel spacing of the resampled images that the annotations and
# segmentations have been performed on. The documentation claims that this is 
# 0.25 microns per pixel, but my calculations revealed that it is actually 0.3
# microns per pixel. This was discussed with the original data providers, who
# agree with my calculations
RESAMPLED_SPACING = 0.0003


def convert_segmentation(
    arrays: list[np.ndarray],
    coords: list[tuple[int, int]],
    source_image: hd.Image,
    segmentation_type: str,
    include_lut: bool,
    container_id: str,
    transfer_syntax_uid: str = ExplicitVRLittleEndian,
    crop_total_pixel_matrix: bool = False,
    bootstrapped: bool = False,
    inner_offset_coords: list[tuple[int, int]] | None = None,
) -> list[tuple[hd.seg.Segmentation, hd.seg.Segmentation, hd.seg.Segmentation]]:
    """Convert input array to DICOM Segmentations.

    Parameters
    ----------
    arrays: list[np.ndarray]
        List of segmentation arrays (height, width, 3), for a single slide,
        read from the PNG files. Each array is the segmentation of a single
        patch drawn from the slide. Each of the three channels encodes a
        different "set" of segmentations. First channel is regions, second is
        nuclei, third is boundaries.
    coords: list[tuple[int, int]]
        List of coordinates associated with each array. Coordinates are in
        pixel units of the source image, and give the (top, left) coordinates
        of the patch that was extracted before it was rescaled and then
        annotated.
    source_image: highdicom.Image
        Dataset (potentially without pixel data) of the whole slide image from
        which the segmentation arrays were derived.
    segmentation_type: str
        Segmentation type to use for the new DICOM segmentation files.
    include_lut: bool
        Whether to include a color palette lookup table for in the new
        segmentation (ignored unless ``segmentation_type`` is "LABELMAP").
    container_id: str
        Container Identifier for the source slide.
    transfer_syntax_uid: str
        Transfer syntax UID for the new segmentations.
    crop_total_pixel_matrix: bool
        If True, limit the size of the total pixel matrix of the output
        segmentation to the area covered by the data. If False, the total pixel
        matrix is defined to cover the same physical area covered by the source
        image. Note this does not affect the frames that are stored, just how
        the total pixel matrix rows and columns are defined, and the positional
        information of the frames.
    bootstrapped: bool
        Whether this is a "bootstrapped" (model-generated), as opposed to a
        manual, segmentation.
    inner_offset_coords: list[tuple[int, int]] | None
        Offset coordinates of the segmentation mask within the larger patch.
        These coordinates are given in pixel units of the resampled patch.
        Format is (top, left). If None, no offset coordinates are used.

    Returns
    -------
    list[tuple[highdicom.seg.Segmentation, highdicom.seg.Segmentation, highdicom.seg.Segmentation]]:
        List of output DICOM segmentations. Each item of the list corresponds
        to one of the regions in the slide, and contains a tuple of three
        Segmentations for (tissues, nuclei, boundaries).

    """
    if bootstrapped:
        algorithm_type = hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC
        region_desc = metadata_config.region_bootstrapped_series_description
        nuclei_desc = metadata_config.nuclei_bootstrapped_series_description
        border_desc = metadata_config.border_bootstrapped_series_description
        algorithm_identification = metadata_config.seg_algorithm_identification
        segment_label_suffix = " (Bootstrapped)"
        type_str = "bootstrapped"
        series_num_start = 200
    else:
        algorithm_type = hd.seg.SegmentAlgorithmTypeValues.MANUAL
        region_desc = metadata_config.region_series_description
        nuclei_desc = metadata_config.nuclei_series_description
        border_desc = metadata_config.border_series_description
        algorithm_identification = None
        segment_label_suffix = " (Manual)"
        type_str = "manual"
        series_num_start = 100

    top, left = coords[0]

    source_geom = cast(hd.VolumeGeometry, source_image.get_volume_geometry())

    if crop_total_pixel_matrix:
        min_top = min(top for (top, _) in coords)
        min_left = min(left for (_, left) in coords)
        seg_position = source_geom.map_indices_to_reference(
            np.array([[0, min_top, min_left]])
        )[0].tolist()
    else:
        seg_position = source_geom.position

    source_spacing = source_geom.pixel_spacing
    slice_spacing = source_geom.spacing_between_slices

    seg_geom = hd.VolumeGeometry.from_attributes(
        image_position=seg_position,
        image_orientation=source_geom.direction_cosines,
        pixel_spacing=[RESAMPLED_SPACING, RESAMPLED_SPACING],
        rows=100,  # not used
        columns=100,  # not used
        number_of_frames=1,
        spacing_between_slices=1.0,
        coordinate_system="SLIDE",
    )

    source_ind_to_seg_ind_transformer = hd.VolumeToVolumeTransformer(
        source_geom,
        seg_geom,
        round_output=False,
    )

    if slice_spacing is None:
        slice_spacing = 0.0

    pixel_measures = hd.PixelMeasuresSequence(
        pixel_spacing=[RESAMPLED_SPACING, RESAMPLED_SPACING],
        spacing_between_slices=slice_spacing,
        slice_thickness=0.0,
    )

    segs = []

    outer_patch_source_pix_indices = np.array(
        [[0, top, left] for (top, left) in coords]
    )

    outer_patch_seg_pix_indices = source_ind_to_seg_ind_transformer(outer_patch_source_pix_indices)
    outer_patch_seg_pix_indices = np.round(outer_patch_seg_pix_indices).astype(np.uint32)

    if inner_offset_coords is not None:
        inner_offsets_array = np.array(
            [[0, top, left] for (top, left) in inner_offset_coords]
        )
        inner_patch_seg_pix_indices = outer_patch_seg_pix_indices + inner_offsets_array
    else:
        inner_patch_seg_pix_indices = outer_patch_seg_pix_indices

    # PlanePositionSequence requires different order convention
    ref_coords = seg_geom.map_indices_to_reference(inner_patch_seg_pix_indices)
    pixel_matrix_positions = np.fliplr(inner_patch_seg_pix_indices[:, 1:]) + 1

    plane_positions = [
        hd.PlanePositionSequence(
            "SLIDE",
            pixel_matrix_position=pix,
            image_position=ref,
        ) for pix, ref in zip(pixel_matrix_positions, ref_coords)
    ]

    scale_y = source_spacing[0] / RESAMPLED_SPACING
    scale_x = source_spacing[1] / RESAMPLED_SPACING

    for patch_num, (patch_array, patch_pos) in enumerate(zip(arrays, plane_positions)):
        patch_segs = []

        # Loop over the three channels. Each creates its own Segmentation instance
        for c, desc, finding_codes in zip(
            range(3),
            [
                region_desc,
                nuclei_desc,
                border_desc,
            ],
            [
                metadata_config.region_finding_codes,
                metadata_config.nuclei_finding_codes,
                metadata_config.border_finding_codes,
            ],
        ):
            # Loop over all the segments within this input channel to create
            # segment descriptions
            segment_descriptions = [
                hd.seg.SegmentDescription(
                    segment_number=number,
                    segment_label=label + segment_label_suffix,
                    segmented_property_category=cat_code,
                    segmented_property_type=prop_code,
                    algorithm_type=algorithm_type,
                    tracking_id=f"{container_id}-{type_str}-ROI{patch_num}-{label}",
                    tracking_uid=hd.UID(),
                    algorithm_identification=algorithm_identification,
                ) for (number, (label, (prop_code, cat_code))) in enumerate(
                    finding_codes.items(),
                    start=1
                )
            ]

            channel_array = patch_array[:, :, c]

            seg = hd.seg.Segmentation(
                source_images=[source_image],
                pixel_array=channel_array,
                segment_descriptions=segment_descriptions,
                series_instance_uid=hd.UID(),
                series_number=series_num_start + c,
                sop_instance_uid=hd.UID(),
                series_description=desc,
                instance_number=1,
                segmentation_type=segmentation_type,
                manufacturer=metadata_config.seg_manufacturer,
                manufacturer_model_name=metadata_config.seg_manufacturer_model_name,
                software_versions=metadata_config.software_versions,
                device_serial_number=metadata_config.device_serial_number,
                transfer_syntax_uid=transfer_syntax_uid,
                dimension_organization_type="TILED_SPARSE",
                plane_positions=[patch_pos],
                pixel_measures=pixel_measures,
            )

            if not crop_total_pixel_matrix:
                # Post-process to set the total pixel matrix of the segmentation to
                # the same physical size as that of the source image. TODO
                # incorporate this option into highdicom
                seg.TotalPixelMatrixRows = int(
                    np.round(source_image.TotalPixelMatrixRows * scale_y)
                )
                seg.TotalPixelMatrixColumns = int(
                    np.round(source_image.TotalPixelMatrixColumns * scale_x)
                )

            patch_segs.append(seg)

        segs.append(tuple(patch_segs))

    return segs


def convert_annotation(
    dataframes: list[pd.DataFrame],
    coords: list[tuple[int, int]],
    source_image: hd.Image,
    bootstrapped: bool,
    inner_offset_coords: list[tuple[int, int]] | None = None,
) -> list[hd.ann.MicroscopyBulkSimpleAnnotations]:
    """Convert vector annotations from a dataframe to DICOM MBSA.

    Parameters
    ----------
    dataframes: list[pandas.DataFrame]
        List of dataframes (read from the original CSVs) containing annotation
        coordinates. Multiple dataframes may be passed, all relating to the
        same source image.
    coords: list[tuple[int, int]
        For each dataframe, a tuple of coordinates giving the coodinates of the
        patch the annotations apply to. Coordinates are in pixel units of the
        source image, and give the (top, left) coordinates of
        the patch that was extracted before it was rescaled and then annotated.
    source_image: highdicom.Image
        Source image that the annotations apply to (pixel data not required).
    bootstrapped: bool
        Whether thses annotations are bootstrapped (created by a model trained
        on manual annotations), or the original manual (False) annotations.
    inner_offset_coords: list[tuple[int, int]] | None
        Offset coordinates of the annotations within the larger patch.
        These coordinates are given in pixel units of the resampled patch.
        Format is (top, left). If None, no offset coordinates are used.
        
    """
    source_geom = cast(hd.VolumeGeometry, source_image.get_volume_geometry())
    source_spacing = source_geom.pixel_spacing
    scale_y = source_spacing[0] / RESAMPLED_SPACING
    scale_x = source_spacing[1] / RESAMPLED_SPACING

    ann_objects = []

    if bootstrapped:
        algorithm_type = hd.ann.AnnotationGroupGenerationTypeValues.AUTOMATIC
        series_description = metadata_config.ann_series_description_boostrapped
        algorithm_identification = metadata_config.ann_algorithm_identification
        series_num_start = 400
    else:
        algorithm_type = hd.ann.AnnotationGroupGenerationTypeValues.MANUAL
        series_description = metadata_config.ann_series_description_manual
        algorithm_identification = None
        series_num_start = 300

    if inner_offset_coords is None:
        inner_offset_coords = [(0, 0)] * len(dataframes)

    for roi_num, (df, patch_offset, inner_offset) in enumerate(
        zip(dataframes, coords, inner_offset_coords)
    ):
        ann_groups = []
        top, left = patch_offset
        inner_top, inner_left = inner_offset

        for n, ((label, ann_type), sub_df) in enumerate(df.groupby(['group', 'type']), 1):

            graphic_type = GRAPHIC_TYPE_MAPPING[ann_type]

            all_graphic_data = []
            for _, row in sub_df.iterrows():

                x = np.array([int(s) for s in row['coords_x'].split(',')])
                y = np.array([int(s) for s in row['coords_y'].split(',')])

                # To find WSI coordinate, add the corner offset and apply the
                # scaling factor
                y = top + (y + inner_top) / scale_y
                x = left + (x + inner_left) / scale_x

                graphic_data = np.stack([x, y]).T

                # graphic data must not duplicate the final point
                # (implicitly closed). Remove duplicates from the end to ensure
                # this. Sometimes there are many duplicates
                while np.array_equal(
                    graphic_data[0],
                    graphic_data[-1]
                ) and len(graphic_data) > 1:
                    graphic_data = graphic_data[:-1]

                all_graphic_data.append(graphic_data)

            property_category, property_type = metadata_config.csv_finding_codes[label]

            ann_groups.append(
                hd.ann.AnnotationGroup(
                    number=n,
                    uid=hd.UID(),
                    label=f"ROI {roi_num + 1}: {label} ({ann_type})",
                    annotated_property_type=property_type,
                    annotated_property_category=property_category,
                    algorithm_type=algorithm_type,
                    graphic_type=graphic_type,
                    graphic_data=all_graphic_data,
                    display_color=metadata_config.ann_color_mapping[label],
                    algorithm_identification=algorithm_identification,
                )
            )

        ann_objects.append(
            hd.ann.MicroscopyBulkSimpleAnnotations(
                source_images=[source_image],
                annotation_coordinate_type=hd.ann.AnnotationCoordinateTypeValues.SCOORD,
                annotation_groups=ann_groups,
                series_instance_uid=hd.UID(),
                series_number=series_num_start + roi_num,
                sop_instance_uid=hd.UID(),
                instance_number=1,
                manufacturer=metadata_config.ann_manufacturer,
                manufacturer_model_name=metadata_config.ann_manufacturer_model_name,
                software_versions=metadata_config.software_versions,
                device_serial_number=metadata_config.device_serial_number,
                content_description=metadata_config.ann_content_description,
                series_description=series_description
            )
        )

    return ann_objects
