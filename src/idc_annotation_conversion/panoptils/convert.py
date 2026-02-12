import highdicom as hd
import numpy as np
import pydicom
from typing import cast

from idc_annotation_conversion.panoptils import metadata_config


def convert_segmentation(
    arrays: list[np.ndarray],
    coords: list[tuple[int, int, int, int]],
    source_image: hd.Image,
    segmentation_type: str,
    include_lut: bool,
    container_id: str,
) -> tuple[hd.seg.Segmentation, hd.seg.Segmentation, hd.seg.Segmentation]:
    """Convert input array to DICOM Segmentations.

    Parameters
    ----------
    arrays: list[np.ndarray]
        List of segmentation arrays (height, width, 3), for a single slide,
        read from the PNG files. Each array is the segmentation of a single
        patch drawn from the slide. Each of the three channels encodes a
        different "set" of segmentations. First channel is regions, second is
        nuclei, third is boundaries.
    coords: list[tuple[int, int, int, int]]
        List of coordinates associated with each array. Item i of the list
        gives the coordinates of the patch for item i of the ``arrays`` list.
        Each item is a tuple containing the (top, bottom, left, right)
        coordinates for the patch.
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

    Returns
    -------
    highdicom.seg.Segmentation:
        Output segmentation containing the tissues.
    highdicom.seg.Segmentation:
        Output segmentation containing the nuclei.
    highdicom.seg.Segmentation:
        Output segmentation containing the boundaries.

    """
    array = np.stack(arrays)

    t, b, l, r = coords[0]
    patch_size = b - t

    min_top = min(t for t, _, _, _ in coords)
    min_left = min(l for _, _, l, _ in coords)

    source_geom = cast(hd.VolumeGeometry, source_image.get_volume_geometry())

    seg_position = source_geom.map_indices_to_reference(
        np.array([[0, min_top, min_left]])
    )[0]

    # Create a transformer that maps total pixel matrix indices to frame of
    # reference coordinate locations
    source_spacing = source_geom.pixel_spacing
    slice_spacing = source_geom.spacing_between_slices

    # The segmentation arrays are resampled at a different pixel spacing from
    # the source image. Find the spacing of pixels in the segmentation array
    # from the scaling factor between the two array sizes
    seg_spacing = [
        source_spacing[0] * (patch_size / array.shape[1]),
        source_spacing[1] * (patch_size / array.shape[2]),
    ]

    seg_geom = hd.VolumeGeometry.from_attributes(
        image_position=seg_position.tolist(),
        image_orientation=source_geom.direction_cosines,
        pixel_spacing=seg_spacing,
        rows=100,  # not used
        columns=100,  # not used
        number_of_frames=1,
        spacing_between_slices=1.0,
        coordinate_system="SLIDE",
    )

    source_ind_to_seg_ind_transformer = hd.VolumeToVolumeTransformer(
        source_geom,
        seg_geom,
        round_output=True
    )

    source_pix_indices_3d = np.array([[0, t, l] for (t, _, l, _) in coords])
    seg_pix_indices = source_ind_to_seg_ind_transformer(source_pix_indices_3d)

    # PlanePositionSequence requires different order convention
    ref_coords = seg_geom.map_indices_to_reference(seg_pix_indices)
    pixel_matrix_positions = np.fliplr(seg_pix_indices[:, 1:]) + 1

    plane_positions = [
        hd.PlanePositionSequence(
            "SLIDE",
            pixel_matrix_position=pix,
            image_position=ref,
        ) for pix, ref in zip(pixel_matrix_positions, ref_coords)
    ]

    if slice_spacing is None:
        slice_spacing = 0.0

    pixel_measures = hd.PixelMeasuresSequence(
        pixel_spacing=seg_spacing,
        spacing_between_slices=slice_spacing,
        slice_thickness=0.0,
    )

    segs = []

    # Loop over the three channels. Each creates its own Segmentation instance
    for c, desc, finding_codes in zip(
        range(3),
        [
            metadata_config.region_series_description,
            metadata_config.nuclei_series_description,
            metadata_config.border_series_description,
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
                segment_label=label,
                segmented_property_category=cat_code,
                segmented_property_type=prop_code,
                algorithm_type=hd.seg.SegmentAlgorithmTypeValues.MANUAL,
                tracking_id=f"{container_id}-{label}",
                tracking_uid=hd.UID(),
            ) for (number, (label, (prop_code, cat_code))) in enumerate(
                finding_codes.items(),
                start=1
            )
        ]

        segs.append(
            hd.seg.Segmentation(
                source_images=[source_image],
                pixel_array=array[:, :, :, c],
                segment_descriptions=segment_descriptions,
                series_instance_uid=hd.UID(),
                series_number=20 + c,
                sop_instance_uid=hd.UID(),
                series_description=desc,
                instance_number=1,
                segmentation_type=segmentation_type,
                manufacturer=metadata_config.seg_manufacturer,
                manufacturer_model_name=metadata_config.seg_manufacturer_model_name,
                software_versions=metadata_config.software_versions,
                device_serial_number=metadata_config.device_serial_number,
                #transfer_syntax_uid=transfer_syntax_uid,
                dimension_organization_type="TILED_SPARSE",
                plane_positions=plane_positions,
                pixel_measures=pixel_measures,
            )
        )

    return tuple(segs)
