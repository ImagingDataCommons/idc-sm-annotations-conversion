import highdicom as hd
import numpy as np
import pydicom

from idc_annotation_conversion.panoptils import metadata_config


def convert_segmentation(
    arrays: list[np.ndarray],
    coords: list[tuple[int, int, int, int]],
    source_image: pydicom.Dataset,
    segmentation_type: str,
    include_lut: bool,
    container_id: str,
) -> tuple[hd.seg.Segmentation, hd.seg.Segmentation, hd.seg.Segmentation]:

    array = np.stack(arrays)

    # TODO fix this
    plane_positions = [
        hd.PlanePositionSequence(
            "SLIDE",
            image_position=(l, t, 0.0),
            pixel_matrix_position=(l, t),
        ) for (t, b, l, r) in coords
    ]

    pixel_measures = hd.PixelMeasuresSequence(
        pixel_spacing=[0.00025, 0.00025],
        spacing_between_slices=0.0,
        slice_thickness=0.0,
    )

    segs = []

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
        segment_descriptions = []
        for (number, (label, (prop_code, cat_code))) in enumerate(
            finding_codes.items(),
            start=1
        ):
            segment_descriptions.append(
                hd.seg.SegmentDescription(
                    segment_number=number,
                    segment_label=label,
                    segmented_property_category=cat_code,
                    segmented_property_type=prop_code,
                    algorithm_type=hd.seg.SegmentAlgorithmTypeValues.MANUAL,
                    tracking_id=f"{container_id}-{label}",
                    tracking_uid=hd.UID(),
                )
            )

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
