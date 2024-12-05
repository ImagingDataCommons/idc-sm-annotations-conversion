"""Functions to convert single annotations into SRs"""
import logging
from time import time
from typing import List, Union
import xml

import numpy as np
import highdicom as hd
import pydicom
from pydicom import Dataset
from pydicom.sr.codedict import codes
from pydicom.uid import (
    JPEGLSLossless,
    ExplicitVRLittleEndian,
    JPEG2000Lossless,
)
from rasterio.features import rasterize
from shapely.geometry.polygon import Polygon

from idc_annotation_conversion.rms import metadata_config


def convert_xml_annotation(
    xml_annotation: xml.etree.ElementTree.Element,
    source_images: List[Dataset],
    use_scoord3d: bool = True,
    include_measurements: bool = False,
    create_segmentation: bool = True,
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str] = hd.seg.SegmentationTypeValues.BINARY,
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str] = hd.DimensionOrganizationTypeValues.TILED_FULL,
    workers: int = 0,
    include_lut: bool = False,
) -> hd.sr.Comprehensive3DSR:
    """Convert an ImageScope XML annotation to a DICOM SR.

    Parameters
    ----------
    xml_annotation: xml.etree.ElementTree.Element
        Pre-loaded root element of the annotation file's XML tree.
    source_images: List[pydicom.Dataset]
        List of dataset of the source images to which the annotation applies.
    use_scoord3d: bool
        Use SCOORD3D coordinates to store points.
    include_measurements: bool
        Include the measurements of length and area from the XML.
    create_segmentation: bool
        Whether to create a (raster) segmentation mask from the (vector)
        contour annotations.
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str], optional
        Segmentation type (BINARY or FRACTIONAL) for the Segmentation Image
        (if any).
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str], optional
        Dimension organization type of the output segmentations.
    workers: int
        Number of workers to use for frame compression in segmentations.
    include_lut: bool, optional
        Whether to include a LUT transformation in the segmentation to store as
        a PALETTE COLOR image. Ignored if segmentation_type is not LABELMAP.

    Returns
    -------
    highdicom.sr.Comprehensive3DSR:
        DICOM SR object encoding the annotation.

    """
    source_image = source_images[0]
    assert xml_annotation.tag == "Annotations"
    microns_per_pixel = float(xml_annotation.attrib["MicronsPerPixel"])

    origin_seq = source_image.TotalPixelMatrixOriginSequence[0]
    origin = (
        origin_seq.XOffsetInSlideCoordinateSystem,
        origin_seq.YOffsetInSlideCoordinateSystem,
        0.0
    )
    pixel_spacing_mm = microns_per_pixel / 1000.0
    transformer = hd.spatial.ImageToReferenceTransformer(
        image_position=origin,
        image_orientation=source_image.ImageOrientationSlide,
        pixel_spacing=(pixel_spacing_mm, pixel_spacing_mm),
    )
    container_id = source_images[0].ContainerIdentifier

    roi_groups = []
    polygons = []

    for annotation in xml_annotation:
        assert annotation.tag == "Annotation"

        regions = annotation[1]
        assert regions.tag == "Regions"

        for region in regions:
            if region.tag == "RegionAttributeHeaders":
                continue
            assert region.tag == "Region"
            vertices = region[1]
            assert vertices.tag == "Vertices"

            graphic_data = np.array(
                [
                    (float(v.attrib["X"]), float(v.attrib["Y"]))
                    for v in vertices
                ]
            )

            # Ensure polygon is closed (this seem to be inconsistent in the
            # source XML files)
            if not np.array_equal(graphic_data[0], graphic_data[-1]):
                graphic_data = np.vstack([graphic_data, graphic_data[0]])

            region_id = f"Region {region.attrib['Id']}: {region.attrib['Text']}"
            tracking_identifier = hd.sr.TrackingIdentifier(hd.UID(), region_id)
            finding_str = region.attrib["Text"].upper()
            finding_type, finding_category = metadata_config.finding_codes[finding_str]

            mask_value = metadata_config.segmentation_channel_order.index(finding_str) + 1
            polygons.append(
                (Polygon(graphic_data), mask_value)
            )

            if use_scoord3d:
                graphic_data_3d = transformer(graphic_data)
                image_region: Union[
                    hd.sr.ImageRegion,
                    hd.sr.ImageRegion3D
                ] = hd.sr.ImageRegion3D(
                    graphic_type=hd.sr.GraphicTypeValues3D.POLYGON,
                    graphic_data=graphic_data_3d,
                    frame_of_reference_uid=source_image.FrameOfReferenceUID,
                )
            else:
                image_region = hd.sr.ImageRegion(
                    graphic_type=hd.sr.GraphicTypeValues.POLYLINE,
                    graphic_data=graphic_data,
                    source_image=hd.sr.SourceImageForRegion.from_source_image(
                        source_images[0]
                    ),
                )

            if include_measurements:
                area_measurement = hd.sr.Measurement(
                    name=codes.SCT.Area,
                    value=float(region.attrib["AreaMicrons"]),
                    unit=codes.UCUM.SquareMicrometer,
                )
                length_measurement = hd.sr.Measurement(
                    name=codes.SCT.Length,
                    value=float(region.attrib["LengthMicrons"]),
                    unit=codes.UCUM.Micrometer,
                )
                measurements = [area_measurement, length_measurement]
            else:
                measurements = None

            roi = hd.sr.PlanarROIMeasurementsAndQualitativeEvaluations(
                tracking_identifier=tracking_identifier,
                referenced_region=image_region,
                finding_type=finding_type,
                finding_category=finding_category,
                measurements=measurements,
            )
            roi_groups.append(roi)

    subject_context = hd.sr.SubjectContext.from_image(
        source_image,
    )

    observation_context = hd.sr.ObservationContext(
        observer_person_context=metadata_config.observer_person_context,
        subject_context=subject_context,
    )

    measurement_report = hd.sr.MeasurementReport(
        observation_context=observation_context,
        procedure_reported=metadata_config.procedure_reported,
        imaging_measurements=roi_groups,
        title=metadata_config.title,
        referenced_images=source_images,
    )

    logging.info("Creating SR.")
    sr = hd.sr.Comprehensive3DSR(
        evidence=source_images,
        content=measurement_report,
        is_complete=metadata_config.is_complete,
        is_final=metadata_config.is_final,
        is_verified=metadata_config.is_verified,
        series_number=1,
        series_instance_uid=hd.UID(),
        sop_instance_uid=hd.UID(),
        instance_number=1,
        series_description=metadata_config.sr_series_description,
        manufacturer=metadata_config.sr_manufacturer,
        manufacturer_model_name=metadata_config.sr_manufacturer_model_name,
        software_versions=metadata_config.software_versions,
        device_serial_number=metadata_config.device_serial_number,
        institution_name=metadata_config.institution_name,
        institutional_department_name=metadata_config.institutional_department_name,
    )
    logging.info("SR created.")

    if create_segmentation:
        im_shape = (
            source_image.TotalPixelMatrixRows,
            source_image.TotalPixelMatrixColumns
        )

        segmentation_mask = rasterize(
            polygons,
            out_shape=im_shape,
            dtype=np.uint8
        )

        segment_descriptions = []
        for (number, label) in enumerate(
            metadata_config.segmentation_channel_order,
            start=1
        ):
            desc = hd.seg.SegmentDescription(
                segment_number=number,
                segment_label=label,
                segmented_property_category=metadata_config.finding_codes[label][1],
                segmented_property_type=metadata_config.finding_codes[label][0],
                algorithm_type=hd.seg.SegmentAlgorithmTypeValues.MANUAL,
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
            hd.seg.SegmentationTypeValues.FRACTIONAL: JPEGLSLossless,
            hd.seg.SegmentationTypeValues.LABELMAP: JPEGLSLossless,
        }[segmentation_type]

        omit_empty_frames = dimension_organization_type.value != "TILED_FULL"

        logging.info("Creating DICOM Segmentations")
        seg_start_time = time()
        segmentations = hd.seg.pyramid.create_segmentation_pyramid(
            source_images=source_images,
            pixel_arrays=[segmentation_mask],
            segmentation_type=segmentation_type,
            segment_descriptions=segment_descriptions,
            series_instance_uid=hd.UID(),
            series_number=25,
            manufacturer=metadata_config.sr_manufacturer,
            manufacturer_model_name=metadata_config.sr_manufacturer_model_name,
            software_versions=metadata_config.software_versions,
            device_serial_number=metadata_config.device_serial_number,
            transfer_syntax_uid=transfer_syntax_uid,
            dimension_organization_type=dimension_organization_type,
            omit_empty_frames=omit_empty_frames,
            workers=workers,
            series_description=metadata_config.manual_segmentation_series_description_by_type[segmentation_type],
            palette_color_lut_transformation=lut_transform,
        )
        seg_time = time() - seg_start_time
        logging.info(f"Created DICOM Segmentations in {seg_time:.1f}s.")
    else:
        segmentations = []


    return sr, segmentations


def convert_segmentation(
    source_images: List[pydicom.Dataset],
    segmentation_array: np.ndarray,
    create_pyramid: bool,
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str],
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str],
    workers: int = 0,
    include_lut: bool = False,
) -> List[hd.seg.Segmentation]:
    """Store segmentation masks as DICOM segmentations.

    Parameters
    ----------
    source_images: Sequence[pydicom.Dataset]
        List of pydicom datasets containing the metadata of the image (already
        converted to DICOM format). Note that the metadata of the image at full
        resolution should appear first in this list. These can be the full image
        datasets, but the PixelData attributes are not required.
    segmentation_array: np.ndarray
        Segmentation output (before thresholding). Should have shape (rows,
        columns, 5), where 5 is the number of classes. Values are in the range
        0 to 1.
    create_pyramid: bool, optional
        Whether to create a full pyramid of segmentations (rather than a single
        segmentation at the highest resolution).
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str], optional
        Segmentation type (BINARY or FRACTIONAL) for the Segmentation Image
        (if any).
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str], optional
        Dimension organization type of the output segmentations.
    workers: int
        Number of workers to use for frame compression.
    include_lut: bool, optional
        Whether to include a LUT transformation in the segmentation to store as
        a PALETTE COLOR image. Ignored if segmentation_type is not LABELMAP.


    Returns
    -------
    segmentations: list[hd.seg.Segmentation]
        DICOM segmentation image(s) encoding the original annotations

    """
    seg_start_time = time()

    container_id = source_images[0].ContainerIdentifier
    segment_descriptions = []
    for (number, label) in enumerate(
        metadata_config.segmentation_channel_order,
        start=1
    ):
        desc = hd.seg.SegmentDescription(
            segment_number=number,
            segment_label=label,
            segmented_property_category=metadata_config.finding_codes[label][1],
            segmented_property_type=metadata_config.finding_codes[label][0],
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
        hd.seg.SegmentationTypeValues.FRACTIONAL: JPEGLSLossless,
        hd.seg.SegmentationTypeValues.LABELMAP: JPEGLSLossless,
    }[segmentation_type]

    omit_empty_frames = dimension_organization_type.value != "TILED_FULL"

    if segmentation_type == hd.seg.SegmentationTypeValues.FRACTIONAL:
        # Add frame axis and remove background class
        mask = segmentation_array[None, :, :, 1:]
    else:
        mask = np.argmax(segmentation_array, axis=2).astype(np.uint8)

    if create_pyramid:
        logging.info("Creating DICOM Segmentations")
        seg_start_time = time()
        segmentations = hd.seg.pyramid.create_segmentation_pyramid(
            source_images=source_images,
            pixel_arrays=[mask],
            segmentation_type=segmentation_type,
            segment_descriptions=segment_descriptions,
            series_instance_uid=hd.UID(),
            series_number=20,
            manufacturer=metadata_config.seg_manufacturer,
            manufacturer_model_name=metadata_config.seg_manufacturer_model_name,
            software_versions=metadata_config.software_versions,
            device_serial_number=metadata_config.device_serial_number,
            transfer_syntax_uid=transfer_syntax_uid,
            dimension_organization_type=dimension_organization_type,
            omit_empty_frames=omit_empty_frames,
            workers=workers,
            series_description=metadata_config.segmentation_series_description_by_type[segmentation_type],
            palette_color_lut_transformation=lut_transform,
        )
        seg_time = time() - seg_start_time
        logging.info(f"Created DICOM Segmentations in {seg_time:.1f}s.")
    else:
        logging.info("Creating DICOM Segmentation")
        seg_start_time = time()
        segmentation = hd.seg.Segmentation(
            source_images=source_images[0:1],
            pixel_array=mask,
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
            tile_pixel_array=True,
            dimension_organization_type=dimension_organization_type,
            omit_empty_frames=omit_empty_frames,
            workers=workers,
            series_description=metadata_config.segmentation_series_description_by_type[segmentation_type],
            palette_color_lut_transformation=lut_transform,
        )
        segmentations = [segmentation]
        seg_time = time() - seg_start_time
        logging.info(f"Created DICOM Segmentation in {seg_time:.1f}s.")

    return segmentations
