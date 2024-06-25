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
from pydicom.uid import JPEGLSLossless, ExplicitVRLittleEndian

from idc_annotation_conversion.rms import metadata_config


def convert_xml_annotation(
    xml_annotation: xml.etree.ElementTree.Element,
    source_images: List[Dataset],
    use_scoord3d: bool = True,
    include_measurements: bool = False,
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

    roi_groups = []

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

    # TODO fold this logic into highdicom library
    subject_context_specimen = hd.sr.SubjectContextSpecimen(
        uid=source_image.SpecimenDescriptionSequence[0].SpecimenUID,
        identifier=source_image.SpecimenDescriptionSequence[0].SpecimenIdentifier,
        container_identifier=source_image.ContainerIdentifier,
    )
    subject_context = hd.sr.SubjectContext(
        subject_class=codes.DCM.Specimen,
        subject_class_specific_context=subject_context_specimen,
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
        series_description=metadata_config.series_description,
        manufacturer=metadata_config.manufacturer,
        manufacturer_model_name=metadata_config.manufacturer_model_name,
        software_versions=metadata_config.software_versions,
        device_serial_number=metadata_config.device_serial_number,
        institution_name=metadata_config.institution_name,
        institutional_department_name=metadata_config.institutional_department_name,
    )

    return sr


def convert_segmentation(
    source_images: List[pydicom.Dataset],
    segmentation_array: np.ndarray,
    create_pyramid: bool,
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str],
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str],
    workers: int = 0,
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

    Returns
    -------
    segmentations: list[hd.seg.Segmentation]
        DICOM segmentation image(s) encoding the original annotations

    """
    seg_start_time = time()

    segment_descriptions = []
    for number, (label, (prop_code, cat_code)) in enumerate(
        metadata_config.finding_codes.items(),
        start=1
    ):
        desc = hd.seg.SegmentDescription(
            segment_number=number,
            segment_label=label,
            segmented_property_category=cat_code,
            segmented_property_type=prop_code,
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
            algorithm_identification=metadata_config.algorithm_identification
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
        hd.seg.SegmentationTypeValues.FRACTIONAL: JPEGLSLossless,
    }[segmentation_type]

    omit_empty_frames = dimension_organization_type.value != "TILED_FULL"

    if segmentation_type == hd.seg.SegmentationTypeValues.FRACTIONAL:
        # Add frame axis and remove background class
        mask = segmentation_array[None, :, :, 1:]
    elif segmentation_type == hd.seg.SegmentationTypeValues.BINARY:
        print("Before", segmentation_array.shape)
        mask = np.argmax(segmentation_array, axis=2).astype(np.uint8)
        print("after", mask.shape)

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
            manufacturer=metadata_config.manufacturer,
            manufacturer_model_name=metadata_config.manufacturer_model_name,
            software_versions=metadata_config.software_versions,
            device_serial_number=metadata_config.device_serial_number,
            transfer_syntax_uid=transfer_syntax_uid,
            max_fractional_value=1,
            dimension_organization_type=dimension_organization_type,
            omit_empty_frames=omit_empty_frames,
            workers=workers,
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
            manufacturer=metadata_config.manufacturer,
            manufacturer_model_name=metadata_config.manufacturer_model_name,
            software_versions=metadata_config.software_versions,
            device_serial_number=metadata_config.device_serial_number,
            transfer_syntax_uid=transfer_syntax_uid,
            max_fractional_value=1,
            tile_pixel_array=True,
            dimension_organization_type=dimension_organization_type,
            omit_empty_frames=omit_empty_frames,
            workers=workers,
        )
        segmentations = [segmentation]
        seg_time = time() - seg_start_time
        logging.info(f"Created DICOM Segmentation in {seg_time:.1f}s.")

    return segmentations
