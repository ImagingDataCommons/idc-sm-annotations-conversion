"""Functions to convert single annotations into SRs"""
import xml

import highdicom as hd
import numpy as np
from pydicom import Dataset
from pydicom.sr.codedict import codes


FINDINGS = {
    "STROMA": (
        hd.sr.CodedConcept(meaning="Stroma", value="stroma", scheme="custom"),
        hd.sr.CodedConcept(meaning="Stroma", value="stroma", scheme="custom"),
    ),
    "ARMS": (
        hd.sr.CodedConcept(meaning="ARMS", value="arms", scheme="custom"),
        hd.sr.CodedConcept(meaning="ARMS", value="arms", scheme="custom"),
    ),
    "ERMS": (
        hd.sr.CodedConcept(meaning="ERMS", value="erms", scheme="custom"),
        hd.sr.CodedConcept(meaning="ERMS", value="erms", scheme="custom"),
    ),
    "NECROSIS": (
        hd.sr.CodedConcept(meaning="Necrosis", value="necrosis", scheme="custom"),
        hd.sr.CodedConcept(meaning="Necrosis", value="necrosis", scheme="custom"),
    ),
}


def convert_xml_annotation(
    xml_annotation: xml.etree.ElementTree.Element,
    source_image: Dataset,
) -> hd.sr.ComprehensiveSR:
    """Convert an ImageScope XML annotation to a DICOM SR.

    Parameters
    ----------
    xml_annotation: xml.etree.ElementTree.Element
        Pre-loaded root element of the annotation file's XML tree.
    source_image: pydicom.Dataset
        Dataset of the source image to which the annotation applies.

    Returns
    -------
    highdicom.sr.ComprehensiveSR:
        DICOM SR object encoding the annotation.

    """
    assert xml_annotation.tag == "Annotations"
    microns_per_pixel = xml_root.attrib["MicronsPerPixel"]

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
                    (float(v.attrib["X"]), float(v.attrib("Y")))
                    for v in vertices
                ]
            )
            # TODO is Id unique or only within an annotation??
            tracking_identifier = hd.sr.TrackingIdentifier(hd.UID(), region.attrib["Id"])
            finding_str = region.attrib["Text"]
            finding_type, finding_category = FINDINGS[finding_str]
            image_region = hd.sr.ImageRegion(
                graphic_type=hd.sr.GraphicTypeValues.POLYLINE,
                graphic_data=graphic_data,
                source_image=hd.sr.SourceImageForRegion.from_source_image(source_image),
            )
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
            roi = hd.sr.PlanarROIMeasurementsAndQualitativeEvaluations(
                tracking_identifier=tracking_identifier,
                referenced_region=image_region,
                finding_type=finding_type,
                finding_category=finding_category,
                measurements=[area_measurement, length_measurement],
            )
            roi_groups.append(roi)
