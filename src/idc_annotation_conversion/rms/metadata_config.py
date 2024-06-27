"""Metadata used for RMS conversions."""
import highdicom as hd
import pydicom
from pydicom.sr.codedict import codes

from idc_annotation_conversion.git_utils import (
    get_git_remote_url,
    get_git_commit_hash,
)


# Dictionary mapping text label found in the XML annotations to tuple of
# (finding_type, finding_category) codes to encode that finding
finding_codes = {
    "STROMA": (
        hd.sr.CodedConcept(
            meaning="Connective tissue",
            value="181769001",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Body substance",
            value="91720002",
            scheme_designator="SCT",
        ),
    ),
    "ARMS": (
        hd.sr.CodedConcept(
            meaning="Alveolar rhabdomyosarcoma",
            value="63449009",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Morphologic abnormality",
            value="49755003",
            scheme_designator="SCT",
        ),
    ),
    "ERMS": (
        hd.sr.CodedConcept(
            meaning="Embryonal rhabdomyosarcoma",
            value="14269005",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Morphologic abnormality",
            value="49755003",
            scheme_designator="SCT",
        ),
    ),
    "NECROSIS": (
        hd.sr.CodedConcept(
            meaning="Necrosis",
            value="6574001",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Morphologic abnormality",
            value="49755003",
            scheme_designator="SCT",
        ),
    ),
}

# Metadata shared between SRs and SEGs
software_versions = get_git_remote_url(simplify=True)
device_serial_number = get_git_commit_hash()
institution_name = None
institutional_department_name = None

# Metadata Specific to XML hand anntations
observer_person_context = hd.sr.ObserverContext(
    observer_type=codes.DCM.Person,
    observer_identifying_attributes=hd.sr.PersonObserverIdentifyingAttributes(
        name=pydicom.valuerep.PersonName.from_named_components(
            family_name="Anonymized",
        )
    )
)
sr_series_description = "Manual region annotations"
sr_manufacturer = "Leica Biosystems"
sr_manufacturer_model_name = "Aperio ImageScope converted with highdicom"
title = codes.DCM.ImagingMeasurementReport
procedure_reported = hd.sr.CodedConcept(
    meaning="Light microscopy",
    value="104157003",
    scheme_designator="SCT",
)
is_complete = True
is_final = True
is_verified = False

# Metadata Specific to model output segmentations
segmentation_channel_order = [
    "NECROSIS",
    "STROMA",
    "ARMS",
    "ERMS",
]
algorithm_identification = hd.AlgorithmIdentificationSequence(
    name="Rhabdomyosarcoma Pathology CNN AI Segmentation Model",
    version="V1.0",
    source="Frederick National Lab for Cancer Research",
    family=codes.cid7162.ArtificialIntelligence,
)
segmentation_series_description = "AI Model Tissue Segmentations"
seg_manufacturer = "NCI/FNLCR"
seg_manufacturer_model_name = "FNLCR_IVG_RMS_iou_0.7343_0.7175_epoch_60"
