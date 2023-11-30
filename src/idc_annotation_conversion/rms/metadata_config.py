"""Metadata used for RMS conversions."""
import highdicom as hd
import pydicom
from pydicom.sr.codedict import codes

from idc_annotation_conversion.git_utils import (
    get_git_remote_url,
    get_git_commit_hash,
)


# Dictionary mapping text label found in the XML annoations to tuple of
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

observer_person_context = hd.sr.ObserverContext(
    observer_type=codes.DCM.Person,
    observer_identifying_attributes=hd.sr.PersonObserverIdentifyingAttributes(
        name=pydicom.valuerep.PersonName.from_named_components(
            family_name="Brown",
            given_name="G",
            middle_name="Tom",
        )
    )
)
observation_context = hd.sr.ObservationContext(
    observer_person_context=observer_person_context,
)

series_description = "Manual region annotations"
manufacturer = "manufacturer"
manufacturer_model_name = "manufacturer_model_name"
software_versions = get_git_remote_url(simplify=True)
device_serial_number = get_git_commit_hash()
institution_name = "institution_name"
institutional_department_name = "institutional_department_name"
title = codes.DCM.ImagingMeasurementReport
procedure_reported = codes.SCT.ImagingProcedure
