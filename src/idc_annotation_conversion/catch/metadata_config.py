"""Metadata for the Catch conversions."""
from highdicom.sr import CodedConcept
from highdicom.base_content import ContributingEquipment
from pydicom.sr.codedict import codes
from idc_annotation_conversion.git_utils import (
    get_git_remote_url,
    get_git_commit_hash,
)

manufacturer = "Manufacturer"
manufacturer_model_name = "Model Name"
series_description = "Region Annotations"
series_number = 201
content_label = "REGIONS"
content_description = "Description"
software_versions = get_git_remote_url(simplify=True)
device_serial_number = get_git_commit_hash()

contributing_equipment = [
    ContributingEquipment(
        manufacturer="Friedrich-Alexander-Universität Erlangen-Nürnberg",
        manufacturer_model_name="SlideRunner",
        purpose_of_reference=codes.DCM.AcquisitionEquipment,
        contribution_description="Annotation tool used to acquire manual annotations."
    ),
]


finding_codes = {
    "Bone": CodedConcept('Bone', 'TMP', 'Bone'),
    "Cartilage": CodedConcept('Cartilage', 'TMP', 'Cartilage'),
    "Dermis": CodedConcept('Dermis', 'TMP', 'Dermis'),
    "Epidermis": CodedConcept('Epidermis', 'TMP', 'Epidermis'),
    "Subcutis": CodedConcept('Subcutis', 'TMP', 'Subcutis'),
    "Inflamm/Necrosis": CodedConcept('Necrosis', 'TMP', 'Necrosis'),
    "Melanoma": CodedConcept('Melanoma', 'TMP', 'Melanoma'),
    "Plasmacytoma": CodedConcept('Plasmacytoma', 'TMP', 'Plasmacytoma'),
    "Mast Cell Tumor": CodedConcept('Mast Cell Tumor', 'TMP', 'Mast Cell Tumor'),
    "PNST": CodedConcept('PNST', 'TMP', 'PNST'),
    "SCC": CodedConcept('SCC', 'TMP', 'SCC'),
    "Trichoblastoma": CodedConcept('Trichoblastoma', 'TMP', 'Trichoblastoma'),
    "Histiocytoma": CodedConcept('Histiocytoma', 'TMP', 'Histiocytoma'),
}


category_codes = {
    "Bone": CodedConcept('Bone', 'TMP', 'Bone'),
    "Cartilage": CodedConcept('Cartilage', 'TMP', 'Cartilage'),
    "Dermis": CodedConcept('Dermis', 'TMP', 'Dermis'),
    "Epidermis": CodedConcept('Epidermis', 'TMP', 'Epidermis'),
    "Subcutis": CodedConcept('Subcutis', 'TMP', 'Subcutis'),
    "Inflamm/Necrosis": CodedConcept('Necrosis', 'TMP', 'Necrosis'),
    "Melanoma": CodedConcept('Melanoma', 'TMP', 'Melanoma'),
    "Plasmacytoma": CodedConcept('Plasmacytoma', 'TMP', 'Plasmacytoma'),
    "Mast Cell Tumor": CodedConcept('Mast Cell Tumor', 'TMP', 'Mast Cell Tumor'),
    "PNST": CodedConcept('PNST', 'TMP', 'PNST'),
    "SCC": CodedConcept('SCC', 'TMP', 'SCC'),
    "Trichoblastoma": CodedConcept('Trichoblastoma', 'TMP', 'Trichoblastoma'),
    "Histiocytoma": CodedConcept('Histiocytoma', 'TMP', 'Histiocytoma'),
}
