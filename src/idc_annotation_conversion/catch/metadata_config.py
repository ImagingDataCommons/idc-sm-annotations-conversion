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
    "Bone": codes.SCT.Bone,
    "Cartilage": codes.SCT.Cartilage,
    "Dermis": CodedConcept('53534000', 'SCT', 'Dermis'),
    "Epidermis": CodedConcept('55988001', 'SCT', 'Epidermis'),
    "Subcutis": CodedConcept('71966008', 'SCT', 'Subcutis'),
    "Inflamm/Necrosis": codes.SCT.Necrosis,
    "Melanoma": codes.SCT.MalignantMelanoma,
    "Plasmacytoma": codes.SCT.Plasmacytoma,
    "Mast Cell Tumor": CodedConcept('89796001', 'SCT', 'Mastocytoma'),
    "PNST": CodedConcept('19897006', 'SCT', 'Malignant Peripheral Nerve Sheath Tumor'),
    "SCC": CodedConcept('1162767002', 'SCT', 'Squamous Cell Carcinoma'),
    "Trichoblastoma": CodedConcept("878881002", "SCT", "Trichoblastoma"),
    "Histiocytoma": CodedConcept('302843004', 'SCT', 'Histiocytoma'),
}


category_codes = {
    "Bone": codes.SCT.AnatomicStructure,
    "Cartilage": codes.SCT.AnatomicStructure,
    "Dermis": codes.SCT.AnatomicStructure,
    "Epidermis": codes.SCT.AnatomicStructure,
    "Subcutis": codes.SCT.AnatomicStructure,
    "Inflamm/Necrosis": codes.SCT.MorphologicallyAbnormalStructure,
    "Melanoma": codes.SCT.MorphologicallyAbnormalStructure,
    "Plasmacytoma": codes.SCT.MorphologicallyAbnormalStructure,
    "Mast Cell Tumor": codes.SCT.MorphologicallyAbnormalStructure,
    "PNST": codes.SCT.MorphologicallyAbnormalStructure,
    "SCC": codes.SCT.MorphologicallyAbnormalStructure,
    "Trichoblastoma": codes.SCT.MorphologicallyAbnormalStructure,
    "Histiocytoma": codes.SCT.MorphologicallyAbnormalStructure,
}
