"""Metadata for the Catch conversions."""
from highdicom.sr import CodedConcept
from highdicom.base_content import ContributingEquipment
import pydicom
from pydicom.sr.codedict import codes
from pydicom.valuerep import PersonName
from idc_annotation_conversion.git_utils import (
    get_git_remote_url,
    get_git_commit_hash,
)

manufacturer = "Friedrich-Alexander-Universität Erlangen-Nürnberg converted by IDC"
manufacturer_model_name = "Annotations"
series_description = "Manual Region Annotations of Tumor and Tissue"
content_creator_name = PersonName.from_named_components(
    family_name="Fragoso",
    given_name="Marco"
)
series_number = 201
content_label = "ANNOTATIONS"
content_description = None
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

clinical_trial_ids_item = pydicom.Dataset()
clinical_trial_ids_item.IssuerOfClinicalTrialProtocolID = "DOI"
clinical_trial_ids_item.ClinicalTrialProtocolID = "doi:10.5281/zenodo.19488392"

finding_codes = {
    "Bone": codes.SCT.Bone,
    "Cartilage": codes.SCT.Cartilage,
    "Dermis": CodedConcept('53534000', 'SCT', 'Dermis'),
    "Epidermis": CodedConcept('55988001', 'SCT', 'Epidermis'),
    "Subcutis": CodedConcept('71966008', 'SCT', 'Subcutis'),
    "Inflamm/Necrosis": CodedConcept('316010', '99MP', 'Inflammation and/or necrosis'),
    "Melanoma": codes.SCT.MalignantMelanoma,
    "Plasmacytoma": codes.SCT.Plasmacytoma,
    "Mast Cell Tumor": CodedConcept('89796001', 'SCT', 'Mastocytoma'),
    "PNST": CodedConcept('19897006', 'SCT', 'Malignant Peripheral Nerve Sheath Tumor'),
    "SCC": CodedConcept('1162767002', 'SCT', 'Squamous Cell Carcinoma'),
    "Trichoblastoma": CodedConcept("878881002", "SCT", "Trichoblastoma"),
    "Histiocytoma": CodedConcept('302843004', 'SCT', 'Histiocytoma'),
}


category_codes = {
    "Bone": codes.SCT.AnatomicalStructure,
    "Cartilage": codes.SCT.AnatomicalStructure,
    "Dermis": codes.SCT.AnatomicalStructure,
    "Epidermis": codes.SCT.AnatomicalStructure,
    "Subcutis": codes.SCT.AnatomicalStructure,
    "Inflamm/Necrosis": codes.SCT.MorphologicallyAbnormalStructure,
    "Melanoma": codes.SCT.MorphologicallyAbnormalStructure,
    "Plasmacytoma": codes.SCT.MorphologicallyAbnormalStructure,
    "Mast Cell Tumor": codes.SCT.MorphologicallyAbnormalStructure,
    "PNST": codes.SCT.MorphologicallyAbnormalStructure,
    "SCC": codes.SCT.MorphologicallyAbnormalStructure,
    "Trichoblastoma": codes.SCT.MorphologicallyAbnormalStructure,
    "Histiocytoma": codes.SCT.MorphologicallyAbnormalStructure,
}
