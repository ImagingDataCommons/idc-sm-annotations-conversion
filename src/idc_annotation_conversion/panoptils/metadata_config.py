import highdicom as hd
from idc_annotation_conversion.git_utils import (
    get_git_remote_url,
    get_git_commit_hash,
)

region_series_description = "PanopTILs Manual Region Segmentations"
nuclei_series_description = "PanopTILs Manual Nuclei Segmentations"
border_series_description = "PanopTILs Manual Border Segmentations"

region_finding_codes = {
    "Cancerous epithelium": (
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
    "Stroma": (
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
    "TILs": (
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
    "Normal epithelium": (
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
    "Junk/Debris": (
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
    "Blood": (
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
    "Other": (
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
    "Whitespace/Empty": (
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

nuclei_finding_codes = {
    "Cancer nucleus": (
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
    "Stromal nucleus": (
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
    "Large stromal nucleus": (
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
    "Lymphocyte nucleus": (
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
    "Plasma cell/large TIL nucleus": (
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
    "Normal epithelial nucleus": (
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
    "Other nucleus": (
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
    "Unknown/ambiguous nucleus": (
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
    "Background (non-nuclear material)": (
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


border_finding_codes = {
    "Border": (
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
}


seg_manufacturer = "Manufacturer"
seg_manufacturer_model_name = "Model Name"
software_versions = get_git_remote_url(simplify=True)
device_serial_number = get_git_commit_hash()
