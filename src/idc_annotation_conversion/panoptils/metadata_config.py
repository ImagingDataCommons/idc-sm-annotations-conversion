from pydicom.sr.codedict import codes
import highdicom as hd
from idc_annotation_conversion.git_utils import (
    get_git_remote_url,
    get_git_commit_hash,
)


# Metadata shared by annotations and segmentations
manufacturer = "Northwestern University converted by IDC"
manufacturer_model_name = "PanopTILs"
algorithm_identification = hd.AlgorithmIdentificationSequence(
    name="MuTILs",
    family=codes.DCM.ArtificialIntelligence,
    version="1.0",
    source="Northwestern University",
)  # only for bootstrapped
software_versions = get_git_remote_url(simplify=True)
device_serial_number = get_git_commit_hash()

# Segmentation-specific metadata
region_series_description = "PanopTILs Manual Region Segmentations"
nuclei_series_description = "PanopTILs Manual Nuclei Segmentations"
border_series_description = "PanopTILs Manual Border Segmentations"
region_bootstrapped_series_description = "PanopTILs Boostrapped Region Segmentations"
nuclei_bootstrapped_series_description = "PanopTILs Boostrapped Nuclei Segmentations"
border_bootstrapped_series_description = "PanopTILs Boostrapped Border Segmentations"

# Annotation-specific metadata
ann_content_description = "Cell type annotations"
ann_series_description_manual = "PanopTILs Manual Cell Type Annotations"
ann_series_description_boostrapped = "PanopTILs Boostrapped Cell Type Annotations"

# From fig 3 in the paper
ann_color_mapping = {
    'ActiveStromalCellNOS': hd.color.CIELabColor.from_string("mediumorchid"),
    'ActiveTILsCell': hd.color.CIELabColor.from_string("royalblue"),
    'CancerEpithelium': hd.color.CIELabColor.from_string("darkmagenta"),
    'NormalEpithelium': hd.color.CIELabColor.from_string("green"),
    'OtherCell': hd.color.CIELabColor.from_string("dimgray"),
    'StromalCellNOS': hd.color.CIELabColor.from_string("plum"),
    'TILsCell': hd.color.CIELabColor.from_string("royalblue"),
    'UnknownOrAmbiguousCell': hd.color.CIELabColor.from_string("gray"),
}


region_finding_codes = {
    "Cancerous epithelium": (
        hd.sr.CodedConcept(
            meaning="Malignant epithelial neoplasm",
            value="1187225007",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Morphologic abnormality",
            value="49755003",
            scheme_designator="SCT",
        ),
    ),
    "Stroma": (
        hd.sr.CodedConcept(
            meaning="Breast stroma",
            value="314375000",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "TILs": (
        hd.sr.CodedConcept(
            meaning="Lymphocyte",
            value="56972008",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "Normal epithelium": (
        hd.sr.CodedConcept(
            meaning="Epithelium",
            value="31610004",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "Junk/Debris": (
        hd.sr.CodedConcept(
            meaning="Debris",
            value="257159000",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Substance",
            value="105590001",
            scheme_designator="SCT",
        ),
    ),
    "Blood": (
        hd.sr.CodedConcept(
            meaning="Blood",
            value="87612001",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Substance",
            value="105590001",
            scheme_designator="SCT",
        ),
    ),
    "Other": (
        hd.sr.CodedConcept(
            meaning="Other",
            value="74964007",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Substance",
            value="105590001",
            scheme_designator="SCT",
        ),
    ),
    "Whitespace/Empty": (
        hd.sr.CodedConcept(
            meaning="Background",
            value="125040",
            scheme_designator="DCM",
        ),
        hd.sr.CodedConcept(
            meaning="Spatial and relational concept",
            value="309825002",
            scheme_designator="SCT",
        ),
    ),
}

nuclei_finding_codes = {
    "Cancer nucleus": (
        hd.sr.CodedConcept(
            meaning="Malignant epithelial neoplasm",
            value="1187225007",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Morphologic abnormality",
            value="49755003",
            scheme_designator="SCT",
        ),
    ),
    "Stromal nucleus": (
        hd.sr.CodedConcept(
            meaning="Breast stroma",
            value="314375000",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "Large stromal nucleus": (
        hd.sr.CodedConcept(
            meaning="Breast stroma",
            value="314375000",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "Lymphocyte nucleus": (
        hd.sr.CodedConcept(
            meaning="Lymphocyte",
            value="56972008",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "Plasma cell/large TIL nucleus": (
        hd.sr.CodedConcept(
            meaning="Plasma cell",
            value="113335003",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "Normal epithelial nucleus": (
        hd.sr.CodedConcept(
            meaning="Epithelial cell",
            value="4212006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "Other nucleus": (
        hd.sr.CodedConcept(
            meaning="Other",
            value="74964007",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "Unknown/ambiguous nucleus": (
        hd.sr.CodedConcept(
            meaning="Unknown",
            value="261665006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    "Background (non-nuclear material)": (
        hd.sr.CodedConcept(
            meaning="Background",
            value="125040",
            scheme_designator="DCM",
        ),
        hd.sr.CodedConcept(
            meaning="Spatial and relational concept",
            value="309825002",
            scheme_designator="SCT",
        ),
    ),
}


border_finding_codes = {
    "Border": (
        hd.sr.CodedConcept(
            meaning="Interface",
            value="112082",
            scheme_designator="DCM",
        ),
        hd.sr.CodedConcept(
            meaning="Spatial and relational concept",
            value="309825002",
            scheme_designator="SCT",
        ),
    ),
}


csv_finding_codes = {
    'ActiveStromalCellNOS': (
        hd.sr.CodedConcept(
            meaning="Breast stroma",
            value="314375000",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    'ActiveTILsCell': (
        hd.sr.CodedConcept(
            meaning="Lymphocyte",
            value="56972008",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    'CancerEpithelium': (
        hd.sr.CodedConcept(
            meaning="Malignant epithelial neoplasm",
            value="1187225007",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Morphologic abnormality",
            value="49755003",
            scheme_designator="SCT",
        ),
    ),
    'NormalEpithelium': (
        hd.sr.CodedConcept(
            meaning="Epithelial cell",
            value="4212006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    'OtherCell': (
        hd.sr.CodedConcept(
            meaning="Other",
            value="74964007",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    'StromalCellNOS': (
        hd.sr.CodedConcept(
            meaning="Breast stroma",
            value="314375000",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    'TILsCell': (
        hd.sr.CodedConcept(
            meaning="Lymphocyte",
            value="56972008",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
    'UnknownOrAmbiguousCell': (
        hd.sr.CodedConcept(
            meaning="Unknown",
            value="261665006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Anatomical Structure",
            value="91723000",
            scheme_designator="SCT",
        ),
    ),
}
