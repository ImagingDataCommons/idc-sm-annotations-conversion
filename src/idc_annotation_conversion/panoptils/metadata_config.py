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
region_bootstrapped_series_description = "PanopTILs Bootstrapped Region Segmentations"
nuclei_bootstrapped_series_description = "PanopTILs Bootstrapped Nuclei Segmentations"
border_bootstrapped_series_description = "PanopTILs Bootstrapped Border Segmentations"

# Annotation-specific metadata
ann_content_description = "Cell type annotations"
ann_series_description_manual = "PanopTILs Manual Cell Type Annotations"
ann_series_description_boostrapped = "PanopTILs Bootstrapped Cell Type Annotations"

### CODES
# Codes used for findings
MALIGNANT_EPITHELIAL_NEOPLASM_CODE = hd.sr.CodedConcept(
    meaning="Malignant epithelial neoplasm",
    value="1187225007",
    scheme_designator="SCT",
)
STROMA_CODE = hd.sr.CodedConcept(
    meaning="Breast stroma",
    value="314375000",
    scheme_designator="SCT",
)
LYMPHOCYTE_CODE = hd.sr.CodedConcept(
    meaning="Lymphocyte",
    value="56972008",
    scheme_designator="SCT",
)
EPITHELIUM_CODE = hd.sr.CodedConcept(
    meaning="Epithelium",
    value="31610004",
    scheme_designator="SCT",
)
DEBRIS_CODE = hd.sr.CodedConcept(
    meaning="Debris",
    value="257159000",
    scheme_designator="SCT",
)
PLASMA_CELL_CODE = hd.sr.CodedConcept(
    meaning="Plasma cell",
    value="113335003",
    scheme_designator="SCT",
)
EPITHELIAL_CELL_CODE = hd.sr.CodedConcept(
    meaning="Epithelial cell",
    value="4212006",
    scheme_designator="SCT",
)
OTHER_CODE = hd.sr.CodedConcept(
    meaning="Other",
    value="74964007",
    scheme_designator="SCT",
)
BACKGROUND_CODE = hd.sr.CodedConcept(
    meaning="Background",
    value="125040",
    scheme_designator="DCM",
)
SPATIAL_CONCEPT_CODE = hd.sr.CodedConcept(
    meaning="Spatial and relational concept",
    value="309825002",
    scheme_designator="SCT",
)
UNKNOWN_CODE = hd.sr.CodedConcept(
    meaning="Unknown",
    value="261665006",
    scheme_designator="SCT",
)
OTHER_CODE = hd.sr.CodedConcept(
    meaning="Other",
    value="74964007",
    scheme_designator="SCT",
)
SUBSTANCE_CODE = hd.sr.CodedConcept(
    meaning="Substance",
    value="105590001",
    scheme_designator="SCT",
)
ANATOMICAL_STRUCTURE_CODE = hd.sr.CodedConcept(
    meaning="Anatomical Structure",
    value="91723000",
    scheme_designator="SCT",
)
ABNORMALITY_CODE = hd.sr.CodedConcept(
    meaning="Morphologic abnormality",
    value="49755003",
    scheme_designator="SCT",
)
NUCLEUS_CODE = hd.sr.CodedConcept(
    meaning='Nucleus',
    value='84640000',
    scheme_designator='SCT',
)
INTERFACE_CODE = hd.sr.CodedConcept(
    meaning="Interface",
    value="112082",
    scheme_designator="DCM",
)
BLOOD_CODE = hd.sr.CodedConcept(
    meaning="Blood",
    value="87612001",
    scheme_designator="SCT",
)

### REGIONS FINDINGS

# Define the label, property type, category, primary anatomic structure, and
# RGB color for each finding
CANCER_EPITHELIUM_REGION_FINDING = (
    "Cancerous epithelium",
    MALIGNANT_EPITHELIAL_NEOPLASM_CODE,
    ABNORMALITY_CODE,
    None,
    [192, 56, 255],
)
STROMA_REGION_FINDING = (
    "Stroma",
    STROMA_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    None,
    [224, 133, 0],
)
TILS_REGION_FINDING = (
    "TILs",
    LYMPHOCYTE_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    None,
    [128, 128, 255],
)
NORMAL_EPITHELIUM_REGION_FINDING = (
    "Normal epithelium",
    EPITHELIUM_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    None,
    [53, 191,11],
)
DEBRIS_REGION_FINDING = (
    "Junk/debris",
    DEBRIS_CODE,
    SUBSTANCE_CODE,
    None,
    [239, 255, 138],
)
BLOOD_REGION_FINDING = (
    "Blood",
    BLOOD_CODE,
    SUBSTANCE_CODE,
    None,
    [255, 33, 33],
)
OTHER_REGION_FINDING = (
    "Other",
    OTHER_CODE,
    SUBSTANCE_CODE,
    None,
    [168, 168, 168],
)
WHITESPACE_REGION_FINDING = (
    "Whitespace/empty",
    BACKGROUND_CODE,
    SPATIAL_CONCEPT_CODE,
    None,
    [219, 219, 219],
)

png_region_findings_list = [
    CANCER_EPITHELIUM_REGION_FINDING,
    STROMA_REGION_FINDING,
    TILS_REGION_FINDING,
    NORMAL_EPITHELIUM_REGION_FINDING,
    DEBRIS_REGION_FINDING,
    BLOOD_REGION_FINDING,
    OTHER_REGION_FINDING,
    WHITESPACE_REGION_FINDING,
]

### NUCLEI FINDINGS

# Define the label, property type, category, primary anatomic structure, and
# RGB color for each finding
CANCER_NUCLEUS_FINDING = (
    "Cancer cell nucleus",
    NUCLEUS_CODE,
    ABNORMALITY_CODE,
    MALIGNANT_EPITHELIAL_NEOPLASM_CODE,
    [132, 0, 194],
)
STROMAL_NUCLEUS_FINDING = (
    "Stromal cell nucleus",
    NUCLEUS_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    STROMA_CODE,
    [209, 77, 255],
)
ACTIVE_STROMAL_NUCLEUS_FINDING = (
    "Active stromal cell nucleus",
    NUCLEUS_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    STROMA_CODE,  # TODO how to differentiate from normal stroms?
    [245, 22, 200],
)
TIL_NUCLEUS_FINDING = (
    "Tumor-infiltrating lymphocyte nucleus",
    NUCLEUS_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    LYMPHOCYTE_CODE,
    [0, 0, 255],
)
ACTIVE_TIL_NUCLEUS_FINDING = (
    "Active tumor-infiltrating lymphocyte nucleus",
    NUCLEUS_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    LYMPHOCYTE_CODE,  # TODO differentiate from normal TIL
    [0, 255, 255],
)
NORMAL_EPITHELIAL_NUCLEUS_FINDING = (
    "Normal epithelial cell nucleus",
    NUCLEUS_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    EPITHELIAL_CELL_CODE,
    [12, 100, 23],
)
OTHER_NUCLEUS_FINDING = (
    "Other cell nucleus",
    NUCLEUS_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    OTHER_CODE,
    [150, 150, 150],
)
UNKNOWN_NUCLEUS_FINDING = (
    "Unknown/ambiguous cell nucleus",
    NUCLEUS_CODE,
    ANATOMICAL_STRUCTURE_CODE,
    UNKNOWN_CODE,
    [100, 100, 100],
)
BACKGROUND_FINDING = (
    "Background",
    BACKGROUND_CODE,
    SPATIAL_CONCEPT_CODE,
    None,
    [30, 30, 30],
)

png_nuclei_findings_list = [
    CANCER_NUCLEUS_FINDING,
    STROMAL_NUCLEUS_FINDING,
    ACTIVE_STROMAL_NUCLEUS_FINDING,
    TIL_NUCLEUS_FINDING,
    ACTIVE_TIL_NUCLEUS_FINDING,
    NORMAL_EPITHELIAL_NUCLEUS_FINDING,
    OTHER_NUCLEUS_FINDING,
    UNKNOWN_NUCLEUS_FINDING,
    BACKGROUND_FINDING,
]

# Findings in the CSV match the nuclei, but are indexed by names, which do not
# match
csv_finding_mapping = {
    'CancerEpithelium': CANCER_NUCLEUS_FINDING,
    'StromalCellNOS': STROMAL_NUCLEUS_FINDING,
    'ActiveStromalCellNOS': ACTIVE_STROMAL_NUCLEUS_FINDING,
    'TILsCell': TIL_NUCLEUS_FINDING,
    'ActiveTILsCell': ACTIVE_TIL_NUCLEUS_FINDING,
    'NormalEpithelium': NORMAL_EPITHELIAL_NUCLEUS_FINDING,
    'OtherCell': OTHER_NUCLEUS_FINDING,
    'UnknownOrAmbiguousCell': UNKNOWN_NUCLEUS_FINDING,
}

### NUCLEI FINDINGS

# Define the label, property type, category, primary anatomic structure, and
# RGB color for each finding
BORDER_FINDING = (
    "Border",
    INTERFACE_CODE,
    SPATIAL_CONCEPT_CODE,
    None,
    [0, 255, 0],
)

png_border_findings_list = [
    BORDER_FINDING,
]
