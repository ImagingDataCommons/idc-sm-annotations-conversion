"""Metadata used for TCGA TIL map conversions."""
import highdicom as hd
from highdicom.color import CIELabColor
import pydicom
from pydicom.sr.codedict import codes

from idc_annotation_conversion.git_utils import (
    get_git_remote_url,
    get_git_commit_hash,
)

# Shared metadata
software_versions = get_git_remote_url(simplify=True)
device_serial_number = get_git_commit_hash()


# 2018 TIL Maps
# =============
# These are the ones from this paper:

# Saltz, Joel, et al. "Spatial organization and molecular correlation of
# tumor-infiltrating lymphocytes using deep learning on pathology images."
# Cell reports 23.1 (2018): 181-193.

# Dictionary mapping text label found in the XML annotations to tuple of
# (finding_category, finding_type) codes to encode that finding
finding_codes_2018 = {
    "TILs Present": (
        hd.sr.CodedConcept(
            meaning="Morphologically abnormal structure",
            value="49755003",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Tumor infiltration by lymphocytes present",
            value="399721002",
            scheme_designator="SCT",
        ),
    ),
    "TILs Absent": (
        hd.sr.CodedConcept(
            meaning="Morphologically abnormal structure",
            value="49755003",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Tumor infiltration by lymphocytes absent",
            value="396396002",
            scheme_designator="SCT",
        ),
    ),
}

segmentation_channel_order_2018 = [
    "TILs Absent",
    "TILs Present",
]
algorithm_identification_2018 = hd.AlgorithmIdentificationSequence(
    name="Stony Brook TIL Segmentation CNN 2018",
    version="1.0",
    source="doi:10.7937/K9/TCIA.2018.Y75F9W1",
    family=codes.cid7162.ArtificialIntelligence,
    parameters={"patch size": "50 x 50 microns"},
)
segmentation_series_description_2018 = "Stony Brook CNN-generated TIL Map"
seg_manufacturer_2018 = "Stony Brook University converted by Imaging Data Commons"
seg_manufacturer_model_name_2018 = "TIL Custom CNN 2018 converted by Imaging Data Commons"

labelmap_lut = hd.PaletteColorLUTTransformation.from_colors(
    ['black', 'blue', 'red'],
    palette_color_lut_uid=hd.UID(),
)

display_colors = {
    "TILs Present": CIELabColor.from_rgb(255, 0, 0),  # red
    "TILs Absent": CIELabColor.from_rgb(0, 0, 255),  # blue
}

# 2022 TIL Maps
# =============
# These are the ones from this paper:

# Abousamra, S., Gupta, R., Hou, L., Batiste, R., Zhao, T., Shankar, A., Saltz,
# J. (2022). Deep learning-based mapping of tumor infiltrating lymphocytes in
# whole slide images of 23 types of cancer. Frontiers in oncology, 11, 806603.

# Dictionary mapping text label found in the XML annotations to tuple of
# (finding_type, finding_category) codes to encode that finding
finding_codes_2022 = {
    "TILs Present": (
        hd.sr.CodedConcept(
            meaning="Morphologically abnormal structure",
            value="49755003",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Tumor infiltration by lymphocytes present",
            value="399721002",
            scheme_designator="SCT",
        ),
    ),
}

segmentation_channel_order_2022 = ["TILs Present"]
algorithm_identification_2022 = hd.AlgorithmIdentificationSequence(
    name="Stony Brook TIL Segmentation Inception-V4 2022",
    version="1.0",
    source="doi:10.3389/fonc.2021.806603",
    family=codes.cid7162.ArtificialIntelligence,
    parameters={"patch size": "50 x 50 microns"},
)
segmentation_series_description_2022_binary = "Stony Brook Inception-V4 Binary TIL Map"
segmentation_series_description_2022_fractional = "Stony Brook Inception-V4 Fractional TIL Map"
seg_manufacturer_2022 = "Stony Brook University converted by Imaging Data Commons"
seg_manufacturer_model_name_2022 = "TIL Inception-V4 2022 converted by Imaging Data Commons"

# DOI of the conversion page in Zenodo for other clinical trial protocol
# These are not yet in pydicom's data dict so we have to set the elements
# manually
# IssuerOfClinicalTrialProtocolID
issuer_tag = pydicom.tag.Tag(0x0012, 0x0022)
# OtherClinicalTrialProtocolIDsSequence
other_trials_seq_tag = pydicom.tag.Tag(0x0012, 0x0023)

clinical_trial_ids_item = pydicom.Dataset()
issuer_value = "DOI"
clinical_trial_ids_item.add(
    pydicom.DataElement(
        issuer_tag,
        "LO",
        issuer_value,
    )
)
clinical_trial_ids_item.ClinicalTrialProtocolID = "doi:10.5281/zenodo.16966285"
other_trials_seq_element = pydicom.DataElement(
    other_trials_seq_tag,
    "SQ",
    [clinical_trial_ids_item],
)
