"""Metadata used for RMS conversions."""
import numpy as np
import highdicom as hd
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
# (finding_type, finding_category) codes to encode that finding
finding_codes_2018 = {
    "TIL Positive": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Tumor infiltration by lymphocytes present",
            value="399721002",
            scheme_designator="SCT",
        ),
    ),
    "TIL Negative": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Tumor infiltration by lymphocytes absent",
            value="396396002",
            scheme_designator="SCT",
        ),
    ),
}

segmentation_channel_order = [
    "TIL Negative",
    "TIL Positive",
]
algorithm_identification_2018 = hd.AlgorithmIdentificationSequence(
    name="TIL CNN",
    version="1.0",
    source="doi:10.7937/K9/TCIA.2018.Y75F9W1",
    family=codes.cid7162.ArtificialIntelligence,
)
segmentation_series_description_2018 = "Automatic TIL Map"
seg_manufacturer_2018 = "Stony Brook University"
seg_manufacturer_model_name_2018 = "TIL CNN"

labelmap_lut = hd.PaletteColorLUTTransformation.from_colors(
    ['black', 'blue', 'red'],
    palette_color_lut_uid=hd.UID(),
)

# 2022 TIL Maps
# =============
# These are the ones from this paper:

# Abousamra, S., Gupta, R., Hou, L., Batiste, R., Zhao, T., Shankar, A., Saltz,
# J. (2022). Deep learning-based mapping of tumor infiltrating lymphocytes in
# whole slide images of 23 types of cancer. Frontiers in oncology, 11, 806603.

# Dictionary mapping text label found in the XML annotations to tuple of
# (finding_type, finding_category) codes to encode that finding
finding_codes_2022 = {
    "TIL Positive": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Tumor infiltration by lymphocytes present",
            value="399721002",
            scheme_designator="SCT",
        ),
    ),
    "TIL Negative": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="Tumor infiltration by lymphocytes absent",
            value="396396002",
            scheme_designator="SCT",
        ),
    ),
}

algorithm_identification_2022 = hd.AlgorithmIdentificationSequence(
    name="Improved TIL CNN",
    version="1.0",
    source="doi:10.3389/fonc.2021.806603",
    family=codes.cid7162.ArtificialIntelligence,
)
segmentation_series_description_2022 = "Improved Automatic TIL Map"
seg_manufacturer_2022 = "Stony Brook University"
seg_manufacturer_model_name_2022 = "Improved TIL CNN"
