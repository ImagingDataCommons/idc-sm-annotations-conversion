"""Metadata used for RMS conversions."""
import numpy as np
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
    "TIL Positive": (
        hd.sr.CodedConcept(
            meaning="Placeholder",
            value="001",
            scheme_designator="X",
        ),
        hd.sr.CodedConcept(
            meaning="Placeholder",
            value="002",
            scheme_designator="X",
        ),
    ),
    "TIL Negative": (
        hd.sr.CodedConcept(
            meaning="Placeholder",
            value="003",
            scheme_designator="X",
        ),
        hd.sr.CodedConcept(
            meaning="Placeholder",
            value="004",
            scheme_designator="X",
        ),
    ),
}

software_versions = get_git_remote_url(simplify=True)
device_serial_number = get_git_commit_hash()
institution_name = None
institutional_department_name = None

segmentation_channel_order = [
    "TIL Negative",
    "TIL Positive",
]
algorithm_identification = hd.AlgorithmIdentificationSequence(
    name="Saltz et al TIL Model",
    version="V1.0",
    source="Stony Brook University",
    family=codes.cid7162.ArtificialIntelligence,
)
segmentation_series_description = "Automatic TIL Map"
seg_manufacturer = "Stony Brook University"
seg_manufacturer_model_name = "Saltz et al TIL Model"

labelmap_lut = hd.PaletteColorLUTTransformation.from_colors(
    ['black', 'blue', 'red'],
    palette_color_lut_uid=hd.UID(),
)
