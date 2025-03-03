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


# Dictionary mapping text label found in the XML annotations to tuple of
# (finding_type, finding_category) codes to encode that finding
finding_codes = {
    "OPClike": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="placeholder",
            value="123",
            scheme_designator="SCT",
        ),
    ),
    "MESlike1": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="placeholder",
            value="123",
            scheme_designator="SCT",
        ),
    ),
    "MESlike2": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="placeholder",
            value="123",
            scheme_designator="SCT",
        ),
    ),
    "NPClike": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="placeholder",
            value="123",
            scheme_designator="SCT",
        ),
    ),
    "AClike": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="placeholder",
            value="123",
            scheme_designator="SCT",
        ),
    ),
    "Normal": (
        hd.sr.CodedConcept(
            meaning="Function",
            value="246464006",
            scheme_designator="SCT",
        ),
        hd.sr.CodedConcept(
            meaning="placeholder",
            value="123",
            scheme_designator="SCT",
        ),
    ),
}

segmentation_channel_order = [
    "OPClike",
    "MESlike1",
    "MESlike2",
    "NPClike",
    "AClike",
    "Normal",
]
algorithm_identification = hd.AlgorithmIdentificationSequence(
    name="Transcriptional Subtypes CNN",
    version="1.0",
    source="doi:10.1038/s41467-023-39933-0",
    family=codes.cid7162.ArtificialIntelligence,
)
segmentation_series_description = "Automated Transcriptional Subtype Map"
seg_manufacturer = "Stanford University"
seg_manufacturer_model_name = "Transcriptional Subtypes CNN"

labelmap_lut = hd.PaletteColorLUTTransformation.from_colors(
    ['black', 'red', 'brown', 'purple', 'blue', 'green', 'white'],
    palette_color_lut_uid=hd.UID(),
)
