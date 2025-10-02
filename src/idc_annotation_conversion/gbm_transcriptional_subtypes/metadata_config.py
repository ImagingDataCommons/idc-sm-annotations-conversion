"""Metadata used for transcriptomic subtype map conversions."""
import highdicom as hd
from highdicom.color import CIELabColor
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
    name="GBM360",
    version="1.0",
    source="https://github.com/gevaertlab/GBM360/",
    family=codes.cid7162.ArtificialIntelligence,
)
segmentation_series_description = "Automated Transcriptional Subtype Map"
seg_manufacturer = "Stanford University"
seg_manufacturer_model_name = "Transcriptional Subtypes CNN"

labelmap_lut = hd.PaletteColorLUTTransformation.from_colors(
    ['black', 'red', 'brown', 'purple', 'blue', 'green', 'white'],
    palette_color_lut_uid=hd.UID(),
)

display_colors = {
    "OPClike": CIELabColor.from_string("red"),
    "MESlike1": CIELabColor.from_string("brown"),
    "MESlike2": CIELabColor.from_string("purple"),
    "NPClike": CIELabColor.from_string("blue"),
    "AClike": CIELabColor.from_string("green"),
    "Normal": CIELabColor.from_string("white"),
}


# Aggressiveness maps
pmap_manufacturer = "Imaging Data Commons"
pmap_manufacturer_model_name = "IDC SM Annotation Conversion"
pmap_software_versions = get_git_remote_url(simplify=True)
pmap_device_serial_number = get_git_commit_hash()
pmap_series_description = "Aggressiveness Map"
pmap_content_label = "AGGRESSIVENESS"
pmap_content_description = "Map of aggressiveness"
pmap_image_flavor = "OTHER"
pmap_real_world_value_mappings = [
    hd.pm.RealWorldValueMapping(
        lut_label="Aggressiveness",
        lut_explanation="Aggressiveness scores calculated by a neural network",
        unit=codes.UCUM.NoUnits,
        value_range=(0.0, 1.0),
        slope=1.0,
        intercept=0.0,
        quantity_definition=hd.sr.CodedConcept(
            meaning='Severity',
            value='246112005',
            scheme_designator='SCT',
        ),
    )
]
# derived_pixel_contrast = ""
