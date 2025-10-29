"""Metadata used for transcriptomic subtype map conversions."""
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
segmentation_series_description = "Automated Transcriptional Subtype Map"
seg_manufacturer = "Gevaert Lab Converted By Imaging Data Commons"
seg_manufacturer_model_name = "GBM360"

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
pmap_manufacturer = "Gevaert Lab Converted By Imaging Data Commons"
pmap_manufacturer_model_name = "GBM360"
pmap_software_versions = get_git_remote_url(simplify=True)
pmap_device_serial_number = get_git_commit_hash()
pmap_series_description = "Aggressiveness Score Map"
pmap_content_label = "AGGRESSIVENESS"
pmap_content_description = "Map of aggressiveness scores"
pmap_image_flavor = "AGGRESSIVENESS"
pmap_real_world_value_mappings = [
    hd.pm.RealWorldValueMapping(
        lut_label="Aggressiveness",
        lut_explanation="Aggressiveness scores calculated by a neural network",
        unit=codes.UCUM.NoUnits,
        value_range=(0.0, 1.0),
        slope=1.0,
        intercept=0.0,
        quantity_definition=hd.sr.CodedConcept(
            meaning='Aggressiveness score',
            value='314001',
            scheme_designator='99PMP',
        ),
    )
]
derived_pixel_contrast = "AI"

pmap_contributing_equipment = [
    hd.ContributingEquipment(
        manufacturer="Gevaert Lab",
        manufacturer_model_name="GBM360",
        institution_name="Stanford University",
        purpose_of_reference=codes.DCM.ProcessingEquipment,
        contribution_description="Generation of aggressiveness map",
        software_versions="https://github.com/gevaertlab/GBM360/",
    ),
]

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
clinical_trial_ids_item.ClinicalTrialProtocolID = "10.5281/zenodo.17470190"
other_trials_seq_element = pydicom.DataElement(
    other_trials_seq_tag,
    "SQ",
    [clinical_trial_ids_item],
)
