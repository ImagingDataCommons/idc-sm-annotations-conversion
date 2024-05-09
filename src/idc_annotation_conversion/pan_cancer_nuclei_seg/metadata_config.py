import highdicom as hd
import pydicom
from pydicom.sr.codedict import codes
from pydicom.sr.coding import Code

from idc_annotation_conversion.git_utils import (
    get_git_remote_url,
    get_git_commit_hash,
)


# Basic Metadata
manufacturer = "Stony Brook University"
manufacturer_model_name = "Pan-Cancer-Nuclei-Seg"
software_versions = get_git_remote_url(simplify=True)
device_serial_number = get_git_commit_hash()

# Label description
label = "Nuclei"
finding_category = Code("91723000", "SCT", "Anatomical Stucture")
finding_type = Code("84640000", "SCT", "Nucleus")

# Algorithm Identification
algorithm_identification = hd.AlgorithmIdentificationSequence(
    name="Pan-Cancer-Nuclei-Seg",
    version="1.0",
    source="https://doi.org/10.7937/TCIA.2019.4A4DKP9U",
    family=codes.cid7162.ArtificialIntelligence,
)

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
clinical_trial_ids_item.ClinicalTrialProtocolID = "doi:xx/xxx"
other_trials_seq_element = pydicom.DataElement(
    other_trials_seq_tag,
    "SQ",
    [clinical_trial_ids_item],
)
