import highdicom as hd
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
