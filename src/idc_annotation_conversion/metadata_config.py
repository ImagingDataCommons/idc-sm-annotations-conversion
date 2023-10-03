import highdicom as hd
from pydicom.sr.codedict import codes
from pydicom.sr.coding import Code

from idc_annotation_conversion.git_utils import (
    get_git_remote_url,
    get_git_commit_hash,
)


def simplify_remote(remote: str) -> str:
    """Simplify a remote URL if necessary to keep it below 64 characters."""
    if len(remote) <= 64:
        return remote

    # Strip from start and end
    if remote.endswith(".git"):
        remote = remote[:-4]
    for start_str in ["http://", "https://", "git@"]:
        if remote.startswith(start_str):
            remote = remote[len(start_str):]

    if len(remote) <= 64:
        return remote

    remote = remote.replace("github.com:", "github:")
    remote = remote.replace("github.com/", "github:")

    if len(remote) > 64:
        raise ValueError(
            "Cannot simplify URL of the remote to be 64 characters or fewer."
        )

    return remote


# Basic Metadata
manufacturer = "Stony Brook University"
manufacturer_model_name = "Pan-Cancer-Nuclei-Seg"
software_versions = simplify_remote(get_git_remote_url())
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
