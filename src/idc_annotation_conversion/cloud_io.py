"""Utilities for reading and writing objects to/from Google Cloud."""
from io import BytesIO

from google.cloud import storage

import pydicom


def read_dataset_from_blob(
    bucket: storage.Bucket,
    blob_name: str,
    stop_before_pixels: bool = False,
) -> pydicom.Dataset:
    """Read a pydicom Dataset from a bucket.

    Parameters
    ----------
    bucket: storage.Bucket
        Bucket object where the blob is stored.
    blob_name: str
        Name of the blob within the bucket.
    stop_before_pixels: bool
        Whether to stop before reading in the pixel data. I.e. return metadata
        only.

    Returns
    -------
    pydicom.Dataset
        Dataset loaded from the specified blob.

    """
    blob = bucket.get_blob(blob_name)
    dcm_bytes = blob.download_as_bytes()
    dcm = pydicom.dcmread(
        BytesIO(dcm_bytes),
        stop_before_pixels=stop_before_pixels,
    )
    return dcm


def write_dataset_to_blob(
    dataset: pydicom.Dataset,
    bucket: storage.Bucket,
    blob_name: str
) -> None:
    """Write a pydicom Dataset to a bucket.

    Parameters
    ----------
    dataset: pydicom.Dataset
        Dataset object to upload.
    bucket: storage.Bucket
        Bucket object where the blob should be stored.
    blob_name: str
        Name of the blob within the bucket. If it already exists, it will be
        overwritten.

    """
    blob = bucket.blob(blob_name)
    with BytesIO() as buf:
        dataset.save_as(buf)
        buf.seek(0)
        blob.upload_from_file(buf)
