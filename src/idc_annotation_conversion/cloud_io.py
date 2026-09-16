"""Utilities for reading and writing objects to/from Google Cloud."""
import hashlib
from io import BytesIO

from google.cloud import storage
from google.cloud.storage.fileio import BlobReader

import pydicom
import numpy as np
from PIL import Image


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
    chunk_size = None
    if stop_before_pixels:
        # Use a much smaller chunk size when just pulling metadata
        chunk_size = 500_000

    blob = bucket.get_blob(blob_name)
    dcm = pydicom.dcmread(
        blob.open(mode="rb", chunk_size=chunk_size),
        stop_before_pixels=stop_before_pixels,
    )
    return dcm


def write_dataset_to_blob(
    dataset: pydicom.Dataset,
    bucket: storage.Bucket,
    blob_name: str
) -> str:
    """Write a pydicom Dataset to a bucket and return hash.

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
        hash = hashlib.md5(buf.getvalue()).hexdigest()
        buf.seek(0)
        blob.upload_from_file(buf)

    return hash


def read_image_from_blob(
    bucket: storage.Bucket,
    blob_name: str,
) -> np.ndarray:
    """Read a non-DICOM image from a bucket.

    Parameters
    ----------
    bucket: storage.Bucket
        Bucket object where the blob is stored. Should be some image file in a
        format supported by PIL (e.g. JPEG, PNG, etc).
    blob_name: str
        Name of the blob within the bucket.

    Returns
    -------
    numpy.ndarray
        Dataset loaded from the specified blob.

    """
    blob = bucket.get_blob(blob_name)
    im = np.array(Image.open(blob.open("rb")))
    return im


def get_blob_uri(blob: storage.Blob):
    """Return a GCS URI to this blob

    Parameters
    ----------
    blob: storage.Blob
        Blob whose URI is sought.

    Returns
    -------
    str:
        URI in the form `gs://...`

    """
    return f"gs://{blob.bucket.name}/{blob.name}"
