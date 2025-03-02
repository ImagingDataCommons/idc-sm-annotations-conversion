"""Utilities for reading and writing objects to/from Google Cloud."""
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
    blob = bucket.get_blob(blob_name)
    dcm = pydicom.dcmread(
        blob.open(mode="rb"),
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
