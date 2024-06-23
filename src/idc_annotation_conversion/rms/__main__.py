"""Main conversion process for the RMS project."""
import datetime
import logging
import os
from io import BytesIO
from itertools import islice
from pathlib import Path
from typing import Optional, List
from xml.etree import ElementTree

import click
from google.cloud import bigquery, storage
import numpy as np
import pandas as pd
from PIL import Image

from idc_annotation_conversion import cloud_config, cloud_io
from idc_annotation_conversion.rms.convert import convert_xml_annotation


# Bucket containing annotations for this project
ANNOTATION_BUCKET_PROJECT = "idc-external-031"
ANNOTATION_BUCKET = "rms_annotation_test_oct_2023"


def find_series(
    container_prefix: str,
    bq_client: bigquery.Client,
) -> pd.DataFrame:
    """Find the DICOM WSI series given the container prefix.

    Parameters
    ----------
    container_prefix: str
        Container prefix.
    bq_client: google.cloud.bigquery.Client
        Existing BigQuery client object.

    Returns
    -------
    pandas.DataFrame:
        DataFrame containing information on each instance in the matching
        series.

    """
    selection_query = f"""
        SELECT
            ContainerIdentifier,
            crdc_instance_uuid,
            crdc_series_uuid,
            LossyImageCompression,
            TotalPixelMatrixRows,
            TotalPixelMatrixColumns,
            Cast(NumberOfFrames AS int) AS NumberOfFrames
        FROM
            bigquery-public-data.idc_v16.dicom_all
        WHERE
            ContainerIdentifier LIKE "{container_prefix}%"
            AND collection_id = "rms_mutation_prediction"
            AND ARRAY_TO_STRING(ImageType,",") LIKE "%VOLUME%"
            AND LossyImageCompression = "00"
        ORDER BY
            NumberOfFrames DESC
    """
    selection_result = bq_client.query(selection_query)
    selection_df = selection_result.result().to_dataframe()
    return selection_df


@click.command()
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path, file_okay=False),
    help="Output directory, default: no output directory",
)
@click.option(
    "--output-bucket",
    "-b",
    help="Output bucket",
    show_default=True,
)
@click.option(
    "--store-bucket/--no-store-bucket",
    "-k/-K",
    help="Whether to store outputs to the bucket in this run.",
    default=True,
    show_default=True,
)
@click.option(
    "--output-prefix",
    "-p",
    help="Prefix for all output blobs. Default no prefix.",
)
@click.option(
    "--store-wsi-dicom/--no-store-wsi-dicom",
    "-d/-D",
    default=False,
    show_default=True,
    help=(
        "Download all WSI DICOM files and store in the output directory "
        "(if any)."
    ),
)
@click.option(
    "--number",
    "-n",
    type=int,
    help="Number to process per collection. All by default.",
)
@click.option(
    "--use-scoord3d/--no-use-scoord3d",
    "-s/-S",
    default=True,
    show_default=True,
    help="Store coordinates as SCOORD3D (versus 2D SCOORD).",
)
@click.option(
    "--include-measurements/--no-include-measurements",
    "-m/-M",
    default=False,
    show_default=True,
    help="Include length and area measurements from the XML.",
)
@click.option(
    "--excluded-cases",
    "-e",
    multiple=True,
    help="Container IDs (before underscore) to skip. Multiple may be provided",
)
def run(
    number: Optional[int],
    output_dir: Optional[Path],
    output_bucket: str,
    output_prefix: Optional[str],
    store_bucket: bool,
    store_wsi_dicom: bool,
    use_scoord3d: bool,
    include_measurements: bool,
    excluded_cases: Optional[List[str]] = None,
):
    """Main process for conversion of RMS XML annotations to DICOM SRs."""
    logging.basicConfig(level=logging.INFO)

    # Suppress highdicom logging (very talkative)
    logging.getLogger("highdicom.base").setLevel(logging.WARNING)
    logging.getLogger("highdicom.seg.sop").setLevel(logging.WARNING)

    # Setup project and authenticate
    os.environ["GCP_PROJECT_ID"] = cloud_config.GCP_PROJECT_ID

    if output_dir is not None:
        output_dir.mkdir(exist_ok=True)

    # Access bucket containing annotations
    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    public_bucket = storage_client.bucket(cloud_config.DICOM_IMAGES_BUCKET)

    ann_storage_client = storage.Client(project=ANNOTATION_BUCKET_PROJECT)
    ann_bucket = ann_storage_client.bucket(ANNOTATION_BUCKET)

    # Setup bigquery client
    bq_client = bigquery.Client(cloud_config.GCP_PROJECT_ID)

    prefix = "RMS-XML-hand-annotations"
    for ann_blob in islice(
        ann_bucket.list_blobs(prefix=prefix),
        1,  # first item is the directory itself
        number + 1 if number is not None else None,
    ):
        if not ann_blob.name.endswith(".xml"):
            continue

        logging.info(f"Processing annotation in {ann_blob.name}.")
        container_prefix, *_ = ann_blob.name.split("/")[1].replace(".xml", "").split("_", maxsplit=1)
        if excluded_cases is not None:
            if container_prefix in excluded_cases:
                logging.info(f"Skipping case {container_prefix}.")
                continue

        selection_df = find_series(container_prefix, bq_client)
        container_id = selection_df.ContainerIdentifier.iloc[0]

        assert selection_df.crdc_series_uuid.nunique() == 1, "Found multiple source series"

        text = ann_blob.download_as_text()
        xml_root = ElementTree.fromstring(text)

        logging.info("Retrieving source images.")
        wsi_dcm = [
            cloud_io.read_dataset_from_blob(
                bucket=public_bucket,
                blob_name=f"{row.crdc_series_uuid}/{row.crdc_instance_uuid}.dcm",
            )
            for _, row in selection_df.iterrows()
        ]

        logging.info("Creating SR.")
        sr_dcm = convert_xml_annotation(
            xml_root,
            wsi_dcm,
            use_scoord3d=use_scoord3d,
            include_measurements=include_measurements,
        )
        logging.info("SR created.")

        # Store objects to bucket
        if store_bucket:
            if output_bucket is None:
                data_str = (datetime.date.today())
                output_bucket = (
                    f"rms_manual_annotation_sr_conversion_{data_str}"
                )
            output_bucket_obj = output_client.bucket(output_bucket)

            if not output_bucket_obj.exists():
                output_bucket_obj.create(
                    location=cloud_config.GCP_DEFAULT_LOCATION
                )

            blob_root = (
                "" if output_prefix is None else f"{output_prefix}/"
            )
            sr_blob_name = (
                f"{blob_root}{container_id}_sr.dcm"
            )

            logging.info(f"Uploading SR to {sr_blob_name}.")
            cloud_io.write_dataset_to_blob(
                sr_dcm,
                output_bucket_obj,
                sr_blob_name,
            )

        # Store objects to filesystem
        if output_dir is not None:
            sr_path = output_dir / f"{container_id}_sr.dcm"

            logging.info(f"Writing sr to {str(sr_path)}.")
            sr_dcm.save_as(sr_path)

            if store_wsi_dicom:
                for i, dcm in enumerate(wsi_dcm):
                    dcm_path = output_dir / f"{container_id}_im{i}.dcm"
                    dcm.save_as(dcm_path)


@click.command()
@click.option(
    "--fractional/--binary",
    "-f/-F",
    help="Store as a fractional (probabilistic) mask.",
    show_default=True,
)
@click.option(
    "--number",
    "-n",
    type=int,
    help="Number to process per collection. All by default.",
)
@click.option(
    "--excluded-cases",
    "-e",
    multiple=True,
    help="Container IDs (before underscore) to skip. Multiple may be provided",
)
def run_seg(
    fractional: bool,
    number: Optional[int],
    excluded_cases: Optional[List[str]] = None,
):

    # Images are so large that they will trigger decompression errors unless
    # you do this...
    Image.MAX_IMAGE_PIXELS = 1000000000

    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    public_bucket = storage_client.bucket(cloud_config.DICOM_IMAGES_BUCKET)

    mask_storage_client = storage.Client(project=ANNOTATION_BUCKET_PROJECT)
    mask_bucket = mask_storage_client.bucket(ANNOTATION_BUCKET)

    # Setup bigquery client
    bq_client = bigquery.Client(cloud_config.GCP_PROJECT_ID)

    prefix = "HyunReferenceModelAllOutputs/predictions"
    for mask_blob in islice(
        mask_bucket.list_blobs(prefix=prefix),
        1,  # first item is the directory itself
        number + 1 if number is not None else None,
    ):
        if not mask_blob.name.endswith(".npy"):
            continue

        logging.info(f"Processing annotation in {mask_blob.name}.")
        container_prefix, *_ = mask_blob.name.split("/")[-1].replace(".npy", "").split("_", maxsplit=1)
        if excluded_cases is not None:
            if container_prefix in excluded_cases:
                logging.info(f"Skipping case {container_prefix}.")
                continue

        selection_df = find_series(container_prefix, bq_client)
        container_id = selection_df.ContainerIdentifier.iloc[0]

        assert selection_df.crdc_series_uuid.nunique() == 1, "Found multiple source series"

        mask_bytes = mask_blob.download_as_bytes()
        mask_im = np.load(BytesIO(mask_bytes))

        print("mask", mask_im.shape, mask_im.dtype)
        print(mask_im.min())
        print(mask_im.max())
        print("images")
        for _, row in selection_df.iterrows():
            print(row.TotalPixelMatrixRows, row.TotalPixelMatrixColumns)
        print()

        logging.info("Retrieving source images.")
        wsi_dcm = [
            cloud_io.read_dataset_from_blob(
                bucket=public_bucket,
                blob_name=f"{row.crdc_series_uuid}/{row.crdc_instance_uuid}.dcm",
            )
            for _, row in selection_df.iterrows()
        ]



if __name__ == "__main__":
    run_seg()
