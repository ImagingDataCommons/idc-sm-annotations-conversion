"""Main conversion process for the TCGA TIL Maps project."""
import datetime
import logging
from pathlib import Path
import os
from time import time
from typing import Any, Optional

import click
from google.cloud import bigquery
from google.cloud import storage
import highdicom as hd
import pandas as pd

from idc_annotation_conversion import cloud_io
from idc_annotation_conversion.tcga_til_maps.convert import (
    convert_segmentation,
)
from idc_annotation_conversion import cloud_config


IMAGE_BUCKET = "idc-open-data"
ANNOTATION_BUCKET = "til-wsi-tcga"
ANNOTATION_PREFIX = 'TIL_maps_after_thres_v1/'


COLLECTIONS = [
    "blca",
    "brca",
    "cesc",
    "coad",
    "luad",
    "lusc",
    "paad",
    "prad",
    "read",
    "skcm",
    "stad",
    "ucec",
    "uvm",
]


@click.command()
@click.option(
    "-c",
    "--collections",
    multiple=True,
    type=click.Choice(COLLECTIONS),
    help="Collections to use, all by default.",
    show_choices=True,
)
@click.option(
    "-f",
    "--blob-filter",
    help=(
        "Only process annotations blobs whose name contains this string."
    ),
)
@click.option(
    "--number",
    "-n",
    type=int,
    help="Number of annotations to process. All by default.",
)
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
    "--keep-existing/--overwrite-existing",
    "-m/-M",
    help="Only process a case if the output does not exist in the bucket.",
    default=False,
    show_default=True,
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
    "--dimension-organization-type",
    "-T",
    type=click.Choice(
        [v.name for v in hd.DimensionOrganizationTypeValues],
        case_sensitive=False,
    ),
    default="TILED_FULL",
    show_default=True,
    help=(
        "Dimension organization type for segmentations. TILED_FULL (default) "
        "or TILED_SPARSE."
    ),
)
@click.option(
    "--segmentation-type",
    "-t",
    type=click.Choice(
        [v.name for v in hd.seg.SegmentationTypeValues],
        case_sensitive=False,
    ),
    default="BINARY",
    show_default=True,
    help="Segmentation type for the Segmentation Image, if any.",
)
@click.option(
    "--include-lut/--no-include-lut",
    "-l/-L",
    default=True,
    show_default=True,
    help="Include a LUT in a labelmap seg.",
)
def run(
    collections: Optional[list[str]],
    number: Optional[int],
    output_dir: Optional[Path],
    output_bucket: str,
    store_bucket: bool,
    store_wsi_dicom: bool,
    segmentation_type: str,
    dimension_organization_type: str,
    blob_filter: Optional[str] = None,
    keep_existing: bool = False,
    include_lut: bool = True,
):
    """Convert TCGA tumor infiltrating lymphocyte (TIL) maps to DICOM
    segmentations.

    Images and annotations are automatically pulled down from cloud buckets as
    required.

    By default, output is to a cloud bucket. It is also possible to output
    to a local directory.

    """
    # Use all collections if none specified
    collections = collections or COLLECTIONS

    logging.basicConfig(level=logging.INFO)

    # Suppress highdicom logging (very talkative)
    logging.getLogger("highdicom.base").setLevel(logging.WARNING)
    logging.getLogger("highdicom.seg.sop").setLevel(logging.WARNING)

    # Setup project and authenticate
    os.environ["GCP_PROJECT_ID"] = cloud_config.GCP_PROJECT_ID
    bq_client = bigquery.Client(cloud_config.GCP_PROJECT_ID)

    if keep_existing and not store_bucket:
        raise ValueError("keep_existing requires store_bucket")

    # Access bucket containing annotations
    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    im_bucket = storage_client.bucket(IMAGE_BUCKET)
    ann_bucket = storage_client.bucket(ANNOTATION_BUCKET)
    output_bucket_obj = None
    if store_bucket:
        if output_bucket is None:
            date_str = (datetime.date.today())
            output_bucket = (
                f"tcga_til_maps_{date_str}"
            )
        output_bucket_obj = output_client.bucket(output_bucket)

        if not output_bucket_obj.exists():
            output_bucket_obj.create(
                location=cloud_config.GCP_DEFAULT_LOCATION
            )

    # Create output directory
    if output_dir is not None:
        output_dir.mkdir(exist_ok=True)

    logging.info("Listing PNGs")

    to_process = []

    # Loop over requested collections
    for collection in collections:
        collection_prefix = f'{ANNOTATION_PREFIX}{collection}/'

        collection_blobs = [
            b for b in ann_bucket.list_blobs(prefix=collection_prefix)
            if b.name.endswith('.png')
        ]

        if blob_filter is not None:
            collection_blobs = [
                b for b in collection_blobs if blob_filter in b.name
            ]

        for blob in collection_blobs:
            container_id = blob.name.split('/')[-1].replace('.png', '')

            # Check whether the output blobs already exist, and skip if
            # they do
            if keep_existing:
                seg_blob_name = f"{container_id}_seg.dcm"
                seg_output_blob = output_bucket_obj.get_blob(seg_blob_name)
                if seg_output_blob.exists():
                    continue

            to_process.append(blob)

    if number is not None:
        to_process = to_process[:number]

    errors = []

    for blob in to_process:
        container_id = blob.name.split('/')[-1].replace('.png', '')
        logging.info(f"Running on container ID {container_id}")

        mask = cloud_io.read_image_from_blob(
            ann_bucket,
            blob.name
        )
        logging.info("Pulled mask")

        selection_query = f"""
            SELECT
                gcs_url,
                CAST(NumberOfFrames AS int) AS NumberOfFrames
            FROM
                bigquery-public-data.idc_current.dicom_all
            WHERE
                ContainerIdentifier='{container_id}' and SamplesPerPixel = 3
            ORDER BY
                NumberOfFrames DESC
        """
        selection_result = bq_client.query(selection_query)
        selection_df = selection_result.result().to_dataframe()
        if len(selection_df) == 0:
            msg = f"ERROR cannot find image: {container_id}"
            errors.append(msg)
            logging.error(msg)
        if len(selection_df) > 1:
            # Check that there aren't multiple options
            if (
                selection_df.iloc[0].NumberOfFrames ==
                selection_df.iloc[1].NumberOfFrames
            ):
                msg = f"WARNING multiple images: {container_id}"
                logging.warning(msg)
                errors.append(msg)

        url = selection_df.iloc[0].gcs_url
        blob_name = "/".join(url.split("/")[3:])

        wsi_dcm = cloud_io.read_dataset_from_blob(im_bucket, blob_name)
        logging.info("Pulled WSI image")

        seg = convert_segmentation(
            mask,
            wsi_dcm,
            dimension_organization_type=dimension_organization_type,
            segmentation_type=segmentation_type,
            include_lut=include_lut,
        )

        # Store objects to filesystem
        if output_dir is not None:
            seg_path = output_dir / f"{container_id}_seg.dcm"

            logging.info(f"Writing segmentation to {str(seg_path)}.")
            seg.save_as(seg_path)

        if store_wsi_dicom:
            dcm_path = output_dir / f"{container_id}_im.dcm"
            wsi_dcm.save_as(dcm_path)

        # Store to bucket
        if output_bucket_obj is not None:
            seg_blob_name = f"{container_id}_seg.dcm"

            logging.info("Writing segmentation to output bucket.")
            cloud_io.write_dataset_to_blob(
                seg,
                output_bucket_obj,
                seg_blob_name,
            )

    logging.info("Printing errors")
    for msg in errors:
        logging.info(msg)


if __name__ == "__main__":
    run()

