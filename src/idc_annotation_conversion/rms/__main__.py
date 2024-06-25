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
import highdicom as hd
import numpy as np
import pandas as pd
from PIL import Image

from idc_annotation_conversion import cloud_config, cloud_io
from idc_annotation_conversion.rms.convert import (
    convert_xml_annotation,
    convert_segmentation,
)


# Bucket containing annotations for this project
ANNOTATION_BUCKET_PROJECT = "idc-external-031"
ANNOTATION_BUCKET = "rms_annotation_test_oct_2023"


@click.group()
def cli():
    pass


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


@cli.command()
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
def convert_xml_annotations(
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
    """Convert RMS XML annotations to DICOM SRs."""
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


@cli.command()
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path, file_okay=False),
    help="Output directory, default: no output directory",
)
@click.option(
    "--output-bucket",
    "-b",
    help="Name of output bucket, if any. Default: no output bucket.",
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
    "--create-pyramid/--no-create-pyramid",
    "-q/-Q",
    default=True,
    show_default=True,
    help="Create a full segmentation pyramid series.",
)
@click.option(
    "--segmentation-type",
    "-t",
    type=click.Choice(
        [v.name for v in hd.seg.SegmentationTypeValues],
        case_sensitive=False,
    ),
    default="FRACTIONAL",
    show_default=True,
    help="Segmentation type for the Segmentation Image, if any.",
)
@click.option(
    "--number",
    "-n",
    type=int,
    help="Number of cases to process. All by default.",
)
@click.option(
    "--workers",
    "-w",
    type=int,
    default=0,
    help="Numbers of workers to use for frame compression.",
)
@click.option(
    "--excluded-cases",
    "-e",
    multiple=True,
    help="Container IDs (before underscore) to skip. Multiple may be provided",
)
def convert_segmentations(
    output_dir: Optional[Path],
    output_bucket: Optional[str],
    number: Optional[int],
    segmentation_type: str,
    dimension_organization_type: str,
    create_pyramid: bool,
    workers: int,
    excluded_cases: Optional[List[str]] = None,
):
    """Convert RMS model segmentation masks to DICOM segmentations."""
    logging.basicConfig(level=logging.INFO)

    # Suppress highdicom logging (very talkative)
    logging.getLogger("highdicom.base").setLevel(logging.WARNING)
    logging.getLogger("highdicom.seg.sop").setLevel(logging.WARNING)

    if output_dir is not None:
        output_dir.mkdir(exist_ok=True)

    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    public_bucket = storage_client.bucket(cloud_config.DICOM_IMAGES_BUCKET)
    if output_bucket is not None:
        output_bucket_obj = output_client.bucket(output_bucket)

        if not output_bucket_obj.exists():
            output_bucket_obj.create(
                location=cloud_config.GCP_DEFAULT_LOCATION
            )
    else:
        output_bucket_obj = None

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
        container_prefix, *_ = (
            mask_blob.name
            .split("/")[-1]
            .replace(".npy", "")
            .split("_", maxsplit=1)
        )
        if excluded_cases is not None:
            if container_prefix in excluded_cases:
                logging.info(f"Skipping case {container_prefix}.")
                continue

        selection_df = find_series(container_prefix, bq_client)
        container_id = selection_df.ContainerIdentifier.iloc[0]

        if selection_df.crdc_series_uuid.nunique() != 1:
            raise RuntimeError("Found multiple source series")

        mask_bytes = mask_blob.download_as_bytes()
        mask = np.load(BytesIO(mask_bytes))

        # The segmentation mask matches the resolution of one level of the
        # source pyramid, but which one varies by case. Finding the
        # matching level and pull only this one and lower resolution
        # images
        for i, row in selection_df.iterrows():
            if mask.shape[:2] == (
                row.TotalPixelMatrixRows, row.TotalPixelMatrixColumns
            ):
                start_index = i
                break
        else:
            raise RuntimeError(
                f"No source image matching mask dimensions for {container_id}."
            )

        logging.info("Retrieving source images.")
        wsi_dcm = [
            cloud_io.read_dataset_from_blob(
                bucket=public_bucket,
                blob_name=f"{row.crdc_series_uuid}/{row.crdc_instance_uuid}.dcm",
            )
            for _, row in selection_df[start_index:].iterrows()
        ]

        segmentations = convert_segmentation(
            source_images=wsi_dcm,
            segmentation_array=mask,
            create_pyramid=create_pyramid,
            segmentation_type=segmentation_type,
            dimension_organization_type=dimension_organization_type,
            workers=workers,
        )

        # Store objects to filesystem
        if output_dir is not None:
            for i, seg in enumerate(segmentations):
                if create_pyramid:
                    seg_path = output_dir / f"{container_id}_seg_{i}.dcm"
                else:
                    seg_path = output_dir / f"{container_id}_seg.dcm"

                logging.info(f"Writing segmentation to {str(seg_path)}.")
                seg.save_as(seg_path)

        # Store to bucket
        if output_bucket_obj is not None:
            for i, seg in enumerate(segmentations):
                if create_pyramid:
                    seg_blob_name = f"{container_id}_seg_{i}.dcm"
                else:
                    seg_blob_name = f"{container_id}_seg.dcm"

                logging.info("Writing segmentation to output bucket.")
                cloud_io.write_dataset_to_blob(
                    seg,
                    output_bucket_obj,
                    seg_blob_name,
                )


if __name__ == "__main__":
    cli()
