"""Main conversion process for GBP transcriptional subtypes."""
import datetime
from concurrent.futures import ProcessPoolExecutor
import logging
from pathlib import Path
import os
from typing import Optional

import click
from google.cloud import bigquery
from google.cloud import storage
import highdicom as hd
import numpy as np
import pandas as pd

from idc_annotation_conversion import cloud_io
from idc_annotation_conversion.gbm_transcriptional_subtypes.convert import (
    convert_segmentation,
    convert_aggressiveness_map,
)
from idc_annotation_conversion import cloud_config
import h5py


SOURCE_BUCKET = "gevaert-cell-subtyping"
AGGRESSIVENESS_MAP_PREFIX = "aggressiveness-maps"


@click.group()
def cli():
    pass


def run_subtype_map_blob(
    slide_df: pd.DataFrame,
    dimension_organization_type: str,
    segmentation_type: str,
    include_lut: bool = True,
    output_dir: Path | None = None,
    store_wsi_dicom: bool = False,
    output_bucket: str | None = None,
) -> str | None:
    """Convert a single PNG blob for the 2018 TIL Maps.

    Parameters
    ----------
    slide_df: pandas.DataFrame,
        Dataframe containing values for this slide.
    dimension_organization_type: str
        Dimension organization of the output segmentation.
    segmentation_type: str
        Segmentation type of the output segmentation.
    include_lut: bool, optional
        Whether to include a palette color LUT in the output segmentation.
        Ignored if segmentation type is not "LABELMAP".
    output_dir: Path | None, optional
        Directory, if any, to store converted segmentations.
    store_wsi_dicom: bool, optional
        Whether to store the original source image after pulling it.
    output_bucket: str | None, optional
        Name of output bucket, if any, to store new segmentations.

    Returns
    -------
    str | None:
        An error/warning message, if any. Otherwise None.

    """
    error: str | None = None
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)

    if output_bucket is not None:
        output_bucket_obj = output_client.bucket(output_bucket)
    else:
        output_bucket_obj = None

    bq_client = bigquery.Client(cloud_config.GCP_PROJECT_ID)
    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    im_bucket = storage_client.bucket(cloud_config.OPEN_DATA_BUCKET)

    container_id = slide_df.slide_id.iloc[0].split('.')[0]
    logging.info(f"Running on container ID {container_id}")

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
        error = f"ERROR cannot find image: {container_id}"
        logging.error(error)
        return error

    if len(selection_df) > 1:
        # Check that there aren't multiple options
        if (
            selection_df.iloc[0].NumberOfFrames ==
            selection_df.iloc[1].NumberOfFrames
        ):
            error = f"WARNING multiple images: {container_id}"
            logging.warning(error)

    url = selection_df.iloc[0].gcs_url
    im_blob_name = "/".join(url.split("/")[3:])

    wsi_dcm = cloud_io.read_dataset_from_blob(
        im_bucket,
        im_blob_name,
        stop_before_pixels=not store_wsi_dicom,
    )
    logging.info("Pulled WSI image")

    seg = convert_segmentation(
        slide_df,
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

    return error


def run_aggressiveness_map_blob(
    input_blob: str,
    dimension_organization_type: str,
    output_dir: Path | None = None,
    store_wsi_dicom: bool = False,
    output_bucket: str | None = None,
) -> str | None:
    """Convert a single h5 blob to an aggressiveness maps.

    Parameters
    ----------
    input_blob: str
        Path of storage blob contaning the aggressiveness map in h5 file
        format.
    dimension_organization_type: str
        Dimension organization of the output parametric map.
    output_dir: Path | None, optional
        Directory, if any, to store converted parametric map.
    store_wsi_dicom: bool, optional
        Whether to store the original source image after pulling it.
    output_bucket: str | None, optional
        Name of output bucket, if any, to store new parametric maps.

    Returns
    -------
    str | None:
        An error/warning message, if any. Otherwise None.

    """
    error: str | None = None
    input_client = storage.Client(project=cloud_config.SOURCE_DATA_PROJECT)
    input_bucket = input_client.bucket(SOURCE_BUCKET)
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)

    if output_bucket is not None:
        output_bucket_obj = output_client.bucket(output_bucket)
    else:
        output_bucket_obj = None

    bq_client = bigquery.Client(cloud_config.GCP_PROJECT_ID)
    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    im_bucket = storage_client.bucket(cloud_config.OPEN_DATA_BUCKET)

    container_id = input_blob.split("/")[1].split(".")[0]
    logging.info(f"Running on container ID {container_id}")

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
        error = f"ERROR cannot find image: {container_id}"
        logging.error(error)
        return error

    if len(selection_df) > 1:
        # Check that there aren't multiple options
        if (
            selection_df.iloc[0].NumberOfFrames ==
            selection_df.iloc[1].NumberOfFrames
        ):
            error = f"WARNING multiple images: {container_id}"
            logging.warning(error)

    url = selection_df.iloc[0].gcs_url
    im_blob_name = "/".join(url.split("/")[3:])

    wsi_dcm = cloud_io.read_dataset_from_blob(
        im_bucket,
        im_blob_name,
        stop_before_pixels=not store_wsi_dicom,
    )
    logging.info("Pulled WSI image")

    # Read the data from the h5 file
    blob_object = input_bucket.get_blob(input_blob)
    with h5py.File(blob_object.open("rb")) as f:
        scores = np.array(f['scores'])
        coords = np.array(
            [(int(s.split('_')[0]), int(s.split('_')[1])) for s in f['coords'].asstr()]
        )

    pmap = convert_aggressiveness_map(
        scores=scores,
        coords=coords,
        source_image=wsi_dcm,
        dimension_organization_type=dimension_organization_type,
    )

    # Store objects to filesystem
    if output_dir is not None:
        out_path = output_dir / f"{container_id}_aggressiveness_pmap.dcm"

        logging.info(f"Writing parametric map to {str(out_path)}.")
        pmap.save_as(out_path)

    if store_wsi_dicom:
        dcm_path = output_dir / f"{container_id}_im.dcm"
        wsi_dcm.save_as(dcm_path)

    # Store to bucket
    if output_bucket_obj is not None:
        out_blob_name = f"{container_id}_aggressiveness_pmap.dcm"

        logging.info("Writing parametric map to output bucket.")
        cloud_io.write_dataset_to_blob(
            pmap,
            output_bucket_obj,
            out_blob_name,
        )

    return error


@cli.command()
@click.option(
    "--number",
    "-n",
    type=int,
    help="Number of annotations to process. All by default.",
)
@click.option(
    "--workers",
    "-w",
    type=int,
    default=0,
    help="Number of worker processes.",
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
    help="Segmentation type for the Segmentation Image",
)
@click.option(
    "--include-lut/--no-include-lut",
    "-l/-L",
    default=True,
    show_default=True,
    help="Include a LUT in a labelmap seg.",
)
def convert_subtype_maps(
    number: Optional[int],
    output_dir: Optional[Path],
    output_bucket: str | None,
    store_bucket: bool,
    store_wsi_dicom: bool,
    segmentation_type: str,
    dimension_organization_type: str,
    keep_existing: bool = False,
    include_lut: bool = True,
    workers: int = 0,
):
    """Convert transcriptional subtype maps to DICOM segmentations.

    The analysis that produced these results is described in:

    Zheng, Yuanning, et al. "Spatial cellular architecture predicts prognosis
    in glioblastoma." Nature Communications 14.1 (2023): 4122.

    By default, output is to a cloud bucket. It is also possible to output
    to a local directory.

    """
    logging.basicConfig(level=logging.INFO)

    # Suppress highdicom logging (very talkative)
    logging.getLogger("highdicom.base").setLevel(logging.WARNING)
    logging.getLogger("highdicom.seg.sop").setLevel(logging.WARNING)

    # Setup project and authenticate
    os.environ["GCP_PROJECT_ID"] = cloud_config.GCP_PROJECT_ID

    if keep_existing and not store_bucket:
        raise ValueError("keep_existing requires store_bucket")

    # The annotations are loaded in from file
    full_dataframe = pd.read_feather("test_cases_table.feather")

    # Access bucket containing annotations
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    output_bucket_obj = None
    if store_bucket:
        if output_bucket is None:
            date_str = (datetime.date.today())
            output_bucket = (
                f"tcga-gbm-transcriptional-subtype-maps-{date_str}"
            )
        output_bucket_obj = output_client.bucket(output_bucket)

        if not output_bucket_obj.exists():
            output_bucket_obj.create(
                location=cloud_config.GCP_DEFAULT_LOCATION
            )
    else:
        output_bucket = None

    # Create output directory
    if output_dir is not None:
        output_dir.mkdir(exist_ok=True)

    logging.info("Listing cases")

    to_process = []

    unique_slides = np.unique(full_dataframe.slide_id.values)

    # Loop over requested collections
    for slide_id in unique_slides:
        # Check whether the output blobs already exist, and skip if
        # they do
        if keep_existing:
            seg_blob_name = f"{slide_id}_seg.dcm"
            seg_output_blob = output_bucket_obj.get_blob(seg_blob_name)
            if seg_output_blob.exists():
                continue

        slide_df = full_dataframe[full_dataframe.slide_id == slide_id].copy()

        to_process.append(slide_df)

    if number is not None:
        to_process = to_process[:number]

    errors = []

    if workers == 0:
        errors = [
            run_subtype_map_blob(
                slide_df=slide_df,
                dimension_organization_type=dimension_organization_type,
                segmentation_type=segmentation_type,
                include_lut=include_lut,
                output_dir=output_dir,
                store_wsi_dicom=store_wsi_dicom,
                output_bucket=output_bucket,
            )
            for slide_df in to_process
        ]
    else:
        with ProcessPoolExecutor(workers) as pool:
            futures = []
            for slide_df in to_process:
                fut = pool.submit(
                    run_subtype_map_blob,
                    slide_df=slide_df,
                    dimension_organization_type=dimension_organization_type,
                    segmentation_type=segmentation_type,
                    include_lut=include_lut,
                    output_dir=output_dir,
                    store_wsi_dicom=store_wsi_dicom,
                    output_bucket=output_bucket,
                )
                futures.append(fut)

            errors = [fut.result() for fut in futures]

    if len(errors) > 0:
        logging.info("Printing errors")
        for msg in errors:
            if msg is not None:
                logging.info(msg)
    else:
        logging.info("No errors!")


@cli.command()
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
    "--workers",
    "-w",
    type=int,
    default=0,
    help="Number of worker processes.",
)
def convert_aggressiveness_maps(
    number: Optional[int],
    output_dir: Optional[Path],
    output_bucket: str | None,
    store_bucket: bool,
    store_wsi_dicom: bool,
    dimension_organization_type: str,
    keep_existing: bool = False,
    workers: int = 0,
):
    """Convert aggressiveness maps to DICOM Parametric Maps.

    The analysis that produced these results is described in:

    Zheng, Yuanning, et al. "Spatial cellular architecture predicts prognosis
    in glioblastoma." Nature Communications 14.1 (2023): 4122.

    By default, output is to a cloud bucket. It is also possible to output
    to a local directory.

    """
    logging.basicConfig(level=logging.INFO)

    # Suppress highdicom logging (very talkative)
    logging.getLogger("highdicom.base").setLevel(logging.WARNING)

    # Setup project and authenticate
    os.environ["GCP_PROJECT_ID"] = cloud_config.GCP_PROJECT_ID

    if keep_existing and not store_bucket:
        raise ValueError("keep_existing requires store_bucket")

    input_client = storage.Client(project=cloud_config.SOURCE_DATA_PROJECT)
    input_bucket = input_client.bucket(SOURCE_BUCKET)

    # Access bucket containing annotations
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    output_bucket_obj = None
    if store_bucket:
        if output_bucket is None:
            date_str = (datetime.date.today())
            output_bucket = (
                f"tcga-gbm-aggressiveness-maps-{date_str}"
            )
        output_bucket_obj = output_client.bucket(output_bucket)

        if not output_bucket_obj.exists():
            output_bucket_obj.create(
                location=cloud_config.GCP_DEFAULT_LOCATION
            )
    else:
        output_bucket = None

    # Create output directory
    if output_dir is not None:
        output_dir.mkdir(exist_ok=True)

    logging.info("Listing cases")
    h5_list = input_bucket.list_blobs(prefix=AGGRESSIVENESS_MAP_PREFIX)

    to_process = []

    # Loop over requested collections
    for input_blob in h5_list:
        if not input_blob.name.endswith("h5"):
            continue

        container_id = input_blob.name.split("/")[1].split(".")[0]

        # Check whether the output blobs already exist, and skip if
        # they do
        if keep_existing:
            pmap_blob_name = f"{container_id}_aggressiveness_pmap.dcm"
            pmap_output_blob = output_bucket_obj.get_blob(pmap_blob_name)
            if pmap_output_blob.exists():
                continue

        to_process.append(input_blob)

    if number is not None:
        to_process = to_process[:number]

    errors = []

    if workers == 0:
        errors = [
            run_aggressiveness_map_blob(
                input_blob=input_blob.name,
                dimension_organization_type=dimension_organization_type,
                output_dir=output_dir,
                store_wsi_dicom=store_wsi_dicom,
                output_bucket=output_bucket,
            )
            for input_blob in to_process
        ]
    else:
        with ProcessPoolExecutor(workers) as pool:
            futures = []
            for input_blob in to_process:
                fut = pool.submit(
                    run_aggressiveness_map_blob,
                    input_blob=input_blob.name,
                    dimension_organization_type=dimension_organization_type,
                    output_dir=output_dir,
                    store_wsi_dicom=store_wsi_dicom,
                    output_bucket=output_bucket,
                )
                futures.append(fut)

            errors = [fut.result() for fut in futures]

    if len(errors) > 0:
        logging.info("Printing errors")
        for msg in errors:
            if msg is not None:
                logging.info(msg)
    else:
        logging.info("No errors!")


if __name__ == "__main__":
    cli()
