"""Main conversion process for the TCGA TIL Maps project."""
import datetime
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
import logging
from pathlib import Path
import os
from typing import Optional

import click
from google.cloud import bigquery
from google.cloud import storage
import highdicom as hd

from idc_annotation_conversion import cloud_io
from idc_annotation_conversion.tcga_til_maps.convert import (
    convert_segmentation,
    convert_txt_file,
)
from idc_annotation_conversion import cloud_config


IMAGE_BUCKET = "idc-open-data"
ANNOTATION_BUCKET_2018 = "til-wsi-tcga"
ANNOTATION_PREFIX_2018 = 'TIL_maps_after_thres_v1'
ANNOTATION_BUCKET_2022 = "til-wsi-tcga-nature-new-results"
ANNOTATION_PREFIX_2022 = 'TCGA_TIL_Maps'

COLLECTIONS_2018 = [
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

COLLECTIONS_2022 = [
    "acc",
    "blca",
    "brca",
    "cesc",
    "coad",
    "esca",
    "hnsc",
    "kirc",
    "lihc",
    "luad",
    "lusc",
    "meso",
    "ov",
    "paad",
    "prad",
    "read",
    "sarc",
    "skcm",
    "stad",
    "tgct",
    "thym",
    "ucec",
    "uvm",
]


@click.group()
def cli():
    pass


def run_blob_2018(
    blob_name: str,
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
    blob_name: str
        Name of the blob in the annotation bucket containing the PNG.
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
    im_bucket = storage_client.bucket(IMAGE_BUCKET)
    ann_bucket = storage_client.bucket(ANNOTATION_BUCKET_2018)

    container_id = blob_name.split('/')[-1].replace('.png', '')
    logging.info(f"Running on container ID {container_id}")

    mask = cloud_io.read_image_from_blob(
        ann_bucket,
        blob_name
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

    return error


@cli.command()
@click.option(
    "-c",
    "--collections",
    multiple=True,
    type=click.Choice(COLLECTIONS_2018),
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
def convert_2018(
    collections: Optional[list[str]],
    number: Optional[int],
    output_dir: Optional[Path],
    output_bucket: str | None,
    store_bucket: bool,
    store_wsi_dicom: bool,
    segmentation_type: str,
    dimension_organization_type: str,
    blob_filter: Optional[str] = None,
    keep_existing: bool = False,
    include_lut: bool = True,
    workers: int = 0,
):
    """Convert TCGA tumor infiltrating lymphocyte (TIL) maps to DICOM
    segmentations.

    This routine is for the earlier set of TIL maps described here:

    https://www.cancerimagingarchive.net/analysis-result/til-wsi-tcga/

    relating to this publication:

    Saltz, Joel, et al. "Spatial organization and molecular correlation of
    tumor-infiltrating lymphocytes using deep learning on pathology images."
    Cell reports 23.1 (2018): 181-193.

    Images and annotations are automatically pulled down from cloud buckets as
    required.

    By default, output is to a cloud bucket. It is also possible to output
    to a local directory.

    """
    # Use all collections if none specified
    collections = collections or COLLECTIONS_2018

    logging.basicConfig(level=logging.INFO)

    # Suppress highdicom logging (very talkative)
    logging.getLogger("highdicom.base").setLevel(logging.WARNING)
    logging.getLogger("highdicom.seg.sop").setLevel(logging.WARNING)

    # Setup project and authenticate
    os.environ["GCP_PROJECT_ID"] = cloud_config.GCP_PROJECT_ID

    if keep_existing and not store_bucket:
        raise ValueError("keep_existing requires store_bucket")

    # Access bucket containing annotations
    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    ann_bucket = storage_client.bucket(ANNOTATION_BUCKET_2018)
    output_bucket_obj = None
    if store_bucket:
        if output_bucket is None:
            date_str = (datetime.date.today())
            output_bucket = (
                f"tcga_til_maps_older_paper_{date_str}"
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

    logging.info("Listing PNGs")

    to_process = []

    # Loop over requested collections
    for collection in collections:
        collection_prefix = f'{ANNOTATION_PREFIX_2018}/{collection}/'

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

    if workers == 0:
        errors = [
            run_blob_2018(
                blob_name=blob.name,
                dimension_organization_type=dimension_organization_type,
                segmentation_type=segmentation_type,
                include_lut=include_lut,
                output_dir=output_dir,
                store_wsi_dicom=store_wsi_dicom,
                output_bucket=output_bucket,
            )
            for blob in to_process
        ]
    else:
        with ProcessPoolExecutor(workers) as pool:
            futures = []
            for blob in to_process:
                fut = pool.submit(
                    run_blob_2018,
                    blob_name=blob.name,
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


def run_blob_2022(
    container_id: str,
    binary_blob_name: str | None,
    fractional_blob_name: str | None,
    segmentation_type: str,
    output_dir: Path | None = None,
    store_wsi_dicom: bool = False,
    output_bucket: str | None = None,
    dimension_organization_type: str = "TILED_FULL",
) -> str | None:
    """Convert segmentations for a single source image for the 2022 dataset.

    Parameters
    ----------
    container_id: str
        Container ID of the source image.
    binary_blob_name: str | None
        Name of the blob in the annotation bucket containing the binary TIL map
        (as a text file). May be omitted.
    fractional_blob_name: str | None
        Name of the blob in the annotation bucket containing the probabilitic
        TIL map (as a text file). May be omitted.
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
    dimension_organization_type: str
        Dimension organization of the output segmentation.

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
    im_bucket = storage_client.bucket(IMAGE_BUCKET)
    ann_bucket = storage_client.bucket(ANNOTATION_BUCKET_2022)

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
    wsi_dcm = hd.Image.from_dataset(wsi_dcm, copy=False)
    logging.info("Pulled WSI image")

    if store_wsi_dicom and output_dir is not None:
        dcm_path = output_dir / f"{container_id}_im.dcm"
        wsi_dcm.save_as(dcm_path)

    if binary_blob_name is not None:
        binary_blob = ann_bucket.get_blob(binary_blob_name)

        with binary_blob.open("r") as tf:
            binary_seg = convert_txt_file(
                tf,
                container_id=container_id,
                segmentation_type=segmentation_type,
                source_image=wsi_dcm,
                is_probability=False,
                dimension_organization_type=dimension_organization_type,
            )

        # Store objects to filesystem
        if output_dir is not None:
            seg_path = output_dir / f"{container_id}_binary_seg.dcm"

            logging.info(f"Writing binary segmentation to {str(seg_path)}.")
            binary_seg.save_as(seg_path)

        # Store to bucket
        if output_bucket_obj is not None:
            seg_blob_name = f"{container_id}_binary_seg.dcm"

            logging.info("Writing binary segmentation to output bucket.")
            cloud_io.write_dataset_to_blob(
                binary_seg,
                output_bucket_obj,
                seg_blob_name,
            )

    if fractional_blob_name is not None:
        fractional_blob = ann_bucket.get_blob(fractional_blob_name)

        with fractional_blob.open("r") as tf:
            fractional_seg = convert_txt_file(
                tf,
                container_id=container_id,
                segmentation_type=hd.seg.SegmentationTypeValues.FRACTIONAL,
                source_image=wsi_dcm,
                is_probability=True,
                dimension_organization_type=dimension_organization_type,
            )

        # Store objects to filesystem
        if output_dir is not None:
            seg_path = output_dir / f"{container_id}_fractional_seg.dcm"

            logging.info(f"Writing fractional segmentation to {str(seg_path)}.")
            fractional_seg.save_as(seg_path)

        # Store to bucket
        if output_bucket_obj is not None:
            seg_blob_name = f"{container_id}_fractional_seg.dcm"

            logging.info("Writing fractional segmentation to output bucket.")
            cloud_io.write_dataset_to_blob(
                fractional_seg,
                output_bucket_obj,
                seg_blob_name,
            )

    return error


@cli.command()
@click.option(
    "-c",
    "--collections",
    multiple=True,
    type=click.Choice(COLLECTIONS_2022),
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
    "--segmentation-type",
    "-t",
    type=click.Choice(
        [v.name for v in hd.seg.SegmentationTypeValues],
        case_sensitive=False,
    ),
    default="BINARY",
    show_default=True,
    help="Segmentation type for the binary segmentations.",
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
    "--include-lut/--no-include-lut",
    "-l/-L",
    default=True,
    show_default=True,
    help="Include a LUT in a labelmap seg.",
)
def convert_2022(
    collections: Optional[list[str]],
    number: Optional[int],
    output_dir: Optional[Path],
    output_bucket: str | None,
    store_bucket: bool,
    store_wsi_dicom: bool,
    segmentation_type: str,
    dimension_organization_type: str,
    blob_filter: Optional[str] = None,
    keep_existing: bool = False,
    include_lut: bool = True,
    workers: int = 0,
):
    """Convert TIL maps from the later 2022 paper:

    Abousamra, S., Gupta, R., Hou, L., Batiste, R., Zhao, T., Shankar, A., Saltz,
    J. (2022). Deep learning-based mapping of tumor infiltrating lymphocytes in
    whole slide images of 23 types of cancer. Frontiers in oncology, 11, 806603.

    to DICOM segmentations.

    Images and annotations are automatically pulled down from cloud buckets as
    required. There are binary and "heatmap" versions of each map. Both are
    converted to segmentations.

    By default, output is to a cloud bucket. It is also possible to output
    to a local directory.

    """
    # Use all collections if none specified
    collections = collections or COLLECTIONS_2022

    logging.basicConfig(level=logging.INFO)

    # Suppress highdicom logging (very talkative)
    logging.getLogger("highdicom.base").setLevel(logging.WARNING)
    logging.getLogger("highdicom.seg.sop").setLevel(logging.WARNING)

    # Setup project and authenticate
    os.environ["GCP_PROJECT_ID"] = cloud_config.GCP_PROJECT_ID

    if keep_existing and not store_bucket:
        raise ValueError("keep_existing requires store_bucket")

    # Access bucket containing annotations
    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    ann_bucket = storage_client.bucket(ANNOTATION_BUCKET_2022)
    output_bucket_obj = None
    if store_bucket:
        if output_bucket is None:
            date_str = (datetime.date.today())
            output_bucket = (
                f"tcga_til_maps_newer_paper_{date_str}"
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

    logging.info("Listing input files")

    # Loop over requested collections
    to_process_pairs = defaultdict(dict)

    for collection in collections:
        collection_prefix = f'{ANNOTATION_PREFIX_2022}/{collection}/'

        for b in ann_bucket.list_blobs(prefix=collection_prefix):
            if not 'prediction' in b.name or b.name.endswith('low_res'):
                continue

            if blob_filter is not None:
                if blob_filter not in b.name:
                    continue

            dir_name = b.name.split('/')[-2]
            is_binary = 'binary' in dir_name

            container_id = b.name.split('/')[-1].replace('prediction-', '').split('.')[0]

            # Check whether the output blobs already exist, and skip if
            # they do
            if keep_existing:
                if is_binary:
                    seg_blob_name = f"{container_id}_binary_seg.dcm"
                else:
                    seg_blob_name = f"{container_id}_fractional_seg.dcm"

                seg_output_blob = output_bucket_obj.get_blob(seg_blob_name)
                if seg_output_blob.exists():
                    continue

            if is_binary:
                to_process_pairs[container_id]['binary'] = b
            else:
                to_process_pairs[container_id]['fractional'] = b

    # Tuple of (container_id, binary_blob, fractional_blob)
    to_process = [
        (k, v.get('binary'), v.get('fractional')) for k, v in to_process_pairs.items()
    ]

    if number is not None:
        to_process = to_process[:number]

    errors = []

    if workers == 0:
        errors = [
            run_blob_2022(
                container_id=container_id,
                binary_blob_name=binary_blob.name if binary_blob is not None else None,
                fractional_blob_name=fractional_blob.name if fractional_blob is not None else None,
                segmentation_type=segmentation_type,
                output_dir=output_dir,
                store_wsi_dicom=store_wsi_dicom,
                output_bucket=output_bucket,
                dimension_organization_type=dimension_organization_type,
            )
            for (container_id, binary_blob, fractional_blob) in to_process
        ]
    else:
        with ProcessPoolExecutor(workers) as pool:
            futures = []
            for (container_id, binary_blob, fractional_blob) in to_process:
                fut = pool.submit(
                    run_blob_2022,
                    container_id=container_id,
                    binary_blob_name=binary_blob.name if binary_blob is not None else None,
                    fractional_blob_name=fractional_blob.name if fractional_blob is not None else None,
                    segmentation_type=segmentation_type,
                    output_dir=output_dir,
                    store_wsi_dicom=store_wsi_dicom,
                    output_bucket=output_bucket,
                    dimension_organization_type=dimension_organization_type,
                )
                futures.append(fut)

            errors = [fut.result() for fut in futures]

    logging.info("Printing errors")
    for msg in errors:
        if msg is not None:
            logging.info(msg)


if __name__ == "__main__":
    cli()
