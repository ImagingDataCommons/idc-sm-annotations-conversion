"""Main conversion process for PanOpTILs."""
import datetime
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import logging
from pathlib import Path
import os
from typing import Optional

import click
from google.cloud import bigquery
from google.cloud import storage
from pydicom.uid import (
    ExplicitVRLittleEndian,
    JPEG2000Lossless,
    JPEGLSLossless,
    RLELossless,
)
import highdicom as hd

from idc_annotation_conversion import cloud_io
from idc_annotation_conversion.panoptils.convert import (
    convert_segmentation,
)
from idc_annotation_conversion import cloud_config, cloud_io


# TODO fix reference coordinates in plane positions (currently fails highdicom check)
# TODO check why multiple images are matching
# TODO crop border cell tiles
# TODO include bootstrapped versions
# TODO metadata
# TODO include lut


SOURCE_BUCKET = "idc-panoptils-manual"
SOURCE_BUCKET_PROJECT = "idc-source-data"
MANUAL_PREFIX = "panoptils-manual/masks"
BOOTSTRAPPED_PREFIX = "panoptils-bootstrapped"


def run_png_blob(
    container_id: str,
    annotation_blobs: list[str],
    output_dir: Path | None = None,
    store_wsi_dicom: bool = False,
    output_bucket: str | None = None,
    segmentation_type: str = "BINARY",
    include_lut: bool = False,
    transfer_syntax_uid: str = ExplicitVRLittleEndian,
    crop_total_pixel_matrix: bool = False,
) -> str | None:
    """Convert a single h5 blob to an aggressiveness maps.

    Parameters
    ----------
    input_blob: str
        Path of storage blob contaning the aggressiveness map in h5 file
        format.
    output_dir: Path | None, optional
        Directory, if any, to store converted parametric map.
    store_wsi_dicom: bool, optional
        Whether to store the original source image after pulling it.
    output_bucket: str | None, optional
        Name of output bucket, if any, to store new parametric maps.
    transfer_syntax_uid: str
        Transfer syntax UID for the new segmentations.
    crop_total_pixel_matrix: bool
        If True, limit the size of the total pixel matrix of the output
        segmentation to the area covered by the data. If False, the total pixel
        matrix is defined to cover the same physical area covered by the source
        image. Note this does not affect the frames that are stored, just how
        the total pixel matrix rows and columns are defined, and the positional
        information of the frames.

    Returns
    -------
    str | None:
        An error/warning message, if any. Otherwise None.

    """
    error: str | None = None

    # Access bucket containing annotations
    ann_storage_client = storage.Client(project=SOURCE_BUCKET_PROJECT)
    ann_bucket = ann_storage_client.bucket(SOURCE_BUCKET)

    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    if output_bucket is not None:
        output_bucket_obj = output_client.bucket(output_bucket)
    else:
        output_bucket_obj = None

    bq_client = bigquery.Client(cloud_config.GCP_PROJECT_ID)
    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    im_bucket = storage_client.bucket(cloud_config.OPEN_DATA_BUCKET)

    logging.info(f"Running on container ID {container_id}")

    # The container IDs used in the file names are not "complete", they miss
    # out some parts, so need to use a partial match
    patient_id, slide_id = container_id.rsplit("-", maxsplit=1)

    selection_query = f"""
        SELECT
            SeriesDescription,
            SeriesInstanceUID,
            LossyImageCompression,
            ContainerIdentifier,
            gcs_url,
            CAST(NumberOfFrames AS int) AS NumberOfFrames
        FROM
            bigquery-public-data.idc_current.dicom_all
        WHERE
            ContainerIdentifier LIKE '{patient_id}-%-{slide_id}'
            AND SamplesPerPixel = 3
        ORDER BY
            NumberOfFrames DESC
    """
    selection_result = bq_client.query(selection_query)
    selection_df = selection_result.result().to_dataframe()
    if selection_df["ContainerIdentifier"].nunique() != 1:
        error = f"ERROR multiple matching container IDs: {container_id}"
        logging.error(error)
        return error

    if len(selection_df) == 0:
        error = f"ERROR cannot find image: {container_id}"
        logging.error(error)
        return error

    if len(selection_df) > 1:
        # Check that there aren't multiple options
        # TODO check why this seems to be happening a lot
        if (
            selection_df.iloc[0].NumberOfFrames ==
            selection_df.iloc[1].NumberOfFrames
        ):
            error = f"WARNING multiple images: {container_id}"
            logging.warning(error)

    url = selection_df.iloc[0].gcs_url
    im_blob_name = "/".join(url.split("/")[3:])

    wsi_ds = cloud_io.read_dataset_from_blob(
        im_bucket,
        im_blob_name,
        stop_before_pixels=not store_wsi_dicom,
    )
    wsi_im = hd.Image.from_dataset(wsi_ds, copy=False)
    logging.info("Pulled WSI image")

    # Read the data from the png files
    arrays = [
        cloud_io.read_image_from_blob(ann_bucket, ann_blob)
        for ann_blob in annotation_blobs
    ]
    coords = [
        (
            int(b.split("top-")[1].split("_")[0]),
            int(b.split("bottom-")[1].split("_")[0]),
            int(b.split("left-")[1].split("_")[0]),
            int(b.replace(".png", "").split("right-")[1].split("_")[0]),
        )
        for b in annotation_blobs
    ]

    region_seg, nuclei_seg, border_seg = convert_segmentation(
        arrays=arrays,
        coords=coords,
        source_image=wsi_im,
        segmentation_type=segmentation_type,
        include_lut=include_lut,
        container_id=selection_df.ContainerIdentifier.iloc[0],
        transfer_syntax_uid=transfer_syntax_uid,
        crop_total_pixel_matrix=crop_total_pixel_matrix,
    )

    # Store objects to filesystem
    if output_dir is not None:

        out_path = output_dir / f"{container_id}_nuclei_seg.dcm"
        logging.info(f"Writing segmentation to {str(out_path)}.")
        nuclei_seg.save_as(out_path)

        out_path = output_dir / f"{container_id}_border_seg.dcm"
        logging.info(f"Writing segmentation to {str(out_path)}.")
        border_seg.save_as(out_path)

        out_path = output_dir / f"{container_id}_region_seg.dcm"
        logging.info(f"Writing segmentation to {str(out_path)}.")
        region_seg.save_as(out_path)

    if store_wsi_dicom:
        dcm_path = output_dir / f"{container_id}_im.dcm"
        wsi_im.save_as(dcm_path)

    # Store to bucket
    if output_bucket_obj is not None:
        logging.info("Writing segmentations to output bucket.")
        cloud_io.write_dataset_to_blob(
            border_seg,
            output_bucket_obj,
            f"{container_id}_border_seg.dcm",
        )
        cloud_io.write_dataset_to_blob(
            nuclei_seg,
            output_bucket_obj,
            f"{container_id}_nuclei_seg.dcm",
        )
        cloud_io.write_dataset_to_blob(
            region_seg,
            output_bucket_obj,
            f"{container_id}_region_seg.dcm",
        )

    return error


@click.command()
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
    help="Segmentation type for the Segmentation Image",
)
@click.option(
    "--transfer-syntax-uid",
    "-x",
    type=click.Choice(
        [
            ExplicitVRLittleEndian,
            JPEG2000Lossless,
            JPEGLSLossless,
            RLELossless,
        ],
        case_sensitive=False,
    ),
    default=ExplicitVRLittleEndian,
    show_default=True,
    help="Transfer syntax for the segmentations",
)
@click.option(
    "--include-lut/--no-include-lut",
    "-l/-L",
    default=True,
    show_default=True,
    help="Include a LUT in a labelmap seg.",
)
@click.option(
    "--crop-total-pixel-matrix/--no-crop-total-pixel-matrix",
    "-c/-C",
    default=False,
    show_default=True,
    help=(
        "Whether to limit the total pixel matrix of the segmentations to "
        "the region containing data."
    ),
)
def main(
    number: Optional[int],
    output_dir: Optional[Path],
    output_bucket: str | None,
    store_bucket: bool,
    store_wsi_dicom: bool,
    segmentation_type: str,
    keep_existing: bool = False,
    include_lut: bool = True,
    workers: int = 0,
    transfer_syntax_uid: str = ExplicitVRLittleEndian,
    crop_total_pixel_matrix: bool = False,
):
    """Convert PNG files to DICOM segmentations.

    Dataset is originally downloaded from:

    https://sites.google.com/view/panoptils/

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

    # Access bucket containing annotations
    ann_storage_client = storage.Client(project=SOURCE_BUCKET_PROJECT)
    ann_bucket = ann_storage_client.bucket(SOURCE_BUCKET)

    # Create output directory
    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    output_bucket_obj = None
    if store_bucket:
        if output_bucket is None:
            date_str = (datetime.date.today())
            output_bucket = (
                f"panoptils-segmentations-{date_str}"
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
    to_process = defaultdict(list)

    for ann_blob in ann_bucket.list_blobs(prefix=MANUAL_PREFIX):
        if not ann_blob.name.endswith(".png"):
            continue

        blob_name = ann_blob.name.split("/")[-1].replace(".png", "")
        container_id = blob_name.split("_")[0]

        # Check whether the output blobs already exist, and skip if
        # they do
        if keep_existing:
            seg_blob_name = f"{container_id}_region_seg.dcm"
            seg_output_blob = output_bucket_obj.get_blob(seg_blob_name)
            if seg_output_blob.exists():
                continue

        if number is not None:
            if (
                container_id not in to_process
                and len(to_process) == number
            ):
                break

        to_process[container_id].append(ann_blob.name)

    errors = []

    if workers == 0:
        errors = [
            run_png_blob(
                container_id=container_id,
                annotation_blobs=blobs,
                segmentation_type=segmentation_type,
                include_lut=include_lut,
                output_dir=output_dir,
                store_wsi_dicom=store_wsi_dicom,
                output_bucket=output_bucket,
                transfer_syntax_uid=transfer_syntax_uid,
                crop_total_pixel_matrix=crop_total_pixel_matrix,
            )
            for container_id, blobs in to_process.items()
        ]
    else:
        with ProcessPoolExecutor(workers) as pool:
            futures = []
            for container_id, blobs in to_process.items():
                fut = pool.submit(
                    run_png_blob,
                    container_id=container_id,
                    annotation_blobs=blobs,
                    segmentation_type=segmentation_type,
                    include_lut=include_lut,
                    output_dir=output_dir,
                    store_wsi_dicom=store_wsi_dicom,
                    output_bucket=output_bucket,
                    crop_total_pixel_matrix=crop_total_pixel_matrix,
                    transfer_syntax_uid=transfer_syntax_uid,
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
    main()
