from io import BytesIO, BufferedReader
from itertools import islice
from getpass import getpass
import logging
from pathlib import Path
import os
import tarfile
from time import time
from typing import List, Generator, Optional

import click
from dicomweb_client import DICOMwebClient
from google.cloud import bigquery
from google.cloud import storage
import highdicom as hd
from oauthlib.oauth2 import BackendApplicationClient
from requests_oauthlib import OAuth2Session

from idc_annotation_conversion import cloud_config, cloud_io
from idc_annotation_conversion.convert import convert_annotations


COLLECTIONS = [
    "blca_polygon",
    "brca_polygon",
    "cesc_polygon",
    "coad_polygon",
    "gbm_polygon",
    "luad_polygon",
    "lusc_polygon",
    "paad_polygon",
    "prad_polygon",
    "read_polygon",
    "skcm_polygon",
    "stad_polygon",
    "ucec_polygon",
    "uvm_polygon",
]


def iter_csvs(ann_blob: storage.Blob) -> Generator[BufferedReader, None, None]:
    """Iterate over individual CSV files in a loaded tar.gz file.

    Parameters:
    -----------
    ann_bytes: google.cloud.storage.Blob
        Blob object for the annotation tar.gz file.

    Yields:
    -------
    io.BufferedReader
        Buffer for each CSV file in the original tar. This buffer can be
        used to load the contents of the CSV from the in-memory tar file.

    """
    # Download the annotation
    ann_bytes = ann_blob.download_as_bytes()

    # Untar in memory and yield each csv as a buffer
    with tarfile.open(fileobj=BytesIO(ann_bytes), mode='r:gz') as tar:
        for member in tar.getmembers():
            if member.isfile():
                yield tar.extractfile(member)


def get_dicom_web_client(
    url: str,
    token_url: Optional[str] = None,
    client_id: Optional[str] = None,
    client_secret: Optional[str] = None,
) -> DICOMwebClient:
    """Create a DICOM Web Client.

    Parameters
    ----------
    url: str
        URL of the archive.
    token_url: Optional[Str], optional
        URL from which to request if using OAuth. If not provided, no
        authentication will be used.
    client_id: Optional[str]
        Client ID to use when requesting OAuth token. Required if OAuth
        authentication is to be used.
    client_secret: Optional[str]
        Client secret to use when requesting OAuth token. Required if OAuth
        authentication is to be used.

    Returns
    -------
    DICOMwebClient
        Client created using the provided parameters.

    """
    if token_url is None:
        # Simple, no authentication
        return DICOMwebClient(url)

    oauth_client = BackendApplicationClient(client_id=client_id)
    oauth = OAuth2Session(client=oauth_client)
    oauth.fetch_token(
        token_url=token_url,
        client_id=client_id,
        client_secret=client_secret,
    )

    return DICOMwebClient(url, session=oauth)


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
    "--number",
    "-n",
    type=int,
    help="Number to process per collection. All by default.",
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
    default=cloud_config.DEFAULT_OUTPUT_BUCKET,
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
    "--store-boundary/--store-centroid",
    "-B/-C",
    default=True,
    show_default=True,
    help=(
        "Store either the full boundary of each nucleus as a polygon "
        "(default), or the just the centroid as a single point."
    ),
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
    "--annotation-coordinate-type",
    "-a",
    type=click.Choice(
        [v.name for v in hd.ann.AnnotationCoordinateTypeValues],
        case_sensitive=False,
    ),
    default="SCOORD",
    show_default=True,
    help=(
        "Coordinate type for points stored in the microscopy annotations. "
        "SCOORD: 2D, SCOORD3D: 3D."
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
    "--with-segmentation/--without-segmentation",
    "-s/-S",
    default=True,
    show_default=True,
    help="Include a segmentation image in the output.",
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
    default="BINARY",
    show_default=True,
    help="Segmentation type for the Segmentation Image, if any.",
)
@click.option(
    "--dicom-archive",
    "-x",
    help="Additionally store outputs to this DICOM archive.",
)
@click.option(
    "--archive-token-url",
    "-u",
    help="URL to use to request an OAuth token to access the archive.",
)
@click.option(
    "--archive-client-id",
    "-i",
    help="Client ID to use for OAuth token request.",
)
@click.option(
    "--archive-client-secret",
    "-y",
    help=(
        "Client secret to use for OAuth token request. If none, user will "
        "be prompted for secret."
    )
)
@click.option(
    "--workers",
    "-w",
    type=int,
    default=0,
    help="Number of subprocesses to use. If 0, the main thread is used."
)
def run(
    collections: Optional[List[str]],
    number: Optional[int],
    output_dir: Optional[Path],
    output_bucket: str,
    output_prefix: Optional[str],
    store_bucket: bool,
    store_boundary: bool,
    store_wsi_dicom: bool,
    annotation_coordinate_type: str,
    with_segmentation: bool,
    segmentation_type: str,
    dimension_organization_type: str,
    create_pyramid: bool,
    dicom_archive: Optional[str] = None,
    archive_token_url: Optional[str] = None,
    archive_client_id: Optional[str] = None,
    archive_client_secret: Optional[str] = None,
    workers: int = 0,
):
    """Convert TCGA cell nuclei annotations to DICOM format.

    Convert CSV-format annotations of cell nuclei to DICOM format. Images and
    annotations are automatically pulled down from cloud buckets as required.

    By default, all annotations for all collections are processed, which will
    take a very long time. This can be controlled with options.

    Bulk Microscopy Simple Annotations are always produced, Segmentation Images
    may optionally be created.

    By default, output is to a cloud bucket. It is also possible to output
    to a local directory.

    """
    # Use all collections if none specified
    collections = collections or COLLECTIONS
    output_prefix = output_prefix or ""

    logging.basicConfig(level=logging.INFO)

    # Suppress highdicom logging (very talkative)
    logging.getLogger("highdicom.base").setLevel(logging.WARNING)
    logging.getLogger("highdicom.seg.sop").setLevel(logging.WARNING)

    # Setup project and authenticate
    os.environ["GCP_PROJECT_ID"] = cloud_config.GCP_PROJECT_ID

    # Access bucket containing annotations
    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    ann_bucket = storage_client.bucket(cloud_config.ANNOTATION_BUCKET)
    public_bucket = storage_client.bucket(cloud_config.DICOM_IMAGES_BUCKET)

    # Setup bigquery client
    bq_client = bigquery.Client(cloud_config.GCP_PROJECT_ID)

    # Create output directory
    if output_dir is not None:
        output_dir.mkdir(exist_ok=True)

    # Setup DICOM archive for outputs
    if dicom_archive is not None:
        if archive_token_url is not None and archive_client_secret is None:
            archive_client_secret = getpass(
                "Enter client secret for DICOM archive: "
            )

        web_client = get_dicom_web_client(
            url=dicom_archive,
            token_url=archive_token_url,
            client_id=archive_client_id,
            client_secret=archive_client_secret,
        )

        try:
            web_client.search_for_studies(limit=1)
        except Exception as e:
            raise RuntimeError(
                "Unsuccessful connecting to requested DICOM archive."
            ) from e

    # Loop over requested collections
    for collection in collections:
        prefix = f'cnn-nuclear-segmentations-2019/data-files/{collection}/'

        if output_dir is not None:
            collection_dir = output_dir / collection
            collection_dir.mkdir(exist_ok=True)

        # Loop over annotations in the bucket for this collection
        for ann_blob in islice(ann_bucket.list_blobs(prefix=prefix), number):
            if not ann_blob.name.endswith('.svs.tar.gz'):
                continue

            image_start_time = time()

            # Massage the blob name to derive container information
            # eg TCGA-05-4244-01Z-00-DX1.d4ff32cd-38cf-40ea-8213-45c2b100ac01
            filename = (
                ann_blob.name
                .replace(prefix, '')
                .split('/')[0]
                .replace('.svs.tar.gz', '')
            )

            # eg TCGA-05-4244-01Z-00-DX1, d4ff32cd-38cf-40ea-8213-45c2b100ac01
            if '.' in filename:
                container_id, _ = filename.split('.')
            else:
                container_id = filename

            logging.info(f"Processing container: {container_id}")

            selection_query = f"""
                SELECT
                    crdc_instance_uuid,
                    Cast(NumberOfFrames AS int) AS NumberOfFrames
                FROM
                    bigquery-public-data.idc_current.dicom_all
                WHERE
                    ContainerIdentifier='{container_id}'
                ORDER BY
                    NumberOfFrames DESC
            """
            selection_result = bq_client.query(selection_query)
            selection_df = selection_result.result().to_dataframe()

            if len(selection_df) == 0:
                # No image found, skip this for now
                logging.error(
                    f"Could not locate image for container {container_id}."
                )
                continue

            source_images = []
            for i, uuid in enumerate(selection_df.crdc_instance_uuid):
                wsi_dcm = cloud_io.read_dataset_from_blob(
                    bucket=public_bucket,
                    blob_name=f"{uuid}.dcm",
                )
                source_images.append(wsi_dcm)

                # Store to disk
                if output_dir is not None:
                    wsi_path = (
                        collection_dir / f"{container_id}_im_{i}.dcm"
                    )
                    wsi_dcm.save_as(wsi_path)

                # Store to DICOM archive
                if dicom_archive is not None:
                    web_client = get_dicom_web_client(
                        url=dicom_archive,
                        token_url=archive_token_url,
                        client_id=archive_client_id,
                        client_secret=archive_client_secret,
                    )
                    web_client.store_instances([wsi_dcm])

            ann_dcm, seg_dcms = convert_annotations(
                annotation_csvs=iter_csvs(ann_blob),
                source_images=source_images,
                include_segmentation=with_segmentation,
                segmentation_type=segmentation_type,
                annotation_coordinate_type=annotation_coordinate_type,
                dimension_organization_type=dimension_organization_type,
                create_pyramid=create_pyramid,
                store_boundary=store_boundary,
                workers=workers,
            )

            # Store objects to bucket
            if store_bucket:
                output_bucket_obj = storage_client.bucket(output_bucket)
                blob_root = (
                    "" if output_prefix is None else f"{output_prefix}/"
                )
                ann_blob_name = (
                    f"{blob_root}{collection}/{container_id}_ann.dcm"
                )

                logging.info(f"Uploading annotation to {ann_blob_name}.")
                cloud_io.write_dataset_to_blob(
                    ann_dcm,
                    output_bucket_obj,
                    ann_blob_name,
                )
                if with_segmentation:
                    for s, seg_dcm in enumerate(seg_dcms):
                        seg_blob_name = (
                            f"{blob_root}{collection}/{container_id}_seg_{s}.dcm"
                        )
                        logging.info(
                            f"Uploading segmentation to {seg_blob_name}."
                        )
                        cloud_io.write_dataset_to_blob(
                            seg_dcm,
                            output_bucket_obj,
                            seg_blob_name,
                        )

            # Store objects to filesystem
            if output_dir is not None:
                ann_path = collection_dir / f"{container_id}_ann.dcm"

                logging.info(f"Writing annotation to {str(ann_path)}.")
                ann_dcm.save_as(ann_path)

                if with_segmentation:
                    for s, seg_dcm in enumerate(seg_dcms):
                        seg_path = collection_dir / f"{container_id}_seg_{s}.dcm"
                        logging.info(f"Writing segmentation to {str(seg_path)}.")
                        seg_dcm.save_as(seg_path)

            # Store objects to DICOM archive
            if dicom_archive is not None:
                # Recreate client each time to deal with token expiration
                web_client = get_dicom_web_client(
                    url=dicom_archive,
                    token_url=archive_token_url,
                    client_id=archive_client_id,
                    client_secret=archive_client_secret,
                )

                logging.info(f"Writing annotation to {dicom_archive}.")
                web_client.store_instances([ann_dcm])

                if with_segmentation:
                    logging.info(f"Writing segmentation(s) to {dicom_archive}.")
                    for seg_dcm in seg_dcms:
                        web_client.store_instances([seg_dcm])

            image_stop_time = time()
            time_for_image = image_stop_time - image_start_time
            logging.info(f"Processed {container_id} in {time_for_image:.2f}s")


if __name__ == "__main__":
    run()
