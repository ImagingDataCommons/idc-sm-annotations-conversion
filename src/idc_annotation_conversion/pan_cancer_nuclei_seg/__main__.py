"""Entrypoint for Pan Cancer Nuclei Segmentation annotation conversion."""
import datetime
from io import BytesIO, BufferedReader
from getpass import getpass
import logging
from pathlib import Path
import os
import tarfile
from time import time
from typing import Any, Generator, Optional

import click
from dicomweb_client import DICOMwebClient
from google.cloud import bigquery
from google.cloud import storage
import highdicom as hd
from oauthlib.oauth2 import BackendApplicationClient
import pandas as pd
from requests_oauthlib import OAuth2Session

from idc_annotation_conversion import cloud_io
from idc_annotation_conversion.pan_cancer_nuclei_seg.convert import (
    get_graphic_data,
    create_bulk_annotations,
    create_segmentations,
)
from idc_annotation_conversion import cloud_config
from idc_annotation_conversion.mp_utils import Pipeline


ANNOTATION_BUCKET = "tcia-nuclei-seg"
ANNOTATION_PREFIX = 'cnn-nuclear-segmentations-2019/data-files/'


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


def get_ann_blob_name(
    output_blob_root: str,
    collection: str,
    container_id: str
) -> str:
    """Get the name of the blob where an annotation should be stored."""
    return f"{output_blob_root}{collection}/{container_id}_ann.dcm"


def get_seg_blob_name(
    output_blob_root: str,
    collection: str,
    container_id: str,
    pyramid_level: int,
) -> str:
    """Get the name of the blob where a segmentation should be stored."""
    return (
        f"{output_blob_root}{collection}/"
        f"{container_id}_seg_{pyramid_level}.dcm"
    )


def split_csv_blob_name(csv_blob: storage.Blob) -> tuple[str, str]:
    """Get collection name and container name from blob name."""
    collection = csv_blob.name.split("/")[-3]
    collection_prefix = f'{ANNOTATION_PREFIX}{collection}/'

    # Massage the blob name to derive container information
    # eg TCGA-05-4244-01Z-00-DX1.d4ff32cd-38cf-40ea-8213-45c2b100ac01
    filename = (
        csv_blob.name
        .replace(collection_prefix, '')
        .split('/')[0]
        .replace('.svs.tar.gz', '')
    )

    # eg TCGA-05-4244-01Z-00-DX1, d4ff32cd-38cf-40ea-8213-45c2b100ac01
    if '.' in filename:
        container_id, _ = filename.split('.')
    else:
        container_id = filename

    return collection, container_id


def iter_csvs(ann_bytes: bytes) -> Generator[BufferedReader, None, None]:
    """Iterate over individual CSV files in a loaded tar.gz file.

    Parameters:
    -----------
    ann_bytes: bytes
        Raw bytes of the downloaded .tar.gz file.

    Yields:
    -------
    io.BufferedReader
        Buffer for each CSV file in the original tar. This buffer can be
        used to load the contents of the CSV from the in-memory tar file.

    """
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


class FileDownloader:

    """Object that finds and downloads files relevant to a case."""

    def __init__(
        self,
        output_dir: Optional[Path],
        dicom_archive: Optional[str] = None,
        archive_token_url: Optional[str] = None,
        archive_client_id: Optional[str] = None,
        archive_client_secret: Optional[str] = None,
    ):
        """

        Parameters
        ----------
        output_dir: Optional[pathlib.Path]
            A local output directory to store a copy of the downloaded files
            in, if required.
        dicom_archive: Optional[str], optional
            Additionally store images to this DICOM archive.
        archive_token_url: Optional[str], optional
            URL to use to request an OAuth token to access the archive.,
        archive_client_id: Optional[str], optional
            Client ID to use for OAuth token request.,
        archive_client_secret: Optional[str], optional
            Client secret to use for OAuth token request. If none, user will
            be prompted for secret.

        """
        self._output_dir = output_dir
        self._dicom_archive = dicom_archive
        self._archive_token_url = archive_token_url
        self._archive_client_id = archive_client_id
        self._archive_client_secret = archive_client_secret

        # Setup bigquery client
        self._bq_client = bigquery.Client(cloud_config.GCP_PROJECT_ID)

        # Public bucket containing IDC images
        self._storage_client = storage.Client(
            project=cloud_config.GCP_PROJECT_ID
        )
        self._public_bucket = self._storage_client.bucket(
            cloud_config.DICOM_IMAGES_BUCKET
        )

        self._errors = []

    def __call__(
        self,
        csv_blob: storage.Blob,
    ) -> dict[str, Any]:
        """Pull all necessary files.

        Parameters
        ----------
        csv_blob: google.cloud.storage.Blob
            CSV blob to process.

        Returns
        -------
        data: dict[str, Any]
            Output data packed into a dict. This will contain the same keys as
            the input dictionary, plus the following additional keys:

            - collection: str
                Collection name for the case.
            - container_id: str
                Container ID for the case.
            - csv_bytes: bytes
                Raw bytes of the downloaded annotation .tar.gz file.
            - source_images: list[pydicom.Dataset]
                Source images for this case.

        """
        start_time = time()
        collection, container_id = split_csv_blob_name(csv_blob)

        try:
            # Download the annotation
            csv_bytes = csv_blob.download_as_bytes()

            logging.info(f"Pulling images for container: {container_id}")

            if self._output_dir is not None:
                collection_dir = self._output_dir / collection
                collection_dir.mkdir(exist_ok=True)

            selection_query = f"""
                SELECT
                    gcs_url,
                    Cast(NumberOfFrames AS int) AS NumberOfFrames
                FROM
                    bigquery-public-data.idc_current.dicom_all
                WHERE
                    ContainerIdentifier='{container_id}'
                ORDER BY
                    NumberOfFrames DESC
            """
            selection_result = self._bq_client.query(selection_query)
            selection_df = selection_result.result().to_dataframe()

            if len(selection_df) == 0:
                # No image found, skip this for now
                msg = f"Could not locate image for container {container_id}."
                raise RuntimeError(msg)

            source_images = []
            for i, url in enumerate(selection_df.gcs_url):
                blob_name = "/".join(url.split("/")[3:])
                wsi_dcm = cloud_io.read_dataset_from_blob(
                    bucket=self._public_bucket,
                    blob_name=blob_name,
                )
                source_images.append(wsi_dcm)

                # Store to disk
                if self._output_dir is not None:
                    wsi_path = (
                        collection_dir / f"{container_id}_im_{i}.dcm"
                    )
                    wsi_dcm.save_as(wsi_path)

                # Store to DICOM archive
                if self._dicom_archive is not None:
                    web_client = get_dicom_web_client(
                        url=self._dicom_archive,
                        token_url=self._archive_token_url,
                        client_id=self._archive_client_id,
                        client_secret=self._archive_client_secret,
                    )
                    web_client.store_instances([wsi_dcm])

            # The pixels are no longer needed, and they take loads of memory.
            # So just delete them
            for wsi_dcm in source_images:
                del wsi_dcm.PixelData

        except Exception as e:
            logging.error(f"Error {str(e)}")
            self._errors.append(
                {
                    "collection": collection,
                    "container_id": container_id,
                    "error_message": str(e),
                    "datetime": str(datetime.datetime.now()),
                }
            )
            errors_df = pd.DataFrame(self._errors)
            errors_df.to_csv("download_error_log.csv")
            return None

        stop_time = time()
        duration = stop_time - start_time
        logging.info(f"Pulled images for {container_id} in {duration:.2f}s")

        return dict(
            collection=collection,
            container_id=container_id,
            csv_bytes=csv_bytes,
            source_images=source_images,
        )


class CSVParser:

    """Class that parses CSV annotations to graphic data."""

    def __init__(
        self,
        annotation_coordinate_type: str,
        graphic_type: str,
        workers: int,
    ):
        """

        Parameters
        ----------
        graphic_type: str
            Graphic type to use to store all nuclei. Note that all but
            'POLYGON' result in simplification and loss of information in the
            annotations.
        annotation_coordinate_type: str
            Store coordinates in the Bulk Microscopy Bulk Simple Annotations in
            the (3D) frame of reference (SCOORD3D), or the (2D) total pixel
            matrix (SCOORD, default).
        workers: int
            Number of subprocess workers to spawn. If 0, all computation will
            use the main thread.

        """
        self._annotation_coordinate_type = annotation_coordinate_type
        self._graphic_type = graphic_type
        self._workers = workers
        self._errors = []

    def __call__(self, data: dict[str, Any]) -> dict[str, Any]:
        """Parse CSV to graphic data.

        Parameters
        ----------
        data: dict[str, Any]
            Input data packed into a dict, including at least:

            - collection: str
                Collection name for the case.
            - container_id: str
                Container ID for the case.
            - csv_bytes: bytes
                Raw bytes of the downloaded annotation .tar.gz file.
            - source_images: list[pydicom.Dataset]
                Source images for this case.

        Returns
        -------
        data: dict[str, Any]
            Output data packed into a dict. This will contain the same keys as
            the input dictionary, plus the following additional keys:

            - polygons: list[shapely.geometry.polygon.Polygon]
                List of polygons. Note that this is always the list of full
                original, 2D polygons regardless of the requested graphic type
                and annotation coordinate type.
            - graphic_data: list[np.ndarray]
                List of graphic data as numpy arrays in the format required for
                the MicroscopyBulkSimpleAnnotations object. These are correctly
                formatted for the requested graphic type and annotation
                coordinate type.
            - areas: list[float]
                Areas for each of the polygons, in square micrometers. Does not
                depend on requested graphic type or coordinate type.

        """
        start_time = time()
        container_id = data["container_id"]
        logging.info(f"Parsing CSV for container: {container_id}")
        try:
            polygons, graphic_data, areas = get_graphic_data(
                annotation_csvs=iter_csvs(data["csv_bytes"]),
                source_image_metadata=data["source_images"][0],
                graphic_type=self._graphic_type,
                annotation_coordinate_type=self._annotation_coordinate_type,
                workers=self._workers,
            )
        except Exception as e:
            logging.error(f"Error {str(e)}")
            self._errors.append(
                {
                    "collection": data["collection"],
                    "container_id": data["container_id"],
                    "error_message": str(e),
                    "datetime": str(datetime.datetime.now()),
                }
            )
            errors_df = pd.DataFrame(self._errors)
            errors_df.to_csv("conversion_error_log.csv")
            return None

        stop_time = time()
        duration = stop_time - start_time
        logging.info(
            f"Processed CSV for container {container_id} in {duration:.2f}s"
        )
        del data["csv_bytes"]  # save some memory
        data["polygons"] = polygons
        data["graphic_data"] = graphic_data
        data["areas"] = areas
        return data


class AnnotationCreator:

    """Class that creates bulk annotations DICOM objects."""

    def __init__(
        self,
        graphic_type: str,
        annotation_coordinate_type: str,
    ):
        """

        Parameters
        ----------
        graphic_type: str
            Graphic type to use to store all nuclei. Note that all but
            'POLYGON' result in simplification and loss of information in the
            annotations.
        annotation_coordinate_type: str
            Store coordinates in the Bulk Microscopy Bulk Simple Annotations in
            the (3D) frame of reference (SCOORD3D), or the (2D) total pixel
            matrix (SCOORD, default).

        """
        self._annotation_coordinate_type = annotation_coordinate_type
        self._graphic_type = graphic_type
        self._errors = []

    def __call__(self, data: dict[str, Any]) -> dict[str, Any]:
        """Parse CSV to graphic data.

        Parameters
        ----------
        data: dict[str, Any]
            Input data packed into a dict, including at least:

            - collection: str
                Collection name for the case.
            - container_id: str
                Container ID for the case.
            - graphic_data: list[np.ndarray]
                List of graphic data as numpy arrays in the format required for
                the MicroscopyBulkSimpleAnnotations object. These are correctly
                formatted for the requested graphic type and annotation
                coordinate type.
            - source_images: list[pydicom.Dataset]
                Source images for this case.
            - areas: list[float]
                Areas for each of the polygons, in square micrometers. Does not
                depend on requested graphic type or coordinate type.

        Returns
        -------
        data: dict[str, Any]
            Output data packed into a dict. This will contain the same keys as
            the input dictionary, plus the following additional keys:

            - ann_dcm: hd.ann.MicroscopyBulkSimpleAnnotations:
                DICOM bulk microscopy annotation encoding the original
                annotations in vector format.

        """
        # Unpack inputs
        collection = data["collection"]
        container_id = data["container_id"]
        source_images = data["source_images"]
        graphic_data = data["graphic_data"]
        areas = data["areas"]

        start_time = time()

        logging.info(f"Creating annotation for container: {container_id}")

        try:
            ann_dcm = create_bulk_annotations(
                source_image_metadata=source_images[0],
                graphic_data=graphic_data,
                areas=areas,
                annotation_coordinate_type=self._annotation_coordinate_type,
                graphic_type=self._graphic_type,
            )
        except Exception as e:
            logging.error(f"Error {str(e)}")
            self._errors.append(
                {
                    "collection": collection,
                    "container_id": container_id,
                    "error_message": str(e),
                    "datetime": str(datetime.datetime.now()),
                }
            )
            errors_df = pd.DataFrame(self._errors)
            errors_df.to_csv("annotation_creator_error_log.csv")
            return None

        stop_time = time()
        duration = stop_time - start_time
        logging.info(
            f"Created annotation for {container_id} in {duration:.2f}s"
        )

        data["ann_dcm"] = ann_dcm

        # Save some memory
        del data["graphic_data"]
        del data["areas"]

        return data


class SegmentationCreator:

    """Class that creates DICOM segmentations from graphic data."""

    def __init__(
        self,
        segmentation_type: str,
        dimension_organization_type: str,
        create_pyramid: bool,
    ):
        """

        Parameters
        ----------
        segmentation_type: str
            Segmentation type (BINARY or FRACTIONAL) for the Segmentation Image
            (if any).
        dimension_organization_type: str
            Dimension organization type of the output segmentations.
        create_pyramid: bool
            Whether to create a full pyramid of segmentations (rather than a
            single segmentation at the highest resolution).

        """
        self._segmentation_type = segmentation_type
        self._dimension_organization_type = dimension_organization_type
        self._create_pyramid = create_pyramid
        self._errors = []

    def __call__(self, data: dict[str, Any]) -> dict[str, Any]:
        """Parse CSV to graphic data.

        Parameters
        ----------
        data: dict[str, Any]
            Input data packed into a dict, including at least:

            - collection: str
                Collection name for the case.
            - container_id: str
                Container ID for the case.
            - graphic_data: list[np.ndarray]
                List of graphic data as numpy arrays in the format required for
                the MicroscopyBulkSimpleAnnotations object. These are correctly
                formatted for the requested graphic type and annotation
                coordinate type.
            - source_images: list[pydicom.Dataset]
                Source images for this case.
            - areas: list[float]
                Areas for each of the polygons, in square micrometers. Does not
                depend on requested graphic type or coordinate type.

        Returns
        -------
        data: dict[str, Any]
            Output data packed into a dict. This will contain the same keys as
            the input dictionary, plus the following additional keys:

            - segmentations: List[hd.seg.Segmentation]
                DICOM segmentation image(s) encoding the original annotations in raster
                format.

        """
        # Unpack inputs
        collection = data["collection"]
        container_id = data["container_id"]
        source_images = data["source_images"]
        polygons = data["polygons"]

        image_start_time = time()

        logging.info(f"Creating segmentation for: {container_id}")

        try:
            seg_dcms = create_segmentations(
                polygons=polygons,
                source_images=source_images,
                dimension_organization_type=self._dimension_organization_type,
                segmentation_type=self._segmentation_type,
                create_pyramid=self._create_pyramid,
            )
        except Exception as e:
            logging.error(f"Error {str(e)}")
            self._errors.append(
                {
                    "collection": collection,
                    "container_id": container_id,
                    "error_message": str(e),
                    "datetime": str(datetime.datetime.now()),
                }
            )
            errors_df = pd.DataFrame(self._errors)
            errors_df.to_csv("segmentation_creator_error_log.csv")
            return None

        image_stop_time = time()
        duration = image_stop_time - image_start_time
        logging.info(
            f"Created segmentations for container {container_id} in "
            f"{duration:.2f}s"
        )

        del data["polygons"]
        del data["source_images"]
        for seg in seg_dcms:
            del seg._db_man  # Prevents pickling...
        data["seg_dcms"] = seg_dcms

        return data


class FileUploader:

    def __init__(
        self,
        output_bucket_obj: Optional[storage.Bucket],
        output_dir: Optional[Path],
        dicom_archive: Optional[str] = None,
        archive_token_url: Optional[str] = None,
        archive_client_id: Optional[str] = None,
        archive_client_secret: Optional[str] = None,
    ):
        """

        Parameters
        ----------
        output_bucket_obj: Optional[google.storage.Bucket]
            Output bucket, if storing in a bucket.
        output_dir: Optional[pathlib.Path]
            A local output directory to store a copy of the downloaded files
            in, if required.
        dicom_archive: Optional[str], optional
            Additionally store images to this DICOM archive.
        archive_token_url: Optional[str], optional
            URL to use to request an OAuth token to access the archive.,
        archive_client_id: Optional[str], optional
            Client ID to use for OAuth token request.,
        archive_client_secret: Optional[str], optional
            Client secret to use for OAuth token request. If none, user will
            be prompted for secret.

        """
        self._output_bucket_obj = output_bucket_obj
        self._output_dir = output_dir
        self._dicom_archive = dicom_archive
        self._archive_token_url = archive_token_url
        self._archive_client_id = archive_client_id
        self._archive_client_secret = archive_client_secret
        self._errors = []

    def __call__(
        self,
        data: dict[str, Any],
    ) -> None:
        """Upload files.

        Parameters
        ----------
        data: dict[str, Any]
            Input data packed into a dictionary, containing at least:

            - collection: str
                Collection name for the case.
            - container_id: str
                Container ID for the case.
            - ann_dcm: highdicom.ann.MicroscopyBulkSimpleAnnotations
                Bulk Annotation DICOM object.
            - seg_dcm: Optional[list[highdicom.seg.Segmentation]]
                Converted segmentations, if required.

        """
        image_start_time = time()

        # Unpack inputs
        collection = data["collection"]
        container_id = data["container_id"]
        ann_dcm = data["ann_dcm"]
        seg_dcms = data.get("seg_dcms")

        logging.info(f"Uploading annotations for {container_id}")

        if self._output_dir is not None:
            collection_dir = self._output_dir / collection
            collection_dir.mkdir(exist_ok=True)

        try:

            # Store objects to bucket
            if self._output_bucket_obj is not None:
                ann_blob_name = get_ann_blob_name(
                    output_blob_root=self._output_blob_root,
                    collection=collection,
                    container_id=container_id,
                )

                logging.info(f"Uploading annotation to {ann_blob_name}.")
                cloud_io.write_dataset_to_blob(
                    ann_dcm,
                    self._output_bucket_obj,
                    ann_blob_name,
                )
                if seg_dcms is not None:
                    for s, seg_dcm in enumerate(seg_dcms):
                        seg_blob_name = get_seg_blob_name(
                            output_blob_root=self._output_blob_root,
                            collection=collection,
                            container_id=container_id,
                            pyramid_level=s,
                        )
                        logging.info(
                            f"Uploading segmentation to {seg_blob_name}."
                        )
                        cloud_io.write_dataset_to_blob(
                            seg_dcm,
                            self._output_bucket_obj,
                            seg_blob_name,
                        )

            # Store objects to filesystem
            if self._output_dir is not None:
                ann_path = collection_dir / f"{container_id}_ann.dcm"

                logging.info(f"Writing annotation to {str(ann_path)}.")
                ann_dcm.save_as(ann_path)

                if seg_dcms is not None:
                    for s, seg_dcm in enumerate(seg_dcms):
                        seg_path = (
                            collection_dir / f"{container_id}_seg_{s}.dcm"
                        )
                        logging.info(
                            f"Writing segmentation to {str(seg_path)}."
                        )
                        seg_dcm.save_as(seg_path)

            # Store objects to DICOM archive
            if self._dicom_archive is not None:
                # Recreate client each time to deal with token expiration
                web_client = get_dicom_web_client(
                    url=self._dicom_archive,
                    token_url=self._archive_token_url,
                    client_id=self._archive_client_id,
                    client_secret=self._archive_client_secret,
                )

                logging.info(f"Writing annotation to {self._dicom_archive}.")
                web_client.store_instances([ann_dcm])

                if seg_dcms is not None:
                    logging.info(
                        f"Writing segmentation(s) to {self._dicom_archive}."
                    )
                    for seg_dcm in seg_dcms:
                        web_client.store_instances([seg_dcm])

            image_stop_time = time()
            time_for_image = image_stop_time - image_start_time
            logging.info(
                f"Uploaded annotations for {container_id} in "
                f"{time_for_image:.2f}s"
            )
        except Exception as e:
            logging.error(f"Error {str(e)}")
            self._errors.append(
                {
                    "collection": collection,
                    "container_id": container_id,
                    "error_message": str(e),
                    "datetime": str(datetime.datetime.now()),
                }
            )
            errors_df = pd.DataFrame(self._errors)
            errors_df.to_csv("upload_error_log.csv")
            return None


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
    "-l",
    "--csv-blob",
    help=(
        "Specify a single CSV blob to process, using its path within "
        "the bucket."
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
    "--output-prefix",
    "-p",
    help="Prefix for all output blobs. Default no prefix.",
)
@click.option(
    "--graphic-type",
    "-g",
    default="POLYGON",
    type=click.Choice(
        [v.name for v in hd.ann.GraphicTypeValues],
        case_sensitive=False,
    ),
    show_default=True,
    help=(
        "Graphic type to use to store all nuclei. Note that all "
        "but 'POLYGON' result in simplification and loss of "
        "information in the annotations."
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
@click.option(
    "--pull-process/--no-pull-process",
    "-r/-R",
    default=True,
    show_default=True,
    help="Use a separate process to pull images.",
)
def run(
    collections: Optional[list[str]],
    number: Optional[int],
    output_dir: Optional[Path],
    output_bucket: str,
    output_prefix: Optional[str],
    store_bucket: bool,
    graphic_type: str,
    store_wsi_dicom: bool,
    annotation_coordinate_type: str,
    with_segmentation: bool,
    segmentation_type: str,
    dimension_organization_type: str,
    create_pyramid: bool,
    csv_blob: Optional[str] = None,
    dicom_archive: Optional[str] = None,
    archive_token_url: Optional[str] = None,
    archive_client_id: Optional[str] = None,
    archive_client_secret: Optional[str] = None,
    keep_existing: bool = False,
    workers: int = 0,
    pull_process: bool = True,
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
    ann_bucket = storage_client.bucket(ANNOTATION_BUCKET)
    output_blob_root = (
        "" if output_prefix is None else f"{output_prefix}/"
    )
    if store_bucket:
        if output_bucket is None:
            date_str = (datetime.date.today())
            output_bucket = (
                f"pan_cancer_nuclei_seg_annotation_conversion_{date_str}"
            )
        output_bucket_obj = output_client.bucket(output_bucket)

        if not output_bucket_obj.exists():
            output_bucket_obj.create(
                location=cloud_config.GCP_DEFAULT_LOCATION
            )

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

    logging.info("Listing CSVs")

    if csv_blob is not None:
        csv_blob_obj = ann_bucket.get_blob(csv_blob)
        if not csv_blob_obj.exists():
            raise RuntimeError(f"No such blob found: {csv_blob}")
        to_process = [csv_blob_obj]
    else:
        to_process = []

        # Loop over requested collections
        for collection in collections:
            collection_prefix = f'{ANNOTATION_PREFIX}{collection}/'

            collection_blobs = [
                b for b in ann_bucket.list_blobs(prefix=collection_prefix)
                if b.name.endswith('.svs.tar.gz')
            ]

            for blob in collection_blobs:

                collection, container_id = split_csv_blob_name(blob)

                # Check whether the output blobs already exist, and skip if
                # they do
                if keep_existing:
                    ann_blob_name = get_ann_blob_name(
                        output_blob_root=output_blob_root,
                        collection=collection,
                        container_id=container_id,
                    )
                    ann_output_blob = output_bucket_obj.get_blob(ann_blob_name)
                    if with_segmentation:
                        seg_blob_name = get_seg_blob_name(
                            output_blob_root=output_blob_root,
                            collection=collection,
                            container_id=container_id,
                            pyramid_level=0,
                        )
                        seg_output_blob = output_bucket_obj.get_blob(
                            seg_blob_name
                        )
                        if (
                            ann_output_blob is not None and
                            seg_output_blob is not None
                        ):
                            continue
                    else:
                        if ann_output_blob is not None:
                            continue

                to_process.append(blob)

        if number is not None:
            to_process = to_process[:number]

    logging.info(f"Found {len(to_process)} CSVs to process")

    operations = []
    pull_kwargs = dict(
        output_dir=output_dir if store_wsi_dicom else None,
        dicom_archive=dicom_archive,
        archive_token_url=archive_token_url,
        archive_client_id=archive_client_id,
        archive_client_secret=archive_client_secret,
    )
    operations.append((FileDownloader, [], pull_kwargs))
    parser_kwargs = dict(
        annotation_coordinate_type=annotation_coordinate_type,
        graphic_type=graphic_type,
        workers=workers,
    )
    operations.append((CSVParser, [], parser_kwargs))
    annotation_creator_kwargs = dict(
        annotation_coordinate_type=annotation_coordinate_type,
        graphic_type=graphic_type,
    )
    operations.append((AnnotationCreator, [], annotation_creator_kwargs))
    if with_segmentation:
        segmentation_creator_kwargs = dict(
            segmentation_type=segmentation_type,
            dimension_organization_type=dimension_organization_type,
            create_pyramid=create_pyramid,
        )
        operations.append(
            (SegmentationCreator, [], segmentation_creator_kwargs)
        )
    upload_kwargs = dict(
        output_dir=output_dir if store_wsi_dicom else None,
        output_bucket_obj=output_bucket_obj if store_bucket else None,
        dicom_archive=dicom_archive,
        archive_token_url=archive_token_url,
        archive_client_id=archive_client_id,
        archive_client_secret=archive_client_secret,
    )
    operations.append((FileUploader, [], upload_kwargs))
    pipeline = Pipeline(
        operations,
        same_process=not pull_process,
    )
    pipeline(to_process)


if __name__ == "__main__":
    run()
