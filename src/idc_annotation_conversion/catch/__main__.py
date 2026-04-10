import datetime
from itertools import islice
import logging
from pathlib import Path
import sqlite3

import click
from google.cloud import storage
import highdicom as hd
from highdicom.sr import MeasurementsAndQualitativeEvaluations
import numpy as np

from idc_annotation_conversion import cloud_io, cloud_config
from idc_annotation_conversion.catch import metadata_config


IMAGES_BUCKET_PROJECT = "idc-converted-data"
IMAGES_BUCKET = "catch-pathology"

ANNOTATIONS_BUCKET_PROJECT = "idc-source-data"
ANNOTATIONS_BUCKET = "catch_pathology_annotations"

SLIDERUNNER_TYPE_GRAPHIC_TYPE_MAP = {
    1: hd.ann.GraphicTypeValues.POINT,
    2: hd.ann.GraphicTypeValues.RECTANGLE,  # guess, would need to check
    3: hd.ann.GraphicTypeValues.POLYGON,
}


@click.command()
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
    "--store-wsi-dicom/--no-store-wsi-dicom",
    "-d/-D",
    default=False,
    show_default=True,
    help=(
        "Download all WSI DICOM files and store in the output directory "
        "(if any)."
    ),
)
def main(
    number: int | None,
    output_dir: Path | None,
    output_bucket: str | None,
    store_bucket: bool,
    store_wsi_dicom: bool,
):
    """Convert manual annotations from the Canine Cutaneous Cancer Histology dataset
    to DICOM Microscopy Bulk Simple Annotations

    TCIA Link: https://doi.org/10.7937/TCIA.2M93-FX66
    Paper: https://www.nature.com/articles/s41597-022-01692-w.pdf

    Annotations are originally in a SQLITE database created by the SlideRunner software.

    """
    image_storage_client = storage.Client(project=IMAGES_BUCKET_PROJECT)
    image_bucket = image_storage_client.bucket(IMAGES_BUCKET)

    annotations_storage_client = storage.Client(project=ANNOTATIONS_BUCKET_PROJECT)
    annotations_bucket = annotations_storage_client.bucket(ANNOTATIONS_BUCKET)

    output_client = storage.Client(project=cloud_config.OUTPUT_GCP_PROJECT_ID)
    if store_bucket:
        if output_bucket is None:
            today = datetime.date.today()
            output_bucket = f"catch-annotations-{today}"

        output_bucket_obj = output_client.bucket(output_bucket)

        if not output_bucket_obj.exists():
            output_bucket_obj.create(
                location=cloud_config.GCP_DEFAULT_LOCATION
            )
    else:
        output_bucket_obj = None

    # Create output directory
    if output_dir is not None:
        output_dir.mkdir(exist_ok=True)

    annotations_db_blob = annotations_bucket.blob("CATCH.sqlite")
    db_file = Path.cwd() / "CATCH.sqlite"

    if not db_file.exists():
        annotations_db_blob.download_to_file(db_file)

    db = sqlite3.connect(db_file)

    slide_query = "SELECT uid, filename FROM Slides;"

    for slide_id, slide_filename in islice(db.execute(slide_query), number):

        slide_name = slide_filename.replace('.svs', '')
        image_blob_name = f"{slide_name}/DCM_0"
        print(image_blob_name)
        image_dcm = cloud_io.read_dataset_from_blob(
            image_bucket,
            image_blob_name,
            stop_before_pixels=not store_wsi_dicom,
        )

        # Each distinct combination of annotation type and class will become
        # its own annotation group in the ANN object
        groups_query = (
            "SELECT DISTINCT type, agreedClass FROM Annotations "
            f"WHERE slide = {slide_id} AND deleted = 0;"
        )

        annotation_groups = []

        for grp_id, (ann_type, ann_class) in enumerate(db.execute(groups_query), 1):

            graphic_type = SLIDERUNNER_TYPE_GRAPHIC_TYPE_MAP[ann_type]

            class_query = (
                f"SELECT name, color FROM Classes WHERE uid = {ann_class}"
            )
            class_name, color = list(db.execute(class_query))[0]

            annotation_query = (
                "SELECT uid FROM Annotations "
                f"WHERE slide = {slide_id} "
                "AND deleted = 0 "
                f"AND type = {ann_type} "
                f"AND agreedClass = {ann_class};"
            )

            graphic_data = []

            for (ann_id, ) in db.execute(annotation_query):

                coordinate_query = (
                    "SELECT coordinateX, coordinateY FROM Annotations_coordinates "
                    f"WHERE annoId = {ann_id} ORDER BY orderIdx;"
                )

                arr = np.array(list(db.execute(coordinate_query)))

                # Remove the last coordinate because the last duplicate
                # coodinate should not be included. Sometimes it seems that
                # there are multiple duplicates at the end
                while (arr[0] == arr[-1]).all():
                    arr = arr[:-1]

                graphic_data.append(arr)

            grp = hd.ann.AnnotationGroup(
                number=grp_id,
                uid=hd.UID(),
                label=class_name,
                annotated_property_category=metadata_config.category_codes[class_name],
                annotated_property_type=metadata_config.finding_codes[class_name],
                graphic_data=graphic_data,
                graphic_type=graphic_type,
                algorithm_type="MANUAL",
                display_color=hd.color.CIELabColor.from_string(color),
            )

            annotation_groups.append(grp)

        ann_dcm = hd.ann.MicroscopyBulkSimpleAnnotations(
            [image_dcm],
            annotation_coordinate_type=hd.ann.AnnotationCoordinateTypeValues.SCOORD,
            annotation_groups=annotation_groups,
            series_description=metadata_config.series_description,
            series_instance_uid=hd.UID(),
            sop_instance_uid=hd.UID(),
            instance_number=1,
            series_number=metadata_config.series_number,
            manufacturer=metadata_config.manufacturer,
            manufacturer_model_name=metadata_config.manufacturer_model_name,
            content_label=metadata_config.content_label,
            device_serial_number=metadata_config.device_serial_number,
            content_description=metadata_config.content_description,
            software_versions=metadata_config.software_versions,
            contributing_equipment=metadata_config.contributing_equipment,
            content_creator_name=metadata_config.content_creator_name,
        )

        # Override collection DOI
        ann_dcm.OtherClinicalTrialProtocolIDsSequence = [
            metadata_config.clinical_trial_ids_item
        ]

        ann_name = slide_filename.replace(".svs", '_ann.dcm')
        im_name = slide_filename.replace(".svs", '_im.dcm')

        # Store objects to filesystem
        if output_dir is not None:
            out_path = output_dir / ann_name
            logging.info(f"Writing annotation to {str(out_path)}.")
            ann_dcm.save_as(out_path)

            if store_wsi_dicom:
                slide_path = output_dir / im_name
                image_dcm.save_as(slide_path)

        # Store to bucket
        if output_bucket_obj is not None:
            logging.info("Writing objects to output bucket.")
            cloud_io.write_dataset_to_blob(
                ann_dcm,
                output_bucket_obj,
                ann_name,
            )
            cloud_io.write_dataset_to_blob(
                image_dcm,
                output_bucket_obj,
                im_name,
            )


if __name__ == "__main__":
    main()
