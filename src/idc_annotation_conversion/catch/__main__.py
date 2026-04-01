import sqlite3

import click
from google.cloud import storage
import highdicom as hd
import numpy as np

from idc_annotation_conversion import cloud_io
from idc_annotation_conversion.catch import metadata_config


# TODO read input sqlite file from bucket
# TODO output to bucket
# TODO metadata


IMAGES_BUCKET_PROJECT = "idc-converted-data"
IMAGES_BUCKET = "catch-pathology"


SLIDERUNNER_TYPE_GRAPHIC_TYPE_MAP = {
    1: hd.ann.GraphicTypeValues.POINT,
    2: hd.ann.GraphicTypeValues.RECTANGLE,  # guess, would need to check
    3: hd.ann.GraphicTypeValues.POLYGON,
}


def main():
    """Convert manual annotations from the Canine Cutaneous Cancer Histology dataset
    to DICOM Microscopy Bulk Simple Annotations

    https://doi.org/10.7937/TCIA.2M93-FX66

    Annotations are originally in a SQLITE database created by the SlideRunner software.

    """
    # TODO download this first
    db = sqlite3.connect("file:///Users/cpb28/Developer/idc-annotation-conversion/CATCH.sqlite")

    slide_query = (
        "SELECT uid, filename FROM Slides;"
    )

    image_storage_client = storage.Client(project=IMAGES_BUCKET_PROJECT)
    image_bucket = image_storage_client.bucket(IMAGES_BUCKET)

    for slide_id, slide_filename in db.execute(slide_query):

        slide_name = slide_filename.replace('.svs', '')
        image_blob_name = f"{slide_name}/DCM_0"
        print(image_blob_name)
        image_dcm = cloud_io.read_dataset_from_blob(
            image_bucket,
            image_blob_name,
            stop_before_pixels=True,
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
        )
        ann_dcm.save_as(slide_filename.replace(".svs", '_ann.dcm'))


if __name__ == "__main__":
    main()
