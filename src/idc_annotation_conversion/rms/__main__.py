"""Main conversion process for the RMS project."""
import os
from itertools import islice
from typing import Optional
from xml.etree import ElementTree

import click
from google.cloud import bigquery, storage
import highdicom as hd
import numpy as np

from idc_annotation_conversion import cloud_config
from idc_annotation_conversion.rms.convert import convert_xml_annotation


# Bucket containing annotations for this project
ANNOTATION_BUCKET_PROJECT = "idc-external-031"
ANNOTATION_BUCKET = "rms_annotation_test_oct_2023"


@click.command()
@click.option(
    "--number",
    "-n",
    type=int,
    help="Number to process per collection. All by default.",
)
def run(
    number: Optional[int],
):
    # Setup project and authenticate
    os.environ["GCP_PROJECT_ID"] = cloud_config.GCP_PROJECT_ID

    # Access bucket containing annotations
    storage_client = storage.Client(project=cloud_config.GCP_PROJECT_ID)
    public_bucket = storage_client.bucket(cloud_config.DICOM_IMAGES_BUCKET)

    ann_storage_client = storage.Client(project=ANNOTATION_BUCKET_PROJECT)
    ann_bucket = ann_storage_client.bucket(ANNOTATION_BUCKET)

    # Setup bigquery client
    bq_client = bigquery.Client(cloud_config.GCP_PROJECT_ID)

    prefix = "RMS-XML-hand-annotations"
    # prefix = "rms_annotation_test_oct_2023/RMS-XML-hand-"
    for ann_blob in islice(ann_bucket.list_blobs(prefix=prefix), number):
        if not ann_blob.name.endswith(".xml"):
            continue

        print(ann_blob.name)
        id1, id2 = ann_blob.name.split("/")[1].replace(".xml", "").split("-", maxsplit=1)
        print(id1, id2)
        selection_query = f"""
            SELECT
                crdc_instance_uuid,
                Cast(NumberOfFrames AS int) AS NumberOfFrames
            FROM
                bigquery-public-data.idc_current.dicom_all
            WHERE
                ContainerIdentifier='{id2}'
            ORDER BY
                NumberOfFrames DESC
        """
        selection_result = bq_client.query(selection_query)
        selection_df = selection_result.result().to_dataframe()
        print(selection_df)
        text = ann_blob.download_as_text()
        xml_root = ElementTree.fromstring(text)

        sr = convert_xml_annotation(xml_root, source_image)





if __name__ == "__main__":
    run()
