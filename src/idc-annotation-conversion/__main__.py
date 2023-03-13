from io import BytesIO
import os
import tarfile

from google.cloud import bigquery
from google.cloud import storage
import pandas as pd
import pydicom


def run():
    # Setup project and authenticate
    project_id = "idc-etl-processing"
    os.environ["GCP_PROJECT_ID"] = project_id

    # Access bucket containing annotations
    storage_client = storage.Client(project=project_id)
    ann_bucket = storage_client.bucket('tcia-nuclei-seg')
    public_bucket = storage_client.bucket('public-datasets-idc')

    # Setup bigquery client
    bq_client = bigquery.Client(project_id)

    subset = 'luad_polygon'
    prefix = f'cnn-nuclear-segmentations-2019/data-files/{subset}/'
    for ann_blob in ann_bucket.list_blobs(prefix=prefix):

        if not ann_blob.name.endswith('.svs.tar.gz'):
            continue

        # e.g. TCGA-05-4244-01Z-00-DX1.d4ff32cd-38cf-40ea-8213-45c2b100ac01
        filename = (
            ann_blob.name
            .replace(prefix, '')
            .split('/')[0]
            .replace('.svs.tar.gz', '')
        )

        # e.g. TCGA-05-4244-01Z-00-DX1, d4ff32cd-38cf-40ea-8213-45c2b100ac01
        container_id, crdc_instance_uuid = filename.split('.')

        print(container_id)
        print(crdc_instance_uuid)

        selection_query = f"""
        SELECT crdc_instance_uuid,
               crdc_study_uuid,
               crdc_series_uuid,
               gcs_url,
               SOPInstanceUID,
               StudyInstanceUID
        FROM
          bigquery-public-data.idc_current.dicom_all
        WHERE
            ContainerIdentifier='{container_id}'
        """
        # AND crdc_instance_uuid='{crdc_instance_uuid}'

        selection_result = bq_client.query(selection_query)
        selection_df = selection_result.result().to_dataframe()

        # Choose an instance uid (TODO do this properly)
        ins_uuid = selection_df.crdc_instance_uuid.iloc[0]

        # Download the annotation
        ann_bytes = ann_blob.download_as_bytes()

        # Download the DICOM file
        dcm_blob = public_bucket.get_blob(f'{ins_uuid}.dcm')
        dcm_bytes = dcm_blob.download_as_bytes()
        dcm = pydicom.dcmread(BytesIO(dcm_bytes))

        # Read the annotations from the tarfile of CSVs
        with tarfile.open(fileobj=BytesIO(ann_bytes), mode='r:gz') as tar:
            for member in tar.getmembers():
                if member.isfile():
                    f = tar.extractfile(member)
                    df = pd.read_csv(f)
                    print(df.columns)

        break


if __name__ == "__main__":
    run()
