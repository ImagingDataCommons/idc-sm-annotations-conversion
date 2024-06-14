# IDC Annotation Conversion

Python project for converting various pathology annotations into DICOM
format for ingestion into the Imaging Data Commons.

The code in this repository is currently under development.

### Installation

This repository is structured to be directly installable as a Python
distribution named `idc-annotation-conversion` via pip. You should be able to
run this command from the root of the cloned repository to install the packages
along with all its dependencies (defined in `pyproject.toml`) in your current
Python environment:

```bash
pip install .
```
Alternatively, you can install the package directly from remote with:

```bash
pip install https://github.com/ImagingDataCommons/idc-sm-annotations-conversion.git
```

### Cloud Authentication

You need to authenticate to the relevant Google cloud buckets to run the code
in this package. Specifically, access to the following resources is required:

- Project `idc-etl-processing`
- Bucket `public-datasets-idc`, the public bucket containing DICOM-format whole
  slide images.
- Bucket `idc-annotation-conversion-outputs`, or any other bucket specified
  as the output bucket, if any.

Depending on the conversion process that you are running, you may also need
access to:

- Bucket `tcia-nuclei-seg`, which contains the original (CSV format)
  segmentations for the `pan_cancer_nuclei_seg` conversion process.
- Project `idc-external-031` and bucket `rms_annotation_test_oct_2023`, which contains the
  original (XML format) annotations for the `rms` conversion process.

If you are using an IDC cloud VM, this should be handled
automatically for you. Otherwise, you should run:

```
gcloud auth application-default login --billing-project idc-etl-processing
```

and then once you are finished:

```
gcloud auth application-default revoke
```

### Use

Each conversion process is implemented as a submodule of the `idc_annotation_conversion`
module, which is installed when you installed this package. Each submodule has an
an entrypoint (a `__main__.py` file), meaning that to run the process once this
package is installed you run:

```bash
python -m idc_annotation_conversion.<module> <args>
```

So for example to run the `pan_cancer_nuclei_seg` conversion process:

```bash
python -m idc_annotation_conversion.pan_cancer_nuclei_seg <args>
```

In each case, the default parameters should be sufficient to run a conversion processon
on the entire collection but there a number of optional arguments to control the process.
You can see the options by running `--help` when calling the submodule. E.g.:

```bash
python -m idc_annotation_conversion.pan_cancer_nuclei_seg --help
```

### Modules

The following modules are currently available:

- `pan_cancer_nuclei_seg`: Conversion of Pan Cancer Nuclei segmentations from
  XML to ANN and SEGs for various TCGA collections.
- `rms`: Conversion of annotations related to the "RMS-Mutation-Prediction"
  collection. Specifically conversion of hand annotated regions to SR, and
  ML generated segmentations to SEG.
