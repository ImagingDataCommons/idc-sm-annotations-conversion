# IDC Annotation Conversion

Python project for converting Nuclei segmentation annotations to the DICOM
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
pip install https://github.com/CPBridge/idc-annotation-conversion.git
```

### Cloud Authentication

You need to authenticate to the relevant Google cloud buckets to run the code
in this package. Specifically, access to the following resources is required:

- Project `idc-etl-processing`
- Bucket `tcia-nuclei-seg`, which contains the original (CSV format)
  segmentations.
- Bucket `public-datasets-idc`, the public bucket containing DICOM-format whole
  slide images.
- Bucket `idc-annotation-conversion-outputs`, or any other bucket specified
  as the output bucket, if any.

If you are using an IDC cloud VM, this should be handled
automatically for you. Otherwise, you should run:

```
gcloud auth application-default login --billing-project idc-etl-processing
```

and then

```
gcloud auth application-default revoke
```

### Running the Process

To run the conversion process, simply execute the `idc-annotation-conversion`
command, which should be installed in your environment by pip.

```bash
idc-annotation-conversion
```

(Alternatively `python -m idc_annotation_conversion` should also work).
Without further arguments, this will convert all collections using the default
parameters. However there are many options to control this process, you can see
them by running

```bash
idc-annotation-conversion --help
```
