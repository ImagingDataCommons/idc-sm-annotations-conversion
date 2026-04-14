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

### Collection Details

The following modules are currently available:

- `pan_cancer_nuclei_seg`: This module implements conversion of Pan Cancer
  Nuclei Segmentations for several collections within TCGA. The original data are supplied
  in a non-standard CSV format giving the image coordinates points on the contours of
  nuclei as segmented by a deep-learning based segmentation model. These data were previously released
  [here](https://www.cancerimagingarchive.net/analysis-result/pan-cancer-nuclei-seg/) as part of
  The Cancer Imaging Archive. These coordinates are converted to DICOM Microscopy Bulk Simple
  Annotation objects, and in addition, the contours are converted to masks and stored as
  a pyramid of binary DICOM Segmentation objects. Since this "raster conversion" takes place at the
  highest resolution, this process is very slow and memory intensive.

- `rms`: Conversion of annotations related to the rhabdomyosarcoma mutation prediction project from the Frederick National Laboratory.
  Both hand annotated regions for tissue type (necrosis, stroma, ARMS, ERMS), used as training data in the project, and model-generated prediction results (for the same tissue classes) are available.
  Hand annotated regions are provided as ImageScope format XML annotations and are converted to DICOM Structured Report objects with the `convert-xml-annotations` sub-command.
  Model-generated probabilistic segmentation maps are provided as serialized NumPy arrays (`.npy` files) and converted to both binary and fractional DICOM Segmentation objects with the `convert-segmentations` sub-command.

- `tcga_til_maps`: There are two versions of this collection, both containing patch-level maps of tumor-infiltrating lymphocytes (TILs) predicted by a neural network for several collections within TCGA. The two versions correspond to two different versions of the model, published in 2018 and 2022 by the same lab at Stony Brook University. Conversion routines for these two versions implemented as two separate sub-commands within this module.

  The 2018 versions covers a smaller subset of the TCGA collections. The algorithm is published in [this paper](https://www.cell.com/cell-reports/pdf/S2211-1247(18)30447-9.pdf) and the source files are available [here](https://stonybrookmedicine.app.box.com/v/cellreportspaper). The collection was also described by TCIA on [this page](https://www.cancerimagingarchive.net/analysis-result/til-wsi-tcga/). The files are supplied as low-resolution PNG images, where each pixel in the PNG corresponds to a 50 micron patch in the original slide and the pixel value indicates the presence of TILs within the patch. The `convert-2018` command converts these to binary DICOM segmentation objects.

  The 2022 versions covers a wider range of TCGA images and additionally has probabilistic segmentations (before thresholding) available in addition to binarized versions. This algorithm is described in [this paper](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2021.806603/full) and the source data are available [here](https://stonybrookmedicine.box.com/v/til-results-new-model). These are supplied as a non-standard text file containing a list of patch coordinates and associated binary or probabilistic pixel values. The `convert-2022` command coverts these to pixel arrays and stores them as DICOM Segmentation objects, giving one binary and one fractional (probabilistic) segmentation object for each slide.

- `gbm_transcriptional_subtypes`: This module relates to a collection of results from [this paper](https://www.nature.com/articles/s41467-023-39933-0) from Stanford University on transcriptional subtypes within glioblastoma. There are two data types of interest here: transcriptional subtype maps classifying an image patch into a set of transcriptional subtypes, and aggressiveness maps giving the aggressiveness of each image patch. While the conversion process for both is implemented in this repository, only the aggressiveness maps have been released at this time. The source data are not publicly available elsewhere. The aggressiveness maps are supplied as arrays of image coordinates and corresponding aggressiveness scores (between 0 and 1) within an h5 format file, with one aggressiveness score for an entire image patch. These are converted to DICOM Parametric Map objects using the `convert-aggressiveness-maps` sub-command of this module.

- `catch`: This module relates to the conversion of pathologist-drawn annotations of tumor and tissue regions in the [CAnine CuTaneous Cancer Histology](https://www.cancerimagingarchive.net/collection/catch/) dataset. Annotations are supplied in a SQLITE database created by the SlideRunner annotation tool. These are converted to DICOM Microscopy Bulk Simple Annotation objects.
