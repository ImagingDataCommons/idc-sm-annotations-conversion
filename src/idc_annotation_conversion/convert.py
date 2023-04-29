"""Utilities for converting annotations. Clear of cloud-specific things."""
import logging
from os import PathLike
import multiprocessing as mp
from tempfile import TemporaryDirectory
from typing import Iterable, Optional, Tuple, Union

import dicomslide
import highdicom as hd
import mahotas as mh
import numpy as np
import pandas as pd
from dicomweb_client import DICOMfileClient
from pydicom import Dataset
from pydicom.sr.codedict import codes
from pydicom.uid import JPEGLSLossless, ExplicitVRLittleEndian
from shapely.geometry.polygon import Polygon
from idc_annotation_conversion import metadata_config


def disassemble_total_pixel_matrix(
    seg_total_pixel_matrix: np.ndarray,
    source_image_metadata: Dataset,
) -> np.ndarray:
    """Disassemble a total pixel matrix into individual tiles.

    Parameters
    ----------
    seg_total_pixel_matrix: numpy.ndarray
        Total pixel matrix of the segmentation as a 2D NumPy array.
    source_metadata: pydicom.Dataset
        DICOM metadata of the source image.

    Returns
    -------
    numpy.ndarray
        Stacked image tiles

    """
    if seg_total_pixel_matrix.ndim != 2:
        raise ValueError(
            "Total pixel matrix has unexpected number of dimensions."
        )

    # Need a client object to work with DICOM slide so just create a dummy one
    with TemporaryDirectory() as tmpdir:
        client = DICOMfileClient(f"file://{tmpdir}")

        im_tpm = dicomslide.TotalPixelMatrix(client, source_image_metadata)

        tile_rows, tile_cols, _ = im_tpm.tile_shape

        return dicomslide.disassemble_total_pixel_matrix(
            seg_total_pixel_matrix,
            im_tpm.tile_positions,
            tile_rows,
            tile_cols,
        )


def process_csv_row(
    csv_row: pd.Series,
    transformer: hd.spatial.ImageToReferenceTransformer,
    store_boundary: bool = True,
    annotation_coordinate_type: hd.ann.AnnotationCoordinateTypeValues = hd.ann.AnnotationCoordinateTypeValues.SCOORD,  # noqa: E501
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Process a single annotation CSV file.

    Parameters
    ----------
    csv_row: pd.Series
        Single row of a loaded annotation CSV.
    transformer: hd.spatial.ImageToReferenceTransformer
        Transformer object to map image coordinates to reference coordinates
        for the image.
    store_boundary: bool, optional
        Store the full nucleus boundary polygon in the Bulk Microscopy Bulk
        Simple Annotations. If False, just the centroid is stored as a single
        point.
    annotation_coordinate_type: Union[hd.ann.AnnotationCoordinateTypeValues, str], optional
        Store coordinates in the Bulk Microscopy Bulk Simple Annotations in the
        (3D) frame of reference (SCOORD3D), or the (2D) total pixel matrix
        (SCOORD, default).

    Returns
    -------
    segmentation_coordinates: np.ndarray
        Numpy array of boundary coordinates in the form used by
        mahotas.polygon.fill_polygon to create the segmentation mask.
    graphic_data: np.ndarray
        Numpy array of coordinates to include in the Bulk Microscopy Simple
        Annotations.
    area: float
        Area measurement of this polygon.

    """  # noqa: E501
    points = np.array(
        csv_row.Polygon[1:-1].split(':'),
        dtype=np.float32
    )
    n = len(points) // 2
    coordinates_image = points.reshape(n, 2)
    if coordinates_image.shape[0] < 3:
        return None, None, None
    polygon_image = Polygon(coordinates_image)
    c, r = polygon_image.exterior.xy

    contour_image = np.stack([r, c]).T.astype(np.int32)

    coordinates_ref = transformer(coordinates_image)
    polygon_ref = Polygon(coordinates_ref)

    if store_boundary:
        # Store the full polygon in graphic data
        if (
            annotation_coordinate_type ==
            hd.ann.AnnotationCoordinateTypeValues.SCOORD3D
        ):
            coords = np.array(polygon_ref.exterior.coords)
        else:
            # 2D total pixel matrix coordinates
            coords = np.array(polygon_image.exterior.coords)

        # Remove the final point (polygon should not be closed)
        coords = coords[:-1, :]
        graphic_data = coords
    else:
        # Store the centroid of the polygon only, as a single point
        if (
            annotation_coordinate_type ==
            hd.ann.AnnotationCoordinateTypeValues.SCOORD3D
        ):
            x, y = polygon_ref.centroid.xy
            centroid = np.array([[x[0], y[0], 0.]])
        else:
            # 2D total pixel matrix coordinates
            x, y = polygon_image.centroid.xy
            centroid = np.array([[x[0], y[0]]])
        graphic_data = centroid

    area = float(polygon_ref.area)

    return contour_image, graphic_data, area


def pool_init(
    transformer: hd.spatial.ImageToReferenceTransformer,
    store_boundary: bool = True,
    annotation_coordinate_type: hd.ann.AnnotationCoordinateTypeValues = hd.ann.AnnotationCoordinateTypeValues.SCOORD,  # noqa: E501
):
    global transformer_global
    transformer_global = transformer
    global store_boundary_global
    store_boundary_global = store_boundary
    global annotation_coordinate_type_global
    annotation_coordinate_type_global = annotation_coordinate_type


def pool_fun(csv_row: pd.Series):
    return process_csv_row(
        csv_row,
        transformer_global,
        store_boundary_global,
        annotation_coordinate_type_global,
    )


def convert_annotations(
    annotation_csvs: Iterable[Union[str, PathLike]],
    source_image_metadata: Dataset,
    *,
    include_segmentation: bool = False,
    store_boundary: bool = True,
    annotation_coordinate_type: Union[
        hd.ann.AnnotationCoordinateTypeValues,
        str
    ] = hd.ann.AnnotationCoordinateTypeValues.SCOORD,
    segmentation_type: Union[
        hd.seg.SegmentationTypeValues,
        str
    ] = hd.seg.SegmentationTypeValues.BINARY,
    workers: int = 0,
) -> Tuple[
    hd.ann.MicroscopyBulkSimpleAnnotations,
    Optional[hd.seg.Segmentation]
]:
    """Convert an annotation into DICOM format.

    Specifically, a Bulk Microscopy Simple Annotation object (vector
    graphics) is created and a Segmentation Image (raster) is optionally
    created.

    Parameters
    ----------
    annotation_csvs: Iterable[Union[str, os.PathLike]]
        Iterable over pathlike objects, each representing the path to a
        CSV-format file containing an annotation for this image.
    source_image_metadata: pydicom.Dataset
        Pydicom dataset containing the metadata of the image (already converted
        to DICOM format). Note that this should be the metadata of the image
        at full resolution. This can be the full image dataset, but the
        PixelData attribute is not required.
    include_segmentation: bool, optional
        Include the segmentation output.
    store_boundary: bool, optional
        Store the full nucleus boundary polygon in the Bulk Microscopy Bulk
        Simple Annotations. If False, just the centroid is stored as a single
        point.
    annotation_coordinate_type: Union[hd.ann.AnnotationCoordinateTypeValues, str], optional
        Store coordinates in the Bulk Microscopy Bulk Simple Annotations in the
        (3D) frame of reference (SCOORD3D), or the (2D) total pixel matrix
        (SCOORD, default).
    segmentation_type: Union[hd.seg.SegmentationTypeValues, str], optional
        Segmentation type (BINARY or FRACTIONAL) for the Segmentation Image
        (if any).
    workers: int
        Number of subprocess workers to spawn. If 0, all computation will use
        the main thread.

    Returns
    -------
    annotation: hd.ann.MicroscopyBulkSimpleAnnotations:
        DICOM bulk microscopy annotation encoding the original annotations in
        vector format.
    segmentation: Optional[hd.seg.Segmentation]:
        DICOM segmentation image encoding the original annotations in raster
        format, if requested. None otherwise.

    """  # noqa: E501
    segmentation_type = hd.seg.SegmentationTypeValues[segmentation_type]
    annotation_coordinate_type = hd.ann.AnnotationCoordinateTypeValues[
        annotation_coordinate_type
    ]

    image_orientation = source_image_metadata.ImageOrientationSlide
    origin = source_image_metadata.TotalPixelMatrixOriginSequence[0]
    image_position = (
        float(origin.XOffsetInSlideCoordinateSystem),
        float(origin.YOffsetInSlideCoordinateSystem),
        0.0,
    )
    pixel_spacing = (
        source_image_metadata.SharedFunctionalGroupsSequence[0]
        .PixelMeasuresSequence[0]
        .PixelSpacing
    )
    transformer = hd.spatial.ImageToReferenceTransformer(
        image_orientation=image_orientation,
        image_position=image_position,
        pixel_spacing=pixel_spacing,
    )

    graphic_type = (
        hd.ann.GraphicTypeValues.POLYGON
        if store_boundary
        else hd.ann.GraphicTypeValues.POINT
    )
    graphic_data = []
    measurements = {
        (codes.SCT.Area, codes.UCUM.SquareMicrometer): [],
    }

    if include_segmentation:
        segmentation_mask = np.zeros(
            (
                source_image_metadata.TotalPixelMatrixRows,
                source_image_metadata.TotalPixelMatrixColumns,
            ),
            dtype=bool
        )
        logging.info(
            f"Total Pixel Matrix {segmentation_mask.shape} "
            f"{segmentation_mask.nbytes:.3g} bytes."
        )

    if workers > 0:
        # Use multiprcocessing
        # Read in all CSVs using threads
        logging.info("Reading CSV files.")
        all_dfs = [pd.read_csv(f) for f in annotation_csvs]

        # Join CSVs and process each row in parallel
        df = pd.concat(all_dfs, ignore_index=True)
        logging.info(f"Found {len(df)} annotations.")

        input_data = [row for _, row in df.iterrows()]

        with mp.Pool(
            workers,
            initializer=pool_init,
            initargs=(transformer, store_boundary, annotation_coordinate_type)
        ) as pool:
            results = pool.map(
                pool_fun,
                input_data,
            )

        for contour_image, graphic_item, area in results:
            if contour_image is None:
                # Failures due to too few points
                continue

            measurements[
                (codes.SCT.Area, codes.UCUM.SquareMicrometer)
            ].append(area)
            graphic_data.append(graphic_item)

            # Can't multithread this bit because each contour is written to the
            # same array
            if include_segmentation:
                mh.polygon.fill_polygon(contour_image, segmentation_mask)

    else:
        # Use the main thread

        for f in annotation_csvs:
            df = pd.read_csv(f)

            for _, csv_row in df.iterrows():
                contour_image, graphic_item, area = process_csv_row(
                    csv_row,
                    transformer,
                    store_boundary,
                    annotation_coordinate_type,
                )

                if contour_image is None:
                    continue

                graphic_data.append(graphic_item)

                measurements[
                    (codes.SCT.Area, codes.UCUM.SquareMicrometer)
                ].append(area)

                if include_segmentation:
                    mh.polygon.fill_polygon(contour_image, segmentation_mask)

    logging.info(f"Parsed {len(graphic_data)} annotations.")

    logging.info("Creating annotation.")
    group = hd.ann.AnnotationGroup(
        number=1,
        uid=hd.UID(),
        label=metadata_config.annotation_label,
        annotated_property_category=metadata_config.finding_category,
        annotated_property_type=metadata_config.finding_type,
        graphic_type=graphic_type,
        graphic_data=graphic_data,
        algorithm_type=hd.ann.AnnotationGroupGenerationTypeValues.AUTOMATIC,
        algorithm_identification=metadata_config.algorithm_identification,
        measurements=[
            hd.ann.Measurements(
                name=name,
                unit=unit,
                values=np.array(values),
            )
            for (name, unit), values in measurements.items()
        ]
    )
    annotations = hd.ann.MicroscopyBulkSimpleAnnotations(
        source_images=[source_image_metadata],
        annotation_coordinate_type=annotation_coordinate_type,
        annotation_groups=[group],
        series_instance_uid=hd.UID(),
        series_number=204,
        sop_instance_uid=hd.UID(),
        instance_number=1,
        manufacturer=metadata_config.manufacturer,
        manufacturer_model_name=metadata_config.manufacturer_model_name,
        software_versions=metadata_config.software_versions,
        device_serial_number=metadata_config.device_serial_number,
    )

    if include_segmentation:
        logging.info("Creating segmentation.")
        segmentation_tiles = disassemble_total_pixel_matrix(
            seg_total_pixel_matrix=segmentation_mask,
            source_image_metadata=source_image_metadata,
        )

        del segmentation_mask  # help free up some memory

        segment_description = hd.seg.SegmentDescription(
            segment_number=1,
            segment_label=metadata_config.label,
            segmented_property_category=metadata_config.finding_category,
            segmented_property_type=metadata_config.finding_type,
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
            algorithm_identification=metadata_config.algorithm_identification
        )

        # Compression method depends on what is possible given the chosen
        # segmentation type
        transfer_syntax_uid = {
            hd.seg.SegmentationTypeValues.BINARY: ExplicitVRLittleEndian,
            hd.seg.SegmentationTypeValues.FRACTIONAL: JPEGLSLossless,
        }[segmentation_type]

        segmentation = hd.seg.Segmentation(
            source_images=[source_image_metadata],
            pixel_array=segmentation_tiles,
            segmentation_type=segmentation_type,
            segment_descriptions=[segment_description],
            series_instance_uid=hd.UID(),
            series_number=20,
            sop_instance_uid=hd.UID(),
            instance_number=1,
            manufacturer=metadata_config.manufacturer,
            manufacturer_model_name=metadata_config.manufacturer_model_name,
            software_versions=metadata_config.software_versions,
            device_serial_number=metadata_config.device_serial_number,
            transfer_syntax_uid=transfer_syntax_uid,
            max_fractional_value=1,
        )
    else:
        segmentation = None

    return annotations, segmentation
