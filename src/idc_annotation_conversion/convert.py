"""Utilities for converting annotations. Clear of cloud-specific things."""
import logging
from os import PathLike
import multiprocessing as mp
from time import time
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import highdicom as hd
import numpy as np
import pandas as pd
from pydicom import Dataset
from pydicom.sr.codedict import codes
from pydicom.uid import JPEGLSLossless, ExplicitVRLittleEndian
from rasterio.features import rasterize
from shapely.geometry.polygon import Polygon


from idc_annotation_conversion import metadata_config


def process_csv_row(
    csv_row: pd.Series,
    transformer: hd.spatial.ImageToReferenceTransformer,
    area_per_pixel_um2: float,
    store_boundary: bool = True,
    annotation_coordinate_type: hd.ann.AnnotationCoordinateTypeValues = hd.ann.AnnotationCoordinateTypeValues.SCOORD,  # noqa: E501
) -> Tuple[Polygon, np.ndarray, float]:
    """Process a single annotation CSV file.

    Parameters
    ----------
    csv_row: pd.Series
        Single row of a loaded annotation CSV.
    transformer: hd.spatial.ImageToReferenceTransformer
        Transformer object to map image coordinates to reference coordinates
        for the image.
    area_per_pixel_um2: float
        Area of each pixel in square micrometers.
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
    polygon_image: shapely.Polygon
        Polygon (in image coordinates) representing the annotation in the CSV.
    graphic_data: np.ndarray
        Numpy array of float32 coordinates to include in the Bulk Microscopy
        Simple Annotations.
    area: float
        Area measurement of this polygon.

    """  # noqa: E501
    points = np.array(
        csv_row.Polygon[1:-1].split(':'),
        dtype=np.float32
    )
    area_pix = float(csv_row.AreaInPixels)
    area_um2 = area_pix * area_per_pixel_um2
    n = len(points) // 2
    coordinates_image = points.reshape(n, 2)
    if coordinates_image.shape[0] < 3:
        return None, None, None
    polygon_image = Polygon(coordinates_image)

    use_3d = (
        annotation_coordinate_type ==
        hd.ann.AnnotationCoordinateTypeValues.SCOORD3D
    )
    if use_3d:
        coordinates_ref = transformer(coordinates_image)
        polygon_ref = Polygon(coordinates_ref)

    if store_boundary:
        # Store the full polygon in graphic data
        if use_3d:
            coords = np.array(polygon_ref.exterior.coords)
        else:
            # 2D total pixel matrix coordinates
            coords = np.array(polygon_image.exterior.coords)

        # Remove the final point (polygon should not be closed)
        coords = coords[:-1, :]
        graphic_data = coords
    else:
        # Store the centroid of the polygon only, as a single point
        if use_3d:
            x, y = polygon_ref.centroid.xy
            centroid = np.array([[x[0], y[0], 0.]])
        else:
            # 2D total pixel matrix coordinates
            x, y = polygon_image.centroid.xy
            centroid = np.array([[x[0], y[0]]])
        graphic_data = centroid

    return polygon_image, graphic_data.astype(np.float32), area_um2


def pool_init(
    transformer: hd.spatial.ImageToReferenceTransformer,
    area_per_pixel_um2: float,
    store_boundary: bool = True,
    annotation_coordinate_type: hd.ann.AnnotationCoordinateTypeValues = hd.ann.AnnotationCoordinateTypeValues.SCOORD,  # noqa: E501
):
    global transformer_global
    transformer_global = transformer
    global area_per_pixel_um2_global
    area_per_pixel_um2_global = area_per_pixel_um2
    global store_boundary_global
    store_boundary_global = store_boundary
    global annotation_coordinate_type_global
    annotation_coordinate_type_global = annotation_coordinate_type


def pool_fun(csv_row: pd.Series):
    return process_csv_row(
        csv_row,
        transformer_global,
        area_per_pixel_um2_global,
        store_boundary_global,
        annotation_coordinate_type_global,
    )


def convert_annotations(
    annotation_csvs: Iterable[Union[str, PathLike]],
    source_images: Sequence[Dataset],
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
    dimension_organization_type: hd.DimensionOrganizationTypeValues = hd.DimensionOrganizationTypeValues.TILED_FULL,  # noqa: E501
    create_pyramid: bool = True,
) -> Tuple[
    hd.ann.MicroscopyBulkSimpleAnnotations,
    Optional[List[hd.seg.Segmentation]]
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
    source_images: Sequence[pydicom.Dataset]
        List of pydicom datasets containing the metadata of the image (already
        converted to DICOM format). Note that the metadata of the image at full
        resolution should appear first in this list. These can be the full image
        datasets, but the PixelData attributes are not required.
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
    dimension_organization_type: Union[hd.DimensionOrganizationTypeValues, str], optional
        Dimension organization type of the output segmentations.
    create_pyramid: bool, optional
        Whether to create a full pyramid of segmentations (rather than a single
        segmentation at the highest resolution).

    Returns
    -------
    annotation: hd.ann.MicroscopyBulkSimpleAnnotations:
        DICOM bulk microscopy annotation encoding the original annotations in
        vector format.
    segmentation: Optional[List[hd.seg.Segmentation]]:
        DICOM segmentation image(s) encoding the original annotations in raster
        format, if requested. None otherwise.

    """  # noqa: E501
    segmentation_type = hd.seg.SegmentationTypeValues[segmentation_type]
    annotation_coordinate_type = hd.ann.AnnotationCoordinateTypeValues[
        annotation_coordinate_type
    ]

    source_image_metadata = source_images[0]
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
    area_per_pixel_um2 = pixel_spacing[0] * pixel_spacing[1] * 1e6

    graphic_type = (
        hd.ann.GraphicTypeValues.POLYGON
        if store_boundary
        else hd.ann.GraphicTypeValues.POINT
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
            initargs=(
                transformer,
                area_per_pixel_um2,
                store_boundary,
                annotation_coordinate_type,
            )
        ) as pool:
            results = pool.map(
                pool_fun,
                input_data,
            )

        results = [r for r in results if r[0] is not None]  # remove failures

        polygons, graphic_data, areas = zip(*results)

        measurements = {
            (codes.SCT.Area, codes.UCUM.SquareMicrometer): areas
         }

    else:
        # Use the main thread
        graphic_data = []
        measurements = {
            (codes.SCT.Area, codes.UCUM.SquareMicrometer): [],
        }
        for f in annotation_csvs:
            df = pd.read_csv(f)

            for _, csv_row in df.iterrows():
                contour_image, graphic_item, area = process_csv_row(
                    csv_row,
                    transformer,
                    area_per_pixel_um2,
                    store_boundary,
                    annotation_coordinate_type,
                )

                if contour_image is None:
                    continue

                graphic_data.append(graphic_item)

                measurements[
                    (codes.SCT.Area, codes.UCUM.SquareMicrometer)
                ].append(area)

    logging.info(f"Parsed {len(graphic_data)} annotations.")

    logging.info("Creating annotation.")
    group = hd.ann.AnnotationGroup(
        number=1,
        uid=hd.UID(),
        label=metadata_config.label,
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
        logging.info("Rasterizing segmentation.")
        mask_shape = (
            source_image_metadata.TotalPixelMatrixRows,
            source_image_metadata.TotalPixelMatrixColumns,
        )
        raster_start_time = time()
        segmentation_mask = rasterize(polygons, mask_shape, dtype=np.uint8)
        raster_time = time() - raster_start_time
        logging.info(f"Completed rasterization in {raster_time:.1f}s.")

        logging.info("Creating DICOM segmentation image.")

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

        omit_empty_frames = dimension_organization_type != "TILED_FULL"

        if create_pyramid:
            seg_start_time = time()
            segmentations = hd.seg.pyramid.create_segmentation_pyramid(
                source_images=source_images,
                pixel_arrays=[segmentation_mask],
                segmentation_type=segmentation_type,
                segment_descriptions=[segment_description],
                series_instance_uid=hd.UID(),
                series_number=20,
                manufacturer=metadata_config.manufacturer,
                manufacturer_model_name=metadata_config.manufacturer_model_name,
                software_versions=metadata_config.software_versions,
                device_serial_number=metadata_config.device_serial_number,
                transfer_syntax_uid=transfer_syntax_uid,
                max_fractional_value=1,
                dimension_organization_type=dimension_organization_type,
                omit_empty_frames=omit_empty_frames,
            )
            seg_time = time() - seg_start_time
            logging.info(f"Created DICOM Segmentations in {seg_time:.1f}s.")
        else:
            seg_start_time = time()
            segmentation = hd.seg.Segmentation(
                source_images=[source_image_metadata],
                pixel_array=segmentation_mask,
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
                tile_pixel_array=True,
                dimension_organization_type=dimension_organization_type,
                omit_empty_frames=omit_empty_frames,
            )
            segmentations = [segmentation]
            seg_time = time() - seg_start_time
            logging.info(f"Created DICOM Segmentation in {seg_time:.1f}s.")
    else:
        segmentations = None

    return annotations, segmentations
