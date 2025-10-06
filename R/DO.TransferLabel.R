#' @title DO.TransferLabel
#' @description Transfers cell-type annotations from a re-annotated subset of a
#' Seurat or SCE object back to the full Seurat or SCE object. This is useful
#' when clusters have been refined or re-labeled in a subset and need to be
#' reflected in the original object.
#'
#' @param sce_object Seurat or SCE object with annotation in meta.data
#' @param Subset_obj subsetted Seurat or SCE object with re-annotated clusters
#' @param annotation_column column name in meta.data with annotation
#' @param subset_annotation column name in meta.data with annotation in the
#' subsetted object
#'
#'
#' @import Seurat
#'
#'
#' @examples
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' sce_data <- DO.TransferLabel(sce_data,
#'     sce_data,
#'     annotation_column = "annotation",
#'     subset_annotation = "annotation"
#' )
#'
#' @return Seurat or SCE Object with transfered labels
#' @export

DO.TransferLabel <- function(sce_object,
    Subset_obj,
    annotation_column,
    subset_annotation) {
    # support for single cell experiment objects
    if (methods::is(sce_object, "SingleCellExperiment")) {
        class_obj <- "SingleCellExperiment"
        sce_object <- as.Seurat(sce_object)
        Subset_obj <- as.Seurat(Subset_obj)
    } else {
        class_obj <- "Seurat"
    }

    # Get annotation with barcodes as rownames
    annotation_pre <- sce_object@meta.data[annotation_column]
    annotation_pre[[annotation_column]] <- as.character(
        annotation_pre[[annotation_column]]
    )
    annotation_subset <- Subset_obj@meta.data[subset_annotation]
    annotation_subset[[subset_annotation]] <- as.character(
        annotation_subset[[subset_annotation]]
    )

    # assign the new labels
    barcodes <- rownames(annotation_pre)[
        rownames(annotation_pre) %in% rownames(annotation_subset)
    ]

    annotation_pre[
        rownames(annotation_pre) %in% barcodes,
    ] <- annotation_subset[[subset_annotation]]

    sce_object@meta.data[[annotation_column]] <-
        factor(annotation_pre$annotation)

    if (class_obj == "SingleCellExperiment") {
        sce_object <- as.SingleCellExperiment(sce_object)
    }

    return(sce_object)
}
