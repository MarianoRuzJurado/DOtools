# SCE Subset function that actually works
#' @author Mariano Ruz Jurado
#' @title DO.Subset
#' @description Creates a subset of a Seurat or SCE object based on either
#' categorical or numeric thresholds in metadata. Allows for subsetting by
#' specifying the ident column, group name, or threshold criteria. Ideal for
#' extracting specific cell populations or clusters based on custom conditions.
#' Returns a new Seurat or SCE object containing only the subsetted cells and
#' does not come with the Seuratv5 subset issue. Please be aware that right now,
#' after using this function the subset might be treated with Seuv5=False in
#' other functions.
#' @param sce_object The seurat or SCE object
#' @param assay assay to subset by
#' @param ident meta data column to subset for
#' @param ident_name name of group of barcodes in ident of subset for
#' @param ident_thresh numeric thresholds as character,
#'  e.g ">5" or c(">5", "<200"), to subset barcodes in ident
#' @return a subsetted Seurat or SCE object
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom SingleCellExperiment reducedDim
#' @examples
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' sce_data_sub <- DO.Subset(
#'     sce_object = sce_data,
#'     ident = "condition",
#'     ident_name = "healthy"
#' )
#'
#' @export
DO.Subset <- function(sce_object,
    assay = "RNA",
    ident,
    ident_name = NULL,
    ident_thresh = NULL) {
    # support for Seurat objects
    if (methods::is(sce_object, "Seurat")) {
        Seu_obj <- sce_object
        class_obj <- "Seurat"
        reduction_names <- names(sce_object@reductions)
        sce_object <- as.SingleCellExperiment(sce_object)
    } else {
        class_obj <- "SingleCellExperiment"
    }

    if (!is.null(ident_name) && !is.null(ident_thresh)) {
        stop(
            "Please provide ident_name for subsetting by a name in the ",
            "column or ident_thresh if it by a threshold"
        )
    }

    # By a name in the provided column
    if (!is.null(ident_name) && is.null(ident_thresh)) {
        .logger("Specified 'ident_name': expecting a categorical variable.")
        sce_object_sub <- sce_object[
            , SingleCellExperiment::colData(sce_object)[, ident] %in% ident_name
        ]
    }

    # By a threshold in the provided column
    if (is.null(ident_name) && !is.null(ident_thresh)) {
        .logger(
            paste0(
                "Specified 'ident_thresh': expecting numeric thresholds ",
                "specified as character, ident_thresh = ",
                paste0(ident_thresh, collapse = " ")
            )
        )

        # Extract the numeric value and operator
        operator <- gsub("[0-9.]", "", ident_thresh)
        threshold <- as.numeric(gsub("[^0-9.]", "", ident_thresh))

        # operator is valid?
        if ("TRUE" %in% !(operator %in% c("<", ">", "<=", ">="))) {
            stop(
                "Invalid threshold operator provided. Use ",
                "one of '<', '>', '<=', '>='"
            )
        }

        # solo case
        if (length(operator) == 1) {
            if (operator == "<") {
                # filtered_cells <- ident_values[ident_values < threshold]
                sce_object_sub <- sce_object[
                    , SingleCellExperiment::colData(sce_object)[, ident] <
                        threshold
                ]
            } else if (operator == ">") {
                sce_object_sub <- sce_object[
                    , SingleCellExperiment::colData(sce_object)[, ident] >
                        threshold
                ]
            } else if (operator == "<=") {
                sce_object_sub <- sce_object[
                    , SingleCellExperiment::colData(sce_object)[, ident] <=
                        threshold
                ]
            } else if (operator == ">=") {
                sce_object_sub <- sce_object[
                    , SingleCellExperiment::colData(sce_object)[, ident] >=
                        threshold
                ]
            }
        }

        # second case
        if (length(operator) == 2) {
            if (paste(operator, collapse = "") == "><") {
                # filtered_cells <- ident_values[ident_values > threshold[1] &
                #                                 ident_values < threshold[2]]
                sce_object_sub <- sce_object[
                    , SingleCellExperiment::colData(sce_object)[, ident] >
                        threshold[1] &
                        SingleCellExperiment::colData(sce_object)[, ident] <
                            threshold[2]
                ]
            } else if (paste(operator, collapse = "") == "<>") {
                sce_object_sub <- sce_object[
                    , SingleCellExperiment::colData(sce_object)[, ident] <
                        threshold[1] &
                        SingleCellExperiment::colData(sce_object)[, ident] >
                            threshold[2]
                ]
            }
        }
    }

    if (ncol(sce_object_sub) == 0) {
        stop("No cells left after subsetting!\n")
    }

    # Seurat support
    if (class_obj == "Seurat") {
        sce_object_sub <- .suppressDeprecationWarnings(
            as.Seurat(sce_object_sub)
        )
        sce_object_sub[[assay]] <- methods::as(
            object = sce_object_sub[[assay]], Class = "Assay5"
        )

        names(sce_object_sub@reductions) <- reduction_names
        sce_object_sub$ident <- NULL
        # some checks
        ncells_interest_prior <- nrow(
            Seu_obj@meta.data[
                Seu_obj@meta.data[[ident]] %in% ident_name,
            ]
        )

        ncells_interest_after <- nrow(
            sce_object_sub@meta.data[
                sce_object_sub@meta.data[[ident]] %in% ident_name,
            ]
        )
        if (ncells_interest_prior != ncells_interest_after) {
            stop(
                "Number of subsetted cell types is not equal in both objects! ",
                "Before: %s; After: %s. Please check your metadata!",
                ncells_interest_prior,
                ncells_interest_after
            )
        }
    }

    return(sce_object_sub)
}
