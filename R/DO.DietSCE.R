#' @title Remove Layers from Seurat or SCE Object by Pattern
#' @author Mariano Ruz Jurado
#' @description This function removes layers from a Seurat or SCE object's RNA
#' assay based on a specified regular expression pattern. It is supposed to
#' remove no longer needed layers from th object.
#' @param sce_object Seurat or SCE object.
#' @param assay Name of the assay from where to remove layers from
#' @param pattern regular expression pattern to match layer names. Default
#' "^scale\\.data\\."
#' @return Seurat or SCE object with specified layers removed.
#'
#' @import SeuratObject
#'
#' @examples
#' sce_data <-
#'   readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' sce_data <- DO.DietSCE(sce_data, pattern = "data")
#'
#' @export
DO.DietSCE <- function(sce_object,
                       assay = "RNA",
                       pattern = "^scale\\.data\\.") {
  .logger(paste("pattern: ", pattern))

  # support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    class_obj <- "SingleCellExperiment"
    sce_object <- as.Seurat(sce_object)
  } else {
    class_obj <- "Seurat"
  }

  layers_to_remove <- grep(pattern, Layers(sce_object), value = TRUE)

  if ("layers" %in% slotNames(sce_object@assays[[assay]])) {
    sce_object@assays[[assay]]@layers[layers_to_remove] <- NULL
    layerNames <- Layers(sce_object)
    .logger(paste(layers_to_remove, "is removed."))
  } else {
    .logger(paste0(
      "Object has no layers, pattern does not need to be removed from ",
      "layers."
    ))
  }

  if (class_obj == "SingleCellExperiment") {
    sce_object <- as.SingleCellExperiment(sce_object)
  }

  return(sce_object)
}
