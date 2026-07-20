#' @title DO.scANVI
#' @description Run the scANVI annotation transfer and latent representation
#' inference using a previously trained scVI model. The function loads an
#' existing scVI model, initializes a scANVI model, trains it, optionally saves
#' the trained scANVI model, and returns the latent embedding together with
#' optional predicted cell labels in the SCE or Seurat object provided.
#'
#' @param sce_object A Seurat or SingleCellExperiment object.
#' @param model_path Path to a previously trained scVI model.
#' @param labels_key Name of the cell annotation column used during scANVI
#' training.
#' @param unlabeled_category Name of the category representing unlabeled cells.
#' @param batch_size Mini-batch size used during scANVI training.
#' @param save_model Optional path where the trained scANVI model should be
#' saved. If `NULL`, the model is not saved.
#' @param annot_predict Logical; if `TRUE`, predicted cell labels are returned
#' and stored in the output object.
#'
#'
#' @import Seurat
#' @import SingleCellExperiment
#' @importFrom basilisk basiliskRun
#'
#' @examples
#' \dontrun{
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' # Run scANVI using a previous saved scVI model and annotation categories
#' sce_data <- DO.scANVI(
#'     sce_data,
#'     model_path = "/path/scVI_model",
#'     labels_key = "annotation",
#'     unlabeled_category = "Unknown")
#' }
#'
#' @return Seurat or SCE Object with dimensionality reduction from scANVI and
#' annotation column from scANVI
#' @export
#'
DO.scANVI <- function(
    sce_object,
    model_path,
    labels_key,
    unlabeled_category,
    batch_size = 128,
    save_model = NULL,
    annot_predict = TRUE) {
      # support for Seurat objects
      if (methods::is(sce_object, "Seurat")) {
        class_obj <- "Seurat"

      } else {
        class_obj <- "SingleCellExperiment"
      }

      #arguments for scANVI
      args <- list(
        # sce_object = sce_object,
        # HVG = HVG,
        model_path = model_path,
        labels_key = labels_key,
        unlabeled_category = unlabeled_category,
        batch_size = as.integer(batch_size),
        save_model = save_model,
        annot_predict = annot_predict
      )

      # basilisk implementation
      res <- basilisk::basiliskRun(
        env = DOtoolsEnv,
        fun = function(args) {

          sc <- reticulate::import("scvi")
          ad <- reticulate::import("anndata")

          #load model
          model <- sc$model$SCVI$load(args$model_path)

          #configurate SCANVI
          model_scanvi <- sc$model$SCANVI$from_scvi_model(
            scvi_model = model,
            unlabeled_category = args$unlabeled_category,
            labels_key = args$labels_key
          )

          model_scanvi$view_anndata_setup()
          model_scanvi$train(batch_size = args$batch_size)

          if (!is.null(args$save_model)) {
            model_scanvi$save(args$save_model,
              overwrite = TRUE,
              save_anndata = TRUE)
          }

          labels_pred <- NULL
          if (args$annot_predict) {
            labels_pred = model_scanvi$predict()
          }

          res <- list(
            embedding = model_scanvi$get_latent_representation(),
            labels_pred = labels_pred
          )

          return(res)

        }, args = args
      )

      #retrieve embedding from results
      scANVI_embedding <- res$embedding
      rownames(scANVI_embedding) <- colnames(sce_object)


      #Assign the embedding and predicitons from scANVI to R objects
      if (class_obj == "Seurat") {
        scANVI_reduction <- .suppressAllWarnings(Seurat::CreateDimReducObject(
          embeddings = scANVI_embedding,
          key = "scANVI_",
          assay = DefaultAssay(sce_object)
        ))
        sce_object@reductions[["scANVI"]] <- scANVI_reduction

        if (annot_predict) {
          sce_object[["scANVI_labels"]] <- res$labels_pred
        }

      } else {
        reducedDim(sce_object, "scANVI") <- scANVI_embedding

        if (annot_predict) {
          sce_object[["scANVI_labels"]] <- res$labels_pred
        }

      }
      return(sce_object)
    }
