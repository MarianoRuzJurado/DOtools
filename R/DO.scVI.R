#' @title DO.scVI
#' @description This function will run the scVI Integration from the scVI python package. It includes all parameters from the actual python package
#' and runs it by using an internal python script. The usage of a gpu is incorporated and highly recommended.
#'
#' @param sce_object Seurat or SCE object with annotation in meta.data
#' @param batch_key meta data column with batch information.
#' @param layer_counts layer with counts. Raw counts are required.
#' @param layer_logcounts layer with log-counts. Log-counts required for calculation of HVG.
#' @param categorical_covariates meta data column names with categorical covariates for scVI inference.
#' @param continuos_covariates meta data  column names with continuous covariates for scVI inference.
#' @param n_hidden number of hidden layers.
#' @param n_latent dimensions of the latent space.
#' @param n_layers number of layers.
#' @param dispersion dispersion mode for scVI.
#' @param gene_likelihood gene likelihood.
#' @param get_model return the trained model.
#'
#'
#' @import Seurat
#' @import zellkonverter
#' @import SingleCellExperiment
#' @importFrom basilisk basiliskRun
#'
#' @examples
#' \dontrun{
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' # Run scVI using the 'orig.ident' column as the batch key
#' sce_data <- DO.scVI(sce_data, batch_key = "orig.ident")
#'
#' }
#'
#' @return Seurat or SCE Object with dimensionality reduction from scVI
#' @export
DO.scVI <- function(sce_object,
                    batch_key,
                    layer_counts="counts",
                    layer_logcounts="logcounts",
                    categorical_covariates=NULL,
                    continuos_covariates=NULL,
                    n_hidden=128,
                    n_latent=30,
                    n_layers=3,
                    dispersion="gene-batch",
                    gene_likelihood="zinb",
                    get_model=FALSE
){

  #support for Seurat objects
  if (is(sce_object, "Seurat")) {
    class_obj <- "Seurat"
    HVG <- VariableFeatures(sce_object)
    #check if Variable Features are defined in the object
    if (length(HVG)==0) {

      .logger("No highly variable genes found in the object! Calculating...")
      HVG <- VariableFeatures(FindVariableFeatures(sce_obj))

      if(length(HVG)==0){
        stop("No highly variable genes found in the object, after automatic calculation!")
      }
    }

    sc_obj <- Seurat::as.SingleCellExperiment(sce_object)
    #Subset the object to the HVG
    sce_object_sub <- sc_obj[rownames(sc_obj) %in% HVG, ]
  } else{
    class_obj <- "SingleCellExperiment"
    Seu_obj <- as.Seurat(sce_object)
    HVG <- VariableFeatures(FindVariableFeatures(Seu_obj))
    sce_object_sub <- sce_object[rownames(sce_object) %in% HVG, ]

  }

  #Make Anndata object
  if (!"counts" %in% names(sce_object_sub@assays)) {
    stop("counts not found in assays of object!")
  }

  #source PATH to python script in install folder
  path_py <- system.file("python", "scVI.py", package = "DOtools")

  #argument list passed to heatmap inside basilisk
  args <- list(
    sce_object = sce_object_sub,
    batch_key = batch_key,
    layer_counts = layer_counts,
    layer_logcounts = layer_logcounts,
    categorical_covariates = categorical_covariates,
    continuos_covariates = continuos_covariates,
    n_hidden = as.integer(n_hidden),
    n_latent = as.integer(n_latent),
    n_layers = as.integer(n_layers),
    dispersion = dispersion,
    gene_likelihood = gene_likelihood,
    get_model = get_model
  )

  #basilisk implementation
  scvi_embedding <- basilisk::basiliskRun(env = DOtoolsEnv, fun=function(args){

    anndata_object <- zellkonverter::SCE2AnnData(args$sce_object)
    anndata_object$layers['counts'] <- anndata_object$X # set

    reticulate::source_python(path_py)

    run_scvi(adata = anndata_object,
             batch_key = args$batch_key,
             layer_counts = args$layer_counts,
             layer_logcounts = args$layer_logcounts,
             categorical_covariates = args$categorical_covariates,
             continuos_covariates = args$continuos_covariates,
             n_hidden = args$n_hidden,
             n_latent = args$n_latent,
             n_layers = args$n_layers,
             dispersion = args$dispersion,
             gene_likelihood = args$gene_likelihood,
             get_model = args$get_model
    )

    return(anndata_object$obsm[["X_scVI"]])
  }, args=args)

  rownames(scvi_embedding) <- colnames(sce_object)

  if (class_obj == "Seurat") {
    scVI_reduction <- suppressWarnings(Seurat::CreateDimReducObject(
      embeddings = scvi_embedding,
      key = "scVI_",
      assay = DefaultAssay(sce_object)
    ))
    sce_object@reductions[["scVI"]] <- scVI_reduction
  } else{
    reducedDim(sce_object, "scVI") <- scvi_embedding
  }
  return(sce_object)
}
