
#' @author Mariano Ruz Jurado
#' @title DO Celltypist
#' @description Runs the CellTypist model on a Seurat object to predict cell type labels, storing the results as metadata.
#' If the number of cells is less than the specified threshold, it returns NAs for the labels.
#' Optionally updates the CellTypist models and returns the probability matrix.
#' Useful for annotating cell types in single-cell RNA sequencing datasets.
#' @param Seu_object The seurat object
#' @param minCellsToRun If the input seurat object has fewer than this many cells, NAs will be added for all expected columns and celltypist will not be run.
#' @param runCelltypistUpdate If true, --update-models will be run for celltypist prior to scoring cells.
#' @param over_clustering Column in metadata in object with clustering assignments for cells, default seurat_clusters
#' @param assay_normalized Assay with log1p normalized expressions
#' @param returnProb will additionally return the probability matrix, return will give a list with the first element beeing the object and second prob matrix
#'
#' @return a seurat
#'
#' @examples
#' \dontrun{
#'
#'
#' DO.CellTypist(
#'   Seu_object = Seurat,
#'   modelName = "Healthy_Adult_Heart.pkl",
#'   runCelltypistUpdate = TRUE,
#'   over_clustering = "seurat_clusters",
#'   SeuV5=T
#' )
#' }
#'
#' @export
DO.CellTypist <- function(Seu_object,
                          modelName = "Healthy_Adult_Heart.pkl",
                          minCellsToRun = 200,
                          runCelltypistUpdate = TRUE,
                          over_clustering = "seurat_clusters",
                          assay_normalized = "RNA",
                          returnProb=FALSE,
                          SeuV5=T) {

  # Make sure R Zellkonverter package is installed
  zk <- system.file(package = "zellkonverter")
  ifelse(nzchar(zk), "", stop("Install zellkonverter R package for Seurat tranformation!"))

  # Make sure R reticulate package is installed
  rt <- system.file(package = "reticulate")
  ifelse(nzchar(rt), "", stop("Install reticulate R package for Python usage in R!"))

  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not correctly configured. Run "usethis::edit_r_environ()" to specify your Python instance'))
  }

  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed in this python environment!')
  }

  if (ncol(Seu_object) < minCellsToRun) {
    warning('Too few cells, will not run celltypist. NAs will be added instead')
    expectedCols <- c('predicted_labels_celltypist')
    Seu_object[[expectedCols]] <- NA
    return(Seu_object)
  }

  message(paste0('Running celltypist using model: ', modelName))

  #Create temporary folder
  outDir <- tempfile(fileext = '')
  if (endsWith(outDir, "/")){
    outDir <- gsub(outDir, pattern = "/$", replacement = "")
  }
  dir.create(outDir)
  message(paste0("Saving celltypist results to temporary folder: ", outDir))

  #Uppercasing gene names
  #zellkonverter h5ad
  DefaultAssay(Seu_object) <- assay_normalized
  if (SeuV5 == TRUE) {
    tmp.assay <- Seu_object
    tmp.assay[["RNA"]] <- as(tmp.assay[["RNA"]], Class = "Assay")
    tmp.sce <- Seurat::as.SingleCellExperiment(tmp.assay, assay = assay_normalized)
    rownames(tmp.sce) <- toupper(rownames(tmp.sce))

  } else{
    tmp.sce <- Seurat::as.SingleCellExperiment(Seu_object, assay = assay_normalized)
    rownames(tmp.sce) <- toupper(rownames(tmp.sce))
  }
  zellkonverter::writeH5AD(tmp.sce, file = paste0(outDir,"/ad.h5ad"), X_name = "logcounts")

  # Ensure models present:
  if (runCelltypistUpdate) {
    system2(reticulate::py_exe(), c("-m", "celltypist.command_line", "--update-models", "--quiet"))
  }

  # Run celltypist:
  args <- c("-m", "celltypist.command_line", "--indata",  paste0(outDir,"/ad.h5ad"), "--model", modelName, "--outdir", outDir,"--majority-voting", "--over-clustering", over_clustering)
  system2(reticulate::py_exe(), args)

  labelFile <- paste0(outDir, "/predicted_labels.csv")
  probFile <- paste0(outDir, "/probability_matrix.csv")
  labels <- utils::read.csv(labelFile, header = T, row.names = 1, stringsAsFactors = FALSE)

  #ad <- import(module = "anndata")
  #ad_obj <- ad$AnnData(X = labels)

  #ct <- import(module = "celltypist")
  #ct$dotplot(ad_obj, use_as_reference = "cell_type", use_as_prediction = "majority_voting")

  probMatrix <- utils::read.csv(probFile, header = T, row.names = 1, stringsAsFactors = FALSE)
  Seu_object@meta.data$predicted_labels_celltypist <- labels$majority_voting
  if (returnProb==TRUE) {
    returnProb <- list(Seu_object, probMatrix)
    names(returnProb) <- c("SeuratObject", "probMatrix")
    return(returnProb)
  } else {
    return(Seu_object)}
}


# Recluster function using FindSubCluster function from Seurat iterative over each cluster -> perfect for fine tuning annotation
#' @author Mariano Ruz Jurado
#' @title DO.FullRecluster
#' @description Performs iterative reclustering on each major cluster found by FindClusters in a Seurat object.
#' It refines the clusters using the FindSubCluster function for better resolution and fine-tuned annotation.
#' The new clustering results are stored in a metadata column called `seurat_Recluster`.
#' Suitable for improving cluster precision and granularity after initial clustering.
#' @param Seu_object The seurat object
#' @param res Resolution for the new clusters, default 0.5
#' @param algorithm Set one of the available algorithms found in FindSubCLuster function, default = 4: leiden
#' @return a Seurat Object with new clustering named seurat_Recluster
#'
#' @import Seurat
#' @import progress
#'
#' @examples
#' \dontrun{
#'
#'
#' DO.FullRecluster(
#'   Seu_object = Seurat
#' )
#' }
#'
#' @export
DO.FullRecluster <- function(Seu_object,
                             res = 0.5,
                             algorithm=4,
                             graph.name="RNA_snn"){

  if (is.null(Seu_object$seurat_clusters)) {
    stop("No seurat clusters defined, please run FindClusters before Reclustering, or fill the slot with a clustering")
  }
  Idents(Seu_object) <- "seurat_clusters"

  Seu_object$seurat_Recluster <- as.vector(Seu_object$seurat_clusters)
  pb <- progres::progress_bar$new(total = length(unique(Seu_object$seurat_clusters)))
  for (cluster in unique(Seu_object$seurat_clusters)) {
    pb$tick()
    Seu_object <- FindSubCluster(Seu_object,
                                   cluster = as.character(cluster),
                                   graph.name = graph.name,
                                   algorithm = algorithm,
                                   resolution = res)

    cluster_cells <- rownames(Seu_object@meta.data)[Seu_object$seurat_clusters == cluster]
    Seu_object$seurat_Recluster[cluster_cells] <- Seu_object$sub.cluster[cluster_cells]
  }
  Seu_object$sub.cluster <- NULL
  return(Seu_object)
}

# Polished UMAP function using Dimplot or FeaturePlot function from Seurat
#' @author Mariano Ruz Jurado
#' @title DO.UMAP
#' @description Creates a polished UMAP plot using Seurat's DimPlot or FeaturePlot functions.
#' It allows customization of colors, labels, and other plot elements for better visualization.
#' The function handles both cluster-based visualizations and gene-based visualizations in a UMAP plot.
#' Ideal for refining UMAP outputs with added flexibility and enhanced presentation.
#' @param Seu_object The seurat object
#' @param FeaturePlot Is it going to be a Dimplot or a FeaturePlot?
#' @param features features for Featureplot
#' @param group.by grouping of plot in DImplot and defines in featureplot the labels
#' @param ... Further arguments passed to DimPlot or FeaturePlot function from Seurat
#' @return Plot with Refined colors and axes
#'
#' @import Seurat
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#'
#'
#' DO.UMAP(
#'   Seu_object = Seurat,
#'   group.by="seurat_clusters"
#' )
#'
#' DO.UMAP(
#'   Seu_object = Seurat,
#'   FeaturePlot=T,
#'   features=c("CDH5","TTN")
#' )
#' }
#'
#' @export
DO.UMAP <- function(Seu_object,
                    FeaturePlot=F,
                    features=NULL,
                    group.by="seurat_clusters",
                    umap_colors=NULL,
                    text_size=14,
                    label=T,
                    order=T,
                    plot.title=T,
                    legend.position="none",
                    ...){

  #Dimplot
  if (FeaturePlot==F) {
    if (is.null(umap_colors)) {
      umap_colors <- rep(c(
        "#1f77b4", "#ff7f0e", "#2ca02c", "tomato2", "#9467bd", "chocolate3","#e377c2", "#ffbb78", "#bcbd22",
        "#17becf","darkgoldenrod2", "#aec7e8", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94","#f7b6d2", "#c7c7c7", "#dbdb8d",
        "#9edae5","sandybrown","moccasin","lightsteelblue","darkorchid","salmon2","forestgreen","bisque"
      ),5)
    }

    p <- DimPlot(Seu_object, group.by = group.by, cols = umap_colors, ...) +
      labs(x="UMAP1",y = "UMAP2")+
      theme(plot.title = element_blank(),
            # text = element_text(face = "bold",size = 20),
            axis.title.x = element_text(size = text_size, family="Helvetica"),
            axis.title.y = element_text(size = text_size, family="Helvetica"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = legend.position,
            legend.text = element_text(face = "bold"))

    if (label==T) {
      p <-LabelClusters(p, id  = group.by, fontface="bold", box = F)
    }
    return(p)
  }

  #FeaturePlot
  if (FeaturePlot==T) {

    if (is.null(features)) {
      stop("Please provide any gene names if using FeaturePlot=T.")
    }

    if (is.null(umap_colors)) {
      umap_colors <- c("lightgrey","red2")
    }

    Idents(Seu_object) <- group.by
    p <- FeaturePlot(Seu_object,
                     features = features,
                     cols = umap_colors,
                     label = label,
                     order = order,
                     ...)&
      labs(x="UMAP1",y = "UMAP2")&
      theme(axis.title.x = element_text(size = 14, family="Helvetica"),
            axis.title.y = element_text(size = 14, family="Helvetica"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = legend.position,
            legend.text = element_text(face = "bold"))

    if (plot.title == F) {
      p <- p & theme(plot.title = element_blank())
    }

    return(p)
  }
}


# Seurat Subset function that actually works
#' @author Mariano Ruz Jurado
#' @title DO.Subset
#' @description Creates a subset of a Seurat object based on either categorical or numeric thresholds in metadata.
#' Allows for subsetting by specifying the ident column, group name, or threshold criteria.
#' Ideal for extracting specific cell populations or clusters based on custom conditions.
#' Returns a new Seurat object containing only the subsetted cells and does not come with the Seuratv5 subset issue
#' Please be aware that right now, after using this function the subset might be treated with Seuv5=False in other functions.
#' @param Seu_object The seurat object
#' @param assay assay to subset by
#' @param ident meta data column to subset for
#' @param ident_name name of group of barcodes in ident of subset for
#' @param ident_thresh numeric thresholds as character, e.g ">5" or c(">5", "<200"), to subset barcodes in ident for
#' @return a subsetted Seurat object
#'
#' @import Seurat
#' @import SeuratObject
#'
#' @examples
#' \dontrun{
#'
#'
#' DO.Subset(
#'   Seu_object = Seurat,
#'   ident="condition",
#'   ident_name="CTRL"
#' )
#'
#' DO.Subset(
#'   Seu_object = Seurat,
#'   ident="nFeature_RNA",
#'   ident_thresh=c(">5", "<200")
#' )
#' }
#'
#'
#' @export
DO.Subset <- function(Seu_object,
                      assay="RNA",
                      ident,
                      ident_name=NULL,
                      ident_thresh=NULL){

  reduction_names <- names(Seu_object@reductions)
  SCE_Object <- as.SingleCellExperiment(Seu_object)

  if (!is.null(ident_name) && !is.null(ident_thresh))  {
    stop("Please provide ident_name for subsetting by a name in the column or ident_thresh if it by a threshold")
  }

  #By a name in the provided column
  if (!is.null(ident_name) && is.null(ident_thresh))  {
    cat("Specified 'ident_name': expecting a categorical variable.\n")
    SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] %in% ident_name]
  }

  #By a threshold in the provided column
  if (is.null(ident_name) && !is.null(ident_thresh))  {
    cat(paste0("Specified 'ident_thresh': expecting numeric thresholds specified as character, ident_thresh = ", paste0(ident_thresh, collapse = " "),"\n"))

    #Extract the numeric value and operator
    operator <- gsub("[0-9.]", "", ident_thresh)
    threshold <- as.numeric(gsub("[^0-9.]", "", ident_thresh))

    #operator is valid?
    if ("TRUE" %in% !(operator %in% c("<", ">", "<=", ">="))) {
      stop("Invalid threshold operator provided. Use one of '<', '>', '<=', '>='")
    }

    #solo case
    if (length(operator) ==1) {
      if (operator == "<") {
        # filtered_cells <- ident_values[ident_values < threshold]
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] < threshold]
      } else if (operator == ">") {
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] > threshold]
      } else if (operator == "<=") {
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] <= threshold]
      } else if (operator == ">=") {
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] >= threshold]
      }
    }

    #second case
    if (length(operator) == 2) {
      if (paste(operator,collapse = "") ==  "><") {
        # filtered_cells <- ident_values[ident_values > threshold[1] & ident_values < threshold[2]]
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] > threshold[1] &
                                       SingleCellExperiment::colData(SCE_Object)[, ident] < threshold[2]]
      } else if (paste(operator,collapse = "") ==  "<>") {
        SCE_Object_sub <- SCE_Object[, SingleCellExperiment::colData(SCE_Object)[, ident] < threshold[1] &
                                       SingleCellExperiment::colData(SCE_Object)[, ident] > threshold[2]]
      }
    }
  }

  if (ncol(SCE_Object_sub) == 0) {
    stop("No cells left after subsetting!\n")
  }

  Seu_object_sub <- as.Seurat(SCE_Object_sub)
  Seu_object_sub[[assay]] <- as(object = Seu_object_sub[[assay]], Class = "Assay5")

  #Identify reductions that are not in uppercase -> These are the old ones, not subsetted REMOVE
  # non_uppercase_reductions <- reduction_names[!grepl("^[A-Z0-9_]+$", reduction_names)]
  # Seu_object_sub@reductions <- Seu_object_sub@reductions[!names(Seu_object_sub@reductions) %in% non_uppercase_reductions]
  names(Seu_object_sub@reductions) <- reduction_names
  Seu_object_sub$ident <- NULL

  #some checks
  ncells_interest_prior <- nrow(Seu_object@meta.data[Seu_object@meta.data[[ident]] %in% ident_name, ])
  ncells_interest_after <- nrow(Seu_object_sub@meta.data[Seu_object_sub@meta.data[[ident]] %in% ident_name, ])
  if (ncells_interest_prior != ncells_interest_after) {
    stop(paste0("Number of subsetted cell types is not equal in both objects! Before: ",ncells_interest_prior,"; After: ", ncells_interest_after))
  }

  return(Seu_object_sub)
}



#' @title Remove Layers from Seurat Object by Pattern
#' @description This function removes layers from a Seurat object's RNA assay based on a specified regular expression pattern.
#' It is supposed to make something similar than the no longer working DietSeurat function, by removing no longer needed layers from th object.
#' @param Seu_object Seurat object.
#' @param pattern regular expression pattern to match layer names. Default "^scale\\.data\\."
#' @return Seurat object with specified layers removed.
#'
#' @import SeuratObject
#'
#' @examples
#' \dontrun{
#'
#'
#' DO.Subset(
#'   Seu_object = Seurat,
#'   ident="condition",
#'   ident_name="CTRL"
#' )
#'
#' DO.Subset(
#'   Seu_object = Seurat,
#'   ident="nFeature_RNA",
#'   ident_thresh=c(">5", "<200")
#' )
#' }
#'
#' @export
DO.DietSeurat <- function(Seu_object, pattern = "^scale\\.data\\.") {
  message(paste("pattern: ", pattern))
  stopifnot("object must be a Seurat object" = inherits(Seu_object, "Seurat"))

  layers_to_remove <- grep(pattern, Layers(Seu_object), value = TRUE)
  Seu_object@assays$RNA@layers[layers_to_remove] <- NULL

  layerNames <- Layers(Seu_object)
  message(paste(layers_to_remove, "is removed."))
  return(Seu_object)
}


umap_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "tomato2", "#9467bd", "chocolate3","#e377c2", "#ffbb78", "#bcbd22",
  "#17becf","darkgoldenrod2", "#aec7e8", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94","#f7b6d2", "#c7c7c7", "#dbdb8d",
  "#9edae5","sandybrown","moccasin","lightsteelblue","darkorchid","salmon2","forestgreen","bisque"
)
