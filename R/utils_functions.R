#' @title DO.Import
#' @author Mariano Ruz Jurado & David John
#' @description
#' Imports and processes single-cell RNA-seq data from various formats (10x Genomics, CellBender, or CSV),
#' performs quality control (QC), filtering, normalization, variable gene selection, and optionally detects doublets.
#' Returns a merged and processed Seurat or SCE object ready for downstream analysis.
#'
#' @param pathways A character vector of paths to directories or files containing raw expression matrices.
#' @param ids A character vector of sample identifiers, matching the order of `pathways`.
#' @param minCellGenes Integer. Minimum number of cells a gene must be expressed in to be retained. Default is 5.
#' @param FilterCells Logical. If `TRUE`, applies QC filtering on cells based on mitochondrial content, counts, and feature thresholds. Default is `TRUE`.
#' @param cut_mt Numeric. Maximum allowed mitochondrial gene proportion per cell. Default is 0.05.
#' @param min_counts Numeric. Minimum UMI count threshold (optional, used only if `low_quantile` is `NULL`).
#' @param max_counts Numeric. Maximum UMI count threshold (optional, used only if `high_quantile` is `NULL`).
#' @param min_genes Numeric. Minimum number of genes detected per cell to retain. Optional.
#' @param max_genes Numeric. Maximum number of genes detected per cell to retain. Optional.
#' @param low_quantile Numeric. Quantile threshold (0 to 1) to filter low UMI cells (used if `min_counts` is `NULL`).
#' @param high_quantile Numeric. Quantile threshold (0 to 1) to filter high UMI cells (used if `max_counts` is `NULL`).
#' @param DeleteDoublets Logical. If `TRUE`, doublets are detected and removed using `scDblFinder`. Default is `TRUE`.
#' @param include_rbs Logical. If `TRUE`, calculates ribosomal gene content in addition to mitochondrial content. Default is `TRUE`.
#' @param Seurat Logical. If `TRUE`, returns Seurat object otherwise SCE object.
#' @param ... Additional arguments passed to `RunPCA()`.
#'
#' @return A merged Seurat or SCE object containing all samples, with normalization, QC, scaling, PCA, and optional doublet removal applied.
#'
#' @import Seurat
#' @import SeuratObject
#' @import ggplot2
#' @import ggpubr
#' @import openxlsx
#' @import dplyr
#' @import magrittr
#'
#'
#' @examples
#' \dontrun{
#' merged_obj <- DO.Import(
#'   pathways = c("path/to/sample1", "path/to/sample2"),
#'   ids = c("sample1", "sample2"),
#'   TenX = TRUE,
#'   CellBender = FALSE,
#'   minCellGenes = 5,
#'   FilterCells = TRUE,
#'   cut_mt = 0.05,
#'   min_counts = 1000,
#'   max_counts = 20000,
#'   min_genes = 200,
#'   max_genes = 6000,
#'   DeleteDoublets = TRUE
#' )
#' }
#'
#' @export
DO.Import <- function(pathways,
                      ids,
                      minCellGenes=5,
                      FilterCells=TRUE,
                      cut_mt=.05,
                      min_counts=NULL,
                      max_counts=NULL,
                      min_genes=NULL,
                      max_genes=NULL,
                      low_quantile=NULL,
                      high_quantile=NULL,
                      DeleteDoublets=TRUE,
                      include_rbs=TRUE,
                      Seurat=TRUE,
                      ...){

  object_list <- list()
  for (i in seq_along(pathways)) {

    id <- ids[i]
    pathway <- pathways[i]

    .logger(paste0("Sample: ",id))

    #cover the case when the user provides just the path to the folder with the .h5
    file_path <- ifelse(dir.exists(pathway),
                        list.files(pathway, pattern = "*filtered.h5",full.names = TRUE),
                        pathway)
    outPath <- ifelse(dir.exists(pathway),
                      pathway,
                      dirname(pathway))

    #Try different methods to read provided path
    .logger("Read matrix")
    tryCatch({
      mtx <- Seurat::Read10X_h5(file_path)
      #Capture the case of multiomics saved in the cellranger results
      if (is.list(mtx)) {
        mtx <- mtx[[grep("gene",names(mtx), value = TRUE, ignore.case = TRUE)]]
      }
    }, error = function(x) {
      warning("Matrix is not produced by Cellranger. Trying CellBender read...")
      tryCatch({
        mtx <- scCustomize::Read_CellBender_h5_Mat(file_path)
      }, error = function(x2) {
        warning("Matrix is not produced by CellBender. Trying simple read...")
        mtx <- read.table(file_path,
                          header = TRUE,
                          sep = ",",
                          dec = ".",
                          row.names = 1)
      })
    })

    #Metrics file
    tdy <- format(Sys.Date(), "%y%m%d")
    met_file <- paste0(tdy, "_metrics_", id, ".xlsx")
    df_met <- data.frame(QC_Step = "Input_Shape",
                         nCells = ncol(mtx),
                         nFeatures = nrow(mtx))

    #Create object + Filtering by minimum Genes per cell
    .logger("Create Single Cell Object")
    sce_object <- SeuratObject::CreateSeuratObject(counts = mtx,
                                  project = id,
                                  min.cells=minCellGenes)

    df_met[2,] <- c("Rm_undetected_genes", ncol(sce_object), nrow(sce_object))

    sce_object$sample <- id #some naming
    .logger(paste0("Setting condition in object to: ", sub("[-|_].*", "", id))) # automatized condition settings
    sce_object$condition <- sub("[-|_].*", "", id)

    #Set the filter on mito/ribo genes
    if (include_rbs == TRUE) {
      .logger("Get Mitochondrial+Ribosomal content")
    } else {
      .logger("Get Mitochondrial content")
    }

    if (include_rbs==TRUE) {
      sel_ribofeatures <- grep("^(RPS|RPL)", rownames(sce_object), value = TRUE, ignore.case = TRUE)
      pt_ribo <- Matrix::colSums(GetAssayData(sce_object, layer = 'counts')[sel_ribofeatures, ]) / Matrix::colSums(GetAssayData(sce_object, layer = 'counts'))
      sce_object$pt_ribo <- pt_ribo
    }
    flt_mitofeatures <- grep("^MT-", rownames(sce_object), value = TRUE, ignore.case = TRUE)
    if (length(flt_mitofeatures)<1) {
      warning("Warning: Could not find MT genes")
    }
    pt_mito <- Matrix::colSums(GetAssayData(sce_object, layer = 'counts')[flt_mitofeatures, ]) / Matrix::colSums(GetAssayData(sce_object, layer = 'counts'))
    sce_object$pt_mito <- pt_mito

    .logger("Create QC images")

    #write QC to file
    prefilter_plot <- .QC_Vlnplot(sce_object = sce_object, id, layer = "counts")
    ggsave(plot = prefilter_plot, filename = paste0(outPath, "/QC_Plots_prefiltered.svg"), width = 10, height = 6)

    if (FilterCells==TRUE) {
      .logger("Start Filtering")

      #check if absolute values are set for counts and quantile is set too
      if(!is.null(min_counts) && !is.null(low_quantile)){
        stop("Both filtering minimum counts by absolute values and quantile is set to Value!")
      }

      if(!is.null(max_counts) && !is.null(high_quantile)){
        stop("Both filtering maximum counts by absolute values and quantile is set to Value!")
      }

      if(is.null(min_counts) && is.null(low_quantile)){
        stop("Both minimum count filtering methods are not set!")
      }

      if(is.null(max_counts) && is.null(high_quantile)){
        stop("Both maximum count filtering methods are not set!")
      }

      #Remove all cells with below gene threshold
      if (!is.null(min_genes)) {
        sce_object <- subset(sce_object, subset = nFeature_RNA > min_genes)
      }

      if (!is.null(max_genes)) {
        sce_object <- subset(sce_object, subset = nFeature_RNA < max_genes)
      }
      if (!is.null(min_genes) || !is.null(max_genes)) {
        df_met[3,] <- c("Rm_cells_based_on_genes", ncol(sce_object), nrow(sce_object))
      }

      #Remove all cells with high MT content
      sce_object <- subset(sce_object, subset = pt_mito < cut_mt)
      df_met[4,] <- c("Rm_cell_high_MT", ncol(sce_object), nrow(sce_object))

      #Remove by absolute values or quantiles
      if (!is.null(min_counts) && is.null(low_quantile)) {
        sce_object <- subset(sce_object, subset = nCount_RNA > min_counts)
      } else{
        #Calculate the Quantile for the lower end
        Quality <- data.frame(UMI=sce_object$nCount_RNA,
                              nGene=sce_object$nFeature_RNA,
                              label = factor(sce_object$sample),
                              pt_mit_rib=sce_object$pt_mito)

        Quantile.low.UMI <- Quality %>% group_by(label) %>%
          summarise(UMI = list(tibble::enframe(quantile(UMI,probs = low_quantile)))) %>%
          tidyr::unnest(cols = c(UMI))

        #Subset
        sce_object <- subset(sce_object, subset = nCount_RNA > Quantile.low.UMI$value)
      }

      #Remove by absolute values
      if (!is.null(max_counts) && is.null(high_quantile)) {
        sce_object <- subset(sce_object, subset = nCount_RNA < max_counts)
      } else{
        #Calculate the Quantile for the lower end
        Quality <- data.frame(UMI=sce_object$nCount_RNA,
                              nGene=sce_object$nFeature_RNA,
                              label = factor(sce_object$sample),
                              pt_mit_rib=sce_object$pt_mito)

        Quantile.high.UMI <- Quality %>% group_by(label) %>%
          summarise(UMI = list(tibble::enframe(quantile(UMI,probs = high_quantile)))) %>%
          tidyr::unnest(cols = c(UMI))

        #Subset
        sce_object <- subset(sce_object, subset = nCount_RNA < Quantile.high.UMI$value)
      }
    }


    #save metric files
    write.xlsx(df_met, file = paste0(outPath,"/", met_file))

    #write QC after filtering to file
    postfilter_plot <- .QC_Vlnplot(sce_object = sce_object, id, layer = "counts")
    ggsave(plot = postfilter_plot, filename = paste0(outPath, "/QC_Plots_postfiltered.svg"), width = 10, height = 6)

    #Preprocess steps Seurat
    .logger("Running Normalisation")
    sce_object <- NormalizeData(object = sce_object,verbose = FALSE)

    .logger("Running Variable Gene Detection")
    sce_object <- FindVariableFeatures(object = sce_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

    if(DeleteDoublets==TRUE){
      .logger("Running scDblFinder")
      SCE_obj <- suppressWarnings(as.SingleCellExperiment(sce_object)) # suppress this scale.data empty warning
      SCE_obj <- scDblFinder::scDblFinder(SCE_obj)
      sce_object$scDblFinder_score <- SCE_obj$scDblFinder.score
      sce_object$scDblFinder_class <- SCE_obj$scDblFinder.class
      sce_object <- subset(sce_object, scDblFinder_class == "singlet")
    }

    object_list[[i]] <- sce_object # capture the final objects
  }

  #concatenate objects
  .logger("Merging objects")
  merged_obj <- Reduce(function(x, y) merge(x, y), object_list) # get rid of the error prone approach
  .logger("Running ScaleData")
  merged_obj <- ScaleData(object = merged_obj)
  .logger("Run PCA")
  merged_obj <- RunPCA(merged_obj, verbose = FALSE, ...)
  #No idea why this is needed, but without the next two lines it doesnt work...
  merged_obj <- JoinLayers(merged_obj)

  if (Seurat == FALSE) {
    merged_obj <- as.SingleCellExperiment(merged_obj)
  } else{
    merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$orig.ident)
  }

  return(merged_obj)
}



#' @author Mariano Ruz Jurado
#' @title DO Celltypist
#' @description Runs the CellTypist model on a Seurat or SCE object to predict cell type labels, storing the results as metadata.
#' If the number of cells is less than the specified threshold, it returns NAs for the labels.
#' Optionally updates the CellTypist models and returns the probability matrix.
#' Useful for annotating cell types in single-cell RNA sequencing datasets.
#' @param sce_object The seurat or sce object
#' @param modelName Specify the model you want to use for celltypist
#' @param minCellsToRun If the input seurat or SCE object has fewer than this many cells, NAs will be added for all expected columns and celltypist will not be run.
#' @param runCelltypistUpdate If true, --update-models will be run for celltypist prior to scoring cells.
#' @param over_clustering Column in metadata in object with clustering assignments for cells, default seurat_clusters
#' @param assay_normalized Assay with log1p normalized expressions
#' @param returnProb will additionally return the probability matrix, return will give a list with the first element beeing the object and second prob matrix
#' @param SeuV5 Specify if the Seurat object is made with Seuratv5
#'
#' @importFrom basilisk basiliskRun
#' @import dplyr
#' @import ggplot2
#'
#' @return a seurat or sce object
#'
#' @examples
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#'
#' sce_data <- DO.CellTypist(
#'   sce_object = sce_data,
#'   modelName = "Healthy_Adult_Heart.pkl",
#'   runCelltypistUpdate = TRUE,
#'   over_clustering = "seurat_clusters",
#'   minCellsToRun=5,
#'   SeuV5=TRUE
#' )
#'
#' @export
DO.CellTypist <- function(sce_object,
                          modelName = "Healthy_Adult_Heart.pkl",
                          minCellsToRun = 200,
                          runCelltypistUpdate = TRUE,
                          over_clustering = "seurat_clusters",
                          assay_normalized = "RNA",
                          returnProb=FALSE,
                          SeuV5=TRUE) {

  # Make sure R Zellkonverter package is installed
  zk <- system.file(package = "zellkonverter")
  ifelse(nzchar(zk), "", stop("Install zellkonverter R package for object tranformation!"))

  # Make sure R reticulate package is installed
  rt <- system.file(package = "reticulate")
  ifelse(nzchar(rt), "", stop("Install reticulate R package for Python usage in R!"))

  # if (!reticulate::py_available(initialize = TRUE)) {
  #   stop(paste0('Python/reticulate not correctly configured. Run "usethis::edit_r_environ()" to specify your Python instance'))
  # }
  #
  # if (!reticulate::py_module_available('celltypist')) {
  #   stop('The celltypist python package has not been installed in this python environment!')
  # }

  if (ncol(sce_object) < minCellsToRun) {
    warning('Too few cells, will not run celltypist. NAs will be added instead')
    expectedCols <- c('predicted_labels_celltypist')
    sce_object[[expectedCols]] <- NA
    return(sce_object)
  }

  .logger(paste0('Running celltypist using model: ', modelName))

  #Create temporary folder
  outDir <- tempfile(fileext = '')
  if (endsWith(outDir, "/")){
    outDir <- gsub(outDir, pattern = "/$", replacement = "")
  }
  dir.create(outDir)
  .logger(paste0("Saving celltypist results to temporary folder: ", outDir))

  #Uppercasing gene names
  #zellkonverter h5ad

  #support for Seurat objects
  if (is(sce_object, "Seurat")) {
    DefaultAssay(sce_object) <- assay_normalized
    if (SeuV5 == TRUE) {
      tmp.assay <- sce_object
      tmp.assay[["RNA"]] <- as(tmp.assay[["RNA"]], Class = "Assay")
      tmp.sce <- Seurat::as.SingleCellExperiment(tmp.assay, assay = assay_normalized)
      rownames(tmp.sce) <- toupper(rownames(tmp.sce))

    } else{
      tmp.sce <- Seurat::as.SingleCellExperiment(sce_object, assay = assay_normalized)
      rownames(tmp.sce) <- toupper(rownames(tmp.sce))
    }
  } else{
    tmp.sce <- sce_object
    rownames(tmp.sce) <- toupper(rownames(tmp.sce))
  }

  #Make Anndata object
  if (!"logcounts" %in% names(tmp.sce@assays)) {
    stop("logcounts not found in assays of object!")
  }

  zellkonverter::writeH5AD(tmp.sce, file = paste0(outDir,"/ad.h5ad"), X_name = "logcounts")

  #basilisk implementation
  results <- basilisk::basiliskRun(env = DOtoolsEnv, fun=function(args){

    # Ensure models present:
    if (runCelltypistUpdate) {
      system2(reticulate::py_exe(), c("-m", "celltypist.command_line", "--update-models", "--quiet"))
    }

    .logger("Running Celltypist")

    # Run celltypist:
    args_cty <- c("-m", "celltypist.command_line", "--indata",  paste0(outDir,"/ad.h5ad"), "--model", modelName, "--outdir", outDir,"--majority-voting", "--over-clustering", over_clustering)
    system2(reticulate::py_exe(), args_cty)

  }, args=args)

  labelFile <- paste0(outDir, "/predicted_labels.csv")
  probFile <- paste0(outDir, "/probability_matrix.csv")
  labels <- utils::read.csv(labelFile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)

  #ad <- import(module = "anndata")
  #ad_obj <- ad$AnnData(X = labels)

  #ct <- import(module = "celltypist")
  #ct$dotplot(ad_obj, use_as_reference = "cell_type", use_as_prediction = "majority_voting")

  probMatrix <- utils::read.csv(probFile, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  sce_object$autoAnnot <- labels$majority_voting

  #Create dotplot with prob
  probMatrix$cluster <- as.character(labels$over_clustering)

  #calculate means
  probMatrix_mean <- probMatrix %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(across(where(is.numeric), mean))

  top_cluster <- probMatrix_mean %>%
    rowwise() %>%
    dplyr::mutate(
      prob=max(c_across(-cluster)),
      label= names(.) [which.max(c_across(-cluster))]
      ) %>%
    select(cluster, label, prob)
  top_cluster$pct.exp <- 100 #since majority voting is set TRUE #TODO make this general if majority voting will become a boolean argument in the future

  top_cluster$label <- factor(top_cluster$label, levels = unique(sort(top_cluster$label, decreasing = TRUE)))
  top_cluster$cluster <- factor(top_cluster$cluster, levels = top_cluster$cluster[order(top_cluster$label, decreasing = TRUE)])

  .logger("Creating probality plot")

  pmain <- ggplot(top_cluster, aes(x = cluster, y = label)) +
    geom_point(aes(size = pct.exp, fill = prob), shape=21, color="black", stroke=0.3) +
    scale_fill_gradient2(low = "royalblue3",mid = "white", high = "firebrick",midpoint=0.5, limits = c(0, 1)) +
    scale_size(range = c(2,10), breaks = c(20,40,60,80,100), limits = c(0,100)) +
    theme_box() +
    theme(plot.margin = ggplot2::margin(t = 1,
                                        r = 1,
                                        b = 1,
                                        l = 1,
                                        unit = "cm"),
          axis.text = ggplot2::element_text(color = "black"),
          legend.direction = "horizontal",
          axis.text.x = element_text(color = "black",angle = 90,hjust = 1,vjust = 0.5, size = 14, family = "Helvetica"),
          axis.text.y = element_text(color = "black", size = 14, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title = element_text(size = 14, color = "black", family = "Helvetica"),
          plot.title = element_text(size = 14, hjust = 0.5,face="bold", family = "Helvetica"),
          plot.subtitle = element_text(size = 14, hjust = 0, family = "Helvetica"),
          axis.line = element_line(color = "black"),
          strip.text.x = element_text(size = 14, color = "black", family = "Helvetica", face = "bold"),
          legend.text = element_text(size = 10, color = "black", family = "Helvetica"),
          legend.title = element_text(size = 10, color = "black", family = "Helvetica", hjust =0, face = "bold"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
    )

  guides.layer <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = "Mean propability",
                                                                 title.position = "top",
                                                                 title.hjust = 0.5,
                                                                 barwidth = unit(3.8,"cm"), # changes the width of the color legend
                                                                 barheight = unit(0.5,"cm"),
                                                                 frame.colour = "black",
                                                                 frame.linewidth = 0.3,
                                                                 ticks.colour = "black",
                                                                 order = 2),
                                  size = ggplot2::guide_legend(title = "Fraction of cells \n in group (%)",
                                                               title.position = "top", title.hjust = 0.5, label.position = "bottom",
                                                               override.aes = list(color = "grey40", fill = "grey50"),
                                                               keywidth = ggplot2::unit(0.5, "cm"), # changes the width of the precentage dots in legend
                                                               order = 1))

  pmain <- pmain + guides.layer
  print(pmain)

  if (returnProb==TRUE) {
    returnProb <- list(sce_object, probMatrix)
    names(returnProb) <- c("SingleCellObject", "probMatrix")
    return(returnProb)
  } else {
    return(sce_object)}
}


# Recluster function using FindSubCluster function from Seurat iterative over each cluster -> perfect for fine tuning annotation
#' @author Mariano Ruz Jurado
#' @title DO.FullRecluster
#' @description Performs iterative reclustering on each major cluster found by FindClusters in a Seurat or SCE object.
#' It refines the clusters using the FindSubCluster function for better resolution and fine-tuned annotation.
#' The new clustering results are stored in a metadata column called `annotation_recluster`.
#' Suitable for improving cluster precision and granularity after initial clustering.
#' @param sce_object The seurat or SCE object
#' @param over_clustering Column in metadata in object with clustering assignments for cells, default seurat_clusters
#' @param res Resolution for the new clusters, default 0.5
#' @param algorithm Set one of the available algorithms found in FindSubCLuster function, default = 4: leiden
#' @param graph.name A builded neirest neighbor graph
#' @return a Seurat or SCE Object with new clustering named annotation_recluster
#'
#' @import Seurat
#' @import progress
#'
#' @examples
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' sce_data <- DO.FullRecluster(
#'   sce_object = sce_data
#' )
#'
#' @export
DO.FullRecluster <- function(sce_object,
                             over_clustering = "seurat_clusters",
                             res = 0.5,
                             algorithm=4,
                             graph.name="RNA_snn"){

  #support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    class_obj <- "SingleCellExperiment"
    sce_object <- as.Seurat(sce_object)
    sce_object <- FindNeighbors(sce_object, reduction = "PCA")
  } else{
    class_obj <- "Seurat"
  }

  #TODO add an argument for a new shared nearest neighbor call so you can use this function with other integration methods in the object otherwise it will use the latest calculated shared nearest neighbors and you need to do it outside of the function
  if (is.null(sce_object@meta.data[[over_clustering]])) {
    stop("No clusters defined, please run e.g. FindClusters before Reclustering, or fill the slot with a clustering")
  }
  Idents(sce_object) <- over_clustering

  sce_object$annotation_recluster <- as.vector(sce_object@meta.data[[over_clustering]])

  pb <- progress::progress_bar$new(
    total = length(unique(sce_object@meta.data[[over_clustering]])),
    format = "  Reclustering [:bar] :percent eta: :eta"
  )

  #newline prints when the function exits (to clean up console)
  on.exit(cat("\n"))

  for (cluster in unique(sce_object@meta.data[[over_clustering]])) {
    pb$tick()
    sce_object <- FindSubCluster(sce_object,
                                 cluster = as.character(cluster),
                                 graph.name = graph.name,
                                 algorithm = algorithm,
                                 resolution = res)

    cluster_cells <- rownames(sce_object@meta.data)[sce_object@meta.data[[over_clustering]] == cluster]
    sce_object$annotation_recluster[cluster_cells] <- sce_object$sub.cluster[cluster_cells]
  }
  sce_object$sub.cluster <- NULL

  if (class_obj == "SingleCellExperiment") {
    sce_object <- as.SingleCellExperiment(sce_object)
  }

  return(sce_object)
}

# Polished UMAP function using Dimplot or FeaturePlot function from Seurat
#' @author Mariano Ruz Jurado
#' @title DO.UMAP
#' @description Creates a polished UMAP plot using Seurat's DimPlot or FeaturePlot functions.
#' It allows customization of colors, labels, and other plot elements for better visualisation.
#' The function handles both cluster-based visualisations and gene-based visualisations in a UMAP plot.
#' Ideal for refining UMAP outputs with added flexibility and enhanced presentation.
#' @param sce_object The seurat or SCE object
#' @param FeaturePlot Is it going to be a Dimplot or a FeaturePlot?
#' @param features features for Featureplot
#' @param group.by grouping of plot in DImplot and defines in featureplot the labels
#' @param umap_colors what colors to use for UMAP, specify as vector
#' @param text_size Size of text
#' @param label label the clusters on the plot by group.by column
#' @param order Boolean determining whether to plot cells in order of expression.
#' @param plot.title title for UMAP
#' @param legend.position specify legend position
#' @param ... Further arguments passed to DimPlot or FeaturePlot function from Seurat
#' @return Plot with Refined colors and axes
#'
#' @import Seurat
#' @import ggplot2
#'
#' @examples
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' DO.UMAP(
#'   sce_object = sce_data,
#'   group.by="seurat_clusters"
#' )
#'
#' DO.UMAP(
#'   sce_object = sce_data,
#'   FeaturePlot=TRUE,
#'   features=c("BAG2","CD74")
#' )
#'
#' @export
DO.UMAP <- function(sce_object,
                    FeaturePlot=FALSE,
                    features=NULL,
                    group.by="seurat_clusters",
                    umap_colors=NULL,
                    text_size=14,
                    label=TRUE,
                    order=TRUE,
                    plot.title=TRUE,
                    legend.position="none",
                    ...){

  #support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    sce_object <- as.Seurat(sce_object)
  }

  #Dimplot
  if (FeaturePlot==FALSE) {
    if (is.null(umap_colors)) {
      umap_colors <- rep(c(
        "#1f77b4", "#ff7f0e", "#2ca02c", "tomato2", "#9467bd", "chocolate3","#e377c2", "#ffbb78", "#bcbd22",
        "#17becf","darkgoldenrod2", "#aec7e8", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94","#f7b6d2", "#c7c7c7", "#dbdb8d",
        "#9edae5","sandybrown","moccasin","lightsteelblue","darkorchid","salmon2","forestgreen","bisque"
      ),5)
    }

    p <- DimPlot(sce_object, group.by = group.by, cols = umap_colors, ...) +
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

    if (label==TRUE) {
      p <-LabelClusters(p, id  = group.by, fontface="bold", box = FALSE)
    }
    return(p)
  }

  #FeaturePlot
  if (FeaturePlot==TRUE) {

    if (is.null(features)) {
      stop("Please provide any gene names if using FeaturePlot=TRUE.")
    }

    if (is.null(umap_colors)) {
      umap_colors <- c("lightgrey","red2")
    }

    Idents(sce_object) <- group.by
    p <- FeaturePlot(sce_object,
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

    if (plot.title == FALSE) {
      p <- p & theme(plot.title = element_blank())
    }

    return(p)
  }
}


# SCE Subset function that actually works
#' @author Mariano Ruz Jurado
#' @title DO.Subset
#' @description Creates a subset of a Seurat or SCE object based on either categorical or numeric thresholds in metadata.
#' Allows for subsetting by specifying the ident column, group name, or threshold criteria.
#' Ideal for extracting specific cell populations or clusters based on custom conditions.
#' Returns a new Seurat or SCE object containing only the subsetted cells and does not come with the Seuratv5 subset issue
#' Please be aware that right now, after using this function the subset might be treated with Seuv5=False in other functions.
#' @param sce_object The seurat or SCE object
#' @param assay assay to subset by
#' @param ident meta data column to subset for
#' @param ident_name name of group of barcodes in ident of subset for
#' @param ident_thresh numeric thresholds as character, e.g ">5" or c(">5", "<200"), to subset barcodes in ident for
#' @return a subsetted Seurat or SCE object
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom SingleCellExperiment reducedDim
#' @examples
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' sce_data_sub <- DO.Subset(
#'   sce_object = sce_data,
#'   ident="condition",
#'   ident_name="healthy"
#' )
#'
#'
#'
#' @export
DO.Subset <- function(sce_object,
                      assay="RNA",
                      ident,
                      ident_name=NULL,
                      ident_thresh=NULL){

  #support for Seurat objects
  if (is(sce_object, "Seurat")) {
    Seu_obj <- sce_object
    class_obj <- "Seurat"
    reduction_names <- names(sce_object@reductions)
    sce_object <- as.SingleCellExperiment(sce_object)
  } else{
    class_obj <- "SingleCellExperiment"
  }

  if (!is.null(ident_name) && !is.null(ident_thresh))  {
    stop("Please provide ident_name for subsetting by a name in the column or ident_thresh if it by a threshold")
  }

  #By a name in the provided column
  if (!is.null(ident_name) && is.null(ident_thresh))  {
    .logger("Specified 'ident_name': expecting a categorical variable.")
    sce_object_sub <- sce_object[, SingleCellExperiment::colData(sce_object)[, ident] %in% ident_name]
  }

  #By a threshold in the provided column
  if (is.null(ident_name) && !is.null(ident_thresh))  {
    .logger(paste0("Specified 'ident_thresh': expecting numeric thresholds specified as character, ident_thresh = ", paste0(ident_thresh, collapse = " ")))

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
        sce_object_sub <- sce_object[, SingleCellExperiment::colData(sce_object)[, ident] < threshold]
      } else if (operator == ">") {
        sce_object_sub <- sce_object[, SingleCellExperiment::colData(sce_object)[, ident] > threshold]
      } else if (operator == "<=") {
        sce_object_sub <- sce_object[, SingleCellExperiment::colData(sce_object)[, ident] <= threshold]
      } else if (operator == ">=") {
        sce_object_sub <- sce_object[, SingleCellExperiment::colData(sce_object)[, ident] >= threshold]
      }
    }

    #second case
    if (length(operator) == 2) {
      if (paste(operator,collapse = "") ==  "><") {
        # filtered_cells <- ident_values[ident_values > threshold[1] & ident_values < threshold[2]]
        sce_object_sub <- sce_object[, SingleCellExperiment::colData(sce_object)[, ident] > threshold[1] &
                                       SingleCellExperiment::colData(sce_object)[, ident] < threshold[2]]
      } else if (paste(operator,collapse = "") ==  "<>") {
        sce_object_sub <- sce_object[, SingleCellExperiment::colData(sce_object)[, ident] < threshold[1] &
                                       SingleCellExperiment::colData(sce_object)[, ident] > threshold[2]]
      }
    }
  }

  if (ncol(sce_object_sub) == 0) {
    stop("No cells left after subsetting!\n")
  }

  #Seurat support
  if (class_obj == "Seurat") {
    sce_object_sub <- as.Seurat(sce_object_sub)
    sce_object_sub[[assay]] <- as(object = sce_object_sub[[assay]], Class = "Assay5")
    # Identify reductions that are not in uppercase -> These are the old ones, not subsetted REMOVE
    # non_uppercase_reductions <- reduction_names[!grepl("^[A-Z0-9_]+$", reduction_names)]
    # sce_object_sub@reductions <- sce_object_sub@reductions[!names(sce_object_sub@reductions) %in% non_uppercase_reductions]
    names(sce_object_sub@reductions) <- reduction_names
    sce_object_sub$ident <- NULL
    #some checks
    ncells_interest_prior <- nrow(Seu_obj@meta.data[Seu_obj@meta.data[[ident]] %in% ident_name, ])
    ncells_interest_after <- nrow(sce_object_sub@meta.data[sce_object_sub@meta.data[[ident]] %in% ident_name, ])
    if (ncells_interest_prior != ncells_interest_after) {
      stop(paste0("Number of subsetted cell types is not equal in both objects! Before: ",ncells_interest_prior,"; After: ", ncells_interest_after,". Please check your metadata!"))
    }
  }

  return(sce_object_sub)
}



#' @title Remove Layers from Seurat or SCE Object by Pattern
#' @author Mariano Ruz Jurado
#' @description This function removes layers from a Seurat or SCE object's RNA assay based on a specified regular expression pattern.
#' It is supposed to remove no longer needed layers from th object.
#' @param sce_object Seurat or SCE object.
#' @param assay Name of the assay from where to remove layers from
#' @param pattern regular expression pattern to match layer names. Default "^scale\\.data\\."
#' @return Seurat or SCE object with specified layers removed.
#'
#' @import SeuratObject
#'
#' @examples
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' sce_data <- DO.DietSCE(sce_data, pattern = "data")
#'
#'
#' @export
DO.DietSCE <- function(sce_object, assay = "RNA", pattern = "^scale\\.data\\.") {
  .logger(paste("pattern: ", pattern))

  #support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    class_obj <- "SingleCellExperiment"
    sce_object <- as.Seurat(sce_object)
  } else{
    class_obj <- "Seurat"
  }

  layers_to_remove <- grep(pattern, Layers(sce_object), value = TRUE)

  if ("layers" %in% slotNames(sce_object@assays[[assay]])) {
    sce_object@assays[[assay]]@layers[layers_to_remove] <- NULL
    layerNames <- Layers(sce_object)
    .logger(paste(layers_to_remove, "is removed."))

  } else {
    .logger("Object has no layers, pattern does not need to be removed from layers.")
  }

  if(class_obj == "SingleCellExperiment"){
    sce_object <- as.SingleCellExperiment(sce_object)
  }

  return(sce_object)
}


#' @title DO.PyEnv
#' @description Sets up or connects to a conda Python environment for use with DOtools.
#' If no environment path is provided, it will create one at `~/.venv/DOtools` and install required Python packages:
#' `scvi-tools`, `celltypist`, and `scanpro`.
#'
#' @param conda_path character string specifying the path to an existing or new conda environment.
#' @return None
#'
#' @examples
#' # Automatically create DOtools environment at ~/.venv/DOtools if it doesn't exist
#' DO.PyEnv()
#'
#' # Use an existing conda environment at a custom location
#' DO.PyEnv(conda_path = "~/miniconda3/envs/my_dotools_env")
#'
#'
#' @export
DO.PyEnv <- function(conda_path = NULL) {

  # Handle conda environment
  if (is.null(conda_path)) {
    conda_path <- file.path(path.expand("~"), ".venv/DOtools")
  } else {
    conda_path <- path.expand(conda_path)
  }

  if (!dir.exists(conda_path)) {
    .logger("Creating conda environment for DOtools")
    conda_args1 <- c("conda","create","-y", "-p", file.path(path.expand("~"), ".venv/DOtools"), "python=3.11")
    conda_args2 <- c("conda","run", "-p", file.path(path.expand("~"), ".venv/DOtools"), "pip", "install", "scvi-tools==1.3.0", "celltypist==1.6.3", "scanpro==0.3.2")
    conda_args3 <- c("conda","run", "-p", file.path(path.expand("~"), ".venv/DOtools"), "pip", "install", "scipy==1.15.3")
    tryCatch({
      system2(conda_args1[1], args = conda_args1[-1], stdout = FALSE, stderr = FALSE)
      system2(conda_args2[1], args = conda_args2[-1], stdout = FALSE, stderr = FALSE)
      system2(conda_args3[1], args = conda_args3[-1], stdout = FALSE, stderr = FALSE)
      }, error = function(e) {
        stop("Failed to create conda environment. Provide a valid environment path or fix installation.")
        })
    } else {
      .logger(sprintf("Using existing conda environment at: %s", conda_path))
      }
  .logger("Python packages ready for DOtools!")
}

#' @title DO.CellBender
#' @description It is supposed to make something similar than the no longer working DietSeurat function, by removing no longer needed layers from th object.
#' This function wraps a system call to a bash script for running CellBender on CellRanger outputs.
#' It ensures required inputs are available and optionally installs CellBender in a conda env.
#'
#' @param cellranger_path Path to folder with CellRanger outputs.
#' @param output_path Output directory for CellBender results.
#' @param samplenames Optional vector of sample names. If NULL, will autodetect folders in `cellranger_path`.
#' @param cuda Logical, whether to use GPU (CUDA).
#' @param cpu_threads Number of CPU threads to use.
#' @param epochs Number of training epochs.
#' @param lr Learning rate.
#' @param estimator_multiple_cpu Use estimator with multiple CPU threads (experimental).
#' @param log Whether to enable logging.
#' @param conda_path Optional path to the conda environment.
#' @param BarcodeRanking Optional Calculation of estimated cells in samples through DropletUtils implementation
#' @param bash_script Path to the bash script that runs CellBender.
#'
#'
#' @import DropletUtils
#' @return None
#'
#' @examples
#' \dontrun{
#' # Define paths
#' cellranger_path <- "/mnt/data/cellranger_outputs"
#' output_path <- "/mnt/data/cellbender_outputs"
#'
#' # Optional: specify sample names if automatic detection is not desired
#' samplenames <- c("Sample_1", "Sample_2")
#'
#' # Run CellBender (uses GPU by default)
#' DO.CellBender(cellranger_path = cellranger_path,
#'               output_path = output_path,
#'               samplenames = samplenames,
#'               cuda = TRUE,
#'               cpu_threads = 8,
#'               epochs = 100,
#'               lr = 0.00001,
#'               estimator_multiple_cpu = FALSE,
#'               log = TRUE)
#' }
#'
#' @export
DO.CellBender <- function(cellranger_path,
                          output_path,
                          samplenames = NULL,
                          cuda = TRUE,
                          cpu_threads = 15,
                          epochs = 150,
                          lr = 0.00001,
                          estimator_multiple_cpu = FALSE,
                          log = TRUE,
                          conda_path = NULL,
                          BarcodeRanking = TRUE,
                          bash_script = system.file("bash", "_run_CellBender.sh", package = "DOtools")) {

  # Check input paths
  stopifnot(file.exists(cellranger_path))
  stopifnot(file.exists(output_path))

  # Warnings and logs
  if (estimator_multiple_cpu)
    .logger("Warning: estimator_multiple_cpu is TRUE. Not recommended for large datasets (above 20 to 30k cells).")
  if (epochs > 150)
    .logger("Warning: Training for more than 150 epochs may lead to overfitting.")
  if (!cuda)
    .logger("Warning: Running without CUDA (GPU) may significantly increase run time.")

  # Handle conda environment
  if (is.null(conda_path)) {
    conda_path <- file.path(path.expand("~"), ".cache/R/basilisk/1.20.0/DOtools/0.99.0/CellBen_env")
  } else {
    conda_path <- path.expand(conda_path)
  }

  if (!dir.exists(conda_path)) {
    .logger("Creating conda environment for CellBender...")
    conda_args1 <- c("conda","create","-y", "-p", file.path(path.expand("~"), ".cache/R/basilisk/1.20.0/DOtools/0.99.0/CellBen_env"), "python=3.7")
    conda_args2 <- c("conda","run", "-p", file.path(path.expand("~"), ".cache/R/basilisk/1.20.0/DOtools/0.99.0/CellBen_env"), "pip", "install", "cellbender", "lxml_html_clean")
    tryCatch({
      system2(conda_args1[1], args = conda_args1[-1], stdout = FALSE, stderr = FALSE)
      system2(conda_args2[1], args = conda_args2[-1], stdout = FALSE, stderr = FALSE)
    }, error = function(e) {
      stop("Failed to create conda environment. Provide a valid environment path or fix installation.")
    })
  } else {
    .logger(sprintf("Using existing conda environment at: %s", conda_path))
  }

  # Detect samples
  if (is.null(samplenames)) {
    samples <- list.dirs(cellranger_path, full.names = FALSE, recursive = FALSE)
  } else {
    samples <- samplenames
  }

  .logger(sprintf("Running CellBender for %d sample(s)", length(samples)))

  # Loop through each sample
  for (sample in samples) {
    # Load HDF5 with reticulate or scanpy for expected_cells, total_droplets

    # Run one by one but sequentially
    # Estimate the number of cells to be used as upper and lower bound
    h5_file <- file.path(cellranger_path, sample, "outs", "raw_feature_bc_matrix.h5")
    tdata <- DropletUtils::read10xCounts(h5_file)

    if (BarcodeRanking == TRUE) {
      result_barcoderanks <- .DO.BarcodeRanks(tdata)

      # Build command
      cmd <- c("conda", "run", "-p", conda_path,
               "bash", bash_script,
               "-i", sample, "-o", output_path,
               "--cellRanger-output", cellranger_path,
               "--cpu-threads", cpu_threads,
               "--epochs", epochs,
               "--lr", lr,
               "--expected-cells", result_barcoderanks[1],
               "--total-droplets", result_barcoderanks[2])
    } else{

      # Build command
      cmd <- c("conda", "run", "-p", conda_path,
               "bash", bash_script,
               "-i", sample, "-o", output_path,
               "--cellRanger-output", cellranger_path,
               "--cpu-threads", cpu_threads,
               "--epochs", epochs,
               "--lr", lr)
    }

    if (cuda) cmd <- c(cmd, "--cuda")
    if (log) cmd <- c(cmd, "--log")
    if (estimator_multiple_cpu) cmd <- c(cmd, "--estimator_multiple_cpu")

    # Run command
    .logger(sprintf("Running CellBender for sample: %s", sample))
    tryCatch({
      system2(cmd[1], args = cmd[-1], stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      .logger(sprintf("Error running CellBender for sample %s: %s", sample, e$message))
    })
  }

  .logger("Finished running CellBender.")
  invisible(NULL)
}



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
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' # Run scVI using the 'orig.ident' column as the batch key
#' sce_data <- DO.scVI(sce_data, batch_key = "orig.ident")
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


#' @title DO.TransferLabel
#' @description Transfers cell-type annotations from a re-annotated subset of a Seurat or SCE object
#' back to the full Seurat or SCE object. This is useful when clusters have been refined
#' or re-labeled in a subset and need to be reflected in the original object.
#'
#' @param sce_object Seurat or SCE object with annotation in meta.data
#' @param Subset_obj subsetted Seurat or SCE object with re-annotated clusters
#' @param annotation_column column name in meta.data with annotation
#' @param subset_annotation column name in meta.data with annotation in the subsetted object
#'
#'
#' @import Seurat
#'
#'
#' @examples
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' sce_data <- DO.TransferLabel(sce_data,
#'                             sce_data,
#'                            annotation_column="annotation",
#'                            subset_annotation="annotation"
#'                           )
#'
#'
#'
#' @return Seurat or SCE Object with transfered labels
#' @export

DO.TransferLabel <- function(sce_object,
                             Subset_obj,
                             annotation_column,
                             subset_annotation){

  #support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    class_obj <- "SingleCellExperiment"
    sce_object <- as.Seurat(sce_object)
    Subset_obj <- as.Seurat(Subset_obj)
  } else{
    class_obj <- "Seurat"
  }

  #Get annotation with barcodes as rownames
  annotation_pre <- sce_object@meta.data[annotation_column]
  annotation_pre[[annotation_column]] <- as.character(annotation_pre[[annotation_column]])
  annotation_subset <- Subset_obj@meta.data[subset_annotation]
  annotation_subset[[subset_annotation]] <- as.character(annotation_subset[[subset_annotation]])

  #assign the new labels
  barcodes <- rownames(annotation_pre)[rownames(annotation_pre) %in% rownames(annotation_subset)]
  annotation_pre[rownames(annotation_pre) %in% barcodes,] <- annotation_subset[[subset_annotation]]

  sce_object@meta.data[[annotation_column]] <- factor(annotation_pre$annotation)

  if (class_obj == "SingleCellExperiment") {
    sce_object <- as.SingleCellExperiment(sce_object)
  }

  return(sce_object)

}

#' @title DO.enrichR
#' @description
#' Performs Gene Ontology enrichment analysis on differentially expressed genes using the EnrichR API.
#' Separately analyzes upregulated and downregulated genes and returns results.
#'
#' @param df_DGE data.frame containing differential gene expression results.
#' @param gene_column column name in `df` with gene symbols.
#' @param pval_column column name in `df` with p-values.
#' @param log2fc_column column name in `df` with log2 fold changes.
#' @param pval_cutoff adjusted p-value threshold for significance (default = 0.05).
#' @param log2fc_cutoff log2 fold change threshold for up/down regulation (default = 0.25).
#' @param path folder path where the output Excel file will be saved. A subfolder "GSA_Tables" will be created.
#' @param filename suffix used in the Excel filename (e.g., "GSA_CellType_MyAnalysis.xlsx").
#' @param species species name for enrichment analysis. Options include "Human", "Mouse", "Yeast", etc. (default = "Mouse").
#' @param go_catgs GO databases to use. Defaults to c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023").
#'
#' @return data.frame with GO enrichment results if `path` is NULL, otherwise writes an Excel file.
#'
#'
#' @import enrichR
#' @import openxlsx
#' @import dplyr
#'
#' @examples
#' library(enrichR)
#'
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#' DGE_result <- DO.MultiDGE(sce_data,
#'                          sample_col = "orig.ident",
#'                          method_sc = "wilcox",
#'                          annotation_col = "annotation",
#'                          ident_ctrl = "healthy")
#'
#' DGE_result <- DGE_result[DGE_result$celltype == "CD4_T_cells",]
#'
#' result_GO <- DO.enrichR(df_DGE = DGE_result,
#'                        gene_column = "gene",
#'                        pval_column = "p_val_SC_wilcox",
#'                        log2fc_column = "avg_log2FC_SC_wilcox",
#'                        pval_cutoff = 0.05,
#'                        log2fc_cutoff = 0.25,
#'                        path = NULL,
#'                        filename = "",
#'                        species = "Human",
#'                        go_catgs = "GO_Biological_Process_2023")
#'
#' @export
DO.enrichR <- function(df_DGE,
                       gene_column,
                       pval_column,
                       log2fc_column,
                       pval_cutoff = 0.05,
                       log2fc_cutoff = 0.25,
                       path = NULL,
                       filename = '',
                       species = 'Human',
                       go_catgs = c('GO_Molecular_Function_2023',
                                    'GO_Cellular_Component_2023',
                                    'GO_Biological_Process_2023')) {


  # Subset up/down-regulated
  df_up <- df_DGE[df_DGE[[pval_column]] < pval_cutoff & df_DGE[[log2fc_column]] > log2fc_cutoff, ]
  df_down <- df_DGE[df_DGE[[pval_column]] < pval_cutoff & df_DGE[[log2fc_column]] < -log2fc_cutoff, ]

  # Set enrichR organism
  enrichR::setEnrichrSite("Enrichr")  # explicitly connect to Enrichr API

  # Run enrichment
  res_up <- enrichR::enrichr(genes = df_up[[gene_column]], databases = go_catgs)
  res_down <- enrichR::enrichr(genes = df_down[[gene_column]], databases = go_catgs)

  # Combine results
  enrich_results <- function(result_list, label) {
    do.call(rbind, lapply(names(result_list), function(db) {
      df <- result_list[[db]]
      df$Database <- db
      df$State <- label
      return(df)
    }))
  }

  combined_up <- enrich_results(res_up, "enriched")
  combined_down <- enrich_results(res_down, "depleted")

  combined_res <- rbind(combined_up, combined_down)

  # Write to Excel or return
  if (!is.null(path)) {
    dir.create(file.path(path, "GO_Tables"), showWarnings = FALSE, recursive = TRUE)
    out_file <- file.path(path, "GO_Tables", paste0("GO_", filename, ".xlsx"))
    openxlsx::write.xlsx(combined_res, out_file, rowNames = FALSE)
    .logger(paste0("Write results to: ", out_file))
    return(invisible(NULL))
  } else {
    return(combined_res)
  }
}

#' @author Mariano Ruz Jurado
#' @title DO.MultiDGE
#' @description
#' Performs differential gene expression analysis using both single-cell and pseudo-bulk approaches across all annotated cell types.
#' The single-cell method uses Seurat's `FindMarkers`, while pseudo-bulk testing uses `DESeq2` on aggregated expression profiles.
#' Outputs a merged data frame with DGE statistics from both methods per condition and cell type.
#'
#' @param sce_object The seurat or SCE object
#' @param assay Specified assay in Seurat or SCE object, default "RNA"
#' @param method_sc method to use for single cell DEG analysis, see FindMarkers from Seurat for options, default "wilcox"
#' @param group_by Column in meta data containing groups used for testing, default "condition"
#' @param annotation_col Column in meta data containing information of cell type annotation
#' @param sample_col Column in meta data containing information of sample annotation, default "orig.ident"
#' @param ident_ctrl Name of the condition in group_by to test against as ctrl, default "ctrl"
#' @param min_pct only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations, default is 0
#' @param logfc_threshold Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells, default is 0.
#' @param only_pos Only return positive markers, default FALSE
#' @param min_cells_group Minimum number of cells in one of the groups, default 3
#' @param ... Additional arguments passed to FindMarkers function
#'
#' @import Seurat
#' @import DESeq2
#' @import tibble
#' @import dplyr
#'
#' @return Dataframe containing statistics for each gene from the single cell and the Pseudobulk DGE approach.
#'
#' @examples
#'
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#' DGE_result <- DO.MultiDGE(sce_data,
#'                          sample_col = "orig.ident",
#'                          method_sc = "wilcox",
#'                          annotation_col = "annotation",
#'                          ident_ctrl = "healthy")
#'
#' @export
DO.MultiDGE <- function(sce_object,
                        assay="RNA",
                        method_sc="wilcox",
                        group_by="condition",
                        annotation_col="annotation",
                        sample_col="orig.ident",
                        ident_ctrl="ctrl",
                        min_pct=0,
                        logfc_threshold=0,
                        only_pos=FALSE,
                        min_cells_group=3,
                        ...
                        ){

  #support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    class_obj <- "SingleCellExperiment"
    sce_object <- as.Seurat(sce_object)
  }

  if (!ident_ctrl %in% sce_object@meta.data[[group_by]]) {
    stop(paste0(ident_ctrl, " was not found in meta data under the specified group_by column: ", group_by))
  }

  DEG_stats_collector_sc <- list() #list for collecting sc method results
  DEG_stats_collector_pb <- list() #list for collecting pseudo_bulk method results

  #Create a PseudoBulk representation of the Seurat single cell expression
  sce_object_PB <- AggregateExpression(sce_object,
                                    assays = assay,
                                    return.seurat = TRUE,
                                    group.by = c(group_by, sample_col, annotation_col))
  #internally AggregateExpression uses interaction() this replaces _ with -, check if the names changed in annotation:
  original_annots <- unique(sce_object@meta.data[[annotation_col]])
  pb_annots <- unique(sce_object_PB@meta.data[[annotation_col]])

  #are all pseudo-bulk annotations in the original set
  if (!all(pb_annots %in% original_annots)) {
    #if it is a hyphen exchange this should fix it
    converted_annots <- gsub("-", "_", pb_annots)

    if (all(converted_annots %in% original_annots)) {
      sce_object_PB@meta.data[[annotation_col]] <- gsub("-", "_", sce_object_PB@meta.data[[annotation_col]])
      .logger("Corrected annotation names in pseudo-bulk object by replacing '-' with '_'.")
    } else {
      .logger("Annotation names do not match even after replacing '-' with '_'. Please inspect manually.")
    }
  } else {
    .logger("Annotation names are consistent between original and pseudo-bulk objects.")
  }

  PB_ident <- paste(annotation_col,group_by, sep = "_")
  sce_object_PB@meta.data[[PB_ident]] <- paste(sce_object_PB@meta.data[[annotation_col]], sce_object_PB@meta.data[[group_by]], sep = "_")
  Idents(sce_object_PB) <- PB_ident

  .logger("Starting DGE single cell method analysis")
  for (ident_con in unique(sce_object@meta.data[[group_by]])) {
    # .logger(ident_con)
    if(ident_con != ident_ctrl){ # skip the ctrl ident
      for (celltype in unique(sce_object@meta.data[[annotation_col]])) {
        .logger(paste0("Comparing ",ident_con, " with ", ident_ctrl," in: ", celltype))
        Seu_celltype <- subset(sce_object, subset = !!sym(annotation_col) == celltype)
        #Check if there are groups with less than 3 cells
        table_cells_sc <- table(Seu_celltype@meta.data[[group_by]])

        count_1_sc <- table_cells_sc[ident_con]
        count_2_sc <- table_cells_sc[ident_ctrl]

        if (is.na(count_1_sc)) {
          count_1_sc <- 0
        }

        if (is.na(count_2_sc)) {
          count_2_sc <- 0
        }

        if (count_1_sc >= 3 && count_2_sc >= 3) {
          DEG_stats_sc <- FindMarkers(object = Seu_celltype,
                                      ident.1 = ident_con,
                                      ident.2 = ident_ctrl,
                                      min.pct = min_pct,
                                      logfc.threshold = logfc_threshold,
                                      assay = assay,
                                      only.pos = only_pos,
                                      test.use = method_sc,
                                      min.cells.group = min_cells_group,
                                      group.by = group_by,
                                      ...)

          DEG_stats_sc <- rownames_to_column(DEG_stats_sc, var="gene")
          DEG_stats_sc[["condition"]] <- ident_con
          DEG_stats_sc[["celltype"]] <- celltype
          # DEG_stats_sc[["testmethod"]] <- method_sc
          DEG_stats_collector_sc <- rbind(DEG_stats_collector_sc, DEG_stats_sc)
        } else {
          .logger(paste0("Skipping ", celltype, " since one comparison has fewer than 3 cells!"))
        }
      }
    }
  }

  .logger("Finished DGE single cell method analysis")
  .logger("Starting DGE pseudo bulk method analysis")

  for (celltype in unique(sce_object_PB@meta.data[[annotation_col]])) {
    # .logger(celltype)
    # running find markers condition specific
    for (ident_con in unique(sce_object_PB@meta.data[[group_by]])) {
      if (ident_con != ident_ctrl){
        .logger(paste0("Comparing ",ident_con, " with ", ident_ctrl," in: ", celltype))
        ident_1 <- paste0(c(celltype,ident_con), collapse = "_")
        ident_2 <- paste0(c(celltype,ident_ctrl), collapse = "_")

        #Check if any of the groups have fewer than 3 cells, test fail if this is the case
        table_cells_pb <- table(Idents(sce_object_PB))

        count_1_pb <- table_cells_pb[ident_1]
        count_2_pb <- table_cells_pb[ident_2]

        if (is.na(count_1_pb)) {
          count_1_pb <- 0
        }

        if (is.na(count_2_pb)) {
          count_2_pb <- 0
        }


        if (count_1_pb >= 3 && count_2_pb >= 3) {
          DEG_stats_pb <- FindMarkers(object = sce_object_PB,
                                      ident.1 = ident_1,
                                      ident.2 = ident_2,
                                      min.pct = min_pct,
                                      logfc.threshold = logfc_threshold,
                                      assay = assay,
                                      only.pos = only_pos,
                                      test.use = "DESeq2",
                                      min.cells.group = min_cells_group,
                                      group.by = PB_ident,
                                      ...)

          DEG_stats_pb <- rownames_to_column(DEG_stats_pb, var = "gene")
          DEG_stats_pb[["condition"]] <- ident_con
          DEG_stats_pb[["celltype"]] <- celltype
          # DEG_stats_pb[["testmethod"]] <- "DESeq2"
          DEG_stats_collector_pb <- rbind(DEG_stats_collector_pb, DEG_stats_pb)
        } else{
          .logger(paste0("Skipping ", celltype, " since one comparison has fewer than 3 cells!"))
        }
      }
    }
  }
  .logger("Finished DGE pseudo bulk method analysis")

  if (!length(DEG_stats_collector_pb) == 0) {
    #combine the two df
    df_pb <- DEG_stats_collector_pb %>%
      rename(
        p_val_PB_DESeq2 = p_val,
        p_val_adj_PB_DESeq2 = p_val_adj,
        avg_log2FC_PB_DESeq2 = avg_log2FC
      ) %>%
      select(gene, celltype, condition, p_val_PB_DESeq2, p_val_adj_PB_DESeq2, avg_log2FC_PB_DESeq2)
  } else{
    .logger("DGE pseudo bulk result is empty...")
    df_pb <- data.frame(gene=NA, celltype=NA, condition=NA)
  }


  rename_map <- rlang::set_names(
    c("p_val", "p_val_adj", "avg_log2FC"),
    paste0(c("p_val_SC_", "p_val_adj_SC_", "avg_log2FC_SC_"), method_sc)
  )

  df_sc <- DEG_stats_collector_sc %>%
    rename(!!!rename_map) # !!! to unquote the strings

  #Sorting
  first_cols <- c("gene","pct.1","pct.2", "celltype", "condition")

  merged_df <- df_sc %>%
    left_join(df_pb, by = c("gene", "celltype", "condition")) %>%
    select(all_of(first_cols), everything()) %>%
    select(all_of(first_cols), sort(setdiff(names(.), first_cols)))

  return(merged_df)
}


#' @title DO.BarcodeRanks
#' @author Mariano Ruz Jurado
#' @description Given a raw count matrix (e.g. from a CellRanger HDF5 file), estimate the number of expected cells and droplets
#' using the knee and inflection points from barcodeRanks.
#'
#' @param SCE_obj A Single cell experiment object.
#'
#' @return A named numeric vector: `c(xpc_cells = ..., total_cells = ...)
#' @rdname dot-DO.BarcodeRanks
#' @keywords internal
.DO.BarcodeRanks <- function(SCE_obj) {
  if (!inherits(SCE_obj, c("SingleCellExperiment"))) {
    stop("Input must be a sparse matrix of class 'dgCMatrix' or 'DelayedMatrix'.")
  }

  result <- DropletUtils::barcodeRanks(SCE_obj)
  metadata <- S4Vectors::metadata(result)

  knee <- metadata$knee
  inflection <- metadata$inflection

  colsum <- colSums(SCE_obj@assays@data$counts)

  xpc_cells <- length(colsum[colsum > knee])
  total_cells <- length(colsum[colsum > inflection])

  return(c(xpc_cells = xpc_cells, total_cells = total_cells))
}


umap_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "tomato2", "#9467bd", "chocolate3","#e377c2", "#ffbb78", "#bcbd22",
  "#17becf","darkgoldenrod2", "#aec7e8", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94","#f7b6d2", "#c7c7c7", "#dbdb8d",
  "#9edae5","sandybrown","moccasin","lightsteelblue","darkorchid","salmon2","forestgreen","bisque"
)


#' @title .QC_Vlnplot
#' @author Mariano Ruz Jurado
#' @description
#' Generates violin plots for common quality control (QC) metrics of single-cell RNA-seq data from a Seurat object.
#' The function displays three violin plots for the number of detected genes per cell (`nFeature_RNA`), total UMI counts per cell (`nCount_RNA`),
#' and mitochondrial gene content percentage (`pt_mito`). Useful for visual inspection of QC thresholds and outliers.
#'
#' @param sce_object A Seurat object containing single-cell RNA-seq data.
#' @param layer A character string specifying the assay layer to use (default is "counts").
#' @param features A character vector of length 3 indicating the feature names to plot. Default is c("nFeature_RNA", "nCount_RNA", "pt_mito").
#'
#' @return A `ggplot` object arranged in a single row showing violin plots for the specified features with overlaid boxplots.
#'
#' @import Seurat
#' @import ggplot2
#' @rdname dot-QC_Vlnplot
#' @keywords internal
.QC_Vlnplot <- function(sce_object, id, layer="counts", features=c("nFeature_RNA","nCount_RNA","pt_mito")){
  p1<- VlnPlot(sce_object,layer = "counts", features = features[1], ncol = 1, pt.size = 0, cols = "grey")+
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha=0.5, color="darkred", size=0.4)+
    theme_light()+
    xlab("") +
    theme(
      plot.title = element_text(face = "bold", color = "black", hjust = 0.5, size = 14),
      axis.title.y = element_text(face = "bold", color = "black", size = 14),
      axis.text.x = element_blank(),
      axis.text.y = element_text(face = "bold", color = "black", hjust = 1, size = 14),
      axis.title.x = element_blank(),
      legend.position = "none"
    )

  p2<- VlnPlot(sce_object,layer = "counts", features = features[2], ncol = 1, pt.size = 0, cols = "grey")+
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha=0.5, color="darkred", size=0.4)+
    theme_light()+
    xlab("") +
    theme(
      plot.title = element_text(face = "bold", color = "black", hjust = 0.5, size = 14),
      axis.title.y = element_text(face = "bold", color = "black", size = 14),
      axis.text.x = element_blank(),
      axis.text.y = element_text(face = "bold", color = "black", hjust = 1, size = 14),
      axis.title.x = element_blank(),
      legend.position = "none"
    )

  p3<-VlnPlot(sce_object,layer = "counts", features = features[3], ncol = 1, pt.size = 0, cols = "grey")+
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha=0.5, color="darkred", size=0.4)+
    theme_light()+
    xlab("") +
    theme(
      plot.title = element_text(face = "bold", color = "black", hjust = 0.5, size = 14),
      axis.title.y = element_text(face = "bold", color = "black", size = 14),
      axis.text.x = element_blank(),
      axis.text.y = element_text(face = "bold", color = "black", hjust = 1, size = 14),
      axis.title.x = element_blank(),
      legend.position = "none"
    )

  gg_plot <- ggarrange(p1,p2,p3, nrow = 1)
  gg_plot <- annotate_figure(gg_plot,top = grid::grobTree(grid::rectGrob(gp = grid::gpar(fill = "white", col = NA), height = unit(1.5, "lines")),
                                                          grid::textGrob(id, gp = grid::gpar(fontface = "bold", col = "darkred", fontsize = 18),hjust = 0.35, vjust=0.6  # align horizontally at 30%
                                                          )
  )
  )

  return(gg_plot)
}


#' @keywords internal
.logger <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, " - ", message, "\n"))
}
