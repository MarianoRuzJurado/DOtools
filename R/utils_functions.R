#' @title DO.Import
#' @author Mariano Ruz Jurado & David John
#' @description
#' Imports and processes single-cell RNA-seq data from various formats (10x Genomics, CellBender, or CSV),
#' performs quality control (QC), filtering, normalization, variable gene selection, and optionally detects doublets.
#' Returns a merged and processed Seurat object ready for downstream analysis.
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
#' @param low_quantile Numeric. Quantile threshold (0–1) to filter low UMI cells (used if `min_counts` is `NULL`).
#' @param high_quantile Numeric. Quantile threshold (0–1) to filter high UMI cells (used if `max_counts` is `NULL`).
#' @param DeleteDoublets Logical. If `TRUE`, doublets are detected and removed using `scDblFinder`. Default is `TRUE`.
#' @param include_rbs Logical. If `TRUE`, calculates ribosomal gene content in addition to mitochondrial content. Default is `TRUE`.
#' @param ... Additional arguments passed to `RunPCA()`.
#'
#' @return A merged Seurat object containing all samples, with normalization, QC, scaling, PCA, and optional doublet removal applied.
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
                      FilterCells=T,
                      cut_mt=.05,
                      min_counts=NULL,
                      max_counts=NULL,
                      min_genes=NULL,
                      max_genes=NULL,
                      low_quantile=NULL,
                      high_quantile=NULL,
                      DeleteDoublets=T,
                      include_rbs=T,
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
        mtx <- mtx[[grep("gene",names(mtx), value = T, ignore.case = T)]]
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
    .logger("Create Seurat")
    Seu_obj <- SeuratObject::CreateSeuratObject(counts = mtx,
                                  project = id,
                                  min.cells=minCellGenes)

    df_met[2,] <- c("Rm_undetected_genes", ncol(Seu_obj), nrow(Seu_obj))

    Seu_obj$sample <- id #some naming
    .logger(paste0("Setting condition in object to: ", sub("[-|_].*", "", id))) # automatized condition settings
    Seu_obj$condition <- sub("[-|_].*", "", id)

    #Set the filter on mito/ribo genes
    if (include_rbs == TRUE) {
      .logger("Get Mitochondrial+Ribosomal content")
    } else {
      .logger("Get Mitochondrial content")
    }

    if (include_rbs==T) {
      sel_ribofeatures <- grep("^(RPS|RPL)", rownames(Seu_obj), value = TRUE, ignore.case = TRUE)
      pt_ribo <- Matrix::colSums(GetAssayData(Seu_obj, layer = 'counts')[sel_ribofeatures, ]) / Matrix::colSums(GetAssayData(Seu_obj, layer = 'counts'))
      Seu_obj$pt_ribo <- pt_ribo
    }
    flt_mitofeatures <- grep("^MT-", rownames(Seu_obj), value = TRUE, ignore.case = TRUE)
    if (length(flt_mitofeatures)<1) {
      warning("Warning: Could not find MT genes")
    }
    pt_mito <- Matrix::colSums(GetAssayData(Seu_obj, layer = 'counts')[flt_mitofeatures, ]) / Matrix::colSums(GetAssayData(Seu_obj, layer = 'counts'))
    Seu_obj$pt_mito <- pt_mito

    .logger("Create QC images")

    #write QC to file
    prefilter_plot <- .QC_Vlnplot(Seu_obj = Seu_obj, id, layer = "counts")
    ggsave(plot = prefilter_plot, filename = paste0(outPath, "/QC_Plots_prefiltered.svg"), width = 10, height = 6)

    if (FilterCells==T) {
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
        Seu_obj <- subset(Seu_obj, subset = nFeature_RNA > min_genes)
      }

      if (!is.null(max_genes)) {
        Seu_obj <- subset(Seu_obj, subset = nFeature_RNA < max_genes)
      }
      if (!is.null(min_genes) || !is.null(max_genes)) {
        df_met[3,] <- c("Rm_cells_based_on_genes", ncol(Seu_obj), nrow(Seu_obj))
      }

      #Remove all cells with high MT content
      Seu_obj <- subset(Seu_obj, subset = pt_mito < cut_mt)
      df_met[4,] <- c("Rm_cell_high_MT", ncol(Seu_obj), nrow(Seu_obj))

      #Remove by absolute values or quantiles
      if (!is.null(min_counts) && is.null(low_quantile)) {
        Seu_obj <- subset(Seu_obj, subset = nCount_RNA > min_counts)
      } else{
        #Calculate the Quantile for the lower end
        Quality <- data.frame(UMI=Seu_obj$nCount_RNA,
                              nGene=Seu_obj$nFeature_RNA,
                              label = factor(Seu_obj$sample),
                              pt_mit_rib=Seu_obj$pt_mito)

        Quantile.low.UMI <- Quality %>% group_by(label) %>%
          summarise(UMI = list(tibble::enframe(quantile(UMI,probs = low_quantile)))) %>%
          tidyr::unnest(cols = c(UMI))

        #Subset
        Seu_obj <- subset(Seu_obj, subset = nCount_RNA > Quantile.low.UMI$value)
      }

      #Remove by absolute values
      if (!is.null(max_counts) && is.null(high_quantile)) {
        Seu_obj <- subset(Seu_obj, subset = nCount_RNA < max_counts)
      } else{
        #Calculate the Quantile for the lower end
        Quality <- data.frame(UMI=Seu_obj$nCount_RNA,
                              nGene=Seu_obj$nFeature_RNA,
                              label = factor(Seu_obj$sample),
                              pt_mit_rib=Seu_obj$pt_mito)

        Quantile.high.UMI <- Quality %>% group_by(label) %>%
          summarise(UMI = list(tibble::enframe(quantile(UMI,probs = high_quantile)))) %>%
          tidyr::unnest(cols = c(UMI))

        #Subset
        Seu_obj <- subset(Seu_obj, subset = nCount_RNA < Quantile.high.UMI$value)
      }
    }


    #save metric files
    write.xlsx(df_met, file = paste0(outPath,"/", met_file))

    #write QC after filtering to file
    postfilter_plot <- .QC_Vlnplot(Seu_obj = Seu_obj, id, layer = "counts")
    ggsave(plot = postfilter_plot, filename = paste0(outPath, "/QC_Plots_postfiltered.svg"), width = 10, height = 6)

    #Preprocess steps Seurat
    .logger("Running Normalisation")
    Seu_obj <- NormalizeData(object = Seu_obj,verbose = FALSE)

    .logger("Running Variable Gene Detection")
    Seu_obj <- FindVariableFeatures(object = Seu_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

    if(DeleteDoublets==TRUE){
      .logger("Running scDblFinder")
      SCE_obj <- suppressWarnings(as.SingleCellExperiment(Seu_obj)) # suppress this scale.data empty warning
      SCE_obj <- scDblFinder::scDblFinder(SCE_obj)
      Seu_obj$scDblFinder_score <- SCE_obj$scDblFinder.score
      Seu_obj$scDblFinder_class <- SCE_obj$scDblFinder.class
      Seu_obj <- subset(Seu_obj, scDblFinder_class == "singlet")
    }

    object_list[[i]] <- Seu_obj # capture the final objects
  }

  #concatenate objects
  .logger("Merging objects")
  merged_obj <- Reduce(function(x, y) merge(x, y), object_list) # get rid of the error prone approach
  .logger("Running ScaleData")
  merged_obj <- ScaleData(object = merged_obj)
  .logger("Run PCA")
  merged_obj <- RunPCA(merged_obj, verbose = F, ...)
  #No idea why this is needed, but without the next two lines it doesnt work...
  merged_obj <- JoinLayers(merged_obj)
  merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$orig.ident)

  return(merged_obj)
}



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
#' @import dplyr
#' @import ggplot2
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

  .logger("Running Celltypist")

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
  Seu_object@meta.data$autoAnnot <- labels$majority_voting

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

  top_cluster$label <- factor(top_cluster$label, levels = unique(sort(top_cluster$label, decreasing = T)))
  top_cluster$cluster <- factor(top_cluster$cluster, levels = top_cluster$cluster[order(top_cluster$label, decreasing = T)])

  .logger("Creating probality plot")

  pmain <- ggplot(top_cluster, aes(x = cluster, y = label)) +
    geom_point(aes(size = pct.exp, fill = prob), shape=21, color="black", stroke=0.3) +
    scale_fill_gradient2(low = "royalblue3",mid = "white", high = "firebrick",midpoint=0.5, limits = c(0, 1)) +
    scale_size(range = c(2,10), breaks = c(20,40,60,80,100), limits = c(0,100)) +
    DOtools:::theme_box() +
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
#' The new clustering results are stored in a metadata column called `annotation_recluster`.
#' Suitable for improving cluster precision and granularity after initial clustering.
#' @param Seu_object The seurat object
#' @param over_clustering Column in metadata in object with clustering assignments for cells, default seurat_clusters
#' @param res Resolution for the new clusters, default 0.5
#' @param algorithm Set one of the available algorithms found in FindSubCLuster function, default = 4: leiden
#' @return a Seurat Object with new clustering named annotation_recluster
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
                             over_clustering = "seurat_clusters",
                             res = 0.5,
                             algorithm=4,
                             graph.name="RNA_snn"){
  #TODO add an argument for a new shared nearest neighbor call so you can use this function with other integration methods in the object otherwise it will use the latest calculated shared nearest neighbors and you need to do it outside of the function
  if (is.null(Seu_object@meta.data[[over_clustering]])) {
    stop("No seurat clusters defined, please run FindClusters before Reclustering, or fill the slot with a clustering")
  }
  Idents(Seu_object) <- over_clustering

  Seu_object$annotation_recluster <- as.vector(Seu_object@meta.data[[over_clustering]])

  pb <- progress::progress_bar$new(
    total = length(unique(Seu_object@meta.data[[over_clustering]])),
    format = "  Reclustering [:bar] :percent eta: :eta"
  )

  #newline prints when the function exits (to clean up console)
  on.exit(cat("\n"))

  for (cluster in unique(Seu_object@meta.data[[over_clustering]])) {
    pb$tick()
    Seu_object <- FindSubCluster(Seu_object,
                                 cluster = as.character(cluster),
                                 graph.name = graph.name,
                                 algorithm = algorithm,
                                 resolution = res)

    cluster_cells <- rownames(Seu_object@meta.data)[Seu_object@meta.data[[over_clustering]] == cluster]
    Seu_object$annotation_recluster[cluster_cells] <- Seu_object$sub.cluster[cluster_cells]
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
#' @author Mariano Ruz Jurado
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


#' @title DO.PyEnv
#' @description Sets up or connects to a conda Python environment for use with DOtools.
#' If no environment path is provided, it will create one at `~/.venv/DOtools` and install required Python packages:
#' `scvi-tools`, `celltypist`, and `scanpro`.
#'
#' @param conda_path character string specifying the path to an existing or new conda environment.
#'
#' @examples
#' \dontrun{
#' # Automatically create DOtools environment at ~/.venv/DOtools if it doesn't exist
#' DO.PyEnv()
#'
#' # Use an existing conda environment at a custom location
#' DO.PyEnv(conda_path = "~/miniconda3/envs/my_dotools_env")
#' }
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
    conda_args2 <- c("conda","run", "-p", file.path(path.expand("~"), ".venv/DOtools"), "pip", "install", "scvi-tools", "celltypist", "scanpro")

    tryCatch({
      system2(conda_args1[1], args = conda_args1[-1], stdout = F, stderr = F)
      system2(conda_args2[1], args = conda_args2[-1], stdout = F, stderr = F)
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
#' @param bash_script Path to the bash script that runs CellBender.
#'
#'
#' @import DropletUtils
#'
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
                           bash_script = system.file("bash", "_run_CellBender.sh", package = "DOtools")) {

  # Check input paths
  stopifnot(file.exists(cellranger_path))
  stopifnot(file.exists(output_path))

  # Warnings and logs
  if (estimator_multiple_cpu)
    message("Warning: estimator_multiple_cpu is TRUE. Not recommended for large datasets (>20–30k cells).")
  if (epochs > 150)
    message("Warning: Training for more than 150 epochs may lead to overfitting.")
  if (!cuda)
    message("Warning: Running without CUDA (GPU) may significantly increase run time.")

  # Handle conda environment
  if (is.null(conda_path)) {
    conda_path <- file.path(path.expand("~"), ".venv/cellbender")
  } else {
    conda_path <- path.expand(conda_path)
  }

  if (!dir.exists(conda_path)) {
    message("Creating conda environment for CellBender...")
    conda_args1 <- c("conda","create","-y", "-p", file.path(path.expand("~"), ".venv/cellbender"), "python=3.7")
    conda_args2 <- c("conda","run", "-p", file.path(path.expand("~"), ".venv/cellbender"), "pip", "install", "cellbender", "lxml_html_clean")
    tryCatch({
      system2(conda_args1[1], args = conda_args1[-1], stdout = F, stderr = F)
      system2(conda_args2[1], args = conda_args2[-1], stdout = F, stderr = F)
    }, error = function(e) {
      stop("Failed to create conda environment. Provide a valid environment path or fix installation.")
    })
  } else {
    message(sprintf("Using existing conda environment at: %s", conda_path))
  }

  # Detect samples
  if (is.null(samplenames)) {
    samples <- list.dirs(cellranger_path, full.names = FALSE, recursive = FALSE)
  } else {
    samples <- samplenames
  }

  message(sprintf("Running CellBender for %d sample(s)", length(samples)))

  # Loop through each sample
  for (sample in samples) {
    # Load HDF5 with reticulate or scanpy for expected_cells, total_droplets

    # Run one by one but sequentially
    # Estimate the number of cells to be used as upper and lower bound
    h5_file <- file.path(cellranger_path, sample, "outs", "raw_feature_bc_matrix.h5")
    tdata <- DropletUtils::read10xCounts(h5_file)

    result_barcoderanks = DO.BarcodeRanks(tdata)

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

    if (cuda) cmd <- c(cmd, "--cuda")
    if (log) cmd <- c(cmd, "--log")
    if (estimator_multiple_cpu) cmd <- c(cmd, "--estimator_multiple_cpu")

    # Run command
    message(sprintf("Running CellBender for sample: %s", sample))
    tryCatch({
      system2(cmd[1], args = cmd[-1], stdout = T, stderr = T)
    }, error = function(e) {
      message(sprintf("Error running CellBender for sample %s: %s", sample, e$message))
    })
  }

  message("Finished running CellBender.")
  invisible(NULL)
}



#' @title DO.scVI
#' @description This function will run the scVI Integration from the scVI python package. It includes all parameters from the actual python package
#' and runs it by using an internal python script. The usage of a gpu is incorporated and highly recommended.
#'
#' @param Seu_object Seurat object with annotation in meta.data
#' @param batch_key: meta data column with batch information.
#' @param layer_counts: layer with counts. Raw counts are required.
#' @param layer_logcounts: layer with log-counts. Log-counts required for calculation of HVG.
#' @param categorical_covariates: meta data column names with categorical covariates for scVI inference.
#' @param continuos_covariates: meta data  column names with continuous covariates for scVI inference.
#' @param n_hidden: number of hidden layers.
#' @param n_latent: dimensions of the latent space.
#' @param n_layers: number of layers.
#' @param dispersion: dispersion mode for scVI.
#' @param gene_likelihood: gene likelihood.
#' @param get_model: return the trained model.
#' @param ... additional arguments for `scvi.model.SCVI`.
#'
#' @import Seurat
#' @import reticulate
#' @import zellkonverter
#'
#'
#' @examples
#' \dontrun{
#' # Run scVI using the 'orig.ident' column as the batch key
#' Seu_object <- DO.scVI(Seu_object, batch_key = "orig.ident")
#' }
#'
#' @return Seurat Object with dimensionality reduction from scVI
#' @export
DO.scVI <- function(Seu_object,
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
                    get_model=FALSE,
                    ...){

  #check if Variable Features are defined in the object
  HVG <- VariableFeatures(Seu_object)
  if (is.null(HVG)) {
    stop("No highly variable Features found in the object! Please run FindVariableFeatures()!")
  }

  #Subset the object to the HVG
  Seu_object_sub <- subset(Seu_object, features = HVG)

  #Conversion of Seu_object to anndata through zellkonverter
  SCE_object <- as.SingleCellExperiment(Seu_object_sub)
  anndata_object <- zellkonverter::SCE2AnnData(SCE_object)
  anndata_object$layers['counts'] <- anndata_object$X # set

  #source PATH to python script in install folder
  path_py <- system.file("python", "scVI.py", package = "DOtools")
  source_python(path_py)

  run_scvi(adata = anndata_object,
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
           get_model = get_model,
           ...)

  scvi_embedding <- anndata_object$obsm[["X_scVI"]]
  rownames(scvi_embedding) <- colnames(Seu_object)

  scVI_reduction <- suppressWarnings(Seurat::CreateDimReducObject(
    embeddings = scvi_embedding,
    key = "scVI_",
    assay = DefaultAssay(Seu_object_sub)
  ))

  Seu_object@reductions[["scVI"]] <- scVI_reduction

  return(Seu_object)
}


#' @title DO.TransferLabel
#' @description Transfers cell-type annotations from a re-annotated subset of a Seurat object
#' back to the full Seurat object. This is useful when clusters have been refined
#' or re-labeled in a subset and need to be reflected in the original object.
#'
#' @param Seu_object Seurat object with annotation in meta.data
#' @param Subset_obj subsetted Seurat object with re-annotated clusters
#' @param annotation_column column name in meta.data with annotation
#' @param subset_annotation column name in meta.data with annotation in the subsetted object
#'
#'
#' @import Seurat
#'
#'
#' @examples
#' \dontrun{
#'
#' Seu_obj <- DO.TransferLabel(Seu_obj,
#'                             Subset_obj,
#'                             annotation_column="annotation",
#'                             subset_annotation="annotation"
#'                            )
#' }
#'
#'
#' @return Seurat Object with transfered labels
#' @export

DO.TransferLabel <- function(Seu_obj,
                             Subset_obj,
                             annotation_column,
                             subset_annotation){

  #Get annotation with barcodes as rownames
  annotation_pre <- Seu_obj@meta.data[annotation_column]
  annotation_pre[[annotation_column]] <- as.character(annotation_pre[[annotation_column]])
  annotation_subset <- Subset_obj@meta.data[subset_annotation]
  annotation_subset[[subset_annotation]] <- as.character(annotation_subset[[subset_annotation]])

  #assign the new labels
  barcodes <- rownames(annotation_pre)[rownames(annotation_pre) %in% rownames(annotation_subset)]
  annotation_pre[rownames(annotation_pre) %in% barcodes,] <- annotation_subset[[subset_annotation]]

  Seu_obj@meta.data[[annotation_column]] <- factor(annotation_pre$annotation)

  return(Seu_obj)

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
#' \dontrun{
#' #GO analysis on differential expression results:
#' results <- go_analysis(df_DGE = deg_results,
#'                        gene_key = "gene",
#'                        pval_key = "adj_pval",
#'                        log2fc_key = "log2FC",
#'                        species = "Human")
#'
#' #Or save the results to a file:
#' go_analysis(df_DGE = deg_results,
#'             gene_key = "gene",
#'             pval_key = "p_val_adj",
#'             log2fc_key = "log2FC",
#'             path = "results/",
#'             filename = "experiment.xlsx",
#'             species = "Mouse")
#' }
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
    .logger("Write results to: ", out_file)
    return(invisible(NULL))
  } else {
    return(combined_res)
  }
}


#' @title DO.BarcodeRanks
#' @author Mariano Ruz Jurado
#' @description Given a raw count matrix (e.g. from a CellRanger HDF5 file), estimate the number of expected cells and droplets
#' using the knee and inflection points from barcodeRanks.
#'
#' @param SCE_object A Single cell experiment object.
#'
#' @return A named numeric vector: `c(xpc_cells = ..., total_cells = ...)`
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
#' @param Seu_obj A Seurat object containing single-cell RNA-seq data.
#' @param layer A character string specifying the assay layer to use (default is "counts").
#' @param features A character vector of length 3 indicating the feature names to plot. Default is c("nFeature_RNA", "nCount_RNA", "pt_mito").
#'
#' @return A `ggplot` object arranged in a single row showing violin plots for the specified features with overlaid boxplots.
#'
#' @import Seurat
#' @import ggplot2
#'
#' @keywords internal
.QC_Vlnplot <- function(Seu_obj, id, layer="counts", features=c("nFeature_RNA","nCount_RNA","pt_mito")){
  p1<- VlnPlot(Seu_obj,layer = "counts", features = features[1], ncol = 1, pt.size = 0, cols = "grey")+
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

  p2<- VlnPlot(Seu_obj,layer = "counts", features = features[2], ncol = 1, pt.size = 0, cols = "grey")+
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

  p3<-VlnPlot(Seu_obj,layer = "counts", features = features[3], ncol = 1, pt.size = 0, cols = "grey")+
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
