#' @title DO.Import
#' @author Mariano Ruz Jurado & David John
#' @description
#' Imports and processes single-cell RNA-seq data from various formats
#' (10x Genomics, CellBender, or CSV), performs quality control (QC), filtering
#' , normalization, variable gene selection, and optionally detects doublets.
#' Returns a merged and processed Seurat or SCE object ready for downstream
#' analysis.
#'
#' @param pathways A character vector of paths to directories or files
#' containing raw expression matrices.
#' @param ids A character vector of sample identifiers, matching the order
#' of `pathways`.
#' @param minCellGenes Integer. Minimum number of cells a gene must be expressed
#'  in to be retained. Default is 5.
#' @param FilterCells Logical. If `TRUE`, applies QC filtering on cells based on
#'  mitochondrial content, counts, and feature thresholds. Default is `TRUE`.
#' @param cut_mt Numeric. Maximum allowed mitochondrial gene proportion per
#' cell. Default is 0.05.
#' @param min_counts Numeric. Minimum UMI count threshold
#' (optional, used only if `low_quantile` is `NULL`).
#' @param max_counts Numeric. Maximum UMI count threshold
#' (optional, used only if `high_quantile` is `NULL`).
#' @param min_genes Numeric. Minimum number of genes detected per cell to
#' retain. Optional.
#' @param max_genes Numeric. Maximum number of genes detected per cell to
#' retain. Optional.
#' @param low_quantile Numeric. Quantile threshold (0 to 1) to filter low UMI
#' cells (used if `min_counts` is `NULL`).
#' @param high_quantile Numeric. Quantile threshold (0 to 1) to filter high UMI
#' cells (used if `max_counts` is `NULL`).
#' @param DeleteDoublets Logical. If `TRUE`, doublets are detected and removed
#' using `scDblFinder`. Default is `TRUE`.
#' @param include_rbs Logical. If `TRUE`, calculates ribosomal gene content in
#' addition to mitochondrial content. Default is `TRUE`.
#' @param Seurat Logical. If `TRUE`, returns Seurat object otherwise SCE object.
#' @param ... Additional arguments passed to `RunPCA()`.
#'
#' @return A merged Seurat or SCE object containing all samples, with
#' normalization, QC, scaling, PCA, and optional doublet removal applied.
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
                      minCellGenes = 5,
                      FilterCells = TRUE,
                      cut_mt = .05,
                      min_counts = NULL,
                      max_counts = NULL,
                      min_genes = NULL,
                      max_genes = NULL,
                      low_quantile = NULL,
                      high_quantile = NULL,
                      DeleteDoublets = TRUE,
                      include_rbs = TRUE,
                      Seurat = TRUE,
                      ...) {
  object_list <- list()
  for (i in seq_along(pathways)) {
    id <- ids[i]
    pathway <- pathways[i]

    .logger(paste0("Sample: ", id))

    # case when the user provides just the path to the folder with the .h5
    file_path <- ifelse(
      dir.exists(pathway),
      list.files(pathway, pattern = "*filtered.h5", full.names = TRUE),
      pathway
    )
    outPath <- ifelse(dir.exists(pathway),
      pathway,
      dirname(pathway)
    )

    # Try different methods to read provided path
    .logger("Read matrix")
    mtx <- tryCatch(
      {
        # Try CellRanger 10X
        tmp <- Seurat::Read10X_h5(file_path)
        if (is.list(tmp)) {
          tmp <- tmp[[grep("gene",
            names(tmp),
            value = TRUE,
            ignore.case = TRUE
          )]]
        }
        tmp
      },
      error = function(x) {
        warning("Matrix not produced by Cellranger. Trying CellBender read...")
        tryCatch(
          {
            scCustomize::Read_CellBender_h5_Mat(file_path)
          },
          error = function(x2) {
            warning(
              "Matrix is not produced by CellBender. Trying simple CSV ",
              "read..."
            )
            tryCatch(
              {
                read.csv(
                  file_path,
                  header = TRUE,
                  row.names = 1,
                  check.names = FALSE
                )
              },
              error = function(e) {
                stop(
                  "Failed to read file as CellRanger, CellBender, or CSV. ",
                  "Check file format."
                )
              }
            )
          }
        )
      }
    )

    # Metrics file
    tdy <- format(Sys.Date(), "%y%m%d")
    met_file <- paste0(tdy, "_metrics_", id, ".xlsx")
    df_met <- data.frame(
      QC_Step = "Input_Shape",
      nCells = ncol(mtx),
      nFeatures = nrow(mtx)
    )

    # Create object + Filtering by minimum Genes per cell
    .logger("Create Single Cell Object")
    sce_object <- SeuratObject::CreateSeuratObject(
      counts = mtx,
      project = id,
      min.cells = minCellGenes
    )

    df_met[2, ] <- c("Rm_undetected_genes", ncol(sce_object), nrow(sce_object))

    sce_object$sample <- id # some naming
    # automatized condition settings
    .logger(paste0("Setting condition in object to: ", sub("[-|_].*", "", id)))
    sce_object$condition <- sub("[-|_].*", "", id)

    # Set the filter on mito/ribo genes
    if (include_rbs == TRUE) {
      .logger("Get Mitochondrial+Ribosomal content")
    } else {
      .logger("Get Mitochondrial content")
    }

    if (include_rbs == TRUE) {
      sel_ribofeatures <- grep("^(RPS|RPL)",
        rownames(sce_object),
        value = TRUE,
        ignore.case = TRUE
      )
      pt_ribo <- Matrix::colSums(
        GetAssayData(sce_object, layer = "counts")[sel_ribofeatures, ]
      ) / Matrix::colSums(
        GetAssayData(sce_object, layer = "counts")
      )
      sce_object$pt_ribo <- pt_ribo
    }
    flt_mitofeatures <- grep("^MT-",
      rownames(sce_object),
      value = TRUE,
      ignore.case = TRUE
    )
    if (length(flt_mitofeatures) < 1) {
      warning("Could not find MT genes")
    }
    pt_mito <- Matrix::colSums(
      GetAssayData(sce_object, layer = "counts")[flt_mitofeatures, ]
    ) / Matrix::colSums(
      GetAssayData(sce_object, layer = "counts")
    )
    sce_object$pt_mito <- pt_mito

    .logger("Create QC images")

    # write QC to file
    prefilter_plot <- .suppressDeprecationWarnings(
      .QC_Vlnplot(sce_object = sce_object, id, layer = "counts")
    )
    ggsave(
      plot = prefilter_plot,
      filename = paste0(outPath, "/QC_Plots_prefiltered.svg"),
      width = 10,
      height = 6
    )

    if (FilterCells == TRUE) {
      .logger("Start Filtering")

      # check if absolute values are set for counts and quantile is set too
      if (!is.null(min_counts) && !is.null(low_quantile)) {
        stop(
          "Both filtering minimum counts by absolute values and quantile is ",
          "set to Value!"
        )
      }

      if (!is.null(max_counts) && !is.null(high_quantile)) {
        stop(
          "Both filtering maximum counts by absolute values and quantile ",
          "is set to Value!"
        )
      }

      if (is.null(min_counts) && is.null(low_quantile)) {
        stop("Both minimum count filtering methods are not set!")
      }

      if (is.null(max_counts) && is.null(high_quantile)) {
        stop("Both maximum count filtering methods are not set!")
      }

      # Remove all cells with below gene threshold
      if (!is.null(min_genes)) {
        sce_object <- subset(sce_object, subset = nFeature_RNA > min_genes)
      }

      if (!is.null(max_genes)) {
        sce_object <- subset(sce_object, subset = nFeature_RNA < max_genes)
      }
      if (!is.null(min_genes) || !is.null(max_genes)) {
        df_met[3, ] <- c(
          "Rm_cells_based_on_genes",
          ncol(sce_object),
          nrow(sce_object)
        )
      }

      # Remove all cells with high MT content
      sce_object <- subset(sce_object, subset = pt_mito < cut_mt)
      df_met[4, ] <- c("Rm_cell_high_MT", ncol(sce_object), nrow(sce_object))

      # Remove by absolute values or quantiles
      if (!is.null(min_counts) && is.null(low_quantile)) {
        sce_object <- subset(sce_object, subset = nCount_RNA > min_counts)
      } else {
        # Calculate the Quantile for the lower end
        Quality <- data.frame(
          UMI = sce_object$nCount_RNA,
          nGene = sce_object$nFeature_RNA,
          label = factor(sce_object$sample),
          pt_mit_rib = sce_object$pt_mito
        )

        Quantile.low.UMI <- Quality %>%
          group_by(label) %>%
          summarise(
            UMI = list(tibble::enframe(quantile(UMI, probs = low_quantile)))
          ) %>%
          tidyr::unnest(cols = c(UMI))

        # Subset
        sce_object <- subset(sce_object,
          subset = nCount_RNA > Quantile.low.UMI$value
        )
      }

      # Remove by absolute values
      if (!is.null(max_counts) && is.null(high_quantile)) {
        sce_object <- subset(sce_object, subset = nCount_RNA < max_counts)
      } else {
        # Calculate the Quantile for the lower end
        Quality <- data.frame(
          UMI = sce_object$nCount_RNA,
          nGene = sce_object$nFeature_RNA,
          label = factor(sce_object$sample),
          pt_mit_rib = sce_object$pt_mito
        )

        Quantile.high.UMI <- Quality %>%
          group_by(label) %>%
          summarise(
            UMI = list(tibble::enframe(quantile(UMI, probs = high_quantile)))
          ) %>%
          tidyr::unnest(cols = c(UMI))

        # Subset
        sce_object <- subset(sce_object,
          subset = nCount_RNA < Quantile.high.UMI$value
        )
      }
    }


    # save metric files
    write.xlsx(df_met, file = paste0(outPath, "/", met_file))

    # write QC after filtering to file
    postfilter_plot <- .QC_Vlnplot(
      sce_object = sce_object,
      id,
      layer = "counts"
    )
    ggsave(
      plot = postfilter_plot,
      filename = paste0(outPath, "/QC_Plots_postfiltered.svg"),
      width = 10,
      height = 6
    )

    # Preprocess steps Seurat
    .logger("Running Normalisation")
    sce_object <- NormalizeData(object = sce_object, verbose = FALSE)

    .logger("Running Variable Gene Detection")
    sce_object <- FindVariableFeatures(
      object = sce_object,
      selection.method = "vst",
      nfeatures = 2000,
      verbose = FALSE
    )

    if (DeleteDoublets == TRUE) {
      .logger("Running scDblFinder")
      SCE_obj <- suppressWarnings(as.SingleCellExperiment(sce_object))
      SCE_obj <- scDblFinder::scDblFinder(SCE_obj)
      sce_object$scDblFinder_score <- SCE_obj$scDblFinder.score
      sce_object$scDblFinder_class <- SCE_obj$scDblFinder.class
      sce_object <- subset(sce_object, scDblFinder_class == "singlet")
    }

    object_list[[i]] <- sce_object # capture the final objects
  }

  # concatenate objects
  .logger("Merging objects")
  merged_obj <- Reduce(function(x, y) merge(x, y), object_list)
  .logger("Running ScaleData")
  merged_obj <- ScaleData(object = merged_obj)
  .logger("Run PCA")
  merged_obj <- RunPCA(merged_obj, verbose = FALSE, ...)
  # No idea why this is needed for seurat,
  # but without the join and split doesnt work...
  merged_obj <- JoinLayers(merged_obj)

  if (Seurat == FALSE) {
    merged_obj <- as.SingleCellExperiment(merged_obj)
  } else {
    merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$orig.ident)
  }

  return(merged_obj)
}

#' @title .QC_Vlnplot
#' @author Mariano Ruz Jurado
#' @description
#' Generates violin plots for common quality control (QC) metrics of single-cell
#' RNA-seq data from a Seurat object. The function displays three violin plots
#' for the number of detected genes per cell (`nFeature_RNA`), total UMI counts
#' per cell (`nCount_RNA`), and mitochondrial gene content percentage
#' (`pt_mito`). Useful for visual inspection of QC thresholds and outliers.
#'
#' @param sce_object A Seurat object containing single-cell RNA-seq data.
#' @param layer A character string specifying the assay layer to use
#' (default is "counts").
#' @param features A character vector of length 3 indicating the feature names
#' to plot. Default is c("nFeature_RNA", "nCount_RNA", "pt_mito").
#'
#' @return A `ggplot` object arranged in a single row showing violin plots for
#' the specified features with overlaid boxplots.
#'
#' @import Seurat
#' @import ggplot2
#' @rdname dot-QC_Vlnplot
#' @keywords internal
.QC_Vlnplot <- function(sce_object,
                        id,
                        layer = "counts",
                        features = c("nFeature_RNA", "nCount_RNA", "pt_mito")) {
  p1 <- VlnPlot(
    sce_object,
    layer = "counts",
    features = features[1],
    ncol = 1,
    pt.size = 0,
    cols = "grey"
  ) +
    geom_boxplot(
      width = 0.25,
      outlier.shape = NA,
      alpha = 0.5,
      color = "darkred",
      size = 0.4
    ) +
    theme_light() +
    xlab("") +
    theme(
      plot.title = element_text(
        face = "bold",
        color = "black",
        hjust = 0.5,
        size = 14
      ),
      axis.title.y = element_text(
        face = "bold",
        color = "black",
        size = 14
      ),
      axis.text.x = element_blank(),
      axis.text.y = element_text(
        face = "bold",
        color = "black",
        hjust = 1,
        size = 14
      ),
      axis.title.x = element_blank(),
      legend.position = "none"
    )

  p2 <- VlnPlot(
    sce_object,
    layer = "counts",
    features = features[2],
    ncol = 1,
    pt.size = 0,
    cols = "grey"
  ) +
    geom_boxplot(
      width = 0.25,
      outlier.shape = NA,
      alpha = 0.5,
      color = "darkred",
      size = 0.4
    ) +
    theme_light() +
    xlab("") +
    theme(
      plot.title = element_text(
        face = "bold",
        color = "black",
        hjust = 0.5,
        size = 14
      ),
      axis.title.y = element_text(
        face = "bold",
        color = "black",
        size = 14
      ),
      axis.text.x = element_blank(),
      axis.text.y = element_text(
        face = "bold",
        color = "black",
        hjust = 1,
        size = 14
      ),
      axis.title.x = element_blank(),
      legend.position = "none"
    )


  p3 <- VlnPlot(
    sce_object,
    layer = "counts",
    features = features[3],
    ncol = 1,
    pt.size = 0,
    cols = "grey"
  ) +
    geom_boxplot(
      width = 0.25,
      outlier.shape = NA,
      alpha = 0.5,
      color = "darkred",
      size = 0.4
    ) +
    theme_light() +
    xlab("") +
    theme(
      plot.title = element_text(
        face = "bold",
        color = "black",
        hjust = 0.5,
        size = 14
      ),
      axis.title.y = element_text(
        face = "bold",
        color = "black",
        size = 14
      ),
      axis.text.x = element_blank(),
      axis.text.y = element_text(
        face = "bold",
        color = "black",
        hjust = 1,
        size = 14
      ),
      axis.title.x = element_blank(),
      legend.position = "none"
    )

  gg_plot <- ggarrange(p1, p2, p3, nrow = 1)
  gg_plot <- annotate_figure(gg_plot, top = grid::grobTree(
    grid::rectGrob(
      gp = grid::gpar(fill = "white", col = NA),
      height = unit(1.5, "lines")
    ),
    grid::textGrob(
      id,
      gp = grid::gpar(
        fontface = "bold",
        col = "darkred",
        fontsize = 18
      ),
      hjust = 0.35,
      vjust = 0.6 # align horizontally at 30%
    )
  ))

  return(gg_plot)
}
