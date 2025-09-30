# Polished UMAP function using Dimplot or FeaturePlot function from Seurat
#' @author Mariano Ruz Jurado
#' @title DO.UMAP
#' @description Creates a polished UMAP plot using Seurat's DimPlot or
#' FeaturePlot functions. It allows customization of colors, labels, and other
#' plot elements for better visualisation. The function handles both
#' cluster-based visualisations and gene-based visualisations in a UMAP plot.
#' Ideal for refining UMAP outputs with added flexibility and enhanced
#' presentation.
#' @param sce_object The seurat or SCE object
#' @param FeaturePlot Is it going to be a Dimplot or a FeaturePlot?
#' @param features features for Featureplot
#' @param group.by grouping of plot in DImplot and defines in featureplot the
#' labels
#' @param umap_colors what colors to use for UMAP, specify as vector
#' @param text_size Size of text
#' @param label label the clusters on the plot by group.by column
#' @param order Boolean determining whether to plot cells in order of
#' expression.
#' @param plot.title title for UMAP
#' @param legend.position specify legend position
#' @param ... Further arguments passed to DimPlot or FeaturePlot function
#' from Seurat
#' @return Plot with Refined colors and axes
#'
#' @import Seurat
#' @import ggplot2
#'
#' @examples
#' sce_data <-
#'   readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' DO.UMAP(
#'   sce_object = sce_data,
#'   group.by = "seurat_clusters"
#' )
#'
#' DO.UMAP(
#'   sce_object = sce_data,
#'   FeaturePlot = TRUE,
#'   features = c("BAG2", "CD74")
#' )
#'
#' @export
DO.UMAP <- function(sce_object,
                    FeaturePlot = FALSE,
                    features = NULL,
                    group.by = "seurat_clusters",
                    umap_colors = NULL,
                    text_size = 14,
                    label = TRUE,
                    order = TRUE,
                    plot.title = TRUE,
                    legend.position = "none",
                    ...) {
  # support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    sce_object <- as.Seurat(sce_object)
  }

  # Dimplot
  if (FeaturePlot == FALSE) {
    if (is.null(umap_colors)) {
      umap_colors <- rep(
        c(
          "#1f77b4",
          "#ff7f0e",
          "#2ca02c",
          "tomato2",
          "#9467bd",
          "chocolate3",
          "#e377c2",
          "#ffbb78",
          "#bcbd22",
          "#17becf",
          "darkgoldenrod2",
          "#aec7e8",
          "#98df8a",
          "#ff9896",
          "#c5b0d5",
          "#c49c94",
          "#f7b6d2",
          "#c7c7c7",
          "#dbdb8d",
          "#9edae5",
          "sandybrown",
          "moccasin",
          "lightsteelblue",
          "darkorchid",
          "salmon2",
          "forestgreen",
          "bisque"
        ),
        5
      )
    }

    p <- DimPlot(sce_object, group.by = group.by, cols = umap_colors, ...) +
      labs(x = "UMAP1", y = "UMAP2") +
      theme(
        plot.title = element_blank(),
        # text = element_text(face = "bold",size = 20),
        axis.title.x = element_text(size = text_size, family = "Helvetica"),
        axis.title.y = element_text(size = text_size, family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = legend.position,
        legend.text = element_text(face = "bold")
      )

    if (label == TRUE) {
      p <- LabelClusters(p, id = group.by, fontface = "bold", box = FALSE)
    }
    return(p)
  }

  # FeaturePlot
  if (FeaturePlot == TRUE) {
    if (is.null(features)) {
      stop("Please provide any gene names if using FeaturePlot=TRUE.")
    }

    if (is.null(umap_colors)) {
      umap_colors <- c("lightgrey", "red2")
    }

    Idents(sce_object) <- group.by
    p <- FeaturePlot(sce_object,
                     features = features,
                     cols = umap_colors,
                     label = label,
                     order = order,
                     ...
    ) &
      labs(x = "UMAP1", y = "UMAP2") &
      theme(
        axis.title.x = element_text(size = 14, family = "Helvetica"),
        axis.title.y = element_text(size = 14, family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = legend.position,
        legend.text = element_text(face = "bold")
      )

    if (plot.title == FALSE) {
      p <- p & theme(plot.title = element_blank())
    }

    return(p)
  }
}
