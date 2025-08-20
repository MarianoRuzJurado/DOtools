#' @author Mariano Ruz Jurado
#' @title DO Correlation Plot for visualizing similarity between categories
#' @description Generates a correlation heatmap from expression data to visualize similarity across sample groups.
#' Allows customization of plot type, correlation method, and color scaling using the ggcorrplot2 and ggplot2 architectures.
#' Ideal for comparing transcriptional profiles between conditions or clusters.
#' @param sce_object Seurat or SCE Object
#' @param group_by Column to aggregate the expression over it, default "orig.ident"
#' @param assay Assay in object to use, default "RNA"
#' @param features What genes to include by default all, default "None"
#' @param method Correlation method, default "spearman"
#' @param plotdesign Plot design, default "circle"
#' @param plottype Show the full plot or only half of it, default "full"
#' @param auto_limits Automatically rescales the colour bar based on the values in the correlation matrix, default "TRUE"
#' @param outline.color the outline color of square or circle. Default value is "white".
#' @param colormap Defines the colormap used in the plot, default c("royalblue4", "royalblue2","firebrick","firebrick4")
#' @param lab_size Size to be used for the correlation coefficient labels. used when lab = TRUE.
#' @param lab logical value. If TRUE, add correlation coefficient on the plot.
#' @param lab_col color to be used for the correlation coefficient labels. used when lab = TRUE.
#' @param axis_size_x Controls x labels size
#' @param axis_size_y Controls y labels size
#' @param ... Additionally arguments passed to ggcorrplot function
#'
#'
#' @return ggplot2
#'
#' @import ggcorrplot
#' @import ggplot2
#' @import Seurat
#'
#' @examples
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' DO.Correlation(
#'   sce_object = sce_data,
#'   group_by = "orig.ident",
#'   assay = "RNA",
#'   features = NULL,
#'   method = "spearman",
#'   plotdesign = "square",
#'   plottype = "full",
#'   auto_limits = TRUE,
#'   outline.color = "white",
#'   colormap = c("royalblue4", "lightsteelblue", "tomato","firebrick4"),
#'   lab_size = 10,
#'   lab = TRUE,
#'   lab_col = "white"
#' )
#'
#' @export
DO.Correlation <- function(sce_object,
                           group_by="orig.ident",
                           assay="RNA",
                           features=NULL,
                           method="spearman",
                           plotdesign="square",
                           plottype="full",
                           auto_limits=TRUE,
                           outline.color="white",
                           colormap=c("royalblue4", "lightsteelblue", "tomato","firebrick4"),
                           lab_size=10,
                           lab=TRUE,
                           lab_col="white",
                           axis_size_x = 12,
                           axis_size_y = 12,
                           ...){

  #support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    sce_object <- .suppressDeprecationWarnings(as.Seurat(sce_object))
  }

  #Aggregate Expression, creating Pseudobulk
  corr_frame <- AggregateExpression(sce_object,
                                    assays = assay,
                                    features = features,
                                    group.by = group_by)[[assay]]

  corr_frame <- as.data.frame(corr_frame)

  #Calculate correlation
  corr_df <- cor(corr_frame, method = method)


  # Auto color limits based on matrix
  if (auto_limits == TRUE) {
    min_val <- min(corr_df, na.rm = TRUE)*0.99
    max_val <- max(corr_df, na.rm = TRUE)*1.01
  } else {
    min_val <- -1
    max_val <- 1
  }

  pmain <- ggcorrplot(corr_df,
                      method = plotdesign,
                      type = plottype,
                      colors = colormap,
                      outline.color = outline.color,
                      lab_col = lab_col,
                      lab_size = lab_size,
                      lab = lab,
                      ...)+
    scale_fill_gradientn(colors = colormap,
                         limits = c(min_val, max_val),
                         name = paste0(tools::toTitleCase(method), "Correlation"))+
    ggplot2::theme(axis.text = ggplot2::element_text(color = "black"),
                   legend.direction = "horizontal",
                   axis.text.x = element_text(color = "black",angle = 0,hjust = 0.5, size = axis_size_x, family = "Helvetica"),
                   axis.text.y = element_text(color = "black",angle=90,,hjust = 0.5, size = axis_size_y, family = "Helvetica"),
                   axis.title = element_text(size = 14, color = "black", family = "Helvetica"),
                   plot.title = element_text(size = 14, hjust = 0.5, face="bold", family = "Helvetica"),
                   plot.subtitle = element_text(size = 14, hjust = 0, family = "Helvetica"),
                   strip.text.x = element_text(size = 14, color = "black", family = "Helvetica", face = "bold"),
                   legend.text = element_text(size = 10, color = "black", family = "Helvetica"),
                   legend.title = element_text(size = 10, color = "black", family = "Helvetica", hjust =0),
                   legend.position = "bottom")

  guides.layer <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = paste0(tools::toTitleCase(method), " \nCorrelation"),
                                                                 title.position = "top",
                                                                 title.hjust = 0.5,
                                                                 barwidth = unit(3.8,"cm"), # changes the width of the color legend
                                                                 barheight = unit(0.5,"cm"),
                                                                 frame.colour = "black",
                                                                 frame.linewidth = 0.3,
                                                                 ticks.colour = "black",
                                                                 order = 1))
  pmain <- pmain + guides.layer

  return(pmain)

}
