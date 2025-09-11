#' @author Mariano Ruz Jurado & David Rodriguez Morales
#' @title DO Heatmap of the mean expression of genes across a groups
#' @description Wrapper around heatmap.py, which generates a heatmap of
#' showing the average nUMI for a set of genes in different groups.
#' Differential gene expression analysis between the different groups can be performed.
#'
#' @param sce_object SCE object or Seurat with meta.data
#' @param assay_normalized Assay with raw counts
#' @param group_by meta data column name with categorical values
#' @param groups_order order for the categories in the group_by
#' @param features gene names or continuous value in meta data
#' @param z_score apply z-score transformation
#' @param path path to save the plot
#' @param filename name of the file
#' @param swap_axes whether to swap the axes or not
#' @param cmap color map
#' @param title title for the main plot
#' @param title_fontprop font properties for the title (e.g., 'weight' and 'size')
#' @param clustering_method clustering method to use when hierarchically clustering the x and y-axis
#' @param clustering_metric metric to use when hierarchically clustering the x and y-axis
#' @param cluster_x_axis hierarchically clustering the x-axis
#' @param cluster_y_axis hierarchically clustering the y-axis
#' @param axs matplotlib axis
#' @param figsize figure size
#' @param linewidth line width for the border of cells
#' @param ticks_fontdict font properties for the x and y ticks (e.g.,  'weight' and 'size')
#' @param xticks_rotation rotation of the x-ticks
#' @param yticks_rotation rotations of the y-ticks
#' @param vmin minimum value
#' @param vcenter center value
#' @param vmax maximum value
#' @param legend_title title for the color bar
#' @param add_stats add statistical annotation, will add a square with an '*' in the center if the expression is significantly different in a group with respect to the others
#' @param df_pvals dataframe with the p-values, should be gene x group or group x gene in case of swap_axes is False
#' @param stats_x_size size of the asterisk
#' @param square_x_size size and thickness of the square percentual, vector
#' @param test test to use for test for significance
#' @param pval_cutoff cutoff for the p-value
#' @param log2fc_cutoff minimum cutoff for the log2FC
#' @param only_pos if set to TRUE, only use positive genes in the condition
#' @param square whether to make the cell square or not
#' @param showP if set to false return a dictionary with the axis
#' @param logcounts whether the input is logcounts or not
#' @return Depending on ``showP``, returns the plot if set to `TRUE` or a dictionary with the axes.
#'
#'
#' @import Seurat
#' @importFrom basilisk basiliskRun
#'
#' @examples#'
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' DO.Heatmap(
#' sce_object = sce_data,
#' assay_normalized = "RNA",
#' group_by="seurat_clusters",
#' features = rownames(sce_data)[1:10],
#' z_score = NULL,
#' path = NULL,
#' filename = "Heatmap.svg",
#' swap_axes = TRUE,
#' cmap = "Reds",
#' title = NULL,
#' title_fontprop = NULL,
#' clustering_method = "complete",
#' clustering_metric = "euclidean",
#' cluster_x_axis = FALSE,
#' cluster_y_axis = FALSE,
#' axs = NULL,
#' figsize = c(5, 6),
#' linewidth = 0.1,
#' ticks_fontdict = NULL,
#' xticks_rotation = 45,
#' yticks_rotation = NULL,
#' vmin = 0.0,
#' vcenter = NULL,
#' vmax = NULL,
#' legend_title = "LogMean(nUMI)\nin group",
#' add_stats = TRUE,
#' df_pvals = NULL,
#' stats_x_size = NULL,
#' square_x_size = NULL,
#' test = "wilcox",
#' pval_cutoff = 0.05,
#' log2fc_cutoff = 0,
#' only_pos = TRUE,
#' square = TRUE,
#' showP = FALSE,
#' logcounts = TRUE
#' )
#'
#'
#' @export
DO.Heatmap <- function(
    sce_object,
    assay_normalized = "RNA",
    group_by = "seurat_clusters",
    groups_order = NULL,
    features,
    z_score = NULL,
    path = NULL,
    filename = "Heatmap.svg",
    swap_axes = TRUE,
    cmap = "Reds",
    title = NULL,
    title_fontprop = NULL,
    clustering_method = "complete",
    clustering_metric = "euclidean",
    cluster_x_axis = FALSE,
    cluster_y_axis = FALSE,
    axs = NULL,
    figsize = c(5, 6),
    linewidth = 0.1,
    ticks_fontdict = NULL,
    xticks_rotation = NULL,
    yticks_rotation = NULL,
    vmin = 0.0,
    vcenter = NULL,
    vmax = NULL,
    legend_title = "LogMean(nUMI)\nin group",
    add_stats = TRUE,
    df_pvals = NULL,
    stats_x_size = NULL,
    square_x_size = NULL,
    test = "wilcox",
    pval_cutoff = 0.05,
    log2fc_cutoff = 0,
    only_pos = TRUE,
    square = TRUE,
    showP = TRUE,
    logcounts = TRUE
) {

  #support for Seurat objects
  if (is(sce_object, "Seurat")) {
    DefaultAssay(sce_object) <- assay_normalized
    sce_object <- Seurat::as.SingleCellExperiment(sce_object, assay = assay_normalized)
  }

  #Make Anndata object
  if (!"logcounts" %in% names(sce_object@assays)) {
    stop("logcounts not found in assays of object!")
  }

  if (add_stats == TRUE) {
    if (is.null(df_pvals)) {
      Seu_obj <- as.Seurat(sce_object)
      df_dge <- FindAllMarkers(Seu_obj,
                               features = features,
                               group.by = group_by,
                               min.pct = 0,
                               test.use = test,
                               logfc.threshold = log2fc_cutoff,
                               only.pos = only_pos)

      df_pvals <- data.frame(matrix(1, nrow = length(features), ncol = length(unique(sce_object[[group_by]]))))
      rownames(df_pvals) <- features
      colnames(df_pvals) <- unique(sce_object[[group_by]])

      # Go through each row in df_dge
      for (i in seq_len(nrow(df_dge))) {
        gene <- df_dge$gene[i]
        cluster <- as.character(df_dge$cluster[i])
        pval_adj <- df_dge$p_val_adj[i]
        df_pvals[gene, cluster] <- pval_adj
      }
    }
  }

  #source PATH to python script in install folder
  path_py <- system.file("python", "heatmap.py", package = "DOtools")


  #argument list passed to heatmap inside basilisk
  args <- list(
    sce_object = sce_object,
    group_by = group_by,
    groups_order = groups_order,
    features = features,
    z_score = z_score,
    path = path,
    filename = filename,
    layer = NULL,
    swap_axes = swap_axes,
    cmap = cmap,
    title = title,
    title_fontprop = title_fontprop,
    clustering_method = clustering_method,
    clustering_metric = clustering_metric,
    cluster_x_axis = cluster_x_axis,
    cluster_y_axis = cluster_y_axis,
    axs = axs,
    figsize = figsize,
    linewidth = linewidth,
    ticks_fontdict = ticks_fontdict,
    xticks_rotation = xticks_rotation,
    yticks_rotation = yticks_rotation,
    vmin = vmin,
    vcenter = vcenter,
    vmax = vmax,
    legend_title = legend_title,
    add_stats = add_stats,
    df_pvals = df_pvals,
    stats_x_size = stats_x_size,
    square_x_size = square_x_size,
    pval_cutoff = pval_cutoff,
    square = square,
    showP = showP,
    logcounts = logcounts
  )


  #basilisk implementation
  results <- basilisk::basiliskRun(env = DOtoolsEnv, fun=function(args){

    AnnData_counts <- zellkonverter::SCE2AnnData(args$sce_object, X_name = "logcounts")

    reticulate::source_python(path_py)

    #Initialize matplot package
    plt <- reticulate::import("matplotlib.pyplot")

    #build dicitonary out of square_x_size if specified
    if (is.numeric(square_x_size) && length(square_x_size) == 2) {
      args$square_x_size <- reticulate::dict(weight = square_x_size[1], size = square_x_size[2])
    }


    heatmap(adata = AnnData_counts,
            group_by = args$group_by,
            features = args$features,
            groups_order = args$groups_order,
            z_score = args$z_score,
            path = args$path,
            filename = args$filename,
            layer = args$layer,
            swap_axes = args$swap_axes,
            cmap = args$cmap,
            title = args$title,
            title_fontprop = args$title_fontprop,
            clustering_method = args$clustering_method,
            clustering_metric = args$clustering_metric,
            cluster_x_axis = args$cluster_x_axis,
            cluster_y_axis = args$cluster_y_axis,
            axs = args$axs,
            figsize = args$figsize,
            linewidth = args$linewidth,
            ticks_fontdict = args$ticks_fontdict,
            xticks_rotation = args$xticks_rotation,
            yticks_rotation = args$yticks_rotation,
            vmin = args$vmin,
            vcenter = args$vcenter,
            vmax = args$vmax,
            legend_title = args$legend_title,
            add_stats = args$add_stats,
            df_pvals = args$df_pvals,
            stats_x_size = args$stats_x_size,
            square_x_size = args$square_x_size,
            pval_cutoff = args$pval_cutoff,
            square = args$square,
            showP = args$showP,
            logcounts = args$logcounts
    )
  },args=args)


}
