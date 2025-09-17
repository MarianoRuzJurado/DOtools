#' @author Mariano Ruz Jurado & David Rodriguez Morales
#' @title DO Heatmap of the mean expression of genes across a groups
#' @description Wrapper around heatmap.py, which generates a heatmap of
#' showing the average nUMI for a set of genes in different groups.
#' Addiional an argumnt can be made to show foldchanges between two conditions.
#' Differential gene expression analysis between the different groups can be performed.
#'
#' @param sce_object SCE object or Seurat with meta.data
#' @param features gene names or continuous value in meta data
#' @param assay_normalized Assay with raw counts
#' @param group_by meta data column name with categorical values
#' @param groups_order order for the categories in the group_by
#' @param value_plot plotted values correspond to expression values or foldchanges
#' @param group_fc if foldchanges specified than the groups must be specified that will be compared
#' @param group_fc_ident_1 Defines the first group in the test
#' @param group_fc_ident_2 Defines the second group in the test
#' @param z_score apply z-score transformation, "group" or "var"
#' @param clip_value Clips the colourscale to the 99th percentile, useful if one gene is driving the colourscale
#' @param max_fc Clips super high foldchanges to this value, so changes can still be appreciated
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
    features,
    assay_normalized = "RNA",
    group_by = "seurat_clusters",
    groups_order = NULL,
    value_plot = "expr",
    group_fc = "condition",
    group_fc_ident_1 = NULL,
    group_fc_ident_2 = NULL,
    clip_value = FALSE,
    max_fc = 5,
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
    sce_object <- .suppressDeprecationWarnings(Seurat::as.SingleCellExperiment(sce_object, assay = assay_normalized))
  }

  #Make Anndata object
  if (!"logcounts" %in% names(sce_object@assays)) {
    stop("logcounts not found in assays of object!")
  }

  if (add_stats == TRUE) {
    Seu_obj <- as.Seurat(sce_object)
    if (is.null(df_pvals) && value_plot == "expr") {
      df_dge <- .suppressDeprecationWarnings(FindAllMarkers(Seu_obj,
                               features = features,
                               group.by = group_by,
                               min.pct = 0,
                               test.use = test,
                               logfc.threshold = log2fc_cutoff,
                               only.pos = only_pos))

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

    #just for the argument if fc is set to expr
    ident_1 <- NULL
    ident_2 <- NULL
    if (is.null(df_pvals) && value_plot == "fc") {

      df_pvals <- data.frame(matrix(1, nrow = length(features), ncol = length(unique(sce_object[[group_by]]))))
      rownames(df_pvals) <- features
      colnames(df_pvals) <- unique(sce_object[[group_by]])

      if (is.null(group_fc_ident_2)) {
        ident_2 <- grepv(pattern = paste(c("CTRL","Cntrl","WT","healthy"),
                                         collapse ="|"),
                         ignore.case = TRUE,
                         unique(sce_object[[group_fc]]))

        .logger(paste0("group_fc_ident_2 is set to NULL, using: ", ident_2))
      } else{
        ident_2 <- group_fc_ident_2
      }

      if (is.null(group_fc_ident_1)) {
        ident_1 <- grepv(pattern = paste(c("CTRL","Cntrl","WT","healthy"),
                                         collapse ="|"),
                         ignore.case = TRUE,
                         invert = TRUE,
                         unique(sce_object[[group_fc]]))[1]

        .logger(paste0("group_fc_ident_1 is set to NULL, using: ", ident_1))
      } else{
        ident_1 <- group_fc_ident_1
      }

      for (grp in unique(sce_object[[group_by]])) {
        Seu_obj_grp <- subset(Seu_obj, !!sym(group_by) == grp)

        #Check if there are groups with less than 3 cells
        table_cells_sc <- table(Seu_obj_grp@meta.data[[group_fc]])

        count_1_sc <- table_cells_sc[ident_1]
        count_2_sc <- table_cells_sc[ident_2]

        if (is.na(count_1_sc)) {
          count_1_sc <- 0
        }

        if (is.na(count_2_sc)) {
          count_2_sc <- 0
        }

        if (count_1_sc >= 3 && count_2_sc >= 3) {

          #calculating statistics on the group_fc level
          df_dge <- .suppressDeprecationWarnings(FindMarkers(Seu_obj_grp,
                                features = features,
                                group.by = group_fc,
                                min.pct = 0,
                                test.use = test,
                                logfc.threshold = log2fc_cutoff,
                                only.pos = only_pos,
                                ident.1 = ident_1,
                                ident.2 = ident_2))
          df_dge <- rownames_to_column(df_dge, var = "gene")

        } else{
          df_dge <- data.frame()
        }

        # Go through each row in df_dge
        for (i in seq_len(nrow(df_dge))) {
          gene <- df_dge$gene[i]
          cluster <- grp
          pval_adj <- df_dge$p_val_adj[i]
          df_pvals[gene, cluster] <- pval_adj
        }
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
    value_plot = value_plot,
    group_fc_ident_1 = ident_1,
    group_fc_ident_2 = ident_2,
    group_fc = group_fc,
    clip_value = clip_value,
    max_fc = max_fc,
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
            value_plot = args$value_plot,
            group_fc_ident_1 = args$group_fc_ident_1,
            group_fc_ident_2 = args$group_fc_ident_2,
            group_fc = args$group_fc,
            clip_value = args$clip_value,
            max_fc = args$max_fc,
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


