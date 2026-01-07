#' @author Mariano Ruz Jurado & David Rodriguez Morales
#' @title DO Heatmap of the mean expression of genes across a groups
#' @description Wrapper around heatmap_foldchange, which generates a heatmap of
#' showing the foldchange for a set of gene expressions between specified
#' groups. Differential gene expression analysis between the different groups
#' can be performed.
#'
#' @param sce_object A SingleCellExperiment or Seurat object containing
#' expression data and metadata.
#' @param features Character vector of gene names or metadata column names
#' to be visualized.
#' @param reference Reference condition used for fold-change calculation.
#' @param assay_normalized Name of the assay containing normalized expression
#' values (default: "RNA").
#' @param group_by Metadata column defining the primary grouping variable
#' (e.g. clusters).
#' @param condition_key Metadata column defining the condition or comparison
#' variable.
#' @param groups_order Optional character vector specifying the order of groups
#' in \code{group_by}.
#' @param conditions_order Optional character vector specifying the order of
#' conditions.
#' @param layer Optional layer name to extract expression values from.
#'
#' @param figsize Numeric vector of length two specifying figure width and
#' height.
#' @param ax Optional matplotlib axis object (for Python backend usage).
#' @param swap_axes Logical; whether to swap x- and y-axes.
#' @param title Optional title for the heatmap.
#' @param title_fontproperties Named list specifying font properties for the
#' title (e.g. size, weight).
#' @param palette Color palette used for the heatmap.
#' @param palette_conditions Color palette used for condition annotations.
#' @param ticks_fontproperties Named list specifying font properties for axis
#' tick labels.
#' @param xticks_rotation Rotation angle for x-axis tick labels.
#' @param yticks_rotation Rotation angle for y-axis tick labels.
#' @param vmin Minimum value for the color scale.
#' @param vcenter Center value for the color scale.
#' @param vmax Maximum value for the color scale.
#' @param colorbar_legend_title Title for the color bar.
#' @param groups_legend_title Title for the group legend.
#' @param group_legend_ncols Number of columns in the group legend.
#'
#' @param path Optional path to save the output figure.
#' @param filename Name of the output file.
#' @param show Logical; whether to display the plot.
#'
#' @param add_stats Logical; whether to add statistical annotations.
#' @param test Statistical test to use (currently "wilcox").
#' @param correction_method Multiple-testing correction method
#' (currently "bonferroni").
#' @param df_pvals Optional data frame containing precomputed p-values
#' (groups x genes or genes x groups depending on axis orientation).
#' @param stats_x_size Size of the statistical annotation symbol.
#' @param square_x_size Size of the square annotation.
#' @param pval_cutoff P-value significance threshold.
#' @param log2fc_cutoff Minimum absolute log2 fold-change cutoff.
#'
#' @param linewidth Line width of heatmap cell borders.
#' @param color_axis_ratio Relative size of the color bar axis.
#' @return Depending on ``showP``, returns the plot if set to `TRUE` or a
#' dictionary with the axes.
#'
#'
#' @import Seurat
#' @importFrom basilisk basiliskRun
#' @importFrom tidyr pivot_longer
#' @examples
#' sce_data <-
#'   readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' DO.HeatmapFC(
#' sce_object = sce_data,
#' features = c("HES4", "ISG15", "TNFRSF18", "TNFRSF4", "MMP23B"),
#' reference = NULL,
#' assay_normalized = "RNA",
#' group_by = "seurat_clusters",
#' condition_key = "condition",
#' groups_order = c("1", "2", "3", "4", "5", "6", "7", "8"),
#' conditions_order = NULL,
#' layer = NULL,
#'
#' # Figure parameters
#' figsize = c(5, 6),
#' ax = NULL,
#' swap_axes = TRUE,
#' title = NULL,
#' title_fontproperties = list(size = NULL, weight = NULL),
#' palette = "RdBu_r",
#' palette_conditions = "tab10",
#' ticks_fontproperties = list(size = NULL, weight = NULL),
#' xticks_rotation = 45,
#' yticks_rotation = NULL,
#' vmin = NULL,
#' vcenter = NULL,
#' vmax = NULL,
#' colorbar_legend_title = "Log2FC",
#' groups_legend_title = "Comparison",
#' group_legend_ncols = 1,
#'
#' # IO
#' path = NULL,
#' filename = "Heatmap.svg",
#' show = FALSE,
#'
#' #   # Statistics
#' add_stats = TRUE,
#' test = c("wilcox"),
#' correction_method = c("bonferroni"),
#' df_pvals = NULL,
#' stats_x_size = NULL,
#' square_x_size = NULL,
#' pval_cutoff = 0.05,
#' log2fc_cutoff = 0.0,
#'
#' # Fx specific
#' linewidth = 0.1,
#' color_axis_ratio = 0.15
#' )
#'
#' @export
DO.HeatmapFC <- function(
    sce_object,
    features,
    reference = NULL,
    assay_normalized = "RNA",
    group_by = "seurat_clusters",
    condition_key = "condition",
    groups_order = NULL,
    conditions_order = NULL,
    layer = NULL,

    # Figure parameters
    figsize = c(5, 6),
    ax = NULL,
    swap_axes = TRUE,
    title = NULL,
    title_fontproperties = list(size = NULL, weight = NULL),
    palette = "RdBu_r",
    palette_conditions = "tab10",
    ticks_fontproperties = list(size = NULL, weight = NULL),
    xticks_rotation = NULL,
    yticks_rotation = NULL,
    vmin = NULL,
    vcenter = NULL,
    vmax = NULL,
    colorbar_legend_title = "Log2FC",
    groups_legend_title = "Comparison",
    group_legend_ncols = 1,

    # IO
    path = NULL,
    filename = "Heatmap.svg",
    showP = TRUE,

    # Statistics
    add_stats = FALSE,
    test = c("wilcox"),
    correction_method = c("bonferroni"),
    df_pvals = NULL,
    stats_x_size = NULL,
    square_x_size = NULL,
    pval_cutoff = 0.05,
    log2fc_cutoff = 0.0,

    # Fx specific
    linewidth = 0.1,
    color_axis_ratio = 0.15
    ) {
    # support for Seurat objects
    if (methods::is(sce_object, "Seurat")) {
        DefaultAssay(sce_object) <- assay_normalized
        sce_object <- .suppressDeprecationWarnings(
            Seurat::as.SingleCellExperiment(sce_object,
            assay = assay_normalized
            )
        )
    }

    # check for logarithmised counts in object
    if (!"logcounts" %in% names(sce_object@assays)) {
        stop("logcounts not found in assays of object!")
    }

    # for the case of adding statistics
    ident_2 <- NULL
    if (add_stats == TRUE) {
        Seu_obj <- as.Seurat(sce_object)

    if (is.null(reference)) {
        ident_2 <- grepv(
            pattern = paste(c("CTRL", "Cntrl", "WT", "healthy"),
                collapse = "|"),
            ignore.case = TRUE,
        unique(sce_object[[condition_key]])
        )

    .logger(paste0("reference is set to NULL, using: ", ident_2))
    } else {
        ident_2 <- reference
    }

    disease_cond <- grepv(
        ident_2, unique(Seu_obj@meta.data[[condition_key]])
        ,invert = TRUE
    )

    # Data frame collecting p_values for all comparisons
    df_pvals_collector <- data.frame()
    for (ident_1 in disease_cond) {
        if (is.null(df_pvals)) {
            df_pvals <- data.frame(
                t(matrix(1,
                    nrow = length(features),
                    ncol = length(unique(sce_object[[group_by]]))
                    )
                )
            )

        colnames(df_pvals) <- features
        df_pvals[[group_by]] <- unique(sce_object[[group_by]])
        df_pvals[[condition_key]] <- condition_key

        for (grp in unique(sce_object[[group_by]])) {
            Seu_obj_grp <- subset(Seu_obj, !!sym(group_by) == grp)

            # Seu_obj_grp <- subset(Seu_obj, seurat_clusters == grp)
            # Check if there are groups with less than 3 cells
            table_cells_sc <- table(Seu_obj_grp@meta.data[[condition_key]])

            count_1_sc <- table_cells_sc[ident_1]
            count_2_sc <- table_cells_sc[ident_2]

            if (is.na(count_1_sc)) {
                count_1_sc <- 0
            }

            if (is.na(count_2_sc)) {
                count_2_sc <- 0
            }

            if (count_1_sc >= 3 && count_2_sc >= 3) {
                # calculating statistics on the condition_key level
                df_dge <- .suppressDeprecationWarnings(
                    FindMarkers(Seu_obj_grp,
                        features = features,
                        group.by = condition_key,
                        min.pct = 0,
                        test.use = test,
                        logfc.threshold = log2fc_cutoff,
                        only.pos = FALSE,
                        ident.1 = ident_1,
                        ident.2 = ident_2
                    )
                )
                df_dge <- rownames_to_column(df_dge, var = "gene")
            } else {
                df_dge <- data.frame()
            }

            # Go through each row in df_dge
            for (i in seq_len(nrow(df_dge))) {
                gene <- df_dge$gene[i]
                cluster <- grp
                pval_adj <- df_dge$p_val_adj[i]
                df_pvals[cluster, gene] <- pval_adj
            }
        }
    }
    df_pvals_collector <- rbind(df_pvals_collector, df_pvals)
    }
    df_pvals_collector <- df_pvals_collector[,
        c(group_by, condition_key,
        setdiff(colnames(df_pvals_collector), c(group_by, condition_key)))
    ]
    df_long <- df_pvals_collector %>%
        pivot_longer(
            cols = -c(!!sym(group_by), !!sym(condition_key)),
            names_to = "genes",
            values_to = "value"
        )
    }



    # source PATH to python script in install folder
    path_py <- system.file("python", "heatmap_fc.py", package = "DOtools")


    # argument list passed to heatmap inside basilisk
    args <- list(
        sce_object = sce_object,
        features = features,
        reference = ident_2, # set reference to ident_2
        assay_normalized = assay_normalized,
        group_by = group_by,
        condition_key = condition_key,
        groups_order = groups_order,
        conditions_order = conditions_order,
        layer = layer,

        # Figure parameters
        figsize = figsize,
        ax = ax,
        swap_axes = swap_axes,
        title = title,
        title_fontproperties = title_fontproperties,
        palette = palette,
        palette_conditions = palette_conditions,
        ticks_fontproperties = ticks_fontproperties,
        xticks_rotation = xticks_rotation,
        yticks_rotation = yticks_rotation,
        vmin = vmin,
        vcenter = vcenter,
        vmax = vmax,
        colorbar_legend_title = colorbar_legend_title,
        groups_legend_title = groups_legend_title,
        group_legend_ncols = group_legend_ncols,

        # IO
        path = path,
        filename = filename,
        showP = showP,

        # Statistics
        add_stats = add_stats,
        test = test,
        correction_method = correction_method,
        df_pvals = df_long, # added here in long format
        stats_x_size = stats_x_size,
        square_x_size = square_x_size,
        pval_cutoff = pval_cutoff,
        log2fc_cutoff = log2fc_cutoff,

        # Fx specific
        linewidth = linewidth,
        color_axis_ratio = color_axis_ratio
    )


    # basilisk implementation
    results <- basilisk::basiliskRun(env = DOtoolsEnv, fun = function(args) {
        AnnData_counts <- zellkonverter::SCE2AnnData(args$sce_object,
            X_name = "logcounts"
        )

        reticulate::source_python(path_py)

        # Initialize matplot package
        plt <- reticulate::import("matplotlib.pyplot")

        # build dicitonary out of square_x_size if specified
        if (is.numeric(square_x_size) && length(square_x_size) == 2) {
            args$square_x_size <- reticulate::dict(
                weight = square_x_size[1],
                size = square_x_size[2]
            )
        }


    heatmap_foldchange(
        adata = AnnData_counts,
        group_by = args$group_by,
        groups_order = args$groups_order,
        features = args$features,
        condition_key = args$condition_key,
        reference = args$reference,
        conditions_order = args$conditions_order,
        layer = NULL,
        figsize = args$figsize,
        ax = args$ax,
        swap_axes = args$swap_axes,
        title = args$title,
        title_fontproperties = args$title_fontproperties,
        palette = args$palette,
        palette_conditions = args$palette_conditions,
        ticks_fontproperties = args$ticks_fontproperties,
        xticks_rotation = args$xticks_rotation,
        yticks_rotation = args$yticks_rotation,

        vmin = args$vmin,
        vcenter = args$vcenter,
        vmax = args$vmax,
        colorbar_legend_title = args$colorbar_legend_title,
        groups_legend_title = args$groups_legend_title,
        group_legend_ncols = args$group_legend_ncols,
        path = args$path,
        filename = args$filename,
        show = args$showP,
        add_stats = args$add_stats,
        test = args$test,
        correction_method = args$correction_method,
        df_pvals = args$df_pvals,
        stats_x_size = args$stats_x_size,
        square_x_size = args$square_x_size,
        pval_cutoff = args$pval_cutoff,
        log2fc_cutoff = args$log2fc_cutoff,

        linewidth = args$linewidth,
        color_axis_ratio = args$color_axis_ratio,
        )
    }, args = args)
}
