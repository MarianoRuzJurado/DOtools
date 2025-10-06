# Dotplot function for multiple genes, using expression values
#' @author Mariano Ruz Jurado
#' @title DO Dot plot
#' @description This function generates a dot plot for multiple genes, comparing
#' expression levels across one or two specified groups. It supports both
#' individual and pseudobulk expression calculations. Highly variable
#' customization options allow control over dot size, color scaling,
#' annotations, and axis orientation. The function integrates seamlessly with
#' SCE objects for single-cell RNA-seq analysis.
#' @param sce_object The SCE object or Seurat
#' @param group.by.x group name to plot on x-axis
#' @param group.by.y group name to look for in meta data
#' @param group.by.y2 second group name to look for in meta data
#' @param across.group.by.x calculate a pseudobulk expression approach for the
#' x-axis categories
#' @param sort_x Vector sorting the xaxis
#' @param across.group.by.y calculate a pseudobulk expression approach for the
#' y-axis categories
#' @param dot.size Vector of dot size
#' @param plot.margin = plot margins
#' @param midpoint midpoint in color gradient
#' @param Feature Genes or DF of interest, Data frame should have columns with
#' gene and annotation information, e.g. output of FindAllMarkers
#' @param limits_colorscale Set manually colorscale limits
#' @param scale_gene If True calculates the Z-score of the average expression
#' per gene
#' @param hide_zero Removes dots for genes with 0 expression
#' @param annotation_x Adds annotation on top of x axis instead on y axis
#' @param annotation_x_position specifies the position for the annotation
#' @param annotation_x_rev reverses the annotations label order
#' @param point_stroke Defines the thickness of the black stroke on the dots
#' @param coord_flip flips the coordinates of the plot with each other
#' @param returnValue return the dataframe behind the plot
#' @param log1p_nUMI log1p the plotted values, boolean
#' @param stats_x Perform statistical test over categories on the xaxis
#' @param stats_y Perform statistical test over categories on the yaxis
#' @param sig_size Control the size of the significance stars in the plot
#' @param nudge_x Control the position of the star on x axis
#' @param nudge_y Control the position of the star on y axis
#' @param ... Further arguments passed to annoSegment function if
#' annotation_x == TRUE
#'
#' @import ggplot2
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @import reshape2
#' @import ggtext
#' @importFrom SeuratObject as.Seurat
#'
#' @return a ggplot
#'
#' @examples
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' DO.Dotplot(
#'     sce_object = sce_data,
#'     Feature = c("NKG7", "IL6", "MALAT1"),
#'     group.by.x = "condition"
#' )
#'
#' @export
DO.Dotplot <- function(sce_object,
    Feature,
    group.by.x = NULL,
    group.by.y = NULL,
    group.by.y2 = NULL,
    across.group.by.x = FALSE,
    across.group.by.y = FALSE,
    sort_x = NULL,
    dot.size = c(1, 6),
    plot.margin = c(1, 1, 1, 1),
    midpoint = 0.5,
    scale_gene = FALSE,
    returnValue = FALSE,
    log1p_nUMI = TRUE,
    hide_zero = TRUE,
    annotation_x = FALSE,
    annotation_x_position = 0.25,
    annotation_x_rev = FALSE,
    point_stroke = 0.2,
    limits_colorscale = NULL,
    coord_flip = FALSE,
    stats_x = FALSE,
    stats_y = TRUE,
    sig_size = 6,
    nudge_x = 0.3,
    nudge_y = 0.2,
    ...) {
    # support for single cell experiment objects
    if (methods::is(sce_object, "SingleCellExperiment")) {
        sce_object <- .suppressAllWarnings(as.Seurat(sce_object))
    }

    if (!is.vector(Feature) && !is.data.frame(Feature)) {
        stop("Feature is not a vector of strings or a data frame!")
    }

    if (across.group.by.x == TRUE && across.group.by.y == TRUE) {
        stop("Both Pseudobulk options are set to true, please define just one!")
    }

    # type of Feature
    FeatureType <- mode(Feature)

    # check if Feature is a vector and annotation specified as true
    if (is.vector(Feature) && annotation_x == TRUE) {
        stop(
            "Feature is a vector, but annotation_x is set to TRUE. If ",
            "annotation on xaxis is wanted with specific cluster names you ",
            "need to provide a dataframe with a column containing cluster ",
            "names for the genes, like in Seurat::FindAllMarkers!"
        )
    }

    # check the input if it is a data frame
    if (!is.vector(Feature)) {
        orig_DF <- Feature # save original df for annotation purposes
        orig_DF$cluster <- as.vector(orig_DF$cluster)
        cluster_name <- grep("cluster|group|annotation|cell",
            colnames(Feature),
            value = TRUE
        )
        cluster <- unique(Feature[[grep(
            "cluster|group|annotation|cell",
            colnames(Feature)
        )]])
        Feature <- unique(Feature[[grep(
            "gene|feature",
            colnames(Feature)
        )]])

        if (is.null(cluster) || is.null(Feature)) {
            stop(
                "Couldn't derive Cluster and Feature information from the ",
                "provided Dataframe. \n Supported names for cluster: ",
                "cluster|group|annotation|cell\n Supported names for Feature: ",
                "gene|feature. \n Please make sure that your colnames in the ",
                "provided Dataframe are supported."
            )
        }
    }

    # Create Feature expression data frame with grouping information
    geneExp <- expm1(FetchData(
        object = sce_object,
        vars = Feature,
        layer = "data"
    )) #

    # catch wrong handling of arguments
    if (is.null(group.by.x) && !is.null(group.by.y) && is.null(group.by.y2)) {
        stop(
            "If you want to make a Marker Plot with just one provided ",
            "group.by then please use group.by.x!"
        )
    }

    geneExp$xaxis <- sce_object@meta.data[[group.by.x]]

    if (!is.null(group.by.y) && is.null(group.by.y2)) {
        geneExp$id <- paste(sce_object@meta.data[[group.by.y]], sep = "")
    } else if (!is.null(group.by.y) && !is.null(group.by.y2)) {
        geneExp$id <- paste(sce_object@meta.data[[group.by.y]], " (",
            sce_object@meta.data[[group.by.y2]], ")",
            sep = ""
        )
    } else if (is.null(group.by.y) && is.null(group.by.y2)) {
        geneExp$id <- paste(sce_object@meta.data[[group.by.x]], sep = "")
    }

    # Include xaxis in the overall grouping
    data.plot <- lapply(
        X = unique(geneExp$id),
        FUN = function(ident) {
            data.use <- geneExp[geneExp$id == ident, ]

            lapply(
                X = unique(data.use$xaxis),
                FUN = function(x_axis) {
                    data.cell <- data.use[data.use$xaxis == x_axis,
                        seq_len(ncol(geneExp) - 2),
                        drop = FALSE
                    ]
                    avg.exp <- apply(
                        X = data.cell,
                        MARGIN = 2,
                        FUN = function(x) {
                            return(mean(x))
                        }
                    )
                    pct.exp <- apply(
                        X = data.cell,
                        MARGIN = 2,
                        FUN = PercentAbove,
                        threshold = 0
                    )

                    res <- data.frame(
                        id = ident,
                        xaxis = x_axis,
                        avg.exp = avg.exp,
                        pct.exp = pct.exp * 100
                    )
                    res$gene <- rownames(res)
                    return(res)
                }
            ) %>% do.call("rbind", .)
        }
    ) %>%
        do.call("rbind", .) %>%
        data.frame()

    data.plot.res <- data.plot

    # add the cluster information to the plot if annotation_x is set to true
    # and a dataframe was provided with cluster annotation
    # TODO Clean this part a bit up
    if (annotation_x == TRUE && !is.null(cluster)) {
        data.plot.res <- purrr::map_df(
            seq_len(nrow(data.plot.res)), function(x) {
                tmp <- data.plot.res[x, ]
                tmp$celltype <- orig_DF[
                    which(
                        orig_DF[[grep("gene|feature", colnames(orig_DF))]] ==
                            tmp[[grep("gene|feature", colnames(tmp))]]
                    ),
                    cluster_name
                ][[1]]
                return(tmp)
            }
        )
    }
    # sort x by provided x-axis
    if (is.null(sort_x) && !is.null(group.by.y)) {
        data.plot.res$xaxis <- factor(data.plot.res$xaxis,
            levels = sort(unique(data.plot.res$xaxis))
        )
    } else if (!is.null(sort_x) && !is.null(group.by.y)) {
        data.plot.res$xaxis <- factor(data.plot.res$xaxis,
            levels = sort_x
        )
    }

    # sort x by provided x-axis for just one group.by.x
    if (is.null(sort_x) && is.null(group.by.y)) {
        data.plot.res$gene <- factor(data.plot.res$gene,
            levels = sort(unique(data.plot.res$gene))
        )
    } else if (!is.null(sort_x) && is.null(group.by.y)) {
        data.plot.res$gene <- factor(data.plot.res$gene,
            levels = sort_x
        )
    }

    # create grouping column for multiple grouping variables on the y-axis
    if (!is.null(group.by.y2)) {
        data.plot.res$group <- vapply(
            strsplit(as.character(data.plot.res$id),
                split = "\\(|\\)"
            ), "[",
            character(1), 2
        )
    }

    if (hide_zero == TRUE) {
        # so fraction 0 is not displayed in plot
        data.plot.res$pct.exp <- ifelse(data.plot.res$pct.exp == 0,
            NA,
            data.plot.res$pct.exp
        )
        # remove empty lines
        data.plot.res <- data.plot.res[stats::complete.cases(
            data.plot.res$pct.exp
            ), ]
    }

    # create bulk expression for group.by.x
    if (across.group.by.x == TRUE) {
        bulk_tmp <- data.plot.res %>%
            dplyr::group_by(id, gene) %>%
            summarise(
                avg.exp = mean(avg.exp),
                pct.exp = mean(pct.exp)
            )
        bulk_tmp$xaxis <- "Pseudobulk"
        data.plot.res <- dplyr::bind_rows(data.plot.res, bulk_tmp)
        data.plot.res$xaxis <- factor(data.plot.res$xaxis,
            levels = c(
                "Pseudobulk",
                setdiff(sort(
                    unique(data.plot.res$xaxis)
                ), "Pseudobulk")
            )
        )
    }

    # create bulk expression for group.by.y, will divide by group.by.y2
    if (across.group.by.y == TRUE) {
        if (is.null(group.by.y2)) {
            bulk_tmp <- data.plot.res %>%
                dplyr::group_by(xaxis, gene) %>%
                summarise(avg.exp = mean(avg.exp), pct.exp = mean(pct.exp))
            bulk_tmp$id <- "Pseudobulk"
            pseudo_levels <- unique(bulk_tmp$id)
            data.plot.res <- dplyr::bind_rows(data.plot.res, bulk_tmp)
            data.plot.res$id <-
                factor(data.plot.res$id,
                    levels = c(
                        pseudo_levels,
                        setdiff(
                            sort(unique(data.plot.res$id)),
                            pseudo_levels
                        )
                    )
                )
        } else {
            bulk_tmp <- data.plot.res %>%
                dplyr::group_by(xaxis, gene, group) %>%
                summarise(avg.exp = mean(avg.exp), pct.exp = mean(pct.exp))
            bulk_tmp$id <- paste0("Pseudobulk (", bulk_tmp$group, ")")
            data.plot.res <- dplyr::bind_rows(data.plot.res, bulk_tmp)
            pseudo_levels <- unique(bulk_tmp$id)
            data.plot.res$id <-
                factor(data.plot.res$id,
                    levels = c(
                        pseudo_levels,
                        setdiff(
                            sort(unique(data.plot.res$id)),
                            pseudo_levels
                        )
                    )
                )
        }
    }

    # get the scale pvalue for plotting
    if (log1p_nUMI == TRUE) {
        data.plot.res$avg.exp.plot <- log1p(data.plot.res$avg.exp)
    } else {
        data.plot.res$avg.exp.plot <- data.plot.res$avg.exp
    }

    # define how expression values are transformed
    if (scale_gene == TRUE) {
        data.plot.res <- data.plot.res %>%
            dplyr::group_by(gene) %>%
            dplyr::mutate(
                z_avg_exp = (avg.exp - mean(avg.exp, na.rm = TRUE)) /
                    stats::sd(avg.exp, na.rm = TRUE)
            ) %>%
            ungroup()
        exp.title <- "Z-score expression \n per gene"
        fill.values <- data.plot.res$z_avg_exp
        ###
    } else if (log1p_nUMI == TRUE) {
        exp.title <- "Mean log(nUMI) \n in group"
        fill.values <- data.plot.res$avg.exp.plot
    } else {
        exp.title <- "Mean nUMI \n in group"
        fill.values <- data.plot.res$avg.exp.plot
    }

    # Define which columns to take for dotplot,
    # it should be able to correctly capture one group.by.x
    if (identical(
        as.vector(data.plot.res$id), as.vector(data.plot.res$xaxis)
    ) &&
        FeatureType == "list") {
        # go over input type
        # get rid of previous factoring to set new one,
        # first alphabetical order on y
        data.plot.res$xaxis <- as.vector(data.plot.res$xaxis)
        data.plot.res$id <-
            factor(data.plot.res$id,
                levels = sort(unique(data.plot.res$id))
            )
        data.plot.res$gene <-
            factor(data.plot.res$gene,
                levels = orig_DF[
                    order(orig_DF$cluster, decreasing = FALSE),
                ]$gene
            )

        if (annotation_x_rev == TRUE) {
            data.plot.res$id <-
                factor(data.plot.res$id,
                    levels = rev(sort(unique(data.plot.res$id)))
                )
        }

        aes_var <- c("gene", "id")
        # there is a second case here for providing just a gene list
        # which need to be adressed with the same aes_var
    } else if (identical(
        as.vector(data.plot.res$id),
        as.vector(data.plot.res$xaxis)
    )) {
        aes_var <- c("gene", "id")
    } else { # all other cases where group.by.y is specified
        data.plot.res$id <- factor(data.plot.res$id,
            levels = rev(sort(unique(data.plot.res$id)))
        )
        aes_var <- c("xaxis", "id")
    }


    # Implemenent statistics in dotplot for x-axis
    if (stats_x & !is.null(group.by.x) & !is.null(group.by.y)) {
        data.plot.res$p_adj <- NA
        for (names_y in as.vector(unique(data.plot.res[[aes_var[2]]]))) {
            # in the case of pseudobulk specified
            if (names_y == "Pseudobulk") {
                stats_test <- .suppressAllWarnings(
                    Seurat::FindAllMarkers(
                        sce_object,
                        features = unique(c(data.plot.res$gene)),
                        min.pct = 0,
                        logfc.threshold =
                            0,
                        group.by = group.by.x,
                        only.pos = TRUE
                    )
                )

                if (!nrow(stats_test) == 0) {
                    stats_test_ren <- stats_test %>%
                        rename(xaxis = cluster)
                } else {
                    stats_test_ren <- data.frame(
                        p_val = character(),
                        avg_log2FC = character(),
                        pct.1 = character(),
                        pct.2 = character(),
                        p_val_adj = character(),
                        xaxis = character(),
                        gene = character()
                    )
                }

                # add p-value if test fails add p-value 1
                data.plot.res <- data.plot.res %>%
                    left_join(
                        stats_test_ren %>% select(gene, xaxis, p_val_adj),
                        by = c("gene", "xaxis")
                    ) %>%
                    mutate(
                        p_adj = case_when(
                            id == names_y & !is.na(p_val_adj) ~ p_val_adj,
                            id == names_y & is.na(p_val_adj) ~ 1,
                            TRUE ~ p_adj
                        )
                    ) %>%
                    select(-p_val_adj)
            } else {
                sce_object_sub <- subset(sce_object,
                    subset = !!sym(group.by.y) == names_y
                )

                # FindMarkers from Seurat in each cluster of group.by.x
                stats_test <- .suppressAllWarnings(
                    Seurat::FindAllMarkers(
                        sce_object_sub,
                        features = unique(c(data.plot.res$gene)),
                        min.pct = 0,
                        logfc.threshold =
                            0,
                        group.by = group.by.x,
                        only.pos = TRUE
                    )
                )

                if (!nrow(stats_test) == 0) {
                    stats_test_ren <- stats_test %>%
                        rename(xaxis = cluster)
                } else {
                    stats_test_ren <- data.frame(
                        p_val = character(),
                        avg_log2FC = character(),
                        pct.1 = character(),
                        pct.2 = character(),
                        p_val_adj = character(),
                        xaxis = character(),
                        gene = character()
                    )
                }

                # add p-value if test fails add p-value 1
                data.plot.res <- data.plot.res %>%
                    left_join(
                        stats_test_ren %>% select(gene, xaxis, p_val_adj),
                        by = c("gene", "xaxis")
                    ) %>%
                    mutate(
                        p_adj = case_when(
                            id == names_y & !is.na(p_val_adj) ~ p_val_adj,
                            id == names_y & is.na(p_val_adj) ~ 1,
                            TRUE ~ p_adj
                        )
                    ) %>%
                    select(-p_val_adj)
            }
        }

        data.plot.res <- data.plot.res %>%
            mutate(
                stars = case_when(
                    p_adj < 0.05 ~ "*",
                    TRUE ~ "",
                )
            )
        data.plot.res$significance <- ifelse(data.plot.res$p_adj < 0.05,
            "p < 0.05",
            "ns"
        )
    }


    # Implemenent statistics in dotplot for y-axis
    if (stats_y & !is.null(group.by.x) & !is.null(group.by.y)) {
        data.plot.res$p_adj <- NA
        for (names_x in as.vector(unique(data.plot.res[[aes_var[1]]]))) {
            # in the case of pseudobulk specified
            if (names_x == "Pseudobulk") {
                stats_test <- .suppressAllWarnings(
                    Seurat::FindAllMarkers(
                        sce_object,
                        features = unique(c(data.plot.res$gene)),
                        min.pct = 0,
                        logfc.threshold =
                            0,
                        group.by = group.by.y,
                        only.pos = TRUE
                    )
                )

                if (!nrow(stats_test) == 0) {
                    stats_test_ren <- stats_test %>%
                        rename(id = cluster)
                } else {
                    stats_test_ren <- data.frame(
                        p_val = character(),
                        avg_log2FC = character(),
                        pct.1 = character(),
                        pct.2 = character(),
                        p_val_adj = character(),
                        id = character(),
                        gene = character()
                    )
                }

                # add p-value if test fails add p-value 1
                data.plot.res <- data.plot.res %>%
                    left_join(stats_test_ren %>% select(gene, id, p_val_adj),
                        by = c("gene", "id")
                    ) %>%
                    mutate(
                        p_adj = case_when(
                            xaxis == names_x &
                                !is.na(
                                    as.numeric(p_val_adj)
                                ) ~ as.numeric(p_val_adj),
                            xaxis == names_x & is.na(as.numeric(p_val_adj)) ~ 1,
                            TRUE ~ as.numeric(p_adj)
                        )
                    ) %>%
                    select(-p_val_adj)
            } else {
                sce_object_sub <-
                    subset(sce_object, subset = !!sym(group.by.x) == names_x)

                # FindMarkers from Seurat in each cluster of group.by.x
                stats_test <- .suppressAllWarnings(
                    Seurat::FindAllMarkers(
                        sce_object_sub,
                        features = unique(c(data.plot.res$gene)),
                        min.pct = 0,
                        logfc.threshold =
                            0,
                        group.by = group.by.y,
                        only.pos = TRUE
                    )
                )
                if (!nrow(stats_test) == 0) {
                    stats_test_ren <- stats_test %>%
                        rename(id = cluster)
                } else {
                    stats_test_ren <- data.frame(
                        p_val = character(),
                        avg_log2FC = character(),
                        pct.1 = character(),
                        pct.2 = character(),
                        p_val_adj = character(),
                        id = character(),
                        gene = character()
                    )
                }

                # add p-value if test fails add p-value 1
                data.plot.res <- data.plot.res %>%
                    left_join(stats_test_ren %>% select(gene, id, p_val_adj),
                        by = c("gene", "id")
                    ) %>%
                    mutate(
                        p_adj = case_when(
                            xaxis == names_x &
                                !is.na(
                                    as.numeric(p_val_adj)
                                ) ~ as.numeric(p_val_adj),
                            xaxis == names_x & is.na(as.numeric(p_val_adj)) ~ 1,
                            TRUE ~ as.numeric(p_adj)
                        )
                    ) %>%
                    select(-p_val_adj)
            }
        }
        data.plot.res <- data.plot.res %>%
            mutate(
                stars = case_when(
                    p_adj < 0.05 ~ "*",
                    TRUE ~ "",
                )
            )
        data.plot.res$significance <- ifelse(data.plot.res$p_adj < 0.05,
            "p < 0.05",
            "ns"
        )
    }

    pmain <- ggplot2::ggplot(data.plot.res, ggplot2::aes(
        x = !!sym(aes_var[1]),
        y = !!sym(aes_var[2])
    )) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::coord_fixed(clip = "off") +
        ggplot2::theme(
            plot.margin = ggplot2::margin(
                t = plot.margin[1],
                r = plot.margin[2],
                b = plot.margin[3],
                l = plot.margin[4],
                unit = "cm"
            ),
            axis.text = ggplot2::element_text(color = "black"),
            legend.direction = "horizontal",
            axis.text.x = element_text(
                color = "black",
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = 14,
                family = "Helvetica"
            ),
            axis.text.y = element_text(
                color = "black",
                size = 14,
                family = "Helvetica"
            ),
            axis.title.x = element_text(
                color = "black",
                size = 14,
                family = "Helvetica"
            ),
            axis.title = element_text(
                size = 14,
                color = "black",
                family = "Helvetica"
            ),
            plot.title = element_text(
                size = 14,
                hjust = 0.5,
                face = "bold",
                family = "Helvetica"
            ),
            plot.subtitle = element_text(
                size = 14,
                hjust = 0,
                family = "Helvetica"
            ),
            axis.line = element_line(color = "black"),
            strip.text.x = element_text(
                size = 14,
                color = "black",
                family = "Helvetica",
                face = "bold"
            ),
            legend.text = element_text(
                size = 10,
                color = "black",
                family = "Helvetica"
            ),
            legend.title = element_text(
                size = 10,
                color = "black",
                family = "Helvetica",
                hjust = 0
            ),
            legend.position = "right",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
        )

    guides.layer <- ggplot2::guides(
        fill = ggplot2::guide_colorbar(
            title = exp.title,
            title.position = "top",
            title.hjust = 0.5,
            barwidth = unit(3.8, "cm"),
            # changes the width of the color legend
            barheight = unit(0.5, "cm"),
            frame.colour = "black",
            frame.linewidth = 0.3,
            ticks.colour = "black",
            order = 2
        ),
        size = ggplot2::guide_legend(
            title = "Fraction of cells \n in group (%)",
            title.position = "top",
            title.hjust = 0.5,
            label.position = "bottom",
            override.aes = list(color = "black", fill = "grey50"),
            keywidth = ggplot2::unit(0.5, "cm"),
            # changes the width of the precentage dots in legend
            order = 1
        )
    )

    # TODO change the scalefillgradient to +n in the else part
    dot.col <- c("#fff5f0", "#990000")
    gradient_colors <- c("#fff5f0", "#fcbba1", "#fc9272", "#fb6a4a", "#990000")
    # "#FFFFFF","#08519C","#BDD7E7" ,"#6BAED6", "#3182BD",
    if (length(dot.col) == 2) {
        breaks <- scales::breaks_extended(n = 5)(range(fill.values))

        if (is.null(limits_colorscale)) {
            limits_colorscale <- c(
                min(range(fill.values)) * .99,
                max(range(fill.values)) * 1.01
            )
        }

        if (max(breaks) > max(limits_colorscale)) {
            limits_colorscale[length(limits_colorscale)] <-
                breaks[length(breaks)]
        }

        pmain <- pmain + ggplot2::scale_fill_gradientn(
            colours = gradient_colors,
            breaks = breaks,
            limits = limits_colorscale
        )
    } else {
        pmain <- pmain + ggplot2::scale_fill_gradient2(
            low = dot.col[1],
            mid = dot.col[2],
            high = dot.col[3],
            midpoint = midpoint,
            name = "Gradient"
        )
    }
    if (across.group.by.x == TRUE) {
        pmain <- pmain +
            ggplot2::geom_point(
                ggplot2::aes(
                    fill = fill.values,
                    size = pct.exp
                ),
                shape = 21,
                stroke = point_stroke
            ) +
            guides.layer +
            facet_grid(
                cols = vars(gene %>% factor(levels = Feature)),
                scales = "fixed"
            ) +
            ggplot2::scale_size(range = c(dot.size[1], dot.size[2])) +
            ggplot2::scale_size_continuous(
                breaks = pretty(round(as.vector(
                    stats::quantile(data.plot.res$pct.exp)
                )), n = 10)[seq(1, 10, by = 2)],
                limits = c(
                    min(data.plot.res$pct.exp) * 1.05,
                    max(data.plot.res$pct.exp) * 1.05
                )
            ) +
            theme(
                panel.spacing = unit(0, "lines"),
                axis.text.x = ggtext::element_markdown(
                    color = "black",
                    angle = 90,
                    hjust = 1,
                    vjust = 0.5,
                    size = 14,
                    family = "Helvetica"
                )
            ) +
            scale_x_discrete(
                labels = function(labels) {
                    labels <- ifelse(labels == "Pseudobulk",
                        paste0("<b>", labels, "</b>"),
                        labels
                    )
                    return(labels)
                }
            )
    } else if (across.group.by.y == TRUE) {
        pmain <- pmain + ggplot2::geom_point(
            ggplot2::aes(
                fill = fill.values,
                size = pct.exp
            ),
            shape = 21,
            stroke = point_stroke
        ) +
            guides.layer + facet_grid(
                cols = vars(gene %>% factor(levels = Feature)),
                scales = "fixed"
            ) +
            ggplot2::scale_size(range = c(dot.size[1], dot.size[2])) +
            ggplot2::scale_size_continuous(
                breaks = pretty(round(as.vector(
                    stats::quantile(data.plot.res$pct.exp)
                )), n = 10)[seq(1, 10, by = 2)],
                limits = c(
                    min(data.plot.res$pct.exp) * 1.05,
                    max(data.plot.res$pct.exp) * 1.05
                )
            ) +
            theme(
                panel.spacing = unit(0, "lines"),
                axis.text.y = ggtext::element_markdown(
                    color = "black",
                    angle = 0,
                    hjust = 1,
                    vjust = 0.5,
                    size = 14,
                    family = "Helvetica"
                )
            ) +
            scale_y_discrete(
                labels = function(labels) {
                    labels <- ifelse(labels %in% pseudo_levels,
                        paste0("<b>", labels, "</b>"),
                        labels
                    )
                    return(labels)
                }
            )
    } else if (identical(
        as.vector(data.plot.res$id),
        as.vector(data.plot.res$xaxis)
    )) {
        pmain <- pmain +
            ggplot2::geom_point(
                ggplot2::aes(
                    fill = fill.values,
                    size = pct.exp
                ),
                shape = 21,
                stroke = point_stroke
            ) +
            guides.layer +
            # facet_wrap(~facet_group, scales="free_x")+
            ggplot2::scale_size(range = c(dot.size[1], dot.size[2])) +
            ggplot2::scale_size_continuous(
                breaks = pretty(round(as.vector(
                    stats::quantile(data.plot.res$pct.exp)
                )), n = 10)[seq(1, 10, by = 2)],
                limits = c(min(pretty(
                    round(as.vector(stats::quantile(
                        data.plot.res$pct.exp
                    ))),
                    n = 10
                )[seq(1, 10, by = 2)]) * .95, max(data.plot.res$pct.exp) * 1.05)
            ) +
            theme(panel.spacing = unit(0, "lines"))

        if (annotation_x == TRUE) {
            plot_max_y <- ggplot_build(pmain)
            plot_max_y <- plot_max_y$layout$panel_params[[1]]$y.range[2] +
                annotation_x_position

            pmain <- .annoSegment(
                object = pmain,
                annoPos = "top",
                aesGroup = TRUE,
                aesGroName = "celltype",
                fontface = "bold",
                fontfamily = "Helvetica",
                pCol = rep("black", length(cluster)),
                textCol = rep("black", length(cluster)),
                addBranch = TRUE,
                branDirection = -1,
                addText = TRUE,
                yPosition = plot_max_y,
                # textSize = 14,
                # hjust = 0.5,
                # vjust = 0,
                # textRot = 0,
                # segWidth = 0.3,
                # lwd = 3
                ...
            )
        }

        if (coord_flip == TRUE) {
            pmain <- pmain + ggplot2::coord_flip()
        }

        if (coord_flip == TRUE && annotation_x == TRUE) {
            warning(
                "Annotation_x and coord_flip set on TRUE, might result in ",
                "unwanted behaviour!"
            )
        }
    } else {
        pmain <- pmain +
            ggplot2::geom_point(
                ggplot2::aes(
                    fill = fill.values,
                    size = pct.exp
                ),
                shape = 21,
                stroke = point_stroke
            ) +
            guides.layer +
            facet_grid(
                cols = vars(gene %>% factor(levels = Feature)),
                scales = "fixed"
            ) +
            ggplot2::scale_size(range = c(dot.size[1], dot.size[2])) +
            ggplot2::scale_size_continuous(
                breaks = pretty(round(as.vector(
                    stats::quantile(data.plot.res$pct.exp)
                )), n = 10)[seq(1, 10, by = 2)],
                limits = c(min(pretty(
                    round(as.vector(stats::quantile(
                        data.plot.res$pct.exp
                    ))),
                    n = 10
                )[seq(1, 10, by = 2)]) * .95, max(data.plot.res$pct.exp) * 1.05)
            ) +
            theme(panel.spacing = unit(0, "lines"))
    }

    if (!is.null(data.plot.res$stars)) {
        pmain <- pmain +
            geom_text(
                aes(label = stars),
                color = "black",
                size = sig_size,
                nudge_x = nudge_x,
                nudge_y = nudge_y,
                fontface = "bold"
            ) +
            geom_point(
                data = subset(data.plot.res, significance == "p < 0.05"),
                aes(shape = significance),
                alpha = 0
            ) +
            scale_shape_manual(
                values = c("p < 0.05" = 8),
                name = "Significance"
            ) +
            guides(shape = guide_legend(override.aes = list(
                alpha = 1,
                size = 4
            )))
    }

    if (returnValue == TRUE) {
        return(data.plot.res)
    }
    return(pmain)
}


# AnnoSegment function modifications
#' @author Mariano Ruz Jurado (edited from: Jun Zhang)
#' @title Annotation modifier for plots
#'
#' @description Used for segment the plot for further annotations
#'
#' @param object ggplot list. Default(NULL).
#' @param relSideDist The relative distance ratio to the y axis range.
#' Default(0.1).
#' @param aesGroup Whether use your group column to add rect annotation.
#' Default("FALSE").
#' @param aesGroName The mapping column name. Default(NULL).
#' @param annoPos The position for the annotation to be added. Default("top").
#' @param xPosition The x axis coordinate for the segment. Default(NULL).
#' @param yPosition The y axis coordinate for the segment. Default(NULL).
#' @param pCol The segment colors. Default(NULL).
#' @param segWidth The relative segment width. Default(1).
#' @param lty The segment line type. Default(NULL).
#' @param lwd The segment line width. Default(NULL).
#' @param alpha The segment color alpha. Default(NULL).
#' @param lineend The segment line end. Default("square").
#' @param annoManual Whether annotate by yourself by supplying with x and y
#' coordinates. Default(FALSE).
#' @param mArrow Whether add segment arrow. Default(FALSE).
#' @param addBranch Whether add segment branch. Default(FALSE).
#' @param bArrow Whether add branch arrow. Default(FALSE).
#' @param branDirection The branch direction. Default(1).
#' @param branRelSegLen The branch relative length to the segment. Default(0.3).
#' @param addText Whether add text label on segment. Default(FALSE).
#' @param textCol The text colors. Default(NULL).
#' @param textSize The text size. Default(NULL).
#' @param fontfamily The text fontfamily. Default(NULL).
#' @param fontface The text fontface. Default(NULL).
#' @param textLabel The text textLabel. Default(NULL).
#' @param textRot The text angle. Default(NULL).
#' @param textHVjust The text distance from the segment. Default(0.2).
#' @param hjust The text hjust. Default(NULL).
#' @param vjust The text vjust. Default(NULL).
#' @param myFacetGrou Your facet group name to be added with annotation when
#' object is a faceted object. Default(NULL).
#' @param aes_x = NULL You should supply the plot X mapping name when annotate
#' a facetd plot. Default(NULL).
#' @param aes_y = NULL You should supply the plot Y mapping name when annotate
#' a facetd plot. Default(NULL).
#'
#' @return ggplot
#'
#' @keywords internal
.annoSegment <- function(object = NULL,
    relSideDist = 0.1,
    aesGroup = FALSE,
    aesGroName = NULL,
    annoPos = "top",
    xPosition = NULL,
    yPosition = NULL,
    pCol = NULL,
    segWidth = 1,
    lty = NULL,
    lwd = 10,
    alpha = NULL,
    lineend = "square",
    annoManual = FALSE,
    mArrow = NULL,
    addBranch = FALSE,
    bArrow = NULL,
    branDirection = 1,
    branRelSegLen = 0.3,
    addText = FALSE,
    textCol = NULL,
    textSize = NULL,
    fontfamily = NULL,
    fontface = NULL,
    textLabel = NULL,
    textRot = 0,
    textHVjust = 0.2,
    hjust = NULL,
    vjust = NULL,
    myFacetGrou = NULL,
    aes_x = NULL,
    aes_y = NULL) {
    facetName <- names(object$facet$params$facets)
    if (is.null(myFacetGrou) & !is.null(facetName)) {
        myFacetGrou <- unique(data[, facetName])[1]
    } else if (!is.null(myFacetGrou) & !is.null(facetName)) {
        myFacetGrou <- myFacetGrou
    } else {
    }
    data <- object$data
    if (is.null(facetName)) {
        aes_x <- ggiraphExtra::getMapping(object$mapping, "x")
        aes_y <- ggiraphExtra::getMapping(object$mapping, "y")
    } else {
        aes_x <- aes_x
        aes_y <- aes_y
    }
    data_x <- data[, c(aes_x)]
    data_y <- data[, c(aes_y)]
    if (annoManual == FALSE) {
        if (annoPos %in% c("top", "botomn")) {
            if (aesGroup == FALSE) {
                nPoints <- length(xPosition)
                xPos <- xPosition
                xmin <- xPos - segWidth / 2
                xmax <- xPos + segWidth / 2
            } else {
                groupInfo <- data %>%
                    dplyr::select(.data[[aes_x]], .data[[aesGroName]]) %>%
                    unique() %>%
                    dplyr::select(.data[[aesGroName]]) %>%
                    table() %>%
                    data.frame()

                # safer indexing
                start <- c(
                    1, groupInfo$Freq[seq_len(length(groupInfo$Freq) - 1)]
                ) %>%
                    cumsum()
                end <- cumsum(groupInfo$Freq)
                xmin <- start - segWidth / 2
                xmax <- end + segWidth / 2
                nPoints <- length(start)
            }
            if (is.null(yPosition)) {
                if (is.numeric(data_y)) {
                    if (annoPos == "top") {
                        ymax <- max(data_y) + relSideDist * max(data_y)
                        ymin <- ymax
                    } else {
                        ymin <- min(data_y) - relSideDist * max(data_y)
                        ymax <- ymin
                    }
                } else {
                    if (annoPos == "top") {
                        ymax <- length(unique(data_y)) + relSideDist *
                            length(unique(data_y))
                        ymin <- ymax
                    } else {
                        ymin <- -relSideDist * length(unique(data_y))
                        ymax <- ymin
                    }
                }
            } else {
                ymax <- yPosition[1]
                ymin <- yPosition[1]
            }
        } else if (annoPos %in% c("left", "right")) {
            if (aesGroup == FALSE) {
                nPoints <- length(yPosition)
                yPos <- yPosition
                ymin <- yPos - segWidth / 2
                ymax <- yPos + segWidth / 2
            } else {
                groupInfo <- data %>%
                    dplyr::select(.data[[aes_y]], .data[[aesGroName]]) %>%
                    unique() %>%
                    dplyr::select(.data[[aesGroName]]) %>%
                    table() %>%
                    data.frame()

                # safer indexing with seq_len()
                start <- c(
                    1, groupInfo$Freq[seq_len(length(groupInfo$Freq) - 1)]
                ) %>%
                    cumsum()
                end <- cumsum(groupInfo$Freq)
                ymin <- start - segWidth / 2
                ymax <- end + segWidth / 2
                nPoints <- length(start)
            }
            if (is.null(xPosition)) {
                if (is.numeric(data_x)) {
                    if (annoPos == "left") {
                        xmin <- min(data_x) - relSideDist * max(data_x)
                        xmax <- xmin
                    } else {
                        xmax <- max(data_x) + relSideDist * max(data_x)
                        xmin <- xmax
                    }
                } else {
                    if (annoPos == "left") {
                        xmin <- -relSideDist * length(unique(data_x))
                        xmax <- xmin
                    } else {
                        xmax <- length(unique(data_x)) + relSideDist *
                            length(unique(data_x))
                        xmin <- xmax
                    }
                }
            } else {
                xmin <- xPosition[1]
                xmax <- xPosition[1]
            }
        }
    } else {
        if (annoPos %in% c("top", "botomn")) {
            xmin <- xPosition[[1]] - segWidth / 2
            xmax <- xPosition[[2]] + segWidth / 2
            ymax <- yPosition[[1]]
            ymin <- yPosition[[1]]
        } else {
            xmin <- xPosition[[1]]
            xmax <- xPosition[[1]]
            ymin <- yPosition[[1]] - segWidth / 2
            ymax <- yPosition[[2]] + segWidth / 2
        }
        nPoints <- max(length(xmin), length(ymin))
    }
    annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf,
                                    ymin = -Inf, ymax = Inf, data) {
        ggplot2::layer(
            data = data,
            stat = StatIdentity,
            position = PositionIdentity,
            geom = ggplot2::GeomCustomAnn,
            inherit.aes = TRUE,
            params = list(
                grob = grob,
                xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax
            )
        )
    }
    if (is.null(pCol)) {
        pCol <- useMyCol("stallion", n = nPoints)
    } else {
        pCol <- pCol
    }
    if (is.null(facetName)) {
        if (annoPos %in% c("top", "botomn")) {
            for (i in seq_len(nPoints)) {
                object <- object + ggplot2::annotation_custom(
                    grob = grid::segmentsGrob(
                        gp = grid::gpar(
                            col = pCol[i],
                            fill = pCol[i],
                            lty = lty,
                            lwd = lwd,
                            lineend = lineend,
                            alpha = alpha
                        ),
                        arrow = mArrow
                    ),
                    xmin = ggplot2::unit(xmin[i], "native"),
                    xmax = ggplot2::unit(xmax[i], "native"),
                    ymin = ggplot2::unit(ymin, "native"),
                    ymax = ggplot2::unit(ymax, "native")
                )
            }
        } else if (annoPos %in% c("left", "right")) {
            for (i in seq_len(nPoints)) {
                object <- object + ggplot2::annotation_custom(
                    grob = grid::segmentsGrob(
                        gp = grid::gpar(
                            col = pCol[i],
                            fill = pCol[i],
                            lty = lty,
                            lwd = lwd,
                            lineend = lineend,
                            alpha = alpha
                        ),
                        arrow = mArrow
                    ),
                    xmin = ggplot2::unit(xmin, "native"),
                    xmax = ggplot2::unit(xmax, "native"),
                    ymin = ggplot2::unit(ymin[i], "native"),
                    ymax = ggplot2::unit(ymax[i], "native")
                )
            }
        } else {
        }
    } else {
        facet_data <- data.frame(myFacetGrou)
        colnames(facet_data) <- facetName
        if (annoPos %in% c("top", "botomn")) {
            for (i in seq_len(nPoints)) {
                object <- object + annotation_custom2(
                    grob = grid::segmentsGrob(
                        gp = grid::gpar(
                            col = pCol[i],
                            fill = pCol[i],
                            lty = lty,
                            lwd = lwd,
                            lineend = lineend,
                            alpha = alpha
                        ),
                        arrow = mArrow
                    ),
                    data = facet_data,
                    xmin = xmin[i],
                    xmax = xmax[i],
                    ymin = ymin,
                    ymax = ymax
                )
            }
        } else if (annoPos %in% c("left", "right")) {
            for (i in seq_len(nPoints)) {
                object <- object + annotation_custom2(
                    grob = grid::segmentsGrob(
                        gp = grid::gpar(
                            col = pCol[i],
                            fill = pCol[i],
                            lty = lty,
                            lwd = lwd,
                            lineend = lineend,
                            alpha = alpha
                        ),
                        arrow = mArrow
                    ),
                    data = facet_data,
                    xmin = xmin,
                    xmax = xmax,
                    ymin = ymin[i],
                    ymax = ymax[i]
                )
            }
        } else {
        }
    }
    if (addBranch == TRUE) {
        if (annoPos %in% c("top", "botomn")) {
            brXmin <- c(xmin, xmax)
            brXmax <- c(xmin, xmax)
            brYmin <- ymax + branRelSegLen * segWidth * branDirection
            brYmax <- ymax
        } else {
            brXmin <- xmax
            brXmax <- xmax + branRelSegLen * segWidth * branDirection
            brYmin <- c(ymin, ymax)
            brYmax <- c(ymin, ymax)
        }
        pCol2 <- rep(pCol, 2)
    }
    if (is.null(facetName)) {
        if (addBranch == TRUE & annoPos %in% c("top", "botomn")) {
            for (i in seq_len(2 * nPoints)) {
                object <- object + ggplot2::annotation_custom(
                    grob = grid::segmentsGrob(
                        gp = grid::gpar(
                            col = pCol2[i],
                            fill = pCol2[i],
                            lty = lty,
                            lwd = lwd,
                            lineend = lineend,
                            alpha = alpha
                        ),
                        arrow = bArrow
                    ),
                    xmin = ggplot2::unit(brXmin[i], "native"),
                    xmax = ggplot2::unit(brXmax[i], "native"),
                    ymin = ggplot2::unit(brYmin, "native"),
                    ymax = ggplot2::unit(brYmax, "native")
                )
            }
        } else if (addBranch == TRUE & annoPos %in% c("left", "right")) {
            for (i in seq_len(2 * nPoints)) {
                object <- object + ggplot2::annotation_custom(
                    grob = grid::segmentsGrob(
                        gp = grid::gpar(
                            col = pCol2[i],
                            fill = pCol2[i],
                            lty = lty,
                            lwd = lwd,
                            lineend = lineend,
                            alpha = alpha
                        ),
                        arrow = bArrow
                    ),
                    xmin = ggplot2::unit(brXmin, "native"),
                    xmax = ggplot2::unit(brXmax, "native"),
                    ymin = ggplot2::unit(brYmin[i], "native"),
                    ymax = ggplot2::unit(brYmax[i], "native")
                )
            }
        } else {
        }
    } else {
        if (addBranch == TRUE & annoPos %in% c("top", "botomn")) {
            for (i in seq_len(2 * nPoints)) {
                object <- object + annotation_custom2(
                    grob = grid::segmentsGrob(
                        gp = grid::gpar(
                            col = pCol2[i],
                            fill = pCol2[i],
                            lty = lty,
                            lwd = lwd,
                            lineend = lineend,
                            alpha = alpha
                        ),
                        arrow = bArrow
                    ),
                    data = facet_data,
                    xmin = brXmin[i],
                    xmax = brXmax[i],
                    ymin = brYmin,
                    ymax = brYmax
                )
            }
        } else if (addBranch == TRUE & annoPos %in% c("left", "right")) {
            for (i in seq_len(2 * nPoints)) {
                object <- object + annotation_custom2(
                    grob = grid::segmentsGrob(
                        gp = grid::gpar(
                            col = pCol2[i],
                            fill = pCol2[i],
                            lty = lty,
                            lwd = lwd,
                            lineend = lineend,
                            alpha = alpha
                        ),
                        arrow = bArrow
                    ),
                    data = facet_data,
                    xmin = brXmin,
                    xmax = brXmax,
                    ymin = brYmin[i],
                    ymax = brYmax[i]
                )
            }
        } else {
        }
    }
    if (is.null(textCol)) {
        textCol <- useMyCol("stallion", n = nPoints)
    } else {
        textCol <- textCol
    }
    if (aesGroup == FALSE) {
        textLabel <- textLabel
    } else {
        textLabel <- groupInfo[, 1]
    }
    if (is.null(facetName)) {
        if (addText == TRUE & annoPos %in% c("top", "botomn")) {
            for (i in seq_len(nPoints)) {
                object <- object + ggplot2::annotation_custom(
                    grob = grid::textGrob(
                        gp = grid::gpar(
                            col = textCol[i],
                            fontsize = textSize,
                            fontfamily = fontfamily,
                            fontface = fontface
                        ),
                        hjust = hjust,
                        vjust = vjust,
                        label = textLabel[i],
                        check.overlap = TRUE,
                        just = "centre",
                        rot = textRot
                    ),
                    xmin = ggplot2::unit(xmin[i], "native"),
                    xmax = ggplot2::unit(xmax[i], "native"),
                    ymin = ggplot2::unit(ymin + textHVjust, "native"),
                    ymax = ggplot2::unit(ymax + textHVjust, "native")
                )
            }
        } else if (addText == TRUE & annoPos %in% c("left", "right")) {
            for (i in seq_len(nPoints)) {
                object <- object + ggplot2::annotation_custom(
                    grob = grid::textGrob(
                        gp = grid::gpar(
                            col = textCol[i],
                            fontsize = textSize,
                            fontfamily = fontfamily,
                            fontface = fontface
                        ),
                        hjust = hjust,
                        vjust = vjust,
                        label = textLabel[i],
                        check.overlap = TRUE,
                        just = "centre",
                        rot = textRot
                    ),
                    xmin = ggplot2::unit(xmin + textHVjust, "native"),
                    xmax = ggplot2::unit(xmax + textHVjust, "native"),
                    ymin = ggplot2::unit(ymin[i], "native"),
                    ymax = ggplot2::unit(ymax[i], "native")
                )
            }
        } else {
        }
    } else {
        if (addText == TRUE & annoPos %in% c("top", "botomn")) {
            for (i in seq_len(nPoints)) {
                object <- object + annotation_custom2(
                    grob = grid::textGrob(
                        gp = grid::gpar(
                            col = textCol[i],
                            fontsize = textSize,
                            fontfamily = fontfamily,
                            fontface = fontface
                        ),
                        hjust = hjust,
                        vjust = vjust,
                        label = textLabel[i],
                        check.overlap = TRUE,
                        just = "centre",
                        rot = textRot
                    ),
                    data = facet_data,
                    xmin = xmin[i],
                    xmax = xmax[i],
                    ymin = ymin + textHVjust,
                    ymax = ymax + textHVjust
                )
            }
        } else if (addText == TRUE & annoPos %in% c("left", "right")) {
            for (i in seq_len(nPoints)) {
                object <- object + annotation_custom2(
                    grob = grid::textGrob(
                        gp = grid::gpar(
                            col = textCol[i],
                            fontsize = textSize,
                            fontfamily = fontfamily,
                            fontface = fontface
                        ),
                        hjust = hjust,
                        vjust = vjust,
                        label = textLabel[i],
                        check.overlap = TRUE,
                        just = "centre",
                        rot = textRot
                    ),
                    data = facet_data,
                    xmin = xmin + textHVjust,
                    xmax = xmin + textHVjust,
                    ymin = ymin[i],
                    ymax = ymax[i]
                )
            }
        } else {
        }
    }
    object
}
