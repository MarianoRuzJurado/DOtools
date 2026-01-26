# perform better Violins
#' @author Mariano Ruz Jurado
#' @title Violin Graph with wilcox test on single cell level
#' @description Creates a violin plot to compare gene expression across
#' different conditions or groups within a Seurat object. It incorporates
#' different tests to evaluate statistical differences between
#' conditions. The plot can be customized with options for data transformation,
#' jitter display, and significance annotations. The function also supports
#' multiple conditions and allows for visualisation of statistical results from
#' different test.
#' @param sce_object combined SCE object or Seurat
#' @param Feature name of the feature
#' @param ListTest List for which conditions wilcox will be performed, if NULL
#' always CTRL group against everything
#' @param returnValues return df.melt.sum data frame containing means and SEM
#' for the set group
#' @param ctrl.condition set your ctrl condition, relevant if running with empty
#'  comparison List
#' @param group.by select the seurat sce_object slot where your conditions can
#' be found, default conditon
#' @param group.by.2 relevant for multiple group testing, e.g. for each
#' cell type the test between each of them in two conditions provided
#' @param geom_jitter_args vector for dots visualisation in vlnplot:
#' size, width, alpha value
#' @param vector_colors specify a minimum number of colours as you have entries
#' in your condition, default 2
#' @param test_use perform one of c(
#' "wilcox", "wilcox_limma", "bimod", "t", "negbinom",
#' "poisson", "LR", "MAST", "DESeq2", "none"
#' ). default "wilcox"
#' @param correction_method correction for p-value calculation. One of
#' c("BH", "bonferroni", "holm", "BY", "fdr", "none"). default "fdr"
#' @param stat_pos_mod value for modifiyng statistics height
#' @param hjust_test value for adjusting height of the text
#' @param vjust_test value for vertical of text
#' @param hjust_test_2 value for adjusting height of the text, with group.by.2
#'  specified
#' @param vjust_test_2 value for vertical of text, with group.by.2 specified
#' @param sign_bar adjusts the sign_bar with group.by.2 specified
#' @param size_test value for size of text of statistical test
#' @param step_mod value for defining the space between one test and the next
#' one
#' @param geom_jitter_args_group_by2 controls the jittering of points if
#' group.by.2 is specified
#' @param y_title specify title on the y axis. default "log(nUMI)"
#' @param p_values Manually providing p-values for plotting, be aware of
#' group size and if necessary make your test return the same amount of values
#' @param random_seed parameter for random state initialisation
#'
#' @import ggplot2
#' @import ggpubr
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#' @import reshape2
#' @importFrom SeuratObject as.Seurat
#'
#' @return a ggplot or a list used data frames
#'
#' @examples
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' ListTest <- list()
#' ListTest[[1]] <- c("healthy", "disease")
#'
#' DO.VlnPlot(
#'     sce_object = sce_data,
#'     Feature = "NKG7",
#'     ListTest = ListTest,
#'     ctrl.condition = "healthy",
#'     group.by = "condition"
#' )
#'
#' @export
DO.VlnPlot <- function(
    sce_object,
    Feature,
    ListTest = NULL,
    returnValues = FALSE,
    ctrl.condition = NULL,
    group.by = "condition",
    group.by.2 = NULL,
    geom_jitter_args = c(0.20, 0.25, 0.25),
    geom_jitter_args_group_by2 = c(0.1, 0.1, 1),
    vector_colors = c(
        "#1f77b4",
        "#ea7e1eff",
        "royalblue4",
        "tomato2",
        "darkgoldenrod",
        "palegreen4",
        "maroon",
        "thistle3"
    ),
    test_use = "wilcox",
    correction_method = "fdr",
    p_values = NULL,
    y_title = "log(nUMI)",
    stat_pos_mod = 1.15,
    hjust_test = 0.8,
    vjust_test = 2.0,
    size_test = 3.33,
    step_mod = 0,
    hjust_test_2 = 0.5,
    vjust_test_2 = 0,
    sign_bar = 0.8,
    random_seed = 42
) {
    # support for single cell experiment objects
    if (methods::is(sce_object, "SingleCellExperiment")) {
        SCE <- TRUE
        sce_object <- .suppressDeprecationWarnings(as.Seurat(sce_object))
    } else {
        SCE <- FALSE
    }

    if (!(Feature %in% rownames(sce_object)) &&
        !(Feature %in% names(sce_object@meta.data))) {
        stop("Feature not found in SCE Object!")
    }
    ## add test and correction methods
    test <- match.arg(
        test_use, c(
            "wilcox", "wilcox_limma", "bimod", "t", "negbinom",
            "poisson", "LR", "MAST", "DESeq2", "none"
        )
    )
    p_method <- match.arg(
        correction_method, c("BH", "bonferroni", "holm", "BY", "fdr", "none")
    )
    ##

    if (test != "none" & is.null(ctrl.condition)) {
        stop("Please specify the ctrl condition as string!")
    }

    ## add similar df expression extraction as in DO.Barplot.R

    # create data frame with conditions from provided sce_object
    vln_df <- data.frame(
        condition = stats::setNames(
            sce_object[[group.by]][, group.by],
            rownames(sce_object[[group.by]])
        ),
        orig.ident = sce_object$orig.ident
    )

    # add values: extracts log1p values or metadata values
    vln_df[, Feature] <- FetchData(sce_object, vars = Feature)

    # add second group for individual splitting and testing in the wilcoxon
    if (!is.null(group.by.2)) {
        vln_df[group.by.2] <- sce_object[[group.by.2]]
        # df[group.by.2] <- sce_object[[group.by.2]]
    }


    vln_df_melt <- melt(vln_df)


    vln_df$group <- factor(
        vln_df$condition,
        levels = c(
            as.character(ctrl.condition),
            levels(factor(vln_df$condition))[
                !(levels(factor(vln_df$condition)) %in% ctrl.condition)
            ]
        )
    )


    # create comparison list for wilcox, always against control
    # ,alternative add your own list as argument
    if (is.null(ListTest)) {
        .logger("ListTest empty, comparing every sample with each other")
        group <- unique(sce_object[[group.by]][, group.by])
        # set automatically ctrl condition if not provided
        if (is.null(ctrl.condition)) {
            ctrl.condition <- group[grep(
                pattern = paste(c(
                    "CTRL",
                    "Ctrl",
                    "ctrl",
                    "WT",
                    "Wt",
                    "wt"
                ), collapse = "|"),
                group
            )[1]]
        }


        # create ListTest
        ListTest <- list()
        for (i in seq_along(group)) {
            cndtn <- as.character(group[i])
            if (cndtn != ctrl.condition) {
                ListTest[[length(ListTest) + 1]] <- c(ctrl.condition, cndtn)
            }
        }
    }
    # delete Null values
    ListTest <- ListTest[!vapply(ListTest, is.null, logical(1))]
    if (!is.null(group.by.2)) {
        indices <- vapply(ListTest, function(x) {
            match(x[2], vln_df[[group.by.2]])
        }, integer(1))
    } else {
        indices <- vapply(ListTest, function(x) {
            match(x[2], vln_df[[group.by]])
        }, integer(1))
    }
    ListTest <- ListTest[order(indices)]

    remove_zeros <- function(lst, df) {
        lst_filtered <- lst
        for (i in seq_along(lst)) {
            elements <- lst[[i]]
            if (all(df[df$condition %in% elements, "Mean"] == 0)) {
                lst_filtered <- lst_filtered[-i]
                warning(sprintf(
                    "Removing Test %s vs %s since both values are 0",
                    elements[1], elements[2]
                ))
            }
        }
        return(lst_filtered)
    }

    # group results and summarize
    if (is.null(group.by.2)) {
        vln_df_sum <- vln_df_melt %>%
            dplyr::group_by(condition, variable) %>%
            dplyr::summarise(Mean = mean(value))
    } else {
        vln_df_sum <- vln_df_melt %>%
            dplyr::group_by(condition, !!sym(group.by.2), variable) %>%
            dplyr::summarise(Mean = mean(value))
    }


    # Remove vectors with both elements having a mean of 0
    ListTest <- remove_zeros(ListTest, vln_df_sum)


    # check there are groups in the data which contain only 0 values
    # and therefore let the test fail
    if (is.null(group.by.2)) {
        group_of_zero <- vln_df_melt %>%
            dplyr::group_by(condition) %>%
            summarise(all_zeros = all(value == 0), .groups = "drop") %>%
            filter(all_zeros)

        if (nrow(group_of_zero) > 0) {
            warning(
                "Some comparisons have no expression in both groups, setting ",
                "expression to minimum value to ensure test does not fail!"
            )
            vln_df_melt <- vln_df_melt %>%
                dplyr::group_by(condition) %>%
                dplyr::mutate(
                    value = if_else(
                        row_number() == 1 & all(value == 0),
                        .Machine$double.xmin,
                        value
                    )
                ) %>%
                ungroup()
        }
    } else {
        group_of_zero <- vln_df_melt %>%
            dplyr::group_by(condition, !!sym(group.by.2)) %>%
            summarise(all_zeros = all(value == 0), .groups = "drop") %>%
            filter(all_zeros)

        # check now the result for multiple entries in group.by.2
        groupby2_check <- group_of_zero %>%
            dplyr::group_by(!!sym(group.by.2)) %>%
            summarise(group_count = n_distinct(condition), .groups = "drop") %>%
            filter(group_count > 1)

        if (nrow(groupby2_check) > 0) {
            warning(
                "Some comparisons have no expression in both groups, setting ",
                "expression to minimum value to ensure test does not fail!"
            )
            vln_df_melt <- vln_df_melt %>%
                dplyr::group_by(condition, !!sym(group.by.2)) %>%
                dplyr::mutate(
                    value = if_else(
                        row_number() == 1 & all(value == 0),
                        .Machine$double.xmin,
                        value
                    )
                ) %>%
                ungroup()
        }
    }
    ### new test, useable with all methods implemented in FindMarkers
    ### and with multiple correction methods
    ### does also take manual p-values as input
    if (test_use != "none" & is.null(group.by.2)) {
        stat.test <- data.frame()
        for (grp in ListTest) {
            # If feature is a gene
            if (Feature %in% rownames(sce_object)) {
                degs <- FindMarkers(sce_object,
                    test.use = test_use,
                    ident.1 = grp[2],
                    ident.2 = grp[1],
                    logfc.threshold = 0,
                    min.pct = 0,
                    min.diff.pct = -Inf,
                    group.by = group.by
                )
                degs$p_val_adj <- p.adjust(degs$p_val, method = p_method)
                degs <- degs[rownames(degs) %in% Feature, ]
            }

            # If feature is a meta data column
            if (!Feature %in% rownames(sce_object)) {
                mat_Feature <- t(as.matrix(sce_object@meta.data[Feature]))
                .suppressAllWarnings(
                    sce_object[[Feature]] <- CreateAssayObject(
                        data = mat_Feature
                    )
                )
                degs <- FindMarkers(sce_object,
                    assay = Feature,
                    test.use = test_use,
                    ident.1 = grp[2], ident.2 = grp[1], logfc.threshold = 0,
                    min.pct = 0, min.diff.pct = -Inf, group.by = group.by
                )
            }

            group_dis <- grp[2]
            group_ctrl <- grp[1]

            test_df <- degs %>%
                rownames_to_column(var = ".y.") %>%
                mutate(
                    group1 = group_dis,
                    group2 = group_ctrl,
                    n1 = sum(sce_object[[group.by]][, group.by] == group1),
                    n2 = sum(sce_object[[group.by]][, group.by] == group2),
                    statistic = NA_real_, # keep this column for consistency
                    p = p_val,
                    p.adj = p_val_adj
                ) %>%
                rstatix::add_significance() %>%
                select(
                    .y., group1, group2, n1, n2,
                    statistic, p, p.adj, p.adj.signif
                )

            if (!is.null(p_values)) {
                # Check for equal numbers of  p-values and comparisons
                if (length(p_values) != length(ListTest)) {
                    stop(sprintf(
                        "Number of provided p-values: %s, does",
                        "not match number of comparisons: %s",
                        length(p_values), length(ListTest)
                    ))
                }
                # Assign new p-values and add significance
                test_df$p.adj <- p_values
                test_df <- test_df %>% rstatix::add_significance()
            }

            # add lowest number to replace 0s and apply scientifc writing
            test_df$p.adj <- ifelse(
                test_df$p.adj == 0,
                sprintf("%.2e", .Machine$double.xmin),
                sprintf("%.2e", test_df$p.adj)
            )
            test_df$p <- ifelse(
                test_df$p == 0,
                sprintf("%.2e", .Machine$double.xmin),
                sprintf("%.2e", test_df$p)
            )
            stat.test <- rbind(stat.test, test_df)
        }
    }

    ## for multiple group argument set
    if (test_use != "none" & !is.null(group.by.2)) {
        stat.test <- data.frame()

        for (grp2 in sort(unique(vln_df_melt[[group.by.2]]))) {
            sce_object_sub <- subset(sce_object, !!sym(group.by.2) == grp2)

            for (grp in ListTest) {
                # If feature is a gene
                if (Feature %in% rownames(sce_object)) {
                    degs <- FindMarkers(sce_object_sub,
                        test.use = test_use,
                        ident.1 = grp[2],
                        ident.2 = grp[1],
                        logfc.threshold = 0,
                        min.pct = 0,
                        min.diff.pct = -Inf,
                        group.by = group.by
                    )
                    degs$p_val_adj <- p.adjust(degs$p_val, method = p_method)
                    degs <- degs[rownames(degs) %in% Feature, ]
                }
                # If feature is a meta data column
                if (!Feature %in% rownames(sce_object)) {
                    mat_Feature <- t(
                        as.matrix(sce_object_sub@meta.data[Feature])
                    )
                    .suppressAllWarnings(
                        sce_object_sub[[Feature]] <- CreateAssayObject(
                            data = mat_Feature
                        )
                    )
                    degs <- FindMarkers(sce_object_sub,
                        assay = Feature,
                        test.use = test_use,
                        ident.1 = grp[2], ident.2 = grp[1], logfc.threshold = 0,
                        min.pct = 0, min.diff.pct = -Inf, group.by = group.by
                    )
                }
                group_dis <- grp[2]
                group_ctrl <- grp[1]

                test_df <- degs %>%
                    rownames_to_column(var = ".y.") %>%
                    mutate(
                        group1 = group_dis,
                        group2 = group_ctrl,
                        n1 = sum(
                            sce_object[[group.by]][, group.by] == group1
                        ),
                        n2 = sum(
                            sce_object[[group.by]][, group.by] == group2
                        ),
                        statistic = NA_real_, # for consistency
                        p = p_val,
                        p.adj = p_val_adj
                    ) %>%
                    rstatix::add_significance() %>%
                    select(
                        .y., group1, group2, n1, n2,
                        statistic, p, p.adj, p.adj.signif
                    )
                if (!is.null(p_values)) {
                    # Check for equal numbers of  p-values and comparisons
                    if (length(p_values) != length(ListTest)) {
                        stop(sprintf(
                            "Number of provided p-values: %s, does",
                            "not match number of comparisons: %s",
                            length(p_values), length(ListTest)
                        ))
                    }
                    # Assign new p-values and add significance
                    test_df$p.adj <- p_values
                    test_df <- test_df %>% rstatix::add_significance()
                }

                # replace 0s and apply scientific writing
                test_df$p.adj <- ifelse(
                    test_df$p.adj == 0,
                    sprintf("%.2e", .Machine$double.xmin),
                    sprintf("%.2e", test_df$p.adj)
                )
                test_df$p <- ifelse(
                    test_df$p == 0,
                    sprintf("%.2e", .Machine$double.xmin),
                    sprintf("%.2e", test_df$p)
                )

                # Add information of group.by.2 to test df
                test_df <- cbind(
                    setNames(data.frame(grp2), group.by.2),
                    test_df
                )
                stat.test <- rbind(stat.test, test_df)
            }
        }
    }
    if (length(unique(vln_df[[group.by]])) > length(vector_colors)) {
        stop(sprintf(
            "Only %s colors provided, but %s needed!",
            length(vector_colors),
            length(unique(vln_df[[group.by]]))
        ))
    }

    #plotting with reproducibility
    set.seed(random_seed)

    # normal violin
    if (is.null(group.by.2)) {
        p <- ggplot(vln_df, aes(x = condition, y = !!sym(Feature))) +
            geom_violin(aes(fill = condition), trim = TRUE, scale = "width", ) +
            geom_jitter(
                size = geom_jitter_args[1],
                width = geom_jitter_args[2],
                alpha = geom_jitter_args[3]
            ) +
            labs(title = Feature, y = y_title) +
            xlab("") +
            # ylab("") +
            theme_classic() +
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
                axis.text.x = element_text(
                    face = "bold",
                    color = "black",
                    angle = 45,
                    hjust = 1,
                    size = 14
                ),
                axis.text.y = element_text(
                    face = "bold",
                    color = "black",
                    hjust = 1,
                    size = 14
                ),
                legend.position = "none"
            ) +
            scale_fill_manual(values = vector_colors)

        if (Feature %in% rownames(sce_object)) {
            p_label <- "p = {p.adj}"
        } else {
            p_label <- "p = {p}"
        }

        if (test_use != "none") {
            p <- p + stat_pvalue_manual(
                stat.test,
                label = p_label,
                y.position = max(vln_df[[Feature]]) * 1.15,
                step.increase = 0.2
            )
        }
    }


    if (!is.null(group.by.2)) {
        # better boxplot alignment
        p <- ggplot(vln_df, aes(
            x = !!sym(group.by.2),
            y = !!sym(Feature),
            fill = !!sym(group.by)
        )) +
            geom_violin(
                aes(fill = !!sym(group.by)),
                trim = TRUE,
                scale = "width",
                position = position_dodge(0.9)
            ) +
            geom_boxplot(
                aes(group = interaction(!!sym(group.by.2), !!sym(group.by))),
                fill = "black",
                color = "grey",
                width = 0.1,
                alpha = 0.7,
                position = position_dodge(0.9),
                outlier.shape = NA
            ) +
            labs(title = Feature, y = y_title) +
            xlab("") +
            theme_classic() +
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
                axis.text.x = element_text(
                    face = "bold",
                    color = "black",
                    angle = 45,
                    hjust = 1,
                    size = 14
                ),
                axis.text.y = element_text(
                    face = "bold",
                    color = "black",
                    hjust = 1,
                    size = 14
                ),
                legend.position = "bottom",
                panel.grid.major = element_line(
                    colour = "grey90",
                    linetype = "dotted"
                ),
                panel.grid.minor = element_line(
                    colour = "grey90",
                    linetype = "dotted"
                ),
                axis.line = element_line(colour = "black"),
                strip.background = element_rect(
                    fill = "lightgrey",
                    colour = "black",
                    linewidth = 1
                ),
                strip.text = element_text(colour = "black", size = 12),
            ) +
            scale_fill_manual(values = vector_colors)

        # Changes: now integrated into the main plot

        if (test_use != "none" & !is.null(group.by.2)) {
            if (Feature %in% rownames(sce_object)) {
                stat.test_plot <- stat.test %>%
                    mutate(
                        y.position = seq(
                            from = max(
                                sce_object@assays$RNA$data[Feature, ][
                                    !is.na(
                                        sce_object@assays$RNA$data[Feature, ]
                                    )
                                ]
                            ) * stat_pos_mod,
                            by = step_mod,
                            length.out = nrow(stat.test)
                        )
                    ) %>%
                    mutate(
                        x = as.numeric(
                            factor(
                                stat.test[[group.by.2]],
                                levels = unique(stat.test[[group.by.2]])
                            )
                        ),
                        xmin = as.numeric(
                            factor(
                                stat.test[[group.by.2]],
                                levels = unique(stat.test[[group.by.2]])
                            )
                        ) - 0.2,
                        xmax = as.numeric(
                            factor(
                                stat.test[[group.by.2]],
                                levels = unique(stat.test[[group.by.2]])
                            )
                        ) + 0.2
                    )
            } else {
                stat.test_plot <- stat.test %>%
                    mutate(
                        y.position = seq(
                            from = max(
                                sce_object[[Feature]][, 1][
                                    !is.na(sce_object[[Feature]][, 1])
                                ]
                            ) * stat_pos_mod,
                            by = step_mod,
                            length.out = nrow(stat.test)
                        )
                    ) %>%
                    mutate(
                        x = as.numeric(
                            factor(
                                stat.test[[group.by.2]],
                                levels = unique(stat.test[[group.by.2]])
                            )
                        ),
                        xmin = as.numeric(
                            factor(
                                stat.test[[group.by.2]],
                                levels = unique(stat.test[[group.by.2]])
                            )
                        ) - 0.2,
                        xmax = as.numeric(
                            factor(
                                stat.test[[group.by.2]],
                                levels = unique(stat.test[[group.by.2]])
                            )
                        ) + 0.2
                    )
            }

            if (Feature %in% rownames(sce_object)) {
                p_label <- "p = {p.adj}"
            } else {
                p_label <- "p = {p}"
            }

            p <- p + stat_pvalue_manual(stat.test_plot,
                label = p_label,
                y.position = "y.position",
                xmin = "xmin",
                xmax = "xmax",
                inherit.aes = FALSE,
                size = size_test,
                angle = 0,
                hjust = hjust_test_2,
                vjust = vjust_test_2,
                tip.length = 0.02,
                bracket.size = sign_bar
            )
        }

        plot_p <- p
    }

    if (returnValues == TRUE) {
        returnList <- list(vln_df, vln_df_melt, stat.test)
        names(returnList) <- c("vln_df", "vln_df_melt", "stat.test")
        return(returnList)
    }

    return(p)
}
