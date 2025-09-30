# perform better Violins
#' @author Mariano Ruz Jurado
#' @title Violin Graph with wilcox test on single cell level
#' @description Creates a violin plot to compare gene expression across
#' different conditions or groups within a Seurat object. It incorporates
#' Wilcoxon rank-sum tests to evaluate statistical differences between
#' conditions. The plot can be customized with options for data transformation,
#' jitter display, and significance annotations. The function also supports
#' multiple conditions and allows for visualisation of statistical results from
#' wilcoxon test.
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
#' @param wilcox_test Bolean if TRUE a bonferoni wilcoxon test will be carried
#' out between ctrl.condition and the rest
#' @param stat_pos_mod value for modifiyng statistics height
#' @param hjust.wilcox value for adjusting height of the text
#' @param vjust.wilcox value for vertical of text
#' @param hjust.wilcox.2 value for adjusting height of the text, with group.by.2
#'  specified
#' @param vjust.wilcox.2 value for vertical of text, with group.by.2 specified
#' @param sign_bar adjusts the sign_bar with group.by.2 specified
#' @param size.wilcox value for size of text of statistical test
#' @param step_mod value for defining the space between one test and the next
#' one
#' @param geom_jitter_args_group_by2 controls the jittering of points if
#' group.by.2 is specified
#' @param SeuV5 Seuratv5 object? (TRUE or FALSE)
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
#'     SeuV5 = TRUE,
#'     Feature = "NKG7",
#'     ListTest = ListTest,
#'     ctrl.condition = "healthy",
#'     group.by = "condition"
#' )
#'
#' @export
DO.VlnPlot <- function(sce_object,
    SeuV5 = TRUE,
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
    wilcox_test = TRUE,
    stat_pos_mod = 1.15,
    hjust.wilcox = 0.8,
    vjust.wilcox = 2.0,
    size.wilcox = 3.33,
    step_mod = 0,
    hjust.wilcox.2 = 0.5,
    vjust.wilcox.2 = 0,
    sign_bar = 0.8) {
    # support for single cell experiment objects
    if (is(sce_object, "SingleCellExperiment")) {
        SCE <- TRUE
        sce_object <- as.Seurat(sce_object)
    } else {
        SCE <- FALSE
    }

    if (!(Feature %in% rownames(sce_object)) &&
        !(Feature %in% names(sce_object@meta.data))) {
        stop("Feature not found in Seurat Object!")
    }

    if (wilcox_test == TRUE) {
        rstat <- system.file(package = "rstatix") # package is installed
        ifelse(nzchar(rstat),
            "",
            stop("Install rstatix R package for wilcox statistic!")
        )
    }

    if (is.null(ctrl.condition)) {
        stop("Please specify the ctrl condition as string!")
    }

    if (SeuV5 == TRUE) {
        rlang::warn(
            "SeuV5 TRUE, if object Seuratv4 or below change SeuV5 to FALSE",
            .frequency = "once",
            .frequency_id = "v5Mean"
        )
        if (Feature %in% rownames(sce_object)) {
            vln.df <- data.frame(
                Feature = sce_object[["RNA"]]$data[Feature, ],
                cluster = sce_object[[group.by]]
            )
        } else {
            vln.df <- data.frame(
                Feature = FetchData(sce_object, vars = Feature)[, 1],
                cluster = sce_object[[group.by]]
            )
        }



        df <- data.frame(
            group = setNames(
                sce_object[[group.by]][, group.by],
                rownames(sce_object[[group.by]])
            ),
            orig.ident = sce_object$orig.ident
        )

        # add second group for individual splitting and testing in the wilcoxon
        if (!is.null(group.by.2)) {
            vln.df[group.by.2] <- sce_object[[group.by.2]]
            df[group.by.2] <- sce_object[[group.by.2]]
        }

        # get expression values for genes from individual cells, add to df
        if (Feature %in% rownames(sce_object)) {
            df[, Feature] <- sce_object@assays$RNA$data[Feature, ]
        } else {
            df[, Feature] <- FetchData(sce_object, vars = Feature)[, 1]
        }
    }

    if (SeuV5 == FALSE) {
        if (Feature %in% rownames(sce_object)) {
            vln.df <- data.frame(
                Feature = sce_object[["RNA"]]@data[Feature, ],
                cluster = sce_object[[group.by]]
            )
        } else {
            vln.df <- data.frame(
                Feature = FetchData(sce_object, vars = Feature)[, 1],
                cluster = sce_object[[group.by]]
            )
        }

        df <- data.frame(
            group = setNames(
                sce_object[[group.by]][, group.by],
                rownames(sce_object[[group.by]])
            ),
            orig.ident = sce_object$orig.ident
        )
        # add second group for individual splitting and testing in the wilcoxon
        if (!is.null(group.by.2)) {
            vln.df[group.by.2] <- sce_object[[group.by.2]]
            df[group.by.2] <- sce_object[[group.by.2]]
        }

        # get expression values for genes from individual cells, add to df
        if (Feature %in% rownames(sce_object)) {
            df[, Feature] <- sce_object@assays$RNA@data[Feature, ]
        } else {
            df[, Feature] <- FetchData(sce_object, vars = Feature)[, 1]
        }
    }


    df.melt <- melt(df)

    vln.df$group <- factor(
        vln.df[[group.by]],
        levels = c(
            as.character(ctrl.condition),
            levels(factor(vln.df[[group.by]]))[
                !(levels(factor(vln.df[[group.by]])) %in% ctrl.condition)
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
            match(x[2], vln.df[[group.by.2]])
        }, integer(1))
    } else {
        indices <- vapply(ListTest, function(x) {
            match(x[2], vln.df[[group.by]])
        }, integer(1))
    }
    ListTest <- ListTest[order(indices)]

    remove_zeros <- function(lst, df) {
        lst_filtered <- lst
        for (i in seq_along(lst)) {
            elements <- lst[[i]]
            if (all(df[df$group %in% elements, "Mean"] == 0)) {
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
        df.melt.sum <- df.melt %>%
            dplyr::group_by(group, variable) %>%
            dplyr::summarise(Mean = mean(value))
    } else {
        df.melt.sum <- df.melt %>%
            dplyr::group_by(group, !!sym(group.by.2), variable) %>%
            dplyr::summarise(Mean = mean(value))
    }


    # Remove vectors with both elements having a mean of 0
    ListTest <- remove_zeros(ListTest, df.melt.sum)

    # check there are groups in the data which contain only 0 values
    # and therefore let the test fail
    if (is.null(group.by.2)) {
        group_of_zero <- df.melt %>%
            dplyr::group_by(group) %>%
            summarise(all_zeros = all(value == 0), .groups = "drop") %>%
            filter(all_zeros)

        if (nrow(group_of_zero) > 0) {
            warning(
                "Some comparisons have no expression in both groups, setting ",
                "expression to minimum value to ensure test does not fail!"
            )
            df.melt <- df.melt %>%
                dplyr::group_by(group) %>%
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
        group_of_zero <- df.melt %>%
            dplyr::group_by(group, !!sym(group.by.2)) %>%
            summarise(all_zeros = all(value == 0), .groups = "drop") %>%
            filter(all_zeros)

        # check now the result for multiple entries in group.by.2
        groupby2_check <- group_of_zero %>%
            dplyr::group_by(!!sym(group.by.2)) %>%
            summarise(group_count = n_distinct(group), .groups = "drop") %>%
            filter(group_count > 1)

        if (nrow(groupby2_check) > 0) {
            warning(
                "Some comparisons have no expression in both groups, setting ",
                "expression to minimum value to ensure test does not fail!"
            )
            df.melt <- df.melt %>%
                dplyr::group_by(group, !!sym(group.by.2)) %>%
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

    # do statistix with rstatix + stats package
    if (wilcox_test == TRUE & is.null(group.by.2)) {
        stat.test <- df.melt %>%
            ungroup() %>%
            rstatix::wilcox_test(value ~ group,
                comparisons = ListTest,
                p.adjust.method = "none"
            ) %>%
            rstatix::add_significance()
        stat.test$p.adj <- stats::p.adjust(stat.test$p,
            method = "bonferroni",
            n = length(rownames(sce_object))
        )
        stat.test$p.adj <- ifelse(
            stat.test$p.adj == 0,
            sprintf("%.2e", .Machine$double.xmin),
            sprintf("%.2e", stat.test$p.adj)
        )
        stat.test$p <- ifelse(
            stat.test$p == 0,
            sprintf("%.2e", .Machine$double.xmin),
            sprintf("%.2e", stat.test$p)
        )
    }

    if (wilcox_test == TRUE & !is.null(group.by.2)) {
        stat.test <- df.melt %>%
            dplyr::group_by(!!sym(group.by.2)) %>%
            rstatix::wilcox_test(value ~ group,
                comparisons = ListTest,
                p.adjust.method = "none"
            ) %>%
            rstatix::add_significance()
        stat.test$p.adj <- stats::p.adjust(stat.test$p,
            method = "bonferroni",
            n = length(rownames(sce_object))
        )
        stat.test$p.adj <- ifelse(
            stat.test$p.adj == 0,
            sprintf("%.2e", .Machine$double.xmin),
            sprintf("%.2e", stat.test$p.adj)
        )
        stat.test$p <- ifelse(
            stat.test$p == 0,
            sprintf("%.2e", .Machine$double.xmin),
            sprintf("%.2e", stat.test$p)
        )
    }

    if (length(unique(vln.df[[group.by]])) > length(vector_colors)) {
        stop(sprintf(
            "Only %s colors provided, but %s needed!",
            length(vector_colors),
            length(unique(vln.df[[group.by]]))
        ))
    }

    # normal violin
    if (is.null(group.by.2)) {
        p <- ggplot(vln.df, aes(x = group, y = Feature)) +
            geom_violin(aes(fill = group), trim = TRUE, scale = "width", ) +
            geom_jitter(
                size = geom_jitter_args[1],
                width = geom_jitter_args[2],
                alpha = geom_jitter_args[3]
            ) +
            labs(title = Feature, y = "Expression Level") +
            xlab("") +
            ylab("") +
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

        if (wilcox_test == TRUE) {
            p <- p + stat_pvalue_manual(
                stat.test,
                label = p_label,
                y.position = max(vln.df$Feature) * 1.15,
                step.increase = 0.2
            )
        }
        return(p)
    }


    if (!is.null(group.by.2)) {
        # plot
        p <- ggplot(vln.df, aes(
            x = !!sym(group.by.2),
            y = Feature,
            fill = !!sym(group.by)
        )) +
            geom_violin(aes(fill = group), trim = TRUE, scale = "width") +
            labs(title = Feature, y = "log(nUMI)") +
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

        p2 <- suppressWarnings(
            ggplot(vln.df, aes(
                x = !!sym(group.by.2),
                y = Feature,
                fill = factor(!!sym(group.by), levels = levels(vln.df$group))
            )) +
                geom_boxplot(
                    width = .1,
                    color = "grey",
                    position = position_dodge(width = 0.9),
                    outlier.shape = NA
                ) +
                xlab("") +
                scale_fill_manual(
                    values = rep("black", length(vector_colors)), name =
                        group.by
                ) +
                theme_classic() +
                theme(
                    plot.title = element_text(
                        face = "bold",
                        color = "transparent",
                        hjust = 0.5,
                        size = 14
                    ),
                    axis.title.y = element_text(
                        face = "bold",
                        color = "transparent",
                        size = 14
                    ),
                    axis.text.x = element_text(
                        face = "bold",
                        color = "transparent",
                        angle = 45,
                        hjust = 1,
                        size = 14
                    ),
                    axis.text.y = element_text(
                        face = "bold",
                        color = "transparent",
                        hjust = 1,
                        size = 14
                    ),
                    legend.position = "bottom",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_blank(),
                    panel.background = element_rect(
                        fill = "transparent",
                        colour = NA
                    ),
                    plot.background = element_rect(
                        fill = "transparent",
                        colour = NA
                    ),
                    strip.background = element_rect(
                        fill = "transparent",
                        colour = NA
                    ),
                    axis.ticks = element_blank(),
                    legend.background = element_rect(
                        fill = "transparent",
                        color = NA
                    ),
                    # Transparent legend background
                    legend.key = element_rect(
                        fill = "transparent",
                        color = NA
                    ),
                    # Transparent legend keys
                    legend.title = element_text(
                        face = "bold",
                        color = "transparent"
                    ),
                    # Legend title styling
                    legend.text = element_text(color = "transparent"),
                    # strip.text = element_blank(),
                )
        )

        if (wilcox_test == TRUE & !is.null(group.by.2)) {
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

            # dplyr::select(x.axis, y.position, p.adj)

            # if Feature not a gene than use the uncorrected p
            if (Feature %in% rownames(sce_object)) {
                p_label <- "p = {p.adj}"
            } else {
                p_label <- "p = {p}"
            }

            p <- p + stat_pvalue_manual(stat.test_plot,
                label = p_label,
                y.position = "y.position",
                # x="x",
                xmin = "xmin",
                xmax = "xmax",
                # xend="xend",
                # step.increase = 0.2,
                inherit.aes = FALSE,
                size = size.wilcox,
                angle = 0,
                hjust = hjust.wilcox.2,
                vjust = vjust.wilcox.2,
                tip.length = 0.02,
                bracket.size = sign_bar
            )

            p2 <- p2 + stat_pvalue_manual(stat.test_plot,
                label = p_label,
                y.position = "y.position",
                # x="x",
                xmin = "xmin",
                xmax = "xmax",
                # xend="xend",
                # step.increase = 0.2,
                inherit.aes = FALSE,
                size = size.wilcox,
                angle = 0,
                hjust = hjust.wilcox.2,
                vjust = vjust.wilcox.2,
                tip.length = 0.02,
                bracket.size = sign_bar,
                color = "transparent"
            )
        }
        plot_p <- cowplot::ggdraw() +
            cowplot::draw_plot(p) +
            cowplot::draw_plot(p2)
        plot_p
    }



    if (returnValues == TRUE) {
        returnList <- list(vln.df, df.melt, stat.test)
        names(returnList) <- c("vln.df", "df.melt", "stat.test")
        return(returnList)
    }

    return(plot_p)
}
