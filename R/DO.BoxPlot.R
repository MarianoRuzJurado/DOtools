# Boxplot function for one or two given groups per gene
#' @author Mariano Ruz Jurado
#' @title Box Graph with wilcox test on single cell level
#' @description Creates a box plot using a pseudo-bulk approach and performs a
#' Wilcoxon test on single-cell level. Allows customization of outlier removal,
#' statistical labels, and color schemes. Supports comparison of conditions
#' with optional second grouping. Useful for visualizing gene expression and
#' statistical differences.
#' @param sce_object The SCE object or Seurat
#' @param Feature name of the feature/gene
#' @param sample.column meta data column containing sample IDs
#' @param ListTest List for which conditions wilcox will be performed, if NULL
#' always CTRL group against everything
#' @param plot_sample Plot individual sample dot in graph
#' @param group.by group name to look for in meta data
#' @param group.by.2 second group name to look for in meta data
#' @param ctrl.condition select condition to compare to
#' @param outlier_removal Outlier calculation
#' @param vector_colors get the colours for the plot
#' @param wilcox_test If you want to have wilcoxon performed, boolean default
#' TRUE
#' @param stat_pos_mod modificator for where the p-value is plotted increase
#' for higher
#' @param hjust.wilcox value for adjusting height of the text
#' @param vjust.wilcox value for vertical of text
#' @param hjust.wilcox.2 value for adjusting height of the text, with
#' group.by.2 specified
#' @param vjust.wilcox.2 value for vertical of text, with group.by.2 specified
#' @param sign_bar adjusts the sign_bar with group.by.2 specified
#' @param size.wilcox value for size of text of statistical test
#' @param step_mod value for defining the space between one test and the
#' next one
#' @param orderAxis vector for xaxis sorting, alphabetically by default
#'
#' @import ggplot2
#' @import ggpubr
#' @import tidyverse
#' @import dplyr
#' @import magrittr
#' @import Seurat
#'
#' @return a ggplot
#'
#' @examples
#' sce_data <-
#'   readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' set.seed(123)
#' sce_data$orig.ident <-
#'   sample(rep(c("A", "B", "C"), length.out = ncol(sce_data)))
#'
#' ListTest <- list()
#' ListTest[[1]] <- c("healthy", "disease")
#'
#' DO.BoxPlot(
#'   sce_object = sce_data,
#'   Feature = "NKG7",
#'   sample.column="orig.ident",
#'   ListTest = ListTest,
#'   group.by="condition",
#'   ctrl.condition = "healthy",
#' )
#'
#' @export
DO.BoxPlot <- function(sce_object,
                       Feature,
                       sample.column = "orig.ident",
                       ListTest=NULL,
                       group.by = "condition",
                       group.by.2 = NULL,
                       ctrl.condition=NULL,
                       outlier_removal = TRUE,
                       plot_sample=TRUE,
                       vector_colors = c("#1f77b4",
                                         "#ea7e1eff",
                                         "royalblue4",
                                         "tomato2",
                                         "darkgoldenrod",
                                         "palegreen4",
                                         "maroon",
                                         "thistle3"),
                       wilcox_test = TRUE,
                       stat_pos_mod = 1.15,
                       step_mod = 0,
                       hjust.wilcox=0.5,
                       vjust.wilcox = 0.25,
                       size.wilcox=3.33,
                       hjust.wilcox.2=0.5,
                       vjust.wilcox.2=0,
                       sign_bar=0.8,
                       orderAxis=NULL
){

  #support for single cell experiment objects
  if (is(sce_object, "SingleCellExperiment")) {
    sce_object <- as.Seurat(sce_object)
  }

  #aggregate expression, pseudobulk to visualize the boxplot
  if (is.null(group.by.2)) {
    pseudo_Seu <- Seurat::AggregateExpression(sce_object,
                                              assays = "RNA",
                                              return.seurat = TRUE,
                                              group.by = c(group.by,
                                                           sample.column),
                                              verbose = FALSE)

    pseudo_Seu$celltype.con <- pseudo_Seu[[group.by]]

  } else{
    pseudo_Seu <- Seurat::AggregateExpression(sce_object,
                                              assays = "RNA",
                                              return.seurat = TRUE,
                                              group.by = c(group.by,
                                                           group.by.2,
                                                           sample.column),
                                              verbose = FALSE)

    #cover the case of subsetted to only have one cell type
    if (length(unique(sce_object@meta.data[[group.by.2]])) == 1) {
      pseudo_Seu@meta.data[[group.by.2]] <-
        unique(sce_object@meta.data[[group.by.2]])
      pseudo_Seu@meta.data[[group.by]] <-
        gsub(paste0(".*(", paste(unique(
          sce_object$condition
        ), collapse = "|"), ").*"), "\\1", pseudo_Seu@meta.data[[group.by]])
    }

    pseudo_Seu$celltype.con <- paste(pseudo_Seu[[group.by]][, 1],
                                     pseudo_Seu[[group.by.2]][, 1], sep = "_")

  }


  if (Feature %in% rownames(sce_object)) {
    df_Feature <- data.frame(
      group = setNames(sce_object[[group.by]][, group.by],
                       rownames(sce_object[[group.by]])),
      Feature = sce_object[["RNA"]]$data[Feature, ],
      cluster = sce_object[[sample.column]][, 1]
    )
    df_Feature[,Feature] <- expm1(sce_object@assays$RNA$data[Feature,])

  }else{
    df_Feature <- data.frame(
      group = setNames(sce_object[[group.by]][, group.by],
                       rownames(sce_object[[group.by]])),
      Feature = FetchData(sce_object, vars = Feature)[, 1],
      cluster = sce_object[[sample.column]][, 1]
    )

    if (is.null(group.by.2)) {
      aggregated_meta <- sce_object@meta.data %>%
        group_by(!!sym(group.by), !!sym(sample.column)) %>%
        summarise(Feature = mean(!!sym(Feature), na.rm = TRUE))

      aggregated_meta$comb <- paste(aggregated_meta[[group.by]],
                                    aggregated_meta[[sample.column]],
                                    sep = "_")

      #check if the matching agrees with the dash names
      if (!any(is.na(
        aggregated_meta$Feature[
          match(
            pseudo_Seu[[sample.column]][, 1],
            aggregated_meta$comb
          )
        ]
      ))) {
        pseudo_Seu[[Feature]] <- aggregated_meta$Feature[
          match(
            pseudo_Seu[[sample.column]][, 1],
            aggregated_meta$comb
          )
        ]
      } else {
        unified_sep_column <- gsub(
          "-",
          "_",
          pseudo_Seu[[sample.column]][, 1]
        )
        pseudo_Seu[[Feature]] <- aggregated_meta$Feature[
          match(
            unified_sep_column,
            aggregated_meta$comb
          )
        ]
      }


    } else{
      #Compute mean for each orig.ident
      aggregated_meta <- sce_object@meta.data %>%
        group_by(!!sym(group.by), !!sym(sample.column), !!sym(group.by.2)) %>%
        summarise(Feature = mean(!!sym(Feature), na.rm = TRUE))

      aggregated_meta$comb <- paste(aggregated_meta[[group.by]],
                                    aggregated_meta[[group.by.2]],
                                    aggregated_meta[[sample.column]],
                                    sep = "_")

      #cover the case of subsetted to only have one cell type
      if (!length(unique(sce_object@meta.data[[group.by.2]])) == 1) {
        pseudo_Seu[[sample.column]][, 1] <-
          gsub("_", "-", pseudo_Seu[[sample.column]][, 1])
      } else{
        pseudo_Seu[[sample.column]][,1] <-
          paste0(gsub("_",
                      "-",
                      pseudo_Seu[[sample.column]][,1]),
                 "-",
                 unique(sce_object@meta.data[[group.by.2]]))
        pseudo_Seu[[sample.column]][,1] <-
          gsub("_",
               "-",
               pseudo_Seu[[sample.column]][,1])
        rownames(pseudo_Seu@meta.data) <- pseudo_Seu[[sample.column]][,1]
      }
      aggregated_meta$comb <- gsub("_", "-", aggregated_meta$comb)
      pseudo_Seu[[Feature]] <- aggregated_meta$Feature[
        match(
          pseudo_Seu[[sample.column]][, 1],
          aggregated_meta$comb
        )
      ]
    }


  }

  #group results and summarize
  if (is.null(group.by.2)) {
    df_melt <- reshape2::melt(df_Feature)
    df_melt_sum <- df_melt %>%
      dplyr::group_by(group, variable) %>%
      dplyr::summarise(Mean = mean(value))
  } else{
    df_Feature[,{group.by.2}] <- setNames(sce_object[[group.by.2]][,group.by.2],
                                          rownames(sce_object[[group.by.2]]))
    df_melt <- reshape2::melt(df_Feature)
    df_melt_sum <- df_melt %>%
      dplyr::group_by(group, !!sym(group.by.2), variable) %>%
      dplyr::summarise(Mean = mean(value))
  }

  #create comparison list for wilcox, always against control
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    .logger("ListTest empty, comparing every sample with each other")
    group <- unique(sce_object[[group.by]][,group.by])
    #set automatically ctrl condition if not provided
    if (is.null(ctrl.condition)) {
      ctrl.condition <- group[grep(pattern = paste(c("CTRL",
                                                     "Ctrl",
                                                     "WT",
                                                     "Wt",
                                                     "wt"),collapse ="|")
                                   ,group)[1]]
    }


    #create ListTest
    ListTest <- list()
    for (i in seq_along(group)) {
      cndtn <- as.character(group[i])
      if (cndtn != ctrl.condition) {
        ListTest[[length(ListTest) + 1]] <- c(ctrl.condition, cndtn)
      }
    }
  }

  ListTest <- ListTest[!vapply(ListTest, is.null, logical(1))]
  indices <- vapply(ListTest,
                    function(x) match(x[2],
                                      df_melt_sum$group),
                    integer(1))
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

  # Remove vectors with both elements having a mean of 0
  ListTest <- remove_zeros(ListTest, df_melt_sum)

  if (is.null(group.by.2)) {
    group_of_zero <- df_melt %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      dplyr::filter(all_zeros)

    if (nrow(group_of_zero) > 0) {
      warning(
        "Some comparisons have no expression in both groups, setting ",
        "expression to minimum value to ensure test does not fail!"
      )
      df_melt <- df_melt %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(value = dplyr::if_else(
          dplyr::row_number() == 1 &
            all(value == 0),
          .Machine$double.xmin,
          value
        )) %>%
        ungroup()
    }
  } else{

    group_of_zero <- df_melt %>%
      dplyr::group_by(group, !!sym(group.by.2)) %>%
      dplyr::summarise(all_zeros = all(value == 0), .groups = "drop") %>%
      dplyr::filter(all_zeros)

    #check now the result for multiple entries in group.by.2
    groupby2_check <- group_of_zero %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      dplyr::summarise(group_count = dplyr::n_distinct(group),
                       .groups = "drop") %>%
      dplyr::filter(group_count > 1)

    if (nrow(groupby2_check) > 0) {
      warning("Some comparisons have no expression in both groups, setting ",
              "expression to minimum value to ensure test does not fail!")
      df_melt <- df_melt %>%
        dplyr::group_by(group, !!sym(group.by.2)) %>%
        dplyr::mutate(value = dplyr::if_else(row_number() == 1 &
                                               all(value == 0),
                                             .Machine$double.xmin, value)) %>%
        dplyr::ungroup()
    }
  }

  #do statistix with rstatix + stats package
  if (wilcox_test == TRUE & is.null(group.by.2)) {
    stat.test <- df_melt %>%
      dplyr::ungroup() %>%
      rstatix::wilcox_test(value ~ group,
                           comparisons = ListTest,
                           p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p,
                                       method = "bonferroni",
                                       n = length(rownames(sce_object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0,
                              sprintf("%.2e",.Machine$double.xmin),
                              sprintf("%.2e", stat.test$p.adj))
    stat.test$p <- ifelse(stat.test$p == 0,
                          sprintf("%.2e",.Machine$double.xmin),
                          sprintf("%.2e", stat.test$p))

  }

  #do statistix with rstatix + stats package add second group
  if (wilcox_test == TRUE & !is.null(group.by.2)) {
    stat.test <- df_melt %>%
      dplyr::group_by(!!sym(group.by.2)) %>%
      rstatix::wilcox_test(value ~ group,
                           comparisons = ListTest,
                           p.adjust.method = "none") %>%
      rstatix::add_significance()
    stat.test$p.adj <- stats::p.adjust(stat.test$p,
                                       method = "bonferroni",
                                       n = length(rownames(sce_object)))
    stat.test$p.adj <- ifelse(stat.test$p.adj == 0,
                              sprintf("%.2e",.Machine$double.xmin),
                              sprintf("%.2e", stat.test$p.adj))
    stat.test$p <- ifelse(stat.test$p == 0,
                          sprintf("%.2e",.Machine$double.xmin),
                          sprintf("%.2e", stat.test$p))

  }


  #pseudobulk boxplot
  theme_box <- function() {
    theme_bw() +
      theme(
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          linewidth = 1
        ),
        panel.grid.major = element_line(colour = "grey90", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey90", linetype = "dotted"),
        axis.line = element_line(colour = "black"),
        #facet_grid colors
        strip.background = element_rect(
          fill = "lightgrey",
          colour = "black",
          linewidth = 1
        ),
        strip.text = element_text(colour = "black", size = 12),
        # legend.background = element_rect(colour = "grey", fill = "white"),
        # legend.box.background = element_rect(colour = "grey", size = 0.5),
      )
  }

  #TODO there might be a problem with the quantiles, recheck
  #SCpubr does not support outlier removal, therefore identify manually
  if (outlier_removal == TRUE) {

    if (Feature %in% rownames(sce_object)) {
      data_matrix <- pseudo_Seu@assays$RNA$data[Feature,]
    } else{
      data_matrix <- pseudo_Seu[[Feature]]
      data_matrix <- setNames(data_matrix[[Feature]], rownames(data_matrix))
    }

    for (grp2 in unique(pseudo_Seu[[group.by.2]][,1])) {
      for(grp in unique(pseudo_Seu[[group.by]][,1])){
        group_cells <- pseudo_Seu@meta.data[[group.by.2]] == grp2 &
          pseudo_Seu@meta.data[[group.by]] == grp
        subset_mat <- data_matrix[group_cells]

        Q1 <- quantile(subset_mat, 0.25)
        Q3 <- quantile(subset_mat, 0.75)
        IQR <- Q3 - Q1  # interquartile range calculation

        lower_bound <- Q1 - 1.5 * IQR
        upper_bound <- Q3 + 1.5 * IQR

        data_matrix_sub <- ifelse(subset_mat >= lower_bound &
                                    subset_mat <= upper_bound,
                                  subset_mat,
                                  NA)

        if (Feature %in% rownames(sce_object)) {
          pseudo_Seu@assays$RNA$data[Feature, group_cells] <- data_matrix_sub
        } else{
          pseudo_Seu@meta.data[group_cells,Feature] <- data_matrix_sub
        }
      }
    }
  }

  if (!is.null(orderAxis)) {
    pseudo_Seu[[group.by]][, 1] <- factor(pseudo_Seu[[group.by]][, 1],
                                          levels = orderAxis)
  }


  if (is.null(group.by.2)) {
    p <- SCpubr::do_BoxPlot(sample = pseudo_Seu,
                            feature = Feature,
                            group.by = group.by,
                            order = FALSE,
                            boxplot.width = 0.8,
                            legend.position = "right")

  } else {
    p <- SCpubr::do_BoxPlot(sample = pseudo_Seu,
                            feature = Feature,
                            group.by = group.by.2,
                            split.by = group.by,
                            boxplot.width = 0.8,
                            order = FALSE)
  }

  if (plot_sample == TRUE) {
    p <- p + geom_point(size=2, alpha=1, position = position_dodge(width = 0.8))
  }

  p <- p +
    scale_fill_manual(values = rep(vector_colors, 2)) +
    theme_box() +
    theme(
      axis.text.x = element_text(
        color = "black",
        angle = 45,
        hjust = 1,
        size = 16,
        family = "Helvetica"
      ),
      axis.text.y = element_text(
        color = "black",
        size = 16,
        family = "Helvetica"
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_text(
        size = 16,
        family = "Helvetica",
        face = "bold"
      ),
      axis.title = element_text(
        size = 16,
        color = "black",
        family = "Helvetica"
      ),
      plot.title = element_text(
        size = 16,
        hjust = 0.5,
        family = "Helvetica"
      ),
      plot.subtitle = element_text(
        size = 16,
        hjust = 0,
        family = "Helvetica"
      ),
      axis.line = element_line(color = "black"),
      strip.text.x = element_text(
        size = 16,
        color = "black",
        family = "Helvetica"
      ),
      legend.text = element_text(
        size = 14,
        color = "black",
        family = "Helvetica"
      ),
      legend.title = element_text(
        size = 14,
        color = "black",
        family = "Helvetica",
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "bottom"
    )

  #p
  # for only one group
  if (wilcox_test == TRUE & is.null(group.by.2)) {

    if (Feature %in% rownames(sce_object)) {
      p_label <- "p = {p.adj}"
    } else{
      p_label <- "p = {p}"
    }


    if (Feature %in% rownames(sce_object)) {
      stat.test_plot <- stat.test %>%
        mutate(
          y.position = seq(
            from = max(
              pseudo_Seu@assays$RNA$data[Feature, ][
                !is.na(pseudo_Seu@assays$RNA$data[Feature, ])
              ]
            ) * stat_pos_mod,
            by = step_mod,
            length.out = nrow(stat.test)
          )
        )
    } else {
      stat.test_plot <- stat.test %>%
        mutate(
          y.position = seq(
            from = max(
              pseudo_Seu@meta.data[, Feature][
                !is.na(pseudo_Seu@meta.data[, Feature])
              ]
            ) * stat_pos_mod,
            by = step_mod,
            length.out = nrow(stat.test)
          )
        )
    }


    p <- p + stat_pvalue_manual(stat.test_plot,
                                label = "p = {p.adj}",
                                y.position = "y.position",
                                # step.increase = 0.2,
                                inherit.aes = FALSE,
                                size = size.wilcox,
                                angle= 0,
                                hjust= hjust.wilcox,
                                vjust = vjust.wilcox)
  }

  if (wilcox_test == TRUE & !is.null(group.by.2)) {


    if (Feature %in% rownames(sce_object)) {
      p_label <- "p = {p.adj}"
    } else{
      p_label <- "p = {p}"
    }


    if (Feature %in% rownames(sce_object)) {
      stat.test_plot <- stat.test %>%
        mutate(
          y.position = seq(
            from = max(
              pseudo_Seu@assays$RNA$data[Feature, ][
                !is.na(pseudo_Seu@assays$RNA$data[Feature, ])
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
              pseudo_Seu@meta.data[, Feature][
                !is.na(pseudo_Seu@meta.data[, Feature])
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
                                angle= 0,
                                hjust= hjust.wilcox.2,
                                vjust = vjust.wilcox.2,
                                tip.length = 0.02,
                                bracket.size = sign_bar)
  }
  return(p)
}
