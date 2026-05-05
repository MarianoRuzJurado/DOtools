#' @author Mariano Ruz Jurado
#' @title DO.MultiDGE
#' @description
#' Performs differential gene expression analysis using both single-cell and
#' pseudo-bulk approaches across all annotated cell types. The single-cell
#' method uses Seurat's `FindMarkers`, while pseudo-bulk testing uses `DESeq2`
#' on aggregated expression profiles. Outputs a merged data frame with DGE
#' statistics from both methods per condition and cell type.
#'
#' @param sce_object The seurat or SCE object
#' @param assay Specified assay in Seurat or SCE object, default "RNA"
#' @param method_sc method to use for single cell DEG analysis, see FindMarkers
#' from Seurat for options, default "wilcox"
#' @param method_pb method to use for pseudobulk DEG analysis, currently
#' supports DESeq2 implemented in Seurat and glmGamPoi
#' @param group_by Column in meta data containing groups used for testing,
#' default "condition"
#' @param annotation_col Column in meta data containing information of cell type
#'  annotation
#' @param sample_col Column in meta data containing information of sample
#' annotation, default "orig.ident"
#' @param ident_ctrl Name of the condition in group_by to test against as ctrl,
#' default "ctrl"
#' @param min_pct only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations, default is 0
#' @param logfc_threshold Limit testing to genes which show, on average, at
#' least X-fold difference (log-scale) between the two groups of cells, default
#' is 0.
#' @param only_pos Only return positive markers, default FALSE
#' @param min_cells_group Minimum number of cells in one of the groups,
#' default 3
#' @param design_fit_glm Design for fitting the glmGamPoi model
#' @param ... Additional arguments passed to FindMarkers function
#'
#' @import Seurat
#' @import DESeq2
#' @import tibble
#' @import dplyr
#'
#' @return Dataframe containing statistics for each gene from the single cell
#' and the Pseudobulk DGE approach.
#'
#' @examples
#'
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#' DGE_result <- DO.MultiDGE(sce_data,
#'     sample_col = "orig.ident",
#'     method_sc = "wilcox",
#'     annotation_col = "annotation",
#'     ident_ctrl = "healthy"
#' )
#'
#' @export
DO.MultiDGE <- function(sce_object,
    assay = "RNA",
    method_sc = "wilcox",
    method_pb = "DESeq2",
    group_by = "condition",
    annotation_col = "annotation",
    sample_col = "orig.ident",
    ident_ctrl = "ctrl",
    min_pct = 0,
    logfc_threshold = 0,
    only_pos = FALSE,
    min_cells_group = 3,
    design_fit_glm = NULL,
    ...) {
    # support for single cell experiment objects
    if (methods::is(sce_object, "SingleCellExperiment")) {
        class_obj <- "SingleCellExperiment"
        sce_object <- .suppressDeprecationWarnings(as.Seurat(sce_object))
    }

    if (!ident_ctrl %in% sce_object@meta.data[[group_by]]) {
        stop(sprintf(
            "%s was not found in meta data under the specified ",
            "group_by column: %s",
            ident_ctrl,
            group_by
        ))
    }

    # data.frame for collecting sc method results
    DEG_stats_collector_sc <- data.frame()
    # data.frame for collecting pseudo_bulk method results
    DEG_stats_collector_pb <- data.frame()

    # Create a PseudoBulk representation of the Seurat single cell expression
    sce_object_PB <- .suppressDeprecationWarnings(
        AggregateExpression(sce_object,
            assays = assay,
            return.seurat = TRUE,
            group.by = c(group_by, sample_col, annotation_col)
        )
    )
    # internally AggregateExpression uses interaction() this replaces _ with -
    original_annots <- unique(sce_object@meta.data[[annotation_col]])
    pb_annots <- unique(sce_object_PB@meta.data[[annotation_col]])

    # are all pseudo-bulk annotations in the original set
    if (!all(pb_annots %in% original_annots)) {
        # if it is a hyphen exchange this should fix it
        converted_annots <- gsub("-", "_", pb_annots)

        if (all(converted_annots %in% original_annots)) {
            sce_object_PB@meta.data[[annotation_col]] <-
                gsub(
                    "-",
                    "_",
                    sce_object_PB@meta.data[[annotation_col]]
                )
            .logger(
                paste0(
                    "Corrected annotation names in pseudo-bulk object ",
                    "by replacing '-' with '_'."
                )
            )
        } else {
            .logger(paste0(
                "Annotation names do not match even after ",
                "replacing '-' with '_'. Please inspect manually."
            ))
        }
    } else {
        .logger(paste0(
            "Annotation names are consistent between original and ",
            "pseudo-bulk objects."
        ))
    }

    PB_ident <- paste(annotation_col, group_by, sep = "_")
    sce_object_PB@meta.data[[PB_ident]] <- paste(
        sce_object_PB@meta.data[[annotation_col]],
        sce_object_PB@meta.data[[group_by]],
        sep = "_"
    )

    Idents(sce_object_PB) <- PB_ident

    if (method_sc != "none") {
      .logger("Starting DGE single cell method analysis")
      for (ident_con in unique(sce_object@meta.data[[group_by]])) {
        # .logger(ident_con)
        if (ident_con != ident_ctrl) { # skip the ctrl ident
          for (celltype in unique(sce_object@meta.data[[annotation_col]])) {
            .logger(
              paste0(
                "Comparing ",
                ident_con,
                " with ",
                ident_ctrl,
                " in: ",
                celltype
              )
            )

            Seu_celltype <- subset(sce_object,
              subset = !!sym(annotation_col) == celltype
            )
            # Check if there are groups with less than 3 cells
            table_cells_sc <- table(Seu_celltype@meta.data[[group_by]])

            count_1_sc <- table_cells_sc[ident_con]
            count_2_sc <- table_cells_sc[ident_ctrl]

            if (is.na(count_1_sc)) {
              count_1_sc <- 0
            }

            if (is.na(count_2_sc)) {
              count_2_sc <- 0
            }

            if (count_1_sc >= 3 && count_2_sc >= 3) {
              DEG_stats_sc <- .suppressDeprecationWarnings(
                FindMarkers(
                  object = Seu_celltype,
                  ident.1 = ident_con,
                  ident.2 = ident_ctrl,
                  min.pct = min_pct,
                  logfc.threshold = logfc_threshold,
                  assay = assay,
                  only.pos = only_pos,
                  test.use = method_sc,
                  min.cells.group = min_cells_group,
                  group.by = group_by,
                  ...
                )
              )
              ###TODO Add the FDR correction
              DEG_stats_sc <- rownames_to_column(DEG_stats_sc,
                var = "gene"
              )
              DEG_stats_sc[["condition"]] <- ident_con
              DEG_stats_sc[["celltype"]] <- celltype
              # DEG_stats_sc[["testmethod"]] <- method_sc
              DEG_stats_collector_sc <- rbind(
                DEG_stats_collector_sc,
                DEG_stats_sc
              )
            } else {
              .logger(paste0(
                "Skipping ",
                celltype,
                " since one comparison has fewer than 3 cells!"
              ))
            }
          }
        }
      }

      .logger("Finished DGE single cell method analysis")
    } else{
      .logger("Skipping DGE single cell method analysis!")
    }

    if (method_pb == "DESeq2") {
      .logger("Starting DGE pseudo bulk method analysis")
      for (celltype in unique(sce_object_PB@meta.data[[annotation_col]])) {
        # .logger(celltype)
        # running find markers condition specific
        for (ident_con in unique(sce_object_PB@meta.data[[group_by]])) {
          if (ident_con != ident_ctrl) {
            .logger(
              paste0(
                "Comparing ", ident_con, " with ",
                ident_ctrl, " in: ", celltype
              )
            )
            ident_1 <- paste0(c(celltype, ident_con), collapse = "_")
            ident_2 <- paste0(c(celltype, ident_ctrl), collapse = "_")

            # Check if any of the groups have fewer than 3 cells
            table_cells_pb <- table(Idents(sce_object_PB))

            count_1_pb <- table_cells_pb[ident_1]
            count_2_pb <- table_cells_pb[ident_2]

            if (is.na(count_1_pb)) {
              count_1_pb <- 0
            }

            if (is.na(count_2_pb)) {
              count_2_pb <- 0
            }


            if (count_1_pb >= 3 && count_2_pb >= 3) {
              DEG_stats_pb <- .suppressDeprecationWarnings(
                FindMarkers(
                  object = sce_object_PB,
                  ident.1 = ident_1,
                  ident.2 = ident_2,
                  min.pct = min_pct,
                  logfc.threshold = logfc_threshold,
                  assay = assay,
                  only.pos = only_pos,
                  test.use = "DESeq2",
                  min.cells.group = min_cells_group,
                  group.by = PB_ident,
                  ...
                )
              )

              DEG_stats_pb <- rownames_to_column(DEG_stats_pb,
                var = "gene"
              )
              DEG_stats_pb[["condition"]] <- ident_con
              DEG_stats_pb[["celltype"]] <- celltype
              # DEG_stats_pb[["testmethod"]] <- "DESeq2"
              DEG_stats_collector_pb <- rbind(
                DEG_stats_collector_pb,
                DEG_stats_pb
              )
            } else {
              .logger(paste0(
                "Skipping ", celltype, " since one comparison has ",
                "fewer than 3 cells!"
              ))
            }
          }
        }
      }
      .logger("Finished DGE pseudo bulk method analysis")
    }

    if (method_pb == "glmGamPoi") {
      .logger("Starting DGE pseudo bulk method analysis")
      # support for Seurat objects
      if (methods::is(sce_object, "Seurat")) {
        DefaultAssay(sce_object) <- assay
        sce_object <- .suppressDeprecationWarnings(
          Seurat::as.SingleCellExperiment(sce_object,
            assay = assay
          )
        )
      }

      # vector for grouping the single cell data for
      group_meta <- c(group_by, sample_col, annotation_col)

      #Check if all supplied columns exist in metadata
      if(length(setdiff(group_meta, names(SummarizedExperiment::colData(sce_object)))) > 0){
        stop("Not all group_by arguments are found in metadata!")
      }

      #catch for default design
      if (is.null(design_fit_glm)) {
        design_fit_glm <- ~ annotation + condition + condition:annotation - 1
      }

      #PB with specified columns
      sce_object_pb <- glmGamPoi::pseudobulk(
        sce_object,
        group_by = glmGamPoi::vars(!!!rlang::syms(group_meta))
      )

      #checking if some cell types are not present in all the conditions
      tab <- table(sce_object_pb[[annotation_col]], sce_object_pb[[group_by]])
      bad_annotations <- rownames(tab)[rowSums(tab > 0) < ncol(tab)]

      if (length(bad_annotations) > 0) {
        .logger(paste0("Removing ", length(bad_annotations),
         " annotation level(s) not present in all conditions: ",
            paste(bad_annotations, collapse = ", "))
        )
      }

      #keep only valid ones
      keep_annotations <- rownames(tab)[rowSums(tab > 0) == ncol(tab)]
      sce_object_pb <- sce_object_pb[, sce_object_pb[[annotation_col]] %in% keep_annotations]

      .logger("Fitting Gamma-Poisson model...")
      #fit Gamma-Poisson model
      fit <- glmGamPoi::glm_gp(
        sce_object_pb, design = design_fit_glm
      )

      comp <- setdiff(unique(sce_object[[group_by]]), ident_ctrl)

      # loop over comparisons and cell types
      DEG_stats_collector_pb <- data.frame()
      for (grp in comp) {
        if (!is.null(annotation_col)){
          for (celltype in keep_annotations) {

            #build the contrast string
            contrast <- paste0(
              "cond(", annotation_col, "='", celltype, "', ",
              group_by, "='", grp, "') - ",
              "cond(", annotation_col, "='", celltype, "', ",
              group_by, "='", ident_ctrl, "')"
            )

            de_res <- glmGamPoi::test_de(fit, contrast = contrast)
            de_res[["celltype"]] <- celltype
            de_res[["condition"]] <- grp

            DEG_stats_collector_pb <- rbind(DEG_stats_collector_pb, de_res)
          }
        } else{
          #build the contrast string
          contrast <- paste0(
            "cond(", group_by, "='", grp, "') - ",
            "cond(", group_by, "='", ident_ctrl, "')"
          )

          de_res <- glmGamPoi::test_de(fit, contrast = contrast)
          de_res[["condition"]] <- grp

          DEG_stats_collector_pb <- rbind(DEG_stats_collector_pb, de_res)
          colnames(DEG_stats_collector_pb) <- c(
            "gene", "p_val", "p_val_adj", colnames(DEG_stats_collector_pb)[4:6],
            "avg_log2FC", colnames(DEG_stats_collector_pb)[8:length(
              colnames(DEG_stats_collector_pb)
              )]
          )
        }

      }
      .logger("Finished DGE pseudo bulk method analysis")
    }

    if (method_pb == "none") {
      .logger("Skipping DGE pseudo-bulk method analysis!")
    }

    if (!length(DEG_stats_collector_pb) == 0) {
        # combine the two df
        df_pb <- DEG_stats_collector_pb %>%
            rename(
                p_val_PB_DESeq2 = p_val,
                p_val_adj_PB_DESeq2 = p_val_adj,
                avg_log2FC_PB_DESeq2 = avg_log2FC
            ) %>%
            select(
                gene,
                celltype,
                condition,
                p_val_PB_DESeq2,
                p_val_adj_PB_DESeq2,
                avg_log2FC_PB_DESeq2
            )
    } else {
        .logger("DGE pseudo bulk result is empty...")
        df_pb <- data.frame(gene = NA, celltype = NA, condition = NA)
    }


    rename_map <- rlang::set_names(
        c("p_val", "p_val_adj", "avg_log2FC"),
        paste0(c("p_val_SC_", "p_val_adj_SC_", "avg_log2FC_SC_"), method_sc)
    )

    df_sc <- DEG_stats_collector_sc %>%
        rename(!!!rename_map) # !!! to unquote the strings

    # Sorting
    first_cols <- c("gene", "pct.1", "pct.2", "celltype", "condition")

    merged_df <- df_sc %>%
        left_join(df_pb, by = c("gene", "celltype", "condition")) %>%
        select(all_of(first_cols), everything()) %>%
        select(all_of(first_cols), sort(setdiff(names(.), first_cols)))

    return(merged_df)
}
