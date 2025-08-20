#' @title DO.enrichR
#' @description
#' Performs Gene Ontology enrichment analysis on differentially expressed genes using the EnrichR API.
#' Separately analyzes upregulated and downregulated genes and returns results.
#'
#' @param df_DGE data.frame containing differential gene expression results.
#' @param gene_column column name in `df` with gene symbols.
#' @param pval_column column name in `df` with p-values.
#' @param log2fc_column column name in `df` with log2 fold changes.
#' @param pval_cutoff adjusted p-value threshold for significance (default = 0.05).
#' @param log2fc_cutoff log2 fold change threshold for up/down regulation (default = 0.25).
#' @param path folder path where the output Excel file will be saved. A subfolder "GSA_Tables" will be created.
#' @param filename suffix used in the Excel filename (e.g., "GSA_CellType_MyAnalysis.xlsx").
#' @param species species name for enrichment analysis. Options include "Human", "Mouse", "Yeast", etc. (default = "Mouse").
#' @param go_catgs GO databases to use. Defaults to c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023").
#'
#' @return data.frame with GO enrichment results if `path` is NULL, otherwise writes an Excel file.
#'
#'
#' @import enrichR
#' @import openxlsx
#' @import dplyr
#'
#' @examples
#' library(enrichR)
#'
#' sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#' DGE_result <- DO.MultiDGE(sce_data,
#'                          sample_col = "orig.ident",
#'                          method_sc = "wilcox",
#'                          annotation_col = "annotation",
#'                          ident_ctrl = "healthy")
#'
#' DGE_result <- DGE_result[DGE_result$celltype == "CD4_T_cells",]
#'
#' result_GO <- DO.enrichR(df_DGE = DGE_result,
#'                        gene_column = "gene",
#'                        pval_column = "p_val_SC_wilcox",
#'                        log2fc_column = "avg_log2FC_SC_wilcox",
#'                        pval_cutoff = 0.05,
#'                        log2fc_cutoff = 0.25,
#'                        path = NULL,
#'                        filename = "",
#'                        species = "Human",
#'                        go_catgs = "GO_Biological_Process_2023")
#'
#' @export
DO.enrichR <- function(df_DGE,
                       gene_column,
                       pval_column,
                       log2fc_column,
                       pval_cutoff = 0.05,
                       log2fc_cutoff = 0.25,
                       path = NULL,
                       filename = '',
                       species = 'Human',
                       go_catgs = c('GO_Molecular_Function_2023',
                                    'GO_Cellular_Component_2023',
                                    'GO_Biological_Process_2023')) {


  # Subset up/down-regulated
  df_up <- df_DGE[df_DGE[[pval_column]] < pval_cutoff & df_DGE[[log2fc_column]] > log2fc_cutoff, ]
  df_down <- df_DGE[df_DGE[[pval_column]] < pval_cutoff & df_DGE[[log2fc_column]] < -log2fc_cutoff, ]

  # Set enrichR organism
  enrichR::setEnrichrSite("Enrichr")  # explicitly connect to Enrichr API

  # Run enrichment
  res_up <- enrichR::enrichr(genes = df_up[[gene_column]], databases = go_catgs)
  res_down <- enrichR::enrichr(genes = df_down[[gene_column]], databases = go_catgs)

  # Combine results
  enrich_results <- function(result_list, label) {
    do.call(rbind, lapply(names(result_list), function(db) {
      df <- result_list[[db]]
      df$Database <- db
      df$State <- label
      return(df)
    }))
  }

  combined_up <- enrich_results(res_up, "enriched")
  combined_down <- enrich_results(res_down, "depleted")

  combined_res <- rbind(combined_up, combined_down)

  # Write to Excel or return
  if (!is.null(path)) {
    dir.create(file.path(path, "GO_Tables"), showWarnings = FALSE, recursive = TRUE)
    out_file <- file.path(path, "GO_Tables", paste0("GO_", filename, ".xlsx"))
    openxlsx::write.xlsx(combined_res, out_file, rowNames = FALSE)
    .logger(paste0("Write results to: ", out_file))
    return(invisible(NULL))
  } else {
    return(combined_res)
  }
}
