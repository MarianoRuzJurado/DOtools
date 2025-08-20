#' @author Mariano Ruz Jurado
#' @title DO Bar plot for GSEA df result
#' @description This function generates a split barplot. This is a plot where the top 10 Go terms
#'are shown, sorted based on a column ('col_split'). Two conditions are shown at the same
#'time. One condition is shown in the positive axis, while the other in the negative one.
#'The condition to be shown as positive is set with 'pos_col'.
#'
#'The GO terms will be shown inside the bars, if the term is too long, using 'cutoff',
#'you can control the maximum number of characters per line.
#'
#' Pre-filter of the dataframe to contain significant Terms is recommended
#'
#' @param df_GSEA dataframe with the results of a gene set enrichment analysis
#' @param term_col column in the dataframe that contains the terms
#' @param col_split column in the dataframe that will be used to sort and split the plot
#' @param cond_col column in the dataframe that contains the condition information
#' @param pos_cond condition that will be shown in the positive side of the plot
#' @param cutoff maximum number of characters per line
#' @param log10_transform if col_split contains values between 0 and 1, assume they are pvals and apply a -log10 transformation
#' @param figsize figure size
#' @param topN how many terms are shown
#' @param path path to save the plot
#' @param filename filename for the plot
#' @param spacing space to add between bars and origin. It is a percentage value, indicating that the bars start at 5 % of the maximum X axis value.
#' @param txt_size size of the go terms text
#' @param alpha_colors alpha value for the colors of the bars
#' @param colors_pairs colors for each condition (1st color --> negative axis; 2nd color --> positive axis)
#' @param title title of the plot
#' @param showP if False, the axis is return
#' @param celltype vector with cell types you want to subset for, use "all" for all celltypes contained in the dataframe column "celltype"
#'
#' @return: None or the axis
#'
#' @importFrom basilisk basiliskRun
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
#' result_GO$celltype <- "CM1"
#'
#' #Run SplitBarGSEA visualisation
#' DO.SplitBarGSEA(df_GSEA = result_GO,
#'                 term_col = "Term",
#'                 col_split = "Combined.Score",
#'                 cond_col = "State",
#'                 pos_cond = "enriched",
#'                 cutoff = 40,
#'                 log10_transform = TRUE,
#'                 figsize = c(12, 8),
#'                 topN = 10,
#'                 colors_pairs = c("sandybrown", "royalblue"),
#'                 alpha_colors = 0.3,
#'                 path = NULL,
#'                 spacing = 5,
#'                 txt_size = 12,
#'                 filename = "SplitBar.svg",
#'                 title = "Top 10 GO Terms in each Condition: ",
#'                 showP = FALSE,
#'                 celltype = "all")
#'
#' @export
DO.SplitBarGSEA <- function(df_GSEA,
                            term_col,
                            col_split,
                            cond_col,
                            pos_cond,
                            cutoff=40,
                            log10_transform=TRUE,
                            figsize=c(12,8),
                            topN=10,
                            colors_pairs=c("sandybrown","royalblue"),
                            alpha_colors=0.3,
                            path=NULL,
                            spacing=5,
                            txt_size=12,
                            filename="SplitBar.svg",
                            title="Top 10 GO Terms in each Condition: ",
                            showP=FALSE,
                            celltype="all")
{

  #Create argument list to pass
  args <- list(
    df_GSEA = df_GSEA,
    term_col = term_col,
    col_split = col_split,
    cond_col = cond_col,
    pos_cond = pos_cond,
    title = paste0(title, celltype),
    cutoff = cutoff,
    log10_transform = log10_transform,
    figsize = figsize,
    topN = topN,
    colors_pairs = colors_pairs,
    alpha_colors = alpha_colors,
    path = path,
    spacing = spacing,
    txt_size = txt_size,
    filename = filename,
    showP = showP
  )

  if (!"celltype" %in% colnames(df_GSEA)) {
    stop("Provided data frame has no column named 'celltype'! Please add a column with cell type information.")
  }

  #source PATH to python script in install folder
  path_py <- system.file("python", "GSEA_split_bar.py", package = "DOtools")

  #Subset by given celltype names
  if (!celltype=="all") {
    df_GSEA <- df_GSEA[df_GSEA$celltype %in% celltype,]
  }

  #basilisk implementation
  results <- basilisk::basiliskRun(env = DOtoolsEnv, fun=function(args){

    reticulate::source_python(path_py)
    #Initialize matplot package
    plt <- reticulate::import("matplotlib.pyplot")

    #Run for each celltype the split_bar_gsea function
    for (celltype in unique(df_GSEA$celltype)) {
      df_GSEA_sub <- df_GSEA[df_GSEA$celltype == celltype,]
      df_GSEA_sub$Term <- sapply(strsplit(df_GSEA_sub$Term, " \\(GO"), function(x)x[1])

      #Conversion of the data.frame
      df_GSEA_sub_pd <- reticulate::r_to_py(df_GSEA_sub)

      #run python function with given arguments
      plot <- split_bar_gsea(df = args$df_GSEA,
                             term_col = args$term_col,
                             col_split = args$col_split,
                             cond_col = args$cond_col,
                             pos_cond = args$pos_cond,
                             title = args$title,
                             cutoff = args$cutoff,
                             log10_transform = args$log10_transform,
                             figsize = args$figsize,
                             topN = args$topN,
                             colors_pairs = args$colors_pairs,
                             alpha_colors = args$alpha_colors,
                             path = args$path,
                             spacing = args$spacing,
                             txt_size = args$txt_size,
                             filename = args$filename,
                             showP = args$showP)

      if (showP == TRUE) {
        #x-title settings
        # plot$set_xlabel("Combined Score")
        plt$show()
      }

      if (!is.null(path)) {
        #Save under provided PATH
        plt$savefig(paste0(path, celltype, "_", filename), dpi=300, bbox_inches="tight")
      }

    }

  }, args = args)
}
