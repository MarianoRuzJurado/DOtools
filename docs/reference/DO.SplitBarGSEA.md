# DO Bar plot for GSEA df result

This function generates a split barplot. This is a plot where the top 10
Go terms are shown, sorted based on a column ('col_split'). Two
conditions are shown at the same time. One condition is shown in the
positive axis, while the other in the negative one. The condition to be
shown as positive is set with 'pos_col'.

The GO terms will be shown inside the bars, if the term is too long,
using 'cutoff', you can control the maximum number of characters per
line.

Pre-filter of the dataframe to contain significant Terms is recommended

## Usage

``` r
DO.SplitBarGSEA(
  df_GSEA,
  term_col,
  col_split,
  cond_col,
  pos_cond,
  cutoff = 40,
  log10_transform = TRUE,
  figsize = c(12, 8),
  topN = 10,
  colors_pairs = c("sandybrown", "royalblue"),
  alpha_colors = 0.3,
  path = NULL,
  spacing = 5,
  txt_size = 12,
  filename = "SplitBar.svg",
  title = "Top 10 GO Terms in each Condition: ",
  showP = FALSE,
  celltype = "all"
)
```

## Arguments

- df_GSEA:

  dataframe with the results of a gene set enrichment analysis

- term_col:

  column in the dataframe that contains the terms

- col_split:

  column in the dataframe that will be used to sort and split the plot

- cond_col:

  column in the dataframe that contains the condition information

- pos_cond:

  condition that will be shown in the positive side of the plot

- cutoff:

  maximum number of characters per line

- log10_transform:

  if col_split contains values between 0 and 1, assume they are pvals
  and apply a -log10 transformation

- figsize:

  figure size

- topN:

  how many terms are shown

- colors_pairs:

  colors for each condition (1st color –\> negative axis; 2nd color –\>
  positive axis)

- alpha_colors:

  alpha value for the colors of the bars

- path:

  path to save the plot

- spacing:

  space to add between bars and origin. It is a percentage value ,
  indicating that the bars start at 5 % of the maximum X axis value.

- txt_size:

  size of the go terms text

- filename:

  filename for the plot

- title:

  title of the plot

- showP:

  if False, the axis is return

- celltype:

  vector with cell types you want to subset for, use "all" for all
  celltypes contained in the dataframe column "celltype"

## Value

: None or the axis

## Author

Mariano Ruz Jurado

## Examples

``` r
library(enrichR)
#> Welcome to enrichR
#> Checking connections ... 
#> Enrichr ... 
#> Connection is Live!
#> FlyEnrichr ... 
#> Connection is Live!
#> WormEnrichr ... 
#> Connection is Live!
#> YeastEnrichr ... 
#> Connection is Live!
#> FishEnrichr ... 
#> Connection is Live!
#> OxEnrichr ... 
#> Connection is Live!

sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
DGE_result <- DO.MultiDGE(sce_data,
    sample_col = "orig.ident",
    method_sc = "wilcox",
    annotation_col = "annotation",
    ident_ctrl = "healthy"
)
#> Centering and scaling data matrix
#> 2026-01-14 13:00:27 - Corrected annotation names in pseudo-bulk object by replacing '-' with '_'.
#> 2026-01-14 13:00:27 - Starting DGE single cell method analysis
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: Monocytes
#> 2026-01-14 13:00:27 - Skipping Monocytes since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: CD4_T_cells
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: NK
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: CD8_T_cells
#> 2026-01-14 13:00:27 - Skipping CD8_T_cells since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: B_cells
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: pDC
#> 2026-01-14 13:00:27 - Finished DGE single cell method analysis
#> 2026-01-14 13:00:27 - Starting DGE pseudo bulk method analysis
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: B_cells
#> 2026-01-14 13:00:27 - Skipping B_cells since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: CD4_T_cells
#> 2026-01-14 13:00:27 - Skipping CD4_T_cells since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: Monocytes
#> 2026-01-14 13:00:27 - Skipping Monocytes since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: NK
#> 2026-01-14 13:00:27 - Skipping NK since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: pDC
#> 2026-01-14 13:00:27 - Skipping pDC since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:27 - Comparing disease with healthy in: CD8_T_cells
#> 2026-01-14 13:00:27 - Skipping CD8_T_cells since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:27 - Finished DGE pseudo bulk method analysis
#> 2026-01-14 13:00:27 - DGE pseudo bulk result is empty...
DGE_result <- DGE_result[DGE_result$celltype == "CD4_T_cells", ]

result_GO <- DO.enrichR(
    df_DGE = DGE_result,
    gene_column = "gene",
    pval_column = "p_val_SC_wilcox",
    log2fc_column = "avg_log2FC_SC_wilcox",
    pval_cutoff = 0.05,
    log2fc_cutoff = 0.25,
    path = NULL,
    filename = "",
    species = "Human",
    go_catgs = "GO_Biological_Process_2023"
)
#> Connection changed to https://maayanlab.cloud/Enrichr/
#> Connection is Live!
#> Uploading data to Enrichr... Done.
#>   Querying GO_Biological_Process_2023... Done.
#> Parsing results... Done.
#> Uploading data to Enrichr... Done.
#>   Querying GO_Biological_Process_2023... Done.
#> Parsing results... Done.

result_GO$celltype <- "CM1"

# Run SplitBarGSEA visualisation
DO.SplitBarGSEA(
    df_GSEA = result_GO,
    term_col = "Term",
    col_split = "Combined.Score",
    cond_col = "State",
    pos_cond = "enriched",
    cutoff = 40,
    log10_transform = TRUE,
    figsize = c(12, 8),
    topN = 10,
    colors_pairs = c("sandybrown", "royalblue"),
    alpha_colors = 0.3,
    path = NULL,
    spacing = 5,
    txt_size = 12,
    filename = "SplitBar.svg",
    title = "Top 10 GO Terms in each Condition: ",
    showP = FALSE,
    celltype = "all"
)
```
