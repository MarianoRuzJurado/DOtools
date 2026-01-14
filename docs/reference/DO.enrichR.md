# DO.enrichR

Performs Gene Ontology enrichment analysis on differentially expressed
genes using the EnrichR API. Separately analyzes upregulated and
downregulated genes and returns results.

## Usage

``` r
DO.enrichR(
  df_DGE,
  gene_column,
  pval_column,
  log2fc_column,
  pval_cutoff = 0.05,
  log2fc_cutoff = 0.25,
  path = NULL,
  filename = "",
  species = "Human",
  go_catgs = c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023",
    "GO_Biological_Process_2023")
)
```

## Arguments

- df_DGE:

  data.frame containing differential gene expression results.

- gene_column:

  column name in `df` with gene symbols.

- pval_column:

  column name in `df` with p-values.

- log2fc_column:

  column name in `df` with log2 fold changes.

- pval_cutoff:

  adjusted p-value threshold for significance (default = 0.05).

- log2fc_cutoff:

  log2 fold change threshold for up/down regulation (default = 0.25).

- path:

  folder path where the output Excel file will be saved. A subfolder
  "GSA_Tables" will be created.

- filename:

  suffix used in the Excel filename (e.g.,
  "GSA_CellType_MyAnalysis.xlsx").

- species:

  species name for enrichment analysis. Options include "Human",
  "Mouse", "Yeast", etc. (default = "Mouse").

- go_catgs:

  GO databases to use. Defaults to c(GO_Biological_Process_2023").

## Value

data.frame with GO enrichment results if `path` is NULL, otherwise
writes an Excel file.

## Examples

``` r
library(enrichR)

sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
DGE_result <- DO.MultiDGE(sce_data,
    sample_col = "orig.ident",
    method_sc = "wilcox",
    annotation_col = "annotation",
    ident_ctrl = "healthy"
)
#> Centering and scaling data matrix
#> 2026-01-14 13:00:43 - Corrected annotation names in pseudo-bulk object by replacing '-' with '_'.
#> 2026-01-14 13:00:43 - Starting DGE single cell method analysis
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: Monocytes
#> 2026-01-14 13:00:43 - Skipping Monocytes since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: CD4_T_cells
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: NK
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: CD8_T_cells
#> 2026-01-14 13:00:43 - Skipping CD8_T_cells since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: B_cells
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: pDC
#> 2026-01-14 13:00:43 - Finished DGE single cell method analysis
#> 2026-01-14 13:00:43 - Starting DGE pseudo bulk method analysis
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: B_cells
#> 2026-01-14 13:00:43 - Skipping B_cells since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: CD4_T_cells
#> 2026-01-14 13:00:43 - Skipping CD4_T_cells since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: Monocytes
#> 2026-01-14 13:00:43 - Skipping Monocytes since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: NK
#> 2026-01-14 13:00:43 - Skipping NK since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: pDC
#> 2026-01-14 13:00:43 - Skipping pDC since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:43 - Comparing disease with healthy in: CD8_T_cells
#> 2026-01-14 13:00:43 - Skipping CD8_T_cells since one comparison has fewer than 3 cells!
#> 2026-01-14 13:00:43 - Finished DGE pseudo bulk method analysis
#> 2026-01-14 13:00:43 - DGE pseudo bulk result is empty...

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
```
