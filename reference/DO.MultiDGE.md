# DO.MultiDGE

Performs differential gene expression analysis using both single-cell
and pseudo-bulk approaches across all annotated cell types. The
single-cell method uses Seurat's `FindMarkers`, while pseudo-bulk
testing uses `DESeq2` on aggregated expression profiles. Outputs a
merged data frame with DGE statistics from both methods per condition
and cell type.

## Usage

``` r
DO.MultiDGE(
  sce_object,
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
  group_glm = NULL,
  ...
)
```

## Arguments

- sce_object:

  The seurat or SCE object

- assay:

  Specified assay in Seurat or SCE object, default "RNA"

- method_sc:

  method to use for single cell DEG analysis, see FindMarkers from
  Seurat for options, default "wilcox", set "none" for skip

- method_pb:

  method to use for pseudobulk DEG analysis, currently supports DESeq2
  implemented in Seurat and glmGamPoi, set "none" for skip

- group_by:

  Column in meta data containing groups used for testing, default
  "condition"

- annotation_col:

  Column in meta data containing information of cell type annotation,
  set to NULL for all cell type analysis

- sample_col:

  Column in meta data containing information of sample annotation,
  default "orig.ident"

- ident_ctrl:

  Name of the condition in group_by to test against as ctrl, default
  "ctrl"

- min_pct:

  only test genes that are detected in a minimum fraction of min.pct
  cells in either of the two populations, default is 0

- logfc_threshold:

  Limit testing to genes which show, on average, at least X-fold
  difference (log-scale) between the two groups of cells, default is 0.

- only_pos:

  Only return positive markers, default FALSE

- min_cells_group:

  Minimum number of cells in one of the groups, default 3

- design_fit_glm:

  Design for fitting the glmGamPoi model

- group_glm:

  grouping vector for glmGamPoi, can have covariates defined

- ...:

  Additional arguments passed to FindMarkers function

## Value

Dataframe containing statistics for each gene from the single cell and
the Pseudobulk DGE approach.

## Author

Mariano Ruz Jurado

## Examples

``` r

sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
DGE_result <- DO.MultiDGE(sce_data,
    sample_col = "orig.ident",
    method_sc = "wilcox",
    annotation_col = "annotation",
    ident_ctrl = "healthy"
)
#> Names of identity class contain underscores ('_'), replacing with dashes ('-')
#> This message is displayed once every 8 hours.
#> Centering and scaling data matrix
#> 2026-07-20 10:10:24 - Corrected annotation names in pseudo-bulk object by replacing '-' with '_'.
#> 2026-07-20 10:10:24 - Starting DGE single cell method analysis
#> 2026-07-20 10:10:24 - Comparing disease with healthy in: Monocytes
#> 2026-07-20 10:10:24 - Skipping Monocytes since one comparison has fewer than 3 cells!
#> 2026-07-20 10:10:24 - Comparing disease with healthy in: CD4_T_cells
#> 2026-07-20 10:10:26 - Comparing disease with healthy in: NK
#> 2026-07-20 10:10:26 - Comparing disease with healthy in: CD8_T_cells
#> 2026-07-20 10:10:26 - Skipping CD8_T_cells since one comparison has fewer than 3 cells!
#> 2026-07-20 10:10:26 - Comparing disease with healthy in: B_cells
#> 2026-07-20 10:10:27 - Comparing disease with healthy in: pDC
#> 2026-07-20 10:10:27 - Finished DGE single cell method analysis
#> 2026-07-20 10:10:27 - Starting DGE pseudo bulk method analysis
#> 2026-07-20 10:10:27 - Comparing disease with healthy in: B_cells
#> 2026-07-20 10:10:27 - Skipping B_cells since one comparison has fewer than 3 cells!
#> 2026-07-20 10:10:27 - Comparing disease with healthy in: CD4_T_cells
#> 2026-07-20 10:10:27 - Skipping CD4_T_cells since one comparison has fewer than 3 cells!
#> 2026-07-20 10:10:27 - Comparing disease with healthy in: Monocytes
#> 2026-07-20 10:10:27 - Skipping Monocytes since one comparison has fewer than 3 cells!
#> 2026-07-20 10:10:27 - Comparing disease with healthy in: NK
#> 2026-07-20 10:10:27 - Skipping NK since one comparison has fewer than 3 cells!
#> 2026-07-20 10:10:27 - Comparing disease with healthy in: pDC
#> 2026-07-20 10:10:27 - Skipping pDC since one comparison has fewer than 3 cells!
#> 2026-07-20 10:10:27 - Comparing disease with healthy in: CD8_T_cells
#> 2026-07-20 10:10:27 - Skipping CD8_T_cells since one comparison has fewer than 3 cells!
#> 2026-07-20 10:10:27 - Finished DGE pseudo bulk method analysis
#> 2026-07-20 10:10:27 - DGE pseudo bulk result is empty...
```
