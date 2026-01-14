# DO.Import

Imports and processes single-cell RNA-seq data from various formats (10x
Genomics, CellBender, or CSV), performs quality control (QC), filtering
, normalization, variable gene selection, and optionally detects
doublets. Returns a merged and processed Seurat or SCE object ready for
downstream analysis.

## Usage

``` r
DO.Import(
  pathways,
  ids,
  minCellGenes = 5,
  FilterCells = TRUE,
  cut_mt = 0.05,
  min_counts = NULL,
  max_counts = NULL,
  min_genes = NULL,
  max_genes = NULL,
  low_quantile = NULL,
  high_quantile = NULL,
  DeleteDoublets = TRUE,
  include_rbs = TRUE,
  Seurat = TRUE,
  ...
)
```

## Arguments

- pathways:

  A character vector of paths to directories or files containing raw
  expression matrices.

- ids:

  A character vector of sample identifiers, matching the order of
  `pathways`.

- minCellGenes:

  Integer. Minimum number of cells a gene must be expressed in to be
  retained. Default is 5.

- FilterCells:

  Logical. If `TRUE`, applies QC filtering on cells based on
  mitochondrial content, counts, and feature thresholds. Default is
  `TRUE`.

- cut_mt:

  Numeric. Maximum allowed mitochondrial gene proportion per cell.
  Default is 0.05.

- min_counts:

  Numeric. Minimum UMI count threshold (optional, used only if
  `low_quantile` is `NULL`).

- max_counts:

  Numeric. Maximum UMI count threshold (optional, used only if
  `high_quantile` is `NULL`).

- min_genes:

  Numeric. Minimum number of genes detected per cell to retain.
  Optional.

- max_genes:

  Numeric. Maximum number of genes detected per cell to retain.
  Optional.

- low_quantile:

  Numeric. Quantile threshold (0 to 1) to filter low UMI cells (used if
  `min_counts` is `NULL`).

- high_quantile:

  Numeric. Quantile threshold (0 to 1) to filter high UMI cells (used if
  `max_counts` is `NULL`).

- DeleteDoublets:

  Logical. If `TRUE`, doublets are detected and removed using
  `scDblFinder`. Default is `TRUE`.

- include_rbs:

  Logical. If `TRUE`, calculates ribosomal gene content in addition to
  mitochondrial content. Default is `TRUE`.

- Seurat:

  Logical. If `TRUE`, returns Seurat object otherwise SCE object.

- ...:

  Additional arguments passed to `RunPCA()`.

## Value

A merged Seurat or SCE object containing all samples, with
normalization, QC, scaling, PCA, and optional doublet removal applied.

## Author

Mariano Ruz Jurado & David John

## Examples

``` r
if (FALSE) { # \dontrun{
merged_obj <- DO.Import(
    pathways = c("path/to/sample1", "path/to/sample2"),
    ids = c("sample1", "sample2"),
    TenX = TRUE,
    CellBender = FALSE,
    minCellGenes = 5,
    FilterCells = TRUE,
    cut_mt = 0.05,
    min_counts = 1000,
    max_counts = 20000,
    min_genes = 200,
    max_genes = 6000,
    DeleteDoublets = TRUE
)
} # }
```
