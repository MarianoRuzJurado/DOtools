# .QC_Vlnplot

Generates violin plots for common quality control (QC) metrics of
single-cell RNA-seq data from a Seurat object. The function displays
three violin plots for the number of detected genes per cell
(`nFeature_RNA`), total UMI counts per cell (`nCount_RNA`), and
mitochondrial gene content percentage (`pt_mito`). Useful for visual
inspection of QC thresholds and outliers.

## Usage

``` r
.QC_Vlnplot(
  sce_object,
  id,
  layer = "counts",
  features = c("nFeature_RNA", "nCount_RNA", "pt_mito")
)
```

## Arguments

- sce_object:

  A Seurat object containing single-cell RNA-seq data.

- layer:

  A character string specifying the assay layer to use (default is
  "counts").

- features:

  A character vector of length 3 indicating the feature names to plot.
  Default is c("nFeature_RNA", "nCount_RNA", "pt_mito").

## Value

A `ggplot` object arranged in a single row showing violin plots for the
specified features with overlaid boxplots.

## Author

Mariano Ruz Jurado
