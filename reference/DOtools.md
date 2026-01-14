# DOtools: A Toolkit for scRNA Data Analysis

The \`DOtools\` package provides a set of functions for advanced data
processing, visualisation, and statistical analysis in Seurat objects.
It includes functions for cell-type prediction, reclustering, creating
polished UMAP plots, subsetting Seurat objects, and various statistical
analyses like Wilcoxon tests and SEM graphs.

## Value

This is a package-level documentation file and does not return a value.

## Details

This package includes the following functions:

- [`DO.BoxPlot`](https://marianoruzjurado.github.io/DOtools/reference/DO.BoxPlot.md):
  A function for creating box plots with Wilcoxon test results.

- [`DO.CellTypist`](https://marianoruzjurado.github.io/DOtools/reference/DO.CellTypist.md):
  A function for running CellTypist on Seurat and SCE objects to predict
  cell types.

- [`DO.DietSCE`](https://marianoruzjurado.github.io/DOtools/reference/DO.DietSCE.md):
  A function for diet-based analysis of Seurat and SCE objects.

- [`DO.Dotplot`](https://marianoruzjurado.github.io/DOtools/reference/DO.Dotplot.md):
  A function for creating dot plots for visualizing gene expression.

- [`DO.FullRecluster`](https://marianoruzjurado.github.io/DOtools/reference/DO.FullRecluster.md):
  A function for fine-grained reclustering of Seurat and SCE objects.

- [`DO.BarplotClustert`](https://marianoruzjurado.github.io/DOtools/reference/DO.BarplotClustert.md):
  A function for generating mean and SEM graphs for cluster-based
  analysis with t-tests.

- [`DO.Barplot`](https://marianoruzjurado.github.io/DOtools/reference/DO.Barplot.md):
  A function for generating mean and SEM graphs with a statistical test
  indicating significance.

- [`DO.Subset`](https://marianoruzjurado.github.io/DOtools/reference/DO.Subset.md):
  A function for subsetting Seurat and SCE objects based on metadata.

- [`DO.UMAP`](https://marianoruzjurado.github.io/DOtools/reference/DO.UMAP.md):
  A function for creating polished UMAP plots using either DimPlot or
  FeaturePlot.

- [`DO.VlnPlot`](https://marianoruzjurado.github.io/DOtools/reference/DO.VlnPlot.md):
  A function for generating violin plots with Wilcoxon test results.

- [`DO.CellComposition`](https://marianoruzjurado.github.io/DOtools/reference/DO.CellComposition.md):
  A function for visualizing and statistically analyzing cell-type
  composition changes across conditions using the Scanpro Python
  package, with support for bootstrapping, proportion plots, and
  customizable output.

- [`DO.Import`](https://marianoruzjurado.github.io/DOtools/reference/DO.Import.md):
  A function for building a merged Seurat and SCE object from 10x
  software output, or directly from provided tables.

- [`DO.Integration`](https://marianoruzjurado.github.io/DOtools/reference/DO.Integration.md):
  A function for integrating SCE objects and Seurat objects with the
  provided method.

- [`DO.CellBender`](https://marianoruzjurado.github.io/DOtools/reference/DO.CellBender.md):
  A function for running CellBender in a virtual conda env with provided
  raw count h5 files.

- [`DO.SplitBarGSEA`](https://marianoruzjurado.github.io/DOtools/reference/DO.SplitBarGSEA.md):
  A function for viusalizing GSEA result from a provided df from e.g.
  metascape

- [`DO.scVI`](https://marianoruzjurado.github.io/DOtools/reference/DO.scVI.md):
  A function for running the scVI Integration implemented in scvi-tools.

- [`DO.TransferLabel`](https://marianoruzjurado.github.io/DOtools/reference/DO.TransferLabel.md):
  A function for transfering annotation from a subseted object to the
  original seurat and SCE object.

- [`DO.PyEnv`](https://marianoruzjurado.github.io/DOtools/reference/DO.PyEnv.md):
  A function for creating a conda envrionment holding all python
  packages needed for some functions.

- [`DO.Correlation`](https://marianoruzjurado.github.io/DOtools/reference/DO.Correlation.md):
  A function for creating a correlation plot between provided samples in
  the category specified.

- [`DO.Heatmap`](https://marianoruzjurado.github.io/DOtools/reference/DO.Heatmap.md):
  A function for generating Heat maps on gene expression data.

- [`DO.HeatmapFC`](https://marianoruzjurado.github.io/DOtools/reference/DO.HeatmapFC.md):
  A function for generating Heat maps showing foldchanges in expression
  between specified conditions.

- [`DO.MultiDGE`](https://marianoruzjurado.github.io/DOtools/reference/DO.MultiDGE.md):
  A function for calculating DEGs on a single cell and speudo bulk
  level.

- `dot-Do.BarcodeRanks`: A function for estimating the number of
  expected cells and droplets.

- `dot-QC.Vlnplot`: A function for estimating the number of expected
  cells and droplets.

## Author

Mariano Ruz Jurado, David Rodriguez Morales

## See also

[`DO.BoxPlot`](https://marianoruzjurado.github.io/DOtools/reference/DO.BoxPlot.md),
[`DO.CellTypist`](https://marianoruzjurado.github.io/DOtools/reference/DO.CellTypist.md),
[`DO.DietSCE`](https://marianoruzjurado.github.io/DOtools/reference/DO.DietSCE.md),
[`DO.Dotplot`](https://marianoruzjurado.github.io/DOtools/reference/DO.Dotplot.md),
[`DO.FullRecluster`](https://marianoruzjurado.github.io/DOtools/reference/DO.FullRecluster.md),
[`DO.BarplotClustert`](https://marianoruzjurado.github.io/DOtools/reference/DO.BarplotClustert.md),
[`DO.Barplot`](https://marianoruzjurado.github.io/DOtools/reference/DO.Barplot.md),
[`DO.Subset`](https://marianoruzjurado.github.io/DOtools/reference/DO.Subset.md),
[`DO.UMAP`](https://marianoruzjurado.github.io/DOtools/reference/DO.UMAP.md),
[`DO.VlnPlot`](https://marianoruzjurado.github.io/DOtools/reference/DO.VlnPlot.md),
[`DO.Import`](https://marianoruzjurado.github.io/DOtools/reference/DO.Import.md),[`DO.Integration`](https://marianoruzjurado.github.io/DOtools/reference/DO.Integration.md),
[`DO.CellBender`](https://marianoruzjurado.github.io/DOtools/reference/DO.CellBender.md),
[`DO.SplitBarGSEA`](https://marianoruzjurado.github.io/DOtools/reference/DO.SplitBarGSEA.md),
[`DO.scVI`](https://marianoruzjurado.github.io/DOtools/reference/DO.scVI.md),
[`DO.TransferLabel`](https://marianoruzjurado.github.io/DOtools/reference/DO.TransferLabel.md),[`DO.Heatmap`](https://marianoruzjurado.github.io/DOtools/reference/DO.Heatmap.md),[`DO.HeatmapFC`](https://marianoruzjurado.github.io/DOtools/reference/DO.HeatmapFC.md),
[`DO.PyEnv`](https://marianoruzjurado.github.io/DOtools/reference/DO.PyEnv.md),
[`DO.Correlation`](https://marianoruzjurado.github.io/DOtools/reference/DO.Correlation.md),
[`DO.MultiDGE`](https://marianoruzjurado.github.io/DOtools/reference/DO.MultiDGE.md),[`DO.TransferLabel`](https://marianoruzjurado.github.io/DOtools/reference/DO.TransferLabel.md),
`dot-Do.BarcodeRanks`, `dot-QC.Vlnplot`
