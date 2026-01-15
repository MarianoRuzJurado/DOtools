# DO.Integration SCE object integration directly

Integrates single-cell RNA-seq data directly from SingleCellExperiment
or Seurat objects. Supports detection of variable genes , scaling, PCA,
neighbor graph construction, clustering, and UMAP embedding, with
multiple integration methods.

## Usage

``` r
DO.Integration(
  sce_object,
  split_key = "orig.ident",
  HVG = FALSE,
  scale = FALSE,
  pca = FALSE,
  neighbors = TRUE,
  neighbors_dim = seq_len(50),
  clusters = TRUE,
  clusters_res = 0.3,
  clusters_algorithm = 4,
  umap = TRUE,
  umap_key = "UMAP",
  umap_dim = seq_len(50),
  integration_method = "CCAIntegration",
  selection_method = "vst",
  loess_span = 0.3,
  clip_max = "auto",
  num_bin = 20,
  binning_method = "equal_width",
  scale_max = 10,
  pca_key = "PCA",
  integration_key = "INTEGRATED.CCA",
  npcs = 50,
  verbose = FALSE
)
```

## Arguments

- sce_object:

  Seurat or SCE Object

- split_key:

  Character. Column in meta data to split the samples by, default
  orig.ident

- HVG:

  Logical. Perform detection of highly variable genes

- scale:

  Logical. Perform scaling of the expression data

- pca:

  Logical. Perform principal component analysis

- neighbors:

  Logical. Perform Nearest-neighbor graph after integration

- neighbors_dim:

  Numeric range. Dimensions of reduction to use as input

- clusters:

  Logical. Perform clustering of cells

- clusters_res:

  Numeric. Value of the resolution parameter, use a value above (below)
  1.0 if you want to obtain a larger (smaller) number of communities.

- clusters_algorithm:

  Numeric. Define the algorithm for clustering, default 4 for "Leiden"

- umap:

  Logical. Runs the Uniform Manifold Approximation and Projection

- umap_key:

  Character name for

- umap_dim:

  Numeric range. Which dimensions to use as input features

- integration_method:

  Character. Define the integration method, please check what versions
  are supported in Seurat::IntegrateLayers function

- selection_method:

  Character. Default "vst". Options: "mean.var.plot", "dispersion"

- loess_span:

  Numeric. Loess span parameter used when fitting the variance-mean
  relationship

- clip_max:

  Character. After standardization values larger than clip.max will be
  set to clip.max; default is 'auto' which sets this value to the square
  root of the number of cells

- num_bin:

  Numeric. Total number of bins to use in the scaled analysis (default
  is 20)

- binning_method:

  Character. “equal_width”: each bin is of equal width along the x-axis
  (default). Options: “equal_frequency”:

- scale_max:

  Numeric. Max value to return for scaled data. The default is 10.

- pca_key:

  Character. Key name to save the pca result in

- integration_key:

  Character. Key name to save the integration result in

- npcs:

  Numeric. Total Number of PCs to compute and store (50 by default)

- verbose:

  Logical. Verbosity for all functions

## Value

integrated sce/seurat object

## Author

Mariano Ruz Jurado

## Examples

``` r
sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

DO.Integration(
    sce_object = sce_data,
    split_key = "orig.ident",
    HVG = TRUE,
    scale = TRUE,
    pca = TRUE,
    integration_method = "CCAIntegration"
)
#> 2026-01-14 16:45:37 - Splitting object for integration with CCAIntegration by orig.ident
#> 2026-01-14 16:45:37 - Calculating highly variable genes
#> 2026-01-14 16:45:37 - Scaling object
#> 2026-01-14 16:45:38 - Running pca, saved in key: PCA
#> Splitting ‘counts’, ‘data’ layers. Not splitting ‘scale.data’. If you would like to split other layers, set in `layers` argument.
#> 2026-01-14 16:45:39 - Running integration, saved in key: INTEGRATED.CCA
#> 2026-01-14 16:45:43 - Running Nearest-neighbor graph construction
#> 2026-01-14 16:45:44 - Running cluster detection
#> 2026-01-14 16:45:45 - Creating UMAP
#> class: SingleCellExperiment 
#> dim: 800 2807 
#> metadata(0):
#> assays(3): counts logcounts scaledata
#> rownames(800): HES4 ISG15 ... SERPINA9 DSG2
#> rowData names(0):
#> colnames(2807): AAACCCAGTGCATTTG-1_1 AAACCCATCTCAACGA-1_1 ...
#>   TTTGGAGCAACTGGTT-1_2 TTTGGAGGTTACCTGA-1_2
#> colData names(15): orig.ident nCount_RNA ... annotation_recluster
#>   leiden0.3
#> reducedDimNames(3): PCA INTEGRATED.CCA UMAP
#> mainExpName: RNA
#> altExpNames(0):
```
