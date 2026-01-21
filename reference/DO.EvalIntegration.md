# Do batch correction metrics for integration

This function calculates different metrics to evaluate the integration
of scRNA expression matrices in a new dimension. Its a wrapper function
around scib batch correction metrics

## Usage

``` r
DO.EvalIntegration(
  sce_object,
  label_key = "annotation",
  batch_key = "orig.ident",
  type_ = "embed",
  pcr_covariate = "orig.ident",
  pcr_n_comps = 30,
  scale = TRUE,
  verbose = FALSE,
  n_cores = 10,
  assay = "RNA",
  integration = "INTEGRATED.CCA",
  kBET = TRUE,
  cells.use = NULL,
  subsample = NULL,
  min_per_batch = NULL,
  all_scores_silhouette = FALSE,
  ...
)
```

## Arguments

- sce_object:

  Seurat or SCE object.

- label_key:

  character, Annotation column

- batch_key:

  character, Sample column

- type\_:

  character, default: "embed"

- pcr_covariate:

  character, covariate column for pcr

- pcr_n_comps:

  integer, number of components for pcr

- scale:

  boolean, default: TRUE

- verbose:

  boolean, defult: FALSE

- n_cores:

  integer, Number of cores used for calculations

- assay:

  character, Name of the assay the integration is saved in

- integration:

  character, Name of the integration to evaluate

- kBET:

  boolean, if kBET should be run

- cells.use:

  vector, named cells to use for kBET subsetting

- subsample:

  float, for starified subsampling,

- min_per_batch:

  integer, minimum number of cells per batch

- all_scores_silhouette:

  boolean, define if all scores of silhouette return

- ...:

  Additionally arguments for kBET

## Value

DataFrame with score for the given integration

## Author

Mariano Ruz Jurado

## Examples

``` r
if (FALSE) { # \dontrun{
sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

DO.EvalIntegration(
    sce_object = sce_data,
    label_key = "annotation",
    batch_key = "orig.ident",
    type_ = "embed",
    pcr_covariate = "orig.ident",
    pcr_n_comps = 30,
    scale = TRUE,
    verbose = FALSE,
    n_cores = 10,
    assay = "RNA",
    integration = "INTEGRATED.CCA",
    kBET = TRUE,
    cells.use = NULL,
    subsample = NULL,
    min_per_batch = NULL,
    all_scores_silhouette = FALSE
)
} # }
```
