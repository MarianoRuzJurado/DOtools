# DO.scVI

This function will run the scVI Integration from the scVI python
package. It includes all parameters from the actual python package and
runs it by using an internal python script. The usage of a gpu is
incorporated and highly recommended.

## Usage

``` r
DO.scVI(
  sce_object,
  batch_key,
  layer_counts = "counts",
  layer_logcounts = "logcounts",
  categorical_covariates = NULL,
  continuos_covariates = NULL,
  n_hidden = 128,
  n_latent = 30,
  n_layers = 3,
  dispersion = "gene-batch",
  gene_likelihood = "zinb",
  get_model = FALSE
)
```

## Arguments

- sce_object:

  Seurat or SCE object with annotation in meta.data

- batch_key:

  meta data column with batch information.

- layer_counts:

  layer with counts. Raw counts are required.

- layer_logcounts:

  layer with log-counts. Log-counts required for calculation of HVG.

- categorical_covariates:

  list of meta data column names with categorical covariates for scVI
  inference.

- continuos_covariates:

  list of meta data column names with continuous covariates for scVI
  inference.

- n_hidden:

  number of hidden layers.

- n_latent:

  dimensions of the latent space.

- n_layers:

  number of layers.

- dispersion:

  dispersion mode for scVI.

- gene_likelihood:

  gene likelihood.

- get_model:

  return the trained model.

## Value

Seurat or SCE Object with dimensionality reduction from scVI

## Examples

``` r
if (FALSE) { # \dontrun{
sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

# Run scVI using the 'orig.ident' column as the batch key
sce_data <- DO.scVI(sce_data, batch_key = "orig.ident")
} # }
```
