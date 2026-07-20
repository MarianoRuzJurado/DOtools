# DO.scANVI

Run the scANVI annotation transfer and latent representation inference
using a previously trained scVI model. The function loads an existing
scVI model, initializes a scANVI model, trains it, optionally saves the
trained scANVI model, and returns the latent embedding together with
optional predicted cell labels in the SCE or Seurat object provided.

## Usage

``` r
DO.scANVI(
  sce_object,
  model_path,
  labels_key,
  unlabeled_category,
  batch_size = 128,
  save_model = NULL,
  annot_predict = TRUE
)
```

## Arguments

- sce_object:

  A Seurat or SingleCellExperiment object.

- model_path:

  Path to a previously trained scVI model.

- labels_key:

  Name of the cell annotation column used during scANVI training.

- unlabeled_category:

  Name of the category representing unlabeled cells.

- batch_size:

  Mini-batch size used during scANVI training.

- save_model:

  Optional path where the trained scANVI model should be saved. If
  `NULL`, the model is not saved.

- annot_predict:

  Logical; if `TRUE`, predicted cell labels are returned and stored in
  the output object.

## Value

Seurat or SCE Object with dimensionality reduction from scANVI and
annotation column from scANVI

## Examples

``` r
if (FALSE) { # \dontrun{
sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

# Run scANVI using a previous saved scVI model and annotation categories
sce_data <- DO.scANVI(
    sce_data,
    model_path = "/path/scVI_model",
    labels_key = "annotation",
    unlabeled_category = "Unknown")
} # }
```
