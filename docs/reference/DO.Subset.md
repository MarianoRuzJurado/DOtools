# DO.Subset

Creates a subset of a Seurat or SCE object based on either categorical
or numeric thresholds in metadata. Allows for subsetting by specifying
the ident column, group name, or threshold criteria. Ideal for
extracting specific cell populations or clusters based on custom
conditions. Returns a new Seurat or SCE object containing only the
subsetted cells and does not come with the Seuratv5 subset issue. Please
be aware that right now, after using this function the subset might be
treated with Seuv5=False in other functions.

## Usage

``` r
DO.Subset(
  sce_object,
  assay = "RNA",
  ident,
  ident_name = NULL,
  ident_thresh = NULL
)
```

## Arguments

- sce_object:

  The seurat or SCE object

- assay:

  assay to subset by

- ident:

  meta data column to subset for

- ident_name:

  name of group of barcodes in ident of subset for

- ident_thresh:

  numeric thresholds as character, e.g "\>5" or c("\>5", "\<200"), to
  subset barcodes in ident

## Value

a subsetted Seurat or SCE object

## Author

Mariano Ruz Jurado

## Examples

``` r
sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

sce_data_sub <- DO.Subset(
    sce_object = sce_data,
    ident = "condition",
    ident_name = "healthy"
)
#> 2026-01-14 13:00:37 - Specified 'ident_name': expecting a categorical variable.
```
