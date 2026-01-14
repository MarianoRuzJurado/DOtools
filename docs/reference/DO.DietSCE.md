# Remove Layers from Seurat or SCE Object by Pattern

This function removes layers from a Seurat or SCE object's RNA assay
based on a specified regular expression pattern. It is supposed to
remove no longer needed layers from th object.

## Usage

``` r
DO.DietSCE(sce_object, assay = "RNA", pattern = "^scale\\.data\\.")
```

## Arguments

- sce_object:

  Seurat or SCE object.

- assay:

  Name of the assay from where to remove layers from

- pattern:

  regular expression pattern to match layer names. Default
  "^scale\\data\\"

## Value

Seurat or SCE object with specified layers removed.

## Author

Mariano Ruz Jurado

## Examples

``` r
sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

sce_data <- DO.DietSCE(sce_data, pattern = "data")
#> 2026-01-14 13:00:02 - pattern:  data
#> 2026-01-14 13:00:02 - Object has no layers, pattern does not need to be removed from layers.
```
