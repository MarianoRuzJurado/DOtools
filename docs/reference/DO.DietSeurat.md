# Remove Layers from Seurat Object by Pattern

This function removes layers from a Seurat object's RNA assay based on a
specified regular expression pattern. It is supposed to make something
similar than the no longer working DietSeurat function, by removing no
longer needed layers from th object.

## Usage

``` r
DO.DietSeurat(Seu_object, pattern = "^scale\\.data\\.")
```

## Arguments

- Seu_object:

  Seurat object.

- pattern:

  regular expression pattern to match layer names. Default
  "^scale\\data\\"

## Value

Seurat object with specified layers removed.

## Author

Mariano Ruz Jurado

## Examples

``` r
sc_data <- readRDS(system.file("extdata", "sc_data.rds", package = "DOtools"))

sc_data <- DO.DietSeurat(sc_data, pattern = "data")
#> pattern:  data
#> data is removed.scale.data is removed.

```
