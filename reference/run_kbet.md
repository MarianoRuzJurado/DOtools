# kBET function

kBET quantifies batch mixing in single-cell data by testing whether
local neighborhood batch composition deviates from the global
distribution, with lower rejection indicating better integration.

## Usage

``` r
run_kbet(
  sce_object,
  embedding,
  batch,
  cells.use = NULL,
  subsample = NULL,
  min_per_batch = NULL,
  ...
)
```

## Arguments

- sce_object:

  Seurat or SCE object.

- embedding:

  Name of the embedding to test

- batch:

  Name of the sample column in meta.data

- cells.use:

  specify a specific set of cells as vector

- subsample:

  set a fraction for stratified subsetting

- min_per_batch:

  should always be higher than or equal k0

## Value

list with summary stats of kBET

## Author

Mariano Ruz Jurado
