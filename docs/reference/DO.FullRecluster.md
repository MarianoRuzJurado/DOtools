# DO.FullRecluster

Performs iterative reclustering on each major cluster found by
FindClusters in a Seurat or SCE object. It refines the clusters using
the FindSubCluster function for better resolution and fine-tuned
annotation. The new clustering results are stored in a metadata column
called `annotation_recluster`. Suitable for improving cluster precision
and granularity after initial clustering.

## Usage

``` r
DO.FullRecluster(
  sce_object,
  over_clustering = "seurat_clusters",
  res = 0.5,
  algorithm = 4,
  graph.name = "RNA_snn"
)
```

## Arguments

- sce_object:

  The seurat or SCE object

- over_clustering:

  Column in metadata in object with clustering assignments for cells,
  default seurat_clusters

- res:

  Resolution for the new clusters, default 0.5

- algorithm:

  Set one of the available algorithms found in FindSubCLuster function,
  default = 4: leiden

- graph.name:

  A builded neirest neighbor graph

## Value

a Seurat or SCE Object with new clustering named annotation_recluster

## Author

Mariano Ruz Jurado

## Examples

``` r
sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

sce_data <- DO.FullRecluster(
    sce_object = sce_data
)
#> Computing nearest neighbor graph
#> Computing SNN
#> 2 singletons identified. 3 final clusters.
#> 1 singletons identified. 3 final clusters.
#> 1 singletons identified. 2 final clusters.
#> 1 singletons identified. 3 final clusters.
#> 
```
