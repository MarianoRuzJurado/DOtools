# DO Celltypist

Runs the CellTypist model on a Seurat or SCE object to predict cell type
labels, storing the results as metadata. If the number of cells is less
than the specified threshold, it returns NAs for the labels. Optionally
updates the CellTypist models and returns the probability matrix. Useful
for annotating cell types in single-cell RNA sequencing datasets.

## Usage

``` r
DO.CellTypist(
  sce_object,
  modelName = "Healthy_Adult_Heart.pkl",
  minCellsToRun = 200,
  runCelltypistUpdate = TRUE,
  over_clustering = "seurat_clusters",
  assay_normalized = "RNA",
  returnAll = FALSE,
  SeuV5 = TRUE
)
```

## Arguments

- sce_object:

  The seurat or sce object

- modelName:

  Specify the model you want to use for celltypist

- minCellsToRun:

  If the input seurat or SCE object has fewer than this many cells, NAs
  will be added for all expected columns and celltypist will not be run.

- runCelltypistUpdate:

  If true, –update-models will be run for celltypist prior to scoring
  cells.

- over_clustering:

  Column in metadata in object with clustering assignments for cells,
  default seurat_clusters

- assay_normalized:

  Assay with log1p normalized expressions

- returnAll:

  will additionally return the probability matrix, return will give a
  list with the first element being the object and second plot and third
  probability matrix

- SeuV5:

  Specify if the Seurat object is made with Seuratv5

## Value

a seurat or sce object

## Author

Mariano Ruz Jurado

## Examples

``` r
sce_data <-
    readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))


sce_data <- DO.CellTypist(
    sce_object = sce_data,
    modelName = "Healthy_Adult_Heart.pkl",
    runCelltypistUpdate = TRUE,
    over_clustering = "seurat_clusters",
    minCellsToRun = 5,
    SeuV5 = TRUE
)
#> 2026-04-16 10:04:06 - Running celltypist using model: Healthy_Adult_Heart.pkl
#> 2026-04-16 10:04:06 - Saving celltypist results to temporary folder: /tmp/RtmpCmogAv/file1e96551f8620
#> Using Python: /home/runner/.pyenv/versions/3.14.0/bin/python3.14
#> Creating virtual environment '/home/runner/.cache/R/basilisk/1.22.0/zellkonverter/1.20.1/zellkonverterAnnDataEnv-0.12.3' ... 
#> + /home/runner/.pyenv/versions/3.14.0/bin/python3.14 -m venv /home/runner/.cache/R/basilisk/1.22.0/zellkonverter/1.20.1/zellkonverterAnnDataEnv-0.12.3
#> Done!
#> Installing packages: pip, wheel, setuptools
#> + /home/runner/.cache/R/basilisk/1.22.0/zellkonverter/1.20.1/zellkonverterAnnDataEnv-0.12.3/bin/python -m pip install --upgrade pip wheel setuptools
#> Installing packages: 'anndata==0.12.3', 'h5py==3.15.1', 'natsort==8.4.0', 'numpy==2.3.4', 'pandas==2.3.3', 'scipy==1.16.2'
#> + /home/runner/.cache/R/basilisk/1.22.0/zellkonverter/1.20.1/zellkonverterAnnDataEnv-0.12.3/bin/python -m pip install --upgrade --no-user 'anndata==0.12.3' 'h5py==3.15.1' 'natsort==8.4.0' 'numpy==2.3.4' 'pandas==2.3.3' 'scipy==1.16.2'
#> Virtual environment '/home/runner/.cache/R/basilisk/1.22.0/zellkonverter/1.20.1/zellkonverterAnnDataEnv-0.12.3' successfully created.
#> 2026-04-16 10:19:01 - Downloading CellTypist model: Healthy_Adult_Heart.pkl
#> 2026-04-16 10:19:05 - Running Celltypist
#> 2026-04-16 10:19:07 - Creating probality plot
```
