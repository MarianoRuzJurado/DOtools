library(testthat)
library(Seurat)
library(DOtools)
library(SingleCellExperiment)

test_that("DO.FullRecluster basic functionality", {
  # Create a larger counts matrix for PCA
  set.seed(123)
  counts <- matrix(rpois(2000, lambda = 5), nrow = 200, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:200)
  colnames(counts) <- paste0("Cell", 1:10)

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, assay = "RNA")

  # Normalize, find variable features, scale
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = FALSE)

  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, npcs = 5, verbose = FALSE)

  # Build neighbors and initial clusters
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:5, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

  # Suppress messages/warnings during reclustering
  reclustered <- suppressMessages(suppressWarnings(
    DO.FullRecluster(
      sce_object = seurat_obj,
      over_clustering = "seurat_clusters",
      res = 0.5,
      algorithm = 4,
      graph.name = "RNA_snn"
    )
  ))

  # Validate output
  expect_true("annotation_recluster" %in% colnames(reclustered@meta.data))
  expect_true(is(reclustered, "Seurat"))
})
