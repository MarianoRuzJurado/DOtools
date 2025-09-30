library(testthat)
library(Seurat)
library(SingleCellExperiment)
library(DOtools)
library(Matrix)

test_that("DO.TransferLabel transfers annotations correctly", {
  set.seed(123)

  # Create mock Seurat object
  counts <- matrix(rpois(200*10, lambda = 5), nrow = 200, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:200)
  colnames(counts) <- paste0("Cell", 1:10)
  seurat_obj <- CreateSeuratObject(counts = counts, assay = "RNA")
  seurat_obj$annotation <- rep(c("A", "B"), each = 5)

  # Subset and re-annotate
  subset_obj <- seurat_obj[, 1:5]
  subset_obj$annotation <- rep("C", 5)

  # Transfer back
  transferred <- DO.TransferLabel(
    sce_object = seurat_obj,
    Subset_obj = subset_obj,
    annotation_column = "annotation",
    subset_annotation = "annotation"
  )

  expect_true(is(transferred, "Seurat"))
  expect_equal(as.character(transferred$annotation[1:5]), rep("C", 5))
  expect_equal(as.character(transferred$annotation[6:10]), rep("B", 5))

  # SingleCellExperiment version
  sce_obj <- as.SingleCellExperiment(seurat_obj)

  # Ensure 'logcounts' exists
  assay(sce_obj, "logcounts") <- assay(sce_obj, "counts")

  subset_sce <- as.SingleCellExperiment(subset_obj)
  assay(subset_sce, "logcounts") <- assay(subset_sce, "counts")

  transferred_sce <- DO.TransferLabel(
    sce_object = sce_obj,
    Subset_obj = subset_sce,
    annotation_column = "annotation",
    subset_annotation = "annotation"
  )

  expect_true(is(transferred_sce, "SingleCellExperiment"))
  expect_equal(as.character(colData(transferred_sce)$annotation[1:5]), rep("C", 5))
  expect_equal(as.character(colData(transferred_sce)$annotation[6:10]), rep("B", 5))
})
