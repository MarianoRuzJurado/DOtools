library(testthat)
library(Seurat)

test_that("DO.Import runs on a small synthetic dataset and returns Seurat", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scDblFinder")

  set.seed(42)

  # Create a toy counts matrix (100 genes × 20 cells)
  mat <- matrix(rpois(2000, 5), nrow = 100, ncol = 20)
  rownames(mat) <- paste0("gene", 1:100)
  colnames(mat) <- paste0("cell", 1:20)

  # Save as CSV to temporary file
  tmpfile <- tempfile(fileext = ".csv")
  write.csv(mat, file = tmpfile, row.names = TRUE)

  # Call DO.Import with reduced PCs
  obj <- DO.Import(
    pathways = c(tmpfile),
    ids = c("sample1"),
    minCellGenes = 1,
    FilterCells = FALSE,
    DeleteDoublets = FALSE,
    Seurat = TRUE,
    npcs = 10   #cap PCs to fit small dataset
  )

  # Basic checks
  expect_s4_class(obj, "Seurat")
  expect_true("sample" %in% colnames(obj@meta.data))
  expect_true("condition" %in% colnames(obj@meta.data))
  expect_gt(ncol(obj), 0)
  expect_gt(nrow(obj), 0)
})

test_that("DO.Import can return SingleCellExperiment", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")

  set.seed(123)

  # Smaller dataset (30 genes × 10 cells)
  mat <- matrix(rpois(300, 5), nrow = 30, ncol = 10)
  rownames(mat) <- paste0("gene", 1:30)
  colnames(mat) <- paste0("cell", 1:10)

  tmpfile <- tempfile(fileext = ".csv")
  write.csv(mat, tmpfile, row.names = TRUE)

  # Call DO.Import with capped PCs
  obj <- DO.Import(
    pathways = c(tmpfile),
    ids = c("sample2"),
    minCellGenes = 1,
    FilterCells = FALSE,
    DeleteDoublets = FALSE,
    Seurat = FALSE,
    npcs = 5    #lower PCs for tiny dataset
  )

  expect_s4_class(obj, "SingleCellExperiment")
  expect_true("sample" %in% colnames(colData(obj)))
  expect_true("condition" %in% colnames(colData(obj)))
})
