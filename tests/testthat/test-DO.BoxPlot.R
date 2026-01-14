library(testthat)
library(Seurat)
library(SingleCellExperiment)
library(DOtools)

test_that("DO.BoxPlot works for basic Seurat object without group.by.2", {
  skip_if_not_installed("Seurat")

  set.seed(42)

  # Create a properly sized dataset with recognizable control group
  mat <- matrix(rpois(10*30, lambda = 5), nrow = 10, ncol = 30)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  # Create Seurat object with proper structure
  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

  # Use recognizable control group names that the function can auto-detect
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  # Test basic functionality with explicit control condition
  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    ctrl.condition = "CTRL"  # Explicitly provide control
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot works with auto-detected control condition", {
  skip_if_not_installed("Seurat")

  set.seed(123)

  mat <- matrix(rpois(10*30, lambda = 5), nrow = 10, ncol = 30)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

  # Use "CTRL" which the function should auto-detect
  seurat_obj$condition <- rep(c("CTRL", "DISEASE"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  # Let function auto-detect control
  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition"
    # ctrl.condition not provided - should auto-detect "CTRL"
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot works with group.by.2", {
  skip_if_not_installed("Seurat")

  set.seed(456)

  # Larger dataset for group.by.2 to ensure enough samples per group
  mat <- matrix(rpois(10*60, lambda = 5), nrow = 10, ncol = 60)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:60)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 30)
  seurat_obj$celltype <- rep(rep(c("Tcell", "Bcell"), each = 15), times = 2)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3", "S4"), each = 15)

  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene2",
    sample.column = "orig.ident",
    group.by = "condition",
    group.by.2 = "celltype",
    ctrl.condition = "CTRL"
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot works with metadata features", {
  skip_if_not_installed("Seurat")

  set.seed(789)

  mat <- matrix(rpois(10*30, lambda = 5), nrow = 10, ncol = 30)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  # Add metadata feature
  seurat_obj$nFeature_RNA <- colSums(mat > 0)  # Real metadata column

  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "nFeature_RNA",  # Use metadata column
    sample.column = "orig.ident",
    group.by = "condition",
    ctrl.condition = "CTRL"
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot works with custom ListTest", {
  skip_if_not_installed("Seurat")

  set.seed(111)

  mat <- matrix(rpois(10*45, lambda = 5), nrow = 10, ncol = 45)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:45)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT1", "TREAT2"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 15)

  # Provide custom comparison list
  ListTest <- list(c("CTRL", "TREAT1"), c("CTRL", "TREAT2"))

  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    ListTest = ListTest  # Use custom comparisons
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot works without statistical test", {
  skip_if_not_installed("Seurat")

  set.seed(222)

  mat <- matrix(rpois(10*30, lambda = 5), nrow = 10, ncol = 30)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    ctrl.condition = "CTRL",
    test_use = "none"  # Disable statistical test
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot handles SingleCellExperiment objects", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")

  set.seed(333)

  mat <- matrix(rpois(10*30, lambda = 5), nrow = 10, ncol = 30)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  # Create a proper SingleCellExperiment object that will convert correctly
  # The function uses as.Seurat() which expects certain structure
  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  # Convert to SingleCellExperiment using Seurat's conversion
  sce <- as.SingleCellExperiment(seurat_obj)

  # Test that conversion works properly
  p <- DO.BoxPlot(
    sce_object = sce,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    ctrl.condition = "CTRL"
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot works with different parameters", {
  skip_if_not_installed("Seurat")

  set.seed(444)

  mat <- matrix(rpois(10*30, lambda = 5), nrow = 10, ncol = 30)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  # Test various parameter combinations
  p1 <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    ctrl.condition = "CTRL",
    plot_sample = FALSE,
    outlier_removal = FALSE
  )

  p2 <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    ctrl.condition = "CTRL",
    orderAxis = c("TREAT", "CTRL")  # Custom order
  )

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("DO.BoxPlot handles edge cases gracefully", {
  skip_if_not_installed("Seurat")

  set.seed(555)

  mat <- matrix(rpois(10*30, lambda = 5), nrow = 10, ncol = 30)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  # Test with empty ListTest (should use default behavior)
  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    ctrl.condition = "CTRL",
    ListTest = NULL  # Empty list
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot works with single celltype in group.by.2", {
  skip_if_not_installed("Seurat")

  set.seed(666)

  mat <- matrix(rpois(10*30, lambda = 5), nrow = 10, ncol = 30)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$celltype <- "Tcell"  # Only one cell type
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    group.by.2 = "celltype",
    ctrl.condition = "CTRL"
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot error handling for invalid inputs", {
  skip_if_not_installed("Seurat")

  set.seed(777)

  mat <- matrix(rpois(10*30, lambda = 5), nrow = 10, ncol = 30)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  # Test with non-existent feature - expect error with the exact message
  expect_error(
    DO.BoxPlot(
      sce_object = seurat_obj,
      Feature = "nonexistent_gene",
      sample.column = "orig.ident",
      group.by = "condition"
    ),
    "Feature not found in SCE Object!"
  )
})

test_that("DO.BoxPlot works with larger dataset", {
  skip_if_not_installed("Seurat")

  set.seed(888)

  # Larger dataset for performance testing
  mat <- matrix(rpois(100*100, lambda = 5), nrow = 100, ncol = 100)
  rownames(mat) <- paste0("gene", 1:100)
  colnames(mat) <- paste0("cell", 1:100)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT1", "TREAT2"), length.out = 100)
  seurat_obj$celltype <- rep(c("Tcell", "Bcell"), length.out = 100)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3", "S4", "S5"), each = 20)

  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    group.by.2 = "celltype",
    ctrl.condition = "CTRL"
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot handles partial zero expression gracefully", {
  skip_if_not_installed("Seurat")

  set.seed(999)

  # Create matrix with partial zero expression (one group has expression, one doesn't)
  mat <- matrix(0, nrow = 10, ncol = 30)
  # Add expression only to the CTRL group
  mat[1:5, 1:15] <- rpois(5*15, lambda = 2)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  # This should handle partial zero expression gracefully
  # Test gene1 which has expression in CTRL but not in TREAT
  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",
    sample.column = "orig.ident",
    group.by = "condition",
    ctrl.condition = "CTRL",
    test_use = "none"
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.BoxPlot works with minimal non-zero expression", {
  skip_if_not_installed("Seurat")

  set.seed(101010)

  # Create matrix with minimal non-zero expression to avoid complete zeros
  mat <- matrix(rpois(10*30, lambda = 0.1), nrow = 10, ncol = 30)  # Very low lambda
  # Ensure at least some cells have non-zero counts
  mat[1, ] <- 1  # Set first gene to have at least 1 count in all cells
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:30)

  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$condition <- rep(c("CTRL", "TREAT"), each = 15)
  seurat_obj$orig.ident <- rep(c("S1", "S2", "S3"), each = 10)

  # Test with gene that has minimal expression
  p <- DO.BoxPlot(
    sce_object = seurat_obj,
    Feature = "gene1",  # This gene has at least some expression
    sample.column = "orig.ident",
    group.by = "condition",
    ctrl.condition = "CTRL"
  )

  expect_s3_class(p, "ggplot")
})
