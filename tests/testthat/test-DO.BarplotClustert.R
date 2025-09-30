library(testthat)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)

# Helper function to create test Seurat object
create_test_seurat <- function() {
  set.seed(42)

  # Toy counts matrix: 10 genes x 12 cells
  mat <- matrix(rpois(120, lambda = 5), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:12)

  # Metadata: four conditions, three cells each
  meta <- data.frame(
    orig.ident = rep(c("sample1", "sample2", "sample3", "sample4"), each = 3),
    condition = rep(c("A", "B", "C", "D"), each = 3),
    cluster = rep(c("cluster1", "cluster2"), times = 6)
  )
  rownames(meta) <- colnames(mat)

  # Create Seurat object - convert matrix to dgCMatrix to avoid warnings
  mat_sparse <- as(mat, "dgCMatrix")
  seurat_obj <- CreateSeuratObject(counts = mat_sparse, meta.data = meta)
  seurat_obj <- NormalizeData(seurat_obj)

  # Add a metadata feature
  seurat_obj$metadata_feature <- rnorm(ncol(seurat_obj))

  return(seurat_obj)
}

# Helper function to create test SCE object
create_test_sce <- function() {
  set.seed(42)

  # Toy counts matrix: 10 genes x 12 cells
  mat <- matrix(rpois(120, lambda = 5), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:12)

  # Create SingleCellExperiment
  sce <- SingleCellExperiment(
    assays = list(counts = mat, logcounts = log1p(mat)),
    colData = data.frame(
      orig.ident = rep(c("sample1", "sample2", "sample3", "sample4"), each = 3),
      condition = rep(c("A", "B", "C", "D"), each = 3),
      cluster = rep(c("cluster1", "cluster2"), times = 6)
    )
  )

  return(sce)
}

test_that("DO.BarplotClustert returns ggplot and data frame correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test 1: returnPlot = TRUE (ggplot)
  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    returnPlot = TRUE,
    returnValues = FALSE,
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p, "ggplot")

  # Test 2: returnValues = TRUE (data frame)
  df <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    returnPlot = FALSE,
    returnValues = TRUE,
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(df, "data.frame")
  expect_true(all(c("condition", "variable", "Mean", "SEM") %in% colnames(df)))
  # Fixed: There are 4 conditions in the test data, not 2
  expect_equal(nrow(df), 4)
})

test_that("DO.BarplotClustert works with SingleCellExperiment objects", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  sce_obj <- create_test_sce()

  # Test with SCE object
  p <- DO.BarplotClustert(
    sce_object = sce_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p, "ggplot")
})

test_that("DO.BarplotClustert works with metadata features", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test with metadata feature
  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "metadata_feature",
    ListTest = list(c("A","B")),
    ctrl.condition = "A",
    group.by = "condition",
    log1p_nUMI = FALSE
  )
  expect_s3_class(p, "ggplot")
})

test_that("DO.BarplotClustert handles automatic ListTest generation", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test with NULL ListTest and explicit ctrl condition
  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = NULL,
    ctrl.condition = "A",  # Explicitly set control
    group.by = "condition"
  )
  expect_s3_class(p, "ggplot")

  # Test with NULL ListTest and automatic ctrl detection
  # Create a test object with a condition that matches the pattern
  set.seed(42)
  mat <- matrix(rpois(60, lambda = 5), nrow = 5, ncol = 12)
  rownames(mat) <- paste0("gene", 1:5)
  colnames(mat) <- paste0("cell", 1:12)

  meta <- data.frame(
    orig.ident = rep(c("sample1", "sample2", "sample3"), each = 4),
    condition = rep(c("CTRL", "disease1", "disease2"), each = 4)
  )
  rownames(meta) <- colnames(mat)

  mat_sparse <- as(mat, "dgCMatrix")
  seurat_ctrl <- CreateSeuratObject(counts = mat_sparse, meta.data = meta)
  seurat_ctrl <- NormalizeData(seurat_ctrl)

  p2 <- DO.BarplotClustert(
    sce_object = seurat_ctrl,
    Feature = "gene1",
    ListTest = NULL,
    ctrl.condition = NULL,  # Let it detect automatically
    group.by = "condition"
  )
  expect_s3_class(p2, "ggplot")
})

test_that("DO.BarplotClustert handles custom visual parameters", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test with custom colors
  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    bar_colours = c("red", "blue", "green", "yellow"),
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p, "ggplot")

  # Test with custom y-limits
  p2 <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    y_limits = c(0, 10),
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p2, "ggplot")

  # Test with x-label rotation
  p3 <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    x_label_rotation = 90,
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p3, "ggplot")
})

test_that("DO.BarplotClustert handles log1p_nUMI parameter", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test with log1p_nUMI = FALSE
  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    log1p_nUMI = FALSE,
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p, "ggplot")

  # Test with log1p_nUMI = TRUE (default)
  p2 <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    log1p_nUMI = TRUE,
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p2, "ggplot")
})

test_that("DO.BarplotClustert handles different group.by parameters", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test with different grouping variable
  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("cluster1", "cluster2")),
    ctrl.condition = "cluster1",
    group.by = "cluster"
  )
  expect_s3_class(p, "ggplot")
})

test_that("DO.BarplotClustert handles single feature only", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test with single feature (the function only supports one feature)
  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",  # Single feature only
    ListTest = list(c("A","B")),
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p, "ggplot")
})

test_that("DO.BarplotClustert handles edge cases", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test with empty ListTest (should generate automatically)
  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(),
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p, "ggplot")

  # Test with step_mod and stat_pos_mod parameters
  p2 <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B"), c("A","C")),
    ctrl.condition = "A",
    group.by = "condition",
    step_mod = 0.3,
    stat_pos_mod = 1.2
  )
  expect_s3_class(p2, "ggplot")
})

test_that("DO.BarplotClustert handles error conditions", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test error when feature not found
  expect_error(
    DO.BarplotClustert(
      sce_object = seurat_obj,
      Feature = "nonexistent_gene",
      ListTest = list(c("A","B")),
      ctrl.condition = "A",
      group.by = "condition"
    ),
    "Feature not found in Seurat Object!"
  )

  # Test error when group.by column not found
  expect_error(
    DO.BarplotClustert(
      sce_object = seurat_obj,
      Feature = "gene1",
      ListTest = list(c("A","B")),
      ctrl.condition = "A",
      group.by = "nonexistent_column"
    )
  )
})

test_that("DO.BarplotClustert handles zero-mean conditions", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  # Create a test object where some conditions have zero expression
  set.seed(42)
  mat <- matrix(0, nrow = 5, ncol = 12)  # All zeros
  rownames(mat) <- paste0("gene", 1:5)
  colnames(mat) <- paste0("cell", 1:12)

  # Set some non-zero values for condition A
  mat[1:3, 1:3] <- 5  # Condition A has expression

  meta <- data.frame(
    orig.ident = rep(c("sample1", "sample2", "sample3", "sample4"), each = 3),
    condition = rep(c("A", "B", "C", "D"), each = 3)
  )
  rownames(meta) <- colnames(mat)

  mat_sparse <- as(mat, "dgCMatrix")
  seurat_obj <- CreateSeuratObject(counts = mat_sparse, meta.data = meta)
  seurat_obj <- NormalizeData(seurat_obj)

  # This should work and remove zero-mean comparisons with a warning
  expect_warning(
    p <- DO.BarplotClustert(
      sce_object = seurat_obj,
      Feature = "gene1",
      ListTest = list(c("A","B"), c("B","C")),  # B vs C both zero
      ctrl.condition = "A",
      group.by = "condition"
    ),
    "Removing Test"
  )
  expect_s3_class(p, "ggplot")
})

test_that("DO.BarplotClustert return values structure is correct", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  df <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    returnPlot = FALSE,
    returnValues = TRUE,
    ctrl.condition = "A",
    group.by = "condition"
  )

  # Check structure of returned data frame
  expect_named(df, c("condition", "variable", "Mean", "SEM"))
  expect_type(df$condition, "character")
  # Fixed: variable is a factor, not character
  expect_s3_class(df$variable, "factor")
  expect_type(df$Mean, "double")
  expect_type(df$SEM, "double")

  # Check that SEM is calculated correctly (non-negative)
  expect_true(all(df$SEM >= 0 | is.na(df$SEM)))
})

test_that("DO.BarplotClustert handles multiple comparisons", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test with multiple comparisons - ensure enough samples per group
  # Create a larger dataset for reliable t-tests
  set.seed(42)
  mat <- matrix(rpois(200, lambda = 5), nrow = 5, ncol = 40)
  rownames(mat) <- paste0("gene", 1:5)
  colnames(mat) <- paste0("cell", 1:40)

  meta <- data.frame(
    orig.ident = rep(c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7", "sample8"), each = 5),
    condition = rep(c("A", "B", "C", "D"), each = 10)
  )
  rownames(meta) <- colnames(mat)

  mat_sparse <- as(mat, "dgCMatrix")
  seurat_large <- CreateSeuratObject(counts = mat_sparse, meta.data = meta)
  seurat_large <- NormalizeData(seurat_large)

  # Test with multiple comparisons
  p <- DO.BarplotClustert(
    sce_object = seurat_large,
    Feature = "gene1",
    ListTest = list(c("A","B"), c("A","C"), c("B","D")),
    ctrl.condition = "A",
    group.by = "condition"
  )
  expect_s3_class(p, "ggplot")
})

test_that("DO.BarplotClustert plot elements are present", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  # Create a dataset with matching dimensions
  set.seed(42)
  mat <- matrix(rpois(120, lambda = 5), nrow = 5, ncol = 24)
  rownames(mat) <- paste0("gene", 1:5)
  colnames(mat) <- paste0("cell", 1:24)

  # Fix: Ensure metadata has exactly 24 rows to match 24 columns in matrix
  meta <- data.frame(
    orig.ident = rep(c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6"), each = 4),
    condition = rep(c("A", "B"), each = 12)
  )
  rownames(meta) <- colnames(mat)  # This ensures exact match

  mat_sparse <- as(mat, "dgCMatrix")
  seurat_obj <- CreateSeuratObject(counts = mat_sparse, meta.data = meta)
  seurat_obj <- NormalizeData(seurat_obj)

  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    ctrl.condition = "A",
    group.by = "condition"
  )

  # Check that essential plot elements are present
  expect_s3_class(p, "ggplot")

  # Extract plot data - handle potential warnings from stat_signif
  suppressWarnings({
    plot_data <- ggplot2::ggplot_build(p)
  })

  # Check that we have the expected layers
  layer_types <- sapply(plot_data$plot$layers, function(layer) class(layer$geom)[1])
  expect_true("GeomCol" %in% layer_types)  # Bar plot
  expect_true("GeomErrorbar" %in% layer_types)  # Error bars
  expect_true("GeomPoint" %in% layer_types)  # Individual points
})

test_that("DO.BarplotClustert works with minimal parameters", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  # Create a test object with a condition that matches the pattern for auto-detection
  set.seed(42)
  mat <- matrix(rpois(60, lambda = 5), nrow = 5, ncol = 12)
  rownames(mat) <- paste0("gene", 1:5)
  colnames(mat) <- paste0("cell", 1:12)

  meta <- data.frame(
    orig.ident = rep(c("sample1", "sample2", "sample3"), each = 4),
    condition = rep(c("CTRL", "disease1", "disease2"), each = 4)
  )
  rownames(meta) <- colnames(mat)

  mat_sparse <- as(mat, "dgCMatrix")
  seurat_obj <- CreateSeuratObject(counts = mat_sparse, meta.data = meta)
  seurat_obj <- NormalizeData(seurat_obj)

  # Test with minimal required parameters - ensure auto-detection works
  p <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1"
  )
  expect_s3_class(p, "ggplot")
})

test_that("DO.BarplotClustert handles various parameter combinations", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ggplot2")

  seurat_obj <- create_test_seurat()

  # Test with different parameter combinations
  p1 <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    ctrl.condition = "A",
    group.by = "condition",
    stat_pos_mod = 1.1,
    step_mod = 0.15,
    x_label_rotation = 0
  )
  expect_s3_class(p1, "ggplot")

  p2 <- DO.BarplotClustert(
    sce_object = seurat_obj,
    Feature = "gene1",
    ListTest = list(c("A","B")),
    ctrl.condition = "A",
    group.by = "condition",
    returnPlot = TRUE,
    returnValues = FALSE
  )
  expect_s3_class(p2, "ggplot")
})
