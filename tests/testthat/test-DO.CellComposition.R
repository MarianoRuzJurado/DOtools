library(testthat)
library(DOtools)
library(dplyr)
library(Seurat)
library(SingleCellExperiment)

SCE_obj <- readRDS(
  system.file("extdata", "sce_data.rds", package = "DOtools")
)

# Create a minimal mock SCE object for testing
create_minimal_mock_sce <- function() {
  set.seed(123)
  n_cells <- 30
  n_genes <- 20

  # Create mock counts
  counts <- matrix(rpois(n_genes * n_cells, lambda = 10),
                   ncol = n_cells, nrow = n_genes)
  rownames(counts) <- paste0("Gene", 1:n_genes)
  colnames(counts) <- paste0("Cell", 1:n_cells)

  col_data <- data.frame(
    annotation = sample(paste0("Cluster", 1:3), n_cells, replace = TRUE),
    orig.ident = sample(paste0("Sample", 1:2), n_cells, replace = TRUE),
    condition = sample(c("Control", "Treatment"), n_cells, replace = TRUE),
    row.names = colnames(counts)
  )

  SingleCellExperiment(
    assays = list(counts = counts),
    colData = col_data
  )
}

test_that("DO.CellComposition basic functionality with SCE object", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("reticulate")

  # Use the actual SCE object from package
  mock_sce <- SCE_obj

  # Create proper mock data that matches the function's expected structure
  mock_basilisk_result <- list(
    df_res = data.frame(
      p_values = c(0.01, 0.1, 0.001),
      adjusted_p_values = c(0.02, 0.15, 0.005),
      row.names = c("Cluster1", "Cluster2", "Cluster3")
    ),
    df_counts = data.frame(
      Cluster1 = c(10, 15, 20),
      Cluster2 = c(5, 8, 12),
      Cluster3 = c(3, 6, 9),
      row.names = c("Sample1", "Sample2", "Sample3")
    ),
    df_design = data.frame(
      Control = c(1, 0, 0),
      Treatment = c(0, 1, 1),
      row.names = c("Sample1", "Sample2", "Sample3")
    ),
    df_merge = data.frame(
      condition = c("Control", "Treatment", "Treatment"),
      Cluster1 = c(0.4, 0.35, 0.45),
      Cluster2 = c(0.3, 0.4, 0.35),
      Cluster3 = c(0.3, 0.25, 0.2)
    )
  )

  # Mock the basiliskRun function
  testthat::local_mocked_bindings(
    basiliskRun = function(env, fun, ...) {
      return(mock_basilisk_result)
    },
    .package = "basilisk"
  )

  result <- DO.CellComposition(
    sce_object = mock_sce,
    cluster_column = "annotation",
    sample_column = "orig.ident",
    condition_column = "condition",
    scanpro_plots = FALSE,
    return_df = FALSE
  )

  # Test basic return types
  expect_true(ggplot2::is.ggplot(result))
  expect_true(!is.null(result))
})

test_that("DO.CellComposition with return_df = TRUE", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("reticulate")

  mock_sce <- SCE_obj

  mock_basilisk_result <- list(
    df_res = data.frame(
      p_values = c(0.01, 0.1, 0.001),
      adjusted_p_values = c(0.02, 0.15, 0.005),
      row.names = c("Cluster1", "Cluster2", "Cluster3")
    ),
    df_counts = data.frame(
      Cluster1 = c(10, 15, 20),
      Cluster2 = c(5, 8, 12),
      Cluster3 = c(3, 6, 9),
      row.names = c("Sample1", "Sample2", "Sample3")
    ),
    df_design = data.frame(
      Control = c(1, 0, 0),
      Treatment = c(0, 1, 1),
      row.names = c("Sample1", "Sample2", "Sample3")
    ),
    df_merge = data.frame(
      condition = c("Control", "Treatment", "Treatment"),
      Cluster1 = c(0.4, 0.35, 0.45),
      Cluster2 = c(0.3, 0.4, 0.35),
      Cluster3 = c(0.3, 0.25, 0.2)
    )
  )

  testthat::local_mocked_bindings(
    basiliskRun = function(env, fun, ...) {
      return(mock_basilisk_result)
    },
    .package = "basilisk"
  )

  result <- DO.CellComposition(
    sce_object = mock_sce,
    cluster_column = "annotation",
    condition_column = "condition",
    return_df = TRUE
  )

  # Test list structure when return_df = TRUE
  expect_type(result, "list")
  expect_length(result, 2)
  expect_true(is.data.frame(result[[1]])) # prop_df
  expect_true(ggplot2::is.ggplot(result[[2]])) # plot_p
})

test_that("DO.CellComposition with sorting and subsetting options", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("reticulate")

  mock_sce <- SCE_obj

  mock_basilisk_result <- list(
    df_res = data.frame(
      p_values = c(0.01, 0.1, 0.001),
      adjusted_p_values = c(0.02, 0.15, 0.005),
      row.names = c("Cluster1", "Cluster2", "Cluster3")
    ),
    df_counts = data.frame(
      Cluster1 = c(10, 15, 20),
      Cluster2 = c(5, 8, 12),
      Cluster3 = c(3, 6, 9),
      row.names = c("Sample1", "Sample2", "Sample3")
    ),
    df_design = data.frame(
      Control = c(1, 0, 0),
      Treatment = c(0, 1, 1),
      row.names = c("Sample1", "Sample2", "Sample3")
    ),
    df_merge = data.frame(
      condition = c("Control", "Treatment", "Treatment"),
      Cluster1 = c(0.4, 0.35, 0.45),
      Cluster2 = c(0.3, 0.4, 0.35),
      Cluster3 = c(0.3, 0.25, 0.2)
    )
  )

  testthat::local_mocked_bindings(
    basiliskRun = function(env, fun, ...) {
      return(mock_basilisk_result)
    },
    .package = "basilisk"
  )

  # Test with sorting and subsetting
  result <- DO.CellComposition(
    sce_object = mock_sce,
    cluster_column = "annotation",
    condition_column = "condition",
    sort_x = c("Treatment", "Control"),
    sub_ident = c("Cluster1", "Cluster2"),
    sort_fill = c("Cluster2", "Cluster1"),
    bar_colors = c(Cluster1 = "red", Cluster2 = "blue", Cluster3 = "green")
  )

  expect_true(ggplot2::is.ggplot(result))
})

test_that("DO.CellComposition with bootstrapping", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("reticulate")

  mock_sce <- SCE_obj

  mock_basilisk_result <- list(
    df_res = data.frame(
      p_values = c(0.01, 0.1, 0.001),
      adjusted_p_values = c(0.02, 0.15, 0.005),
      row.names = c("Cluster1", "Cluster2", "Cluster3")
    ),
    df_counts = data.frame(
      Cluster1 = c(10, 15, 20),
      Cluster2 = c(5, 8, 12),
      Cluster3 = c(3, 6, 9),
      row.names = c("Sample1", "Sample2", "Sample3")
    ),
    df_design = data.frame(
      Control = c(1, 0, 0),
      Treatment = c(0, 1, 1),
      row.names = c("Sample1", "Sample2", "Sample3")
    ),
    df_merge = data.frame(
      condition = c("Control", "Treatment", "Treatment"),
      Cluster1 = c(0.4, 0.35, 0.45),
      Cluster2 = c(0.3, 0.4, 0.35),
      Cluster3 = c(0.3, 0.25, 0.2)
    )
  )

  testthat::local_mocked_bindings(
    basiliskRun = function(env, fun, ...) {
      return(mock_basilisk_result)
    },
    .package = "basilisk"
  )

  # Test with bootstrapping
  result <- DO.CellComposition(
    sce_object = mock_sce,
    cluster_column = "annotation",
    condition_column = "condition",
    n_reps = 5
  )

  expect_true(ggplot2::is.ggplot(result))
})

test_that("DO.CellComposition error conditions", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("reticulate")

  # Create a minimal SCE object without counts assay
  mock_sce_no_counts <- create_minimal_mock_sce()
  # Rename the assay to something other than "counts"
  names(assays(mock_sce_no_counts)) <- "not_counts"

  # Mock basiliskRun to not be called
  testthat::local_mocked_bindings(
    basiliskRun = function(env, fun, ...) {
      stop("basiliskRun should not be called in this test")
    },
    .package = "basilisk"
  )

  expect_error(
    DO.CellComposition(
      sce_object = mock_sce_no_counts,
      cluster_column = "annotation",
      condition_column = "condition"
    ),
    "counts not found in assays of object"
  )
})

test_that("DO.CellComposition with custom plot parameters", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("reticulate")

  mock_sce <- SCE_obj

  mock_basilisk_result <- list(
    df_res = data.frame(
      p_values = c(0.01, 0.1),
      adjusted_p_values = c(0.02, 0.15),
      row.names = c("Cluster1", "Cluster2")
    ),
    df_counts = data.frame(
      Cluster1 = c(10, 15),
      Cluster2 = c(5, 8),
      row.names = c("Sample1", "Sample2")
    ),
    df_design = data.frame(
      Control = c(1, 0),
      Treatment = c(0, 1),
      row.names = c("Sample1", "Sample2")
    ),
    df_merge = data.frame(
      condition = c("Control", "Treatment"),
      Cluster1 = c(0.4, 0.35),
      Cluster2 = c(0.6, 0.65)
    )
  )

  testthat::local_mocked_bindings(
    basiliskRun = function(env, fun, ...) {
      return(mock_basilisk_result)
    },
    .package = "basilisk"
  )

  # Test with custom plot parameters
  result <- DO.CellComposition(
    sce_object = mock_sce,
    cluster_column = "annotation",
    condition_column = "condition",
    legend.pos.x = 0.6,
    legend.pos.y = -0.1,
    cowplot_width = 0.8,
    cowlegend_width = 0.8
  )

  expect_true(ggplot2::is.ggplot(result))
})

test_that("DO.CellComposition with different column names", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("reticulate")

  mock_sce <- SCE_obj

  mock_basilisk_result <- list(
    df_res = data.frame(
      p_values = c(0.01, 0.1),
      adjusted_p_values = c(0.02, 0.15),
      row.names = c("Cluster1", "Cluster2")
    ),
    df_counts = data.frame(
      Cluster1 = c(10, 15),
      Cluster2 = c(5, 8),
      row.names = c("Sample1", "Sample2")
    ),
    df_design = data.frame(
      Control = c(1, 0),
      Treatment = c(0, 1),
      row.names = c("Sample1", "Sample2")
    ),
    df_merge = data.frame(
      condition = c("Control", "Treatment"),
      Cluster1 = c(0.4, 0.35),
      Cluster2 = c(0.6, 0.65)
    )
  )

  testthat::local_mocked_bindings(
    basiliskRun = function(env, fun, ...) {
      return(mock_basilisk_result)
    },
    .package = "basilisk"
  )

  # Test with different column mappings
  result <- DO.CellComposition(
    sce_object = mock_sce,
    cluster_column = "annotation",
    sample_column = "orig.ident",
    condition_column = "condition",
    scanpro_plots = FALSE
  )

  expect_true(ggplot2::is.ggplot(result))
})



test_that("DO.CellComposition validates required packages", {
  # This test just ensures the function contains the dependency checks

  # Check that the function source code contains the expected dependency checks
  func_source <- capture.output(DO.CellComposition)
  expect_true(any(grepl("zellkonverter", func_source)))
  expect_true(any(grepl("reticulate", func_source)))
})

# Additional edge case tests

test_that("DO.CellComposition handles scanpro_plots = TRUE gracefully", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("reticulate")

  mock_sce <- SCE_obj

  mock_basilisk_result <- list(
    df_res = data.frame(
      p_values = c(0.01, 0.1),
      adjusted_p_values = c(0.02, 0.15),
      row.names = c("Cluster1", "Cluster2")
    ),
    df_counts = data.frame(
      Cluster1 = c(10, 15),
      Cluster2 = c(5, 8),
      row.names = c("Sample1", "Sample2")
    ),
    df_design = data.frame(
      Control = c(1, 0),
      Treatment = c(0, 1),
      row.names = c("Sample1", "Sample2")
    ),
    df_merge = data.frame(
      condition = c("Control", "Treatment"),
      Cluster1 = c(0.4, 0.35),
      Cluster2 = c(0.6, 0.65)
    )
  )

  testthat::local_mocked_bindings(
    basiliskRun = function(env, fun, ...) {
      return(mock_basilisk_result)
    },
    .package = "basilisk"
  )

  # Test with scanpro_plots = TRUE but no outputFolder (should use getwd)
  result <- DO.CellComposition(
    sce_object = mock_sce,
    cluster_column = "annotation",
    condition_column = "condition",
    scanpro_plots = TRUE,
    outputFolder = tempdir()  # Use tempdir to avoid permission issues
  )

  expect_true(ggplot2::is.ggplot(result))
})

test_that("DO.CellComposition handles scanpro_group parameter", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("reticulate")

  mock_sce <- SCE_obj

  mock_basilisk_result <- list(
    df_res = data.frame(
      p_values = c(0.01, 0.1),
      adjusted_p_values = c(0.02, 0.15),
      row.names = c("Cluster1", "Cluster2")
    ),
    df_counts = data.frame(
      Cluster1 = c(10, 15),
      Cluster2 = c(5, 8),
      row.names = c("Sample1", "Sample2")
    ),
    df_design = data.frame(
      Control = c(1, 0),
      Treatment = c(0, 1),
      row.names = c("Sample1", "Sample2")
    ),
    df_merge = data.frame(
      condition = c("Control", "Treatment"),
      Cluster1 = c(0.4, 0.35),
      Cluster2 = c(0.6, 0.65)
    )
  )

  testthat::local_mocked_bindings(
    basiliskRun = function(env, fun, ...) {
      return(mock_basilisk_result)
    },
    .package = "basilisk"
  )

  # Test with scanpro_group
  result <- DO.CellComposition(
    sce_object = mock_sce,
    cluster_column = "annotation",
    condition_column = "condition",
    scanpro_plots = TRUE,
    scanpro_group = "Cluster1",
    outputFolder = tempdir()
  )

  expect_true(ggplot2::is.ggplot(result))
})
