library(testthat)
library(mockery)
library(Seurat)
library(SingleCellExperiment)
library(Matrix)

SCE_obj <- readRDS(
  system.file("extdata",
              "sce_data.rds",
              package = "DOtools")
)

# ------------------------------
# Improved Helper functions for test objects
# ------------------------------
setup_minimal_seurat <- function() {
  set.seed(42)
  counts_matrix <- as(
    matrix(rpois(20*10, lambda = 10), nrow = 20, ncol = 10), "dgCMatrix"
    )
  rownames(counts_matrix) <- paste0("Gene", 1:20)
  colnames(counts_matrix) <- paste0("Cell", 1:10)

  seurat_obj <- CreateSeuratObject(counts = counts_matrix)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

  seurat_obj$orig.ident <- rep(c("A", "B"), each = 5)
  VariableFeatures(seurat_obj) <- rownames(seurat_obj)[1:10]
  seurat_obj
}

setup_minimal_sce <- function() {
  set.seed(42)
  counts_matrix <- as(
    matrix(rpois(20*10, lambda = 10), nrow = 20, ncol = 10), "dgCMatrix"
    )
  rownames(counts_matrix) <- paste0("Gene", 1:20)
  colnames(counts_matrix) <- paste0("Cell", 1:10)

  sce <- SingleCellExperiment(
    assays = list(
      counts = counts_matrix,
      logcounts = log1p(counts_matrix)
    )
  )
  colData(sce)$orig.ident <- rep(c("A", "B"), each = 5)
  sce
}

# Global mock for basiliskRun to avoid Python execution
mock_basilisk_run <- function(env, fun, args) {
  # Return a simple matrix as the scVI embedding
  embedding <- matrix(
    runif(ncol(args$sce_object) * 5), nrow = ncol(args$sce_object), ncol = 5
    )
  rownames(embedding) <- colnames(args$sce_object)
  embedding
}

# ------------------------------
# Enhanced Tests with Real Data
# ------------------------------

test_that("DO.scVI works with real SCE data", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  # Ensure the object has the required structure
  expect_s4_class(SCE_obj, "SingleCellExperiment")
  expect_true("counts" %in% assayNames(SCE_obj))

  # Add batch information if not present
  if (!"orig.ident" %in% colnames(colData(SCE_obj))) {
    colData(SCE_obj)$orig.ident <- rep(c("A", "B"), length.out = ncol(SCE_obj))
  }

  mock_as_seurat <- function(x) {
    # Create a minimal Seurat object that mimics conversion from real SCE
    counts_matrix <- assay(x, "counts")
    seurat_obj <- CreateSeuratObject(counts = counts_matrix)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    VariableFeatures(seurat_obj) <-
      rownames(seurat_obj)[1:min(2000, nrow(seurat_obj))]
    seurat_obj$orig.ident <- colData(x)$orig.ident
    seurat_obj
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      result <- DO.scVI(SCE_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "SingleCellExperiment")
  expect_true("scVI" %in% reducedDimNames(result))
})

test_that("DO.scVI handles real data with custom batch keys", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  # Create a custom batch column
  colData(SCE_obj)$custom_batch <- sample(c("Batch1", "Batch2", "Batch3"),
                                          ncol(SCE_obj), replace = TRUE)

  mock_as_seurat <- function(x) {
    counts_matrix <- assay(x, "counts")
    seurat_obj <- CreateSeuratObject(counts = counts_matrix)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    VariableFeatures(seurat_obj) <-
      rownames(seurat_obj)[1:min(2000, nrow(seurat_obj))]
    seurat_obj$custom_batch <- colData(x)$custom_batch
    seurat_obj
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      result <- DO.scVI(SCE_obj, batch_key = "custom_batch")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "SingleCellExperiment")
})

test_that("DO.scVI handles real data with covariates", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  # Add some mock covariates
  colData(SCE_obj)$batch <- rep(c("A", "B"), length.out = ncol(SCE_obj))
  colData(SCE_obj)$condition <- sample(
    c("Ctrl", "Treat"), ncol(SCE_obj), replace = TRUE
    )
  colData(SCE_obj)$n_genes <- rnorm(ncol(SCE_obj))

  mock_as_seurat <- function(x) {
    counts_matrix <- assay(x, "counts")
    seurat_obj <- CreateSeuratObject(counts = counts_matrix)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    VariableFeatures(seurat_obj) <-
      rownames(seurat_obj)[1:min(2000, nrow(seurat_obj))]
    seurat_obj$batch <- colData(x)$batch
    seurat_obj$condition <- colData(x)$condition
    seurat_obj$n_genes <- colData(x)$n_genes
    seurat_obj
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      result <- DO.scVI(SCE_obj,
                        batch_key = "batch",
                        categorical_covariates = "condition",
                        continuos_covariates = "n_genes")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "SingleCellExperiment")
})

test_that("DO.scVI handles real data with different parameters", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  colData(SCE_obj)$batch <- rep(c("A", "B"), length.out = ncol(SCE_obj))

  mock_as_seurat <- function(x) {
    counts_matrix <- assay(x, "counts")
    seurat_obj <- CreateSeuratObject(counts = counts_matrix)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

    VariableFeatures(seurat_obj) <-
      rownames(seurat_obj)[1:min(2000, nrow(seurat_obj))]

    seurat_obj$batch <- colData(x)$batch
    seurat_obj
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      result <- DO.scVI(SCE_obj,
                        batch_key = "batch",
                        n_latent = 20,
                        n_hidden = 256,
                        n_layers = 2,
                        dispersion = "gene",
                        gene_likelihood = "nb")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "SingleCellExperiment")
})

test_that("DO.scVI converts real Seurat data correctly", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  # Convert SCE to Seurat for testing
  counts_matrix <- assay(SCE_obj, "counts")
  seurat_real <- CreateSeuratObject(counts = counts_matrix)
  seurat_real <- NormalizeData(seurat_real, verbose = FALSE)
  seurat_real <- ScaleData(seurat_real, verbose = FALSE)

  # Add batch information
  seurat_real$batch <- rep(c("A", "B"), length.out = ncol(seurat_real))

  VariableFeatures(seurat_real) <-
    rownames(seurat_real)[1:min(2000, nrow(seurat_real))]

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      result <- DO.scVI(seurat_real, batch_key = "batch")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "Seurat")
  expect_true("scVI" %in% names(result@reductions))
})

# Test custom layer parameters properly
test_that("DO.scVI handles custom layer parameters", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  colData(SCE_obj)$batch <- rep(c("A", "B"), length.out = ncol(SCE_obj))

  mock_as_seurat <- function(x) {
    counts_matrix <- assay(x, "counts")
    seurat_obj <- CreateSeuratObject(counts = counts_matrix)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

    VariableFeatures(seurat_obj) <-
      rownames(seurat_obj)[1:min(2000, nrow(seurat_obj))]

    seurat_obj$batch <- colData(x)$batch
    seurat_obj
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      # Test that the function accepts custom layer parameters without error
      result <- DO.scVI(SCE_obj,
                        batch_key = "batch",
                        layer_counts = "counts",  # Use the actual counts layer
                        # Use the actual logcounts layer
                        layer_logcounts = "logcounts")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "SingleCellExperiment")
})

# Test edge cases with real data
test_that("DO.scVI handles real data edge cases", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  # Test with a subset of the real data
  small_sce <- SCE_obj[1:min(100, nrow(SCE_obj)), 1:min(50, ncol(SCE_obj))]
  colData(small_sce)$batch <- rep(c("A", "B"), length.out = ncol(small_sce))

  mock_as_seurat <- function(x) {
    counts_matrix <- assay(x, "counts")
    seurat_obj <- CreateSeuratObject(counts = counts_matrix)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

    VariableFeatures(seurat_obj) <-
      rownames(seurat_obj)[1:min(10, nrow(seurat_obj))]

    seurat_obj$batch <- colData(x)$batch
    seurat_obj
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      result <- DO.scVI(small_sce, batch_key = "batch")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "SingleCellExperiment")
})

# Test error conditions with real data
test_that("DO.scVI handles real data with missing counts", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  # Create a copy without counts to test error handling
  sce_no_counts <- SCE_obj

  # Remove counts assay properly
  if ("counts" %in% assayNames(sce_no_counts)) {
    # Keep only non-counts assays
    non_counts_assays <- setdiff(assayNames(sce_no_counts), "counts")
    if (length(non_counts_assays) > 0) {
      assays(sce_no_counts) <- assays(sce_no_counts)[non_counts_assays]
    } else {
      # If no other assays, create a dummy assay without counts
      dummy_matrix <- matrix(1, nrow = nrow(sce_no_counts),
                             ncol = ncol(sce_no_counts))
      rownames(dummy_matrix) <- rownames(sce_no_counts)
      colnames(dummy_matrix) <- colnames(sce_no_counts)
      assays(sce_no_counts) <- list(dummy = dummy_matrix)
    }
  }

  colData(sce_no_counts)$batch <- rep(c("A", "B"),
                                      length.out = ncol(sce_no_counts))

  mock_as_seurat <- function(x) {
    seurat_obj <- setup_minimal_seurat()
    seurat_obj$batch <- colData(x)$batch
    seurat_obj
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      # This should trigger the "counts not found" error in function
      expect_error(DO.scVI(sce_no_counts, batch_key = "batch"),
                   "counts not found in assays of object!")
    },
    .package = "basilisk"
  )
})

# Test that DO.scVI preserves cell identities with real data
test_that("DO.scVI preserves cell identities with real data", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  colData(SCE_obj)$batch <- rep(c("A", "B"), length.out = ncol(SCE_obj))
  original_cells <- colnames(SCE_obj)

  mock_as_seurat <- function(x) {
    counts_matrix <- assay(x, "counts")
    seurat_obj <- CreateSeuratObject(counts = counts_matrix)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    VariableFeatures(seurat_obj) <-
      rownames(seurat_obj)[1:min(2000, nrow(seurat_obj))]
    seurat_obj$batch <- colData(x)$batch
    seurat_obj
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      result <- DO.scVI(SCE_obj, batch_key = "batch")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "SingleCellExperiment")
  expect_equal(colnames(result), original_cells)
})

# Test that DO.scVI handles real data with specific HVG subsets
test_that("DO.scVI handles real data with specific HVG", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  colData(SCE_obj)$batch <- rep(c("A", "B"), length.out = ncol(SCE_obj))

  # Convert to Seurat to control HVG
  counts_matrix <- assay(SCE_obj, "counts")
  seurat_real <- CreateSeuratObject(counts = counts_matrix)
  seurat_real <- NormalizeData(seurat_real, verbose = FALSE)
  seurat_real <- ScaleData(seurat_real, verbose = FALSE)
  seurat_real$batch <- colData(SCE_obj)$batch

  # Set specific HVG
  n_genes <- min(500, nrow(seurat_real))
  VariableFeatures(seurat_real) <- rownames(seurat_real)[1:n_genes]

  captured_n_genes <- NULL
  hvg_mock <- function(env, fun, args) {
    captured_n_genes <<- nrow(args$sce_object)
    embedding <- matrix(
      runif(ncol(args$sce_object) * 5), nrow = ncol(args$sce_object), ncol = 5
      )
    rownames(embedding) <- colnames(args$sce_object)
    embedding
  }

  testthat::with_mocked_bindings(
    basiliskRun = hvg_mock,
    code = {
      result <- DO.scVI(seurat_real, batch_key = "batch")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "Seurat")
  expect_equal(captured_n_genes, n_genes)
})

# Test the actual logger function in code
test_that("DO.scVI logger works with real data", {
  skip_if_not(exists("SCE_obj"), "Real SCE data not available")

  colData(SCE_obj)$batch <- rep(c("A", "B"), length.out = ncol(SCE_obj))

  mock_as_seurat <- function(x) {
    counts_matrix <- assay(x, "counts")
    seurat_obj <- CreateSeuratObject(counts = counts_matrix)
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    # Don't set VariableFeatures to trigger HVG calculation message
    seurat_obj$batch <- colData(x)$batch
    seurat_obj
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      result <- DO.scVI(SCE_obj, batch_key = "batch")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "SingleCellExperiment")
})

# Test that DO.scVI correctly uses the environment
test_that("DO.scVI uses correct environment", {
  seurat_obj <- setup_minimal_seurat()

  captured_env <- NULL
  env_mock <- function(env, fun, args) {
    captured_env <<- env
    embedding <- matrix(
      runif(ncol(args$sce_object) * 5), nrow = ncol(args$sce_object), ncol = 5
      )
    rownames(embedding) <- colnames(args$sce_object)
    embedding
  }

  testthat::with_mocked_bindings(
    basiliskRun = env_mock,
    code = {
      result <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_true(!is.null(captured_env))
})

# ------------------------------
# Existing Tests
# ------------------------------

test_that("DO.scVI runs on Seurat object and adds reduction", {
  seurat_obj <- setup_minimal_seurat()

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      result <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "Seurat")
  expect_true("scVI" %in% names(result@reductions))
  expect_equal(ncol(result@reductions$scVI@cell.embeddings), 5)
})

test_that("DO.scVI runs on SingleCellExperiment object", {
  sce_obj <- setup_minimal_sce()

  mock_as_seurat <- function(x) {
    setup_minimal_seurat()
  }

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      result <- DO.scVI(sce_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "SingleCellExperiment")
  expect_true("scVI" %in% reducedDimNames(result))
})

test_that("DO.scVI calculates HVG if none found in Seurat object", {
  seurat_obj <- setup_minimal_seurat()
  VariableFeatures(seurat_obj) <- character(0)

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      result <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "Seurat")
  expect_true("scVI" %in% names(result@reductions))
  expect_true(length(VariableFeatures(result)) > 0)
})

test_that("DO.scVI correctly subsets to HVG", {
  seurat_obj <- setup_minimal_seurat()
  VariableFeatures(seurat_obj) <- rownames(seurat_obj)[1:5]

  captured_sce <- NULL
  custom_mock <- function(env, fun, args) {
    captured_sce <<- args$sce_object
    embedding <- matrix(
      runif(ncol(args$sce_object) * 5), nrow = ncol(args$sce_object), ncol = 5
      )
    rownames(embedding) <- colnames(args$sce_object)
    embedding
  }

  testthat::with_mocked_bindings(
    basiliskRun = custom_mock,
    code = {
      result <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_equal(nrow(captured_sce), 5)
  expect_true(all(rownames(captured_sce) %in% VariableFeatures(seurat_obj)))
})

test_that("DO.scVI preserves original object structure", {
  seurat_obj <- setup_minimal_seurat()
  seurat_obj$custom_meta <- sample(
    c("X", "Y"), ncol(seurat_obj), replace = TRUE
    )

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      result <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "Seurat")
  expect_true("custom_meta" %in% colnames(result@meta.data))
  expect_true("scVI" %in% names(result@reductions))
})

test_that("DO.scVI creates correct reduction keys", {
  seurat_obj <- setup_minimal_seurat()

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      result <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_equal(result@reductions$scVI@key, "scVI_")
  expect_equal(DefaultAssay(result@reductions$scVI), "RNA")
})

test_that("DO.scVI handles different parameter values", {
  seurat_obj <- setup_minimal_seurat()

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      result1 <- DO.scVI(seurat_obj, batch_key = "orig.ident", n_latent = 10)
      result2 <- DO.scVI(seurat_obj, batch_key = "orig.ident", n_hidden = 64)
    },
    .package = "basilisk"
  )

  expect_s4_class(result1, "Seurat")
  expect_s4_class(result2, "Seurat")
})

test_that("DO.scVI correctly handles object class detection", {
  seurat_obj <- setup_minimal_seurat()

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      result_seurat <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )
  expect_s4_class(result_seurat, "Seurat")

  sce_obj <- setup_minimal_sce()

  mock_as_seurat <- function(x) setup_minimal_seurat()

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      mockery::stub(DO.scVI, "as.Seurat", mock_as_seurat)
      result_sce <- DO.scVI(sce_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )
  expect_s4_class(result_sce, "SingleCellExperiment")
})

test_that("DO.scVI handles small dataset", {
  set.seed(42)
  counts_matrix <- as(
    matrix(rpois(5*5, lambda = 10), nrow = 5, ncol = 5), "dgCMatrix"
    )
  rownames(counts_matrix) <- paste0("Gene", 1:5)
  colnames(counts_matrix) <- paste0("Cell", 1:5)

  seurat_obj <- CreateSeuratObject(counts = counts_matrix)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj$orig.ident <- c("A", "B", "A", "B", "A")
  VariableFeatures(seurat_obj) <- rownames(seurat_obj)[1:3]

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      result <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "Seurat")
  expect_true("scVI" %in% names(result@reductions))
})

test_that("DO.scVI basic functionality works", {
  seurat_obj <- setup_minimal_seurat()

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      result <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_s4_class(result, "Seurat")
  expect_true("scVI" %in% names(result@reductions))
})

test_that("DO.scVI passes parameters correctly", {
  seurat_obj <- setup_minimal_seurat()

  captured_args <- NULL
  param_mock <- function(env, fun, args) {
    captured_args <<- args
    embedding <- matrix(
      runif(ncol(args$sce_object) * 5), nrow = ncol(args$sce_object), ncol = 5
      )
    rownames(embedding) <- colnames(args$sce_object)
    embedding
  }

  testthat::with_mocked_bindings(
    basiliskRun = param_mock,
    code = {
      result <- DO.scVI(
        seurat_obj,
        batch_key = "test_batch",
        n_latent = 15,
        n_hidden = 64,
        n_layers = 2
      )
    },
    .package = "basilisk"
  )

  expect_equal(captured_args$batch_key, "test_batch")
  expect_equal(captured_args$n_latent, 15L)
  expect_equal(captured_args$n_hidden, 64L)
  expect_equal(captured_args$n_layers, 2L)
})

test_that("DO.scVI subsets to HVG", {
  seurat_obj <- setup_minimal_seurat()
  VariableFeatures(seurat_obj) <- rownames(seurat_obj)[1:5]

  captured_nrow <- NULL
  subset_mock <- function(env, fun, args) {
    captured_nrow <<- nrow(args$sce_object)
    embedding <- matrix(
      runif(ncol(args$sce_object) * 5), nrow = ncol(args$sce_object), ncol = 5
      )
    rownames(embedding) <- colnames(args$sce_object)
    embedding
  }

  testthat::with_mocked_bindings(
    basiliskRun = subset_mock,
    code = {
      result <- DO.scVI(seurat_obj, batch_key = "orig.ident")
    },
    .package = "basilisk"
  )

  expect_equal(captured_nrow, 5)
})

test_that("DO.scVI handles missing batch_key", {
  seurat_obj <- setup_minimal_seurat()

  testthat::with_mocked_bindings(
    basiliskRun = mock_basilisk_run,
    code = {
      expect_error(DO.scVI(seurat_obj), "argument \"batch_key\" is missing")
    },
    .package = "basilisk"
  )
})

test_that("DO.scVI handles invalid object types", {
  invalid_obj <- data.frame(x = 1:10)

  expect_error(DO.scVI(invalid_obj, batch_key = "test"),
               "no applicable method for 'as.Seurat'")
})
