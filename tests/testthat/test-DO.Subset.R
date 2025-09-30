# tests/testthat/test-DO.Subset-extended.R
library(testthat)
library(DOtools)
library(Seurat)
library(SingleCellExperiment)
library(Matrix)

# load example objects
sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
sce_seurat <- as.Seurat(sce_data)

# ensure data slot exists
if (!"data" %in% slotNames(sce_seurat[["RNA"]])) {
  sce_seurat[["RNA"]]@data <- as.matrix(GetAssayData(sce_seurat, layer = "counts"))
}

# ensure condition present
if (!("condition" %in% colnames(sce_seurat@meta.data))) {
  sce_seurat@meta.data$condition <- sample(c("healthy", "disease"), ncol(sce_seurat), replace = TRUE)
}

ident_vec_meta <- sce_seurat@meta.data$nCount_RNA
names(ident_vec_meta) <- colnames(sce_seurat)

# helper to compute expected names for threshold tests
compute_expected <- function(ident_vec, ident_thresh) {
  operator <- gsub("[0-9.]", "", ident_thresh)
  threshold <- as.numeric(gsub("[^0-9.]", "", ident_thresh))
  if (length(operator) == 1) {
    if (operator == "<") {
      idx <- which(ident_vec < threshold)
    } else if (operator == ">") {
      idx <- which(ident_vec > threshold)
    } else if (operator == "<=") {
      idx <- which(ident_vec <= threshold)
    } else if (operator == ">=") {
      idx <- which(ident_vec >= threshold)
    } else {
      stop("Invalid operator in test helper")
    }
  } else if (length(operator) == 2) {
    op_pair <- paste0(operator, collapse = "")
    if (op_pair == "><") {
      idx <- which(ident_vec > threshold[1] & ident_vec < threshold[2])
    } else if (op_pair == "<>") {
      idx <- which(ident_vec < threshold[1] & ident_vec > threshold[2])
    } else {
      stop("Invalid operator pair in test helper")
    }
  } else {
    stop("Unsupported operator format in test helper")
  }
  names(ident_vec)[idx]
}

test_that("DO.Subset subsets by ident_name correctly (existing test kept)", {
  ident <- "condition"
  ident_name <- unique(sce_seurat@meta.data[[ident]])[1]
  sub_obj <- DO.Subset(
    sce_object = sce_seurat,
    ident = ident,
    ident_name = ident_name
  )
  expected_cells <- colnames(sce_seurat)[sce_seurat@meta.data[[ident]] %in% ident_name]
  expect_s4_class(sub_obj, "Seurat")
  expect_equal(sort(colnames(sub_obj)), sort(expected_cells))
})

test_that("DO.Subset subsets by numeric threshold '>' correctly", {
  ident <- "nCount_RNA"
  thresh_val <- floor(median(ident_vec_meta, na.rm = TRUE))
  ident_thresh <- paste0(">", thresh_val)

  expected_cells <- compute_expected(ident_vec_meta, ident_thresh)

  sub_obj <- DO.Subset(
    sce_object = sce_seurat,
    ident = ident,
    ident_thresh = ident_thresh
  )

  expect_s4_class(sub_obj, "Seurat")
  expect_equal(sort(colnames(sub_obj)), sort(expected_cells))
})

test_that("DO.Subset subsets by numeric threshold '<' correctly", {
  ident <- "nCount_RNA"
  thresh_val <- ceiling(median(ident_vec_meta, na.rm = TRUE))
  ident_thresh <- paste0("<", thresh_val)

  expected_cells <- compute_expected(ident_vec_meta, ident_thresh)

  sub_obj <- DO.Subset(
    sce_object = sce_seurat,
    ident = ident,
    ident_thresh = ident_thresh
  )

  expect_s4_class(sub_obj, "Seurat")
  expect_equal(sort(colnames(sub_obj)), sort(expected_cells))
})

test_that("DO.Subset supports '<=' and '>=' operators", {
  ident <- "nCount_RNA"
  low <- quantile(ident_vec_meta, probs = 0.25, na.rm = TRUE)
  high <- quantile(ident_vec_meta, probs = 0.75, na.rm = TRUE)

  # <=
  ident_thresh_le <- paste0("<=", floor(low))
  expected_le <- compute_expected(ident_vec_meta, ident_thresh_le)
  sub_le <- DO.Subset(sce_object = sce_seurat, ident = ident, ident_thresh = ident_thresh_le)
  expect_equal(sort(colnames(sub_le)), sort(expected_le))

  # >=
  ident_thresh_ge <- paste0(">=", ceiling(high))
  expected_ge <- compute_expected(ident_vec_meta, ident_thresh_ge)
  sub_ge <- DO.Subset(sce_object = sce_seurat, ident = ident, ident_thresh = ident_thresh_ge)
  expect_equal(sort(colnames(sub_ge)), sort(expected_ge))
})

test_that("DO.Subset supports two-operator range '><' (greater than then less than)", {
  ident <- "nCount_RNA"
  low <- as.numeric(quantile(ident_vec_meta, 0.25, na.rm = TRUE))
  high <- as.numeric(quantile(ident_vec_meta, 0.75, na.rm = TRUE))
  ident_thresh <- c(paste0(">", floor(low)), paste0("<", ceiling(high)))

  expected_cells <- compute_expected(ident_vec_meta, ident_thresh)

  sub_obj <- DO.Subset(
    sce_object = sce_seurat,
    ident = ident,
    ident_thresh = ident_thresh
  )
  expect_s4_class(sub_obj, "Seurat")
  expect_equal(sort(colnames(sub_obj)), sort(expected_cells))
})

test_that("DO.Subset supports two-operator '<>' (less than then greater than)", {
  ident <- "nCount_RNA"
  low <- as.numeric(quantile(ident_vec_meta, 0.1, na.rm = TRUE))
  high <- as.numeric(quantile(ident_vec_meta, 0.9, na.rm = TRUE))
  # choose thresholds where some cells satisfy < high & > low
  ident_thresh <- c(paste0("<", ceiling(high)), paste0(">", floor(low)))

  expected_cells <- compute_expected(ident_vec_meta, ident_thresh)

  sub_obj <- DO.Subset(
    sce_object = sce_seurat,
    ident = ident,
    ident_thresh = ident_thresh
  )
  expect_s4_class(sub_obj, "Seurat")
  expect_equal(sort(colnames(sub_obj)), sort(expected_cells))
})

test_that("DO.Subset throws error when both ident_name and ident_thresh are provided", {
  expect_error(
    DO.Subset(
      sce_object = sce_seurat,
      ident = "condition",
      ident_name = unique(sce_seurat@meta.data$condition)[1],
      ident_thresh = ">500"
    ),
    "Please provide ident_name for subsetting by a name in the column or ident_thresh if it by a threshold"
  )
})

test_that("DO.Subset throws error when no cells left", {
  expect_error(
    DO.Subset(
      sce_object = sce_seurat,
      ident = "nCount_RNA",
      ident_thresh = ">999999"
    ),
    "No cells left after subsetting!"
  )
})

test_that("DO.Subset errors on invalid operator", {
  expect_error(
    DO.Subset(
      sce_object = sce_seurat,
      ident = "nCount_RNA",
      ident_thresh = "!=5"
    ),
    "Invalid threshold operator provided"
  )
})

test_that("DO.Subset works with SingleCellExperiment input and returns SCE", {
  ident <- "condition"
  ident_name <- unique(colData(sce_data)[[ident]])[1]
  sub_sce <- DO.Subset(
    sce_object = sce_data,
    ident = ident,
    ident_name = ident_name
  )
  expect_s4_class(sub_sce, "SingleCellExperiment")
  expect_true(all(colData(sub_sce)[[ident]] %in% ident_name))
})

test_that("DO.Subset preserves reduction names when input is Seurat (reduction name transfer)", {
  # Add a fake reduction to original seurat
  so <- sce_seurat
  emb <- matrix(rnorm(ncol(so) * 2), ncol = 2)
  rownames(emb) <- colnames(so)
  so@reductions[["MYRED"]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "MYR_")
  original_names <- names(so@reductions)

  ident <- "condition"
  ident_name <- unique(so@meta.data[[ident]])[1]
  sub_obj <- DO.Subset(
    sce_object = so,
    ident = ident,
    ident_name = ident_name
  )

  # After subsetting the function sets names(sce_object_sub@reductions) <- reduction_names
  expect_equal(names(sub_obj@reductions), original_names)
})
