library(testthat)
library(DOtools)
library(SingleCellExperiment)
library(mockery)
library(dplyr)
library(ggplot2)

# ------------------------------
# SCE object
# ------------------------------
counts <- matrix(sample(0:5, 8, replace = TRUE), nrow = 4)
rownames(counts) <- paste0("Gene", 1:4)

coldata <- data.frame(
  orig.ident = c("S1", "S2"),
  condition = c("A", "B"),
  annotation = c("Tcell", "Bcell")
)
rownames(coldata) <- paste0("Cell", 1:2)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(RNA = counts, logcounts = counts),
  colData = coldata
)

# ------------------------------
# Unit test for DO.CellTypist
# ------------------------------
test_that("DO.CellTypist handles too few cells", {

  # Force minCellsToRun higher than number of cells
  result <- DO.CellTypist(
    sce_object = sce,
    minCellsToRun = 3
  )

  # Check that NA column was added in colData
  expect_true("predicted_labels_celltypist" %in% colnames(colData(result)))
  expect_true(all(is.na(colData(result)$predicted_labels_celltypist)))
})

test_that("DO.CellTypist returns object and prob matrix when returnProb=TRUE", {

  # Mock reading labels and probabilities
  mock_labels <- data.frame(
    majority_voting = c("Tcell", "Bcell"),
    over_clustering = c("C1", "C2"),
    row.names = colnames(sce),
    stringsAsFactors = FALSE
  )

  mock_probMatrix <- data.frame(
    cluster = c("C1","C2"),
    Tcell = c(0.8, 0.1),
    Bcell = c(0.2, 0.9)
  )

  # Stub basiliskRun and CSV reads
  stub(DO.CellTypist, 'basilisk::basiliskRun', NULL)
  stub(DO.CellTypist, 'utils::read.csv', function(file, ...) {
    if (grepl("predicted_labels.csv", file)) {
      return(mock_labels)
    } else if (grepl("probability_matrix.csv", file)) {
      return(mock_probMatrix)
    } else {
      stop("Unexpected file")
    }
  })

  result <- DO.CellTypist(
    sce_object = sce,
    minCellsToRun = 1,
    returnProb = TRUE
  )

  # Check that result is a list
  expect_type(result, "list")
  expect_named(result, c("SingleCellObject", "probMatrix"))

  # Check SCE object inside list
  expect_s4_class(result$SingleCellObject, "SingleCellExperiment")

  # Check probMatrix
  expect_s3_class(result$probMatrix, "data.frame")
  expect_true(all(c("cluster","Tcell","Bcell") %in% colnames(result$probMatrix)))
})
