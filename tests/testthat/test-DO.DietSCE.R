library(testthat)
library(DOtools)
library(SingleCellExperiment)
library(Seurat)

# ------------------------------
# Seurat object for testing
# ------------------------------
counts <- matrix(
  sample(1:5, 12, replace = TRUE),
  nrow = 4,
  dimnames = list(paste0("Gene", 1:4), paste0("Cell", 1:3))
)

coldata <- data.frame(
  orig.ident = c("S1","S2","S3"),
  condition = c("A","B","A"),
  row.names = paste0("Cell", 1:3)
)

sce <- SingleCellExperiment(
  assays = list(
    counts = counts,
    logcounts = counts
  ),
  colData = coldata
)

# ------------------------------
# Create Seurat object safely without touching layers
# ------------------------------
seurat_obj <- CreateSeuratObject(counts = counts)

# ------------------------------
# Unit tests for DO.DietSCE
# ------------------------------
test_that("DO.DietSCE returns Seurat object unchanged when no layers exist", {
  obj <- DO.DietSCE(seurat_obj, assay = "RNA", pattern = "^scale\\.data\\.")
  expect_s4_class(obj, "Seurat")
})

test_that("DO.DietSCE works with SCE input and no matching layers", {
  obj <- DO.DietSCE(sce, assay = "counts", pattern = "^none$")
  expect_s4_class(obj, "SingleCellExperiment")
})
