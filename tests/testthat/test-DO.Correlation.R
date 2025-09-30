library(testthat)
library(DOtools)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)

# ------------------------------
# SCE object for testing
# ------------------------------
counts <- matrix(
  sample(1:5, 12, replace = TRUE),
  nrow = 4
)
rownames(counts) <- paste0("Gene", 1:4)

coldata <- data.frame(
  orig.ident = c("S1","S2","S3"),
  condition = c("A","B","A"),
  row.names = paste0("Cell", 1:3)
)

sce <- SingleCellExperiment(
  assays = list(
    originalexp = counts,
    logcounts = counts
  ),
  colData = coldata
)

# ------------------------------
# Convert SCE to Seurat
# ------------------------------
seurat_obj <- as.Seurat(sce, counts = "originalexp", data = "logcounts")
DefaultAssay(seurat_obj) <- "originalexp"

# ------------------------------
# Unit test for DO.Correlation
# ------------------------------
test_that("DO.Correlation returns a ggplot object", {
  p <- DO.Correlation(
    sce_object = seurat_obj,
    group_by = "orig.ident",
    assay = "originalexp",
    features = NULL,
    method = "spearman",
    plotdesign = "square",
    plottype = "full",
    auto_limits = TRUE,
    outline.color = "white",
    colormap = c("royalblue4", "lightsteelblue", "tomato","firebrick4"),
    lab_size = 10,
    lab = TRUE,
    lab_col = "white"
  )

  expect_s3_class(p, "ggplot")
})

test_that("DO.Correlation works with minimal features and SCE input", {
  p <- DO.Correlation(
    sce_object = seurat_obj,
    group_by = "orig.ident",
    assay = "originalexp",
    features = rownames(sce)[1:2] # subset of genes
  )

  expect_s3_class(p, "ggplot")
})
