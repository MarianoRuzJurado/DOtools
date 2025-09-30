library(testthat)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

# Load example Seurat/SCE data
sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))

test_that("DO.UMAP returns a ggplot object with DimPlot", {
  p <- DO.UMAP(sce_object = sce_data, group.by = "seurat_clusters")
  expect_s3_class(p, "gg")
  expect_true(inherits(p, "ggplot"))
})

test_that("DO.UMAP returns a ggplot object with FeaturePlot", {
  genes <- c("BAG2", "CD74")
  p <- DO.UMAP(sce_object = sce_data, FeaturePlot = TRUE, features = genes)
  expect_s3_class(p, "gg")
  expect_true(inherits(p, "ggplot"))
})

test_that("DO.UMAP throws error if features are missing for FeaturePlot", {
  expect_error(DO.UMAP(sce_object = sce_data, FeaturePlot = TRUE),
               "Please provide any gene names if using FeaturePlot=TRUE")
})

test_that("DO.UMAP handles custom colors", {
  # Repeat colors to ensure enough for all clusters
  n_clusters <- length(unique(sce_data$seurat_clusters))
  custom_colors <- rep(c("blue", "forestgreen", "firebrick", "purple", "orange", "cyan", "pink", "moccasin"), length.out = n_clusters)
  p <- DO.UMAP(sce_object = sce_data, group.by = "seurat_clusters", umap_colors = custom_colors)
  expect_s3_class(p, "gg")
})

test_that("DO.UMAP respects label and plot.title arguments", {
  p <- DO.UMAP(sce_object = sce_data, label = FALSE, plot.title = FALSE)
  # Ensure no geom_text layer is present (no labels)
  expect_false(any(sapply(p$layers, function(x) inherits(x$geom, "GeomText"))))
  # Check plot title is blank
  plot_title <- ggplot_build(p)$layout$plot$title
  expect_true(is.null(plot_title) || plot_title == "")
})

test_that("DO.UMAP converts SCE to Seurat object", {
  # Ensure SCE object has UMAP reduction for plotting
  sce_obj <- as(sce_data, "SingleCellExperiment")
  p <- DO.UMAP(sce_object = sce_obj, group.by = "seurat_clusters")
  expect_s3_class(p, "gg")
})
