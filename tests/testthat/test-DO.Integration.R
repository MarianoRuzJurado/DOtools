# tests/testthat/test-DO.Integration.R
library(testthat)
library(mockery)
library(Seurat)
library(SingleCellExperiment)
library(Matrix)

# small helper to create a Seurat object
make_small_seurat <- function(n_genes = 50, n_cells = 20, n_samples = 2) {
  set.seed(1)
  counts <- matrix(rpois(n_genes * n_cells, lambda = 5), nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("G", seq_len(n_genes))
  colnames(counts) <- paste0("C", seq_len(n_cells))
  so <- CreateSeuratObject(counts = counts, assay = "RNA")
  so$orig.ident <- rep(paste0("S", seq_len(n_samples)), length.out = n_cells)
  so
}

# Stubs for Seurat heavy-lifting used inside DO.Integration
stub_FindVariableFeatures <- function(object, ...) {
  VariableFeatures(object) <- head(rownames(object), 5)
  object
}
stub_ScaleData <- function(object, ...) {
  object
}
stub_RunPCA <- function(object, reduction.name = "PCA", npcs = 50, ...) {
  ncell <- ncol(object)
  emb <- matrix(rnorm(ncell * 2), nrow = ncell, ncol = 2)
  rownames(emb) <- colnames(object)
  colnames(emb) <- paste0("PC", 1:2)
  object@reductions[[reduction.name]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "PC_")
  object
}
stub_JoinLayers <- function(object, ...) object
stub_IntegrateLayers <- function(object, method, orig.reduction, new.reduction, verbose = FALSE, ...) {
  ncell <- ncol(object)
  emb <- matrix(rnorm(ncell * 2), nrow = ncell, ncol = 2)
  rownames(emb) <- colnames(object)
  colnames(emb) <- paste0("INTPC", 1:2)
  object@reductions[[new.reduction]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "INTPC_")
  object
}
stub_FindNeighbors <- function(object, reduction, dims, verbose = FALSE, ...) object
stub_FindClusters <- function(object, resolution, algorithm, random.seed, cluster.name, verbose = FALSE, ...) {
  object[[cluster.name]] <- sample(0:1, size = ncol(object), replace = TRUE)
  object
}
stub_RunUMAP <- function(object, reduction, reduction.name, dims, verbose = FALSE, ...) {
  ncell <- ncol(object)
  emb <- matrix(rnorm(ncell * 2), nrow = ncell, ncol = 2)
  rownames(emb) <- colnames(object)
  colnames(emb) <- paste0("UMAP", 1:2)
  object@reductions[[reduction.name]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "UMAP_")
  object
}

# A safe split stub used inside DO.Integration to avoid calling split.Assay5
stub_split_assay <- function(assay_obj, f, ...) {
  assay_obj
}

test_that("DO.Integration: Seurat input runs through HVG/scale/pca/integrate/neighbors/clusters/umap (stubs)", {
  so <- make_small_seurat(n_genes = 30, n_cells = 24, n_samples = 2)

  stub(DOtools::DO.Integration, "FindVariableFeatures", stub_FindVariableFeatures)
  stub(DOtools::DO.Integration, "ScaleData", stub_ScaleData)
  stub(DOtools::DO.Integration, "RunPCA", stub_RunPCA)
  stub(DOtools::DO.Integration, "JoinLayers", stub_JoinLayers)
  stub(DOtools::DO.Integration, "IntegrateLayers", stub_IntegrateLayers)
  stub(DOtools::DO.Integration, "FindNeighbors", stub_FindNeighbors)
  stub(DOtools::DO.Integration, "FindClusters", stub_FindClusters)
  stub(DOtools::DO.Integration, "RunUMAP", stub_RunUMAP)
  stub(DOtools::DO.Integration, "split", stub_split_assay)      # avoid split.Assay5 complexity
  stub(DOtools::DO.Integration, "Layers", function(x) character(0)) # force the "no layers" branch

  out <- DOtools::DO.Integration(
    sce_object = so,
    split_key = "orig.ident",
    HVG = TRUE,
    scale = TRUE,
    pca = TRUE,
    npcs = 10,
    neighbors = TRUE,
    neighbors_dim = seq_len(5),
    clusters = TRUE,
    clusters_res = 0.2,
    clusters_algorithm = 4,
    umap = TRUE,
    umap_key = "UMAPTEST",
    umap_dim = seq_len(5),
    integration_method = "CCAIntegration",
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_gt(length(VariableFeatures(out)), 0)
  expect_true("PCA" %in% names(out@reductions))
  expect_true("INTEGRATED.CCA" %in% names(out@reductions))
  expect_true("UMAPTEST" %in% names(out@reductions))
  expected_cl_name <- paste0("leiden", 0.2)
  expect_true(expected_cl_name %in% colnames(out@meta.data))
  expect_false("seurat_clusters" %in% colnames(out@meta.data))
})

test_that("DO.Integration: accepts SingleCellExperiment input and returns SCE (round-trip)", {
  # create SCE with a 'logcounts' assay so as.Seurat conversion is safe
  set.seed(2)
  mat <- matrix(rpois(20 * 6, lambda = 4), nrow = 20, ncol = 6)
  rownames(mat) <- paste0("g", seq_len(20))
  colnames(mat) <- paste0("c", seq_len(6))
  sce <- SingleCellExperiment(assays = list(counts = mat, logcounts = log1p(mat)))
  sce$orig.ident <- rep(c("A", "B"), length.out = ncol(sce))

  # Build a simple Seurat object explicitly (ensures an "RNA" assay exists)
  seurat_from_sce <- CreateSeuratObject(counts = assay(sce, "counts"), assay = "RNA")
  seurat_from_sce$orig.ident <- sce$orig.ident

  # stub conversion helper inside DO.Integration to return our pre-built Seurat
  stub(DOtools::DO.Integration, ".suppressDeprecationWarnings", function(x) seurat_from_sce)

  # stub other heavy operations and split to be safe
  stub(DOtools::DO.Integration, "FindVariableFeatures", stub_FindVariableFeatures)
  stub(DOtools::DO.Integration, "ScaleData", stub_ScaleData)
  stub(DOtools::DO.Integration, "RunPCA", stub_RunPCA)
  stub(DOtools::DO.Integration, "JoinLayers", stub_JoinLayers)
  stub(DOtools::DO.Integration, "IntegrateLayers", stub_IntegrateLayers)
  stub(DOtools::DO.Integration, "FindNeighbors", stub_FindNeighbors)
  stub(DOtools::DO.Integration, "FindClusters", stub_FindClusters)
  stub(DOtools::DO.Integration, "RunUMAP", stub_RunUMAP)
  stub(DOtools::DO.Integration, "split", stub_split_assay)
  stub(DOtools::DO.Integration, "Layers", function(x) character(0))

  out_sce <- DOtools::DO.Integration(
    sce_object = sce,
    split_key = "orig.ident",
    HVG = TRUE,
    scale = TRUE,
    pca = TRUE,
    neighbors = FALSE,
    clusters = FALSE,
    umap = FALSE,
    verbose = FALSE
  )

  expect_s4_class(out_sce, "SingleCellExperiment")
  expect_true("counts" %in% names(assays(out_sce)))
  expect_true("orig.ident" %in% colnames(colData(out_sce)))
})

test_that("DO.Integration: skip split branch when Layers indicates existing layers", {
  so <- make_small_seurat(n_genes = 40, n_cells = 12, n_samples = 3)

  # stub Layers to return a matching-layer name so split is skipped
  stub(DOtools::DO.Integration, "Layers", function(x) c("counts.S1"))
  stub(DOtools::DO.Integration, "JoinLayers", stub_JoinLayers)
  stub(DOtools::DO.Integration, "IntegrateLayers", stub_IntegrateLayers)
  stub(DOtools::DO.Integration, "split", stub_split_assay)

  out <- DOtools::DO.Integration(
    sce_object = so,
    split_key = "orig.ident",
    HVG = FALSE,
    scale = FALSE,
    pca = FALSE,
    neighbors = FALSE,
    clusters = FALSE,
    umap = FALSE,
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_true("INTEGRATED.CCA" %in% names(out@reductions))
})

test_that("DO.Integration: selection_method/integration_method options run (smoke)", {
  so <- make_small_seurat(n_genes = 25, n_cells = 10, n_samples = 2)

  stub(DOtools::DO.Integration, "FindVariableFeatures", stub_FindVariableFeatures)
  stub(DOtools::DO.Integration, "IntegrateLayers", stub_IntegrateLayers)
  stub(DOtools::DO.Integration, "JoinLayers", stub_JoinLayers)
  stub(DOtools::DO.Integration, "Layers", function(x) character(0))
  stub(DOtools::DO.Integration, "split", stub_split_assay)

  out <- DOtools::DO.Integration(
    sce_object = so,
    split_key = "orig.ident",
    HVG = TRUE,
    selection_method = "mean.var.plot",
    integration_method = "CCAIntegration",
    scale = FALSE,
    pca = FALSE,
    neighbors = FALSE,
    clusters = FALSE,
    umap = FALSE,
    verbose = FALSE
  )

  expect_s4_class(out, "Seurat")
  expect_true(length(VariableFeatures(out)) > 0)
  expect_true("INTEGRATED.CCA" %in% names(out@reductions))
})
