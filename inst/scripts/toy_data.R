#Example data was created by using h5 files from 10x public datasets using the following links:
# healthy_link1 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5"
# healthy_link2 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_raw_feature_bc_matrix.h5"
# disease_link1 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/malt_10k_protein_v3/malt_10k_protein_v3_filtered_feature_bc_matrix.h5"
# disease_link2 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/malt_10k_protein_v3/malt_10k_protein_v3_raw_feature_bc_matrix.h5"

#Data was downloaded using .example_10x() function implemented in the DOtools package and processed with DO.Import
library(Seurat)
library(DOtools)

base <- DOtools:::.example_10x()

paths = c(file.path(base, "healthy/outs/filtered_feature_bc_matrix.h5"),
          file.path(base, "disease/outs/filtered_feature_bc_matrix.h5"))

SCE_obj <- DO.Import(pathways = paths,
                     ids = c("healthy-1", "disease-1"),
                     DeleteDoublets = TRUE,
                     cut_mt = .05,
                     min_counts = 500,
                     min_genes = 300,
                     high_quantile = .95,
                     Seurat=FALSE) # Set to TRUE for Seurat object

Seu_obj <- as.Seurat(SCE_obj)
Seu_obj[["RNA"]] <- split(Seu_obj[["RNA"]], f = Seu_obj$orig.ident)
Seu_obj <- FindVariableFeatures(Seu_obj)
Seu_obj <- ScaleData(object = Seu_obj)
Seu_obj <- RunPCA(Seu_obj, verbose = FALSE, reduction.name = "PCA")
Seu_obj <- JoinLayers(Seu_obj)
Seu_obj[["RNA"]] <- split(Seu_obj[["RNA"]], f = Seu_obj$orig.ident)

#Integration through Seurat
Seu_obj <- IntegrateLayers(object = Seu_obj,
                           method = CCAIntegration,
                           orig.reduction = "PCA",
                           new.reduction = "INTEGRATED.CCA",
                           verbose = TRUE)

Seu_obj <- JoinLayers(Seu_obj)

#afterwards object was filtered by highlyvariable features and the subsetted by 800 random genes to decrease size, including markers for immune cells
HVG <- VariableFeatures(Seu_obj)
immune_markers <- c("PTPRC", "CD79A", "BANK1", "MS4A1", "CD3E", "CD4", "IL7R", "NKG7", "KLRD1","CD68", "CD14","ITGAM", "LILRA4", "CLEC4C", "LRRC26")
comb_genes <- unique(c(HVG, immune_markers))

# Keep only genes that exist in the object
comb_genes <- intersect(comb_genes, rownames(Seu_obj))

set.seed(42)

# Limit to 800 genes (including required ones)
# ensures immune_markers are included
if (length(comb_genes) > 800) {
  # Keep immune_markers + most variable others up to 800
  extra_needed <- 800 - length(immune_markers)
  other_genes <- setdiff(HVG, immune_markers)
  selected_genes <- unique(c(immune_markers, head(other_genes, extra_needed)))
} else {
  selected_genes <- comb_genes
}

SCE_obj <- as.SingleCellExperiment(subset(Seu_obj, features = selected_genes))

saveRDS(SCE_obj, file = "/mnt/mariano/scStorage/Mariano/Rpackage/DOtools/inst/extdata", compress = "xz")

# Subset the Seurat object
sc_data <- subset(Seu_obj, features = random_genes)
