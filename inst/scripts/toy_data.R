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

Seu_obj <- DO.Import(pathways = paths,
                     ids = c("healthy-1", "disease-1"),
                     DeleteDoublets = TRUE,
                     cut_mt = .05,
                     min_counts = 500,
                     min_genes = 300,
                     high_quantile = .95)

#Integration through Seurat
Seu_obj <- IntegrateLayers(object = Seu_obj,
                           method = CCAIntegration,
                           orig.reduction = "pca",
                           new.reduction = "integrated.cca",
                           verbose = TRUE)

Seu_obj <- JoinLayers(Seu_obj)

#afterwards object was filtered by highlyvariable features and the subsetted by 800 random genes to decrease size
HVG <- VariableFeatures(Seu_obj)
Seu_obj <- subset(Seu_obj, features = HVG)

#all gene names in the Seurat object
all_genes <- rownames(Seu_obj)

#sample 800 genes
set.seed(42)
n_genes <- 800
random_genes <- sample(all_genes, n_genes)

# Subset the Seurat object
sc_data <- subset(Seu_obj, features = random_genes)
