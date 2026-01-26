#' @author Mariano Ruz Jurado
#' @title DO.Integration SCE object integration directly
#' @description Integrates single-cell RNA-seq data directly from
#' SingleCellExperiment or Seurat objects. Supports detection of variable genes
#' , scaling, PCA, neighbor graph construction, clustering, and UMAP embedding,
#' with multiple integration methods.
#' @param sce_object Seurat or SCE Object
#' @param split_key Character. Column in meta data to split the samples by,
#' default orig.ident
#' @param HVG Logical. Perform detection of highly variable genes
#' @param scale Logical. Perform scaling of the expression data
#' @param pca Logical. Perform principal component analysis
#' @param neighbors Logical. Perform Nearest-neighbor graph after integration
#' @param neighbors_dim Numeric range. Dimensions of reduction to use as input
#' @param clusters Logical. Perform clustering of cells
#' @param clusters_res Numeric. Value of the resolution parameter, use a value
#' above (below) 1.0 if you want to obtain a larger (smaller) number of
#' communities.
#' @param clusters_algorithm Numeric. Define the algorithm for clustering,
#' default 4 for "Leiden"
#' @param umap Logical. Runs the Uniform Manifold Approximation and Projection
#' @param umap_key Character name for
#' @param umap_dim Numeric range. Which dimensions to use as input features
#' @param pca Logical. Perform principal component analysis
#' @param integration_method Character. Define the integration method, please
#' check what versions are supported in Seurat::IntegrateLayers function
#' @param selection_method Character. Default "vst". Options: "mean.var.plot",
#' "dispersion"
#' @param loess_span Numeric. Loess span parameter used when fitting the
#' variance-mean relationship
#' @param clip_max Character. After standardization values larger than clip.max
#' will be set to clip.max; default is 'auto' which sets this value to the
#' square root of the number of cells
#' @param num_bin Numeric. Total number of bins to use in the scaled analysis
#' (default is 20)
#' @param binning_method Character. “equal_width”: each bin is of equal width
#' along the x-axis (default). Options: “equal_frequency”:
#' @param scale_max Numeric. Max value to return for scaled data. The default
#' is 10.
#' @param pca_key Character. Key name to save the pca result in
#' @param integration_key Character. Key name to save the integration result in
#' @param npcs Numeric. Total Number of PCs to compute and store (50 by default)
#' @param verbose Logical. Verbosity for all functions
#' @param random_seed parameter for random state initialisation
#'
#' @return integrated sce/seurat object
#'
#' @import Seurat
#'
#' @examples
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' DO.Integration(
#'     sce_object = sce_data,
#'     split_key = "orig.ident",
#'     HVG = TRUE,
#'     scale = TRUE,
#'     pca = TRUE,
#'     integration_method = "CCAIntegration"
#' )
#'
#' @export
DO.Integration <- function(sce_object,
    split_key = "orig.ident",
    HVG = FALSE,
    scale = FALSE,
    pca = FALSE,
    neighbors = TRUE,
    neighbors_dim = seq_len(50),
    clusters = TRUE,
    clusters_res = 0.3,
    clusters_algorithm = 4,
    umap = TRUE,
    umap_key = "UMAP",
    umap_dim = seq_len(50),
    integration_method = "CCAIntegration",
    selection_method = "vst",
    loess_span = 0.3,
    clip_max = "auto",
    num_bin = 20,
    binning_method = "equal_width",
    scale_max = 10,
    pca_key = "PCA",
    integration_key = "INTEGRATED.CCA",
    npcs = 50,
    verbose = FALSE,
    random_seed = 42) {
    # Cover SCE and Seurat
    if (methods::is(sce_object, "SingleCellExperiment")) {
        SCE <- TRUE
        sce_object <- .suppressDeprecationWarnings(as.Seurat(sce_object))
    } else {
        SCE <- FALSE
    }

    if (!any(grepl("^(counts|data|scale\\.data)\\..+", Layers(sce_object)))) {
        .logger(
            paste0(
                "Splitting object for integration with ",
                integration_method,
                " by ",
                split_key
            )
        )
        .suppressAllWarnings(
                sce_object[["RNA"]] <-
                    split(sce_object[["RNA"]],
                        f = sce_object@meta.data[[split_key]]
            )
        )
    }

    if (HVG) {
        .logger("Calculating highly variable genes")
        sce_object <- FindVariableFeatures(sce_object,
            selection.method = selection_method,
            loess.span = loess_span,
            clip.max = clip_max,
            num.bin = num_bin,
            binning.method = binning_method,
            verbose = verbose
        )
    }

    if (scale) {
        .logger("Scaling object")
        sce_object <- ScaleData(
            object = sce_object,
            scale.max = scale_max,
            verbose = verbose
        )
    }

    if (pca) {
        .logger(paste0("Running pca, saved in key: ", pca_key))
        sce_object <- RunPCA(sce_object,
            verbose = verbose,
            reduction.name = pca_key,
            npcs = npcs,
            seed.use = random_seed
        )
    }


    sce_object <- JoinLayers(sce_object)
    .suppressAllWarnings(
        sce_object[["RNA"]] <-
                split(sce_object[["RNA"]],
                    f = sce_object@meta.data[[split_key]]
            )
        )

    .logger(paste0("Running integration, saved in key: ", integration_key))
    # Integration of Seurat covers many different methods in one call
    sce_object <- .suppressDeprecationWarnings(
        IntegrateLayers(
            object = sce_object,
            method = integration_method,
            orig.reduction = pca_key,
            new.reduction = integration_key,
            verbose = verbose
        )
    )

    # After Integration we join the layers
    sce_object <- JoinLayers(sce_object)

    if (neighbors) {
        .logger("Running Nearest-neighbor graph construction")
        sce_object <- FindNeighbors(
            object = sce_object,
            reduction = integration_key,
            dims = neighbors_dim,
            verbose = verbose
        )
    }

    if (clusters) {
        algo_names <- c(
            "1" = "louvainOri",
            "2" = "louvainRef",
            "3" = "SLM",
            "4" = "leiden"
        )

        cl_name <- paste0(
            algo_names[as.character(clusters_algorithm)], clusters_res
        )

        .logger("Running cluster detection")
        sce_object <- FindClusters(
            object = sce_object,
            resolution = clusters_res,
            algorithm = clusters_algorithm,
            random.seed = 42,
            cluster.name = cl_name,
            verbose = verbose
        )
        # clean up
        if (cl_name != "seurat_clusters") {
            sce_object$seurat_clusters <- NULL
        }
    }

    if (umap) {
        .logger("Creating UMAP")
        options(Seurat.warn.umap.uwot = FALSE)
        sce_object <- RunUMAP(
            object = sce_object,
            reduction = integration_key,
            reduction.name = umap_key,
            dims = umap_dim,
            verbose = verbose
        )
    }

    if (SCE) {
        sce_object <- as.SingleCellExperiment(sce_object)
    }

    return(sce_object)
}
