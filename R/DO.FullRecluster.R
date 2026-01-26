# Recluster function using FindSubCluster function from Seurat
#' @author Mariano Ruz Jurado
#' @title DO.FullRecluster
#' @description Performs iterative reclustering on each major cluster found by
#' FindClusters in a Seurat or SCE object. It refines the clusters using the
#' FindSubCluster function for better resolution and fine-tuned annotation. The
#' new clustering results are stored in a metadata column called
#' `annotation_recluster`. Suitable for improving cluster precision and
#' granularity after initial clustering.
#' @param sce_object The seurat or SCE object
#' @param over_clustering Column in metadata in object with clustering
#' assignments for cells, default seurat_clusters
#' @param res Resolution for the new clusters, default 0.5
#' @param algorithm Set one of the available algorithms found in FindSubCLuster
#' function, default = 4: leiden
#' @param graph.name A builded neirest neighbor graph
#' @param random_seed parameter for random state initialisation
#' @return a Seurat or SCE Object with new clustering named annotation_recluster
#'
#' @import Seurat
#' @import progress
#'
#' @examples
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' sce_data <- DO.FullRecluster(
#'     sce_object = sce_data
#' )
#'
#' @export
DO.FullRecluster <- function(sce_object,
    over_clustering = "seurat_clusters",
    res = 0.5,
    algorithm = 4,
    graph.name = "RNA_snn",
    random_seed = 42) {

    set.seed(random_seed)
    # support for single cell experiment objects
    if (methods::is(sce_object, "SingleCellExperiment")) {
        class_obj <- "SingleCellExperiment"
        sce_object <- .suppressAllWarnings(as.Seurat(sce_object))
        sce_object <- FindNeighbors(sce_object, reduction = "PCA")
    } else {
        class_obj <- "Seurat"
    }

    # TODO add an argument for a new shared nearest neighbor call so you can use
    # this function with other integration methods in the object otherwise it
    # will use the latest calculated shared nearest neighbors and you need to do
    # it outside of the function
    if (is.null(sce_object@meta.data[[over_clustering]])) {
        stop(
            "No clusters defined, please run a cluster algorithm before ",
            "reclustering and save it to the object"
        )
    }
    Idents(sce_object) <- over_clustering

    sce_object$annotation_recluster <-
        as.vector(sce_object@meta.data[[over_clustering]])

    pb <- progress::progress_bar$new(
        total = length(unique(sce_object@meta.data[[over_clustering]])),
        format = "  Reclustering [:bar] :percent eta: :eta"
    )

    # newline prints when the function exits (to clean up console)
    on.exit(message(""))

    for (cluster in unique(sce_object@meta.data[[over_clustering]])) {
        pb$tick()
        sce_object <- .suppressAllWarnings(FindSubCluster(sce_object,
            cluster = as.character(cluster),
            graph.name = graph.name,
            algorithm = algorithm,
            resolution = res
        ))

        cluster_cells <- rownames(sce_object@meta.data)[
            sce_object@meta.data[[over_clustering]] == cluster
        ]
        sce_object$annotation_recluster[cluster_cells] <-
            sce_object$sub.cluster[cluster_cells]
    }
    sce_object$sub.cluster <- NULL

    if (class_obj == "SingleCellExperiment") {
        sce_object <- as.SingleCellExperiment(sce_object)
    }

    return(sce_object)
}
