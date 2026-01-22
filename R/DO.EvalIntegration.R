#' @title Do batch correction metrics for integration
#' @author Mariano Ruz Jurado
#' @description This function calculates different metrics to evaluate the
#' integration of scRNA expression matrices in a new dimension. Its a wrapper
#' function around scib batch correction metrics
#' @param sce_object Seurat or SCE object.
#' @param assay character, Name of the assay the integration is saved in
#' @param label_key character, Annotation column
#' @param batch_key character, Sample column
#' @param type_ character, default: "embed"
#' @param pcr_covariate character, covariate column for pcr
#' @param pcr_n_comps integer, number of components for pcr
#' @param scale boolean, default: TRUE
#' @param verbose boolean, defult: FALSE
#' @param n_cores integer, Number of cores used for calculations
#' @param integration character, Name of the integration to evaluate
#' @param kBET boolean, if kBET should be run
#' @param cells.use vector, named cells to use for kBET subsetting
#' @param subsample float, for starified subsampling,
#' @param min_per_batch integer, minimum number of cells per batch
#' @param all_scores_silhouette boolean,
#' define if all scores of silhouette return
#' @param ... Additionally arguments for kBET
#' @return DataFrame with score for the given integration
#'
#' @import SeuratObject
#'
#' @examples
#' \dontrun{
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#' DO.EvalIntegration(
#'     sce_object = sce_data,
#'     label_key = "annotation",
#'     batch_key = "orig.ident",
#'     type_ = "embed",
#'     pcr_covariate = "orig.ident",
#'     pcr_n_comps = 30,
#'     scale = TRUE,
#'     verbose = FALSE,
#'     n_cores = 10,
#'     assay = "RNA",
#'     integration = "INTEGRATED.CCA",
#'     kBET = TRUE,
#'     cells.use = NULL,
#'     subsample = NULL,
#'     min_per_batch = NULL,
#'     all_scores_silhouette = FALSE
#' )
#' }
#' @export
DO.EvalIntegration <- function(
    sce_object,
    label_key = "annotation",
    batch_key = "orig.ident",
    type_ = "embed",
    pcr_covariate = "orig.ident",
    pcr_n_comps = 30,
    scale = TRUE,
    verbose = FALSE,
    n_cores = 10,
    assay = "RNA",
    integration = "INTEGRATED.CCA",
    kBET = TRUE,
    cells.use = NULL,
    subsample = NULL,
    min_per_batch = NULL,
    all_scores_silhouette = FALSE,
    ...
) {
    # #Testing
    # sce_object=sce_data
    # label_key = "annotation"
    # batch_key = "orig.ident"
    # type_ = "embed"
    # pcr_covariate = "orig.ident"
    # pcr_n_comps = 30
    # scale = TRUE
    # verbose = FALSE
    # n_cores = 10
    # assay = "RNA"
    # integration = "INTEGRATED.CCA"
    # kBET = TRUE
    # cells.use = NULL
    # subsample = NULL
    # min_per_batch = NULL
    # all_scores_silhouette = FALSE


    .logger(paste("Evaluating integration:", integration))

    # support for Seurat objects
    if (methods::is(sce_object, "Seurat")) {
        DefaultAssay(sce_object) <- assay
        sce_object <- .suppressDeprecationWarnings(
            Seurat::as.SingleCellExperiment(sce_object,
                assay = assay
            )
        )
    }

    if (kBET) {
        res_kbet <- run_kbet(
            sce_object = sce_object,
            embedding = integration,
            batch = batch_key,
            subsample = subsample,
            cells.use = cells.use,
            ...
        )

        # Take the mean from kBET scoring
        mean_kBET <- mean(res_kbet$stats$kBET.observed)
    }


    # source PATH to python script in install folder
    path_py <- system.file("python", "eval_metrics.py", package = "DOtools")

    # arguments for scib call in python
    args <- list(
        sce_object = sce_object,
        label_key = label_key,
        batch_key = batch_key,
        type_ = type_,
        embed = integration,
        use_rep = integration,
        covariate = pcr_covariate,
        n_comps = pcr_n_comps,
        scale = scale,
        verbose = verbose,
        n_cores = n_cores,
        return_all_silhouette = all_scores_silhouette
    )


    # basilisk implementation
    # this part runs from scib package: graph connectivity,
    # pcr_comparison and silhouette_batch
    results <- basilisk::basiliskRun(env = DOtoolsEnv, fun = function(args) {
        AnnData_counts <- zellkonverter::SCE2AnnData(args$sce_object,
            X_name = "logcounts"
        )

        reticulate::source_python(path_py)

        # python code call
        eval_integration(
            adata = AnnData_counts,
            label_key = args$label,
            batch_key = args$batch_key,
            embed = args$embed,
            use_rep = args$use_rep,
            covariate = args$covariate,
            type_ = args$type_,
            n_comps = as.integer(args$n_comps),
            n_cores = as.integer(args$n_cores),
            scale = args$scale,
            verbose = args$verbose,
            return_all_silhouette = args$all_scores_silhouette
        )
    }, args = args)

    if (condition) {}
    kBET_df <- data.frame(
        setNames(list(mean_kBET), integration),
        row.names = "kBET"
    )
    results <- rbind(results, kBET_df)

    return(results)
}


#' @title kBET function
#' @author Mariano Ruz Jurado
#' @description kBET quantifies batch mixing in single-cell data by testing
#' whether local neighborhood batch composition deviates from the global
#' distribution, with lower rejection indicating better integration.
#' @param sce_object Seurat or SCE object.
#' @param embedding Name of the embedding to test
#' @param batch Name of the sample column in meta.data
#' @param cells.use specify a specific set of cells as vector
#' @param subsample set a fraction for stratified subsetting
#' @param min_per_batch should always be higher than or equal k0
#' @return list with summary stats of kBET
#'
#'
#' @keywords internal
run_kbet <- function(
    sce_object,
    embedding,
    batch,
    cells.use = NULL,
    subsample = NULL,
    min_per_batch = NULL,
    ...
) {
    # kBET needs to be installed
    if (!requireNamespace("kBET", quietly = TRUE)) {
        stop(
            "kBET is not installed.\n",
            "Install it with:\n",
            "  remotes::install_github('theislab/kBET')",
            call. = FALSE
        )
    }

    # extract embedding from SCE Objects and Seurat
    if (inherits(sce_object, "Seurat")) {
        if (!embedding %in% names(sce_object@reductions)) {
            stop("Embedding '", embedding, "' not found in Seurat object.")
        }
        X <- Embeddings(sce_object, reduction = embedding)

        meta <- sce_object@meta.data
    } else if (inherits(sce_object, "SingleCellExperiment")) {
        if (!embedding %in% SingleCellExperiment::reducedDimNames(sce_object)) {
            stop(
                "Embedding '",
                embedding, "' not found in SingleCellExperiment."
                )
        }
        X <- SingleCellExperiment::reducedDim(sce_object, embedding)

        meta <- as.data.frame(SummarizedExperiment::colData(sce_object))
    } else {
        stop("Unsupported object type. Use Seurat or SingleCellExperiment.")
    }


    # check if batch available in meta data
    if (!batch %in% colnames(meta)) {
        stop("Batch column '", batch, "' not found in metadata.")
    }

    # subsetting by stratified fraction
    if (!is.null(subsample)) {
        batches <- unique(meta[[batch]])

        cells.sub <- unlist(lapply(batches, function(b) {
            cells <- rownames(meta)[meta[[batch]] == b]
            n <- floor(length(cells) * subsample)

            if (!is.null(min_per_batch)) {
                n <- max(n, min_per_batch)
            }

            if (n > length(cells)) {
                stop("Not enough cells in batch ",
                    b,
                    " for requested subsample.")
            }

            sample(cells, n, replace = FALSE)
        }))

        cells.use <- cells.sub
    }

    # subsetting if specified
    if (!is.null(cells.use)) {
        X <- X[cells.use, , drop = FALSE]
        meta <- meta[cells.use, , drop = FALSE]
    }

    batch_vec <- as.factor(meta[[batch]])

    # kBET call with additional parameters
    res <- kBET::kBET(
        df = as.matrix(X),
        batch = batch_vec,
        plot = FALSE,
        ...
    )
    return(res)
}
