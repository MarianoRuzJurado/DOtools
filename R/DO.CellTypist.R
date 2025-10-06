#' @author Mariano Ruz Jurado
#' @title DO Celltypist
#' @description Runs the CellTypist model on a Seurat or SCE object to predict
#' cell type labels, storing the results as metadata. If the number of cells is
#' less than the specified threshold, it returns NAs for the labels. Optionally
#' updates the CellTypist models and returns the probability matrix. Useful for
#' annotating cell types in single-cell RNA sequencing datasets.
#' @param sce_object The seurat or sce object
#' @param modelName Specify the model you want to use for celltypist
#' @param minCellsToRun If the input seurat or SCE object has fewer than this
#' many cells, NAs will be added for all expected columns and celltypist will
#' not be run.
#' @param runCelltypistUpdate If true, --update-models will be run for
#' celltypist prior to scoring cells.
#' @param over_clustering Column in metadata in object with clustering
#' assignments for cells, default seurat_clusters
#' @param assay_normalized Assay with log1p normalized expressions
#' @param returnProb will additionally return the probability matrix, return
#' will give a list with the first element beeing the object and second prob
#' matrix
#' @param SeuV5 Specify if the Seurat object is made with Seuratv5
#'
#' @importFrom basilisk basiliskRun
#' @import dplyr
#' @import ggplot2
#'
#' @return a seurat or sce object
#'
#' @examples
#' sce_data <-
#'     readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))
#'
#'
#' sce_data <- DO.CellTypist(
#'     sce_object = sce_data,
#'     modelName = "Healthy_Adult_Heart.pkl",
#'     runCelltypistUpdate = TRUE,
#'     over_clustering = "seurat_clusters",
#'     minCellsToRun = 5,
#'     SeuV5 = TRUE
#' )
#'
#' @export
DO.CellTypist <- function(sce_object,
    modelName = "Healthy_Adult_Heart.pkl",
    minCellsToRun = 200,
    runCelltypistUpdate = TRUE,
    over_clustering = "seurat_clusters",
    assay_normalized = "RNA",
    returnProb = FALSE,
    SeuV5 = TRUE) {
    # Make sure R Zellkonverter package is installed
    zk <- system.file(package = "zellkonverter")
    ifelse(nzchar(zk), "", stop(
        "Install zellkonverter R package for object ",
        "tranformation!"
    ))

    # Make sure R reticulate package is installed
    rt <- system.file(package = "reticulate")
    ifelse(nzchar(rt), "", stop(
        "Install reticulate R package for Python usage ",
        "in R!"
    ))


    if (ncol(sce_object) < minCellsToRun) {
        warning(
            "Too few cells, will not run celltypist. NAs will be added instead"
            )
        expectedCols <- c("predicted_labels_celltypist")
        sce_object[[expectedCols]] <- NA
        return(sce_object)
    }

    .logger(paste0("Running celltypist using model: ", modelName))

    # Create temporary folder
    outDir <- tempfile(fileext = "")
    if (endsWith(outDir, "/")) {
        outDir <- gsub(outDir, pattern = "/$", replacement = "")
    }
    dir.create(outDir)
    .logger(paste0("Saving celltypist results to temporary folder: ", outDir))

    # Uppercasing gene names
    # zellkonverter h5ad

    # support for Seurat objects
    if (methods::is(sce_object, "Seurat")) {
        DefaultAssay(sce_object) <- assay_normalized
        if (SeuV5 == TRUE) {
            tmp.assay <- sce_object
            tmp.assay[["RNA"]] <- methods::as(tmp.assay[["RNA"]],
                Class = "Assay")
            tmp.sce <- Seurat::as.SingleCellExperiment(tmp.assay,
                assay = assay_normalized
            )
            rownames(tmp.sce) <- toupper(rownames(tmp.sce))
        } else {
            tmp.sce <- Seurat::as.SingleCellExperiment(sce_object,
                assay = assay_normalized
            )
            rownames(tmp.sce) <- toupper(rownames(tmp.sce))
        }
    } else {
        tmp.sce <- sce_object
        rownames(tmp.sce) <- toupper(rownames(tmp.sce))
    }

    # Make Anndata object
    if (!"logcounts" %in% names(tmp.sce@assays)) {
        stop("logcounts not found in assays of object!")
    }

    zellkonverter::writeH5AD(tmp.sce,
        file = paste0(outDir, "/ad.h5ad"),
        X_name = "logcounts"
    )

    # basilisk implementation
    results <- basilisk::basiliskRun(env = DOtoolsEnv, fun = function(args) {
        # Ensure models present:
        # if (runCelltypistUpdate) {
        #     system2(
        #     reticulate::py_exe(),
        #     c("-m", "celltypist.command_line", "--update-models", "--quiet")
        #     )
        # }
        if (runCelltypistUpdate) {
            ct <- reticulate::import("celltypist")

            models_path <- ct$models$models_path
            model_file <- file.path(models_path, modelName)

            if (!file.exists(model_file)) {
                .logger(paste0("Downloading CellTypist model: ", modelName))
                # download only the specific model
                ct$models$download_models(
                    model = modelName, force_update = FALSE
                )
            } else {
                .logger(paste0("Model already present: ", model_file))
                }
            }
        .logger("Running Celltypist")

        # Run celltypist:
        args_cty <- c(
            "-m",
            "celltypist.command_line",
            "--indata",
            paste0(outDir, "/ad.h5ad"),
            "--model",
            modelName,
            "--outdir",
            outDir,
            "--majority-voting",
            "--over-clustering",
            over_clustering
        )
        system2(reticulate::py_exe(), args_cty)
    }, args = args)

    labelFile <- paste0(outDir, "/predicted_labels.csv")
    probFile <- paste0(outDir, "/probability_matrix.csv")
    labels <- utils::read.csv(
        labelFile,
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
    )

    # ad <- import(module = "anndata")
    # ad_obj <- ad$AnnData(X = labels)

    # ct <- import(module = "celltypist")

    probMatrix <- utils::read.csv(
        probFile,
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
    )
    sce_object$autoAnnot <- labels$majority_voting

    # Create dotplot with prob
    probMatrix$cluster <- as.character(labels$over_clustering)

    # calculate means
    probMatrix_mean <- probMatrix %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(across(where(is.numeric), mean))

    top_cluster <- probMatrix_mean %>%
        rowwise() %>%
        dplyr::mutate(
            prob = max(c_across(-cluster)),
            label = names(.)[which.max(c_across(-cluster))]
        ) %>%
        select(cluster, label, prob)
    top_cluster$pct.exp <- 100
    # since majority voting is set TRUE
    # TODO make this general if majority voting will
    # become a boolean argument in the future

    top_cluster$label <-
        factor(top_cluster$label,
            levels = unique(sort(top_cluster$label,
                decreasing = TRUE
            ))
        )
    top_cluster$cluster <-
        factor(top_cluster$cluster,
            levels = top_cluster$cluster[order(top_cluster$label,
                decreasing = TRUE
            )]
        )

    .logger("Creating probality plot")

    pmain <- ggplot(top_cluster, aes(x = cluster, y = label)) +
        geom_point(
            aes(size = pct.exp, fill = prob),
            shape = 21,
            color = "black",
            stroke = 0.3
        ) +
        scale_fill_gradient2(
            low = "royalblue3",
            mid = "white",
            high = "firebrick",
            midpoint = 0.5,
            limits = c(0, 1)
        ) +
        scale_size(
            range = c(2, 10),
            breaks = c(20, 40, 60, 80, 100),
            limits = c(0, 100)
        ) +
        theme_box() +
        theme(
            plot.margin = ggplot2::margin(
                t = 1,
                r = 1,
                b = 1,
                l = 1,
                unit = "cm"
            ),
            axis.text = ggplot2::element_text(color = "black"),
            legend.direction = "horizontal",
            axis.text.x = element_text(
                color = "black",
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = 14,
                family = "Helvetica"
            ),
            axis.text.y = element_text(
                color = "black",
                size = 14,
                family = "Helvetica"
            ),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.title = element_text(
                size = 14,
                color = "black",
                family = "Helvetica"
            ),
            plot.title = element_text(
                size = 14,
                hjust = 0.5,
                face = "bold",
                family = "Helvetica"
            ),
            plot.subtitle = element_text(
                size = 14,
                hjust = 0,
                family = "Helvetica"
            ),
            axis.line = element_line(color = "black"),
            strip.text.x = element_text(
                size = 14,
                color = "black",
                family = "Helvetica",
                face = "bold"
            ),
            legend.text = element_text(
                size = 10,
                color = "black",
                family = "Helvetica"
            ),
            legend.title = element_text(
                size = 10,
                color = "black",
                family = "Helvetica",
                hjust = 0,
                face = "bold"
            ),
            legend.position = "right",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
        )

    guides.layer <- ggplot2::guides(
        fill = ggplot2::guide_colorbar(
            title = "Mean propability",
            title.position = "top",
            title.hjust = 0.5,
            barwidth = unit(3.8, "cm"),
            # changes the width of the color legend
            barheight = unit(0.5, "cm"),
            frame.colour = "black",
            frame.linewidth = 0.3,
            ticks.colour = "black",
            order = 2
        ),
        size = ggplot2::guide_legend(
            title = "Fraction of cells \n in group (%)",
            title.position = "top",
            title.hjust = 0.5,
            label.position = "bottom",
            override.aes = list(color = "grey40", fill = "grey50"),
            keywidth = ggplot2::unit(0.5, "cm"),
            # changes the width of the percentage dots in legend
            order = 1
        )
    )

    pmain <- pmain + guides.layer
    pmain

    if (returnProb == TRUE) {
        returnProb <- list(sce_object, probMatrix)
        names(returnProb) <- c("SingleCellObject", "probMatrix")
        return(returnProb)
    } else {
        return(sce_object)
    }
}
