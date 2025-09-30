#' @title DO.CellBender
#' @description This function wraps a system call to a bash script for running
#' CellBender on CellRanger outputs. It ensures required inputs are available
#' and optionally installs CellBender in a conda env.
#'
#' @param cellranger_path Path to folder with CellRanger outputs.
#' @param output_path Output directory for CellBender results.
#' @param samplenames Optional vector of sample names. If NULL, will
#' autodetect folders in `cellranger_path`.
#' @param cuda Logical, whether to use GPU (CUDA).
#' @param cpu_threads Number of CPU threads to use.
#' @param epochs Number of training epochs.
#' @param lr Learning rate.
#' @param estimator_multiple_cpu Use estimator with multiple CPU threads.
#' @param log Whether to enable logging.
#' @param conda_path Optional path to the conda environment.
#' @param BarcodeRanking Optional Calculation of estimated cells in samples
#' through DropletUtils implementation
#' @param bash_script Path to the bash script that runs CellBender.
#'
#'
#' @import DropletUtils
#' @return None
#'
#' @examples
#' \dontrun{
#' # Define paths
#' cellranger_path <- "/mnt/data/cellranger_outputs"
#' output_path <- "/mnt/data/cellbender_outputs"
#'
#' # Optional: specify sample names if automatic detection is not desired
#' samplenames <- c("Sample_1", "Sample_2")
#'
#' # Run CellBender (uses GPU by default)
#' DO.CellBender(
#'     cellranger_path = cellranger_path,
#'     output_path = output_path,
#'     samplenames = samplenames,
#'     cuda = TRUE,
#'     cpu_threads = 8,
#'     epochs = 100,
#'     lr = 0.00001,
#'     estimator_multiple_cpu = FALSE,
#'     log = TRUE
#' )
#' }
#'
#' @export
DO.CellBender <- function(cellranger_path,
    output_path,
    samplenames = NULL,
    cuda = TRUE,
    cpu_threads = 15,
    epochs = 150,
    lr = 0.00001,
    estimator_multiple_cpu = FALSE,
    log = TRUE,
    conda_path = NULL,
    BarcodeRanking = TRUE,
    bash_script = system.file("bash",
        "_run_CellBender.sh",
        package = "DOtools"
    )) {
    # Check input paths
    stopifnot(file.exists(cellranger_path))
    stopifnot(file.exists(output_path))

    # Warnings and logs
    if (estimator_multiple_cpu) {
        .logger(paste0(
            "Warning: estimator_multiple_cpu is TRUE. Not recommended for ",
            "large datasets (above 20 to 30k cells)."
        ))
    }
    if (epochs > 150) {
        .logger(paste0(
            "Warning: Training for more than 150 epochs may lead to ",
            "overfitting."
        ))
    }
    if (!cuda) {
        .logger(paste0(
            "Warning: Running without CUDA (GPU) may significantly increase ",
            "run time."
        ))
    }

    # Handle conda environment
    if (is.null(conda_path)) {
        conda_path <-
            file.path(
                path.expand("~"),
                ".cache/R/basilisk/1.20.0/DOtools/0.99.0/CellBen_env"
            )
    } else {
        conda_path <- path.expand(conda_path)
    }

    if (!dir.exists(conda_path)) {
        .logger("Creating conda environment for CellBender...")

        conda_args1 <-
            c(
                "conda",
                "create",
                "-y",
                "-p",
                file.path(
                    path.expand("~"),
                    ".cache/R/basilisk/1.20.0/DOtools/0.99.0/CellBen_env"
                ),
                "python=3.7"
            )

        conda_args2 <-
            c(
                "conda",
                "run",
                "-p",
                file.path(
                    path.expand("~"),
                    ".cache/R/basilisk/1.20.0/DOtools/0.99.0/CellBen_env"
                ),
                "pip",
                "install",
                "cellbender",
                "lxml_html_clean"
            )

        tryCatch(
            {
                system2(conda_args1[1],
                    args = conda_args1[-1],
                    stdout = FALSE,
                    stderr = FALSE
                )
                system2(conda_args2[1],
                    args = conda_args2[-1],
                    stdout = FALSE,
                    stderr = FALSE
                )
            },
            error = function(e) {
                stop(
                    "Failed to create conda environment. Provide a valid ",
                    "environment path or fix installation."
                )
            }
        )
    } else {
        .logger(sprintf("Using existing conda environment at: %s", conda_path))
    }

    # Detect samples
    if (is.null(samplenames)) {
        samples <- list.dirs(cellranger_path,
            full.names = FALSE,
            recursive = FALSE
        )
    } else {
        samples <- samplenames
    }

    .logger(sprintf("Running CellBender for %d sample(s)", length(samples)))

    # Loop through each sample
    for (sample in samples) {
        # Load HDF5 with reticulate or scanpy for expected_cells, total_droplets

        # Run one by one but sequentially
        # Estimate the number of cells to be used as upper and lower bound
        h5_file <- file.path(
            cellranger_path,
            sample,
            "outs",
            "raw_feature_bc_matrix.h5"
        )

        tdata <- DropletUtils::read10xCounts(h5_file)

        if (BarcodeRanking == TRUE) {
            result_barcoderanks <- .DO.BarcodeRanks(tdata)

            # Build command
            cmd <- c(
                "conda", "run", "-p", conda_path,
                "bash", bash_script,
                "-i", sample, "-o", output_path,
                "--cellRanger-output", cellranger_path,
                "--cpu-threads", cpu_threads,
                "--epochs", epochs,
                "--lr", lr,
                "--expected-cells", result_barcoderanks[1],
                "--total-droplets", result_barcoderanks[2]
            )
        } else {
            # Build command
            cmd <- c(
                "conda", "run", "-p", conda_path,
                "bash", bash_script,
                "-i", sample, "-o", output_path,
                "--cellRanger-output", cellranger_path,
                "--cpu-threads", cpu_threads,
                "--epochs", epochs,
                "--lr", lr
            )
        }

        if (cuda) cmd <- c(cmd, "--cuda")
        if (log) cmd <- c(cmd, "--log")
        if (estimator_multiple_cpu) cmd <- c(cmd, "--estimator_multiple_cpu")

        # Run command
        .logger(sprintf("Running CellBender for sample: %s", sample))
        tryCatch(
            {
                system2(cmd[1], args = cmd[-1], stdout = TRUE, stderr = TRUE)
            },
            error = function(e) {
                .logger(
                    sprintf(
                        "Error running CellBender for sample %s: %s",
                        sample,
                        e$message
                    )
                )
            }
        )
    }

    .logger("Finished running CellBender.")
    invisible(NULL)
}

#' @title DO.BarcodeRanks
#' @author Mariano Ruz Jurado
#' @description Given a raw count matrix (e.g. from a CellRanger HDF5 file),
#' estimate the number of expected cells and droplets using the knee and
#' inflection points from barcodeRanks.
#'
#' @param SCE_obj A Single cell experiment object.
#'
#' @return A named numeric vector: `c(xpc_cells = ..., total_cells = ...)
#' @rdname dot-DO.BarcodeRanks
#' @keywords internal
.DO.BarcodeRanks <- function(SCE_obj) {
    if (!inherits(SCE_obj, c("SingleCellExperiment"))) {
        stop(
            "Input must be a sparse matrix of class 'dgCMatrix' ",
            "or 'DelayedMatrix'."
        )
    }

    result <- DropletUtils::barcodeRanks(SCE_obj)
    metadata <- S4Vectors::metadata(result)

    knee <- metadata$knee
    inflection <- metadata$inflection

    colsum <- colSums(SCE_obj@assays@data$counts)

    xpc_cells <- length(colsum[colsum > knee])
    total_cells <- length(colsum[colsum > inflection])

    return(c(xpc_cells = xpc_cells, total_cells = total_cells))
}
