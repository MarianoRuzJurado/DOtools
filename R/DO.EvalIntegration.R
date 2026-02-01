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
        res_kbet <- .run_kbet(
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
.run_kbet <- function(
    sce_object,
    embedding,
    batch,
    cells.use = NULL,
    subsample = NULL,
    min_per_batch = NULL,
    ...
) {
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
                stop(
                    "Not enough cells in batch ",
                    b,
                    " for requested subsample."
                )
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
    res <- .kBET_fct(
        df = as.matrix(X),
        batch = batch_vec,
        plot = FALSE,
        ...
    )
    return(res)
}

#' kBET - k-nearest neighbour batch effect test
#' @author Mariano Ruz Jurado (edited from: Maren Buettner)
#'
#' @description \code{kBET} runs a chi square test to evaluate
#' the probability of a batch effect.
#'
#' @param df dataset (rows: cells, columns: features)
#' @param batch batch id for each cell or a data frame with
#' both condition and replicates
#' @param k0 number of nearest neighbours to test on (neighbourhood size)
#' @param knn an n x k matrix of nearest neighbours for each cell (optional)
#' @param testSize number of data points to test,
#' (10 percent sample size default, but at least 25)
#' @param do.pca perform a pca prior to knn search? (defaults to TRUE)
#' @param dim.pca if do.pca=TRUE, choose the number of dimensions
#' to consider (defaults to 50)
#' @param heuristic compute an optimal neighbourhood size k
#' (defaults to TRUE)
#' @param n_repeat to create a statistics on batch estimates,
#' evaluate 'n_repeat' subsets
#' @param alpha significance level
#' @param adapt In some cases, a number of cells do not contribute
#' to any neighbourhood
#' and this may cause an imbalance in observed and expected batch
#' label frequencies.
#' Frequencies will be adapted if adapt=TRUE (default).
#' @param addTest perform an LRT-approximation to the multinomial
#' test AND a multinomial exact test (if appropriate)
#' @param plot if stats > 10, then a boxplot of the resulting
#' rejection rates is created
#' @param verbose displays stages of current computation (defaults to FALSE)
#' @return list object
#'    \enumerate{
#'    \item \code{summary} - a rejection rate for the data,
#'         an expected rejection rate for random
#'         labeling and the significance for the observed result
#'    \item \code{results} - detailed list for each tested cells;
#'         p-values for expected and observed label distribution
#'    \item \code{average.pval} - significance level over the averaged
#'         batch label distribution in all neighbourhoods
#'    \item \code{stats} - extended test summary for every sample
#'    \item \code{params} - list of input parameters and adapted parameters,
#'    respectively
#'    \item \code{outsider} - only shown if \code{adapt=TRUE}. List of samples
#'         without mutual nearest neighbour: \itemize{
#'     \item \code{index} - index of each outsider sample)
#'     \item \code{categories} - tabularised labels of outsiders
#'     \item \code{p.val} - Significance level of outsider batch label
#'     distribution vs expected frequencies.
#'     If the significance level is lower than \code{alpha},
#'     expected frequencies will be adapted}
#'    }
#' @return If the optimal neighborhood size (k0) is smaller than 10,
#' NA is returned.
#' @import ggplot2
#' @importFrom stats quantile pchisq
#' @importFrom utils data
#' @importFrom FNN get.knn
#'
#' @keywords internal
.kBET_fct <- function(
    df,
    batch,
    k0 = NULL,
    knn = NULL,
    testSize = NULL,
    do.pca = TRUE,
    dim.pca = 50,
    heuristic = TRUE,
    n_repeat = 100,
    alpha = 0.05,
    addTest = FALSE,
    verbose = FALSE,
    plot = TRUE,
    adapt = TRUE
) {
    # preliminaries:
    dof <- length(unique(batch)) - 1 # degrees of freedom
    if (is.factor(batch)) {
        batch <- droplevels(batch)
    }

    frequencies <- table(batch) / length(batch)
    # get 3 different permutations of the batch label
    batch.shuff <- replicate(3, batch[sample.int(length(batch))])


    class.frequency <- data.frame(
        class = names(frequencies),
        freq = as.numeric(frequencies)
    )
    dataset <- df
    dim.dataset <- dim(dataset)
    # check the feasibility of data input
    if (dim.dataset[1] != length(batch) && dim.dataset[2] != length(batch)) {
        stop("Input matrix and batch information do not match. Execution halted.")
    }

    if (dim.dataset[2] == length(batch) && dim.dataset[1] != length(batch)) {
        if (verbose) {
            sprintf("kBET needs samples as rows. Transposing matrix...")
        }
        dataset <- t(dataset)
        dim.dataset <- dim(dataset)
    }
    # check if the dataset is too small per se
    if (dim.dataset[1] <= 10) {
        if (verbose) {
            sprintf("Your dataset has less than 10 samples. Abort and return NA.")
        }
        return(NA)
    }

    stopifnot(is.numeric(n_repeat), n_repeat > 0)


    ####
    chi_batch_test <- function(knn.set, class.freq, batch, df) {
        # knn.set: indices of nearest neighbours
        # empirical frequencies in nn-environment (sample 1)
        if (all(is.na(knn.set))) { # if all values of a neighbourhood are NA
            return(NA)
        } else {
            freq.env <- table(batch[knn.set[!is.na(knn.set)]])
            full.classes <- rep(0, length(class.freq$class))
            full.classes[class.freq$class %in% names(freq.env)] <- freq.env
            exp.freqs <- class.freq$freq * length(knn.set)
            # compute chi-square test statistics
            chi.sq.value <- sum((full.classes - exp.freqs)^2 / exp.freqs)
            result <- 1 - pchisq(chi.sq.value, df) # p-value for the result
            if (is.na(result)) { # I actually would like to now when 'NA' arises.
                return(NA)
            } else {
                result
            }
        }
    }
    ####

    do_heuristic <- FALSE
    if (is.null(k0) || k0 >= dim.dataset[1]) {
        do_heuristic <- heuristic
        if (!heuristic) {
            # default environment size: quarter the size of the largest batch
            k0 <- floor(mean(class.frequency$freq) * dim.dataset[1] / 4)
        } else {
            # default environment size: three quarter the size of the largest batch
            k0 <- floor(mean(class.frequency$freq) * dim.dataset[1] * 0.75)
            if (k0 < 10) {
                if (verbose) {
                    warning(
                        "Your dataset has too few samples to run a heuristic.\n",
                        "Return NA.\n",
                        "Please assign k0 and set heuristic=FALSE."
                    )
                }

                return(NA)
            }
        }
        if (verbose) {
            sprintf(paste0("Initial neighbourhood size is set to ", k0, ".\n"))
        }
    }

    # if k0 was set by the user and is too small & we do not operate on a
    # knn graph, abort
    # the reason is that if we want to test kBET on knn graph data integration
    # methods we usually face small numbers of nearest neighbours.
    if (k0 < 10 & is.null(knn)) {
        if (verbose) {
            warning(
                "Your dataset has too few samples to run a heuristic.\n",
                "Return NA.\n",
                "Please assign k0 and set heuristic=FALSE."
            )
        }
        return(NA)
    }
    # find KNNs
    if (is.null(knn)) {
        if (!do.pca) {
            if (verbose) {
                sprintf("finding knns...")
                tic <- proc.time()
            }
            # use the nearest neighbour index directly for further use in the package
            knn <- get.knn(dataset, k = k0, algorithm = "cover_tree")$nn.index
        } else {
            dim.comp <- min(dim.pca, dim.dataset[2])
            if (verbose) {
                sprintf("reducing dimensions with svd first...")
            }
            data.pca <- svd(x = dataset, nu = dim.comp, nv = 0)
            if (verbose) {
                sprintf("finding knns...")
                tic <- proc.time()
            }
            knn <- get.knn(data.pca$u, k = k0, algorithm = "cover_tree")
        }
        if (verbose) {
            .logger("Done")
        }
    }

    # backward compatibility for knn-graph
    if (is(knn, "list")) {
        knn <- knn$nn.index
        if (verbose) {
            sprintf("KNN input is a list, extracting nearest neighbour index.")
        }
    }

    # set number of tests
    if (is.null(testSize) || (floor(testSize) < 1 || dim.dataset[1] < testSize)) {
        test.frac <- 0.1
        testSize <- ceiling(dim.dataset[1] * test.frac)
        if (testSize < 25 && dim.dataset[1] > 25) {
            testSize <- 25
        }
        if (verbose) {
            sprintf("Number of kBET tests is set to %s", testSize)
        }
    }
    # decide to adapt general frequencies
    if (adapt) {
        # idx.run <- sample.int(dim.dataset[1], size = min(2*testSize, dim.dataset[1]))
        outsider <- which(!(seq_len(dim.dataset[1]) %in% knn[, seq_len(k0 - 1)]))
        is.imbalanced <- FALSE # initialisation
        p.out <- 1
        # avoid unwanted things happen if length(outsider) == 0
        if (length(outsider) > 0) {
            p.out <- chi_batch_test(outsider, class.frequency, batch, dof)
            if (!is.na(p.out)) {
                is.imbalanced <- p.out < alpha
                if (is.imbalanced) {
                    new.frequencies <- table(batch[-outsider]) / length(batch[-outsider])
                    new.class.frequency <- data.frame(
                        class = names(new.frequencies),
                        freq = as.numeric(new.frequencies)
                    )
                    if (verbose) {
                        outs_percent <-
                            round(length(outsider) / length(batch) * 100, 3)
                        message(
                            sprintf(
                                "There are %s cells (%s%%) that do not appear in any neighbourhood.\n%s\n%s",
                                length(outsider),
                                outs_percent,
                                "The expected frequencies for each category have been adapted.",
                                "Cell indexes are saved to result list."
                            )
                        )
                    }
                } else {
                    if (verbose) {
                        sprintf("No outsiders found.")
                    }
                }
            } else {
                if (verbose) {
                    sprintf("No outsiders found.")
                }
            }
        }
    }

    residual_score_batch <- function(knn.set, class.freq, batch) {
        # knn.set: indices of nearest neighbours
        # empirical frequencies in nn-environment (sample 1)
        # ignore NA entries (which may arise from subsampling a knn-graph)
        if (all(is.na(knn.set))) { # if all values of a neighbourhood are NA
            return(NA)
        } else {
            freq.env <- table(batch[knn.set[!is.na(knn.set)]]) / length(!is.na(knn.set))
            full.classes <- rep(0, length(class.freq$class))
            full.classes[class.freq$class %in% names(freq.env)] <- freq.env
            exp.freqs <- class.freq$freq
            # compute chi-square test statistics
            sum((full.classes - exp.freqs)^2 / exp.freqs)
        }
    }

    ptnorm <- function(x, mu, sd, a = 0, b = 1, alpha = 0.05, verbose = FALSE) {
        # this is the cumulative density of the truncated normal distribution
        # x ~ N(mu, sd^2), but we condition on a <= x <= b
        if (!is.na(x)) {
            if (a > b) {
                warning("Lower and upper bound are interchanged.")
                tmp <- a
                a <- b
                b <- tmp
            }

            if (sd <= 0 || is.na(sd)) {
                if (verbose) {
                    warning("Standard deviation must be positive.")
                }
                if (alpha <= 0) {
                    stop("False positive rate alpha must be positive.")
                }
                sd <- alpha
            }
            if (x < a || x > b) {
                warning("x out of bounds.")
                cdf <- as.numeric(x > a)
            } else {
                alp <- pnorm((a - mu) / sd)
                bet <- pnorm((b - mu) / sd)
                zet <- pnorm((x - mu) / sd)
                cdf <- (zet - alp) / (bet - alp)
            }
            cdf
        } else {
            return(NA)
        }
    }

    if (do_heuristic) {
        # btw, when we bisect here that returns some interval with
        # the optimal neihbourhood size
        if (verbose) {
            sprintf("Determining optimal neighbourhood size ...")
        }
        ###
        bisect <- function(foo, bounds, known = NULL, ..., tolx = 5, toly = 0.01) {
            if (is.null(known)) {
                evalFoo <- vapply(bounds, foo, FUN.VALUE = numeric(1), ...)

                if (diff(evalFoo) < -toly && diff(bounds) > tolx) {
                    bounds[2] <- round(sum(bounds) / 2, 0)
                    known <- c(evalFoo[1], 0)

                    bisect(foo, bounds, known, ..., tolx = tolx, toly = toly)
                } else if (diff(evalFoo) > toly && diff(bounds) > tolx) {
                    bounds[1] <- round(sum(bounds) / 2, 0)
                    known <- c(0, evalFoo[2])

                    bisect(foo, bounds, known, ..., tolx = tolx, toly = toly)
                } else if (dist(evalFoo) < toly) {
                    mid <- round(sum(bounds) / 2, 0)

                    center <- vapply(mid, foo, FUN.VALUE = numeric(1), ...)

                    dist_center <- vapply(
                        evalFoo,
                        function(x, y) dist(c(x, y)),
                        FUN.VALUE = numeric(1),
                        y = center
                    )

                    if (max(abs(dist_center)) < toly || dist(bounds) < tolx) {
                        bounds
                    } else if (dist_center[1] > toly) {
                        bounds[2] <- mid
                        known <- c(evalFoo[1], 0)

                        bisect(foo, bounds, known, ..., tolx = tolx, toly = toly)
                    } else if (dist_center[2] > toly) {
                        bounds[1] <- mid
                        known <- c(0, evalFoo[2])

                        bisect(foo, bounds, known, ..., tolx = tolx, toly = toly)
                    }
                }
            } else {
                new.eval <- which(known == 0)
                old.eval <- which(known != 0)

                evalFoo <- vapply(
                    bounds[new.eval],
                    foo,
                    FUN.VALUE = numeric(1),
                    ...
                )

                result <- numeric(length(known))
                result[new.eval] <- evalFoo
                result[old.eval] <- known[old.eval]

                if (diff(result) < -toly && diff(bounds) > tolx) {
                    bounds[2] <- round(sum(bounds) / 2, 0)
                    known <- c(evalFoo[1], 0)

                    bisect(foo, bounds, known, ..., tolx = tolx, toly = toly)
                } else if (diff(result) > toly && diff(bounds) > tolx) {
                    bounds[1] <- round(sum(bounds) / 2, 0)
                    known <- c(0, evalFoo[1])

                    bisect(foo, bounds, known, ..., tolx = tolx, toly = toly)
                } else if (dist(result) < toly || dist(bounds) < tolx) {
                    bounds
                }
            }
        }

        scan_nb <- function(x, df, batch, knn) {
            res <- .kBET_fct(
                df = df, batch = batch, k0 = x, knn = knn, testSize = NULL,
                heuristic = FALSE, n_repeat = 10, alpha = 0.05,
                addTest = FALSE, plot = FALSE, verbose = FALSE, adapt = FALSE
            )
            result <- res$summary
            result$kBET.observed[1]
        }

        ###
        opt.k <- bisect(scan_nb, bounds = c(10, k0), known = NULL, dataset, batch, knn)
        # result
        if (length(opt.k) > 1) {
            k0 <- opt.k[2]
            if (verbose) {
                sprintf("done. New size of neighbourhood is set to %s", k0)
            }
        } else {
            if (verbose) {
                message(
                    sprintf(
                        "done.\nHeuristic did not change the neighbourhood.\nIf results appear inconclusive, change k0=%s.",
                        k0
                    )
                )
            }
        }
    }

    # initialise result list
    rejection <- list()
    rejection$summary <- data.frame(
        kBET.expected = numeric(4),
        kBET.observed = numeric(4),
        kBET.signif = numeric(4)
    )

    rejection$results <- data.frame(
        tested = numeric(dim.dataset[1]),
        kBET.pvalue.test = rep(0, dim.dataset[1]),
        kBET.pvalue.null = rep(0, dim.dataset[1])
    )

    # get average residual score
    env <- as.vector(cbind(knn[, seq_len(k0 - 1)], seq_len(dim.dataset[1])))
    cf <- if (adapt && is.imbalanced) new.class.frequency else class.frequency
    rejection$average.pval <- 1 - pchisq(k0 * residual_score_batch(env, cf, batch), dof)


    # initialise intermediates
    kBET.expected <- numeric(n_repeat)
    kBET.observed <- numeric(n_repeat)
    kBET.signif <- numeric(n_repeat)


    if (addTest) {
        # initialize result list
        rejection$summary$lrt.expected <- numeric(4)
        rejection$summary$lrt.observed <- numeric(4)

        rejection$results$lrt.pvalue.test <- rep(0, dim.dataset[1])
        rejection$results$lrt.pvalue.null <- rep(0, dim.dataset[1])

        lrt.expected <- numeric(n_repeat)
        lrt.observed <- numeric(n_repeat)
        lrt.signif <- numeric(n_repeat)
        # decide to perform exact test or not
        if (choose(k0 + dof, dof) < 5e5 && k0 <= min(table(batch))) {
            exact.expected <- numeric(n_repeat)
            exact.observed <- numeric(n_repeat)
            exact.signif <- numeric(n_repeat)

            rejection$summary$exact.expected <- numeric(4)
            rejection$summary$exact.observed <- numeric(4)
            rejection$results$exact.pvalue.test <- rep(0, dim.dataset[1])
            rejection$results$exact.pvalue.null <- rep(0, dim.dataset[1])
        }

        for (i in seq_len(n_repeat)) {
            # choose a random sample from dataset (rows: samples, columns: features)
            idx.runs <- sample.int(dim.dataset[1], size = testSize)
            env <- cbind(knn[idx.runs, seq_len(k0 - 1)], idx.runs)
            # env.rand <- t(sapply(rep(dim.dataset[1],testSize),  sample.int, k0))

            # perform test
            if (adapt && is.imbalanced) {
                p.val.test <- apply(env, 1,
                    FUN = chi_batch_test,
                    new.class.frequency, batch, dof
                )
            } else {
                p.val.test <- apply(env, 1,
                    FUN = chi_batch_test,
                    class.frequency, batch, dof
                )
            }

            is.rejected <- p.val.test < alpha


            # p.val.test <- apply(env, 1, FUN = chi_batch_test, class.frequency,
            # batch,  dof)
            p.val.test.null <- apply(apply(
                batch.shuff, 2,
                function(x, freq, dof, envir) {
                    apply(envir, 1, FUN = chi_batch_test, freq, x, dof)
                },
                class.frequency, dof, env
            ), 1, mean, na.rm = TRUE)
            # p.val.test.null <- apply(env.rand, 1, FUN = chi_batch_test,
            # class.frequency, batch, dof)

            # summarise test results
            kBET.expected[i] <- sum(p.val.test.null < alpha,
                na.rm = TRUE
            ) / sum(!is.na(p.val.test.null))
            kBET.observed[i] <- sum(is.rejected, na.rm = TRUE) / sum(!is.na(p.val.test))

            # compute significance
            kBET.signif[i] <-
                1 - ptnorm(kBET.observed[i],
                    mu = kBET.expected[i],
                    sd = sqrt(kBET.expected[i] * (1 - kBET.expected[i]) / testSize),
                    alpha = alpha
                )

            # assign results to result table
            rejection$results$tested[idx.runs] <- 1
            rejection$results$kBET.pvalue.test[idx.runs] <- p.val.test
            rejection$results$kBET.pvalue.null[idx.runs] <- p.val.test.null

            ###
            lrt_approximation <- function(knn.set, class.freq, batch, df) {
                # knn.set: indices of nearest neighbours
                # empirical frequencies in nn-environment (sample 1)
                if (all(is.na(knn.set))) { # if all values of a neighbourhood are NA
                    return(NA)
                } else {
                    obs.env <- table(batch[knn.set[!is.na(knn.set)]]) # observed realisations of each category
                    freq.env <- obs.env / sum(obs.env) # observed 'probabilities'
                    full.classes <- rep(0, length(class.freq$class))
                    obs.classes <- class.freq$class %in% names(freq.env)
                    # for stability issues (to avoid the secret division by 0): introduce
                    # another alternative model where the observed probability
                    # is either the empirical frequency or 1/(sample size) at minimum
                    if (length(full.classes) > sum(obs.classes)) {
                        dummy.count <- length(full.classes) - sum(obs.classes)
                        full.classes[obs.classes] <- obs.env / (sum(obs.env) + dummy.count)
                        pmin <- 1 / (sum(obs.env) + dummy.count)
                        full.classes[!obs.classes] <- pmin
                    } else {
                        full.classes[obs.classes] <- freq.env
                    }
                    exp.freqs <- class.freq$freq # expected 'probabilities'
                    # compute likelihood ratio of null and alternative hypothesis,
                    # test statistics converges to chi-square distribution
                    full.obs <- rep(0, length(class.freq$class))
                    full.obs[obs.classes] <- obs.env

                    lrt.value <- -2 * sum(full.obs * log(exp.freqs / full.classes))

                    result <- 1 - pchisq(lrt.value, df) # p-value for the result
                    if (is.na(result)) { # I actually would like to now when 'NA' arises.
                        return(NA)
                    } else {
                        result
                    }
                }
            }
            ###

            # compute likelihood-ratio test (approximation for multinomial exact test)
            cf <- if (adapt && is.imbalanced) new.class.frequency else class.frequency
            p.val.test.lrt <- apply(env, 1, FUN = lrt_approximation, cf, batch, dof)
            p.val.test.lrt.null <- apply(apply(
                batch.shuff, 2,
                function(x, freq, dof, envir) {
                    apply(envir, 1, FUN = lrt_approximation, freq, x, dof)
                },
                class.frequency, dof, env
            ), 1, mean, na.rm = TRUE)

            lrt.expected[i] <-
                sum(p.val.test.lrt.null < alpha, na.rm = TRUE) / sum(!is.na(p.val.test.lrt.null))
            lrt.observed[i] <-
                sum(p.val.test.lrt < alpha, na.rm = TRUE) / sum(!is.na(p.val.test.lrt))

            lrt.signif[i] <-
                1 - ptnorm(lrt.observed[i],
                    mu = lrt.expected[i],
                    sd = sqrt(lrt.expected[i] * (1 - lrt.expected[i]) / testSize),
                    alpha = alpha
                )

            rejection$results$lrt.pvalue.test[idx.runs] <- p.val.test.lrt
            rejection$results$lrt.pvalue.null[idx.runs] <- p.val.test.lrt.null


            # exact result test consume a fairly high amount of computation time,
            # as the exact test computes the probability of ALL
            # possible configurations (under the assumption that all batches
            # are large enough to 'imitate' sampling with replacement)
            # For example: k0=33 and dof=5 yields 501942 possible
            # choices and a computation time of several seconds (on a 16GB RAM machine)
            if (exists(x = "exact.observed")) {
                if (adapt && is.imbalanced) {
                    p.val.test.exact <- apply(
                        env, 1, multiNom,
                        new.class.frequency$freq, batch
                    )
                } else {
                    p.val.test.exact <- apply(
                        env, 1, multiNom,
                        class.frequency$freq, batch
                    )
                }
                p.val.test.exact.null <- apply(apply(
                    batch.shuff, 2,
                    function(x, freq, envir) {
                        apply(envir, 1, FUN = multiNom, freq, x)
                    },
                    class.frequency$freq, env
                ), 1, mean, na.rm = TRUE)
                # apply(env, 1, multiNom, class.frequency$freq, batch.shuff)

                exact.expected[i] <- sum(p.val.test.exact.null < alpha, na.rm = TRUE) / testSize
                exact.observed[i] <- sum(p.val.test.exact < alpha, na.rm = TRUE) / testSize
                # compute the significance level for the number of rejected data points
                exact.signif[i] <-
                    1 - ptnorm(exact.observed[i],
                        mu = exact.expected[i],
                        sd = sqrt(exact.expected[i] * (1 - exact.expected[i]) / testSize),
                        alpha = alpha
                    )
                # p-value distribution
                rejection$results$exact.pvalue.test[idx.runs] <- p.val.test.exact
                rejection$results$exact.pvalue.null[idx.runs] <- p.val.test.exact.null
            }
        }


        if (n_repeat > 1) {
            # summarize chi2-results
            CI95 <- c(0.025, 0.5, 0.975)
            rejection$summary$kBET.expected <- c(
                mean(kBET.expected, na.rm = TRUE),
                quantile(kBET.expected, CI95,
                    na.rm = TRUE
                )
            )
            rownames(rejection$summary) <- c("mean", "2.5%", "50%", "97.5%")
            rejection$summary$kBET.observed <- c(
                mean(kBET.observed, na.rm = TRUE),
                quantile(kBET.observed, CI95,
                    na.rm = TRUE
                )
            )
            rejection$summary$kBET.signif <- c(
                mean(kBET.signif, na.rm = TRUE),
                quantile(kBET.signif, CI95,
                    na.rm = TRUE
                )
            )
            # summarize lrt-results
            rejection$summary$lrt.expected <- c(
                mean(lrt.expected, na.rm = TRUE),
                quantile(lrt.expected, CI95,
                    na.rm = TRUE
                )
            )
            rejection$summary$lrt.observed <- c(
                mean(lrt.observed, na.rm = TRUE),
                quantile(lrt.observed, CI95,
                    na.rm = TRUE
                )
            )
            rejection$summary$lrt.signif <- c(
                mean(lrt.signif, na.rm = TRUE),
                quantile(lrt.signif, CI95,
                    na.rm = TRUE
                )
            )
            # summarize exact test results
            if (exists(x = "exact.observed")) {
                rejection$summary$exact.expected <- c(
                    mean(exact.expected, na.rm = TRUE),
                    quantile(exact.expected, CI95,
                        na.rm = TRUE
                    )
                )
                rejection$summary$exact.observed <- c(
                    mean(exact.observed, na.rm = TRUE),
                    quantile(exact.observed, CI95,
                        na.rm = TRUE
                    )
                )
                rejection$summary$exact.signif <- c(
                    mean(exact.signif, na.rm = TRUE),
                    quantile(exact.signif, CI95,
                        na.rm = TRUE
                    )
                )
            }

            if (n_repeat < 10) {
                sprintf(
                    "The quantile computation for %s, subset is not meaningful.",
                    n_repeat
                )
            }

            if (plot && exists(x = "exact.observed")) {
                plot.data <- data.frame(
                    class = rep(
                        c(
                            "kBET", "kBET (random)",
                            "lrt", "lrt (random)",
                            "exact", "exact (random)"
                        ),
                        each = n_repeat
                    ),
                    data = c(
                        kBET.observed, kBET.expected,
                        lrt.observed, lrt.expected,
                        exact.observed, exact.expected
                    )
                )
                g <- ggplot(plot.data, aes(class, data)) +
                    geom_boxplot() +
                    theme_bw() +
                    labs(x = "Test", y = "Rejection rate") +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                g
            }
            if (plot && !exists(x = "exact.observed")) {
                plot.data <- data.frame(
                    class = rep(
                        c(
                            "kBET", "kBET (random)",
                            "lrt", "lrt (random)"
                        ),
                        each = n_repeat
                    ),
                    data = c(
                        kBET.observed, kBET.expected,
                        lrt.observed, lrt.expected
                    )
                )
                g <- ggplot(plot.data, aes(class, data)) +
                    geom_boxplot() +
                    theme_bw() +
                    labs(x = "Test", y = "Rejection rate") +
                    scale_y_continuous(limits = c(0, 1)) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                g
            }
        } else { # i.e. no n_repeat
            rejection$summary$kBET.expected <- kBET.expected
            rejection$summary$kBET.observed <- kBET.observed
            rejection$summary$kBET.signif <- kBET.signif

            if (addTest) {
                rejection$summary$lrt.expected <- lrt.expected
                rejection$summary$lrt.observed <- lrt.observed
                rejection$summary$lrt.signif <- lrt.signif
                if (exists(x = "exact.observed")) {
                    rejection$summary$exact.expected <- exact.expected
                    rejection$summary$exact.observed <- exact.observed
                    rejection$summary$exact.signif <- exact.signif
                }
            }
        }
    } else { # kBET only
        for (i in seq_len(n_repeat)) {
            # choose a random sample from dataset
            # (rows: samples, columns: parameters)
            idx.runs <- sample.int(dim.dataset[1], size = testSize)
            env <- cbind(knn[idx.runs, seq_len(k0 - 1)], idx.runs)

            # perform test
            cf <- if (adapt && is.imbalanced) new.class.frequency else class.frequency
            p.val.test <- apply(env, 1, chi_batch_test, cf, batch, dof)

            is.rejected <- p.val.test < alpha

            p.val.test.null <- apply(
                batch.shuff, 2,
                function(x) apply(env, 1, chi_batch_test, class.frequency, x, dof)
            )
            # p.val.test.null <- apply(env, 1, FUN = chi_batch_test,
            # class.frequency, batch.shuff, dof)

            # summarise test results
            # kBET.expected[i] <- sum(p.val.test.null < alpha) / length(p.val.test.null)
            kBET.expected[i] <- mean(apply(
                p.val.test.null, 2,
                function(x) sum(x < alpha, na.rm = TRUE) / sum(!is.na(x))
            ))

            kBET.observed[i] <- sum(is.rejected, na.rm = TRUE) / sum(!is.na(p.val.test))

            # compute significance
            kBET.signif[i] <- 1 - ptnorm(
                kBET.observed[i],
                mu = kBET.expected[i],
                sd = sqrt(kBET.expected[i] * (1 - kBET.expected[i]) / testSize),
                alpha = alpha
            )
            # assign results to result table
            rejection$results$tested[idx.runs] <- 1
            rejection$results$kBET.pvalue.test[idx.runs] <- p.val.test
            rejection$results$kBET.pvalue.null[idx.runs] <- rowMeans(p.val.test.null, na.rm = TRUE)
        }

        if (n_repeat > 1) {
            # summarize chi2-results
            CI95 <- c(0.025, 0.5, 0.975)
            rejection$summary$kBET.expected <- c(
                mean(kBET.expected, na.rm = TRUE),
                quantile(kBET.expected,
                    CI95,
                    na.rm = TRUE
                )
            )
            rownames(rejection$summary) <- c("mean", "2.5%", "50%", "97.5%")
            rejection$summary$kBET.observed <- c(
                mean(kBET.observed, na.rm = TRUE),
                quantile(kBET.observed, CI95,
                    na.rm = TRUE
                )
            )
            rejection$summary$kBET.signif <- c(
                mean(kBET.signif, na.rm = TRUE),
                quantile(kBET.signif, CI95,
                    na.rm = TRUE
                )
            )

            # return also n_repeat
            rejection$stats$kBET.expected <- kBET.expected
            rejection$stats$kBET.observed <- kBET.observed
            rejection$stats$kBET.signif <- kBET.signif

            if (n_repeat < 10) {
                sprintf(
                    "Quantile computation for %s subset results is not meaningful.",
                    n_repeat
                )
            }
            if (plot) {
                plot.data <- data.frame(
                    class = rep(
                        c(
                            "observed(kBET)",
                            "expected(random)"
                        ),
                        each = n_repeat
                    ),
                    data = c(kBET.observed, kBET.expected)
                )
                g <- ggplot(plot.data, aes(class, data)) +
                    geom_boxplot() +
                    theme_bw() +
                    labs(x = "Test", y = "Rejection rate") +
                    scale_y_continuous(limits = c(0, 1))
                g
            }
        } else {
            rejection$summary$kBET.expected <- kBET.expected[1]
            rejection$summary$kBET.observed <- kBET.observed[1]
            rejection$summary$kBET.signif <- kBET.signif[1]
        }
    }
    # collect parameters
    rejection$params <- list()
    rejection$params$k0 <- k0
    rejection$params$testSize <- testSize
    rejection$params$do.pca <- do.pca
    rejection$params$dim.pca <- dim.pca
    rejection$params$heuristic <- heuristic
    rejection$params$n_repeat <- n_repeat
    rejection$params$alpha <- alpha
    rejection$params$addTest <- addTest
    rejection$params$verbose <- verbose
    rejection$params$plot <- plot

    # add outsiders
    if (adapt) {
        rejection$outsider <- list()
        rejection$outsider$index <- outsider
        rejection$outsider$categories <- table(batch[outsider])
        rejection$outsider$p.val <- p.out
    }
    rejection
}
