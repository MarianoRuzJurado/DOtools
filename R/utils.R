#' @title Download example dataset 10x
#' @import cli
#' @import curl
#'
#' @return directory path where the data was saved
#'
#' @keywords internal
.example_10x <- function(base_dir = tempfile("dotools_datasets_")) {
    healthy_path <- file.path(base_dir, "healthy", "outs")
    disease_path <- file.path(base_dir, "disease", "outs")

    dir.create(healthy_path, recursive = TRUE, showWarnings = FALSE)
    dir.create(disease_path, recursive = TRUE, showWarnings = FALSE)

    healthy_link1 <-
        paste0(
            "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k",
            "_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5"
        )

    healthy_link2 <-
        paste0(
            "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k",
            "_protein_v3/pbmc_10k_protein_v3_raw_feature_bc_matrix.h5"
        )

    disease_link1 <-
        paste0(
            "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/malt_10k",
            "_protein_v3/malt_10k_protein_v3_filtered_feature_bc_matrix.h5"
        )

    disease_link2 <- paste0(
        "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/malt_10k",
        "_protein_v3/malt_10k_protein_v3_raw_feature_bc_matrix.h5"
    )
    links <- list(
        "healthy filtered" = healthy_link1,
        "healthy raw" = healthy_link2,
        "disease filtered" = disease_link1,
        "disease raw" = disease_link2
    )

    message("Downloading data to ", base_dir)

    for (name in names(links)) {
        link <- links[[name]]
        filename <- sub(".*10k_protein_v3_", "", link)
        dest_dir <- if (grepl("healthy", name)) healthy_path else disease_path
        dest_file <- file.path(dest_dir, filename)

        message(sprintf("Downloading %s to %s", name, dest_file))
        curl_download(url = link, destfile = dest_file, mode = "wb")
    }

    return(base_dir)
}


#' @keywords internal
.logger <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(timestamp, " - ", message)
}

#' @keywords internal
.suppressDeprecationWarnings <-
    function(expr,
                pattern = "PackageCheck\\(\\) was deprecated") {
        withCallingHandlers(
            expr,
            warning = function(w) {
                # Muffle lifecycle deprecation warnings
                if (inherits(w, "lifecycle_warning_deprecated")) {
                    rlang::cnd_muffle(w)
                }
                # Muffle warnings that match a provided pattern
                else if (!is.null(pattern) &&
                    grepl(pattern, conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                }
            }
        )
    }

# used for seurat warnings
#' @keywords internal
.suppressAllWarnings <- function(expr) {
    withCallingHandlers(
        expr,
        warning = function(w) {
            invokeRestart("muffleWarning")
            }
        )
    }

#

#' @keywords internal
.suppressPatternWarning <- function(expr, patterns) {
    withCallingHandlers(
    expr,
    warning = function(w) {
        # If warning message matches ANY supplied pattern → suppress it
        if (any(vapply(patterns, grepl, logical(1), x = conditionMessage(w)))) {
            invokeRestart("muffleWarning")
            }
        }
    )
}

theme_box <- function() {
    theme_bw() +
        theme(
            panel.border = element_rect(
                colour = "black",
                fill = NA,
                linewidth = 1
            ),
            panel.grid.major = element_line(
                colour = "grey90",
                linetype = "dotted"
            ),
            panel.grid.minor = element_line(
                colour = "grey90",
                linetype = "dotted"
            ),
            axis.line = element_line(colour = "black"),
            # facet_grid colors
            strip.background = element_rect(
                fill = "lightgrey",
                colour = "black",
                linewidth = 1
            ),
            strip.text = element_text(colour = "black", size = 12),
            # legend.background = element_rect(colour = "grey", fill = "white"),
            # legend.box.background = element_rect(colour = "grey", size = 0.5),
        )
}

#' @keywords internal
glmGamPoi_test <- function(
    sce_object,
    assay_normalized = "RNA",
    group_by = c("orig.ident","condition","annotation"),
    design_fit = NULL,
    annotation_key = "annotation",
    condition_key = "condition",
    reference = "ctrl"
){


  # support for Seurat objects
  if (methods::is(sce_object, "Seurat")) {
    DefaultAssay(sce_object) <- assay_normalized
    sce_object <- DOtools:::.suppressDeprecationWarnings(
      Seurat::as.SingleCellExperiment(sce_object,
        assay = assay_normalized
      )
    )
  }

  #Check if all supplied columns exist in metadata
  if(length(setdiff(group_by, names(SummarizedExperiment::colData(sce_object)))) > 0){
    stop("Not all group_by arguments are found in metadata!")
  }

  #catch for default design
  if (is.null(design_fit)) {
    design_fit <- ~ annotation + condition + condition:annotation - 1
  }

  #PB with specified columns
  sce_object_pb <- glmGamPoi::pseudobulk(
    sce_object,
    group_by = glmGamPoi::vars(!!!rlang::syms(group_by))
  )

  #checking if some cell types are not present in all the conditions
  tab <- table(sce_object_pb[[annotation_col]], sce_object_pb[[group_by]])
  bad_annotations <- rownames(tab)[rowSums(tab > 0) < ncol(tab)]

  if (length(bad_annotations) > 0) {
    .logger(paste0("Removing ", length(bad_annotations),
      " annotation level(s) not present in all conditions: ",
      paste(bad_annotations, collapse = ", "))
    )
  }

  #keep only valid ones
  keep_annotations <- rownames(tab)[rowSums(tab > 0) == ncol(tab)]
  sce_object_pb <- sce_object_pb[, sce_object_pb[[annotation_col]] %in% keep_annotations]

  DOtools:::.logger("Fitting Gamma-Poisson model...")
  #fit Gamma-Poisson model
  fit <- glmGamPoi::glm_gp(
    sce_object_pb, design = design_fit
  )

  comp <- setdiff(unique(sce_object[[condition_key]]), reference)

  # loop over comparisons and cell types
  de_collector <- data.frame()
  for (grp in comp) {
    if (!is.null(annotation_key)){
      for (celltype in unique(sce_object[[annotation_key]])) {

        #build the contrast string
        contrast <- paste0(
          "cond(", annotation_key, "='", celltype, "', ",
          condition_key, "='", grp, "') - ",
          "cond(", annotation_key, "='", celltype, "', ",
          condition_key, "='", reference, "')"
        )

        de_res <- glmGamPoi::test_de(fit, contrast = contrast)
        de_res[["celltype"]] <- celltype
        de_res[["condition"]] <- grp

        de_collector <- rbind(de_collector, de_res)
      }
    } else{
      #build the contrast string
      contrast <- paste0(
        "cond(", condition_key, "='", grp, "') - ",
        "cond(", condition_key, "='", reference, "')"
      )

      de_res <- glmGamPoi::test_de(fit, contrast = contrast)
      de_res[["condition"]] <- grp

      de_collector <- rbind(de_collector, de_res)

    }

  }

  return(de_collector)
}




umap_colors <- c(
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "tomato2",
    "#9467bd",
    "chocolate3",
    "#e377c2",
    "#ffbb78",
    "#bcbd22",
    "#17becf",
    "darkgoldenrod2",
    "#aec7e8",
    "#98df8a",
    "#ff9896",
    "#c5b0d5",
    "#c49c94",
    "#f7b6d2",
    "#c7c7c7",
    "#dbdb8d",
    "#9edae5",
    "sandybrown",
    "moccasin",
    "lightsteelblue",
    "darkorchid",
    "salmon2",
    "forestgreen",
    "bisque"
)
