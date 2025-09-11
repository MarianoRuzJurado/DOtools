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

  healthy_link1 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5"
  healthy_link2 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_raw_feature_bc_matrix.h5"
  disease_link1 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/malt_10k_protein_v3/malt_10k_protein_v3_filtered_feature_bc_matrix.h5"
  disease_link2 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/malt_10k_protein_v3/malt_10k_protein_v3_raw_feature_bc_matrix.h5"

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
  message(paste0(timestamp, " - ", message, "\n"))
}

#' @keywords internal
.suppressDeprecationWarnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (inherits(w, "lifecycle_warning_deprecated")) {
        rlang::cnd_muffle(w)
      }
    }
  )
}


theme_box <- function(){
  theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(colour = "grey90", linetype = "dotted"),
      panel.grid.minor = element_line(colour = "grey90", linetype = "dotted"),
      axis.line = element_line(colour = "black"),
      #facet_grid colors
      strip.background = element_rect(fill = "lightgrey", colour = "black", linewidth = 1),
      strip.text = element_text(colour = "black", size = 12),
      # legend.background = element_rect(colour = "grey", fill = "white"),
      # legend.box.background = element_rect(colour = "grey", size = 0.5),
    )
}


umap_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "tomato2", "#9467bd", "chocolate3","#e377c2", "#ffbb78", "#bcbd22",
  "#17becf","darkgoldenrod2", "#aec7e8", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94","#f7b6d2", "#c7c7c7", "#dbdb8d",
  "#9edae5","sandybrown","moccasin","lightsteelblue","darkorchid","salmon2","forestgreen","bisque"
)
