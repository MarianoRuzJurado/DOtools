#' @title Download example dataset 10x
#' @import cli
#' @import curl
#'
#' @keywords internal
.example_10x <- function() {
  dir.create("/tmp/dotools_datasets/healthy/outs", recursive = TRUE, showWarnings = FALSE)
  dir.create("/tmp/dotools_datasets/disease/outs", recursive = TRUE, showWarnings = FALSE)

  healthy_path <- "/tmp/dotools_datasets/healthy/outs/"
  disease_path <- "/tmp/dotools_datasets/disease/outs/"

  healthy_link1 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.tar.gz"
  healthy_link2 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_raw_feature_bc_matrix.tar.gz"
  disease_link1 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/malt_10k_protein_v3/malt_10k_protein_v3_filtered_feature_bc_matrix.tar.gz"
  disease_link2 <- "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/malt_10k_protein_v3/malt_10k_protein_v3_raw_feature_bc_matrix.tar.gz"

  links <- list(
    "healthy filtered" = healthy_link1,
    "healthy raw" = healthy_link2,
    "disease filtered" = disease_link1,
    "disease raw" = disease_link2
  )

  message("📥 Downloading data to /tmp/dotools_datasets/")

  for (name in names(links)) {
    link <- links[[name]]
    filename <- sub(".*10k_protein_v3_", "", link)
    dest_dir <- if (grepl("healthy", name)) healthy_path else disease_path
    dest_file <- file.path(dest_dir, filename)

    message(sprintf("⬇️  Downloading %s to %s", name, dest_file))
    curl_download(url = link, destfile = dest_file, mode = "wb")
  }

  invisible(NULL)
}
