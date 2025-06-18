#' @title Download example dataset 10x
#' @import cli
#' @import curl
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

  message("📥 Downloading data to ", base_dir)

  for (name in names(links)) {
    link <- links[[name]]
    filename <- sub(".*10k_protein_v3_", "", link)
    dest_dir <- if (grepl("healthy", name)) healthy_path else disease_path
    dest_file <- file.path(dest_dir, filename)

    message(sprintf("⬇️  Downloading %s to %s", name, dest_file))
    curl_download(url = link, destfile = dest_file, mode = "wb")
  }

  return(base_dir)
}
