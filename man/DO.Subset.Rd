% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_functions.R
\name{DO.Subset}
\alias{DO.Subset}
\title{DO.Subset}
\usage{
DO.Subset(
  Seu_object,
  assay = "RNA",
  ident,
  ident_name = NULL,
  ident_thresh = NULL
)
}
\arguments{
\item{Seu_object}{The seurat object}

\item{assay}{assay to subset by}

\item{ident}{meta data column to subset for}

\item{ident_name}{name of group of barcodes in ident of subset for}

\item{ident_thresh}{numeric thresholds as character, e.g ">5" or c(">5", "<200"), to subset barcodes in ident for}
}
\value{
a subsetted Seurat object
}
\description{
Creates a subset of a Seurat object based on either categorical or numeric thresholds in metadata.
Allows for subsetting by specifying the ident column, group name, or threshold criteria.
Ideal for extracting specific cell populations or clusters based on custom conditions.
Returns a new Seurat object containing only the subsetted cells and does not come with the Seuratv5 subset issue
Please be aware that right now, after using this function the subset might be treated with Seuv5=False in other functions.
}
\examples{
sc_data <- readRDS(system.file("extdata", "sc_data.rds", package = "DOtools"))

sc_data_sub <- DO.Subset(
  Seu_object = sc_data,
  ident="condition",
  ident_name="ctrl"
)



}
\author{
Mariano Ruz Jurado
}
