% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_functions.R
\name{DO.DietSeurat}
\alias{DO.DietSeurat}
\title{Remove Layers from Seurat Object by Pattern}
\usage{
DO.DietSeurat(Seu_object, pattern = "^scale\\\\.data\\\\.")
}
\arguments{
\item{Seu_object}{Seurat object.}

\item{pattern}{regular expression pattern to match layer names. Default "^scale\\.data\\."}
}
\value{
Seurat object with specified layers removed.
}
\description{
This function removes layers from a Seurat object's RNA assay based on a specified regular expression pattern.
It is supposed to make something similar than the no longer working DietSeurat function, by removing no longer needed layers from th object.
}
\examples{
sc_data <- readRDS(system.file("extdata", "sc_data.rds", package = "DOtools"))

sc_data <- DO.DietSeurat(sc_data, pattern = "data")


}
\author{
Mariano Ruz Jurado
}
