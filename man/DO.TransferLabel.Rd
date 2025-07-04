% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_functions.R
\name{DO.TransferLabel}
\alias{DO.TransferLabel}
\title{DO.TransferLabel}
\usage{
DO.TransferLabel(Seu_obj, Subset_obj, annotation_column, subset_annotation)
}
\arguments{
\item{Seu_obj}{Seurat object with annotation in meta.data}

\item{Subset_obj}{subsetted Seurat object with re-annotated clusters}

\item{annotation_column}{column name in meta.data with annotation}

\item{subset_annotation}{column name in meta.data with annotation in the subsetted object}
}
\value{
Seurat Object with transfered labels
}
\description{
Transfers cell-type annotations from a re-annotated subset of a Seurat object
back to the full Seurat object. This is useful when clusters have been refined
or re-labeled in a subset and need to be reflected in the original object.
}
\examples{
sc_data <- readRDS(system.file("extdata", "sc_data.rds", package = "DOtools"))

sc_data <- DO.TransferLabel(sc_data,
                            sc_data,
                            annotation_column="annotation_refined",
                            subset_annotation="annotation_refined"
                           )



}
