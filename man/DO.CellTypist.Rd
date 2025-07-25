% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_functions.R
\name{DO.CellTypist}
\alias{DO.CellTypist}
\title{DO Celltypist}
\usage{
DO.CellTypist(
  sce_object,
  modelName = "Healthy_Adult_Heart.pkl",
  minCellsToRun = 200,
  runCelltypistUpdate = TRUE,
  over_clustering = "seurat_clusters",
  assay_normalized = "RNA",
  returnProb = FALSE,
  SeuV5 = TRUE
)
}
\arguments{
\item{sce_object}{The seurat or sce object}

\item{modelName}{Specify the model you want to use for celltypist}

\item{minCellsToRun}{If the input seurat or SCE object has fewer than this many cells, NAs will be added for all expected columns and celltypist will not be run.}

\item{runCelltypistUpdate}{If true, --update-models will be run for celltypist prior to scoring cells.}

\item{over_clustering}{Column in metadata in object with clustering assignments for cells, default seurat_clusters}

\item{assay_normalized}{Assay with log1p normalized expressions}

\item{returnProb}{will additionally return the probability matrix, return will give a list with the first element beeing the object and second prob matrix}

\item{SeuV5}{Specify if the Seurat object is made with Seuratv5}
}
\value{
a seurat or sce object
}
\description{
Runs the CellTypist model on a Seurat or SCE object to predict cell type labels, storing the results as metadata.
If the number of cells is less than the specified threshold, it returns NAs for the labels.
Optionally updates the CellTypist models and returns the probability matrix.
Useful for annotating cell types in single-cell RNA sequencing datasets.
}
\examples{
sce_data <- readRDS(system.file("extdata", "sce_data.rds", package = "DOtools"))


sce_data <- DO.CellTypist(
  sce_object = sce_data,
  modelName = "Healthy_Adult_Heart.pkl",
  runCelltypistUpdate = TRUE,
  over_clustering = "seurat_clusters",
  minCellsToRun=5,
  SeuV5=TRUE
)

}
\author{
Mariano Ruz Jurado
}
