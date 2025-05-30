% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_functions.R
\name{DO.CellTypist}
\alias{DO.CellTypist}
\title{DO Celltypist}
\usage{
DO.CellTypist(
  Seu_object,
  modelName = "Healthy_Adult_Heart.pkl",
  minCellsToRun = 200,
  runCelltypistUpdate = TRUE,
  over_clustering = "seurat_clusters",
  assay_normalized = "RNA",
  returnProb = FALSE,
  SeuV5 = T
)
}
\arguments{
\item{Seu_object}{The seurat object}

\item{minCellsToRun}{If the input seurat object has fewer than this many cells, NAs will be added for all expected columns and celltypist will not be run.}

\item{runCelltypistUpdate}{If true, --update-models will be run for celltypist prior to scoring cells.}

\item{over_clustering}{Column in metadata in object with clustering assignments for cells, default seurat_clusters}

\item{assay_normalized}{Assay with log1p normalized expressions}

\item{returnProb}{will additionally return the probability matrix, return will give a list with the first element beeing the object and second prob matrix}
}
\value{
a seurat
}
\description{
Runs the CellTypist model on a Seurat object to predict cell type labels, storing the results as metadata.
If the number of cells is less than the specified threshold, it returns NAs for the labels.
Optionally updates the CellTypist models and returns the probability matrix.
Useful for annotating cell types in single-cell RNA sequencing datasets.
}
\examples{
\dontrun{


DO.CellTypist(
  Seu_object = Seurat,
  modelName = "Healthy_Adult_Heart.pkl",
  runCelltypistUpdate = TRUE,
  over_clustering = "seurat_clusters",
  SeuV5=T
)
}

}
\author{
Mariano Ruz Jurado
}
