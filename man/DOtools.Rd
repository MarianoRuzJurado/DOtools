\name{DOtools}
\alias{DOtools}
\title{DOtools: A Toolkit for Seurat Data Analysis}
\description{
    The `DOtools` package provides a set of functions for advanced data processing, visualization, and statistical analysis in Seurat objects. It includes functions for cell-type prediction, reclustering, creating polished UMAP plots, subsetting Seurat objects, and various statistical analyses like Wilcoxon tests and SEM graphs.
}
\details{
    This package includes the following functions:
    \itemize{
        \item \code{\link{DO.Box.Plot.wilcox}}: A function for creating box plots with Wilcoxon test results.
        \item \code{\link{DO.CellTypist}}: A function for running CellTypist on Seurat objects to predict cell types.
        \item \code{\link{DO.DietSeurat}}: A function for diet-based analysis of Seurat objects.
        \item \code{\link{DO.Dotplot}}: A function for creating dot plots for visualizing gene expression.
        \item \code{\link{DO.FullRecluster}}: A function for fine-grained reclustering of Seurat objects.
        \item \code{\link{DO.Mean.SEM.Graphs.cluster.t}}: A function for generating mean and SEM graphs for cluster-based analysis with t-tests.
        \item \code{\link{DO.Mean.SEM.Graphs.wilcox}}: A function for generating mean and SEM graphs with Wilcoxon test results.
        \item \code{\link{DO.Subset}}: A function for subsetting Seurat objects based on metadata.
        \item \code{\link{DO.UMAP}}: A function for creating polished UMAP plots using either DimPlot or FeaturePlot.
        \item \code{\link{DO.Vln.Plot.wilcox}}: A function for generating violin plots with Wilcoxon test results.
        \item \code{\link{dot-annoSegment}}: A function for adding annotation segments to dot plots.
    }
}
\author{Mariano Ruz Jurado}
\seealso{
    \code{\link{DO.Box.Plot.wilcox}}, \code{\link{DO.CellTypist}}, \code{\link{DO.DietSeurat}}, \code{\link{DO.Dotplot}},
    \code{\link{DO.FullRecluster}}, \code{\link{DO.Mean.SEM.Graphs.cluster.t}}, \code{\link{DO.Mean.SEM.Graphs.wilcox}},
    \code{\link{DO.Subset}}, \code{\link{DO.UMAP}}, \code{\link{DO.Vln.Plot.wilcox}}, \code{\link{dot-annoSegment}}
}
