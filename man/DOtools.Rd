\name{DOtools}
\alias{DOtools}
\title{DOtools: A Toolkit for scRNA Data Analysis}
\description{
    The `DOtools` package provides a set of functions for advanced data processing, visualisation, and statistical analysis in Seurat objects. It includes functions for cell-type prediction, reclustering, creating polished UMAP plots, subsetting Seurat objects, and various statistical analyses like Wilcoxon tests and SEM graphs.
}
\details{
    This package includes the following functions:
    \itemize{
        \item \code{\link{DO.BoxPlot}}: A function for creating box plots with Wilcoxon test results.
        \item \code{\link{DO.CellTypist}}: A function for running CellTypist on Seurat objects to predict cell types.
        \item \code{\link{DO.DietSeurat}}: A function for diet-based analysis of Seurat objects.
        \item \code{\link{DO.Dotplot}}: A function for creating dot plots for visualizing gene expression.
        \item \code{\link{DO.FullRecluster}}: A function for fine-grained reclustering of Seurat objects.
        \item \code{\link{DO.BarplotClustert}}: A function for generating mean and SEM graphs for cluster-based analysis with t-tests.
        \item \code{\link{DO.BarplotWilcox}}: A function for generating mean and SEM graphs with Wilcoxon test results.
        \item \code{\link{DO.Subset}}: A function for subsetting Seurat objects based on metadata.
        \item \code{\link{DO.UMAP}}: A function for creating polished UMAP plots using either DimPlot or FeaturePlot.
        \item \code{\link{DO.VlnPlot}}: A function for generating violin plots with Wilcoxon test results.
        \item \code{\link{DO.CellComposition}}: A function for visualizing and statistically analyzing cell-type composition changes across conditions using the Scanpro Python package, with support for bootstrapping, proportion plots, and customizable output.
        \item \code{\link{DO.Import}}: A function for building a merged Seurat object from 10x software output, or directly from provided tables.
        \item \code{\link{DO.CellBender}}: A function for running CellBender in a virtual conda env with provided raw count h5 files.
        \item \code{\link{DO.SplitBarGSEA}}: A function for viusalizing GSEA result from a provided df from e.g. metascape
        \item \code{\link{DO.scVI}}: A function for running the scVI Integration implemented in scvi-tools.
        \item \code{\link{DO.TransferLabel}}: A function for transfering annotation from a subseted object to the original seurat object.
        \item \code{\link{DO.PyEnv}}: A function for creating a conda envrionment holding all python packages needed for some functions.
        \item \code{\link{DO.Correlation}}: A function for creating a correlation plot between provided samples in the category specified.
        \item \code{\link{DO.MultiDGE}}: A function for calculating DEGs on a single cell and speudo bulk level.
        \item \code{\link{dot-Do.BarcodeRanks}}: A function for estimating the number of expected cells and droplets.
        \item \code{\link{dot-QC.Vlnplot}}: A function for estimating the number of expected cells and droplets.
    }
}
\author{Mariano Ruz Jurado, David Rodriguez Morales}
\seealso{
    \code{\link{DO.Box.Plot.wilcox}}, \code{\link{DO.CellTypist}}, \code{\link{DO.DietSeurat}}, \code{\link{DO.Dotplot}},
    \code{\link{DO.FullRecluster}}, \code{\link{DO.Mean.SEM.Graphs.cluster.t}}, \code{\link{DO.Mean.SEM.Graphs.wilcox}},
    \code{\link{DO.Subset}}, \code{\link{DO.UMAP}}, \code{\link{DO.Vln.Plot.wilcox}}, \code{\link{DO.Import}}, \code{\link{DO.CellBender}},
    \code{\link{DO.SplitBarGSEA}}, \code{\link{DO.scVI}}, \code{\link{DO.TransferLabel}}, \code{\link{DO.PyEnv}}, \code{\link{DO.Correlation}},
    \code{\link{DO.MultiDGE}}, \code{\link{dot-Do.BarcodeRanks}}, \code{\link{dot-QC.Vlnplot}}
}
