# DOtools ![](reference/figures/LogoDoTools.png)

## Overview

DOtools is a user-friendly R package designed to streamline common
workflows in single-cell RNA sequencing (scRNA-seq) data analysis using
the Seurat ecosystem and third-party tools such as scVI, CellTypist, and
CellBender. It provides high-level wrappers and visualisation utilities
to help efficiently preprocess, analyze, and interpret single-cell data.

# **Installation**

DOtools is currently available through Bioconductor. To install the
package, start R and run:

``` r
install.packages("BiocManager")
BiocManager::install("DOtools")
```

## **Requirements**

Some functions in this package depend on Python packages or scripts.
These functions use
[basilisk](https://www.bioconductor.org/packages/release/bioc/html/basilisk.html)
to create isolated environments, which are built using conda. As a
result, a working conda installation is required. The DoTools pipeline
was tested on Linux and MacOS and an example dataset with ~65k cells can
be processed with a machine with at least 16GBs of RAM memomy and 5
CPUs.

## **Python package**

If you prefer python over R for data analysis, we also provide a python
version of the DOtools package. Please refer for the python version to:
[DOtools_py](https://github.com/davidrm-bio/DOTools_py)

## **Contribution Guidelines**

Raising up an issue in this Github repository might be the fastest way
of submitting suggestions and bugs. Alternatively you can write an
email: <ruzjurado@med.uni-frankfurt.de>

## **Citation**

tba
