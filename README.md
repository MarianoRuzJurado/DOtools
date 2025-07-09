# DOtools <img src="man/figures/LogoDoTools.png" align="right" width="240"/>

<!-- badges: start -->

<!-- [![BioC status](https://www.bioconductor.org/shields/build/release/bioc/DOtools.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/DOtools) -->

<!-- [![BioC dev status](https://www.bioconductor.org/shields/build/devel/bioc/DOtools.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/DOtools) -->

[![Issues](https://img.shields.io/github/issues/MarianoRuzJurado/DOtools)](https://github.com/MarianoRuzJurado/DOtools/issues) [![Stars](https://img.shields.io/github/stars/MarianoRuzJurado/DOtools?style=flat&logo=github&color=yellow)](https://github.com/MarianoRuzJurado/DOtools/stargazers)

<!-- badges: end -->

## Overview

DOtools is a user-friendly R package designed to streamline common workflows in single-cell RNA sequencing (scRNA-seq) data analysis using the Seurat ecosystem and third-party tools such as scVI, CellTypist, and CellBender.
It provides high-level wrappers and visualisation utilities to help efficiently preprocess, analyze, and interpret single-cell data.

# <b> Installation </b>

DOtools is currently available through github.
To install the package, start R and run:

``` ruby
install.packages("devtools") # if not installed already
devtools::install_github("MarianoRuzJurado/DOtools")
```

We highly recommend having [conda](https://www.anaconda.com/docs/getting-started/miniconda/main) installed for python environments.
Some functions in this package rely on python packages/scripts, please make sure you specify with reticulate a useable python.
We recommend to use our `DO.PyEnv()` function for easy-to-use solution.

## <b> Python package </b>

If you prefer python over R for data analysis, we also provide a python version of the DOtools package.
Please refer for the python version to: [DOtools_py](https://github.com/davidrm-bio/DOTools_py)

## <b> Contribution Guidelines </b>

Raising up an issue in this Github repository might be the fastest way of submitting suggestions and bugs.
Alternatively you can write an email: [ruzjurado\@med.uni-frankfurt.de](mailto:ruzjurado@med.uni-frankfurt.de)

## <b> Citation </b>

tba
