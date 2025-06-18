## ----chunk_setup, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)



## ----vignette_setup, echo=FALSE, message=FALSE, warning = FALSE---------------
# Track time spent on making the vignette.
start_time <- Sys.time()

# Bib setup.
library(RefManageR)

# Write bibliography information
bib <- c(
    DOtools = citation("DOtools")[1],
    scDbl = citation("scDblFinder")[2],
    Seurat = citation("Seurat")[3]
)

## ----bioconductor_install, eval=FALSE-----------------------------------------
# install.packages("BiocManager") # WORK iN PROGRESS
# BiocManager::install("DOtools")

## ----github_install, eval=FALSE-----------------------------------------------
# install.packages("devtools")
# devtools::install_github("MarianoRuzJurado/DOtools")

## ----load_library, message=FALSE----------------------------------------------
library(DOtools)

#Additional packages
library(Seurat)
library(plyr)
library(dplyr)
library(tibble)
library(enrichR)

#Python installation set up
DO.PyEnv()
reticulate::use_python("~/.venv/DOtools/bin/python")

