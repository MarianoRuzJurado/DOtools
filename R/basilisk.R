#' Cellbender Python Environment
#'
#' A `BasiliskEnvironment` containing Cellbender and Python 3.7, used in the DOtools package for running Cellbender-related functionality.
#'
#' @return A [basilisk::BasiliskEnvironment()] object preconfigured with Cellbender.
#' @export
#TODO write this function correctly for creating the Cellbender env maybe use that function in cellbender call and if its not there create it
cellbender_env <- basilisk::BasiliskEnvironment(
  envname = "cellbender",
  pkgname = "DOtools",
  packages = c("python=3.7.16",
               "cellbender=0.3.2")
)

#' DOtools Python Environment
#'
#' A `BasiliskEnvironment` containing general-purpose Python packages used in DOtools, including celltypist and scanpro, with Python 3.11.
#'
#' @return A [basilisk::BasiliskEnvironment()] object preconfigured with DOtools-related packages.
#' @export
#TODO write this function correctly for creating the remaining python function env maybe use that function in cellbender call and if its not there create it
dotools_env <- basilisk::BasiliskEnvironment(
  envname = "dotools",
  pkgname = "DOtools",
  packages = c("python=3.11.11",
               "celltypist=1.6.3",
               "scanpro=0.3.2")
)
