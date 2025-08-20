#' @title DO.PyEnv
#' @description Sets up or connects to a conda Python environment for use with DOtools.
#' If no environment path is provided, it will create one at `~/.venv/DOtools` and install required Python packages:
#' `scvi-tools`, `celltypist`, and `scanpro`.
#'
#' @param conda_path character string specifying the path to an existing or new conda environment.
#' @return None
#'
#' @examples
#' # Automatically create DOtools environment at ~/.venv/DOtools if it doesn't exist
#' DO.PyEnv()
#'
#' # Use an existing conda environment at a custom location
#' DO.PyEnv(conda_path = "~/miniconda3/envs/my_dotools_env")
#'
#'
#' @export
DO.PyEnv <- function(conda_path = NULL) {

  # Handle conda environment
  if (is.null(conda_path)) {
    conda_path <- file.path(path.expand("~"), ".venv/DOtools")
  } else {
    conda_path <- path.expand(conda_path)
  }

  if (!dir.exists(conda_path)) {
    .logger("Creating conda environment for DOtools")
    conda_args1 <- c("conda","create","-y", "-p", file.path(path.expand("~"), ".venv/DOtools"), "python=3.11")
    conda_args2 <- c("conda","run", "-p", file.path(path.expand("~"), ".venv/DOtools"), "pip", "install", "scvi-tools==1.3.0", "celltypist==1.6.3", "scanpro==0.3.2")
    conda_args3 <- c("conda","run", "-p", file.path(path.expand("~"), ".venv/DOtools"), "pip", "install", "scipy==1.15.3")
    tryCatch({
      system2(conda_args1[1], args = conda_args1[-1], stdout = FALSE, stderr = FALSE)
      system2(conda_args2[1], args = conda_args2[-1], stdout = FALSE, stderr = FALSE)
      system2(conda_args3[1], args = conda_args3[-1], stdout = FALSE, stderr = FALSE)
    }, error = function(e) {
      stop("Failed to create conda environment. Provide a valid environment path or fix installation.")
    })
  } else {
    .logger(sprintf("Using existing conda environment at: %s", conda_path))
  }
  .logger("Python packages ready for DOtools!")
}
