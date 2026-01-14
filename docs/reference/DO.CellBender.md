# DO.CellBender

This function wraps a system call to a bash script for running
CellBender on CellRanger outputs. It ensures required inputs are
available and optionally installs CellBender in a conda env.

## Usage

``` r
DO.CellBender(
  cellranger_path,
  output_path,
  samplenames = NULL,
  cuda = TRUE,
  cpu_threads = 15,
  epochs = 150,
  lr = 1e-05,
  estimator_multiple_cpu = FALSE,
  log = TRUE,
  conda_path = NULL,
  BarcodeRanking = TRUE,
  bash_script = system.file("bash", "_run_CellBender.sh", package = "DOtools")
)
```

## Arguments

- cellranger_path:

  Path to folder with CellRanger outputs.

- output_path:

  Output directory for CellBender results.

- samplenames:

  Optional vector of sample names. If NULL, will autodetect folders in
  `cellranger_path`.

- cuda:

  Logical, whether to use GPU (CUDA).

- cpu_threads:

  Number of CPU threads to use.

- epochs:

  Number of training epochs.

- lr:

  Learning rate.

- estimator_multiple_cpu:

  Use estimator with multiple CPU threads.

- log:

  Whether to enable logging.

- conda_path:

  Optional path to the conda environment.

- BarcodeRanking:

  Optional Calculation of estimated cells in samples through
  DropletUtils implementation

- bash_script:

  Path to the bash script that runs CellBender.

## Value

None

## Examples

``` r
if (FALSE) { # \dontrun{
# Define paths
cellranger_path <- "/mnt/data/cellranger_outputs"
output_path <- "/mnt/data/cellbender_outputs"

# Optional: specify sample names if automatic detection is not desired
samplenames <- c("Sample_1", "Sample_2")

# Run CellBender (uses GPU by default)
DO.CellBender(
    cellranger_path = cellranger_path,
    output_path = output_path,
    samplenames = samplenames,
    cuda = TRUE,
    cpu_threads = 8,
    epochs = 100,
    lr = 0.00001,
    estimator_multiple_cpu = FALSE,
    log = TRUE
)
} # }
```
