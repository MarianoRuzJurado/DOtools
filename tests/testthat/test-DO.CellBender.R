# tests/testthat/test-DO.CellBender.R
library(testthat)
library(mockery)
library(SingleCellExperiment)
library(S4Vectors)

# ---- Test .DO.BarcodeRanks using a real barcodeRanks-friendly SCE ----
test_that(".DO.BarcodeRanks returns correct names and numeric values", {
  # Create SCE with column sums that are distinct so barcodeRanks can compute knee/inflection
  # Use 1 gene (row) with different counts per column (cell)
  counts <- matrix(c(1000, 900, 800, 50, 20, 5), nrow = 1)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

  res <- DOtools:::.DO.BarcodeRanks(sce)

  expect_named(res, c("xpc_cells", "total_cells"))
  # Accept integer or double (both are numeric in R)
  expect_true(is.numeric(res["xpc_cells"]))
  expect_true(is.numeric(res["total_cells"]))
  # Values should be non-negative
  expect_true(as.numeric(res["xpc_cells"]) >= 0)
  expect_true(as.numeric(res["total_cells"]) >= 0)
})

# ---- Test DO.CellBender behavior for missing paths ----
test_that("DO.CellBender errors if paths missing (stopifnot)", {
  # Passing non-existent paths should trigger stopifnot -> expect_error
  expect_error(
    DOtools::DO.CellBender(cellranger_path = tempfile(), output_path = tempfile())
  )
})

# ---- Test DO.CellBender builds correct command with BarcodeRanking = TRUE ----
test_that("DO.CellBender builds commands (BarcodeRanking=TRUE) and includes expected args", {
  # Create a fake CellRanger folder with one sample subdir and dummy h5 file
  base_dir <- tempfile("cellr_")
  dir.create(base_dir)
  sample_dir <- file.path(base_dir, "SampleA")
  dir.create(file.path(sample_dir, "outs"), recursive = TRUE)
  h5_path <- file.path(sample_dir, "outs", "raw_feature_bc_matrix.h5")
  file.create(h5_path)

  output_dir <- tempfile("out_")
  dir.create(output_dir)

  bash_script <- tempfile("cb_sh_"); file.create(bash_script)

  # Prepare a place to capture system2 calls
  calls <- list()
  fake_system2 <- function(cmd, args, stdout = TRUE, stderr = TRUE) {
    calls[[length(calls) + 1]] <<- list(cmd = cmd, args = args)
    # mimic success
    return(0)
  }

  # Stub inside DO.CellBender:
  stub(DOtools::DO.CellBender, "system2", fake_system2)
  stub(DOtools::DO.CellBender, "DropletUtils::read10xCounts",
       function(h5) SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, nrow = 3, ncol = 5))))
  stub(DOtools::DO.CellBender, ".DO.BarcodeRanks", function(x) c(10, 20))
  stub(DOtools::DO.CellBender, ".logger", function(...) invisible(NULL))

  # choose a conda path that does not exist to exercise the "create conda env" branch
  tmp_conda <- file.path(tempdir(), "fake_conda_env_nonexistent")
  if (dir.exists(tmp_conda)) unlink(tmp_conda, recursive = TRUE, force = TRUE)

  # Run the function (no external commands executed because of stubs)
  DOtools::DO.CellBender(
    cellranger_path = base_dir,
    output_path = output_dir,
    samplenames = NULL,
    cuda = TRUE,
    cpu_threads = 4,
    epochs = 10,
    lr = 1e-5,
    estimator_multiple_cpu = TRUE,
    log = TRUE,
    conda_path = tmp_conda,
    BarcodeRanking = TRUE,
    bash_script = bash_script
  )

  # There should be at least one system2 call that attempts to run the bash script
  run_calls <- Filter(function(x) any(grepl(bash_script, as.character(x$args), fixed = TRUE)), calls)
  expect_true(length(run_calls) >= 1, info = "At least one system2 call should include the bash script")

  # Check that the run call contains expected extra flags
  found_expected <- FALSE
  for (rc in run_calls) {
    args_chr <- as.character(rc$args)
    if (any(grepl("--expected-cells", args_chr, fixed = TRUE))) {
      expect_true(any(grepl("10", args_chr)), info = "expected-cells value should be present")
      expect_true(any(grepl("--total-droplets", args_chr, fixed = TRUE)), info = "total-droplets flag present")
      expect_true(any(grepl("20", args_chr, fixed = TRUE)), info = "total-droplets value present")
      expect_true(any(grepl("--cuda", args_chr, fixed = TRUE)), info = "cuda flag present")
      expect_true(any(grepl("--log", args_chr, fixed = TRUE)), info = "log flag present")
      expect_true(any(grepl("--estimator_multiple_cpu", args_chr, fixed = TRUE)), info = "estimator_multiple_cpu flag present")
      found_expected <- TRUE
    }
  }
  expect_true(found_expected, info = "Run command with expected flags was not found in system2 calls")
})

# ---- Test DO.CellBender builds command when BarcodeRanking = FALSE and uses existing conda path ----
test_that("DO.CellBender without BarcodeRanking builds a simpler command and respects existing conda_path", {
  base_dir <- tempfile("cellr2_"); dir.create(base_dir)
  sample_dir <- file.path(base_dir, "SampleB"); dir.create(file.path(sample_dir, "outs"), recursive = TRUE)
  file.create(file.path(sample_dir, "outs", "raw_feature_bc_matrix.h5"))
  output_dir <- tempfile("out2_"); dir.create(output_dir)
  bash_script <- tempfile("cb_sh2_"); file.create(bash_script)

  existing_conda <- tempfile("conda_exist_"); dir.create(existing_conda)

  calls2 <- list()
  fake_system2_b <- function(cmd, args, stdout = TRUE, stderr = TRUE) {
    calls2[[length(calls2) + 1]] <<- list(cmd = cmd, args = args)
    return(0)
  }

  stub(DOtools::DO.CellBender, "system2", fake_system2_b)
  stub(DOtools::DO.CellBender, "DropletUtils::read10xCounts",
       function(h5) SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix(1, nrow = 2, ncol = 3))))
  stub(DOtools::DO.CellBender, ".DO.BarcodeRanks", function(x) c(5, 8))
  stub(DOtools::DO.CellBender, ".logger", function(...) invisible(NULL))

  DOtools::DO.CellBender(
    cellranger_path = base_dir,
    output_path = output_dir,
    samplenames = NULL,
    cuda = TRUE,
    cpu_threads = 2,
    epochs = 20,
    lr = 1e-6,
    estimator_multiple_cpu = FALSE,
    log = FALSE,
    conda_path = existing_conda,
    BarcodeRanking = FALSE,
    bash_script = bash_script
  )

  run_calls <- Filter(function(x) any(grepl(bash_script, as.character(x$args), fixed = TRUE)), calls2)
  expect_true(length(run_calls) >= 1)
  for (rc in run_calls) {
    args_chr <- as.character(rc$args)
    expect_false(any(grepl("--expected-cells", args_chr)))
    expect_false(any(grepl("--total-droplets", args_chr)))
    expect_true(any(grepl("--cuda", args_chr)))
  }
})
