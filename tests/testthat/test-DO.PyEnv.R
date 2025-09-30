library(testthat)
library(mockery)  # for mocking functions

test_that("DO.PyEnv sets default path when conda_path is NULL", {
  # Mock dir.exists to always return FALSE
  stub(DO.PyEnv, "dir.exists", function(path) FALSE)
  # Mock system2 to just record calls
  calls <- list()
  stub(DO.PyEnv, "system2", function(command, args, stdout, stderr) {
    calls <<- append(calls, list(list(command = command, args = args)))
    return(0)  # pretend success
  })

  DO.PyEnv()

  # Check that system2 was called three times
  expect_equal(length(calls), 3)

  # Check that the first system2 call contains "create"
  expect_true(any(grepl("create", calls[[1]]$args)))
})

test_that("DO.PyEnv expands provided path", {
  test_path <- "~/my_env_test"

  stub(DO.PyEnv, "dir.exists", function(path) TRUE)  # pretend env exists
  stub(DO.PyEnv, ".logger", function(msg) msg)       # avoid printing

  expect_silent(DO.PyEnv(conda_path = test_path))
})

test_that("DO.PyEnv throws error if system2 fails", {
  stub(DO.PyEnv, "dir.exists", function(path) FALSE)
  stub(DO.PyEnv, "system2", function(command, args, stdout, stderr) stop("system failure"))

  expect_error(DO.PyEnv(), "Failed to create conda environment")
})
