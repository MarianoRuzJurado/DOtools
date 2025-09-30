library(testthat)
library(mockery)
library(ggplot2)
library(curl)

test_that(".logger prints a timestamped message", {
  msg <- "Hello world"
  output <- capture_messages(.logger(msg))
  expect_true(any(grepl(msg, output)))
  expect_true(any(grepl("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}", output)))
})

test_that(".suppressDeprecationWarnings muffles warnings", {
  f <- function() warning("PackageCheck() was deprecated", call. = FALSE, class = "lifecycle_warning_deprecated")
  expect_silent(.suppressDeprecationWarnings(f()))
})

test_that("theme_box returns a ggplot2 theme", {
  th <- theme_box()
  expect_s3_class(th, "theme")
  expect_true("panel.border" %in% names(th))
  expect_true("strip.background" %in% names(th))
})

test_that("umap_colors vector has expected length and valid colors", {
  expect_true(length(umap_colors) >= 20)
  valid <- sapply(umap_colors, function(x) {
    tryCatch({
      col2rgb(x)
      TRUE
    }, error = function(e) FALSE)
  })
  expect_true(all(valid))
})

test_that(".example_10x creates directories and calls curl_download", {
  tmpdir <- tempfile("dotools_test_")
  stub(.example_10x, "curl_download", function(url, destfile, mode) {
    cat("Mock download", url, "->", destfile, "\n")
  })

  dir_created <- .example_10x(base_dir = tmpdir)

  expect_true(grepl(tmpdir, dir_created))
  expect_true(dir.exists(file.path(tmpdir, "healthy", "outs")))
  expect_true(dir.exists(file.path(tmpdir, "disease", "outs")))
})
