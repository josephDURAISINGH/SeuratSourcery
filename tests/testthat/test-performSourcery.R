context("make sure pipeline works")
library(SeuratSourcery)

test_that("performSourcery runs full pipeline", {
  runes <- list(A = demo_small(), B = demo_small())

  out <- performSourcery(datasets = runes, visualize = FALSE, verbose = FALSE)

  expect_type(out, "list")
  expect_length(out, 2)
  expect_true(all(vapply(out, inherits, logical(1), what = "Seurat")))
})

# [END]
