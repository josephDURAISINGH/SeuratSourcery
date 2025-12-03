context("Check loadRune appends a single dataset properly")
library(SeuratSourcery)
test_that("loadRune appends a single dataset", {
  obj1 <- demo_small()
  obj2 <- demo_small()

  initial <- list(A = obj1)
  out <- loadRune("dummy", data = initial)

  expect_length(out, 2)
  expect_true(all(vapply(out, inherits, logical(1), what = "Seurat")))
})

# [END]
