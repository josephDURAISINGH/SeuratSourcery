context("Check loadRune appends a single dataset properly")
library(SeuratSourcery)
test_that("loadRune appends a single dataset", {
  td <- tempdir()
  obj1 <- demo_rune_data()
  obj2 <- demo_rune_data()

  f1 <- file.path(td, "A.rds")
  f2 <- file.path(td, "B.rds")
  saveRDS(obj1, f1)
  saveRDS(obj2, f2)

  out <- loadRune(f1)
  out <- loadRune(f2, data = out)

  # Tests
  expect_length(out, 2)
  expect_true(all(vapply(out, inherits, logical(1), what = "Seurat")))
})

# [END]
