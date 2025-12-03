context("Testing data loading")
library(SeuratSourcery)

test_that("summonData loads multiple RDS Seurat objects", {
  td <- tempdir()
  obj1 <- demo_rune_data()
  obj2 <- demo_rune_data()

  saveRDS(obj1, file.path(td, "A.rds"))
  saveRDS(obj2, file.path(td, "B.rds"))

  out <- summonData(td)

  expect_type(out, "list")
  expect_length(out, 2)
  expect_true(all(vapply(out, inherits, logical(1), what = "Seurat")))
})

# [END]
