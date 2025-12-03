context("Check activateRune actually modifies the seurat objects during harmonization")
library(SeuratSourcery)

test_that("activateRune harmonizes basic Seurat objects", {
  runes <- list(A = demo_rune_data(), B = demo_rune_data())
  out <- activateRune(runes)

  expect_type(out, "list")
  expect_length(out, 2)
  expect_true(all(vapply(out, inherits, logical(1), what = "Seurat")))

  # Each object should now have misc$harmonized = TRUE
  expect_true(all(vapply(out, function(x) isTRUE(x@misc$harmonized), logical(1))))
})

# [END]
