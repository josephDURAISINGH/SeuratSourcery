context("Valud metric outputs")
library(SeuratSourcery)

test_that("runeInspection computes metrics correctly", {
  runes <- list(A = demo_rune_data(), B = demo_rune_data())
  out <- runeInspection(runes)

  expect_true(is.list(out))
  expect_true("summary" %in% names(out))
  expect_true("overlap" %in% names(out))
  expect_true("common_meta" %in% names(out))
  expect_true("datasets" %in% names(out))

  expect_equal(nrow(out$summary), 2)
  expect_equal(ncol(out$overlap), 2)
})

# [END]
