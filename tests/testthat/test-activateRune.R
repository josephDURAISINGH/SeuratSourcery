test_that("activateRune harmonizes and adds layers", {
  runes <- list(sample1 = demo_rune_data())
  harmonized <- activateRune(runes)

  expect_true("harmonized" %in% names(harmonized[[1]]@assays$RNA@layers))
  expect_true(all(harmonized[[1]]@misc$harmonized))
})
