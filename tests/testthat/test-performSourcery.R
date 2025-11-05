test_that("performSourcery runs full pipeline", {
  tmp <- tempfile()
  saveRDS(demo_rune_data(), file.path(tmp, "sample1.rds"))
  saveRDS(demo_rune_data(), file.path(tmp, "sample2.rds"))

  integrated <- performSourcery(path = tmp, method = "standard", dims = 1:5, visualize = FALSE)
  expect_s3_class(integrated, "Seurat")
})
