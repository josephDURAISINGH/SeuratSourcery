test_that("basicIntegration integrates without error", {
  runes <- list(a = demo_rune_data(), b = demo_rune_data())
  integrated <- basicIntegration(runes, dims = 1:5)

  expect_s3_class(integrated, "Seurat")
  expect_true("integrated" %in% Seurat::Layers(integrated))
})
