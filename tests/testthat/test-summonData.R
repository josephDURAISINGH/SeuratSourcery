test_that("summonData handles basic input", {
  tmpdir <- tempdir()
  tmpfile <- file.path(tmpdir, "sample1.rds")
  saveRDS(demo_rune_data(), tmpfile)

  datasets <- summonData(path = tmpdir, pattern = "\\.rds$")
  expect_type(datasets, "list")
  expect_s3_class(datasets[[1]], "Seurat")
})
