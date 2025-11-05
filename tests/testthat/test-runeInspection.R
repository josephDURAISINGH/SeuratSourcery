test_that("runeInspection returns correct structure", {
  runes <- list(a = demo_rune_data(), b = demo_rune_data())
  report <- runeInspection(runes)

  expect_true("summary" %in% names(report))
  expect_s3_class(report$summary, "tbl_df")
})
