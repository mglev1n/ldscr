test_that("ldsc_rg works with two traits", {
  rg_res <- ldsc_rg(
    munged_sumstats = list(
      "APOB" = sumstats_munged_example(example = "APOB"),
      "LDL" = sumstats_munged_example(example = "LDL")
    ),
    ancestry = "EUR"
  )

  expect_s3_class(rg_res, "ldscr_list")
  expect_s3_class(rg_res$rg, "data.frame")
})

test_that("ldsc_rg works with three traits", {
  rg_res <- ldsc_rg(
    munged_sumstats = list(
      "BMI" = sumstats_munged_example(example = "BMI"),
      "LDL" = sumstats_munged_example(example = "LDL", dataframe = TRUE),
      "BMI2" = sumstats_munged_example(example = "BMI", dataframe = TRUE)
    ),
    ancestry = "EUR"
  )

  expect_type(rg_res, "list")
  expect_equal(nrow(rg_res$rg), 3)
})
