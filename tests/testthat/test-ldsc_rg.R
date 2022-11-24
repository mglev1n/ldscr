test_that("ldsc_rg works", {
  rg_res <- ldsc_rg(
    munged_sumstats = list(
      "EUR" = sumstats_munged_example(),
      "EUR2" = sumstats_munged_example()
    ),
    ancestry = "EUR"
  )

  expect_type(rg_res, "list")
})
