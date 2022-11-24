test_that("heritability estimation using built-in data works", {
  ldsc_res <- ldsc_h2("extdata/BMI-sumstats-munged.txt.gz", ancestry = "EUR")
  expect_equal(dim(ldsc_res), c(1, 9))
})
