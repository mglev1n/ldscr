test_that("Heritability estimation using supplied ld data works", {
  ldsc_res <- ldsc_h2(sumstats_munged_example(example = "BMI", dataframe = FALSE),
    ld = unique(fs::path_dir(ldscore_files("EUR"))),
    wld = unique(fs::path_dir(ldscore_files("EUR")))
  )
  expect_equal(dim(ldsc_res), c(1, 10))
})

test_that("Heritability estimation using built-in ld data works", {
  ldsc_res <- ldsc_h2(sumstats_munged_example(example = "BMI", dataframe = FALSE),
    ancestry = "EUR"
  )
  expect_equal(dim(ldsc_res), c(1, 10))
})

test_that("Heritability estimation using munged dataframe works", {
  ldsc_res <- ldsc_h2(sumstats_munged_example(example = "BMI", dataframe = TRUE), ancestry = "EUR")
  expect_equal(dim(ldsc_res), c(1, 10))
})

test_that("Heritability estimation on liability scale works", {
  ldsc_res <- ldsc_h2(sumstats_munged_example(example = "BMI", dataframe = TRUE), ancestry = "EUR", sample_prev = 0.2, population_prev = 0.01)
  expect_equal(dim(ldsc_res), c(1, 12))
})

test_that("Conversion from observed to liability heritability works", {
  expect_type(h2_liability(0.28, 0.1, 0.05), "double")
})
