test_that("Heritability estimation using supplied ld data works", {
  ldsc_res <- ldsc_h2(sumstats_munged_example(dataframe = FALSE),
                      ld = unique(fs::path_dir(ldscore_files("EUR"))),
                      wld = unique(fs::path_dir(ldscore_files("EUR"))))
  expect_equal(dim(ldsc_res), c(1, 9))
})

test_that("Heritability estimation using built-in ld data works", {
  ldsc_res <- ldsc_h2(sumstats_munged_example(dataframe = FALSE),
                      ancestry = "EUR")
  expect_equal(dim(ldsc_res), c(1, 9))
})

test_that("Heritability estimation using munged dataframe works", {
  ldsc_res <- ldsc_h2(sumstats_munged_example(dataframe = TRUE), ancestry = "EUR")
  expect_equal(dim(ldsc_res), c(1, 9))
})

test_that("Conversion from observed to liability heritability works", {
  expect_type(h2_liability(0.28, 0.1, 0.05), "double")
})

