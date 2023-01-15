test_that("read_ld works on internal data", {
  x <- read_ld(ancestry = "EUR")
  expect_s3_class(x, "data.frame")
})

test_that("read_wld works on internal data", {
  w <- read_wld(ancestry = "EUR")
  expect_s3_class(w, "data.frame")
})

test_that("read_m works on internal data", {
  m <- read_m(ancestry = "EUR")
  expect_s3_class(m, "data.frame")
})

test_that("read_sumstats works on internal data", {
  sumstats_df <- read_sumstats(sumstats_munged_example(example = "BMI", dataframe = TRUE))
  expect_s3_class(sumstats_df, "data.frame")
})

test_that("merge_sumstats works on internal data", {
  x <- read_ld(ancestry = "EUR")
  w <- x %>% dplyr::mutate(wLD = L2)
  sumstats_df <- read_sumstats(sumstats_munged_example(example = "BMI", dataframe = TRUE))
  sumstats_df <- merge_sumstats(sumstats_df, w, x, chr_filter = seq(2, 22, 2))
  expect_s3_class(sumstats_df, "data.frame")
})
