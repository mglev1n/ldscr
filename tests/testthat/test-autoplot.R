test_that("Correlation heatmap autoplot works", {
  rg_res <- ldsc_rg(
    munged_sumstats = list(
      "APOB" = sumstats_munged_example(example = "APOB"),
      "LDL" = sumstats_munged_example(example = "LDL")
    ),
    ancestry = "EUR"
  )

  autoplot.ldscr_list(rg_res)

  expect_error(.plot <- ggplot2::autoplot(rg_res), NA)
  expect_s3_class(.plot, "gg")
})

