#' Estimate heritability
#'
#' @description
#'
#' `ldsc_h2()` uses ldscore regression to estimate the heritability of a trait from GWAS summary statistics and reference LD information.
#'
#'
#' @param munged_sumstats Either a dataframe, or a path to a file containing munged summary statistics. Must contain at least columns named `SNP` (rsid), `A1` (effect allele), `A2` (non-effect allele), `N` (total sample size) and `Z` (Z-score)
#' @param ancestry (character) One of "AFR", "AMR", "CSA", "EAS", "EUR", or "MID", which will utilize the appropriate built-in `ld` and `wld` files from Pan-UK Biobank. If empty or `NULL`, the user must specify paths to `ld` and `wld` files.
#' @param sample_prev (numeric) For binary traits, this should be the prevalence of cases in the current sample, used for conversion from observed heritability to liability-scale heritability. The default is `NA`, which is appropriate for quantitative traits or estimating heritability on the observed scale.
#' @param population_prev (numeric) For binary traits, this should be the population prevalence of the trait, used for conversion from observed heritability to liability-scale heritability. The default is `NA`, which is appropriate for quantitative traits or estimating heritability on the observed scale.
#' @param ld (character) Path to directory containing ld score files, ending in `*.l2.ldscore.gz`. Default is `NA`, which will utilize the built-in ld score files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @param wld (character) Path to directory containing weight files. Default is `NA`, which will utilize the built-in weight files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @param n_blocks (numeric) Number of blocks used to produce block jackknife standard errors. Default is `200`
#' @param chisq_max (numeric) Maximum value of Z^2 for SNPs to be included in LD-score regression. Default is to set `chisq_max` to the maximum of 80 and N*0.001.
#'
#' @return A [tibble][tibble::tibble-package] containing heritability information. If `sample_prev` and `population_prev` were provided, the heritability estimate will also be returned on the liability scale.
#' @export
#'
#' @import dtplyr
#' @import data.table
#'
#' @examples
#' \donttest{
#' ldsc_h2(sumstats_munged_example(example = "BMI", dataframe = TRUE), ancestry = "EUR")
#' }
#'
ldsc_h2 <- function(munged_sumstats, ancestry, sample_prev = NA, population_prev = NA, ld, wld, n_blocks = 200, chisq_max = NA) {
  # Check function arguments
  if (missing(ancestry)) {
    cli::cli_progress_step("No ancestry specified, checking for user-specified `ld` and `wld`")
    checkmate::assert_directory_exists(ld)
    checkmate::assert_directory_exists(wld)
  } else {
    checkmate::assert_choice(ancestry, c("AFR", "AMR", "CSA", "EAS", "EUR", "MID"), null.ok = FALSE)
    cli::cli_progress_step("Using {ancestry} reference from Pan-UKB")
  }

  # Dimensions
  n.traits <- 1
  n.V <- 1

  # Storage:
  cov <- matrix(NA, nrow = n.traits, ncol = n.traits)
  V.hold <- matrix(NA, nrow = n_blocks, ncol = n.V)
  N.vec <- matrix(NA, nrow = 1, ncol = n.V)
  Liab.S <- rep(1, n.traits)
  I <- matrix(NA, nrow = n.traits, ncol = n.traits)


  # READ LD SCORES:
  cli::cli_progress_step("Reading LD Scores")
  x <- read_ld(ancestry, ld)
  x$CM <- x$MAF <- NULL


  # READ weights:
  cli::cli_progress_step("Reading weights")
  w <- read_wld(ancestry, wld)
  w$CM <- w$MAF <- NULL
  colnames(w)[ncol(w)] <- "wLD"

  # READ M
  cli::cli_progress_step("Reading M")
  m <- read_m(ancestry, ld)
  M.tot <- sum(m)
  m <- M.tot

  ### READ ALL CHI2 + MERGE WITH LDSC FILES
  s <- 0
  cli::cli_progress_step("Reading summary statistics")
  sumstats_df <- read_sumstats(munged_sumstats)

  cli::cli_progress_step("Merging summary statistics with LD-score files")
  merged <- merge_sumstats(sumstats_df, w, x)

  cli::cli_alert_info(glue::glue("{nrow(merged)}/{nrow(sumstats_df)} SNPs remain after merging with LD-score files"))

  ## REMOVE SNPS with excess chi-square:
  if (is.na(chisq_max)) {
    chisq_max <- max(0.001 * max(merged$N), 80)
  }
  rm <- (merged$Z^2 > chisq_max)
  merged <- merged[!rm, ]

  cli::cli_alert_info(glue::glue("Removed {sum(rm)} SNPs with Chi^2 > {chisq_max}; {nrow(merged)} SNPs remain"))

  ## ESTIMATE Heritability
  cli::cli_progress_step("Estimating heritability")

  merged$chi1 <- merged$Z^2

  n.snps <- nrow(merged)

  ## ADD INTERCEPT:
  merged$intercept <- 1
  merged$x.tot <- merged$L2
  merged$x.tot.intercept <- 1

  ## MAKE WEIGHTS:
  initial.w <- make_weights(chi1 = merged$chi1, L2 = merged$L2, wLD = merged$wLD, N = merged$N, M.tot)

  # return(list(merged, initial.w))

  merged$weights <- initial.w / sum(initial.w)

  N.bar <- mean(merged$N)

  ## Preweight LD and chi:

  weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) * merged$weights)
  weighted.chi <- as.matrix(merged$chi1 * merged$weights)


  ## Perform analysis:
  analysis_res <- perform_analysis(n.blocks = n_blocks, n.snps, weighted.LD, weighted.chi, N.bar, m)

  lambda.gc <- median(merged$chi1) / qchisq(0.5, df = 1)
  mean.Chi <- mean(merged$chi1)
  ratio <- (analysis_res$intercept - 1) / (mean.Chi - 1)
  ratio.se <- analysis_res$intercept.se / (mean.Chi - 1)

  if (is.na(population_prev) == F & is.na(sample_prev) == F) {
    # conversion.factor <- (population_prev^2 * (1 - population_prev)^2) / (sample_prev * (1 - sample_prev) * dnorm(qnorm(1 - population_prev))^2)
    # Liab.S <- conversion.factor
    h2_lia <- h2_liability(h2 = analysis_res$reg.tot, sample_prev, population_prev)

    h2_res <- tibble(
      mean_chisq = mean.Chi,
      lambda_gc = lambda.gc,
      intercept = analysis_res$intercept,
      intercept_se = analysis_res$intercept.se,
      ratio = ratio,
      ratio_se = ratio.se,
      h2_observed = analysis_res$reg.tot,
      h2_observed_se = analysis_res$tot.se,
      h2_Z = analysis_res$reg.tot / analysis_res$tot.se,
      h2_p = 2 * pnorm(abs(h2_Z), lower.tail = FALSE),
      h2_liability = h2_lia,
      h2_liability_se = h2_lia / h2_Z
    )
  } else {
    h2_res <- tibble(
      mean_chisq = mean.Chi,
      lambda_gc = lambda.gc,
      intercept = analysis_res$intercept,
      intercept_se = analysis_res$intercept.se,
      ratio = ratio,
      ratio_se = ratio.se,
      h2_observed = analysis_res$reg.tot,
      h2_observed_se = analysis_res$tot.se,
      h2_Z = analysis_res$reg.tot / analysis_res$tot.se,
      h2_p = 2 * pnorm(abs(h2_Z), lower.tail = FALSE)
    )
  }

  return(h2_res)
}


#' Convert Heritability to Liability Scale
#'
#' @description
#'
#' `h2_liability()` converts heritability estimates from the observed to liability scale.
#'
#'
#' @param h2 (numeric) Estimate of observed-scale heritability
#' @param sample_prev (numeric) Proportion of cases in the current sample
#' @param population_prev (numeric) Population prevalence of trait
#' @return A list containing liability-scale heritability, standard error, and p-value
#' @export

#' @examples
#' h2_liability(0.28, 0.1, 0.05)
#'
h2_liability <- function(h2, sample_prev, population_prev) {
  checkmate::assert_double(h2, lower = 0, upper = 1)
  checkmate::assert_double(sample_prev, lower = 0, upper = 1)
  checkmate::assert_double(population_prev, lower = 0, upper = 1)

  # From equation 23 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059431/ Estimating Missing Heritability for Disease from Genome-wide Association Studies
  K <- population_prev
  P <- sample_prev
  zv <- dnorm(qnorm(K))

  h2_liab <- h2 * K^2 * (1 - K)^2 / P / (1 - P) / zv^2

  return(h2_liab)
}



