#' Estimate cross-trait genetic correlations
#'
#' @description
#'
#' `ldsc_rg()` uses ldscore regression to estimate the pairwise genetic correlations between traits. The function relies on named lists of traits, sample prevalences, and populatin prevalences. The name of each trait should be consistent across each argument.
#'
#'
#' @param munged_sumstats (list) A named list of dataframes, or paths to files containing munged summary statistics. Each set of munged summary statistics contain at least columns named `SNP` (rsid), `A1` (effect allele), `A2` (non-effect allele), `N` (total sample size) and `Z` (Z-score)
#' @param ancestry (character) One of `AFR`, `AMR`, `CSA`, `EAS`, `EUR`, or `MID`, which will utilize the appropriate built-in `ld` and `wld` files from Pan-UK Biobank. If empty or `NULL`, the user must specify paths to `ld` and `wld` files.
#' @param sample_prev (list) A named list containing the prevalence of cases in the current sample, used for conversion from observed heritability to liability-scale heritability. The default is `NA`, which is appropriate for quantitative traits or estimating heritability on the observed scale.
#' @param population_prev (list) A named list containing the population prevalence of the trait, used for conversion from observed heritability to liability-scale heritability. The default is `NA`, which is appropriate for quantitative traits or estimating heritability on the observed scale.
#' @param ld (character) Path to directory containing ld score files, ending in `*.l2.ldscore.gz`. Default is `NA`, which will utilize the built-in ld score files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @param wld (character) Path to directory containing weight files. Default is `NA`, which will utilize the built-in weight files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @param n_blocks (numeric) Number of blocks used to produce block jackknife standard errors. Default is `200`
#' @param chisq_max (numeric) Maximum value of Z^2 for SNPs to be included in LD-score regression. Default is to set `chisq_max` to the maximum of 80 and N*0.001.
#'
#' @return A list of heritablilty and genetic correlation information
#'  - `rg` = [tibble][tibble::tibble-package] containing pairwise genetic correlations information.
#'  - `h2` = [tibble][tibble::tibble-package] containing heritability information for each trait. If `sample_prev` and `population_prev` were provided, the heritability estimates will also be returned on the liability scale.
#' @export
#'

ldsc_rg <- function(munged_sumstats, ancestry, sample_prev = NA, population_prev = NA, ld, wld, n_blocks = 200, chisq_max = NA) {
  # Check function arguments
  if (missing(ancestry)) {
    cli::cli_progress_step("No ancestry specified, checking for user-specified `ld` and `wld`")
    checkmate::assert_file_exists(ld)
    checkmate::assert_file_exists(wld)
  } else {
    checkmate::assert_choice(ancestry, c("AFR", "AMR", "CSA", "EAS", "EUR", "MID"), null.ok = FALSE)
    cli::cli_progress_step("Using {ancestry} reference from Pan-UKB")
  }

  checkmate::assert_list(munged_sumstats)

  if (missing(sample_prev)) {
    cli::cli_alert_info("No sample prevalence data provided. Estimating heritabilities on the observed scale.")
  }

  # Dimensions
  n.traits <- length(munged_sumstats)
  n.V <- n.traits * (n.traits + 1) / 2

  # Set number of blocks
  if (n.traits > 18) {
    n.blocks <- (((n.traits + 1) * (n.traits + 2)) / 2) + 1
    cli::cli_alert_info("Setting the number of blocks used to perform the block jacknife used to estimate the sampling covariance matrix (V) to {n.blocks}")
    if (n.blocks > 1000) {
      cli_alert_warning("The number of blocks needed to estimate V is > 1000, which may result in sampling dependencies across the blocks used to estimate standard errors and can bias results.")
    }
  }

  # Storage:
  cov <- matrix(NA, nrow = n.traits, ncol = n.traits)
  V.hold <- matrix(NA, nrow = n.blocks, ncol = n.V)
  N.vec <- matrix(NA, nrow = 1, ncol = n.V)
  Liab.S <- rep(1, n.traits)
  I <- matrix(NA, nrow = n.traits, ncol = n.traits)

  # READ LD SCORES:
  cli::cli_progress_step("Reading LD Scores")

  if (missing(ancestry)) {
    x <- fs::dir_ls(ld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  } else {
    x <- ldscore_files(ancestry, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  }

  x$CM <- NULL
  x$MAF <- NULL


  # READ weights:
  cli::cli_progress_step("Reading weights")
  if (missing(ancestry)) {
    w <- fs::dir_ls(wld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  } else {
    w <- ldscore_files(ancestry, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  }

  w$CM <- NULL
  w$MAF <- NULL

  colnames(w)[ncol(w)] <- "wLD"

  # READ M
  cli::cli_progress_step("Reading M")
  if (missing(ancestry)) {
    m <- fs::dir_ls(ld, glob = "*.l2.M_5_50") %>%
      vroom::vroom(col_types = vroom::cols(), col_names = FALSE, delim = "\t")
  } else {
    m <- ldscore_files(ancestry, glob = "*.l2.M_5_50") %>%
      vroom::vroom(col_types = vroom::cols(), col_names = FALSE, delim = "\t")
  }
}
