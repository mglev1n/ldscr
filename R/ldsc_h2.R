#' Estimate heritability
#'
#' @description
#'
#' `ldsc_h2()` uses ldscore regression to estimate the heritability of a trait from GWAS summary statistics and reference LD information.
#'
#'
#' @param munged_sumstats Either a dataframe, or a path to a file containing munged summary statistics. Must contain at least columns named `SNP` (rsid), `A1` (effect allele), `A2` (non-effect allele), `N` (total sample size) and `Z` (Z-score)
#' @param ancestry (character) One of `AFR`, `AMR`, `CSA`, `EAS`, `EUR`, or `MID`, which will utilize the appropriate built-in `ld` and `wld` files from Pan-UK Biobank. If empty or `NULL`, the user must specify paths to `ld` and `wld` files.
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

  if (missing(ancestry)) {
    x <- fs::dir_ls(ld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  } else {
    x <- ldscore_files(ancestry, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  }

  x$CM <- x$MAF <- NULL


  # READ weights:
  cli::cli_progress_step("Reading weights")
  if (missing(ancestry)) {
    w <- fs::dir_ls(wld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  } else {
    w <- ldscore_files(ancestry, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  }

  w$CM <- w$MAF <- NULL

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

  M.tot <- sum(m)
  m <- M.tot

  ### READ ALL CHI2 + MERGE WITH LDSC FILES
  s <- 0
  cli::cli_progress_step("Reading summary statistics")
  if (is.character(munged_sumstats)) {
    cli::cli_progress_step("Reading summary statistics from {munged_sumstats}")
    sumstats_df <- vroom::vroom(munged_sumstats, col_types = vroom::cols())
  } else {
    cli::cli_progress_step("Reading summary statistics from dataframe")
    sumstats_df <- munged_sumstats
  }

  cli::cli_progress_step("Merging summary statistics with LD-score files")
  merged <- sumstats_df %>%
    dtplyr::lazy_dt() %>%
    dplyr::select(SNP, N, Z, A1) %>%
    dplyr::inner_join(w[, c("SNP", "wLD")], by = c("SNP")) %>%
    dplyr::inner_join(x, by = c("SNP")) %>%
    dplyr::arrange(CHR, BP) %>%
    na.omit() %>%
    unique() %>%
    tibble::as_tibble()

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


#' Paths to LD-score files from Pan-UKB
#'
#' @param ancestry One of `AFR`, `AMR`, `CSA`, `EAS`, `EUR`, or `MID`
#' @param ... arguments passed to `fs::dir_ls()`
#'
#' @return a list of paths
#'
#' @export
#' @examples
#' ldscore_files("AFR")
#'
ldscore_files <- function(ancestry, ...) {
  fs::dir_ls(fs::path(fs::path_package("extdata", package = "ldscr"), ancestry), ...)
}

#' Example munged dataframe
#'
#' @param dataframe (logical) If `TRUE` (default), return an example munged dataframe. If `FALSE`, return path to the file on disk.
#' @param example (character) One of `BMI` or `LDL` which have been included as example traits.
#' @return either a [tibble][tibble::tibble-package] containing a munged dataframe, or a path to the file on disk.
#'
#' @export
#' @examples
#' sumstats_munged_example(example = "BMI", dataframe = TRUE)
sumstats_munged_example <- function(example, dataframe = TRUE) {
  if (dataframe) {
    vroom::vroom(fs::path(fs::path_package("extdata", paste0(example, "-sumstats-munged.txt.gz"), package = "ldscr")), col_types = vroom::cols())
  } else {
    fs::path(fs::path_package("extdata", paste0(example, "-sumstats-munged.txt.gz"), package = "ldscr"))
  }
}

# Internal Function to make weights:
make_weights <- function(chi1, L2, wLD, N, M.tot) {
  tot.agg <- (M.tot * (mean(chi1) - 1)) / mean(L2 * N)
  tot.agg <- max(tot.agg, 0)
  tot.agg <- min(tot.agg, 1)
  ld <- pmax(L2, 1)
  w.ld <- pmax(wLD, 1)
  c <- tot.agg * N / M.tot
  het.w <- 1 / (2 * (1 + (c * ld))^2)
  oc.w <- 1 / w.ld
  w <- het.w * oc.w
  initial.w <- sqrt(w)

  return(initial.w)
}

# Internal function to perform LDSC heritability/covariance analysis
perform_analysis <- function(n.blocks, n.snps, weighted.LD, weighted.chi, N.bar, m) {
  n.annot <- 1

  select.from <- floor(seq(from = 1, to = n.snps, length.out = (n.blocks + 1)))
  select.to <- c(select.from[2:n.blocks] - 1, n.snps)

  xty.block.values <- matrix(data = NA, nrow = n.blocks, ncol = (n.annot + 1))
  xtx.block.values <- matrix(data = NA, nrow = ((n.annot + 1) * n.blocks), ncol = (n.annot + 1))
  colnames(xty.block.values) <- colnames(xtx.block.values) <- colnames(weighted.LD)
  replace.from <- seq(from = 1, to = nrow(xtx.block.values), by = (n.annot + 1))
  replace.to <- seq(from = (n.annot + 1), to = nrow(xtx.block.values), by = (n.annot + 1))
  for (i in 1:n.blocks) {
    xty.block.values[i, ] <- t(t(weighted.LD[select.from[i]:select.to[i], ]) %*% weighted.chi[select.from[i]:select.to[i], ])
    xtx.block.values[replace.from[i]:replace.to[i], ] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i], ]) %*% weighted.LD[select.from[i]:select.to[i], ])
  }
  xty <- as.matrix(colSums(xty.block.values))
  xtx <- matrix(data = NA, nrow = (n.annot + 1), ncol = (n.annot + 1))
  colnames(xtx) <- colnames(weighted.LD)
  for (i in 1:nrow(xtx)) {
    xtx[i, ] <- t(colSums(xtx.block.values[seq(from = i, to = nrow(xtx.block.values), by = ncol(weighted.LD)), ]))
  }

  reg <- solve(xtx, xty)
  intercept <- reg[2]
  coefs <- reg[1] / N.bar
  reg.tot <- coefs * m

  delete.from <- seq(from = 1, to = nrow(xtx.block.values), by = ncol(xtx.block.values))
  delete.to <- seq(from = ncol(xtx.block.values), to = nrow(xtx.block.values), by = ncol(xtx.block.values))
  delete.values <- matrix(data = NA, nrow = n.blocks, ncol = (n.annot + 1))
  colnames(delete.values) <- colnames(weighted.LD)
  for (i in 1:n.blocks) {
    xty.delete <- xty - xty.block.values[i, ]
    xtx.delete <- xtx - xtx.block.values[delete.from[i]:delete.to[i], ]
    delete.values[i, ] <- solve(xtx.delete, xty.delete)
  }

  tot.delete.values <- delete.values[, 1:n.annot]
  pseudo.values <- matrix(data = NA, nrow = n.blocks, ncol = length(reg))
  colnames(pseudo.values) <- colnames(weighted.LD)
  for (i in 1:n.blocks) {
    pseudo.values[i, ] <- (n.blocks * reg) - ((n.blocks - 1) * delete.values[i, ])
  }

  jackknife.cov <- cov(pseudo.values) / n.blocks
  jackknife.se <- sqrt(diag(jackknife.cov))
  intercept.se <- jackknife.se[length(jackknife.se)]
  coef.cov <- jackknife.cov[1:n.annot, 1:n.annot] / (N.bar^2)
  cat.cov <- coef.cov * (m %*% t(m))
  tot.cov <- sum(cat.cov)
  tot.se <- sqrt(tot.cov)

  return(
    list(
      reg.tot = reg.tot,
      tot.se = tot.se,
      intercept = intercept,
      intercept.se = intercept.se,
      pseudo.values = pseudo.values[, 1],
      N.bar = N.bar
    )
  )
}
