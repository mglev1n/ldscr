#' Estimate heritability
#'
#' @param munged_sumstats Either a dataframe, or a path to a file containing munged summary statistics. Must contain at least columns named `SNP` (rsid), `A1` (effect allele), `A2` (non-effect allele), `N` (total sample size) and `Z` (Z-score)
#' @param sample_prev (numeric) For binary traits, this should be the prevalence of cases in the current sample, used for conversion from observed heritability to liability-scale heritability. The default is `NA`, which is appropriate for quantitative traits or estimating heritability on the observed scale.
#' @param ancestry (character) One of `AFR`, `AMR`, `CSA`, `EAS`, `EUR`, or `MID`, which will utilize the appropriate built-in `ld` and `wld` files from Pan-UK Biobank. If empty or `NULL`, the user must specify paths to `ld` and `wld` files.
#' @param population_prev (numeric) For binary traits, this should be the population prevalence of the trait, used for conversion from observed heritability to liability-scale heritability. The default is `NA`, which is appropriate for quantitative traits or estimating heritability on the observed scale.
#' @param ld (character) Path to directory containing ld score files, ending in `*.l2.ldscore.gz`. Default is `NA`, which will utilize the built-in ld score files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @param wld (character) Path to directory containing weight files. Default is `NA`, which will utilize the built-in weight files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @param n_blocks (numeric) Number of blocks used to produce block jackknife standard errors. Default is `200`
#' @param chisq_max (numeric) Maximum value of Z^2 for SNPs to be included in LD-score regression. Default is to set `chisq_max` to the maximum of 80 and N*0.001.
#'
#' @return a [tibble][tibble::tibble-package] of heritability information
#'
#' @export
#'

ldsc_h2 <- function(munged_sumstats, ancestry, sample_prev = NA, population_prev = NA, ld, wld, n_blocks = 200, chisq_max = NA) {
  # Check function arguments
  # if (!missing(ancestry)) {
  #   checkmate::assert_choice(ancestry, c("AFR", "AMR", "CSA", "EAS", "EUR", "MID"), null.ok = FALSE)
  #   cli::cli_progress_step("Using {ancestry} reference from Pan-UKB")
  # } else {
  #   cli::cli_progress_step("No ancestry specified, checking for user-specified `ld` and `wld`")
  #   checkmate::assert_character(ld)
  #   checkmate::assert_character(wld)
  # }

  # Dimensions
  n.traits <- 1
  n.V <- 1

  # Storage:
  cov <- matrix(NA, nrow = n.traits, ncol = n.traits)
  V.hold <- matrix(NA, nrow = n_blocks, ncol = n.V)
  N.vec <- matrix(NA, nrow = 1, ncol = n.V)
  Liab.S <- rep(1, n.traits)
  I <- matrix(NA, nrow = n.traits, ncol = n.traits)


  #########  READ LD SCORES:
  cli::cli_progress_step("Reading LD Scores")

  x <- fs::dir_ls(ld, glob = "*.l2.ldscore.gz") %>%
    vroom::vroom(col_types = vroom::cols())

  x$CM <- NULL
  x$MAF <- NULL


  ######### READ weights:
  cli::cli_progress_step("Reading weights")
  w <- fs::dir_ls(wld, glob = "*.l2.ldscore.gz") %>%
    vroom::vroom(col_types = vroom::cols())

  w$CM <- NULL
  w$MAF <- NULL

  colnames(w)[ncol(w)] <- "wLD"

  ### READ M
  cli::cli_progress_step("Reading M")
  m <- fs::dir_ls(ld, glob = "*.l2.M_5_50") %>%
    vroom::vroom(col_types = vroom::cols(), col_names = FALSE, delim = "\t")

  M.tot <- sum(m)
  m <- M.tot

  ### READ ALL CHI2 + MERGE WITH LDSC FILES
  cli::cli_progress_step("Reading summary statistics")
  s <- 0

  if (checkmate::check_character(munged_sumstats)) {
    sumstats_df <- vroom::vroom(munged_sumstats, col_types = vroom::cols())
  } else if (checkmate::check_data_frame(munged_sumstats)) {
    sumstats_df <- munged_sumstats
  }

  cli::cli_progress_step("Merging summary statistics with LD-score files")
  merged <- dplyr::inner_join(sumstats_df, w %>% dplyr::select(SNP, wLD), by = "SNP") %>%
    dplyr::inner_join(x, by = "SNP") %>%
    dplyr::arrange(CHR, BP)

  cli::cli_alert_info(glue::glue("{nrow(merged)}/{nrow(sumstats_df)} SNPs remain after merging with LD-score files"))

  ## REMOVE SNPS with excess chi-square:
  if (is.na(chisq_max)) {
    chisq_max <- max(0.001 * max(merged$N), 80)
  }
  rm <- (merged$Z^2 > chisq_max)
  merged <- merged[!rm, ]

  cli::cli_alert_info(glue::glue("Removed {sum(rm)} SNPs with Chi^2 > {chisq_max}"))

  ## ESTIMATE Heritability
  cli::cli_progress_step("Estimating heritability")

  merged$chi1 <- merged$Z^2

  samp.prev <- sample_prev
  pop.prev <- population_prev

  n.snps <- nrow(merged)

  ## ADD INTERCEPT:
  merged$intercept <- 1
  merged$x.tot <- merged$L2
  merged$x.tot.intercept <- 1

  #### MAKE WEIGHTS:
  tot.agg <- (M.tot * (mean(merged$chi1) - 1)) / mean(merged$L2 * merged$N)
  tot.agg <- max(tot.agg, 0)
  tot.agg <- min(tot.agg, 1)
  merged$ld <- pmax(merged$L2, 1)
  merged$w.ld <- pmax(merged$wLD, 1)
  merged$c <- tot.agg * merged$N / M.tot
  merged$het.w <- 1 / (2 * (1 + (merged$c * merged$ld))^2)
  merged$oc.w <- 1 / merged$w.ld
  merged$w <- merged$het.w * merged$oc.w
  merged$initial.w <- sqrt(merged$w)
  merged$weights <- merged$initial.w / sum(merged$initial.w)

  N.bar <- mean(merged$N)

  ## preweight LD and chi:

  weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) * merged$weights)
  weighted.chi <- as.matrix(merged$chi1 * merged$weights)


  ## Perfrom analysis:

  n.annot <- 1

  select.from <- floor(seq(from = 1, to = n.snps, length.out = (n_blocks + 1)))
  select.to <- c(select.from[2:n_blocks] - 1, n.snps)

  xty.block.values <- matrix(data = NA, nrow = n_blocks, ncol = (n.annot + 1))
  xtx.block.values <- matrix(data = NA, nrow = ((n.annot + 1) * n_blocks), ncol = (n.annot + 1))
  colnames(xty.block.values) <- colnames(xtx.block.values) <- colnames(weighted.LD)
  replace.from <- seq(from = 1, to = nrow(xtx.block.values), by = (n.annot + 1))
  replace.to <- seq(from = (n.annot + 1), to = nrow(xtx.block.values), by = (n.annot + 1))
  for (i in 1:n_blocks) {
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
  delete.values <- matrix(data = NA, nrow = n_blocks, ncol = (n.annot + 1))
  colnames(delete.values) <- colnames(weighted.LD)
  for (i in 1:n_blocks) {
    xty.delete <- xty - xty.block.values[i, ]
    xtx.delete <- xtx - xtx.block.values[delete.from[i]:delete.to[i], ]
    delete.values[i, ] <- solve(xtx.delete, xty.delete)
  }

  tot.delete.values <- delete.values[, 1:n.annot]
  pseudo.values <- matrix(data = NA, nrow = n_blocks, ncol = length(reg))
  colnames(pseudo.values) <- colnames(weighted.LD)
  for (i in 1:n_blocks) {
    pseudo.values[i, ] <- (n_blocks * reg) - ((n_blocks - 1) * delete.values[i, ])
  }

  jackknife.cov <- cov(pseudo.values) / n_blocks
  jackknife.se <- sqrt(diag(jackknife.cov))
  intercept.se <- jackknife.se[length(jackknife.se)]
  coef.cov <- jackknife.cov[1:n.annot, 1:n.annot] / (N.bar^2)

  cat.cov <- coef.cov * (m %*% t(m))
  tot.cov <- sum(cat.cov)
  tot.se <- sqrt(tot.cov)

  V.hold[, s] <- pseudo.values[, 1]
  N.vec[1, s] <- N.bar

  if (is.na(pop.prev) == F & is.na(samp.prev) == F) {
    conversion.factor <- (pop.prev^2 * (1 - pop.prev)^2) / (samp.prev * (1 - samp.prev) * dnorm(qnorm(1 - pop.prev))^2)
    Liab.S <- conversion.factor
  }

  # cov[j, j] <- reg.tot
  # I[j, j] <- intercept

  lambda.gc <- median(merged$chi1) / qchisq(0.5, df = 1)
  mean.Chi <- mean(merged$chi1)
  ratio <- (intercept - 1) / (mean.Chi - 1)
  ratio.se <- intercept.se / (mean.Chi - 1)

  res <- tibble(
    mean_chisq = mean.Chi,
    lambda_gc = lambda.gc,
    intercept = intercept,
    intercept_se = intercept.se,
    ratio = ratio,
    ratio_se = ratio.se,
    h2_observed = reg.tot,
    h2_observed_se = tot.se,
    h2_Z = reg.tot / tot.se
  )

  return(res)
}
