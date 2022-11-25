#' Estimate cross-trait genetic correlations
#'
#' @description
#'
#' `ldsc_rg()` uses ldscore regression to estimate the pairwise genetic correlations between traits. The function relies on named lists of traits, sample prevalences, and populatin prevalences. The name of each trait should be consistent across each argument.
#'
#'
#' @param munged_sumstats (list) A named list of dataframes, or paths to files containing munged summary statistics. Each set of munged summary statistics contain at least columns named `SNP` (rsid), `A1` (effect allele), `A2` (non-effect allele), `N` (total sample size) and `Z` (Z-score)
#' @param ancestry (character) One of "AFR", "AMR", "CSA", "EAS", "EUR", or "MID", which will utilize the appropriate built-in `ld` and `wld` files from Pan-UK Biobank. If empty or `NULL`, the user must specify paths to `ld` and `wld` files.
#' @param sample_prev (list) A named list containing the prevalence of cases in the current sample, used for conversion from observed heritability to liability-scale heritability. The default is `NA`, which is appropriate for quantitative traits or estimating heritability on the observed scale.
#' @param population_prev (list) A named list containing the population prevalence of the trait, used for conversion from observed heritability to liability-scale heritability. The default is `NA`, which is appropriate for quantitative traits or estimating heritability on the observed scale.
#' @param ld (character) Path to directory containing ld score files, ending in `*.l2.ldscore.gz`. Default is `NA`, which will utilize the built-in ld score files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @param wld (character) Path to directory containing weight files. Default is `NA`, which will utilize the built-in weight files from Pan-UK Biobank for the ancestry specified in `ancestry`.
#' @param n_blocks (numeric) Number of blocks used to produce block jackknife standard errors. Default is `200`
#' @param chisq_max (numeric) Maximum value of Z^2 for SNPs to be included in LD-score regression. Default is to set `chisq_max` to the maximum of 80 and N*0.001.
#'
#' @return A list of heritablilty and genetic correlation information
#'  - `h2` = [tibble][tibble::tibble-package] containing heritability information for each trait. If `sample_prev` and `population_prev` were provided, the heritability estimates will also be returned on the liability scale.
#'  - `rg` = [tibble][tibble::tibble-package] containing pairwise genetic correlations information.
#'  - `raw` = A list of correlation/covariance matrices
#'
#' @export
#'
#' @import dtplyr
#' @import data.table
#'

ldsc_rg <- function(munged_sumstats, ancestry, sample_prev = NA, population_prev = NA, ld, wld, n_blocks = 200, chisq_max = NA) {
  # Check function arguments
  if (missing(ancestry)) {
    cli::cli_progress_step("No ancestry specified, checking for user-specified `ld` and `wld`")
    checkmate::assert_directory_exists(ld)
    checkmate::assert_directory_exists(wld)
  } else {
    checkmate::assert_choice(ancestry, c("AFR", "AMR", "CSA", "EAS", "EUR", "MID"), null.ok = FALSE)
    cli::cli_progress_step("Using {ancestry} reference from Pan-UKB")
  }

  checkmate::assert_list(munged_sumstats)

  if (missing(sample_prev)) {
    cli::cli_alert_info("No sample prevalence data provided. Estimating heritabilities on the observed scale.")
  }

  checkmate::assert_number(n_blocks)
  n.blocks <- n_blocks

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
  h2_res <- tibble()

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

  M.tot <- sum(m)
  m <- M.tot

  cli::cli_progress_step("Reading summary statistics")
  # READ summary statistics

  all_y <- purrr::imap(munged_sumstats, ~ {
    if (is.character(.x)) {
      cli::cli_progress_step("Reading summary statistics for '{.y}' from {.x}")
      sumstats_df <- vroom::vroom(.x, col_types = vroom::cols())
    } else {
      cli::cli_progress_step("Reading summary statistics for '{.y}' from dataframe")
      sumstats_df <- .x
    }
    sumstats_df <- na.omit(sumstats_df)

    cli::cli_progress_step("Merging '{.y}' with LD-score files")
    merged <- sumstats_df %>%
      dtplyr::lazy_dt() %>%
      dplyr::select(SNP, N, Z, A1) %>%
      dplyr::inner_join(w[, c("SNP", "wLD")], by = c("SNP")) %>%
      dplyr::inner_join(x, by = c("SNP")) %>%
      dplyr::arrange(CHR, BP) %>%
      na.omit() %>%
      unique() %>%
      tibble::as_tibble()

    cli::cli_alert_info(glue::glue("{nrow(merged)}/{nrow(sumstats_df)} SNPs remain after merging '{.y}' with LD-score files"))

    ## REMOVE SNPS with excess chi-square:
    if (is.na(chisq_max)) {
      chisq_max <- max(0.001 * max(merged$N), 80)
    }
    rm <- (merged$Z^2 > chisq_max)
    merged <- merged[!rm, ]

    cli::cli_alert_info(glue::glue("Removed {sum(rm)} SNPs with Chi^2 > {chisq_max} from '{.y}'; {nrow(merged)} SNPs remain"))

    return(merged)
  })

  # count the total nummer of runs, both loops
  s <- 1

  for (j in 1:n.traits) {
    # chi1 <- traits[j]

    y1 <- all_y[[j]]
    y1$chi1 <- y1$Z^2

    for (k in j:n.traits) {
      ##### HERITABILITY code
      if (j == k) {
        trait <- names(munged_sumstats[j])
        cli::cli_progress_step("Estimating heritability for '{trait}'")

        samp.prev <- sample_prev[j]
        pop.prev <- population_prev[j]

        merged <- y1
        n.snps <- nrow(merged)

        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1

        initial.w <- make_weights(chi1 = merged$chi1, L2 = merged$L2, wLD = merged$wLD, N = merged$N, M.tot)

        merged$weights <- initial.w / sum(initial.w)

        N.bar <- mean(merged$N)

        ## preweight LD and chi:
        weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) * merged$weights)
        weighted.chi <- as.matrix(merged$chi1 * merged$weights)

        ## Perform analysis:
        analysis_res <- perform_analysis(n.blocks, n.snps, weighted.LD, weighted.chi, N.bar, m)

        V.hold[, s] <- analysis_res$pseudo.values
        N.vec[1, s] <- analysis_res$N.bar

        lambda.gc <- median(merged$chi1) / qchisq(0.5, df = 1)
        mean.Chi <- mean(merged$chi1)
        ratio <- (analysis_res$intercept - 1) / (mean.Chi - 1)
        ratio.se <- analysis_res$intercept.se / (mean.Chi - 1)

        if (is.na(population_prev) == F & is.na(sample_prev) == F) {
          # conversion.factor <- (population_prev^2 * (1 - population_prev)^2) / (sample_prev * (1 - sample_prev) * dnorm(qnorm(1 - population_prev))^2)
          # Liab.S <- conversion.factor
          h2_lia <- h2_liability(h2 = analysis_res$reg.tot, sample_prev, population_prev)

          h2_res <- h2_res %>%
            bind_rows(
              tibble(
                trait = trait,
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
            )
        } else {
          h2_res <- h2_res %>%
            bind_rows(
              tibble(
                trait = trait,
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
            )
        }

        cov[j, j] <- analysis_res$reg.tot
        I[j, j] <- analysis_res$intercept
      }


      ##### GENETIC COVARIANCE code

      if (j != k) {
        trait1 <- names(munged_sumstats[j])
        trait2 <- names(munged_sumstats[k])
        cli::cli_progress_step("Estimating genetic covariance for for '{trait1}' and '{trait2}'")

        # Reuse the data read in for heritability
        y2 <- all_y[[k]]
        y <- merge(y1, y2[, c("SNP", "N", "Z", "A1")], by = "SNP", sort = FALSE)

        y$Z.x <- ifelse(y$A1.y == y$A1.x, y$Z.x, -y$Z.x)
        y$ZZ <- y$Z.y * y$Z.x
        y$chi2 <- y$Z.y^2
        merged <- na.omit(y)
        n.snps <- nrow(merged)

        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1


        #### MAKE WEIGHTS:
        initial.w <- make_weights(chi1 = merged$chi1, L2 = merged$L2, wLD = merged$wLD, N = merged$N.x, M.tot)
        initial.w2 <- make_weights(chi1 = merged$chi2, L2 = merged$L2, wLD = merged$wLD, N = merged$N.y, M.tot)

        merged$weights_cov <- (initial.w + initial.w2) / sum(initial.w + initial.w2)

        N.bar <- sqrt(mean(merged$N.x) * mean(merged$N.y))

        ## preweight LD and chi:

        weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) * merged$weights)
        weighted.chi <- as.matrix(merged$ZZ * merged$weights_cov)

        ## Perform analysis:
        covariance_res <- perform_analysis(n.blocks, n.snps, weighted.LD, weighted.chi, N.bar, m)

        V.hold[, s] <- covariance_res$pseudo.values
        N.vec[1, s] <- covariance_res$N.bar

        cov[k, j] <- cov[j, k] <- covariance_res$reg.tot
        I[k, j] <- I[j, k] <- covariance_res$intercept

      }

      ### Total count
      s <- s + 1
    }
  }


  ## Scale V to N per study (assume m constant)
  # /!\ crossprod instead of tcrossprod because N.vec is a one-row matrix
  v.out <- cov(V.hold) / crossprod(N.vec * (sqrt(n.blocks) / m))

  ### Scale S and V to liability:
  ratio <- tcrossprod(sqrt(Liab.S))
  S <- cov * ratio

  # calculate the ratio of the rescaled and original S matrices
  scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)

  # rescale the sampling correlation matrix by the appropriate diagonals
  V <- v.out * tcrossprod(scaleO)


  # name traits according to the names of the input summart statistics
  # use general format of V1-VX if no names provided
  colnames(S) <- if (is.null(names(munged_sumstats))) paste0("V", 1:ncol(S)) else names(munged_sumstats)
  rownames(S) <- if (is.null(names(munged_sumstats))) paste0("V", 1:ncol(S)) else names(munged_sumstats)

  if (mean(Liab.S) != 1) {
    r <- nrow(S)
    SE <- matrix(0, r, r)
    SE[lower.tri(SE, diag = TRUE)] <- sqrt(diag(V))

    colnames(SE) <- colnames(S)
    rownames(SE) <- rownames(S)
  }


  if (all(diag(S) > 0)) {
    ## calculate standardized results to print genetic correlations to log and screen
    ratio <- tcrossprod(1 / sqrt(diag(S)))
    S_Stand <- S * ratio

    # calculate the ratio of the rescaled and original S matrices
    scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)

    ## Make sure that if ratio in NaN (devision by zero) we put the zero back in
    # -> not possible because of 'all(diag(S) > 0)'
    # scaleO[is.nan(scaleO)] <- 0

    # rescale the sampling correlation matrix by the appropriate diagonals
    V_Stand <- V * tcrossprod(scaleO)

    # enter SEs from diagonal of standardized V
    r <- nrow(S)
    SE_Stand <- matrix(0, r, r)
    SE_Stand[lower.tri(SE_Stand, diag = TRUE)] <- sqrt(diag(V_Stand))

    colnames(SE_Stand) <- colnames(S)
    rownames(SE_Stand) <- rownames(S)

  } else {
    cli::cli_alert_warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
  }


  ind <- which(lower.tri(S, diag = F), arr.ind = TRUE)

  rg_res <- tibble(
    trait1 = dimnames(S_Stand)[[2]][ind[, 2]],
    trait2 = dimnames(S_Stand)[[1]][ind[, 1]],
    rg = S_Stand[ind],
    rg_se = SE_Stand[ind],
    rg_p = 2*pnorm(abs(rg/rg_se), lower.tail = FALSE)
  )

  return(
    list(
      h2 = h2_res,
      rg = rg_res,
      raw = list(V = V, S = S, I = I, N = N.vec, m = m, V_Stand = V_Stand, S_Stand = S_Stand, SE_Stand = SE_Stand)
    )
  )
}
