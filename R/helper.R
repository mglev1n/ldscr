# Read ld from either internal or external file
read_ld <- function(ancestry, ld) {
  if (missing(ancestry)) {
    x <- fs::dir_ls(ld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  } else {
    x <- ldscore_files(ancestry, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  }
  return(x)
}

# Read wld from either internal or external file
read_wld <- function(ancestry, wld) {
  if (missing(ancestry)) {
    w <- fs::dir_ls(wld, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  } else {
    w <- ldscore_files(ancestry, glob = "*.l2.ldscore.gz") %>%
      vroom::vroom(col_types = vroom::cols())
  }
  return(w)
}

# Read M from either internal or external file
read_m <- function(ancestry, ld) {
  if (missing(ancestry)) {
    m <- fs::dir_ls(ld, glob = "*.l2.M_5_50") %>%
      vroom::vroom(col_types = vroom::cols(), col_names = FALSE, delim = "\t")
  } else {
    m <- ldscore_files(ancestry, glob = "*.l2.M_5_50") %>%
      vroom::vroom(col_types = vroom::cols(), col_names = FALSE, delim = "\t")
  }
  return(m)
}


# Read summary statistics from either internal or external file
read_sumstats <- function(munged_sumstats, name) {
  # Check if name is present to enable more informative logging
  if (missing(name)) {
    if (is.character(munged_sumstats)) {
      cli::cli_progress_step("Reading summary statistics from {munged_sumstats}")
      sumstats_df <- vroom::vroom(munged_sumstats, col_types = vroom::cols())
    } else {
      cli::cli_progress_step("Reading summary statistics from dataframe")
      sumstats_df <- munged_sumstats
    }
  } else {
    if (is.character(munged_sumstats)) {
      cli::cli_progress_step("Reading summary statistics for '{name}' from {munged_sumstats}")
      sumstats_df <- vroom::vroom(.x, col_types = vroom::cols())
    } else {
      cli::cli_progress_step("Reading summary statistics for '{name}' from dataframe")
      sumstats_df <- .x
    }
    sumstats_df <- na.omit(sumstats_df)

    return(sumstats_df)
  }
}

merge_sumstats <- function(sumstats_df, w, x, chr_filter) {
  merged <- sumstats_df %>%
    dtplyr::lazy_dt() %>%
    dplyr::select(SNP, N, Z, A1) %>%
    dplyr::inner_join(w[, c("SNP", "wLD")], by = c("SNP")) %>%
    dplyr::inner_join(x, by = c("SNP")) %>%
    dplyr::arrange(CHR, BP) %>%
    dplyr::filter(CHR %in% chr_filter) %>%
    na.omit() %>%
    unique() %>%
    tibble::as_tibble()

  return(merged)
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
