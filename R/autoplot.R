autoplot.ldscr_list <- function(object, remove_na = TRUE, show_na_info = TRUE, ...) {
  result <- object$raw

  r <- nrow(result$S)

  # Handle case where S_Stand doesn't exist (shouldn't happen with robust version, but safety check)
  if (is.null(result$S_Stand)) {
    cli::cli_alert_warning("Standardized correlation matrix not found. Cannot create plot.")
    return(NULL)
  }

  # Calculate p-values, handling NAs
  result$Pval_Stand <- matrix(NA, nrow = r, ncol = r)
  for (i in 1:r) {
    for (j in 1:r) {
      if (!is.na(result$S_Stand[i, j]) && !is.na(result$SE_Stand[i, j]) && result$SE_Stand[i, j] > 0) {
        result$Pval_Stand[i, j] <- 2 * pnorm(abs(result$S_Stand[i, j] / result$SE_Stand[i, j]), lower.tail = FALSE)
      }
    }
  }

  # Ensure matrices are symmetric
  result$Pval_Stand <- as.matrix(Matrix::forceSymmetric(result$Pval_Stand, uplo = "L"))
  result$S_Stand <- as.matrix(Matrix::forceSymmetric(result$S_Stand, uplo = "L"))

  # Set trait names
  result$trait_names <- colnames(result$S_Stand)
  if (is.null(result$trait_names)) {
    result$trait_names <- paste0("Trait", 1:r)
    colnames(result$S_Stand) <- result$trait_names
    rownames(result$S_Stand) <- result$trait_names
  }

  rownames(result$S_Stand) <- result$trait_names
  rownames(result$Pval_Stand) <- result$trait_names
  colnames(result$Pval_Stand) <- result$trait_names
  rownames(result$SE_Stand) <- result$trait_names
  colnames(result$SE_Stand) <- result$trait_names

  # Create the plotting dataframe
  plot_df <- corrr::as_cordf(result$S_Stand) %>%
    corrr::stretch() %>%
    dplyr::left_join(
      corrr::as_cordf(result$Pval_Stand) %>%
        corrr::stretch() %>%
        dplyr::rename(pval = r),
      by = c("x", "y")
    ) %>%
    dplyr::mutate(
      x = stringr::str_replace(x, "_", " "),
      y = stringr::str_replace(y, "_", " ")
    )

  # Count NAs before potential removal
  na_count <- sum(is.na(plot_df$r) & !is.na(plot_df$x) & !is.na(plot_df$y))
  total_pairs <- nrow(plot_df[!is.na(plot_df$x) & !is.na(plot_df$y), ])

  # Handle NAs based on remove_na parameter
  if (remove_na) {
    # Keep rows where correlation is NOT NA, or where it's the diagonal
    plot_df <- plot_df %>%
      dplyr::filter(!is.na(r) | x == y)

    # Show info about removed correlations
    if (show_na_info && na_count > 0) {
      cli::cli_alert_info("Removed {na_count} correlation{?s} with NA values from plot (out of {total_pairs} total)")
    }
  } else {
    # Keep NAs but mark them for special handling in the plot
    plot_df <- plot_df %>%
      dplyr::mutate(is_na = is.na(r))
  }

  # Set diagonal to 1 for traits that aren't NA
  plot_df <- plot_df %>%
    dplyr::mutate(r = dplyr::case_when(
      x == y & !is.na(r) ~ 1,
      x == y & is.na(r) ~ NA_real_,
      TRUE ~ r
    )) %>%
    dplyr::mutate(
      x = forcats::fct_inorder(x),
      y = forcats::fct_inorder(y)
    )

  # Calculate Bonferroni correction based on valid (non-NA) correlations
  valid_correlations <- plot_df %>%
    dplyr::filter(!is.na(r), x != y) %>%
    nrow()

  bonferroni_threshold <- if (valid_correlations > 0) {
    0.05 / valid_correlations
  } else {
    0.05
  }

  plot_df <- plot_df %>%
    dplyr::mutate(
      pval_bonferroni = !is.na(pval) & pval < bonferroni_threshold,
      label = dplyr::case_when(
        pval_bonferroni ~ as.character(glue::glue("P < {signif(bonferroni_threshold, 2)}")),
        TRUE ~ NA_character_
      )
    )

  # Create the base plot
  p <- plot_df %>%
    ggplot(aes(x, y, fill = r))

  if (remove_na) {
    # Standard heatmap without NA handling
    p <- p +
      geom_tile() +
      geom_point(aes(shape = label)) +
      scale_fill_gradient2(
        low = "#2c477a",
        mid = "white",
        high = "#ad171c",
        name = "*r<sub>g</sub>*",
        na.value = "grey90"
      )
  } else {
    # Heatmap with NA visualization
    p <- p +
      geom_tile(aes(alpha = ifelse(is.na(r), 0.3, 1))) +
      geom_point(aes(shape = label)) +
      scale_fill_gradient2(
        low = "#2c477a",
        mid = "white",
        high = "#ad171c",
        name = "*r<sub>g</sub>*",
        na.value = "grey50"
      ) +
      scale_alpha_identity() +
      # Add text for NA cells
      geom_text(
        data = plot_df %>% dplyr::filter(is.na(r), x != y),
        aes(label = "NA"),
        color = "darkgrey",
        size = 3,
        fill = NULL
      )
  }

  # Complete the plot
  p <- p +
    scale_shape_manual(values = c(8), na.translate = FALSE, name = NULL) +
    guides(fill = guide_colourbar(order = 1)) +
    labs(
      x = NULL,
      y = NULL,
      caption = if (!remove_na && na_count > 0) {
        glue::glue("Grey cells indicate correlations that could not be estimated (n = {na_count})")
      } else {
        NULL
      }
    ) +
    coord_equal() +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.title = ggtext::element_markdown(),
      plot.caption = element_text(size = 10, color = "grey60")
    )

  return(p)
}
