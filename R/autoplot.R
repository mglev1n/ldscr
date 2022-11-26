autoplot.ldscr_list <- function(object, ...){
  result <- object$raw

  r <- nrow(result$S)
  result$Pval_Stand <- 2*pnorm(abs(result$S_Stand/result$SE_Stand), lower.tail = FALSE)
  result$Pval_Stand <- as.matrix(Matrix::forceSymmetric(result$Pval_Stand,uplo = "L"))
  result$trait_names <- colnames(result$S_Stand)
  rownames(result$S_Stand) <- colnames(result$S_Stand)
  ldsc_res <- result
  rownames(ldsc_res$Pval_Stand) <- ldsc_res$trait_names
  colnames(ldsc_res$Pval_Stand) <- ldsc_res$trait_names
  rownames(ldsc_res$SE_Stand) <- ldsc_res$trait_names
  colnames(ldsc_res$SE_Stand) <- ldsc_res$trait_names

  plot_df <- corrr::as_cordf(ldsc_res$S_Stand) %>%
    corrr::rearrange(method = "HC") %>%
    corrr::stretch() %>%
    left_join(
      corrr::as_cordf(ldsc_res$Pval_Stand) %>%
        corrr::stretch() %>%
        rename(pval = r),
      by = c("x", "y")
    ) %>%
    mutate(
      x = stringr::str_replace(x, "_", " "),
      y = stringr::str_replace(y, "_", " ")
    ) %>%
    mutate(
      x = forcats::fct_inorder(x),
      y = forcats::fct_inorder(y)
    ) %>%
    mutate(r = dplyr::case_when(is.na(r) ~ 1, TRUE ~ r)) %>%
    mutate(pval_bonferroni = pval < 0.05 / ((dplyr::n_distinct(.$x)^2 - dplyr::n_distinct(.$x)) / 2)) %>%
    mutate(label = dplyr::case_when(
      pval_bonferroni ~ as.character(glue::glue("P < {signif(0.05 / ((dplyr::n_distinct(.$x)^2 - dplyr::n_distinct(.$x))/2), 2)}")),
      TRUE ~ NA_character_
    ))


  plot_df %>%
    ggplot(aes(x, y, fill = r)) +
    geom_tile() +
    geom_point(aes(shape = label)) +
    scale_fill_gradient2(low = "#2c477a", mid = "white", high = "#ad171c", name = "*r<sub>g</sub>*") +
    scale_shape_manual(values = c(8), na.translate = F, name = NULL) +
    guides(fill = guide_colourbar(order = 1)) +
    labs(
      x = NULL,
      y = NULL
    ) +
    coord_equal() +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.title = ggtext::element_markdown()
    )
}
