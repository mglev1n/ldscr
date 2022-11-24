#' Estimate cross-trait genetic correlations
#'
#' @description
#'
#' `ldsc_rg()` uses ldscore regression to estimate the pairwise genetic correlations between traits.
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

}
