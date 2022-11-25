#' Paths to LD-score files from Pan-UKB
#'
#' @param ancestry (character) One of "AFR", "AMR", "CSA", "EAS", "EUR", or "MID"
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
#' @param example (character) One of "BMI" or "LDL" which have been included as example traits.
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
