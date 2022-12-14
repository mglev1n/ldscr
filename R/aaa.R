# nocov start

# Onload -----------------------------------------------------------------------

## Taken from https://github.com/tidyverse/dplyr/blob/d310ad1cef1c14d770c94e1a9a4c79c888f46af6/R/zzz.r#L2-L9

.onLoad <- function(libname, pkgname) {
  # dynamically register autoplot methods
  s3_register("ggplot2::autoplot", "ldscr_list")

  invisible()
}

# Dynamic reg helper -----------------------------------------------------------

# vctrs/register-s3.R
# https://github.com/r-lib/vctrs/blob/master/R/register-s3.R
s3_register <- function(generic, class, method = NULL) {
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)

  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]

  if (is.null(method)) {
    method <- get(paste0(generic, ".", class), envir = parent.frame())
  }
  stopifnot(is.function(method))

  if (package %in% loadedNamespaces()) {
    registerS3method(generic, class, method, envir = asNamespace(package))
  }

  # Always register hook in case package is later unloaded & reloaded
  setHook(
    packageEvent(package, "onLoad"),
    function(...) {
      registerS3method(generic, class, method, envir = asNamespace(package))
    }
  )
}

# nocov end
