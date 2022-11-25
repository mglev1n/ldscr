
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ldscr

<!-- badges: start -->

[![R-CMD-check](https://github.com/mglev1n/ldscr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mglev1n/ldscr/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/mglev1n/ldscr/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/mglev1n/ldscr/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

The goal of ldscr is to provide functionality to estimate genetic
heritability and cross-trait genetic correlations from GWAS summary
statistics using LD score regression within R.

## Installation

You can install the development version of ldscr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mglev1n/ldscr")
```

## Usage

`ldsc_h2()` can be used to estimate heritability. Sample GWAS data is
provided in `sumstats_munged_example()`. Users can utilize built-in LD
reference data, or provide their own.

``` r
library(ldscr)
df <- sumstats_munged_example(example = "BMI")
h2_res <- ldsc_h2(munged_sumstats = df, ancestry = "EUR")
```

`ldsc_rg()` can be used to estimate cross-trait genetic correlations
between two or more traits.

``` r
rg_res <- ldsc_rg(
  munged_sumstats = list(
    "APOB" = sumstats_munged_example(example = "APOB"),
    "LDL" = sumstats_munged_example(example = "LDL")
  ),
  ancestry = "EUR"
)
```
