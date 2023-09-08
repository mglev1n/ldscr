
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ldscr

<!-- badges: start -->

[![R-CMD-check](https://github.com/mglev1n/ldscr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mglev1n/ldscr/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/mglev1n/ldscr/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/mglev1n/ldscr/actions/workflows/test-coverage.yaml)

<!-- badges: end -->

The goal of ldscr is to provide functionality to estimate genetic
heritability and cross-trait genetic correlations from GWAS summary
statistics using LD score regression within R. Details of LD score
regression for estimation of heritabliity and genetic correlations have
been previously published.<sup>1,2</sup> This package adapts code and
functionality originally implemented in
[GenomicSEM](https://github.com/GenomicSEM/GenomicSEM).<sup>3</sup>

## Installation

You can install the development version of ldscr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")

devtools::install_github("mglev1n/ldscr")
devtools::install_github("skoyamamd/ldscr")
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

## References

<div id="refs" class="references csl-bib-body" line-spacing="2">

<div id="ref-Bulik-Sullivan2015" class="csl-entry">

<span class="csl-left-margin">1. </span><span
class="csl-right-inline">Bulik-Sullivan, B. *et al.* [An atlas of
genetic correlations across human diseases and
traits](https://doi.org/10.1038/ng.3406). *Nature Genetics* **47**,
1236–1241 (2015).</span>

</div>

<div id="ref-Bulik-Sullivan2015a" class="csl-entry">

<span class="csl-left-margin">2. </span><span
class="csl-right-inline">Bulik-Sullivan, B. K. *et al.* [LD Score
regression distinguishes confounding from polygenicity in genome-wide
association studies](https://doi.org/10.1038/ng.3211). *Nature Genetics*
**47**, 291–295 (2015).</span>

</div>

<div id="ref-Grotzinger2019" class="csl-entry">

<span class="csl-left-margin">3. </span><span
class="csl-right-inline">Grotzinger, A. D. *et al.* [Genomic structural
equation modelling provides insights into the multivariate genetic
architecture of complex
traits](https://doi.org/10.1038/s41562-019-0566-x). *Nature Human
Behaviour* **3**, 513–525 (2019).</span>

</div>

</div>
