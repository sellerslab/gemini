
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![DOI](https://zenodo.org/badge/175870293.svg)](https://zenodo.org/badge/latestdoi/175870293)

[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

gemini
======

The `gemini` package allows users to analyze combinatorial CRISPR screens as described in Zamanighomi et al. ?.

Installation
------------

To install the latest development version of `gemini`, use the [`devtools`](%22https://github.com/r-lib/devtools%22) package as follows:

``` r
devtools::install_github(repo = "sellerslab/gemini", build_vignettes = TRUE)
```

You can (eventually) install the stable release version of gemini from [Bioconductor](https://www.bioconductor.org/) with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("gemini")
```

Details
-------

See the vignette for usage instructions:

``` r
vignette("gemini-quickstart", package = "gemini")
#> Warning: vignette 'gemini-quickstart' not found
```
