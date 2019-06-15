
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![DOI](https://zenodo.org/badge/175870293.svg)](https://zenodo.org/badge/latestdoi/175870293)

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
source("https://www.bioconductor.org/biocLite.R")
biocLite("gemini")
```

Quickstart
----------

Given a counts matrix and the relevant annotations, GEMINI uses the following workflow:

counts (with replicate annotation and guide annotations) → `gemini.input` → `gemini.model` → `gemini.score`

Each of these steps can be seen here:

``` r
library(gemini)
data("example-bigPapi", package = "gemini")

# COUNTS -> gemini.input
gemini.input <- gemini_create_input(counts.matrix = counts,
                                    sample.replicate.annotation = sample.replicate.annotation,
                                    guide.annotation = guide.annotation,
                                    gene.column.names = c("U6.gene", "H1.gene"),
                                    sample.column.name = "samplename")

gemini.input %<>% gemini_calculate_lfc()

# gemini.input -> gemini.model
gemini.model <- gemini_initialize(Input = gemini.input,
                                  nc_gene = "CD81",
                                  verbose = F)
gemini.model %<>% gemini_inference(verbose = F)

# gemini.model -> gemini.score 
# Here, we use EEF2, a known essential gene, as the positive control.
# Additionally, a set of essential genes can be used instead.
gemini.score <- gemini_score(gemini.model, pc_gene = "EEF2")

# To calculate significance, identify a set of non-interacting genes - here, we use genes paired with other negative controls (HPRT intron, 6T)
noninteracting_pairs <- rownames(gemini.model$s)[grep("6T|HPRT",rownames(gemini.model$s))]
gemini.score <- gemini_score(gemini.model, pc_gene = "EEF2", nc_pairs = noninteracting_pairs)
```

Visualization
-------------

#### Convergence

After running `gemini_inference`, the resulting convergence rate can be visualized. If the model is divergent, alternative parameters for the priors should be selected until convergence is achieved. This can be done through a form of cross-validation. See the manual for more details.

``` r
gemini_plot_mae(gemini.model)
```

#### Interactions

If a genetic interaction is identified, check the interaction using the `gemini_boxplot` function.

``` r
gemini_boxplot(Model = gemini.model, gene.column.names = c("U6.gene", "H1.gene"), g = "BRCA2", h = "PARP1", nc_gene = "CD81", sample = "A549", show_inference = T, identify_guides = T)

gemini_boxplot(Model = gemini.model, gene.column.names = c("U6.gene", "H1.gene"), g = "BRCA2", h = "PARP1", nc_gene = "CD81", sample = "A549", show_inference = T, color_x = T)
```
