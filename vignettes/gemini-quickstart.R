## ------------------------------------------------------------------------
library("gemini")
require("magrittr")
data("example-bigPapi", package = "gemini")

## ---- echo = F-----------------------------------------------------------
knitr::kable(head(counts[,1:5]), caption = "Counts matrix", align = 'l')

## ---- echo = F-----------------------------------------------------------
knitr::kable(head(sample.replicate.annotation), caption = "Sample/replicate annotations")
knitr::kable(head(guide.annotation), caption = "Guide/gene annotation")

## ------------------------------------------------------------------------
Input <- gemini_create_input(counts.matrix = counts,
                    sample.replicate.annotation = sample.replicate.annotation,
                    guide.annotation = guide.annotation,
                    ETP.column = 'pDNA', # this can also be specified by column index
                    gene.column.names = c("U6.gene", "H1.gene"),
                    sample.column.name = "samplename",
                    verbose = T)

## ------------------------------------------------------------------------
Input %<>% gemini_calculate_lfc(normalize = T, 
                                CONSTANT = 32)

## ------------------------------------------------------------------------
Model <- gemini_initialize(Input = Input, 
                  nc_gene = "HPRT intron", 
                  pattern_join = '-',
                  pattern_split = ';', 
                  cores = 1,
                  verbose = T)

## ------------------------------------------------------------------------
Model %<>% gemini_inference(cores = 1,
                            verbose = T,
							save_iterations = F)

## ------------------------------------------------------------------------
gemini_plot_mae(Model)

## ------------------------------------------------------------------------
# Use other negative controls or known non-interacting gene pairs as nc_pairs
nc_pairs <- grep("6T|CD81", rownames(Model$s), value = T)
# Some nc_pairs...
head(nc_pairs, n = 5)

Score <- gemini_score(Model = Model,
             pc_gene = "EEF2",
             nc_pairs = nc_pairs)

## ---- echo = F-----------------------------------------------------------
knitr::kable(Score$fdr_strong[order(Score$fdr_strong[,"A549"], decreasing = F)[1:10],], caption = "FDRs for top 10 interactions for A549")

## ------------------------------------------------------------------------
gemini_boxplot(Model = Model, 
               gene.column.names = c("U6.gene", "H1.gene"), 
               g = "BRCA2",
               h = "BRCA1",
               nc_gene = "HPRT intron",
               sample = "A549",
			   show_inference = T, 
			   identify_guides = T
			   )

