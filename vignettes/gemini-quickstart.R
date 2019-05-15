## ------------------------------------------------------------------------
library("gemini")
data("counts", "guide.annotation", "sample.replicate.annotation", package = "gemini")

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(head(counts[,1:5]), caption = "Counts matrix", align = 'l')

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(head(sample.replicate.annotation), caption = "Sample/replicate annotations")
knitr::kable(head(guide.annotation[,1:3]), caption = "Guide/gene annotation")

## ------------------------------------------------------------------------
Input <- gemini_create_input(counts.matrix = counts,
                    sample.replicate.annotation = sample.replicate.annotation,
                    guide.annotation = guide.annotation,
                    ETP.column = 'pDNA', 
                    gene.column.names = c("U6.gene", "H1.gene"),
                    sample.column.name = "samplename",
                    verbose = TRUE)

# Note: ETP column can also be specified by column index 
# (e.g. ETP.column = c(1))

## ------------------------------------------------------------------------
Input %<>% gemini_calculate_lfc(normalize = TRUE, 
                                CONSTANT = 32)

## ------------------------------------------------------------------------
Model <- gemini_initialize(Input = Input, 
                  nc_gene = "CD81", 
                  pattern_join = ';',
                  pattern_split = ';', 
                  cores = 1,
                  verbose = TRUE)

## ------------------------------------------------------------------------
Model %<>% gemini_inference(cores = 1,
                            verbose = FALSE)

## ---- eval=T-------------------------------------------------------------
gemini_plot_mae(Model)

## ------------------------------------------------------------------------
# Use non-interacting gene pairs as nc_pairs. A caveat here is that this set is constructed only using negative controls! This probably leads to biased sampling of the null distribution, thereby overestimating the number of significant hits, but still is useful in this case.
nc_pairs <- grep("6T|HPRT", rownames(Model$s), value = TRUE)

# An example of some nc_pairs...
head(nc_pairs, n = 5)

Score <- gemini_score(Model = Model,
             pc_gene = "EEF2",
             nc_pairs = nc_pairs)

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(Score$strong[order(Score$strong[,"A549"], decreasing = TRUE)[1:10],], caption = "Strong scores for top 10 interactions from A549 in all samples")

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(Score$fdr_strong[order(Score$fdr_strong[,"A549"], decreasing = FALSE)[1:10],], caption = "FDRs for top 10 interactions in A549")

## ---- eval = F-----------------------------------------------------------
#  gemini_boxplot(Model = Model,
#                 g = "BRCA2",
#                 h = "PARP1",
#                 nc_gene = "CD81",
#                 sample = "A549",
#  			   show_inference = TRUE,
#  			   identify_guides = TRUE
#  			   )

## ---- eval = F-----------------------------------------------------------
#  gemini_boxplot(Model = Model,
#                 g = "BRCA2",
#                 h = "PARP1",
#                 nc_gene = "CD81",
#                 sample = "A549",
#  			   show_inference = TRUE,
#  			   color_x = TRUE
#  			   )

