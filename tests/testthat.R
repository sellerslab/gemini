library(testthat)
library(gemini)

test_check("gemini")

testthat::test_that("data can be imported", code = ({
	Input <- gemini_create_input(counts.matrix = gemini::counts, 
						sample.replicate.annotation = gemini::sample.replicate.annotation,
						guide.annotation = gemini::guide.annotation,
						sample.column.name = 'samplename',
						gene.column.names = c("U6.gene", "H1.gene"),
						ETP.column = 1,
						verbose = TRUE)}))