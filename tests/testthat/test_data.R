library("gemini")

testthat::test_that("data is correct", code = ({
        data("counts", package = "gemini")
        data("guide.annotation", package = "gemini")
        data("sample.replicate.annotation", package = "gemini")
        data("Input", package = "gemini")
        data("Model", package = "gemini")
        
        # Make sure that dimensions match up
        expect_equal(nrow(counts), nrow(guide.annotation))
        expect_equal(rownames(counts), guide.annotation[,1])
        
        expect_equal(ncol(counts), nrow(sample.replicate.annotation))
        expect_equal(colnames(counts), sample.replicate.annotation[,1])
}))

testthat::test_that("Input object is reproducible", code = ({
    data("counts", package = "gemini")
    data("guide.annotation", package = "gemini")
    data("sample.replicate.annotation", package = "gemini")
    Input.new <- gemini_create_input(counts.matrix = counts, 
                                     sample.replicate.annotation = sample.replicate.annotation,
                                     guide.annotation = guide.annotation,
                                     sample.column.name = 'samplename',
                                     gene.column.names = c("U6.gene", "H1.gene"),
                                     ETP.column = 1,
                                     verbose = TRUE)
    Input.new %<>% gemini_calculate_lfc()
    data("Input", package = "gemini")
    expect_equal(Input.new, Input)
})
)
