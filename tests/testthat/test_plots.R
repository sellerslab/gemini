library("gemini")

test_that("plots are working", {
    data("Model", package = "gemini")
    gg = gemini_plot_mae(Model)
    testthat::expect_is(object = gg, class = "ggplot")
    
    gg2 = gemini_boxplot(Model, g = "BRCA1", h = "BRCA2", nc_gene = "CD81", sample = "A549")
    testthat::expect_is(object = gg2, class = "ggplot")
})