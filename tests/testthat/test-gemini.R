context("loading")

test_that("package loads successfully",
          require('gemini')
          )

test_that("package data loads successfully",
          {
            expect_true(object = is.matrix(counts))
          })