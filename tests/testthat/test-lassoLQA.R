context("output result")
library(lassoLQA)

test_that("lassoLQA", {
  expect_error(lassoLQA(matrix(rnorm(90),30,3), matrix(rnorm(90),30,3)%*%c(2,1,0)+rnorm(30,0,4)), NA)
  expect_error(lassoLQA(matrix(rnorm(1200),200,6), matrix(rnorm(1200),200,6)%*%c(3,2,1,0,0,0)+rnorm(200,0,4)),  NA)
})
