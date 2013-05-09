library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))

context("W-correlations")

test_that("wcor method returns proper matrix", {
  set.seed(1)
  N <- 48
  L <- 24
  v <- rnorm(N)
  ss <- ssa(v, L = L, svd.method = "eigen")
  w <- wcor(ss)

  expect_true(all(abs(w) <= 1))
  expect_true(all(diag(w) == 1))
  expect_equal(w, t(w))
})
