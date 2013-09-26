library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))
context("1dSSA")

test_that("Lcov works correctly for marginal L values", {
  Ns <- c(1005, 1500, 2000, 5000)

  set.seed(1)
  for (N in Ns) {
    for (L in c(1, N)) {
      F <- rcauchy(N)
      C.exact <- tcrossprod(hankel(F, L))
      C.fast <- Lcov.matrix(F, L = L)
      expect_equal(C.fast, C.exact,
                   info = sprintf("L = %d, N = %d", L, N))
    }
  }
})
