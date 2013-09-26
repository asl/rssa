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

test_that("1dSSA works correctly for marginal L values", {
  Ns <- c(15, 150, 200)

  set.seed(1)
  for (N in Ns) {
    for (L in c(1, N)) {
      for (svd.method in c("propack", "eigen", "svd")) {
        F <- rcauchy(N)
        ss <- ssa(F, L = L, svd.method = svd.method)

        expect_equal(wnorm(ss), sqrt(sum(F^2)),
                     info = sprintf("L = %d, N = %d, svd.method = %s", L, N, svd.method))
        expect_equal(reconstruct(ss, 1)$F1, F,
                     info = sprintf("L = %d, N = %d, svd.method = %s", L, N, svd.method))
      }
    }
  }
})

test_that("Toeplitz SSA works correctly for marginal L values", {
  Ns <- c(15, 150, 200)

  set.seed(1)
  for (N in Ns) {
    for (L in 1) {
      for (svd.method in c("propack", "eigen", "svd")) {
        F <- rcauchy(N)
        ss <- ssa(F, kind = "toeplitz-ssa", L = L, svd.method = svd.method)

        expect_equal(wnorm(ss), sqrt(sum(F^2)),
                     info = sprintf("L = %d, N = %d, svd.method = %s", L, N, svd.method))
        expect_equal(reconstruct(ss, 1)$F1, F,
                     info = sprintf("L = %d, N = %d, svd.method = %s", L, N, svd.method))
      }
    }
  }
})
