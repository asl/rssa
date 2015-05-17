library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))
context("1dSSA")

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

test_that("Marginal case for shaped 2dSSA", {
  s <- ssa(rbind(rep(1, 20), c(rep(1, 10), rep(NA, 10))), L = c(2, 10),  kind = "2d-ssa")
  expect_equal(dim(s$fmask), c(1, 11))
  expect_equal(s$lambda[1:5]^2, c(20, 0, 0, 0, 0))


  s <- ssa(cbind(rep(1, 20), c(rep(1, 10), rep(NA, 10))), L = c(10, 2),  kind = "2d-ssa")
  expect_equal(dim(s$fmask), c(11, 1))
  expect_equal(s$lambda[1:5]^2, c(20, 0, 0, 0, 0))
})
