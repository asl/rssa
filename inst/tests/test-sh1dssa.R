library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))

context("Shaped 1dSSA")

test_that("Shaped 1d SSA works correctly", {
  set.seed(1)
  N <- 200
  L <- 40
  rank <- 6
  F <- rnorm(N)
  F[50] <- F[101] <- F[143] <- NA

  wmask0 <- wmask1 <- wmask2 <- rep(TRUE, L)
  wmask1[2] <- wmask1[30] <- wmask2[23] <- FALSE
  wmasks <- list(wmask0, wmask1, wmask2)
  for (wmask in wmasks) {
    circulars <- c(FALSE, TRUE)
    for (circular in circulars) {
      s2 <- ssa(F, L = c(L, 1), wmask = as.matrix(wmask),
                kind = "2d-ssa", circular = c(circular, FALSE), neig = rank + 1)

      svd.methods <- c("svd", "nutrlan", "propack", "eigen")
      for (svd.method in svd.methods) {
        s1 <- ssa(F, L = L, wmask = wmask,
                  svd.method = svd.method, neig = rank + 1, circular = circular)

        for (r in seq_len(rank)) {
          expect_equal(reconstruct(s1, r)$F1, reconstruct(s2, r)$F1,
                       info = sprintf("component = %d, svd.method = %s", r, svd.method))
        }
      }
    }
  }
})
