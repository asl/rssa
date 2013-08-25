library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))
source(system.file("extdata", "mssa.batching.R", package = "Rssa"))
context("2dSSA")

test_that("All kinds of 1d-forecasting and reconstruction work correctly for finite rank series", {
set.seed(1)
N <- 71
sigma <- 0
Ls <- c(12, 24, 36, 48, 60)
len <- 400

examples.all <- list(list(xs = list(cosine(N + len, 30, 12), cosine(N + len, 20, 12, pi / 4)), ssa.dim = 2, mssa.dim = 2, cssa.dim = 2),
                     list(xs = list(cosine(N + len, 30, 12), cosine(N + len, 30, 12, pi / 2)), ssa.dim = 2, mssa.dim = 2, cssa.dim = 1),
                     list(xs = list(cosine(N + len, 30, 12), cosine(N + len, 20, 8,  pi / 4)), ssa.dim = 2, mssa.dim = 4, cssa.dim = 4))

examples.mssa.dlen <- list(list(xs = list(cosine(N + len + 10, 30, 12), cosine(N + len, 20, 12, pi / 4)), mssa.dim = 2),
                           list(xs = list(cosine(N + len + 10, 30, 12), cosine(N + len, 30, 12, pi / 2)), mssa.dim = 2),
                           list(xs = list(cosine(N + len + 10, 30, 12), cosine(N + len, 20, 8,  pi / 4)), mssa.dim = 4))

for (svd.method in c("svd", "eigen", "nutrlan", "propack")) {
  rec.all <- all.MSE(examples.all,
                     N = 1,
                     sigma = sigma,
                     Ls = Ls,
                     type = "reconstruct",
                     len = len,
                     svd.method = svd.method,
                     eval.sd = FALSE)


  expect_true(all(abs(unlist(rec.all)) < .Machine$double.eps^.5),
              label = sprintf("Reconstruction, series of same length, svd.method = %s", svd.method))

  fore.all <- all.MSE(examples.all,
                      N = 1,
                      sigma = sigma,
                      Ls = Ls,
                      type = "forecast",
                      len = len,
                      svd.method = svd.method,
                      eval.sd = FALSE)

  expect_true(all(abs(unlist(fore.all)) < .Machine$double.eps^.5),
              label = sprintf("Forecast, series of same length, svd.method = %s", svd.method))

  expect_true(all(abs(unlist(fore.all)) < .Machine$double.eps^.5),
              label = sprintf("MSSA reconstruction, series of different lengths, svd.method = %s", svd.method))


  fore.mssa.dlen <- all.MSE(examples.mssa.dlen,
                            kinds = c("r-mssa-row", "r-mssa-column", "v-mssa-row", "v-mssa-column"),
                            N = 1,
                            sigma = sigma,
                            Ls = Ls,
                            type = "forecast",
                            len = len,
                            svd.method = svd.method,
                            eval.sd = FALSE)
  expect_true(all(abs(unlist(fore.all)) < .Machine$double.eps^.5),
              label = sprintf("MSSA forecast, series of different lengths, svd.method = %s", svd.method))
}
})
