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

test_that("wcor method returns proper matrix for 2dSSA", {
  set.seed(1)
  N <- c(48, 49)
  L <- c(24, 25)
  mx <- matrix(rnorm(prod(N)), N[1], N[2])
  ss <- ssa(mx, L = L, kind = "2d-ssa", svd.method = "propack", neig = 10)
  w <- wcor(ss)

  expect_true(all(abs(w) <= 1))
  expect_true(all(diag(w) == 1))
  expect_equal(w, t(w))
})

test_that("wcor method return correct matrix for 2dSSA case", {
  env <- new.env()
  load(file = system.file("extdata", "wcor.2dssa.testdata.rda", package = "Rssa"), envir = env)

  set.seed(1)
  mx <- outer(1:50, 1:50,
              function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
        rnorm(50^2, sd = 0.1)

  for (svd.method in c("nutrlan", "propack")) {
    s <- ssa(mx, kind = "2d-ssa", svd.method = svd.method)
    w <- wcor(s, groups = 1:12)

    expect_equal(w, env$w,
                 info = sprintf("wcor.2d.ssa.%s", svd.method))
  }
})

test_that("Hankel weights computed correctly for marginal cases", {
expect_equal(.hweights.default(1, 1), 1)
expect_equal(.hweights.default(10, 1), rep(1, 10))
expect_equal(.hweights.default(10, 10), rep(1, 10))
})

test_that("Hankel weights computed correctly for common case", {
expect_equal(.hweights.default(5, 2), c(1, 2, 2, 2, 1))
expect_equal(.hweights.default(5, 3), c(1, 2, 3, 2, 1))
expect_equal(.hweights.default(5, 4), c(1, 2, 2, 2, 1))
})

test_that("`wnorm' works correctly for MSSA", {
Nss <- list(20,
            c(17, 17),
            14,
            c(14, 32, 36, 36, 31, 37))

set.seed(1)
for (Ns in Nss) {
  f <- lapply(Ns, rnorm)

  sss <- lapply(f, ssa, kind = "1d-ssa", L = 13,
                force.decompose = FALSE)
  ss <- ssa(f, kind = "mssa", L = 13,
            force.decompose = FALSE)

  w1 <- sum(sapply(sss, wnorm) ^ 2)
  w2 <- wnorm(ss) ^ 2
  expect_equal(w2, w1, label = sprintf("lengths: %s",
                                       paste0(Ns, collapse = ", ")))
}
})
