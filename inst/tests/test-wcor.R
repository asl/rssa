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
