library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))
context("Cadzow")

test_that("Cadzow limit is a series of finite rank", {
  s <- ssa(co2)
  eps <- sqrt(.Machine$double.eps)

  ranks <- 1:5
  for (rank in ranks) {
    cz <- cadzow(s, rank = rank)
    expect_true(high.rank.rate(cz, rank = rank, ssaobj = s) < eps)
  }
})

.series.wsqdistance <- function(F1, F2, weights = 1) {
  mask <- weights > 0

  weights <- weights[mask]
  F1 <- as.vector(unlist(F1))[mask]
  F2 <- as.vector(unlist(F2))[mask]

  sum(weights * (F1-F2)^2)
}

test_that("Cadzow correction really works", {
  s <- ssa(co2)
  eps <- sqrt(.Machine$double.eps)
  delta <- 0.0001
  w <- .hweights(s)

  ranks <- 1:5
  for (rank in ranks) {
    cz <- cadzow(s, rank = rank, correct = TRUE)
    expect_true(.series.wsqdistance(cz, .F(s), w) < .series.wsqdistance((1 + delta) * cz, .F(s), w))
    expect_true(.series.wsqdistance(cz, .F(s), w) < .series.wsqdistance((1 - delta) * cz, .F(s), w))
  }
})

test_that("Cadzow for Complex SSA", {
  set.seed(1)
  N <- 100
  v <- rnorm(N) + 1i * rnorm(N)
  s <- ssa(v, kind = "cssa")
  eps <- sqrt(.Machine$double.eps)

  ranks <- 1:5
  for (rank in ranks) {
    cz <- cadzow(s, rank = rank)
    expect_true(high.rank.rate(cz, rank = rank, ssaobj = s) < eps)
  }
})
