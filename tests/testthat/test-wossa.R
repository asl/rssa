library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))
context("Weighted Oblique SSA and Cadzow")

test_that("Weighted Oblique Cadzow limit is a series of finite rank", {
  L <- 228
  row.oblique <- c(rep(1, 22), rep(0.12, 197), rep(1, 22))
  column.oblique <- rep(1, L)
  s <- ssa(co2, L = L, row.oblique = row.oblique, column.oblique = column.oblique)
  eps <- sqrt(.Machine$double.eps)

  ranks <- 1:5
  for (rank in ranks) {
    cz <- cadzow(s, rank = rank)
    expect_true(high.rank.rate(cz, rank = rank, ssaobj = s) < eps)
  }
})

.series.wsqdistance <- function(F1, F2, weights = rep(1, length(F1))) {
  mask <- weights > 0

  weights <- weights[mask]
  F1 <- as.vector(unlist(F1))[mask]
  F2 <- as.vector(unlist(F2))[mask]

  sum(weights * (F1-F2)^2)
}

test_that("Weighted Oblique Cadzow gives closer answer than Basic Cadzow", {
  L <- 228
  row.oblique <- c(rep(1, 22), rep(0.12, 197), rep(1, 22))
  column.oblique <- rep(1, L)
  s_wo <- ssa(co2, L = L, row.oblique = row.oblique, column.oblique = column.oblique)

  s_b <- ssa(co2, L = L)

  ranks <- 1:5
  for (rank in ranks) {
    cz_wo <- cadzow(s_wo, rank = rank)
    cz_b <- cadzow(s_b, rank = rank)
    
    expect_true(.series.wsqdistance(cz_wo, .F(s_wo)) <= 
      .series.wsqdistance(cz_b, .F(s_wo)))
  }
})
