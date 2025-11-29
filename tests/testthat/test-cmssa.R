library(testthat)
library(Rssa)
context("CMSSA")

test.cases <- c("linear", "periodic")

for (case in test.cases) {
test_that(sprintf("Complex MSSA reconstruction is correct, %s case", case), {
  N <- 100; P <- 5; coefs <- seq(from = 1, to = P)
  v <- switch(case,
              linear = (1 + 1i) * seq(N) %o% coefs,
              periodic = exp(1i * seq(N) / 12) %o% coefs)
  s <- ssa(v, kind = "cmssa")
  eps <- sqrt(.Machine$double.eps)

  rank <- switch(case,
                  linear = 2,
                  periodic = 1)

  rec <- reconstruct(s, groups = list(1:rank))$F1
  expect_true(sqrt(mean(abs(rec - v)^2)) < eps)
})

test_that(sprintf("Complex MSSA recurrent column forecast is correct, %s case", case), {
  N <- 110; P <- 5; coefs <- seq(from = 1, to = P); len <- 10
  v <- switch(case,
              linear = (1 + 1i) * seq(N) %o% coefs,
              periodic = exp(1i * seq(N) / 12) %o% coefs)
  s <- ssa(v[1:(N - len)], kind = "cmssa")
  eps <- sqrt(.Machine$double.eps)

  rank <- switch(case,
                 linear = 2,
                 periodic = 1)

  pred <- rforecast(s, groups = list(1:rank), direction = "column", len = len)
  expect_true(sqrt(mean(abs(pred - v[(N - len + 1):N])^2)) < eps)
})
}

test_that("Built-in SVD and Primme SVD yield same reconstructions in Complex MSSA", {
  N <- 100; P <- 5;
  v <- matrix(rnorm(N * P) + 1i * rnorm(N * P), ncol = P)
  s_default <- ssa(v, kind = "cmssa")
  s_primme <- ssa(v, kind = "cmssa", svd = "primme")
  eps <- sqrt(.Machine$double.eps)

  ranks <- 1:5

  for (rank in ranks) {
    rec_default <- reconstruct(s_default, groups = list(1:rank))$F1
    rec_primme <- reconstruct(s_primme, groups = list(1:rank))$F1
    expect_true(sqrt(mean(abs(rec_default - rec_primme)^2)) < eps)
  }

})