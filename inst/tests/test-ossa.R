library(testthat)
library(Rssa)
context("OSSA")

test_that("I-OSSA separates 3 sines exactly", {
  N <- 150
  L <- 70

  omega1 <- 0.05
  omega2 <- 0.06
  omega3 <- 0.07

  F1.real <- 4*sin(2*pi*omega1*(1:N))
  F2.real <- 2*sin(2*pi*omega2*(1:N))
  F3.real <- sin(2*pi*omega3*(1:N))
  F <- F1.real + F2.real + F3.real
  ss <- ssa(F, L)
  ioss <- iossa(ss, nested.groups = list(1:2, 3:4, 5:6), maxIter = 200, tol = 1e-8, kappa = NULL)

  rec <- reconstruct(ioss, groups = ioss$iossa.groups)
  expect_equal(rec$F1, F1.real, tolerance = 1e-6)
  expect_equal(rec$F2, F2.real, tolerance = 1e-6)
  expect_equal(rec$F3, F3.real, tolerance = 1e-6)
})

test_that("I-OSSA and F-OSSA", {
  N <- 200
  L <- 100
  omega1 <- 0.07
  omega2 <- 0.06

  F1.real <- 2*sin(2*pi*omega1*(1:N))
  F2.real <- 2*sin(2*pi*omega2*(1:N))
  ss <- ssa(F1.real + F2.real, L, svd.method = "eigen", neig = 28)
  fss <- fossa(ss, nested.groups = list(c(1,2), c(3,4)), kappa = 100)
  ioss <- iossa(fss, nested.groups = list(c(1,2), c(3,4)), maxIter = 1000, kappa = 2, tol = 1e-8)

  rec <- reconstruct(ioss, groups = ioss$iossa.groups)
  expect_equal(rec$F1, F1.real, tolerance = 1e-6)
  expect_equal(rec$F2, F2.real, tolerance = 1e-6)

  wc <- owcor(ioss, groups = list(1:2, 3:4))
  expect_equivalent(wc[,], diag(2))

  expect_true(ioss$iossa.result$conv)
})

test_that("FOSSA", {
  N <- 150
  L <- 70
  omega1 <- 1/5
  omega2 <- 1/10

  F1.real <- 2*sin(2*pi*omega1*(1:N))
  F2.real <- 2*sin(2*pi*omega2*(1:N))
  v <- F1.real + F2.real
  ss <- ssa(v, L, svd.method = "eigen")
  fss <- fossa(ss,  nested.groups = list(1:2, 3:4), gamma = 100.5)
  wc <- wcor(fss, groups = list(1:2, 3:4))

  expect_equivalent(wc[,], diag(2))

  rec <- reconstruct(fss, groups = list(1:2, 3:4))
  expect_equal(rec$F1, F1.real, tolerance = 1e-6)
  expect_equal(rec$F2, F2.real, tolerance = 1e-6)
})

test_that ("OSSA + PSSA forecast is correct", {
  N <- 100
  len <- 20
  tt <- seq_len(N + len)
  F <- 0.01 * tt^2 + 10 * sin(2*pi * tt / 10)
  pss <- ssa(F[seq_len(N)], row.projector = "centering", column.projector = "centering")
  ios <- iossa(pss, nested.groups = list(c(1:2), c(3:5)))
  fos <- fossa(ios, nested.groups = ios$iossa.groups, gamma = 1000)

  rforec.ios <- rforecast(ios, groups = list(1:5), len = len, only.new = FALSE)
  vforec.ios <- vforecast(ios, groups = list(1:5), len = len, only.new = FALSE)
  expect_equal(rforec.ios, F)
  expect_equal(vforec.ios, F)

  rforec.fos <- rforecast(fos, groups = list(1:5), len = len, only.new = FALSE)
  vforec.fos <- vforecast(fos, groups = list(1:5), len = len, only.new = FALSE)
  expect_equal(rforec.fos, F)
  expect_equal(vforec.fos, F)
})
