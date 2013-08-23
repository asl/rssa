library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))

context("Shaped 2dSSA")

test_that("new.hbhmat returns matrix with proper dimension", {
  data(Mars)
  s <- ssa(Mars, mask = Mars != 0, wmask = circle(50), neig = 50, kind = "2d-ssa", force.decompose = FALSE)
  h <- .get.or.create.hbhmat(s)
  expect_equal(sum(s$wmask), hbhrows(h))
  expect_equal(sum(s$fmask), hbhcols(h))

  s <- ssa(matrix(1, 100, 100), wmask = circle(10), kind = "2d-ssa", force.decompose = FALSE)
  h <- .get.or.create.hbhmat(s)
  expect_equal(sum(s$wmask), hbhrows(h))
  expect_equal(prod(s$length - s$window + 1), hbhcols(h))
})

.sh2hbhmat <- function(x) {
  N <- x$length
  L <- x$window
  K <- N - L + 1
  wmask <- .get(x, "wmask")
  if (is.null(wmask))
    wmask <- matrix(TRUE, L[1], L[2])

  fmask <- .get(x, "fmask")
  if (is.null(fmask))
    fmask <- matrix(TRUE, K[1], K[2])

  F <- .get(x, "F")
  weights <- .get(x, "weights")
  if (!is.null(weights)) {
    mask <- weights > 0
    F[!mask] <- mean(F[mask]) # Improve FFT stability & remove NAs
  } else {
    weights <- tcrossprod(.hweights.default(N[1], L[1]),
                          .hweights.default(N[2], L[2]))
  }

  matmul <- function(v) fmatmul(F, v, wmask, fmask)
  tmatmul <- function(u) fmatmul(F, u, wmask, fmask, transposed = TRUE)

  extmat(matmul, tmatmul, sum(wmask), sum(fmask))
}

fmatmul <- function(field, X, wmask, fmask, transposed = FALSE) {
  if (transposed) {
    tmp <- wmask; wmask <- fmask; fmask <- tmp
  }

  stopifnot(length(X) == sum(fmask))

  x <- fmask
  x[fmask] <- X
  y <- convolve2.filter(field, x)
  Y <- y[wmask]

  Y
}

.hankelize.one.shaped2d.ssa <- function(x, U, V) {
  N <- x$length
  L <- x$window
  K <- N - L + 1
  wmask <- .get(x, "wmask")
  if (is.null(wmask))
    wmask <- matrix(TRUE, L[1], L[2])

  fmask <- .get(x, "fmask")
  if (is.null(fmask))
    fmask <- matrix(TRUE, K[1], K[2])

  weights <- .get(x, "weights")
  if (is.null(weights)) {
    weights <- tcrossprod(.hweights.default(N[1], L[1]),
                          .hweights.default(N[2], L[2]))
  }
  weights <- as.integer(round(weights))

  x <- wmask
  x[wmask] <- U

  y <- fmask
  y[fmask] <- V

  res <- convolve2.open(x, y)
  res <- res / weights
  res[weights == 0] <- NA

  res
}

test_that("Shaped SSA is common case of 2dSSA", {
  F0 <- matrix(c(0, 1, 0,
                 1, 1, 1),
               2, 3)
  F <- rbind(cbind(F0, F0), cbind(F0, F0))

  # Test simple matrix multiplication
  vtest <- as.vector(F0)
  ht_mul_v <- c(4, 2, 4, 3, 2, 3, 3, 2, 3, 4, 2, 4)
  utest <- c(0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0)
  h_mul_u <- c(3, 2, 1, 3, 1, 2)
  Frec <- tcrossprod(c(1, 3, 9, 27), c(1, 2, 4, 8, 16, 32))


  # Test 2D-SSA
  s <- ssa(F, kind = "2d-ssa", force.decompose = FALSE)
  h.CPP <- .get.or.create.hbhmat(s)

  # Test multiplication
  expect_equal(hbhmatmul(h.CPP, vtest, transposed = TRUE), ht_mul_v,
               label = "tmatmul, 2D-SSA")
  expect_equal(hbhmatmul(h.CPP, utest), h_mul_u,
               label = "matmul, 2D-SSA")

  # Test Hankelization
  Frec_res <- .hankelize.one.2d.ssa(s,
                                    as.vector(Frec[1:2, 1:3]),
                                    as.vector(Frec[1:3, 1:4]))
  dim(Frec_res) <- dim(Frec)
  expect_equal(Frec_res, Frec, label = "hankelize, 2D-SSA")


  # Test R-SH-SSA
  h.R <- .sh2hbhmat(s)

  # Test multiplication
  expect_equal(ematmul(h.R, vtest, transposed = TRUE), ht_mul_v,
               label = "tmatmul, SH-SSA")
  expect_equal(ematmul(h.R, utest), h_mul_u,
               label = "matmul, SH-SSA")

  # Test Hankelization
  Frec_res <- .hankelize.one.shaped2d.ssa(s,
                                          as.vector(Frec[1:2, 1:3]),
                                          as.vector(Frec[1:3, 1:4]))
  dim(Frec_res) <- dim(Frec)
  expect_equal(Frec_res, Frec, label = "hankelize, SH-SSA")


  # Test 2D-SH-SSA (hbhmat with rectangular wmask)
  s3 <- ssa(F, kind = "2d-ssa", wmask = matrix(TRUE, 2, 3), force.decompose = FALSE)
  h3 <- .get.or.create.hbhmat(s3)

  # Test multiplication
  expect_equal(hbhmatmul(h3, vtest, transposed = TRUE), ht_mul_v,
               label = "tmatmul, 2D-SSA")
  expect_equal(hbhmatmul(h3, utest), h_mul_u,
               label = "matmul, 2D-SH-SSA")

  # Test Hankelization
  Frec_res <- .hankelize.one.2d.ssa(s3,
                                    as.vector(Frec[1:2, 1:3]),
                                    as.vector(Frec[1:3, 1:4]))
  dim(Frec_res) <- dim(Frec)
  expect_equal(Frec_res, Frec,
               label = "hankelize, 2D-SH-SSA")
})


test_that("Shaped 2D-SSA test", {
  F0 <- matrix(c(0, 1, 0,
                 1, 1, 1),
               2, 3)
  F <- rbind(cbind(F0, F0), cbind(F0, F0))
  F <- rbind(cbind(F, F), cbind(F, F))

  mask <- matrix(TRUE, nrow(F), ncol(F))
  mask[nrow(F), 1:2] <- FALSE

  wmask <- matrix(c(FALSE, TRUE,  TRUE,  TRUE,
                    TRUE,  FALSE, TRUE,  FALSE,
                    FALSE, FALSE, FALSE, TRUE),
                  3, 4)

  utest <- rnorm(6)
  vtest <- rnorm(52)

  s <- ssa(F, kind = "2d-ssa", mask = mask, wmask = wmask, force.decompose = FALSE)
  h.CPP <- .get.or.create.hbhmat(s)
  h.R <- .sh2hbhmat(s)

  # Test Shaped matrices multiplication
  expect_equal(ematmul(h.R, utest, transposed = TRUE), hbhmatmul(h.CPP, utest, transposed = TRUE),
               label = "tmatmul, SH-SSA <-> 2D-SH-SSA",
               tolerance = 1e-5)
  expect_equal(ematmul(h.R, vtest), hbhmatmul(h.CPP, vtest),
               label = "matmul, SH-SSA <-> 2D-SH-SSA",
               tolerance = 1e-5)

  # Test hankelization
  rec.R <- .hankelize.one.shaped2d.ssa(s, utest, vtest)
  rec.CPP <- .hankelize.one.2d.ssa(s, utest, vtest)
  dim(rec.CPP) <- dim(F)
  expect_equal(is.na(rec.R), is.na(rec.CPP))
  mask <- !is.na(rec.R)
  expect_equal(rec.R[mask], rec.CPP[mask])
})
