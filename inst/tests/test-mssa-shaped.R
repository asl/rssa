library(testthat)
library(Rssa)
context("MSSA")

new.hmat.striped.old <- function(F, L) {
  N <- sapply(F, length); K <- N - L + 1

  h <- lapply(seq_along(N),
              function(idx) new.hmat(F[[idx]], L = L))
  b <- c(0, cumsum(K))
  matmul <- function(v) {
    res <- numeric(L)
    for (idx in seq_along(h)) {
      res <- res + hmatmul(h[[idx]], v[(b[idx]+1):b[idx+1]], transposed = FALSE)
    }
    res
  }
  tmatmul <- function(v) unlist(lapply(h, hmatmul, v = v, transposed = TRUE))

  extmat(matmul, tmatmul, nrow = L, ncol = sum(K))
}

as.matrix.extmat <- function (x) {
  apply(diag(extmat.ncol(x)), 2, ematmul, emat = x)
}

test_that("new.hmat.striped and new.hmat.striped.old produce equal matrices", {
  set.seed(1)
  Ns <- list(10, 11, c(113, 113), 100, c(10, 14, 13, 11, 17))
  for (N in Ns) {
    Ls <- c(2, 4, 7, 9)
    F <- lapply(N, rnorm)

    for (L in Ls) {
      new <- .hmat.striped(ssa(F, L = L, kind = "mssa", force.decompose = FALSE))
      old <- new.hmat.striped.old(F, L)

      expect_equal(hbhrows(new), extmat.nrow(old))
      expect_equal(hbhcols(new), extmat.ncol(old))


      for (i in 1:42) {
        v <- rnorm(hbhcols(new))
        u <- rnorm(hbhrows(new))

        expect_equal(hbhmatmul(new, v, transposed = FALSE),
                     ematmul(old, v, transposed = FALSE))
        expect_equal(hbhmatmul(new, u, transposed = TRUE),
                     ematmul(old, u, transposed = TRUE))
      }
    }
  }
})

test_that("shaped MSSA works correct with gaps", {
  v1 <- 1:100
  v2 <- 1:200
  v1[60:63] <- NA
  v2[50:55] <- NA
  v1[1:3] <- NA
  v2[198:200] <- NA

  v <- list(v1, v2)
  r <- 2

  ss <- ssa(v, L = 10, kind = "mssa")
  expect_equal(reconstruct(ss, groups = list(1:r))$F1, v)

  w.exp <- structure(c(1, 0.0372400664289056, 0.0372400664289056, 1),
                     .Dim = c(2L, 2L),
                     .Dimnames = list(c("F1", "F2"), c("F1", "F2")),
                     class = "wcor.matrix")
  expect_equal(wcor(ss, 1:r), w.exp)
})

test_that("shaped MSSA works correct with gaps and uncovered points", {
  v1 <- 1:100
  v2 <- 1:200
  v1[60:63] <- NA
  v2[50:55] <- NA
  v1[3:5] <- NA
  v2[197:199] <- NA

  L <- 10
  v <- list(v1, v2)
  r <- 2
  v.res <- v
  v.res[[1]][1:2] <- NA
  v.res[[2]][200] <- NA

  ss <- ssa(v, L = L, kind = "mssa")
  expect_equal(reconstruct(ss, groups = list(1:r))$F1, v.res)

  w.exp <- structure(c(1, 0.0365310054106304, 0.0365310054106304, 1),
                     .Dim = c(2L, 2L),
                     .Dimnames = list(c("F1", "F2"), c("F1", "F2")),
                     class = "wcor.matrix")
  expect_equal(wcor(ss, 1:r), w.exp)
})
