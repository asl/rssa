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
