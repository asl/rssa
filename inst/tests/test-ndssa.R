library(testthat)
library(Rssa)
context("n-dimensional SSA")

mouter <- function(xs, FUN) {
  grid <- do.call(expand.grid, xs)
  args <- as.list(grid)
  names(args) <- names(xs)
  res <- do.call(FUN, args)
  dim(res) <- sapply(xs, length)

  res
}


test_that("nd-SSA works for 3d arrays of the finite rank", {
  sssin <- function(x, y, z) sin(x + y + z)
  r <- 2

  x <- mouter(list(1:10, 1:9, 1:12), sssin)
  ss <- ssa(x, kind = "nd-ssa")

  rec <- reconstruct(ss, groups = list(1:r))$F1

  expect_equal(rec, x)
})

test_that("nd-SSA works for 4d arrays of the finite rank", {
  ssssin <- function(x, y, z, t) sin(x + y + z + t)
  r <- 2

  x <- mouter(list(1:10, 1:9, 1:12, 1:5), ssssin)
  ss <- ssa(x)

  rec <- reconstruct(ss, groups = list(1:r))$F1

  expect_equal(rec, x)
})
