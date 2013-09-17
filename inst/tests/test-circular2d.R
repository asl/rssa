library(testthat)
library(Rssa)

context("Circular 2d SSA")
test_that("Circular SSA works correct for sines with finite circular rank", {
N <- 32
M <- 36
ii <- 3 * 2*pi * (1:N) / N
jj <- 4 * 2*pi * (1:M) / M

mx <- outer(ii, jj, function(i, j) sin(i) + sin(j))

ss <- ssa(mx, circular = TRUE, kind = "2d-ssa")

expect_true(sum(ss$lambda[-(1:4)]) < .Machine$double.eps^.25)
expect_equal(reconstruct(ss, groups = list(1:4))$F1, mx)
})

test_that("Half-circular SSA works correct for sines with finite circular rank", {
N <- 32
M <- 36
ii <- 3 * 2*pi * (1:N) / N
jj <- (1:M) / M

mx <- outer(ii, jj, function(i, j) sin(i) * exp(j))

ss <- ssa(mx, circular = c(TRUE, FALSE), kind = "2d-ssa")

expect_true(sum(ss$lambda[-(1:2)]) < .Machine$double.eps^.25)
expect_equal(reconstruct(ss, groups = list(1:2))$F1, mx)
})
