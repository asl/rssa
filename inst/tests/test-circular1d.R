library(testthat)
library(Rssa)

test_that("Circular Toeplitz SSA works correct for sines with finite circular rank", {
N <- 640
ii <- 3 * 2*pi * (1:N) / N

Ls <- c(3, 10, 17, 50, 371)

F <- sin(ii)
for (L in Ls) {
  for (svd.method in c("eigen", "svd", "nutrlan", "propack")) {
    if (identical(svd.method, "nutrlan") && L < 10)
      svd.method <- "propack"
    ss <- ssa(F, L = L, circular = TRUE, kind = "toeplitz-ssa",
              svd.method = svd.method,
              neig = 3)

    expect_true(sum(ss$lambda[-(1:2)]) < .Machine$double.eps^.25)
    expect_equal(reconstruct(ss, groups = list(1:2))$F1, F)
  }
}
})

test_that("Lcor computation works correctly for circular case", {
  Ls <- c(3, 10, 17, 50, 371, 500, 1000)
  Ns <- c(1005, 1500, 2000, 5000)

  set.seed(1)
  for (N in Ns) {
    for (L in Ls) {
      F <- rcauchy(N)
      C.exact <- convolve(F, F, conj = TRUE)[1:L] / N
      C <- Lcor(F, L = L, circular = TRUE)
      expect_equal(C, C.exact,
                   info = sprintf("L = %d, N = %d", L, N))
    }
  }
})
