library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("1dSSA");

test_that("1dSSA reconstruct test", {
  env <- new.env();
  load(system.file("extdata", "1dssa.testdata.rda", package = "Rssa"), envir = env);
  #names <- c("co2.td", "fr50.td", "fr1k.td", "fr50k.td", "fr50.nz.td", "fr1k.nz.td", "fr50k.nz.td");
  names <- c("co2.td", "fr50.td", "fr1k.td", "fr50.nz.td", "fr1k.nz.td");
  for (name in names) {
    test.test.data(what = "reconstruct",
                   test.data = env[[name]]);
  }
});

test_that("Fast Lcov matrix' computation works correctly", {
  Ls <- c(3, 10, 17, 50, 371, 500, 1000);
  Ns <- c(1005, 1500, 2000, 5000);

  set.seed(1);
  for (N in Ns) {
    for (L in Ls) {
      F <- rcauchy(N);
      C.exact <- tcrossprod(hankel(F, L));
      C.fast <- Lcov.matrix(F, L = L);
      expect_equal(C.fast, C.exact,
                   info = sprintf("L = %d, N = %d", L, N));
    }
  }
});
