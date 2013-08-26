library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("2dSSA");


test_that("simple 2s-ssa test", {
  load(system.file("extdata", "2dssa.testdata.rda", package = "Rssa"));

  for (svd.method in c("nutrlan", "propack")) {
    ss <- ssa(field, kind = "2d-ssa", L = L, neig = 20, svd.method = svd.method);
    cur.rec <- reconstruct(ss, groups = groups);

    expect_equal(cur.rec, expected.reconstruction,
                 label = sprintf("%s.2d.ssa reconstruction", svd.method));
  }
});

test_that("2d SSA works correctly with finite rank fields", {
# Artificial field for 2dSSA
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20))
for (svd.method in c("eigen", "svd", "nutrlan", "propack")) {
  # Decompose
  s <- ssa(mx, kind = "2d-ssa", neig = 5, svd.method = svd.method)
  # Reconstruct
  r <- reconstruct(s, groups = list(1:5))$F1
  expect_equal(r, mx, label = sprintf("svd.method = %s", svd.method))
}
})
