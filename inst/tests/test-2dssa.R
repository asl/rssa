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
