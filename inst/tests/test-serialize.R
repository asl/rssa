library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("Serialization");

test_that("Serialization works correctly", {
  Ls <- c(17, 100, 222, 234);
  kinds <- c("1d-ssa", "toeplitz-ssa");
  svd.methods <- c("eigen", "nutrlan", "propack", "svd");
  groups <- list(1, 1:2, 3:5, 1:5);
  len <- 100;
  neig = 15;

  for (kind in kinds) for (svd.method in svd.methods) {
    if (identical(kind, "toeplitz-ssa") && identical(svd.method, "svd"))
      next;

    for (L in Ls) {
      suppressWarnings(ss <- ssa(co2, L = L, kind = kind, svd.method = svd.method, neig = neig));

      # Serialize ssa-object to raw vector
      rw <- serialize(ss, connection = NULL);

      # Unserialize ssa-object
      ss.uns <- unserialize(rw);

      expect_equal(ss.uns$U, ss$U);
      expect_equal(ss.uns$V, ss$V);
      expect_equal(ss.uns$lambda, ss$lambda);

      expect_equal(reconstruct(ss.uns, groups = groups), reconstruct(ss, groups = groups));
      expect_equal(rforecast(ss.uns, groups = groups, len = len, base = "original"),
                   rforecast(ss, groups = groups, len = len, base = "original"));

      expect_equal(rforecast(ss.uns, groups = groups, len = len, base = "reconstructed"),
                   rforecast(ss, groups = groups, len = len, base = "reconstructed"));

      expect_equal(vforecast(ss.uns, groups = groups, len = len),
                   vforecast(ss, groups = groups, len = len));
    }
  }
});
