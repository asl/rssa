library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("Signal parameters' estimation");

# test_that("parestimate.esprit works correctly for polynomial trends", {
# for (d in 0:3) {
#   N <- 40;
#   h <- t(hankel((1:N) ^ d, d + 1));
#   U <- svd(h)$u
#   par <- parestimate.esprit(U);
#
#   mu <- par$moduli * exp(pi * 2i / par$periods);
#   # Sad but true. Multiple eigenroots are evaluated with bad precision
#   expect_true(is_multisets_approx_equal(mu, rep(1, d + 1), tol = 0.001),
#               label = sprintf("Estimated characteristical roots for (1:N)^%d are correct", d));
#   }
# });

test_that("parestimate.esprit works correctly for two sines", {
  r <- 4;
  N <- 40;
  T1 <- 6;
  T2 <- 11;
  v <- sin(2 * pi * (1:N) / T1) + sin(2 * pi * (1:N) / T2);
  h <- t(hankel(v, r));
  U <- qr.Q(qr(h));
  par <- parestimate.esprit(U);
  mu <- par$moduli * exp(pi * 2i / par$periods);
  expectred.mu <- exp(pi * 2i / c(T1, -T1, T2, -T2));

  expect_true(is_multisets_approx_equal(mu, expectred.mu),
              label = sprintf("Est. ch. roots for sum of sines with periods %3.1f and %3.1f are correct", T1, T2));
});
