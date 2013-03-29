library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("FFT plan");

test_that("is.fft.plan-check works correctly", {
  expect_true(is.fft.plan(fft.plan.1d(42)));

  expect_false(is.fft.plan(42));
  expect_false(is.fft.plan(new.hmat(1:42)));

  regexp <- "pointer provided is not a fft plan";

  no.fft.plan <- 1;
  expect_error(new.hmat(1:42, L = 10, no.fft.plan), regexp = regexp);
  expect_error(Lcov.matrix(1:42, L = 10, no.fft.plan), regexp = regexp);
  expect_error(.Call("hankelize_one_fft", 1:10, 1:10, no.fft.plan), regexp = regexp);
  expect_error(.hankelize.multi.hankel(as.matrix(1:10), as.matrix(1:10), no.fft.plan), regexp = regexp);

  expect_error(new.tmat(1:42, L = 10, no.fft.plan), regexp = regexp);

  expect_error(dim.fft.plan(no.fft.plan), regexp = regexp);
});

test_that("dim.fft.plan works correctly", {
  expect_equal(dim.fft.plan(fft.plan.1d(42)), 42);
});
