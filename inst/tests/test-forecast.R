library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("SSA forecast");

test_that("1dSSA forecast test", {
  env <- new.env();
  load(system.file("extdata", "1dssa.testdata.rda", package = "Rssa"), envir = env);
  #names <- c("co2.td", "fr50.td", "fr1k.td", "fr50k.td", "fr50.nz.td", "fr1k.nz.td", "fr50k.nz.td");
  names <- c("co2.td", "fr50.td", "fr1k.td", "fr50.nz.td", "fr1k.nz.td");
  for (name in names) {
    test.test.data(what = c("rforecast", "vforecast"),
                   test.data = env[[name]]);
  }
});

test_that("toeplitz SSA forecast test", {
  env <- new.env();
  load(system.file("extdata", "toeplitz.testdata.rda", package = "Rssa"), envir = env);
  #names <- c("co2.td", "fr50.td", "fr1k.td", "fr50k.td", "fr50.nz.td", "fr1k.nz.td", "fr50k.nz.td");
  names <- c("co2.td", "fr50.td", "fr1k.td", "fr50.nz.td", "fr1k.nz.td");
  for (name in names) {
    test.test.data(what = c("rforecast", "vforecast"),
                   test.data = env[[name]]);
  }
});
