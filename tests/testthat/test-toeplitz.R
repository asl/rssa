library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("Toeplitz SSA");

test_that("toeplitz SSA reconstruct test", {
  env <- new.env();
  load(system.file("extdata", "toeplitz.testdata.rda", package = "Rssa"), envir = env);
  #names <- c("co2.td", "fr50.td", "fr1k.td", "fr50k.td", "fr50.nz.td", "fr1k.nz.td", "fr50k.nz.td");
  names <- c("co2.td", "fr50.td", "fr1k.td", "fr50.nz.td", "fr1k.nz.td");
  for (name in names) {
    test.test.data(what = "reconstruct",
                   test.data = env[[name]]);
  }
});
