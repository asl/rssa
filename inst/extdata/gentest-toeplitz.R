library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));

co2.td <- make.test.data(series = co2,
                         Ls = c(17, 234, 235, 300),
                         Ls.forecast = c(17, 100, 222, 234),
                         groups = as.list(1:5),
                         groups.forecast = list(1, 1:2, 3:5, 1:5),
                         len = 100,
                         kind = "toeplitz-ssa",
                         svd.method = "e",
                         neig = 10);
test.test.data(test.data = co2.td);

finite.rank.r5ex1 <- function(N) {
  tt <- 1:N;
  cos(2*pi*(1:N) / 7) + sin(2*pi*(1:N) / 17) * exp(tt / N * 1.5) + exp(-tt / N * 1.2);
}

fr50 <- finite.rank.r5ex1(50);
fr1k <- finite.rank.r5ex1(1000);
fr50k <- finite.rank.r5ex1(50000);

fr50.td <- make.test.data(series = fr50,
                          Ls = c(17, 25, 40),
                          Ls.forecast = c(17, 24, 25),
                          groups = as.list(1:5),
                          groups.forecast = list(1, 1:2, 3:5, 1:5, 5),
                          len = 100,
                          kind = "toeplitz-ssa",
                          svd.method = "e",
                          neig = 5);
test.test.data(test.data = fr50.td);

fr1k.td <- make.test.data(series = fr1k,
                          Ls = c(17, 493, 499, 500, 670),
                          Ls.forecast = c(17, 493, 499, 500),
                          groups = as.list(1:5),
                          groups.forecast = list(1, 1:2, 3:5, 1:5, 5),
                          len = 100,
                          kind = "toeplitz-ssa",
                          svd.method = "e",
                          neig = 5);
test.test.data(test.data = fr1k.td);

#fr50k.td <- make.test.data(series = fr50k,
#                           Ls = c(17, 493, 23800, 25000, 40000),
#                           Ls.forecast = c(17, 493, 23000, 23800, 25000),
#                           groups = as.list(1:5),
#                           groups.forecast = list(1, 1:2, 3:5, 1:5, 5),
#                           len = 100,
#                           kind = "toeplitz-ssa",
#                           svd.method = "p",
#                           neig = 5,
#                           tolerance = 1e-6,
#                           svd.methods = c("p", "n"));
#test.test.data(test.data = fr50k.td);

set.seed(1);
fr50.nz.td <- make.test.data(series = fr50 + rnorm(fr50),
                             name = "fr50.nz",
                             Ls = c(17, 25, 40),
                             Ls.forecast = c(17, 24, 25),
                             groups = as.list(1:10),
                             groups.forecast = list(1, 1:2, 3:5, 1:5, c(1, 3, 6, 10), 1:10),
                             len = 100,
                             kind = "toeplitz-ssa",
                             svd.method = "e",
                             neig = 15);
test.test.data(test.data = fr50.nz.td);

set.seed(1);
fr1k.nz.td <- make.test.data(series = fr1k + rnorm(fr1k),
                             name = "fr1k.nz",
                             Ls = c(17, 493, 499, 500, 670),
                             Ls.forecast = c(17, 493, 499, 500),
                             groups = as.list(1:10),
                             groups.forecast = list(1, 1:2, 3:5, 1:5, c(1, 3, 6, 10), 1:10),
                             len = 100,
                             kind = "toeplitz-ssa",
                             svd.method = "e",
                             neig = 15);
test.test.data(test.data = fr1k.nz.td);

#set.seed(1);
#fr50k.nz.td <- make.test.data(series = fr50k + rnorm(fr50k),
#                              name = "fr50k.nz",
#                              Ls = c(17, 493, 23800, 25000, 40000),
#                              Ls.forecast = c(17, 493, 23000, 23800, 25000),
#                              groups = as.list(1:10),
#                              groups.forecast = list(1, 1:2, 3:5, 1:5, c(1, 3, 6, 10), 1:10),
#                              len = 100,
#                              kind = "toeplitz-ssa",
#                              svd.method = "p",
#                              neig = 15,
#                              tolerance = 1e-6,
#                              svd.methods = c("p", "n"));
#test.test.data(test.data = fr50k.nz.td);

#save(co2.td, fr50.td, fr1k.td, fr50k.td, fr50.nz.td, fr1k.nz.td, fr50k.nz.td,
save(co2.td, fr50.td, fr1k.td, fr50.nz.td, fr1k.nz.td,
     file = system.file("extdata", "toeplitz.testdata.rda", package = "Rssa"),
     compress = "xz", compression_level = 9);
