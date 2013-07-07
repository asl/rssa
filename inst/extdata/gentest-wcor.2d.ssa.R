library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))


set.seed(1)
mx <- outer(1:50, 1:50,
            function(i, j) sin(2*pi * i/17) * cos(2*pi * j/7) + exp(i/25 - j/20)) +
      rnorm(50^2, sd = 0.1)
# Decompose 'mx' with default parameters
s <- ssa(mx, kind = "2d-ssa", svd.method = "nutrlan")
w <- wcor(s, groups = 1:12)

save(w,
     file = system.file("extdata", "wcor.2dssa.testdata.rda", package = "Rssa"),
     compress = "xz", compression_level = 9);
