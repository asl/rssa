library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));

N <- c(110, 117);
L <- c(55, 53);
groups <- as.list(1:10);

set.seed(1);
field <- matrix(rnorm(prod(N)), N[1], N[2]);

ss <- ssa(field, kind = "2d-ssa", L = L, neig = 20);
expected.reconstruction <- reconstruct(ss, groups = groups);

save(L, groups, field, expected.reconstruction,
     file = system.file("extdata", "2dssa.testdata.rda", package = "Rssa"),
     compress = "xz", compression_level = 9);
