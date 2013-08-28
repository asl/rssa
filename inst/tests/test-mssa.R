library(testthat)
library(Rssa)

context("Multichannel SSA AustralianWine multitest")

load(system.file("extdata", "mssa.testdata.rda", package = "Rssa"))

for (testcase.name in names(testcases)) {{{
testcase <- testcases[[testcase.name]]

ssa.data <- testcase$ssa.data
L <- testcase$L
len <- testcase$len

test_that(sprintf("Real data MSSA reconstruction test (testcase: %s)", testcase.name), {
for (svd.method in c("svd", "eigen", "nutrlan", "propack")) {
  L <- testcase$L
  s <- ssa(testcase$ssa.data, L = L, kind = "mssa",
           neig = 15,
           svd.method = svd.method)

  rec <- reconstruct(s,
                     groups = list(Trend = c(1, 6),
                                   Seasonality = c(2:5, 7:12)))
  expect_equal(rec, testcase$rec,
               label = sprintf("%s.mssa reconstruction", svd.method))
}
})

test_that(sprintf("Real data MSSA forecast test (testcase: %s)", testcase.name), {
s <- ssa(ssa.data, L = L, kind = "mssa",
         neig = 15,
         svd.method = "propack")
vrfore <- vforecast(s,
                    groups = list(1, 1:12),
                    direction = "row",
                    len = len, only.new = FALSE)
expect_equal(vrfore, testcase$vrfore,
             label = "mssa row-vector forecast")

vcfore <- vforecast(s,
                    groups = list(1, 1:12),
                    direction = "column",
                    len = len, only.new = FALSE)
expect_equal(vcfore, testcase$vcfore,
             label = "mssa column-vector forecast")

rrofore <- rforecast(s,
                     groups = list(1, 1:12),
                     direction = "row",
                     base = "original",
                     len = len, only.new = FALSE)
expect_equal(rrofore, testcase$rrofore,
             label = "mssa row-reccurent forecast (base = original)")

rrrfore <- rforecast(s,
                     groups = list(1, 1:12),
                     direction = "row",
                     base = "reconstructed",
                     len = len, only.new = FALSE)
expect_equal(rrrfore, testcase$rrrfore,
             label = "mssa row-reccurent forecast (base = reconstructed)")

rcofore <- rforecast(s,
                     groups = list(1, 1:12),
                     direction = "column",
                     base = "original",
                     len = len, only.new = FALSE)
expect_equal(rcofore, testcase$rcofore,
             label = "mssa column-reccurent forecast (base = original)")

rcrfore <- rforecast(s,
                     groups = list(1, 1:12),
                     direction = "column",
                     base = "reconstructed",
                     len = len, only.new = FALSE)
expect_equal(rcrfore, testcase$rcrfore,
             label = "mssa column-reccurent forecast (base = reconstructed)")
})
}}}
