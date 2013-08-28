library(testthat)
library(Rssa)

data(AustralianWine)


params <- list(same.length.mx = list(ssa.data = AustralianWine[1:176, c("Total", "Drywhite", "Fortified")],
                                     L = 84,
                                     len = 60),
               diff.lengths.mx = list(ssa.data = AustralianWine[, c("Total", "Drywhite", "Fortified")],
                                      L = 84,
                                      len = 60),
               same.length.list = list(ssa.data = list(AustralianWine[1:176, "Total"],
                                                       AustralianWine[1:176, "Drywhite"],
                                                       AustralianWine[1:176, "Fortified"]),
                                       L = 84,
                                       len = 60),
               diff.lengths.list = list(ssa.data = list(AustralianWine[1:176, "Total"],
                                                        AustralianWine[1:187, "Drywhite"],
                                                        AustralianWine[1:184, "Fortified"]),
                                        L = 84,
                                        len = 60),
               diff.lengths.list.NA = list(ssa.data = list(c(rep(NA, 1002), AustralianWine[21:176, "Total"]),
                                                           c(AustralianWine[1:187, "Drywhite"], rep(NA, 101)),
                                                           c(rep(NA, 503), as.vector(AustralianWine[1:184, "Fortified"]), rep(NA, 511))),
                                           L = 84,
                                           len = 60),
               dataframe = list(ssa.data = as.data.frame(AustralianWine[, c("Total", "Drywhite", "Fortified")]),
                                L = 84,
                                len = 60))

testcases <- lapply(params, function(param) {
ssa.data <- param$ssa.data; L <- param$L; len <- param$len
s <- ssa(ssa.data, L = L, kind = "mssa", neig = 15)
rec <- reconstruct(s,
                   groups = list(Trend = c(1, 6),
                                 Seasonality = c(2:5, 7:12)))

vrfore <- vforecast(s,
                    groups = list(1, 1:12),
                    direction = "row",
                    len = len, only.new = FALSE)

vcfore <- vforecast(s,
                    groups = list(1, 1:12),
                    direction = "column",
                    len = len, only.new = FALSE)

rrofore <- rforecast(s,
                     groups = list(1, 1:12),
                     direction = "row",
                     base = "original",
                     len = len, only.new = FALSE)

rrrfore <- rforecast(s,
                     groups = list(1, 1:12),
                     direction = "row",
                     base = "reconstructed",
                     len = len, only.new = FALSE)

rcofore <- rforecast(s,
                     groups = list(1, 1:12),
                     direction = "column",
                     base = "original",
                     len = len, only.new = FALSE)

rcrfore <- rforecast(s,
                     groups = list(1, 1:12),
                     direction = "column",
                     base = "reconstructed",
                     len = len, only.new = FALSE)

list(L = L, len = len,
     ssa.data = ssa.data,
     rec = rec,
     vrfore = vrfore, rrofore = rrofore, rrrfore = rrrfore,
     vcfore = vcfore, rcofore = rcofore, rcrfore = rcrfore)
})

save(testcases,
     file = system.file("extdata", "mssa.testdata.rda", package = "Rssa"),
     compress = "xz", compression_level = 9)
