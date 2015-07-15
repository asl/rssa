library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))
context("ProjectionSSA")

test_that("PSSA with no any projections is just SSA", {
  set.seed(1)
  Ns <- c(100, 200, 500, 1000)
  Ls <- c(70, 100, 200, 500)
  svd.methods <- c("propack", "nutrlan", "svd", "eigen")
  neig <- 10
  groups <- list(1, 2:3, 2:10)
  len <- 50

  for (i in seq_along(Ns)) {
    N <- Ns[i]
    L <- Ls[i]

    v <- rnorm(N)
    for (svd.method in svd.methods) {
      ss <- ssa(v, L = L, kind = "1d-ssa", svd.method = svd.method, neig = neig + 1)
      pss <- ssa(v, L = L, row.projector = 0, column.projector = 0,
                 svd.method = svd.method, neig = neig + 1)

      expect_equal(pss$sigma[seq_len(neig)], ss$sigma[seq_len(neig)],
                   label = sprintf("N = %d, L = %d, svd.method = %s",
                                   N, L, svd.method))

      expect_equal(reconstruct(pss, groups = groups)[], reconstruct(ss, groups = groups)[],
                   label = sprintf("N = %d, L = %d, svd.method = %s",
                                   N, L, svd.method))

      expect_equal(rforecast(pss, groups = groups, len = len, base = "original")[],
                   rforecast(ss, groups = groups, len = len, base = "original")[],
                   tolerance = 1e-6,
                   label = sprintf("N = %d, L = %d, svd.method = %s",
                                   N, L, svd.method))
      expect_equal(rforecast(pss, groups = groups, len = len, base = "reconstructed")[],
                   rforecast(ss, groups = groups, len = len, base = "reconstructed")[],
                   tolerance = 1e-6,
                   label = sprintf("N = %d, L = %d, svd.method = %s",
                                   N, L, svd.method))

      expect_equal(vforecast(pss, groups = groups, len = len)[],
                   vforecast(ss, groups = groups, len = len)[],
                   tolerance = 1e-6,
                   label = sprintf("N = %d, L = %d, svd.method = %s",
                                   N, L, svd.method))
    }
  }
})

test_that("All svd.methods provides the same decomposition", {
  set.seed(1)
  Ns <- c(100, 200, 500)
  Ls <- c(70, 100, 200)
  svd.methods <- c("propack", "nutrlan", "eigen")

  row.projectors <- 0:5
  column.projectors <- 0:5

  neig <- 20
  groups <- as.list(1:neig)

  for (i in seq_along(Ns)) {
    N <- Ns[i]
    L <- Ls[i]
    v <- rnorm(N)

    for (row.projector in row.projectors) {
      for (column.projector in column.projectors) {
        pss.svd <- ssa(v, L = L,
                       row.projector = row.projector, column.projector = column.projector,
                       svd.method = "svd", neig = neig + 1)

        for (svd.method in svd.methods) {
          pss <- ssa(v, L = L,
                     row.projector = row.projector, column.projector = column.projector,
                     svd.method = svd.method, neig = neig + 1)

          expect_equal(pss$sigma[seq_len(neig)], pss.svd$sigma[(seq_len(neig))],
                       label = sprintf("N = %d, L = %d, cp = %d, rp = %d, svd.method = %s",
                                       N, L, column.projector, row.projector, svd.method))

          expect_equal(reconstruct(pss, groups = groups)[], reconstruct(pss.svd, groups = groups)[],
                       label = sprintf("N = %d, L = %d, cp = %d, rp = %d, svd.method = %s",
                                       N, L, column.projector, row.projector, svd.method))
        }
      }
    }
  }
})

test_that("PSSA reconstruct and predict finite rank series exactly", {
  svd.methods <- c("propack", "nutrlan", "svd", "eigen")
  N <- 500
  len <- 50
  tt <- seq_len(N + len) / sqrt(N + len)
  vvs <- list(tt^4 + 2 * tt - 13,
              sin(2 * pi * tt/13),
              cos(tt) * exp(tt/N),
              sin(2 * pi * tt/11) + tt^2,
              tt^3 + 1)
  pranks <- c(5, 0, 0, 3, 4)
  npranks <- c(0, 2, 2, 2, 0)
  ranks <- npranks + pranks

  for (vvi in seq_along(vvs)) {
    vv <- vvs[[vvi]]
    nprank <- npranks[vvi]
    prank <- pranks[vvi]

    v <- vv[seq_len(N)]

    rc.projectors <- expand.grid(row.projector = 0:pranks[vvi], column.projector = 0:pranks[vvi])

    for (rci in seq_len(nrow(rc.projectors))) {
      row.projector <- rc.projectors$row.projector[rci]
      column.projector <- rc.projectors$column.projector[rci]
      rank <- max(prank, column.projector + row.projector) + nprank
      neig <- rank - (column.projector + row.projector) + 1

      for (svd.method in svd.methods) {
        pss <- ssa(v, row.projector = row.projector, column.projector = column.projector,
                   svd.method = svd.method, neig = neig)

        expect_equal(reconstruct(pss, groups = list(all = 1:rank))$all, v,
                     tolerance = 1e-6)

        expect_equal(rforecast(pss,
                               groups = list(all = 1:rank),
                               len = 50,
                               base = "original",
                               only.new = FALSE),
                     vv,
                     tolerance = 1e-6)
        expect_equal(rforecast(pss,
                               groups = list(all = 1:rank),
                               len = 50,
                               base = "reconstructed",
                               only.new = FALSE),
                     vv,
                     tolerance = 1e-6)

        expect_equal(vforecast(pss,
                               groups = list(all = 1:rank),
                               len = 50,
                               only.new = FALSE),
                     vv,
                     tolerance = 1e-6)
      }
    }
  }
})

test_that("PSSA (double centering) reconstruct test", {
  env <- new.env()
  load(system.file("extdata", "pssa.testdata.rda", package = "Rssa"), envir = env)
  #names <- c("co2.td", "fr50.td", "fr1k.td", "fr50k.td", "fr50.nz.td", "fr1k.nz.td", "fr50k.nz.td")
  names <- c("co2.td", "fr50.td", "fr1k.td", "fr50.nz.td", "fr1k.nz.td")
  for (name in names) {
    test.test.data(what = "reconstruct",
                   test.data = env[[name]],
                   kind = "1d-ssa",
                   row.projector = "centering",
                   column.projector = "centering")
  }
})

test_that("PSSA (double centerging) forecast test", {
  env <- new.env()
  load(system.file("extdata", "pssa.testdata.rda", package = "Rssa"), envir = env)
  #names <- c("co2.td", "fr50.td", "fr1k.td", "fr50k.td", "fr50.nz.td", "fr1k.nz.td", "fr50k.nz.td")
  names <- c("co2.td", "fr50.td", "fr1k.td", "fr50.nz.td", "fr1k.nz.td")
  for (name in names) {
    test.test.data(what = c("rforecast", "vforecast"),
                   test.data = env[[name]],
                   kind = "1d-ssa",
                   row.projector = "centering",
                   column.projector = "centering")
  }
})

test_that("PSSA backward predict finite rank series exactly", {
  svd.methods <- c("svd")
  N <- 500
  len <- 50
  tt <- seq_len(N + len)
  vvs <- list(tt^4 + 2 * tt - 13,
              sin(2 * pi * tt/13),
              cos(tt) * exp(tt/N),
              sin(2 * pi * tt/11) + tt^2,
              tt^3 + 1)
  vvs <- lapply(vvs, rev)
  pranks <- c(5, 0, 0, 3, 4)
  npranks <- c(0, 2, 2, 2, 0)
  ranks <- npranks + pranks

  for (vvi in seq_along(vvs)) {
    vv <- vvs[[vvi]]
    nprank <- npranks[vvi]
    prank <- pranks[vvi]

    v <- vv[len + seq_len(N)]

    rc.projectors <- expand.grid(row.projector = 0:pranks[vvi], column.projector = 0:pranks[vvi])

    for (rci in seq_len(nrow(rc.projectors))) {
      row.projector <- rc.projectors$row.projector[rci]
      column.projector <- rc.projectors$column.projector[rci]
      rank <- max(prank, column.projector + row.projector) + nprank
      neig <- rank - (column.projector + row.projector) + 1

      for (svd.method in svd.methods) {
        pss <- ssa(v, row.projector = row.projector, column.projector = column.projector,
                   svd.method = svd.method, neig = neig)

        expect_equal(reconstruct(pss, groups = list(all = 1:rank))$all, v,
                     tolerance = 1e-6)

        expect_equal(rforecast(pss,
                               groups = list(all = 1:rank),
                               reverse = TRUE,
                               len = len,
                               base = "original",
                               only.new = FALSE),
                     vv,
                     tolerance = 1e-6)
        expect_equal(rforecast(pss,
                               groups = list(all = 1:rank),
                               reverse = TRUE,
                               len = len,
                               base = "reconstructed",
                               only.new = FALSE),
                     vv,
                     tolerance = 1e-6)

        # expect_equal(vforecast(pss,
        #                        groups = list(all = 1:rank),
        #                        len = 50,
        #                        only.new = FALSE),
        #              vv,
        #              tolerance = 1e-6)
      }
    }
  }
})
