library(Rssa)

cosine <- function(N, A, T, phase = 0) {
  k <- 0 : (N-1)
  A * cos(2 * pi * k / T + phase)
}

MSE <- function(x, y) {
  mean((x - y) ^ 2)
}

ssa.wrap <- function(xs, L, dim, svd.method = "eigen") {
  res <- lapply(xs, function(x) {
      ss <- ssa(x, L = L, kind = "1d-ssa", svd.method = svd.method, neig = dim + 2)
      reconstruct(ss, groups = list(1:dim))[[1]]
    })

  names(res) <- paste("x", seq_along(res), sep = "")

  res
}

cssa.wrap <- function(xs, L, dim, svd.method = "svd") {
  svd.method <- "svd" # TODO fix SVD methods in package

  stopifnot(length(xs) == 2)
  ss <- ssa(xs[[1]] + 1i * xs[[2]], L = L, kind = "cssa", svd.method = svd.method)
  rec <- reconstruct(ss, groups = list(1:dim))[[1]]

  list(x1 = Re(rec), x2 = Im(rec))
}

mssa.wrap <- function(xs, L, dim, svd.method = "auto") {
  ss <- ssa(xs, L = L, kind = "mssa", svd.method = svd.method, neig = dim + 2)
  res <- reconstruct(ss, groups = list(1:dim))[[1]]

  names(res) <- paste("x", seq_along(res), sep = "")

  res
}

ssa.wrap.forecast <- function(xs, L, dim, len,
                              ftype = c("reccurent", "vector"),
                              svd.method = "eigen",
                              ...) {
  ftype <- match.arg(ftype)

  forecast <- switch(ftype,
                     reccurent = rforecast,
                     vector = vforecast)

  res <- lapply(xs, function(x) {
      ss <- ssa(x, L = L, kind = "1d-ssa", svd.method = svd.method, neig = dim + 2)
      forecast(ss, groups = list(1:dim), len = len, drop = FALSE, only.new = TRUE)[[1]]
    })

  names(res) <- paste("x", seq_along(res), sep = "")

  res
}

cssa.wrap.forecast <- function(xs, L, dim, len,
                               ftype = c("reccurent", "vector"),
                               svd.method = "svd",
                               ...) {
  svd.method <- "svd" # TODO fix SVD methods in package
  stopifnot(length(xs) == 2)

  ftype <- match.arg(ftype)

  forecast <- switch(ftype,
                     reccurent = rforecast,
                     vector = vforecast)

  ss <- ssa(xs[[1]] + 1i * xs[[2]], L = L, kind = "cssa", svd.method = svd.method)
  fore <- forecast(ss, groups = list(1:dim), len = len, drop = FALSE, only.new = TRUE)[[1]]

  list(x1 = Re(fore), x2 = Im(fore))
}

mssa.wrap.forecast <- function(xs, L, dim, len,
                               direction = c("row", "column"),
                               ftype = c("reccurent", "vector"),
                               svd.method = "auto") {
  direction <- match.arg(direction)
  ftype <- match.arg(ftype)

  forecast <- switch(ftype,
                     reccurent = rforecast,
                     vector = vforecast)

  ss <- ssa(xs, L = L, kind = "mssa", svd.method = svd.method, neig = dim + 2)
  res <- forecast(ss, groups = list(1:dim), len = len, drop = FALSE, only.new = TRUE, direction = direction)[[1]]

  names(res) <- paste("x", seq_along(res), sep = "")

  res
}

MSE.ssa <- function(xs, xns, L, dim,
                    type = c("reconstruction", "forecast"),
                    kind = c("ssa", "r-ssa", "v-ssa",
                             "mssa", "mssa-row", "mssa-column", "r-mssa-row", "r-mssa-column", "v-mssa-row", "v-mssa-column",
                             "cssa", "r-cssa", "v-cssa"),
                    len, mse.idx = NULL,
                    svd.method = "auto") {
  if (is.null(mse.idx)) mse.idx <- seq_along(xs)

  kind <- match.arg(kind)
  if (identical(kind, "mssa")) kind <- "r-mssa-row"
  if (identical(kind, "mssa-row")) kind <- "r-mssa-row"
  if (identical(kind, "mssa-column")) kind <- "r-mssa-column"
  if (identical(kind, "ssa")) kind <- "r-ssa"
  if (identical(kind, "cssa")) kind <- "r-cssa"
  type <- match.arg(type)

  if (identical(type, "reconstruction")) {
    wrap <- switch(kind,
                   `r-ssa` = , `v-ssa` = ssa.wrap,
                   `r-mssa-row` = , `r-mssa-column` =, `v-mssa-row` = , `v-mssa-column` = mssa.wrap,
                   `r-cssa` = , `v-cssa` = cssa.wrap)

    xs <- lapply(xs, function(x) x[seq_len(length(x) - len)])
    xns <- lapply(xns, function(x) x[seq_len(length(x) - len)])
    xrs <- wrap(xns, L = L, dim = dim, svd.method = svd.method)
  } else if (identical(type, "forecast")) {
    wrap.forecast <- switch(kind,
                            `r-ssa` = , `v-ssa` = ssa.wrap.forecast,
                            `r-mssa-row` = , `r-mssa-column` =, `v-mssa-row` = , `v-mssa-column` = mssa.wrap.forecast,
                            `r-cssa` = , `v-cssa` = cssa.wrap.forecast)

   ftype <- switch(kind,
                   `r-ssa` = , `r-cssa` = , `r-mssa-row` = , `r-mssa-column` = "reccurent",
                   `v-ssa` = , `v-cssa` = , `v-mssa-row` = , `v-mssa-column` = "vector")

    xs <- lapply(xs, function(x) x[(length(x) - len + 1) : length(x)])
    xns <- lapply(xns, function(x) x[seq_len(length(x) - len)])
    xrs <- wrap.forecast(xns, L = L, dim = dim, len = len,
                         direction = switch(kind,
                                            `r-mssa-row` = , `v-mssa-row` = "row",
                                            `r-mssa-column` = , `v-mssa-column` = "column"),
                         ftype = ftype,
                         svd.method = svd.method)
  }

  mean(sapply(mse.idx, function(i) MSE(xs[[i]], xrs[[i]])))
}

tune.forecast <- function(XS,
                          base = XS,
                          Ls, dims,
                          kind,
                          len,
                          W = N,
                          mse.idx = NULL,
                          svd.method = "auto") {
  grid <- expand.grid(L = Ls, dim = dims)

  N <- length(XS[[1]])
  for (i in 1:(N - W + 1)) {
    xs <- lapply(XS, function(x) x[i : (i + W - 1)])
    xns <- lapply(base, function(x) x[i : (i + W - 1)])

    grid[[paste("MSE", i, sep = "")]] <- sapply(seq_len(nrow(grid)), function(i) {
      MSE.ssa(xs = xs, xns = xns,
              L = grid$L[i], dim = grid$dim[i],
              type = "forecast", kind = kind, len = len, mse.idx = mse.idx, svd.method = svd.method)
    })
  }

  MSE <- rowMeans(as.matrix(grid[, -(1:2), drop = FALSE]))
  grid <- grid[, 1:2]
  grid$MSE <- MSE
  opt.i <- which.min(grid$MSE)

  list(L = grid$L[opt.i], dim = grid$dim[opt.i], grid = grid)
}

all.MSE <- function(examples, N, sigma, Ls, len,
                    type = c("reconstruction", "forecast"),
                    kinds,
                    mse.idx = NULL,
                    svd.method = svd.method,
                    eval.sd = TRUE) {
  type <- match.arg(type)

  if (missing(kinds)) {
    kinds <- switch(type,
                    reconstruction = c("ssa", "mssa", "cssa"),
                    forecast = c("r-ssa", "v-ssa",
                                 "r-mssa-row", "r-mssa-column", "v-mssa-row", "v-mssa-column",
                                 "r-cssa", "v-cssa"))
  }

  get.dim <- function(example, kind) {
    name <- if (grepl("mssa", kind)) {
        "mssa"
      } else if (grepl("cssa", kind)) {
        "cssa"
      } else {
        "ssa"
      }

    example[[paste(name, "dim", sep = ".")]]
  }

  lapply(examples, function(example) {
    sapply(Ls, function(L) {
      xs <- example$xs
      base <- example$base
      if (is.null(base)) base <- xs
      mses <- replicate(N, {
          xns <- lapply(xs, function(x) x + rnorm(length(x), sd = sigma))


          sapply(kinds,
                 function(kind) MSE.ssa(if (identical(base, "original")) xns else base,
                                        xns,
                                        L = L,
                                        dim = get.dim(example, kind),
                                        kind = kind,
                                        len = len,
                                        type = type, mse.idx = mse.idx, svd.method = svd.method))
        })

      if (!is.matrix(mses)) mses <- matrix(mses, nrow = 1, ncol = length(mses))
      means <- rowMeans(mses)
      names(means) <- kinds
      if (eval.sd) {
        sds <- apply(mses, 1, sd) / sqrt(N)
        names(sds) <- paste(kinds, "sd", sep = ".")

        c(means, sds)
      } else {
        means
      }
    })
  })
}
